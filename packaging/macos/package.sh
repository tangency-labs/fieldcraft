#!/usr/bin/env bash
# Build, sign and (with credentials) notarize FieldCraft for macOS:
#
#   $DIST/fieldcraft-<version>-macos-<arch>.dmg        FieldCraft.app on a drag-to-Applications DMG
#   $DIST/fieldcraft-cli-<version>-macos-<arch>.zip    fieldcraft-cli with the licence files
#
# Usage: packaging/macos/package.sh [--arch universal|aarch64|x86_64] [--skip-build]
#
# Signing (env):
#   MACOS_SIGN_IDENTITY   codesign identity (name or SHA-1). Default "-": ad-hoc, which Gatekeeper
#                         blocks on other Macs until the user allows it. In CI, import-cert.sh sets
#                         it from the APPLE_CERTIFICATE secret.
#   MACOS_KEYCHAIN        keychain holding the identity (optional)
# Notarization (env; all three, and a real identity): APPLE_ID, APPLE_PASSWORD (an app-specific
# password), APPLE_TEAM_ID.
#
# Writes `signing=notarized|signed|adhoc` to $GITHUB_OUTPUT when run in GitHub Actions.
# Adapted from SolidCraft's and VectorCraft's packaging/macos/package.sh (MIT OR Apache-2.0; see
# ATTRIBUTION.md).
set -euo pipefail
# shellcheck source=../env.sh
. "$(dirname "${BASH_SOURCE[0]}")/../env.sh"
HERE="$ROOT/packaging/macos"
ICONS="$ROOT/assets/app-icon"

ARCH=universal
SKIP_BUILD=0
while [ $# -gt 0 ]; do
  case "$1" in
    --arch) ARCH="$2"; shift 2 ;;
    --skip-build) SKIP_BUILD=1; shift ;;
    -h | --help) sed -n '2,17p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
case "$ARCH" in
  universal) TARGETS=(aarch64-apple-darwin x86_64-apple-darwin); WANT=(arm64 x86_64) ;;
  aarch64) TARGETS=(aarch64-apple-darwin); WANT=(arm64) ;;
  x86_64) TARGETS=(x86_64-apple-darwin); WANT=(x86_64) ;;
  *) echo "unknown --arch $ARCH" >&2; exit 2 ;;
esac

# The oldest macOS the binaries run on; also the Info.plist's LSMinimumSystemVersion.
export MACOSX_DEPLOYMENT_TARGET=11.0
IDENTITY="${MACOS_SIGN_IDENTITY:--}"
SHORT_VERSION="${VERSION%%-*}"
WORK="$CARGO_TARGET_DIR/macos-package"
APP="$WORK/FieldCraft.app"
DMG="$DIST/fieldcraft-$VERSION-macos-$ARCH.dmg"
CLI_ZIP="$DIST/fieldcraft-cli-$VERSION-macos-$ARCH.zip"

NOTARIZE=0
SIGNING=adhoc
if [ "$IDENTITY" = "-" ]; then
  warn "macOS: ad-hoc signing (no MACOS_SIGN_IDENTITY); the app is not notarized"
elif [ -n "${APPLE_ID:-}" ] && [ -n "${APPLE_PASSWORD:-}" ] && [ -n "${APPLE_TEAM_ID:-}" ]; then
  NOTARIZE=1
  SIGNING=notarized
else
  SIGNING=signed
  warn "macOS: APPLE_ID / APPLE_PASSWORD / APPLE_TEAM_ID incomplete; signed but not notarized"
fi

echo "==> FieldCraft $VERSION for macOS ($ARCH), identity: $IDENTITY, notarize: $NOTARIZE"

# ---- build -------------------------------------------------------------------------------------
if [ "$SKIP_BUILD" = 0 ]; then
  args=()
  for t in "${TARGETS[@]}"; do args+=(--target "$t"); done
  (cd "$ROOT" && cargo build --release --locked -p fieldcraft -p fieldcraft-cli "${args[@]}")
fi

rm -rf "$WORK"
mkdir -p "$WORK/bin"
for bin in fieldcraft fieldcraft-cli; do
  inputs=()
  for t in "${TARGETS[@]}"; do inputs+=("$CARGO_TARGET_DIR/$t/release/$bin"); done
  lipo -create -output "$WORK/bin/$bin" "${inputs[@]}"
  lipo -info "$WORK/bin/$bin"
done

# ---- signing helpers ---------------------------------------------------------------------------
sign() {
  # Hardened runtime and a secure timestamp for a real identity; ad-hoc can't be timestamped.
  local ts=(--timestamp)
  [ "$IDENTITY" = "-" ] && ts=(--timestamp=none)
  local kc=()
  [ -n "${MACOS_KEYCHAIN:-}" ] && kc=(--keychain "$MACOS_KEYCHAIN")
  codesign --force --sign "$IDENTITY" ${kc[@]+"${kc[@]}"} "${ts[@]}" "$@"
}

notarize() {
  local file="$1" out id status
  echo "==> notarizing $(basename "$file") (this can take a few minutes)"
  out="$(xcrun notarytool submit "$file" --apple-id "$APPLE_ID" --password "$APPLE_PASSWORD" \
    --team-id "$APPLE_TEAM_ID" --wait --timeout 1h --output-format json)" || true
  echo "$out"
  id="$(printf '%s' "$out" | plutil -extract id raw -o - - 2>/dev/null || true)"
  status="$(printf '%s' "$out" | plutil -extract status raw -o - - 2>/dev/null || true)"
  if [ "$status" != "Accepted" ]; then
    if [ -n "$id" ]; then
      xcrun notarytool log "$id" --apple-id "$APPLE_ID" --password "$APPLE_PASSWORD" --team-id "$APPLE_TEAM_ID" || true
    fi
    echo "error: notarization of $(basename "$file") failed (status: ${status:-unknown})" >&2
    exit 1
  fi
}

# hdiutil attaches a device while it works, which fails now and then on CI runners ("Resource
# busy"). Usage: hdiutil_retry <hdiutil args…>; five tries with growing pauses.
hdiutil_retry() {
  local what="$1" attempt
  for attempt in 1 2 3 4 5; do
    if hdiutil "$@"; then return 0; fi
    [ "$attempt" = 5 ] && { echo "hdiutil $what failed 5 times" >&2; return 1; }
    echo "hdiutil $what failed (attempt $attempt of 5); retrying in $((attempt * 10))s" >&2
    sleep $((attempt * 10))
  done
}

# ---- FieldCraft.app ----------------------------------------------------------------------------
echo "==> assembling $APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/Licenses"
# The executable and the icons carry the display name (CFBundleExecutable, CFBundleIconFile,
# CFBundleIconName). Assets.car holds the layered icon for macOS 26 and later, the .icns the
# icon for earlier versions (assets/app-icon/README.md).
cp "$WORK/bin/fieldcraft" "$APP/Contents/MacOS/FieldCraft"
cp "$WORK/bin/fieldcraft-cli" "$APP/Contents/MacOS/fieldcraft-cli"
cp "$ICONS/icon.icns" "$APP/Contents/Resources/FieldCraft.icns"
cp "$ICONS/Assets.car" "$APP/Contents/Resources/Assets.car"
copy_docs "$APP/Contents/Resources/Licenses"
sed -e "s/@APP_ID@/$APP_ID/g" -e "s/@VERSION@/$VERSION/g" -e "s/@SHORT_VERSION@/$SHORT_VERSION/g" \
  -e "s/@BUILD_SHA@/$FIELDCRAFT_BUILD_SHA/g" -e "s/@MIN_MACOS@/$MACOSX_DEPLOYMENT_TARGET/g" "$HERE/Info.plist.in" >"$APP/Contents/Info.plist"
plutil -lint "$APP/Contents/Info.plist"
printf 'APPL????' >"$APP/Contents/PkgInfo"

# Sign inside-out: both executables, then the bundle (no --deep on the bundle's own signature).
sign --options runtime --identifier "$APP_ID.cli" "$APP/Contents/MacOS/fieldcraft-cli"
sign --options runtime --entitlements "$HERE/entitlements.plist" "$APP/Contents/MacOS/FieldCraft"
sign --options runtime --entitlements "$HERE/entitlements.plist" "$APP"
codesign --verify --strict --deep --verbose=2 "$APP"
codesign -dv "$APP"

if [ "$NOTARIZE" = 1 ]; then
  ditto -c -k --keepParent "$APP" "$WORK/FieldCraft-notarize.zip"
  notarize "$WORK/FieldCraft-notarize.zip"
  xcrun stapler staple "$APP"
  xcrun stapler validate "$APP"
  spctl --assess --type execute -vvv "$APP"
fi

# ---- DMG ---------------------------------------------------------------------------------------
echo "==> building $DMG"
STAGE="$WORK/dmg"
mkdir -p "$STAGE"
ditto "$APP" "$STAGE/FieldCraft.app"
ln -s /Applications "$STAGE/Applications"
cp "$HERE/README.txt" "$STAGE/README.txt"
rm -f "$DMG"
# `create -srcfolder` writes the files as they are; `makehybrid -hfs` gives every file a
# com.apple.FinderInfo, which `codesign --verify --strict` rejects (VectorCraft found this).
hdiutil_retry create -srcfolder "$STAGE" -volname "FieldCraft" -fs HFS+ -format UDZO \
  -imagekey zlib-level=9 -ov "$DMG"
sign "$DMG"
codesign --verify --strict --verbose=2 "$DMG"

# Check the image as a user gets it: mount it and look inside.
MNT="$WORK/dmg-check"
rm -rf "$MNT"
mkdir -p "$MNT"
hdiutil_retry attach -nobrowse -readonly -mountpoint "$MNT" "$DMG" >/dev/null
# Detach on the way out, also when a check below fails and ends the script.
trap 'hdiutil detach "$MNT" >/dev/null 2>&1 || true' EXIT
M_APP="$MNT/FieldCraft.app"
M_PLIST="$M_APP/Contents/Info.plist"
test -L "$MNT/Applications"
test "$(readlink "$MNT/Applications")" = /Applications
cmp "$HERE/README.txt" "$MNT/README.txt"
plutil -lint "$M_PLIST"
for key in CFBundleName CFBundleDisplayName CFBundleExecutable CFBundleIconFile CFBundleIconName; do
  test "$(/usr/libexec/PlistBuddy -c "Print :$key" "$M_PLIST")" = FieldCraft
done
test "$(/usr/libexec/PlistBuddy -c 'Print :CFBundleIdentifier' "$M_PLIST")" = "$APP_ID"
test "$(/usr/libexec/PlistBuddy -c 'Print :FieldCraftVersion' "$M_PLIST")" = "$VERSION"
cmp "$ICONS/icon.icns" "$M_APP/Contents/Resources/FieldCraft.icns"
cmp "$ICONS/Assets.car" "$M_APP/Contents/Resources/Assets.car"
for f in "${DOC_FILES[@]}"; do
  test -s "$M_APP/Contents/Resources/Licenses/$f"
done
for bin in FieldCraft fieldcraft-cli; do
  test -x "$M_APP/Contents/MacOS/$bin"
  for arch in "${WANT[@]}"; do lipo "$M_APP/Contents/MacOS/$bin" -verify_arch "$arch"; done
  codesign --verify --strict --verbose=2 "$M_APP/Contents/MacOS/$bin"
done
# Both binaries report the version they were packaged as (a nightly's suffix included).
test "$("$M_APP/Contents/MacOS/FieldCraft" --version)" = "fieldcraft $VERSION"
test "$("$M_APP/Contents/MacOS/fieldcraft-cli" --version)" = "fieldcraft-cli $VERSION"
codesign --verify --strict --deep --verbose=2 "$M_APP"
"$M_APP/Contents/MacOS/fieldcraft-cli" commands >/dev/null
hdiutil_retry detach "$MNT" >/dev/null
trap - EXIT
if [ "$NOTARIZE" = 1 ]; then
  notarize "$DMG"
  xcrun stapler staple "$DMG"
  xcrun stapler validate "$DMG"
  spctl --assess --type open --context context:primary-signature -vvv "$DMG"
fi

# ---- CLI ---------------------------------------------------------------------------------------
echo "==> building $CLI_ZIP"
CLI_DIR="$WORK/fieldcraft-cli-$VERSION-macos-$ARCH"
mkdir -p "$CLI_DIR"
cp "$WORK/bin/fieldcraft-cli" "$CLI_DIR/"
copy_docs "$CLI_DIR"
sign --options runtime --identifier "$APP_ID.cli" "$CLI_DIR/fieldcraft-cli"
codesign --verify --strict --verbose=2 "$CLI_DIR/fieldcraft-cli"
rm -f "$CLI_ZIP"
ditto -c -k --keepParent "$CLI_DIR" "$CLI_ZIP"
# A bare Mach-O can't carry a stapled ticket; Gatekeeper looks the notarization up online.
if [ "$NOTARIZE" = 1 ]; then notarize "$CLI_ZIP"; fi

test "$("$CLI_DIR/fieldcraft-cli" --version)" = "fieldcraft-cli $VERSION"
if [ -n "${GITHUB_OUTPUT:-}" ]; then echo "signing=$SIGNING" >>"$GITHUB_OUTPUT"; fi
echo "==> done (signing: $SIGNING)"
ls -lh "$DMG" "$CLI_ZIP"
