#!/usr/bin/env bash
# Build and package FieldCraft for Linux x86_64:
#
#   $DIST/fieldcraft-<version>-linux-x86_64.tar.gz    bin/fieldcraft, bin/fieldcraft-cli, the
#                                                     .desktop entry, icons and licences
#   $DIST/fieldcraft-<version>-linux-x86_64.AppImage  the app as one file (chmod +x, run)
#   $DIST/fieldcraft_<version>_amd64.deb              Debian/Ubuntu package (cargo-deb)
#   $DIST/fieldcraft-<version>.x86_64.rpm             Fedora/RHEL/openSUSE package (cargo-generate-rpm)
#
# Usage: packaging/linux/package.sh [--skip-build] [--formats "tar appimage deb rpm"]
#
# Build on an old distribution (CI: Ubuntu 22.04, glibc 2.35) so the binaries run on newer ones.
# The AppImage uses appimagetool and the type 2 runtime at pinned releases, checked against their
# SHA-256 (downloaded into $CARGO_TARGET_DIR). The .deb and .rpm need cargo-deb and
# cargo-generate-rpm (the versions release.yml installs); their metadata and file lists are in
# apps/fieldcraft/Cargo.toml. FIELDCRAFT_MAINTAINER sets the .deb's Maintainer field.
# Optional: desktop-file-validate.
# Adapted from SolidCraft's packaging/linux/package.sh (MIT OR Apache-2.0; see ATTRIBUTION.md).
set -euo pipefail
# shellcheck source=../env.sh
. "$(dirname "${BASH_SOURCE[0]}")/../env.sh"
HERE="$ROOT/packaging/linux"
ICONSET="$ROOT/assets/app-icon/icon.iconset"

SKIP_BUILD=0
FORMATS="tar appimage deb rpm"
while [ $# -gt 0 ]; do
  case "$1" in
    --skip-build) SKIP_BUILD=1; shift ;;
    --formats) FORMATS="$2"; shift 2 ;;
    -h | --help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
has() { case " $FORMATS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

ARCH="$(uname -m)"
# The AppImage tools, pinned (github.com/AppImage/appimagetool, github.com/AppImage/type2-runtime;
# MIT). A version bump changes the SHA-256 below too.
APPIMAGETOOL_VERSION=1.9.1
RUNTIME_VERSION=20251108
case "$ARCH" in
  x86_64)
    APPIMAGETOOL_SHA256=ed4ce84f0d9caff66f50bcca6ff6f35aae54ce8135408b3fa33abfc3cb384eb0
    RUNTIME_SHA256=2fca8b443c92510f1483a883f60061ad09b46b978b2631c807cd873a47ec260d
    ;;
  *) echo "unsupported architecture $ARCH (release builds are x86_64)" >&2; exit 2 ;;
esac
APPIMAGETOOL_URL="https://github.com/AppImage/appimagetool/releases/download/$APPIMAGETOOL_VERSION/appimagetool-$ARCH.AppImage"
RUNTIME_URL="https://github.com/AppImage/type2-runtime/releases/download/$RUNTIME_VERSION/runtime-$ARCH"
BASENAME="fieldcraft-$VERSION-linux-$ARCH"

echo "==> FieldCraft $VERSION for Linux $ARCH ($FORMATS)"
if [ "$SKIP_BUILD" = 0 ]; then
  (cd "$ROOT" && cargo build --release --locked -p fieldcraft -p fieldcraft-cli)
fi
BIN="$CARGO_TARGET_DIR/release"
# Every format ships stripped binaries; cargo-generate-rpm packages target/release as it is.
strip "$BIN/fieldcraft" "$BIN/fieldcraft-cli"
WORK="$CARGO_TARGET_DIR/linux-package"
STAGE="$WORK/root/usr"
rm -rf "$WORK"

# ---- an FHS tree, shared by both formats --------------------------------------------------------
install -Dm755 "$BIN/fieldcraft" "$STAGE/bin/fieldcraft"
install -Dm755 "$BIN/fieldcraft-cli" "$STAGE/bin/fieldcraft-cli"
install_share "$STAGE/share"
if command -v desktop-file-validate >/dev/null; then
  desktop-file-validate "$STAGE/share/applications/$APP_ID.desktop"
fi

# ---- .tar.gz ------------------------------------------------------------------------------------
if has tar; then
  mkdir -p "$WORK/tar"
  cp -R "$STAGE" "$WORK/tar/$BASENAME"
  tar -C "$WORK/tar" -czf "$DIST/$BASENAME.tar.gz" "$BASENAME"
  echo "wrote $DIST/$BASENAME.tar.gz"
fi

# ---- .deb and .rpm ------------------------------------------------------------------------------
# From the binaries built and stripped above (--no-build). Semver pre-releases (0.2.0-rc.1) become
# 0.2.0~rc.1 inside the packages so they sort before the release; the file names keep $VERSION.
PKG_VERSION="$(printf '%s' "$VERSION" | sed 's/-/~/')"
if has deb; then
  case "$ARCH" in x86_64) DEB_ARCH=amd64 ;; *) DEB_ARCH="$ARCH" ;; esac
  OUT="$DIST/fieldcraft_${VERSION}_$DEB_ARCH.deb"
  (cd "$ROOT" && cargo deb --locked -p fieldcraft --no-build --deb-version "$PKG_VERSION" \
    --maintainer "${FIELDCRAFT_MAINTAINER:-FieldCraft contributors <jyviko@users.noreply.github.com>}" -o "$OUT")
  echo "wrote $OUT"
fi
if has rpm; then
  OUT="$DIST/fieldcraft-$VERSION.$ARCH.rpm"
  (cd "$ROOT" && cargo generate-rpm -p apps/fieldcraft -s "version = \"$PKG_VERSION\"" -o "$OUT")
  echo "wrote $OUT"
fi

# ---- AppImage -----------------------------------------------------------------------------------
# Usage: fetch <url> <sha256> <file>; downloads once and checks the hash every time.
fetch() {
  if [ ! -s "$3" ]; then curl -fsSL -o "$3.part" "$1" && mv "$3.part" "$3"; fi
  if ! echo "$2  $3" | sha256sum --check --status; then
    echo "error: $3 does not match its pinned SHA-256 ($1)" >&2
    rm -f "$3"
    exit 1
  fi
}
if has appimage; then
  TOOL="$CARGO_TARGET_DIR/appimagetool-$APPIMAGETOOL_VERSION-$ARCH.AppImage"
  RUNTIME="$CARGO_TARGET_DIR/appimage-runtime-$RUNTIME_VERSION-$ARCH"
  fetch "$APPIMAGETOOL_URL" "$APPIMAGETOOL_SHA256" "$TOOL"
  fetch "$RUNTIME_URL" "$RUNTIME_SHA256" "$RUNTIME"
  chmod +x "$TOOL"

  APPDIR="$WORK/FieldCraft.AppDir"
  mkdir -p "$APPDIR"
  cp -R "$STAGE" "$APPDIR/usr"
  install -Dm755 "$HERE/AppRun" "$APPDIR/AppRun"
  cp "$HERE/$APP_ID.desktop" "$APPDIR/$APP_ID.desktop"
  cp "$ICONSET/icon_256x256.png" "$APPDIR/$APP_ID.png"
  ln -s "$APP_ID.png" "$APPDIR/.DirIcon"
  OUT="$DIST/$BASENAME.AppImage"
  rm -f "$OUT"
  # Extract-and-run: appimagetool is itself an AppImage, and CI runners have no FUSE.
  ARCH="$ARCH" APPIMAGE_EXTRACT_AND_RUN=1 "$TOOL" --no-appstream --runtime-file "$RUNTIME" "$APPDIR" "$OUT"
  echo "wrote $OUT"
fi

# Both binaries report the version they were packaged as (a nightly's suffix included).
test "$("$STAGE/bin/fieldcraft" --version)" = "fieldcraft $VERSION"
test "$("$STAGE/bin/fieldcraft-cli" --version)" = "fieldcraft-cli $VERSION"
echo "==> done"
ls -lh "$DIST"
