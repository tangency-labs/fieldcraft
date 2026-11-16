#!/usr/bin/env bash
# Build FieldCraft from source as a single-file Flatpak bundle:
#
#   $DIST/fieldcraft-<version>-linux-x86_64.flatpak
#
# Usage: packaging/linux/flatpak-bundle.sh [--no-test]
#
# Builds the committed tree (HEAD; uncommitted changes are not included) with flatpak-builder
# from packaging/linux/flatpak/io.github.tangencylabs.fieldcraft.yml. The crates are fetched here with
# `cargo vendor --locked`, so the build inside flatpak-builder runs offline. Needs flatpak,
# flatpak-builder, git, cargo and THIRD-PARTY-LICENSES.txt (`cargo xtask licences`); the
# freedesktop runtime, SDK and Rust extension named in the manifest are installed per user from
# Flathub. FIELDCRAFT_POSTHOG_KEY, when set, is built in as in the other packages
# (docs/privacy.md). Unless --no-test, the bundle is installed per user and both programs print
# their version.
# After SolidCraft's packaging/linux/flatpak-bundle.sh (MIT OR Apache-2.0; see ATTRIBUTION.md).
set -euo pipefail
# shellcheck source=../env.sh
. "$(dirname "${BASH_SOURCE[0]}")/../env.sh"
HERE="$ROOT/packaging/linux"

TEST=1
while [ $# -gt 0 ]; do
  case "$1" in
    --no-test) TEST=0; shift ;;
    -h | --help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

ARCH="$(uname -m)"
[ "$ARCH" = x86_64 ] || { echo "unsupported architecture $ARCH (release builds are x86_64)" >&2; exit 2; }
for tool in flatpak flatpak-builder git cargo; do
  command -v "$tool" >/dev/null || { echo "error: $tool not found (apt install flatpak flatpak-builder)" >&2; exit 1; }
done
[ -s "$ROOT/THIRD-PARTY-LICENSES.txt" ] || { echo "error: THIRD-PARTY-LICENSES.txt is missing (cargo xtask licences)" >&2; exit 1; }

echo "==> FieldCraft $VERSION Flatpak bundle for $ARCH from $(git -C "$ROOT" rev-parse --short HEAD)"
WORK="$CARGO_TARGET_DIR/flatpak-bundle"
SRC="$WORK/source"
rm -rf "$SRC" "$WORK/build" "$WORK/repo"
mkdir -p "$SRC"
git -C "$ROOT" archive --format=tar HEAD | tar -x -C "$SRC"
cp "$ROOT/THIRD-PARTY-LICENSES.txt" "$SRC/"
# The vendored crates, and the source replacement that makes cargo use them (`cargo vendor`
# prints it; its relative path resolves against the folder holding .cargo/).
(cd "$SRC" && cargo vendor --locked vendor) >>"$SRC/.cargo/config.toml"
{
  printf 'export FIELDCRAFT_VERSION=%q\n' "$VERSION"
  printf 'export FIELDCRAFT_RELEASE_DATE=%q\n' "$(git -C "$ROOT" log -1 --format=%cs)"
  printf 'export FIELDCRAFT_BUILD_SHA=%q\n' "$FIELDCRAFT_BUILD_SHA"
  if [ -n "${FIELDCRAFT_POSTHOG_KEY:-}" ]; then printf 'export FIELDCRAFT_POSTHOG_KEY=%q\n' "$FIELDCRAFT_POSTHOG_KEY"; fi
} >"$SRC/flatpak-build-env.sh"
cp "$HERE/flatpak/$APP_ID.yml" "$WORK/$APP_ID.yml"

flatpak remote-add --user --if-not-exists flathub "$FLATHUB"
# --disable-rofiles-fuse: no FUSE in containers and CI runners. The state folder keeps the
# build cache between local runs.
flatpak-builder --user --install-deps-from=flathub --disable-rofiles-fuse --force-clean \
  --default-branch=stable --repo="$WORK/repo" --state-dir="$WORK/state" "$WORK/build" "$WORK/$APP_ID.yml"
# The analytics key was in the staged source only.
rm -rf "$SRC"

OUT="$DIST/fieldcraft-$VERSION-linux-$ARCH.flatpak"
# --runtime-repo lets `flatpak install` fetch the freedesktop runtime from Flathub if missing.
flatpak build-bundle --runtime-repo="$FLATHUB" "$WORK/repo" "$OUT" "$APP_ID" stable
echo "wrote $OUT"

if [ "$TEST" = 1 ]; then
  flatpak install --user -y --noninteractive --reinstall "$OUT"
  flatpak info --user "$APP_ID"
  test "$(flatpak run --command=fieldcraft-cli "$APP_ID" --version)" = "fieldcraft-cli $VERSION"
  test "$(flatpak run "$APP_ID" --version)" = "fieldcraft $VERSION"
  flatpak run --command=sh "$APP_ID" -c "ls /app/share/applications /app/share/metainfo /app/share/icons/hicolor/256x256/apps /app/share/doc/fieldcraft"
fi
echo "==> done"
ls -lh "$OUT"
