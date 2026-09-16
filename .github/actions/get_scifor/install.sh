#!/usr/bin/env bash
set -euo pipefail

repo=SciFortran/SciFortran
case "${SCIFOR_RUNNER_OS}/${SCIFOR_RUNNER_ARCH}" in
  Linux/X64)
    platform=ubuntu-24.04-x86_64
    if [[ "$(. /etc/os-release; printf '%s' "$ID:$VERSION_ID")" != ubuntu:24.04 ]]; then
      echo "SciFortran release $platform requires an ubuntu-24.04 runner" >&2
      exit 1
    fi
    ;;
  macOS/ARM64)
    platform=macos-15-arm64
    if [[ "$(sw_vers -productVersion)" != 15.* ]]; then
      echo "SciFortran release $platform requires a macos-15 runner" >&2
      exit 1
    fi
    ;;
  *)
    echo "No SciFortran binary release for ${SCIFOR_RUNNER_OS}/${SCIFOR_RUNNER_ARCH}" >&2
    exit 1
    ;;
esac

release=${SCIFOR_REQUESTED_RELEASE:-latest}
if [[ "$release" == latest ]]; then
  release=$(gh release list -R "$repo" --exclude-drafts --limit 100 \
    --json tagName,publishedAt \
    --jq '[.[] | select(.tagName | startswith("scifor-"))] | sort_by(.publishedAt) | last | .tagName // empty')
fi
if [[ "$release" != scifor-* ]]; then
  echo "No SciFortran release found; expected a scifor-* tag, got: $release" >&2
  exit 1
fi

asset="${release}-${platform}.tar.gz"
download_dir="$RUNNER_TEMP/scifor"
mkdir -p "$download_dir"
gh release download "$release" -R "$repo" --pattern "$asset" --dir "$download_dir" --clobber
tar -C "$download_dir" -xzf "$download_dir/$asset"
root="$download_dir/${asset%.tar.gz}"
test -s "$root/lib/libscifor.a"
test -s "$root/lib/pkgconfig/scifor.pc"
test -s "$root/include/scifor.mod"

export PKG_CONFIG_PATH="$root/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --exists scifor
cflags=$(pkg-config --cflags scifor)
libs=$(pkg-config --libs scifor)

{
  echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"
  echo "SCIFOR_ROOT=$root"
  echo "SCIFOR_RELEASE=$release"
  echo "GLOB_INC=$cflags"
  echo "GLOB_LIB=$libs"
} >> "$GITHUB_ENV"
{
  echo "release=$release"
  echo "root=$root"
} >> "$GITHUB_OUTPUT"

echo "Installed $release from $asset at $root"
