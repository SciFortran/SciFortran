#!/usr/bin/env bash
# Build a relocatable, CI-oriented archive from a completed CMake build.
set -euo pipefail

if (( $# != 3 )); then
  echo "usage: $0 BUILD_DIR OUTPUT_DIR PLATFORM" >&2
  exit 2
fi

build_dir=$(cd "$1" && pwd)
output_dir=$2
platform=$3
repo_dir=$(cd "$(dirname "$0")/.." && pwd)
sha=${GITHUB_SHA:-$(git -C "$repo_dir" rev-parse HEAD)}
version=$(git -C "$repo_dir" describe --tags --exclude='ci-*' --abbrev=0)
short_sha=${sha:0:8}
name="scifor-${version}-g${short_sha}-${platform}"
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
prefix="$stage/$name"
mkdir -p "$prefix/lib/pkgconfig" "$prefix/include"

test -s "$build_dir/libscifor.a"
test -s "$build_dir/etc/scifor.pc"
test -s "$build_dir/include/scifor.mod"
cp "$build_dir/libscifor.a" "$prefix/lib/"
cp "$build_dir/include/"*.mod "$prefix/include/"

# pcfiledir is resolved by pkg-config at the archive's final location.
# Preserve the dependency flags selected by the build's CMake configuration.
sed 's|^scifor_dir=.*|scifor_dir=${pcfiledir}/../..|' \
  "$build_dir/etc/scifor.pc" > "$prefix/lib/pkgconfig/scifor.pc"

{
  echo "commit=$sha"
  echo "version_tag=$version"
  echo "platform=$platform"
  echo "fortran_compiler=$(${FC:-gfortran} --version | head -n 1)"
  echo "mpi_compiler=$(mpif90 --version | head -n 1)"
  echo "pkg_config_libs=$(PKG_CONFIG_PATH="$prefix/lib/pkgconfig" pkg-config --libs scifor)"
} > "$prefix/BUILD-INFO"

tar -C "$stage" -czf "$output_dir/$name.tar.gz" "$name"
echo "$output_dir/$name.tar.gz"
