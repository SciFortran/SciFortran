#!/usr/bin/env bash
# Extract elsewhere and compile a consumer without using the build tree.
set -euo pipefail

if (( $# != 1 )); then
  echo "usage: $0 PACKAGE.tar.gz" >&2
  exit 2
fi

archive=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
tar -C "$stage" -xzf "$archive"
prefix=$(find "$stage" -mindepth 1 -maxdepth 1 -type d -print -quit)
test -n "$prefix"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --validate scifor

cat > "$stage/consumer.f90" <<'EOF'
program consumer
  use sf_constants, only: pi, timestamp
  implicit none
  if (abs(pi - acos(-1.0d0)) > 1.0d-12) stop 1
  call timestamp()
end program consumer
EOF

# Deliberately run from the extracted archive, after moving it away from
# the build tree. This checks the .mod files and pkg-config's resolved paths.
${FC:-mpif90} $(pkg-config --cflags scifor) "$stage/consumer.f90" \
  -o "$stage/consumer" $(pkg-config --libs scifor)
"$stage/consumer"
