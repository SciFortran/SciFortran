# Binary releases for GitHub CI

Pushes to `master-release` build two prerelease assets:
`ubuntu-24.04-x86_64` and `macos-15-arm64`. Each release is tagged
`ci-<full commit SHA>` and each asset is named
`scifor-<full commit SHA>-<platform>.tar.gz`. Push the desired commit from
`master` to `master-release` when a binary release is wanted.
The workflow only publishes after a push to `master-release`.

The archive contains `lib/libscifor.a`, Fortran `.mod` files under `include`,
`lib/pkgconfig/scifor.pc`, and `BUILD-INFO`. It can be extracted into any
directory. It does not contain the compiler, MPI, BLAS or LAPACK libraries.
Install those dependencies with the same setup action used by the build.

In a downstream GitHub Actions job on a matching runner, pin the release
commit and download the corresponding asset. For example on Ubuntu:

```bash
SHA=<full-commit-sha-from-master-release>
PLATFORM=ubuntu-24.04-x86_64
gh release download "ci-${SHA}" --repo SciFortran/SciFortran \
  --pattern "scifor-${SHA}-${PLATFORM}.tar.gz" --dir /tmp/scifor-download
mkdir -p "$HOME/opt"
tar -C "$HOME/opt" -xzf "/tmp/scifor-download/scifor-${SHA}-${PLATFORM}.tar.gz"
export PKG_CONFIG_PATH="$HOME/opt/scifor-${SHA}-${PLATFORM}/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --cflags --libs scifor
```

Compile Fortran consumers using a compatible GNU Fortran version and the
same MPI implementation. The asset's `BUILD-INFO` records the build compiler
and link flags. The `.mod` files are compiler-dependent; a package built on
one runner is not a universal Linux or macOS binary package.
