# Binary releases for GitHub CI

Pushes to `master-release` build two prerelease assets:
`ubuntu-24.04-x86_64` and `macos-15-arm64`. The workflow reads the nearest
version tag reachable from the pushed commit and the first eight characters
of that commit's SHA. Both the GitHub release tag and title have the form
`scifor-<version tag>-<sha8>`; for example, `scifor-4.23.13-1234abcd`.
Each asset is named `scifor-<version tag>-<sha8>-<platform>.tar.gz`.
The full commit SHA is recorded in the release notes and `BUILD-INFO`.
Push the desired commit from `master` to `master-release` when a binary
release is wanted.
The workflow only publishes after a push to `master-release`.

The archive contains `lib/libscifor.a`, Fortran `.mod` files under `include`,
`lib/pkgconfig/scifor.pc`, and `BUILD-INFO`. It can be extracted into any
directory. It does not contain the compiler, MPI, BLAS or LAPACK libraries.
Install those dependencies with the same setup action used by the build.

In a downstream GitHub Actions job, use the `SciFortran/checkout` action. It
installs GNU Fortran, Open MPI, BLAS/LAPACK and `pkg-config`, then selects and
extracts the newest published binary release for an `ubuntu-24.04` or `macos-15` arm64
runner. To pin a particular release:

```yaml
permissions:
  contents: read

jobs:
  build:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4
      - id: scifor
        uses: SciFortran/checkout@main
        with:
          release: scifor-4.23.13-1234abcd
```

Follow the action with the consuming project's own configure and build steps.
For a CMake project, run `cmake -S . -B build` before `cmake --build build`.

Omit `with: release` to use the most recently published `scifor-*` release.
The action sets `PKG_CONFIG_PATH`, `SFROOT`, `SCIFOR_ROOT`, `SCIFOR_RELEASE`,
`LIBRARY_PATH`, `LD_LIBRARY_PATH`, `INCLUDE_PATH`, `FC`, `GLOB_INC` and
`GLOB_LIB` for later steps in the same job, and provides `release` and `root`
outputs.
The `uses:` reference chooses the action version, while `with: release`
chooses the binary package. To pin both, use a commit SHA in `uses:` and a
release tag in `with:`.

For installation outside GitHub Actions, download the corresponding asset.
For example on Ubuntu, after installing the same dependencies:

```bash
RELEASE=scifor-4.23.13-1234abcd  # Replace with a published release tag.
PLATFORM=ubuntu-24.04-x86_64
ASSET="${RELEASE}-${PLATFORM}.tar.gz"
gh release download "$RELEASE" --repo SciFortran/SciFortran \
  --pattern "$ASSET" --dir /tmp/scifor-download
mkdir -p "$HOME/opt"
tar -C "$HOME/opt" -xzf "/tmp/scifor-download/$ASSET"
export PKG_CONFIG_PATH="$HOME/opt/${ASSET%.tar.gz}/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --cflags --libs scifor
```

Choose `macos-15-arm64` on a matching macOS runner. Use an explicit tag with
`gh release download`: these CI builds are prereleases and are not marked as
GitHub's latest stable release.

Compile Fortran consumers using a compatible GNU Fortran version and the
same MPI implementation. The asset's `BUILD-INFO` records the build compiler
and link flags. The `.mod` files are compiler-dependent; a package built on
one runner is not a universal Linux or macOS binary package.
