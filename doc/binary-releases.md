# Binary releases for GitHub CI

Pushes to `master-release` build two prerelease assets:
`ubuntu-24.04-x86_64` and `macos-15-arm64`. Each release is tagged
`ci-<version tag>-g<8-character commit SHA>`; for example,
`ci-4.23.13-g1234abcd`. Each asset is named
`scifor-<version tag>-g<8-character commit SHA>-<platform>.tar.gz`.
The full commit SHA is recorded in the release notes and `BUILD-INFO`.
Push the desired commit from `master` to `master-release` when a binary
release is wanted.
The workflow only publishes after a push to `master-release`.

The archive contains `lib/libscifor.a`, Fortran `.mod` files under `include`,
`lib/pkgconfig/scifor.pc`, and `BUILD-INFO`. It can be extracted into any
directory. It does not contain the compiler, MPI, BLAS or LAPACK libraries.
Install those dependencies with the same setup action used by the build.

In a downstream GitHub Actions job on a matching runner, pin the release
tag and download the corresponding asset. For example on Ubuntu, after
installing GNU Fortran, Open MPI, BLAS/LAPACK and `pkg-config`:

```bash
RELEASE=ci-4.23.13-g1234abcd  # Replace with a published release tag.
PLATFORM=ubuntu-24.04-x86_64
ASSET="scifor-${RELEASE#ci-}-${PLATFORM}.tar.gz"
gh release download "$RELEASE" --repo SciFortran/SciFortran \
  --pattern "$ASSET" --dir /tmp/scifor-download
mkdir -p "$HOME/opt"
tar -C "$HOME/opt" -xzf "/tmp/scifor-download/$ASSET"
export PKG_CONFIG_PATH="$HOME/opt/${ASSET%.tar.gz}/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --cflags --libs scifor
```

In GitHub Actions, set `GH_TOKEN: ${{ github.token }}` on the download step
and write `PKG_CONFIG_PATH` to `$GITHUB_ENV` so later steps can use it.
Choose `macos-15-arm64` instead on a matching macOS runner. Use an explicit
tag with `gh release download`: these CI builds are prereleases and are not
marked as GitHub's latest stable release.

For a downstream repository, a complete Ubuntu setup step looks like this
(replace the example tag in both places):

```yaml
permissions:
  contents: read

jobs:
  build:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4
      - uses: SciFortran/SciFortran/.github/actions@ci-4.23.13-g1234abcd
        with:
          pack-type: open
      - name: Download SciFortran
        env:
          GH_TOKEN: ${{ github.token }}
          RELEASE: ci-4.23.13-g1234abcd
        run: |
          platform=ubuntu-24.04-x86_64
          asset="scifor-${RELEASE#ci-}-${platform}.tar.gz"
          mkdir -p "$RUNNER_TEMP/scifor"
          gh release download "$RELEASE" -R SciFortran/SciFortran \
            -p "$asset" -D "$RUNNER_TEMP/scifor"
          tar -C "$RUNNER_TEMP/scifor" -xzf "$RUNNER_TEMP/scifor/$asset"
          prefix="$RUNNER_TEMP/scifor/${asset%.tar.gz}"
          echo "PKG_CONFIG_PATH=$prefix/lib/pkgconfig" >> "$GITHUB_ENV"
      - name: Build dependent project
        run: |
          pkg-config --cflags --libs scifor
          # Invoke your project's build with FC=mpif90 and the flags above.
```

Compile Fortran consumers using a compatible GNU Fortran version and the
same MPI implementation. The asset's `BUILD-INFO` records the build compiler
and link flags. The `.mod` files are compiler-dependent; a package built on
one runner is not a universal Linux or macOS binary package.
