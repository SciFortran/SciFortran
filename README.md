# SciFortran

[![Ubuntu](https://img.shields.io/github/actions/workflow/status/QcmPlab/SciFortran/Ubuntu_Scheduled.yml?label=Ubuntu&logo=ubuntu&style=flat-square)](https://github.com/SciFortran/SciFortran/actions/workflows/Scheduled.yml) 
[![MacOS](https://img.shields.io/github/actions/workflow/status/QcmPlab/SciFortran/MacOS_Scheduled.yml?label=macOS&logo=apple&style=flat-square)](https://github.com/SciFortran/SciFortran/actions/workflows/Scheduled.yml) 
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://SciFortran.github.io/SciFortran)


## An open-source Fortran library for mathematics, science and engineering.

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

#### Dependencies

* [GNU Fortran (`gfortran`)](https://gcc.gnu.org/fortran/) > 5.0 **OR** [Intel Fortran Compiler Classic (`ifort`)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)  > 13.0
* [CMake](https://cmake.org/) ≥ 3.5 
* [Make](https://www.gnu.org/software/make/) **OR** [Ninja](https://ninja-build.org/) ≥ 1.10 

See documentation for further details:
[SciFortran.github.io/SciFortran](https://scifortran.github.io/SciFortran/dependencies.html)


## BUILD & INSTALL 

Detailed instructions for building and installing `SciFor` please read the documentation:
[SciFortran.github.io/SciFortran](https://scifortran.github.io/SciFortran/installation.html)

## GitHub workflows and binary releases

The three workflow branches have separate purposes. Push a commit from local
`master` to the branch whose workflow you want to run:

| Branch | Workflow | Result |
| --- | --- | --- |
| `master-doc` | `documentation.yml` | Build the documentation and publish it to GitHub Pages. |
| `master-ci` | `CI.yml` | Compile SciFortran on Ubuntu and macOS; no tests or installation. |
| `master-release` | `Release.yml` | Build and publish Ubuntu x86_64 and macOS arm64 binary archives. |

For example, `git push origin master:master-ci` runs the compilation workflow,
while `git push origin master:master-release` creates a binary release if both
platform builds and package checks pass. The release workflow reads the nearest
version tag reachable from the pushed commit and its eight-character commit
SHA. A release tag and title look like `scifor-4.23.13-1234abcd`; its assets
are named `scifor-4.23.13-1234abcd-ubuntu-24.04-x86_64.tar.gz` and
`scifor-4.23.13-1234abcd-macos-15-arm64.tar.gz`. These per-push builds are
marked as prereleases. Find them under [Releases](https://github.com/SciFortran/SciFortran/releases)
and download the platform archive from **Assets**. The automatically generated
"Source code" archives do not contain the compiled library.

In another repository, the `SciFortran/checkout` composite action installs the
dependencies and the newest published SciFortran release for the runner:

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
```

Add the consuming project's own configure and build steps after the action.
For a CMake project, configure a build directory first and then use
`cmake --build` on that directory.

The action selects the newest published `scifor-*` prerelease, downloads the
asset matching `ubuntu-24.04` or `macos-15` (arm64), and makes
`PKG_CONFIG_PATH`, `SFROOT`, `SCIFOR_ROOT`, `SCIFOR_RELEASE`, `LIBRARY_PATH`,
`LD_LIBRARY_PATH`, `INCLUDE_PATH`, `FC`, `GLOB_INC` and `GLOB_LIB`
available to later steps in the same job. It also exposes `release` and `root`
outputs. To reproduce a build, set
`with: release: scifor-4.23.13-1234abcd` on the action.
That pins the binary release; use a commit SHA instead of `main` in `uses:`
if you also want to pin the action code. See
[binary release details](doc/binary-releases.md).
If a new `master-release` workflow is still running, this selects the previous
successfully published release.
Because CI releases are prereleases, the action selects a tag explicitly.

The setup action installs the external compiler, MPI, BLAS and LAPACK
dependencies. The archive supplies `libscifor.a`, Fortran `.mod` files and
`scifor.pc`. Use a compatible GNU Fortran and MPI setup on the consuming runner;
the `.mod` files are compiler-dependent. For macOS, use `macos-15` and the
`macos-15-arm64` asset. See [binary release details](doc/binary-releases.md).



## AUTHORS
[Adriano Amaricci](https://github.com/aamaricci)  
[Lorenzo Crippa](https://github.com/lcrippa)  
[Samuele Giuli](https://github.com/SamueleGiuli)  
[Gabriele Bellomia](https://github.com/beddalumia)  
[Giacomo Mazza](https://github.com/GiacMazza)


If you encounter bugs or difficulties, please [file an issue](https://github.com/SciFortran/SciFortran/issues/new/choose). For any other communication, please reach out any of the contributors or developers:         


--

***LICENSE***  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LGPL for more details.

You should have received a copy of the GNU LGPL along with this program.  If not, see <http://www.gnu.org/licenses/>.
