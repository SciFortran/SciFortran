Install
#########################


SciFortran is available in the form of a static Fortran library
`libscifor.a` and the related Fortran module `SCIFOR`.
Our build system relies on CMake. Experimental support for linking
Intel's MKL is provided, although it may fail on some systems.


Building SciFortran
======================

Clone the repo:

.. code-block:: bash
		
   git clone https://github.com/SciFortran/SciFortran scifor

   
Optionally [1]_ define the fortran compiler:

.. code-block:: bash
		
   export FC=mpif90/gfortran/ifort


From the repository directory (`cd scifor`) make a standard
out-of-source CMake compilation:


**GNU Make**

Using GNU `make` is the default CMake workflow, with widest version
support (CMake > 3.0). Note that parallel `make` execution is tested
and working.

.. code-block:: bash
		
   mkdir build 
   cd build  
   cmake .. 
   make -j

The `CMake` compilation can be customized using the following
additional variables (default values between `< >`, optional in `[ ]`):   

.. code-block:: bash

   -DLONG_PREFIX = <yes>/no  #set a long or short prefix for the default installation directory  
       
   -DCMAKE_INSTALL_PREFIX    #specify custom installation prefix  
   
   -DUSE_MPI = <yes>/no        #set use of MPI 

   -DVERBOSE = yes/<no>

   -DBUILD_TYPE = <RELEASE>/TESTING/DEBUG/AGGRESSIVE  #compilation options. TESTING:mild or no optimization, DEBUG:relevant debugging options, AGGRESSIVE: all debug options of (might not compile on  some systems) 

   -DWITH_BLAS_LAPACK = yes/<no>   # skip search of preinstalled linear algebra libraries and enforce compilation from local source

   -DWITH_SCALAPACK = <yes>/no       #search and link to available ScaLAPACK library


Install SciFortran
======================
System-wide installation is completed after the build step using
either:

.. code-block:: bash

   make install


To actually link the library to any of your project we provide
different solutions:

* A generated `environment module <https://github.com/cea-hpc/modules>`_, installed to `~/.modules.d/scifor/<PLAT>/<VERSION>`  
* A generated `bash` script at `<PREFIX>/bin/configvars.sh`, to be sourced for permanent loading.
* A generated `pkg-config
  <https://github.com/freedesktop/pkg-config>`_ file to, installed to
  `~/.pkgconfig.d/scifor.pc`
  
which you can choose among by following the instructions printed on screen.


Binary releases for GitHub CI
===============================

Pushes to ``master-release`` build two prerelease assets:
``ubuntu-24.04-x86_64`` and ``macos-15-arm64``. The workflow reads the nearest
version tag reachable from the pushed commit and the first eight characters
of that commit's SHA. Both the GitHub release tag and title have the form
``scifor-<version tag>-<sha8>``; for example, ``scifor-4.23.13-1234abcd``.
Each asset is named ``scifor-<version tag>-<sha8>-<platform>.tar.gz``.
The full commit SHA is recorded in the release notes and ``BUILD-INFO``.
Push the desired commit from ``master`` to ``master-release`` when a binary
release is wanted.
The workflow only publishes after a push to ``master-release``.

The archive contains ``lib/libscifor.a``, Fortran ``.mod`` files under ``include``,
``lib/pkgconfig/scifor.pc``, and ``BUILD-INFO``. It can be extracted into any
directory. It does not contain the compiler, MPI, BLAS or LAPACK libraries.
Install those dependencies with the same setup action used by the build.

In a downstream GitHub Actions job, use the ``SciFortran/checkout`` action. It
installs GNU Fortran, Open MPI, BLAS/LAPACK and ``pkg-config``, then selects and
extracts the newest published binary release for an ``ubuntu-24.04`` or ``macos-15`` arm64
runner. To pin a particular release:

.. code-block:: yaml

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

Follow the action with the consuming project's own configure and build steps.
For a CMake project, run ``cmake -S . -B build`` before ``cmake --build build``.

Omit ``with: release`` to use the most recently published ``scifor-*`` release.
The action sets ``PKG_CONFIG_PATH``, ``SFROOT``, ``SCIFOR_ROOT``, ``SCIFOR_RELEASE``,
``LIBRARY_PATH``, ``LD_LIBRARY_PATH``, ``INCLUDE_PATH``, ``FC``, ``GLOB_INC`` and
``GLOB_LIB`` for later steps in the same job, and provides ``release`` and ``root``
outputs.
The ``uses:`` reference chooses the action version, while ``with: release``
chooses the binary package. To pin both, use a commit SHA in ``uses:`` and a
release tag in ``with:``.

For installation outside GitHub Actions, download the corresponding asset.
For example on Ubuntu, after installing the same dependencies:

.. code-block:: bash

   RELEASE=scifor-4.23.13-1234abcd  # Replace with a published release tag.
   PLATFORM=ubuntu-24.04-x86_64
   ASSET="${RELEASE}-${PLATFORM}.tar.gz"
   gh release download "$RELEASE" --repo SciFortran/SciFortran \
     --pattern "$ASSET" --dir /tmp/scifor-download
   mkdir -p "$HOME/opt"
   tar -C "$HOME/opt" -xzf "/tmp/scifor-download/$ASSET"
   export PKG_CONFIG_PATH="$HOME/opt/${ASSET%.tar.gz}/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
   pkg-config --cflags --libs scifor

Choose ``macos-15-arm64`` on a matching macOS runner. Use an explicit tag with
``gh release download``: these CI builds are prereleases and are not marked as
GitHub's latest stable release.

Compile Fortran consumers using a compatible GNU Fortran version and the
same MPI implementation. The asset's ``BUILD-INFO`` records the build compiler
and link flags. The ``.mod`` files are compiler-dependent; a package built on
one runner is not a universal Linux or macOS binary package.


Uninstall
===================

CMake does not officially provide uninstall procedures in the
generated Make files. Hence SciFortran supplies a homebrew
method to remove the generated files by calling (from the relevant
build folder):

.. code-block:: bash
		
   make uninstall


Known issues
======================
`SciFortran` has been tested with success on several Unix/Linux
platforms. Support for Windows, through `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_, is still experimental, although few people reported successful installation with minimal efforts. 

Some have reported issues concerning the wrong setup for the library `pkg-config` file, contained in  `$PREFIX/<PLAT>/<VERSION>/etc/scifor.pc`. The variable `Libs=-L${libdir} -lscifor <blas/lapack/scalapack>` produced by `cmake` during the configuration and installation process can be not properly defined for the part corresponding to third parties libraries such as Blas/Lapack/Scalapack. This breaks compilation against `scifor` whenever `pkg-config` is used to generate the linking options. 


.. tip::

   Corrupted `pkg-config` file can be fixed manually. Edit the
   `scifor.pc` file  overwriting the definition of the variable
   `Libs`, as appropriate for your system.
   







.. rubric:: footnotes
      
.. [1] In some cases CMake fails to find the MPI fortran compiler,
       even if it is effectively installed and loaded into the
       system. An easy fix is to setup and export the `FC=mpif90`
       environment variable before invoking the `cmake <options> ..`
       command.

