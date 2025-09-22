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

