Dependencies
#########################


Installation of `SciFortran` requires a number of third party software and, quite
obviously, a modern Fortran compiler. A restricted number of
this software is mandatory while another part can be considered
optional, though it is still strongly recommended to get maximum
performance and access all `SciFortran` features. 


* `GNU Fortran`_ > 5.0 *OR* `Intel Fortran Compiler`_ > 13.0
  
*  CMake_ ≥ 3.5
   
*  Make_ *OR* Ninja_ ≥ 1.10


**Optional**

*  MPI_
   
* `BLAS/LAPACK`_
  
*  ScaLAPACK_


.. tip::
   
   Should Lapack/Blas not be available in your system, SciFortran
   will compile internal copies of such libraries. This option,
   however, in most cases leads to noticeable performance
   degradation, with respect to optimized versions of the same
   libraries, such as `Intel  MKL`_, `Apple vecLib`_, OpenBLAS_, etc.
   Support for linking MKL is offered via a custom `CMake
   macro`_   which should however be considered in beta development,
   as it is not universal and may end up in wrong linking
   directives. Apple's vecLib is instead automatically recognized by
   CMake, in any standard macOS system.
   
   
.. note::
   
   If any of the required libraries is not available in your system, or
   the version requirement is not satisfied, please install/upgrade
   them. We generally advice for pre-packaged versions, as provided by
   either apt_, pip_, brew_, conda_ or spack_.
   The latter may provide the best options for HPC environments
   (trivial installation without sudo, easy management of
   interdependent multiple versions, automatic generation of
   environment modules, etc.), but choose freely according to your
   needs.



.. _GNU Fortran: https://gcc.gnu.org/fortran/
.. _Intel Fortran Compiler: https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html
.. _CMake: https://cmake.org
.. _Make: https://www.gnu.org/software/make
.. _Ninja: https://ninja-build.org
.. _apt: https://en.wikipedia.org/wiki/APT_(software)
.. _pip: https://pypi.org/project/pip/
.. _brew: https://formulae.brew.sh/
.. _conda: https://docs.conda.io/en/latest/
.. _spack: https://spack.io/
.. _Intel  MKL: https://en.wikipedia.org/wiki/Math_Kernel_Library
.. _Apple vecLib: https://developer.apple.com/documentation/accelerate/veclib
.. _OpenBLAS: https://www.openblas.net
.. _CMake macro:  https://github.com/SciFortran/SciFortran/blob/master/cmake/FindMKL.cmake
.. _MPI: https://github.com/open-mpi/ompi
.. _BLAS/LAPACK: https://github.com/Reference-LAPACK/lapack 
.. _ScaLAPACK: https://github.com/Reference-ScaLAPACK/scalapack


       

  

