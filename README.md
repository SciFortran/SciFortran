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

