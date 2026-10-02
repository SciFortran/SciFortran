Structure
###################

**SciFortran** is organized as a collection of modules containing
types, constants and routines pertaining to different tasks.

Specific modules or members of a module can be loaded via

.. code-block:: fortran

   USE SF_ARRAYS
   USE SF_LINALG, only: inv,eigh,eye

While the whole library can be imported via

.. code-block:: fortran

   USE SCIFOR
   
**SciFortran** contains the following modules:


.. f:automodule::   scifor

.. toctree::
   :caption: Structure
   :maxdepth: 2
   :glob:
   :hidden:

   modules/*

