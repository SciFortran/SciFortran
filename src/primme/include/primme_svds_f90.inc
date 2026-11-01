!*******************************************************************************
!  Copyright (c) 2018, College of William & Mary                                   
!  All rights reserved.                                                            
!                                                                                  
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:     
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.             
!      * Redistributions in binary form must reproduce the above copyright         
!        notice, this list of conditions and the following disclaimer in the       
!        documentation and/or other materials provided with the distribution.      
!      * Neither the name of the College of William & Mary nor the
!        names of its contributors may be used to endorse or promote products      
!        derived from this software without specific prior written permission.     
!                                                                                  
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
!  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
!  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!  PRIMME: https://github.com/primme/primme
!  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
!*******************************************************************************
!  File: primme_svds_f90.inc
!  
!  Purpose - Main header with the PRIMME SVDS F2008 interface functions.
!  
!******************************************************************************

!-------------------------------------------------------
!     Defining easy to remember labels for setting the 
!     method in primme_svds_set_method from Fortran
!-------------------------------------------------------

integer, parameter :: primme_svds_default = 0
integer, parameter :: primme_svds_hybrid = 1
integer, parameter :: primme_svds_normalequations = 2
integer, parameter :: primme_svds_augmented = 3

!-------------------------------------------------------
!     Defining easy to remember labels for setting the 
!     members of the primme_svds structure from Fortran
!-------------------------------------------------------
 
integer, parameter ::  PRIMME_SVDS_primme                       = 1
integer, parameter ::  PRIMME_SVDS_primmeStage2                 = 2
integer, parameter ::  PRIMME_SVDS_m                            = 3
integer, parameter ::  PRIMME_SVDS_n                            = 4
integer, parameter ::  PRIMME_SVDS_matrixMatvec                 = 5
integer, parameter ::  PRIMME_SVDS_matrixMatvec_type            = 6
integer, parameter ::  PRIMME_SVDS_applyPreconditioner          = 7
integer, parameter ::  PRIMME_SVDS_applyPreconditioner_type     = 8
integer, parameter ::  PRIMME_SVDS_numProcs                     = 9
integer, parameter ::  PRIMME_SVDS_procID                       = 10
integer, parameter ::  PRIMME_SVDS_mLocal                       = 11
integer, parameter ::  PRIMME_SVDS_nLocal                       = 12
integer, parameter ::  PRIMME_SVDS_commInfo                     = 13
integer, parameter ::  PRIMME_SVDS_globalSumReal                = 14
integer, parameter ::  PRIMME_SVDS_globalSumReal_type           = 15
integer, parameter ::  PRIMME_SVDS_broadcastReal                = 16
integer, parameter ::  PRIMME_SVDS_broadcastReal_type           = 17
integer, parameter ::  PRIMME_SVDS_numSvals                     = 18
integer, parameter ::  PRIMME_SVDS_target                       = 19
integer, parameter ::  PRIMME_SVDS_numTargetShifts              = 20
integer, parameter ::  PRIMME_SVDS_targetShifts                 = 21
integer, parameter ::  PRIMME_SVDS_method                       = 22
integer, parameter ::  PRIMME_SVDS_methodStage2                 = 23
integer, parameter ::  PRIMME_SVDS_matrix                       = 24
integer, parameter ::  PRIMME_SVDS_preconditioner               = 25
integer, parameter ::  PRIMME_SVDS_locking                      = 26
integer, parameter ::  PRIMME_SVDS_numOrthoConst                = 27
integer, parameter ::  PRIMME_SVDS_aNorm                        = 28
integer, parameter ::  PRIMME_SVDS_eps                          = 29
integer, parameter ::  PRIMME_SVDS_precondition                 = 30
integer, parameter ::  PRIMME_SVDS_initSize                     = 31
integer, parameter ::  PRIMME_SVDS_maxBasisSize                 = 32
integer, parameter ::  PRIMME_SVDS_maxBlockSize                 = 33
integer, parameter ::  PRIMME_SVDS_maxMatvecs                   = 34
integer, parameter ::  PRIMME_SVDS_iseed                        = 35
integer, parameter ::  PRIMME_SVDS_printLevel                   = 36
integer, parameter ::  PRIMME_SVDS_internalPrecision            = 37
integer, parameter ::  PRIMME_SVDS_outputFile                   = 38
integer, parameter ::  PRIMME_SVDS_stats_numOuterIterations     = 39
integer, parameter ::  PRIMME_SVDS_stats_numRestarts            = 40
integer, parameter ::  PRIMME_SVDS_stats_numMatvecs             = 41
integer, parameter ::  PRIMME_SVDS_stats_numPreconds            = 42
integer, parameter ::  PRIMME_SVDS_stats_numGlobalSum           = 43
integer, parameter ::  PRIMME_SVDS_stats_volumeGlobalSum        = 44
integer, parameter ::  PRIMME_SVDS_stats_numBroadcast           = 45
integer, parameter ::  PRIMME_SVDS_stats_volumeBroadcast        = 46
integer, parameter ::  PRIMME_SVDS_stats_numOrthoInnerProds     = 47
integer, parameter ::  PRIMME_SVDS_stats_elapsedTime            = 48
integer, parameter ::  PRIMME_SVDS_stats_timeMatvec             = 49
integer, parameter ::  PRIMME_SVDS_stats_timePrecond            = 50
integer, parameter ::  PRIMME_SVDS_stats_timeOrtho              = 51
integer, parameter ::  PRIMME_SVDS_stats_timeGlobalSum          = 52
integer, parameter ::  PRIMME_SVDS_stats_timeBroadcast          = 53
integer, parameter ::  PRIMME_SVDS_stats_lockingIssue           = 54
integer, parameter ::  PRIMME_SVDS_convTestFun                  = 55
integer, parameter ::  PRIMME_SVDS_convTestFun_type             = 56
integer, parameter ::  PRIMME_SVDS_convtest                     = 57
integer, parameter ::  PRIMME_SVDS_monitorFun                   = 58
integer, parameter ::  PRIMME_SVDS_monitorFun_type              = 59
integer, parameter ::  PRIMME_SVDS_monitor                      = 60
integer, parameter ::  PRIMME_SVDS_queue                        = 61
integer, parameter ::  PRIMME_SVDS_profile                      = 62 

!-------------------------------------------------------
!    Defining easy to remember labels for setting the 
!    enum members for targeting and operator
!-------------------------------------------------------

integer(kind=c_int64_t), parameter ::  primme_svds_largest = 0
integer(kind=c_int64_t), parameter ::  primme_svds_smallest = 1
integer(kind=c_int64_t), parameter ::  primme_svds_closest_abs = 2
integer(kind=c_int64_t), parameter ::  primme_svds_op_none = 0
integer(kind=c_int64_t), parameter ::  primme_svds_op_AtA = 1
integer(kind=c_int64_t), parameter ::  primme_svds_op_AAt = 2
integer(kind=c_int64_t), parameter ::  primme_svds_op_augmented = 3


!-------------------------------------------------------
! Declare interface
!-------------------------------------------------------

abstract interface
   subroutine primme_svds_matvec_double(x, ldx, y, ldy, blockSize, mode, primme_svds, ierr) bind(c)
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      integer(c_int64_t) :: ldx, ldy
      real(c_double) :: x(ldx,*), y(ldy,*)
      type(c_ptr), value :: primme_svds
      integer(c_int), intent(in) :: blockSize, mode
      integer(c_int) :: ierr
   end subroutine

   subroutine primme_svds_matvec_complexdouble(x, ldx, y, ldy, blockSize, mode, primme_svds, ierr) bind(c)
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      integer(c_int64_t) :: ldx, ldy
      complex(c_double) :: x(ldx,*), y(ldy,*)
      type(c_ptr), value :: primme_svds
      integer(c_int), intent(in) :: blockSize, mode
      integer(c_int) :: ierr
   end subroutine

   subroutine primme_svds_matvec_single(x, ldx, y, ldy, blockSize, mode, primme_svds, ierr) bind(c)
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_float
      integer(c_int64_t) :: ldx, ldy
      real(c_float) :: x(ldx,*), y(ldy,*)
      type(c_ptr), value :: primme_svds
      integer(c_int), intent(in) :: blockSize, mode
      integer(c_int) :: ierr
   end subroutine

   subroutine primme_svds_matvec_complexsingle(x, ldx, y, ldy, blockSize, mode, primme_svds, ierr) bind(c)
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_float
      integer(c_int64_t) :: ldx, ldy
      complex(c_float) :: x(ldx,*), y(ldy,*)
      type(c_ptr), value :: primme_svds
      integer(c_int), intent(in) :: blockSize, mode
      integer(c_int) :: ierr
   end subroutine
end interface

interface
   function primme_svds_params_create() bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr) :: primme_svds_params_create
   end function

   function primme_svds_params_destroy(primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), intent(in), value :: primme_svds
      integer(c_int) :: primme_svds_params_destroy
   end function
end interface

interface primme_svds_get_member
   function primme_svds_get_member_int(primme_svds, label, value) bind(c, name="primme_svds_get_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t
      type(c_ptr), intent(in), value :: primme_svds
      integer(c_int), intent(in), value :: label
      integer(kind=c_int64_t), intent(out) :: value
      integer(c_int) :: primme_svds_get_member_int
   end function

   function primme_svds_get_member_double(primme_svds, label, value) bind(c, name="primme_svds_get_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      real(kind=c_double), intent(out) :: value
      integer(c_int) :: primme_svds_get_member_double
   end function

   function primme_svds_get_member_ptr(primme_svds, label, value) bind(c, name="primme_svds_get_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      type(c_ptr), intent(out) :: value
      integer(c_int) :: primme_svds_get_member_ptr
   end function
end interface primme_svds_get_member
 
interface primme_svds_set_member
   function primme_svds_set_member_int(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t
      type(c_ptr), intent(in), value :: primme_svds
      integer(c_int), intent(in), value :: label
      integer(kind=c_int64_t), intent(in) :: value
      integer(c_int) :: primme_svds_set_member_int
   end function

   function primme_svds_set_member_double(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      real(kind=c_double), intent(in) :: value
      integer(c_int) :: primme_svds_set_member_double
   end function

   function primme_svds_set_member_doubles(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      real(kind=c_double), dimension(*) :: value
      integer(c_int) :: primme_svds_set_member_doubles
   end function

   function primme_svds_set_member_ptr(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      type(c_ptr), intent(in), value :: value
      integer(c_int) :: primme_svds_set_member_ptr
   end function

   function primme_svds_set_member_fun(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      type(c_funptr), intent(in), value :: value
      integer(c_int) :: primme_svds_set_member_fun
   end function
end interface primme_svds_set_member

! NOTE: The following functions will be part of the interface primme_set_member as soon as
!       gfortran supports it
interface
    function primme_svds_set_member_matvec_double(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t
      import :: primme_svds_matvec_double
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      procedure(primme_svds_matvec_double) :: value
      integer(c_int) :: primme_svds_set_member_matvec_double
   end function

   function primme_svds_set_member_matvec_complexdouble(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t
      import :: primme_svds_matvec_complexdouble
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      procedure(primme_svds_matvec_complexdouble) :: value
      integer(c_int) :: primme_svds_set_member_matvec_complexdouble
   end function
    function primme_svds_set_member_matvec_single(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t
      import :: primme_svds_matvec_single
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      procedure(primme_svds_matvec_single) :: value
      integer(c_int) :: primme_svds_set_member_matvec_single
   end function
    function primme_svds_set_member_matvec_complexsingle(primme_svds, label, value) bind(c, name="primme_svds_set_member")
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t
      import :: primme_svds_matvec_complexsingle
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: label
      procedure(primme_svds_matvec_complexsingle) :: value
      integer(c_int) :: primme_svds_set_member_matvec_complexsingle
   end function
end interface

interface
   function primme_svds_set_method(method, methodStage1, methodStage2, primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value :: primme_svds
      integer(c_int), value :: method, methodStage1, methodStage2
      integer(c_int) :: primme_svds_set_method
   end function
end interface

interface xprimme_svds
   function dprimme_svds(svals, svecs, rnorms, primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t, c_double
      real(c_double), dimension(*), intent(out) :: svals, rnorms
      real(c_double), dimension(*) :: svecs
      type(c_ptr), value :: primme_svds
      integer(c_int) :: dprimme_svds
   end function

   function zprimme_svds(svals, svecs, rnorms, primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t, c_double
      real(c_double), dimension(*), intent(out) :: svals, rnorms
      complex(c_double), dimension(*) :: svecs
      type(c_ptr), value :: primme_svds
      integer(c_int) :: zprimme_svds
   end function

   function sprimme_svds(svals, svecs, rnorms, primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t, c_float
      real(c_float), dimension(*), intent(out) :: svals, rnorms
      real(c_float), dimension(*) :: svecs
      type(c_ptr), value :: primme_svds
      integer(c_int) :: sprimme_svds
   end function

   function cprimme_svds(svals, svecs, rnorms, primme_svds) bind(c)
      use iso_c_binding, only: c_ptr, c_funptr, c_int, c_int64_t, c_float
      real(c_float), dimension(*), intent(out) :: svals, rnorms
      complex(c_float), dimension(*) :: svecs
      type(c_ptr), value :: primme_svds
      integer(c_int) :: cprimme_svds
   end function
end interface xprimme_svds
