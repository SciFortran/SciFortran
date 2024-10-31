module SF_SP_LINALG
!SciFortran module for Arpack/Lanczos diagonalization
  USE SF_MISC,   only: assert_shape
  USE SF_RANDOM, only: mt_random
  USE SF_LINALG, only: eye,eigh
#ifdef _MPI
  USE SF_MPI
#endif
  implicit none
#ifdef _MPI
  include 'mpif.h'
#endif
  private
  

  interface sp_eigh
     module procedure :: lanczos_arpack_d
     module procedure :: lanczos_arpack_c
#ifdef _MPI
     module procedure :: lanczos_parpack_d
     module procedure :: lanczos_parpack_c
#endif
  end interface sp_eigh


  interface sp_lanc_eigh
     module procedure :: lanczos_eigh_d
     module procedure :: lanczos_eigh_c
#ifdef _MPI
     module procedure :: mpi_lanczos_eigh_d
     module procedure :: mpi_lanczos_eigh_c
#endif
  end interface sp_lanc_eigh


  interface sp_lanc_tridiag
     module procedure :: lanczos_tridiag_d
     module procedure :: lanczos_tridiag_c
#ifdef _MPI
     module procedure :: mpi_lanczos_tridiag_d
     module procedure :: mpi_lanczos_tridiag_c
#endif
  end interface sp_lanc_tridiag


  interface sp_dvdson_eigh
     module procedure :: dvdson_eigh_d
  end interface sp_dvdson_eigh


  complex(8),parameter              :: zero=(0d0,0d0)
  complex(8),parameter              :: one=(1d0,0d0)
  complex(8),parameter              :: xi=(0d0,1d0)
  integer,allocatable               :: seed_random(:)
  integer                           :: nrandom
  logical                           :: verb=.false.
  real(8)                           :: threshold_=1.d-12
  integer                           :: ncheck_=10



  !****************************************************************************************
  !                                      PUBLIC 
  !****************************************************************************************
  public :: sp_eigh
  public :: sp_lanc_eigh
  public :: sp_lanc_tridiag
  public :: sp_dvdson_eigh
  !****************************************************************************************




contains


  !##################################################################
  ! ARPACK METHOD for LOWEST part of the spectrum of a of a SPARSE
  !   MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  !        [COMM    MPI  Communicator for the processor grid.  (INPUT)]
  include "arpack_d.f90"
  include "arpack_c.f90"
#ifdef _MPI
  include "parpack_d.f90"
  include "parpack_c.f90"
#endif


  !##################################################################
  ! LANCZOS METHOD for LOWEST EigenSolution OR tri-diagonalization
  !    of a SPARSE MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  include "lanczos_d.f90"
  include "lanczos_c.f90"
#ifdef _MPI
  include "mpi_lanczos_d.f90"
  include "mpi_lanczos_c.f90"
#endif


  !##################################################################
  ! DAVIDSON METHOD for LOWEST EigenSolution of a SPARSE MATRIX (defined via H*v)
  include "dvdson_serial.f90"



end module SF_SP_LINALG
























! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !---------------------------------------------------------------------
! ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
! !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
! !    tridiagonal matrix by the QL method.  The eigenvectors of a full
! !    symmetric matrix can also be found if TRED2 has been used to reduce this
! !    full matrix to tridiagonal form.
! !  Parameters:
! !    Input, integer ( kind = 4 ) N, the order of the matrix.
! !
! !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
! !    the matrix.  On output, the eigenvalues in ascending order.  If an error
! !    exit is made, the eigenvalues are correct but unordered for indices
! !    1,2,...,IERR-1.
! !
! !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
! !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
! !    On output, E has been destroyed.
! !
! !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
! !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
! !    the tridiagonal matrix are desired, Z must contain the identity matrix.
! !    On output, Z contains the orthonormal eigenvectors of the symmetric
! !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
! !    the eigenvectors associated with the stored eigenvalues.
! !
! !    Output, integer ( kind = 4 ) IERR, error flag.
! !    0, normal return,
! !    J, if the J-th eigenvalue has not been determined after
! !    30 iterations.
! !
! !---------------------------------------------------------------------
! subroutine tql2 ( n, d, e, z, ierr )
!   integer :: n
!   real(8) :: c
!   real(8) :: c2
!   real(8) :: c3
!   real(8) :: d(n)
!   real(8) :: dl1
!   real(8) :: e(n)
!   real(8) :: el1
!   real(8) :: f
!   real(8) :: g
!   real(8) :: h
!   integer ( kind = 4 ) i
!   integer ( kind = 4 ) ierr
!   integer ( kind = 4 ) ii
!   integer ( kind = 4 ) j
!   integer ( kind = 4 ) k
!   integer ( kind = 4 ) l
!   integer ( kind = 4 ) l1
!   integer ( kind = 4 ) l2
!   integer ( kind = 4 ) m
!   integer ( kind = 4 ) mml
!   real(8) :: p
!   real(8) :: r
!   real(8) :: s
!   real(8) :: s2
!   real(8) :: tst1
!   real(8) :: tst2
!   real(8) :: z(n,n)
!   ierr = 0
!   if ( n == 1 ) then
!      return
!   end if
!   do i = 2, n
!      e(i-1) = e(i)
!   end do
!   f = 0.0D+00
!   tst1 = 0.0D+00
!   e(n) = 0.0D+00
!   do l = 1, n
!      j = 0
!      h = abs ( d(l) ) + abs ( e(l) )
!      tst1 = max ( tst1, h )
!      !
!      !  Look for a small sub-diagonal element.
!      !
!      do m = l, n
!         tst2 = tst1 + abs ( e(m) )
!         if ( tst2 == tst1 ) then
!            exit
!         end if
!      end do
!      if ( m == l ) then
!         go to 220
!      end if
! 130  continue
!      if ( 30 <= j ) then
!         ierr = l
!         return
!      end if
!      j = j + 1
!      !
!      !  Form shift.
!      !
!      l1 = l + 1
!      l2 = l1 + 1
!      g = d(l)
!      p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
!      r = pythag ( p, 1.0D+00 )
!      d(l) = e(l) / ( p + sign ( r, p ) )
!      d(l1) = e(l) * ( p + sign ( r, p ) )
!      dl1 = d(l1)
!      h = g - d(l)
!      d(l2:n) = d(l2:n) - h
!      f = f + h
!      !
!      !  QL transformation.
!      !
!      p = d(m)
!      c = 1.0D+00
!      c2 = c
!      el1 = e(l1)
!      s = 0.0D+00
!      mml = m - l
!      do ii = 1, mml
!         c3 = c2
!         c2 = c
!         s2 = s
!         i = m - ii
!         g = c * e(i)
!         h = c * p
!         r = pythag ( p, e(i) )
!         e(i+1) = s * r
!         s = e(i) / r
!         c = p / r
!         p = c * d(i) - s * g
!         d(i+1) = h + s * ( c * g + s * d(i) )
!         !
!         !  Form vector.
!         !
!         do k = 1, n
!            h = z(k,i+1)
!            z(k,i+1) = s * z(k,i) + c * h
!            z(k,i) = c * z(k,i) - s * h
!         end do
!      end do
!      p = - s * s2 * c3 * el1 * e(l) / dl1
!      e(l) = s * p
!      d(l) = c * p
!      tst2 = tst1 + abs ( e(l) )
!      if ( tst2 > tst1 ) then
!         go to 130
!      end if
! 220  continue
!      d(l) = d(l) + f
!   end do
!   !
!   !  Order eigenvalues and eigenvectors.
!   !
!   do ii = 2, n
!      i = ii - 1
!      k = i
!      p = d(i)
!      do j = ii, n
!         if ( d(j) < p ) then
!            k = j
!            p = d(j)
!         end if
!      end do
!      if ( k /= i ) then
!         d(k) = d(i)
!         d(i) = p
!         do j = 1, n
!            call r8_swap ( z(j,i), z(j,k) )
!         end do
!      end if
!   end do
!   return
! end subroutine tql2


! !---------------------------------------------------------------------
! ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
! !    The formula
! !    PYTHAG = sqrt ( A * A + B * B )
! !    is reasonably accurate, but can fail if, for example, A**2 is larger
! !    than the machine overflow.  The formula can lose most of its accuracy
! !    if the sum of the squares is very large or very small.
! !  Parameters:
! !    Input, real(8) :: A, B, the two legs of a right triangle.
! !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
! !---------------------------------------------------------------------
! function pythag ( a, b )
!   implicit none
!   real(8) :: a
!   real(8) :: b
!   real(8) :: p
!   real(8) :: pythag
!   real(8) :: r
!   real(8) :: s
!   real(8) :: t
!   real(8) :: u
!   p = max ( abs ( a ), abs ( b ) )
!   if ( p /= 0.0D+00 ) then
!      r = ( min ( abs ( a ), abs ( b ) ) / p )**2
!      do
!         t = 4.0D+00 + r
!         if ( t == 4.0D+00 ) then
!            exit
!         end if
!         s = r / t
!         u = 1.0D+00 + 2.0D+00 * s
!         p = u * p
!         r = ( s / u )**2 * r
!      end do
!   end if
!   pythag = p
!   return
! end function pythag

! !---------------------------------------------------------------------
! ! PURPOSE: swaps two R8's.
! !  Parameters:
! !    Input/output, real(8) :: X, Y.  On output, the values of X and
! !    Y have been interchanged.
! !---------------------------------------------------------------------
! subroutine r8_swap ( x, y )
!   real(8) :: x
!   real(8) :: y
!   real(8) :: z
!   z = x
!   x = y
!   y = z
!   return
! end subroutine r8_swap



