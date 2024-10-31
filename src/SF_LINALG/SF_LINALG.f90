module SF_LINALG
!SciFortran module for linear algebra
  USE SF_BLACS
  implicit none
  private

  !COMMONLY USED PARAMETERS
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)

  !>EIGENVALUE PROBLEM:
  !Eigenvalue/-vector problem for general matrices:
  public :: eig
  !Eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  public :: eigh
#ifdef _SCALAPACK
  public :: p_eigh
#endif
  !Eigenvalue/-vector problem for real symmetric/complex hermitian matrices using Jacobi method (unsorted out):
  public :: eigh_jacobi
  !Eigenvalues for general matrices:
  public :: eigvals
  !Eigenvalues for symmetric/hermitian matrices:
  public :: eigvalsh
  !
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig
  !
  interface eigh
     module procedure deigh_generalized
     module procedure zeigh_generalized
     module procedure deigh_simple
     module procedure zeigh_simple
     module procedure deigh_tridiag
  end interface eigh
#ifdef _SCALAPACK
  interface p_eigh
     module procedure p_deigh_simple
     module procedure p_zeigh_simple
  end interface p_eigh
#endif
  !
  interface eigh_jacobi
     module procedure d_jacobi
     module procedure c_jacobi
  end interface eigh_jacobi
  !
  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals
  !
  interface eigvalsh
     module procedure deigvalsh
     module procedure zeigvalsh
  end interface eigvalsh





  !>SVD DECOMPOSITION:
  !singular values of real/complex matrices:
  public :: svdvals
  !singular value decomposition of real/complex matrices:
  public :: svd
  !
  interface svdvals
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals
  !
  interface svd
     module procedure dsvd
     module procedure zsvd
  end interface svd





  !>MATRIX INVERSION:
  !Matrix inversion for real/complex matrices:
  public :: inv
#ifdef _SCALAPACK
  public :: p_inv
#endif
  !Matrix inversion for real/complex symmetric matrices:
  public :: inv_sym
  ! matrix inversion for complex hermitian matrices:
  public :: inv_her
  ! matrix inversion for real/complex triangular matrices:
  public :: inv_triang
  ! matrix inversion for real/complex  matrices using Gauss-Jordan elimination:
  public :: inv_gj
  !
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv
#ifdef _SCALAPACK
  interface p_inv
     module procedure p_dinv
     module procedure p_zinv
  end interface p_inv
#endif
  !
  interface inv_sym
     module procedure dinv_sym
     module procedure zinv_sym
  end interface inv_sym
  !
  interface inv_her
     module procedure zinv_her
  end interface inv_her
  !
  interface inv_triang
     module procedure dinv_triang
     module procedure zinv_triang
  end interface inv_triang
  !
  interface inv_gj
     module procedure dinv_gj
     module procedure zinv_gj
  end interface inv_gj





  !>LINEAR SYSTEM SOLUTION:
  ! solution to linear systems of equation with real/complex coefficients:
  public :: solve
  !least square solutions the real/complex systems of equations of possibly non-square shape:
  public :: lstsq
  !
  interface solve
     module procedure dsolve_1rhs
     module procedure zsolve_1rhs
     module procedure dsolve_Mrhs
     module procedure zsolve_Mrhs
  end interface solve
  !
  interface lstsq
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq





  !>TRIDIAGONAL MATRICES BUILD, CHECK and INVERT:
  !invert a (block) tridiagonal matrix using the iterative algorithm (diagonal elements only):
  public :: inv_tridiag
  !Returns the real/complex block identity matrix of size n x n .
  public :: deye_tridiag, zeye_tridiag
  !check the matrix is actually (block) tridiagonal
  public :: check_tridiag
  !get the main (block) diagonals from a (block) tridiagonal matrix
  public :: get_tridiag
  !build a (block) tridigonal matrix from the main (block) diagonals
  public :: build_tridiag
  !
  interface inv_tridiag
     module procedure d_invert_tridiag_matrix
     module procedure c_invert_tridiag_matrix
     module procedure d_invert_tridiag_block_matrix
     module procedure c_invert_tridiag_block_matrix
     module procedure d_invert_tridiag_matrix_mat
     module procedure c_invert_tridiag_matrix_mat
     module procedure d_invert_tridiag_block_matrix_mat
     module procedure c_invert_tridiag_block_matrix_mat
  end interface inv_tridiag
  !
  interface eye_tridiag
     module procedure deye_tridiag
  end interface eye_tridiag
  !
  interface check_tridiag
     module procedure d_check_tridiag
     module procedure c_check_tridiag
     module procedure d_check_tridiag_block
     module procedure c_check_tridiag_block
  end interface check_tridiag
  !
  interface get_tridiag
     module procedure d_get_tridiag
     module procedure c_get_tridiag
     module procedure d_get_tridiag_block
     module procedure c_get_tridiag_block
  end interface get_tridiag
  !
  interface build_tridiag
     module procedure d_build_tridiag
     module procedure c_build_tridiag
     module procedure d_build_tridiag_block
     module procedure c_build_tridiag_block
  end interface build_tridiag








  !>AUXILIARY:
  ! determinants of real/complex square matrices:
  public :: det
  !Returns the real/complex identity matrix of size n x n .
  public :: eye, deye, zeye
  !Returns a matrix of zeros or ones of specified size
  public :: zeros, ones
  !construction of square matrices from the diagonal elements:
  public :: diag
  !get the diagonal from a matrix:
  public :: diagonal
  !trace of real/complex matrices:
  public :: trace
  !
  interface det
     module procedure ddet
     module procedure zdet
  end interface det
  !
  interface deye
     module procedure deye_matrix
     module procedure deye_indices
  end interface deye
  !
  interface zeye
     module procedure zeye_matrix
     module procedure zeye_indices
  end interface zeye
  !
  interface eye
     module procedure deye_matrix
     module procedure deye_indices
  end interface eye
  !
  interface diag
     module procedure ddiag
     module procedure zdiag
  end interface diag
  !
  interface diagonal
     module procedure d_diagonal
     module procedure z_diagonal
  end interface diagonal
  !
  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace
  !
  interface zeros
     module procedure zzeros_1
     module procedure zzeros_2
     module procedure zzeros_3
     module procedure zzeros_4
     module procedure zzeros_5
     module procedure zzeros_6
     module procedure zzeros_7
  end interface zeros
  !
  interface ones
     module procedure zones_1
     module procedure zones_2
     module procedure zones_3
     module procedure zones_4
     module procedure zones_5
     module procedure zones_6
     module procedure zones_7
  end interface ones


  !>EXTERNAL PRODUCTS
  !Kroenecker product of matrices
  public :: kron
  public :: kronecker_product
  public :: operator(.kx.)
  !outer product of two 1d arrays to form a matrix
  public :: outerprod
  public :: cross_product
  public :: s3_product
  !
  interface kron
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface kron
  !
  interface kronecker_product
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface kronecker_product
  !
  interface operator(.kx.)
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface operator(.kx.)
  !
  interface outerprod
     module procedure outerprod_d,outerprod_c
  end interface outerprod
  !
  interface cross_product
     module procedure :: cross_3d_d
     module procedure :: cross_3d_c
  end interface cross_product
  !
  interface s3_product
     module procedure :: s3_product_d
     module procedure :: s3_product_c
  end interface s3_product




  !>BLAS INTERFACE
  !Matrix-matrix product
  public :: mat_product
  public :: operator(.x.)
#ifdef _SCALAPACK
  public :: p_mat_product
  public :: operator(.Px.)
#endif
  !
  interface mat_product
     module procedure :: d_matmul
     module procedure :: z_matmul
  end interface mat_product
  interface operator (.x.)
     module procedure :: d_matmul_
     module procedure :: z_matmul_
  end interface operator (.x.)
#ifdef _SCALAPACK
  interface p_mat_product
     module procedure :: p_d_matmul
     module procedure :: p_z_matmul
  end interface p_mat_product
  interface operator (.Px.)
     module procedure :: p_d_matmul_f
     module procedure :: p_z_matmul_f
  end interface operator (.Px.)
#endif  



  !NOT PUBLIC:
  !Assert shape of matrices:
  interface assert_shape
     module procedure dassert_shape
     module procedure zassert_shape
  end interface assert_shape

  !Swap two elements A and B (from nr90)
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_z,swap_zv,swap_zm
  end interface swap

  !Interface to Lapack function ilaenv
  interface
     integer function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
       character*(*) name, opts
       integer       ispec, n1, n2, n3, n4
     end function ilaenv
  end interface

#ifdef _MPI
#  ifdef _SCALAPACK
  interface Distribute_BLACS
     module procedure :: D_Distribute_BLACS
     module procedure :: Z_Distribute_BLACS
  end interface Distribute_BLACS

  interface Gather_BLACS
     module procedure :: D_Gather_BLACS
     module procedure :: Z_Gather_BLACS
  end interface Gather_BLACS
#  endif
#endif



contains


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  SOLUTION TO EIGEN PROBLEMS
  ! - general matrices (in general non-symmetric or non-complex-hermitian matrices).
  ! - real symmetric/complex hermitian matrices
  ! - Jacobi method
  ! - N-by-N real/complex nonsymmetric matrix A eigenvalues and, optionally,the left/right eigenvectors.
  !-------------------------------------------------------------------------------------------
  include "linalg_eig.f90"
  include "linalg_eigh.f90"
  include "linalg_eigh_jacobi.f90"
  include "linalg_eigvals.f90"
  include "linalg_eigvalsh.f90"
#ifdef _SCALAPACK
  include "linalg_p_eigh.f90"
#endif

  !-------------------------------------------------------------------------------------------
  !PURPOSE: compute singular values s_i of a real/complex m x n matrix A
  !PURPOSE: compute the singular value decomposition A = U sigma Vtransp / U sigma V^H of 
  ! real/complex matrix A
  !-------------------------------------------------------------------------------------------
  include "linalg_svdvals.f90"
  include "linalg_svd.f90"


  !-------------------------------------------------------------------------------------------
  !PURPOSE: INVERSION OF A MATRIX USING LAPACK LIBRARY
  ! - General M*N (D,C)
  ! - Symmetric N*N (D,C)
  ! - Hermitial (C)
  ! - Triangular N*N (D,C)
  ! - Gauss-Jordan
  ! note: M is destroyed and replaces by its inverse M^-1
  !-------------------------------------------------------------------------------------------
  include "linalg_inv.f90"
  include "linalg_inv_sym.f90"
  include "linalg_inv_her.f90"
  include "linalg_inv_triang.f90"
  include "linalg_inv_gj.f90"
#ifdef _SCALAPACK
  include "linalg_p_inv.f90"
#endif

  !+-----------------------------------------------------------------+
  !PROGRAM  : SOLVE LINEAR SYSTEM using LAPACK algorithms
  ! - direct solution of A*x=b
  ! -least square solution to A x = b for real A, b
  !+-----------------------------------------------------------------+
  include "linalg_solve.f90"
  include "linalg_lstsq.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE: wrap BLAS 1,2,3 operations
  !-------------------------------------------------------------------------------------------
  include "linalg_blas.f90"
#ifdef _SCALAPACK
  include "linalg_p_blas.f90"
#endif

  !+-----------------------------------------------------------------+
  !PROGRAM  : BUILD,CHECK & INVERT TRIDIAGONAL MATRICES
  !+-----------------------------------------------------------------+
  include "linalg_inv_tridiag.f90"
  include "linalg_check_tridiag.f90"
  include "linalg_get_tridiag.f90"
  include "linalg_build_tridiag.f90"
  function deye_tridiag(Nblock,N) result(eye_block)
    integer                       :: Nblock
    integer                       :: N
    real(8),dimension(Nblock,N,N) :: eye_block
    integer                       :: iblock
    do iblock=1,Nblock
       eye_block(iblock,:,:) = eye(N)
    enddo
  end function deye_tridiag
  !
  function zeye_tridiag(Nblock,N) result(eye_block)
    integer                          :: Nblock
    integer                          :: N
    complex(8),dimension(Nblock,N,N) :: eye_block
    integer                          :: iblock
    do iblock=1,Nblock
       eye_block(iblock,:,:) = zeye(N)
    enddo
  end function zeye_tridiag






  !+-----------------------------------------------------------------+
  !PROGRAM  : AUXILIARY AND COMPUTATIONAL ROUTINES
  ! - det: compute determinant of a real/complex matrix
  ! - diag: construct real matrix from diagonal elements
  ! - trace: return trace along the main diagonal
  ! - Xeye: returns the identity matrix of size n x n
  !+-----------------------------------------------------------------+
  include "linalg_auxiliary.f90"




  !+-----------------------------------------------------------------+
  !PROGRAM  : EXTERNAL PRODUCTS ROUTINES
  ! - kronecker:  compute the tensor product (M1_kp_M2) of 
  ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
  ! the number of rows and columns of the Matrix M1 and M2
  ! - outerprod: Form a matrix A(:,:) from the outerproduct of two 1d arrays:
  ! A(i,j) = a_i*b_j
  ! - cross: cross or vector product for 2d and 3d vectors.
  ! - s3_product: evaluate the S3 product A.(BxC) for 3d vectors
  !+-----------------------------------------------------------------+
  include "linalg_external_products.f90"











  !##################################################################
  !                  OTHER COMPUTATIONAL ROUTINES
  !##################################################################
#ifdef _MPI
#  ifdef _SCALAPACK
  include "linalg_blacs_aux.f90"
#  endif
#endif


  !-------------------------------------------------------------------------------------------
  !PURPOSE: Asser the correct shape of a matrix
  !-------------------------------------------------------------------------------------------
  subroutine dassert_shape(A, shap, routine, matname)
    real(8),intent(in)  :: A(:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
       print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print*, "Shape should be ", shap
       stop "Aborting due to illegal matrix operation"
    end if
  end subroutine dassert_shape
  subroutine zassert_shape(A, shap, routine, matname)
    complex(8),intent(in) :: A(:,:)
    integer,intent(in)    :: shap(:)
    character(len=*)      :: routine, matname
    if(any(shape(A) /= shap)) then
       print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print*, "Shape should be ", shap
       stop "Aborting due to illegal matrix operation"
    end if
  end subroutine zassert_shape




  function upper_triangle(j,k,extra)
    integer,intent(in)           :: j,k
    integer,optional, intent(in) :: extra
    logical,dimension(j,k)       :: upper_triangle
    integer                      :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth(1,1,j),arth(1,1,k)) < n)
  end function upper_triangle

  function outerdiff(a,b)
    integer,dimension(:), intent(in)   :: a,b
    integer,dimension(size(a),size(b)) :: outerdiff
    outerdiff = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiff

  function arth(first,increment,n) result(arth_i)
    integer,parameter    :: npar_arth=16,npar2_arth=8
    integer,intent(in)   :: first,increment,n
    integer,dimension(n) :: arth_i
    integer              :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= npar_arth) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,npar2_arth
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*npar2_arth
       k=npar2_arth
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  end function arth


  function shift_dw(x_in) result(x)
    integer, intent(in) :: x_in
    integer             :: x
    x=x_in
    x = ior(x,rshift(x, 1))
    x = ior(x,rshift(x, 2))
    x = ior(x,rshift(x, 4))
    x = ior(x,rshift(x, 8))
    x = ior(x,rshift(x, 16))
    x = ior(x,rshift(x, 32))
    x = x + 1
    x = x/2
  end function shift_dw



  function free_unit(n) result(unit_)
    integer,optional :: n
    integer          :: unit_,ios
    logical          :: opened
    unit_=100
    do 
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=opened,iostat=ios)
       if(.not.opened.AND.ios==0)exit 
       if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
    enddo
    if(present(n))n=unit_
  end function free_unit


end module SF_LINALG


