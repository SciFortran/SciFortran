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
  !This subroutine computes the eigenvalues and the right eigenvectors of a general (non-symmetric real or non-hermitian
  !complex) square matrix :f:var:`A`, :math:`A v_j = \lambda_j v_j`, using the LAPACK routines :f:func_inline:`dgeev` (real,
  !:f:func_inline:`deig`) and :f:func_inline:`zgeev` (complex, :f:func_inline:`zeig`).
  !
  !The input matrix is not modified. The eigenvalues :f:var:`Eval` are complex and are returned in the order given by LAPACK.
  !The eigenvectors are stored column-wise in :f:var:`Evec`, :code:`Evec(:,j)` being the eigenvector associated to
  !:code:`Eval(j)`. The right eigenvectors are computed unless :code:`jobvr='N'`. The left eigenvectors are not available: the
  !program stops if :code:`jobvl='V'`, if :f:var:`A` is not square, or if :f:var:`Evec` does not have shape :code:`[n,n]`.
  !
     module procedure deig
     module procedure zeig
  end interface eig
  !
  interface eigh
  !This subroutine solves the eigenvalue problem of real symmetric or complex hermitian matrices. The eigenvalues are real and
  !are returned in ascending order. The specific procedures cover
  !
  !* the generalized problem :math:`A c = \lambda B c`, with :f:var:`Am` symmetric/hermitian and :f:var:`Bm` positive definite
  !  (:f:func_inline:`dsygvd`, :f:func_inline:`zhegvd`), returning the eigenvalues :f:var:`lam` and the eigenvectors :f:var:`c`,
  !  stored column-wise and normalized as :math:`c^H B c = 1`. Only the lower triangles of :f:var:`Am` and :f:var:`Bm` are
  !  referenced and the matrices are not modified
  !* the simple problem :math:`A v = E v`, solved in place with the LAPACK driver chosen by :f:var:`method`:
  !  :f:func_inline:`dsyevd` (default), :f:func_inline:`dsyevr`, :f:func_inline:`dsyev` or :f:func_inline:`dsyevx`, and the
  !  :code:`zhe...` counterparts for hermitian matrices. If :code:`jobz='V'` (default) the eigenvectors overwrite the columns of
  !  :f:var:`A`, otherwise :f:var:`A` is destroyed. The triangle :f:var:`uplo` (default :code:`'U'`) is used. A subset of the
  !  spectrum can be selected, with the :code:`r` and :code:`x` drivers only, by value (:f:var:`vl`, :f:var:`vu`) or by index
  !  (:f:var:`il`, :f:var:`iu`)
  !* a real symmetric tridiagonal matrix, given by its main diagonal :f:var:`D` and off-diagonal :f:var:`U` and diagonalized
  !  with :f:func_inline:`dstevr`. :f:var:`D` and :f:var:`U` are overwritten, and :f:var:`D` contains the eigenvalues on output.
  !  The eigenvectors are computed only if :f:var:`Ev` is present, and a subset of the spectrum can be selected by index with
  !  :f:var:`Irange` or by value with :f:var:`Vrange`
  !
  !The program stops if the shapes of the arguments are inconsistent, or if the LAPACK routine returns an error code.
  !
     module procedure deigh_generalized
     module procedure zeigh_generalized
     module procedure deigh_simple
     module procedure zeigh_simple
     module procedure deigh_tridiag
  end interface eigh
#ifdef _SCALAPACK
  interface p_eigh
  !This subroutine is the distributed-memory counterpart of :f:func_inline:`eigh` for the simple problem, for real symmetric
  !(:f:func_inline:`p_deigh_simple`) and complex hermitian (:f:func_inline:`p_zeigh_simple`) matrices, using ScaLAPACK. It is
  !available only when compiling with :code:`_SCALAPACK`.
  !
  !The matrix :f:var:`A`, a full copy of which is expected on every process, is distributed on the BLACS process grid of the
  !module in a block-cyclic layout with square blocks of size :f:var:`Nblock`, and diagonalized with the driver chosen by
  !:f:var:`method`: :f:func_inline:`PDSYEVR` (default), :f:func_inline:`PDSYEV`, :f:func_inline:`PDSYEVD` or
  !:f:func_inline:`PDSYEVX`. The eigenvalues are returned in ascending order in :f:var:`W`, and if :code:`jobz='V'` the
  !eigenvectors are gathered back into :f:var:`A`, column-wise.
  !
     module procedure p_deigh_simple
     module procedure p_zeigh_simple
  end interface p_eigh
#endif
  !
  interface eigh_jacobi
  !This subroutine computes all the eigenvalues and eigenvectors of a real symmetric (:f:func_inline:`d_jacobi`) or complex
  !hermitian (:f:func_inline:`c_jacobi`) matrix with the Jacobi method of plane rotations. Only the upper triangle of the input
  !matrix is referenced and the matrix is destroyed. The eigenvalues are NOT sorted.
  !
  !* real case: the eigenvectors are stored column-wise in :f:var:`v`, :code:`v(:,j)` being associated to :code:`d(j)`. The
  !  number of rotations performed is returned in :f:var:`nrot`
  !* complex case: the eigenvectors are stored in the ROWS of :f:var:`U`, conjugated, so that
  !  :math:`A = U^\dagger \mathrm{diag}(D) U`. The sweep at which convergence was detected is returned in :f:var:`sweep`
  !
  !The program stops if the matrix is not square or if convergence is not reached within 50 sweeps.
  !
     module procedure d_jacobi
     module procedure c_jacobi
  end interface eigh_jacobi
  !
  interface eigvals
  !This function returns the eigenvalues of a general real or complex square matrix :f:var:`A`, computed with the LAPACK
  !routines :f:func_inline:`dgeev` (:f:func_inline:`deigvals`) and :f:func_inline:`zgeev` (:f:func_inline:`zeigvals`) without
  !computing the eigenvectors. The input matrix is not modified. The result is a complex array of size :code:`n`, in the order
  !given by LAPACK (not sorted). The program stops if :f:var:`A` is not square or if the LAPACK routine returns an error code.
  !
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals
  !
  interface eigvalsh
  !This function returns the eigenvalues of a real symmetric (:f:func_inline:`deigvalsh`) or complex hermitian
  !(:f:func_inline:`zeigvalsh`) matrix :f:var:`A`, computed with the LAPACK routines :f:func_inline:`dsyevd` and
  !:f:func_inline:`zheevd` without computing the eigenvectors. Only the upper triangle of :f:var:`A` is referenced and the input
  !matrix is not modified. The result is a real array, in ascending order. The program stops if :f:var:`A` is not square or if
  !the LAPACK routine returns an error code.
  !
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
  !This function returns the singular values of a real (:f:func_inline:`dsvdvals`) or complex (:f:func_inline:`zsvdvals`)
  !:math:`m\times n` matrix :f:var:`A`, in descending order, computed with the LAPACK routine :f:func_inline:`dgesvd`
  !(:f:func_inline:`zgesvd`) without computing the singular vectors. The input matrix is not modified. The result is a real
  !array allocated with size :code:`min(m,n)`. The program stops if the LAPACK routine returns an error code.
  !
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals
  !
  interface svd
  !This subroutine computes the singular value decomposition :math:`A = U \Sigma V^T` of a real (:f:func_inline:`dsvd`) or
  !:math:`A = U \Sigma V^H` of a complex (:f:func_inline:`zsvd`) :math:`m\times n` matrix :f:var:`A`, with the LAPACK routines
  !:f:func_inline:`dgesvd` and :f:func_inline:`zgesvd`. The input matrix is not modified.
  !
  !The singular values :f:var:`s` are in descending order, and :math:`\Sigma` is the :math:`m\times n` matrix with
  !:math:`\Sigma_{ii} = s_i`. The matrices :f:var:`U` and :f:var:`Vtransp` are full, with shapes :code:`[m,m]` and
  !:code:`[n,n]`. Note that the routine returns the transpose (conjugate transpose in the complex case) of :math:`V`, not
  !:math:`V`. The program stops if :f:var:`U` or :f:var:`Vtransp` do not have these shapes, or if the LAPACK routine returns an
  !error code.
  !
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
  !This subroutine inverts in place a general real (:f:func_inline:`dinv`) or complex (:f:func_inline:`zinv`) square matrix
  !:f:var:`Am` using the LU factorization with partial pivoting, :f:func_inline:`dgetrf` + :f:func_inline:`dgetri`
  !(:f:func_inline:`zgetrf` + :f:func_inline:`zgetri`). On output :f:var:`Am` is replaced by its inverse. The program stops if
  !:f:var:`Am` is not square, if it is singular, or if the LAPACK routines return an error code.
  !
     module procedure dinv
     module procedure zinv
  end interface inv
#ifdef _SCALAPACK
  interface p_inv
  !This subroutine is the distributed-memory counterpart of :f:func_inline:`inv`, for real (:f:func_inline:`p_dinv`) and complex
  !(:f:func_inline:`p_zinv`) square matrices, using ScaLAPACK. It is available only when compiling with :code:`_SCALAPACK`.
  !
  !The matrix :f:var:`A`, a full copy of which is expected on every process, is distributed on the BLACS process grid of the
  !module in a block-cyclic layout with square blocks of size :f:var:`Nblock`, and inverted in place with the LU factorization,
  !:f:func_inline:`PDGETRF` + :f:func_inline:`PDGETRI` (:f:func_inline:`PZGETRF` + :f:func_inline:`PZGETRI`). The inverse is
  !gathered back into :f:var:`A`. The program stops if :f:var:`A` is not square or if a ScaLAPACK routine fails.
  !
     module procedure p_dinv
     module procedure p_zinv
  end interface p_inv
#endif
  !
  interface inv_sym
  !This subroutine inverts in place a real symmetric (:f:func_inline:`dinv_sym`) or complex symmetric
  !(:f:func_inline:`zinv_sym`, NOT hermitian, see :f:func_inline:`inv_her`) matrix :f:var:`A` using the Bunch-Kaufman
  !factorization, :f:func_inline:`dsytrf` + :f:func_inline:`dsytri` (:f:func_inline:`zsytrf` + :f:func_inline:`zsytri`). Only
  !the triangle selected by :f:var:`uplo` (default :code:`'U'`) is referenced on input. On output :f:var:`A` contains the full
  !inverse, the computed triangle being mirrored onto the other one. The program stops if a LAPACK routine returns an error code.
  !
     module procedure dinv_sym
     module procedure zinv_sym
  end interface inv_sym
  !
  interface inv_her
  !This subroutine inverts in place a complex hermitian matrix :f:var:`A` using the Bunch-Kaufman factorization,
  !:f:func_inline:`zhetrf` + :f:func_inline:`zhetri`. The hermiticity of :f:var:`A` is tested first. Only the triangle selected
  !by :f:var:`uplo` (default :code:`'U'`) is referenced by LAPACK, and on output :f:var:`A` contains the full inverse, the other
  !triangle being the complex conjugate of the computed one. The program stops if :f:var:`A` is not hermitian or if a LAPACK
  !routine returns an error code.
  !
     module procedure zinv_her
  end interface inv_her
  !
  interface inv_triang
  !This subroutine inverts in place a real (:f:func_inline:`dinv_triang`) or complex (:f:func_inline:`zinv_triang`) triangular
  !matrix :f:var:`A` using the LAPACK routine :f:func_inline:`dtrtri` (:f:func_inline:`ztrtri`). On output the triangle selected
  !by :f:var:`uplo` (default :code:`'U'`) contains the triangle of the inverse, the other one is not referenced. If
  !:code:`diag='U'` the matrix is assumed unit triangular and its diagonal is not referenced, the default being :code:`'N'`. The
  !program stops if the LAPACK routine returns an error code, e.g. for a singular matrix.
  !
     module procedure dinv_triang
     module procedure zinv_triang
  end interface inv_triang
  !
  interface inv_gj
  !This subroutine inverts in place a real (:f:func_inline:`dinv_gj`) or complex (:f:func_inline:`zinv_gj`) square matrix
  !:f:var:`a` with the Gauss-Jordan elimination with full pivoting. On output :f:var:`a` is replaced by its inverse. The program
  !stops if the matrix is singular.
  !
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
  !This subroutine solves the real or complex linear system :math:`A x = b` using the LU factorization, :f:func_inline:`dgetrf`
  !+ :f:func_inline:`dgetrs` (:f:func_inline:`zgetrf` + :f:func_inline:`zgetrs`). The specific procedures cover a single
  !right-hand side (:f:func_inline:`dsolve_1rhs`, :f:func_inline:`zsolve_1rhs`) and several right-hand sides stored as the
  !columns of :f:var:`b` (:f:func_inline:`dsolve_mrhs`, :f:func_inline:`zsolve_mrhs`).
  !
  !The solution overwrites :f:var:`b`. Note that :f:var:`A` is factorized in place: despite the :code:`intent(in)` it contains
  !the LU factors on exit, so a copy must be passed if the original matrix is needed afterwards. If :f:var:`trans` is
  !:code:`'T'` or :code:`'C'` the transposed system :math:`A^T x = b` is solved instead, the default being :code:`'N'`. The
  !program stops if a LAPACK routine returns an error code, e.g. for a singular matrix.
  !
     module procedure dsolve_1rhs
     module procedure zsolve_1rhs
     module procedure dsolve_Mrhs
     module procedure zsolve_Mrhs
  end interface solve
  !
  interface lstsq
  !This function returns the least-squares solution :math:`x = \mathrm{argmin}_x \| A x - b \|_2` of the real
  !(:f:func_inline:`dlstsq`) or complex (:f:func_inline:`zlstsq`) linear system :math:`A x = b`, with :math:`A` an
  !:math:`m\times n` matrix of possibly non-square shape. It uses the LAPACK routine :f:func_inline:`dgelsy`
  !(:f:func_inline:`zgelsy`), a complete orthogonal factorization with column pivoting, which is suitable also for
  !rank-deficient matrices. The input arrays are not modified, and the solution is allocated with size :code:`n`. The program
  !stops if the LAPACK routine returns an error code.
  !
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
  !This subroutine returns the diagonal elements, or the diagonal blocks, of the inverse of a real or complex (block)
  !tridiagonal matrix, without building the full inverse, using a downward and an upward recursion. The specific procedures cover
  !
  !* a tridiagonal matrix given by its sub-diagonal :f:var:`sub_`, main diagonal :f:var:`diag_` and over-diagonal
  !  :f:var:`over_`, of size :f:var:`N`: the :code:`N` diagonal elements of the inverse are returned in :f:var:`Inv`
  !* a block tridiagonal matrix given by its block diagonals, made of :f:var:`Nb` blocks of size :math:`N\times N`: the
  !  :code:`Nb` diagonal blocks of the inverse are returned in :f:var:`Ainv`
  !* a full tridiagonal matrix :f:var:`Amat`, or a full block tridiagonal matrix :f:var:`Amat` made of :code:`Nb` blocks of size
  !  :math:`N\times N`, treated in place: the diagonals are extracted with :f:func_inline:`get_tridiag` and on output
  !  :f:var:`Amat` is a diagonal (block diagonal) matrix containing the diagonal elements (blocks) of the inverse
  !
  !The tridiagonal structure of the input is not checked, see :f:func_inline:`check_tridiag`. The program stops if a zero pivot
  !is found, i.e. if the matrix is ill-conditioned, or if :f:var:`Amat` is not square.
  !
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
  !This function checks whether a real or complex matrix :f:var:`Amat` is tridiagonal, i.e. whether all the elements outside the
  !main, sub- and over-diagonals are exactly zero. The specific procedures cover a plain matrix and a matrix made of
  !:f:var:`Nblock` blocks of size :math:`N_{size}\times N_{size}`, to be checked for block tridiagonal structure. The result is
  !:code:`.true.` for a (block) tridiagonal matrix. The program stops if a plain matrix is not square.
  !
     module procedure d_check_tridiag
     module procedure c_check_tridiag
     module procedure d_check_tridiag_block
     module procedure c_check_tridiag_block
  end interface check_tridiag
  !
  interface get_tridiag
  !This subroutine extracts the three diagonals of a real or complex tridiagonal matrix :f:var:`Amat`: the sub-diagonal
  !:code:`sub(i) = A(i+1,i)`, the main diagonal :code:`diag(i) = A(i,i)` and, optionally, the over-diagonal
  !:code:`over(i) = A(i,i+1)`. The block version, for a matrix made of :f:var:`Nblock` blocks of size
  !:math:`N_{size}\times N_{size}`, extracts the block diagonals :code:`sub(i,:,:) = A(i+1,i)`, :code:`diag(i,:,:) = A(i,i)` and
  !:code:`over(i,:,:) = A(i,i+1)`. The input matrix is not modified and its structure is not checked, see
  !:f:func_inline:`check_tridiag`. The inverse operation is :f:func_inline:`build_tridiag`.
  !
     module procedure d_get_tridiag
     module procedure c_get_tridiag
     module procedure d_get_tridiag_block
     module procedure c_get_tridiag_block
  end interface get_tridiag
  !
  interface build_tridiag
  !This function builds a real or complex tridiagonal matrix from its three diagonals: the sub-diagonal
  !:code:`A(i+1,i) = sub(i)`, the main diagonal :code:`A(i,i) = diag(i)` and the over-diagonal :code:`A(i,i+1) = over(i)`. If
  !:f:var:`over` is not present the matrix is symmetric, :code:`over = sub`. The block version, with :f:var:`Nblock` blocks of
  !size :math:`N_{size}\times N_{size}`, builds a block tridiagonal matrix from the block diagonals :code:`sub(i,:,:)`,
  !:code:`diag(i,:,:)` and :code:`over(i,:,:)`. The inverse operation is :f:func_inline:`get_tridiag`.
  !
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
  !This function returns the determinant of a real (:f:func_inline:`ddet`) or complex (:f:func_inline:`zdet`) square matrix
  !:f:var:`A`, computed from its LU factorization (:f:func_inline:`dgetrf`, :f:func_inline:`zgetrf`) as the product of the
  !diagonal elements of :math:`U` with the sign changes due to the row permutations. The input matrix is not modified. The
  !program stops if :f:var:`A` is not square, or if the LAPACK routine returns an error code, which includes the case of an
  !exactly singular matrix.
  !
     module procedure ddet
     module procedure zdet
  end interface det
  !
  interface deye
  !This function returns real identity elements: :code:`deye(n)` is the :math:`n\times n` identity matrix
  !(:f:func_inline:`deye_matrix`) and :code:`deye(i,j)` is the Kronecker delta :math:`\delta_{ij}`, 1 if :code:`i==j` and 0
  !otherwise (:f:func_inline:`deye_indices`). It is the same as :f:func_inline:`eye`, and the complex version is
  !:f:func_inline:`zeye`.
  !
     module procedure deye_matrix
     module procedure deye_indices
  end interface deye
  !
  interface zeye
  !This function returns complex identity elements: :code:`zeye(n)` is the :math:`n\times n` identity matrix
  !(:f:func_inline:`zeye_matrix`) and :code:`zeye(i,j)` is the Kronecker delta :math:`\delta_{ij}`, 1 if :code:`i==j` and 0
  !otherwise (:f:func_inline:`zeye_indices`). The real version is :f:func_inline:`deye`.
  !
     module procedure zeye_matrix
     module procedure zeye_indices
  end interface zeye
  !
  interface eye
  !This function returns real identity elements: :code:`eye(n)` is the :math:`n\times n` identity matrix
  !(:f:func_inline:`deye_matrix`) and :code:`eye(i,j)` is the Kronecker delta :math:`\delta_{ij}`, 1 if :code:`i==j` and 0
  !otherwise (:f:func_inline:`deye_indices`). It is the same as :f:func_inline:`deye`. The complex version is
  !:f:func_inline:`zeye`.
  !
     module procedure deye_matrix
     module procedure deye_indices
  end interface eye
  !
  interface diag
  !This function builds a real (:f:func_inline:`ddiag`) or complex (:f:func_inline:`zdiag`) square matrix with the elements of
  !the array :f:var:`x` on the main diagonal and zeros elsewhere. The inverse operation is :f:func_inline:`diagonal`.
  !
     module procedure ddiag
     module procedure zdiag
  end interface diag
  !
  interface diagonal
  !This function returns the main diagonal of a real (:f:func_inline:`d_diagonal`) or complex (:f:func_inline:`z_diagonal`)
  !matrix :f:var:`A`, expected to be square: the size of the result is :code:`size(A,1)`. The inverse operation is
  !:f:func_inline:`diag`.
  !
     module procedure d_diagonal
     module procedure z_diagonal
  end interface diagonal
  !
  interface trace
  !This function returns the trace of a real (:f:func_inline:`dtrace`) or complex (:f:func_inline:`ztrace`) matrix :f:var:`A`,
  !i.e. the sum of the elements along the main diagonal. For a non-square matrix the sum runs up to :code:`min(m,n)`.
  !
     module procedure dtrace
     module procedure ztrace
  end interface trace
  !
  interface zeros
  !This function returns a complex(8) array with all the elements set to 0. The shape is given by one integer extent for each
  !dimension, :code:`zeros(n1,...,nk)`, and the ranks from 1 to 7 are covered.
  !
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
  !This function returns a complex(8) array with all the elements set to 1. The shape is given by one integer extent for each
  !dimension, :code:`ones(n1,...,nk)`, and the ranks from 1 to 7 are covered.
  !
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
  !This function returns the Kronecker (tensor) product :math:`A \otimes B` of two matrices :f:var:`A` and :f:var:`B`, i.e. the
  !block matrix whose :math:`(i,j)` block is :math:`A_{ij} B`. The result has shape
  !:code:`[size(A,1)*size(B,1),size(A,2)*size(B,2)]`. The specific procedures cover integer, real, real :math:`\otimes` complex,
  !complex :math:`\otimes` real and complex matrices, the result being complex if any of the two is complex. It is the same as
  !:f:func_inline:`kronecker_product` and as the operator :code:`.kx.`.
  !
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface kron
  !
  interface kronecker_product
  !This function returns the Kronecker (tensor) product :math:`A \otimes B` of two matrices :f:var:`A` and :f:var:`B`, i.e. the
  !block matrix whose :math:`(i,j)` block is :math:`A_{ij} B`. The result has shape
  !:code:`[size(A,1)*size(B,1),size(A,2)*size(B,2)]`. The specific procedures cover integer, real, real :math:`\otimes` complex,
  !complex :math:`\otimes` real and complex matrices, the result being complex if any of the two is complex. It is the same as
  !:f:func_inline:`kron` and as the operator :code:`.kx.`.
  !
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface kronecker_product
  !
  interface operator(.kx.)
  !This operator returns the Kronecker (tensor) product :math:`A \otimes B` of two matrices :f:var:`A` and :f:var:`B`, i.e. the
  !block matrix whose :math:`(i,j)` block is :math:`A_{ij} B`. The result has shape
  !:code:`[size(A,1)*size(B,1),size(A,2)*size(B,2)]`. The specific procedures cover integer, real, real :math:`\otimes` complex,
  !complex :math:`\otimes` real and complex matrices, the result being complex if any of the two is complex. It is used as
  !:code:`AxB = A .kx. B`, and it is the same as :f:func_inline:`kron`.
  !
     module procedure :: i_kronecker_product
     module procedure :: d_kronecker_product
     module procedure :: dc_kronecker_product
     module procedure :: cd_kronecker_product
     module procedure :: c_kronecker_product
  end interface operator(.kx.)
  !
  interface outerprod
  !This function returns the outer product of two real or complex one-dimensional arrays :f:var:`a` and :f:var:`b`, i.e. the
  !matrix :math:`A_{ij} = a_i b_j`, of shape :code:`[size(a),size(b)]`. No complex conjugation is applied.
  !
     module procedure outerprod_d,outerprod_c
  end interface outerprod
  !
  interface cross_product
  !This function returns the cross (vector) product :math:`c = a \times b` of two real or complex vectors :f:var:`a` and
  !:f:var:`b` with 3 components, :math:`c_1 = a_2 b_3 - a_3 b_2`, :math:`c_2 = a_3 b_1 - a_1 b_3`,
  !:math:`c_3 = a_1 b_2 - a_2 b_1`. Only the 3-dimensional version is available through this interface.
  !
     module procedure :: cross_3d_d
     module procedure :: cross_3d_c
  end interface cross_product
  !
  interface s3_product
  !This function returns the triple product :math:`a \cdot (b \times c)` of three real or complex vectors :f:var:`a`, :f:var:`b`
  !and :f:var:`c` with 3 components, evaluated as :code:`dot_product(a,cross_product(b,c))`. The result is a real number: in the
  !complex case :f:var:`a` is conjugated by :code:`dot_product` and the imaginary part of the result is discarded.
  !
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
  !This subroutine computes the real (:f:func_inline:`d_matmul`) or complex (:f:func_inline:`z_matmul`) matrix-matrix product
  !:math:`C = \alpha A B + \beta C` with the BLAS routine :f:func_inline:`dgemm` (:f:func_inline:`zgemm`), :math:`A` being
  !:math:`N\times K`, :math:`B` :math:`K\times M` and :math:`C` :math:`N\times M`. By default :math:`\alpha=1` and
  !:math:`\beta=0`, so that :f:var:`C` contains the plain product :math:`A B` on output. :f:var:`A` and :f:var:`B` are not
  !modified. The program stops if :f:var:`B` or :f:var:`C` do not have the shapes :code:`[K,M]` and :code:`[N,M]`. The function
  !form is the operator :code:`.x.`.
  !
     module procedure :: d_matmul
     module procedure :: z_matmul
  end interface mat_product
  interface operator(.x.)
  !This operator returns the real or complex matrix-matrix product :math:`C = A B`, computed with the BLAS routine
  !:f:func_inline:`dgemm` (:f:func_inline:`zgemm`), :math:`A` being :math:`N\times K` and :math:`B` :math:`K\times M`. It is
  !used as :code:`C = A .x. B`, through the functions :f:func_inline:`d_matmul_` and :f:func_inline:`z_matmul_`. The program
  !stops if :f:var:`B` does not have the shape :code:`[K,M]`. The subroutine form is :f:func_inline:`mat_product`.
  !
     module procedure :: d_matmul_
     module procedure :: z_matmul_
  end interface operator(.x.)
#ifdef _SCALAPACK
  interface p_mat_product
  !This subroutine is the distributed-memory counterpart of :f:func_inline:`mat_product`: it computes the real
  !(:f:func_inline:`p_d_matmul`) or complex (:f:func_inline:`p_z_matmul`) matrix-matrix product :math:`C = \alpha A B + \beta C`
  !with the ScaLAPACK routine :f:func_inline:`PDGEMM` (:f:func_inline:`PZGEMM`). It is available only when compiling with
  !:code:`_SCALAPACK`.
  !
  !The matrices :f:var:`A` and :f:var:`B`, full copies of which are expected on every process, are distributed on the BLACS
  !process grid of the module in a block-cyclic layout with square blocks of size :f:var:`Nblock`, and the product is gathered
  !back into :f:var:`C`. By default :math:`\alpha=1` and :math:`\beta=0`. The program stops if :f:var:`B` or :f:var:`C` do not
  !have the shapes :code:`[K,M]` and :code:`[N,M]`. The function form is the operator :code:`.Px.`.
  !
     module procedure :: p_d_matmul
     module procedure :: p_z_matmul
  end interface p_mat_product
  interface operator(.Px.)
  !This operator returns the real or complex matrix-matrix product :math:`C = A B`, computed in parallel with the ScaLAPACK
  !routine :f:func_inline:`PDGEMM` (:f:func_inline:`PZGEMM`), :math:`A` being :math:`N\times K` and :math:`B` :math:`K\times M`.
  !It is used as :code:`C = A .Px. B`, through the functions :f:func_inline:`p_d_matmul_f` and :f:func_inline:`p_z_matmul_f`,
  !and it is available only when compiling with :code:`_SCALAPACK`. The block size of the block-cyclic distribution is chosen
  !automatically, as the largest power of 2 smaller than :code:`min(N,K,M)`, up to 64. The program stops if :f:var:`B` does not
  !have the shape :code:`[K,M]`. The subroutine form is :f:func_inline:`p_mat_product`.
  !
     module procedure :: p_d_matmul_f
     module procedure :: p_z_matmul_f
  end interface operator(.Px.)
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
  !This function returns :code:`Nblock` real :math:`N\times N` identity matrices stacked in a rank-3 array, i.e. the
  !blocks of a block identity matrix. It can be used to build block tridiagonal matrices, see
  !:f:func_inline:`build_tridiag`.
  !
    integer                       :: Nblock     ! number of blocks
    integer                       :: N          ! size of each block
    real(8),dimension(Nblock,N,N) :: eye_block  ! identity blocks, eye_block(i,:,:) = 1, [Nblock,N,N]
    integer                       :: iblock
    do iblock=1,Nblock
       eye_block(iblock,:,:) = eye(N)
    enddo
  end function deye_tridiag
  !
  function zeye_tridiag(Nblock,N) result(eye_block)
  !This function is the complex counterpart of :f:func_inline:`deye_tridiag`: it returns :code:`Nblock` complex
  !:math:`N\times N` identity matrices stacked in a rank-3 array.
  !
    integer                          :: Nblock     ! number of blocks
    integer                          :: N          ! size of each block
    complex(8),dimension(Nblock,N,N) :: eye_block  ! identity blocks, eye_block(i,:,:) = 1, [Nblock,N,N]
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
