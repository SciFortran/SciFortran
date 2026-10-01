subroutine dinv(Am)
  !This subroutine inverts in place a general real square matrix :f:var:`Am` using the LU factorization with partial
  !pivoting, :f:func_inline:`dgetrf` + :f:func_inline:`dgetri`. The complex version :f:func_inline:`zinv` uses
  !:f:func_inline:`zgetrf` + :f:func_inline:`zgetri`. Both are instances of the generic interface :f:func_inline:`inv`.
  !
  !On output :f:var:`Am` is replaced by its inverse. The program stops if :f:var:`Am` is not square, if it is singular,
  !or if the LAPACK routines return an error code.
  !
  real(8), intent(inout) :: Am(:,:)  ! in: square matrix to be inverted; out: its inverse
  real(8), allocatable   :: Amt(:,:),work(:)  ! temporary work arrays
  integer                :: info,lda,n,lwork,nb
  integer, allocatable   :: ipiv(:)
  n = size(Am(1, :))
  call assert_shape(Am, [n, n], "inv", "Am")
  lda = n
  nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
  lwork = n*nb
  if (nb < 1) nb = max(1, n)
  allocate(Amt(n,n), work(lwork), ipiv(n))
  Amt = Am
  call dgetrf(n, n, Amt, lda, ipiv, info)
  if(info /= 0) then
     print *, "dgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop ' D_mat_invert: dgetrf error'
  end if
  call dgetri(n, Amt, n, ipiv, work, lwork, info)
  if (info /= 0) then
     print *, "dgetri returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
        print *, "singular and its inverse could not be computed."
     end if
     stop ' D_mat_invert: dgetri error'
  end if
  Am = Amt
end subroutine Dinv
!
subroutine Zinv(Am)
  ! Inverts the general complex matrix Am
  complex(8), intent(inout) :: Am(:,:)   ! Matrix to be inverted
  complex(8), allocatable   :: Amt(:,:), work(:)
  integer                   :: n, nb
  integer                   :: lwork, info
  integer, allocatable      :: ipiv(:)
  n = size(Am, 1)
  call assert_shape(Am, [n, n], "inv", "Am")
  nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
  if (nb < 1) nb = max(1, n)
  lwork = n*nb
  allocate(Amt(n,n), ipiv(n), work(lwork))
  Amt = Am
  call zgetrf(n, n, Amt, n, ipiv, info)
  if (info /= 0) then
     print *, "zgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'Z_mat_invert: zgetrf error'
  end if
  call zgetri(n, Amt, n, ipiv, work, lwork, info)
  if (info /= 0) then
     print *, "zgetri returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
        print *, "singular and its inverse could not be computed."
     end if
     stop 'Z_mat_invert: zgetri error'
  end if
  Am = Amt
end subroutine Zinv
