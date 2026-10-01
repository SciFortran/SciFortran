function dsvdvals(A) result(s)
  !This function returns the singular values of a real :math:`m\times n` matrix :f:var:`A`, in descending order,
  !computed with the LAPACK routine :f:func_inline:`dgesvd` without singular vectors. The complex version
  !:f:func_inline:`zsvdvals` uses :f:func_inline:`zgesvd`. Both are instances of the generic interface
  !:f:func_inline:`svdvals`.
  !
  !The input matrix is copied and left untouched. The result is allocated with size :code:`min(m,n)`; the program
  !stops if the LAPACK routine returns an error code.
  !
  real(8), intent(in)  :: A(:,:)  ! real matrix [m,n], not modified
  real(8), allocatable :: s(:)    ! singular values, descending, size min(m,n)
  integer              :: info, lwork, m, n
  real(8), allocatable :: work(:), At(:,:)
  real(8)              :: u(1,1), vt(1,1)  ! not used if only s is to be computed
  m = size(A(:,1))  ! = lda
  n = size(A(1,:))
  allocate(At(m,n), s(min(m,n)))
  At(:,:) = A(:, :)  ! A is overwritten in dgesvd
  ! query optimal lwork and allocate workspace:
  allocate(work(1))
  call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, info)
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, info)
  if(info /= 0) then
     print *, "dgesvd returned info = ", info
     if(info < 0) then
        print *, "the ", -info, "-th argument had an illegal value"
     else
        print *, "DBDSQR did not converge, there are ", info
        print *, "superdiagonals of an intermediate bidiagonal form B"
        print *, "did not converge to zero. See the description of WORK"
        print *, "in DGESVD's man page for details."
     endif
     stop 'dsvdvals error: dgesvd '
  endif
end function dsvdvals

function zsvdvals(A) result(s)
  complex(8), intent(in)  :: A(:,:)
  real(8), allocatable    :: s(:)
  integer                 :: info, lwork, m, n, lrwork
  complex(8), allocatable :: work(:), At(:,:)
  real(8), allocatable    :: rwork(:)
  complex(8)              :: u(1,1), vt(1,1)  ! not used if only s is to be computed
  m = size(A(:,1))  ! = lda
  n = size(A(1,:))
  lrwork = 5*min(m,n)
  allocate(At(m,n), s(min(m,n)), rwork(lrwork))
  At(:,:) = A(:,:)  ! A is overwritten in zgesvd!
  ! query optimal lwork and allocate workspace:
  allocate(work(1))
  call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, rwork, info)
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, rwork, info)
  if(info /= 0) then
     print *, "zgesvd returned info = ", info
     if(info < 0) then
        print *, "the ", -info, "-th argument had an illegal value"
     else
        print *, "ZBDSQR did not converge, there are ", info
        print *, "superdiagonals of an intermediate bidiagonal form B"
        print *, "did not converge to zero. See the description of RWORK"
        print *, "in ZGESVD's man page for details."
     endif
     stop 'zsvdvals error: zgesvd'
  endif
end function zsvdvals
