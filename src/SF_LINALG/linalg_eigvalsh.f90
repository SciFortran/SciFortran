function deigvalsh(A) result(lam)
  !This function returns the eigenvalues of a real symmetric matrix :f:var:`A`, computed with the LAPACK routine
  !:f:func_inline:`dsyevd` without eigenvectors. The complex version :f:func_inline:`zeigvalsh` works on hermitian
  !matrices with :f:func_inline:`zheevd`. Both are instances of the generic interface :f:func_inline:`eigvalsh`.
  !
  !Only the upper triangle of :f:var:`A` is referenced; the input matrix is copied and left untouched. The eigenvalues
  !are real (also for hermitian matrices) and are returned in ascending order.
  !
  !The program stops if :f:var:`A` is not square, or if the LAPACK routine returns an error code.
  !
  real(8), intent(in)  :: A(:, :)         ! real symmetric matrix [n,n], not modified
  real(8)              :: lam(size(A,1))  ! eigenvalues of A, ascending, size n
  real(8), allocatable :: At(:,:),work(:),iwork(:)
  real(8)              :: lwork_guess(1),liwork_guess(1)
  integer              :: info,lda,lwork,liwork,n
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  allocate(At(n,n))
  !Copy the input Matrix
  At    = A
  !1st Call: Query the right size for the working array.
  call dsyevd('N','U', n, At, lda, lam, lwork_guess, -1, liwork_guess, -1, info )
  if(info /= 0) then
     print*, "dsyevd returned info = ",info
     stop
  endif
  lwork  = int(lwork_guess(1))  ;allocate(work(lwork))
  liwork = int(liwork_guess(1)) ;allocate(iwork(liwork))
  !2nd Call: Actual solution of the eigenproblem
  call dsyevd('N','U', n, At, lda, lam, work, lwork, iwork, liwork, info )
  if(info /= 0) then
     print *, "ssyevd returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*,"the algorithm failed to converge"
        print*,"the",info,"off-diagonal elements of an intermediate"
        print*,"tridiagonal form did not converge to zero"
     end if
     stop 'deigvalsh error: 2nd call dsyevd'
  end if
end function deigvalsh

function zeigvalsh(A) result(lam)
  complex(8),intent(in) :: A(:, :)         ! matrix for eigenvalue compuation
  real(8)               :: lam(size(A,1))  ! eigenvalues: A c = lam c
  real(8), allocatable  :: At(:,:),work(:),rwork(:),iwork(:)
  real(8)               :: lwork_guess(1),lrwork_guess(1),liwork_guess(1)
  integer               :: info, lda,lwork,lrwork,liwork,n
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  allocate(At(n,n))
  !Copy the input Matrix
  At    = A
  !1st Call: Query the right size for the working array.
  call zheevd('N','U', n, At, lda, lam, lwork_guess,-1, lrwork_guess,-1, liwork_guess,-1, info )
  if(info /= 0) then
     print*, "zsyevd returned info = ",info
     stop
  endif
  lwork  = int(lwork_guess(1))  ;allocate(work(lwork))
  rwork  = int(lrwork_guess(1)) ;allocate(rwork(lrwork))
  liwork = int(liwork_guess(1)) ;allocate(iwork(liwork))
  !2nd Call: Actual solution of the eigenproblem
  call zheevd('N','U', n, At, lda, lam, work,lwork, rwork,lrwork, iwork,liwork, info )
  if(info /= 0) then
     print *, "zheevd returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*,"the algorithm failed to converge"
        print*,"the",info,"off-diagonal elements of an intermediate"
        print*,"tridiagonal form did not converge to zero"
     end if
     stop 'zeigvalsh error: 2nd call zheevd'
  end if
end function zeigvalsh
