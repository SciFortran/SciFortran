function deigvals(A) result(lam)
  !This function returns the eigenvalues of a general real matrix :f:var:`A`, computed with the LAPACK routine
  !:f:func_inline:`dgeev` without eigenvectors. The complex version :f:func_inline:`zeigvals` uses
  !:f:func_inline:`zgeev`. Both are instances of the generic interface :f:func_inline:`eigvals`.
  !
  !The input matrix is copied and left untouched. The eigenvalues are complex and are returned in the order given
  !by LAPACK (not sorted).
  !
  !The program stops if :f:var:`A` is not square, or if the LAPACK routine returns an error code.
  !
  real(8), intent(in)     :: A(:, :)         ! real square matrix [n,n], not modified
  complex(8)              :: lam(size(A,1))  ! eigenvalues of A, A c = lam c, size n
  real(8), allocatable    :: At(:,:),vl(:,: ),vr(:,:),wi(:),work(:),wr(:)
  real(8)                 :: lwork_guess(1)
  integer                 :: info, lda, ldvl, ldvr, lwork, n
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  ldvl  = n
  ldvr  = n
  allocate(At(n,n),wr(n),wi(n),vl(ldvl,n),vr(ldvr,n))
  !Copy the input Matrix
  At    = A
  !1st Call: Query the right size for the working array.
  call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr,  lwork_guess, -1, info)
  if(info /= 0) then
     print*, "dgeev returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*, "the QR algorithm failed to compute all the"
        print*, "eigenvalues, and no eigenvectors have been computed;"
        print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
        print*, "have converged."
     end if
     stop 'deig error: 1st call dgeev'
  endif
  lwork = int(lwork_guess(1)) ;allocate(work(lwork))
  !2nd Call: Actual solution of the eigenproblem
  call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
  if(info /= 0) then
     print *, "dgeev returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*, "the QR algorithm failed to compute all the"
        print*, "eigenvalues, and no eigenvectors have been computed;"
        print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
        print*, "have converged."
     end if
     stop 'deig error: 2nd call dgeev'
  end if
  lam = wr + xi*wi
  deallocate(At,wr,wi,vl,vr)
end function deigvals

function zeigvals(A) result(lam)
  complex(8), intent(in)  :: A(:, :)  ! matrix to solve eigenproblem for
  complex(8)              :: lam(size(A,1))  ! eigenvalues: A c = lam c
  integer                 :: info, lda, ldvl, ldvr, lwork, n, lrwork
  real(8), allocatable    :: rwork(:)
  complex(8), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)
  complex(8)              :: lwork_guess(1)
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  ldvl  = n
  ldvr  = n
  lrwork= 2*n
  allocate(At(n,n),vl(ldvl,n),vr(ldvr,n),rwork(lrwork))
  !Copy the input Matrix
  At  = A
  !1st Call: Query the right size for the working array.
  call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, lwork_guess, -1, rwork, info)
  if(info /= 0) then
     print*, "zgeev returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*, "the QR algorithm failed to compute all the"
        print*, "eigenvalues, and no eigenvectors have been computed;"
        print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
        print*, "have converged."
     end if
     stop 'deig error: 1st call zgeev'
  endif
  lwork = int(lwork_guess(1)) ;allocate(work(lwork))
  !2nd Call: Actual solution of the eigenproblem
  call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
  if(info /= 0) then
     print *, "zgeev returned info = ",info
     if (info < 0) then
        print*, "the",-info,"-th argument had an illegal value"
     else
        print*, "the QR algorithm failed to compute all the"
        print*, "eigenvalues, and no eigenvectors have been computed;"
        print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
        print*, "have converged."
     end if
     stop 'deig error: 2nd call zgeev'
  end if
  deallocate(At,vl,vr,rwork)
end function zeigvals
