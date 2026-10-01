!<TODO  add optional switch for left or right eigenvectors in deig() and zeig()?
!>TODO
subroutine deig(A,Eval,Evec,jobvl,jobvr)
  !This subroutine computes the eigenvalues and the right eigenvectors of a general (non-symmetric) real matrix
  !:f:var:`A`, :math:`A v_j = \lambda_j v_j`, using the LAPACK routine :f:func_inline:`dgeev`. The complex version
  !:f:func_inline:`zeig` uses :f:func_inline:`zgeev`. Both are instances of the generic interface :f:func_inline:`eig`.
  !
  !The input matrix is copied and left untouched. The eigenvalues are returned as complex numbers, in the order
  !given by LAPACK, and the eigenvectors are stored column-wise, i.e. :f:var:`Evec(i,j)` is the :math:`i`-th component
  !of the eigenvector associated to :f:var:`Eval(j)`. For real input the eigenvectors of complex-conjugate pairs of
  !eigenvalues are rebuilt from the two real columns returned by :f:func_inline:`dgeev`.
  !
  !The left eigenvectors are not available. The program stops if :code:`jobvl='V'`, if :f:var:`A` is not square,
  !if :f:var:`Evec` does not have shape :code:`[n,n]`, or if the LAPACK routine returns an error code.
  !
  real(8),intent(in)                 :: A(:,:)     ! real matrix [n,n], not modified
  complex(8),intent(out)             :: Eval(:)    ! eigenvalues, size n
  complex(8),intent(out)             :: Evec(:,:)  ! right eigenvectors, [n,n]; Evec(:,j) <-> Eval(j)
  character(len=1),optional          :: jobvl      ! optional, 'N' (default); 'V' is not supported
  character(len=1),optional          :: jobvr      ! optional, 'V' (default) computes the right eigenvectors, 'N' does not
  character(len=1)                   :: jobvl_,jobvr_
  real(8),dimension(:,:),allocatable :: At,vl,vr
  real(8),dimension(:),allocatable   :: wi,wr
  real(8),dimension(:),allocatable   :: work
  real(8),dimension(1)               :: lwork_guess
  integer                            :: info,lda,ldvl,ldvr,lwork,n,i
  jobvl_='N';if(present(jobvl))jobvl_=jobvl
  jobvr_='V';if(present(jobvr))jobvr_=jobvr
  if(jobvl_=='V')stop "deig error: jobvl = V is not supported yet."
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  call assert_shape(Evec,[n,n],"solve","Evec")
  ldvl  = n
  ldvr  = n
  allocate(At(n,n),wr(n),wi(n),vl(ldvl,n),vr(ldvr,n))
  !Copy the input Matrix
  At    = A
  !1st Call: Query the right size for the working array.
  call dgeev(jobvl_, jobvr_, n, At, lda, wr, wi, vl, ldvl, vr, ldvr, lwork_guess, -1, info)
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
  call dgeev(jobvl_, jobvr_, n, At, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
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
  Eval = wr + xi*wi
  ! as DGEEV has a rather complicated way of returning the eigenvectors,
  ! it is necessary to build the complex array of eigenvectors from
  ! two real arrays:
  do i = 1,n
     if(wi(i) > 0d0) then  ! first of two conjugate eigenvalues
        Evec(:, i) = vr(:, i) + xi*vr(:, i+1)
     elseif(wi(i) < 0d0) then  ! second of two conjugate eigenvalues
        Evec(:, i) = vr(:, i-1) - xi*vr(:, i)
     else
        Evec(:, i) = vr(:, i)
     end if
  end do
end subroutine deig

subroutine zeig(A,Eval,Evec,jobvl,jobvr)
  !IN matrix - OUT eigenvectors: A.evec = eval.evec ; evec(i,j) = ith component of jth vec.
  complex(8),intent(in)                 :: A(:,:)
  complex(8),intent(out)                :: Eval(:)
  complex(8),intent(out)                :: Evec(:,:)
  character(len=1),optional             :: jobvl,jobvr
  character(len=1)                      :: jobvl_,jobvr_
  real(8),dimension(:),allocatable      :: rwork
  complex(8),dimension(:),allocatable   :: work
  complex(8),dimension(1)               :: lwork_guess
  complex(8),dimension(:,:),allocatable :: vl,vr
  integer                               :: info,lda,ldvl,ldvr,lwork,n,i,lrwork

  jobvl_='N';if(present(jobvl))jobvl_=jobvl
  jobvr_='V';if(present(jobvr))jobvr_=jobvr
  if(jobvl_=='V')stop "deig error: jobvl = V is not supported yet."
  lda   = size(A(:,1))
  n     = size(A(1,:))
  call assert_shape(A,[n,n],"solve","A")
  call assert_shape(Evec,[n,n],"solve","Evec")
  ldvl  = n
  ldvr  = n
  lrwork= 2*n
  allocate(vl(ldvl,n),vr(ldvr,n),rwork(lrwork))
  !Copy the input Matrix
  Evec  = A
  !1st Call: Query the right size for the working array.
  call zgeev(jobvl_, jobvr_, n, Evec, lda, Eval, vl, ldvl, vr, ldvr, lwork_guess, -1, rwork, info)
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
  call zgeev(jobvl_, jobvr_, n, Evec, lda, Eval, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
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
  Evec = vr
end subroutine zeig
