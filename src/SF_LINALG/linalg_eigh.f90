!-------------------------------------------------------------------------------------------
!PURPOSE:  eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
!-------------------------------------------------------------------------------------------
subroutine deigh_generalized(Am, Bm, lam, c)
  !This subroutine solves the generalized eigenvalue problem :math:`A c = \lambda B c` for all eigenvalues and
  !eigenvectors, where :f:var:`Am` is real symmetric and :f:var:`Bm` is real symmetric positive definite, using the
  !LAPACK routine :f:func_inline:`dsygvd`. The complex version :f:func_inline:`zeigh_generalized` (hermitian :f:var:`Am`
  !and :f:var:`Bm`, real eigenvalues) uses :f:func_inline:`zhegvd`. Both are instances of the generic interface
  !:f:func_inline:`eigh`.
  !
  !Only the lower triangular parts of :f:var:`Am` and :f:var:`Bm` are referenced; both matrices are copied and left
  !untouched. The eigenvalues are returned in ascending order and the eigenvectors are stored column-wise,
  !normalized such that :math:`c^T B c = 1` (:math:`c^H B c = 1` in the complex case).
  !
  !The program stops if :f:var:`Am`, :f:var:`Bm` or :f:var:`c` do not have shape :code:`[n,n]`, if :f:var:`Bm` is not
  !positive definite, or if the LAPACK routine fails to converge.
  !
  real(8), intent(in)  :: Am(:,:)  ! LHS matrix of Am c = lam Bm c, [n,n]
  real(8), intent(in)  :: Bm(:,:)  ! RHS matrix of Am c = lam Bm c, positive definite, [n,n]
  real(8), intent(out) :: lam(:)   ! eigenvalues, ascending, size n
  real(8), intent(out) :: c(:,:)   ! eigenvectors, c(i,j) = i-th component of j-th vector, [n,n]
  integer              :: n
  ! lapack variables
  integer              :: lwork, liwork, info
  integer, allocatable :: iwork(:)
  real(8), allocatable :: Bmt(:,:), work(:)
  ! solve
  n = size(Am,1)
  call assert_shape(Am, [n, n], "eigh", "Am")
  call assert_shape(Bm, [n, n], "eigh", "B")
  call assert_shape(c, [n, n], "eigh", "c")
  lwork = 1 + 6*n + 2*n**2
  liwork = 3 + 5*n
  allocate(Bmt(n,n), work(lwork), iwork(liwork))
  c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
  call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
  if (info /= 0) then
     print *, "dsygvd returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else if (info <= n) then
        print *, "the algorithm failed to compute an eigenvalue while working"
        print *, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
        print *, "through", mod(info, n+1)
     else
        print *, "The leading minor of order ", info-n, &
             "of B is not positive definite. The factorization of B could ", &
             "not be completed and no eigenvalues or eigenvectors were computed."
     end if
     stop 'deigh_generalized error: dsygvd'
  end if
end subroutine deigh_generalized
!
subroutine zeigh_generalized(Am, Bm, lam, c)
  ! solves generalized eigen value problem for all eigenvalues and eigenvectors
  ! Am must by hermitian, Bm hermitian positive definite.
  ! Only the lower triangular part of Am and Bm is used.
  complex(8), intent(in)  :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
  complex(8), intent(in)  :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
  real(8), intent(out)    :: lam(:)      ! eigenvalues: Am c = lam Bm c
  complex(8), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
  ! lapack variables
  integer                 :: info, liwork, lrwork, lwork, n
  integer, allocatable    :: iwork(:)
  real(8), allocatable    :: rwork(:)
  complex(8), allocatable :: Bmt(:,:), work(:)
  n = size(Am,1)
  call assert_shape(Am, [n, n], "eigh", "Am")
  call assert_shape(Bm, [n, n], "eigh", "Bm")
  call assert_shape(c, [n, n], "eigh", "c")
  lwork = 2*n + n**2
  lrwork = 1 + 5*N + 2*n**2
  liwork = 3 + 5*n
  allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
  c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd
  call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
  if (info /= 0) then
     print *, "zhegvd returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else if (info <= n) then
        print *, "the algorithm failed to compute an eigenvalue while working"
        print *, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
        print *, "through", mod(info, n+1)
     else
        print *, "The leading minor of order ", info-n, &
             "of B is not positive definite. The factorization of B could ", &
             "not be completed and no eigenvalues or eigenvectors were computed."
     end if
     stop 'deigh_generalized error: zhegvd'
  end if
end subroutine zeigh_generalized







subroutine deigh_simple(A,E,method,jobz,uplo,vl,vu,il,iu,tol)
  !This subroutine computes the eigenvalues and, optionally, the eigenvectors of a real symmetric matrix :f:var:`A`,
  !using the LAPACK driver selected by :f:var:`method`: :f:func_inline:`dsyevd` (default), :f:func_inline:`dsyevr`,
  !:f:func_inline:`dsyev` or :f:func_inline:`dsyevx`. An unrecognized method name falls back to :f:func_inline:`dsyevd`.
  !The complex version :f:func_inline:`zeigh_simple` works on hermitian matrices with :f:func_inline:`zheevd` (default),
  !:f:func_inline:`zheevr`, :f:func_inline:`zheev` and :f:func_inline:`zheevx`; the eigenvalues are real. Both are
  !instances of the generic interface :f:func_inline:`eigh`.
  !
  !The matrix is diagonalized in place: if :code:`jobz='V'` the first :code:`mE` columns of :f:var:`A` contain on
  !output the eigenvectors, stored column-wise, otherwise :f:var:`A` is destroyed. The eigenvalues are returned in
  !:f:var:`E` in ascending order. With the :code:`r` and :code:`x` drivers a subset of the spectrum can be selected
  !either by value, passing both :f:var:`vl` and :f:var:`vu` (eigenvalues in :math:`(v_l,v_u]`), or by index,
  !passing :f:var:`il` and/or :f:var:`iu` (a missing one defaults to 1); the other drivers ignore these arguments
  !and compute the full spectrum.
  !
  !The program stops if both a value range and an index range are given, if :f:var:`jobz` or :f:var:`uplo` have an
  !illegal value, or if :f:var:`A` is not square.
  !
  real(8),dimension(:,:),intent(inout)       :: A       ! in: symmetric matrix; out: eigenvectors (columns) if jobz=V
  real(8),dimension(size(A,2)),intent(inout) :: E       ! eigenvalues, ascending, size n (zero beyond the mE found)
  character(len=*),optional                  :: method  ! optional: 'dsyevd' (default), 'dsyevr', 'dsyev', 'dsyevx'
  character(len=1),optional                  :: jobz    ! optional: 'V' (default) eigenvectors too, 'N' eigenvalues only
  character(len=1),optional                  :: uplo    ! optional: 'U' (default) or 'L', triangle of A to be used
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl      ! optional: lower bound of the eigenvalue interval (vl,vu]
  real(8),optional                           :: vu      ! optional: upper bound of the eigenvalue interval (vl,vu]
  real(8),optional                           :: tol     ! optional: absolute eigenvalue tolerance (dsyevr/dsyevx), default 0
  integer,optional                           :: il      ! optional: index of the first eigenvalue wanted
  integer,optional                           :: iu      ! optional: index of the last eigenvalue wanted
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: i,j,info    
  integer                                    :: lwork,liwork,mE
  integer                                    :: guess_liwork(1)
  real(8)                                    :: guess_lwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  real(8),dimension(:),allocatable           :: work,rwork
  integer,dimension(:),allocatable           :: Isuppz,Ifail
  integer,dimension(:),allocatable           :: Iwork
  !
  method_='dsyevd';if(present(method))method_=trim(method)
  jobz_='V'  ;if(present(jobz))jobz_=jobz
  uplo_='U'  ;if(present(uplo))uplo_=uplo
  vl_  = 1d0 ;if(present(vL))vL_=vL
  vu_  = 1d0 ;if(present(vU))vU_=vU
  iL_  = 1   ;if(present(iL))iL_=iL
  iU_  = 1   ;if(present(iU))iU_=iU
  tol_ = 0d0 ;if(present(tol))tol_=tol
  !
  E=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  if(jobz_/='V'.AND.jobz_/='N')stop "deigh_simple error: jobz has illegal value"
  if(uplo_/='U'.AND.uplo_/='L')stop "deigh_simple error: uplo has illegal value"
  !
  Ns = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "deigh_simple error: A has illegal shape"
  !
  select case(method_)
  case ("dsyevr")
     allocate(Isuppz(2*Ns))
     allocate(Z(Ns,Ns))
     call dsyevr(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,Isuppz,guess_lwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call dsyevr(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,Isuppz,work,lwork,iwork,liwork,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
     !
  case ("dsyev")
     call dsyev(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call dsyev(jobz_,uplo_,Ns,A,Ns,E,work,lwork,info)
     !
  case default
     ! case("dsyevd")
     call dsyevd(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call dsyevd(jobz_,uplo_,Ns,A,Ns,E,work,lwork,iwork,liwork,info)
     !
  case ("dsyevx")
     allocate(Z(Ns,Ns))
     allocate(Ifail(Ns))
     allocate(rwork(7*Ns))
     allocate(iwork(5*Ns))
     call dsyevx(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,guess_lwork,-1,rwork,iwork,ifail,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call dsyevx(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,work,lwork,rwork,iwork,ifail,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
  end select
  return
end subroutine deigh_simple


subroutine zeigh_simple(A,E,method,jobz,uplo,vl,vu,il,iu,tol)
  complex(8),dimension(:,:),intent(inout)       :: A ! M v = E v/v(i,j) = ith cmpt of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: E ! eigenvalues
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: i,j,info    
  integer                                    :: lwork,liwork,lrwork,mE
  complex(8)                                 :: guess_lwork(1)
  real(8)                                    :: guess_lrwork(1)
  integer                                    :: guess_liwork(1)
  logical                                    :: boolV,boolI
  complex(8),dimension(:,:),allocatable      :: Z
  complex(8),dimension(:),allocatable        :: Work
  real(8),dimension(:),allocatable           :: Rwork
  integer,dimension(:),allocatable           :: Iwork
  integer,dimension(:),allocatable           :: Isuppz
  integer,dimension(:),allocatable           :: Ifail
  !
  method_='zheevd';if(present(method))method_=trim(method)
  jobz_='V'  ;if(present(jobz))jobz_=jobz
  uplo_='U'  ;if(present(uplo))uplo_=uplo
  vl_  = 1d0 ;if(present(vL))vL_=vL
  vu_  = 1d0 ;if(present(vU))vU_=vU
  iL_  = 1   ;if(present(iL))iL_=iL
  iU_  = 1   ;if(present(iU))iU_=iU
  tol_ = 0d0 ;if(present(tol))tol_=tol
  !
  E=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  if(jobz_/='V'.AND.jobz_/='N')stop "zeigh_simple error: jobz has illegal value"
  if(uplo_/='U'.AND.uplo_/='L')stop "zeigh_simple error: uplo has illegal value"
  !
  Ns = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "zeigh_simple error: A has illegal shape"
  !
  mE = Ns
  select case(method_)
  case ("zheevr")
     allocate(Isuppz(2*Ns))
     allocate(Z(Ns,Ns))
     call zheevr(jobz_,range,uplo_,&
          Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,&
          ME,E,Z,Ns,&
          Isuppz,guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1)
     lrwork= guess_lrwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call zheevr(jobz_,range,uplo_,&
          Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,&
          ME,E,Z,Ns,&
          Isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
     !<copy the Evecs from Z to the input matrix A
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)    
     !
  case ("zheev")
     allocate(rwork(max(1,3*Ns)))
     call zheev(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,rwork,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call zheev(jobz_,uplo_,Ns,A,Ns,E,work,lwork,rwork,info)
     !
  case default
     ! case("zheevd")
     call zheevd(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1) ; lrwork= guess_lrwork(1) ; liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call zheevd(jobz_,uplo_,Ns,A,Ns,E,work,lwork,rwork,lrwork,iwork,liwork,info)
     !
  case ("zheevx")
     allocate(Z(Ns,Ns))
     allocate(Ifail(Ns))
     allocate(rwork(7*Ns))
     allocate(iwork(5*Ns))
     call zheevx(jobz_,range,uplo_,Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,guess_lwork,-1,rwork,iwork,ifail,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call zheevx(jobz_,range,uplo_,Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,work,lwork,rwork,iwork,ifail,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
  end select
  return
end subroutine zeigh_simple






! subroutine deigh_simple(A,W,jobz,uplo,vl,vu,il,iu,tol)
!   real(8),dimension(:,:),intent(inout)       :: A ! M v = E v/v(i,j) = ith component of jth vec.
!   real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
!   character(len=1),optional                  :: jobz,uplo
!   character(len=1)                           :: jobz_,uplo_,range
!   real(8),optional                           :: vl,vu,tol
!   integer,optional                           :: il,iu
!   real(8)                                    :: vL_,vU_,tol_
!   integer                                    :: iL_,iU_
!   integer                                    :: i,j,n,lda,info,ldz
!   integer                                    :: lwork,liwork,m
!   integer                                    :: guess_iwork(1)
!   real(8)                                    :: guess_work(1)
!   logical                                    :: boolV,boolI
!   real(8),dimension(:,:),allocatable         :: Z
!   real(8),dimension(:),allocatable           :: work
!   integer,dimension(:),allocatable           :: Isuppz
!   integer,dimension(:),allocatable           :: Iwork
!   !
!   !
!   jobz_='V'  ;if(present(jobz))jobz_=jobz
!   uplo_='U'  ;if(present(uplo))uplo_=uplo
!   vl_  = 1d0 ;if(present(vL))vL_=vL
!   vu_  = 1d0 ;if(present(vU))vU_=vU
!   iL_  = 1   ;if(present(iL))iL_=iL
!   iU_  = 1   ;if(present(iU))iU_=iU
!   tol_ = 0d0 ;if(present(tol))tol_=tol
!   !
!   W=0d0
!   !
!   range='A'
!   boolV=present(vL).AND.present(vU)
!   boolI=present(iL).OR.present(iU)
!   if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
!   if(boolV)range='V'
!   if(boolI)range='I'
!   !
!   N     = max(1,size(A,1))
!   LDA   = max(1,size(A,2))
!   ! M     = N ; if(range=="I")M = iU_ - iL_ + 1
!   LDZ   = N
!   if(any(shape(A)/=[N,N]))stop "my_eighD error: A has illegal shape"
!   !
!   allocate(Z(N,N)) !Supplying N columns is always safe.
!   !
!   allocate(Isuppz(2*N))
!   !
!   call dsyevr(jobz_,range,uplo_,N,A,N,vl_,vu_,iL_,iU_,tol_,M,W,Z,N,Isuppz,guess_work,-1,guess_iwork,-1,info)
!   lwork = guess_work(1)
!   liwork= guess_iwork(1)
!   allocate(work(lwork))
!   allocate(iwork(liwork))
!   call dsyevr(jobz_,range,uplo_,N,A,N,vl_,vu_,iL_,iU_,tol_,M,W,Z,N,Isuppz,work,lwork,iwork,liwork,info)
!   !<copy the Evecs from Z to the input matrix A
!   !print*,range,M
!   A = 0d0
!   A(:,1:M) = Z(:,1:M)    
! end subroutine deigh_simple

! subroutine zeigh_simple(A,W,jobz,uplo,vl,vu,il,iu,tol)
!   complex(8),dimension(:,:),intent(inout)    :: A ! M v = E v/v(i,j) = ith component of jth vec.
!   real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
!   character(len=1),optional                  :: jobz,uplo
!   character(len=1)                           :: jobz_,uplo_,range
!   real(8),optional                           :: vl,vu,tol
!   integer,optional                           :: il,iu
!   real(8)                                    :: vL_,vU_,tol_
!   integer                                    :: iL_,iU_
!   integer                                    :: i,j,n,lda,info,ldz
!   integer                                    :: lwork,liwork,lrwork,m
!   integer                                    :: guess_iwork(1)
!   real(8)                                    :: guess_rwork(1)
!   complex(8)                                 :: guess_work(1)
!   logical                                    :: boolV,boolI
!   complex(8),dimension(:,:),allocatable      :: Z
!   complex(8),dimension(:),allocatable        :: work
!   complex(8),dimension(:),allocatable        :: rwork
!   integer,dimension(:),allocatable           :: Isuppz
!   integer,dimension(:),allocatable           :: Iwork
!   !
!   !
!   jobz_='V'  ;if(present(jobz))jobz_=jobz
!   uplo_='U'  ;if(present(uplo))uplo_=uplo
!   vl_  = 1d0 ;if(present(vL))vL_=vL
!   vu_  = 1d0 ;if(present(vU))vU_=vU
!   iL_  = 1   ;if(present(iL))iL_=iL
!   iU_  = 1   ;if(present(iU))iU_=iU
!   tol_ = 0d0 ;if(present(tol))tol_=tol
!   !
!   W=0d0
!   !
!   range='A'
!   boolV=present(vL).AND.present(vU)
!   boolI=present(iL).OR.present(iU)
!   if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
!   if(boolV)range='V'
!   if(boolI)range='I'
!   !
!   N     = max(1,size(A,1))
!   LDA   = max(1,size(A,2))
!   ! M     = N ; if(range=="I")M = iU_ - iL_ + 1
!   LDZ   = N
!   if(any(shape(A)/=[N,N]))stop "my_eighD error: A has illegal shape"
!   !
!   allocate(Z(N,N)) !Supplying N columns is always safe.
!   !
!   allocate(Isuppz(2*N))
!   !
!   call zheevr(jobz_,range,uplo_,N,A,N,vl_,vu_,iL_,iU_,tol_,M,W,Z,N,Isuppz,guess_work,-1,guess_rwork,-1,guess_iwork,-1,info)
!   lwork = guess_work(1)
!   lrwork= guess_rwork(1)
!   liwork= guess_iwork(1)
!   allocate(work(lwork))
!   allocate(rwork(lrwork))
!   allocate(iwork(liwork))
!   call zheevr(jobz_,range,uplo_,N,A,N,vl_,vu_,iL_,iU_,tol_,M,W,Z,N,Isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
!   !<copy the Evecs from Z to the input matrix A
!   A = dcmplx(0d0,0d0)
!   A(:,1:M) = Z(:,1:M)    
! end subroutine zeigh_simple


! subroutine deigh_simple(M,E,jobz,uplo)
!   real(8),dimension(:,:),intent(inout) :: M ! M v = E v/v(i,j) = ith component of jth vec.
!   real(8),dimension(:),intent(inout)   :: E ! eigenvalues
!   character(len=1),optional            :: jobz,uplo
!   character(len=1)                     :: jobz_,uplo_
!   integer                              :: i,j,n,lda,info,lwork
!   real(8),dimension(:),allocatable     :: work
!   real(8),dimension(1)                 :: lwork_guess
!   jobz_='V';if(present(jobz))jobz_=jobz
!   uplo_='U';if(present(uplo))uplo_=uplo
!   lda = max(1,size(M,1))
!   n   = size(M,2)
!   call assert_shape(M,[n,n],"eigh","M")
!   Call dsyev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,info)
!   if (info /= 0) then
!      print*, "dsyevd returned info =", info
!      if (info < 0) then
!         print*, "the", -info, "-th argument had an illegal value"
!      else
!         print*, "the algorithm failed to compute an eigenvalue while working"
!         print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
!         print*, "through", mod(info, n+1)
!      end if
!      stop 'error deigh: 1st call dsyev'
!   end if
!   lwork=lwork_guess(1)
!   allocate(work(lwork))
!   call dsyev(jobz_,uplo_,n,M,lda,E,work,lwork,info)
!   if (info /= 0) then
!      print*, "dsyevd returned info =", info
!      if (info < 0) then
!         print*, "the", -info, "-th argument had an illegal value"
!      else
!         print*, "the algorithm failed to compute an eigenvalue while working"
!         print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
!         print*, "through", mod(info, n+1)
!      end if
!      stop 'error deigh: 2nd call dsyev'
!   end if
!   deallocate(work)
! end subroutine deigh_simple
! !-----------------------------
! subroutine zeigh_simple(M,E,jobz,uplo)
!   complex(8),dimension(:,:),intent(inout) :: M! M v = E v/v(i,j) = ith component of jth vec.
!   real(8),dimension(:),intent(inout)      :: E! eigenvalues
!   character(len=1),optional               :: jobz,uplo
!   character(len=1)                        :: jobz_,uplo_
!   integer                                 :: i,j,n,lda,info,lwork
!   complex(8),dimension(1)                 :: lwork_guess
!   complex(8),dimension(:),allocatable     :: work
!   real(8),dimension(:),allocatable        :: rwork
!   !write(*,*)"matrix_diagonalization called with: jobz="//jobz_//" uplo="//uplo_
!   jobz_='V';if(present(jobz))jobz_=jobz
!   uplo_='U';if(present(uplo))uplo_=uplo
!   lda = max(1,size(M,1))
!   n   = size(M,2)
!   call assert_shape(M,[n,n],"eigh","M")
!   allocate(rwork(max(1,3*N-2)))
!   call zheev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,rwork,info)
!   if(info/=0) then
!      print*, "zheev returned info =", info
!      if (info < 0) then
!         print*, "the", -info, "-th argument had an illegal value"
!      else
!         print*, "the algorithm failed to compute an eigenvalue while working"
!         print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
!         print*, "through", mod(info, n+1)
!      end if
!      stop 'error zeigh: 1st call zheev'
!   end if
!   lwork=lwork_guess(1) ; allocate(work(lwork))
!   call zheev(jobz_,uplo_,n,M,lda,E,work,lwork,rwork,info)
!   if(info/=0) then
!      print*, "zheev returned info =", info
!      if (info < 0) then
!         print*, "the", -info, "-th argument had an illegal value"
!      else
!         print*, "the algorithm failed to compute an eigenvalue while working"
!         print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
!         print*, "through", mod(info, n+1)
!      end if
!      stop 'error zeigh: 2nd call zheev'
!   end if
!   deallocate(work,rwork)
! end subroutine zeigh_simple




subroutine deigh_tridiag(D,U,Ev,Irange,Vrange)
  !This subroutine computes the eigenvalues and, optionally, the eigenvectors of a real symmetric tridiagonal matrix
  !using the LAPACK routine :f:func_inline:`dstevr`. The matrix is defined by its main diagonal :f:var:`D` and its
  !off-diagonal :f:var:`U`, :code:`U(i) = A(i+1,i) = A(i,i+1)`. It is an instance of the generic interface
  !:f:func_inline:`eigh`.
  !
  !Both :f:var:`D` and :f:var:`U` are overwritten. On output the first :code:`m` elements of :f:var:`D` contain the
  !eigenvalues in ascending order, with :code:`m=n` by default, :code:`m=Irange(2)-Irange(1)+1` if :f:var:`Irange` is
  !given, or the number of eigenvalues found if :f:var:`Vrange` is given. The eigenvectors are computed only if
  !:f:var:`Ev` is present, and are stored column-wise. A subset of the spectrum is selected either by index with
  !:f:var:`Irange` = :code:`[il,iu]`, or by value with :f:var:`Vrange` = :code:`[vl,vu]` (interval :math:`(v_l,v_u]`).
  !
  !The program stops if both :f:var:`Irange` and :f:var:`Vrange` are present, if :f:var:`Ev` does not have shape
  !:code:`[n,m]` (:code:`m=n` for :f:var:`Vrange`), or if :f:func_inline:`dstevr` returns an error code.
  !
  real(8),dimension(:)                :: d       ! in: main diagonal, size n; out: eigenvalues (first m valid)
  real(8),dimension(max(1,size(d)-1)) :: u       ! off-diagonal, size n-1; destroyed on exit
  real(8),dimension(:,:),optional     :: Ev      ! optional out: eigenvectors, [n,m]
  integer,dimension(2),optional       :: Irange  ! optional: [il,iu] indices of the first and last eigenvalue wanted
  integer,dimension(2),optional       :: Vrange  ! optional: [vl,vu] eigenvalue interval (vl,vu]
  !
  integer                             :: n
  ! = 'N':  Compute eigenvalues only;
  ! = 'V':  Compute eigenvalues and eigenvectors.
  character                           :: jobz_
  ! = 'A': all eigenvalues will be found
  ! = 'V': all eigenvalues in the half-open interval (VL,VU] will be found.
  ! = 'I': the IL-th through IU-th eigenvalues will be found.
  ! For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
  ! DSTEIN are called
  character                           :: range_
  real(8)                             :: vl,vu 
  integer                             :: il,iu
  real(8),dimension(size(d))          :: w
  real(8),dimension(:,:),allocatable  :: z
  integer,dimension(:),allocatable    :: isuppz
  integer                             :: lwork,liwork
  real(8),dimension(:),allocatable    :: work
  integer,dimension(:),allocatable    :: iwork
  real(8)                             :: abstol = 0d0 !tolerance on approximation error of eigenvalues
  integer                             :: m
  integer                             :: ldz
  integer                             :: info
  !
  n = size(d)
  !
  jobz_ = 'N'; if(present(Ev))jobz_='V'
  !
  range_= 'A'
  m     = n
  if(present(Irange).AND.present(Vrange))stop "EIGH_T: Irange & Vrange both present"
  if(present(Irange))then
     range_ = 'I'
     il = Irange(1)
     iu = Irange(2)
     m  = iu-il+1
  endif
  if(present(Vrange))then
     range_ = 'V'
     vl     = Vrange(1)
     vu     = Vrange(2)
     m      = n
  endif
  !
  if(present(Ev))then
     if(any(shape(Ev)/=[n,m]))stop "EIGH_T ERROR: wrong dimension in Ev"
  endif
  !
  ldz = n
  !
  allocate(z(n,m), isuppz(2*max(1,m)))
  allocate(work(1),iwork(1))
  call dstevr(jobz_,range_,n,d,u,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,-1,iwork,-1,info)
  if(info/=0) then
     print*, "dstevr returned info =", info
     if (info < 0) then
        print*, "the", -info, "-th argument had an illegal value"
     else
        print*, "Internal error"
     end if
     stop 'error dstevr: 1st call dstevr'
  end if
  !
  lwork  = work(1)
  liwork = iwork(1)
  deallocate(work,iwork);allocate(work(lwork),iwork(liwork))
  !
  call dstevr(jobz_,range_,n,d,u,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
  if(info/=0) then
     print*, "dstevr returned info =", info
     if (info < 0) then
        print*, "the", -info, "-th argument had an illegal value"
     else
        print*, "Internal error"
     end if
     stop 'error dstevr: 2nd call dstevr'
  end if
  !
  d  = w
  if(present(Ev)) Ev = z
  !
  deallocate(work,iwork,isuppz,z)
  !
end subroutine deigh_tridiag
