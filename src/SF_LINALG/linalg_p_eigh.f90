subroutine p_deigh_simple(A,W,Nblock,method,jobz,uplo,vl,vu,il,iu,tol)
  !This subroutine is the distributed-memory (ScaLAPACK) counterpart of :f:func_inline:`deigh_simple`: it computes
  !the eigenvalues and, optionally, the eigenvectors of a real symmetric matrix :f:var:`A`. The complex version
  !:f:func_inline:`p_zeigh_simple` works on hermitian matrices. Both are instances of the generic interface
  !:f:func_inline:`p_eigh`, available only when compiling with :code:`_SCALAPACK`.
  !
  !The matrix :f:var:`A` (a full copy is expected on every process) is distributed on the BLACS process grid of the
  !module (:code:`p_context`, :code:`p_Nx` x :code:`p_Ny`) in a block-cyclic layout with square blocks of size
  !:f:var:`Nblock`, and diagonalized with the driver selected by :f:var:`method`: :f:func_inline:`PDSYEVR` (default),
  !:f:func_inline:`PDSYEV`, :f:func_inline:`PDSYEVD` or :f:func_inline:`PDSYEVX`. If :code:`jobz='V'` the eigenvectors
  !are gathered back into :f:var:`A`, column-wise. The eigenvalues are returned in ascending order in :f:var:`W`.
  !The selection of a subset of the spectrum by :f:var:`vl`, :f:var:`vu`, :f:var:`il`, :f:var:`iu` works as in
  !:f:func_inline:`deigh_simple` (:code:`PDSYEVR` and :code:`PDSYEVX` only). Timings are written by the master
  !process to the file :code:`p_eigh.info`.
  !
  !The program stops if both a value range and an index range are given, or if :f:var:`A` is not square.
  !
  real(8),dimension(:,:),intent(inout)       :: A       ! in: symmetric matrix; out: eigenvectors (columns) if jobz=V
  real(8),dimension(size(A,2)),intent(inout) :: W       ! eigenvalues, ascending, size n
  integer                                    :: Nblock  ! block size of the block-cyclic distribution
  integer                                    :: Nb
  character(len=*),optional                  :: method  ! optional: 'dsyevr' (default), 'dsyev', 'dsyevd', 'dsyevx'
  character(len=1),optional                  :: jobz    ! optional: 'V' (default) eigenvectors too, 'N' eigenvalues only
  character(len=1),optional                  :: uplo    ! optional: 'L' (default) or 'U', triangle of A to be used
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl      ! optional: lower bound of the eigenvalue interval (vl,vu]
  real(8),optional                           :: vu      ! optional: upper bound of the eigenvalue interval (vl,vu]
  real(8),optional                           :: tol     ! optional: accepted but currently not used
  integer,optional                           :: il      ! optional: index of the first eigenvalue wanted
  integer,optional                           :: iu      ! optional: index of the last eigenvalue wanted
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: Qrows,Qcols
  integer                                    :: i,j,lda,info,ldz
  integer                                    :: lwork,liwork,mW,mZ
  real(8),dimension(:),allocatable           :: work,Rwork,Gap
  integer                                    :: guess_liwork(1)
  real(8)                                    :: guess_lwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  integer,dimension(:),allocatable           :: Isuppz,Ifail,Iclustr
  integer,dimension(:),allocatable           :: Iwork
  !
  integer,external                           :: numroc,indxG2L,indxG2P,indxL2G
  real(8),external                           :: dlamch      
  !
  real(8),dimension(:,:),allocatable         :: A_loc,Z_loc
  integer                                    :: rankX,rankY
  integer                                    :: Nrow,Ncol
  integer                                    :: sendR,sendC
  integer                                    :: myi,myj,unit
  integer,dimension(9)                       :: descA,descAloc,descZloc
  real(8)                                    :: t_stop,t_start
  logical                                    :: master
  !
  method_='dsyevr';if(present(method))method_=trim(method)
  jobz_='V'      ;if(present(jobz))jobz_=jobz
  uplo_='L'      ;if(present(uplo))uplo_=uplo
  vl_  = 1d0     ;if(present(vL))vL_=vL
  vu_  = 1d0     ;if(present(vU))vU_=vU
  iL_  = 1       ;if(present(iL))iL_=iL
  iU_  = 1       ;if(present(iU))iU_=iU
  tol_ = dlamch('s')     ;if(present(tol))tol_=tol
  !
  W=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 100
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_eigh.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Qrows x Qcols=",Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols))
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  !< Distribute A
  call Distribute_BLACS(A,A_loc,descAloc,unit)
  !
  !< Allocate distributed eigenvector matrix
  allocate(Z_loc(Qrows,Qcols));Z_loc=0d0
  call descinit( descZloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)write(unit,"(A20,A21)")"Using Method:",method_
  if(master)call cpu_time(t_start)
  select case(method_)
  case default
     call  PDSYEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call  PDSYEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          work, lwork, iwork, liwork,info)
     !
  case ("dsyev")
     call PDSYEV(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call PDSYEV( jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,info)
     !
     !
  case ("dsyevd")
     call PDSYEVD(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_liwork,-1,INFO)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call PDSYEVD(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,iwork,liwork,info)
     !
  case ("dsyevx")
     allocate(Ifail(Ns))
     allocate(Iclustr(2*p_Nx*p_Ny))
     allocate(Gap(p_Nx*p_Ny))
     call PDSYEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW, mZ,W,-1d0,Z_loc,1,1,descAloc,&
          guess_lwork,-1,guess_liwork,-1,ifail,iclustr,gap,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call PDSYEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW,mZ,W,-1d0,Z_loc,1,1,descZloc,&
          work,lwork,iwork,liwork,ifail,iclustr,gap,info)
  end select
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time diag:",t_stop-t_start
  !
  if(jobz_=='V')then
     A=0d0
     call Gather_BLACS(Z_loc,A,descA,unit)
  endif
  if(master)close(unit)
100 continue
  return
  !
end subroutine p_deigh_simple




subroutine p_zeigh_simple(A,W,Nblock,method,jobz,uplo,vl,vu,il,iu,tol)
  complex(8),dimension(:,:),intent(inout)    :: A ! M v = E v/v(i,j) = ith component of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
  integer                                    :: Nblock
  integer                                    :: Nb
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: Qrows,Qcols
  integer                                    :: i,j,lda,info,ldz
  integer                                    :: lwork,liwork,lrwork,mW,mZ
  complex(8),dimension(:),allocatable        :: Work,RRwork
  real(8),dimension(:),allocatable           :: Rwork
  integer,dimension(:),allocatable           :: Iwork
  complex(8)                                 :: guess_lwork(1)
  complex(8)                                 :: guess_lrrwork(1)
  real(8)                                    :: guess_lrwork(1)
  integer                                    :: guess_liwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  integer,dimension(:),allocatable           :: Isuppz
  integer,dimension(:),allocatable           :: Ifail
  integer,dimension(:),allocatable           :: Iclustr
  real(8),dimension(:),allocatable           :: Gap
  !
  integer,external                           :: numroc,indxG2L,indxL2G
  real(8),external                           :: dlamch      
  !
  complex(8),dimension(:,:),allocatable      :: A_loc,Z_loc
  integer                                    :: rankX,rankY
  integer                                    :: Nrow,Ncol
  integer                                    :: sendR,sendC
  integer                                    :: myi,myj,unit
  integer,dimension(9)                       :: descA,descAloc,descZloc
  real(8)                                    :: t_stop,t_start
  logical                                    :: master
  !
  method_='zheevr';if(present(method))method_=trim(method)
  !
  jobz_='V'      ;if(present(jobz))jobz_=jobz
  uplo_='L'      ;if(present(uplo))uplo_=uplo
  vl_  = 1d0     ;if(present(vL))vL_=vL
  vu_  = 1d0     ;if(present(vU))vU_=vU
  iL_  = 1       ;if(present(iL))iL_=iL
  iU_  = 1       ;if(present(iU))iU_=iU
  tol_ = dlamch('s')     ;if(present(tol))tol_=tol
  !
  !>DEBUG
  if(method_=="zheev" .OR. method_=="zheevd") then
     write(0,"(A)")"WARNING p_eigh: pzheev & pzheevd are bugged. Fall on pzheevr."
     method_='zheevr'
  endif
  if(uplo_=='U') then
     write(0,"(A)")"WARNING p_eigh: pzheevr OR pzheevx are bugged. Only work with with uplo=L"
     uplo_  ='L'
  endif
  !<DEBUG
  !
  W=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 200
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_eigh.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Qrows x Qcols=",Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols));A_loc=zero
  !
  !< Distribute A
  call Distribute_BLACS(A,A_loc,descAloc,unit)
  !
  !< Allocate distributed eigenvector matrix
  allocate(Z_loc(Qrows,Qcols));Z_loc=zero
  call descinit( descZloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)write(unit,"(A20,A21)")"Using Method:",method_
  if(master)call cpu_time(t_start)
  select case(method_)
  case default
     call  PZHEEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1) ; lrwork= guess_lrwork(1)  ; liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call  PZHEEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          work, lwork, rwork, lrwork, iwork, liwork,info)
     !
  case ("zheevx")
     allocate(Ifail(Ns))
     allocate(Iclustr(2*p_Nx*p_Ny))
     allocate(Gap(p_Nx*p_Ny))
     call PZHEEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW, mZ,W,-1d0,Z_loc,1,1,descAloc,&
          guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,ifail,iclustr,gap,info)
     lwork = guess_lwork(1)
     lrwork= guess_lrwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call PZHEEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW,mZ,W,-1d0,Z_loc,1,1,descZloc,&
          work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)
     !
     !
     ! >>>ACTHUNG<< BUGGED DO NOT USE: --> zheevr, w UPLO='L'
  case ("zheev")              !Bugged
     call PZHEEV(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_lrwork,-1,info)
     lwork = guess_lwork(1) ; lrwork = guess_lrwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     call PZHEEV(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,rwork,lrwork,info)
     !
     !
     ! >>>ACTHUNG<< BUGGED DO NOT USE: --> zheevr, w UPLO='L'
  case ("zheevd")
     call PZHEEVD('V',uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,INFO)
     lwork = guess_lwork(1) ; lrwork= guess_lrwork(1) ; liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call PZHEEVD('V',uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,rwork,lrwork,iwork,liwork,info)
  end select
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time diag:",t_stop-t_start
  !
  if(jobz_=='V')then
     A=zero
     call Gather_BLACS(Z_loc,A,descA,unit)
  endif
  !
  if(master)close(unit)
200 continue
  return
end subroutine p_zeigh_simple

