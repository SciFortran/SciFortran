subroutine lanczos_parpack_d(MpiComm,MatVec,eval,evec,Nblock,Nitermax,v0,tol,iverbose,vrandom)
  !Arguments
  integer                   :: MpiComm
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(nchunk,vin,vout)
       integer                   :: nchunk
       real(8),dimension(nchunk) :: vin,vout
     end subroutine MatVec
  end interface
  !Arguments
  real(8)                   :: eval(:)![Neigen]
  real(8)                   :: evec(:,:)![Ns,Neigen]
  integer,optional          :: Nblock
  integer,optional          :: Nitermax
  ! character(len=2),optional :: which
  real(8),optional          :: v0(size(evec,1))
  real(8),optional          :: tol
  logical,optional          :: iverbose
  logical,optional          :: vrandom
  !Dimensions:
  integer                   :: Ns
  integer                   :: Neigen
  integer                   :: maxn,maxnev,maxncv,ldv,maxncv_
  integer                   :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable       :: evec_tmp(:) ![Nloc] see above
  real(8),allocatable       :: ax(:),d(:,:)
  real(8),allocatable       :: resid(:),vec(:)
  real(8),allocatable       :: workl(:),workd(:)
  real(8),allocatable       :: v(:,:)
  logical,allocatable       :: select(:)
  integer                   :: iparam(11)
  integer                   :: ipntr(11)
  !Control Vars:
  integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                   :: verb,vran
  integer                   :: i
  real(8)                   :: sigma,norm,norm_tmp
  real(8)                   :: tol_
  character                 :: bmat  
  character(len=2)          :: which_
  real(8),external          :: dnrm2
  !MPI
  logical                   :: mpi_master
  !
  integer ::  logfil, ndigit, mgetv0,&
       msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
       mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
       mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
  common /debug/logfil, ndigit, mgetv0,&
       msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
       mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
       mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
  !
  if(MpiComm == MPI_COMM_NULL)return
  !
  !MPI setup:
  mpi_master= Get_master_MPI(MpiComm)
  !
  Ns     = size(evec,1)
  Neigen = size(eval)
  call assert_shape(Evec,[Ns,Neigen],"P_arpack_d","Evec")
  !
  maxn   = Ns
  maxnev = Neigen
  !  
  maxncv = 10*Neigen ; if(present(Nblock))maxncv = Nblock
  maxitr = 512       ; if(present(Nitermax))maxitr = Nitermax
  which_ = 'SA'     ! ; if(present(which))which_=which
  tol_   = 0d0       ; if(present(tol))tol_=tol
  verb   = .false.   ; if(present(iverbose))verb=iverbose
  vran   = .true.    ; if(present(vrandom))vran=vrandom


  if(verb)then
     ndigit=-4
     logfil = 6
     mcaupd=1;mnaupd=1
     mcaup2=1;mnaup2=1
     mceupd=4;mneupd=4
  endif
  !
  maxncv_=maxncv                !every rank can have a different maxncv
  if(maxncv_>Ns)then
     maxncv_=Ns                 !some rank may have ncv > Ns (Ns=N/#mpi)
     if(mpi_master.and.verb)then
       print*,"PARPACK WARNING Ncv > Ns: reset block size to ",Ns
     endif
  endif
  !
  !BUG FIX FOR THE BLOCK RESIZE STUCK BEHAVIOR, revised from Xuanyu Long
  maxncv=0
  call MPI_ALLREDUCE(maxncv_,maxncv,1,MPI_INTEGER,MPI_MIN,MpiComm,ierr)
  !call MPI_ALLREDUCE(MPI_IN_PLACE,maxncv,1,MPI_INTEGER,MPI_MIN,MpiComm,ierr)
  !
  !
  !=========================================================================
  ldv    = maxn
  n      = maxn
  nev    = maxnev
  ncv    = maxncv
  bmat   = 'I'
  !
  !=========================================================================
  ! Setup distribution of data to nodes:
  if ( ldv > maxn ) then
     stop ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
  else if ( nev > maxnev ) then
     stop ' ERROR with _SDRV1: NEV is greater than MAXNEV '
  else if ( ncv > maxncv ) then
     stop ' ERROR with _SDRV1: NCV is greater than MAXNCV '
  end if
  !
  allocate(ax(ldv))
  allocate(resid(ldv))
  allocate(workd(3*ldv))
  allocate(v(ldv,ncv))
  allocate(d(ncv,2))
  allocate(workl(ncv*(ncv+8)))
  allocate(select(ncv))
  !
  ax     =0.d0
  d      =0.d0
  resid  =0.d0
  workl  =0.d0
  workd  =0.d0
  v      =0.d0
  lworkl = ncv*(ncv+8)+10
  info   = 1
  ido    = 0
  select = .true.
  !
  if(present(v0))then
     resid = v0
  else
     allocate(vec(ldv))
     if(vran)then
        call mt_random(vec)
     else
        vec = 1d0              !start with unitary vector 1/sqrt(Ndim)
     endif
     norm_tmp = dot_product(vec,vec)
     norm = 0d0
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     resid=vec/sqrt(norm)
     deallocate(vec)
  endif
  !
  ishfts    = 1
  mode1     = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
  !
  !  MAIN LOOP (Reverse communication loop)
  do
     call pdsaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
          ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info )
     !
     if(ido/=-1.AND.ido/=1)exit
     !  Perform matrix vector multiplication:
     !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
  end do
  !
  !POST PROCESSING:
  if(info/=0)then
     if(mpi_master)then
        write(*,'(a,i6)')'Warning/Error in PDSAUPD, info = ', info
        include "error_msg_arpack.h90"
     endif
  else
     !
     call pdseupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,bmat,&
          ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)        
     do j=1,neigen
        eval(j)=d(j,1)
     enddo
     !
     evec=0d0
     do j=1,neigen
        do i=1,ldv
           evec(i,j)=v(i,j)
        enddo
     enddo
     !
     !  Compute the residual norm ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and eigenvectors.
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of PDSEUPD.'
        stop
     else
        nconv =  iparam(5)
     end if
     !
     if(mpi_master)then
        include "info_msg_arpack.h90"
     endif
     !
     ! if(mpi_master.and.nconv==0.and.verb)stop "None of the required values was found."
  endif
  deallocate(ax,resid,workd,v,d,workl,select)
end subroutine lanczos_parpack_d
