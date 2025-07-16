subroutine lanczos_arpack_d(MatVec,eval,evec,Nblock,Nitermax,bmat,v0,tol,iverbose,iexit)
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  !Arguments
  real(8)                   :: eval(:)![Neigen]
  real(8)                   :: evec(:,:)![Ns,Neigen]
  integer,optional          :: Nblock
  integer,optional          :: Nitermax
  ! character(len=2),optional :: which
  character(len=1),optional :: bmat
  real(8),optional          :: v0(size(evec,1))!(ns)
  real(8),optional          :: tol
  logical,optional          :: iverbose
  integer,optional          :: iexit
  !Dimensions:
  integer                   :: Ns
  integer                   :: Neigen
  integer                   :: maxn,maxnev,maxncv,ldv
  integer                   :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable       :: ax(:),d(:,:)
  real(8),allocatable       :: resid(:)
  real(8),allocatable       :: workl(:),workd(:)
  real(8),allocatable       :: v(:,:)
  logical,allocatable       :: select(:)
  integer                   :: iparam(11)
  integer                   :: ipntr(11)
  !Control Vars:
  integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                   :: rvec,verb
  integer                   :: i
  real(8)                   :: sigma
  real(8)                   :: tol_
  character                 :: bmat_
  character(len=2)          :: which_
  real(8),external          :: dnrm2
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
  Ns     = size(evec,1)
  Neigen = size(eval)
  call assert_shape(Evec,[Ns,Neigen],"Arpack_d","Evec")
  !
  maxn   = Ns
  maxnev = Neigen
  !
  maxncv = 10*Neigen ; if(present(Nblock))maxncv = Nblock
  maxitr = 512       ; if(present(Nitermax))maxitr = Nitermax
  bmat_  = 'I'       ; if(present(bmat))bmat_=bmat
  which_='SA'        !; if(present(which))which_=which
  tol_  = 0d0        ; if(present(tol))tol_=tol
  verb  =.false.     ; if(present(iverbose))verb=iverbose
  !
  if(bmat_/='I'.AND.bmat_/='G')stop "ARPACK: selected *bmat* is wrong. can be [I, G]"
  ! if(&
  !      which_/="LA" .OR. &
  !      which_/="SA" .OR. &
  !      which_/="LM" .OR. &
  !      which_/="SM" .OR. &
  !      which_/="BE")then
  !    stop "ARPACK: selected *which_* is wrong. can be [LA, SA, LM, SM, BE]"
  ! endif
  if(verb)then
     ndigit=-4
     logfil = 6
     mcaupd=1;mnaupd=1
     mcaup2=1;mnaup2=1
     mceupd=4;mneupd=4
  endif
  !
  ldv    = Ns  
  maxncv = max(maxncv,maxnev+2)
  if(maxncv>Ns)maxncv=Ns
  n      = maxn
  nev    = maxnev
  ncv    = maxncv
  ! 
  allocate(ax(n))
  allocate(d(ncv,2))
  allocate(resid(n))
  allocate(workl(ncv*(ncv+8)))
  allocate(workd(3*n))
  allocate(v(ldv,ncv))
  allocate(select(ncv))
  ax        = 0d0
  d         = 0d0
  resid     = 0d0
  workl     = 0d0
  workd     = 0d0
  v         = 0d0
  lworkl    = ncv*(ncv+8)
  info      = 1
  ido       = 0
  ishfts    = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  mode1     = 1
  iparam(7) = mode1
  if(present(v0))then
     resid=v0
  else
     call mt_random(resid)
  endif
  resid=resid/sqrt(dot_product(resid,resid))
  !
  !MAIN LOOP, REVERSE COMMUNICATION
  do
     call dsaupd(ido,bmat_,n,which_,nev,tol_,resid,ncv,v,ldv,&
          iparam,ipntr,workd,workl,lworkl,info)
     if(ido/=-1.AND.ido/=1)then
        exit
     end if
     !  Perform matrix vector multiplication:
     !  y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
  end do
  !
  if(present(iexit))iexit=info
  !
  !POST PROCESSING:
  if(info>=0)then
     if (info > 0)then
        write(*,'(a,i6)')'Soft failure in DSAUPD, info = ',info
        include "error_msg_arpack.h90"
     endif
     rvec = .true.
     call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat_,n,which_,&
          nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     do j=1,neigen
        eval(j)=d(j,1)
        do i=1,ns
           evec(i,j)=v(i,j)
        enddo
     enddo
     !
     !  Compute the residual norm ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and eigenvectors.
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with DSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of DSEUPD.'
        stop
     else
        nconv =  iparam(5)
     end if
     !
     include "info_msg_arpack.h90"
     !
     !if(nconv == 0) stop "ARPACK:: no converged eigenvalues have been found."
     !
  else
     write(*,'(a,i6)')'Hard failure in DSAUPD, info = ',info
     include "error_msg_arpack.h90"
  endif

  deallocate(ax,d,resid,workl,workd,v,select)
end subroutine lanczos_arpack_d
