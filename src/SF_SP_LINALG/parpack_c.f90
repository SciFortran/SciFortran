subroutine lanczos_parpack_c(MpiComm,MatVec,eval,evec,Nblock,Nitermax,v0,tol,iverbose,vrandom)
  !Arguments
  integer                    :: MpiComm
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(nchunk,vin,vout)
       integer                      :: nchunk
       complex(8),dimension(nchunk) :: vin,vout
     end subroutine MatVec
  end interface
  !Arguments
  real(8)                   :: eval(:)   ![Neigen]
  complex(8)                :: evec(:,:) ![Nloc,Neigen]
  integer,optional          :: Nblock
  integer,optional          :: Nitermax
  ! character(len=2),optional :: which
  complex(8),optional       :: v0(size(evec,1))
  real(8),optional          :: tol
  logical,optional          :: iverbose
  logical,optional          :: vrandom
  !Dimensions:
  integer                   :: Ns
  integer                   :: Neigen
  integer                   :: maxn,maxnev,maxncv,ldv,nloc
  integer                   :: n,nconv,ncv,nev
  !Arrays:
  complex(8),allocatable    :: evec_tmp(:) ![Nloc] see above
  complex(8),allocatable    :: ax(:),d(:)
  complex(8),allocatable    :: v(:,:)
  complex(8),allocatable    :: workl(:),workd(:),workev(:)
  complex(8),allocatable    :: resid(:),vec(:)
  real(8),allocatable       :: rwork(:),rd(:,:)
  logical,allocatable       :: select(:)
  integer                   :: iparam(11)
  integer                   :: ipntr(14)
  !Control Vars:
  integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                   :: verb,vran
  integer                   :: i
  real(8)                   :: sigma,norm,norm_tmp
  real(8)                   :: tol_
  character                 :: bmat  
  character(len=2)          :: which_
  real(8),external          :: dznrm2,dlapy2
  real(8),allocatable       :: reV(:),imV(:)
  integer,allocatable       :: Eorder(:)
  !MPI
  logical                   :: mpi_master

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
  which_ = 'SR'      !; if(present(which))which_=which
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
  !
  if(maxncv>Ns)then
     maxncv=Ns
     print*,"PARPACK WARNING Ncv > Ns: reset block size to ",Ns
  endif
  !BUG FIX FOR THE BLOCK RESIZE STUCK BEHAVIOR, from Xuanyu Long
  call MPI_ALLREDUCE(MPI_IN_PLACE,maxncv,1,MPI_INTEGER,MPI_MIN,MpiComm,ierr)
  !
  !
  n      = maxn
  nev    = maxnev
  ncv    = maxncv
  bmat   = 'I'
  !
  !=========================================================================
  ! Setup distribution of data to nodes:
  ldv = size(evec,1)             !ldv is the SMALL dimension
  ! if(ldv>0.AND.ldv < ncv)stop "LANCZOS_PARPACK_C error: ldv < maxNblock. Unstable! Find a way to increase ldv... (less cpu?)"
  if ( ldv > maxn ) then
     stop ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
  else if ( nev > maxnev ) then
     stop ' ERROR with _SDRV1: NEV is greater than MAXNEV '
  else if ( ncv > maxncv ) then
     stop ' ERROR with _SDRV1: NCV is greater than MAXNCV '
  end if
  !
  allocate(ax(ldv))
  allocate(d(ncv))
  allocate(resid(ldv))
  allocate(v(ldv,ncv))
  allocate(workd(3*ldv))
  allocate(workev(3*ncv))
  allocate(workl(ncv*(3*ncv+5) + 10))
  allocate(rwork(ncv))
  allocate(rd(ncv,3))
  allocate(select(ncv))
  !
  ax     = zero
  d      = zero
  v      = zero
  workl  = zero
  workd  = zero
  resid  = zero
  workev = zero
  rwork  = zero
  rd     = zero
  lworkl = ncv*(3*ncv+5) + 10
  info   = 1
  ido    = 0
  select = .true.
  !
  !
  if(present(v0))then
     resid = v0
  else
     if(vran)then
        call mt_random(resid)
     else
        resid = one               !start with unitary vector 1/sqrt(Ndim)
     endif
     norm_tmp = dot_product(resid,resid)
     norm = 0d0
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     resid=resid/sqrt(norm)
  endif
  !
  ishfts    = 1
  mode1     = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
  !
  !
  !  MAIN LOOP (Reverse communication loop)
  do
     call pznaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
          ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)
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
        write(*,'(a,i6)')'Warning/Error in PZNAUPD, info = ', info
        include "error_msg_arpack.h90"
     endif
  else
     !
     call pzneupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,workev,bmat,&
          ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
     !
     do j=1,neigen
        eval(j)=dreal(d(j))
     enddo
     allocate(Eorder(Neigen))
     call sort_array(Eval,Eorder)
     !
     evec=zero
     do j=1,neigen
        do i=1,ldv
           evec(i,j)=v(i,Eorder(j))
        enddo
     enddo
     !
     !  Compute the residual norm ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and eigenvectors.
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with PZNEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of PZNEUPD.'
     else
        nconv =  iparam(5)
     end if
     !
     if(mpi_master)then
        include "info_msg_arpack.h90"
     end if
     !
     !
     ! if(mpi_master.and.nconv==0.and.verb)stop "None of the required values was found."
  endif
  deallocate(ax,d,resid,v,workd,workev,workl,rwork,rd,select)
  !
  !
contains
  !
  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order2,no_touch_array)
    implicit none
    real(8),dimension(:)                    :: array
    real(8),dimension(size(array))          :: backup
    integer,dimension(size(array))          :: order
    integer,dimension(size(array)),optional :: order2
    integer                                 :: i
    logical,optional                        :: no_touch_array
    do i=1,size(order)
       order(i)=i
    enddo
    call qsort_sort( array, order, 1, size(array) )
    if(.not.present(no_touch_array))then
       do i=1,size(order)
          backup(i)=array(order(i))
       enddo
       array=backup
    endif
    if(present(order2)) order2=order
  end subroutine sort_array
  !---------------------------------------------!  
  recursive subroutine qsort_sort( array, order, left, right )
    implicit none
    real(8), dimension(:)                 :: array
    integer, dimension(:)                 :: order
    integer                               :: left
    integer                               :: right
    integer                               :: i
    integer                               :: last
    if ( left .ge. right ) return
    call qsort_swap( order, left, qsort_rand(left,right) )
    last = left
    do i = left+1, right
       if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
          last = last + 1
          call qsort_swap( order, last, i )
       endif
    enddo
    call qsort_swap( order, left, last )
    call qsort_sort( array, order, left, last-1 )
    call qsort_sort( array, order, last+1, right )
  end subroutine qsort_sort
  !---------------------------------------------!
  subroutine qsort_swap( order, first, second )
    implicit none
    integer, dimension(:)                 :: order
    integer                               :: first, second
    integer                               :: tmp
    tmp           = order(first)
    order(first)  = order(second)
    order(second) = tmp
  end subroutine qsort_swap
  !---------------------------------------------!
  integer function qsort_rand( lower, upper )
    implicit none
    integer                               :: lower, upper
    real(4)                               :: r
    r=drand()
    qsort_rand =  lower + nint(r * (upper-lower))
  end function qsort_rand
  !---------------------------------------------!
  function drand()
    implicit none
    real(8)                               :: drand
    real(4)                               :: r
    call random_number(r)
    drand=dble(r)
  end function drand
  !---------------------------------------------!
  function compare(f,g)
    implicit none
    real(8)                               :: f,g
    integer                               :: compare
    if(f<g) then
       compare=-1
    else
       compare=1
    endif
  end function compare
  !
  !
end subroutine lanczos_parpack_c









