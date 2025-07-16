subroutine lanczos_arpack_c(MatVec,eval,evec,Nblock,Nitermax,bmat,v0,tol,iverbose)
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: vin
       complex(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  !Arguments
  real(8)                      :: eval(:)![Neigen]
  complex(8)                   :: evec(:,:)![Ns,Neigen]
  integer,optional             :: Nblock
  integer,optional             :: Nitermax
  ! character(len=2),optional    :: which
  character(len=1),optional    :: bmat
  complex(8),optional          :: v0(size(evec,1))
  real(8),optional             :: tol
  logical,optional             :: iverbose
  !Dimensions:
  integer                      :: Ns
  integer                      :: Neigen
  integer                      :: maxn,maxnev,maxncv,ldv
  integer                      :: n,nconv,ncv,nev
  !Arrays:
  complex(8),allocatable       :: ax(:),d(:)
  complex(8),allocatable       :: v(:,:)
  complex(8),allocatable       :: workl(:),workd(:),workev(:)
  complex(8),allocatable       :: resid(:)
  real(8),allocatable          :: rwork(:),rd(:,:)
  logical,allocatable          :: select(:)
  integer                      :: iparam(11)
  integer                      :: ipntr(14)
  !Control Vars:
  integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                      :: rvec,verb
  integer                      :: i
  real(8)                      :: sigma
  real(8)                      :: tol_
  character                    :: bmat_
  character(len=2)             :: which_
  real(8),external             :: dznrm2,dlapy2
  real(8),allocatable          :: reV(:),imV(:)
  integer,allocatable          :: Eorder(:)

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
  call assert_shape(Evec,[Ns,Neigen],"Arpack_c","Evec")
  !
  maxn   = Ns
  maxnev = Neigen
  maxncv = 10*Neigen ; if(present(Nblock))maxncv = Nblock
  maxitr = 512       ; if(present(Nitermax))maxitr = Nitermax
  bmat_  = 'I'       ; if(present(bmat))bmat_=bmat
  which_='SR'       ! ; if(present(which))which_=which
  tol_  = 0d0        ; if(present(tol))tol_=tol
  verb  =.false.     ; if(present(iverbose))verb=iverbose
  !
  if(bmat_/='I'.AND.bmat_/='G')stop "ARPACK: selected *bmat* is wrong. can be [I, G]"
  ! which_/="LM" .OR. &
  ! which_/="SM" .OR. &       
  ! which_/="LI" .OR. &
  ! which_/="SI"
  ! if(which_/="LR" .OR. which_/="SR")then
  !    stop "ARPACK: selected *which_* is wrong. can be [LM, SM, LR, SR, LI, SI]"
  ! endif
  if(verb)then
     ndigit=-4
     logfil = 6
     mcaupd=1;mnaupd=1
     mcaup2=1;mnaup2=1
     mceupd=4;mneupd=4
  endif
  !
  ldv    = Ns        ; if(maxncv>Ns)maxncv=Ns
  n      = maxn
  nev    = maxnev
  ncv    = maxncv
  ! 
  allocate(ax(n))
  allocate(d(ncv))
  allocate(resid(n))
  allocate(v(ldv,ncv))
  allocate(workd(3*n))
  allocate(workev(3*ncv))
  allocate(workl(ncv*(3*ncv+5) + 10))
  allocate(rwork(ncv))
  allocate(rd(ncv,3))
  allocate(select(ncv))
  ax     = zero
  d      = zero
  v      = zero
  workl  = zero
  workd  = zero
  resid  = zero
  workev = zero
  rwork  = 0d0
  rd     = 0d0
  lworkl = ncv*(3*ncv+5) + 10
  info   = 1
  ido    = 0
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
     call znaupd(ido,bmat_,n,which_,nev,tol_,resid,ncv,v,ldv,&
          iparam,ipntr,workd,workl,lworkl,rwork,info)
     if(ido/=-1.AND.ido/=1)exit
     !  Perform matrix vector multiplication:
     !  y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
  end do
  !
  !POST PROCESSING:
  if(info>=0)then
     if (info > 0)then
        write(*,'(a,i6)')'Soft failure in ZNAUPD, info = ',info
        include "error_msg_arpack.h90"
     endif
     rvec = .true.
     call zneupd  (rvec,'A',select,d,v,ldv,sigma,workev,bmat_,n,which_,&
          nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
     !  Eigenvalues are returned in the first column of the two dimensional 
     !  array D and the corresponding eigenvectors are returned in the first 
     !  NCONV (=IPARAM(5)) columns of the two dimensional array V if requested.
     !  Otherwise, an orthogonal basis for the invariant subspace corresponding 
     !  to the eigenvalues in D is returned in V.
     do j=1,neigen
        eval(j)=dreal(d(j))
     end do
     allocate(Eorder(Neigen))
     call sort_array(Eval,Eorder)
     do j=1,neigen
        evec(:,j)=v(:,Eorder(j))
     enddo
     !
     !  Compute the residual norm ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and eigenvectors.
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with ZNEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of ZNEUPD.'
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
    write(*,'(a,i6)')'Hard failure in ZNAUPD, info = ',info
    include "error_msg_arpack.h90"
    STOP
  endif
  deallocate(ax,d,resid,v,workd,workev,workl,rwork,rd,select)
  !
contains

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
end subroutine lanczos_arpack_c
