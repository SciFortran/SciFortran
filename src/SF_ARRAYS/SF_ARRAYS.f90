module SF_ARRAYS
!SciFortran module for array creation and manipulation
  implicit none

  !ADD BUILDING MATRICES ROUTINES

  !GRIDS:
  public :: linspace
  public :: logspace
  public :: arange
  public :: powspace
  public :: upmspace
  public :: upminterval


contains


  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function linspace(start,stop,num,istart,iend,mesh) result(array)
    real(8)          :: start,stop,step,array(num)
    integer          :: num,i
    logical,optional :: istart,iend
    logical          :: startpoint_,endpoint_
    real(8),optional :: mesh
    !
    if(num<0)stop "linspace: N<0, abort."
    !
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.;if(present(iend))endpoint_=iend
    !
    if(startpoint_.AND.endpoint_)then
       if(num<2)stop "linspace: N<2 with both start and end points"
       step = (stop-start)/(dble(num)-1d0)
       forall(i=1:num)array(i)=start + (dble(i)-1d0)*step
    elseif(startpoint_.AND.(.not.endpoint_))then
       step = (stop-start)/dble(num)
       forall(i=1:num)array(i)=start + (dble(i)-1d0)*step
    elseif(.not.startpoint_.AND.endpoint_)then
       step = (stop-start)/dble(num)
       forall(i=1:num)array(i)=start + dble(i)*step
    else
       step = (stop-start)/(dble(num)+1d0)
       forall(i=1:num)array(i)=start + dble(i)*step
    endif
    if(present(mesh))mesh=step
  end function linspace



  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function logspace(start,stop,num,base) result(array)
    real(8)          :: start,stop,array(num)
    integer          :: num,i
    ! logical,optional :: iend
    ! logical          :: endpoint_
    real(8),optional :: base
    real(8)          :: base_
    real(8)          :: A,B
    base_= 10.d0;if(present(base))base_=base
    if(num<0)stop "logspace: N<0, abort."
    A=start;if(start==0.d0)A=1.d-12
    B=stop;if(stop==0.d0)B=1.d-12
    A=log(A)/log(base_) ; B=log(B)/log(base_)
    array = linspace(A,B,num=num,iend=.true.)
    array = base_**array
  end function logspace



  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function arange(start,num,iend) result(array)
    integer          :: start,array(num)
    integer          :: num,i
    logical,optional :: iend
    logical          :: endpoint_
    if(num<0)stop "arange: N<0, abort."
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(endpoint_)then
       forall(i=1:num)array(i)=start+i-1
    else
       forall(i=1:num-1)array(i)=start+i-1
    end if
  end function arange





  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function upminterval(start,stop,midpoint,p,q,type,base,mesh) result(array)
    integer  :: i,p,q,N,Nhalf
    real(8)  :: start,stop,midpoint,array(2*p*q+1)
    real(8),optional :: base,mesh(2*P*Q+1)
    real(8)          :: base_
    integer,optional :: type
    integer          :: type_
    type_= 0          ;if(present(type))type_=type
    base_= 2.d0       ;if(present(base))base_=base
    N=2*p*q+1
    Nhalf=p*q
    if(type_==0)then
       array(1:Nhalf+1)   = upmspace(start,midpoint,p,q,Nhalf+1,base=base_)
       array(N:Nhalf+2:-1)= upmspace(stop,midpoint,p,q,Nhalf,base=base_,iend=.false.)
    else
       array(Nhalf+1:1:-1) = upmspace(midpoint,start,p,q,Nhalf+1,base=base_)
       array(Nhalf+2:N)    = upmspace(midpoint,stop,p,q,Nhalf,base=base_,istart=.false.)
    endif
    if(present(mesh))then
       do i=1,N-1
          mesh(i)=(array(i+1)-array(i))
       enddo
       Mesh(N)=(array(N)-array(N-1))
    endif
  end function upminterval


  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function upmspace(start,stop,p,u,ndim,base,istart,iend,mesh) result(aout)
    integer          :: p,u,ndim,pindex,uindex,pa,pb
    real(8)          :: start,stop,step,array(p*u+1),aout(ndim)
    real(8),optional :: mesh(ndim)
    real(8)          :: ustart,ustop
    integer          :: i,j
    logical,optional :: iend,istart
    logical          :: endpoint_,startpoint_,check
    real(8),optional :: base
    real(8)          :: base_
    ! real(8),optional :: mesh(p*u+1)
    if(ndim<0)stop "upmspace: N<0, abort."
    check=(ndim==(p*u)).OR.(ndim==(p*u+1))
    if(.not.check)stop "upmspace: wrong Ndim, abort."
    base_= 2.d0       ;if(present(base))base_=base
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.  ;if(present(iend))endpoint_=iend
    check=startpoint_.AND.endpoint_
    pindex=1
    array(pindex) = start
    do i=1,p
       pindex=1+i*u               !index of the next p-mesh point
       pa=1+(i-1)*u               !index of the previous p-mesh point
       pb=1+i*u                   !
       array(pindex)=start + (stop-start)*base_**(-p+i) !create p-mesh
       ustart = array(pa)         !u-interval start
       ustop  = array(pb)         !u-interval stop
       step   = (ustop-ustart)/dble(u) !u-interval step
       do j=1,u-1
          uindex=pa+j    !u-mesh points
          array(uindex)=ustart + dble(j)*step
       enddo
    enddo
    if(check)then
       aout(1:ndim)=array
    elseif(.not.endpoint_)then
       aout(1:ndim)=array(1:p*u)
    elseif(.not.startpoint_)then
       aout(1:ndim)=array(2:)   
    endif
    if(present(mesh))then
       do i=1,ndim-1
          mesh(i)=abs(aout(i+1)-aout(i))
       enddo
       mesh(ndim)=abs(aout(ndim)-aout(ndim-1))
    endif
  end function upmspace



  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function powspace(start,stop,num,base) result(array)
    real(8)          :: start,stop,step,array(num)
    integer          :: num,i
    real(8),optional :: base
    real(8)          :: base_
    if(num<0)stop "powspace: N<0, abort."
    base_= 2.d0;if(present(base))base_=base
    array(1) = start
    forall(i=2:num)array(i)=start + (stop-start)*base_**(-num+i)
  end function powspace




end module SF_ARRAYS
