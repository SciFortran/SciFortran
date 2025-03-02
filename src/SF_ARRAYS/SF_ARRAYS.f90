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
  !
  !Returns an array of evenly spaced numbers over a specified interval.
  !Returns :f:var:`num` evenly spaced samples, calculated over the interval [:f:var:`start`, :f:var:`stop`].
  !The start and end points of the interval can optionally be excluded.
  !
    real(8)          :: start           !Starting value of the sequence
    real(8)          :: stop            !End value of the sequence
    integer          :: num             !Number of samples to generate
    logical,optional :: istart          !If :code:`.true.`, :f:var:`start` is included in the resulting array. Default :code:`.true.`
    logical,optional :: iend            !If :code:`.true.`, :f:var:`stop` is included in the resulting array. Default :code:`.true.`
    real(8),optional :: mesh            !If present, the step is saved in this variable
    real(8)          :: array(num)      !Contains :f:var:`num` equally spaced samples in the interval 
                                        ![:f:var:`start`, :f:var:`stop`], left/right open or closed depending on 
                                        !:f:var:`istart` and :f:var:`iend`
    integer          :: i               !
    real(8)          :: step            
    logical          :: startpoint_,endpoint_
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
  !
  !Returns numbers spaced evenly on a log scale.
  !In linear space, the sequence starts at :f:var:`start` and ends with :f:var:`stop` (differently from numpy).
  !
    real(8)          :: start      ! The starting value of the sequence. Must be positive. If set to :code:`0`, 
                                   ! it will be reshifted to :code:`1e-12`
    real(8)          :: stop       ! The end value of the sequence. Must be positive. If set to :code:`0`, 
                                   ! it will be reshifted to :code:`1e-12`
    integer          :: num        ! Number of samples to generate
    real(8),optional :: base       ! The base of the exponential. Default :code:`10`
    real(8)          :: array(num) ! Contains :f:var:`num` samples, equally spaced on a log scale 
                                   ! in the closed interval [start, stop]
    integer          :: i
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
  function arange(start,num) result(array)
  !
  !Returns an array of :f:var:`num` integers starting with :f:var:`start`
  !
    integer          :: start         !First element of the array
    integer          :: num           !Length of the array
    integer          :: array(num)    !Contains the integer numbers in [:f:var:`start`, 
                                      !:f:var:`start`:code:`+`:f:var:`num`:code:`-1`]
    integer          :: i
    if(num<0)stop "arange: N<0, abort."
    forall(i=1:num)array(i)=start+i-1
  end function arange



  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function upmspace(start,stop,p,u,ndim,base,istart,iend,mesh) result(aout)
    !
    !Returns an array of number spaced linearly between exponentially spaced checkpoints.
    !The interval [:f:var:`start`, :f:var:`stop`] is first divided into :f:var:`p`
    !coarse regions, the width of which is exponentially increasing so that the 
    !i-th checkpoint is (:f:var:`stop` - :f:var:`start`) :math:`\cdot` :f:var:`base` ^ ( - :f:var:`p` + :code:`i` ).
    !Each of the coarse intervals is then linearly divided in :f:var:`u` subintervals.
    !
    real(8)          :: start  !First element of the array
    real(8)          :: stop   !Last element of the array
    integer          :: p      !Number of coarse subdivisions
    integer          :: u      !Number of fine subdivisions
    integer          :: ndim   !Length of the out array, must be :math:`p \cdot u` or  :math:`p \cdot u + 1`
    real(8),optional :: base   !Base of the exponential spacing (default :code:`2`)
    logical,optional :: istart !If :code:`.true.`, :f:var:`start` is included in the resulting array. Default :code:`.true.`
    logical,optional :: iend   !If :code:`.true.`, :f:var:`stop` is included in the resulting array. Default :code:`.true.`
    real(8),optional :: mesh(ndim) !If present, contains the distances between consecutive points in :f:var:`aout`. 
                                   !The last element is :code:`aout(ndim) - aout(ndim-1)`
    real(8)          :: aout(ndim) !Contains :f:var:`p` coarse exponentially-spaced checkpoints, 
                                   !each two of which separated by :f:var:`u` linearly spaced points
    real(8)          :: step,array(p*u+1)
    integer          :: pindex,uindex,pa,pb
    real(8)          :: ustart,ustop
    integer          :: i,j
    logical          :: endpoint_,startpoint_,check
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
  function upminterval(start,stop,midpoint,p,q,type,base,mesh) result(array)
    !
    !Returns an array of number spaced linearly between exponentially spaced checkpoints,
    !specularly around a middle point between the interval extrema. This is achieved by
    !calling :f:func:`upmspace` twice in a mirror-like fashion. The exponential
    !thickening of the mesh can be around the middle ( :f:var:`type` :code:`!=0` ) or the 
    !boundaries( :f:var:`type` :code:`=0` )
    !
    real(8)  :: start                   !First element of the array
    real(8)  :: stop                    !Last element of the array
    real(8)  :: midpoint                !Middle point of the array. The distance between the coarse checkpoints 
                                        !behaves specularly on the two sides
    integer  :: p                       !Number of coarse subdivisions
    integer  :: q                       !Number of fine subdivisions
    integer,optional :: type            !If :code:`=0`, mesh is thicker around :f:var:`start` and :f:var:`stop`. 
                                        !If :code:`=1`, mesh is thicker around :f:var:`midpoint`. Default :code:`0`
    real(8),optional :: base            !Base of the exponential spacing (default :code:`2`)
    real(8),optional :: mesh(2*P*Q+1)   !If present, contains the distances between consecutive points in :f:var:`array`. 
                                        !The last element is :code:`array(2 · p · q+1) - array(2 · p · q)`
    real(8)          :: array(2*p*q+1)  !Contains :math:`2 \cdot p \cdot q+1` points
    real(8)          :: base_
    integer          :: i
    integer          :: N
    integer          :: Nhalf
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
  function powspace(start,stop,num,base) result(array)
  !
  !Returns numbers spaced evenly on an exponential  scale.
  !The :code:`i`-th number is :f:var:`start` +
  !(:f:var:`stop` - :f:var:`start`) :math:`\cdot` :f:var:`base` ^ ( - :f:var:`num` + :code:`i` )
  !
    real(8)          :: start !first element of the array
    real(8)          :: stop  !last element of the array
    integer          :: num   !number of elements of the array
    real(8),optional :: base  !base of the exponential (defaule :code:`2`)
    real(8)          :: array(num) !contains :f:var:`num` exponentially spaced reals
    integer          :: i
    real(8)          :: base_
    if(num<0)stop "powspace: N<0, abort."
    base_= 2.d0;if(present(base))base_=base
    array(1) = start
    forall(i=2:num)array(i)=start + (stop-start)*base_**(-num+i)
  end function powspace




end module SF_ARRAYS
