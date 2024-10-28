subroutine bfgs_with_grad(func,grad,x,l,u,nbd,factr,pgtol,iprint,nloop)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8) :: func
     end function func
  end interface
  interface
     function grad(x)
       real(8),dimension(:) :: x
       real(8),dimension(size(x)) :: grad
     end function grad
  end interface
  integer                                   :: n,m = 5, iprint_ = -1
  integer,optional                          :: iprint,nloop
  real(8),optional                          :: factr, pgtol
  real(8)                                   :: factr_,pgtol_
  character(len=60)                         :: task, csave
  logical                                   :: lsave(4)
  integer                                   :: isave(44)
  real(8)                                   :: dsave(29),f
  integer,dimension(:),allocatable          :: iwa,nbd_
  integer,dimension(:),allocatable,optional :: nbd
  real(8),dimension(:),allocatable          :: x,l_,u_,g,wa
  real(8),dimension(:),allocatable,optional :: l,u
  !
  n=size(x)
  factr_  = 1.0d+7
  pgtol_  = 1.0d-5
  iprint_ = -1
  !
  allocate( g(n))
  allocate ( iwa(3*n) )
  allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
  task = 'START'
  if(present(factr))factr_=factr
  if(present(pgtol))pgtol_=pgtol
  if(present(iprint))iprint_=iprint
  if(present(nbd))then
     nbd_=nbd
     l_=l
     u_=u
  else
     allocate ( nbd_(n), l_(n), u_(n))
     nbd_=0
     l_=0.d0
     u_=0.d0
  endif
  !     The beginning of the loop
  do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
     !    
     !     This is the call to the L-BFGS-B code.
     call setulb ( n, m, x, l_, u_, nbd_, f, g, factr_, pgtol_, &
          wa, iwa, task, iprint_,&
          csave, lsave, isave, dsave )
     if (task(1:2) .eq. 'FG') then  
        f=func(x)
        g=grad(x)
     endif
     if(present(nloop))then
        if(isave(30) .ge. nloop)task='STOP: TOTAL NO. OF LOOPS EXCEEDS LIMIT'
     endif
  end do
  deallocate (nbd_,l_,u_,iwa,wa,g)
end subroutine bfgs_with_grad



subroutine bfgs_no_grad(func,x,l,u,nbd,factr,pgtol,iprint,nloop)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8) :: func
     end function func
  end interface
  integer                                   :: n,m = 5, iprint_ = -1
  integer,optional                          :: iprint,nloop
  real(8),optional                          :: factr, pgtol
  real(8)                                   :: factr_,pgtol_
  character(len=60)                         :: task, csave
  logical                                   :: lsave(4)
  integer                                   :: isave(44)
  real(8)                                   :: dsave(29),f
  integer,dimension(:),allocatable          :: iwa,nbd_
  integer,dimension(:),allocatable,optional :: nbd
  real(8),dimension(:),allocatable          :: x,l_,u_,g,wa
  real(8),dimension(:),allocatable,optional :: l,u
  !
  n=size(x)
  factr_  = 1.0d+7
  pgtol_  = 1.0d-5
  iprint_ = -1
  !
  allocate( g(n))
  allocate ( iwa(3*n) )
  allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
  task = 'START'
  if(present(factr))factr_=factr
  if(present(pgtol))pgtol_=pgtol
  if(present(iprint))iprint_=iprint
  if(present(nbd))then
     nbd_=nbd
     l_=l
     u_=u
  else
     allocate ( nbd_(n), l_(n), u_(n))
     nbd_=0
     l_=0.d0
     u_=0.d0
  endif
  !     The beginning of the loop
  do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
     !     This is the call to the L-BFGS-B code.
     call setulb ( n, m, x, l_, u_, nbd_, f, g, factr_, pgtol_, &
          wa, iwa, task, iprint_,&
          csave, lsave, isave, dsave )
     if (task(1:2) .eq. 'FG') then  
        f=func(x)
        g=func_func_jacobian(func,size(x),x)
     endif
     if(present(nloop))then
        if(isave(30) .ge. nloop)task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
     endif
  end do
  deallocate (nbd_,l_,u_,iwa,wa,g)
end subroutine bfgs_no_grad
