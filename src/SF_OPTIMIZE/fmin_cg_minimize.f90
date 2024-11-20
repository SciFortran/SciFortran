!+-------------------------------------------------------------------+
!     PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
!     Adapted from unkown minimize.f routine.
!     don't worry it works...
!+-------------------------------------------------------------------+
subroutine fmin_cgminimize_func(p,fcn,iter,fret,ftol,itmax,iverbose,mode,new_version,hh_par)
  real(8),dimension(:),intent(inout) :: p
  procedure(cgfit_func)              :: fcn
  integer                            :: iter
  real(8)                            :: fret
  real(8),optional                   :: ftol
  real(8)                            :: ftol_
  integer, optional                  :: itmax,mode
  integer                            :: itmax_,mode_,iprint_
  integer                            :: n
  real(8)                            :: f
  real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
  real(8)                            :: dfn,deps,hh
  integer                            :: iexit,itn
  logical,optional                   :: new_version
  logical                            :: new_version_
  logical,optional                   :: iverbose
  logical                            :: iverbose_
  real(8),optional                   :: hh_par
  !
  if(associated(func))nullify(func) ; func=>fcn
  !    
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  new_version_=.false.;if(present(new_version))new_version_=new_version
  iprint_=0;if(iverbose_)iprint_=1
  !
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iverbose_)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
  endif
  itmax_=1000
  if(present(itmax))then
     itmax_=itmax
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
  endif
  mode_ =1
  if(present(mode))then
     mode_=mode_
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
  endif
  !
  N=size(p)
  !
  !itmax is actually the number of calls to fcn
  itmax_ = 2*N*itmax_
  !
  allocate(x(N),g(N),h(N*N),w(100*N),xprmt(N))
  dfn=-0.5d0
  hh = 1d-5 ; if(present(hh_par))hh = hh_par
  iexit=0
  !set initial point
  x=p
  xprmt=abs(p)+1.d-15
  select case(new_version_)
  case (.true.)
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: using New version:"
     call minimize_sascha_(fcn_,n,x,f,g,h,w,&
          dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
  case(.false.)
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: using Old version:"
     call minimize_krauth_(fcn_,n,x,f,g,h,w,&
          dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
  end select
  !
  if( iexit == 0 )then
     write(*,"(A,I5)")"CG-minimize: IEXIT = 0, hessian matrix is not positive definite"
  elseif( iexit == 2 )then
    write(*,"(A,I5)")"CG-minimize: IEXIT = 2, G DX > 0"
  elseif( iexit == 3 )then
    write(*,"(A,I5)")"CG-minimize: IEXIT = 3, maximum number of calls to function exceeded"
  endif
  !set output variables
  iter=itn
  fret=f
  p=x
  deallocate(x,g,h,w,xprmt)
end subroutine fmin_cgminimize_func
subroutine fcn_(n,x,f)
  integer :: n
  real(8) :: x(n)
  real(8) :: f
  f=func(x)
end subroutine fcn_
!
subroutine fmin_cgminimize_sub(p,fcn,iter,fret,ftol,itmax,iverbose,mode,new_version,hh_par)
  real(8),dimension(:),intent(inout) :: p
  interface 
     subroutine fcn(n,x_,f_)
       integer                       :: n
       real(8),dimension(n)          :: x_
       real(8)                       :: f_
     end subroutine fcn
  end interface
  integer                            :: iter
  real(8)                            :: fret
  real(8),optional                   :: ftol
  real(8)                            :: ftol_
  integer, optional                  :: itmax,mode
  integer                            :: itmax_,mode_,iprint_
  integer                            :: n
  real(8)                            :: f
  real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
  real(8)                            :: dfn,deps,hh
  integer                            :: iexit,itn
  logical,optional                   :: new_version
  logical                            :: new_version_
  logical,optional                   :: iverbose
  logical                            :: iverbose_
  real(8),optional                   :: hh_par
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  new_version_=.false.;if(present(new_version))new_version_=new_version
  iprint_=0;if(iverbose_)iprint_=1
  !
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iverbose_)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
  endif
  itmax_=1000
  if(present(itmax))then
     itmax_=itmax
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
  endif
  mode_ =1
  if(present(mode))then
     mode_=mode_
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
  endif
  n=size(p)   
  !
  !itmax is actually the number of calls to fcn
  itmax_ = 2*n*itmax_
  !
  allocate(x(n),g(n),h(n*n),w(100*n),xprmt(n))
  dfn=-0.5d0
  hh = 1d-5 ; if(present(hh_par))hh = hh_par
  iexit=0
  !set initial point
  x=p
  xprmt=abs(p)+1.d-15
  select case(new_version_)
  case (.true.)
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: using New version:"
     call minimize_sascha_(fcn,n,x,f,g,h,w,&
          dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
  case (.false.)
     if(iverbose_)write(*,"(A,I5)")"CG-minimize: using Old version:"
     call minimize_krauth_(fcn,n,x,f,g,h,w,&
          dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
  end select
  !
  if( iexit == 0 )then
     write(*,"(A,I5)")"CG-minimize: IEXIT = 0, hessian matrix is not positive definite"
  elseif( iexit == 2 )then
    write(*,"(A,I5)")"CG-minimize: IEXIT = 2, G DX > 0"
  elseif( iexit == 3 )then
    write(*,"(A,I5)")"CG-minimize: IEXIT = 3, maximum number of calls to function exceeded"
  endif
  !set output variables
  iter=itn
  fret=f
  p=x
  deallocate(x,g,h,w,xprmt)
end subroutine fmin_cgminimize_sub
