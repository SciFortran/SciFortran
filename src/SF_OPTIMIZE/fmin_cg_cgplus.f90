subroutine fmin_cgplus_df(p,func,fjac,iter,fret,ftol,itmax,imethod,iverb1,iverb2)
  real(8),dimension(:),intent(inout) :: p
  integer                            :: N,i
  interface 
     function func(a)
       real(8),dimension(:)          ::  a
       real(8)                       ::  func
     end function func
     function fjac(a)
       real(8),dimension(:)          :: a
       real(8),dimension(size(a))    :: fjac
     end function fjac
  end interface
  integer,intent(out)                :: iter
  real(8)                            :: fret
  real(8),optional                   :: ftol
  real(8)                            :: ftol_
  integer, optional                  :: itmax,iverb1,iverb2,imethod
  integer                            :: itmax_
  real(8),allocatable,dimension(:)   :: x,g,d,gold,w
  real(8)                            :: f,eps,tlev
  logical                            :: finish
  integer                            :: iprint(2),iflag,method
  integer                            :: nfun,irest
  iprint(1)= -1;if(present(iverb1))iprint(1)=iverb1
  iprint(2)= 0;if(present(iverb2))iprint(2)=iverb2
  method   = 2;if(present(imethod))method=imethod
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG+: ftol updated to:",ftol
  endif
  itmax_=1000
  if(present(itmax))then
     itmax_=itmax
     if(iprint(1)>=0)write(*,"(A,I5)")"CG+: itmax updated to:",itmax
  endif
  n     = size(p)
  finish= .false. 
  irest = 1
  allocate(x(n),g(n),d(n),gold(n),w(n))
  x     = p
  iflag = 0
  fgloop: do
     !calculate the function and gradient values here
     f = func(x)
     g = fjac(x)
     cgloop: do
        !call the CG code
        !iflag= 0 : successful termination
        !       1 : return to evaluate f and g
        !       2 : return with a new iterate, try termination test
        !      -i : error
        call cgfam(n,x,f,g,d,gold,iprint,ftol_,w,iflag,irest,method,finish,iter,nfun)
        if(iflag <= 0 .OR. iter > itmax_) exit fgloop
        if(iflag == 1) cycle fgloop
        if(iflag == 2) then
           ! termination test.  the user may replace it by some other test. however,
           ! the parameter 'finish' must be set to 'true' when the test is satisfied.
           tlev= ftol_*(1.d0 + dabs(f))
           i=0
           iloop: do 
              i=i+1
              if(i > n) then
                 finish = .true. 
                 cycle cgloop
              endif
              if(dabs(g(i)) > tlev) then
                 cycle cgloop
              else
                 cycle iloop
              endif
           enddo iloop
        endif
     enddo cgloop
  enddo fgloop
  p=x
  fret=f
  if(iflag<0)stop "CG+ error: iflag < 0"
  if(iprint(1)>=0.AND.iter>=itmax_)write(0,*)"CG+ exit with iter >= itmax"
end subroutine fmin_cgplus_df

subroutine fmin_cgplus_f(p,fcn,iter,fret,ftol,itmax,imethod,deps,iverb1,iverb2)
  real(8),dimension(:),intent(inout) :: p
  integer                            :: N,i
  ! interface 
  !    function fcn(a)
  !      real(8),dimension(:)          ::  a
  !      real(8)                       ::  fcn
  !    end function fcn
  ! end interface
  procedure(cgfit_func)              :: fcn
  integer,intent(out)                :: iter
  real(8)                            :: fret
  real(8),optional                   :: ftol,deps
  real(8)                            :: ftol_,deps_
  integer, optional                  :: itmax,iverb1,iverb2,imethod
  integer                            :: itmax_
  real(8),allocatable,dimension(:)   :: x,g,d,gold,w
  real(8)                            :: f,eps,tlev
  logical                            :: finish
  integer                            :: iprint(2),iflag,method
  integer                            :: nfun,irest
  if(associated(func))nullify(func) ; func=>fcn
  if(associated(fjac))nullify(fjac) ; fjac=>dfcn
  iprint(1)= -1;if(present(iverb1))iprint(1)=iverb1
  iprint(2)= 0;if(present(iverb2))iprint(2)=iverb2
  method   = 2;if(present(imethod))method=imethod
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG+: ftol updated to:",ftol
  endif
  itmax_=1000
  if(present(itmax))then
     itmax_=itmax
     if(iprint(1)>=0)write(*,"(A,I5)")"CG+: itmax updated to:",itmax
  endif
  deps_=0d0
  if(present(deps))then
     deps_=deps
     if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG: deps update to:",deps
  endif
  df_eps = deps_
  !
  n     = size(p)
  finish= .false. 
  irest = 1
  allocate(x(n),g(n),d(n),gold(n),w(n))
  x     = p
  iflag = 0
  fgloop: do
     !calculate the function and gradient values here
     f = func(x)
     g = fjac(x)
     cgloop: do
        !call the CG code
        !iflag= 0 : successful termination
        !       1 : return to evaluate f and g
        !       2 : return with a new iterate, try termination test
        !      -i : error
        call cgfam(n,x,f,g,d,gold,iprint,ftol_,w,iflag,irest,method,finish,iter,nfun)
        if(iflag <= 0 .OR. iter > itmax_) exit fgloop
        if(iflag == 1) cycle fgloop
        if(iflag == 2) then
           ! termination test.  the user may replace it by some other test. however,
           ! the parameter 'finish' must be set to 'true' when the test is satisfied.
           tlev= ftol_*(1.d0 + dabs(f))
           i=0
           iloop: do 
              i=i+1
              if(i > n) then
                 finish = .true. 
                 cycle cgloop
              endif
              if(dabs(g(i)) > tlev) then
                 cycle cgloop
              else
                 cycle iloop
              endif
           enddo iloop
        endif
     enddo cgloop
  enddo fgloop
  p=x
  fret=f
  if(iflag<0)stop "CG+ error: iflag < 0"
  if(iprint(1)>=0.AND.iter>=itmax_)write(0,*)"CG+ exit with iter >= itmax"
  nullify(func)
  nullify(fjac)
  return
end subroutine fmin_cgplus_f
!
function dfcn(p) 
  real(8),dimension(:)       :: p
  real(8),dimension(size(p)) :: dfcn
  dfcn=func_func_jacobian(func,size(p),p)
end function dfcn
