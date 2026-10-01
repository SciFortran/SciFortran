
!+-------------------------------------------------------------------+
!PURPOSE  : 
! Given a function f, and given a bracketing triplet of abscissas
! ax, bx, cx (such that bx is between ax and cx, and f(bx) is less
! than both f(ax) and f(cx)), this routine isolates the minimum to a
! fractional precision of about tol using Brent’s method. The abscissa
! of the minimum is returned as xmin, and the minimum function value
! is returned as brent, the returned function value.
! Parameters: Maximum allowed number of iterations; golden ratio;
! and a small number that protects against trying to achieve
! fractional accuracy for a minimum that happens to be exactly zero.
!+-------------------------------------------------------------------+
subroutine brent(func,xmin,brack,tol,niter)
  !This subroutine finds the minimum of a function of one variable, :f:var:`func`, with the Brent method (golden section
  !search and parabolic interpolation). The abscissa of the minimum is returned in :f:var:`xmin`. The minimum must be
  !bracketed: :f:var:`brack` contains either a bracketing triplet :code:`(a,x,b)`, with :code:`a<x<b` and :code:`func(x)`
  !lower than both :code:`func(a)` and :code:`func(b)`, or two starting points, from which a bracketing triplet is searched
  !downhill with :f:func_inline:`bracket`. If :f:var:`brack` is absent the points 0 and 1 are used as starting points. The
  !program stops if :f:var:`brack` has a single element.
  !
  !@param func: function whose minimum is searched, real(8) function func(x) with real(8) x
  !@param xmin: abscissa of the minimum
  !@param brack: optional: 2 starting points or 3 bracketing points, default [0,1]
  !@param tol: optional: fractional precision of the minimum, default 1d-9
  !@param niter: optional: maximum number of iterations, default 200
  !
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
  end interface
  real(8),intent(inout)         :: xmin
  real(8),dimension(:),optional :: brack
  real(8),optional              :: tol
  integer,optional              :: niter
  real(8)                       :: tol_
  integer                       :: niter_
  integer                       :: iter
  real(8)                       :: ax,xx,bx,fa,fx,fb,fret
  !
  tol_=1d-9;if(present(tol))tol_=tol
  Niter_=200;if(present(Niter))Niter_=Niter
  !
  if(present(brack))then
     select case(size(brack))
     case(1)
        stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
     case(2)
        ax = brack(1)
        xx = brack(2)
        call bracket(ax,xx,bx,fa,fx,fb,func)
     case (3)
        ax = brack(1)
        xx = brack(2)
        bx = brack(3)
     end select
  else
     ax=0d0
     xx=1d0
     call bracket(ax,xx,bx,fa,fx,fb,func)
  endif
  fret=brent_optimize(ax,xx,bx,func,tol_,niter_,xmin)
end subroutine brent
!

!
function brent_optimize(ax,bx,cx,func,tol,itmax,xmin)
  real(8), intent(in)  :: ax,bx,cx,tol
  real(8), intent(out) :: xmin
  real(8)              :: brent_optimize
  integer              :: itmax
  real(8), parameter   :: cgold=0.3819660d0,zeps=1.0d-3*epsilon(ax)
  integer              :: iter
  real(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
  end interface
  a=min(ax,cx)
  b=max(ax,cx)
  v=bx
  w=v
  x=v
  e=0.d0
  fx=func(x)
  fv=fx
  fw=fx
  do iter=1,itmax
     xm=0.5d0*(a+b)
     tol1=tol*abs(x)+zeps
     tol2=2.0*tol1
     if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
        xmin=x
        brent_optimize=fx
        return
     end if
     if (abs(e) > tol1) then
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.d0*(q-r)
        if (q > 0.d0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5d0*q*etemp) .or. &
             p <= q*(a-x) .or. p >= q*(b-x)) then
           e=merge(a-x,b-x, x >= xm )
           d=cgold*e
        else
           d=p/q
           u=x+d
           if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
        end if
     else
        e=merge(a-x,b-x, x >= xm )
        d=cgold*e
     end if
     u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
     fu=func(u)
     if (fu <= fx) then
        if (u >= x) then
           a=x
        else
           b=x
        end if
        call shft(v,w,x,u)
        call shft(fv,fw,fx,fu)
     else
        if (u < x) then
           a=u
        else
           b=u
        end if
        if (fu <= fw .or. w == x) then
           v=w
           fv=fw
           w=u
           fw=fu
        else if (fu <= fv .or. v == x .or. v == w) then
           v=u
           fv=fu
        end if
     end if
  end do
  !pause 'brent: exceed maximum iterations'
contains
  subroutine shft(a,b,c,d)
    real(8), intent(out) :: a
    real(8), intent(inout) :: b,c
    real(8), intent(in) :: d
    a=b
    b=c
    c=d
  end subroutine shft
end function brent_optimize






!+-------------------------------------------------------------------+
!PURPOSE  : 
!  Given a function f and its derivative function df, and given a
!  bracketing triplet of abscissas ax, bx, cx [such that bx is between
!  ax and cx, and f(bx) is less than both f(ax) and f(cx)], this
!  routine isolates the minimum to a fractional precision of about
!  tol using a modification of Brent’s method that uses derivatives.
!  The abscissa of the minimum is returned as xmin, and the minimum
!  function value is returned as dbrent, the returned function value.
!+-------------------------------------------------------------------+ 
subroutine dbrent_wgrad(func,dfunc,xmin,brack,tol,niter)
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
     function dfunc(x)
       real(8) :: x
       real(8) :: dfunc
     end function dfunc
  end interface
  real(8),intent(inout)         :: xmin   ! abscissa of the minimum
  real(8),dimension(:),optional :: brack  ! optional: 2 starting points or 3 bracketing points, default [0,1]
  real(8),optional              :: tol    ! optional: fractional precision of the minimum, default 1d-9
  integer,optional              :: niter  ! optional: maximum number of iterations, default 200
  real(8)                       :: tol_
  integer                       :: niter_
  integer                       :: iter
  real(8)                       :: ax,xx,bx,fa,fx,fb,fret
  !
  tol_=1d-9;if(present(tol))tol_=tol
  Niter_=200;if(present(Niter))Niter_=Niter
  !
  if(present(brack))then
     select case(size(brack))
     case(1)
        stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
     case(2)
        ax = brack(1)
        xx = brack(2)
        call bracket(ax,xx,bx,fa,fx,fb,func)
     case (3)
        ax = brack(1)
        xx = brack(2)
        bx = brack(3)
     end select
  else
     ax=0d0
     xx=1d0
     call bracket(ax,xx,bx,fa,fx,fb,func)
  endif
  fret=dbrent_optimize(ax,xx,bx,func,dfunc,tol_,niter_,xmin)
end subroutine dbrent_wgrad
!
subroutine dbrent_nograd(func,xmin,brack,tol,niter)
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
  end interface
  real(8),intent(inout)         :: xmin   ! abscissa of the minimum
  real(8),dimension(:),optional :: brack  ! optional: 2 starting points or 3 bracketing points, default [0,1]
  real(8),optional              :: tol    ! optional: fractional precision of the minimum, default 1d-9
  integer,optional              :: niter  ! optional: maximum number of iterations, default 200
  real(8)                       :: tol_
  integer                       :: niter_
  integer                       :: iter
  real(8)                       :: ax,xx,bx,fa,fx,fb,fret
  !
  tol_=1d-9;if(present(tol))tol_=tol
  Niter_=200;if(present(Niter))Niter_=Niter
  !
  if(present(brack))then
     select case(size(brack))
     case(1)
        stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
     case(2)
        ax = brack(1)
        xx = brack(2)
        call bracket(ax,xx,bx,fa,fx,fb,func)
     case (3)
        ax = brack(1)
        xx = brack(2)
        bx = brack(3)
     end select
  else
     ax=0d0
     xx=1d0
     call bracket(ax,xx,bx,fa,fx,fb,func)
  endif
  fret=dbrent_optimize(ax,xx,bx,func,dfunc,tol_,niter_,xmin)
contains
  function dfunc(x)
    real(8) :: x
    real(8) :: dfunc
    call fgradient_func(func,x,dfunc)
  end function dfunc
  !
  subroutine fgradient_func(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       function funcv(x)
         real(8) :: x
         real(8) :: funcv
       end function funcv
    end interface
    integer          ::  n
    real(8),intent(in) ::  x
    real(8)            ::  x_
    real(8)          ::  fvec
    real(8)          ::  fjac
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    x_ = x
    eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    fvec = funcv(x_)
    temp = x_
    h    = eps*abs(temp)
    if(h==0d0) h = eps
    x_   = temp + h
    wa1  = funcv(x_)
    x_   = temp
    fjac = (wa1 - fvec)/h
  end subroutine fgradient_func
end subroutine dbrent_nograd



function dbrent_optimize(ax,bx,cx,func,fjac,tol,itmax,xmin) result(dbrent)
  real(8),intent(in)  :: ax,bx,cx,tol
  real(8),intent(out) :: xmin
  real(8)             :: dbrent
  integer             :: itmax
  real(8), parameter  :: zeps=1.d-3*epsilon(ax)
  integer             :: iter
  real(8)             :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2
  real(8)             :: u,u1,u2,v,w,x,xm
  logical             :: ok1,ok2
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
     function fjac(x)
       real(8) :: x
       real(8) :: fjac
     end function fjac
  end interface
  a=min(ax,cx)
  b=max(ax,cx)
  v=bx
  w=v
  x=v
  e=0.d0
  fx=func(x)
  fv=fx
  fw=fx
  dx=fjac(x)
  dv=dx
  dw=dx
  do iter=1,ITMAX
     xm=0.5d0*(a+b)
     tol1=tol*abs(x)+ZEPS
     tol2=2.0d0*tol1
     if (abs(x-xm) <= (tol2-0.5d0*(b-a))) exit
     if (abs(e) > tol1) then
        d1=2.0d0*(b-a)
        d2=d1
        if (dw /= dx) d1=(w-x)*dx/(dx-dw)
        if (dv /= dx) d2=(v-x)*dx/(dx-dv)
        u1=x+d1
        u2=x+d2
        ok1=((a-u1)*(u1-b) > 0.d0) .and. (dx*d1 <= 0.d0)
        ok2=((a-u2)*(u2-b) > 0.d0) .and. (dx*d2 <= 0.d0)
        olde=e
        e=d
        if (ok1 .or. ok2) then
           if (ok1 .and. ok2) then
              d=merge(d1,d2, abs(d1) < abs(d2))
           else
              d=merge(d1,d2,ok1)
           end if
           if (abs(d) <= abs(0.5d0*olde)) then
              u=x+d
              if (u-a < tol2 .or. b-u < tol2) &
                   d=sign(tol1,xm-x)
           else
              e=merge(a,b, dx >= 0.d0)-x
              d=0.5d0*e
           end if
        else
           e=merge(a,b, dx >= 0.d0)-x
           d=0.5d0*e
        end if
     else
        e=merge(a,b, dx >= 0.d0)-x
        d=0.5d0*e
     end if
     if (abs(d) >= tol1) then
        u=x+d
        fu=func(u)
     else
        u=x+sign(tol1,d)
        fu=func(u)
        if (fu > fx) exit
     end if
     du=fjac(u)
     if (fu <= fx) then
        if (u >= x) then
           a=x
        else
           b=x
        end if
        call mov3(v,fv,dv,w,fw,dw)
        call mov3(w,fw,dw,x,fx,dx)
        call mov3(x,fx,dx,u,fu,du)
     else
        if (u < x) then
           a=u
        else
           b=u
        end if
        if (fu <= fw .or. w == x) then
           call mov3(v,fv,dv,w,fw,dw)
           call mov3(w,fw,dw,u,fu,du)
        else if (fu <= fv .or. v == x .or. v == w) then
           call mov3(v,fv,dv,u,fu,du)
        end if
     end if
  end do
  if (iter > ITMAX) stop 'dbrent: exceeded maximum iterations'
  xmin=x
  dbrent=fx
contains
  !bl
  subroutine mov3(a,b,c,d,e,f)
    real(8), intent(in) :: d,e,f
    real(8), intent(out) :: a,b,c
    a=d
    b=e
    c=f
  end subroutine mov3
end function dbrent_optimize




!+-------------------------------------------------------------------+
!PURPOSE  : 
!     Given a function FUNC, and given distinct initial points AX and BX,
!     this routine searches in the downhill direction (defined by the 
!     function as evaluated at the initial points) and returns new points
!     AX, BX, CX which bracket a minimum of the function.  
!     Also returned are the function values at the three points, 
!     FA, FB, and FC.
!+-------------------------------------------------------------------+
subroutine bracket(ax,bx,cx,fa,fb,fc,func)
  !Given a function :f:var:`func` and two distinct initial points :f:var:`ax` and :f:var:`bx`, this subroutine searches in the
  !downhill direction, defined by the function evaluated at the initial points, and returns three points :f:var:`ax`,
  !:f:var:`bx`, :f:var:`cx` which bracket a minimum of the function, i.e. :code:`bx` is between :code:`ax` and :code:`cx` and
  !:code:`func(bx)` is lower than both :code:`func(ax)` and :code:`func(cx)`. The values of the function at the three points
  !are returned too.
  !
  !@param ax: in: first initial point; out: first point of the bracketing triplet
  !@param bx: in: second initial point; out: middle point of the bracketing triplet
  !@param cx: out: last point of the bracketing triplet
  !@param fa: out: func(ax)
  !@param fb: out: func(bx)
  !@param fc: out: func(cx)
  !@param func: function whose minimum is bracketed, real(8) func(x) with real(8) x
  !
  real(8), intent(inout) :: ax,bx
  real(8), intent(out) :: cx,fa,fb,fc
  !...the first parameter is the default ratio by which successive intervals
  !   are magnified; the second is the maximum magnification allowed for a
  !   parabolic-fit step
  real(8), parameter :: gold=1.618034d0,glimit=100.d0,tiny=1.d-20
  real(8) :: fu,q,r,u,ulim
  interface
     function func(x)
       real(8) :: x
       real(8) :: func
     end function func
  end interface
  fa=func(ax)
  fb=func(bx)
  if (fb > fa) then
     call swap(ax,bx)
     call swap(fa,fb)
  end if
  cx=bx+GOLD*(bx-ax)
  fc=func(cx)
  do
     if (fb < fc) RETURN
     r=(bx-ax)*(fb-fc)
     q=(bx-cx)*(fb-fa)
     u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(max(abs(q-r),TINY),q-r))
     ulim=bx+GLIMIT*(cx-bx)
     if ((bx-u)*(u-cx) > 0.d0) then
        fu=func(u)
        if (fu < fc) then
           ax=bx
           fa=fb
           bx=u
           fb=fu
           RETURN
        else if (fu > fb) then
           cx=u
           fc=fu
           RETURN
        end if
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     else if ((cx-u)*(u-ulim) > 0.d0) then
        fu=func(u)
        if (fu < fc) then
           bx=cx
           cx=u
           u=cx+GOLD*(cx-bx)
           call shft(fb,fc,fu,func(u))
        end if
     else if ((u-ulim)*(ulim-cx) >= 0.d0) then
        u=ulim
        fu=func(u)
     else
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     end if
     call shft(ax,bx,cx,u)
     call shft(fa,fb,fc,fu)
  end do
contains
  subroutine swap(a,b)
    real(8), intent(inout) :: a,b
    real(8) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap
  !-------------------
  subroutine shft(a,b,c,d)
    real(8), intent(out) :: a
    real(8), intent(inout) :: b,c
    real(8), intent(in) :: d
    a=b
    b=c
    c=d
  end subroutine shft
end subroutine bracket
