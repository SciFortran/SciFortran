!ENTIRELY BASED ON THE WORK OF 
! Jacob Williams https://github.com/jacobwilliams/quadrature-fortran
!SEE There for a more object-oriented approach and further information.
!LICENSE CAN BE FOUND HERE:
!https://github.com/jacobwilliams/quadrature-fortran/blob/master/LICENSE
!
!LIBRARY EXTENDED TO DEAL WITH:
! 1D, 2D sampled functions (dble arrays)
! vectorial functions of vectorial variables F:R^n-->R^m
! this required to slightly change *dgauss_generic routine and gX=6,8,10,12,14
! to handle multiple functions simultaneously.
MODULE GAUSS_QUADRATURE
  implicit none
  private
  real(8),parameter :: zero      = 0.0d0
  real(8),parameter :: one_half  = 0.5d0
  real(8),parameter :: one       = 1.0d0
  real(8),parameter :: two       = 2.0d0
  real(8),parameter :: three     = 3.0d0
  real(8),parameter :: four      = 4.0d0

  type :: integration_type
     real(8)                                         :: val
     real(8)                                         :: tol
     real(8)                                         :: a
     real(8)                                         :: b
     integer                                         :: dim=1
     !
     procedure(gauss_func_method),pointer            :: g  => null()
     procedure(f_x),pointer,nopass                   :: fx => null()
     procedure(f_xvec),pointer,nopass                :: fxvec => null()
     !
     type(integration_type),dimension(:),allocatable :: ivec
  end type integration_type


  abstract interface
     !Fvec(x): R-->R^m integrand
     function f_x(x,m) result(f)
       real(8) :: x
       integer :: m
       real(8) :: f(m)
     end function f_x

     !F(xvec): R^n-->R^m integrand  
     function f_xvec(x,m) result(f)
       real(8) :: x(:)          !n
       integer :: m
       real(8) :: f(m)
     end function f_xvec

     !Gaussian quadrature method to be used in integration
     function gauss_func_method(self, x, h, m) result(f)
       import :: integration_type
       class(integration_type),intent(inout) :: self
       real(8), intent(in)                   :: x
       real(8), intent(in)                   :: h
       integer                               :: m
       real(8)                               :: f(m)
     end function gauss_func_method
  end interface

  interface gauss_quad
     module procedure :: integrate_1d_func_main
     module procedure :: integrate_nd_func_main
     module procedure :: integrate_1d_func_1
     module procedure :: integrate_nd_func_1
     module procedure :: integrate_1d_sample
     module procedure :: integrate_2d_sample
  end interface gauss_quad

  interface integrate
     module procedure :: integrate_1d_func_main
     module procedure :: integrate_nd_func_main
     module procedure :: integrate_1d_func_1
     module procedure :: integrate_nd_func_1
     module procedure :: integrate_1d_sample
     module procedure :: integrate_2d_sample
  end interface integrate

  public :: gauss_quad
  public :: integrate



  !FUNCTION INTERPOLATION TYPE: used to integrate sampled functions
  type finter1d_type
     real(8),allocatable,dimension(:)   :: X
     real(8),allocatable,dimension(:)   :: F
     real(8),allocatable,dimension(:)   :: G
     integer                            :: N=0
     integer                            :: Imin=0,Imax=0
     logical                            :: status=.false.
  end type finter1d_type

  type finter2d_type
     real(8),allocatable,dimension(:)   :: X
     real(8),allocatable,dimension(:)   :: y
     real(8),allocatable,dimension(:,:) :: F 
     integer                            :: N=0
     integer                            :: Imin=0,Imax=0
     integer                            :: Jmin=0,Jmax=0
     logical                            :: status=.false.
  end type finter2d_type



contains



  subroutine init_type(self,xl,xu,tol,method)
    type(integration_type)    :: self     !! for the 1d integration
    real(8),intent(in)        :: xl       !! x integration lower bound
    real(8),intent(in)        :: xu       !! x integration upper bound
    real(8),intent(in)        :: tol      !! error tolerance for dx integration
    integer,intent(in)        :: method   !! quadrature method to use for x
    !
    select case (method)
    case(6);  self%g => g6
    case(8);  self%g => g8
    case(10); self%g => g10
    case(12); self%g => g12
    case(14); self%g => g14
    case default
       error stop 'invalid quadrature method in initialize integration: method=6,8,10,12,14'
    end select
    !
    self%tol    =  tol             !tolerance
    self%a      =  xl              !lower bound
    self%b      =  xu              !upper bound
  end subroutine init_type




  !##################################################################
  !##################################################################
  !               FUNCTIONS: nD R^n --> R^m
  !              m=1 is available as a simplified interface
  !##################################################################
  !##################################################################

  subroutine integrate_1d_func_main(m,fx,xl,xu,ans,tol,method,ierr,err)
    integer                            :: m
    procedure(f_x)                     :: fx       !! 1d function: f(x)
    real(8),intent(in)                 :: xl       !! x integration lower bound
    real(8),intent(in)                 :: xu       !! x integration upper bound
    real(8),intent(inout),dimension(m) :: ans
    !
    real(8),intent(in),optional        :: tol     !! error tolerance for dx integration
    integer,intent(in),optional        :: method  !! quadrature method to use for x
    integer,intent(out),optional       :: ierr
    real(8),intent(out),optional       :: err
    !
    type(integration_type)             :: self     !! for the 1d integration
    !
    real(8)                            :: tol_
    integer                            :: method_
    integer                            :: ierr_
    real(8)                            :: err_
    !
    tol_   = 1d-9 ; if(present(tol))tol_ = tol
    method_= 10   ; if(present(method))method_=method
    !
    self%fx => fx
    call init_type(self,xl,xu,tol_,method)
    call dgauss_generic(self,ans, ierr_,err_,m)
    select case(ierr_)
    case(-1)
       stop "Lower Bound ~ Upper Bound => Ans set to zero"
    case(2)
       stop "Dgauss_Generic ABNORMAL EXIT CODE: ANS does not meet requested tolerance"
    case default
       continue
    end select
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
  end subroutine integrate_1d_func_main


  subroutine integrate_1d_func_1(func,xl,xu,ans,tol,method,ierr,err)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    real(8),intent(in)                 :: xl       !! x integration lower bound
    real(8),intent(in)                 :: xu       !! x integration upper bound
    real(8),intent(inout)              :: ans
    !
    real(8),intent(in),optional        :: tol     !! error tolerance for dx integration
    integer,intent(in),optional        :: method  !! quadrature method to use for x
    integer,intent(out),optional       :: ierr
    real(8),intent(out),optional       :: err
    !
    real(8)                            :: tol_,ans_(1)
    integer                            :: method_
    integer                            :: ierr_
    real(8)                            :: err_
    tol_   = 1d-9 ; if(present(tol))tol_ = tol
    method_= 10   ; if(present(method))method_=method
    call integrate_1d_func_main(1,fx,xl,xu,ans_,tol_,method_,ierr_,err_)
    ans = ans_(1)
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
  contains
    function fx(x,m)
      real(8) :: x
      integer :: m
      real(8) :: fx(m)
      fx(1) = func(x)
    end function fx
  end subroutine integrate_1d_func_1


  subroutine integrate_nd_func_main(m,fxvec,xl,xu,ans,methods,method,tols,tol,ierr,err)
    integer                                         :: m
    procedure(f_xvec)                               :: fxvec    !! 2d function: f(x,y)
    real(8),dimension(:),intent(in)                 :: xl       !! integration lower bounds
    real(8),dimension(size(xl)),intent(in)          :: xu       !! integration upper bounds
    real(8),intent(inout),dimension(m)              :: ans
    !
    integer,dimension(size(xl)),intent(in),optional :: methods  !! quadrature methods nD-->nD to use
    integer,intent(in),optional                     :: method   !! quadrature method   1-->nD to use
    real(8),dimension(size(xl)),intent(in),optional :: tols     !! error tolerances nD-->nD for dx integration
    real(8),intent(in),optional                     :: tol     !! error tolerance   1-->nD for dx integration
    integer,intent(out),optional                    :: ierr
    real(8),intent(out),optional                    :: err
    !
    real(8),dimension(size(xl))                     :: tol_
    integer,dimension(size(xl))                     :: method_
    integer                                         :: ierr_
    real(8)                                         :: err_
    !
    integer                                         :: i
    type(integration_type)                          :: self     !! for the 1d integration
    !
    tol_   = 1d-9 ; if(present(tol))tol_ = tol; if(present(tols))tol_ = tols
    method_= 10   ; if(present(method))method_=method; if(present(methods))method_=methods
    !
    self%fxvec => fxvec  !the user-defined f(x,y) function to integrate
    self%dim   = size(xl)
    allocate(self%ivec(self%dim))
    select case (self%dim)
    case default
       stop "INTEGRATE_ND ERROR: calling multi-dimensional quadrature for dimension N>6. Ask developer or DIY."
    case(1)
       stop "INTEGRATE_ND ERROR: calling multi-dimensional quadrature for 1d problem. Use scalar version."
    case (2)
       self%ivec(1)%fx => f_of_x1
       self%ivec(2)%fx => f_of_x2
    case (3)
       self%ivec(1)%fx => f_of_x1
       self%ivec(2)%fx => f_of_x2
       self%ivec(3)%fx => f_of_x3
    case (4)
       self%ivec(1)%fx => f_of_x1
       self%ivec(2)%fx => f_of_x2
       self%ivec(3)%fx => f_of_x3
       self%ivec(4)%fx => f_of_x4
    case (5)
       self%ivec(1)%fx => f_of_x1
       self%ivec(2)%fx => f_of_x2
       self%ivec(3)%fx => f_of_x3
       self%ivec(4)%fx => f_of_x4
       self%ivec(5)%fx => f_of_x5
    case (6)
       self%ivec(1)%fx => f_of_x1
       self%ivec(2)%fx => f_of_x2
       self%ivec(3)%fx => f_of_x3
       self%ivec(4)%fx => f_of_x4
       self%ivec(5)%fx => f_of_x5
       self%ivec(6)%fx => f_of_x6
    end select
    !
    do i=1,self%dim
       call init_type(self%ivec(i),xl(i),xu(i),tol_(i),method_(i))
    enddo
    !
    call dgauss_generic(self%ivec(self%dim),ans,ierr_,err_,m)
    select case(ierr_)
    case(-1)
       stop "Lower Bound ~ Upper Bound => Ans set to zero"
    case(2)
       stop "Dgauss_Generic ABNORMAL EXIT CODE: ANS does not meet requested tolerance"
    case default
       continue
    end select
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
    !
  contains
    !
    function f_of_x1(x1,m) result(f)
      real(8)                       :: x1
      integer                       :: m
      real(8)                       :: f(m)
      real(8),dimension(size(xl)-1) :: xvec
      integer                       :: i
      xvec = (/( self%ivec(i)%val, i=2,size(xl) )/)
      f = self%fxvec([x1,xvec],m)
    end function f_of_x1
    !
    function f_of_x2(x2,m) result(f)
      real(8) :: x2
      integer :: m
      real(8) :: f(m)
      self%ivec(2)%val = x2
      call dgauss_generic(self%ivec(1),f, ierr_, err_,m)
    end function f_of_x2
    !
    function f_of_x3(x3,m) result(f)
      real(8) :: x3
      integer :: m
      real(8) :: f(m)
      self%ivec(3)%val = x3
      call dgauss_generic(self%ivec(2),f, ierr_, err_,m)
    end function f_of_x3
    !
    function f_of_x4(x4,m) result(f)
      real(8) :: x4
      integer :: m
      real(8) :: f(m)
      self%ivec(4)%val = x4
      call dgauss_generic(self%ivec(3),f, ierr_, err_,m)
    end function f_of_x4
    !
    function f_of_x5(x5,m) result(f)
      real(8) :: x5
      integer :: m
      real(8) :: f(m)
      self%ivec(5)%val = x5
      call dgauss_generic(self%ivec(4),f, ierr_, err_,m)
    end function f_of_x5
    !
    function f_of_x6(x6,m) result(f)
      real(8) :: x6
      integer :: m
      real(8) :: f(m)
      self%ivec(6)%val = x6
      call dgauss_generic(self%ivec(5),f, ierr_, err_,m)
    end function f_of_x6
    !
  end subroutine integrate_nd_func_main



  subroutine integrate_nd_func_1(func,xl,xu,ans,methods,method,tols,tol,ierr,err)
    interface
       function func(x)
         real(8),dimension(:) :: x
         real(8)              :: func
       end function func
    end interface
    real(8),dimension(:),intent(in)                 :: xl       !! integration lower bounds
    real(8),dimension(size(xl)),intent(in)          :: xu       !! integration upper bounds
    real(8),intent(inout)                           :: ans
    integer,dimension(size(xl)),intent(in),optional :: methods  !! quadrature methods nD-->nD to use
    integer,intent(in),optional                     :: method   !! quadrature method   1-->nD to use
    real(8),dimension(size(xl)),intent(in),optional :: tols     !! error tolerances nD-->nD for dx integration
    real(8),intent(in),optional                     :: tol     !! error tolerance   1-->nD for dx integration
    integer,intent(out),optional                    :: ierr
    real(8),intent(out),optional                    :: err
    real(8),dimension(size(xl))                     :: tol_,ans_
    integer,dimension(size(xl))                     :: method_
    integer                                         :: ierr_
    real(8)                                         :: err_
    !
    tol_   = 1d-9 ; if(present(tol))tol_ = tol; if(present(tols))tol_ = tols
    method_= 10   ; if(present(method))method_=method; if(present(methods))method_=methods
    call integrate_nd_func_main(1,fxvec,xl,xu,ans_,methods=method_,tols=tol_,ierr=ierr_,err=err_)
    ans = ans_(1)
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
    !
  contains
    !
    function fxvec(x,m)
      real(8),dimension(:) :: x
      integer              :: m
      real(8),dimension(m) :: fxvec
      fxvec(1) = func(x)
    end function fxvec
    !
  end subroutine integrate_nd_func_1






  !##################################################################
  !##################################################################
  !               SAMPLED FUNCTIONS: 1D + 2D R^n --> R
  !##################################################################
  !##################################################################
  subroutine integrate_1d_sample(fsample,xl,xu,ans,tol,method,ierr,err)
    real(8),dimension(:)             :: fsample       !! 1d array: f(x)
    real(8),intent(in)               :: xl       !! x integration lower bound
    real(8),intent(in)               :: xu       !! x integration upper bound
    real(8),intent(inout)            :: ans
    !
    real(8),intent(in),optional      :: tol     !! error tolerance for dx integration
    integer,intent(in),optional      :: method  !! quadrature method to use for x
    integer,intent(out),optional     :: ierr
    real(8),intent(out),optional     :: err
    !
    real(8)                          :: tol_,ans_(1)
    integer                          :: method_
    integer                          :: ierr_
    real(8)                          :: err_
    !
    type(finter1d_type)              :: Finterp
    real(8),dimension(size(fsample)) :: Xsample
    integer                          :: Lsample
    integer                          :: Ninterp
    !
    tol_   = 1d-9 ; if(present(tol))tol_ = tol
    method_= 10   ; if(present(method))method_=method
    !
    Lsample = size(fsample)
    Xsample = gauss_quad_linspace(xl,xu,Lsample)
    call init_finter_1d(Finterp,Xsample,fsample,method_)
    call integrate_1d_func_main(1,fx,xl,xu,ans_,tol_,method_,ierr_,err_)
    ans = ans_(1)
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
    !
  contains
    !
    function fx(x,m)
      real(8) :: x
      integer :: m
      real(8) :: fx(m)
      real(8) :: dy
      integer :: j,k,k0,k1
      integer :: n
      fx(1) = 0d0
      !
      N = Finterp%N    !order of polynomial interpolation
      j = locate(finterp%X(finterp%Imin:finterp%Imax),x)
      !
      k=max(j-(N-1)/2,1)
      k0=k
      if(k0 < finterp%Imin)k0=finterp%Imin
      k1=k+N+1
      if(k1 > finterp%Imax)then
         k1=finterp%Imax
         k0=k1-N-1
      endif
      call polint(finterp%X(k0:k1),finterp%F(k0:k1),x,fx(1),dy)
    end function fx
  end subroutine integrate_1d_sample


  subroutine integrate_2d_sample(fsample,xl,xu,ans,methods,method,tols,tol,ierr,err)
    real(8),dimension(:,:),intent(in)               :: fsample  !! 2d function: f(x,y)
    real(8),dimension(:),intent(in)                 :: xl       !! integration lower bounds
    real(8),dimension(size(xl)),intent(in)          :: xu       !! integration upper bounds
    real(8),intent(inout)                           :: ans
    !
    integer,dimension(size(xl)),intent(in),optional :: methods  !! quadrature methods nD-->nD to use
    integer,intent(in),optional                     :: method   !! quadrature method   1-->nD to use
    real(8),dimension(size(xl)),intent(in),optional :: tols     !! error tolerances nD-->nD for dx integration
    real(8),intent(in),optional                     :: tol     !! error tolerance   1-->nD for dx integration
    integer,intent(out),optional                    :: ierr
    real(8),intent(out),optional                    :: err
    !
    real(8),dimension(size(xl))                     :: tol_,ans_
    integer,dimension(size(xl))                     :: method_
    integer                                         :: ierr_
    real(8)                                         :: err_
    !
    type(finter2d_type)                             :: Finterp
    real(8),dimension(size(fsample,1))              :: Xsample
    real(8),dimension(size(fsample,2))              :: Ysample
    integer                                         :: Lsample1,Lsample2
    integer                                         :: Ninterp
    !
    tol_   = 1d-9 ; if(present(tol))tol_ = tol; if(present(tols))tol_ = tols
    method_= 10   ; if(present(method))method_=method; if(present(methods))method_=methods
    !
    Lsample1 = size(fsample,1)
    Lsample2 = size(fsample,2)
    Xsample  = gauss_quad_linspace(xl(1),xu(1),Lsample1)
    Ysample  = gauss_quad_linspace(xl(2),xu(2),Lsample2)
    call init_finter_2d(Finterp,Xsample,Ysample,Fsample,method_(1))
    call integrate_nd_func_main(1,fxvec,xl,xu,ans_,methods=method_,tols=tol_,ierr=ierr_,err=err_)
    ans = ans_(1)
    if(present(ierr))ierr=ierr_
    if(present(err))err=err_
    !
  contains
    !
    function fxvec(xvec,m)
      real(8),dimension(:) :: xvec
      integer              :: m
      real(8)              :: x,y
      real(8)              :: fxvec(m)
      real(8)              :: df
      integer              :: itmp,jtmp,kx,ky,k0x,k0y,k1x,k1y
      integer              :: n
      !
      x = xvec(1)
      y = xvec(2)
      !
      N=Finterp%N    !order of polynomial interpolation
      itmp=locate(Finterp%X(Finterp%Imin:Finterp%Imax),x)
      jtmp=locate(Finterp%Y(Finterp%Jmin:Finterp%Jmax),y)
      kx=max(itmp-(N-1)/2,1)
      ky=max(jtmp-(N-1)/2,1)
      k0x = kx ; if(k0x < Finterp%Imin)k0x=Finterp%Imin
      k0y = ky ; if(k0y < Finterp%Jmin)k0y=Finterp%Jmin
      k1x = kx+N+1
      if(k1x > Finterp%Imax)then         
         k1x=Finterp%Imax
         k0x=k1x-N-1
      endif
      k1y = ky+N+1
      if(k1y > Finterp%Jmax)then
         k1y=Finterp%Jmax
         k0y=k1y-N-1
      endif
      call polin2(Finterp%X(k0x:k1x),Finterp%Y(k0y:k1y),Finterp%F(k0x:k1x,k0y:k1y),x,y,fxvec(1),df)
    end function fxvec
    !
  end subroutine integrate_2d_sample




















  !  Integrate a real function of one variable over a finite
  !  interval using the specified adaptive algorithm.
  !  Intended primarily for high accuracy
  !  integration or integration of smooth functions.
  !  * SLATEC is public domain software: http://www.netlib.org/slatec/guide
  !  * Original sourcecode from: http://www.netlib.org/slatec/src/dgaus8.f
  !  * Jones, R. E., (SNLA) -- Original SLATEC code.
  !  * Jacob Williams : 1/20/2020 : refactored to modern Fortran and generalized.
  !@note This function is recursive. [It can call itself indirectly during double integration]
  !! pick a value for abs(error_tol) so that
  !! dtol < abs(error_tol) <= 1.0e-3 where dtol is the larger
  !! of 1.0e-18 and the real(8) unit roundoff d1mach(4).
  !! ans will normally have no more error than abs(error_tol)
  !! times the integral of the absolute value of fun(x).  usually,
  !! smaller values of error_tol yield more accuracy and require
  !! more function evaluations.
  recursive subroutine dgauss_generic(self, ans, ierr, err, m)
    type(integration_type),intent(inout) :: self
    real(8)                              :: lb         !! lower bound of the integration
    real(8)                              :: ub         !! upper bound of the integration
    real(8)                              :: error_tol  !! is a requested pseudorelative error tolerance.
    integer                              :: m
    real(8),intent(out),dimension(m)     :: ans        !! computed value of integral
    integer,intent(out)                  :: ierr       !! status code:
    real(8),intent(out)                  :: err        !! an estimate of the absolute error in `ans`.
    !
    real(8),parameter                    :: sq2      = sqrt(two)
    real(8),parameter                    :: ln2      = log(two)
    integer,parameter                    :: nlmn     = 1                   !! ??
    integer,parameter                    :: kmx      = 5000                !! ??
    integer,parameter                    :: kml      = 6                   !! ??
    real(8),parameter                    :: magic    = 0.30102000d0       !! ??
    integer,parameter                    :: iwork    = 60                  !! size of the work arrays. ?? Why 60 ??
    real(8),parameter                    :: bb       = radix(one)          !! machine constant
    real(8),parameter                    :: d1mach4  = bb**(1-digits(one)) !! machine constant
    real(8),parameter                    :: d1mach5  = log10(bb)           !! machine constant
    integer                              :: k,l,lmn,lmx,mxl,nbits,nib,nlmx
    real(8)                              :: ae(m),anib,area(m),c,ee(m),ef,eps,est(m),gl(m),glr(m),tol,merr(m)
    real(8),dimension(iwork)             :: aa,hh
    real(8),dimension(iwork,m)           :: gr,vl
    integer,dimension(iwork)             :: lr
    !
    lb = self%a
    ub = self%b
    error_tol = self%tol
    !
    ans = zero
    ierr = 1
    err = zero
    if (lb == ub) return
    aa = zero
    hh = zero
    vl = zero
    gr = zero
    lr = 0
    merr = zero
    k = digits(one)
    anib = d1mach5*k/magic
    nbits = anib
    nlmx = min(60,(nbits*5)/8)         ! ... is this the same 60 as iwork???
    lmx = nlmx
    lmn = nlmn
    if (ub /= zero) then
       if (sign(one,ub)*lb > zero) then
          c = abs(one-lb/ub)
          if (c <= 0.1d0) then
             if (c <= zero) return
             anib = one_half - log(c)/ln2
             nib = anib
             lmx = min(nlmx,nbits-nib-7)
             if (lmx < 1) then
                ! lb and ub are too nearly equal to allow
                ! normal integration [ans is set to zero]
                ierr = -1
                return
             end if
             lmn = min(lmn,lmx)
          end if
       end if
    end if
    tol = max(abs(error_tol),two**(5-nbits))/two
    if (error_tol == zero) tol = sqrt(d1mach4)
    eps = tol
    hh(1) = (ub-lb)/four
    aa(1) = lb
    lr(1) = 1
    l = 1
    est = self%g(aa(l)+two*hh(l),two*hh(l),m)
    k = 8
    area = abs(est)
    ef = one_half
    mxl = 0
    !compute refined estimates, estimate the error, etc.
    main : do
       gl = self%g(aa(l)+hh(l),hh(l),m)
       gr(l,:) = self%g(aa(l)+three*hh(l),hh(l),m)
       k = k + 16
       area = area + (abs(gl)+abs(gr(l,:))-abs(est))
       glr = gl + gr(l,:)
       ee = abs(est-glr)*ef
       ae = max(eps*area,tol*abs(glr))
       if(any(ee-ae > zero)) then
          !consider the left half of this level
          if (k > kmx) lmx = kml
          if (l >= lmx) then
             mxl = 1
          else
             l = l + 1
             eps = eps*one_half
             ef = ef/sq2
             hh(l) = hh(l-1)*one_half
             lr(l) = -1
             aa(l) = aa(l-1)
             est = gl
             cycle main
          end if
       end if
       merr = merr + (est-glr)
       if (lr(l) > 0) then
          !return one level
          ans = glr
          do
             if (l <= 1) exit main ! finished
             l = l - 1
             eps = eps*two
             ef = ef*sq2
             if (lr(l) <= 0) then
                vl(l,:) = vl(l+1,:) + ans
                est = gr(l-1,:)
                lr(l) = 1
                aa(l) = aa(l) + four*hh(l)
                cycle main
             end if
             ans = vl(l+1,:) + ans
          end do
       else
          !proceed to right half at this level
          vl(l,:) = glr
          est = gr(l-1,:)
          lr(l) = 1
          aa(l) = aa(l) + four*hh(l)
          cycle main
       end if
    end do main
    err = maxval(abs(merr))
    if ((mxl/=0) .and. any(abs(merr)-two*tol*area > 0)) ierr = 2 ! ans is probably insufficiently accurate
  end subroutine dgauss_generic




  !>  !### See also
  !  * Coefficients from:
  !    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php
  !  6-point method.
  recursive function g6(self, x, h, m) result(f)
    class(integration_type),intent(inout) :: self
    real(8), intent(in)                   :: x
    real(8), intent(in)                   :: h
    integer                               :: m
    real(8)                               :: f(m)
    !> abscissae:
    real(8),dimension(3),parameter ::  a = [   0.6612093864662645136613995950199053470064485643&
         &951700708145267058521834966071431009442864037464&
         &614564298883716392751466795573467722253804381723&
         &198010093367423918538864300079016299442625145884&
         &9024557188219703863032236201173523213570221879361&
         &8906974301231555871064213101639896769013566165126&
         &1150514997832d0,&
         &0.2386191860831969086305017216807119354186106301&
         &400213501813951645742749342756398422492244272573&
         &491316090722230970106872029554530350772051352628&
         &872175189982985139866216812636229030578298770859&
         &440976999298617585739469216136216592222334626416&
         &400139367778945327871453246721518889993399000945&
         &408150514997832d0,&
         &0.9324695142031520278123015544939946091347657377&
         &122898248725496165266135008442001962762887399219&
         &259850478636797265728341065879713795116384041921&
         &786180750210169211578452038930846310372961174632&
         &524612619760497437974074226320896716211721783852&
         &305051047442772222093863676553669179038880252326&
         &771150514997832d0 ]
    !> weights:
    real(8),dimension(3),parameter ::  w = [   0.36076157304813860756983351383771611166152189274&
         &674548228973924023714003783726171832096220198881&
         &934794311720914037079858987989027836432107077678&
         &721140858189221145027225257577711260007323688285&
         &916316028951118005174081368554707448247248610118&
         &325993144981721640242558677752676819993095031068&
         &73150514997832d0,&
         0.46791393457269104738987034398955099481165560576&
         &921053531162531996391420162039812703111009258479&
         &198230476626878975479710092836255417350295459356&
         &355927338665933648259263825590180302812735635025&
         &362417046193182590009975698709590053347408007463&
         &437682443180817320636917410341626176534629278889&
         &17150514997832d0,&
         0.17132449237917034504029614217273289352682250148&
         &404398239863543979894576054234015464792770542638&
         &866975211652206987440430919174716746217597462964&
         &922931803144845206713510916832108437179940676688&
         &721266924855699404815942932735702498405343382418&
         &236324411837461039120523911904421970357029774978&
         &12150514997832d0 ]
    f = h * ( w(1)*( self%fx(x-a(1)*h,m) + self%fx(x+a(1)*h,m) ) + &
         w(2)*( self%fx(x-a(2)*h,m) + self%fx(x+a(2)*h,m) ) + &
         w(3)*( self%fx(x-a(3)*h,m) + self%fx(x+a(3)*h,m) ) )
  end function g6


  !  This is the 8-point formula from the original SLATEC routine
  recursive function g8(self, x, h,m ) result(f)
    class(integration_type),intent(inout) :: self
    real(8), intent(in)                   :: x
    real(8), intent(in)                   :: h
    integer                               :: m
    real(8)                               :: f(m)
    !> abscissae:
    real(8),parameter ::   x1 = 0.18343464249564980493947614236018398066675781291297378231718847&
         &369920447422154211411606822371112335374526765876&
         &428676660891960125238768656837885699951606635681&
         &044755516171385019663858107642055323708826547494&
         &928123149612477646193635627706457164566131594051&
         &34052985058171969174306064445289638150514997832d0
    real(8),parameter ::   x2 = 0.52553240991632898581773904918924634904196424312039285775085709&
         &927245482076856127252396140019363198206190968292&
         &482526085071087937666387799398053953036682536311&
         &190182730324023600607174700061279014795875767562&
         &412888953366196435283308256242634705401842246036&
         &88817537938539658502113876953598879150514997832d0
    real(8),parameter ::   x3 = 0.79666647741362673959155393647583043683717173161596483207017029&
         &503921730567647309214715192729572593901919745345&
         &309730926536564949170108596027725620746216896761&
         &539350162903423256455826342053015458560600957273&
         &426035574157612651404288519573419337108037227831&
         &36113628137267630651413319993338002150514997832d0
    real(8),parameter ::   x4 = 0.96028985649753623168356086856947299042823523430145203827163977&
         &737242489774341928443943895926331226831042439281&
         &729417621023895815521712854793736422049096997004&
         &339826183266373468087812635533469278673596634808&
         &705975425476039293185338665681328688426134748962&
         &8923208763998895240977248938732425615051499783203d0
    !> weights:
    real(8),parameter ::   w1 = 0.36268378337836198296515044927719561219414603989433054052482306&
         &756668673472390667732436604208482850955025876992&
         &629670655292582155698951738449955760078620768427&
         &783503828625463057710075533732697147148942683287&
         &804318227790778467229655355481996014024877675059&
         &28976560993309027632737537826127502150514997832d0
    real(8),parameter ::   w2 = 0.31370664587788728733796220198660131326032899900273493769026394&
         &507495627194217349696169807623392855604942757464&
         &107780861624724683226556160568906242764697589946&
         &225031187765625594632872220215204316264677947216&
         &038226012952768986525097231851579983531560624197&
         &51736972560423953923732838789657919150514997832d0
    real(8),parameter ::   w3 = 0.22238103445337447054435599442624088443013087005124956472590928&
         &929361681457044904085365314237719792784215926610&
         &121221812311143757985257224193818266745320905779&
         &086132895368404027893986488760043856972021574820&
         &632532471955902286315706513199655897335454406059&
         &52819880671616779621183704306688233150514997832d0
    real(8),parameter ::   w4 = 0.10122853629037625915253135430996219011539409105168495705900369&
         &806474017876347078486028273930404500655815438933&
         &141326670771549403089234876787319730411360735846&
         &905332088240507319763065757292054679614357794675&
         &524923287300550259929540899466768105108107294683&
         &66466585774650346143712142008566866150514997832d0
    f = h * ( w1*( self%fx(x-x1*h,m) + self%fx(x+x1*h,m)) + &
         w2*( self%fx(x-x2*h,m) + self%fx(x+x2*h,m)) + &
         w3*( self%fx(x-x3*h,m) + self%fx(x+x3*h,m)) + &
         w4*( self%fx(x-x4*h,m) + self%fx(x+x4*h,m)) )
  end function g8




  !  10-point method.
  recursive function g10(self, x, h, m) result(f)
    class(integration_type),intent(inout) :: self
    real(8), intent(in)                   :: x
    real(8), intent(in)                   :: h
    integer                               :: m
    real(8)                               :: f(m)
    !> abscissae:
    real(8),dimension(5),parameter ::  a = [   0.14887433898163121088482600112971998461756485942&
         &069169570798925351590361735566852137117762979946&
         &369123003116080525533882610289018186437654023167&
         &619699680909130507378277203710590709424758594227&
         &432498371771742473462169148529029429290031934666&
         &590824338380943550759968335702300050038372806343&
         &51d0,&
         0.43339539412924719079926594316578416220007183765&
         &624649650270151314376698907770350122510275795011&
         &772122368293504099893794727422475772324920512677&
         &410328220862009523192709334620320113283203876915&
         &840634111498011298231414887874432043247664144215&
         &76788807708483879452488118549797039287926964254222d0,&
         0.67940956829902440623432736511487357576929471183&
         &480946766481718895255857539507492461507857357048&
         &037949983390204739931506083674084257663009076827&
         &417182029235431978528469774097183691437120135529&
         &628377331531086791269325449548547293413247272116&
         &80274268486617121011712030227181051010718804444161d0,&
         0.86506336668898451073209668842349304852754301496&
         &533045252195973184537475513805556135679072894604&
         &577069440463108641176516867830016149345356373927&
         &293968909500115713496898930516120724357604809009&
         &797259233179237955357392905958797769568324277022&
         &36942765911483643714816923781701572597289139322313d0,&
         0.97390652851717172007796401208445205342826994669&
         &238211923121206669659520323463615962572356495626&
         &855625823304251877421121502216860143447777992054&
         &095872599424367044136957648812587991466331435107&
         &587371198778752105670674524353687136830338609093&
         &88311646653581707125686970668737259229449284383797d0 ]
    !> weights:
    real(8),dimension(5),parameter ::  w = [   0.29552422471475287017389299465133832942104671702&
         &685360135430802975599593821715232927035659579375&
         &421672271716440125255838681849078955200582600193&
         &634249418696660956271864888416804323130506153586&
         &740908305127066386528748390174687472659751595445&
         &0775158914556548308329986393605934912382356670244d0,&
         0.26926671930999635509122692156946935285975993846&
         &088379580056327624215343231917927676422663670925&
         &276075559581145036869830869292346938114524155646&
         &588466344237116560144322599601417290445280303444&
         &112979029770671425375348062846083992765750069116&
         &86749842814086288868533208042150419508881916391898d0,&
         0.21908636251598204399553493422816319245877187052&
         &267708988095654363519991065295128124268399317720&
         &219278659121687281288763476662690806694756883092&
         &118433166566771052699153220775367726528266710278&
         &782468510102088321733200642734832547562506684158&
         &85349420711613410227291565477768928313300688702802d0,&
         0.14945134915058059314577633965769733240255663966&
         &942736783547726875323865472663001094594726463473&
         &195191400575256104543633823445170674549760147137&
         &160119371095287981348288651187709535664396393337&
         &739399092016902046490838156187791575225783003434&
         &27785361756927642128792412282970150172590842897331d0,&
         0.06667134430868813759356880989333179285786483432&
         &015814512869488161341206408408710177678550968505&
         &887782109005471452041933148750712625440376213930&
         &498731699404163449536370640018701124231550439352&
         &624245062983271819871864748056604411786208647844&
         &9236378557180717569208295026105115288152794421677d0 ]
    f = h * ( w(1)*(  self%fx(x-a(1)*h,m)   +  self%fx(x+a(1)*h,m) ) + &
         w(2)*(  self%fx(x-a(2)*h,m)   +  self%fx(x+a(2)*h,m) ) + &
         w(3)*(  self%fx(x-a(3)*h,m)   +  self%fx(x+a(3)*h,m) ) + &
         w(4)*(  self%fx(x-a(4)*h,m)   +  self%fx(x+a(4)*h,m) ) + &
         w(5)*(  self%fx(x-a(5)*h,m)   +  self%fx(x+a(5)*h,m) )   )
  end function g10


  !  12-point method.
  recursive function g12(self, x, h, m) result(f)
    class(integration_type),intent(inout) :: self
    real(8), intent(in)                   :: x
    real(8), intent(in)                   :: h
    integer                               :: m
    real(8)                               :: f(m)
    !> abscissae:
    real(8),dimension(6),parameter ::  a = [   0.12523340851146891547244136946385312998339691630&
         &544427321292175474846205624138968874286829846949&
         &135959410459879132051097315159969664463407959720&
         &578930281363427149751877364610797786290401085851&
         &749803458163536009061915338533985792224380950454&
         &5097342064247739686883799517760948964137522919201d0,&
         0.36783149899818019375269153664371756125636014133&
         &540962131179987950408992951678787387873442850054&
         &657723463312639597714521513515217932743935324199&
         &163774275382871320389664162274303718284470963188&
         &934547884841822611461227526979609371629600504639&
         &62319787423676668046033025242558536362617894366679d0,&
         0.58731795428661744729670241894053428036909851404&
         &805248151027087966734069937589526243571076498874&
         &820190960155999292889267723106959108867175142499&
         &189843704151965799654931521792486834699342245746&
         &542270556959107871794349154143635139191674285545&
         &96877940491139756923177447689738849120865435563147d0,&
         0.76990267419430468703689383321281807598492575001&
         &893163766441906424911654310847122401642499922342&
         &191061761754045422185620704016285265354759491942&
         &035158754711514435184626896570143367857869960707&
         &068262822102488760216156789235759062543109515384&
         &10899341797549230707021382467596975621464477134163d0,&
         0.90411725637047485667846586611909619253759670921&
         &329754655407576068123479572923579048696942782373&
         &326781186038289641042234889971981954299601063524&
         &901258268291998347354448614206140899100247009682&
         &576258221693446448698746167580757842398074380920&
         &64065954540171679180850205196702894963912359448494d0,&
         0.98156063424671925069054909014928082296015519981&
         &373151046268212180779324431825398222525726789045&
         &223578555649237284127318524545703044707716708276&
         &967488752886112565550184482662910041202137201539&
         &996961235882788466302337187351583920530374414763&
         &9383170419389543470920618543180673569225988370568d0]
    !> weights:
    real(8),dimension(6),parameter ::  w = [   0.24914704581340278500056243604295121083046090256&
         &961883139535100311627942745728804303115680061804&
         &235306483347611787718583305851107360364968803964&
         &210377008542941502737221091728257019684301646591&
         &924021619820796255207324340857766137885796625403&
         &29347837170742904111565650371846972323325015720931d0,&
         0.23349253653835480876084989892487805625940997219&
         &975487473052349782149200007941167528067902650856&
         &369046673875643970886883389854278840891609661975&
         &038847380753533248145179488750388812162792803042&
         &489598308782293577290791644231030018795306547073&
         &15375809270840669989018891281956753131165193423269d0,&
         0.20316742672306592174906445580979837650651814727&
         &459014639859456579764563251047284379514439506460&
         &523243116042933686325996496137135190210132907910&
         &420189599423685656890245260738280276852445703846&
         &681240064758134063899875305215261728059344541572&
         &2327927963339557545261423500783899286052850767594d0,&
         0.16007832854334622633465252954335907187201173049&
         &086417790989954415795422517329115068165655263705&
         &773052707487709681280262724376386088264904467503&
         &100243409511213679869020659979278560098046375913&
         &998387244938872586336059061667742863824552984448&
         &70458396283884610940466728874776625823124924247387d0,&
         0.10693932599531843096025471819399622421457017347&
         &032488000512604210281899362749757654053731809631&
         &645741357635933314114416117033051696355084484800&
         &865232691960050114390447649204829355515358576079&
         &107052492180710337954704248957128309309678064675&
         &98358517298903536137451828089012822811396588037254d0,&
         0.04717533638651182719461596148501706031702907399&
         &484708956050534700380972115203871067082590707541&
         &453609661610167559673857966748040823913299673846&
         &365109909808575797967885849598965975687054894525&
         &799700269519193179311245399071070942125321236826&
         &63180160342232703368882666374567833050364187887189d0]
    f = h * ( w(1)*(  self%fx(x-a(1)*h,m)   +  self%fx(x+a(1)*h,m) ) + &
         w(2)*(  self%fx(x-a(2)*h,m)   +  self%fx(x+a(2)*h,m) ) + &
         w(3)*(  self%fx(x-a(3)*h,m)   +  self%fx(x+a(3)*h,m) ) + &
         w(4)*(  self%fx(x-a(4)*h,m)   +  self%fx(x+a(4)*h,m) ) + &
         w(5)*(  self%fx(x-a(5)*h,m)   +  self%fx(x+a(5)*h,m) ) + &
         w(6)*(  self%fx(x-a(6)*h,m)   +  self%fx(x+a(6)*h,m) ) )
  end function g12



  !>
  !  14-point method.
  recursive function g14(self, x, h, m) result(f)
    class(integration_type),intent(inout) :: self
    real(8), intent(in)                   :: x
    real(8), intent(in)                   :: h
    integer                               :: m
    real(8)                               :: f(m)
    !> abscissae:
    real(8),dimension(7),parameter ::  a = [  0.10805494870734366206624465021983474761195160547&
         &4237557040821061308013529011730007130100688176689&
         &3672374502026424466474638099232632258191427567218&
         &1973150409752806137273842265069487944308775321508&
         &8445556391329819060204836416480024319739665907101&
         &2506161702814425014635643221773541001328892761d0,&
         &0.31911236892788976043567182416847546683426120353&
         &3843956596650187257333440512792783164933705421346&
         &4131802793151826090394496145640578710017716508863&
         &2222396245608012120993128542172348808287716458637&
         &8479374239121304478425121768114783511643536777896&
         &2949997448460558214759676525644841351801594858d0,&
         &0.51524863635815409196529071855118866230888528256&
         &9306036951504769092784951832055660452072020350772&
         &8923922907932905090138695274035571340047593918260&
         &5653057211011637652073200342580823038532041784020&
         &3436173906624491224801618641571038235567674745455&
         &3979637438627635490786064892912451481973721288d0,&
         &0.68729290481168547014801980301933413753840121274&
         &7170675619266488628184896183133256947373070505211&
         &8384106603630216790054729627432715418501010682124&
         &6881727389082952662885443589912839338608106959371&
         &4595904926885388784713769175169784875289055161406&
         &7877996475717650653147982694804026342351254071d0,&
         &0.82720131506976499318979474265039496103970110147&
         &5081181560709054241479830810028873570426390137889&
         &5453991241406273986535333275661226737816179582645&
         &1069907936808669317564778014567859855078251147291&
         &5830426696849656086721489336979443959282673643228&
         &6425172143208924251106624044295037127737490111d0,&
         &0.92843488366357351733639113937787426447703921040&
         &9837618717962447482131093544359853111413905683657&
         &5176363551261559882603607008578010786539258018984&
         &5400440650494157888098179531161147719130825235345&
         &8596605653673043686690855550898698329741248613224&
         &5749388483890945436457404705549484348178721002d0,&
         &0.98628380869681233884159726670405280167609140723&
         &9225881644070811777749554132491637910646239665151&
         &7527602612562941358578689852603067447974494119727&
         &0324710898207170072955675048180261687970555989447&
         &5396929426197069500447181272675429908986256542893&
         &3676463914802477677291745002965827767360741735d0 ]
    !> weights:
    real(8),dimension(7),parameter ::  w = [   0.21526385346315779019587644331626003527499755805&
         &4128800219776392543618787353994604001024441410819&
         &5782372566723324367709929481659764649301890356019&
         &0805098142804175780269156508228762641736544919294&
         &6281203662033345376460522564310634412912654698349&
         &487266562730897512393716549425155133887783267d0,&
         &0.20519846372129560396592406566121805571033906130&
         &9419451716897290283367144825249720339431839991890&
         &8957243692694424494287284534856133850644865918702&
         &3021403166714178733299347482783913811132568481282&
         &5439676020905052976535424973123755325146919285189&
         &8072394707049964721031773292256965337005468577d0,&
         &0.18553839747793781374171659012515703624892260293&
         &7331659020034925069098350263525444425552731146712&
         &2229825611215057289188990778964974252160895085525&
         &2415283643607286404060027232379697141385075345609&
         &3331227890449938852384485366393922617921879824760&
         &6150274514935557012909889503067356410067833406d0,&
         &0.15720316715819353456960193862384215660566803733&
         &7323374969317043874768176369608298513958093362418&
         &0762768531519990811885018854374920646576267489242&
         &9103726460198700102219564745910784232280561068611&
         &6907713218466935160138377442838502265889923868443&
         &9084685022864905124096570215866733146092008329d0,&
         &0.12151857068790318468941480907247662595666934569&
         &0074672291075392543159743892526492318819906270375&
         &0071489155506530592569942811574313408868548096421&
         &2571445460802891854106154207862005646754562932960&
         &2540610239636717985405900755004972904989241013019&
         &1072357341821083329663867464821867539341968434d0,&
         &0.08015808715976020980563327706285430958369778539&
         &4594765201399065489571474457287169863536190819137&
         &7559686225015908038847487953091382572604434376755&
         &1198447409477973877237005366105771785226539545491&
         &2313554662497115946457665357652160093748935412771&
         &0937535198838649279475628473516378736712929573d0,&
         &0.03511946033175186303183287613819178061970560927&
         &7127276581499890196416322837808270537676796998646&
         &4636614217324764405511345585478510619843098677334&
         &0884595716394793248808744456729064741484147706750&
         &3186014306010893702617623540676052379390445897465&
         &9810087587180865408885105556219147609526200925d0 ]
    f = h * ( w(1)*(  self%fx(x-a(1)*h,m)   +  self%fx(x+a(1)*h,m) ) + &
         w(2)*(  self%fx(x-a(2)*h,m)   +  self%fx(x+a(2)*h,m) ) + &
         w(3)*(  self%fx(x-a(3)*h,m)   +  self%fx(x+a(3)*h,m) ) + &
         w(4)*(  self%fx(x-a(4)*h,m)   +  self%fx(x+a(4)*h,m) ) + &
         w(5)*(  self%fx(x-a(5)*h,m)   +  self%fx(x+a(5)*h,m) ) + &
         w(6)*(  self%fx(x-a(6)*h,m)   +  self%fx(x+a(6)*h,m) ) + &
         w(7)*(  self%fx(x-a(7)*h,m)   +  self%fx(x+a(7)*h,m) ) )
  end function g14








  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






  subroutine init_finter_1d(self,Xin,Fin,N)
    type(finter1d_type) :: self
    real(8)           :: xin(:)
    real(8)           :: fin(size(xin))
    integer           :: N,Lin
    if(self%status)call delete_finter_1d(self)
    Lin=size(xin)
    allocate(self%x(Lin),self%f(Lin))
    self%X    = Xin
    self%F    = Fin
    self%Imax = Lin
    self%Imin = 1
    self%N    = N
    self%status=.true.
  end subroutine init_finter_1d

  subroutine init_finter_2d(self,Xin,Yin,Fin,N)
    type(finter2d_type) :: self
    real(8)       :: xin(:),yin(:)
    real(8)       :: fin(size(xin),size(yin))
    integer       :: N,Lx,Ly
    if(self%status)call delete_finter_2d(self)
    Lx=size(xin) ; Ly=size(yin)
    allocate(self%x(Lx),self%y(Ly),self%f(Lx,Ly))
    self%X    = Xin
    self%Y    = Yin
    self%F    = Fin
    self%Imin = 1
    self%Jmin = 1
    self%Imax = Lx
    self%Jmax = Ly
    self%N    = N
    self%status=.true.
  end subroutine init_finter_2d

  subroutine delete_finter_1d(self)
    type(finter1d_type) :: self
    if(allocated(self%X))deallocate(self%X)
    if(allocated(self%F))deallocate(self%F)
    self%imax=0
    self%imin=0
    self%N   =0
    self%status=.false.
  end subroutine delete_finter_1d

  subroutine delete_finter_2d(self)
    type(finter2d_type) :: self
    if(allocated(self%x))deallocate(self%x)
    if(allocated(self%y))deallocate(self%y)
    if(allocated(self%f))deallocate(self%f)
    self%imin=0
    self%jmin=0
    self%imax=0
    self%jmax=0
    self%N   =0
    self%status=.false.
  end subroutine delete_finter_2d


  function locate(xx,x)
    real(8), dimension(:), intent(in) :: xx
    real(8), intent(in)               :: x
    integer                           :: locate
    integer                           :: n,jl,jm,ju
    logical                           :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  end function locate


  subroutine polint(xa,ya,x,y,dy)
    real(8), dimension(:), intent(in) :: xa,ya
    real(8), intent(in)               :: x
    real(8), intent(out)              :: y,dy
    integer                           :: m,n,ns
    real(8), dimension(size(xa))      :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.d0))then
          print*,'polint: calculation failure'
          stop
       endif
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint


  subroutine polin2(x1a,x2a,ya,x1,x2,y,dy)
    real(8), dimension(:), intent(in)   :: x1a,x2a
    real(8), dimension(:,:), intent(in) :: ya
    real(8), intent(in)                 :: x1,x2
    real(8), intent(out)                :: y,dy
    integer                             :: j,m,ndum
    real(8), dimension(size(x1a))       :: ymtmp
    real(8), dimension(size(x2a))       :: yntmp
    m = size(x1a);if(m/=size(ya,1))stop "POLINT: wrong dimensions m"
    ndum=size(x2a);if(ndum/=size(ya,2))stop "POLINT: wrong dimensions ndum"
    do j=1,m
       yntmp=ya(j,:)
       call polint(x2a,yntmp,x2,ymtmp(j),dy)
    end do
    call polint(x1a,ymtmp,x1,y,dy)
  end subroutine polin2


  function iminloc(arr)
    real(8), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc

  function assert_eq(n1,n2,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq
    if (n1 == n2) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq'
    end if
  end function assert_eq




  function gauss_quad_linspace(start,stop,num,istart,iend,mesh) result(array)
    integer          :: num,i
    real(8)          :: start,stop,step,array(num)
    logical,optional :: istart,iend
    logical          :: startpoint_,endpoint_
    real(8),optional :: mesh
    if(num<0)stop "gauss_quad_linspace: N<0, abort."
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(startpoint_.AND.endpoint_)then
       if(num<2)stop "gauss_quad_linspace: N<2 with both start and end points"
       step = (stop-start)/real(num-1,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(startpoint_.AND.(.not.endpoint_))then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(.not.startpoint_.AND.endpoint_)then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    else
       step = (stop-start)/real(num+1,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    endif
    if(present(mesh))mesh=step
  end function gauss_quad_linspace


END MODULE GAUSS_QUADRATURE








