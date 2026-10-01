
!LMDIF INTERFACE:
!solve M nonlinear equations in N unknowns with M>N
!so f(x)=0 can NOT be solved.
!This look for a solution x so that the norm
! transpose(f(x))*f(x) is minimized.
subroutine leastsq_lmdif_func(func,a,m,tol,info)
  interface
     function func(a,m)
       real(8),dimension(:) :: a
       integer              :: m
       real(8),dimension(m) :: func
     end function func
  end interface
  real(8),dimension(:) :: a     ! in: starting guess; out: solution, size n
  integer              :: m     ! number of functions, m >= n
  real(8),optional     :: tol   ! optional: tolerance, default 1d-15
  integer,optional     :: info  ! optional out: status of the MINPACK routine
  real(8)              :: tol_
  integer              :: info_
  integer              :: n
  real(8),dimension(m) :: fvec
  tol_ = 1.d-15;if(present(tol))tol_=tol
  n=size(a)
  call lmdif1(leastsq_lmdif1_func2sub,m,n,a,fvec,tol_,info_)
  if(present(info))info=info_
  include "leastsq_error.h90"
contains
  subroutine leastsq_lmdif1_func2sub(m,n,a,fvec,iflag)
    integer ::  m
    integer ::  n
    real(8) ::  a(n)
    real(8) ::  fvec(m)
    integer ::  iflag
    fvec = func(a,m)
    if(iflag<0)stop "LEASTSQ_lmdif1_func2sub ERROR: iflag < 0 "
  end subroutine leastsq_lmdif1_func2sub
end subroutine leastsq_lmdif_func

subroutine leastsq_lmdif_sub(func,a,m,tol,info)
  interface
     subroutine func(a,m,f)
       real(8),dimension(:) :: a
       integer              :: m
       real(8),dimension(m) :: f
     end subroutine func
  end interface
  !
  real(8),dimension(:) :: a     ! in: starting guess; out: solution, size n
  integer              :: m     ! number of functions, m >= n
  real(8),optional     :: tol   ! optional: tolerance, default 1d-15
  integer,optional     :: info  ! optional out: status of the MINPACK routine
  real(8)              :: tol_
  integer              :: info_
  integer              :: n
  real(8),dimension(m) :: fvec
  tol_ = 1.d-15;if(present(tol))tol_=tol
  n=size(a)
  call lmdif1(leastsq_lmdif1_sub2sub,m,n,a,fvec,tol_,info_)
  if(present(info))info=info_
  include "leastsq_error.h90"
contains
  subroutine leastsq_lmdif1_sub2sub(m,n,a,fvec,iflag)
    integer ::  m
    integer ::  n
    real(8) ::  a(n)
    real(8) ::  fvec(m)
    integer ::  iflag
    call func(a,m,fvec)
    if(iflag<0)stop "LEASTSQ_LMDIF1_sub2sub ERROR: iflag < 0 "
  end subroutine leastsq_lmdif1_sub2sub
end subroutine leastsq_lmdif_sub






!LMDER INTERFACE:
!solve M nonlinear equations in N unknowns with M>N
!so f(x)=0 can NOT be solved.
!This look for a solution x so that the norm
! transpose(f(x))*f(x) is minimized.
subroutine leastsq_lmder_func(func,dfunc,a,m,tol,info)
  interface
     function func(a,m)
       real(8),dimension(:) :: a
       integer              :: m
       real(8),dimension(m) :: func
     end function func
     !
     function dfunc(a,m)
       real(8),dimension(:)         :: a
       integer                      :: m
       real(8),dimension(m,size(a)) :: dfunc
     end function dfunc
  end interface
  real(8),dimension(:)         :: a     ! in: starting guess; out: solution, size n
  integer                      :: m     ! number of functions, m >= n
  real(8),optional             :: tol   ! optional: tolerance, default 1d-15
  integer,optional             :: info  ! optional out: status of the MINPACK routine
  real(8)                      :: tol_
  integer                      :: info_
  integer                      :: n
  real(8),dimension(m)         :: fvec
  real(8),dimension(m,size(a)) :: fjac
  tol_ = 1.d-15;if(present(tol))tol_=tol
  n=size(a)
  call lmder1(leastsq_lmder1_func2sub,m,n,a,fvec,fjac,m,tol_,info_)
  if(present(info))info=info_
  include "leastsq_error.h90"

contains
  subroutine leastsq_lmder1_func2sub(m,n,a,fvec,fjac,ldfjac,iflag)
    integer ::  m
    integer ::  n
    integer ::  ldfjac
    real(8) ::  a(n)
    real(8) ::  fvec(m)
    real(8) ::  fjac(ldfjac,n)
    integer ::  iflag
    if(iflag==1)then
       fvec = func(a,m)
    elseif(iflag==2)then
       fjac = dfunc(a,m)
    endif
    if(iflag<0)stop "LEASTSQ_LMDER1_func2sub ERROR: iflag < 0 "
  end subroutine leastsq_lmder1_func2sub
end subroutine leastsq_lmder_func

subroutine leastsq_lmder_sub(func,dfunc,a,m,tol,info)
  interface
     subroutine func(a,m,f)
       real(8),dimension(:) :: a
       integer              :: m
       real(8),dimension(m) :: f
     end subroutine func
     !
     subroutine dfunc(a,m,df)
       real(8),dimension(:)         :: a
       integer                      :: m
       real(8),dimension(m,size(a)) :: df
     end subroutine dfunc
  end interface
  real(8),dimension(:)         :: a     ! in: starting guess; out: solution, size n
  integer                      :: m     ! number of functions, m >= n
  real(8),optional             :: tol   ! optional: tolerance, default 1d-15
  integer,optional             :: info  ! optional out: status of the MINPACK routine
  real(8)                      :: tol_
  integer                      :: info_
  integer                      :: n
  real(8),dimension(m)         :: fvec
  real(8),dimension(m,size(a)) :: fjac
  tol_ = 1.d-15;if(present(tol))tol_=tol
  n=size(a)
  call lmder1(leastsq_lmder1_sub2sub,m,n,a,fvec,fjac,m,tol_,info_)
  if(present(info))info=info_
  include "leastsq_error.h90"
contains
  subroutine leastsq_lmder1_sub2sub(m,n,a,fvec,fjac,ldfjac,iflag)
    integer ::  m
    integer ::  n
    integer ::  ldfjac
    real(8) ::  a(n)
    real(8) ::  fvec(m)
    real(8) ::  fjac(ldfjac,n)
    integer ::  iflag
    if(iflag==1)then
       call func(a,m,fvec)
    elseif(iflag==2)then
       call dfunc(a,m,fjac)
    endif
    if(iflag<0)stop "LEASTSQ_LMDER1_sub2sub ERROR: iflag < 0 "
  end subroutine leastsq_lmder1_sub2sub
end subroutine leastsq_lmder_sub





