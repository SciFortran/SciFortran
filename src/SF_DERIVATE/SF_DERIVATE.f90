MODULE SF_DERIVATE
  implicit none
  private

  complex(8),parameter :: one= (1.d0,0.d0)
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi = (0.d0,1.d0)
  real(8),parameter    :: pi=3.14159265358979323846264338327950288419716939937510D0


  !CENTRAL FINITE DIFFERENCE PARAMETERS:
  real(8),dimension(-1:1),parameter :: d_ccoeff_n2=[-1d0/2d0, 0d0, 1d0/2d0]
  real(8),dimension(-2:2),parameter :: d_ccoeff_n4=[1d0/12d0, -2d0/3d0, 0d0, 2d0/3d0, -1d0/12d0]
  real(8),dimension(-3:3),parameter :: d_ccoeff_n6=[-1d0/60d0, 3d0/20d0, -3d0/4d0, 0d0, 3d0/4d0, -3d0/20d0, 1d0/60d0] 
  real(8),dimension(-4:4),parameter :: d_ccoeff_n8=[1d0/280d0, -4d0/105d0, 1d0/5d0, -4d0/5d0, 0d0, 4d0/5d0, -1d0/5d0, 4d0/105d0, -1d0/280d0]


  real(8),dimension(-1:1),parameter :: d2_ccoeff_n2=[1d0, -2d0, 1d0]
  real(8),dimension(-2:2),parameter :: d2_ccoeff_n4=[-1d0/12d0, 4d0/3d0, -5d0/2d0, 4d0/3d0, -1d0/12d0]
  real(8),dimension(-3:3),parameter :: d2_ccoeff_n6=[1d0/90d0, -3d0/20d0, 3d0/2d0, -49d0/18d0, 3d0/2d0, -3d0/20d0, 1d0/90d0]
  real(8),dimension(-4:4),parameter :: d2_ccoeff_n8=[-1d0/560d0, 8d0/315d0, -1d0/5d0, 8d0/5d0, -205d0/72d0, 8d0/5d0, -1d0/5d0, 8d0/315d0, -1d0/560d0]


  real(8),dimension(-2:2),parameter :: d3_ccoeff_n2=[-1d0/2d0, 1d0, 0d0, -1d0, 1d0/2d0]
  real(8),dimension(-3:3),parameter :: d3_ccoeff_n4=[1d0/8d0, -1d0, 13d0/8d0, 0d0, -13d0/8d0, 1d0, -1d0/8d0]
  real(8),dimension(-4:4),parameter :: d3_ccoeff_n6=[-7d0/240d0, 3d0/10d0, -169d0/120d0, 61d0/30d0, 0d0, -61d0/30d0, 169d0/120d0, -3d0/10d0, 7d0/240d0]


  real(8),dimension(-2:2),parameter :: d4_ccoeff_n2=[1d0, -4d0, 6d0, -4d0, 1d0]
  real(8),dimension(-3:3),parameter :: d4_ccoeff_n4=[-1d0/6d0, 2d0, -13d0/2d0, 28d0/3d0, -13d0/2d0, 2d0, -1d0/6d0]
  real(8),dimension(-4:4),parameter :: d4_ccoeff_n6=[7d0/240d0, -2d0/5d0, 169d0/60d0, -122d0/15d0, 91d0/8d0, -122d0/15d0, 169d0/60d0, -2d0/5d0, 7d0/240d0]



  !FORWARD FINITE DIFFERENCE PARAMETERS:
  !BACKWARD = +FORWARD for N=even derivative order
  !BACKWARD = -FORWARD for N=odd derivative order
  real(8),dimension(0:1),parameter :: d_fcoeff_n1=[-1d0,1d0]
  real(8),dimension(0:2),parameter :: d_fcoeff_n2=[-3d0/2d0, 2d0, -1d0/2d0]
  real(8),dimension(0:3),parameter :: d_fcoeff_n3=[-11d0/6d0, 3d0, -3d0/2d0, 1d0/3d0]
  real(8),dimension(0:4),parameter :: d_fcoeff_n4=[-25d0/12d0, 4d0, -3d0, 4d0/3d0, -1d0/4d0]
  real(8),dimension(0:5),parameter :: d_fcoeff_n5=[-137d0/60d0, 5d0, -5d0, 10d0/3d0, -5d0/4d0, 1d0/5d0]
  real(8),dimension(0:6),parameter :: d_fcoeff_n6=[-49d0/20d0, 6d0, -15d0/2d0, 20d0/3d0, -15d0/4d0, 6d0/5d0, -1d0/6d0]

  real(8),dimension(0:2),parameter :: d2_fcoeff_n1=[1d0, -2d0, 1d0]
  real(8),dimension(0:3),parameter :: d2_fcoeff_n2=[2d0, -5d0, 4d0, -1d0]
  real(8),dimension(0:4),parameter :: d2_fcoeff_n3=[35d0/12d0, -26d0/3d0, 19d0/2d0, -14d0/3d0, 11d0/12d0]
  real(8),dimension(0:5),parameter :: d2_fcoeff_n4=[15d0/4d0, -77d0/6d0, 107d0/6d0, -13d0, 61d0/12d0, -5d0/6d0]
  real(8),dimension(0:6),parameter :: d2_fcoeff_n5=[203d0/45d0, -87d0/5d0, 117d0/4d0, -254d0/9d0, 33d0/2d0, -27d0/5d0, 137d0/180d0]
  real(8),dimension(0:7),parameter :: d2_fcoeff_n6=[469d0/90d0, -223d0/10d0, 879d0/20d0, -949d0/18d0, 41d0, -201d0/10d0, 1019d0/180d0, -7d0/10d0]

  real(8),dimension(0:3),parameter :: d3_fcoeff_n1=[-1d0, 3d0, -3d0, 1d0]
  real(8),dimension(0:4),parameter :: d3_fcoeff_n2=[-5d0/2d0, 9d0, -12d0, 7d0, -3d0/2d0]
  real(8),dimension(0:5),parameter :: d3_fcoeff_n3=[-17d0/4d0, 71d0/4d0, -59d0/2d0, 49d0/2d0, -41d0/4d0, 7d0/4d0]
  real(8),dimension(0:6),parameter :: d3_fcoeff_n4=[-49d0/8d0, 29d0, -461d0/8d0, 62d0, -307d0/8d0, 13d0, -15d0/8d0]
  real(8),dimension(0:7),parameter :: d3_fcoeff_n5=[-967d0/120d0, 638d0/15d0, -3929d0/40d0, 389d0/3d0, -2545d0/24d0, 268d0/5d0, -1849d0/120d0, 29d0/15d0]
  real(8),dimension(0:8),parameter :: d3_fcoeff_n6=[-801d0/80d0, 349d0/6d0, -18353d0/120d0, 2391d0/10d0, -1457d0/6d0, 4891d0/30d0, -561d0/8d0, 527d0/30d0, -469d0/240d0]


  real(8),dimension(0:4),parameter :: d4_fcoeff_n1=[1d0, -4d0, 6d0, -4d0, 1d0]
  real(8),dimension(0:5),parameter :: d4_fcoeff_n2=[3d0, -14d0, 26d0, -24d0, 11d0, -2d0]
  real(8),dimension(0:6),parameter :: d4_fcoeff_n3=[35d0/6d0, -31d0, 137d0/2d0, -242d0/3d0, 107d0/2d0, -19d0, 17d0/6d0]
  real(8),dimension(0:7),parameter :: d4_fcoeff_n4=[28d0/3d0, -111d0/2d0, 142d0, -1219d0/6d0, 176d0, -185d0/2d0, 82d0/3d0, -7d0/2d0]
  real(8),dimension(0:8),parameter :: d4_fcoeff_n5=[1069d0/80d0, -1316d0/15d0, 15289d0/60d0, -2144d0/5d0, 10993d0/24d0, -4772d0/15d0, 2803d0/20d0, -536d0/15d0, 967d0/240d0]



  interface djacobian
     module procedure fdjac_nn_func, fdjac_nn_sub, fdjac_mn_func,fdjac_mn_sub
  end interface djacobian

  interface dgradient
     module procedure fdjac_1n_func, fdjac_1n_sub
  end interface dgradient

  interface f_djacobian
     module procedure f_jac_nn_func, f_jac_nn_sub,  f_jac_mn_func, f_jac_mn_sub
  end interface f_djacobian

  interface f_dgradient
     module procedure f_jac_1n_func, f_jac_1n_sub
  end interface f_dgradient


  interface cjacobian
     module procedure c_fdjac_nn_func, c_fdjac_nn_sub, c_fdjac_mn_func, c_fdjac_mn_sub
  end interface cjacobian

  interface cgradient
     module procedure c_fdjac_1n_func,c_fdjac_1n_sub
  end interface cgradient

  interface f_cjacobian
     module procedure c_f_jac_nn_func , c_f_jac_nn_sub , c_f_jac_mn_func , c_f_jac_mn_sub
  end interface f_cjacobian

  interface f_cgradient
     module procedure c_f_jac_1n_func , c_f_jac_1n_sub
  end interface f_cgradient


  public :: djacobian
  public :: dgradient
  public :: f_djacobian
  public :: f_dgradient
  !
  public :: cjacobian
  public :: cgradient
  public :: f_cjacobian
  public :: f_cgradient
  !
  public :: deriv
  public :: derivative
  public :: derivative2
  public :: derivative3
  public :: derivative4
  public :: derivativeN

  
contains

  function deriv(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1)= (f(2)-f(1))/dh
    do i=2,L-1
       df(i) = (f(i+1)-f(i-1))/(2.d0*dh)
    enddo
    df(L)= (f(L)-f(L-1))/dh
  end function deriv










  ! Purpose: estimates an N/M/1 by N jacobian matrix using forward differences.
  !   computes a forward-difference approximation
  !   to the N by N jacobian matrix associated with a specified
  !   problem of N functions in N variables. If the jacobian has
  !   a banded form, then function evaluations are saved by only
  !   approximating the nonzero terms.
  !    -Input FCN: the name of the user-supplied function which
  !    calculates the functions. The routines accept functions/routines.
  !    Check out their explicit form in the interfaces of the subroutines
  !    -Input, real(8) X(:) the point where the jacobian is evaluated.
  !    -Output, real(8) FJAC(N,N) the approximate jacobian matrix.
  !    -Input,optional integer ML, MU, specify the number of subdiagonals and
  !    superdiagonals within the band of the jacobian matrix. Default values N-1
  !    -Input,optional real(8) EPSFCN, is used in determining a suitable step
  !    length for the forward-difference approximation.  This approximation
  !    assumes that the relative errors in the functions are of the order of
  !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
  !    the relative errors in the functions are of the order of the machine
  !    precision.
  include 'derivate_fjacobian_d.f90'
  include 'derivate_fjacobian_c.f90'




  function derivative(f,dh,order)  result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    integer,intent(in),optional     :: order
    real(8),dimension(size(f))      :: df
    integer                         :: i,L,order_
    L=size(f)
    order_=4;if(present(order))order_=order
    if(L < order_ + 1) stop "derivative: L < order+1."
    select case(order_)
    case(1)
       df = derivF_n121(f,dh)
    case(2)
       df = derivF_n222(f,dh)
    case(4)
       df = derivF_n444(f,dh)
    case(6)
       df = derivF_n666(f,dh)
    case default
       stop "derivative: order is not 1,2,4 or 6"
    end select
  end function derivative

  function derivF_n121(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d_fcoeff_n1,f(1:2))/dh
    do i=2,L-1
       df(i) = dot_product(d_ccoeff_n2(-1:1),f(i-1:i+1))/dh
    enddo
    df(L) = dot_product(-d_fcoeff_n1,f(L:L-1:-1))/dh
  end function derivF_n121

  function derivF_n222(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d_fcoeff_n2,f(1:3))/dh
    do i=2,L-1
       df(i) = dot_product(d_ccoeff_n2(-1:1),f(i-1:i+1))/dh
    enddo
    df(L) = dot_product(-d_fcoeff_n2,f(L:L-2:-1))/dh
  end function derivF_n222

  function derivF_n444(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d_fcoeff_n4,f(1:5))/dh
    df(2) = dot_product(d_fcoeff_n4,f(2:6))/dh
    do i=3,L-2
       df(i) = dot_product(d_ccoeff_n4(-2:2),f(i-2:i+2))/dh
    enddo
    df(L-1) = dot_product(-d_fcoeff_n4,f(L-1:L-5:-1))/dh
    df(L) = dot_product(-d_fcoeff_n4,f(L:L-4:-1))/dh
  end function derivF_n444

  function derivF_n666(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d_fcoeff_n6(0:6),f(1:7))/dh
    df(2) = dot_product(d_fcoeff_n6(0:6),f(2:8))/dh
    df(3) = dot_product(d_fcoeff_n6(0:6),f(3:9))/dh
    do i=4,L-3
       df(i) = dot_product(d_ccoeff_n6(-3:3),f(i-3:i+3))/dh
    enddo
    df(L-2) = dot_product(-d_fcoeff_n6(0:6),f(L-2:L-8:-1))/dh
    df(L-1) = dot_product(-d_fcoeff_n6(0:6),f(L-1:L-7:-1))/dh
    df(L) = dot_product(-d_fcoeff_n6(0:6),f(L:L-6:-1))/dh
  end function derivF_n666






  function derivative2(f,dh,order)  result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    integer,intent(in),optional     :: order
    real(8),dimension(size(f))      :: df
    integer                         :: i,L,order_
    L=size(f)
    order_=4;if(present(order))order_=order
    if(L < order_ + 1) stop "derivative2: L < order+1."
    select case(order_)
    case(2)
       df = derivF2_n222(f,dh)
    case(4)
       df = derivF2_n444(f,dh)
    case(6)
       df = derivF2_n666(f,dh)
    case default
       stop "derivative2: order is not 2,4 or 6"
    end select
  end function derivative2

  function derivF2_n222(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d2_fcoeff_n2,f(1:4))/dh**2
    df(2) = dot_product(d2_fcoeff_n2,f(2:5))/dh**2
    do i=3,L-2
       df(i) = dot_product(d2_ccoeff_n2(-1:1),f(i-1:i+1))/dh**2
    enddo
    !foo = f(L-1:L-4:-1)
    df(L-1) = dot_product(d2_fcoeff_n2,f(L-1:L-4:-1))/dh**2
    !foo = f(L:L-3:-1)
    df(L) = dot_product(d2_fcoeff_n2,f(L:L-3:-1))/dh**2
  end function derivF2_n222

  function derivF2_n444(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d2_fcoeff_n4,f(1:6))/dh**2
    df(2) = dot_product(d2_fcoeff_n4,f(2:7))/dh**2
    do i=3,L-2
       df(i) = dot_product(d2_ccoeff_n4(-2:2),f(i-2:i+2))/dh**2
    enddo
    df(L-1) = dot_product(d2_fcoeff_n4,f(L-1:L-6:-1))/dh**2
    df(L) = dot_product(d2_fcoeff_n4,f(L:L-5:-1))/dh**2
  end function derivF2_n444

  function derivF2_n666(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d2_fcoeff_n6(0:7),f(1:8))/dh**2
    df(2) = dot_product(d2_fcoeff_n6(0:7),f(2:9))/dh**2
    df(3) = dot_product(d2_fcoeff_n6(0:7),f(3:10))/dh**2
    do i=4,L-3
       df(i) = dot_product(d2_ccoeff_n6(-3:3),f(i-3:i+3))/dh**2
    enddo
    df(L-2) = dot_product(d2_fcoeff_n6,f(L-2:L-9:-1))/dh**2
    df(L-1) = dot_product(d2_fcoeff_n6,f(L-1:L-8:-1))/dh**2
    df(L) = dot_product(d2_fcoeff_n6,f(L:L-7:-1))/dh**2
  end function derivF2_n666







  function derivative3(f,dh,order)  result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    integer,intent(in),optional     :: order
    real(8),dimension(size(f))      :: df
    integer                         :: i,L,order_
    L=size(f)
    order_=4;if(present(order))order_=order
    if(L < order_ + 1) stop "derivative3: L < order+1."
    select case(order_)
    case(2)
       df = derivF3_n222(f,dh)
    case(4)
       df = derivF3_n444(f,dh)
    case(6)
       df = derivF3_n666(f,dh)
    case default
       stop "derivative3: order is not 2,4 or 6"
    end select
  end function derivative3

  function derivF3_n222(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,j,L
    L=size(f)
    df=0.d0
    df(1) = dot_product(d3_fcoeff_n2,f(1:5))/dh**3
    df(2) = dot_product(d3_fcoeff_n2,f(2:6))/dh**3
    do i=3,L-2
       df(i) = dot_product(d3_ccoeff_n2(-2:2),f(i-2:i+2))/dh**3
    enddo
    df(L-1) = dot_product(-d3_fcoeff_n2,f(L-1:L-5:-1))/dh**3
    df(L) = dot_product(-d3_fcoeff_n2,f(L:L-4:-1))/dh**3
  end function derivF3_n222

  function derivF3_n444(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d3_fcoeff_n4,f(1:7))/dh**3
    df(2) = dot_product(d3_fcoeff_n4,f(2:8))/dh**3
    df(3) = dot_product(d3_fcoeff_n4,f(3:9))/dh**3
    do i=4,L-3
       df(i) = dot_product(d3_ccoeff_n4(-3:3),f(i-3:i+3))/dh**3
    enddo
    df(L-2) = dot_product(-d3_fcoeff_n4,f(L-2:L-8:-1))/dh**3
    df(L-1) = dot_product(-d3_fcoeff_n4,f(L-1:L-7:-1))/dh**3
    df(L) = dot_product(-d3_fcoeff_n4,f(L:L-6:-1))/dh**3
  end function derivF3_n444

  function derivF3_n666(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d3_fcoeff_n6,f(1:9))/dh**3
    df(2) = dot_product(d3_fcoeff_n6,f(2:10))/dh**3
    df(3) = dot_product(d3_fcoeff_n6,f(3:11))/dh**3
    df(4) = dot_product(d3_fcoeff_n6,f(4:12))/dh**3
    do i=5,L-4
       df(i) = dot_product(d3_ccoeff_n6(-4:4),f(i-4:i+4))/dh**3
    enddo
    df(L-3) = dot_product(-d3_fcoeff_n6,f(L-3:L-11:-1))/dh**3
    df(L-2) = dot_product(-d3_fcoeff_n6,f(L-2:L-10:-1))/dh**3
    df(L-1) = dot_product(-d3_fcoeff_n6,f(L-1:L-9:-1))/dh**3
    df(L) = dot_product(-d3_fcoeff_n6,f(L:L-8:-1))/dh**3
  end function derivF3_n666








  function derivative4(f,dh,order)  result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    integer,intent(in),optional     :: order
    real(8),dimension(size(f))      :: df
    integer                         :: i,L,order_
    L=size(f)
    order_=4;if(present(order))order_=order
    if(L < order_ + 1) stop "derivative4: L < order+1."
    select case(order_)
    case(2)
       df = derivF4_n222(f,dh)
    case(4)
       df = derivF4_n444(f,dh)
    case(6)
       df = derivF4_n666(f,dh)
    case default
       stop "derivative4: order is not 2,4 or 6"
    end select
  end function derivative4

  function derivF4_n222(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,j,L
    L=size(f)
    df(1) = dot_product(d4_fcoeff_n2,f(1:6))/dh**4
    df(2) = dot_product(d4_fcoeff_n2,f(2:7))/dh**4
    do i=3,L-2
       df(i) = dot_product(d4_ccoeff_n2(-2:2),f(i-2:i+2))/dh**4
    enddo
    df(L-1) = dot_product(d4_fcoeff_n2,f(L-1:L-6:-1))/dh**4
    df(L) = dot_product(d4_fcoeff_n2,f(L:L-5:-1))/dh**4
  end function derivF4_n222


  function derivF4_n444(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d4_fcoeff_n4,f(1:8))/dh**4
    df(2) = dot_product(d4_fcoeff_n4,f(2:9))/dh**4
    df(3) = dot_product(d4_fcoeff_n4,f(3:10))/dh**4
    do i=4,L-3
       df(i) = dot_product(d4_ccoeff_n4(-3:3),f(i-3:i+3))/dh**4
    enddo
    df(L-2) = dot_product(d4_fcoeff_n4,f(L-2:L-9:-1))/dh**4
    df(L-1) = dot_product(d4_fcoeff_n4,f(L-1:L-8:-1))/dh**4
    df(L) = dot_product(d4_fcoeff_n4,f(L:L-7:-1))/dh**4
  end function derivF4_n444

  function derivF4_n666(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1) = dot_product(d4_fcoeff_n5,f(1:9))/dh**4
    df(2) = dot_product(d4_fcoeff_n5,f(2:10))/dh**4
    df(3) = dot_product(d4_fcoeff_n5,f(3:11))/dh**4
    df(4) = dot_product(d4_fcoeff_n5,f(4:12))/dh**4
    do i=5,L-4
       df(i) = dot_product(d4_ccoeff_n6(-4:4),f(i-4:i+4))/dh**4
    enddo
    df(L-3) = dot_product(d4_fcoeff_n5,f(L-3:L-11:-1))/dh**4
    df(L-2) = dot_product(d4_fcoeff_n5,f(L-2:L-10:-1))/dh**4
    df(L-1) = dot_product(d4_fcoeff_n5,f(L-1:L-9:-1))/dh**4
    df(L) = dot_product(d4_fcoeff_n5,f(L:L-8:-1))/dh**4
  end function derivF4_n666




  function derivativeN(f,dh,n)  result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    integer,intent(in)              :: n
    real(8),dimension(size(f))      :: df,tmp,dtmp
    integer                         :: i,L
    L=size(f)
    if(L < n + 1) stop "derivative4: L < order+1."
    tmp=f
    do i=1,n
       dtmp = derivF_n666(tmp,dh)
       tmp=dtmp
    enddo
    df=dtmp
  end function derivativeN


END MODULE SF_DERIVATE
