!PURPOSE:
! Trapezoidal rule for 1d data function integration between a and b
! or with respect to a given dh
!
! + _ab: given a,b and f(:) integrate f(:)
! + _dh: given dh=(b-a)/L-1/2 integrate f(:)
! + _nonlin: integrate f(:) using given x(:)
function d_trapz_ab_sample(f,a,b) result(sum)
  real(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8) :: a,b   !integration limits a (lower), b (upper); default 0,1 for functions
  real(8) :: dh
  real(8) :: sum
  integer :: i,L
  L=size(f)
  dh=(b-a)/(L-1)/2d0
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh
  enddo
end function d_trapz_ab_sample
function c_trapz_ab_sample(f,a,b) result(sum)
  complex(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8)    :: a,b   !integration limits a (lower), b (upper); default 0,1 for functions
  real(8)    :: dh
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  dh=(b-a)/(L-1)/2.d0
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh
  enddo
end function c_trapz_ab_sample




function d_trapz_dh_sample(f,dh) result(sum)
  real(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8) :: dh    !uniform grid spacing
  real(8) :: sum
  integer :: i,L
  L=size(f)
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh/2d0
  enddo
end function d_trapz_dh_sample
function c_trapz_dh_sample(f,dh) result(sum)
  complex(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8)    :: dh    !uniform grid spacing
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  sum=0.d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh/2.d0
  enddo
end function c_trapz_dh_sample





function d_trapz_nonlin_sample(f,x) result(sum)
  real(8) :: f(:)        !integrand: function f(x), or array of its samples f(x_i)
  real(8) :: x(size(f))  !non-uniform grid points x_i
  real(8) :: a,b,dh
  real(8) :: sum
  integer :: i,L
  L=size(f)
  a=minval(x)
  b=maxval(x)
  sum=0.d0
  do i=1,L-1
     dh  = (x(i+1)-x(i))/2.d0
     sum = sum + (f(i+1)+f(i))*dh
  enddo
end function d_trapz_nonlin_sample
function c_trapz_nonlin_sample(f,x) result(sum)
  complex(8) :: f(:)        !integrand: function f(x), or array of its samples f(x_i)
  real(8)    :: x(size(f))  !non-uniform grid points x_i
  real(8)    :: a,b,dh
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  a=minval(x)
  b=maxval(x)
  sum=0.d0
  do i=1,L-1
     dh  = (x(i+1)-x(i))/2.d0
     sum = sum + (f(i+1)+f(i))*dh
  enddo
end function c_trapz_nonlin_sample









!+-----------------------------------------------------------------+
!PURPOSE: Simpson rule for data function integration between extrema
! [a,b] or with respect to a given dh
!
! + _ab: given a,b and f(:) integrate f(:)
! + _dh: fiven dh=(b-a)/L-1/2 integrate f(:)
! + _nonlin: integrate f(:) using given x(:)
!+-----------------------------------------------------------------+
function d_simpson_ab_sample(f,a,b) result(sum)
  real(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8) :: a,b   !integration limits a (lower), b (upper); default 0,1 for functions
  real(8) :: dh
  real(8) :: sum
  integer :: L
  L=size(f)
  dh=(b-a)/dble(L-1)
  sum = d_simpson_dh_sample(f,dh)
end function d_simpson_ab_sample
function c_simpson_ab_sample(f,a,b) result(sum)
  complex(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  real(8)    :: a,b   !integration limits a (lower), b (upper); default 0,1 for functions
  real(8)    :: dh
  complex(8) :: sum
  integer    :: L
  L=size(f)
  dh=(b-a)/dble(L-1)
  sum = c_simpson_dh_sample(f,dh)
end function c_simpson_ab_sample


function d_simpson_dh_sample(f,dh) result(sum)
  real(8) :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  integer :: n
  real(8) :: dh    !uniform grid spacing
  real(8) :: sum,sum1,sum2,int1,int2
  integer :: i,p,m,mm,mmm
  N=size(f)
  sum=0d0
  if(N==1)return
  sum1=0d0
  sum2=0d0
  int1=0d0
  int2=0d0
  if(mod(n-1,2)==0)then                !if n-1 is even:
     do i=2,N-1,2
        sum1 = sum1 + f(i)
     enddo
     do i=3,N-2,2
        sum2 = sum2 + f(i)
     enddo
     sum = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n))*dh/3d0
  else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
     if (N>=6) then
        do i=2,N-4,2
           sum1 = sum1 + f(i)
        enddo
        do i=3,N-5,2
           sum2 = sum2 + f(i)
        enddo
        int1 = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n-3))*dh/3d0
     endif
     int2 = (f(n-3)+3d0*f(n-2)+3d0*f(n-1)+f(n))*dh*3d0/8d0
     sum  = int1 + int2
  end if
end function d_simpson_dh_sample
function c_simpson_dh_sample(f,dh) result(sum)
  complex(8)           :: f(:)  !integrand: function f(x), or array of its samples f(x_i)
  integer              :: n
  real(8)              :: dh    !uniform grid spacing
  complex(8)           :: sum,sum1,sum2,int1,int2
  integer              :: i,p,m
  complex(8),parameter :: zero=cmplx(0d0,0d0,8)
  N=size(f)
  sum=zero  
  if(N==1)return
  sum1=zero
  sum2=zero
  int1=zero
  int2=zero
  if(mod(n-1,2)==0)then                !if n-1 is even:
     do i=2,N-1,2
        sum1 = sum1 + f(i)
     enddo
     do i=3,N-2,2
        sum2 = sum2 + f(i)
     enddo
     sum = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n))*dh/3d0
  else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
     if (N>=6) then
        do i=2,N-4,2
           sum1 = sum1 + f(i)
        enddo
        do i=3,N-5,2
           sum2 = sum2 + f(i)
        enddo
        int1 = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n-3))*dh/3d0
     endif
     int2 = (f(n-3)+3d0*f(n-2)+3d0*f(n-1)+f(n))*dh*3d0/8d0
     sum  = int1 + int2
  end if
end function c_simpson_dh_sample




function d_simpson_nonlin_sample(f,x) result(sum)
  real(8) :: f(:)        !integrand: function f(x), or array of its samples f(x_i)
  real(8) :: x(size(f))  !non-uniform grid points x_i
  real(8) :: dx(size(f)+1)
  real(8) :: sum,sum1,sum2,sum3
  real(8) :: a,b,dh
  integer :: i,n,m
  n=size(f)
  m=n-1
  dx=0d0
  forall(i=1:m)dx(i)=x(i+1)-x(i)
  sum1=0d0
  sum2=0d0
  sum3=0d0
  i=0
  do while(i<n)
     !Simpson's 3/8 rule
     if((dx(i)==dx(i+1)).AND.(dx(i)==dx(i+2)))then
        sum1=sum1+(3.d0*dx(i)*(f(i)+&
             3.d0*(f(i+1)+f(i+2))+f(i+3)))/8.d0
        i=i+3
        !Simpson's 1/3 rule
     elseif(dx(i)==dx(i+1))then
        sum2=sum2+(2.d0*dx(i)*(f(i)+&
             4.d0*f(i+1)+f(i+2)))/6.d0
        i=i+2
        !trapezoidal rule
     elseif(dx(i)/=dx(i+1)) then
        sum3=sum3+dx(i)*(f(i)+f(i+1))/2.d0
        i = i + 1
     endif
  enddo
  sum = sum1+sum2+sum3
end function d_simpson_nonlin_sample
function c_simpson_nonlin_sample(f,x) result(sum)
  complex(8) :: f(:)       !integrand: function f(x), or array of its samples f(x_i)
  complex(8) :: sum
  real(8)    :: x(size(f))  !non-uniform grid points x_i
  real(8)    :: rsum,isum
  rsum=d_simpson_nonlin_sample(dreal(f),x)
  isum=d_simpson_nonlin_sample(dimag(f),x)
  sum  = dcmplx(rsum,isum)
end function c_simpson_nonlin_sample








