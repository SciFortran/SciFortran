!PURPOSE: Trapezoidal rule for 2d data function integration
! + _std: integrate f(:) between a and b. a,b default are 0,1
function d_trapz2d_sample(func,dhx,dhy,xrange,yrange) result(int)
  real(8),dimension(:,:)        :: func           !integrand: function f(x) with x=[x,y], or 2d array of samples f(x_i,y_j)
  real(8),optional              :: dhx,dhy        !grid spacings along x and y (default: from xrange,yrange and grid size)
  real(8),dimension(2),optional :: xrange,yrange  !integration limits [lower,upper] along x and y (default [0,1] for samples)
  real(8)                       :: hx,hy
  real(8),dimension(2)          :: xrange_,yrange_
  integer                       :: Nx,Ny,i,j
  real(8)                       :: int
  Nx=size(func,1)
  Ny=size(func,2)
  xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
  yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
  hx=xrange_(2)-xrange_(1);hx=hx/Nx
  if(present(dhx))hx=dhx
  hy=yrange_(2)-yrange_(1);hy=hy/Ny
  if(present(dhy))hy=dhy
  !
  int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
  do i=2,Nx
     do j=2,Ny
        int=int+4d0*func(i,j)
     enddo
  enddo
  do j=2,Ny
     int=int+2d0*( func(1,j) + func(Nx,j) )
  enddo
  do i=2,Nx
     int=int+2d0*( func(i,1) + func(i,Ny) )
  enddo
  int=int*hx*hy/4d0
end function d_trapz2d_sample

function c_trapz2d_sample(func,dhx,dhy,xrange,yrange) result(int)
  complex(8),dimension(:,:)     :: func           !integrand: function f(x) with x=[x,y], or 2d array of samples f(x_i,y_j)
  real(8),optional              :: dhx,dhy        !grid spacings along x and y (default: from xrange,yrange and grid size)
  real(8),dimension(2),optional :: xrange,yrange  !integration limits [lower,upper] along x and y (default [0,1] for samples)
  real(8)                       :: hx,hy
  real(8),dimension(2)          :: xrange_,yrange_
  integer                       :: Nx,Ny,i,j
  complex(8)                    :: int
  Nx=size(func,1)
  Ny=size(func,2)
  xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
  yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
  hx=xrange_(2)-xrange_(1);hx=hx/Nx
  if(present(dhx))hx=dhx
  hy=yrange_(2)-yrange_(1);hy=hy/Ny
  if(present(dhy))hy=dhy
  !
  int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
  do i=2,Nx
     do j=2,Ny
        int=int+4d0*func(i,j)
     enddo
  enddo
  do j=2,Ny
     int=int+2d0*( func(1,j) + func(Nx,j) )
  enddo
  do i=2,Nx
     int=int+2d0*( func(i,1) + func(i,Ny) )
  enddo
  int=int*hx*hy/4d0
end function c_trapz2d_sample




!+-----------------------------------------------------------------+
!PURPOSE: Simpson rule for 2d data function integration
!+-----------------------------------------------------------------+
function d_simps2d_sample(func,dhx,dhy,xrange,yrange) result(int)
  real(8),dimension(:,:)        :: func           !integrand: function f(x) with x=[x,y], or 2d array of samples f(x_i,y_j)
  real(8),optional              :: dhx,dhy        !grid spacings along x and y (default: from xrange,yrange and grid size)
  real(8),dimension(2),optional :: xrange,yrange  !integration limits [lower,upper] along x and y (default [0,1] for samples)
  real(8)                       :: hx,hy
  real(8),dimension(2)          :: xrange_,yrange_
  integer                       :: Nx,Ny,i,j
  real(8)                       :: int
  Nx=size(func,1)
  Ny=size(func,2)
  xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
  yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
  hx=xrange_(2)-xrange_(1);hx=hx/Nx/2
  if(present(dhx))hx=dhx
  hy=yrange_(2)-yrange_(1);hy=hy/Ny/2
  if(present(dhy))hy=dhy
  !
  int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
  !
  do j=1,Ny
     int=int+4d0*( func(1,2*j-1) + func(Nx,2*j-1) )
  enddo
  do j=1,Ny-1
     int=int+2d0*( func(1,2*j) + func(Nx,2*j) )
  enddo
  !
  do i=1,Nx
     int=int+4d0*( func(2*i-1,1) + func(2*i-1,Ny) )
  enddo
  do i=1,Nx-1
     int=int+2d0*( func(2*i,1) + func(2*i,Ny) )
  enddo
  !
  do j=1,Ny
     do i=1,Nx
        int=int+16d0*func(2*i-1,2*j-1)
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx
        int=int+8d0*func(2*i-1,2*j)
     enddo
  enddo
  !
  do j=1,Ny
     do i=1,Nx-1
        int=int+8d0*func(2*i,2*j-1)
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx-1
        int=int+4d0*func(2*i,2*j)
     enddo
  enddo
  int=int*hx*hy/9d0
end function d_simps2d_sample

function c_simps2d_sample(func,dhx,dhy,xrange,yrange) result(int)
  complex(8),dimension(:,:)     :: func           !integrand: function f(x) with x=[x,y], or 2d array of samples f(x_i,y_j)
  real(8),optional              :: dhx,dhy        !grid spacings along x and y (default: from xrange,yrange and grid size)
  real(8),dimension(2),optional :: xrange,yrange  !integration limits [lower,upper] along x and y (default [0,1] for samples)
  real(8)                       :: hx,hy
  real(8),dimension(2)          :: xrange_,yrange_
  integer                       :: Nx,Ny,i,j
  complex(8)                    :: int
  Nx=size(func,1)
  Ny=size(func,2)
  xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
  yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
  hx=xrange_(2)-xrange_(1);hx=hx/Nx/2
  if(present(dhx))hx=dhx
  hy=yrange_(2)-yrange_(1);hy=hy/Ny/2
  if(present(dhy))hy=dhy
  !
  int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
  !
  do j=1,Ny
     int=int+4d0*( func(1,2*j-1) + func(Nx,2*j-1) )
  enddo
  do j=1,Ny-1
     int=int+2d0*( func(1,2*j) + func(Nx,2*j) )
  enddo
  !
  do i=1,Nx
     int=int+4d0*( func(2*i-1,1) + func(2*i-1,Ny) )
  enddo
  do i=1,Nx-1
     int=int+2d0*( func(2*i,1) + func(2*i,Ny) )
  enddo
  !
  do j=1,Ny
     do i=1,Nx
        int=int+16d0*func(2*i-1,2*j-1)
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx
        int=int+8d0*func(2*i-1,2*j)
     enddo
  enddo
  !
  do j=1,Ny
     do i=1,Nx-1
        int=int+8d0*func(2*i,2*j-1)
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx-1
        int=int+4d0*func(2*i,2*j)
     enddo
  enddo
  int=int*hx*hy/9d0
end function c_simps2d_sample
