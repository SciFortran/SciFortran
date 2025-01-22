function d_trapz2d_func(func,xrange,yrange,Nx,Ny) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8)              :: func
     end function func
  end interface
  real(8),dimension(2) :: xrange,yrange
  integer              :: Nx,Ny,i,j
  real(8)              :: xx(Nx),yy(Ny)
  real(8)              :: hx,hy
  real(8)              :: int
  hx=xrange(2)-xrange(1)
  hx=hx/Nx
  hy=yrange(2)-yrange(1)
  hy=hy/Ny
  int=&
       func([xrange(1),yrange(1)])+&
       func([xrange(1),yrange(2)])+&
       func([xrange(2),yrange(1)])+&
       func([xrange(2),yrange(2)])
  xx=sf_integrate_linspace(xrange(1),xrange(2),Nx)
  yy=sf_integrate_linspace(yrange(1),yrange(2),Ny)
  do i=2,Nx
     do j=2,Ny
        int=int+4d0*func([xx(i),yy(j)])
     enddo
  enddo
  do j=2,Ny
     int=int+2d0*( func([xrange(1),yy(j)]) + func([xrange(2),yy(j)]) )
  enddo
  do i=2,Nx
     int=int+2d0*( func([xx(i),yrange(1)]) + func([xx(i),yrange(2)]) )
  enddo
  int=int*hx*hy/4d0
end function d_trapz2d_func

function c_trapz2d_func(func,xrange,yrange,Nx,Ny) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       complex(8)           :: func
     end function func
  end interface
  real(8),dimension(2) :: xrange,yrange
  integer              :: Nx,Ny,i,j
  real(8)              :: xx(Nx),yy(Ny)
  real(8)              :: hx,hy
  complex(8)           :: int
  hx=xrange(2)-xrange(1)
  hx=hx/Nx
  hy=yrange(2)-yrange(1)
  hy=hy/Ny
  int=&
       func([xrange(1),yrange(1)])+&
       func([xrange(1),yrange(2)])+&
       func([xrange(2),yrange(1)])+&
       func([xrange(2),yrange(2)])
  xx=sf_integrate_linspace(xrange(1),xrange(2),Nx)
  yy=sf_integrate_linspace(yrange(1),yrange(2),Ny)
  do i=2,Nx
     do j=2,Ny
        int=int+4d0*func([xx(i),yy(j)])
     enddo
  enddo
  do j=2,Ny
     int=int+2d0*( func([xrange(1),yy(j)]) + func([xrange(2),yy(j)]) )
  enddo
  do i=2,Nx
     int=int+2d0*( func([xx(i),yrange(1)]) + func([xx(i),yrange(2)]) )
  enddo
  int=int*hx*hy/4d0
end function c_trapz2d_func

!RECURSIVE VERSIONS:
function d_trapz2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8)              :: func
     end function func
  end interface
  real(8),dimension(2)      :: xrange,yrange
  integer                   :: N,icount
  real(8)                   :: int,eps,int0
  integer,optional          :: N0
  integer                   :: N0_
  logical,optional          :: iterative
  logical                   :: iterative_
  real(8),optional          :: threshold
  real(8)                   :: threshold_
  iterative_=.false.;if(present(iterative))iterative_=iterative
  N0_=51;if(present(N0))N0_=N0
  threshold_=1d0;if(iterative_)threshold_=5.d-3
  if(present(threshold))threshold_=threshold
  N=N0_
  eps=1d0
  icount=1
  int=d_trapz2d_func(func,xrange,yrange,N,N)
  do while (eps>threshold_)
     icount=icount+1
     int0=int
     N=2*N-10
     int=d_trapz2d_func(func,xrange,yrange,N,N)
     eps=abs(int-int0)/abs(int)
  enddo
end function d_trapz2d_func_recursive

function c_trapz2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       complex(8)           :: func
     end function func
  end interface
  real(8),dimension(2)      :: xrange,yrange
  integer                   :: N,icount
  real(8)                   :: eps,int0
  integer,optional          :: N0
  integer                   :: N0_
  logical,optional          :: iterative
  logical                   :: iterative_
  real(8),optional          :: threshold
  real(8)                   :: threshold_
  complex(8)                :: int
  iterative_=.false.;if(present(iterative))iterative_=iterative
  N0_=51;if(present(N0))N0_=N0
  threshold_=1d0;if(iterative_)threshold_=5.d-3
  if(present(threshold))threshold_=threshold
  N=N0_
  eps=1d0
  icount=1
  int=c_trapz2d_func(func,xrange,yrange,N,N)
  do while (eps>threshold_)
     icount=icount+1
     int0=int
     N=2*N-10
     int=c_trapz2d_func(func,xrange,yrange,N,N)
     eps=abs(int-int0)/abs(int)
  enddo
end function c_trapz2d_func_recursive







function d_simps2d_func(func,xrange,yrange,Nx,Ny) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8)              :: func
     end function func
  end interface
  real(8),dimension(2) :: xrange,yrange
  integer :: Nx,Ny,i,j
  real(8) :: xx(2*Nx),yy(2*Ny)
  real(8) :: hx,hy
  real(8) :: int
  hx=xrange(2)-xrange(1)
  hx=hx/Nx/2
  hy=yrange(2)-yrange(1)
  hy=hy/Ny/2
  xx=sf_integrate_linspace(xrange(1),xrange(2),2*Nx)
  yy=sf_integrate_linspace(yrange(1),yrange(2),2*Ny)
  !
  int=&
       func([xrange(1),yrange(1)])+&
       func([xrange(1),yrange(2)])+&
       func([xrange(2),yrange(1)])+&
       func([xrange(2),yrange(2)])
  !
  do j=1,Ny
     int=int+4d0*(func([xrange(1),yy(2*j-1)])+func([xrange(2),yy(2*j-1)]))
  enddo
  do j=1,Ny-1
     int=int+2d0*(func([xrange(1),yy(2*j)])+func([xrange(2),yy(2*j)]))
  enddo
  !
  do i=1,Nx
     int=int+4d0*(func([xx(2*i-1),yrange(1)])+func([xx(2*i-1),yrange(2)]))
  enddo
  do i=1,Nx-1
     int=int+2d0*(func([xx(2*i),yrange(1)])+func([xx(2*i),yrange(2)]))
  enddo
  !
  do j=1,Ny
     do i=1,Nx
        int=int+16d0*func([xx(2*i-1),yy(2*j-1)])
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx
        int=int+8d0*func([xx(2*i-1),yy(2*j)])
     enddo
  enddo
  !
  do j=1,Ny
     do i=1,Nx-1
        int=int+8d0*func([xx(2*i),yy(2*j-1)])
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx-1
        int=int+4d0*func([xx(2*i),yy(2*j)])
     enddo
  enddo
  int=int*hx*hy/9d0
end function d_simps2d_func

function c_simps2d_func(func,xrange,yrange,Nx,Ny) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       complex(8)           :: func
     end function func
  end interface
  real(8),dimension(2) :: xrange,yrange
  integer              :: Nx,Ny,i,j
  real(8)              :: xx(2*Nx),yy(2*Ny)
  real(8)              :: hx,hy
  complex(8)           :: int
  hx=xrange(2)-xrange(1)
  hx=hx/Nx/2
  hy=yrange(2)-yrange(1)
  hy=hy/Ny/2
  xx=sf_integrate_linspace(xrange(1),xrange(2),2*Nx)
  yy=sf_integrate_linspace(yrange(1),yrange(2),2*Ny)
  !
  int=&
       func([xrange(1),yrange(1)])+&
       func([xrange(1),yrange(2)])+&
       func([xrange(2),yrange(1)])+&
       func([xrange(2),yrange(2)])
  !
  do j=1,Ny
     int=int+4d0*(func([xrange(1),yy(2*j-1)])+func([xrange(2),yy(2*j-1)]))
  enddo
  do j=1,Ny-1
     int=int+2d0*(func([xrange(1),yy(2*j)])+func([xrange(2),yy(2*j)]))
  enddo
  !
  do i=1,Nx
     int=int+4d0*(func([xx(2*i-1),yrange(1)])+func([xx(2*i-1),yrange(2)]))
  enddo
  do i=1,Nx-1
     int=int+2d0*(func([xx(2*i),yrange(1)])+func([xx(2*i),yrange(2)]))
  enddo
  !
  do j=1,Ny
     do i=1,Nx
        int=int+16d0*func([xx(2*i-1),yy(2*j-1)])
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx
        int=int+8d0*func([xx(2*i-1),yy(2*j)])
     enddo
  enddo
  !
  do j=1,Ny
     do i=1,Nx-1
        int=int+8d0*func([xx(2*i),yy(2*j-1)])
     enddo
  enddo
  do j=1,Ny-1
     do i=1,Nx-1
        int=int+4d0*func([xx(2*i),yy(2*j)])
     enddo
  enddo
  int=int*hx*hy/9d0
end function c_simps2d_func

!RECURSIVE VERSIONS:
function d_simps2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       real(8)              :: func
     end function func
  end interface
  real(8),dimension(2)      :: xrange,yrange
  integer                   :: N,icount
  real(8)                   :: eps,int0
  integer,optional          :: N0
  integer                   :: N0_
  logical,optional          :: iterative
  logical                   :: iterative_
  real(8),optional          :: threshold
  real(8)                   :: threshold_
  real(8)                   :: int
  iterative_=.false.;if(present(iterative))iterative_=iterative
  N0_=51;if(present(N0))N0_=N0
  threshold_=1d0;if(iterative_)threshold_=1.d-3
  if(present(threshold))threshold_=threshold
  N=N0_
  eps=1d0
  icount=1
  int=d_simps2d_func(func,xrange,yrange,N,N)
  do while (eps>threshold_)
     icount=icount+1
     int0=int
     N=2*N-10
     int=d_simps2d_func(func,xrange,yrange,N,N)
     eps=abs(int-int0)/abs(int)
  enddo
end function d_simps2d_func_recursive

function c_simps2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
  interface
     function func(x)
       real(8),dimension(:) :: x
       complex(8)           :: func
     end function func
  end interface
  real(8),dimension(2)      :: xrange,yrange
  integer                   :: N,icount
  real(8)                   :: eps,int0
  integer,optional          :: N0
  integer                   :: N0_
  logical,optional          :: iterative
  logical                   :: iterative_
  real(8),optional          :: threshold
  real(8)                   :: threshold_
  complex(8)                :: int
  iterative_=.false.;if(present(iterative))iterative_=iterative
  N0_=51;if(present(N0))N0_=N0
  threshold_=1d0;if(iterative_)threshold_=1.d-3
  if(present(threshold))threshold_=threshold
  N=N0_
  eps=1d0
  icount=1
  int=c_simps2d_func(func,xrange,yrange,N,N)
  do while (eps>threshold_)
     icount=icount+1
     int0=int
     N=2*N-10
     int=c_simps2d_func(func,xrange,yrange,N,N)
     eps=abs(int-int0)/abs(int)
  enddo
end function c_simps2d_func_recursive
