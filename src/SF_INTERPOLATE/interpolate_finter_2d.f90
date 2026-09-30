subroutine init_finter2d(func,Xin,Yin,Fin,N)
  !This subroutine initializes a :f:type:`finter2d_type` object with the tabulated data of a
  !real function of two variables on the regular grid :f:var:`Xin` :math:`\times` :f:var:`Yin`, 
  !so that it can be interpolated with :f:func_inline:`finter2d`. The data are copied, 
  !:f:var:`func` does not keep a reference to :f:var:`Xin`, :f:var:`Yin` and :f:var:`Fin`. If 
  !:f:var:`func` was already initialized its memory is first released.
  type(finter2d_type) :: func                !finter2d_type object to initialize (previous content is released)
  real(8)       :: xin(:),yin(:)             !grid points along x and y, strictly increasing
  real(8)       :: fin(size(xin),size(yin))  !function values on the grid, Fin(i,j)=F(Xin(i),Yin(j))
  integer       :: N                         !order of the local polynomial interpolation (window of (N+2)x(N+2) grid points)
  integer       :: Lx,Ly
  if(func%status)deallocate(func%x,func%y,func%f)
  Lx=size(xin) ; Ly=size(yin)
  allocate(func%x(Lx),func%y(Ly),func%f(Lx,Ly))
  func%X    = Xin
  func%Y    = Yin
  func%F    = Fin
  func%Imin = 1
  func%Jmin = 1
  func%Imax = Lx
  func%Jmax = Ly
  func%N    = N
  func%status=.true.
end subroutine init_finter2d





subroutine delete_finter2d(func)
  !This subroutine releases the memory of the :f:type:`finter2d_type` object :f:var:`func` and
  !resets it to the uninitialized state: :code:`Imin=Imax=Jmin=Jmax=N=0` and :code:`status=.false.`.
  type(finter2d_type) :: func  !finter2d_type object to release
  if(allocated(func%x))deallocate(func%x)
  if(allocated(func%y))deallocate(func%y)
  if(allocated(func%f))deallocate(func%f)
  func%imin=0
  func%jmin=0
  func%imax=0
  func%jmax=0
  func%N   =0
  func%status=.false.
end subroutine delete_finter2d








function finter2d(func,x,y)
  !This function evaluates at the point (:f:var:`x`, :f:var:`y`) the function stored in the
  !:f:type:`finter2d_type` object :f:var:`func`. It finds the grid cell containing the point 
  !with :f:func_inline:`locate` along each direction and interpolates with :f:func_inline:`polin2` 
  !on a window of :math:`(N+2) \times (N+2)` grid points around it, with :code:`N=func%N`. 
  !The window is shifted to stay inside the grid, so points outside the grid are extrapolated 
  !with the polynomial of the first or last window.
  real(8)         :: x,y       !coordinates x,y of the point where the function is interpolated
  type(finter2d_type) :: func  !finter2d_type object, initialized with init_finter2d
  real(8)         :: finter2d
  real(8)         :: f,df
  integer         :: itmp,jtmp,kx,ky,k0x,k0y,k1x,k1y
  integer         :: n
  N=func%N    !order of polynomial interpolation
  finter2d=0.d0
  itmp=locate(func%X(func%Imin:func%Imax),x)
  jtmp=locate(func%Y(func%Jmin:func%Jmax),y)
  kx=max(itmp-(N-1)/2,1)
  ky=max(jtmp-(N-1)/2,1)
  k0x = kx ; if(k0x < func%Imin)k0x=func%Imin
  k0y = ky ; if(k0y < func%Jmin)k0y=func%Jmin
  k1x = kx+N+1
  if(k1x > func%Imax)then         
     k1x=func%Imax
     k0x=k1x-N-1
  endif
  k1y = ky+N+1
  if(k1y > func%Jmax)then
     k1y=func%Jmax
     k0y=k1y-N-1
  endif
  call polin2(func%X(k0x:k1x),func%Y(k0y:k1y),func%F(k0x:k1x,k0y:k1y),x,y,f,df)
  finter2d=f
  return
end function finter2d

