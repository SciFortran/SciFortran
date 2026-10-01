subroutine pdf_allocate_2d(self,Nvec)
  type(pdf_kernel_2d) :: self     ! pdf object, not allocated
  integer             :: Nvec(2)  ! number of grid points along x and y
  if(self%status)stop "PDF_ALLOCATE: PDF already allocated"
  self%N     = Nvec
  self%Ndata = 0
  allocate(self%x(self%N(1)))
  allocate(self%y(self%N(2)))
  allocate(self%pdf(self%N(1),self%N(2)))
  self%status=.true.
  self%pdf   = 0d0
  self%sigma = 0d0
end subroutine pdf_allocate_2d

subroutine pdf_deallocate_2d(self)
  type(pdf_kernel_2d) :: self  ! pdf object
  if(.not.self%status)stop "PDF_DEALLOCATE: PDF not allocated"
  self%xmin = 0d0
  self%xmax = 0d0
  self%dx   = 0d0
  self%Sigma= 0d0
  self%N    = 0
  self%Ndata= 0
  deallocate(self%x)
  deallocate(self%y)
  deallocate(self%pdf)
  self%status=.false.
end subroutine pdf_deallocate_2d


subroutine pdf_save_2d(self,pfile)
  type(pdf_kernel_2d) :: self  ! pdf object
  character(len=*) :: pfile    ! name of the output file
  integer          :: i,unit
  if(.not.self%status)stop "PDF_SAVE: PDF not allocated"
  open(free_unit(unit),file=trim(pfile))
  write(unit,*)self%N
  write(unit,*)self%xmin
  write(unit,*)self%xmax
  write(unit,*)self%dx
  write(unit,*)self%Ndata
  write(unit,*)self%x
  write(unit,*)self%y
  write(unit,*)self%pdf
  write(unit,*)self%sigma
  write(unit,*)self%status
  write(unit,*)self%variance
  write(unit,*)self%rescale
  close(unit)
end subroutine pdf_save_2d


subroutine pdf_read_2d(self,pfile)
  type(pdf_kernel_2d) :: self  ! pdf object, not allocated; allocated and set on output
  character(len=*) :: pfile    ! name of the file written by pdf_save
  integer          :: i,N1,N2,unit
  if(.not.self%status)then
     print*,"PDF_READ: PDF not allocated"
  endif
  open(free_unit(unit),file=trim(pfile))
  read(unit,*)N1,N2
  call pdf_allocate(self,[N1,N2])
  read(unit,*)self%xmin
  read(unit,*)self%xmax
  read(unit,*)self%dx
  read(unit,*)self%Ndata   
  read(unit,*)self%x
  read(unit,*)self%y
  read(unit,*)self%pdf
  read(unit,*)self%sigma
  read(unit,*)self%status
  read(unit,*)self%variance
  read(unit,*)self%rescale
  close(unit)
end subroutine pdf_read_2d



subroutine pdf_set_range_2d(self,a,b)
  type(pdf_kernel_2d) :: self  ! pdf object
  real(8)             :: a(2)  ! lower bounds of the grid, (x,y)
  real(8)             :: b(2)  ! upper bounds of the grid, (x,y)
  if(.not.self%status)stop "PDF_SET_RANGE: PDF not allocated"
  self%xmin=a
  self%xmax=b
  self%x = linspace(a(1),b(1),self%N(1),mesh=self%dx(1))
  self%y = linspace(a(2),b(2),self%N(2),mesh=self%dx(2))
end subroutine pdf_set_range_2d




subroutine pdf_push_sigma_2d(self,sigma)
  type(pdf_kernel_2d) :: self        ! pdf object
  real(8)             :: sigma(2,2)  ! covariance matrix of the Gaussian kernel
  if(.not.self%status)stop "PDF_PUSH_SIGMA: PDF not allocated"
  self%sigma = sigma
  self%variance=.true.
end subroutine pdf_push_sigma_2d




subroutine pdf_get_sigma_2d(self,sigma)
  type(pdf_kernel_2d) :: self        ! pdf object
  real(8)             :: sigma(2,2)  ! covariance matrix of the Gaussian kernel stored in self
  if(.not.self%status)stop "PDF_GET_SIGMA: PDF not allocated"
  sigma = self%sigma
end subroutine pdf_get_sigma_2d




subroutine pdf_sigma_data_2d(self,data,h)
  type(pdf_kernel_2d)    :: self    ! pdf object
  real(8),dimension(:,:) :: data    ! sample of data, [2,L]
  real(8)                :: h(2,2)  ! covariance matrix of the kernel, Silverman rule
  integer                :: L
  if(.not.self%status)stop "PDF_=SIGMA: PDF not allocated"
  L = size(data,2)
  if(any(shape(data)/=[2,L]))stop "PDF_SET_SIGMA: DATA wrong dimensions"
  h = 0d0
  h(1,1) = get_var(data(1,:))
  h(2,2) = get_var(data(2,:))
  h(1,1) = (1d0/size(data,2))**(2/6d0)*h(1,1) !Silverman's rule of thumb.
  h(2,2) = (1d0/size(data,2))**(2/6d0)*h(2,2) !Silverman's rule of thumb.
end subroutine pdf_sigma_data_2d

subroutine pdf_sigma_sdev_2d(self,sdev,Nvec,h)
  type(pdf_kernel_2d) :: self     ! pdf object
  real(8)             :: sdev(2)  ! standard deviations of the data, (x,y)
  integer             :: Nvec(2)  ! number of data, (x,y)
  real(8)             :: h(2,2)   ! covariance matrix of the kernel, Silverman rule
  if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
  h      = 0d0
  h(1,1) = (1d0/Nvec(1))**(2/6d0)*sdev(1)**2 !Silverman's rule of thumb.
  h(2,2) = (1d0/Nvec(2))**(2/6d0)*sdev(2)**2 !Silverman's rule of thumb.
end subroutine pdf_sigma_sdev_2d




subroutine pdf_accumulate_s_2d(self,data,sigma)
  type(pdf_kernel_2d) :: self        ! pdf object
  real(8)             :: data(2)     ! data point (x,y)
  real(8),optional    :: sigma(2,2)  ! optional: kernel covariance matrix, if not stored in self
  real(8)             :: sigma_(2,2)
  integer             :: i,j
  if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
  if(self%variance)then
     sigma_ = self%sigma
  else
     if(present(sigma))then
        sigma_ = sigma
     else
        stop "PDF_ACCUMULATE: PDF sigma not set or passed"
     endif
  endif
  !
  self%pdf   = self%pdf + gaussian_kernel_2d(self%x,self%y,data,sigma_)
  self%Ndata = self%Ndata + 1
end subroutine pdf_accumulate_s_2d



function gaussian_kernel_2d(x,y,mean,sigma) result(gaussian_kernel)
  !This function returns the bivariate Gaussian kernel :math:`e^{-\frac{1}{2}(r-m)^T\Sigma^{-1}(r-m)}/2\pi\sqrt{\det\Sigma}`
  !with mean :f:var:`mean` and covariance matrix :f:var:`sigma`, evaluated on the grid of the points :math:`r=(x_i,y_j)`.
  !
  real(8),intent(in) :: x(:)                              ! grid along x
  real(8),intent(in) :: y(:)                              ! grid along y
  real(8),intent(in) :: mean(2)                           ! center of the kernel, (x,y)
  real(8),intent(in) :: sigma(2,2)                        ! covariance matrix of the kernel, [2,2]
  real(8)            :: gaussian_kernel(size(x),size(y))  ! value of the kernel on the grid, [size(x),size(y)]
  real(8),parameter  :: pi2=2d0*acos(-1d0)
  real(8)            :: detSigma
  real(8)            :: InvSigma(2,2)
  real(8)            :: gaussian_factor
  real(8)            :: arg,xvec(2)
  integer            :: i,j
  !
  detSigma = det(sigma)
  gaussian_factor=1/pi2/sqrt(detSigma)
  !
  InvSigma = Sigma
  call inv(InvSigma)
  do i=1,size(x)
     do j=1,size(y)
        xvec = [x(i),y(j)] - mean
        arg  = dot_product(xvec,matmul(InvSigma,xvec))
        !
        gaussian_kernel(i,j)=gaussian_factor*exp(-0.5d0*arg)
        !
     enddo
  enddo
end function gaussian_kernel_2d




subroutine pdf_normalize_2d(self)
  type(pdf_kernel_2d) :: self  ! pdf object
  if(.not.self%status)stop "PDF_NORMALIZE: PDF not allocated"
  if(.not.self%rescale)self%pdf = self%pdf/self%Ndata
  self%rescale = .true.
  self%pdf = self%pdf/sum(self%pdf)/product(self%dx)
end subroutine pdf_normalize_2d



subroutine pdf_print_pfile_2d(self,pfile,normalize)
  type(pdf_kernel_2d) :: self       ! pdf object
  character(len=*)    :: pfile      ! name of the output file
  logical,optional    :: normalize  ! optional: normalize first, default .true.
  logical             :: normalize_
  integer             :: i,unit
  if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
  normalize_ = .true.; if(present(normalize))normalize_=normalize
  if(normalize_)call pdf_normalize(self)
  call splot3d(pfile,self%x,self%y,self%pdf)
end subroutine pdf_print_pfile_2d

