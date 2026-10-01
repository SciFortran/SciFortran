subroutine pdf_allocate_1d(self,N)
  type(pdf_kernel) :: self  ! pdf object, not allocated
  integer          :: N     ! number of grid points
  if(self%status)stop "PDF_ALLOCATE: PDF already allocated"
  self%N    = N
  self%Ndata= 0
  allocate(self%x(N))
  allocate(self%pdf(N))
  self%status=.true.
  self%pdf = 0d0
  self%sigma = 0d0
end subroutine pdf_allocate_1d

subroutine pdf_deallocate_1d(self)
  type(pdf_kernel) :: self  ! pdf object
  if(.not.self%status)stop "PDF_DEALLOCATE: PDF not allocated"
  self%xmin=0d0
  self%xmax=0d0
  self%dx  =0d0
  self%sigma = 0d0
  self%N=0
  self%Ndata= 0
  deallocate(self%x)
  deallocate(self%pdf)
  self%status=.false.
end subroutine pdf_deallocate_1d



subroutine pdf_save_1d(self,pfile)
  type(pdf_kernel) :: self   ! pdf object
  character(len=*) :: pfile  ! name of the output file
  integer          :: i,unit
  if(.not.self%status)stop "PDF_SAVE: PDF not allocated"
  open(free_unit(unit),file=trim(pfile))
  write(unit,*)self%N
  write(unit,*)self%xmin
  write(unit,*)self%xmax
  write(unit,*)self%dx
  write(unit,*)self%Ndata
  write(unit,*)self%x
  write(unit,*)self%pdf
  write(unit,*)self%sigma
  write(unit,*)self%status
  write(unit,*)self%variance
  write(unit,*)self%rescale
  close(unit)
end subroutine pdf_save_1d


subroutine pdf_read_1d(self,pfile)
  type(pdf_kernel) :: self   ! pdf object, not allocated; allocated and set on output
  character(len=*) :: pfile  ! name of the file written by pdf_save
  integer          :: i,N,unit
  if(.not.self%status)write(*,"(A)")"PDF_READ: PDF not allocated"
  open(free_unit(unit),file=trim(pfile))
  read(unit,*)N
  call pdf_allocate(self,N)
  read(unit,*)self%xmin
  read(unit,*)self%xmax
  read(unit,*)self%dx
  read(unit,*)self%Ndata   
  read(unit,*)self%x
  read(unit,*)self%pdf
  read(unit,*)self%sigma
  read(unit,*)self%status
  read(unit,*)self%variance
  read(unit,*)self%rescale
  close(unit)
end subroutine pdf_read_1d



subroutine pdf_set_range_1d(self,a,b)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: a     ! lower bound of the grid
  real(8)          :: b     ! upper bound of the grid
  if(.not.self%status)stop "PDF_SET_RANGE: PDF not allocated"
  self%xmin=a
  self%xmax=b
  self%x = linspace(a,b,self%N,mesh=self%dx)
end subroutine pdf_set_range_1d




subroutine pdf_push_sigma_1d(self,sigma)
  type(pdf_kernel) :: self   ! pdf object
  real(8)          :: sigma  ! width (standard deviation) of the Gaussian kernel
  if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
  self%sigma = sigma
  self%variance=.true.
end subroutine pdf_push_sigma_1d



subroutine pdf_get_sigma_1d(self,sigma)
  type(pdf_kernel) :: self    
  real(8)          :: sigma
  if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
  sigma = self%sigma
end subroutine pdf_get_sigma_1d




subroutine pdf_sigma_data_1d(self,data,h)
  type(pdf_kernel)     :: self  ! pdf object
  real(8),dimension(:) :: data  ! sample of data
  real(8)              :: h     ! width of the kernel, Silverman rule
  if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
  h = get_sd(data)
  h = (4d0/3d0/size(data))**(1d0/5d0)*h !Silverman's rule of thumb.
end subroutine pdf_sigma_data_1d

subroutine pdf_sigma_sdev_1d(self,sdev,N,h)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: sdev  ! standard deviation of the data
  integer          :: N     ! number of data
  real(8)          :: h     ! width of the kernel, Silverman rule
  if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
  h = (4d0/3d0/N)**(1d0/5d0)*sdev !Silverman's rule of thumb.
end subroutine pdf_sigma_sdev_1d





subroutine pdf_accumulate_s_1d(self,data,sigma)
  type(pdf_kernel) :: self   ! pdf object
  real(8)          :: sigma_
  real(8)          :: data   ! data point
  real(8),optional :: sigma  ! optional: kernel width, if not stored in self
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
  self%pdf   = self%pdf + gaussian_kernel_1d(self%x,data,sigma_)
  self%Ndata = self%Ndata + 1
end subroutine pdf_accumulate_s_1d

subroutine pdf_accumulate_v_1d(self,data,sigma)
  type(pdf_kernel) :: self     ! pdf object
  real(8)          :: sigma_
  real(8)          :: data(:)  ! data points
  real(8),optional :: sigma    ! optional: kernel width, if not stored in self
  integer          :: i
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
  do i=1,size(data)
     self%pdf   = self%pdf + gaussian_kernel_1d(self%x,data(i),sigma_)
     self%Ndata = self%Ndata + 1
  enddo
end subroutine pdf_accumulate_v_1d




elemental function gaussian_kernel_1d(x,mean,sigma) result(gaussian_kernel)
  real(8),intent(in) :: x,mean,sigma
  real(8)            :: gaussian_kernel
  real(8),parameter  :: pi2=2d0*acos(-1d0)
  gaussian_kernel=exp(-0.5d0*((x-mean)/sigma)**2)/sqrt(pi2)/sigma
end function gaussian_kernel_1d


subroutine pdf_normalize_1d(self)
  type(pdf_kernel) :: self  ! pdf object
  if(.not.self%status)stop "PDF_NORMALIZE: PDF not allocated"
  if(.not.self%rescale)self%pdf = self%pdf/self%Ndata
  self%rescale = .true.
  self%pdf = self%pdf/sum(self%pdf)/self%dx
end subroutine pdf_normalize_1d



subroutine pdf_print_pfile_1d(self,pfile,normalize)
  type(pdf_kernel) :: self       ! pdf object
  character(len=*) :: pfile      ! name of the output file (appended)
  logical,optional :: normalize  ! optional: normalize first, default .true.
  logical          :: normalize_
  integer          :: i,unit
  if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
  normalize_ = .true.; if(present(normalize))normalize_=normalize
  if(normalize_)call pdf_normalize_1d(self)
  open(free_unit(unit),file=trim(pfile),access='append')
  do i=1,self%N
     write(unit,*)self%x(i),self%pdf(i)
  enddo
  write(unit,*)
  close(unit)
end subroutine pdf_print_pfile_1d




function pdf_mean_1d(self) result(mean)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: mean  ! mean of the pdf
  mean = simps(self%x*self%pdf,self%xmin,self%xmax)
end function pdf_mean_1d

function pdf_var_1d(self) result(var)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: mean
  real(8)          :: var   ! variance of the pdf
  mean = pdf_mean_1d(self)
  var  = simps((self%x-mean)**2*self%pdf,self%xmin,self%xmax)
end function pdf_var_1d

function pdf_sdev_1d(self) result(sdev)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: var
  real(8)          :: sdev  ! standard deviation of the pdf
  var  = pdf_var_1d(self)
  sdev = sqrt(var)
end function pdf_sdev_1d

function pdf_moment_1d(self,n,mu) result(mom)
  type(pdf_kernel) :: self  ! pdf object
  integer          :: n     ! order of the moment
  real(8),optional :: mu    ! optional: center of the moment, default 0
  real(8)          :: mu_
  real(8)          :: mom   ! moment of order n of the pdf
  mu_ = 0d0;if(present(mu))mu_=mu
  mom = simps((self%x-mu_)**n*self%pdf,self%xmin,self%xmax)
end function pdf_moment_1d

function pdf_skew_1d(self) result(skew)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: mean,sdev
  real(8)          :: skew  ! skewness of the pdf
  mean = pdf_mean_1d(self)
  sdev = pdf_sdev_1d(self)
  skew = pdf_moment_1d(self,3)
  skew = skew - 3*mean*sdev**2 - mean**3
  skew = skew/sdev**3
end function pdf_skew_1d

function pdf_curt_1d(self) result(curt)
  type(pdf_kernel) :: self  ! pdf object
  real(8)          :: mean,var
  real(8)          :: curt  ! kurtosis of the pdf (not excess)
  mean = pdf_mean_1d(self)
  var  = pdf_var_1d(self)
  curt = pdf_moment_1d(self,4,mean)
  curt = curt/var**2
end function pdf_curt_1d

subroutine pdf_print_moments_pfile_1d(self,pfile)
  type(pdf_kernel) :: self   ! pdf object
  character(len=*) :: pfile  ! name of the output file
  integer          :: i,unit
  if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
  open(free_unit(unit),file=trim(pfile))
  write(unit,*)pdf_mean_1d(self),pdf_sdev_1d(self),pdf_skew_1d(self),pdf_curt_1d(self)
  close(unit)
end subroutine pdf_print_moments_pfile_1d
