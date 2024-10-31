MODULE SF_STAT
!SciFortran module for statistical estimators and visualization
  USE SF_ARRAYS, only: linspace
  USE SF_INTEGRATE, only: simps
  USE SF_IOTOOLS, only: free_unit,splot3d
  USE SF_LINALG, only: det,inv
  implicit none
  private

  type,public :: histogram
     integer                       :: n=0
     real(8),dimension(:),pointer  :: range
     real(8),dimension(:),pointer  :: bin
  end type histogram

  type,public :: pdf_kernel
     integer                          :: N=0
     real(8)                          :: xmin
     real(8)                          :: xmax
     real(8)                          :: dx
     integer                          :: Ndata=0     
     real(8),dimension(:),allocatable :: x
     real(8),dimension(:),allocatable :: pdf
     real(8)                          :: sigma
     logical                          :: status=.false.
     logical                          :: variance=.false.
     logical                          :: rescale=.false.
  end type pdf_kernel

  type,public :: pdf_kernel_2d
     integer,dimension(2)               :: N
     real(8),dimension(2)               :: xmin
     real(8),dimension(2)               :: xmax
     real(8),dimension(2)               :: dx
     integer                            :: Ndata=0   
     real(8),dimension(:),allocatable   :: x,y
     real(8),dimension(:,:),allocatable :: pdf
     real(8),dimension(2,2)             :: Sigma
     logical                            :: status=.false.
     logical                            :: variance=.false.
     logical                            :: rescale=.false.
  end type pdf_kernel_2d


  interface pdf_allocate
     module procedure :: pdf_allocate_1d
     module procedure :: pdf_allocate_2d
  end interface pdf_allocate

  interface pdf_deallocate
     module procedure :: pdf_deallocate_1d
     module procedure :: pdf_deallocate_2d
  end interface pdf_deallocate


  interface pdf_save
     module procedure :: pdf_save_1d
     module procedure :: pdf_save_2d
  end interface pdf_save


  interface pdf_read
     module procedure :: pdf_read_1d
     module procedure :: pdf_read_2d
  end interface pdf_read


  interface pdf_set_range
     module procedure :: pdf_set_range_1d
     module procedure :: pdf_set_range_2d
  end interface pdf_set_range


  interface pdf_push_sigma
     module procedure :: pdf_push_sigma_1d
     module procedure :: pdf_push_sigma_2d
  end interface pdf_push_sigma


  interface pdf_get_sigma
     module procedure :: pdf_get_sigma_1d
     module procedure :: pdf_get_sigma_2d
  end interface pdf_get_sigma


  interface pdf_sigma
     module procedure :: pdf_sigma_data_1d
     module procedure :: pdf_sigma_sdev_1d
     module procedure :: pdf_sigma_data_2d
     module procedure :: pdf_sigma_sdev_2d
  end interface pdf_sigma


  interface pdf_accumulate
     module procedure :: pdf_accumulate_s_1d
     module procedure :: pdf_accumulate_v_1d
     module procedure :: pdf_accumulate_s_2d
  end interface pdf_accumulate



  interface pdf_normalize
     module procedure :: pdf_normalize_1d
     module procedure :: pdf_normalize_2d
  end interface pdf_normalize


  interface pdf_print
     module procedure :: pdf_print_pfile_1d
     module procedure :: pdf_print_pfile_2d
  end interface pdf_print

  interface pdf_write
     module procedure :: pdf_print_pfile_1d
     module procedure :: pdf_print_pfile_2d
  end interface pdf_write



  interface pdf_mean
     module procedure :: pdf_mean_1d
  end interface pdf_mean

  interface pdf_var
     module procedure :: pdf_var_1d
  end interface pdf_var

  interface pdf_sdev
     module procedure :: pdf_sdev_1d
  end interface pdf_sdev

  interface pdf_moment
     module procedure :: pdf_moment_1d
  end interface pdf_moment

  interface pdf_skew
     module procedure :: pdf_skew_1d
  end interface pdf_skew

  interface pdf_curt
     module procedure :: pdf_curt_1d
  end interface pdf_curt

  interface pdf_print_moments
     module procedure :: pdf_print_moments_pfile_1d
  end interface pdf_print_moments



  public :: get_moments
  public :: get_mean,get_sd,get_var,get_skew,get_curt
  public :: get_covariance

  public :: pdf_allocate
  public :: pdf_deallocate
  public :: pdf_set_range
  public :: pdf_sigma
  public :: pdf_push_sigma
  public :: pdf_accumulate
  public :: pdf_normalize
  public :: pdf_save
  public :: pdf_read
  public :: pdf_print
  public :: pdf_print_moments
  public :: pdf_mean
  public :: pdf_var
  public :: pdf_sdev
  public :: pdf_skew
  public :: pdf_curt
  public :: pdf_moment


  public :: histogram_allocate
  public :: histogram_deallocate
  public :: histogram_set_range_uniform
  public :: histogram_accumulate
  public :: histogram_get_range
  public :: histogram_get_value
  public :: histogram_print
  public :: histogram_reset


contains


  subroutine get_moments(data,ave,sdev,var,skew,curt)
    real(8), intent(out),optional     :: ave,sdev,var,skew,curt
    real(8)                           :: ave_,sdev_,var_,skew_,curt_
    real(8), dimension(:), intent(in) :: data
    integer                           :: n
    real(8)                           :: ep,adev
    real(8), dimension(size(data))    :: p,s
    n=size(data)
    if (n <= 1) then
       print*,'data.size must be at least 2'
       stop
    endif
    ave_=sum(data(:))/n
    s(:)=data(:)-ave_
    ep=sum(s(:))
    adev=sum(abs(s(:)))/n
    p(:)=s(:)*s(:)
    var_=sum(p(:))
    p(:)=p(:)*s(:)
    skew_=sum(p(:))
    p(:)=p(:)*s(:)
    curt_=sum(p(:))
    var_=(var_-ep**2/n)/(n-1)
    sdev_=sqrt(var_)
    if (var_ /= 0.0) then
       skew_=skew_/(n*sdev_**3)
       curt_=curt_/(n*var_**2)-3.0d0
    else
       skew_=0.d0
       curt_=0.d0
    endif
    if(present(ave))ave=ave_
    if(present(sdev))sdev=sdev_
    if(present(var))var=var_
    if(present(skew))skew=skew_
    if(present(curt))curt=curt_
  end subroutine get_moments


  function get_mean(data) result(mean)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: mean
    call get_moments(data,ave=mean)
  end function get_mean


  function get_sd(data) result(sd)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: sd
    call get_moments(data,sdev=sd)
  end function get_sd


  function get_var(data) result(var)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: var
    call get_moments(data,var=var)
  end function get_var


  function get_skew(data) result(skew)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: skew
    call get_moments(data,skew=skew)
  end function get_skew


  function get_curt(data) result(curt)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: curt
    call get_moments(data,curt=curt)
  end function get_curt


  function get_covariance(data,mean) result(covariance)
    real(8),dimension(:,:),intent(in)            :: data
    real(8),dimension(size(data,1)),intent(in)   :: mean
    real(8),dimension(size(data,1),size(data,1)) :: covariance
    real(8),dimension(size(data,2))              :: Xi,Xj
    integer                                      :: i,j,ncol,L
    ncol=size(data,1)
    L   =size(data,2)
    covariance=0.d0
    do i=1,ncol
       Xi(:)=data(i,:)-mean(i)
       do j=1,ncol
          Xj(:)=data(j,:)-mean(j)
          covariance(i,j) = sum(Xi(:)*Xj(:))/real(L-1,8)
       enddo
    enddo
  end function get_covariance




  !### PDF
  include "kernel_density_1d.f90"
  include "kernel_density_2d.f90"



  !### HISTOGRAMS
  include "histogram.f90"


END MODULE SF_STAT

