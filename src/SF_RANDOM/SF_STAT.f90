MODULE SF_STAT
!SciFortran module for statistical estimators and visualization
  USE SF_ARRAYS, only: linspace
  USE SF_INTEGRATE, only: simps
  USE SF_IOTOOLS, only: free_unit,splot3d
  USE SF_LINALG, only: det,inv
  implicit none
  private

  type histogram
     !This derived type represents a histogram of :f:var:`n` bins. The bin :math:`i`, with :math:`i=0,\dots,n-1`, covers the
     !interval :code:`[range(i),range(i+1))` and contains in :code:`bin(i)` the sum of the weights of the values accumulated
     !in it. It is created by :f:func_inline:`histogram_allocate` and handled by the :code:`histogram_*` routines.
     !
     integer                       :: n=0    ! number of bins
     real(8),dimension(:),pointer  :: range  ! edges of the bins, range(0:n); bin i is [range(i),range(i+1))
     real(8),dimension(:),pointer  :: bin    ! content of the bins, bin(0:n); only bin(0:n-1) is used
  end type histogram

  type pdf_kernel
     !This derived type holds the probability density function (pdf) of a one-dimensional variable, estimated with the
     !Gaussian kernel density method. The pdf is represented on a uniform grid :f:var:`x` of :f:var:`N` points between
     !:f:var:`xmin` and :f:var:`xmax`: every accumulated data point adds to it a Gaussian of width :f:var:`sigma` centered on
     !the data point. See :f:func_inline:`pdf_allocate`, :f:func_inline:`pdf_set_range`, :f:func_inline:`pdf_accumulate` and
     !:f:func_inline:`pdf_normalize`.
     !
     integer                          :: N=0               ! number of grid points
     real(8)                          :: xmin              ! lower bound of the grid
     real(8)                          :: xmax              ! upper bound of the grid
     real(8)                          :: dx                ! grid step
     integer                          :: Ndata=0           ! number of data points accumulated
     real(8),dimension(:),allocatable :: x                 ! grid, size N
     real(8),dimension(:),allocatable :: pdf               ! pdf on the grid, size N
     real(8)                          :: sigma             ! width (standard deviation) of the Gaussian kernel
     logical                          :: status=.false.    ! .true. if the object is allocated
     logical                          :: variance=.false.  ! .true. if sigma has been stored with pdf_push_sigma
     logical                          :: rescale=.false.   ! .true. once the pdf has been normalized
  end type pdf_kernel

  type pdf_kernel_2d
     !This derived type holds the probability density function (pdf) of a two-dimensional variable, estimated with the
     !Gaussian kernel density method. The pdf is represented on a uniform grid of :code:`N(1)` x :code:`N(2)` points, defined
     !by the grids :f:var:`x` and :f:var:`y` between :f:var:`xmin` and :f:var:`xmax`: every accumulated data point adds to it
     !a bivariate Gaussian with covariance matrix :f:var:`Sigma` centered on the data point. See
     !:f:func_inline:`pdf_allocate`, :f:func_inline:`pdf_set_range`, :f:func_inline:`pdf_accumulate` and
     !:f:func_inline:`pdf_normalize`.
     !
     integer,dimension(2)               :: N                 ! number of grid points along x and y
     real(8),dimension(2)               :: xmin              ! lower bounds of the grid, (x,y)
     real(8),dimension(2)               :: xmax              ! upper bounds of the grid, (x,y)
     real(8),dimension(2)               :: dx                ! grid steps, (x,y)
     integer                            :: Ndata=0           ! number of data points accumulated
     real(8),dimension(:),allocatable   :: x                 ! grid along x, size N(1)
     real(8),dimension(:),allocatable   :: y                 ! grid along y, size N(2)
     real(8),dimension(:,:),allocatable :: pdf               ! pdf on the grid, [N(1),N(2)]
     real(8),dimension(2,2)             :: Sigma             ! covariance matrix of the Gaussian kernel
     logical                            :: status=.false.    ! .true. if the object is allocated
     logical                            :: variance=.false.  ! .true. if Sigma has been stored with pdf_push_sigma
     logical                            :: rescale=.false.   ! .true. once the pdf has been normalized
  end type pdf_kernel_2d


  interface pdf_allocate
     !Generic interface to allocate a pdf object: a :f:var:`pdf_kernel` with :f:var:`N` grid points
     !(:f:func_inline:`pdf_allocate_1d`), or a :f:var:`pdf_kernel_2d` with a grid of :code:`Nvec(1)` x :code:`Nvec(2)` points
     !(:f:func_inline:`pdf_allocate_2d`). The pdf and the kernel width are set to zero, no data is accumulated, and the object
     !is flagged as allocated. The grid is not defined until :f:func_inline:`pdf_set_range` is called. The program stops if
     !the object is already allocated.
     !
     module procedure :: pdf_allocate_1d
     module procedure :: pdf_allocate_2d
  end interface pdf_allocate

  interface pdf_deallocate
     !Generic interface to deallocate a pdf object, :f:func_inline:`pdf_deallocate_1d` or :f:func_inline:`pdf_deallocate_2d`:
     !the arrays are freed and the ranges, the kernel width and the number of data are reset to zero. The program stops if the
     !object is not allocated.
     !
     module procedure :: pdf_deallocate_1d
     module procedure :: pdf_deallocate_2d
  end interface pdf_deallocate


  interface pdf_save
     !Generic interface to save a pdf object, :f:func_inline:`pdf_save_1d` or :f:func_inline:`pdf_save_2d`, to the formatted
     !file :f:var:`pfile`. The file is overwritten, and it contains one record per item: the number of grid points, the
     !ranges, the grid steps, the number of data, the grids, the pdf, the kernel width and the three flags of the object. It
     !can be read back with :f:func_inline:`pdf_read`. The program stops if the object is not allocated.
     !
     module procedure :: pdf_save_1d
     module procedure :: pdf_save_2d
  end interface pdf_save


  interface pdf_read
     !Generic interface to read a pdf object, :f:func_inline:`pdf_read_1d` or :f:func_inline:`pdf_read_2d`, from a file
     !written by :f:func_inline:`pdf_save`. The object is allocated by the routine, therefore it must not be allocated on
     !input, otherwise the program stops.
     !
     module procedure :: pdf_read_1d
     module procedure :: pdf_read_2d
  end interface pdf_read


  interface pdf_set_range
     !Generic interface to define the grid of a pdf object, as a uniform mesh between the lower bound :f:var:`a` and the upper
     !bound :f:var:`b`: a scalar pair for :f:func_inline:`pdf_set_range_1d`, vectors with the bounds along x and y for
     !:f:func_inline:`pdf_set_range_2d`. It sets the ranges and the grid steps of the object. The program stops if the object
     !is not allocated.
     !
     module procedure :: pdf_set_range_1d
     module procedure :: pdf_set_range_2d
  end interface pdf_set_range


  interface pdf_push_sigma
     !Generic interface to store the width of the Gaussian kernel in the pdf object: the standard deviation :f:var:`sigma` in
     !one dimension (:f:func_inline:`pdf_push_sigma_1d`), the :math:`2\times 2` covariance matrix in two dimensions
     !(:f:func_inline:`pdf_push_sigma_2d`). Once stored, it is used by :f:func_inline:`pdf_accumulate`, and the argument
     !:code:`sigma` of the latter is ignored. The program stops if the object is not allocated.
     !
     module procedure :: pdf_push_sigma_1d
     module procedure :: pdf_push_sigma_2d
  end interface pdf_push_sigma


  interface pdf_get_sigma
     module procedure :: pdf_get_sigma_1d
     module procedure :: pdf_get_sigma_2d
  end interface pdf_get_sigma


  interface pdf_sigma
     !Generic interface to estimate the width of the Gaussian kernel with Silverman's rule of thumb. The result is returned in
     !:f:var:`h`, the pdf object is not modified (but it must be allocated). In one dimension :math:`h = (4/3N)^{1/5}\sigma`,
     !with :math:`N` the number of data and :math:`\sigma` their standard deviation; in two dimensions the result is the
     !diagonal covariance matrix :math:`h_{ii} = N^{-1/3}\sigma_i^2`.
     !
     !The standard deviation is computed from the data (:f:func_inline:`pdf_sigma_data_1d`, and
     !:f:func_inline:`pdf_sigma_data_2d` for data of shape :code:`[2,L]`) or it is given together with the number of data
     !(:f:func_inline:`pdf_sigma_sdev_1d`, :f:func_inline:`pdf_sigma_sdev_2d`). The program stops if the object is not
     !allocated or, in two dimensions, if the data do not have shape :code:`[2,L]`.
     !
     module procedure :: pdf_sigma_data_1d
     module procedure :: pdf_sigma_sdev_1d
     module procedure :: pdf_sigma_data_2d
     module procedure :: pdf_sigma_sdev_2d
  end interface pdf_sigma


  interface pdf_accumulate
     !Generic interface to accumulate data in a pdf object. Every data point adds to the pdf on the grid a Gaussian kernel
     !centered on it, and the number of data of the object is increased. In one dimension the kernel is
     !:math:`e^{-(x-m)^2/2\sigma^2}/\sqrt{2\pi}\sigma` and the data is a scalar (:f:func_inline:`pdf_accumulate_s_1d`) or a
     !vector, whose elements are accumulated one by one (:f:func_inline:`pdf_accumulate_v_1d`). In two dimensions it is
     !:math:`e^{-\frac{1}{2}(x-m)^T\Sigma^{-1}(x-m)}/2\pi\sqrt{\det\Sigma}` and the data is a point :code:`(x,y)`
     !(:f:func_inline:`pdf_accumulate_s_2d`).
     !
     !The width of the kernel is the one stored with :f:func_inline:`pdf_push_sigma`; if none has been stored it must be
     !passed in the optional argument :code:`sigma`, otherwise the program stops. The program stops also if the object is not
     !allocated.
     !
     module procedure :: pdf_accumulate_s_1d
     module procedure :: pdf_accumulate_v_1d
     module procedure :: pdf_accumulate_s_2d
  end interface pdf_accumulate



  interface pdf_normalize
     !Generic interface to normalize a pdf object, :f:func_inline:`pdf_normalize_1d` or :f:func_inline:`pdf_normalize_2d`. The
     !pdf is divided by the number of accumulated data (the first time only) and then rescaled so that its integral on the
     !grid is one, i.e. :math:`\sum_i p_i \Delta x = 1` (:math:`\sum_{ij} p_{ij} \Delta x \Delta y = 1` in two dimensions). It
     !is called by :f:func_inline:`pdf_print` unless otherwise requested. The program stops if the object is not allocated.
     !
     module procedure :: pdf_normalize_1d
     module procedure :: pdf_normalize_2d
  end interface pdf_normalize


  interface pdf_print
     !Generic interface to write a pdf object to the file :f:var:`pfile`. The pdf is normalized first, with
     !:f:func_inline:`pdf_normalize`, unless :f:var:`normalize` is :code:`.false.`. In one dimension the two columns
     !:code:`x pdf` are appended to the file, followed by a blank line (:f:func_inline:`pdf_print_pfile_1d`); in two
     !dimensions the pdf is written with :f:func_inline:`splot3d` (:f:func_inline:`pdf_print_pfile_2d`). The program stops if
     !the object is not allocated.
     !
     module procedure :: pdf_print_pfile_1d
     module procedure :: pdf_print_pfile_2d
  end interface pdf_print

  interface pdf_write
     module procedure :: pdf_print_pfile_1d
     module procedure :: pdf_print_pfile_2d
  end interface pdf_write



  interface pdf_mean
     !Generic interface for the mean of a one-dimensional pdf, :math:`\int x\,p(x)\,dx`, integrated with the Simpson's rule on
     !the grid (:f:func_inline:`pdf_mean_1d`).
     !
     module procedure :: pdf_mean_1d
  end interface pdf_mean

  interface pdf_var
     !Generic interface for the variance of a one-dimensional pdf, :math:`\int (x-\mu)^2\,p(x)\,dx`, integrated with the
     !Simpson's rule on the grid (:f:func_inline:`pdf_var_1d`).
     !
     module procedure :: pdf_var_1d
  end interface pdf_var

  interface pdf_sdev
     !Generic interface for the standard deviation of a one-dimensional pdf, the square root of :f:func_inline:`pdf_var`
     !(:f:func_inline:`pdf_sdev_1d`).
     !
     module procedure :: pdf_sdev_1d
  end interface pdf_sdev

  interface pdf_moment
     !Generic interface for the moment of order :math:`n` of a one-dimensional pdf around a point :math:`\mu`,
     !:math:`\int (x-\mu)^n\,p(x)\,dx`, integrated with the Simpson's rule on the grid (:f:func_inline:`pdf_moment_1d`). The
     !point :math:`\mu` is optional and it is zero by default, i.e. the moment is the raw one.
     !
     module procedure :: pdf_moment_1d
  end interface pdf_moment

  interface pdf_skew
     !Generic interface for the skewness of a one-dimensional pdf, the third standardized moment
     !:math:`\int (x-\mu)^3 p(x)\,dx/\sigma^3` (:f:func_inline:`pdf_skew_1d`).
     !
     module procedure :: pdf_skew_1d
  end interface pdf_skew

  interface pdf_curt
     !Generic interface for the kurtosis of a one-dimensional pdf, the fourth standardized moment
     !:math:`\int (x-\mu)^4 p(x)\,dx/\sigma^4`, NOT the excess kurtosis (:f:func_inline:`pdf_curt_1d`).
     !
     module procedure :: pdf_curt_1d
  end interface pdf_curt

  interface pdf_print_moments
     !Generic interface to write to the file :f:var:`pfile` the moments of a one-dimensional pdf, on a single line: mean,
     !standard deviation, skewness and kurtosis (:f:func_inline:`pdf_print_moments_pfile_1d`).
     !
     module procedure :: pdf_print_moments_pfile_1d
  end interface pdf_print_moments


  public :: histogram
  public :: pdf_kernel
  public :: pdf_kernel_2d

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
  !This subroutine computes the mean, the standard deviation, the variance, the skewness and the kurtosis of the real array
  !:f:var:`data`. All the outputs are optional: only the requested ones are returned.
  !
  !The variance is the unbiased estimator, with :math:`n-1` in the denominator, computed with the corrected two-pass
  !algorithm, and the standard deviation is its square root. The skewness is the third central moment divided by
  !:math:`n\sigma^3`. The kurtosis is the EXCESS kurtosis, the fourth central moment divided by :math:`n\sigma^4` minus 3,
  !which is zero for a normal distribution. If the variance is zero the skewness and the kurtosis are set to zero. The program
  !stops if :f:var:`data` has less than two elements.
  !
    real(8), intent(out),optional     :: ave   ! mean
    real(8), intent(out),optional     :: sdev  ! standard deviation (unbiased)
    real(8), intent(out),optional     :: var   ! variance (unbiased)
    real(8), intent(out),optional     :: skew  ! skewness
    real(8), intent(out),optional     :: curt  ! excess kurtosis
    real(8)                           :: ave_,sdev_,var_,skew_,curt_
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
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
  !This function returns the mean of the real array :f:var:`data`, computed with :f:func_inline:`get_moments`. The program
  !stops if :f:var:`data` has less than two elements.
  !
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
    real(8)                           :: mean  ! mean
    call get_moments(data,ave=mean)
  end function get_mean


  function get_sd(data) result(sd)
  !This function returns the standard deviation, with :math:`n-1` in the denominator, of the real array :f:var:`data`,
  !computed with :f:func_inline:`get_moments`. The program stops if :f:var:`data` has less than two elements.
  !
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
    real(8)                           :: sd    ! standard deviation
    call get_moments(data,sdev=sd)
  end function get_sd


  function get_var(data) result(var)
  !This function returns the variance, with :math:`n-1` in the denominator, of the real array :f:var:`data`, computed with
  !:f:func_inline:`get_moments`. The program stops if :f:var:`data` has less than two elements.
  !
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
    real(8)                           :: var   ! variance
    call get_moments(data,var=var)
  end function get_var


  function get_skew(data) result(skew)
  !This function returns the skewness of the real array :f:var:`data`, computed with :f:func_inline:`get_moments`. The program
  !stops if :f:var:`data` has less than two elements.
  !
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
    real(8)                           :: skew  ! skewness
    call get_moments(data,skew=skew)
  end function get_skew


  function get_curt(data) result(curt)
  !This function returns the excess kurtosis of the real array :f:var:`data`, computed with :f:func_inline:`get_moments`. The
  !program stops if :f:var:`data` has less than two elements.
  !
    real(8), dimension(:), intent(in) :: data  ! sample, size n >= 2
    real(8)                           :: curt  ! excess kurtosis
    call get_moments(data,curt=curt)
  end function get_curt


  function get_covariance(data,mean) result(covariance)
  !This function returns the covariance matrix of the variables stored as the rows of :f:var:`data`, for :math:`L` samples
  !stored as the columns: :math:`C_{ij} = \frac{1}{L-1}\sum_k (x_{ik}-\bar{x}_i)(x_{jk}-\bar{x}_j)`, with :math:`\bar{x}_i`
  !the means given in :f:var:`mean`.
  !
    real(8),dimension(:,:),intent(in)            :: data        ! samples, [ncol,L]: ncol variables (rows), L samples (columns)
    real(8),dimension(size(data,1)),intent(in)   :: mean        ! mean of each variable, size ncol
    real(8),dimension(size(data,1),size(data,1)) :: covariance  ! covariance matrix, [ncol,ncol]
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

