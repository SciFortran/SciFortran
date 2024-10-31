module SF_RANDOM
!SciFortran module for RNG
  implicit none
  private

  !MT parameters:
  integer,parameter :: defaultsd = 4357 !Default seed  
  integer,parameter :: N=624, N1=N+1    !Period parameters
  integer,save      :: mt(0:n-1)        !the array for the state vector
  integer,save      :: mti = n1         !mti==N+1 means mt[N] is not initialized

  integer           :: i,j,k,D
  real(8),parameter :: pi    = 3.14159265358979d0
  real(8),parameter :: pi2   = 6.28318530717959d0
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0

  !adam Miller parameters:
  integer,parameter :: dp = SELECTED_REAL_KIND(12, 60)
  real              :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, vsmall = TINY(1.0), vlarge = HUGE(1.0)


  !MT interface:
  interface mersenne_init
     module procedure :: init_genrand
  end interface mersenne_init

  interface mt_init
     module procedure :: init_genrand
  end interface mt_init

  interface mersenne
     module procedure :: grnd
  end interface mersenne

  interface mt_random
     module procedure :: d_grnd_1
     module procedure :: d_grnd_2
     module procedure :: d_grnd_3
     module procedure :: d_grnd_4
     module procedure :: d_grnd_5
     module procedure :: d_grnd_6
     module procedure :: d_grnd_7
     !
     module procedure :: c_grnd_1
     module procedure :: c_grnd_2
     module procedure :: c_grnd_3
     module procedure :: c_grnd_4
     module procedure :: c_grnd_5
     module procedure :: c_grnd_6
     module procedure :: c_grnd_7
  end interface mt_random

  interface mt_uniform
     module procedure :: igrnd
     module procedure :: dgrnd_uniform
  end interface mt_uniform

  interface mt_normal
     module procedure :: gaussrnd
     module procedure :: normalrnd
  end interface mt_normal

  interface mt_exponential
     module procedure :: exponentialrnd
  end interface mt_exponential

  interface mt_gamma
     module procedure :: gammarnd
  end interface mt_gamma

  interface mt_chi_square
     module procedure :: chi_squarernd
  end interface mt_chi_square

  interface mt_inverse_gamma
     module procedure :: inverse_gammarnd
  end interface mt_inverse_gamma

  interface mt_weibull
     module procedure :: weibullrnd
  end interface mt_weibull

  interface mt_cauchy
     module procedure :: cauchyrnd
  end interface mt_cauchy

  interface mt_student_t
     module procedure :: student_trnd
  end interface mt_student_t

  interface mt_laplace
     module procedure :: laplacernd
  end interface mt_laplace

  interface mt_log_normal
     module procedure :: log_normalrnd
  end interface mt_log_normal

  interface mt_beta
     module procedure :: betarnd
  end interface mt_beta

  ! Overload procedures for saving and getting mt state
  interface mt_save
     module procedure :: mtsavef
     module procedure :: mtsaveu
  end interface mt_save

  interface mt_get
     module procedure :: mtgetf
     module procedure :: mtgetu
  end interface mt_get



  public :: mersenne
  public :: mersenne_init
  !
  public :: mt_random
  public :: mt_uniform
  public :: mt_normal
  public :: mt_exponential
  public :: mt_gamma
  public :: mt_chi_square
  public :: mt_inverse_gamma
  public :: mt_weibull
  public :: mt_cauchy
  public :: mt_student_t
  public :: mt_laplace
  public :: mt_log_normal
  public :: mt_beta
  !
  public :: mt_save
  public :: mt_get
  public :: mt_init
  !
  public :: random_number_init
  public :: random_number_seed
  !
  public :: random_normal
  public :: random_gamma
  public :: random_gamma1
  public :: random_gamma2
  public :: random_chisq
  public :: random_exponential
  public :: random_Weibull
  public :: random_beta
  public :: random_inv_gauss
  public :: random_Poisson
  public :: random_binomial1
  public :: bin_prob
  public :: lngamma
  public :: random_binomial2
  public :: random_neg_binomial
  public :: random_von_Mises
  public :: random_Cauchy
  !
  public :: nrand
  public :: random_order


contains




  !+-----------------------------------------------------------------+
  !PURPOSE  : INTRINSIC RNG initialization  
  !+-----------------------------------------------------------------+
  subroutine random_number_init(shift)
    integer,optional                 :: shift
    integer                          :: i, n, clock
    integer,dimension(:),allocatable :: seed
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    call SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    if(present(shift))seed=seed+shift
    call RANDOM_SEED(PUT = seed)
    deallocate(seed)
  end subroutine random_number_init





  !------------------------------------------------------!
  !  Get the RNG seed from /dev/urandom device.          !
  !                                                      !
  !  In order to get positive seed the most              !
  !  significant bit in the number read from the         !
  !  device is cleared (by anding it with LMASK).        !
  !                                                      !
  !  NOTE: Routine uses the default integer type.        !
  !                                                      !
  !  If the device can not be opened or read routine     !
  !  falls back to calculating seed from the current     !
  !  time.                                               !
  !                                                      !
  !  Note that stream i/o is used which is a Fortran     !
  !  2003 feature.                                       !
  !                                                      !
  !  Input parameters:                                   !
  !    info : integer, if /=0 print info to stdout       !
  !    file : integer, 0: use /dev/urandom               !
  !                    1: use /dev/random                !
  !                                                      !
  !  Both parameters are optional, so that the simplest  !
  !  way to call the function is 'getseed()'.            !
  !                                                      !
  !  Generating a large amount of random numbers using   !
  !  /dev/random may be slow because quality of random   !
  !  bits from this device is guaranteed and system may  !
  !  have to wait while enough 'entropy' is collected    !
  !  from network traffic, keyboard etc.                 !
  !                                                      !
  !  A.Kuronen, antti.kuronen@helsinki.fi, 2008-2014     !
  !------------------------------------------------------!
  integer function random_number_seed(info,file)
    integer,optional,intent(in) :: info,file
    integer                     :: t(8),rn,is
    integer,parameter           :: LMASK=huge(rn) ! = 0111...111
    integer,parameter           :: LUN=676769
    character (len=80)          :: rdev0='/dev/urandom',rdev1='/dev/random',rdev
    logical                     :: openok,readok,printinfo
    openok=.true.
    readok=.true.
    if (present(file)) then
       if (file==0) then
          rdev=rdev0
       else
          rdev=rdev1
       end if
    else
       rdev=rdev0
    end if
    if (present(info)) then
       printinfo=(info/=0)
    else
       printinfo=.false.
    end if
    open(LUN,file=rdev,form='unformatted',access='stream',action='read',iostat=is)
    if (is/=0) then
       openok=.false.
       print *,'open',is
    else
       read(LUN,iostat=is) rn
       if (is/=0) then
          readok=.false.
       end if
    end if
    if (openok) close(LUN)
    if (openok.and.readok) then
       rn=iand(rn,LMASK) ! Make it positive, i.e. zero the leftmost bit
       if (printinfo) write(6,'(a,a,a,i0)') 'Seed from ',trim(rdev),': ',rn
    else
       call date_and_time(values=t)
       rn=t(7)+60*(t(6)+60*(t(5)+24*(t(3)-1+31*(t(2)-1+12*t(1)))))+t(8)
       if (printinfo) write(6,'(a,i12)') 'Seed from time:',rn
    end if
    random_number_seed=rn
    return
  end function random_number_seed






  !+-----------------------------------------------------------------+
  !PURPOSE  : MERSENNE TWISTER RNG
  !+-----------------------------------------------------------------+
  include "random_mt.f90"










  !+-----------------------------------------------------------------+
  !purpose  : rng library
  !+-----------------------------------------------------------------+
  ! A module for random number generation from the following distributions:
  !
  !     Distribution                    Function/subroutine name
  !
  !     Normal (Gaussian)               random_normal
  !     Gamma                           random_gamma
  !     Chi-squared                     random_chisq
  !     Exponential                     random_exponential
  !     Weibull                         random_Weibull
  !     Beta                            random_beta
  !     t                               random_t
  !     Multivariate normal             random_mvnorm
  !     Generalized inverse Gaussian    random_inv_gauss
  !     Poisson                         random_Poisson
  !     Binomial                        random_binomial1   *
  !                                     random_binomial2   *
  !     Negative binomial               random_neg_binomial
  !     von Mises                       random_von_Mises
  !     Cauchy                          random_Cauchy
  !
  !  Generate a random ordering of the integers 1 .. N
  !                                     random_order
  !     Initialize (seed) the uniform random number generator for ANY compiler
  !                                     seed_random_number
  !     Lognormal - see note below.
  !  ** Two functions are provided for the binomial distribution.
  !  If the parameter values remain constant, it is recommended that the
  !  first function is used (random_binomial1).   If one or both of the
  !  parameters change, use the second function (random_binomial2).
  ! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
  ! is used to provide a source of uniformly distributed random numbers.
  ! N.B. At this stage, only one random number is generated at each call to
  !      one of the functions above.
  ! The module uses the following functions which are included here:
  ! bin_prob to calculate a single binomial probability
  ! lngamma  to calculate the logarithm to base e of the gamma function
  ! Some of the code is adapted from Dagpunar's book:
  !     Dagpunar, J. 'Principles of random variate generation'
  !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
  !
  ! In most of Dagpunar's routines, there is a test to see whether the value
  ! of one or two floating-point parameters has changed since the last call.
  ! These tests have been replaced by using a logical variable FIRST.
  ! This should be set to .TRUE. on the first call using new values of the
  ! parameters, and .FALSE. if the parameter values are the same as for the
  ! previous call.
  ! Lognormal distribution
  ! If X has a lognormal distribution, then log(X) is normally distributed.
  ! Here the logarithm is the natural logarithm, that is to base e, sometimes
  ! denoted as ln.  To generate random variates from this distribution, generate
  ! a random deviate from the normal distribution with mean and variance equal
  ! to the mean and variance of the logarithms of X, then take its exponential.
  ! Relationship between the mean & variance of log(X) and the mean & variance
  ! of X, when X has a lognormal distribution.
  ! Let m = mean of log(X), and s^2 = variance of log(X)
  ! Then
  ! mean of X     = exp(m + 0.5s^2)
  ! variance of X = (mean(X))^2.[exp(s^2) - 1]
  ! In the reverse direction (rarely used)
  ! variance of log(X) = log[1 + var(X)/(mean(X))^2]
  ! mean of log(X)     = log(mean(X) - 0.5var(log(X))
  ! N.B. The above formulae relate to population parameters; they will only be
  !      approximate if applied to sample values.
  ! Version 1.13, 2 October 2000
  ! Changes from version 1.01
  ! 1. The random_order, random_Poisson & random_binomial routines have been
  !    replaced with more efficient routines.
  ! 2. A routine, seed_random_number, has been added to seed the uniform random
  !    number generator.   This requires input of the required number of seeds
  !    for the particular compiler from a specified I/O unit such as a keyboard.
  ! 3. Made compatible with Lahey's ELF90.
  ! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
  ! 5. INTENT for array f corrected in random_mvnorm.
  include "random_routines.f90"




  !+-----------------------------------------------------------------+
  !PURPOSE:  Numerical Recipes
  !+-----------------------------------------------------------------+
  real(8) function nrand(dseed_)
    integer,optional        :: dseed_
    integer                 :: dseed
    integer,parameter       :: IM1=2147483563, IM2=2147483399, IMM1=IM1-1, IA1=40014, &
         & IA2=40692, IQ1=53668, IQ2=52774, IR1=12211, IR2=3791,  &
         & NTAB=32, NDIV=1+IMM1/NTAB
    real(kind=8), parameter :: AM=1.0d0/IM1, EPS=1.2e-7, RNMX=1.-EPS
    integer                 :: dseed2, j, k, iv(NTAB), iy
    save iv, iy, dseed2
    data dseed2/123456789/, iv/NTAB*0/, iy/0/
    dseed = 123456 ;if(present(dseed_))dseed=dseed_
    if(dseed .le. 0) then
       dseed = max(-dseed,1)
       dseed2 = dseed
       do j=NTAB+8, 1, -1
          k = dseed/IQ1
          dseed = IA1*(dseed-k*IQ1)-k*IR1
          if(dseed .lt. 0) dseed = dseed+IM1
          if(j .le. NTAB) iv(j) = dseed
       enddo
       iy=iv(1)
    endif
    k = dseed/IQ1
    dseed = IA1*(dseed-k*IQ1)-k*IR1
    if(dseed .lt. 0) dseed = dseed+IM1
    k = dseed2/IQ2
    dseed2 = IA2*(dseed2-k*IQ2)-k*IR2
    if(dseed2 .lt. 0) dseed2 = dseed2+IM2
    j = 1+iy/NDIV
    iy = iv(j)-dseed2
    iv(j) = dseed
    if(iy .lt. 1) iy = iy+IMM1
    nrand = min(AM*iy,RNMX)
  end function nrand





  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  SUBROUTINE random_order(order, n)
    !     generate a random ordering of the integers 1 ... n.
    integer, intent(in)  :: n
    integer, intent(out) :: order(n)
    !     local variables
    integer :: i, j, k
    real(8) :: wk
    do i = 1, n
       order(i) = i
    end do
    !     starting at the end, swap the current last indicator with one
    !     randomly chosen from those preceeding it.
    do i = n, 2, -1
       call random_number(wk)
       j = 1 + i * wk
       if (j < i) then
          k = order(i)
          order(i) = order(j)
          order(j) = k
       end if
    end do
    return
  end subroutine random_order


end module SF_RANDOM
