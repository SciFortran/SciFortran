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
  real(8),parameter :: pi    = 3.14159265358979d0           ! pi
  real(8),parameter :: pi2   = 6.28318530717959d0           ! 2 pi
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0  ! square root of 2
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0  ! square root of 3
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0  ! square root of 6

  !adam Miller parameters:
  integer,parameter :: dp = SELECTED_REAL_KIND(12, 60)      ! kind of the double precision reals of the Alan Miller routines
  real              :: zero = 0.0                           ! constant 0, single precision
  real              :: half = 0.5                           ! constant 0.5, single precision
  real              :: one = 1.0                            ! constant 1, single precision
  real              :: two = 2.0                            ! constant 2, single precision
  real              :: vsmall = TINY(1.0)                   ! smallest positive normalized real, tiny(1.0)
  real              :: vlarge = HUGE(1.0)                   ! largest real, huge(1.0)


  !MT interface:
  interface mersenne_init
     !Generic interface to initialize the Mersenne Twister random number generator with the integer :f:var:`seed`
     !(:f:func_inline:`init_genrand`), same as :f:func_inline:`mt_init`. The 624 words of the state are filled from the seed
     !with the linear recurrence of Knuth (multiplier 1812433253), assuming that the integer overflow does not stop the
     !program. If the generator is not initialized, the first call to :f:func_inline:`mersenne` uses the default seed 4357.
     !
     module procedure :: init_genrand
  end interface mersenne_init

  interface mt_init
     !Generic interface to initialize the Mersenne Twister random number generator with the integer :f:var:`seed`
     !(:f:func_inline:`init_genrand`), same as :f:func_inline:`mersenne_init`. See :f:func_inline:`mersenne_init` for the
     !details.
     !
     module procedure :: init_genrand
  end interface mt_init

  interface mersenne
     !Generic interface to the Mersenne Twister random number generator MT19937 (:f:func_inline:`grnd`). Every call returns a
     !real(8) random number uniformly distributed between 0 and 1, obtained dividing a 32 bit random integer by
     !:math:`2^{32}-1`. The 624 words of the state are regenerated every 624 calls. The generator is initialized with
     !:f:func_inline:`mersenne_init`, and its state can be saved and read back with :f:func_inline:`mt_save` and
     !:f:func_inline:`mt_get`. All the :code:`mt_*` generators use it as the source of uniform numbers.
     !
     module procedure :: grnd
  end interface mersenne

  interface mt_random
     !Generic interface to fill an array of rank 1 to 7 with uniform random numbers between 0 and 1, generated with
     !:f:func_inline:`mersenne`. The array is real(8) (:f:func_inline:`d_grnd_1` to :f:func_inline:`d_grnd_7`) or complex(8)
     !(:code:`c_grnd_1` to :code:`c_grnd_7`), in which case the real and the imaginary parts are two independent numbers.
     !
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
     !Generic interface for uniform random numbers generated with :f:func_inline:`mersenne`: an integer between :f:var:`l` and
     !:f:var:`h`, both included (:f:func_inline:`igrnd`), or a real(8) number between :f:var:`a` and :f:var:`b`
     !(:f:func_inline:`dgrnd_uniform`).
     !
     module procedure :: igrnd
     module procedure :: dgrnd_uniform
  end interface mt_uniform

  interface mt_normal
     !Generic interface for normally distributed random numbers: a standard normal deviate, with zero mean and unit standard
     !deviation (:f:func_inline:`gaussrnd`), or a normal deviate with the given mean and standard deviation
     !(:f:func_inline:`normalrnd`). The former uses the polar form of the Box-Muller method (Numerical Recipes), generating
     !two deviates at a time and returning the second one at the following call; the latter uses the basic Box-Muller method,
     !one deviate per call.
     !
     module procedure :: gaussrnd
     module procedure :: normalrnd
  end interface mt_normal

  interface mt_exponential
     !Generic interface for random numbers from an exponential distribution with mean :f:var:`mean`, positive, obtained by
     !inversion (:f:func_inline:`exponentialrnd`). If the mean is not positive a message is printed and the result is
     !undefined.
     !
     module procedure :: exponentialrnd
  end interface mt_exponential

  interface mt_gamma
     !Generic interface for random numbers from a gamma distribution with parameters :f:var:`shape` and :f:var:`scale`, both
     !positive, i.e. with mean :math:`shape\cdot scale` (:f:func_inline:`gammarnd`). It uses the method of Marsaglia and Tsang
     !(ACM TOMS 26, 363, 2000), and for :math:`shape<1` the relation with a gamma deviate of shape :math:`shape+1`. If a
     !parameter is not positive a message is printed and the result is undefined.
     !
     module procedure :: gammarnd
  end interface mt_gamma

  interface mt_chi_square
     !Generic interface for random numbers from a chi-square distribution with :f:var:`dof` degrees of freedom
     !(:f:func_inline:`chi_squarernd`), obtained from :f:func_inline:`mt_gamma`.
     !
     module procedure :: chi_squarernd
  end interface mt_chi_square

  interface mt_inverse_gamma
     !Generic interface for random numbers from an inverse gamma distribution with parameters :f:var:`shape` and
     !:f:var:`scale` (:f:func_inline:`inverse_gammarnd`), obtained as the inverse of a gamma deviate,
     !:f:func_inline:`mt_gamma`, with shape :f:var:`shape` and scale :math:`1/scale`.
     !
     module procedure :: inverse_gammarnd
  end interface mt_inverse_gamma

  interface mt_weibull
     !Generic interface for random numbers from a Weibull distribution with parameters :f:var:`shape` and :f:var:`scale`, both
     !positive (:f:func_inline:`weibullrnd`), obtained by inversion, :math:`scale\,(-\ln u)^{1/shape}`. If a parameter is not
     !positive a message is printed and the result is undefined.
     !
     module procedure :: weibullrnd
  end interface mt_weibull

  interface mt_cauchy
     !Generic interface for random numbers from a Cauchy distribution with median :f:var:`median` and scale :f:var:`scale`,
     !positive (:f:func_inline:`cauchyrnd`), obtained by inversion. If the scale is not positive a message is printed and the
     !result is undefined.
     !
     module procedure :: cauchyrnd
  end interface mt_cauchy

  interface mt_student_t
     !Generic interface for random numbers from a Student t distribution with :f:var:`dof` degrees of freedom, positive
     !(:f:func_inline:`student_trnd`), obtained as the ratio of a standard normal deviate, :f:func_inline:`mt_normal`, and the
     !square root of a chi-square deviate, :f:func_inline:`mt_chi_square`, divided by the degrees of freedom. If :f:var:`dof`
     !is not positive a message is printed and the result is undefined.
     !
     module procedure :: student_trnd
  end interface mt_student_t

  interface mt_laplace
     !Generic interface for random numbers from a Laplace (double exponential) distribution with mean :f:var:`mean` and scale
     !:f:var:`scale`, positive (:f:func_inline:`laplacernd`), obtained by inversion. If the scale is not positive a message is
     !printed and the result is undefined.
     !
     module procedure :: laplacernd
  end interface mt_laplace

  interface mt_log_normal
     !Generic interface for random numbers from a log-normal distribution (:f:func_inline:`log_normalrnd`), i.e. the
     !exponential of a normal deviate, :f:func_inline:`mt_normal`, with mean :f:var:`mu` and standard deviation
     !:f:var:`sigma`. Note that :f:var:`mu` and :f:var:`sigma` are the mean and standard deviation of the logarithm of the
     !variable, not of the variable.
     !
     module procedure :: log_normalrnd
  end interface mt_log_normal

  interface mt_beta
     !Generic interface for random numbers from a beta distribution with positive shape parameters :f:var:`a` and :f:var:`b`
     !(:f:func_inline:`betarnd`), obtained as :math:`u/(u+v)` with :math:`u` and :math:`v` gamma deviates,
     !:f:func_inline:`mt_gamma`, of shapes :f:var:`a` and :f:var:`b`. If a parameter is not positive a message is printed and
     !the result is undefined.
     !
     module procedure :: betarnd
  end interface mt_beta

  ! Overload procedures for saving and getting mt state
  interface mt_save
     !Generic interface to save the state of the Mersenne Twister generator, i.e. the index and the 624 words of the state, to
     !a file (:f:func_inline:`mtsavef`) or to an open unit (:f:func_inline:`mtsaveu`). With :f:var:`forma` equal to
     !:code:`'u'` or :code:`'U'` the state is written unformatted, otherwise formatted. NOTE: a file is not overwritten, the
     !state is APPENDED to it, and the unit 10 is used. The state can be read back with :f:func_inline:`mt_get`.
     !
     module procedure :: mtsavef
     module procedure :: mtsaveu
  end interface mt_save

  interface mt_get
     !Generic interface to read the state of the Mersenne Twister generator, saved with :f:func_inline:`mt_save`, from a file
     !(:f:func_inline:`mtgetf`) or from an open unit (:f:func_inline:`mtgetu`). With :f:var:`forma` equal to :code:`'u'` or
     !:code:`'U'` the state is read unformatted, otherwise formatted. The file must exist and the first saved state is read;
     !the unit 10 is used for the file.
     !
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
  !This subroutine initializes the random number generator of the compiler, the one used by the intrinsic
  !:code:`random_number`, with a seed array built from the system clock, :code:`seed(i) = clock + 37*(i-1)`. The optional
  !argument :f:var:`shift` is added to every element of the seed, and it can be used to make different the seeds of processes
  !started at the same time.
  !
    integer,optional                 :: shift  ! optional: integer added to every element of the seed
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
  !This function returns a positive integer to be used as the seed of a random number generator, read from the device
  !:code:`/dev/urandom` (default) or :code:`/dev/random`, with stream i/o. The most significant bit is cleared to get a
  !positive number. If the device cannot be opened or read, the seed is computed from the current date and time. The simplest
  !way to call it is :code:`random_number_seed()`. Note that :code:`/dev/random` may be slow, since it can wait for enough
  !entropy.
  !
    integer,optional,intent(in) :: info  ! optional: if /= 0 print the origin of the seed to the standard output
    integer,optional,intent(in) :: file  ! optional: 0 use /dev/urandom (default), /= 0 use /dev/random
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
  !This function returns a uniform random number in the open interval (0,1), generated with the algorithm :code:`ran2` of
  !Numerical Recipes, i.e. the combination of two congruential generators of L'Ecuyer with a Bays-Durham shuffle. The optional
  !seed :f:var:`dseed_` initializes the generator if it is not positive, with the seed :code:`-dseed_` (or 1). Only the second
  !generator and the shuffle table are saved between calls, the state of the first one is not.
  !
    integer,optional        :: dseed_  ! optional: a value <= 0 initializes the generator with seed :code:`-dseed_`
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
  !This subroutine generates a random permutation of the integers 1,...,:f:var:`n`, with the Fisher-Yates shuffle: starting
  !from the end, the last element is swapped with one randomly chosen among the preceding ones. The random numbers are
  !generated with the intrinsic :code:`random_number`.
  !
    !     generate a random ordering of the integers 1 ... n.
    integer, intent(in)  :: n         ! number of integers
    integer, intent(out) :: order(n)  ! random permutation of 1,...,n, size n
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
