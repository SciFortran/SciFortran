!---------------------------------------------------------------!
! Initialization subroutine                                     !
!---------------------------------------------------------------!
subroutine sgrnd(seed)
  !This subroutine initializes the state of the Mersenne Twister generator with the congruential generator of Knuth
  !(multiplier 69069). It is used by :f:func_inline:`grnd` with the default seed 4357 if the generator has not been
  !initialized.
  !
  ! Setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
  ! [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]
  integer,intent(in) :: seed  ! seed of the generator
  mt(0) = iand(seed,-1)
  do mti=1,N-1
     mt(mti) = iand(69069 * mt(mti-1),-1)
  enddo
  return
end subroutine sgrnd

subroutine init_genrand(seed)
  ! This initialization is based upon the multiplier given on p.106 of the
  ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
  ! This version assumes that integer overflow does NOT cause a crash.
  integer,intent(in) :: seed  ! seed of the generator
  integer             :: latest
  mt(0) = seed
  latest = seed
  do mti = 1,n-1
     latest = IEOR( latest,ISHFT( latest,-30 ) )
     latest = latest * 1812433253 + mti
     mt(mti) = latest
  end do
  return
end subroutine init_genrand





!---------------------------------------------------------------!
! Random number generator: [0,1[                                !
!---------------------------------------------------------------!
function grnd()
  real(8)            :: grnd  ! random number in [0,1]
  integer,parameter :: M = 397, MATA  = -1727483681 ! constant vector a
  integer,parameter :: LMASK =  2147483647          ! least significant r bits
  integer,parameter :: UMASK = -LMASK - 1           ! most significant w-r bits
  integer,parameter :: TMASKB= -1658038656, TMASKC= -272236544    
  integer,save       :: mag01(0:1)=[0,MATA] ! mag01(x) = x * MATA for x=0,1
  integer            :: kk,y
  ! These statement functions have been replaced with separate functions
  ! tshftu(y) = ISHFT(y,-11)
  ! tshfts(y) = ISHFT(y,7)
  ! tshftt(y) = ISHFT(y,15)
  ! tshftl(y) = ISHFT(y,-18)
  if (mti>=N) then               ! generate N words at one time
     if (mti==N+1) then          ! if sgrnd() has not been called,
        call sgrnd( defaultsd ) ! a default initial seed is used
     endif
     do kk=0,N-M-1
        y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
        mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
     enddo
     do kk=N-M,N-2
        y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
        mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
     enddo
     y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
     mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
     mti = 0
  endif
  y=mt(mti)
  mti = mti + 1 
  y=ieor(y,TSHFTU(y))
  y=ieor(y,iand(TSHFTS(y),TMASKB))
  y=ieor(y,iand(TSHFTT(y),TMASKC))
  y=ieor(y,TSHFTL(y))    
  if (y<0) then
     grnd=(dble(y)+2d0**32)/(2d0**32-1d0)
  else
     grnd=dble(y)/(2d0**32-1d0)
  endif
  return
contains
  integer function TSHFTU(y)
    integer,intent(in) :: y
    TSHFTU=ishft(y,-11)
    return
  end function TSHFTU
  integer function TSHFTS(y)
    integer,intent(in) :: y
    TSHFTS=ishft(y,7)
    return
  end function TSHFTS
  integer function TSHFTT(y)
    integer,intent(in) :: y
    TSHFTT=ishft(y,15)
    return
  end function TSHFTT
  integer function TSHFTL(y)
    integer,intent(in) :: y
    TSHFTL=ishft(y,-18)
    return
  end function TSHFTL
end function grnd





subroutine d_grnd_1(A)
  real(8),dimension(:) :: A  ! real(8) array of rank 1, filled with random numbers
  integer              :: i
  do i=1,size(A)
     A(i) = mersenne()
  enddo
end subroutine d_grnd_1

subroutine d_grnd_2(A)
  real(8),dimension(:,:) :: A  ! real(8) array of rank 2, filled with random numbers
  integer              :: i1,i2
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        A(i1,i2) = mersenne()
     enddo
  enddo
end subroutine d_grnd_2

subroutine d_grnd_3(A)
  real(8),dimension(:,:,:) :: A  ! real(8) array of rank 3, filled with random numbers
  integer                  :: i1,i2,i3
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           A(i1,i2,i3) = mersenne()
        enddo
     enddo
  enddo
end subroutine d_grnd_3

subroutine d_grnd_4(A)
  real(8),dimension(:,:,:,:) :: A  ! real(8) array of rank 4, filled with random numbers
  integer                  :: i1,i2,i3,i4
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              A(i1,i2,i3,i4) = mersenne()
           enddo
        enddo
     enddo
  enddo
end subroutine d_grnd_4

subroutine d_grnd_5(A)
  real(8),dimension(:,:,:,:,:) :: A  ! real(8) array of rank 5, filled with random numbers
  integer                  :: i1,i2,i3,i4,i5
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 A(i1,i2,i3,i4,i5) = mersenne()
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine d_grnd_5

subroutine d_grnd_6(A)
  real(8),dimension(:,:,:,:,:,:) :: A  ! real(8) array of rank 6, filled with random numbers
  integer                  :: i1,i2,i3,i4,i5,i6
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 do i6=1,size(A,6)
                    A(i1,i2,i3,i4,i5,i6) = mersenne()
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine d_grnd_6


subroutine d_grnd_7(A)
  real(8),dimension(:,:,:,:,:,:,:) :: A  ! real(8) array of rank 7, filled with random numbers
  integer                  :: i1,i2,i3,i4,i5,i6,i7
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 do i6=1,size(A,6)
                    do i7=1,size(A,7)
                       A(i1,i2,i3,i4,i5,i6,i7) = mersenne()
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine d_grnd_7










subroutine c_grnd_1(A)
  complex(8),dimension(:) :: A  ! complex(8) array of rank 1, filled with random numbers
  integer                 :: i
  do i=1,size(A)
     A(i) = dcmplx(mersenne(),mersenne())
  enddo
end subroutine c_grnd_1

subroutine c_grnd_2(A)
  complex(8),dimension(:,:) :: A  ! complex(8) array of rank 2, filled with random numbers
  integer                   :: i1,i2
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        A(i1,i2) = dcmplx(mersenne(),mersenne())
     enddo
  enddo
end subroutine c_grnd_2

subroutine c_grnd_3(A)
  complex(8),dimension(:,:,:) :: A  ! complex(8) array of rank 3, filled with random numbers
  integer                     :: i1,i2,i3
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           A(i1,i2,i3) = dcmplx(mersenne(),mersenne())
        enddo
     enddo
  enddo
end subroutine c_grnd_3

subroutine c_grnd_4(A)
  complex(8),dimension(:,:,:,:) :: A  ! complex(8) array of rank 4, filled with random numbers
  integer                       :: i1,i2,i3,i4
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              A(i1,i2,i3,i4) = dcmplx(mersenne(),mersenne())
           enddo
        enddo
     enddo
  enddo
end subroutine c_grnd_4

subroutine c_grnd_5(A)
  complex(8),dimension(:,:,:,:,:) :: A  ! complex(8) array of rank 5, filled with random numbers
  integer                         :: i1,i2,i3,i4,i5
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 A(i1,i2,i3,i4,i5) = dcmplx(mersenne(),mersenne())
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine c_grnd_5

subroutine c_grnd_6(A)
  complex(8),dimension(:,:,:,:,:,:) :: A  ! complex(8) array of rank 6, filled with random numbers
  integer                           :: i1,i2,i3,i4,i5,i6
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 do i6=1,size(A,6)
                    A(i1,i2,i3,i4,i5,i6) = dcmplx(mersenne(),mersenne())
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine c_grnd_6


subroutine c_grnd_7(A)
  complex(8),dimension(:,:,:,:,:,:,:) :: A  ! complex(8) array of rank 7, filled with random numbers
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  do i1=1,size(A,1)
     do i2=1,size(A,2)
        do i3=1,size(A,3)
           do i4=1,size(A,4)
              do i5=1,size(A,5)
                 do i6=1,size(A,6)
                    do i7=1,size(A,7)
                       A(i1,i2,i3,i4,i5,i6,i7) = dcmplx(mersenne(),mersenne())
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine c_grnd_7








!---------------------------------------------------------------!
! Integer random number generator:                              !
! return a random integer between in [l,h]                      !
! (boundaries l and h included)                                 !
!                                                               !
! A.Kuronen, 2009                                               !
!---------------------------------------------------------------!
integer function igrnd(l,h)
  integer,intent(in) :: l  ! lower bound, included
  integer,intent(in) :: h  ! upper bound, included
  real(8)            :: u,r
  u=grnd()
  r=(h-l+1)*u+l
  igrnd=int(r)
  return
end function igrnd


function dgrnd_uniform(a,b) result(c)
  real(8) :: a  ! lower bound
  real(8) :: b  ! upper bound
  real(8) :: c  ! random number in [a,b]
  c = a + (b-a)*mersenne()
end function dgrnd_uniform






!---------------------------------------------------------------!
! Random numbers with normal (Gaussian) distribution.           !
! Mean is 0 and standard deviation is 1                         !
! See W.H.Press et al., Numerical Recipes 1st ed., page 203     !
!                                                               !
! A.Kuronen, 2009                                               !
!---------------------------------------------------------------!  
function gaussrnd()
  real(8)       :: gaussrnd  ! standard normal deviate
  real(8)       :: fac,v1,v2,r
  real(8),save :: gset
  integer,save :: iset=0
  if (iset==0) then ! Create a new RN
     r=100.0
     do while (r>1.0)
        v1 = 2.0*grnd()-1.0
        v2 = 2.0*grnd()-1.0
        r = v1*v1+v2*v2
     end do
     fac = sqrt(-2.0*log(r)/r)
     gset = v1*fac
     gaussrnd = v2*fac
     iset = 1
  else ! Use the 2nd NR from the previous call
     gaussrnd = gset
     iset = 0
  endif
  return
end function gaussrnd


! Random Sample from normal (Gaussian) distribution
!
function normalrnd(mean,stdev) result(c)
  real(8)           :: r,theta
  real(8)           :: mean   ! mean
  real(8)           :: stdev  ! standard deviation, positive
  real(8)           :: c      ! normal deviate
  real(8),parameter :: pi=acos(-1d0)
  if(stdev <= 0d0) then
     Write(*,*) "NORMALRND: Standard Deviation must be +ve"
  else
     r     = (-2d0*log(mersenne()))**0.5d0
     theta = 2d0*pi*mersenne()
     c     = mean+stdev*r*sin(theta)
  end if
end function normalrnd




!  Random smaple from an exponential distribution 
!
function exponentialrnd(mean) result(c)
  real(8) :: mean  ! mean of the distribution, positive
  real(8) :: c     ! exponential deviate
  if (mean <= 0d0) then
     write(*,*) "EXPONENTIALRND: mean must be positive"
  else
     c=-mean*log(mersenne())
  end if
end function exponentialrnd



!  Return a random sample from a gamma distribution
!
recursive function gammarnd(shape,scale) result(ans)
  real(8) :: u,w,d,c,x,xsq,g,v
  real(8) :: shape  ! shape parameter, positive
  real(8) :: scale  ! scale parameter, positive
  real(8) :: ans    ! gamma deviate, mean shape*scale
  if (shape <= 0d0) then
     write(*,*) "GAMMARND: Shape parameter must be positive"
  end if
  if (scale <= 0d0) then
     write(*,*) "GAMMARND: Scale parameter must be positive"
  end if
  !    ## Implementation based on "A Simple Method for Generating Gamma Variables"
  !    ## by George Marsaglia and Wai Wan Tsang.  
  !    ## ACM Transactions on Mathematical Software
  !    ## Vol 26, No 3, September 2000, pages 363-372.
  if (shape >= 1d0) then
     d = shape - 1d0/3d0
     c = 1d0/(9d0*d)**0.5d0
     do while (.true.)
        x = normalrnd(0d0,1d0)
        v = 1.d0 + c*x
        do while (v <= 0d0) 
           x = normalrnd(0d0,1d0)
           v = 1d0 + c*x
        end do
        v = v*v*v
        u = mersenne()
        xsq = x*x
        if ((u < 1d0 -.0331d0*xsq*xsq) .or. (log(u) < 0.5d0*xsq + d*(1d0 - v + log(v))) )then
           ans=scale*d*v
           return 
        end if
     end do
  else
     g = gammarnd(shape+1d0,1d0)
     w = mersenne()
     ans=scale*g*(w)**(1d0/shape)
     return 
  end if
end function gammarnd



! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
function chi_squarernd(dof) result(ans)
  real(8) :: ans  ! chi-square deviate
  real(8) :: dof  ! degrees of freedom
  ans=gammarnd(0.5d0,2d0*dof)
end function chi_squarernd



! ## return a random sample from an inverse gamma random variable
!
function inverse_gammarnd(shape,scale) result(ans)
  real(8) :: shape  ! shape parameter, positive
  real(8) :: scale  ! scale parameter, positive
  real(8) :: ans    ! inverse gamma deviate
  ! If X is gamma(shape, scale) then
  ! 1/Y is inverse gamma(shape, 1/scale)
  ans= 1d0/gammarnd(shape,1d0/scale)
end function inverse_gammarnd



!## return a sample from a Weibull distribution
!
function weibullrnd(shape,scale) result(ans)
  real(8) :: temp
  real(8) :: shape  ! shape parameter, positive
  real(8) :: scale  ! scale parameter, positive
  real(8) :: ans    ! Weibull deviate
  if (shape <= 0d0) then
     write(*,*) "WEIBULLRND: Shape parameter must be positive"
  end if
  if (scale <= 0d0) then
     write(*,*) "WEIBULLRND: Scale parameter must be positive"
  end if
  ans= scale * (-log(mersenne()))**(1.d0/shape)
end function weibullrnd

!## return a random sample from a Cauchy distribution
!
function cauchyrnd(median,scale) result(ans)
  real(8) :: p
  real(8) :: ans     ! Cauchy deviate
  real(8) :: median  ! median of the distribution
  real(8) :: scale   ! scale parameter, positive
  if (scale <= 0d0) then
     write(*,*) "CAUCHYRND: Scale parameter must be positive"
  end if
  ans = median + scale*tan(PI*(mersenne() - 0.5))
end function cauchyrnd

!
!## return a random sample from a Student t distribution
!
function student_trnd(dof) result(ans)
  real(8) :: y1,y2
  real(8) :: ans  ! Student t deviate
  real(8) :: dof  ! degrees of freedom, positive
  if (dof <= 0.d0) then
     write(*,*) "STUDENT_TRND: Degrees of freedom must be positive"
  end if
  !
  ! ## See Seminumerical Algorithms by Knuth
  y1 = normalrnd(0d0,1d0)
  y2 = chi_squarernd(dof)
  ans= y1 / (y2 / dof)**0.5d0
  !
end function student_trnd

!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
function laplacernd(mean,scale)  result(ans)
  real(8) :: u
  real(8) :: ans    ! Laplace deviate
  real(8) :: mean   ! mean of the distribution
  real(8) :: scale  ! scale parameter, positive
  if (scale <= 0d0) then
     write(*,*) "LAPLACERND: Scale parameter must be positive"
  end if
  u = mersenne()
  if (u < 0.5d0) then
     ans = mean + scale*log(2d0*u) 
  else
     ans = mean - scale*log(2*(1-u))
  end if

end function laplacernd

!
! ## return a random sample from a log-normal distribution
!
function log_normalrnd(mu,sigma) result(ans)
  real(8) :: ans    ! log-normal deviate
  real(8) :: mu     ! mean of the logarithm of the variable
  real(8) :: sigma  ! standard deviation of the logarithm
  ans= exp(normalrnd(mu,sigma))
end function log_normalrnd

!
! ## return a random sample from a beta distribution
!
function betarnd(a,b) result(ans)
  real(8) :: u,v
  real(8) :: a    ! first shape parameter, positive
  real(8) :: b    ! second shape parameter, positive
  real(8) :: ans  ! beta deviate, in [0,1]
  if ((a <= 0d0) .or. (b <= 0d0)) then
     write(*,*) "BETARND: Beta parameters must be positive"
  end if
  !    ## There are more efficient methods for generating beta samples.
  !    ## However such methods are a little more efficient and much more complicated.
  !    ## For an explanation of why the following method works, see
  !    ## http://www.johndcook.com/distribution_chart.html#gamma_beta
  u = gammarnd(a,1d0)
  v = gammarnd(b,1d0)
  ans = u / (u + v)
end function betarnd



!---------------------------------------------------------------!
! State saving subroutines.                                     !
!                                                               !
!  Usage:  call mtsave( file_name, format_character )           !
!     or   call mtsave( unit_number, format_character )         !
!  where   format_character = 'u' or 'U' will save in           !
!          unformatted form, otherwise state information will   !
!          be written in formatted form.                        !
!---------------------------------------------------------------!
subroutine mtsavef( fname,forma )
  !NOTE: This subroutine APPENDS to the end of the file "fname".
  character(*),intent(in) :: fname  ! name of the file, the state is appended
  character,intent(in)    :: forma  ! 'u' or 'U': unformatted, otherwise formatted
  select case (forma)
  case('u','U')
     open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED',position='APPEND')
     write(10) mti
     write(10) mt       
  case default
     open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED',position='APPEND')
     write(10,*) mti
     write(10,*) mt       
  end select
  close(10)    
  return
end subroutine mtsavef

subroutine mtsaveu(unum,forma)
  integer,intent(in)    :: unum   ! open unit where the state is written
  character,intent(in)  :: forma  ! 'u' or 'U': unformatted, otherwise formatted
  select case (forma)
  case('u','U')
     write(unum) mti
     write(unum) mt       
  case default
     write(unum,*) mti
     write(unum,*) mt       
  end select
  return
end subroutine mtsaveu



!---------------------------------------------------------------!
! State getting subroutines.                                    !
!                                                               !
!  Usage:  call mtget( file_name, format_character )            !
!     or   call mtget( unit_number, format_character )          !
!  where   format_character = 'u' or 'U' will read in           !
!          unformatted form, otherwise state information will   !
!          be read in formatted form.                           !
!---------------------------------------------------------------!
subroutine mtgetf(fname,forma)
  character(*),intent(in) :: fname  ! name of the file, must exist
  character,intent(in)    :: forma  ! 'u' or 'U': unformatted, otherwise formatted
  select case (forma)
  case('u','U')
     open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
     read(10) mti
     read(10) mt
  case default
     open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
     read(10,*) mti
     read(10,*) mt
  end select
  close(10)
  return
end subroutine mtgetf

subroutine mtgetu(unum,forma)
  integer,intent(in)    :: unum   ! open unit from which the state is read
  character,intent(in)  :: forma  ! 'u' or 'U': unformatted, otherwise formatted
  select case (forma)
  case('u','U')
     read(unum) mti
     read(unum) mt
  case default
     read(unum,*) mti
     read(unum,*) mt
  end select
  return
end subroutine mtgetu
