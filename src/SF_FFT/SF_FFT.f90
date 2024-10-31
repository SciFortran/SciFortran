MODULE SF_FFT_FFTPACK
!SciFortran module for Fourier transform
  USE SF_INTEGRATE, only: simps
  USE SF_ARRAYS, only:linspace
  USE SF_CONSTANTS, only: pi2,xi,pi
  implicit none
  private

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms of time/frequency signal
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface FT_direct
     module procedure :: d_FT_direct
     module procedure :: c_FT_direct
  end interface FT_direct

  interface FT_inverse
     module procedure :: d_FT_inverse
     module procedure :: c_FT_inverse
  end interface FT_inverse

  interface FFT_signal
     module procedure :: d_FFT_signal
     module procedure :: c_FFT_signal
  end interface FFT_signal

  interface iFFT_signal
     module procedure :: d_iFFT_signal
     module procedure :: c_iFFT_signal
  end interface iFFT_signal

  interface tfft
     module procedure d_tfft,c_tfft
  end interface tfft

  interface itfft
     module procedure d_itfft,c_itfft
  end interface itfft

  public :: FT_direct
  public :: FT_inverse
  public :: FFT_signal
  public :: iFFT_signal
  public :: tfft
  public :: itfft


  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fft
     module procedure rfft_1d_forward,cfft_1d_forward
  end interface fft
  public :: fft
  public :: rfft_1d_forward
  public :: cfft_1d_forward

  interface ifft
     module procedure rfft_1d_backward,cfft_1d_backward
  end interface ifft
  public :: ifft
  public :: rfft_1d_backward
  public :: cfft_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fft2
     module procedure rfft_2d_forward,cfft_2d_forward
  end interface fft2
  public :: fft2
  public :: rfft_2d_forward
  public :: cfft_2d_forward

  interface ifft2
     module procedure rfft_2d_backward,cfft_2d_backward
  end interface ifft2
  public :: ifft2
  public :: rfft_2d_backward
  public :: cfft_2d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fftn
     module procedure rfft_nd_forward,cfft_nd_forward
  end interface fftn
  public :: fftn
  public :: rfft_nd_forward
  public :: cfft_nd_forward

  interface ifftn
     module procedure rfft_nd_backward,cfft_nd_backward
  end interface ifftn
  public :: ifftn
  public :: rfft_nd_backward
  public :: cfft_nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosft
     module procedure cost_1d_forward
  end interface cosft
  public :: cosft
  public :: cost_1d_forward

  interface icosft
     module procedure cost_1d_backward
  end interface icosft
  public :: icosft
  public :: cost_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosftn
     module procedure cost_Nd_forward
  end interface cosftn
  public :: cosftn
  public :: cost_nd_forward

  interface icosftn
     module procedure cost_Nd_backward
  end interface icosftn
  public :: icosftn
  public :: cost_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinft
     module procedure sint_1d_forward
  end interface sinft
  public :: sinft
  public :: sint_1d_forward

  interface isinft
     module procedure sint_1d_backward
  end interface isinft
  public :: isinft
  public :: sint_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinftn
     module procedure sint_Nd_forward
  end interface sinftn
  public :: sinftn
  public :: sint_nd_forward

  interface isinftn
     module procedure sint_Nd_backward
  end interface isinftn
  public :: isinftn
  public :: sint_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!





  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !HELPER FUNCTIONS:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fftshift
     module procedure rfft_1d_shift,cfft_1d_shift
  end interface fftshift
  public :: fftshift
  public :: rfft_1d_shift
  public :: cfft_1d_shift

  interface ifftshift
     module procedure rfft_1d_ishift,cfft_1d_ishift
  end interface ifftshift
  public :: ifftshift
  public :: rfft_1d_ishift
  public :: cfft_1d_ishift

  interface fftex
     module procedure rfft_1d_ex,cfft_1d_ex
  end interface fftex
  public :: fftex
  public :: rfft_1d_ex
  public :: cfft_1d_ex

  public :: fft_tmax
  public :: fft_fmax
  public :: fft_tarray
  public :: fft_farray
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  

contains







  !*********************************************************************
  !               TIME <==> FREQUENCY DOMAIN FFT:
  !*********************************************************************
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the simple FT of a function from time to frequency. 
  !+-------------------------------------------------------------------+
  function d_FT_direct(ft,t,w) result(fw)
    real(8),dimension(:),intent(in)        :: ft
    real(8),dimension(size(ft)),intent(in) :: t
    real(8),dimension(:),intent(in)        :: w
    real(8),dimension(size(w))             :: fw
    real(8)                                :: a,b
    integer                                :: i
    a = t(1);b = t(size(t))
    do i=1,size(w)
       fw(i)= simps(ft*exp(-xi*pi2*w(i)*t),a,b)
    enddo
  end function D_FT_Direct
  function c_FT_direct(ft,t,w) result(fw)
    complex(8),dimension(:),intent(in)     :: ft
    real(8),dimension(size(ft)),intent(in) :: t
    real(8),dimension(:),intent(in)        :: w
    complex(8),dimension(size(w))          :: fw
    real(8)                                :: a,b
    integer                                :: i
    a = t(1);b = t(size(t))
    do i=1,size(w)
       fw(i)= simps(ft*exp(-xi*pi2*w(i)*t),a,b)
    enddo
  end function C_FT_Direct


  function d_FT_inverse(fw,t,w) result(ft)
    real(8),dimension(:),intent(in)        :: fw
    real(8),dimension(:),intent(in)        :: t
    real(8),dimension(size(fw)),intent(in) :: w
    real(8),dimension(size(t))             :: ft
    real(8)                                :: a,b
    integer                                :: i
    a = w(1) ; b = w(size(w))
    do i=1,size(t)
       ft(i)= simps(fw*exp(xi*pi2*w*t(i)),a,b)
    enddo
  end function D_FT_Inverse
  function c_FT_inverse(fw,t,w) result(ft)
    complex(8),dimension(:),intent(in)     :: fw
    real(8),dimension(:),intent(in)        :: t
    real(8),dimension(size(fw)),intent(in) :: w
    complex(8),dimension(size(t))          :: ft
    real(8)                                :: a,b
    integer                                :: i
    a = w(1) ; b = w(size(w))
    do i=1,size(t)
       ft(i)= simps(fw*exp(xi*pi2*w*t(i)),a,b)
    enddo
  end function C_FT_Inverse



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a function from time to frequency. 
  !+-------------------------------------------------------------------+
  function d_FFT_signal(ft,dt) result(fw)
    real(8),dimension(:)        :: ft
    real(8)                     :: dt
    real(8),dimension(size(ft)) :: fw
    call tfft(ft)
    fw = ft*dt
  end function d_FFT_signal
  function c_FFT_signal(ft,dt) result(fw)
    complex(8),dimension(:)        :: ft
    real(8)                        :: dt
    complex(8),dimension(size(ft)) :: fw
    call tfft(ft)
    fw = ft*dt
  end function c_FFT_signal


  function d_iFFT_signal(fw,dt) result(ft)
    real(8),dimension(:)        :: fw
    real(8)                     :: dt
    real(8),dimension(size(fw)) :: ft
    call itfft(fw)
    ft=fw/dt
  end function d_iFFT_signal
  function c_iFFT_signal(fw,dt) result(ft)
    complex(8),dimension(:)        :: fw
    real(8)                        :: dt
    complex(8),dimension(size(fw)) :: ft
    call itfft(fw)
    ft=fw/dt
  end function c_iFFT_signal





  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a function from time to frequency. 
  !+-------------------------------------------------------------------+
  subroutine d_tfft(func_in,func_out)
    real(8),dimension(:)                         :: func_in
    real(8),dimension(size(func_in)),optional    :: func_out
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)*size(ftmp)
    else
       func_in  = fftshift(ftmp)*size(ftmp)
    endif
  end subroutine d_tfft
  !
  subroutine c_tfft(func_in,func_out)
    complex(8),dimension(:)                         :: func_in
    complex(8),dimension(size(func_in)),optional    :: func_out
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)*size(ftmp)
    else
       func_in  = fftshift(ftmp)*size(ftmp)
    endif
  end subroutine c_tfft

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a function from frequency to time.
  !+-------------------------------------------------------------------+
  subroutine d_itfft(func_in,func_out)
    real(8),dimension(:)                      :: func_in
    real(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)/size(ftmp)
    else
       func_in  = ifftshift(ftmp)/size(ftmp)
    endif
  end subroutine d_itfft
  !
  subroutine c_itfft(func_in,func_out)
    complex(8),dimension(:)                      :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)/size(ftmp)
    else
       func_in  = ifftshift(ftmp)/size(ftmp)
    endif
  end subroutine c_itfft







  !*********************************************************************
  !             FAST FOURIER TRANSFORM FUNCTIONS:
  !*********************************************************************  
  !                       FORWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM FFT using FFTPACK 5.1 routines.
  !out is :
  ! [y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]   if N is even
  ! [y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]   if N is odd
  !where:
  ! y(j) = sum[k=0..N-1] x[k] * exp(-xi*j*k* 2*pi/N), j = 0..N-1
  ! For index J*INC+1 where J=0,...,N-1 (that is, for the Jth
  !          element of the sequence):
  !               N-1
  ! C(J*INC+1) =  SUM C(K*INC+1)*EXP(-XI*J*K*2*PI/N)
  !               K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_forward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lenwrk = N
    lensav = N + int( log(dble(N))/log(2.d0) ) + 4
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_1d_forward: LENSAV not big enough"
    call rfft1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_1d_forward: LENR not big enough"
    case (2)
       stop "rfft_1d_forward: LENSAV not big enough"
    case (3)
       stop "rfft_1d_forward: LENWRK not big enough"
    case (20)
       stop "rfft_1d_forward: input error returned by lower level routine"
    end select
  end subroutine rfft_1d_forward
  !
  subroutine cfft_1d_forward(func)
    complex(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: N,lenwrk,lensav,lenc,inc,ier
    N      = size(func)
    lenwrk = 2*N
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenc   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_1d_forward: LENSAV not big enough"
    call cfft1f(N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_1d_forward: LENC not big enough"
    case (2)
       stop "cfft_1d_forward: LENSAV not big enough"
    case (3)
       stop "cfft_1d_forward: LENWRK not big enough"
    case (20)
       stop "cfft_1d_forward: input error returned by lower level routine"
    end select
  end subroutine cfft_1d_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 2-DIM FFT using FFTPACK 5.1 routines.
  ! The full complex transform of r(i,j) is given by:
  !                    L-1  M-1
  ! C(I,J) =   1/(L*M)*SUM  SUM  C(L1,M1)*EXP(-XI*2*PI*(I*L1/L + J*M1/M))
  !                    L1=0 M1=0
  !
  !+-------------------------------------------------------------------+
  subroutine rfft_2d_forward(func)
    real(8),dimension(:,:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,M,lenwrk,lensav,ldim,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = M*(L+1)
    lensav = L+3*M  + int(log(dble(L))/log(2.d0))+2*int(log(dble(M))/log(2.d0)) + 12
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "rfft_2d_forward: LENSAV not big enough"
    if(ier==20)stop "rfft_2d_forward: input error returned by lower level routine"
    call rfft2f(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (6)
       stop "rfft_2d_forward: LDIM is less than 2*INT((L+1)/2)"
    case (2)
       stop "rfft_2d_forward: LENSAV not big enough"
    case (3)
       stop "rfft_2d_forward: LENWRK not big enough"
    case (20)
       stop "rfft_2d_forward: input error returned by lower level routine"
    end select
  end subroutine rfft_2d_forward
  !
  subroutine cfft_2d_forward(func)
    complex(8),dimension(:,:),intent(inout) :: func
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,M,lenwrk,lensav,ldim,inc,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = 2*L*M
    lensav = 2*(L+M) + int(log(dble(L))/log(2.d0))+int(log(dble(M))/log(2.d0)) + 8
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "cfft_2d_forward: LENSAV not big enough"
    if(ier==20)stop "cfft_2d_forward: input error returned by lower level routine"
    call cfft2f(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (5)
       stop "cfft_2d_forward: L > LDIM"
    case (2)
       stop "cfft_2d_forward: LENSAV not big enough"
    case (3)
       stop "cfft_2d_forward: LENWRK not big enough"
    case (20)
       stop "cfft_2d_forward: input error returned by lower level routine"
    end select
  end subroutine cfft_2d_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward N-DIM FFT using FFTPACK 5.1 routines.
  ! For index L*JUMP + J*INC +1 where J=0,...,N-1 and 
  ! L=0,...,LOT-1, (that is, for the Jth element of the Lth sequence),
  !                       N-1
  !  C(L*JUMP+J*INC+1) = SUM C(L*JUMP+K*INC+1)*EXP(-XI*J*K*2*PI/N)
  !                       K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_nd_forward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "rfft_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = N*lot
    lensav = N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_Nd_forward: LENSAV not big enough"
    call rfftmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_Nd_forward: LENR   not big enough"
    case (2)
       stop "rfft_Nd_forward: LENSAV not big enough"
    case (3)
       stop "rfft_Nd_forward: LENWRK not big enough"
    case (4)
       stop "rfft_Nd_forward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine rfft_nd_forward
  !
  subroutine cfft_nd_forward(func,n,lot)
    complex(8),dimension(:),intent(inout) :: func
    integer,intent(in)                    :: n,lot
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,lenwrk,lensav,lenc,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cfft_Nd_forward: incommensurate values of parameters." 
    lenc   = N*lot
    lenwrk = 2*N*lot
    lensav = 2*N + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_Nd_forward: LENSAV not big enough"
    call cfftmf(lot,jump,N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_Nd_forward: LENR   not big enough"
    case (2)
       stop "cfft_Nd_forward: LENSAV not big enough"
    case (3)
       stop "cfft_Nd_forward: LENWRK not big enough"
    case (4)
       stop "cfft_Nd_forward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine cfft_nd_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_1d_forward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = N-1
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cost1i(N,wsave,lensav,ier)
    if(ier==2)stop "cost_1d_forward: LENSAV not big enough"
    if(ier==20)stop "cost_1d_forward: error returned by lower level routine"
    call cost1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_1d_forward: LENR not big enough"
    case (2)
       stop "cost_1d_forward: LENSAV not big enough"
    case (3)
       stop "cost_1d_forward: LENWRK not big enough"
    case (20)
       stop "cost_1d_forward: input error returned by lower level routine"
    end select
  end subroutine cost_1d_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_1d_forward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = N/2 + N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = 2*(N+1)
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call sint1i(N,wsave,lensav,ier)
    if(ier==2)stop "sint_1d_forward: LENSAV not big enough"
    if(ier==20)stop "sint_1d_forward: error returned by lower level routine"
    call sint1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_1d_forward: LENR not big enough"
    case (2)
       stop "sint_1d_forward: LENSAV not big enough"
    case (3)
       stop "sint_1d_forward: LENWRK not big enough"
    case (20)
       stop "sint_1d_forward: input error returned by lower level routine"
    end select
  end subroutine sint_1d_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward M-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_nd_forward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*(N+1)
    lensav = 2*N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call costmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_forward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_forward: error returned by lower level routine"
    call costmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_forward: LENR not big enough"
    case (2)
       stop "cost_Nd_forward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_forward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_forward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_forward: input error returned by lower level routine"
    end select
  end subroutine cost_nd_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward M-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_nd_forward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*2*(N+2)
    lensav = N/2 + N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call sintmi(N,wsave,lensav,ier)
    if(ier==2)stop "sint_Nd_forward: LENSAV not big enough"
    if(ier==20)stop "sint_Nd_forward: error returned by lower level routine"
    call sintmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_Nd_forward: LENR not big enough"
    case (2)
       stop "sint_Nd_forward: LENSAV not big enough"
    case (3)
       stop "sint_Nd_forward: LENWRK not big enough"
    case (4)
       stop "sint_Nd_forward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "sint_Nd_forward: input error returned by lower level routine"
    end select
  end subroutine sint_nd_forward



  !*********************************************************************
  !                       BACKWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate 1-DIM backward FFT using FFTPACK 5.1 routines.
  ! The returned real/complex array contains y(0), y(1),..., y(n-1) where
  ! y(j) = sum_{k=-N/2,...,N/2-1}(x(k) * exp(2*pi*XI*j*k/N))
  ! For index J*INC+1 where J=0,...,N-1,
  !              N-1
  ! C(J*INC+1) = SUM C(K*INC+1)*EXP(XI*J*K*2*PI/N)
  !              K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_backward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lenwrk = N
    lensav = N + int( log(dble(N))/log(2.d0) ) + 4
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    call rfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_1d_backward: LENSAV not big enough"
    allocate(work(lenwrk))
    call rfft1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_1d_backward: LENR not big enough"
    case (2)
       stop "rfft_1d_backward: LENSAV not big enough"
    case (3)
       stop "rfft_1d_backward: LENWRK not big enough"
    case (20)
       stop "rfft_1d_backward: input error returned by lower level routine"
    end select
  end subroutine rfft_1d_backward
  !
  subroutine cfft_1d_backward(func)
    complex(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: N,lenwrk,lensav,lenc,inc,ier
    N      = size(func)
    lenwrk = 2*N
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenc   = N
    inc    = 1
    allocate(wsave(lensav))
    call cfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_1d_backward: LENSAV not big enough"
    allocate(work(lenwrk))
    call cfft1b(N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_1d_backward: LENC not big enough"
    case (2)
       stop "cfft_1d_backward: LENSAV not big enough"
    case (3)
       stop "cfft_1d_backward: LENWRK not big enough"
    case (20)
       stop "cfft_1d_backward: input error returned by lower level routine"
    end select
  end subroutine cfft_1d_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 2-DIM FFT using FFTPACK 5.1 routines.
  ! For purposes of exposition, assume the index ranges of array C are defined by
  !  C(0:L-1,0:M-1).
  ! For I=0,...,L-1 and J=0,...,M-1, the C(I,J)'s are given in the traditional 
  ! aliased form by
  !           L-1  M-1
  !  C(I,J) = SUM  SUM  C(L1,M1)*EXP(XI*2*PI*(I*L1/L + J*M1/M))
  !           L1=0 M1=0
  !+-------------------------------------------------------------------+
  subroutine rfft_2d_backward(func)
    real(8),dimension(:,:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,M,lenwrk,lensav,ldim,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = M*(L+1)
    lensav = L+3*M  + int(log(dble(L))/log(2.d0))+2*int(log(dble(M))/log(2.d0)) + 12
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "rfft_2d_backward: LENSAV not big enough"
    if(ier==20)stop "rfft_2d_backward: input error returned by lower level routine"
    call rfft2b(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (6)
       stop "rfft_2d_backward: LDIM is less than 2*INT((L+1)/2)"
    case (2)
       stop "rfft_2d_backward: LENSAV not big enough"
    case (3)
       stop "rfft_2d_backward: LENWRK not big enough"
    case (20)
       stop "rfft_2d_backward: input error returned by lower level routine"
    end select
  end subroutine rfft_2d_backward
  !
  subroutine cfft_2d_backward(func)
    complex(8),dimension(:,:),intent(inout) :: func
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,M,lenwrk,lensav,ldim,inc,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = 2*L*M
    lensav = 2*(L+M) + int(log(dble(L))/log(2.d0))+int(log(dble(M))/log(2.d0)) + 8
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "cfft_2d_backward: LENSAV not big enough"
    if(ier==20)stop "cfft_2d_backward: input error returned by lower level routine"
    call cfft2b(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (5)
       stop "cfft_2d_backward: L > LDIM"
    case (2)
       stop "cfft_2d_backward: LENSAV not big enough"
    case (3)
       stop "cfft_2d_backward: LENWRK not big enough"
    case (20)
       stop "cfft_2d_backward: input error returned by lower level routine"
    end select
  end subroutine cfft_2d_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward N-DIM FFT using FFTPACK 5.1 routines.
  ! For index L*JUMP+J*INC+1 where J=0,...,N-1 and 
  ! L=0,...,LOT-1, (that is, for the Jth element of the Lth sequence),
  !                      N-1
  !  C(L*JUMP+J*INC+1) = SUM C(L*JUMP+K*INC+1)*EXP(XI*J*K*2*PI/N)
  !                      K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_nd_backward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "rfft_Nd_backward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = N*lot
    lensav = N  + int(log(dble(N))/log(2.d0)) + 4
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_Nd_backward: LENSAV not big enough"
    jump   = N
    inc    = 1
    call rfftmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_Nd_backward: LENR   not big enough"
    case (2)
       stop "rfft_Nd_backward: LENSAV not big enough"
    case (3)
       stop "rfft_Nd_backward: LENWRK not big enough"
    case (4)
       stop "rfft_Nd_backward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine rfft_nd_backward
  !
  subroutine cfft_nd_backward(func,n,lot)
    complex(8),dimension(:),intent(inout) :: func
    integer,intent(in)                    :: n,lot
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,lenwrk,lensav,lenc,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cfft_Nd_backward: incommensurate values of parameters." 
    !
    lenc   = N*lot
    lenwrk = 2*N*lot
    lensav = 2*N + int(log(dble(N))/log(2.d0)) + 4
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_Nd_backward: LENSAV not big enough"
    !
    jump   = N
    inc    = 1
    call cfftmb(lot,jump,N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_Nd_backward: LENR   not big enough"
    case (2)
       stop "cfft_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cfft_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cfft_Nd_backward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine cfft_nd_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 1-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_1d_backward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = N-1
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cost1i(N,wsave,lensav,ier)
    if(ier==2)stop "cost_1d_backward: LENSAV not big enough"
    if(ier==20)stop "cost_1d_backward: error returned by lower level routine"
    call cost1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_1d_backward: LENR not big enough"
    case (2)
       stop "cost_1d_backward: LENSAV not big enough"
    case (3)
       stop "cost_1d_backward: LENWRK not big enough"
    case (20)
       stop "cost_1d_backward: input error returned by lower level routine"
    end select
  end subroutine cost_1d_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 1-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_1d_backward(func)
    real(8),dimension(:),intent(inout) :: func
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = N/2 + N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = 2*(N+1)
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call sint1i(N,wsave,lensav,ier)
    if(ier==2)stop "sint_1d_backward: LENSAV not big enough"
    if(ier==20)stop "sint_1d_backward: error returned by lower level routine"
    call sint1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_1d_backward: LENR not big enough"
    case (2)
       stop "sint_1d_backward: LENSAV not big enough"
    case (3)
       stop "sint_1d_backward: LENWRK not big enough"
    case (20)
       stop "sint_1d_backward: input error returned by lower level routine"
    end select
  end subroutine sint_1d_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward M-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_nd_backward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_backward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*(N+1)
    lensav = 2*N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call costmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_backward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_backward: error returned by lower level routine"
    call costmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_backward: LENR not big enough"
    case (2)
       stop "cost_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_backward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_backward: input error returned by lower level routine"
    end select
  end subroutine cost_nd_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward M-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_nd_backward(func,n,lot)
    real(8),dimension(:),intent(inout) :: func
    integer,intent(in)                 :: N,lot
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_backward: incommensurate values of parameters."
    !
    lenr   = N*lot
    lenwrk = 2*lot*(N+2)
    lensav = N/2 + N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call sintmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_backward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_backward: error returned by lower level routine"
    call sintmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_backward: LENR not big enough"
    case (2)
       stop "cost_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_backward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_backward: input error returned by lower level routine"
    end select
  end subroutine sint_nd_backward














  !*********************************************************************
  !                           HELPER FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Shift the zero-frequency to the center of the 1-DIM array
  ! output of the forward-FFT is:
  ! [y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]   if N is even
  ! [y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]   if N is odd
  ! using *shift produces:
  ! [y(-N/2+1),   ..., y(-1), y(0), y(1), ..., y(N/2)]       if N is even
  ! [y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]   if N is even
  !+-------------------------------------------------------------------+
  function rfft_1d_shift(fin) result(fout)
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_shift
  !
  function cfft_1d_shift(fin) result(fout)
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_shift


  !+-------------------------------------------------------------------+
  !PURPOSE  : The inverse of _shift
  !+-------------------------------------------------------------------+
  function rfft_1d_ishift(fin) result(fout)
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_ishift
  !
  function cfft_1d_ishift(fin) result(fout)
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_ishift


  !+-------------------------------------------------------------------+
  !PURPOSE  : Alternate sign of the output of a FFT:
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_ex(func)
    real(8),dimension(:)    :: func
    real(8)                 :: ex
    integer                 :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine rfft_1d_ex
  !
  subroutine cfft_1d_ex(func)
    complex(8),dimension(:) :: func
    real(8)                 :: ex
    integer                 :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine cfft_1d_ex






  ! TIME-FREQUENCY SIGNALS
  function fft_tmax(L,dt)
    integer :: L
    real(8) :: dt
    real(8) :: fft_tmax
    fft_tmax = L*dt/2
  end function fft_tmax

  function fft_fmax(L,dt)
    integer :: L
    real(8) :: dt
    real(8) :: fft_fmax
    fft_fmax = pi/dt
  end function fft_fmax




  function fft_tarray(L,dt) result(time)
    real(8)              :: dt
    integer              :: L
    real(8)              :: tmax
    real(8),dimension(L) :: time
    tmax = fft_tmax(L,dt)
    time = linspace(-tmax,tmax,L)
  end function fft_tarray

  function fft_farray(L,dt,df) result(freq)
    integer              :: L
    real(8)              :: dt
    real(8),optional     :: df
    real(8)              :: wmax,dw
    real(8),dimension(L) :: freq
    wmax = fft_fmax(L,dt)
    freq = linspace(-wmax,wmax,L,iend=.false.,mesh=dw)
    if(present(df))df=dw
  end function fft_farray


END MODULE SF_FFT_FFTPACK
