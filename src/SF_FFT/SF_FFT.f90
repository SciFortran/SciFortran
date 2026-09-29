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
!This function evaluates the direct Fourier transform of a discretized function 
!from time to frequency domain. Takes as input a discretized function :math:`ft_{i}` on a 
!set of time points :math:`t_{i}` and a set of frequencies :math:`w_{j}` and returns 
!a discretized function :math:`fw_{j}` = :f:func_inline:`simps` 
!:math:`(ft_{i} \cdot e^{-i 2\pi w_{j}t_{i}}, t_{1}, t_{N})` , where :code:`N=size(t)`
     module procedure :: d_FT_direct
     module procedure :: c_FT_direct
  end interface FT_direct

  interface FT_inverse
!This function evaluates the inverse Fourier transform of a discretized function 
!from frequency to time domain. Takes as input a discretized function :math:`fw_{i}` on a 
!set of time points :math:`w_{i}` and a set of frequencies :math:`t_{j}` and returns 
!a discretized function :math:`ft_{j}` = :f:func_inline:`simps` 
!:math:`(fw_{i} \cdot e^{i 2\pi t_{j}w_{i}}, w_{1}, w_{N})` , where :code:`N=size(w)`
     module procedure :: d_FT_inverse
     module procedure :: c_FT_inverse
  end interface FT_inverse

  interface FFT_signal
!This function evaluates the fast Fourier transform of a discretized function :f:var:`ft` .
!It returns :f:var:`fw` = :f:var:`dt` · 
!:f:func_inline:`tfft` ( :f:var:`ft` )
!
     module procedure :: d_FFT_signal
     module procedure :: c_FFT_signal
  end interface FFT_signal

  interface iFFT_signal
!This function evaluates the inverse fast Fourier transform of a discretized function 
!:f:var:`fw` . It returns :f:var:`ft` = :f:func_inline:`itfft` ( :f:var:`fw` ) / :f:var:`ft`
!
     module procedure :: d_iFFT_signal
     module procedure :: c_iFFT_signal
  end interface iFFT_signal

  interface tfft
  !This subroutine evaluates the fast Fourier transform of a real or complex array
  !:f:var:`func_in` of length :code:`N`, arranged with the zero of the domain at the
  !center of the array. It performs the sequence
  !
  !* :f:func_inline:`ifftshift` of the input
  !* :f:func_inline:`fft` (forward transform)
  !* :f:func_inline:`fftshift` of the result
  !* multiplication by :code:`N`
  !
  !If the optional argument :f:var:`func_out` is present it contains the result and 
  !:f:var:`func_in` is left untouched, otherwise :f:var:`func_in` is overwritten. 
  !In the real case the imaginary part of the transform is discarded.
  !
     module procedure d_tfft,c_tfft
  end interface tfft

  interface itfft
  !This subroutine evaluates the inverse fast Fourier transform of a real or complex 
  !array :f:var:`func_in` of length :code:`N`, and is the counterpart of :f:func_inline:`tfft`. 
  !It performs the sequence
  !
  !* :f:func_inline:`ifft` (backward transform)
  !* :f:func_inline:`fftex` (alternate the sign of the elements, shifting the output by half a period)
  !* :f:func_inline:`ifftshift` of the result
  !* division by :code:`N`
  !
  !If the optional argument :f:var:`func_out` is present it contains the result and 
  !:f:var:`func_in` is left untouched, otherwise :f:var:`func_in` is overwritten. 
  !In the real case the imaginary part of the transform is discarded.
  !
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
  !This subroutine evaluates the forward 1-dimensional FFT of an array of lenght :code:`N` 
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + :f:func_inline:`rfft1f`
  !for the real case and :f:func_inline:`cfft1i` + :f:func_inline:`cfft1f` for the 
  !complex case. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  ! 
  !The subroutine modifies the input array to contain
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !where
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`  
  !
  !for :math:`j \in [0,N-1]`
  !
     module procedure rfft_1d_forward,cfft_1d_forward
  end interface fft
  public :: fft
  public :: rfft_1d_forward
  public :: cfft_1d_forward

  interface ifft
  !This subroutine evaluates the backward 1-dimensional FFT of an array of lenght :code:`N` 
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + :f:func_inline:`rfft1b`
  !for the real case and :f:func_inline:`cfft1i` + :f:func_inline:`cfft1b` for the 
  !complex case. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  ! 
  !The subroutine modifies the input array to contain
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !where
  !
  !:math:`y(j) = \sum_{k=-N/2}^{N/2-1} x(k) \cdot a^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`
  !
     module procedure rfft_1d_backward,cfft_1d_backward
  end interface ifft
  public :: ifft
  public :: rfft_1d_backward
  public :: cfft_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fft2
  !This subroutine evaluates the forward 2-dimensional FFT of an array of size :code:`(L,M)`
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + :f:func_inline:`rfft2f`
  !for the real case and :f:func_inline:`cfft2i` + :f:func_inline:`cfft2f` for the
  !complex case. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(N))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12` (real)
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(N))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8` (complex)
  !* :code:`lenwrk` = :code:`M*(L+1)` (real), :code:`2*L*M` (complex)
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The real subroutine returns the
  !FFTPACK packed form of the transform of a real array.
  !
     module procedure rfft_2d_forward,cfft_2d_forward
  end interface fft2
  public :: fft2
  public :: rfft_2d_forward
  public :: cfft_2d_forward

  interface ifft2
  !This subroutine evaluates the backward 2-dimensional FFT of an array of size :code:`(L,M)`
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + :f:func_inline:`rfft2b`
  !for the real case and :f:func_inline:`cfft2i` + :f:func_inline:`cfft2b` for the
  !complex case. These are called with the same parameters as :f:func_inline:`fft2`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`.
  !
     module procedure rfft_2d_backward,cfft_2d_backward
  end interface ifft2
  public :: ifft2
  public :: rfft_2d_backward
  public :: cfft_2d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fftn
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmf`
  !for the real case and :f:func_inline:`cfftmi` + :f:func_inline:`cfftmf` for the
  !complex case. These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4` (real), :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4` (complex)
  !* :code:`lenwrk` = :code:`N*lot` (real), :code:`2*N*lot` (complex)
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine modifies the input array to contain
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure rfft_nd_forward,cfft_nd_forward
  end interface fftn
  public :: fftn
  public :: rfft_nd_forward
  public :: cfft_nd_forward

  interface ifftn
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmb`
  !for the real case and :f:func_inline:`cfftmi` + :f:func_inline:`cfftmb` for the
  !complex case. These are called with the same parameters as :f:func_inline:`fftn`.
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine modifies the input array to contain
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
     module procedure rfft_nd_backward,cfft_nd_backward
  end interface ifftn
  public :: ifftn
  public :: rfft_nd_backward
  public :: cfft_nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosft
  !This subroutine evaluates the forward 1-dimensional cosine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + :f:func_inline:`cost1f`.
  !These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The subroutine modifies the input array in place. The transform is the inverse of
  !:f:func_inline:`icosft`, i.e. the normalization is carried by the forward transform.
  !
     module procedure cost_1d_forward
  end interface cosft
  public :: cosft
  public :: cost_1d_forward

  interface icosft
  !This subroutine evaluates the backward 1-dimensional cosine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + :f:func_inline:`cost1b`.
  !These are called with the same parameters as :f:func_inline:`cosft`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
     module procedure cost_1d_backward
  end interface icosft
  public :: icosft
  public :: cost_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosftn
  !This subroutine evaluates :code:`lot` forward cosine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure cost_Nd_forward
  end interface cosftn
  public :: cosftn
  public :: cost_nd_forward

  interface icosftn
  !This subroutine evaluates :code:`lot` backward cosine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmb`.
  !These are called with the same parameters as :f:func_inline:`cosftn`.
  !
     module procedure cost_Nd_backward
  end interface icosftn
  public :: icosftn
  public :: cost_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinft
  !This subroutine evaluates the forward 1-dimensional sine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + :f:func_inline:`sint1f`.
  !These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The subroutine modifies the input array in place. The transform is the inverse of
  !:f:func_inline:`isinft`, i.e. the normalization is carried by the forward transform.
  !
     module procedure sint_1d_forward
  end interface sinft
  public :: sinft
  public :: sint_1d_forward

  interface isinft
  !This subroutine evaluates the backward 1-dimensional sine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + :f:func_inline:`sint1b`.
  !These are called with the same parameters as :f:func_inline:`sinft`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
     module procedure sint_1d_backward
  end interface isinft
  public :: isinft
  public :: sint_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinftn
  !This subroutine evaluates :code:`lot` forward sine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure sint_Nd_forward
  end interface sinftn
  public :: sinftn
  public :: sint_nd_forward

  interface isinftn
  !This subroutine evaluates :code:`lot` backward sine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmb`.
  !These are called with the same parameters as :f:func_inline:`sinftn`.
  !
     module procedure sint_Nd_backward
  end interface isinftn
  public :: isinftn
  public :: sint_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!





  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !HELPER FUNCTIONS:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fftshift
  !This function shifts the zero-frequency component of the output of :f:func_inline:`fft`
  !to the center of the array. The input
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
     module procedure rfft_1d_shift,cfft_1d_shift
  end interface fftshift
  public :: fftshift
  public :: rfft_1d_shift
  public :: cfft_1d_shift

  interface ifftshift
  !This function is the inverse of :f:func_inline:`fftshift`: it moves the
  !zero-frequency component from the center of the array back to the first position.
  !
     module procedure rfft_1d_ishift,cfft_1d_ishift
  end interface ifftshift
  public :: ifftshift
  public :: rfft_1d_ishift
  public :: cfft_1d_ishift

  interface fftex
  !This subroutine alternates the sign of the elements of an array,
  !:math:`f(i) \to (-1)^{i-1} f(i)`. It is used by :f:func_inline:`itfft` to shift the output of
  !:f:func_inline:`ifft` by half a period.
  !
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
  function d_FT_direct(ft,t,w) result(fw)
    real(8),dimension(:),intent(in)        :: ft !Dicretized function to transform (time domain)
    real(8),dimension(size(ft)),intent(in) :: t  !Discretized time points
    real(8),dimension(:),intent(in)        :: w  !Discretized frequency points
    real(8),dimension(size(w))             :: fw !Discretized Fourier-transform function (frequency domain)
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
    real(8),dimension(:),intent(in)        :: fw !Dicretized function to transform (frequency domain)
    real(8),dimension(:),intent(in)        :: t  !Discretized time points
    real(8),dimension(size(fw)),intent(in) :: w  !Discretized frequency points
    real(8),dimension(size(t))             :: ft !Discretized Fourier-transform function (time domain)
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



  function d_FFT_signal(ft,dt) result(fw)
    real(8),dimension(:)        :: ft !Time-domain function
    real(8)                     :: dt !Time step
    real(8),dimension(size(ft)) :: fw !Frequency-domain function
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
    real(8),dimension(:)        :: fw !Frequency-domain function
    real(8)                     :: dt !Time step
    real(8),dimension(size(fw)) :: ft !Time-domain function
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

  subroutine rfft_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + 
  !:f:func_inline:`rfft1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized. Since the transform of a real sequence is 
  !Hermitian, the subroutine returns the FFTPACK packed form of
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !that is, the array contains
  !
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y(N/2-1), Im y(N/2-1), Re y(N/2)]`   if `N` is even
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y((N-1)/2), Im y((N-1)/2)]`           if `N` is odd
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func      !Input/output array of size :code:`N`
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

  subroutine cfft_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional FFT of a complex array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft1i` + 
  !:f:func_inline:`cfft1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N`
  !* :code:`lenc`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized. The subroutine modifies the input array to contain
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, i.e. the array is returned in the aliased order
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine rfft_2d_forward(func)
  !This subroutine evaluates the forward 2-dimensional FFT of a real array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + 
  !:f:func_inline:`rfft2f`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(L))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12`
  !* :code:`lenwrk` = :code:`M*(L+1)`
  !
  !The subroutine computes
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The transform of a real array is 
  !Hermitian, so the array is overwritten with the FFTPACK packed form of :math:`y`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine cfft_2d_forward(func)
  !This subroutine evaluates the forward 2-dimensional FFT of a complex array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft2i` + 
  !:f:func_inline:`cfft2f`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(L))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8`
  !* :code:`lenwrk` = :code:`2*L*M`
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine rfft_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transforms are not normalized. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`,
  !the subroutine computes
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. Since the transform of a real sequence is Hermitian, each 
  !sequence is overwritten with the FFTPACK packed form of :math:`y`, as in 
  !:f:func_inline:`rfft_1d_forward`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the complex sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`cfftmi` + :f:func_inline:`cfftmf`.
  !These are called with the parameters
  !
  !* :code:`lenc` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transforms are not normalized. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`,
  !the subroutine modifies the input array to contain
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. Each sequence is returned in the aliased order described 
  !in :f:func_inline:`cfft_1d_forward`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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


  subroutine cost_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional cosine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + 
  !:f:func_inline:`cost1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`cost_1d_backward`, i.e. the 
  !normalization is carried by the forward transform. The subroutine modifies the 
  !input array to contain the coefficients :math:`y(k)` such that the backward 
  !transform :f:func_inline:`cost_1d_backward` reconstructs
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine sint_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional sine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + 
  !:f:func_inline:`sint1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`sint_1d_backward`, i.e. the 
  !normalization is carried by the forward transform. The subroutine modifies the 
  !input array to contain the coefficients :math:`y(k)` such that the backward 
  !transform :f:func_inline:`sint_1d_backward` reconstructs
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine cost_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional cosine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`cost_1d_forward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`cost_nd_backward` 
  !and the normalization is carried by the forward transform.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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


  subroutine sint_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional sine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`sint_1d_forward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`sint_nd_backward` 
  !and the normalization is carried by the forward transform.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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


  subroutine rfft_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + 
  !:f:func_inline:`rfft1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized, so applying it after :f:func_inline:`rfft_1d_forward`
  !returns the original array multiplied by :code:`N`. The input array is expected in the 
  !FFTPACK packed form produced by :f:func_inline:`rfft_1d_forward`, that is
  !
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y(N/2-1), Im y(N/2-1), Re y(N/2)]`   if `N` is even
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y((N-1)/2), Im y((N-1)/2)]`           if `N` is odd
  !
  !and the subroutine modifies the input array to contain the real sequence
  !
  !:math:`x(j) = \sum_{k=0}^{N-1} y(k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, where :math:`y(k)` is the Hermitian sequence encoded by the packed input.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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
  !This subroutine evaluates the backward 1-dimensional FFT of a complex array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft1i` + 
  !:f:func_inline:`cfft1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N`
  !* :code:`lenc`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized, so applying it after :f:func_inline:`cfft_1d_forward`
  !returns the original array multiplied by :code:`N`. The input array is expected in the 
  !aliased order produced by :f:func_inline:`cfft_1d_forward`, that is
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !and the subroutine modifies the input array to contain
  !
  !:math:`x(j) = \sum_{k=0}^{N-1} y(k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine rfft_2d_backward(func)
  !This subroutine evaluates the backward 2-dimensional FFT of a real array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + 
  !:f:func_inline:`rfft2b`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(L))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12`
  !* :code:`lenwrk` = :code:`M*(L+1)`
  !
  !The input array is expected in the FFTPACK packed form produced by 
  !:f:func_inline:`rfft_2d_forward`, and the subroutine modifies it to contain the real array
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`, where :math:`y(l,m)` is the Hermitian 
  !array encoded by the packed input. The transform carries no normalization factor, 
  !because the :math:`1/(LM)` factor is carried by :f:func_inline:`rfft_2d_forward`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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
  !This subroutine evaluates the backward 2-dimensional FFT of a complex array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft2i` + 
  !:f:func_inline:`cfft2b`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(L))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8`
  !* :code:`lenwrk` = :code:`2*L*M`
  !
  !The input array is expected in the form produced by :f:func_inline:`cfft_2d_forward`,
  !and the subroutine modifies it to contain
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The transform carries no normalization 
  !factor, because the :math:`1/(LM)` factor is carried by :f:func_inline:`cfft_2d_forward`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine rfft_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !Each sequence is expected in the FFTPACK packed form produced by 
  !:f:func_inline:`rfft_nd_forward`. The transforms are not normalized, so applying 
  !them after :f:func_inline:`rfft_nd_forward` returns the original array multiplied 
  !by :code:`N`. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine 
  !modifies the input array to contain the real sequence
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, where :math:`y(l N + k)` is the Hermitian sequence encoded 
  !by the packed input.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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



  subroutine cfft_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the complex sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`cfftmi` + :f:func_inline:`cfftmb`.
  !These are called with the parameters
  !
  !* :code:`lenc` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !Each sequence is expected in the aliased order produced by 
  !:f:func_inline:`cfft_nd_forward`. The transforms are not normalized, so applying 
  !them after :f:func_inline:`cfft_nd_forward` returns the original array multiplied 
  !by :code:`N`. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine 
  !modifies the input array to contain
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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



  subroutine cost_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional cosine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + 
  !:f:func_inline:`cost1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`cost_1d_forward`, i.e. the 
  !normalization is carried by the forward transform, so no scaling is applied here. 
  !The input array contains the coefficients :math:`y(k)` produced by 
  !:f:func_inline:`cost_1d_forward`, and the subroutine modifies it to contain
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine sint_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional sine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + 
  !:f:func_inline:`sint1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`sint_1d_forward`, i.e. the 
  !normalization is carried by the forward transform, so no scaling is applied here. 
  !The input array contains the coefficients :math:`y(k)` produced by 
  !:f:func_inline:`sint_1d_forward`, and the subroutine modifies it to contain
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
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


  subroutine cost_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional cosine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`cost_1d_backward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`cost_nd_forward` 
  !and the normalization is carried by the forward transform, so no scaling is applied here.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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


  subroutine sint_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional sine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`sint_1d_backward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`sint_nd_forward` 
  !and the normalization is carried by the forward transform, so no scaling is applied here.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
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

  function rfft_1d_shift(fin) result(fout)
  !This function shifts the zero-frequency component of a real array :f:var:`fin`
  !of length :code:`N` to the center of the array, and returns the result as an array 
  !of the same size. It is a circular rotation by :code:`floor((N+1)/2)` elements, 
  !so the input
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
  !where :code:`y(N/2)` and :code:`y(-N/2)` are the same element for even :code:`N`.
  !
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_shift
  !
  function cfft_1d_shift(fin) result(fout)
  !This function shifts the zero-frequency component of a complex array :f:var:`fin`
  !of length :code:`N` to the center of the array, and returns the result as an array 
  !of the same size. It acts as :f:func_inline:`rfft_1d_shift` on complex data.
  !
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_shift

  function rfft_1d_ishift(fin) result(fout)
  !This function is the inverse of :f:func_inline:`rfft_1d_shift`. It moves the
  !zero-frequency component of a real array :f:var:`fin` of length :code:`N` from the 
  !center of the array back to the first position, and returns the result as an array 
  !of the same size. It is a circular rotation by :code:`N-floor((N+1)/2)` elements, 
  !so the input
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(0), y(1), ..., y(N/2-1), y(-N/2), ..., y(-1)]`          if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`    if `N` is odd
  !
  !that is, the aliased order returned by :f:func_inline:`fft`.
  !
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_ishift

  function cfft_1d_ishift(fin) result(fout)
  !This function is the inverse of :f:func_inline:`cfft_1d_shift`. It moves the
  !zero-frequency component of a complex array :f:var:`fin` of length :code:`N` from the 
  !center of the array back to the first position, and returns the result as an array 
  !of the same size. It acts as :f:func_inline:`rfft_1d_ishift` on complex data.
  !
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_ishift

  subroutine rfft_1d_ex(func)
  !This subroutine alternates the sign of the elements of a real array :f:var:`func`
  !in place,
  !
  !:math:`f(i) \to (-1)^{i-1} f(i)`
  !
  !for :math:`i \in [1,N]`, where :code:`N=size(func)`. Applying it twice returns 
  !the original array.
  !
    real(8),dimension(:)    :: func
    real(8)                 :: ex
    integer                 :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine rfft_1d_ex

  subroutine cfft_1d_ex(func)
  !This subroutine alternates the sign of the elements of a complex array :f:var:`func`
  !in place,
  !
  !:math:`f(i) \to (-1)^{i-1} f(i)`
  !
  !for :math:`i \in [1,N]`, where :code:`N=size(func)`. Applying it twice returns 
  !the original array.
  !
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
  !This function returns the maximum time of the time grid used with :f:func_inline:`FFT_signal`,
  !for an array of :f:var:`L` points with time step :f:var:`dt`:
  !
  !:math:`t_{max} = L \cdot dt / 2`
  !
    integer :: L
    real(8) :: dt
    real(8) :: fft_tmax
    fft_tmax = L*dt/2
  end function fft_tmax

  function fft_fmax(L,dt)
  !This function returns the maximum frequency of the frequency grid used with :f:func_inline:`FFT_signal`,
  !for a time step :f:var:`dt`:
  !
  !:math:`\omega_{max} = \pi / dt`
  !
  !The argument :f:var:`L` is not used, and is kept for symmetry with :f:func_inline:`fft_tmax`.
  !
    integer :: L
    real(8) :: dt
    real(8) :: fft_fmax
    fft_fmax = pi/dt
  end function fft_fmax




  function fft_tarray(L,dt) result(time)
  !This function returns the time array of :code:`L` equally spaced points 
  !in :math:`[-t_{max}, t_{max}]`, with :math:`t_{max}` = :f:func_inline:`fft_tmax` ( :f:var:`L` , :f:var:`dt` ) .
  !The endpoint is included, so the actual spacing is :math:`2 t_{max}/(L-1)`.
  !
    real(8)              :: dt
    integer              :: L
    real(8)              :: tmax
    real(8),dimension(L) :: time
    tmax = fft_tmax(L,dt)
    time = linspace(-tmax,tmax,L)
  end function fft_tarray

  function fft_farray(L,dt,df) result(freq)
  !This function returns the frequency array of :code:`L` equally spaced points 
  !in :math:`[-\omega_{max}, \omega_{max})`, with :math:`\omega_{max}` = :f:func_inline:`fft_fmax` ( :f:var:`L` , :f:var:`dt` ) .
  !The upper endpoint is excluded, so the spacing is :math:`2 \omega_{max}/L`.
  !If the optional argument :f:var:`df` is present, it returns this spacing.
  !
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
