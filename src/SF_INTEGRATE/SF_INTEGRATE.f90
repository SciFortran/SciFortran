! ABOUT QUADPACK:
! How to decide what routine to use, if your integration region is finite:
!
! If you can factor the integrand as F(X)=W(X)*G(X), where G is smooth on [A,B] and W(X)=COS(OMEGA*X) or SIN(OMEGA*X) then use QAWO.
! Otherwise, if you can factor F(X)=W(X)*G(X) where G is smooth and W(X)=(X-A)**ALFA * (B-X)**BETA * (LOG(X-A))**L * (LOG(B-X))**K with K, L = 0 or 1, and ALFA, BETA greater than -1, then use QAWS.
!
! Otherwise, if you can factor F(X)=W(X)*G(X) where G is smooth and W(X)=1/(X-C) for some constant C, use QAWC.
! Otherwise, if you do not care too much about possible inefficient use of computer time, and do not want to further analyze the problem, use QAGS.
! Otherwise, if the integrand is smooth, use QNG or QAG.
! Otherwise, if there are discontinuities or singularities of the integrand or of its derivative, and you know where they are, split the integration range at these points and analyze each subinterval. You can also use QAGP, which is to be provided with the x-locations of the singularities or discontinuities.
! Otherwise, if the integrand has end point singularities, use QAGS.
! Otherwise, if the integrand has an oscillatory behavior of nonspecific type, and no singularities, use QAG with KEY=6.
! Otherwise, use QAGS.
!  
! Routines for an infinite region:
! If the integrand decays rapidly to zero, truncate the interval and use the finite interval decision tree.
! Otherwise, if the integrand oscillates over the entire infinite range, and the integral is a Fourier transform, use QAWF.
! Or, if the integrand oscillates over the entire infinite range, but is not a Fourier transform, then sum the successive positive and negative contributions by integrating between the zeroes of the integrand. Apply convergence acceleration with QELG.
! Otherwise, if you are not constrained by computer time, and do not wish to analyze the problem further, use QAGI.
! Otherwise, if the integrand has a non-smooth behavior in the range, and you know where it occurs, split off these regions and use the appropriate finite range routines to integrate over them. Then begin this tree again to handle the remainder of the region.
! Otherwise, truncation of the interval, or application of a suitable transformation for reducing the problem to a finite range may be possible. And you may also call QAGI.

module SF_INTEGRATE
!SciFortran module for function integration
  USE GAUSS_QUADRATURE
  implicit none
  private


  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)
  real(8),parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510d0


  !TRAPEZIODAL RULE:
  interface trapz
  !This function evaluates the integral of a discretized real or complex function
  !using the trapezoidal rule (2nd order), or of a function passed as argument.
  !The specific procedures cover
  !
  !* a sampled function on a uniform grid defined by the limits :code:`a`, :code:`b`
  !* a sampled function on a uniform grid defined by the step :code:`dh`
  !* a sampled function on a non-uniform grid
  !* a function evaluated between the limits :code:`a`, :code:`b`
  !* a function evaluated on a non-uniform grid
  !
  !The quadrature weights are provided by :f:func_inline:`get_quadrature_weights`
  !with :code:`nrk=2`.
  !
     module procedure :: d_trapz_ab_sample
     module procedure :: c_trapz_ab_sample
     module procedure :: d_trapz_dh_sample
     module procedure :: c_trapz_dh_sample
     module procedure :: d_trapz_nonlin_sample
     module procedure :: c_trapz_nonlin_sample
     module procedure :: d_trapz_ab_func
     module procedure :: c_trapz_ab_func
     module procedure :: d_trapz_nonlin_func
     module procedure :: c_trapz_nonlin_func
  end interface trapz

  
  interface trapz2d
  !This function evaluates the 2-dimensional integral of a real or complex function
  !using the trapezoidal rule (2nd order). The specific procedures cover a function
  !of two variables, a function evaluated recursively as a nested 1-dimensional
  !integral, and a sampled function on a 2-dimensional grid.
  !
     module procedure :: d_trapz2d_func
     module procedure :: c_trapz2d_func
     module procedure :: d_trapz2d_func_recursive
     module procedure :: c_trapz2d_func_recursive
     module procedure :: d_trapz2d_sample
     module procedure :: c_trapz2d_sample
  end interface trapz2d





  !SIMPSON'S RULE
  interface simps
  !This function evaluates the integral of a discretized real or complex function
  !using Simpson's rule (4th order), or of a function passed as argument.
  !The specific procedures cover
  !
  !* a sampled function on a uniform grid defined by the step :code:`dh`
  !* a sampled function on a uniform grid defined by the limits :code:`a`, :code:`b`
  !* a sampled function on a non-uniform grid
  !* a function evaluated between the limits :code:`a`, :code:`b`
  !* a function evaluated on a non-uniform grid
  !
  !The quadrature weights are provided by :f:func_inline:`get_quadrature_weights`
  !with :code:`nrk=4`.
  !
     module procedure :: d_simpson_dh_sample
     module procedure :: c_simpson_dh_sample
     module procedure :: d_simpson_ab_sample
     module procedure :: c_simpson_ab_sample
     module procedure :: d_simpson_nonlin_sample
     module procedure :: c_simpson_nonlin_sample
     !FUNCTION:
     module procedure :: d_simps_ab_func
     module procedure :: c_simps_ab_func
     module procedure :: d_simps_nonlin_func
     module procedure :: c_simps_nonlin_func
  end interface simps
  !
  interface simps2d
  !This function evaluates the 2-dimensional integral of a real or complex function
  !using Simpson's rule (4th order). The specific procedures cover a function
  !of two variables, a function evaluated recursively as a nested 1-dimensional
  !integral, and a sampled function on a 2-dimensional grid.
  !
     module procedure :: d_simps2d_func
     module procedure :: c_simps2d_func
     module procedure :: d_simps2d_func_recursive
     module procedure :: c_simps2d_func_recursive
     module procedure :: d_simps2d_sample
     module procedure :: c_simps2d_sample
  end interface simps2d
  !


  !QUADPACK global interface
  interface quad
  !This function evaluates the 1-dimensional integral of a function, or of a sampled
  !function, using the adaptive routines of the QUADPACK library. The routine
  !is chosen according to the type of integrand and integration region, as described
  !in the guide at the top of the module (for example :code:`QAGS`, :code:`QAGI`,
  !:code:`QAWO`, :code:`QAWF`, :code:`QAWC`, :code:`QAWS`).
  !The optional argument :f:var:`weight_func` selects the weight function :math:`W(x)`
  !of the integrand :math:`W(x) f(x)` in the routines that integrate with a weight.
  !
  !For :code:`QAWO` and :code:`QAWF`, which also require :f:var:`omega`, the accepted values are
  !
  !* :code:`weight_func=1`: :math:`W(x) = \cos(\omega x)`
  !* :code:`weight_func=2`: :math:`W(x) = \sin(\omega x)`
  !
  !For :code:`QAWS`, which also requires :f:var:`alfa` and :f:var:`beta`, the accepted values are
  !
  !* :code:`weight_func=1`: :math:`W(x) = (x-a)^{\alpha} (b-x)^{\beta}`
  !* :code:`weight_func=2`: :math:`W(x) = (x-a)^{\alpha} (b-x)^{\beta} \ln(x-a)`
  !* :code:`weight_func=3`: :math:`W(x) = (x-a)^{\alpha} (b-x)^{\beta} \ln(b-x)`
  !* :code:`weight_func=4`: :math:`W(x) = (x-a)^{\alpha} (b-x)^{\beta} \ln(x-a) \ln(b-x)`
  !
  !with :math:`\alpha, \beta > -1`. The program stops if :f:var:`weight_func` is outside
  !these ranges. It is ignored by the other routines.
  !
     module procedure :: quad_func
     module procedure :: quad_sample
  end interface quad



  !1D Adaptive QUADPACK
  public :: quad
  !nD Adapative GAUSS RULE 6-14
  public :: gauss_quad
  public :: integrate

  !1D simple
  public :: trapz
  public :: simps
  !2D simple
  public :: trapz2d
  public :: simps2d

  !KRAMERS-KRONIG:
  public :: kronig

  !AUX:
  public :: get_quadrature_weights



  !<TODO
  ! add cubature methods for arbitrary domains
  ! add Montecarlo base 1d/2d/3d integrals
  !>TODO



contains



  !PURPOSE:
  ! evaluate 1D and 2D integrals using trapz (2nd order) and simps (4th order) rule.
  include "integrate_func_1d.f90"
  include "integrate_sample_1d.f90"
  include "integrate_func_2d.f90"
  include "integrate_sample_2d.f90"




  !PURPOSE: 
  !generic interface to QUADPACK routines:
  !non-automatic integrator:
  !estimate the integral on [a,b] using 15, 21,..., 61
  !point rule and return an error estimate.
  ! public :: QNG 
  ! public :: QK15
  ! public :: QK21
  ! public :: QK31
  ! public :: QK41
  ! public :: QK51
  ! public :: QK61
  !automatic integrator:
  !without weights:
  ! public :: QAG  !globally adaptive integrator
  ! public :: QAGI !integration over infinite intervals
  ! public :: QAGS !globally adaptive interval subdivisio with epsilon extrapolation
  ! public :: QAGP !serves the same purposes as qags, but user-supplied information about singularities
  !with weights:
  ! public :: QAWO !integration of cos(omega*x)*f(x) sin(omega*x)*f(x) over a finite interval (a,b).
  ! public :: QAWF !Fourier cosine or fourier sine transform of f(x)
  ! public :: QAWC !computes the cauchy principal value of f(x)/(x-c) over a finite interval (a,b)
  ! public :: QAWS !integrates w(x)*f(x) over (a,b) with a < b finite,
  ! !              ! and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
  ! !              ! where v(x) = 1 or log(x-a) or log(b-x)
  ! !              !              or log(x-a)*log(b-x)
  ! !              ! and   -1 < alfa, -1 < beta.
  include "integrate_quad_func.f90"
  include "integrate_quad_sample.f90"



  !PURPOSE  : Perform a fast Kramers-K\"onig integration: 
  function kronig(fi,wr,M) result(fr)
  !This function evaluates the Kramers-Kronig transform of a function of a
  !real variable. Given the values :f:var:`fi` of the imaginary part on the
  !uniform grid :f:var:`wr` of :f:var:`M` points, it returns the real part :f:var:`fr`
  !on the same grid,
  !
  !:math:`f_r(w_i) = \frac{1}{\pi} \, \mathcal{P} \int_{w_1}^{w_M} \frac{f_i(w')}{w'-w_i} \, dw'`
  !
  !The principal value is evaluated by subtracting the singularity,
  !
  !:math:`\mathcal{P} \int \frac{f_i(w')}{w'-w_i} dw' = \int \frac{f_i(w')-f_i(w_i)}{w'-w_i} dw' + f_i(w_i) \, \mathcal{P} \int \frac{dw'}{w'-w_i}`
  !
  !where the first integral is a sum with weight :code:`dh=wr(2)-wr(1)`, using
  !the finite-difference derivative of :math:`f_i` at the point :math:`w'=w_i`,
  !and the second is known analytically, :math:`\ln\left(\frac{w_M-w_i}{w_i-w_1}\right)`.
  !This term is set to zero at the end points :math:`i=1,M`, where it diverges.
  !The cost scales as :math:`M^2`.
  !
    integer :: i,j,M
    real(8),dimension(M) :: fi,wr,fr
    real(8),dimension(M) :: logo,deriv
    real(8) :: dh,sum
    dh=wr(2)-wr(1)
    logo=0.d0
    do i=2,M-1
       logo(i) = log( (wr(M)-wr(i))/(wr(i)-wr(1)) )
    enddo
    deriv(1)= (fi(2)-fi(1))/dh
    deriv(M)= (fi(M)-fi(M-1))/dh
    do i=2,M-1
       deriv(i) = (fi(i+1)-fi(i-1))/(2*dh)
    enddo
    fr=0.d0
    do i=1,M
       sum=0.d0
       do j=1,M
          if(i/=j)then
             sum=sum+(fi(j)-fi(i))*dh/(wr(j)-wr(i))
          else
             sum=sum+deriv(i)*dh
          endif
       enddo
       fr(i) = (sum + fi(i)*logo(i))/pi
    enddo
    return
  end function kronig














  !PURPOSE: obtain quadrature weights for higher order integration (2,4)
  subroutine get_quadrature_weights(wt,nrk)
  !This subroutine returns the quadrature weights :f:var:`wt` for the integration
  !of a function sampled on :code:`N=size(wt)` points with unit step, such that
  !
  !:math:`\int f(x) dx \approx dh \sum_{i=1}^{N} wt(i) \cdot f(x_i)`
  !
  !The order of the rule is selected by the optional argument :f:var:`nrk`
  !
  !* :code:`nrk=2`: trapezoidal rule (2nd order), :math:`wt = [1/2, 1, ..., 1, 1/2]`
  !* :code:`nrk=4`: Simpson's rule (4th order), the default
  !
  !For :code:`nrk=4` the weights depend on :code:`N`
  !
  !* :code:`N=1`: :math:`wt = 1`
  !* :code:`N=2`: :math:`wt = [1/2, 1/2]`
  !* :code:`N=3`: Simpson's rule :math:`[1/3, 4/3, 1/3]`
  !* :code:`N=4`: Simpson's 3/8 rule :math:`[3/8, 9/8, 9/8, 3/8]`
  !* :code:`N>=5` odd: composite Simpson's rule :math:`[1/3, 4/3, 2/3, ..., 4/3, 1/3]`
  !* :code:`N>=6` even: composite Simpson's rule, closed by a Simpson's 3/8 rule on the last four points
  !
  !The program stops if :f:var:`nrk` is neither 2 nor 4.
  !
    real(8),dimension(:) :: wt    !Quadrature weights
    integer,optional     :: nrk   !Order of the rule
    integer              :: nrk_
    integer              :: N
    nrk_=4;if(present(nrk))nrk_=nrk
    N=size(wt)
    if(nrk_==4)then
       select case(n)           !n>=3
       case (1)
          wt = 1.d0
       case (2)
          wt = 0.5d0
       case (3)                 !simpson's rule
          wt(1)=1.d0/3.d0
          wt(2)=4.d0/3.d0
          wt(3)=1.d0/3.d0
       case(4)                  !simpson's 3/8 rule
          wt(1)=3.d0/8.d0
          wt(2)=9.d0/8.d0
          wt(3)=9.d0/8.d0
          wt(4)=3.d0/8.d0
       case(5)                  !Simpson's rule (E,O n)
          wt(1)=1.d0/3.d0
          wt(2)=4.d0/3.d0
          wt(3)=2.d0/3.d0
          wt(4)=4.d0/3.d0
          wt(5)=1.d0/3.d0
       case default            !Simpson's rule n>=6
          if(mod(n-1,2)==0)then
             wt(1)=1.d0/3.d0
             wt(n)=1.d0/3.d0
             wt(2:n-1:2)=4.d0/3.d0
             wt(3:n-2:2)=2.d0/3.d0
          else
             wt(1)=1.d0/3.d0
             wt(2:n-4:2)=4.d0/3.d0
             wt(3:n-5:2)=2.d0/3.d0
             wt(n-3)=17.d0/24.d0
             wt(n-2)=9.d0/8.d0
             wt(n-1)=9.d0/8.d0
             wt(n)=3.d0/8.d0
          endif
          ! case default             !Simpson's rule n>=6
          !    wt(1)=3.d0/8.d0
          !    wt(2)=7.d0/6.d0
          !    wt(3)=23.d0/24.d0
          !    wt(4:n-3)=1.d0
          !    wt(n-2)=23.d0/24.d0
          !    wt(n-1)=7.d0/6.d0
          !    wt(n)=3.d0/8.d0
       end select
    elseif(nrk_==2)then
       wt(1) = 0.5d0
       wt(2:n-1)=1.d0
       wt(n) = 0.5d0
    else
       stop "error in +get_quadrature_weights: nrk != 2,4" 
    end if
  end subroutine get_quadrature_weights





  function sf_integrate_linspace(start,stop,num,istart,iend,mesh) result(array)
    integer          :: num,i
    real(8)          :: start,stop,step,array(num)
    logical,optional :: istart,iend
    logical          :: startpoint_,endpoint_
    real(8),optional :: mesh
    if(num<0)stop "linspace: N<0, abort."
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(startpoint_.AND.endpoint_)then
       if(num<2)stop "linspace: N<2 with both start and end points"
       step = (stop-start)/real(num-1,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(startpoint_.AND.(.not.endpoint_))then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(.not.startpoint_.AND.endpoint_)then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    else
       step = (stop-start)/real(num+1,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    endif
    if(present(mesh))mesh=step
  end function sf_integrate_linspace
















  subroutine polint(xa,ya,x,y,dy)
    real(8), dimension(:), intent(in) :: xa,ya
    real(8), intent(in)          :: x
    real(8), intent(out)         :: y,dy
    integer                      :: m,n,ns
    real(8), dimension(size(xa)) :: c,d,den,ho
    n=assert_eq2(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0))then
          print*,'polint: calculation failure'
          stop
       endif
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint
  function locate(xx,x)
    real(8), dimension(:), intent(in) :: xx
    real(8), intent(in) :: x
    integer :: locate
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  end function locate
  function iminloc(arr)
    real(8), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc
  function assert_eq2(n1,n2,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2

end module SF_INTEGRATE





















