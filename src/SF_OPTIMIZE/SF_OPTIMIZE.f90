MODULE SF_OPTIMIZE
!SciFortran module for minimization and root finding
  USE CGFIT_ROUTINES
  USE BROYDEN_ROUTINES
  !
  USE SF_CONSTANTS
  USE SF_LINALG, only: inv_sym
  private


  interface fmin_cg
     module procedure fmin_cg_df,fmin_cg_f
  end interface fmin_cg

  interface fmin_cgplus
     module procedure fmin_cgplus_df,fmin_cgplus_f
  end interface fmin_cgplus

  interface fmin_cgminimize
     module procedure fmin_cgminimize_func,fmin_cgminimize_sub
  end interface fmin_cgminimize

  interface leastsq
     module procedure :: leastsq_lmdif_func
     module procedure :: leastsq_lmdif_sub
     module procedure :: leastsq_lmder_func
     module procedure :: leastsq_lmder_sub
  end interface leastsq

  interface curvefit
     module procedure :: curvefit_lmdif_func
     module procedure :: curvefit_lmdif_sub
     module procedure :: curvefit_lmder_func
     module procedure :: curvefit_lmder_sub
  end interface curvefit

  interface dbrent
     module procedure :: dbrent_wgrad
     module procedure :: dbrent_nograd
  end interface dbrent

  interface fmin_bfgs
     module procedure :: bfgs_with_grad
     module procedure :: bfgs_no_grad
  end interface fmin_bfgs

  interface linear_mix
     module procedure :: d_linear_mix_1
     module procedure :: d_linear_mix_2
     module procedure :: d_linear_mix_3
     module procedure :: d_linear_mix_4
     module procedure :: d_linear_mix_5
     module procedure :: d_linear_mix_6
     module procedure :: d_linear_mix_7
     module procedure :: c_linear_mix_1
     module procedure :: c_linear_mix_2
     module procedure :: c_linear_mix_3
     module procedure :: c_linear_mix_4
     module procedure :: c_linear_mix_5
     module procedure :: c_linear_mix_6
     module procedure :: c_linear_mix_7
  end interface linear_mix


  interface adaptive_mix
     module procedure :: d_adaptive_mix
     module procedure :: c_adaptive_mix
  end interface adaptive_mix

  interface broyden_mix
     module procedure :: d_broyden_mix
     module procedure :: c_broyden_mix
  end interface broyden_mix


  interface fsolve
     module procedure :: fsolve_hybrd_func
     module procedure :: fsolve_hybrd_sub
     !
     module procedure :: fsolve_hybrj_func
     module procedure :: fsolve_hybrj_sub
  end interface fsolve




  !OPTIMIZATION:
  public   :: brent         !minimize a given a function of one-variable with a possible bracketing interval without using derivative information
  public   :: dbrent        !minimize a given a function of one-variable with a possible bracketing interval  using derivative information
  public   :: bracket       !Bracket the minimum of the function.
  !General purpose
  public   :: fmin                !Minimize a function using the Nelder-Mead downhill simplex algorithm.
  public   :: fmin_cg             !Conjugate-Gradient 1
  public   :: fmin_cgplus         !Conjugate-Gradient 2
  public   :: fmin_cgminimize     !Conjugate-Gradient 3 (very old f77)
  !Constrained (multivariate)
  public   :: fmin_bfgs    !Minimize a function using the BFGS algorithm.
  public   :: leastsq      !Minimize the sum of squares of a set of equations. Wrap MINPACK: lmdif/lmder
  public   :: curvefit     !Use non-linear least squares to fit a function, f, to data.
  !
  !> TODO:
  ! public :: fmin_powell  !Minimize a function using modified Powell’s method. This method
  ! public :: fmin_ncg     !Unconstrained minimization of a function using the Newton-CG method.
  ! public :: anneal       !Minimize a function using simulated annealing.
  ! public :: basinhopping ! Find the global minimum of a function using the basin-hopping algorithm ..



  !ROOT FINDING:
  public :: brentq
  public :: bisect
  public :: newton
  public :: fzero
  !Multidimensional
  !General nonlinear solvers:
  public :: fsolve              !
  public :: broyden1
  !Large-scale nonlinear solvers:

  !Fixed points accelerators:
  public :: linear_mix
  public :: adaptive_mix
  public :: broyden_mix
  ! public :: broyden2 !Find a root of a function, using Broyden’s second Jacobian approximation.
  ! public :: newton_krylov !Find a root of a function, using Krylov approximation for inverse Jacobian.
  ! public :: anderson !Find a root of a function, using (extended) Anderson mixing.


  real(8)                         :: df_eps=tiny(1d0)
  ! procedure(hybrd_func),pointer :: hybrd_funcv
  real(8), dimension(:),pointer   :: fmin_fvecp

contains

  ! Brent methods, including bracket
  include "brent.f90"


  ! INTERFACES TO MINPACK lmder/lmdif 
  include "leastsq.f90" 
  include "curvefit.f90"


  ! Minimizes a function using the Nelder-Mead algorithm.
  !    This routine seeks the minimum value of a user-specified function.
  !    Simplex function minimisation procedure due to Nelder and Mead (1965),
  !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
  !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
  !    25, 97) and Hill(1978, 27, 380-2)
  include "fmin_Nelder_Mead.f90" 


  ! Minimize the Chi^2 distance using conjugate gradient
  !     Adapted by FRPRM subroutine from NumRec (10.6),, 
  !     the Fletcher-Reeves-Polak-Ribiere minimisation is performed 
  include "fmin_cg.f90"


  !  Minimize the Chi^2 distance using conjugate gradient
  !     Adapted from unkown minimize.f routine.
  include "fmin_cg_minimize.f90"

  ! Conjugate Gradient methods for solving unconstrained nonlinear
  !  optimization problems:
  ! Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
  ! of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
  ! pp. 21-42. 
  include "fmin_cg_cgplus.f90"



  ! Constrained BFGS (L-BFGS_B) optimization problems:
  ! Ciyou Zhu , Richard H. Byrd , Peihuang Lu and Jorge Nocedal: "L-BFGS-B: 
  ! FORTRAN SUBROUTINES FOR LARGE-SCALE BOUND CONSTRAINED OPTIMIZATION"
  include "fmin_bfgs.f90"




  ! Mixing and Acceleration:
  include "linear_mix.f90"
  include "adaptive_mix.f90"
  include "broyden_mix.f90"


  ! Interface to MINPACK hybrd/hybrj: FSOLVE
  include "fsolve.f90"


  ! Broyden root finding method:
  include "broyden1.f90"


  ! Find root of scalar functions
  include "froot_scalar.f90"








  !           AUXILIARY JACOBIAN/GRADIENT CALCULATIONS
  !
  !          1 x N Jacobian (df_i/dx_j for i=1;j=1,...,N)
  subroutine sub_func_jacobian(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    n=size(x)
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    fvec = funcv(x)
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       wa1  = funcv(x)
       x(j) = temp
       fjac(j) = (wa1 - fvec)/h
    enddo
  end subroutine sub_func_jacobian

  ! formerly f_jac_1n_func(funcv,n,x)
  function func_func_jacobian(funcv,n,x) result(df)
    interface
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer               :: n
    real(8), dimension(n) :: x
    real(8), dimension(n) :: df
    call sub_func_jacobian(funcv,x,df)
  end function func_func_jacobian











END MODULE SF_OPTIMIZE
