MODULE SF_OPTIMIZE
!SciFortran module for minimization and root finding
  USE CGFIT_ROUTINES
  USE BROYDEN_ROUTINES
  !
  USE SF_CONSTANTS
  USE SF_LINALG, only: inv_sym
  private


  interface fmin_cg
     !Generic interface to minimize a function of several variables with the Fletcher-Reeves-Polak-Ribiere conjugate gradient
     !method (Numerical Recipes 10.6), with the line minimizations done with the Brent method. The gradient is supplied by the
     !user (:f:func_inline:`fmin_cg_df`) or computed numerically (:f:func_inline:`fmin_cg_f`). The location of the minimum is
     !returned in :f:var:`p`, the number of iterations in :f:var:`iter` and the value of the function in :f:var:`fret`.
     !
     !The convergence is controlled by :f:var:`istop` and :f:var:`ftol`: with :math:`a = |F_n-F_{n-1}|/(1+|F_n|)` and
     !:math:`b = \|p_n-p_{n-1}\|^2/(1+\|p_n\|^2)`, the iteration stops when :math:`a<ftol` and :math:`b<ftol`
     !(:code:`istop=0`), when :math:`a<ftol` (:code:`istop=1`) or when :math:`b<ftol` (:code:`istop=2`). The last values of
     !:math:`a` and :math:`b` are returned in :f:var:`err`.
     !
     module procedure :: fmin_cg_df
     module procedure :: fmin_cg_f
  end interface fmin_cg

  interface fmin_cgplus
     !Generic interface to minimize a function of several variables with the nonlinear conjugate gradient method of the CG+
     !package (Gilbert and Nocedal), a very old Fortran 77 code. The gradient is supplied by the user
     !(:f:func_inline:`fmin_cgplus_df`) or computed numerically (:f:func_inline:`fmin_cgplus_f`). The variant of the method is
     !selected by :f:var:`imethod`, and the iteration stops when every component of the gradient is smaller than
     !:math:`ftol(1+|f|)`. The location of the minimum is returned in :f:var:`p`.
     !
     module procedure :: fmin_cgplus_df
     module procedure :: fmin_cgplus_f
  end interface fmin_cgplus

  interface fmin_cgminimize
     !Generic interface to minimize a function of several variables with a conjugate gradient routine adapted from a very old
     !Fortran 77 code, a quasi-Newton minimizer with numerical gradient. The function is a function
     !(:f:func_inline:`fmin_cgminimize_func`) or a subroutine (:f:func_inline:`fmin_cgminimize_sub`). Two versions of the
     !minimizer are available, selected by :f:var:`new_version`. The location of the minimum is returned in :f:var:`p`.
     !
     module procedure :: fmin_cgminimize_func
     module procedure :: fmin_cgminimize_sub
  end interface fmin_cgminimize

  interface leastsq
     !Generic interface to minimize the sum of squares of :math:`m` nonlinear functions of :math:`n` unknowns, :math:`m\ge n`,
     !i.e. to find :f:var:`a` so that :math:`f(a)^T f(a)` is minimum, with the Levenberg-Marquardt routines of MINPACK. The
     !Jacobian is computed numerically (:code:`lmdif`, :f:func_inline:`leastsq_lmdif_func` and
     !:f:func_inline:`leastsq_lmdif_sub`) or supplied by the user (:code:`lmder`, :f:func_inline:`leastsq_lmder_func` and
     !:f:func_inline:`leastsq_lmder_sub`). The functions are a function :code:`func(a,m)` returning a vector of size
     !:math:`m`, or a subroutine :code:`func(a,m,f)`; the Jacobian :code:`dfunc(a,m)` is a matrix of shape :math:`[m,n]`.
     !
     !The program stops, writing the reason to the file :code:`LEASTSQ_ERROR.err`, if the status returned by MINPACK is not 1.
     !
     module procedure :: leastsq_lmdif_func
     module procedure :: leastsq_lmdif_sub
     module procedure :: leastsq_lmder_func
     module procedure :: leastsq_lmder_sub
  end interface leastsq

  interface curvefit
     !Generic interface to fit a model function to data with the non-linear least squares Levenberg-Marquardt routines of
     !MINPACK: the parameters :f:var:`a` of :code:`model_func(x,a)` are adjusted to minimize the sum of the squares of
     !:code:`model_func(xdata,a)-ydata`. The Jacobian is computed numerically (:code:`lmdif`,
     !:f:func_inline:`curvefit_lmdif_func` and :f:func_inline:`curvefit_lmdif_sub`) or supplied by the user (:code:`lmder`,
     !:f:func_inline:`curvefit_lmder_func` and :f:func_inline:`curvefit_lmder_sub`). The model is a function
     !:code:`model_func(x,a)` returning a vector with the size of :code:`x`, or a subroutine :code:`model_func(x,a,f)`; the
     !derivative :code:`model_dfunc(x,a)` is a matrix of shape :math:`[size(x),size(a)]`. The model is evaluated on the whole
     !array :code:`xdata` at once.
     !
     module procedure :: curvefit_lmdif_func
     module procedure :: curvefit_lmdif_sub
     module procedure :: curvefit_lmder_func
     module procedure :: curvefit_lmder_sub
  end interface curvefit

  interface dbrent
     !Generic interface to find the minimum of a function of one variable, :f:var:`func`, with the Brent method using the
     !derivative: the user supplies the derivative :f:var:`dfunc` (:f:func_inline:`dbrent_wgrad`), or it is computed
     !numerically (:f:func_inline:`dbrent_nograd`). The abscissa of the minimum is returned in :f:var:`xmin`. The bracketing
     !of the minimum, :f:var:`brack`, works as in :f:func_inline:`brent`: a bracketing triplet, or two starting points from
     !which it is searched with :f:func_inline:`bracket`, with default points 0 and 1. The program stops if :f:var:`brack` has
     !a single element.
     !
     module procedure :: dbrent_wgrad
     module procedure :: dbrent_nograd
  end interface dbrent

  interface fmin_bfgs
     !Generic interface to minimize a function of several variables with the limited-memory BFGS algorithm with bounds,
     !L-BFGS-B (Zhu, Byrd, Lu and Nocedal), subject to the optional bounds :math:`l\le x\le u`. The gradient is supplied by
     !the user (:f:func_inline:`bfgs_with_grad`) or computed numerically (:f:func_inline:`bfgs_no_grad`). The type of the
     !bounds is given by :f:var:`nbd`, the bounds are ignored if :f:var:`nbd` is absent. The location of the minimum is
     !returned in :f:var:`x`.
     !
     module procedure :: bfgs_with_grad
     module procedure :: bfgs_no_grad
  end interface fmin_bfgs

  interface linear_mix
     !Generic interface for the linear mixing of an array of rank 1 to 7, real (:f:func_inline:`d_linear_mix_1` to
     !:f:func_inline:`d_linear_mix_7`) or complex (:code:`c_linear_mix_1` to :code:`c_linear_mix_7`):
     !:math:`x \leftarrow x + \alpha F_x`, with :math:`F_x = x_{out}-x_{in}` the residual of the iteration.
     !
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
     !Generic interface for the adaptive linear mixing, real (:f:func_inline:`d_adaptive_mix`) or complex
     !(:f:func_inline:`c_adaptive_mix`): :math:`x \leftarrow x + \beta F_x`, with an element-wise mixing parameter
     !:math:`\beta_j`, which starts from :f:var:`alpha`, is increased by :f:var:`alpha` at every iteration in which the
     !residual :math:`F_{x,j}` keeps its sign, up to 1, and is reset to :f:var:`alpha` if the sign changes (oscillations). The
     !previous residual is stored between calls: it is reset when :f:var:`iter` is 1.
     !
     module procedure :: d_adaptive_mix
     module procedure :: c_adaptive_mix
  end interface adaptive_mix

  interface broyden_mix
     !Generic interface for the modified Broyden mixing of Johnson (Phys. Rev. B 38, 12807, 1988), real
     !(:f:func_inline:`d_broyden_mix`) or complex (:f:func_inline:`c_broyden_mix`): the mixed vector is built from the
     !residual :math:`F_x = x_{out}-x_{in}` and from the history of the last :f:var:`M` iterations, stored between calls and
     !reset when :f:var:`iter` is 1. The first iteration is a linear mixing with :f:var:`alpha`, as are all the iterations if
     !:code:`M=0`. The weight of the initial iteration is :f:var:`w0`.
     !
     module procedure :: d_broyden_mix
     module procedure :: c_broyden_mix
  end interface broyden_mix


  interface fsolve
     !Generic interface to solve a system of :math:`n` nonlinear equations in :math:`n` unknowns, :math:`f(x)=0`, with the
     !hybrid Powell method of MINPACK. The Jacobian is computed numerically (:code:`hybrd`, :f:func_inline:`fsolve_hybrd_func`
     !and :f:func_inline:`fsolve_hybrd_sub`) or supplied by the user (:code:`hybrj`, :f:func_inline:`fsolve_hybrj_func` and
     !:f:func_inline:`fsolve_hybrj_sub`). The system is a function :code:`func(x)` returning a vector with the size of
     !:code:`x`, or a subroutine :code:`func(x,f)`; the Jacobian :code:`dfunc(x)` is a matrix of shape :math:`[n,n]`. If
     !:f:var:`icheck` is true, the default, the program stops, writing the reason to a file, when the status returned by
     !MINPACK is not a success.
     !
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
