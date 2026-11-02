module SF_PRIMME
  USE ISO_C_BINDING, only: c_ptr,c_int,c_int64_t,c_double,c_funloc,c_f_pointer,c_associated
#ifdef _MPI
  use MPI
#endif
  implicit none
  private


  !Include PRIMME header for F90 interface
  include 'primme_f90.inc'


  interface primme_set
    module procedure :: primme_set_int
    module procedure :: primme_set_double
  end interface primme_set


  !Abstract interfaces to MatVec and GlobalSumReal in PRIMME. These
  !require C bindings
  abstract interface
    subroutine primme_matvec_d(n,vin,vout)
      import c_double
      integer, intent(in)         :: n
      real(c_double), intent(in)  :: vin(n)
      real(c_double), intent(out) :: vout(n)
    end subroutine primme_matvec_d
    !
    subroutine primme_matvec_z(n,vin,vout)
      import c_double
      integer, intent(in)            :: n
      complex(c_double), intent(in)  :: vin(n)
      complex(c_double), intent(out) :: vout(n)
    end subroutine primme_matvec_z
    !
    subroutine primme_global_sum_real(sendBuf,recvBuf,count,primme,ierr) bind(c)
      import c_ptr,c_int
      type(c_ptr), value         :: sendBuf,recvBuf,primme
      integer(c_int), intent(in) :: count
      integer(c_int)             :: ierr
    end subroutine primme_global_sum_real
  end interface
  !
  !Local pointers to MatVec (Dble and Zcmplx)
  procedure(primme_matvec_d), pointer, save :: user_matvec_d => null()
  procedure(primme_matvec_z), pointer, save :: user_matvec_z => null()
  !
  !
  !
  !
  !Public interfaces
  interface primme_eigsh
    module procedure sp_eigh_primme_d
    module procedure sp_eigh_primme_z
#ifdef _MPI
    module procedure sp_eigh_primme_mpi_d
    module procedure sp_eigh_primme_mpi_z
#endif
  end interface primme_eigsh
  !
  interface sp_eigh_primme
    module procedure sp_eigh_primme_d
    module procedure sp_eigh_primme_z
#ifdef _MPI
    module procedure sp_eigh_primme_mpi_d
    module procedure sp_eigh_primme_mpi_z
#endif
  end interface sp_eigh_primme
  !
  !
  public :: primme_matvec_d
  public :: primme_matvec_z
  public :: primme_global_sum_real
  !
  public :: primme_eigsh
  public :: sp_eigh_primme
  !
  public :: PRIMME_DYNAMIC
  public :: PRIMME_DEFAULT_MIN_TIME
  public :: PRIMME_DEFAULT_MIN_MATVECS
  !
  !
  !Local variables
  logical, save :: active_mpi = .false.
#ifdef _MPI
  integer, save :: active_comm = MPI_COMM_NULL
#endif


contains

  subroutine sp_eigh_primme_d(MatVec,eval,evec,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm)
    procedure(primme_matvec_d), pointer, intent(in) :: MatVec
    real(c_double), intent(out)                     :: eval(:)
    real(c_double), intent(inout)                   :: evec(:,:)
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,solve_status,chosen_method
    integer(c_int64_t)                              :: n,nconv_local,nmat,niter_,ninit
    integer                                         :: j
    !
    n=int(size(evec,1),c_int64_t)
    if(size(evec,2)<size(eval)) error stop 'sp_eigh_primme_d: evec has too few columns'
    !
    chosen_method=PRIMME_DYNAMIC; if(present(method))chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize))ninit=initSize
    !
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t) .OR. ninit>int(size(evec,2),c_int64_t))&
      error stop 'sp_eigh_primme_d: invalid initSize'
    !
    !Set the MatVec pointer
    user_matvec_d=>MatVec
    !
    active_mpi=.false.
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params)) error stop 'sp_eigh_primme_d: primme_params_create failed'
    !
    !
    !Set PRIMME internal parameters
    call primme_set(params,PRIMME_n,n,'PRIMME_n')
    call primme_set(params,PRIMME_numEvals,int(size(eval),c_int64_t),'PRIMME_numEvals')
    call primme_set(params,PRIMME_target,primme_smallest,'PRIMME_target')
    call primme_set(params,PRIMME_initSize,ninit,'PRIMME_initSize')
    if(present(Nitermax))call primme_set(params,PRIMME_maxOuterIterations,Nitermax,'PRIMME_maxOuterIterations')
    if(present(Nblock))call primme_set(params,PRIMME_maxBasisSize,Nblock,'PRIMME_maxBasisSize')
    if(present(tol))call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !
    !Set the internal PRIMME MatVec callback dble precision
    status=primme_set_member_matvec_double(params,PRIMME_matrixMatvec,mv_d)
    call check(status,'PRIMME_matrixMatvec')
    !
    !Set the internal PRIMME pre-condition method: PRIMME_DEFAULT, PRIMME_MIN_TIME, PRIMME_MIN_MATVECS
    status=primme_set_method(chosen_method,params)
    call check(status,'primme_set_method')
    !
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=dprimme(eval,evec(:,1:size(eval)),rnorms,params)
    !
    !
    call primme_collect_stats(params,nconv_local,nmat,niter_)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_d)
  end subroutine sp_eigh_primme_d




  subroutine sp_eigh_primme_z(MatVec,eval,evec,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm)
    procedure(primme_matvec_z), pointer, intent(in) :: MatVec
    real(c_double), intent(out)                     :: eval(:)
    complex(c_double), intent(inout)                :: evec(:,:)
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,chosen_method
    integer(c_int64_t)                              :: n,nconv_local,nmat,niter_,ninit
    !
    n=int(size(evec,1),c_int64_t)
    if(size(evec,2)<size(eval)) error stop 'sp_eigh_primme_z: evec has too few columns'
    chosen_method=PRIMME_DYNAMIC; if(present(method)) chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize)) ninit=initSize
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t) .OR. ninit>int(size(evec,2),c_int64_t))&
      error stop 'sp_eigh_primme_z: invalid initSize'
    !
    !Set the MatVec pointer
    user_matvec_z=>MatVec
    !
    active_mpi=.false.
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params)) error stop 'sp_eigh_primme_z: primme_params_create failed'
    !
    !Set PRIMME internal parameters
    call primme_set(params,PRIMME_n,n,'PRIMME_n')
    call primme_set(params,PRIMME_numEvals,int(size(eval),c_int64_t),'PRIMME_numEvals')
    call primme_set(params,PRIMME_target,primme_smallest,'PRIMME_target')
    call primme_set(params,PRIMME_initSize,ninit,'PRIMME_initSize')
    if(present(Nitermax)) call primme_set(params,PRIMME_maxOuterIterations,Nitermax,'PRIMME_maxOuterIterations')
    if(present(Nblock))   call primme_set(params,PRIMME_maxBasisSize,Nblock,'PRIMME_maxBasisSize')
    if(present(tol))      call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !Set the internal PRIMME MatVec callback complex precision
    status=primme_set_member_matvec_complexdouble(params,PRIMME_matrixMatvec,mv_z)
    call check(status,'PRIMME_matrixMatvec')
    !
    !Set the internal PRIMME method
    status=primme_set_method(chosen_method,params)
    call check(status,'primme_set_method')
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=zprimme(eval,evec(:,1:size(eval)),rnorms,params)
    !
    call primme_collect_stats(params,nconv_local,nmat,niter_)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_z)
  end subroutine sp_eigh_primme_z




#ifdef _MPI
  subroutine sp_eigh_primme_mpi_d(MpiComm,MatVec,eval,evec,GlobalSumReal,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm)
    integer, intent(in)                             :: MpiComm
    procedure(primme_matvec_d), pointer, intent(in) :: MatVec
    real(c_double), intent(out)                     :: eval(:)
    real(c_double), intent(inout)                   :: evec(:,:)
    procedure(primme_global_sum_real), optional     :: GlobalSumReal
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,mpierr,chosen_method
    integer(c_int64_t)                              :: n,nlocal,nconv_local,nmat,niter_,ninit,nprocs,rank
    integer                                         :: ierr
    !
    !Set MPI global dimensions
    nlocal=int(size(evec,1),c_int64_t)
    call MPI_Allreduce(nlocal,n,1,MPI_INTEGER8,MPI_SUM,MpiComm,mpierr)
    if(mpierr/=MPI_SUCCESS)error stop 'sp_eigh_primme_mpi_d: MPI_Allreduce failed'
    call MPI_Comm_size(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_size')
    nprocs=int(ierr,c_int64_t)
    call MPI_Comm_rank(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_rank')
    rank=int(ierr,c_int64_t)
    !
    if(size(evec,2)<size(eval))error stop 'sp_eigh_primme_mpi_d: evec has too few columns'
    !
    chosen_method=PRIMME_DYNAMIC; if(present(method)) chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize)) ninit=initSize
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t)) &
      error stop 'sp_eigh_primme_mpi_d: invalid initSize'
    !
    !Set the MatVec pointer
    user_matvec_d=>MatVec
    active_mpi=.true.; active_comm=MpiComm
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params))error stop 'sp_eigh_primme_mpi_d: primme_params_create failed'
    !
    !Set PRIMME internal parameters
    call primme_set(params,PRIMME_n,n,'PRIMME_n')
    call primme_set(params,PRIMME_nLocal,nlocal,'PRIMME_nLocal')
    call primme_set(params,PRIMME_numProcs,nprocs,'PRIMME_numProcs')
    call primme_set(params,PRIMME_procID,rank,'PRIMME_procID')
    call primme_set(params,PRIMME_numEvals,int(size(eval),c_int64_t),'PRIMME_numEvals')
    call primme_set(params,PRIMME_target,primme_smallest,'PRIMME_target')
    call primme_set(params,PRIMME_initSize,ninit,'PRIMME_initSize')
    if(present(Nitermax)) call primme_set(params,PRIMME_maxOuterIterations,Nitermax,'PRIMME_maxOuterIterations')
    if(present(Nblock))   call primme_set(params,PRIMME_maxBasisSize,Nblock,'PRIMME_maxBasisSize')
    if(present(tol))      call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !Set the internal PRIMME MatVec callback double precision
    status=primme_set_member_matvec_double(params,PRIMME_matrixMatvec,mv_d)
    call check(status,'PRIMME_matrixMatvec')
    !
    !Set the internal PRIMME GlobalSumReal callback
    call primme_set(params,PRIMME_globalSumReal_type,primme_op_double,'PRIMME_globalSumReal_type')
    if(present(GlobalSumReal))then
      status=primme_set_member(params,PRIMME_globalSumReal,c_funloc(GlobalSumReal))
    else
      status=primme_set_member(params,PRIMME_globalSumReal,c_funloc(default_global_sum_real))
    endif
    call check(status,'PRIMME_globalSumReal')
    !
    !Set the internal PRIMME method
    status=primme_set_method(chosen_method,params)
    call check(status,'primme_set_method')
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=dprimme(eval,evec(:,1:size(eval)),rnorms,params)
    call primme_collect_stats(params,nconv_local,nmat,niter_)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_d); active_mpi=.false.; active_comm=MPI_COMM_NULL
  end subroutine sp_eigh_primme_mpi_d

  subroutine sp_eigh_primme_mpi_z(MpiComm,MatVec,eval,evec,GlobalSumReal,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm)
    integer, intent(in)                             :: MpiComm
    procedure(primme_matvec_z), pointer, intent(in) :: MatVec
    real(c_double), intent(out)                     :: eval(:)
    complex(c_double), intent(inout)                :: evec(:,:)
    procedure(primme_global_sum_real), optional     :: GlobalSumReal
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,mpierr,chosen_method
    integer(c_int64_t)                              :: n,nlocal,nconv_local,nmat,niter_,ninit,nprocs,rank
    integer                                         :: ierr
    !
    !Set MPI global dimensions
    nlocal=int(size(evec,1),c_int64_t)
    call MPI_Allreduce(nlocal,n,1,MPI_INTEGER8,MPI_SUM,MpiComm,mpierr)
    if(mpierr/=MPI_SUCCESS)error stop 'sp_eigh_primme_mpi_z: MPI_Allreduce failed'
    call MPI_Comm_size(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_size')
    nprocs=int(ierr,c_int64_t)
    call MPI_Comm_rank(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_rank')
    rank=int(ierr,c_int64_t)
    if(size(evec,2)<size(eval))error stop 'sp_eigh_primme_mpi_z: evec has too few columns'
    !
    chosen_method=PRIMME_DYNAMIC; if(present(method)) chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize)) ninit=initSize
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t)) &
      error stop 'sp_eigh_primme_mpi_z: invalid initSize'
    !
    !Set the MatVec pointer
    user_matvec_z=>MatVec
    active_mpi=.true.; active_comm=MpiComm
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params))error stop 'sp_eigh_primme_mpi_z: primme_params_create failed'
    !
    !Set PRIMME internal parameters
    call primme_set(params,PRIMME_n,n,'PRIMME_n')
    call primme_set(params,PRIMME_nLocal,nlocal,'PRIMME_nLocal')
    call primme_set(params,PRIMME_numProcs,nprocs,'PRIMME_numProcs')
    call primme_set(params,PRIMME_procID,rank,'PRIMME_procID')
    call primme_set(params,PRIMME_numEvals,int(size(eval),c_int64_t),'PRIMME_numEvals')
    call primme_set(params,PRIMME_target,primme_smallest,'PRIMME_target')
    call primme_set(params,PRIMME_initSize,ninit,'PRIMME_initSize')
    if(present(Nitermax)) call primme_set(params,PRIMME_maxOuterIterations,Nitermax,'PRIMME_maxOuterIterations')
    if(present(Nblock))   call primme_set(params,PRIMME_maxBasisSize,Nblock,'PRIMME_maxBasisSize')
    if(present(tol))      call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !Set the internal PRIMME MatVec callback complex precision
    status=primme_set_member_matvec_complexdouble(params,PRIMME_matrixMatvec,mv_z)
    call check(status,'PRIMME_matrixMatvec')
    !
    !Set the internal PRIMME GlobalSumReal callback
    call primme_set(params,PRIMME_globalSumReal_type,primme_op_double,'PRIMME_globalSumReal_type')
    if(present(GlobalSumReal))then
      status=primme_set_member(params,PRIMME_globalSumReal,c_funloc(GlobalSumReal))
    else
      status=primme_set_member(params,PRIMME_globalSumReal,c_funloc(default_global_sum_real))
    endif
    call check(status,'PRIMME_globalSumReal')
    !
    !Set the internal PRIMME method
    status=primme_set_method(chosen_method,params)
    call check(status,'primme_set_method')
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=zprimme(eval,evec(:,1:size(eval)),rnorms,params)
    call primme_collect_stats(params,nconv_local,nmat,niter_)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_z); active_mpi=.false.; active_comm=MPI_COMM_NULL
  end subroutine sp_eigh_primme_mpi_z
#endif











!######################################################################
!           MATRIX-VECTOR PRODUCT CALLBACKS FOR PRIMME
!######################################################################
  subroutine mv_d(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    real(c_double)             :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    integer                    :: j
    ierr=0
    do j=1,int(blockSize)
      call user_matvec_d(int(ldx),x(1:ldx,j),y(1:ldy,j))
    enddo
  end subroutine mv_d

  subroutine mv_z(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    complex(c_double)          :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    integer                    :: j
    ierr=0
    do j=1,int(blockSize)
      call user_matvec_z(int(ldx),x(1:ldx,j),y(1:ldy,j))
    enddo
  end subroutine mv_z









  !######################################################################
  !           MPI DEFAULT GLOBAL_SUM_REAL + AUX  ROUTINES
  !######################################################################
#ifdef _MPI
  subroutine default_global_sum_real(sendBuf,recvBuf,count,primme,ierr) bind(c)
    type(c_ptr), value         :: sendBuf,recvBuf,primme
    integer(c_int), intent(in) :: count
    integer(c_int)             :: ierr
    real(c_double), pointer    :: send(:),recv(:)
    integer                    :: mpierr
    ierr=0
    call c_f_pointer(recvBuf,recv,[int(count)])
    if(c_associated(sendBuf,recvBuf))then
      call MPI_Allreduce(MPI_IN_PLACE,recv,int(count),MPI_DOUBLE_PRECISION,MPI_SUM,active_comm,mpierr)
    else
      call c_f_pointer(sendBuf,send,[int(count)])
      call MPI_Allreduce(send,recv,int(count),MPI_DOUBLE_PRECISION,MPI_SUM,active_comm,mpierr)
    endif
    ierr=int(mpierr,c_int)
  end subroutine default_global_sum_real

  subroutine check_mpi(status,where)
    integer,intent(in)      :: status
    character(*),intent(in) :: where
    if(status/=MPI_SUCCESS)then
      write(*,'(A,A)') 'SF_PRIMME MPI failure in ',trim(where)
      error stop
    endif
  end subroutine check_mpi
#endif




  !######################################################################
  !                        AUXILIARY ROUTINES
  !######################################################################
  subroutine primme_set_int(params,label,value,name)
    type(c_ptr),intent(in)        :: params
    integer(c_int),intent(in)     :: label
    integer(c_int64_t),intent(in) :: value
    character(*),intent(in)       :: name
    integer(c_int)                :: status
    status=primme_set_member(params,label,value); call check(status,name)
  end subroutine primme_set_int

  subroutine primme_set_double(params,label,value,name)
    type(c_ptr),intent(in)    :: params
    integer(c_int),intent(in) :: label
    real(c_double),intent(in) :: value
    character(*),intent(in)   :: name
    integer(c_int)            :: status
    status=primme_set_member(params,label,value); call check(status,name)
  end subroutine primme_set_double

  subroutine primme_collect_stats(params,nconv_local,nmat,niter)
    type(c_ptr),intent(in)         :: params
    integer(c_int64_t),intent(out) :: nconv_local,nmat,niter
    integer(c_int)                 :: status
    status=primme_get_member(params,PRIMME_initSize,nconv_local);           call check(status,'PRIMME_initSize')
    status=primme_get_member(params,PRIMME_stats_numMatvecs,nmat);          call check(status,'PRIMME_stats_numMatvecs')
    status=primme_get_member(params,PRIMME_stats_numOuterIterations,niter); call check(status,'PRIMME_stats_numOuterIterations')
  end subroutine primme_collect_stats

  subroutine check(status,where)
    integer(c_int),intent(in) :: status
    character(*),intent(in)   :: where
    if(status/=0)then
      write(*,'(A,A,A,I0)') 'SF_PRIMME: error setting ',trim(where),'; code=',status
      error stop
    endif
  end subroutine check




end module SF_PRIMME
