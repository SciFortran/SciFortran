  ! MPI complex Hermitian eigensolver. MatVec receives the local vector
  ! segment; eval is replicated and evec stores local rows. GlobalSumReal
  ! defaults to the module MPI reduction.
  subroutine sp_eigh_primme_mpi_z(MpiComm,MatVec,eval,evec,GlobalSumReal,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm,Preconditioner,NumPrec,TimePrec,NumMatVecCalls,MaxBlockSize)
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
    procedure(primme_preconditioner_z), optional    :: Preconditioner
    integer(c_int64_t), optional, intent(out)       :: NumPrec,NumMatVecCalls
    real(c_double), optional, intent(out)           :: TimePrec
    integer(c_int64_t), optional, intent(in)        :: MaxBlockSize
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,mpierr,chosen_method
    integer(c_int64_t)                              :: n,nlocal,nconv_local,nmat,niter_,ninit,nprocs,rank,nprec
    real(c_double)                                  :: tprec
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
    nullify(user_matvec_z)
    user_matvec_z=>MatVec
    user_matvec_calls_z=0_c_int64_t
    nullify(user_preconditioner_z)
    if(present(Preconditioner))user_preconditioner_z=>Preconditioner
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
    if(present(MaxBlockSize))call primme_set(params,PRIMME_maxBlockSize,MaxBlockSize,'PRIMME_maxBlockSize')
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
    if(present(Preconditioner))then
      status=primme_set_member(params,PRIMME_applyPreconditioner,c_funloc(pre_z))
      call check(status,'PRIMME_applyPreconditioner')
      call primme_set(params,PRIMME_applyPreconditioner_type,primme_op_double,'PRIMME_applyPreconditioner_type')
      call primme_set(params,PRIMME_correctionParams_precondition,1_c_int64_t,'PRIMME_correctionParams_precondition')
    else
      call primme_set(params,PRIMME_correctionParams_precondition,0_c_int64_t,'PRIMME_correctionParams_precondition')
    endif
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=zprimme(eval,evec(:,1:size(eval)),rnorms,params)
    call primme_collect_stats(params,nconv_local,nmat,niter_,nprec,tprec)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(NumPrec))NumPrec=nprec
    if(present(TimePrec))TimePrec=tprec
    if(present(NumMatVecCalls))NumMatVecCalls=user_matvec_calls_z
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme_mpi_z: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_z); nullify(user_preconditioner_z)
    active_mpi=.false.; active_comm=MPI_COMM_NULL
  end subroutine sp_eigh_primme_mpi_z


  ! MPI complex Hermitian MatMat eigensolver. MatMat receives local vector
  ! segments as columns; eigenvalues are replicated and evec is row-local.
  subroutine sp_eigh_primme_block_mpi_z(MpiComm,eval,evec,MatMat,GlobalSumReal,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm,Preconditioner,NumPrec,TimePrec,NumMatVecCalls,MaxBlockSize)
    integer, intent(in)                             :: MpiComm
    procedure(primme_matmat_z), pointer, intent(in) :: MatMat
    real(c_double), intent(out)                     :: eval(:)
    complex(c_double), intent(inout)                :: evec(:,:)
    procedure(primme_global_sum_real), optional     :: GlobalSumReal
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    procedure(primme_preconditioner_z), optional    :: Preconditioner
    integer(c_int64_t), optional, intent(out)       :: NumPrec,NumMatVecCalls
    real(c_double), optional, intent(out)           :: TimePrec
    integer(c_int64_t), optional, intent(in)        :: MaxBlockSize
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,mpierr,chosen_method
    integer(c_int64_t)                              :: n,nlocal,nconv_local,nmat,niter_,ninit,nprocs,rank,nprec
    real(c_double)                                  :: tprec
    integer                                         :: ierr
    !
    !Set MPI global dimensions
    nlocal=int(size(evec,1),c_int64_t)
    call MPI_Allreduce(nlocal,n,1,MPI_INTEGER8,MPI_SUM,MpiComm,mpierr)
    if(mpierr/=MPI_SUCCESS)error stop 'sp_eigh_primme_block_mpi_z: MPI_Allreduce failed'
    call MPI_Comm_size(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_size')
    nprocs=int(ierr,c_int64_t)
    call MPI_Comm_rank(MpiComm,ierr,mpierr); call check_mpi(mpierr,'MPI_Comm_rank')
    rank=int(ierr,c_int64_t)
    if(size(evec,2)<size(eval))error stop 'sp_eigh_primme_block_mpi_z: evec has too few columns'
    !
    chosen_method=PRIMME_DYNAMIC; if(present(method)) chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize)) ninit=initSize
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t)) &
      error stop 'sp_eigh_primme_block_mpi_z: invalid initSize'
    !
    !Set the MatMat pointer
    nullify(user_matmat_z)
    user_matmat_z=>MatMat
    user_matvec_calls_z=0_c_int64_t
    nullify(user_preconditioner_z)
    if(present(Preconditioner))user_preconditioner_z=>Preconditioner
    active_mpi=.true.; active_comm=MpiComm
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params))error stop 'sp_eigh_primme_block_mpi_z: primme_params_create failed'
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
    if(present(MaxBlockSize))call primme_set(params,PRIMME_maxBlockSize,MaxBlockSize,'PRIMME_maxBlockSize')
    if(present(tol))      call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !Set the internal PRIMME MatMat callback complex precision
    status=primme_set_member_matvec_complexdouble(params,PRIMME_matrixMatvec,mm_z)
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
    if(present(Preconditioner))then
      status=primme_set_member(params,PRIMME_applyPreconditioner,c_funloc(pre_z))
      call check(status,'PRIMME_applyPreconditioner')
      call primme_set(params,PRIMME_applyPreconditioner_type,primme_op_double,'PRIMME_applyPreconditioner_type')
      call primme_set(params,PRIMME_correctionParams_precondition,1_c_int64_t,'PRIMME_correctionParams_precondition')
    else
      call primme_set(params,PRIMME_correctionParams_precondition,0_c_int64_t,'PRIMME_correctionParams_precondition')
    endif
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=zprimme(eval,evec(:,1:size(eval)),rnorms,params)
    call primme_collect_stats(params,nconv_local,nmat,niter_,nprec,tprec)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(NumPrec))NumPrec=nprec
    if(present(TimePrec))TimePrec=tprec
    if(present(NumMatVecCalls))NumMatVecCalls=user_matvec_calls_z
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme_block_mpi_z: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matmat_z); nullify(user_preconditioner_z)
    active_mpi=.false.; active_comm=MPI_COMM_NULL
  end subroutine sp_eigh_primme_block_mpi_z

