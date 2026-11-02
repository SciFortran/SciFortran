  ! Serial real symmetric eigensolver. MatVec applies H to one vector.
  ! Eigenvectors are returned as columns of evec; optional solver settings
  ! retain PRIMME defaults when omitted. Preconditioner is applied per vector.
  subroutine sp_eigh_primme_d(MatVec,eval,evec,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm,Preconditioner,NumPrec,TimePrec,NumMatVecCalls,MaxBlockSize)
    procedure(primme_matvec_d), pointer, intent(in) :: MatVec
    real(c_double), intent(out)                     :: eval(:)
    real(c_double), intent(inout)                   :: evec(:,:)
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    procedure(primme_preconditioner_d), optional    :: Preconditioner
    integer(c_int64_t), optional, intent(out)       :: NumPrec,NumMatVecCalls
    real(c_double), optional, intent(out)           :: TimePrec
    integer(c_int64_t), optional, intent(in)        :: MaxBlockSize
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,solve_status,chosen_method
    integer(c_int64_t)                              :: n,nconv_local,nmat,niter_,ninit,nprec
    real(c_double)                                  :: tprec
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
    nullify(user_matvec_d)
    user_matvec_d=>MatVec
    user_matvec_calls_d=0_c_int64_t
    nullify(user_preconditioner_d)
    if(present(Preconditioner))user_preconditioner_d=>Preconditioner
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
    if(present(MaxBlockSize))call primme_set(params,PRIMME_maxBlockSize,MaxBlockSize,'PRIMME_maxBlockSize')
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
    if(present(Preconditioner))then
      status=primme_set_member(params,PRIMME_applyPreconditioner,c_funloc(pre_d))
      call check(status,'PRIMME_applyPreconditioner')
      call primme_set(params,PRIMME_applyPreconditioner_type,primme_op_double,'PRIMME_applyPreconditioner_type')
      call primme_set(params,PRIMME_correctionParams_precondition,1_c_int64_t,'PRIMME_correctionParams_precondition')
    else
      call primme_set(params,PRIMME_correctionParams_precondition,0_c_int64_t,'PRIMME_correctionParams_precondition')
    endif
    !
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=dprimme(eval,evec(:,1:size(eval)),rnorms,params)
    !
    !
    call primme_collect_stats(params,nconv_local,nmat,niter_,nprec,tprec)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(NumPrec))NumPrec=nprec
    if(present(TimePrec))TimePrec=tprec
    if(present(NumMatVecCalls))NumMatVecCalls=user_matvec_calls_d
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme_d: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matvec_d); nullify(user_preconditioner_d)
  end subroutine sp_eigh_primme_d





  ! Serial real symmetric eigensolver with a MatMat callback. MatMat receives
  ! all local vectors in a PRIMME callback as columns of one matrix.
  subroutine sp_eigh_primme_block_d(eval,evec,MatMat,method,tol,Nitermax,Nblock,initSize,iexit,NumOp,Niter,Nconv,rnorm,Preconditioner,NumPrec,TimePrec,NumMatVecCalls,MaxBlockSize)
    procedure(primme_matmat_d), pointer, intent(in) :: MatMat
    real(c_double), intent(out)                     :: eval(:)
    real(c_double), intent(inout)                   :: evec(:,:)
    integer(c_int), optional, intent(in)            :: method
    real(c_double), optional, intent(in)            :: tol
    integer(c_int64_t), optional, intent(in)        :: Nitermax,Nblock,initSize
    integer(c_int), optional, intent(out)           :: iexit
    integer(c_int64_t), optional, intent(out)       :: NumOp,Niter,Nconv
    real(c_double), optional, intent(out)           :: rnorm(:)
    procedure(primme_preconditioner_d), optional    :: Preconditioner
    integer(c_int64_t), optional, intent(out)       :: NumPrec,NumMatVecCalls
    real(c_double), optional, intent(out)           :: TimePrec
    integer(c_int64_t), optional, intent(in)        :: MaxBlockSize
    type(c_ptr)                                     :: params
    real(c_double), allocatable                     :: rnorms(:)
    integer(c_int)                                  :: status,solve_status,chosen_method
    integer(c_int64_t)                              :: n,nconv_local,nmat,niter_,ninit,nprec
    real(c_double)                                  :: tprec
    integer                                         :: j
    !
    n=int(size(evec,1),c_int64_t)
    if(size(evec,2)<size(eval)) error stop 'sp_eigh_primme_block_d: evec has too few columns'
    !
    chosen_method=PRIMME_DYNAMIC; if(present(method))chosen_method=method
    ninit=0_c_int64_t;            if(present(initSize))ninit=initSize
    !
    if(ninit<0 .OR. ninit>int(size(eval),c_int64_t) .OR. ninit>int(size(evec,2),c_int64_t))&
      error stop 'sp_eigh_primme_block_d: invalid initSize'
    !
    !Set the MatMat pointer
    nullify(user_matmat_d)
    user_matmat_d=>MatMat
    user_matvec_calls_d=0_c_int64_t
    nullify(user_preconditioner_d)
    if(present(Preconditioner))user_preconditioner_d=>Preconditioner
    !
    active_mpi=.false.
    !
    !Initialize PRIMME
    params=primme_params_create()
    if(.not.c_associated(params)) error stop 'sp_eigh_primme_block_d: primme_params_create failed'
    !
    !
    !Set PRIMME internal parameters
    call primme_set(params,PRIMME_n,n,'PRIMME_n')
    call primme_set(params,PRIMME_numEvals,int(size(eval),c_int64_t),'PRIMME_numEvals')
    call primme_set(params,PRIMME_target,primme_smallest,'PRIMME_target')
    call primme_set(params,PRIMME_initSize,ninit,'PRIMME_initSize')
    if(present(Nitermax))call primme_set(params,PRIMME_maxOuterIterations,Nitermax,'PRIMME_maxOuterIterations')
    if(present(Nblock))call primme_set(params,PRIMME_maxBasisSize,Nblock,'PRIMME_maxBasisSize')
    if(present(MaxBlockSize))call primme_set(params,PRIMME_maxBlockSize,MaxBlockSize,'PRIMME_maxBlockSize')
    if(present(tol))call primme_set(params,PRIMME_eps,tol,'PRIMME_eps')
    !
    !
    !Set the internal PRIMME MatMat callback dble precision
    status=primme_set_member_matvec_double(params,PRIMME_matrixMatvec,mm_d)
    call check(status,'PRIMME_matrixMatvec')
    !
    !Set the internal PRIMME pre-condition method: PRIMME_DEFAULT, PRIMME_MIN_TIME, PRIMME_MIN_MATVECS
    status=primme_set_method(chosen_method,params)
    call check(status,'primme_set_method')
    if(present(Preconditioner))then
      status=primme_set_member(params,PRIMME_applyPreconditioner,c_funloc(pre_d))
      call check(status,'PRIMME_applyPreconditioner')
      call primme_set(params,PRIMME_applyPreconditioner_type,primme_op_double,'PRIMME_applyPreconditioner_type')
      call primme_set(params,PRIMME_correctionParams_precondition,1_c_int64_t,'PRIMME_correctionParams_precondition')
    else
      call primme_set(params,PRIMME_correctionParams_precondition,0_c_int64_t,'PRIMME_correctionParams_precondition')
    endif
    !
    !
    allocate(rnorms(size(eval))) ; rnorms=0
    !
    !CALL PRIMME:
    status=dprimme(eval,evec(:,1:size(eval)),rnorms,params)
    !
    !
    call primme_collect_stats(params,nconv_local,nmat,niter_,nprec,tprec)
    if(present(iexit))iexit=status
    if(present(NumOp))NumOp=nmat
    if(present(Niter))Niter=niter_
    if(present(Nconv))Nconv=nconv_local
    if(present(NumPrec))NumPrec=nprec
    if(present(TimePrec))TimePrec=tprec
    if(present(NumMatVecCalls))NumMatVecCalls=user_matvec_calls_d
    if(present(rnorm))then
      if(size(rnorm)<size(rnorms)) error stop 'sp_eigh_primme_block_d: rnorm array too small'
      rnorm(1:size(rnorms))=rnorms
    endif
    deallocate(rnorms)
    !
    !Destroy PRIMME
    status=primme_params_destroy(params)
    call check(status,'primme_params_destroy')
    nullify(user_matmat_d); nullify(user_preconditioner_d)
  end subroutine sp_eigh_primme_block_d




