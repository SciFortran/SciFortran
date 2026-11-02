module SF_PRIMME
  USE ISO_C_BINDING, only: c_ptr,c_int,c_int64_t,c_double,c_funloc,c_f_pointer,c_associated,c_null_ptr
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


  !Callback contracts. MatVec receives one local vector; MatMat receives a
  !local multivector with one vector per column. Preconditioner receives one
  !correction vector and its associated Ritz shift.
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
    subroutine primme_matmat_d(n,nvec,vin,vout)
      import c_double
      integer, intent(in)         :: n,nvec
      real(c_double), intent(in)  :: vin(n,nvec)
      real(c_double), intent(out) :: vout(n,nvec)
    end subroutine primme_matmat_d
    !
    subroutine primme_matmat_z(n,nvec,vin,vout)
      import c_double
      integer, intent(in)            :: n,nvec
      complex(c_double), intent(in)  :: vin(n,nvec)
      complex(c_double), intent(out) :: vout(n,nvec)
    end subroutine primme_matmat_z
    !
    subroutine primme_preconditioner_d(n,vin,vout,shift)
      import c_double
      integer, intent(in)         :: n
      real(c_double), intent(in)  :: vin(n)
      real(c_double), intent(out) :: vout(n)
      real(c_double), intent(in)  :: shift
    end subroutine primme_preconditioner_d
    !
    subroutine primme_preconditioner_z(n,vin,vout,shift)
      import c_double
      integer, intent(in)            :: n
      complex(c_double), intent(in)  :: vin(n)
      complex(c_double), intent(out) :: vout(n)
      real(c_double), intent(in)     :: shift
    end subroutine primme_preconditioner_z
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
  procedure(primme_matvec_d), pointer, save         :: user_matvec_d => null()
  procedure(primme_matvec_z), pointer, save         :: user_matvec_z => null()
  procedure(primme_matmat_d), pointer, save         :: user_matmat_d => null()
  procedure(primme_matmat_z), pointer, save         :: user_matmat_z => null()
  procedure(primme_preconditioner_d), pointer, save :: user_preconditioner_d => null()
  procedure(primme_preconditioner_z), pointer, save :: user_preconditioner_z => null()
  integer(c_int64_t), save                          :: user_matvec_calls_d=0_c_int64_t
  integer(c_int64_t), save                          :: user_matvec_calls_z=0_c_int64_t
  !
  public :: primme_matvec_d
  public :: primme_matvec_z
  public :: primme_matmat_d
  public :: primme_matmat_z
  public :: primme_preconditioner_d
  public :: primme_preconditioner_z
  public :: primme_global_sum_real
  !
  !
  !
  !Public interfaces
  interface primme_eigsh
    module procedure sp_eigh_primme_d
    module procedure sp_eigh_primme_z
    module procedure sp_eigh_primme_block_d
    module procedure sp_eigh_primme_block_z
#ifdef _MPI
    module procedure sp_eigh_primme_mpi_d
    module procedure sp_eigh_primme_mpi_z
    module procedure sp_eigh_primme_block_mpi_d
    module procedure sp_eigh_primme_block_mpi_z
#endif
  end interface primme_eigsh
  !
  interface sp_eigh_primme
    module procedure sp_eigh_primme_d
    module procedure sp_eigh_primme_z
    module procedure sp_eigh_primme_block_d
    module procedure sp_eigh_primme_block_z
#ifdef _MPI
    module procedure sp_eigh_primme_mpi_d
    module procedure sp_eigh_primme_mpi_z
    module procedure sp_eigh_primme_block_mpi_d
    module procedure sp_eigh_primme_block_mpi_z
#endif
  end interface sp_eigh_primme
  !
  !
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

  !Precision-specific interfaces are kept in separate include files.
  include 'prime_eigsh_d.f90'
  include 'prime_eigsh_z.f90'
#ifdef _MPI
  include 'prime_eigsh_mpi_d.f90'
  include 'prime_eigsh_mpi_z.f90'
#endif


!######################################################################
!       SINGLE-VECTOR AND BLOCK MATRIX-PRODUCT CALLBACKS FOR PRIMME
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
      user_matvec_calls_d=user_matvec_calls_d+1_c_int64_t
      call user_matvec_d(int(ldx),x(1:ldx,j),y(1:ldy,j))
    enddo
  end subroutine mv_d

  subroutine mm_d(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    real(c_double)             :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    ierr=0
    user_matvec_calls_d=user_matvec_calls_d+1_c_int64_t
    call user_matmat_d(int(ldx),int(blockSize),x(1:ldx,1:blockSize),y(1:ldy,1:blockSize))
  end subroutine mm_d

  subroutine mv_z(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    complex(c_double)          :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    integer                    :: j
    ierr=0
    do j=1,int(blockSize)
      user_matvec_calls_z=user_matvec_calls_z+1_c_int64_t
      call user_matvec_z(int(ldx),x(1:ldx,j),y(1:ldy,j))
    enddo
  end subroutine mv_z

  subroutine mm_z(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    complex(c_double)          :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    ierr=0
    user_matvec_calls_z=user_matvec_calls_z+1_c_int64_t
    call user_matmat_z(int(ldx),int(blockSize),x(1:ldx,1:blockSize),y(1:ldy,1:blockSize))
  end subroutine mm_z






  !######################################################################
  !                    PRE-CONDITIONER CALLBACK
  !######################################################################
  subroutine pre_d(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    real(c_double)             :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    type(c_ptr)                :: shifts_ptr
    real(c_double), pointer    :: shifts(:)
    integer                    :: j
    ierr=0
    shifts_ptr=c_null_ptr
    ierr=primme_get_member(primme,PRIMME_ShiftsForPreconditioner,shifts_ptr)
    if(ierr/=0)return
    if(.not.c_associated(shifts_ptr))then
      ierr=1
      return
    endif
    call c_f_pointer(shifts_ptr,shifts,[int(blockSize)])
    do j=1,int(blockSize)
      call user_preconditioner_d(int(ldx),x(1:ldx,j),y(1:ldy,j),shifts(j))
    enddo
  end subroutine pre_d

  subroutine pre_z(x,ldx,y,ldy,blockSize,primme,ierr) bind(c)
    integer(c_int64_t)         :: ldx,ldy
    complex(c_double)          :: x(ldx,*),y(ldy,*)
    integer(c_int), intent(in) :: blockSize
    type(c_ptr), value         :: primme
    integer(c_int)             :: ierr
    type(c_ptr)                :: shifts_ptr
    real(c_double), pointer    :: shifts(:)
    integer                    :: j
    ierr=0
    shifts_ptr=c_null_ptr
    ierr=primme_get_member(primme,PRIMME_ShiftsForPreconditioner,shifts_ptr)
    if(ierr/=0)return
    if(.not.c_associated(shifts_ptr))then
      ierr=1
      return
    endif
    call c_f_pointer(shifts_ptr,shifts,[int(blockSize)])
    do j=1,int(blockSize)
      call user_preconditioner_z(int(ldx),x(1:ldx,j),y(1:ldy,j),shifts(j))
    enddo
  end subroutine pre_z










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

  subroutine primme_collect_stats(params,nconv_local,nmat,niter,nprec,tprec)
    type(c_ptr),intent(in)         :: params
    integer(c_int64_t),intent(out) :: nconv_local,nmat,niter
    integer(c_int64_t),intent(out),optional :: nprec
    real(c_double),intent(out),optional     :: tprec
    integer(c_int)                 :: status
    status=primme_get_member(params,PRIMME_initSize,nconv_local);
    call check(status,'PRIMME_initSize')
    status=primme_get_member(params,PRIMME_stats_numMatvecs,nmat);
    call check(status,'PRIMME_stats_numMatvecs')
    status=primme_get_member(params,PRIMME_stats_numOuterIterations,niter);
    call check(status,'PRIMME_stats_numOuterIterations')
    if(present(nprec))then
      status=primme_get_member(params,PRIMME_stats_numPreconds,nprec);
      call check(status,'PRIMME_stats_numPreconds')
    endif
    if(present(tprec))then
      status=primme_get_member(params,PRIMME_stats_timePrecond,tprec);
      call check(status,'PRIMME_stats_timePrecond')
    endif
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
