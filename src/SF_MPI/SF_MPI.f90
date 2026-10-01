MODULE SF_MPI
  !SciFortran MPI interface
  implicit none
#ifdef _MPI
  ! USE MPI
  include 'mpif.h'
#endif


  private


#ifdef _MPI
  interface Bcast_MPI
  !This subroutine broadcasts the logical, integer, real(8) or complex(8) scalar or array :f:var:`data` from the process
  !:f:var:`root` to all the processes of the communicator :f:var:`comm`, using :code:`MPI_BCAST`. The specific procedures cover
  !scalars and arrays up to rank 7 (8 with gfortran newer than 8). Despite the :code:`intent(in)`, :f:var:`data` is overwritten
  !on the processes other than :f:var:`root`. Nothing is done if :f:var:`comm` is :code:`MPI_COMM_NULL`. MPI errors are reported
  !by :f:func_inline:`Error_MPI`.
  !
     module procedure :: MPI_Bcast_Bool_0
     module procedure :: MPI_Bcast_Bool_1
     module procedure :: MPI_Bcast_Bool_2
     module procedure :: MPI_Bcast_Bool_3
     module procedure :: MPI_Bcast_Bool_4
     module procedure :: MPI_Bcast_Bool_5
     module procedure :: MPI_Bcast_Bool_6
     module procedure :: MPI_Bcast_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_Bcast_Bool_8
#endif
     !
     module procedure :: MPI_Bcast_Int_0
     module procedure :: MPI_Bcast_Int_1
     module procedure :: MPI_Bcast_Int_2
     module procedure :: MPI_Bcast_Int_3
     module procedure :: MPI_Bcast_Int_4
     module procedure :: MPI_Bcast_Int_5
     module procedure :: MPI_Bcast_Int_6
     module procedure :: MPI_Bcast_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_Bcast_Int_8
#endif
     !
     module procedure :: MPI_Bcast_Dble_0
     module procedure :: MPI_Bcast_Dble_1
     module procedure :: MPI_Bcast_Dble_2
     module procedure :: MPI_Bcast_Dble_3
     module procedure :: MPI_Bcast_Dble_4
     module procedure :: MPI_Bcast_Dble_5
     module procedure :: MPI_Bcast_Dble_6
     module procedure :: MPI_Bcast_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_Bcast_Dble_8
#endif
     !
     module procedure :: MPI_Bcast_Cmplx_0
     module procedure :: MPI_Bcast_Cmplx_1
     module procedure :: MPI_Bcast_Cmplx_2
     module procedure :: MPI_Bcast_Cmplx_3
     module procedure :: MPI_Bcast_Cmplx_4
     module procedure :: MPI_Bcast_Cmplx_5
     module procedure :: MPI_Bcast_Cmplx_6
     module procedure :: MPI_Bcast_Cmplx_7     
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_Bcast_Cmplx_8
#endif
  end interface Bcast_MPI



  interface AllGather_MPI
  !This subroutine gathers the logical, integer, real(8) or complex(8) scalar or array :f:var:`send` of every process of the
  !communicator :f:var:`comm` into :f:var:`data`, on all the processes, in the order of the ranks, using :code:`MPI_ALLGATHER`.
  !The specific procedures cover scalars and arrays up to rank 7 (8 with gfortran newer than 8). :f:var:`data` has to be large
  !enough to hold the contributions of all the processes. Nothing is done if :f:var:`comm` is :code:`MPI_COMM_NULL`.
  !
     module procedure :: MPI_AllGather_Bool_0
     module procedure :: MPI_AllGather_Bool_1
     module procedure :: MPI_AllGather_Bool_2
     module procedure :: MPI_AllGather_Bool_3
     module procedure :: MPI_AllGather_Bool_4
     module procedure :: MPI_AllGather_Bool_5
     module procedure :: MPI_AllGather_Bool_6
     module procedure :: MPI_AllGather_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllGather_Bool_8
#endif
     !
     module procedure :: MPI_AllGather_Int_0
     module procedure :: MPI_AllGather_Int_1
     module procedure :: MPI_AllGather_Int_2
     module procedure :: MPI_AllGather_Int_3
     module procedure :: MPI_AllGather_Int_4
     module procedure :: MPI_AllGather_Int_5
     module procedure :: MPI_AllGather_Int_6
     module procedure :: MPI_AllGather_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllGather_Int_8
#endif
     !
     module procedure :: MPI_AllGather_Dble_0
     module procedure :: MPI_AllGather_Dble_1
     module procedure :: MPI_AllGather_Dble_2
     module procedure :: MPI_AllGather_Dble_3
     module procedure :: MPI_AllGather_Dble_4
     module procedure :: MPI_AllGather_Dble_5
     module procedure :: MPI_AllGather_Dble_6
     module procedure :: MPI_AllGather_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllGather_Dble_8
#endif
     !
     module procedure :: MPI_AllGather_Cmplx_0
     module procedure :: MPI_AllGather_Cmplx_1
     module procedure :: MPI_AllGather_Cmplx_2
     module procedure :: MPI_AllGather_Cmplx_3
     module procedure :: MPI_AllGather_Cmplx_4
     module procedure :: MPI_AllGather_Cmplx_5
     module procedure :: MPI_AllGather_Cmplx_6
     module procedure :: MPI_AllGather_Cmplx_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllGather_Cmplx_8
#endif
  end interface AllGather_MPI




  interface AllReduce_MPI
  !This subroutine sums the logical, integer, real(8) or complex(8) scalar or array :f:var:`send` over all the processes of the
  !communicator :f:var:`comm`, and returns the sum in :f:var:`data` on every process, using :code:`MPI_ALLREDUCE` with
  !:code:`MPI_SUM`. The specific procedures cover scalars and arrays up to rank 7 (8 with gfortran newer than 8). Nothing is
  !done if :f:var:`comm` is :code:`MPI_COMM_NULL`. For the maximum and the sum on a single process see :f:func_inline:`Max_MPI`
  !and :f:func_inline:`Sum_MPI`. Note that the MPI standard does not define :code:`MPI_SUM` for logical data, so the logical
  !procedures may fail or depend on the MPI implementation.
  !
     module procedure :: MPI_AllReduce_Bool_0
     module procedure :: MPI_AllReduce_Bool_1
     module procedure :: MPI_AllReduce_Bool_2
     module procedure :: MPI_AllReduce_Bool_3
     module procedure :: MPI_AllReduce_Bool_4
     module procedure :: MPI_AllReduce_Bool_5
     module procedure :: MPI_AllReduce_Bool_6
     module procedure :: MPI_AllReduce_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllReduce_Bool_8
#endif
     !
     module procedure :: MPI_AllReduce_Int_0
     module procedure :: MPI_AllReduce_Int_1
     module procedure :: MPI_AllReduce_Int_2
     module procedure :: MPI_AllReduce_Int_3
     module procedure :: MPI_AllReduce_Int_4
     module procedure :: MPI_AllReduce_Int_5
     module procedure :: MPI_AllReduce_Int_6
     module procedure :: MPI_AllReduce_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllReduce_Int_8
#endif
     !
     module procedure :: MPI_AllReduce_Dble_0
     module procedure :: MPI_AllReduce_Dble_1
     module procedure :: MPI_AllReduce_Dble_2
     module procedure :: MPI_AllReduce_Dble_3
     module procedure :: MPI_AllReduce_Dble_4
     module procedure :: MPI_AllReduce_Dble_5
     module procedure :: MPI_AllReduce_Dble_6
     module procedure :: MPI_AllReduce_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllReduce_Dble_8
#endif
     !
     module procedure :: MPI_AllReduce_Cmplx_0
     module procedure :: MPI_AllReduce_Cmplx_1
     module procedure :: MPI_AllReduce_Cmplx_2
     module procedure :: MPI_AllReduce_Cmplx_3
     module procedure :: MPI_AllReduce_Cmplx_4
     module procedure :: MPI_AllReduce_Cmplx_5
     module procedure :: MPI_AllReduce_Cmplx_6
     module procedure :: MPI_AllReduce_Cmplx_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_AllReduce_Cmplx_8
#endif
  end interface AllReduce_MPI




  interface Reduce_MPI
  !This subroutine is the same as :f:func_inline:`Sum_MPI`: it sums the logical, integer, real(8) or complex(8) scalar or array
  !:f:var:`data` over all the processes of the communicator :f:var:`comm` on the process :f:var:`root` (default 0), using
  !:code:`MPI_REDUCE` with :code:`MPI_SUM`. The specific procedures cover scalars and arrays up to rank 7 (8 with gfortran newer
  !than 8).
  !
  !If :f:var:`recv` is present the sum is returned in :f:var:`recv` and :f:var:`data` is left unchanged. Otherwise the sum
  !overwrites :f:var:`data` on the process of rank 0, which is the root for the default :f:var:`root`, and :f:var:`data` is
  !unchanged on the other processes. Nothing is done if :f:var:`comm` is :code:`MPI_COMM_NULL`. Note that the MPI standard does
  !not define :code:`MPI_SUM` for logical data, so the logical procedures may fail or depend on the MPI implementation.
  !
     module procedure :: MPI_ReduceSum_Bool_0
     module procedure :: MPI_ReduceSum_Bool_1
     module procedure :: MPI_ReduceSum_Bool_2
     module procedure :: MPI_ReduceSum_Bool_3
     module procedure :: MPI_ReduceSum_Bool_4
     module procedure :: MPI_ReduceSum_Bool_5
     module procedure :: MPI_ReduceSum_Bool_6
     module procedure :: MPI_ReduceSum_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Bool_8
#endif
     !
     module procedure :: MPI_ReduceSum_Int_0
     module procedure :: MPI_ReduceSum_Int_1
     module procedure :: MPI_ReduceSum_Int_2
     module procedure :: MPI_ReduceSum_Int_3
     module procedure :: MPI_ReduceSum_Int_4
     module procedure :: MPI_ReduceSum_Int_5
     module procedure :: MPI_ReduceSum_Int_6
     module procedure :: MPI_ReduceSum_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Int_8
#endif
     !
     module procedure :: MPI_ReduceSum_Dble_0
     module procedure :: MPI_ReduceSum_Dble_1
     module procedure :: MPI_ReduceSum_Dble_2
     module procedure :: MPI_ReduceSum_Dble_3
     module procedure :: MPI_ReduceSum_Dble_4
     module procedure :: MPI_ReduceSum_Dble_5
     module procedure :: MPI_ReduceSum_Dble_6
     module procedure :: MPI_ReduceSum_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Dble_8
#endif
     !
     module procedure :: MPI_ReduceSum_Cmplx_0
     module procedure :: MPI_ReduceSum_Cmplx_1
     module procedure :: MPI_ReduceSum_Cmplx_2
     module procedure :: MPI_ReduceSum_Cmplx_3
     module procedure :: MPI_ReduceSum_Cmplx_4
     module procedure :: MPI_ReduceSum_Cmplx_5
     module procedure :: MPI_ReduceSum_Cmplx_6
     module procedure :: MPI_ReduceSum_Cmplx_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Cmplx_8
#endif
  end interface Reduce_MPI





  interface Sum_MPI
  !This subroutine sums the logical, integer, real(8) or complex(8) scalar or array :f:var:`data` over all the processes of the
  !communicator :f:var:`comm` on the process :f:var:`root` (default 0), using :code:`MPI_REDUCE` with :code:`MPI_SUM`. The
  !specific procedures cover scalars and arrays up to rank 7 (8 with gfortran newer than 8). It is the same as
  !:f:func_inline:`Reduce_MPI`.
  !
  !If :f:var:`recv` is present the sum is returned in :f:var:`recv` and :f:var:`data` is left unchanged. Otherwise the sum
  !overwrites :f:var:`data` on the process of rank 0, which is the root for the default :f:var:`root`, and :f:var:`data` is
  !unchanged on the other processes. Nothing is done if :f:var:`comm` is :code:`MPI_COMM_NULL`. The sum on every process is
  !available with :f:func_inline:`AllReduce_MPI`. Note that the MPI standard does not define :code:`MPI_SUM` for logical data,
  !so the logical procedures may fail or depend on the MPI implementation.
  !
     module procedure :: MPI_ReduceSum_Bool_0
     module procedure :: MPI_ReduceSum_Bool_1
     module procedure :: MPI_ReduceSum_Bool_2
     module procedure :: MPI_ReduceSum_Bool_3
     module procedure :: MPI_ReduceSum_Bool_4
     module procedure :: MPI_ReduceSum_Bool_5
     module procedure :: MPI_ReduceSum_Bool_6
     module procedure :: MPI_ReduceSum_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Bool_8
#endif
     !
     module procedure :: MPI_ReduceSum_Int_0
     module procedure :: MPI_ReduceSum_Int_1
     module procedure :: MPI_ReduceSum_Int_2
     module procedure :: MPI_ReduceSum_Int_3
     module procedure :: MPI_ReduceSum_Int_4
     module procedure :: MPI_ReduceSum_Int_5
     module procedure :: MPI_ReduceSum_Int_6
     module procedure :: MPI_ReduceSum_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Int_8
#endif
     !
     module procedure :: MPI_ReduceSum_Dble_0
     module procedure :: MPI_ReduceSum_Dble_1
     module procedure :: MPI_ReduceSum_Dble_2
     module procedure :: MPI_ReduceSum_Dble_3
     module procedure :: MPI_ReduceSum_Dble_4
     module procedure :: MPI_ReduceSum_Dble_5
     module procedure :: MPI_ReduceSum_Dble_6
     module procedure :: MPI_ReduceSum_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Dble_8
#endif
     !
     module procedure :: MPI_ReduceSum_Cmplx_0
     module procedure :: MPI_ReduceSum_Cmplx_1
     module procedure :: MPI_ReduceSum_Cmplx_2
     module procedure :: MPI_ReduceSum_Cmplx_3
     module procedure :: MPI_ReduceSum_Cmplx_4
     module procedure :: MPI_ReduceSum_Cmplx_5
     module procedure :: MPI_ReduceSum_Cmplx_6
     module procedure :: MPI_ReduceSum_Cmplx_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceSum_Cmplx_8
#endif
  end interface Sum_MPI



  interface Max_MPI
  !This subroutine evaluates the element-wise maximum of the logical, integer, real(8) or complex(8) scalar or array
  !:f:var:`data` over all the processes of the communicator :f:var:`comm` on the process :f:var:`root` (default 0), using
  !:code:`MPI_REDUCE` with :code:`MPI_MAX`. The specific procedures cover scalars and arrays up to rank 7 (8 with gfortran newer
  !than 8).
  !
  !If :f:var:`recv` is present the maximum is returned in :f:var:`recv` and :f:var:`data` is left unchanged. Otherwise the
  !maximum overwrites :f:var:`data` on the process of rank 0, which is the root for the default :f:var:`root`, and :f:var:`data`
  !is unchanged on the other processes. Nothing is done if :f:var:`comm` is :code:`MPI_COMM_NULL`. Note that the MPI standard
  !does not define :code:`MPI_MAX` for logical and complex data, so these procedures may fail or depend on the MPI
  !implementation.
  !
     module procedure :: MPI_ReduceMax_Bool_0
     module procedure :: MPI_ReduceMax_Bool_1
     module procedure :: MPI_ReduceMax_Bool_2
     module procedure :: MPI_ReduceMax_Bool_3
     module procedure :: MPI_ReduceMax_Bool_4
     module procedure :: MPI_ReduceMax_Bool_5
     module procedure :: MPI_ReduceMax_Bool_6
     module procedure :: MPI_ReduceMax_Bool_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceMax_Bool_8
#endif
     !
     module procedure :: MPI_ReduceMax_Int_0
     module procedure :: MPI_ReduceMax_Int_1
     module procedure :: MPI_ReduceMax_Int_2
     module procedure :: MPI_ReduceMax_Int_3
     module procedure :: MPI_ReduceMax_Int_4
     module procedure :: MPI_ReduceMax_Int_5
     module procedure :: MPI_ReduceMax_Int_6
     module procedure :: MPI_ReduceMax_Int_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceMax_Int_8
#endif
     !
     module procedure :: MPI_ReduceMax_Dble_0
     module procedure :: MPI_ReduceMax_Dble_1
     module procedure :: MPI_ReduceMax_Dble_2
     module procedure :: MPI_ReduceMax_Dble_3
     module procedure :: MPI_ReduceMax_Dble_4
     module procedure :: MPI_ReduceMax_Dble_5
     module procedure :: MPI_ReduceMax_Dble_6
     module procedure :: MPI_ReduceMax_Dble_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceMax_Dble_8
#endif
     !
     module procedure :: MPI_ReduceMax_Cmplx_0
     module procedure :: MPI_ReduceMax_Cmplx_1
     module procedure :: MPI_ReduceMax_Cmplx_2
     module procedure :: MPI_ReduceMax_Cmplx_3
     module procedure :: MPI_ReduceMax_Cmplx_4
     module procedure :: MPI_ReduceMax_Cmplx_5
     module procedure :: MPI_ReduceMax_Cmplx_6
     module procedure :: MPI_ReduceMax_Cmplx_7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: MPI_ReduceMax_Cmplx_8
#endif
  end interface Max_MPI





  public :: Init_MPI
  public :: Finalize_MPI
  public :: StartMsg_MPI
  public :: Barrier_MPI
  public :: Check_MPI
  !
  public :: Get_Size_MPI
  public :: Get_Rank_MPI
  public :: Get_Master_MPI
  public :: Get_Last_MPI
  public :: cpu_time_MPI
  public :: Error_MPI
  !
  public :: Bcast_MPI
  public :: Reduce_MPI
  public :: AllGather_MPI
  public :: AllReduce_MPI
  public :: Sum_MPI
  public :: Max_MPI


  integer :: ierr
  integer :: rank
  logical :: master
  integer :: root_

contains


  !****************************************
  !              MPI START/STOP
  !****************************************
  subroutine Init_MPI(comm,msg)
    !This subroutine initializes the MPI environment with :code:`MPI_Init`. If :f:var:`comm` is present it is set to
    !:code:`MPI_COMM_WORLD`, and if :f:var:`msg` is present and true the start message of :f:func_inline:`StartMsg_MPI` is
    !printed. MPI errors are reported by :f:func_inline:`Error_MPI`, which does not stop the program.
    !
    integer,optional :: comm
    logical,optional :: msg  !if present and T print the start message, see StartMsg_MPI
    call MPI_Init(ierr)
    call Error_MPI(ierr,"MPI_Start")
    if(present(comm))comm=MPI_COMM_WORLD
    if(present(msg))then
       if(msg)call StartMsg_MPI(MPI_COMM_WORLD)
    endif
  end subroutine Init_MPI

  subroutine Finalize_MPI(comm)
    !This subroutine terminates the MPI environment with :code:`MPI_Finalize`. If :f:var:`comm` is present it is set to
    !:code:`MPI_COMM_NULL`. MPI errors are reported by :f:func_inline:`Error_MPI`, which does not stop the program.
    !
    integer,optional :: comm
    call MPI_Finalize(ierr)
    call Error_MPI(ierr,"MPI_Stop")
    if(present(comm))comm=MPI_COMM_NULL
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    !This subroutine prints on the standard output a start message made of a header line, one line :code:`rank: i of N alive`
    !for each process of the communicator :f:var:`comm` in the order of the ranks, and a closing line. The printing is
    !synchronized with barriers. The default communicator is :code:`MPI_COMM_WORLD`, and nothing is printed for
    !:code:`MPI_COMM_NULL`.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_,size
    integer          :: i
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    if(comm_ /= Mpi_Comm_Null)then
       rank = Get_Rank_MPI(comm_)
       size = Get_Size_MPI(comm_)
       if(rank==0)write(*,'(a)')"---------------MPI----------------"
       do i=0,size-1
          if(rank==i)write(*,"(A,I6,A,I6,A)")"rank:",rank," of ",size," alive"
          call MPI_Barrier(comm_,ierr)
       enddo
       call MPI_Barrier(comm_,ierr)
       if(rank==0)write(*,'(a)')"----------------------------------"
       if(rank==0)write(*,'(a)')""
    endif
  end subroutine StartMsg_MPI

  subroutine Barrier_MPI(comm)
    !This subroutine blocks the calling process until all the processes of the communicator :f:var:`comm` have reached the
    !barrier, using :code:`MPI_Barrier`. The default communicator is :code:`MPI_COMM_WORLD`, and nothing is done for
    !:code:`MPI_COMM_NULL`.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    if(comm_/=Mpi_Comm_Null)then
       call MPI_Barrier(comm_,ierr)
       call Error_MPI(ierr,"Barrier_MPI")
    endif
  end subroutine Barrier_MPI


  !****************************************
  !              MPI TOOLS
  !****************************************
  function check_MPI() result(bool)
    !This function returns :code:`.true.` if the MPI environment has been initialized, as given by :code:`MPI_Initialized`.
    !
    logical          :: bool    
    call MPI_Initialized(bool,ierr)
  end function check_MPI


  function get_size_MPI(comm) result(size)
    !This function returns the number of processes of the communicator :f:var:`comm` (default :code:`MPI_COMM_WORLD`), using
    !:code:`MPI_Comm_size`. For :code:`MPI_COMM_NULL` the result is not set.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_
    integer          :: size
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm    
    if(comm_/=Mpi_Comm_Null)then
       call MPI_Comm_size(comm_,size,ierr)
       call Error_MPI(ierr,"Get_Size_MPI")
    else
       return
    endif
  end function get_size_MPI

  function Get_rank_MPI(comm) result(rank)
    !This function returns the rank of the calling process in the communicator :f:var:`comm` (default :code:`MPI_COMM_WORLD`),
    !using :code:`MPI_Comm_rank`. For :code:`MPI_COMM_NULL` the result is not set.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_
    integer          :: rank
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    if(comm_/=Mpi_Comm_Null)then
       call MPI_Comm_rank(comm_,rank,ierr)
       call Error_MPI(ierr,"Get_Rank_MPI")
    else
       return
    endif
  end function Get_rank_MPI

  function Get_master_MPI(comm) result(master)
    !This function returns :code:`.true.` on the master process, i.e. the process of rank 0, of the communicator :f:var:`comm`
    !(default :code:`MPI_COMM_WORLD`), and :code:`.false.` on the other processes and for :code:`MPI_COMM_NULL`.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_
    integer          :: rank
    logical          :: master
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    if(comm_/=Mpi_Comm_Null)then    
       call MPI_Comm_rank(comm_,rank,ierr)
       call Error_MPI(ierr,"Get_Master_MPI")
       master=.false.
       if(rank==0)master=.true.
    else
       master=.false.
    endif
  end function Get_master_MPI

  function Get_last_MPI(comm) result(last)
    !This function returns :code:`.true.` on the last process, i.e. the process of rank :code:`size-1`, of the communicator
    !:f:var:`comm` (default :code:`MPI_COMM_WORLD`), and :code:`.false.` on the other processes and for :code:`MPI_COMM_NULL`.
    !
    integer,optional :: comm  !MPI communicator (default: MPI_COMM_WORLD)
    integer          :: comm_
    integer          :: size
    integer          :: rank
    logical          :: last
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    if(comm_/=Mpi_Comm_Null)then
       call MPI_Comm_rank(comm_,rank,ierr)
       call MPI_Comm_size(comm_,size,ierr)    
       last=.false.
       if(rank==size-1)last=.true.
    else
       last=.false.
    endif
  end function Get_last_MPI


  !returns an elapsed time on the calling processor
  function cpu_time_MPI() result(time)
    !This function returns the elapsed time in seconds on the calling processor, as given by :code:`MPI_WTIME`. The time is
    !measured from an arbitrary point in the past, so only the difference between two calls is meaningful.
    !
    real(8) :: time
    time = MPI_WTIME()
  end function Cpu_Time_MPI


  function Get_Processor_MPI() result(workstation)
    integer                               :: istat
    character(len=MPI_MAX_PROCESSOR_NAME) :: workstation
    call MPI_GET_PROCESSOR_NAME(workstation,istat,ierr)
    call Error_MPI(ierr,"Get_Processor_MPI")
  end function Get_Processor_MPI













  !****************************************
  !              MPI BROADCAST
  !****************************************
  !!Bool
  subroutine MPI_Bcast_Bool_0(comm,data,root)
    integer,intent(in)          :: comm  !MPI communicator
    logical,intent(in)          :: data  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root  !rank of the root process (default 0)
    logical,dimension(1)        :: data_
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_0')
  end subroutine MPI_Bcast_Bool_0
  !
  subroutine MPI_Bcast_Bool_1(comm,data,root)
    integer,intent(in)          :: comm     !MPI communicator
    logical,intent(in)          :: data(:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_1')
  end subroutine MPI_Bcast_Bool_1
  !
  subroutine MPI_Bcast_Bool_2(comm,data,root)
    integer,intent(in)          :: comm       !MPI communicator
    logical,intent(in)          :: data(:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_2')
  end subroutine MPI_Bcast_Bool_2
  !
  subroutine MPI_Bcast_Bool_3(comm,data,root)
    integer,intent(in)          :: comm         !MPI communicator
    logical,intent(in)          :: data(:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_3')
  end subroutine MPI_Bcast_Bool_3
  !
  subroutine MPI_Bcast_Bool_4(comm,data,root)
    integer,intent(in)          :: comm           !MPI communicator
    logical,intent(in)          :: data(:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_4')
  end subroutine MPI_Bcast_Bool_4
  !
  subroutine MPI_Bcast_Bool_5(comm,data,root)
    integer,intent(in)          :: comm             !MPI communicator
    logical,intent(in)          :: data(:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_5')
  end subroutine MPI_Bcast_Bool_5
  !
  subroutine MPI_Bcast_Bool_6(comm,data,root)
    integer,intent(in)          :: comm               !MPI communicator
    logical,intent(in)          :: data(:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_6')
  end subroutine MPI_Bcast_Bool_6
  !
  subroutine MPI_Bcast_Bool_7(comm,data,root)
    integer,intent(in)          :: comm                 !MPI communicator
    logical,intent(in)          :: data(:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_7')
  end subroutine MPI_Bcast_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Bool_8(comm,data,root)
    integer,intent(in)          :: comm                   !MPI communicator
    logical,intent(in)          :: data(:,:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_8')
  end subroutine MPI_Bcast_Bool_8
#endif




  !! INTEGER
  subroutine MPI_Bcast_Int_0(comm,data,root)
    integer,intent(in)          :: comm  !MPI communicator
    integer,intent(in)          :: data  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root  !rank of the root process (default 0)
    integer,dimension(1)        :: data_
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_0')
  end subroutine MPI_Bcast_Int_0
  !
  subroutine MPI_Bcast_Int_1(comm,data,root)
    integer,intent(in)          :: comm     !MPI communicator
    integer,intent(in)          :: data(:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_1')
  end subroutine MPI_Bcast_Int_1
  !
  subroutine MPI_Bcast_Int_2(comm,data,root)
    integer,intent(in)          :: comm       !MPI communicator
    integer,intent(in)          :: data(:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_2')
  end subroutine MPI_Bcast_Int_2
  !
  subroutine MPI_Bcast_Int_3(comm,data,root)
    integer,intent(in)          :: comm         !MPI communicator
    integer,intent(in)          :: data(:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_3')
  end subroutine MPI_Bcast_Int_3
  !
  subroutine MPI_Bcast_Int_4(comm,data,root)
    integer,intent(in)          :: comm           !MPI communicator
    integer,intent(in)          :: data(:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_4')
  end subroutine MPI_Bcast_Int_4
  !
  subroutine MPI_Bcast_Int_5(comm,data,root)
    integer,intent(in)          :: comm             !MPI communicator
    integer,intent(in)          :: data(:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_5')
  end subroutine MPI_Bcast_Int_5
  !
  subroutine MPI_Bcast_Int_6(comm,data,root)
    integer,intent(in)          :: comm               !MPI communicator
    integer,intent(in)          :: data(:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_6')
  end subroutine MPI_Bcast_Int_6
  !
  subroutine MPI_Bcast_Int_7(comm,data,root)
    integer,intent(in)          :: comm                 !MPI communicator
    integer,intent(in)          :: data(:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_7')
  end subroutine MPI_Bcast_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Int_8(comm,data,root)
    integer,intent(in)          :: comm                   !MPI communicator
    integer,intent(in)          :: data(:,:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_7')
  end subroutine MPI_Bcast_Int_8
#endif



  !! REAL8
  subroutine MPI_Bcast_Dble_0(comm,data,root)
    integer,intent(in)          :: comm  !MPI communicator
    real(8),intent(in)          :: data  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root  !rank of the root process (default 0)
    real(8),dimension(1)        :: data_
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_0')
  end subroutine MPI_Bcast_Dble_0
  !
  subroutine MPI_Bcast_Dble_1(comm,data,root)
    integer,intent(in)          :: comm     !MPI communicator
    real(8),intent(in)          :: data(:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_1')
  end subroutine MPI_Bcast_Dble_1
  !
  subroutine MPI_Bcast_Dble_2(comm,data,root)
    integer,intent(in)          :: comm       !MPI communicator
    real(8),intent(in)          :: data(:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_2')
  end subroutine MPI_Bcast_Dble_2
  !
  subroutine MPI_Bcast_Dble_3(comm,data,root)
    integer,intent(in)          :: comm         !MPI communicator
    real(8),intent(in)          :: data(:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_3')
  end subroutine MPI_Bcast_Dble_3
  !
  subroutine MPI_Bcast_Dble_4(comm,data,root)
    integer,intent(in)          :: comm           !MPI communicator
    real(8),intent(in)          :: data(:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_4')
  end subroutine MPI_Bcast_Dble_4
  !
  subroutine MPI_Bcast_Dble_5(comm,data,root)
    integer,intent(in)          :: comm             !MPI communicator
    real(8),intent(in)          :: data(:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_5')
  end subroutine MPI_Bcast_Dble_5
  !
  subroutine MPI_Bcast_Dble_6(comm,data,root)
    integer,intent(in)          :: comm               !MPI communicator
    real(8),intent(in)          :: data(:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_6')
  end subroutine MPI_Bcast_Dble_6
  !
  subroutine MPI_Bcast_Dble_7(comm,data,root)
    integer,intent(in)          :: comm                 !MPI communicator
    real(8),intent(in)          :: data(:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_7')
  end subroutine MPI_Bcast_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Dble_8(comm,data,root)
    integer,intent(in)          :: comm                   !MPI communicator
    real(8),intent(in)          :: data(:,:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_8')
  end subroutine MPI_Bcast_Dble_8
#endif




  !!CMPLX8
  subroutine MPI_Bcast_Cmplx_0(comm,data,root)
    integer,intent(in)          :: comm  !MPI communicator
    complex(8),intent(in)       :: data  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root  !rank of the root process (default 0)
    complex(8),dimension(1)        :: data_
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_0')
  end subroutine MPI_Bcast_Cmplx_0
  !
  subroutine MPI_Bcast_Cmplx_1(comm,data,root)
    integer,intent(in)          :: comm     !MPI communicator
    complex(8),intent(in)       :: data(:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_1')
  end subroutine MPI_Bcast_Cmplx_1
  !
  subroutine MPI_Bcast_Cmplx_2(comm,data,root)
    integer,intent(in)          :: comm       !MPI communicator
    complex(8),intent(in)       :: data(:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_2')
  end subroutine MPI_Bcast_Cmplx_2
  !
  subroutine MPI_Bcast_Cmplx_3(comm,data,root)
    integer,intent(in)          :: comm         !MPI communicator
    complex(8),intent(in)       :: data(:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_3')
  end subroutine MPI_Bcast_Cmplx_3
  !
  subroutine MPI_Bcast_Cmplx_4(comm,data,root)
    integer,intent(in)          :: comm           !MPI communicator
    complex(8),intent(in)       :: data(:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_4')
  end subroutine MPI_Bcast_Cmplx_4
  !
  subroutine MPI_Bcast_Cmplx_5(comm,data,root)
    integer,intent(in)          :: comm             !MPI communicator
    complex(8),intent(in)       :: data(:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_5')
  end subroutine MPI_Bcast_Cmplx_5
  !
  subroutine MPI_Bcast_Cmplx_6(comm,data,root)
    integer,intent(in)          :: comm               !MPI communicator
    complex(8),intent(in)       :: data(:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_6')
  end subroutine MPI_Bcast_Cmplx_6
  !
  subroutine MPI_Bcast_Cmplx_7(comm,data,root)
    integer,intent(in)          :: comm                 !MPI communicator
    complex(8),intent(in)       :: data(:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_7')
  end subroutine MPI_Bcast_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Cmplx_8(comm,data,root)
    integer,intent(in)          :: comm                   !MPI communicator
    complex(8),intent(in)       :: data(:,:,:,:,:,:,:,:)  !data to broadcast (overwritten except on root)
    integer,intent(in),optional :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,root_,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_8')
  end subroutine MPI_Bcast_Cmplx_8
#endif















  !****************************************
  !              MPI ALLGATHER
  !****************************************
  !!BOOL
  subroutine MPI_Allgather_Bool_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    logical,intent(inout)       :: data  !gathered data of all processes, by rank
    logical,intent(in)          :: send  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,1,MPI_LOGICAL,data,1,MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_0')
  end subroutine MPI_Allgather_Bool_0
  !
  subroutine MPI_Allgather_Bool_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    logical,intent(inout)       :: data(:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_1')
  end subroutine MPI_Allgather_Bool_1
  !
  subroutine MPI_Allgather_Bool_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    logical,intent(inout)       :: data(:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_2')
  end subroutine MPI_Allgather_Bool_2
  !
  subroutine MPI_Allgather_Bool_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    logical,intent(inout)       :: data(:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_3')
  end subroutine MPI_Allgather_Bool_3
  !
  subroutine MPI_Allgather_Bool_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_4')
  end subroutine MPI_Allgather_Bool_4
  !
  subroutine MPI_Allgather_Bool_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_5')
  end subroutine MPI_Allgather_Bool_5
  !
  subroutine MPI_Allgather_Bool_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_6')
  end subroutine MPI_Allgather_Bool_6
  !
  subroutine MPI_Allgather_Bool_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_7')
  end subroutine MPI_Allgather_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Bool_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    logical,intent(in)          :: send(:,:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_8')
  end subroutine MPI_Allgather_Bool_8
#endif


  !!INTEGER
  subroutine MPI_Allgather_Int_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    integer,intent(inout)       :: data  !gathered data of all processes, by rank
    integer,intent(in)          :: send  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,1,MPI_INTEGER,data,1,MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_0')
  end subroutine MPI_Allgather_Int_0
  !
  subroutine MPI_Allgather_Int_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    integer,intent(inout)       :: data(:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_1')
  end subroutine MPI_Allgather_Int_1
  !
  subroutine MPI_Allgather_Int_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    integer,intent(inout)       :: data(:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_2')
  end subroutine MPI_Allgather_Int_2
  !
  subroutine MPI_Allgather_Int_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    integer,intent(inout)       :: data(:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_3')
  end subroutine MPI_Allgather_Int_3
  !
  subroutine MPI_Allgather_Int_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_4')
  end subroutine MPI_Allgather_Int_4
  !
  subroutine MPI_Allgather_Int_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_5')
  end subroutine MPI_Allgather_Int_5
  !
  subroutine MPI_Allgather_Int_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_6')
  end subroutine MPI_Allgather_Int_6
  !
  subroutine MPI_Allgather_Int_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_7')
  end subroutine MPI_Allgather_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Int_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    integer,intent(in)          :: send(:,:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_8')
  end subroutine MPI_Allgather_Int_8
#endif




  !!REAL8
  subroutine MPI_Allgather_Dble_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    real(8),intent(inout)       :: data  !gathered data of all processes, by rank
    real(8),intent(in)          :: send  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,1,MPI_DOUBLE_PRECISION,data,1,MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_0')
  end subroutine MPI_Allgather_Dble_0
  !
  subroutine MPI_Allgather_Dble_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    real(8),intent(inout)       :: data(:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_1')
  end subroutine MPI_Allgather_Dble_1
  !
  subroutine MPI_Allgather_Dble_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    real(8),intent(inout)       :: data(:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_2')
  end subroutine MPI_Allgather_Dble_2
  !
  subroutine MPI_Allgather_Dble_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    real(8),intent(inout)       :: data(:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_3')
  end subroutine MPI_Allgather_Dble_3
  !
  subroutine MPI_Allgather_Dble_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_4')
  end subroutine MPI_Allgather_Dble_4
  !
  subroutine MPI_Allgather_Dble_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_5')
  end subroutine MPI_Allgather_Dble_5
  !
  subroutine MPI_Allgather_Dble_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_6')
  end subroutine MPI_Allgather_Dble_6
  !
  subroutine MPI_Allgather_Dble_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_7')
  end subroutine MPI_Allgather_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Dble_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    real(8),intent(in)          :: send(:,:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_8')
  end subroutine MPI_Allgather_Dble_8
#endif


  !!CMPLX8
  subroutine MPI_Allgather_Cmplx_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    complex(8),intent(inout)    :: data  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,1,MPI_DOUBLE_COMPLEX,data,1,MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_0')
  end subroutine MPI_Allgather_Cmplx_0
  !
  subroutine MPI_Allgather_Cmplx_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    complex(8),intent(inout)    :: data(:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_1')
  end subroutine MPI_Allgather_Cmplx_1
  !
  subroutine MPI_Allgather_Cmplx_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    complex(8),intent(inout)    :: data(:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_2')
  end subroutine MPI_Allgather_Cmplx_2
  !
  subroutine MPI_Allgather_Cmplx_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_3')
  end subroutine MPI_Allgather_Cmplx_3
  !
  subroutine MPI_Allgather_Cmplx_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_4')
  end subroutine MPI_Allgather_Cmplx_4
  !
  subroutine MPI_Allgather_Cmplx_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_5')
  end subroutine MPI_Allgather_Cmplx_5
  !
  subroutine MPI_Allgather_Cmplx_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_6')
  end subroutine MPI_Allgather_Cmplx_6
  !
  subroutine MPI_Allgather_Cmplx_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_7')
  end subroutine MPI_Allgather_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Cmplx_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:,:)  !gathered data of all processes, by rank
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:,:)  !data sent by each process
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_8')
  end subroutine MPI_Allgather_Cmplx_8
#endif
































  !****************************************
  !              MPI ALLREDUCE
  !****************************************
  !!BOOL
  subroutine MPI_Allreduce_Bool_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    logical,intent(inout)       :: data  !sum of send over the processes
    logical,intent(in)          :: send  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,1,MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_0')
  end subroutine MPI_Allreduce_Bool_0
  !
  subroutine MPI_Allreduce_Bool_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    logical,intent(inout)       :: data(:)  !sum of send over the processes
    logical,intent(in)          :: send(:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_1')
  end subroutine MPI_Allreduce_Bool_1
  !
  subroutine MPI_Allreduce_Bool_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    logical,intent(inout)       :: data(:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_2')
  end subroutine MPI_Allreduce_Bool_2
  !
  subroutine MPI_Allreduce_Bool_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    logical,intent(inout)       :: data(:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_3')
  end subroutine MPI_Allreduce_Bool_3
  !
  subroutine MPI_Allreduce_Bool_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_4')
  end subroutine MPI_Allreduce_Bool_4
  !
  subroutine MPI_Allreduce_Bool_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_5')
  end subroutine MPI_Allreduce_Bool_5
  !
  subroutine MPI_Allreduce_Bool_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_6')
  end subroutine MPI_Allreduce_Bool_6
  !
  subroutine MPI_Allreduce_Bool_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_7')
  end subroutine MPI_Allreduce_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Bool_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    logical,intent(inout)       :: data(:,:,:,:,:,:,:,:)  !sum of send over the processes
    logical,intent(in)          :: send(:,:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_8')
  end subroutine MPI_Allreduce_Bool_8
#endif






  !!INTEGER
  subroutine MPI_Allreduce_Int_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    integer,intent(inout)       :: data  !sum of send over the processes
    integer,intent(in)          :: send  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_0')
  end subroutine MPI_Allreduce_Int_0
  !
  subroutine MPI_Allreduce_Int_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    integer,intent(inout)       :: data(:)  !sum of send over the processes
    integer,intent(in)          :: send(:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_1')
  end subroutine MPI_Allreduce_Int_1
  !
  subroutine MPI_Allreduce_Int_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    integer,intent(inout)       :: data(:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_2')
  end subroutine MPI_Allreduce_Int_2
  !
  subroutine MPI_Allreduce_Int_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    integer,intent(inout)       :: data(:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_3')
  end subroutine MPI_Allreduce_Int_3
  !
  subroutine MPI_Allreduce_Int_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_4')
  end subroutine MPI_Allreduce_Int_4
  !
  subroutine MPI_Allreduce_Int_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_5')
  end subroutine MPI_Allreduce_Int_5
  !
  subroutine MPI_Allreduce_Int_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_6')
  end subroutine MPI_Allreduce_Int_6
  !
  subroutine MPI_Allreduce_Int_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_7')
  end subroutine MPI_Allreduce_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Int_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    integer,intent(inout)       :: data(:,:,:,:,:,:,:,:)  !sum of send over the processes
    integer,intent(in)          :: send(:,:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_8')
  end subroutine MPI_Allreduce_Int_8
#endif




  !!REAL8
  subroutine MPI_Allreduce_Dble_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    real(8),intent(inout)       :: data  !sum of send over the processes
    real(8),intent(in)          :: send  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_0')
  end subroutine MPI_Allreduce_Dble_0
  !
  subroutine MPI_Allreduce_Dble_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    real(8),intent(inout)       :: data(:)  !sum of send over the processes
    real(8),intent(in)          :: send(:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_1')
  end subroutine MPI_Allreduce_Dble_1
  !
  subroutine MPI_Allreduce_Dble_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    real(8),intent(inout)       :: data(:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_2')
  end subroutine MPI_Allreduce_Dble_2
  !
  subroutine MPI_Allreduce_Dble_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    real(8),intent(inout)       :: data(:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_3')
  end subroutine MPI_Allreduce_Dble_3
  !
  subroutine MPI_Allreduce_Dble_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_4')
  end subroutine MPI_Allreduce_Dble_4
  !
  subroutine MPI_Allreduce_Dble_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_5')
  end subroutine MPI_Allreduce_Dble_5
  !
  subroutine MPI_Allreduce_Dble_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_6')
  end subroutine MPI_Allreduce_Dble_6
  !
  subroutine MPI_Allreduce_Dble_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_7')
  end subroutine MPI_Allreduce_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Dble_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:,:)  !sum of send over the processes
    real(8),intent(in)          :: send(:,:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_8')
  end subroutine MPI_Allreduce_Dble_8
#endif


  !!CMPLX8
  subroutine MPI_Allreduce_Cmplx_0(comm,send,data)
    integer,intent(in)          :: comm  !MPI communicator
    complex(8),intent(inout)    :: data  !sum of send over the processes
    complex(8),intent(in)       :: send  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_0')
  end subroutine MPI_Allreduce_Cmplx_0
  !
  subroutine MPI_Allreduce_Cmplx_1(comm,send,data)
    integer,intent(in)          :: comm     !MPI communicator
    complex(8),intent(inout)    :: data(:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_1')
  end subroutine MPI_Allreduce_Cmplx_1
  !
  subroutine MPI_Allreduce_Cmplx_2(comm,send,data)
    integer,intent(in)          :: comm       !MPI communicator
    complex(8),intent(inout)    :: data(:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_2')
  end subroutine MPI_Allreduce_Cmplx_2
  !
  subroutine MPI_Allreduce_Cmplx_3(comm,send,data)
    integer,intent(in)          :: comm         !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_3')
  end subroutine MPI_Allreduce_Cmplx_3
  !
  subroutine MPI_Allreduce_Cmplx_4(comm,send,data)
    integer,intent(in)          :: comm           !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_4')
  end subroutine MPI_Allreduce_Cmplx_4
  !
  subroutine MPI_Allreduce_Cmplx_5(comm,send,data)
    integer,intent(in)          :: comm             !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_5')
  end subroutine MPI_Allreduce_Cmplx_5
  !
  subroutine MPI_Allreduce_Cmplx_6(comm,send,data)
    integer,intent(in)          :: comm               !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_6')
  end subroutine MPI_Allreduce_Cmplx_6
  !
  subroutine MPI_Allreduce_Cmplx_7(comm,send,data)
    integer,intent(in)          :: comm                 !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_7')
  end subroutine MPI_Allreduce_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Cmplx_8(comm,send,data)
    integer,intent(in)          :: comm                   !MPI communicator
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:,:)  !sum of send over the processes
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:,:)  !data to sum over the processes
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_8')
  end subroutine MPI_Allreduce_Cmplx_8
#endif



























  !****************************************
  !    MPI REDUCE IN_PLACE: MPI_SUM
  !****************************************
  !!BOOL
  subroutine MPI_ReduceSum_Bool_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    logical,intent(inout)          :: data  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv  !receives the sum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm) 
    if(present(recv))then
       call MPI_Reduce(data,recv,1,MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_0')
  end subroutine MPI_ReduceSum_Bool_0
  !
  subroutine MPI_ReduceSum_Bool_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    logical,intent(inout)          :: data(:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:)  !receives the sum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_1')
  end subroutine MPI_ReduceSum_Bool_1
  !
  subroutine MPI_ReduceSum_Bool_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    logical,intent(inout)          :: data(:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_2')
  end subroutine MPI_ReduceSum_Bool_2
  !
  subroutine MPI_ReduceSum_Bool_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    logical,intent(inout)          :: data(:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_3')
  end subroutine MPI_ReduceSum_Bool_3
  !
  subroutine MPI_ReduceSum_Bool_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_4')
  end subroutine MPI_ReduceSum_Bool_4
  !
  subroutine MPI_ReduceSum_Bool_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_5')
  end subroutine MPI_ReduceSum_Bool_5
  !
  subroutine MPI_ReduceSum_Bool_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_6')
  end subroutine MPI_ReduceSum_Bool_6
  !
  subroutine MPI_ReduceSum_Bool_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_7')
  end subroutine MPI_ReduceSum_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceSum_Bool_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_Bool_8')
  end subroutine MPI_ReduceSum_Bool_8
#endif




  !INTEGER
  subroutine MPI_ReduceSum_int_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    integer,intent(inout)          :: data  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv  !receives the sum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,1,MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_0')
  end subroutine MPI_ReduceSum_int_0
  !
  subroutine MPI_ReduceSum_int_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    integer,intent(inout)          :: data(:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:)  !receives the sum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_1')
  end subroutine MPI_ReduceSum_int_1
  !
  subroutine MPI_ReduceSum_int_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    integer,intent(inout)          :: data(:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_2')
  end subroutine MPI_ReduceSum_int_2
  !
  subroutine MPI_ReduceSum_int_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    integer,intent(inout)          :: data(:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_3')
  end subroutine MPI_ReduceSum_int_3
  !
  subroutine MPI_ReduceSum_int_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_4')
  end subroutine MPI_ReduceSum_int_4
  !
  subroutine MPI_ReduceSum_int_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_5')
  end subroutine MPI_ReduceSum_int_5
  !
  subroutine MPI_ReduceSum_int_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_6')
  end subroutine MPI_ReduceSum_int_6
  !
  subroutine MPI_ReduceSum_int_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_7')
  end subroutine MPI_ReduceSum_int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceSum_int_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_int_8')
  end subroutine MPI_ReduceSum_int_8
#endif

  !REAL8
  subroutine MPI_ReduceSum_dble_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    real(8),intent(inout)          :: data  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv  !receives the sum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,1,MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_0')
  end subroutine MPI_ReduceSum_dble_0
  !
  subroutine MPI_ReduceSum_dble_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    real(8),intent(inout)          :: data(:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:)  !receives the sum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_1')
  end subroutine MPI_ReduceSum_dble_1
  !
  subroutine MPI_ReduceSum_dble_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    real(8),intent(inout)          :: data(:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_2')
  end subroutine MPI_ReduceSum_dble_2
  !
  subroutine MPI_ReduceSum_dble_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    real(8),intent(inout)          :: data(:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_3')
  end subroutine MPI_ReduceSum_dble_3
  !
  subroutine MPI_ReduceSum_dble_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_4')
  end subroutine MPI_ReduceSum_dble_4
  !
  subroutine MPI_ReduceSum_dble_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_5')
  end subroutine MPI_ReduceSum_dble_5
  !
  subroutine MPI_ReduceSum_dble_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_6')
  end subroutine MPI_ReduceSum_dble_6
  !
  subroutine MPI_ReduceSum_dble_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_7')
  end subroutine MPI_ReduceSum_dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceSum_dble_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_dble_8')
  end subroutine MPI_ReduceSum_dble_8
#endif



  !CMPLX8
  subroutine MPI_ReduceSum_cmplx_0(comm,data,recv,root)
    integer,intent(in)                :: comm  !MPI communicator
    complex(8),intent(inout)          :: data  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv  !receives the sum on root if present
    integer,intent(in),optional       :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,1,MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_0')
  end subroutine MPI_ReduceSum_cmplx_0
  !
  subroutine MPI_ReduceSum_cmplx_1(comm,data,recv,root)
    integer,intent(in)                :: comm     !MPI communicator
    complex(8),intent(inout)          :: data(:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:)  !receives the sum on root if present
    integer,intent(in),optional       :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_1')
  end subroutine MPI_ReduceSum_cmplx_1
  !
  subroutine MPI_ReduceSum_cmplx_2(comm,data,recv,root)
    integer,intent(in)                :: comm       !MPI communicator
    complex(8),intent(inout)          :: data(:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_2')
  end subroutine MPI_ReduceSum_cmplx_2
  !
  subroutine MPI_ReduceSum_cmplx_3(comm,data,recv,root)
    integer,intent(in)                :: comm         !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_3')
  end subroutine MPI_ReduceSum_cmplx_3
  !
  subroutine MPI_ReduceSum_cmplx_4(comm,data,recv,root)
    integer,intent(in)                :: comm           !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_4')
  end subroutine MPI_ReduceSum_cmplx_4
  !
  subroutine MPI_ReduceSum_cmplx_5(comm,data,recv,root)
    integer,intent(in)                :: comm             !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_5')
  end subroutine MPI_ReduceSum_cmplx_5
  !
  subroutine MPI_ReduceSum_cmplx_6(comm,data,recv,root)
    integer,intent(in)                :: comm               !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_6')
  end subroutine MPI_ReduceSum_cmplx_6
  !
  subroutine MPI_ReduceSum_cmplx_7(comm,data,recv,root)
    integer,intent(in)                :: comm                 !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_7')
  end subroutine MPI_ReduceSum_cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceSum_cmplx_8(comm,data,recv,root)
    integer,intent(in)                :: comm                   !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to sum (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the sum on root if present
    integer,intent(in),optional       :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceSum_cmplx_8')
  end subroutine MPI_ReduceSum_cmplx_8
#endif
















  !****************************************
  !    MPI REDUCE IN_PLACE: MPI_MAX
  !****************************************
  !!BOOL
  subroutine MPI_ReduceMax_Bool_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    logical,intent(inout)          :: data  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv  !receives the maximum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm) 
    if(present(recv))then
       call MPI_Reduce(data,recv,1,MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_0')
  end subroutine MPI_ReduceMax_Bool_0
  !
  subroutine MPI_ReduceMax_Bool_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    logical,intent(inout)          :: data(:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_1')
  end subroutine MPI_ReduceMax_Bool_1
  !
  subroutine MPI_ReduceMax_Bool_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    logical,intent(inout)          :: data(:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_2')
  end subroutine MPI_ReduceMax_Bool_2
  !
  subroutine MPI_ReduceMax_Bool_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    logical,intent(inout)          :: data(:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_3')
  end subroutine MPI_ReduceMax_Bool_3
  !
  subroutine MPI_ReduceMax_Bool_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_4')
  end subroutine MPI_ReduceMax_Bool_4
  !
  subroutine MPI_ReduceMax_Bool_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_5')
  end subroutine MPI_ReduceMax_Bool_5
  !
  subroutine MPI_ReduceMax_Bool_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_6')
  end subroutine MPI_ReduceMax_Bool_6
  !
  subroutine MPI_ReduceMax_Bool_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_7')
  end subroutine MPI_ReduceMax_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceMax_Bool_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    logical,intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    logical,intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       call MPI_Reduce(data,recv,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_LOGICAL,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_Bool_8')
  end subroutine MPI_ReduceMax_Bool_8
#endif




  !INTEGER
  subroutine MPI_ReduceMax_int_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    integer,intent(inout)          :: data  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv  !receives the maximum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,1,MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_0')
  end subroutine MPI_ReduceMax_int_0
  !
  subroutine MPI_ReduceMax_int_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    integer,intent(inout)          :: data(:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_1')
  end subroutine MPI_ReduceMax_int_1
  !
  subroutine MPI_ReduceMax_int_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    integer,intent(inout)          :: data(:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_2')
  end subroutine MPI_ReduceMax_int_2
  !
  subroutine MPI_ReduceMax_int_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    integer,intent(inout)          :: data(:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_3')
  end subroutine MPI_ReduceMax_int_3
  !
  subroutine MPI_ReduceMax_int_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_4')
  end subroutine MPI_ReduceMax_int_4
  !
  subroutine MPI_ReduceMax_int_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_5')
  end subroutine MPI_ReduceMax_int_5
  !
  subroutine MPI_ReduceMax_int_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_6')
  end subroutine MPI_ReduceMax_int_6
  !
  subroutine MPI_ReduceMax_int_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_7')
  end subroutine MPI_ReduceMax_int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceMax_int_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    integer,intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    integer,intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0
       call MPI_Reduce(data,recv,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_INTEGER,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_int_8')
  end subroutine MPI_ReduceMax_int_8
#endif

  !REAL8
  subroutine MPI_ReduceMax_dble_0(comm,data,recv,root)
    integer,intent(in)             :: comm  !MPI communicator
    real(8),intent(inout)          :: data  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv  !receives the maximum on root if present
    integer,intent(in),optional    :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,1,MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_0')
  end subroutine MPI_ReduceMax_dble_0
  !
  subroutine MPI_ReduceMax_dble_1(comm,data,recv,root)
    integer,intent(in)             :: comm     !MPI communicator
    real(8),intent(inout)          :: data(:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_1')
  end subroutine MPI_ReduceMax_dble_1
  !
  subroutine MPI_ReduceMax_dble_2(comm,data,recv,root)
    integer,intent(in)             :: comm       !MPI communicator
    real(8),intent(inout)          :: data(:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_2')
  end subroutine MPI_ReduceMax_dble_2
  !
  subroutine MPI_ReduceMax_dble_3(comm,data,recv,root)
    integer,intent(in)             :: comm         !MPI communicator
    real(8),intent(inout)          :: data(:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_3')
  end subroutine MPI_ReduceMax_dble_3
  !
  subroutine MPI_ReduceMax_dble_4(comm,data,recv,root)
    integer,intent(in)             :: comm           !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_4')
  end subroutine MPI_ReduceMax_dble_4
  !
  subroutine MPI_ReduceMax_dble_5(comm,data,recv,root)
    integer,intent(in)             :: comm             !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_5')
  end subroutine MPI_ReduceMax_dble_5
  !
  subroutine MPI_ReduceMax_dble_6(comm,data,recv,root)
    integer,intent(in)             :: comm               !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_6')
  end subroutine MPI_ReduceMax_dble_6
  !
  subroutine MPI_ReduceMax_dble_7(comm,data,recv,root)
    integer,intent(in)             :: comm                 !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_7')
  end subroutine MPI_ReduceMax_dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceMax_dble_8(comm,data,recv,root)
    integer,intent(in)             :: comm                   !MPI communicator
    real(8),intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    real(8),intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional    :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=0d0
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_PRECISION,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_dble_8')
  end subroutine MPI_ReduceMax_dble_8
#endif



  !CMPLX8
  subroutine MPI_ReduceMax_cmplx_0(comm,data,recv,root)
    integer,intent(in)                :: comm  !MPI communicator
    complex(8),intent(inout)          :: data  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv  !receives the maximum on root if present
    integer,intent(in),optional       :: root  !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,1,MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,1,MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,1,MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_0')
  end subroutine MPI_ReduceMax_cmplx_0
  !
  subroutine MPI_ReduceMax_cmplx_1(comm,data,recv,root)
    integer,intent(in)                :: comm     !MPI communicator
    complex(8),intent(inout)          :: data(:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root     !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_1')
  end subroutine MPI_ReduceMax_cmplx_1
  !
  subroutine MPI_ReduceMax_cmplx_2(comm,data,recv,root)
    integer,intent(in)                :: comm       !MPI communicator
    complex(8),intent(inout)          :: data(:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root       !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_2')
  end subroutine MPI_ReduceMax_cmplx_2
  !
  subroutine MPI_ReduceMax_cmplx_3(comm,data,recv,root)
    integer,intent(in)                :: comm         !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root         !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_3')
  end subroutine MPI_ReduceMax_cmplx_3
  !
  subroutine MPI_ReduceMax_cmplx_4(comm,data,recv,root)
    integer,intent(in)                :: comm           !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root           !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_4')
  end subroutine MPI_ReduceMax_cmplx_4
  !
  subroutine MPI_ReduceMax_cmplx_5(comm,data,recv,root)
    integer,intent(in)                :: comm             !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root             !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_5')
  end subroutine MPI_ReduceMax_cmplx_5
  !
  subroutine MPI_ReduceMax_cmplx_6(comm,data,recv,root)
    integer,intent(in)                :: comm               !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root               !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_6')
  end subroutine MPI_ReduceMax_cmplx_6
  !
  subroutine MPI_ReduceMax_cmplx_7(comm,data,recv,root)
    integer,intent(in)                :: comm                 !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root                 !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_7')
  end subroutine MPI_ReduceMax_cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_ReduceMax_cmplx_8(comm,data,recv,root)
    integer,intent(in)                :: comm                   !MPI communicator
    complex(8),intent(inout)          :: data(:,:,:,:,:,:,:,:)  !data to maximize (overwritten on rank 0 if no recv)
    complex(8),intent(inout),optional :: recv(:,:,:,:,:,:,:,:)  !receives the maximum on root if present
    integer,intent(in),optional       :: root                   !rank of the root process (default 0)
    root_=0;if(present(root))root_=root
    if(comm==MPI_COMM_NULL)return
    master=Get_master_MPI(comm)
    if(present(recv))then
       recv=dcmplx(0d0,0d0)
       call MPI_Reduce(data,recv,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
    else     
       if(master)then
          call MPI_Reduce(Mpi_In_Place,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       else
          call MPI_Reduce(data,data,size(data),MPI_DOUBLE_COMPLEX,MPI_MAX,root_,comm,ierr)
       endif
    endif
    call Error_MPI(sub='MPI_ReduceMax_cmplx_8')
  end subroutine MPI_ReduceMax_cmplx_8
#endif











# define STR_ERR_COMM      'invalid communicator in mpi call.'
# define STR_ERR_COUNT     'invalid count in mpi call.'
# define STR_ERR_TYPE      'invalid datatype in mpi call.'
# define STR_ERR_BUFFER    'invalid buffer in mpi call.'
# define STR_ERR_ROOT      'invalid root in mpi call.'
# define STR_ERR_ARG       'invalid argument in mpi call.'
# define STR_ERR_TAG       'invalid tag in mpi call.'
# define STR_ERR_RANK      'invalid rank in mpi call.'
# define STR_ERR_GROUP     'null group passed to mpi call.'
# define STR_ERR_OP        'invalid operation in mpi call.'
# define STR_ERR_TOPOLOGY  'invalid topology in mpi call.'
# define STR_ERR_DIMS      'illegal dimension argument in mpi call.'
# define STR_ERR_UNKNOWN   'unknown error in mpi call.'
# define STR_ERR_TRUNCATE  'message truncated on receive in mpi call.'
# define STR_ERR_OTHER     'other error in mpi call.'
# define STR_ERR_INTERN    'internal error code in mpi call.'
# define STR_ERR_IN_STATUS 'look in status for error value.'
# define STR_ERR_PENDING   'pending request in mpi call.'
# define STR_ERR_REQUEST   'illegal mpi_request handle in mpi call.'
# define STR_ERR_LASTCODE  'last error code in mpi call.'
  subroutine Error_MPI(err,sub)
    !This subroutine prints on the standard output the message associated to the MPI error code :f:var:`err`, preceded by the
    !label :f:var:`sub`. The default code is the one stored by the module after the last MPI call, and the default label is
    !:code:`MPI_Get_Error:`. Nothing is printed for :code:`MPI_SUCCESS` or for an unrecognized code, and the program is never
    !stopped.
    !
    integer,optional,intent(in)          :: err  !MPI error code (default: the last code stored by the module)
    character(len=*),optional,intent(in) :: sub  !label of the calling routine, printed before the message
    integer                              :: err_
    character(len=128)                   :: sub_
    err_=ierr            ; if(present(err))err_=err
    sub_="MPI_Get_Error:"; if(present(sub))sub_=sub
    select case (err_)
    case (MPI_SUCCESS)
       return
       !
    case (MPI_ERR_COMM)
       write(*,'(2A)')  trim(sub_),STR_ERR_COMM
       !
    case (MPI_ERR_COUNT)
       write(*,'(2A)')  trim(sub_),STR_ERR_COUNT
       !
    case (MPI_ERR_TYPE)
       write(*,'(2A)')  trim(sub_),STR_ERR_TYPE
       !
    case (MPI_ERR_BUFFER)
       write(*,'(2A)')  trim(sub_),STR_ERR_BUFFER
       !
    case (MPI_ERR_ROOT)
       write(*,'(2A)')  trim(sub_),STR_ERR_ROOT
       !
    case (MPI_ERR_ARG)
       write(*,'(2A)')  trim(sub_),STR_ERR_ARG
       !
    case (MPI_ERR_TAG)
       write(*,'(2A)')  trim(sub_),STR_ERR_TAG
       !
    case (MPI_ERR_RANK)
       write(*,'(2A)')  trim(sub_),STR_ERR_RANK
       !
    case (MPI_ERR_GROUP)
       write(*,'(2A)')  trim(sub_),STR_ERR_GROUP
       !
    case (MPI_ERR_OP)
       write(*,'(2A)')  trim(sub_),STR_ERR_OP
       !
    case (MPI_ERR_TOPOLOGY)
       write(*,'(2A)')  trim(sub_),STR_ERR_TOPOLOGY
       !
    case (MPI_ERR_DIMS)
       write(*,'(2A)')  trim(sub_),STR_ERR_DIMS
       !
    case (MPI_ERR_UNKNOWN)
       write(*,'(2A)')  trim(sub_),STR_ERR_UNKNOWN
       !
    case (MPI_ERR_TRUNCATE)
       write(*,'(2A)')  trim(sub_),STR_ERR_TRUNCATE
       !
    case (MPI_ERR_OTHER)
       write(*,'(2A)')  trim(sub_),STR_ERR_OTHER
       !
    case (MPI_ERR_INTERN)
       write(*,'(2A)')  trim(sub_),STR_ERR_INTERN
       !
    case (MPI_ERR_IN_STATUS)
       write(*,'(2A)')  trim(sub_),STR_ERR_IN_STATUS
       !
    case (MPI_ERR_PENDING)
       write(*,'(2A)')  trim(sub_),STR_ERR_PENDING
       !
    case (MPI_ERR_REQUEST)
       write(*,'(2A)')  trim(sub_),STR_ERR_REQUEST
       !
    case (MPI_ERR_LASTCODE)
       write(*,'(2A)')  trim(sub_),STR_ERR_LASTCODE
       !
    case default
       return
       !
    end select
  end subroutine Error_MPI



#else



  public :: Init_MPI
  public :: Finalize_MPI
  public :: StartMsg_MPI
  !
  public :: Check_MPI
  public :: Get_Size_MPI
  public :: Get_Rank_MPI
  public :: Get_Master_MPI
  public :: Get_Last_MPI
  !

  integer :: size
  integer :: rank
  integer :: ierr

contains


  !****************************************
  !              MPI START/STOP
  !****************************************
  subroutine Init_MPI()
    !This is the serial version, used when compiling without :code:`_MPI`. The subroutine does nothing.
    !
    return
  end subroutine Init_MPI

  subroutine Finalize_MPI()
    !This is the serial version, used when compiling without :code:`_MPI`. The subroutine does nothing.
    !
    return
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    !This is the serial version, used when compiling without :code:`_MPI`. The subroutine does nothing.
    !
    integer :: comm  !MPI communicator, not used in the serial version
    return
  end subroutine StartMsg_MPI



  !****************************************
  !              MPI TOOLS
  !****************************************
  function Check_MPI() result(bool)
    !This is the serial version, used when compiling without :code:`_MPI`. The function always returns :code:`.false.`, as MPI
    !is never initialized.
    !
    logical :: bool
    bool=.false.
  end function Check_MPI

  function Get_size_MPI(comm) result(size)
    !This is the serial version, used when compiling without :code:`_MPI`. The function always returns 1.
    !
    integer :: comm  !MPI communicator, not used in the serial version
    integer :: size
    size=1
  end function Get_size_MPI

  function Get_rank_MPI(comm) result(rank)
    !This is the serial version, used when compiling without :code:`_MPI`. The function always returns 0.
    !
    integer :: comm  !MPI communicator, not used in the serial version
    integer :: rank
    rank=0
  end function Get_rank_MPI

  function Get_master_MPI(comm) result(master)
    !This is the serial version, used when compiling without :code:`_MPI`. The function always returns :code:`.true.`.
    !
    integer :: comm  !MPI communicator, not used in the serial version
    logical :: master
    master=.true.
  end function Get_master_MPI

  function Get_last_MPI(comm) result(last)
    !This is the serial version, used when compiling without :code:`_MPI`. The function always returns :code:`.true.`.
    !
    integer :: comm  !MPI communicator, not used in the serial version
    logical :: last
    last=.true.
  end function Get_last_MPI


#endif

END MODULE SF_MPI







! function Get_Q_MPI(comm,N) result(mpiQ)
!   integer :: comm
!   integer :: N
!   integer :: size
!   integer :: rank
!   integer :: mpiQ
!   size = Get_size_MPI(comm)
!   mpiQ = N/size
! end function Get_Q_MPI

! function Get_R_MPI(comm,N) result(mpiR)
!   integer :: comm
!   integer :: N
!   integer :: size
!   integer :: rank
!   integer :: mpiR
!   logical :: last
!   size = Get_size_MPI(comm)
!   last = Get_last_MPI(comm)
!   mpiR=0
!   if(last)mpiR = mod(N,size)
! end function Get_R_MPI

! function Get_Chunk_MPI(comm,N) result(Nchunk)
!   integer :: comm
!   integer :: N
!   integer :: Nchunk
!   Nchunk = Get_Q_MPI(comm,N)+Get_R_MPI(comm,N)
! end function Get_Chunk_MPI



! function Get_Q_MPI(comm,N) result(mpiQ)
!   integer :: comm
!   integer :: N
!   integer :: mpiQ
!   mpiQ = N
! end function Get_Q_MPI

! function Get_R_MPI(comm,N) result(mpiR)
!   integer :: comm
!   integer :: N
!   integer :: mpiR
!   mpiR=0
! end function Get_R_MPI

! function Get_Chunk_MPI(comm,N) result(Nchunk)
!   integer :: comm
!   integer :: N
!   integer :: Nchunk
!   Nchunk = N
! end function Get_Chunk_MPI
