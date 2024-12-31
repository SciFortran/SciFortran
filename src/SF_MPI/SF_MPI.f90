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
  public :: AllGather_MPI
  public :: AllReduce_MPI


  integer :: ierr
  integer :: rank

contains


  !****************************************
  !              MPI START/STOP
  !****************************************
  subroutine Init_MPI(comm,msg)
    integer,optional :: comm
    logical,optional :: msg
    call MPI_Init(ierr)
    call Error_MPI(ierr,"MPI_Start")
    if(present(comm))comm=MPI_COMM_WORLD
    if(present(msg))then
       if(msg)call StartMsg_MPI(MPI_COMM_WORLD)
    endif
  end subroutine Init_MPI

  subroutine Finalize_MPI(comm)
    integer,optional :: comm
    call MPI_Finalize(ierr)
    call Error_MPI(ierr,"MPI_Stop")
    if(present(comm))comm=MPI_COMM_NULL
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    integer,optional :: comm
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
    integer,optional :: comm
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
    logical          :: bool    
    call MPI_Initialized(bool,ierr)
  end function check_MPI


  function get_size_MPI(comm) result(size)
    integer,optional :: comm
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
    integer,optional :: comm
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
    integer,optional :: comm
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
    integer,optional :: comm
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
    integer,intent(in)          :: comm
    logical,intent(in)          :: data
    integer,intent(in),optional :: root
    logical,dimension(1)        :: data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_0')
  end subroutine MPI_Bcast_Bool_0
  !
  subroutine MPI_Bcast_Bool_1(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_1')
  end subroutine MPI_Bcast_Bool_1
  !
  subroutine MPI_Bcast_Bool_2(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_2')
  end subroutine MPI_Bcast_Bool_2
  !
  subroutine MPI_Bcast_Bool_3(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_3')
  end subroutine MPI_Bcast_Bool_3
  !
  subroutine MPI_Bcast_Bool_4(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_4')
  end subroutine MPI_Bcast_Bool_4
  !
  subroutine MPI_Bcast_Bool_5(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_5')
  end subroutine MPI_Bcast_Bool_5
  !
  subroutine MPI_Bcast_Bool_6(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_6')
  end subroutine MPI_Bcast_Bool_6
  !
  subroutine MPI_Bcast_Bool_7(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_7')
  end subroutine MPI_Bcast_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Bool_8(comm,data,root)
    integer,intent(in)          :: comm
    logical,intent(in)          :: data(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Bool_8')
  end subroutine MPI_Bcast_Bool_8
#endif




  !! INTEGER
  subroutine MPI_Bcast_Int_0(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data
    integer,intent(in),optional :: root
    integer,dimension(1)        :: data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_0')
  end subroutine MPI_Bcast_Int_0
  !
  subroutine MPI_Bcast_Int_1(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_1')
  end subroutine MPI_Bcast_Int_1
  !
  subroutine MPI_Bcast_Int_2(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_2')
  end subroutine MPI_Bcast_Int_2
  !
  subroutine MPI_Bcast_Int_3(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_3')
  end subroutine MPI_Bcast_Int_3
  !
  subroutine MPI_Bcast_Int_4(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_4')
  end subroutine MPI_Bcast_Int_4
  !
  subroutine MPI_Bcast_Int_5(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_5')
  end subroutine MPI_Bcast_Int_5
  !
  subroutine MPI_Bcast_Int_6(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_6')
  end subroutine MPI_Bcast_Int_6
  !
  subroutine MPI_Bcast_Int_7(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_7')
  end subroutine MPI_Bcast_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Int_8(comm,data,root)
    integer,intent(in)          :: comm
    integer,intent(in)          :: data(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Int_7')
  end subroutine MPI_Bcast_Int_8
#endif



  !! REAL8
  subroutine MPI_Bcast_Dble_0(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data
    integer,intent(in),optional :: root
    real(8),dimension(1)        :: data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_0')
  end subroutine MPI_Bcast_Dble_0
  !
  subroutine MPI_Bcast_Dble_1(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_1')
  end subroutine MPI_Bcast_Dble_1
  !
  subroutine MPI_Bcast_Dble_2(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_2')
  end subroutine MPI_Bcast_Dble_2
  !
  subroutine MPI_Bcast_Dble_3(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_3')
  end subroutine MPI_Bcast_Dble_3
  !
  subroutine MPI_Bcast_Dble_4(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_4')
  end subroutine MPI_Bcast_Dble_4
  !
  subroutine MPI_Bcast_Dble_5(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_5')
  end subroutine MPI_Bcast_Dble_5
  !
  subroutine MPI_Bcast_Dble_6(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_6')
  end subroutine MPI_Bcast_Dble_6
  !
  subroutine MPI_Bcast_Dble_7(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_7')
  end subroutine MPI_Bcast_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Dble_8(comm,data,root)
    integer,intent(in)          :: comm
    real(8),intent(in)          :: data(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Dble_8')
  end subroutine MPI_Bcast_Dble_8
#endif




  !!CMPLX8
  subroutine MPI_Bcast_Cmplx_0(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data
    integer,intent(in),optional :: root
    complex(8),dimension(1)        :: data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !data_(1) = data
    call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_0')
  end subroutine MPI_Bcast_Cmplx_0
  !
  subroutine MPI_Bcast_Cmplx_1(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_1')
  end subroutine MPI_Bcast_Cmplx_1
  !
  subroutine MPI_Bcast_Cmplx_2(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_2')
  end subroutine MPI_Bcast_Cmplx_2
  !
  subroutine MPI_Bcast_Cmplx_3(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_3')
  end subroutine MPI_Bcast_Cmplx_3
  !
  subroutine MPI_Bcast_Cmplx_4(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_4')
  end subroutine MPI_Bcast_Cmplx_4
  !
  subroutine MPI_Bcast_Cmplx_5(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_5')
  end subroutine MPI_Bcast_Cmplx_5
  !
  subroutine MPI_Bcast_Cmplx_6(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_6')
  end subroutine MPI_Bcast_Cmplx_6
  !
  subroutine MPI_Bcast_Cmplx_7(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_7')
  end subroutine MPI_Bcast_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Bcast_Cmplx_8(comm,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(in)       :: data(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
    call Error_MPI(sub='MPI_Bcast_Cmplx_8')
  end subroutine MPI_Bcast_Cmplx_8
#endif















  !****************************************
  !              MPI ALLGATHER
  !****************************************
  !!BOOL
  subroutine MPI_Allgather_Bool_0(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data
    logical,intent(in)          :: send
    integer,intent(in),optional :: root
    logical,dimension(1)        :: send_,data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLGATHER(send,1,MPI_LOGICAL,data,1,MPI_LOGICAL,comm,ierr)
    !data = data_(1)
    call Error_MPI(sub='MPI_Allgather_Bool_0')
  end subroutine MPI_Allgather_Bool_0
  !
  subroutine MPI_Allgather_Bool_1(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:)
    logical,intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_1')
  end subroutine MPI_Allgather_Bool_1
  !
  subroutine MPI_Allgather_Bool_2(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:)
    logical,intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_2')
  end subroutine MPI_Allgather_Bool_2
  !
  subroutine MPI_Allgather_Bool_3(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:)
    logical,intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_3')
  end subroutine MPI_Allgather_Bool_3
  !
  subroutine MPI_Allgather_Bool_4(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_4')
  end subroutine MPI_Allgather_Bool_4
  !
  subroutine MPI_Allgather_Bool_5(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_5')
  end subroutine MPI_Allgather_Bool_5
  !
  subroutine MPI_Allgather_Bool_6(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_6')
  end subroutine MPI_Allgather_Bool_6
  !
  subroutine MPI_Allgather_Bool_7(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_7')
  end subroutine MPI_Allgather_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Bool_8(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Bool_8')
  end subroutine MPI_Allgather_Bool_8
#endif



  !!INTEGER
  subroutine MPI_Allgather_Int_0(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data
    integer,intent(in)          :: send
    integer,intent(in),optional :: root
    integer,dimension(1)        :: send_,data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLGATHER(send,1,MPI_INTEGER,data,1,MPI_INTEGER,comm,ierr)
    !data  = data_(1)
    call Error_MPI(sub='MPI_Allgather_Int_0')
  end subroutine MPI_Allgather_Int_0
  !
  subroutine MPI_Allgather_Int_1(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:)
    integer,intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_1')
  end subroutine MPI_Allgather_Int_1
  !
  subroutine MPI_Allgather_Int_2(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:)
    integer,intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_2')
  end subroutine MPI_Allgather_Int_2
  !
  subroutine MPI_Allgather_Int_3(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:)
    integer,intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_3')
  end subroutine MPI_Allgather_Int_3
  !
  subroutine MPI_Allgather_Int_4(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_4')
  end subroutine MPI_Allgather_Int_4
  !
  subroutine MPI_Allgather_Int_5(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_5')
  end subroutine MPI_Allgather_Int_5
  !
  subroutine MPI_Allgather_Int_6(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_6')
  end subroutine MPI_Allgather_Int_6
  !
  subroutine MPI_Allgather_Int_7(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_7')
  end subroutine MPI_Allgather_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Int_8(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Int_8')
  end subroutine MPI_Allgather_Int_8
#endif





  !!REAL8
  subroutine MPI_Allgather_Dble_0(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data
    real(8),intent(in)          :: send
    integer,intent(in),optional :: root
    real(8),dimension(1)        :: send_,data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLGATHER(send,1,MPI_DOUBLE_PRECISION,data,1,MPI_DOUBLE_PRECISION,comm,ierr)
    !data = data_(1)
    call Error_MPI(sub='MPI_Allgather_Dble_0')
  end subroutine MPI_Allgather_Dble_0
  !
  subroutine MPI_Allgather_Dble_1(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:)
    real(8),intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_1')
  end subroutine MPI_Allgather_Dble_1
  !
  subroutine MPI_Allgather_Dble_2(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:)
    real(8),intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_2')
  end subroutine MPI_Allgather_Dble_2
  !
  subroutine MPI_Allgather_Dble_3(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:)
    real(8),intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_3')
  end subroutine MPI_Allgather_Dble_3
  !
  subroutine MPI_Allgather_Dble_4(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_4')
  end subroutine MPI_Allgather_Dble_4
  !
  subroutine MPI_Allgather_Dble_5(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_5')
  end subroutine MPI_Allgather_Dble_5
  !
  subroutine MPI_Allgather_Dble_6(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_6')
  end subroutine MPI_Allgather_Dble_6
  !
  subroutine MPI_Allgather_Dble_7(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_7')
  end subroutine MPI_Allgather_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Dble_8(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Dble_8')
  end subroutine MPI_Allgather_Dble_8
#endif



  !!CMPLX8
  subroutine MPI_Allgather_Cmplx_0(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data
    complex(8),intent(in)       :: send
    integer,intent(in),optional :: root
    complex,dimension(1)        :: send_,data_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLGATHER(send,1,MPI_DOUBLE_COMPLEX,data,1,MPI_DOUBLE_COMPLEX,comm,ierr)
    !data = data_(1)
    call Error_MPI(sub='MPI_Allgather_Cmplx_0')
  end subroutine MPI_Allgather_Cmplx_0
  !
  subroutine MPI_Allgather_Cmplx_1(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:)
    complex(8),intent(in)       :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_1')
  end subroutine MPI_Allgather_Cmplx_1
  !
  subroutine MPI_Allgather_Cmplx_2(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:)
    complex(8),intent(in)       :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_2')
  end subroutine MPI_Allgather_Cmplx_2
  !
  subroutine MPI_Allgather_Cmplx_3(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:)
    complex(8),intent(in)       :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_3')
  end subroutine MPI_Allgather_Cmplx_3
  !
  subroutine MPI_Allgather_Cmplx_4(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_4')
  end subroutine MPI_Allgather_Cmplx_4
  !
  subroutine MPI_Allgather_Cmplx_5(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_5')
  end subroutine MPI_Allgather_Cmplx_5
  !
  subroutine MPI_Allgather_Cmplx_6(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_6')
  end subroutine MPI_Allgather_Cmplx_6
  !
  subroutine MPI_Allgather_Cmplx_7(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_7')
  end subroutine MPI_Allgather_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allgather_Cmplx_8(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
    call Error_MPI(sub='MPI_Allgather_Cmplx_8')
  end subroutine MPI_Allgather_Cmplx_8
#endif




























  !****************************************
  !              MPI ALLREDUCE
  !****************************************
  !!BOOL
  subroutine MPI_Allreduce_Bool_0(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data
    logical,intent(in)          :: send
    integer,intent(in),optional :: root
    logical,dimension(1)       :: data_,send_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLREDUCE(send,data,1,MPI_LOGICAL,MPI_SUM,comm,ierr)
    !data     = data_(1)
    call Error_MPI(sub='MPI_Allreduce_Bool_0')
  end subroutine MPI_Allreduce_Bool_0
  !
  subroutine MPI_Allreduce_Bool_1(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:)
    logical,intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_1')
  end subroutine MPI_Allreduce_Bool_1
  !
  subroutine MPI_Allreduce_Bool_2(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:)
    logical,intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_2')
  end subroutine MPI_Allreduce_Bool_2
  !
  subroutine MPI_Allreduce_Bool_3(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:)
    logical,intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_3')
  end subroutine MPI_Allreduce_Bool_3
  !
  subroutine MPI_Allreduce_Bool_4(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_4')
  end subroutine MPI_Allreduce_Bool_4
  !
  subroutine MPI_Allreduce_Bool_5(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_5')
  end subroutine MPI_Allreduce_Bool_5
  !
  subroutine MPI_Allreduce_Bool_6(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_6')
  end subroutine MPI_Allreduce_Bool_6
  !
  subroutine MPI_Allreduce_Bool_7(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_7')
  end subroutine MPI_Allreduce_Bool_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Bool_8(comm,send,data,root)
    integer,intent(in)          :: comm
    logical,intent(inout)       :: data(:,:,:,:,:,:,:,:)
    logical,intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Bool_8')
  end subroutine MPI_Allreduce_Bool_8
#endif






  !!INTEGER
  subroutine MPI_Allreduce_Int_0(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data
    integer,intent(in)          :: send
    integer,intent(in),optional :: root
    integer,dimension(1)       :: data_,send_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLREDUCE(send,data,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    !data     = data_(1)
    call Error_MPI(sub='MPI_Allreduce_Int_0')
  end subroutine MPI_Allreduce_Int_0
  !
  subroutine MPI_Allreduce_Int_1(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:)
    integer,intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_1')
  end subroutine MPI_Allreduce_Int_1
  !
  subroutine MPI_Allreduce_Int_2(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:)
    integer,intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_2')
  end subroutine MPI_Allreduce_Int_2
  !
  subroutine MPI_Allreduce_Int_3(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:)
    integer,intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_3')
  end subroutine MPI_Allreduce_Int_3
  !
  subroutine MPI_Allreduce_Int_4(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_4')
  end subroutine MPI_Allreduce_Int_4
  !
  subroutine MPI_Allreduce_Int_5(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_5')
  end subroutine MPI_Allreduce_Int_5
  !
  subroutine MPI_Allreduce_Int_6(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_6')
  end subroutine MPI_Allreduce_Int_6
  !
  subroutine MPI_Allreduce_Int_7(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_7')
  end subroutine MPI_Allreduce_Int_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Int_8(comm,send,data,root)
    integer,intent(in)          :: comm
    integer,intent(inout)       :: data(:,:,:,:,:,:,:,:)
    integer,intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Int_8')
  end subroutine MPI_Allreduce_Int_8
#endif




  !!REAL8
  subroutine MPI_Allreduce_Dble_0(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data
    real(8),intent(in)          :: send
    integer,intent(in),optional :: root
    real(8),dimension(1)        :: data_,send_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    !data     = data_(1)
    call Error_MPI(sub='MPI_Allreduce_Dble_0')
  end subroutine MPI_Allreduce_Dble_0
  !
  subroutine MPI_Allreduce_Dble_1(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:)
    real(8),intent(in)          :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_1')
  end subroutine MPI_Allreduce_Dble_1
  !
  subroutine MPI_Allreduce_Dble_2(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:)
    real(8),intent(in)          :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_2')
  end subroutine MPI_Allreduce_Dble_2
  !
  subroutine MPI_Allreduce_Dble_3(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:)
    real(8),intent(in)          :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_3')
  end subroutine MPI_Allreduce_Dble_3
  !
  subroutine MPI_Allreduce_Dble_4(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_4')
  end subroutine MPI_Allreduce_Dble_4
  !
  subroutine MPI_Allreduce_Dble_5(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_5')
  end subroutine MPI_Allreduce_Dble_5
  !
  subroutine MPI_Allreduce_Dble_6(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_6')
  end subroutine MPI_Allreduce_Dble_6
  !
  subroutine MPI_Allreduce_Dble_7(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_7')
  end subroutine MPI_Allreduce_Dble_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Dble_8(comm,send,data,root)
    integer,intent(in)          :: comm
    real(8),intent(inout)       :: data(:,:,:,:,:,:,:,:)
    real(8),intent(in)          :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Dble_8')
  end subroutine MPI_Allreduce_Dble_8
#endif



  !!CMPLX8
  subroutine MPI_Allreduce_Cmplx_0(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data
    complex(8),intent(in)       :: send
    integer,intent(in),optional :: root
    complex(8),dimension(1)     :: data_,send_
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    !send_(1) = send
    call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    !data     = data_(1)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_0')
  end subroutine MPI_Allreduce_Cmplx_0
  !
  subroutine MPI_Allreduce_Cmplx_1(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:)
    complex(8),intent(in)       :: send(:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_1')
  end subroutine MPI_Allreduce_Cmplx_1
  !
  subroutine MPI_Allreduce_Cmplx_2(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:)
    complex(8),intent(in)       :: send(:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_2')
  end subroutine MPI_Allreduce_Cmplx_2
  !
  subroutine MPI_Allreduce_Cmplx_3(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:)
    complex(8),intent(in)       :: send(:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_3')
  end subroutine MPI_Allreduce_Cmplx_3
  !
  subroutine MPI_Allreduce_Cmplx_4(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_4')
  end subroutine MPI_Allreduce_Cmplx_4
  !
  subroutine MPI_Allreduce_Cmplx_5(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_5')
  end subroutine MPI_Allreduce_Cmplx_5
  !
  subroutine MPI_Allreduce_Cmplx_6(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_6')
  end subroutine MPI_Allreduce_Cmplx_6
  !
  subroutine MPI_Allreduce_Cmplx_7(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_7')
  end subroutine MPI_Allreduce_Cmplx_7
  !
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine MPI_Allreduce_Cmplx_8(comm,send,data,root)
    integer,intent(in)          :: comm
    complex(8),intent(inout)    :: data(:,:,:,:,:,:,:,:)
    complex(8),intent(in)       :: send(:,:,:,:,:,:,:,:)
    integer,intent(in),optional :: root
    rank=0;if(present(root))rank=root
    if(comm==MPI_COMM_NULL)return
    call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call Error_MPI(sub='MPI_Allreduce_Cmplx_8')
  end subroutine MPI_Allreduce_Cmplx_8
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
    integer,optional,intent(in)          :: err
    character(len=*),optional,intent(in) :: sub
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
    return
  end subroutine Init_MPI

  subroutine Finalize_MPI()
    return
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    integer :: comm
    return
  end subroutine StartMsg_MPI



  !****************************************
  !              MPI TOOLS
  !****************************************
  function Check_MPI() result(bool)
    logical :: bool
    bool=.false.
  end function Check_MPI

  function Get_size_MPI(comm) result(size)
    integer :: comm
    integer :: size
    size=1
  end function Get_size_MPI

  function Get_rank_MPI(comm) result(rank)
    integer :: comm
    integer :: rank
    rank=0
  end function Get_rank_MPI

  function Get_master_MPI(comm) result(master)
    integer :: comm
    logical :: master
    master=.true.
  end function Get_master_MPI

  function Get_last_MPI(comm) result(last)
    integer :: comm
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
