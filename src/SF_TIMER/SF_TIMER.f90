module SF_TIMER
  !SciFortran module for timers
#ifdef _MPI
  USE SF_MPI
#endif
  implicit none
#ifdef _MPI
  include "mpif.h"
#endif
  private
  !
  !The months of the year:
  character(len=9),parameter,dimension(12) :: month = (/ &
       'January  ', &
       'February ', &
       'March    ', &
       'April    ', &
       'May      ', &
       'June     ', &
       'July     ', &
       'August   ', &
       'September', &
       'October  ', &
       'November ', &
       'December ' /)
  !from secs to day, hour, min
  integer,parameter                        :: secs_in_one_day=86400
  integer,parameter                        :: secs_in_one_hour=3600
  integer,parameter                        :: secs_in_one_min=60
  !
  integer(4),dimension(8),save             :: data
  integer(4),save                          :: year
  integer(4),save                          :: mese
  integer(4),save                          :: day
  integer(4),save                          :: h
  integer(4),save                          :: m
  integer(4),save                          :: s
  integer(4),save                          :: ms
  !Counter of the started timers:
  integer,save                             :: Tindex=0,Tprof=0
  integer,parameter                        :: Tmax=1000!
  !
  real,save                                :: time
  real,save                                :: old_time
  real,save                                :: dtime
  real,save                                :: elapsed_time
  real,save                                :: eta_time
  integer(8),dimension(Tmax)               :: Funit
  integer(8)                               :: st_rate=1,pr_rate=1
  !
  !Using Date and Time:
  integer,dimension(Tmax,8),save           :: dt_T_start=0
  integer,dimension(Tmax,8),save           :: dt_T_stop=0
  integer,dimension(Tmax,8),save           :: dt_T0=00
  integer,dimension(Tmax,8),save           :: dt_T1=0
  !Using Cpu_time
  real(8),dimension(Tmax),save                :: ct_T_start=0d0
  real(8),dimension(Tmax),save                :: ct_T_stop=0d0
  real(8),dimension(Tmax),save                :: ct_T0=0d0
  real(8),dimension(Tmax),save                :: ct_T1=0d0
  !Using System_clock
  real(8),dimension(Tmax),save          :: st_T_start=0
  real(8),dimension(Tmax),save          :: st_T_stop=0
  real(8),dimension(Tmax),save          :: st_T0=0
  real(8),dimension(Tmax),save          :: st_T1=0
  !
  !Using System_clock
  real(8),dimension(Tmax),save          :: pr_T_start=0
  real(8),dimension(Tmax),save          :: pr_T_stop=0
  real(8),dimension(Tmax),save          :: pr_T0=0
  real(8),dimension(Tmax),save          :: pr_T1=0
  !MPI
  integer                                  :: mpi_id
  integer                                  :: mpi_size
  logical                                  :: mpi_master

  interface start_progress
     !Generic interface to start a timer, to be stopped with :f:func_inline:`stop_progress`. It is an alias of
     !:f:func_inline:`start_timer`, whose arguments are documented there.
     !
     !At the moment of the call the timer records three clocks: the wall clock (:code:`system_clock`, or :code:`MPI_Wtime`
     !when compiled with MPI and MPI is initialized), the CPU time (:code:`cpu_time`) and the date and time
     !(:code:`date_and_time`). It also resets the state used by :f:func_inline:`eta` and :f:func_inline:`progress`, which need
     !therefore a started timer. Timers can be nested, the last started being the first to be stopped, up to 1000 at the same
     !time: the program stops otherwise.
     !
     !All the output of the timer, i.e. the :f:var:`title`, the messages of :f:func_inline:`eta` and
     !:f:func_inline:`progress`, and the elapsed time printed by :f:func_inline:`stop_progress`, is written to the unit
     !:f:var:`unit`. The default is the standard output (unit 6), plus the rank of the process when compiled with MPI, so that
     !each process writes to its own unit; unit 5 is replaced by 6.
     !
     module procedure :: start_timer
  end interface start_progress

  interface stop_progress
     !Generic interface to stop the most recently started timer, see :f:func_inline:`start_progress`. It is an alias of
     !:f:func_inline:`stop_timer`, whose argument is documented there.
     !
     !The elapsed CPU time and, in parentheses, the elapsed wall-clock time are printed on the unit of the timer, both as
     !:code:`hhh:mm:ss.mmm`, followed by :f:var:`msg` if present:
     !:code:`CPU time (wall time): hhh:mm:ss.mmm (hhh:mm:ss.mmm) : msg`. The timer is then released.
     !
     module procedure :: stop_timer
  end interface stop_progress

  public :: t_start,t_stop
  public :: start_timer,start_progress
  public :: stop_timer ,stop_progress
  public :: eta
  public :: progress



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : start a timer to measure elapsed time between two call
  !+-------------------------------------------------------------------+
  function t_start() result(pr_T)
  !This function starts a timer, to be read with :f:func_inline:`t_stop`, to measure the time elapsed between two calls, e.g.
  !for profiling. It is independent of :f:func_inline:`start_timer`. Timers can be nested, the last started being the first to
  !be stopped, up to 1000 at the same time: the program stops otherwise. The wall clock is :code:`MPI_Wtime` when compiled
  !with MPI and MPI is initialized, :code:`system_clock` otherwise. The returned value is always zero.
  !
    real(8)    :: pr_T  ! always 0
    integer(8) :: T_
    !
    Tprof=Tprof+1
    if(Tprof>Tmax)stop "start_timer error: too many timers started"
    !
#ifdef _MPI
    if(check_MPI())then
       pr_T_start(Tprof) = MPI_Wtime()
    else
       call system_clock(count_rate=pr_rate)
       call system_clock(count=T_)
       pr_T_start(Tprof) = dble(T_)
    endif
#else
    call system_clock(count_rate=pr_rate)
    call system_clock(count=T_)
    pr_T_start(Tprof) = dble(T_)
#endif
    pr_T=0d0
  end function t_start


  !+-------------------------------------------------------------------+
  !PURPOSE  : stop the timer and get the partial time
  !+-------------------------------------------------------------------+
  function t_stop() result(pr_T)
  !This function stops the timer most recently started with :f:func_inline:`t_start`, releases it, and returns the time
  !elapsed since its start, in seconds.
  !
    real(8)    :: pr_T  ! time elapsed since the matching t_start, in seconds
    integer(8) :: T_
    !
#ifdef _MPI    
    if(check_MPI())then
       pr_T_stop(Tprof) = MPI_Wtime()
       pr_T = pr_T_stop(Tprof)-pr_T_start(Tprof)
    else
       call system_clock(count=T_)
       pr_T_stop(Tprof)=dble(T_)
       pr_T = dble(pr_T_stop(Tprof)-pr_T_start(Tprof))/st_rate
    endif
#else
    call system_clock(count=T_)
    pr_T_stop(Tprof)=dble(T_)
    pr_T = dble(pr_T_stop(Tprof)-pr_T_start(Tprof))/st_rate
#endif
    !
    pr_T_start(Tprof)   = 0d0
    pr_T_stop(Tprof)    = 0d0
    !
    if(Tprof>1)then
       Tprof=Tprof-1
    else
       Tprof=0
    endif
  end function t_stop







  !+-------------------------------------------------------------------+
  !PURPOSE  : start a timer to measure elapsed time between two call
  !+-------------------------------------------------------------------+
  subroutine start_timer(title,unit)
    character(len=*),optional :: title  ! optional: title, printed on the unit of the timer
    integer,optional          :: unit   ! optional: output unit of the timer, default 6 (+ MPI rank)
    integer(8)                :: T_
    real :: cT_
    !
    Tindex=Tindex+1
    if(Tindex>Tmax)stop "start_timer error: too many timers started"
    !
    mpi_id = 0
#ifdef _MPI    
    if(check_MPI())mpi_id = get_Rank_MPI()
#endif
    funit(Tindex)=6+mpi_id;if(present(unit))funit(Tindex)=unit
    if(funit(Tindex)==5)funit(Tindex)=6
    !
    if(present(title))write(funit(Tindex),"(1x,A)")trim(title)
    !
#ifdef _MPI    
    if(check_MPI())then
       st_T_start(Tindex) = MPI_Wtime()
    else
       call system_clock(count_rate=st_rate)
       call system_clock(count=T_)
       st_T_start(Tindex)=dble(T_)
    endif
#else
    call system_clock(count_rate=st_rate)
    call system_clock(count=T_)
    st_T_start(Tindex)=dble(T_)
#endif
    call cpu_time(cT_)
    ct_T_start(Tindex)=dble(cT_)
    call date_and_time(values=dt_T_start(Tindex,:))
    st_T0 = st_T_start
    ct_T0 = ct_T_start
    dt_T0 = dt_T_start
    !
    !init variables for ETA:
    elapsed_time = 0.0
    old_time     = 0.0
    time         = 0.0
  end subroutine start_timer



  !+-------------------------------------------------------------------+
  !PURPOSE  : stop the timer and get the partial time
  !+-------------------------------------------------------------------+
  subroutine stop_timer(msg)
    character(len=*),optional :: msg  ! optional: message, printed after the elapsed time
    real(8)              :: ct_T
    real(8)              :: st_T
    integer,dimension(8) :: dt_T
    integer(8)           :: T_
    real :: cT_
    !
    !
#ifdef _MPI    
    if(check_MPI())then
       st_T_stop(Tindex) = MPI_Wtime()
    else
       call system_clock(count=T_)
       st_T_stop(Tindex)=dble(T_)
    endif
#else
    call system_clock(count=T_)
    st_T_stop(Tindex)=dble(T_)
#endif
    call cpu_time(cT_)
    ct_T_stop(Tindex)=dble(cT_)
    call date_and_time(values=dt_T_stop(Tindex,:))
    ct_T = ct_time_difference(ct_T_stop(Tindex),ct_T_start(Tindex))
    st_T = st_time_difference(st_T_stop(Tindex),st_T_start(Tindex))
    dt_T = dt_time_difference(dt_T_stop(Tindex,:),dt_T_start(Tindex,:))
    !
    if(present(msg))then
       call print_total_time(ct_T,st_T,dt_T,msg)
    else
       call print_total_time(ct_T,st_T,dt_T)
    endif
    ct_T_start(Tindex)   = 0d0
    st_T_start(Tindex)   = 0d0
    dt_T_start(Tindex,:) = 0
    ct_T_stop(Tindex)    = 0d0
    st_T_stop(Tindex)    = 0d0
    dt_T_stop(Tindex,:)  = 0
    !
    funit(Tindex)        = 6
    !
    if(Tindex>1)then
       Tindex=Tindex-1
    else
       Tindex=0
    endif
  end subroutine stop_timer







  !+-------------------------------------------------------------------+
  !PURPOSE  : get Expected Time of Arrival
  !+-------------------------------------------------------------------+
  subroutine eta(i,L,step,method)
  !This subroutine prints the Expected Time of Arrival (ETA) of a loop of :f:var:`L` iterations, to be called at every
  !iteration :f:var:`i`. It requires a timer started with :f:func_inline:`start_timer`, and it writes on the unit of that
  !timer.
  !
  !The ETA is the average time per iteration spent so far, times the number of iterations, minus the elapsed time. It is
  !printed only when the completion percentage reaches a multiple of :f:var:`step`, as :code:`NNN% |ETA: hhh:mm:ss.mmm`; at
  !10%, 50% and 100% the date and time are appended. The internal state is reset when :code:`i==L`, so that the routine can be
  !used for the next loop. The clock used is selected by the first character of :f:var:`method`: :code:`c` (default) for
  !:code:`cpu_time`, :code:`d` for :code:`date_and_time`, anything else, e.g. :code:`s`, for :code:`system_clock` (wall time).
  !
    integer                   :: i       ! current iteration
    integer                   :: L       ! total number of iterations
    integer,optional          :: step    ! optional: print every step percent, default 10
    character(len=*),optional :: method  ! optional: clock, 'c' cpu_time (default), 'd' date_and_time, 's' system_clock
    !
    integer,save              :: mod_print
    integer                   :: percent,iprint
    integer,save              :: older=0,oldiprint=0
    logical                   :: esc,fullprint
    integer(4),dimension(8)   :: dummy
    character(len=1)          :: method_
    character(len=80)         :: message
    logical,save              :: lentry=.true.
    !
    method_='c';if(present(method))method_=trim(method) !s=system_clock(wall_time),c=cpu_time,d=date_and_time
    !
    !Preambolo:
    if(lentry)then
       mod_print=10;if(present(step))mod_print=step
       if(funit(Tindex)/=6)then
          write(*,"(2x,A,I4)")"+ETA --> unit: ",funit(Tindex)
       endif
       lentry=.false.
    endif
    !
    if(i==L)lentry=.true.
    !
    !avoid repetition of percentage (within the error)
    percent=int(100.0*i/L)
    if(percent==0)return
    if(percent==older)return
    if(percent<mod_print)return
    older=percent
    !
    !set step for printing:
    esc=.true.
    iprint=percent/mod_print
    if(iprint/=oldiprint)esc=.false.
    if(esc)return
    oldiprint=iprint
    !
    !check if fullprint (date) is needed
    fullprint=.false.;if(any(percent==[10,50,100]))fullprint=.true.
    !
    old_time     = time
    time         = total_time(method_) !time since the start of the clock
    !Get ETA:
    dtime        = time-old_time     
    elapsed_time = elapsed_time + dtime
    dtime        = elapsed_time/real(i,4)
    eta_time     = dtime*L - elapsed_time
    !
    ms           = int(fraction(eta_time)*1000.0)
    h            = int(eta_time/secs_in_one_hour)
    m            = int((eta_time - h*secs_in_one_hour)/secs_in_one_min)
    s            = int(eta_time - h*secs_in_one_hour - m*secs_in_one_min)
    !

    write(funit(Tindex),"(1x,1i3,1a7,i3.3,a1,i2.2,a1,i2.2,a1,i3.3)",advance='no')percent,"% |ETA: ",h,":",m,":",s,".",ms
    if(fullprint)then
       call date_and_time(values=dummy)
       write(funit(Tindex),"(a2,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2)")&
            " @",dummy(3),trim(month(dummy(2))),dummy(1),dummy(5),':',dummy(6),':',dummy(7)
    else
       write(funit(Tindex),"(1X)")
    endif
  end subroutine eta






  !+-------------------------------------------------------------------+
  !TYPE     : Subroutine
  !+-------------------------------------------------------------------+
  subroutine progress(i,imax,method)
  !This subroutine prints a progress bar with the Expected Time of Arrival (ETA), to be called at every iteration :f:var:`i`
  !of a loop of :f:var:`imax` iterations. It requires a timer started with :f:func_inline:`start_timer`, and it writes on the
  !unit of that timer.
  !
  !The bar has 50 characters, and it is updated in place on the same line, followed by the ETA computed as in
  !:f:func_inline:`eta`: :code:`NNN% |*****     | ETA h:mm:ss.mmm`. The line is terminated by the last call, :code:`i==imax`.
  !The clock used is selected by the first character of :f:var:`method`: :code:`c` (default) for :code:`cpu_time`, :code:`d`
  !for :code:`date_and_time`, anything else, e.g. :code:`s`, for :code:`system_clock` (wall time).
  !
    integer                   :: i       ! current iteration
    integer                   :: imax    ! total number of iterations
    character(len=*),optional :: method  ! optional: clock, 'c' cpu_time (default), 'd' date_and_time, 's' system_clock
    integer                   :: k,jmax
    character(len=1)          :: method_
    character(len=62)         :: bar="???% |                                                  | ETA "
    !                            bar="100% |-------------50char-------------------------------| ETA "
    !
    method_='c';if(present(method))method_=trim(method) !s=system_clock(wall_time),c=cpu_time,d=date_and_time
    !
    old_time     = time
    time         = total_time(method_)
    !Get ETA
    dtime        = time-old_time     
    elapsed_time = elapsed_time + dtime
    dtime        = elapsed_time/real(i,4)
    eta_time     = dtime*imax - elapsed_time
    !
    ms           = int(fraction(eta_time)*1000.0)
    h            = int(eta_time/secs_in_one_hour)
    m            = int((eta_time - h*secs_in_one_hour)/secs_in_one_min)
    s            = int(eta_time - h*secs_in_one_hour - m*secs_in_one_min)
    !
    write(unit=bar(1:3),fmt="(I3,$)")100*i/imax
    jmax=50*i/imax
    do k=1,jmax
       bar(6+k:6+k)="*"
    enddo
    write(unit=funit(Tindex),fmt="(A1,A1,A62,I2,A1,I2.2,A1,I2.2,A1,I3.3,$)")'+',char(13), bar,h,":",m,":",s,".",ms
    if(i==imax)write(funit(Tindex),*)
  end subroutine progress











  !+-------------------------------------------------------------------+
  !PURPOSE  : a total_time
  !+-------------------------------------------------------------------+
  function total_time(method)
    character(len=*)        :: method
    real(8)                 :: total_time
    integer(4),dimension(8) :: dummy
    integer(8) :: T_
    real(8) :: cT_
    select case(trim(method))
    case default
#ifdef _MPI    
       if(check_MPI())then
          st_T1(Tindex) = MPI_Wtime()
       else
          call system_clock(count=T_)
          st_T1(Tindex)=T_
       endif
#else
       call system_clock(count=T_)
       st_T1(Tindex)=T_
#endif              
       total_time = st_time_difference(st_T1(Tindex),st_T0(Tindex))
    case("c","ct","cpu","cpu_time")
       call cpu_time(cT_)
       ct_T1(Tindex)=dble(cT_)
       total_time = ct_time_difference(ct_T1(Tindex),ct_T0(Tindex))
    case("d","dt","date_time","date","time","date_and_time")
       call date_and_time(values=dt_T1(Tindex,:))
       dummy= dt_time_difference(dt_T1(Tindex,:),dt_T0(Tindex,:))
       year = dummy(1)
       mese = dummy(2)
       day  = dummy(3)
       h    = dummy(5)
       m    = dummy(6)
       s    = dummy(7)
       ms   = dummy(8)
       total_time= dble(ms)/1000        &
            + dble(s)                   &
            + dble(m)*secs_in_one_min   &
            + dble(h)*secs_in_one_hour  &
            + dble(day)*secs_in_one_day
    end select
  end function total_time




  !+-------------------------------------------------------------------+
  !PURPOSE  : get time difference between two events
  !+-------------------------------------------------------------------+
  function ct_time_difference(data1,data0) result(time_difference)
    real(8) :: data1,data0
    real(8) :: time_difference
    time_difference =dble(data1-data0)
  end function ct_time_difference

  function st_time_difference(data1,data0) result(time_difference)
    real(8) :: data1,data0
    real(8) :: time_difference
    time_difference = dble(data1-data0)/st_rate
  end function st_time_difference

  function dt_time_difference(data1,data0) result(time_difference)
    integer(4),dimension(8) :: data1,data0,dummy,time_difference
    dummy =data1-data0
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    if(h<0)then
       day=day-1
       h=24+h
    endif
    if(m<0)then
       h=h-1
       m=60+m
    endif
    if(s<0)then
       m=m-1
       s=60+s
    endif
    if(ms<0)then
       s=s-1
       ms=1000+ms
    endif
    time_difference(1)=year
    time_difference(2)=mese
    time_difference(3)=day
    time_difference(5)=h
    time_difference(6)=m
    time_difference(7)=s
    time_difference(8)=ms
  end function dt_time_difference





  !+-------------------------------------------------------------------+
  !PURPOSE  : print total time
  !+-------------------------------------------------------------------+
  subroutine print_total_time(ct_T,st_T,dt_T,title)
    real(8)                   :: ct_T
    real(8)                   :: st_T
    integer                   :: dt_T(8)
    character(len=*),optional :: title
    !
    write(funit(Tindex),"(a)",advance="no")"CPU time (wall time): "
    ms=int(fraction(ct_T)*1000.0)
    h =int(ct_T/secs_in_one_hour)
    m =int((ct_T - h*secs_in_one_hour)/secs_in_one_min)
    s =int(ct_T  - h*secs_in_one_hour - m*secs_in_one_min)
    write(funit(Tindex),"(i3.3,a1,i2.2,a1,i2.2,a1,i3.3,3x)",advance="no")h,":",m,":",s,".",ms
    ms=int(fraction(real(st_T))*1000.0)
    h =int(st_T/secs_in_one_hour)
    m =int((st_T - h*secs_in_one_hour)/secs_in_one_min)
    s =int(st_T  - h*secs_in_one_hour - m*secs_in_one_min)
    write(funit(Tindex),"(a1,i3.3,a1,i2.2,a1,i2.2,a1,i3.3,a1)",advance="no")"(",h,":",m,":",s,".",ms,")"
    if(present(title))then
       write(funit(Tindex),"(A3,A)",advance="yes")" : ",trim(title)
    else
       write(funit(Tindex),"(1x)",advance="yes")
    endif
  end subroutine print_total_time



end module SF_TIMER






! !+-------------------------------------------------------------------+
! !TYPE     : Subroutine
! !+-------------------------------------------------------------------+
! subroutine progress(i,imax)
!   integer            :: i,imax,k,jmax
!   character(len=7)  :: bar="> ???% "
!   write(unit=bar(3:5),fmt="(I3,$)")100*i/imax
!   write(unit=funit,fmt="(A1,A1,A7,$)")'+',char(13), bar
!   if(i==imax)write(funit,*)
! end subroutine progress



! !+-------------------------------------------------------------------+
! !TYPE     : Subroutine
! !+-------------------------------------------------------------------+
! subroutine progress_bar(i,imax)
!   integer            :: i,imax,k,jmax
!   character(len=57)  :: bar="???% |                                                  |"
!   write(unit=bar(1:3),fmt="(I3,$)")100*i/imax
!   jmax=50*i/imax
!   do k=1,jmax
!      bar(6+k:6+k)="*"
!   enddo
!   write(unit=funit,fmt="(A1,A1,A57,$)")'+',char(13), bar
!   if(i==imax)write(funit,*)
! end subroutine progress_bar
