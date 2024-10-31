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
  integer,save                             :: Tindex=0
  integer,parameter                        :: Tmax=100!
  !
  real,save                                :: time
  real,save                                :: old_time
  real,save                                :: dtime
  real,save                                :: elapsed_time
  real,save                                :: eta_time
  integer(8),dimension(Tmax)               :: Funit
  integer(8)                               :: st_rate=1
  !
  !Using Date and Time:
  integer,dimension(Tmax,8),save           :: dt_T_start, dt_T_stop, dt_T0, dt_T1
  !Using Cpu_time
  real,dimension(Tmax),save                :: ct_T_start, ct_T_stop, ct_T0, ct_T1
  !Using System_clock
  integer(8),dimension(Tmax),save          :: st_T_start, st_T_stop, st_T0, st_T1

  integer                                  :: mpi_id
  integer                                  :: mpi_size
  logical                                  :: mpi_master

  interface start_progress
     module procedure :: start_timer
  end interface start_progress

  interface stop_progress
     module procedure :: stop_timer
  end interface stop_progress

  public :: start_timer,start_progress
  public :: stop_timer ,stop_progress
  public :: eta
  public :: progress



contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : start a timer to measure elapsed time between two call
  !+-------------------------------------------------------------------+
  subroutine start_timer(title,unit)
    character(len=*),optional :: title
    integer,optional          :: unit
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
       call system_clock(count=st_T_start(Tindex))
    endif
#else
    call system_clock(count_rate=st_rate)
    call system_clock(count=st_T_start(Tindex))
#endif
    call cpu_time(ct_T_start(Tindex))
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
    character(len=*),optional :: msg
    real                      :: ct_T
    real                      :: st_T
    integer,dimension(8)      :: dt_T
    !
    !
#ifdef _MPI    
    if(check_MPI())then
       st_T_stop(Tindex) = MPI_Wtime()
    else
       call system_clock(count=st_T_stop(Tindex))
    endif
#else
    call system_clock(count=st_T_stop(Tindex))
#endif
    call cpu_time(ct_T_stop(Tindex))
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
    ct_T_start(Tindex)   = 0
    st_T_start(Tindex)   = 0
    dt_T_start(Tindex,:) = 0
    ct_T_stop(Tindex)    = 0
    st_T_stop(Tindex)    = 0
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
    integer                   :: i
    integer                   :: L
    integer,optional          :: step
    character(len=*),optional :: method
    !
    integer,save              :: mod_print
    integer                   :: percent,iprint
    integer,save              :: older=0,oldiprint=0
    logical                   :: esc,fullprint
    integer(4),dimension(8)   :: dummy
    character(len=16)         :: string
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
          write(string,"(I4)")funit(Tindex)
          write(*,"(2x,A,I3)")"+ETA --> unit:"//trim(adjustl(trim(string)))
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
    integer                   :: i,imax
    character(len=*),optional :: method
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
    select case(trim(method))
    case default
#ifdef _MPI    
       if(check_MPI())then
          st_T1(Tindex) = MPI_Wtime()
       else
          call system_clock(count=st_T1(Tindex))
       endif
#else
       call system_clock(count=st_T1(Tindex))
#endif              
       total_time = st_time_difference(st_T1(Tindex),st_T0(Tindex))
    case("c","ct","cpu","cpu_time")
       call cpu_time(ct_T1(Tindex))
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
    real    :: data1,data0
    real(8) :: time_difference
    time_difference =dble(data1-data0)
  end function ct_time_difference

  function st_time_difference(data1,data0) result(time_difference)
    integer(8):: data1,data0
    real(8)   :: time_difference
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
    real                      :: ct_T
    real                      :: st_T
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
