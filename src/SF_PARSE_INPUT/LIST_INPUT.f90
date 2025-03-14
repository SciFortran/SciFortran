module LIST_INPUT
  implicit none
  private

  type input_var
     integer,pointer             :: i
     real(8),pointer             :: d
     logical,pointer             :: l
     character(len=:),pointer    :: ch
  end type input_var

  type input_node
     type(input_var),dimension(:),allocatable :: var
     character(len=3)                         :: type
     character(len=100)                       :: name
     character(len=512)                       :: comment
     type(input_node),pointer                 :: next !link to next box
  end type input_node

  type input_list
     logical                   :: status=.false.
     integer                   :: size 
     type(input_node),pointer  :: root
  end type input_list


  interface append_to_input_list
     module procedure i_append_to_input_list
     module procedure d_append_to_input_list
     module procedure l_append_to_input_list
     module procedure iv_append_to_input_list
     module procedure dv_append_to_input_list
     module procedure lv_append_to_input_list
     module procedure ch_append_to_input_list
  end interface append_to_input_list


  public :: input_list
  public :: init_input_list
  public :: delete_input_list
  public :: size_input_list
  public :: append_to_input_list
  public :: print_input_list


  type(input_list)   :: default_list

  character(len=255) :: p_buffer
  character(len=7)   :: file_status
  integer,parameter  :: pos_comment=46 !72

  !LOCAL VERSION OF TXTFY//STR
  interface txtfy
     module procedure i_to_ch,r_to_ch,c_to_ch,l_to_ch
  end interface txtfy



contains  





  !+------------------------------------------------------------------+
  !PURPOSE: init the input list
  !+------------------------------------------------------------------+
  subroutine init_input_list(list)
    type(input_list),optional :: list
    if(present(list))then
       allocate(list%root)    
       list%size=0
       list%status=.true.
       list%root%next=>null()
    else
       allocate(default_list%root)    
       default_list%size=0
       default_list%status=.true.
       default_list%root%next=>null()
    endif
  end subroutine init_input_list





  !+------------------------------------------------------------------+
  !PURPOSE: delete the list
  !+------------------------------------------------------------------+
  subroutine delete_input_list(list)
    type(input_list),optional :: list
    type(input_node),pointer  :: p,c
    integer :: i
    if(present(list))then
       do
          p => list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          do i=1,size(c%var)
             nullify(c%var(i)%i)
             nullify(c%var(i)%d)
             nullify(c%var(i)%l)
             nullify(c%var(i)%ch)
          enddo
          deallocate(c%var)
          deallocate(c)
       end do
       list%status=.false.
       deallocate(list%root)
    else
       do
          p => default_list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          do i=1,size(c%var)
             nullify(c%var(i)%i)
             nullify(c%var(i)%d)
             nullify(c%var(i)%l)
             nullify(c%var(i)%ch)
          enddo
          deallocate(c%var)
          deallocate(c)
       end do
       default_list%status=.false.
       deallocate(default_list%root)
    endif
    p => null()
    c => null()
  end subroutine delete_input_list




  !+------------------------------------------------------------------+
  !PURPOSE: get list size
  !+------------------------------------------------------------------+
  function size_input_list(list) result(size)
    type(input_list),optional :: list
    integer                   :: size
    size=default_list%size
    if(present(list))size=list%size
  end function size_input_list





  !+------------------------------------------------------------------+
  !PURPOSE: !Append input data to the list:
  !+------------------------------------------------------------------+
  !========================SCALAR==================================
  subroutine i_append_to_input_list(variable,name,comment)
    integer,target            :: variable
    character(len=*)          :: name
    character(len=*),optional :: comment
    type(input_node),pointer  :: p,c
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    allocate(p%next%var(1))
    p%next%var(1)%i=>variable
    p%next%name= name
    p%next%type='i'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine i_append_to_input_list

  subroutine d_append_to_input_list(variable,name,comment)
    real(8),target           :: variable
    character(len=*)         :: name
    character(len=*),optional:: comment
    type(input_node),pointer :: p,c
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    allocate(p%next%var(1))
    p%next%var(1)%d=>variable
    !
    p%next%name= name
    p%next%type='d'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine d_append_to_input_list

  subroutine l_append_to_input_list(variable,name,comment)
    logical,target           :: variable
    character(len=*)         :: name
    character(len=*),optional:: comment
    type(input_node),pointer :: p,c
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    ! allocate(p%next%l(1))
    ! p%next%l(1) = variable
    !>NEW
    allocate(p%next%var(1))
    p%next%var(1)%l=>variable
    !<
    p%next%name= name
    p%next%type='l'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine l_append_to_input_list


  !========================VECTOR==================================
  subroutine iv_append_to_input_list(variable,name,comment)
    integer,dimension(:),target :: variable
    character(len=*)            :: name
    character(len=*),optional   :: comment
    type(input_node),pointer    :: p,c
    integer :: i
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    ! allocate(p%next%i(size(variable)))
    ! p%next%i   = variable
    !>NEW
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%i=>variable(i)
    enddo
    !<
    p%next%name= name
    p%next%type='i'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine iv_append_to_input_list

  subroutine dv_append_to_input_list(variable,name,comment)
    real(8),dimension(:),target :: variable
    character(len=*)            :: name
    character(len=*),optional   :: comment
    type(input_node),pointer    :: p,c
    integer :: i
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    ! allocate(p%next%d(size(variable)))
    ! p%next%d   = variable
    !>NEW
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%d=>variable(i)
    enddo
    !<
    p%next%name= name
    p%next%type='d'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine dv_append_to_input_list

  subroutine lv_append_to_input_list(variable,name,comment)
    logical,dimension(:),target :: variable
    character(len=*)            :: name
    character(len=*),optional   :: comment
    type(input_node),pointer    :: p,c
    integer :: i
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    ! allocate(p%next%l(size(variable)))
    ! p%next%l   = variable
    !>NEW
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%l=>variable(i)
    enddo
    !<
    p%next%name= name
    p%next%type='l'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine lv_append_to_input_list



  !========================STRING==================================
  subroutine ch_append_to_input_list(variable,name,comment)
    character(len=*),target  :: variable
    character(len=*)         :: name
    character(len=*),optional:: comment
    type(input_node),pointer :: p,c
    if(.not.default_list%status)call init_input_list()
    p => default_list%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    !
    ! allocate(p%next%ch(1))
    ! p%next%ch(1) = variable
    !>NEW
    allocate(p%next%var(1))
    nullify(p%next%var(1)%ch)
    p%next%var(1)%ch=> variable
    !<
    p%next%name= name
    p%next%type='ch'
    p%next%comment=""
    if(present(comment))p%next%comment=trim(comment)
    !
    default_list%size=default_list%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine ch_append_to_input_list







  !+------------------------------------------------------------------+
  !PURPOSE: print the list to file
  !+------------------------------------------------------------------+
  subroutine print_input_list(file,list)
    character(len=*),optional :: file
    type(input_list),optional :: list
    integer                   :: i,counter,size
    type(input_node),pointer  :: c
    if(present(list))then
       c => list%root%next
    else
       c => default_list%root%next
    endif
    counter = 0 
    file_status='replace'
    size=default_list%size
    if(present(list))size=list%size
    if(size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(present(file))then
             call print_input_node(c,file)
          else
             call print_input_node(c)
          endif
          c => c%next
       enddo
    else
       write(*,*)"input list: empty"
       return
    endif
    file_status='replace'
    c => null()
  end subroutine print_input_list
  !---------------------------------------------------------------------
  subroutine print_input_node(c,file)
    type(input_node)          :: c
    character(len=*),optional :: file
    character(len=255)        :: blank=""
    integer                   :: clen
    integer                   :: unit,i
    !
    call s_blank_delete(c%name)
    select case(c%type)
    case('ch')
       p_buffer=trim(c%name)//"="//trim(adjustl(trim(c%var(1)%ch)))

    case('i')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%i)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%i))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%i))
       endif

    case('d')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%d)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%d))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%d))
       endif

    case('l')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%l)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%l))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%l))
       endif
    end select
    !
    call s_blank_delete(p_buffer)
    clen=pos_comment-len(trim(p_buffer))
    if(clen<=0)clen=1
    p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(c%comment)
    !
    ! write(*,"(1x,A)")trim(p_buffer)
    if(present(file))then
       unit=free_unit()
       open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
       write(unit,"(1x,A)")trim(p_buffer)
       close(unit)
    else
       write(*,"(1x,A)")trim(p_buffer)
    endif
    p_buffer=""
  end subroutine print_input_node






















  !+------------------------------------------------------------------+
  !PURPOSE: ANCILLARY routines
  !+------------------------------------------------------------------+
  !Auxiliary routines:
  subroutine upper_case(s)
    character              ch
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       ch = s(i:i)
       call ch_cap ( ch )
       s(i:i) = ch
    end do
  end subroutine upper_case

  subroutine lower_case(s)
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       call ch_low ( s(i:i) )
    end do
  end subroutine lower_case

  subroutine ch_cap(ch)
    character              ch
    integer   ( kind = 4 ) itemp
    itemp = iachar ( ch )
    if ( 97 <= itemp .and. itemp <= 122 ) then
       ch = achar ( itemp - 32 )
    end if
  end subroutine ch_cap

  subroutine ch_low ( ch )
    character ch
    integer ( kind = 4 ) i
    i = iachar ( ch )
    if ( 65 <= i .and. i <= 90 ) then
       ch = achar ( i + 32 )
    end if
  end subroutine ch_low


  subroutine s_blank_delete ( s )
    !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
    !    All TAB characters are also removed.
    !    Input/output, character ( len = * ) S, the string to be transformed.
    implicit none
    character              ch
    integer   ( kind = 4 ) get
    integer   ( kind = 4 ) put
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    character, parameter :: tab = achar ( 9 )
    put = 0
    s_length = len_trim ( s )
    do get = 1, s_length
       ch = s(get:get)
       if ( ch /= ' ' .and. ch /= tab ) then
          put = put + 1
          s(put:put) = ch
       end if
    end do
    s(put+1:s_length) = ' '
    return
  end subroutine s_blank_delete


  function i_to_ch(i4) result(string)
    character(len=32) :: string
    integer           :: i4
    call i4_to_s_left(i4,string)
  end function i_to_ch

  function r_to_ch(r8) result(string)
    character(len=32) :: string
    character(len=16) :: string_
    real(8)           :: r8
    call r8_to_s_left(r8,string_)
    string=adjustl(string_)
  end function r_to_ch

  function c_to_ch(c) result(string)
    character(len=32+3) :: string
    character(len=16) :: sre,sim
    complex(8)        :: c
    real(8)           :: re,im
    re=real(c,8);im=aimag(c)
    call r8_to_s_left(re,sre)
    call r8_to_s_left(im,sim)
    string="("//trim(sre)//","//trim(sim)//")"
  end function c_to_ch

  function l_to_ch(bool) result(string)
    logical :: bool
    character(len=1) :: string
    string='F'
    if(bool)string='T'
  end function l_to_ch

  subroutine i4_to_s_left ( i4, s )
    character :: c
    integer   :: i
    integer   :: i4
    integer   :: idig
    integer   :: ihi
    integer   :: ilo
    integer   :: ipos
    integer   :: ival
    character(len=*) ::  s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !  Make a copy of the integer.
    ival = i4
    !  Handle the negative sign.
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !  The absolute value of the integer goes into S(ILO:IHI).
    ipos = ihi
    !  Find the last digit of IVAL, strip it off, and stick it into the string.
    do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
          exit
       end if
    end do
    !  Shift the string to the left.
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '
  end subroutine i4_to_s_left

  subroutine r8_to_s_left ( r8, s )
    integer :: i
    real(8) :: r8
    integer :: s_length
    character(len=*) ::  s
    character(len=16) :: s2
    s_length = len ( s )
    if ( s_length < 16 ) then
       do i = 1, s_length
          s(i:i) = '*'
       end do
    else if ( r8 == 0.0D+00 ) then
       s(1:16) = '     0.d0     '
    else
       if(abs(r8) < 1.d0)then
          write ( s2, '(ES16.9)' ) r8
       else
          write ( s2, '(F16.9)' ) r8
       endif
       s(1:16) = s2
    end if
    !  Shift the string left.
    s = adjustl ( s )
  end subroutine r8_to_s_left


  subroutine digit_to_ch(digit,ch)
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch


  function free_unit() result(unit_)
    integer :: unit_,ios
    logical :: is_it_opened
    unit_=100
    do 
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=is_it_opened,iostat=ios)
       if(.not.is_it_opened.AND.ios==0)return 
       if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
    enddo
  end function free_unit



end module LIST_INPUT















! !========================POLYMORPHIC==================================
! subroutine poly_append_to_input_list(variable,name,comment)
!   class(*),target           :: variable
!   character(len=*)          :: name
!   character(len=*),optional :: comment
!   type(input_node),pointer  :: p,c
!   if(.not.default_list%status)call init_input_list()
!   p => default_list%root
!   c => p%next
!   do                            !traverse the list until obj < value (ordered list)
!      if(.not.associated(c))exit !empty list or beginning of the list
!      p => c
!      c => c%next
!   end do
!   allocate(p%next)                !Create a new element in the list
!   !
!   allocate(p%next%var(1))
!   p%next%var(1)%item  => variable
!   p%next%name= name
!   p%next%type='i'
!   p%next%comment=""
!   if(present(comment))p%next%comment=trim(comment)
!   !
!   default_list%size=default_list%size+1
!   if(.not.associated(c))then !end of the list special case (current=>current%next)
!      p%next%next  => null()
!   else
!      p%next%next  => c      !the %next of the new node come to current
!   end if
!   p=>null()
!   c=>null()
! end subroutine poly_append_to_input_list

!            
! DISABLED 
!
! !+------------------------------------------------------------------+
! !PURPOSE:   !Get input variable from the list:
! !+------------------------------------------------------------------+
! !========================0-dimension==================================
! subroutine i_get_input_variable(variable,name,list)
!   integer                   :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%i(1)
!            !>NEW
!            variable = c%var(1)%i
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine i_get_input_variable

! subroutine d_get_input_variable(variable,name,list)
!   real(8)                   :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%d(1)
!            !>NEW
!            variable = c%var(1)%d
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine d_get_input_variable

! subroutine l_get_input_variable(variable,name,list)
!   logical                   :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%l(1)
!            !>NEW
!            variable = c%var(1)%l
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine l_get_input_variable


! !========================1-dimension==================================
! subroutine iv_get_input_variable(variable,name,list)
!   integer,dimension(:)      :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%i(1:size(variable))
!            if(size(variable)/=size(c%var))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
!            !>NEW
!            do i=1,size(variable)
!               variable(i) = c%var(i)%i
!            enddo
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine iv_get_input_variable

! subroutine dv_get_input_variable(variable,name,list)
!   real(8),dimension(:)      :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%d(1:size(variable))
!            if(size(variable)/=size(c%var))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
!            !>NEW
!            do i=1,size(variable)
!               variable(i) = c%var(i)%d
!            enddo
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine dv_get_input_variable

! subroutine lv_get_input_variable(variable,name,list)
!   logical,dimension(:)      :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%l(1:size(variable))
!            if(size(variable)/=size(c%var))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
!            !>NEW
!            do i=1,size(variable)
!               variable(i) = c%var(i)%l
!            enddo
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine lv_get_input_variable


! !========================STRING==================================
! subroutine ch_get_input_variable(variable,name,list)
!   character(len=*)          :: variable
!   character(len=*)          :: name
!   type(input_list),optional :: list
!   integer                   :: i,counter,unit,size_
!   type(input_node),pointer  :: c
!   logical                   :: bool
!   character(len=len(name)) :: name_
!   name_=name;call upper_case(name_)
!   if(present(list))then
!      c => list%root%next
!   else
!      c => default_list%root%next
!   endif
!   counter = 0 
!   unit=free_unit()
!   size_=default_list%size
!   if(present(list))size_=list%size
!   if(size_>0)then
!      do
!         if(.not.associated(c))exit
!         counter=counter+1
!         if(trim(c%name)==trim(name_))then
!            ! variable=c%ch(1)
!            !>NEW
!            variable = c%var(1)%ch
!            !<
!            c=>null()
!            return
!         endif
!         c => c%next
!      enddo
!      write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
!   else
!      write(*,"(A)")"input list: empty"
!      return
!   endif
!   c => null()
! end subroutine ch_get_input_variable
