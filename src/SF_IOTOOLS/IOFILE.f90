module IOFILE
  implicit none
  private

  !file size to be stored automagically (in Kb)
  integer,save :: store_size=2048

  interface str
     module procedure str_i_to_ch
     module procedure str_i_to_ch_pad
     module procedure str_r_to_ch
     module procedure str_c_to_ch
     module procedure str_l_to_ch
     module procedure str_ch_to_ch
  end interface str

  interface txtfy
     module procedure str_i_to_ch
     module procedure str_i_to_ch_pad
     module procedure str_r_to_ch
     module procedure str_c_to_ch
     module procedure str_l_to_ch
     module procedure str_ch_to_ch
  end interface txtfy

  interface reg
     module procedure reg_filename
  end interface reg

  interface create_dir
     module procedure create_data_dir
  end interface create_dir

  interface newunit
     module procedure free_unit
  end interface newunit

  interface print_matrix
     module procedure print_array_d
     module procedure print_array_c
  end interface print_matrix


  public :: set_store_size
  !
  public :: str
  public :: txtfy !obsolete
  public :: reg
  public :: to_lower
  public :: to_upper
  !
  public :: file_size
  public :: file_length
  public :: file_info
  !
  public :: file_gzip           !data_store
  public :: file_bzip           !data_store
  public :: file_xz           !data_store
  public :: file_targz
  public :: file_tarbz2
  !
  public :: file_gunzip         !data_open
  public :: file_bunzip         !data_open
  public :: file_unxz         !data_open
  public :: file_untargz
  public :: file_untarbz2
  !
  public :: newunit
  public :: free_unit
  public :: free_units
  !
  public :: create_dir
  !
  public :: get_filename
  public :: get_filepath
  !
  public :: print_matrix

contains




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  pure function reg_filename(file) result(reg)
    character(len=*),intent(in)                        :: file    
    character(len=len_trim(trim(adjustl(trim(file))))) :: reg
    reg=trim(adjustl(trim(file)))
  end function reg_filename

  function get_filename(string) result(fname)
    character(len=*) :: string
    ! character(len=len_trim(string)) :: fname
    character(len=:),allocatable :: fname
    integer :: i,slen
    slen=len_trim(string)
    do i=slen,1,-1
       if(string(i:i)== '/')exit
    enddo
    fname=trim(adjustl(trim(string(i+1:slen))))
  end function get_filename

  function get_filepath(string) result(pname)
    character(len=*) :: string
    ! character(len=len_trim(string)) :: pname
    character(len=:),allocatable :: pname
    integer :: i,slen
    slen=len_trim(string)
    do i=slen,1,-1
       if(string(i:i)== '/')exit
    enddo
    pname=trim(adjustl(trim(string(1:i))))
  end function get_filepath


  function free_unit(n) result(unit_)
    integer,optional :: n
    integer          :: unit_,ios
    logical          :: opened
    unit_=100
    do 
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=opened,iostat=ios)
       if(.not.opened.AND.ios==0)exit 
       if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
    enddo
    if(present(n))n=unit_
  end function free_unit

  function free_units(n) result(unit)
    integer :: n
    integer :: unit(n)
    integer :: i,unit_,ios
    logical :: is_it_opened
    unit_=100
    do i=1,n
       search: do 
          unit_=unit_+1
          inquire(unit=unit_,OPENED=is_it_opened,iostat=ios)
          if(.not.is_it_opened.AND.ios==0)exit search
          if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
       enddo search
       unit(i)=unit_
    enddo
  end function free_units




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_size(file,printf) result(size)
    integer               :: size,status
    character(len=*)      :: file
    integer,dimension(13) :: buff
    logical,optional      :: printf
    logical               :: control,printf_
    printf_=.false.;if(present(printf))printf_=printf
    inquire(file=reg(file),exist=control)
    if(.not.control)then
       write(*,*) 'Cannot read '//reg(file)//'. Skip file_size'
       return
    endif
    open(10,file=reg(file))
    call fstat(10,buff,status)
    size=nint(dble(buff(8))/dble(1024))
    if(printf_)write(*,"(A,A,A,f9.6,A)")"file: **",reg(file),"** is ",size," Kb"
  end function file_size



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_info(file)
    integer               :: file_info
    character(len=*)      :: file
    integer,dimension(13) :: buff
    integer               :: status,ifile
    logical               :: IOfile
    inquire(file=reg(file),exist=IOfile)
    if(.not.IOfile)then
       print*,'Cannot read ',reg(file),': skip file_size'
       file_info=0
       return
    endif
    open(10,file=reg(file))
    call fstat(10,buff,status)
    if(status == 0)then
       WRITE (*, FMT="('Device ID:',               T30, I19)") buff(1)
       WRITE (*, FMT="('Inode number:',            T30, I19)") buff(2)
       WRITE (*, FMT="('File mode (octal):',       T30, O19)") buff(3)
       WRITE (*, FMT="('Number of links:',         T30, I19)") buff(4)
       WRITE (*, FMT="('Owner''s uid:',            T30, I19)") buff(5)
       WRITE (*, FMT="('Owner''s gid:',            T30, I19)") buff(6)
       WRITE (*, FMT="('Device where located:',    T30, I19)") buff(7)
       WRITE (*, FMT="('File size:',               T30, I19)") buff(8)
       WRITE (*, FMT="('Last access time:',        T30, A19)") buff(9)
       WRITE (*, FMT="('Last modification time',   T30, A19)") buff(10)
       WRITE (*, FMT="('Last status change time:', T30, A19)") buff(11)
       WRITE (*, FMT="('Preferred block size:',    T30, I19)") buff(12)
       WRITE (*, FMT="('No. of blocks allocated:', T30, I19)") buff(13)
    endif
    close(10)
  end function file_info



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_length(file,verbose) result(lines)
    integer           :: lines
    logical,optional  :: verbose
    character(len=*)  :: file
    integer           :: ifile,ierr,pos
    logical           :: IOfile,bool,bool1,bool2,verbose_
    character(len=256)::buffer
    verbose_=.true.;if(present(verbose))verbose_=verbose
    inquire(file=reg(file),exist=IOfile)
    if(.not.IOfile)then
       inquire(file=reg(file)//".gz",exist=IOfile)
       if(IOfile)call file_gunzip(reg(file))
    endif
    if(.not.IOfile)then
       inquire(file=reg(file)//".bz2",exist=IOfile)
       if(IOfile)call file_bunzip(reg(file))
    endif
    if(.not.IOfile)then
       inquire(file=reg(file)//".xz",exist=IOfile)
       if(IOfile)call file_unxz(reg(file))
    endif
    lines=0
    if(.not.IOfile)then
       write(*,*) 'Cannot read '//reg(file)//'. Skip file_length'
       return
    endif
    open(99,file=reg(file))
    ierr=0
    do while(ierr==0)
       lines=lines+1
       read(99,*,iostat=ierr)buffer
       bool1=scan(buffer,"#").ne.0
       bool2=len_trim(buffer).eq.0       
       if(bool1 .OR. bool2)lines=lines-1
    enddo
    lines=lines-1
    if(verbose_)write(*,'(A,I9,A)') 'there are', lines,' lines in '//reg(file)
    rewind(99)
    close(99)
  end function file_length



  !******************************************************************
  !******************************************************************
  !******************************************************************







  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine set_store_size(size)
    integer :: size
    store_size=size
    write(*,"(A)")"store size ="//trim(txtfy(size))//"Kb"
  end subroutine set_store_size




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_gzip(file,size)
    character(len=*)  :: file
    integer,optional  :: size
    logical           :: control
    character(len=9)  :: csize 
    integer           :: cstatus,fsize,unit,len
    fsize=store_size;if(present(size))fsize=size
    write(*,"(A)") "Storing "//file
    inquire(file=reg(file),exist=control)
    if(control)then
       len=file_size(reg(file))
       if(len>fsize)call system("gzip -fv --best --rsyncable"//reg(file))
    endif
    inquire(file=reg(file),opened=control,number=unit)
    if(control)close(unit)
  end subroutine file_gzip





  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_gunzip(filename)
    character(len=*)             :: filename
    character(len=3),parameter   :: type='.gz'
    logical                      :: iexist,iopen
    integer                      :: unit
    !
    inquire(file=reg(filename),exist=iexist)
    if(iexist)return           !nothing to be done:
    !
    inquire(file=reg(filename)//type,exist=iexist)
    if(.not.iexist)then
       write(*,"(A)")"file "//reg(filename)//" not found, not even with .gz extension"
       stop
    endif
    write(*,"(A)")"deflate "//reg(filename)//type
    call system("gunzip -l"//reg(filename)//type)
    inquire(file=reg(filename),opened=iopen,number=unit)
    if(iopen)close(unit)
  end subroutine file_gunzip









  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_bzip(file,size)
    character(len=*)  :: file
    integer,optional  :: size
    logical           :: control
    character(len=9)  :: csize 
    integer           :: cstatus,fsize,unit,len
    fsize=store_size;if(present(size))fsize=size
    write(*,"(A)") "Storing "//file
    inquire(file=reg(file),exist=control)
    if(control)then
       len=file_size(reg(file))
       if(len>fsize)call system("bzip2 -zfv "//reg(file))
    endif
    inquire(file=reg(file),opened=control,number=unit)
    if(control)close(unit)
  end subroutine file_bzip


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_bunzip(filename)
    character(len=*)             :: filename
    character(len=4),parameter   :: type='.bz2'
    logical                      :: iexist,iopen
    integer                      :: unit
    !
    inquire(file=reg(filename),exist=iexist)
    if(iexist)return           !nothing to be done:
    !
    inquire(file=reg(filename)//type,exist=iexist)
    if(.not.iexist)then
       write(*,"(A)")"file "//reg(filename)//" not found, not even with .bz2 extension"
       stop
    endif
    write(*,"(A)")"deflate "//reg(filename)//type
    call system("bzip2 -dv "//reg(filename)//type)
    inquire(file=reg(filename),opened=iopen,number=unit)
    if(iopen)close(unit)
  end subroutine file_bunzip








  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_xz(file,size)
    character(len=*)  :: file
    integer,optional  :: size
    logical           :: control
    character(len=9)  :: csize 
    integer           :: cstatus,fsize,unit,len
    fsize=store_size;if(present(size))fsize=size
    write(*,"(A)") "Storing "//file
    inquire(file=reg(file),exist=control)
    if(control)then
       len=file_size(reg(file))
       if(len>fsize)call system("xz -zfv "//reg(file))
    endif
    inquire(file=reg(file),opened=control,number=unit)
    if(control)close(unit)
  end subroutine file_xz


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_unxz(filename)
    character(len=*)             :: filename
    character(len=3),parameter   :: type='.xz'
    logical                      :: iexist,iopen
    integer                      :: unit
    !
    inquire(file=reg(filename),exist=iexist)
    if(iexist)return           !nothing to be done:
    !
    inquire(file=reg(filename)//type,exist=iexist)
    if(.not.iexist)then
       write(*,"(A)")"file "//reg(filename)//" not found, not even with .xz extension"
       stop
    endif
    write(*,"(A)")"deflate "//reg(filename)//type
    call system("xz -dv "//reg(filename)//type)
    inquire(file=reg(filename),opened=iopen,number=unit)
    if(iopen)close(unit)
  end subroutine file_unxz








  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_targz(tarball,pattern,size)
    character(len=*)           :: tarball
    character(len=*)           :: pattern
    integer,optional           :: size
    character(len=4),parameter :: type='.tgz'
    integer                    :: control,fsize
    character(len=100)         :: cmsg
    fsize=store_size;if(present(size))fsize=size
    write(*,"(A)") "Store "//str(pattern)//" in "//str(tarball)//type
    call execute_command_line("tar -czf "//str(tarball)//type//" "//str(pattern),CMDSTAT=control,CMDMSG=cmsg)
    if(control>0)then
       write(*,"(A)")"Command tar -czf failed with error: "//str(cmsg)
    elseif(control<0)then
       write(*,"(A)")"Command tar -czf not supported"
    else
       call execute_command_line("rm -f "//str(pattern))
    endif
  end subroutine file_targz


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_untargz(tarball)
    character(len=*)           :: tarball
    integer                    :: control
    logical                    :: iexist,iopen
    integer                    :: unit
    character(len=4),parameter :: type='.tgz'
    character(len=100)         :: cmsg
    inquire(file=reg(tarball)//type,exist=iexist)
    if(.not.iexist)then
       write(*,"(A)")"Tarball "//str(tarball)//type//" not present: skip"
       return           !nothing to be done:
    endif
    write(*,"(A)")"deflate "//reg(tarball)//type
    call execute_command_line("tar -xzf "//str(tarball)//type//" ",CMDSTAT=control,CMDMSG=cmsg)
    if(control>0)then
       write(*,"(A)")"Command tar -xzf failed with error: "//str(cmsg)
    elseif(control<0)then
       write(*,"(A)")"Command tar -xzf not supported"
    else
       call execute_command_line("rm -f "//str(tarball)//type)
    endif
  end subroutine file_untargz










  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_tarbz2(tarball,pattern,size)
    character(len=*)           :: tarball
    character(len=*)           :: pattern
    integer,optional           :: size
    character(len=8),parameter :: type='.tar.bz2'
    integer                    :: control,fsize
    character(len=100)         :: cmsg
    fsize=store_size;if(present(size))fsize=size
    write(*,"(A)") "Store "//str(pattern)//" in "//str(tarball)//type
    call execute_command_line("tar -cjSf "//str(tarball)//type//" "//str(pattern),CMDSTAT=control,CMDMSG=cmsg)
    if(control>0)then
       write(*,"(A)")"Command tar -cjSf failed with error: "//str(cmsg)
    elseif(control<0)then
       write(*,"(A)")"Command tar -cjSf not supported"
    else
       call execute_command_line("rm -f "//str(pattern))
    endif
  end subroutine file_tarbz2


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine file_untarbz2(tarball)
    character(len=*)           :: tarball
    integer                    :: control
    logical                    :: iexist,iopen
    integer                    :: unit
    character(len=8),parameter :: type='.tar.bz2'
    character(len=100)         :: cmsg
    inquire(file=reg(tarball)//type,exist=iexist)
    if(.not.iexist)then
       write(*,"(A)")"Tarball "//str(tarball)//type//" not present: skip"
       return           !nothing to be done:
    endif
    write(*,"(A)")"deflate "//reg(tarball)//type
    call execute_command_line("tar -xjf "//str(tarball)//type//" ",CMDSTAT=control,CMDMSG=cmsg)
    if(control>0)then
       write(*,"(A)")"Command tar -xjf failed with error: "//str(cmsg)
    elseif(control<0)then
       write(*,"(A)")"Command tar -xjf not supported"
    else
       call execute_command_line("rm -f "//str(tarball)//type)
    endif
  end subroutine file_untarbz2





  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine create_data_dir(dir_name)
    character(len=*) :: dir_name
    logical          :: control
    call system("mkdir -v "//reg(dir_name))
  end subroutine create_data_dir


  !******************************************************************
  !******************************************************************
  !******************************************************************



  function to_upper(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("a":"z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))-32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_upper

  function to_lower(StrIn) result(StrOut)
    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer :: i
    do i = 1,len(StrIn)
       select case(StrIn(i:i))
       case("A":"Z")
          StrOut(i:i) = achar(iachar(StrIn(i:i))+32)
       case default
          StrOut(i:i) = StrIn(i:i)
       end select
    end do
  end function to_lower




  ! function i_to_ch(i4) result(string)
  !   character(len=32) :: string
  !   integer           :: i4
  !   call i4_to_s_left(i4,string)
  ! end function i_to_ch
  ! function i_to_ch_pad(i4,Npad) result(string)
  !   integer             :: Npad
  !   character(len=Npad) :: string
  !   integer             :: i4
  !   call i4_to_s_zero(i4,string)
  ! end function i_to_ch_pad
  ! function r_to_ch(r8) result(string)
  !   character(len=32) :: string
  !   character(len=16) :: string_
  !   real(8)           :: r8
  !   call r8_to_s_left(r8,string_)
  !   string=adjustl(string_)
  ! end function r_to_ch
  ! function c_to_ch(c) result(string)
  !   character(len=32+3) :: string
  !   character(len=16) :: sre,sim
  !   complex(8)        :: c
  !   real(8)           :: re,im
  !   re=real(c,8);im=aimag(c)
  !   call r8_to_s_left(re,sre)
  !   call r8_to_s_left(im,sim)
  !   string="("//trim(sre)//","//trim(sim)//")"
  ! end function c_to_ch
  ! function l_to_ch(bool) result(string)
  !   logical :: bool
  !   character(len=1) :: string
  !   string='F'
  !   if(bool)string='T'
  ! end function l_to_ch
  ! function ch_to_ch(txt) result(string)
  !   character(len=*)              :: txt
  !   character(len=len(trim(txt))) :: string
  !   string=trim(txt)
  ! end function ch_to_ch



  function str_i_to_ch(i4) result(string)
    integer                      :: i4
    character(len=:),allocatable :: string
    character(len=16)            :: string_
    call i4_to_s_left(i4,string_)
    string=trim(adjustl(trim(string_)))
  end function str_i_to_ch

  function str_i_to_ch_pad(i4,Npad) result(string)
    integer                      :: i4
    integer                      :: Npad
    character(len=:),allocatable :: string
    character(len=Npad)          :: string_pad
    call i4_to_s_zero(i4,string_pad)
    string=trim(adjustl(trim(string_pad)))
  end function str_i_to_ch_pad

  function str_r_to_ch(r8,d,lead) result(string)
    real(8)                      :: r8
    integer,optional             :: d,lead
    integer                      :: w_,d_,lead_
    character(len=:),allocatable :: string
    character(len=:),allocatable :: string_
    d_    =6 ;if(present(d))d_=d
    lead_ =1 ;if(present(lead))lead_=lead
    w_ = get_w_(r8,d_,lead_)
    allocate(character(len=w_) :: string_)
    call r8_to_s_left(r8,string_,d_,lead_)
    string=trim(adjustl(trim(string_)))
  end function str_r_to_ch


  function str_c_to_ch(c,d,lead) result(string)
    complex(8)                   :: c
    integer,optional             :: d,lead
    integer                      :: w_,d_,lead_
    character(len=:),allocatable :: string
    character(len=:),allocatable :: sre,sim
    real(8)                      :: re,im
    d_    =6 ;if(present(d))d_=d
    lead_ =1 ;if(present(lead))lead_=lead
    re=dreal(c)
    w_ = get_w_(re,d_,lead_)
    allocate(character(len=w_) :: sre)
    call r8_to_s_left(re,sre,d_,lead_)
    !
    im=dimag(c)
    w_ = get_w_(im,d_,lead_)
    allocate(character(len=w_) :: sim)
    call r8_to_s_left(im,sim,d_,lead_)
    string="("//trim(adjustl(trim(sre)))//","//trim(adjustl(trim(sim)))//")"
  end function str_c_to_ch


  function str_l_to_ch(bool) result(string)
    logical          :: bool
    character(len=1) :: string
    string='F'
    if(bool)string='T'
  end function str_l_to_ch

  function str_ch_to_ch(txt) result(string)
    character(len=*)                             :: txt
    character(len=:),allocatable :: string
    string=trim(adjustl(trim(txt)))
  end function str_ch_to_ch


  function get_w_(r8,d,lead) result(w)
    real(8) :: r8
    integer :: d,lead
    integer :: w
    if(r8==0d0)then
       w=d+4
    else
       w=floor(log10(abs(r8)))
       if(w < -lead)then
          w = d + 4 + 4
       else
          w = abs(w) + d + 4 + lead
       endif
    endif
  end function get_w_




  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine i4_to_s_left ( i4, s )
    !! I4_TO_S_LEFT converts an I4 to a left-justified string.
    !  Example:
    !    Assume that S is 6 characters long:
    !        I4  S
    !         1  1
    !        -1  -1
    !         0  0
    !      1952  1952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Parameters:
    !    Input, integer ( kind = 4 ) I4, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be left-justified.  If there is not enough space,
    !    the string will be filled with stars.
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


  subroutine r8_to_s_left( r8, s, digits, lead)
    !! R8_TO_S_LEFT writes an R8 into a left justified string.
    !    An R8 is a real ( kind = 8 ) value.
    !    A 'F<len(s)>.DIGITS' format is used with a WRITE statement.
    character(len=12)   :: fmt
    integer             :: i
    real(8)             :: r8
    character(len=*)    :: s
    integer             :: s_length,w_
    integer             :: digits,lead
    s_length = len ( s )
    write(fmt,"(A2,I0,A1,I0,A1)")"(F",s_length,".",digits,")"
    if(r8/=0d0)then
       w_=floor(log10(abs(r8)))
       if(w_<-lead)then
          write(fmt,"(A3,I0,A1,I0,A1)")"(ES",s_length,".",digits,")"
       else
          write(fmt,"(A2,I0,A1,I0,A1)")"(F",s_length,".",digits+lead,")"
       endif
    endif
    write ( s, fmt ) r8
    s = trim(adjustl(trim( s )))
  end subroutine r8_to_s_left


  subroutine digit_to_ch(digit,ch)
    !! DIGIT_TO_CH returns the character representation of a decimal digit.
    !    Instead of CHAR, we now use the ACHAR function, which
    !    guarantees the ASCII collating sequence.
    !  Example:
    !    DIGIT   CH
    !    -----  ---
    !      0    '0'
    !      1    '1'
    !    ...    ...
    !      9    '9'
    !     17    '*'
    !  Parameters:
    !    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
    !    Output, character CH, the corresponding character.
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch


  subroutine i4_to_s_zero ( intval, s )
    !! I4_TO_S_ZERO converts an I4 to a string, with zero padding.
    !    An I4 is an integer ( kind = 4 ).
    !  Example:
    !    Assume that S is 6 characters long:
    !    INTVAL  S
    !         1  000001
    !        -1  -00001
    !         0  000000
    !      1952  001952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Parameters:
    !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be right justified, and zero padded.
    !    If there is not enough space, the string will be filled with stars.
    implicit none
    character c
    integer ( kind = 4 ) i
    integer ( kind = 4 ) idig
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) intval
    integer ( kind = 4 ) ipos
    integer ( kind = 4 ) ival
    character ( len = * ) s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !
    !  Make a copy of the integer.
    !
    ival = intval
    !
    !  Handle the negative sign.
    !
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !
    !  Working from right to left, strip off the digits of the integer
    !  and place them into S(ILO:IHI).
    !
    ipos = ihi
    do while ( ival /= 0 .or. ipos == ihi )
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
    end do
    !
    !  Fill the empties with zeroes.
    !
    do i = ilo, ipos
       s(i:i) = '0'
    end do
    return
  end subroutine i4_to_s_zero






  subroutine print_array_d(M,file,w,d)
    real(8),dimension(:,:)    :: M
    character(len=*),optional :: file
    integer,optional          :: w,d
    integer                   :: w_,d_
    integer                   :: i,j
    integer                   :: unit
    character(len=128)        :: fmt
    unit=6
    if(present(file))open(free_unit(unit),file=reg(file))
    w_=5;if(present(w))w_=w
    d_=2;if(present(d))d_=d
    write(fmt,"(A1,I0,A2,I0,A1,I0,A5)")"(",size(M,2),"(F",w_,".",d_,",1x))"
    do i=1,size(M,1)
       write(unit,fmt)(M(i,j),j=1,size(M,2))
    enddo
  end subroutine print_array_d

  subroutine print_array_c(M,file,w,d)
    complex(8),dimension(:,:) :: M
    character(len=*),optional :: file
    integer,optional          :: w,d
    integer                   :: w_,d_
    integer                   :: i,j
    integer                   :: unit
    character(len=128)        :: fmt
    unit=6
    if(present(file))open(free_unit(unit),file=reg(file))
    w_=5;if(present(w))w_=w
    d_=2;if(present(d))d_=d
    write(fmt,"(A1,I0,A5,I0,A1,I0,A5,I0,A1,I0,A7)")"(",size(M,2),"(A1,F",w_,".",d_,",A1,F",w_,".",d_,"A1,1x))"
    do i=1,size(M,1)
       write(unit,fmt)("(",dreal(M(i,j)),",",dimag(M(i,j)),")",j=1,size(M,2))
    enddo
  end subroutine print_array_c


end module IOFILE




