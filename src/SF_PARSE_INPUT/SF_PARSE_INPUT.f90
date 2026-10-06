module SF_PARSE_INPUT
!SciFortran module for input file parsing
  USE LIST_INPUT, delete_input=>delete_input_list
  implicit none
  private

  !cmd line variables:
  !=========================================================
  type input_variable
     !This derived type holds one entry :code:`NAME=value` parsed from a line of an input file, by
     !:f:func_inline:`scan_input_variable`, or from a command line argument, by :f:func_inline:`scan_cmd_variable`.
     !It is used internally by the parsing routines and it is not public.
     !
     character(len=64)                          :: name   ! name of the variable, converted to upper case
     character(len=256)                         :: value  ! text after the :code:`=` (vectors: comma separated list)
     character(len=20),dimension(:),allocatable :: args   ! fields of the list in :f:var:`value`, allocated for vectors only
  end type input_variable
  
  interface parse_cmd_variable
     !Generic interface to update a variable from the command line. It accepts integer, real(8), logical and character scalars
     !(:f:func_inline:`i_parse_cmd_variable`, :f:func_inline:`d_parse_cmd_variable`, :f:func_inline:`l_parse_cmd_variable`,
     !:f:func_inline:`ch_parse_cmd_variable`), and integer, real(8) and logical vectors
     !(:f:func_inline:`iv_parse_cmd_variable`, :f:func_inline:`dv_parse_cmd_variable`,
     !:f:func_inline:`lv_parse_cmd_variable`). The arguments are documented in :f:func_inline:`i_parse_cmd_variable` and
     !:f:func_inline:`iv_parse_cmd_variable`.
     !
     !First :f:var:`variable` is set to :f:var:`default`, if present. Then every command line argument of the form
     !:code:`NAME=value`, with :code:`NAME` matching :f:var:`name` (not case sensitive), is read with list-directed input into
     !:f:var:`variable`; if the argument is repeated the last one wins. A vector is given as a comma separated list without
     !spaces, :code:`NAME=v1,v2,...,vn`, with :code:`n=size(variable)` entries, and the program stops if the number of entries
     !differs from :code:`size(variable)`.
     !
     !A message :code:`Variable NAME updated to value` is written to standard error for each update. Unlike
     !:f:func_inline:`parse_input_variable`, the variable is not registered in the input list.
     !
     module procedure :: i_parse_cmd_variable
     module procedure :: d_parse_cmd_variable
     module procedure :: l_parse_cmd_variable
     module procedure :: iv_parse_cmd_variable
     module procedure :: dv_parse_cmd_variable
     module procedure :: lv_parse_cmd_variable
     module procedure :: ch_parse_cmd_variable
  end interface parse_cmd_variable


  interface parse_input_variable
     !Generic interface to read a variable from an input file and from the command line, and to register it in the input list.
     !It accepts integer, real(8), logical and character scalars (:f:func_inline:`i_parse_input`,
     !:f:func_inline:`d_parse_input`, :f:func_inline:`l_parse_input`, :f:func_inline:`ch_parse_input`), and integer, real(8)
     !and logical vectors (:f:func_inline:`iv_parse_input`, :f:func_inline:`dv_parse_input`, :f:func_inline:`lv_parse_input`).
     !The arguments are documented in :f:func_inline:`i_parse_input` and :f:func_inline:`iv_parse_input`.
     !
     !First :f:var:`variable` is set to :f:var:`default`, if present. Then, if :f:var:`file` exists, it is scanned line by
     !line for the first entry :code:`NAME=value` matching :f:var:`name` (not case sensitive), blanks and tabs being ignored
     !and anything following a comment character (:code:`!`, :code:`#` or :code:`%`) discarded. The value is read with
     !list-directed input; if the entry is not found the variable keeps its default value. A vector is given as a comma
     !separated list without spaces, :code:`NAME=v1,v2,...,vn`, with :code:`n=size(variable)` entries; a single trailing comma
     !is ignored.
     !
     !Finally, the command line is parsed with :f:func_inline:`parse_cmd_variable`, so that an argument :code:`NAME=value`
     !overrides the value read from the file, and the variable is registered, with its :f:var:`comment`, in the input list
     !used by :f:func_inline:`save_input` and :f:func_inline:`print_input`.
     !
     !A missing input file is not an error: the default is used and the module remembers it, so that
     !:f:func_inline:`save_input` dumps a default input file and stops. The program stops if the value can not be converted to
     !the type of :f:var:`variable`, or if the number of entries of a vector differs from :code:`size(variable)`.
     !
     module procedure :: i_parse_input
     module procedure :: d_parse_input
     module procedure :: l_parse_input
     module procedure :: iv_parse_input
     module procedure :: dv_parse_input
     module procedure :: lv_parse_input
     module procedure :: ch_parse_input
  end interface parse_input_variable


  interface save_input
     !Generic interface to dump the input list, i.e. all the variables registered so far by
     !:f:func_inline:`parse_input_variable`, to a file with the same format as the input files, :code:`NAME=value` followed by
     !the comment of the variable. The only instance is :f:func_inline:`save_input_file`, whose argument :f:var:`file` is the
     !base name of the output file.
     !
     !The output file is named :code:`used.`:f:var:`file` (see :f:func_inline:`print_input`), and it contains the final values
     !of the variables, i.e. after the command line has been applied. If any of the input files requested by
     !:f:func_inline:`parse_input_variable` could not be found, the file contains the *default* values and the program stops
     !after writing it, which provides a way to generate a template input file.
     !
     module procedure :: save_input_file
  end interface save_input

  interface print_input
     !Generic interface to print the input list, i.e. all the variables registered so far by
     !:f:func_inline:`parse_input_variable`. The only instance is :f:func_inline:`print_input_list`, whose arguments are the
     !optional base name :f:var:`file` of the output file and the optional :f:var:`list` to be printed (the default list if
     !absent).
     !
     !Each variable is printed, in the order it was registered, on a line :code:`NAME=value`, with vectors as comma separated
     !lists, followed by its comment preceded by :code:`!`, aligned at column 46 when possible: the format is the same as that
     !of the input files. If :f:var:`file` is present the output is written to the file :code:`used.`:f:var:`file`, replacing
     !it if it exists; otherwise it is written to the standard output. An empty list prints :code:`input list: empty`.
     !
     module procedure :: print_input_list
  end interface print_input

  public  :: parse_cmd_variable
  public  :: parse_input_variable
  public  :: check_input_variable
  public  :: save_input_file
  public  :: save_input
  public  :: print_input
  public  :: delete_input

  logical :: IOinput=.true.


contains

   function check_input_variable(name,file) result(check_input)
      !This subroutine checks if the variable "name" is present in the input list.
      !The output is a boolean: :code:`True` if the variable is present, :code:`False` if it is not present. 
      logical                  :: check_input
      character(len=*)         :: name      ! name of the variable in file/command line (case insens.)
      character(len=*)         :: file      ! input file name
      character(len=len(name)) :: name_
      type(input_variable)     :: var
      integer                  :: i,unit,pos
      integer                  :: status
      logical                  :: bool
      character(len=255)       :: buffer
      name_=name;call upper_case(name_)
      check_input=.false.
      inquire(file=file,exist=bool)
      if(bool)then
         unit=free_unit()
         open(unit,file=file)
         status=0
         var_search: do while(status>=0)
            read(unit,"(A255)",iostat=status)buffer
            pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
            var = scan_input_variable(trim(buffer))
            if(var%name==name_)then
               check_input=.true.
               exit var_search
            endif
         enddo var_search
         close(unit)
      endif
      return
   end function check_input_variable





  !---------------------------------------------------------------------
  !PURPOSE: PARSE INPUT  VARIABLES from file/CMD Line
  !---------------------------------------------------------------------
  !=====================SCALAR=====================================
  subroutine i_parse_input(variable,name,file,default,comment)
    integer                  :: variable  ! variable to be set; keeps its value if not found
    integer,optional         :: default   ! optional: default value, used if not in file
    character(len=*)         :: name      ! name of the variable in file/command line (case insens.)
    character(len=*),optional:: comment   ! optional: description, printed in the saved input file
    character(len=*)         :: file      ! input file name
    character(len=len(name)) :: name_
    type(input_variable)     :: var
    integer                  :: i,unit,pos
    integer                  :: status
    logical                  :: bool,idefault=.true.
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             read(var%value,*)variable
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine i_parse_input

  subroutine d_parse_input(variable,name,file,default,comment)
    real(8)                  :: variable
    real(8),optional         :: default
    character(len=*)         :: name
    character(len=*),optional:: comment
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  ::i,unit,pos
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             read(var%value,*)variable
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine d_parse_input

  subroutine l_parse_input(variable,name,file,default,comment)
    logical               :: variable
    logical,optional         :: default
    character(len=*)         :: name
    character(len=*),optional:: comment
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  ::i,unit,pos
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             read(var%value,*)variable
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine l_parse_input


  !=====================VECTOR=====================================
  subroutine iv_parse_input(variable,name,file,default,comment)
    integer,dimension(:)                       :: variable  ! vector to be set, size n
    integer,dimension(size(variable)),optional :: default   ! optional: default vector, size n, used if not in file
    character(len=*)                           :: name      ! name of the variable in file/command line (case insens.)
    character(len=*),optional:: comment                     ! optional: description, printed in the saved input file
    character(len=*)                           :: file      ! input file name
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             nargs=check_cmd_vector_size(ndim,var)
             allocate(var%args(nargs))
             iarg=0
             pos0=0
             do j=1,len(var%value)
                if(var%value(j:j)==",")then
                   iarg=iarg+1
                   var%args(iarg)=var%value(pos0+1:j-1)
                   pos0=j
                endif
             enddo
             var%args(nargs)=var%value(pos0+1:)
             do iarg=1,nargs
                read(var%args(iarg),*)variable(iarg)
             enddo
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine iv_parse_input

  subroutine dv_parse_input(variable,name,file,default,comment)
    real(8),dimension(:)                       :: variable
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*),optional:: comment
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                        :: var
    integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             nargs=check_cmd_vector_size(ndim,var)
             allocate(var%args(nargs))
             iarg=0
             pos0=0
             do j=1,len(var%value)
                if(var%value(j:j)==",")then
                   iarg=iarg+1
                   var%args(iarg)=var%value(pos0+1:j-1)
                   pos0=j
                endif
             enddo
             var%args(nargs)=var%value(pos0+1:)
             do iarg=1,nargs
                read(var%args(iarg),*)variable(iarg)
             enddo
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine dv_parse_input

  subroutine lv_parse_input(variable,name,file,default,comment)
    logical,dimension(:)                       :: variable
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*),optional:: comment
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             nargs=check_cmd_vector_size(ndim,var)
             allocate(var%args(nargs))
             iarg=0
             pos0=0
             do j=1,len(var%value)
                if(var%value(j:j)==",")then
                   iarg=iarg+1
                   var%args(iarg)=var%value(pos0+1:j-1)
                   pos0=j
                endif
             enddo
             var%args(nargs)=var%value(pos0+1:)
             do iarg=1,nargs
                read(var%args(iarg),*)variable(iarg)
             enddo
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine lv_parse_input


  !=====================STRING=====================================
  subroutine ch_parse_input(variable,name,file,default,comment)
    character(len=*)                  :: variable
    character(len=*),optional         :: default
    character(len=*)         :: name
    character(len=*),optional:: comment
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  ::i,unit,pos
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)IOinput=.false.
    if(bool)then
       unit=free_unit()
       open(unit,file=file)
       status=0
       var_search: do while(status>=0)
          read(unit,"(A255)",iostat=status)buffer
          pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
          var = scan_input_variable(trim(buffer))
          if(var%name==name_)then
             read(var%value,*)variable
             exit var_search
          endif
       enddo var_search
       close(unit)
    endif
    call parse_cmd_variable(variable,name_)
    if(present(comment))then
       call append_to_input_list(variable,name_,comment)
    else
       call append_to_input_list(variable,name_)
    endif
    return
  end subroutine ch_parse_input









  !---------------------------------------------------------------------
  !PURPOSE: PARSE INPUT  VARIABLES from CMD LINE
  ! - from file/CMD Line
  !---------------------------------------------------------------------
  !=====================SCALAR=====================================
  subroutine i_parse_cmd_variable(variable,name,default)
    integer                   :: variable  ! variable to be set; keeps its value if not given
    integer,optional          :: default   ! optional: default value
    character(len=*)          :: name      ! name of the variable on the command line (case insens.)
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    If(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine i_parse_cmd_variable

  subroutine d_parse_cmd_variable(variable,name,default)
    real(8)                   :: variable
    real(8),optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine d_parse_cmd_variable

  subroutine l_parse_cmd_variable(variable,name,default)
    logical                   :: variable
    logical,optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine l_parse_cmd_variable


  !=====================VECTOR=====================================
  subroutine iv_parse_cmd_variable(variable,name,default)
    integer,dimension(:)                       :: variable  ! vector to be set, size n
    integer,dimension(size(variable)),optional :: default   ! optional: default vector, size n
    character(len=*)                           :: name      ! name of the variable on the command line (case insens.)
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(0,"(A,100I6)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine iv_parse_cmd_variable

  subroutine dv_parse_cmd_variable(variable,name,default)
    real(8),dimension(:)                       :: variable
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(0,"(A,100F18.9)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine dv_parse_cmd_variable

  subroutine lv_parse_cmd_variable(variable,name,default)
    logical,dimension(:)                       :: variable
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(0,"(A,100L3)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine lv_parse_cmd_variable


  !=====================STRING=====================================
  subroutine ch_parse_cmd_variable(variable,name,default)
    character(len=*)          :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = scan_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine ch_parse_cmd_variable

















  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine save_input_file(file)
  !This subroutine dumps the input list, i.e. all the variables registered so far by
  !:f:func_inline:`parse_input_variable`, to a file with the same format as the input files, :code:`NAME=value` followed
  !by the comment of the variable. It is the public routine behind the generic name :f:func_inline:`save_input`.
  !
  !The output file is named :code:`used.`:f:var:`file` (see :f:func_inline:`print_input_list`), and it contains the
  !final values of the variables, i.e. after the command line has been applied. If any of the input files requested
  !by :f:func_inline:`parse_input_variable` could not be found, the file contains the *default* values and the program
  !stops after writing it, which provides a way to generate a template input file.
  !
    character(len=*)   :: file  ! base name of the output file, written as used.<file>
    if(IOinput)then
       call print_input_list(trim(file))
    else
       call print_input_list(trim(file))
       write(*,*)"input file can not be found. Dumped a *default* version in: used."//trim(file)
       stop
    endif
  end subroutine save_input_file

























  !---------------------------------------------------------------------
  !PURPOSE: AUXILIARY ROUTINES
  !---------------------------------------------------------------------
  function scan_comment(buffer) result(pos)
    character(len=255),intent(in)           :: buffer
    integer                                 :: pos
    character(len=1),dimension(3),parameter :: comments=["!","#","%"]
    integer :: i
    pos=0
    do i=1,3
       pos=scan(buffer,comments(i))
       if(pos/=0)return
    end do
  end function scan_comment

  function scan_cmd_variable(i)  result(var)
    integer            :: i,ncount,pos
    type(input_variable) :: var
    character(len=512) :: buffer
    ncount=command_argument_count()
    if(i>ncount)then
       write(*,*)"get_cmd_variable: i > cmdline_argument_count!"
       return
    endif
    call get_command_argument(i,buffer)
    pos      = scan(buffer,"=")
    var%name = buffer(1:pos-1);call upper_case(var%name)
    var%value= buffer(pos+1:)
  end function scan_cmd_variable

  function scan_input_variable(buffer)  result(var)
    character(len=*)   :: buffer
    integer            :: i,pos,iloc,len
    type(input_variable) :: var
    call s_blank_delete(buffer)
    pos      = scan(buffer,"=")
    var%name = buffer(1:pos-1)
    call upper_case(var%name)
    len=len_trim(buffer)
    if(buffer(len:len)==',')then
       var%value= buffer(pos+1:len-1)
    else
       var%value= buffer(pos+1:)
    endif
  end function scan_input_variable

  function check_cmd_vector_size(ndim,var) result(nargs)
    integer            :: ndim,ncount,j,nargs
    type(input_variable) :: var
    logical            :: iscalar
    
    iscalar=(scan(var%value,",")==0)
    if(ndim==1 .and. iscalar)then
       nargs=ndim
       return
    endif
    ncount=0
    do j=1,len(var%value)
       if(var%value(j:j)==",")ncount=ncount+1
    enddo
    nargs=ncount+1
    if(nargs/=ndim)then
       print*,"wrong dimensions parsing variable:   ",trim(var%name)
       print*,"expecting a comma separated list of: ",ndim," with no tail commas"
       print*,"gotten ",nargs, "argument fields (including after any tail commas)."
       STOP
    endif
  end function check_cmd_vector_size

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

end module SF_PARSE_INPUT


