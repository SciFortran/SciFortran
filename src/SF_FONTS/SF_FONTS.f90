MODULE SF_FONTS
!SciFortran module for font customization
  implicit none
  private

  public :: bold
  public :: underline
  public :: highlight
  public :: erased
  !
  public :: font_red
  public :: font_green
  public :: font_yellow
  public :: font_blue
  !
  public :: bold_red
  public :: bold_green
  public :: bold_yellow
  public :: bold_blue
  !
  public :: bg_red
  public :: bg_green
  public :: bg_yellow
  public :: bg_blue


contains

  function bold(text) result(textout)
  !This function returns the string :f:var:`text` in bold face. The text is enclosed 
  !between the ANSI escape sequences :code:`ESC[1m` and :code:`ESC[0m`, so the 
  !returned string has length :code:`8+len(text)`.
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[1m"//text//achar(27)//"[0m"
  end function bold

  function underline(text) result(textout)
  !This function returns the string :f:var:`text` underlined, using the ANSI escape 
  !sequence :code:`ESC[4m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`8+len(text)`.
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[4m"//text//achar(27)//"[0m"
  end function underline

  function highlight(text) result(textout)
  !This function returns the string :f:var:`text` highlighted, using the ANSI escape 
  !sequence :code:`ESC[7m` (reverse video: foreground and background colors are 
  !swapped) and resetting with :code:`ESC[0m`. The returned string has length 
  !:code:`8+len(text)`.
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[7m"//text//achar(27)//"[0m"
  end function highlight

  function erased(text) result(textout)
  !This function returns the string :f:var:`text` struck through, using the ANSI escape 
  !sequence :code:`ESC[9m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`8+len(text)`.
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[9m"//text//achar(27)//"[0m"
  end function erased

  function font_red(text) result(textout)
  !This function returns the string :f:var:`text` in red font, using the ANSI escape 
  !sequence :code:`ESC[91m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[91m"//text//achar(27)//"[0m"
  end function font_red

  function font_green(text) result(textout)
  !This function returns the string :f:var:`text` in green font, using the ANSI escape 
  !sequence :code:`ESC[92m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[92m"//text//achar(27)//"[0m"
  end function font_green

  function font_yellow(text) result(textout)
  !This function returns the string :f:var:`text` in yellow font, using the ANSI escape 
  !sequence :code:`ESC[93m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[93m"//text//achar(27)//"[0m"
  end function font_yellow

  function font_blue(text) result(textout)
  !This function returns the string :f:var:`text` in blue font, using the ANSI escape 
  !sequence :code:`ESC[94m` and resetting with :code:`ESC[0m`. The returned string has 
  !length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[94m"//text//achar(27)//"[0m"
  end function font_blue

  function bold_red(text) result(textout)
  !This function returns the string :f:var:`text` in bold red font, using the ANSI escape 
  !sequence :code:`ESC[1;91m` and resetting with :code:`ESC[0m`. It combines 
  !:f:func_inline:`bold` and :f:func_inline:`font_red`. The returned string has length 
  !:code:`11+len(text)`.
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;91m"//text//achar(27)//"[0m"
  end function bold_red

  function bold_green(text) result(textout)
  !This function returns the string :f:var:`text` in bold green font, using the ANSI escape 
  !sequence :code:`ESC[1;92m` and resetting with :code:`ESC[0m`. It combines 
  !:f:func_inline:`bold` and :f:func_inline:`font_green`. The returned string has length 
  !:code:`11+len(text)`.
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;92m"//text//achar(27)//"[0m"
  end function bold_green

  function bold_yellow(text) result(textout)
  !This function returns the string :f:var:`text` in bold yellow font, using the ANSI escape 
  !sequence :code:`ESC[1;93m` and resetting with :code:`ESC[0m`. It combines 
  !:f:func_inline:`bold` and :f:func_inline:`font_yellow`. The returned string has length 
  !:code:`11+len(text)`.
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;93m"//text//achar(27)//"[0m"
  end function bold_yellow

  function bold_blue(text) result(textout)
  !This function returns the string :f:var:`text` in bold blue font, using the ANSI escape 
  !sequence :code:`ESC[1;94m` and resetting with :code:`ESC[0m`. It combines 
  !:f:func_inline:`bold` and :f:func_inline:`font_blue`. The returned string has length 
  !:code:`11+len(text)`.
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;94m"//text//achar(27)//"[0m"
  end function bold_blue

  function bg_red(text) result(textout)
  !This function returns the string :f:var:`text` on a red background, using the ANSI 
  !escape sequence :code:`ESC[41m` and resetting with :code:`ESC[0m`. The returned 
  !string has length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[41m"//text//achar(27)//"[0m"
  end function bg_red

  function bg_green(text) result(textout)
  !This function returns the string :f:var:`text` on a green background, using the ANSI 
  !escape sequence :code:`ESC[42m` and resetting with :code:`ESC[0m`. The returned 
  !string has length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[42m"//text//achar(27)//"[0m"
  end function bg_green

  function bg_yellow(text) result(textout)
  !This function returns the string :f:var:`text` on a yellow background, using the ANSI 
  !escape sequence :code:`ESC[43m` and resetting with :code:`ESC[0m`. The returned 
  !string has length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[43m"//text//achar(27)//"[0m"
  end function bg_yellow

  function bg_blue(text) result(textout)
  !This function returns the string :f:var:`text` on a blue background, using the ANSI 
  !escape sequence :code:`ESC[44m` and resetting with :code:`ESC[0m`. The returned 
  !string has length :code:`9+len(text)`.
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[44m"//text//achar(27)//"[0m"
  end function bg_blue

END MODULE SF_FONTS
