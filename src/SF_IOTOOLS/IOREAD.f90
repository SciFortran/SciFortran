module IOREAD
  !Contains procedures for reading files
  USE IOFILE
  implicit none
  private

  integer            :: unit
  character(len=128) :: fmt
  logical            :: control

  interface sread
  !This subroutine reads a real or complex array of rank 1 to 7 from the file :f:var:`pname`, together with the abscissa
  !:code:`X`,
  !and is the counterpart of :f:func_inline:`splot`. The last dimension of :code:`Y1` runs along :code:`X`, and each line of the
  !file
  !contains :code:`X(k)  Y1(...,k)` for a real array and :code:`X(k)  Im Y1(...,k)  Re Y1(...,k)` for a complex array. The blank
  !lines between the curves are skipped. If the file, or its :code:`.gz` version, is not found a message is printed, the
  !program waits 5 seconds, and then stops with an end-of-file error when reading.
  !
     module procedure :: sreadA1_RR
     module procedure :: sreadA1_RC
     module procedure :: sreadA2_RR
     module procedure :: sreadA2_RC
     module procedure :: sreadA3_RR
     module procedure :: sreadA3_RC
     module procedure :: sreadA4_RR
     module procedure :: sreadA4_RC
     module procedure :: sreadA5_RR
     module procedure :: sreadA5_RC
     module procedure :: sreadA6_RR
     module procedure :: sreadA6_RC
     module procedure :: sreadA7_RR
     module procedure :: sreadA7_RC
  end interface sread


  interface read_array
  !This subroutine reads a real or complex scalar, or array of rank 1 to 7, :code:`Y1` from the file :f:var:`pname`, and is the
  !counterpart of :f:func_inline:`save_array`. If :f:var:`pname` does not exist the file :code:`pname.bz2` is first uncompressed
  !with :f:func_inline:`file_bunzip`, and after reading the file is compressed again with :f:func_inline:`file_bzip` if it is
  !larger
  !than the store size, see :f:func_inline:`set_store_size`. If the file, or its :code:`.gz` version, is not found a message is
  !printed,
  !the program waits 5 seconds, and then stops with an end-of-file error when reading.
  !
  !For rank 2 to 7 the order of the elements is set by :f:var:`order`. With :code:`order="R"` the last index varies fastest, with
  !:code:`order="C"` the first index varies fastest. The program stops if :f:var:`order` does not start with :code:`R` or
  !:code:`C`.
  !The blank lines separating the runs are skipped, so :f:var:`wspace` has no effect when reading.
  !
     module procedure :: data_readA0_R
     module procedure :: data_readA0_C
     module procedure :: data_readA1_R
     module procedure :: data_readA1_C
     module procedure :: data_readA2_R
     module procedure :: data_readA2_C
     module procedure :: data_readA3_R
     module procedure :: data_readA3_C
     module procedure :: data_readA4_R
     module procedure :: data_readA4_C
     module procedure :: data_readA5_R
     module procedure :: data_readA5_C
     module procedure :: data_readA6_R
     module procedure :: data_readA6_C
     module procedure :: data_readA7_R
     module procedure :: data_readA7_C
  end interface read_array


  public :: sread
  public :: read_array


contains


  ! SPLOT arrays (1--7)
  include "ioread_sread.f90"

  ! READ_ARRAY arrays (0--7)
  include "ioread_read_array.f90"

  ! READ safety utility
  include "ioread_control.f90"


end module IOREAD
