MODULE SF_IOTOOLS
  !SciFortran module for reading and writing data
  USE IOFILE
  USE IOPLOT
  USE IOREAD
  private

  public :: splot
  public :: splot3d
  public :: save_array
  !
  public :: sread
  public :: read_array
  !
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

END MODULE SF_IOTOOLS
