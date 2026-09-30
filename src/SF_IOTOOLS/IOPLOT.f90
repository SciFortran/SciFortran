module IOPLOT
  USE IOFILE
  !Contains procedures for generating plotting scripts and data files
  implicit none
  private

  interface splot
  !This subroutine writes a real or complex array :f:var:`Y1` of rank 1 to 7 to the file :f:var:`pname`, together with the
  !abscissa :f:var:`X`, in a format that can be plotted with gnuplot. The last dimension of :f:var:`Y1` runs along :f:var:`X`.
  !Each line of the file contains
  !
  !* :code:`X(k)  Y1(...,k)` for a real array
  !* :code:`X(k)  Im Y1(...,k)  Re Y1(...,k)` for a complex array, so the imaginary part is in column 2 and the real part in
  !  column 3
  !
  !For rank 2 to 7 a blank line follows each set of :code:`size(X)` lines, so that each one-dimensional slice is a
  !separate curve. If :f:var:`append` is true the data are added at the end of the file, after a blank line if the file
  !already exists, otherwise an existing file is overwritten. The data can be read back with :f:func_inline:`sread`.
  !
     module procedure :: splotA1_RR
     module procedure :: splotA1_RC
     module procedure :: splotA2_RR
     module procedure :: splotA2_RC
     module procedure :: splotA3_RR
     module procedure :: splotA3_RC
     module procedure :: splotA4_RR
     module procedure :: splotA4_RC
     module procedure :: splotA5_RR
     module procedure :: splotA5_RC
     module procedure :: splotA6_RR
     module procedure :: splotA6_RC
     module procedure :: splotA7_RR
     module procedure :: splotA7_RC
  end interface splot


  interface splot3d
  !This subroutine writes a function of two variables, given on the grid :f:var:`X1` :math:`\times` :f:var:`X2`, to files that
  !can be plotted with gnuplot as a color map and as a surface. The specific procedures cover a real or complex array
  !:f:var:`Y` of shape :code:`(size(X1),size(X2))`, and a real or complex sequence of :code:`Nt` frames :f:var:`Y` of shape
  !:code:`(size(X1),size(X2),Nt)` for an animated map.
  !
  !For a real :f:var:`Y` the data file :f:var:`pname` contains the lines :code:`X1(i)  X2(j)  Y(i,j)`, with a blank line after
  !each value of :code:`i`. The gnuplot scripts :code:`pname_map.gp` (color map) and, unless :f:var:`nosurface` is true,
  !:code:`pname_surface.gp` (surface view) are written next to it. If :f:var:`wlines` is present, whatever its value, the file
  !:code:`pname_withlines` is also written with every :f:var:`nlines`-th value of :code:`i`, to be plotted as lines over the
  !surface.
  !
  !For a complex :f:var:`Y` the real and imaginary parts are written in the files :code:`re_name` and :code:`im_name`, in the
  !directory
  !of :f:var:`pname`, with the scripts :code:`pname_re_map.gp`, :code:`pname_im_map.gp`, :code:`pname_re_surface.gp` and
  !:code:`pname_im_surface.gp`.
  !
  !The animate procedures write all the frames in the data file, one block per frame, and only the map script(s), which loop over
  !the frames with a color range common to all of them. They stop if the first two dimensions of :f:var:`Y` do not match
  !:f:var:`X1` and :f:var:`X2`. The plot ranges default to the extrema of :f:var:`X1` and :f:var:`X2`, see :f:var:`xmin`.
  !
     module procedure :: d_splot3d
     module procedure :: c_splot3d
     module procedure :: d_splot3d_animate
     module procedure :: c_splot3d_animate
  end interface splot3d

  interface save_array
  !This subroutine writes a real or complex scalar, or array of rank 1 to 7, :f:var:`Y1` to the file :f:var:`pname` using
  !list-directed
  !output, one element per line, complex elements as :code:`(re,im)`. After writing, the file is compressed with
  !:f:func_inline:`file_bzip` if it is larger than the store size, see :f:func_inline:`set_store_size`.
  !
  !For rank 2 to 7 the order of the elements is set by :f:var:`order`. With :code:`order="R"` the last index varies fastest, with
  !:code:`order="C"` the first index varies fastest. If :f:var:`wspace` is true a blank line is written after each run of the
  !fastest
  !index. The program stops if :f:var:`order` does not start with :code:`R` or :code:`C`. The data can be read back with
  !:f:func_inline:`read_array`.
  !
     module procedure :: data_saveA0_R
     module procedure :: data_saveA0_C
     module procedure :: data_saveA1_R
     module procedure :: data_saveA1_C
     module procedure :: data_saveA2_R
     module procedure :: data_saveA2_C
     module procedure :: data_saveA3_R
     module procedure :: data_saveA3_C
     module procedure :: data_saveA4_R
     module procedure :: data_saveA4_C
     module procedure :: data_saveA5_R
     module procedure :: data_saveA5_C
     module procedure :: data_saveA6_R
     module procedure :: data_saveA6_C
     module procedure :: data_saveA7_R
     module procedure :: data_saveA7_C
  end interface save_array




  public :: splot
  public :: splot3d
  public :: save_array

  integer            :: unit
  character(len=128) :: fmt


contains


  ! SPLOT ararys (1--7)
  include "ioplot_splot.f90"


  ! SPLOT 3D:
  include "ioplot_splot3d.f90"


  ! SAVE_ARRAY arrays (0--7)
  include "ioplot_save_array.f90"


end module IOPLOT
