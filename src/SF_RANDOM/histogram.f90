function histogram_allocate(n) result(h)
  !This function allocates a histogram of :f:var:`n` bins, with all the edges and the bin contents set to zero. The edges can
  !be defined with :f:func_inline:`histogram_set_range_uniform`. The program stops if :f:var:`n` is not positive.
  !
  integer,intent(in) :: n  ! number of bins
  type(histogram)    :: h  ! histogram, allocated
  if(n<=0)then
     print*,"histogram length must be positive integer. n=",n
     stop
  endif
  allocate(h%range(0:n));h%range=0.d0
  allocate(h%bin(0:n))    ;h%bin=0.d0
  h%n=n
end function histogram_allocate


subroutine histogram_deallocate(h)
  !This subroutine deallocates the histogram :f:var:`h` and sets its number of bins to zero.
  !
  type(histogram)    :: h  ! histogram to be deallocated
  deallocate(h%range)
  deallocate(h%bin)
  h%n=0
end subroutine  histogram_deallocate


subroutine histogram_reset(h)
  !This subroutine sets to zero both the edges and the contents of the bins of the histogram :f:var:`h`, which is kept
  !allocated. The edges must therefore be defined again with :f:func_inline:`histogram_set_range_uniform`.
  !
  type(histogram)    :: h  ! histogram to be reset
  h%range=0.d0
  h%bin=0.d0
end subroutine histogram_reset


subroutine histogram_set_range_uniform(h,xmin,xmax)
  !This subroutine defines the edges of the bins of the histogram :f:var:`h` as :code:`n+1` equally spaced points between
  !:f:var:`xmin` and :f:var:`xmax`, and it resets the contents of the bins to zero. The program stops if :code:`xmin>=xmax`.
  !
  type(histogram),intent(inout) :: h     ! histogram, allocated
  real(8),intent(in)            :: xmin  ! lower edge of the first bin
  real(8),intent(in)            :: xmax  ! upper edge of the last bin
  integer                       :: i,n
  real(8)                       :: f1,f2
  if(xmin>=xmax)then
     print*,"histogram%range: xmin must be less than xmax:",xmin,xmax
     stop
  endif
  n=h%n
  do i=0,n
     f1= real(n-i,8)/real(n,8)
     f2= real(i,8)/real(n,8)
     h%range(i) = f1*xmin + f2*xmax
  enddo
  h%bin=0.d0
end subroutine histogram_set_range_uniform


subroutine histogram_accumulate(h,x,w)
  !This subroutine adds the weight :f:var:`w` to the bin of the histogram :f:var:`h` containing the value :f:var:`x`, found
  !with a bisection search. If :f:var:`x` is outside the range of the histogram the message :code:`X out of range!` is printed
  !and, as the routine is written, the weight is added to the first bin.
  !
  type(histogram),intent(inout) :: h  ! histogram, with the edges defined
  real(8),intent(in)            :: x  ! value to be accumulated
  real(8),intent(in)            :: w  ! weight of the value
  integer                       :: i,index
  index=0
  call find_index(h%n,h%range,x,index)
  if(index>=h%n)then
     print*,"index lies outside valid range of 0 .. n - 1"
     stop
  endif
  h%bin(index)=h%bin(index)+w
end subroutine histogram_accumulate


subroutine find_index(n,range,x,index)
  integer,intent(in)                :: n
  real(8),dimension(0:n),intent(in) :: range
  real(8),intent(in)                :: x
  integer,intent(out)               :: index
  integer                           :: i,upper,lower,mid
  if((x<range(0)) .OR. (x>range(n)))then
     print*,"X out of range!"
     return
  endif
  upper=n
  lower=0
  do while((upper-lower>1))
     mid = (upper+lower)/2    !int(dble(upper + lower)/2.d0)
     if( x >= range(mid))then
        lower=mid
     else
        upper=mid
     endif
  enddo
  index=lower
  if(x<range(lower) .OR. x>range(lower+1))then
     print*,"error: x not found within range!"
     stop
  endif
end subroutine find_index


subroutine histogram_get_range(h,index,lower,upper)
  !This subroutine returns the edges of the bin number :f:var:`index` of the histogram :f:var:`h`: the bins are numbered from
  !0 to :code:`n-1`. The program stops if :f:var:`index` is not smaller than :code:`n`.
  !
  type(histogram),intent(in) :: h      ! histogram
  integer,intent(in)         :: index  ! number of the bin, from 0 to n-1
  real(8),intent(out)        :: lower  ! lower edge of the bin
  real(8),intent(out)        :: upper  ! upper edge of the bin
  if(index>=h%n)then
     print*,"error: *i lies outside valid range=0...n-1"
     stop
  endif
  lower=h%range(index)
  upper=h%range(index+1)
end subroutine histogram_get_range


function histogram_get_value(h,index) result(value)
  !This function returns the content of the bin number :f:var:`index` of the histogram :f:var:`h`: the bins are numbered from
  !0 to :code:`n-1`. The program stops if :f:var:`index` is not smaller than :code:`n`.
  !
  type(histogram),intent(in) :: h      ! histogram
  integer,intent(in)         :: index  ! number of the bin, from 0 to n-1
  real(8)                    :: value  ! content of the bin
  if(index>=h%n)then
     print*,"error: *index lies outside valid range=0...n-1"
     stop
  endif
  value=h%bin(index)
end function histogram_get_value


subroutine histogram_print(h,unit)
  !This subroutine writes the histogram :f:var:`h` on the unit :f:var:`unit`, as a step profile ready to be plotted: pairs
  !:code:`edge content` with the two edges of each bin, preceded and followed by a point of zero content, in the format
  !:code:`2F12.7`.
  !
  type(histogram),intent(in) :: h     ! histogram
  integer,intent(in)         :: unit  ! output unit, already open
  integer                    :: i,n
  real(8)                    :: lower,upper,bin_value
  n=h%n
  call histogram_get_range(h,0,lower,upper)
  write(unit,"(2F12.7)")lower,0.d0
  do i=0,n-1
     call histogram_get_range(h,i,lower,upper)
     bin_value = histogram_get_value(h,i)
     write(unit,"(2F12.7)")lower,bin_value
     write(unit,"(2F12.7)")upper,bin_value
  enddo
  write(unit,"(2F12.7)")upper,0.d0
end subroutine histogram_print

