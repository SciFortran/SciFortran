module SF_MISC
!SciFortran module for sundry functions
  implicit none
  private
  real(8),parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510d0

  !UNIINV = Merge-sort inverse ranking of an array, with removal of duplicate entries. 
  interface uniinv
     module procedure d_uniinv, r_uniinv, i_uniinv
  end interface uniinv

  !UNISTA = (stable unique) removes duplicates from an array, leaving unique entries
  ! in the order of their first appearance in the initial set.
  interface unista
     module procedure d_unista, r_unista, i_unista
  end interface unista


  interface uniq_array
     module procedure :: i_uniq
     module procedure :: d_uniq
  end interface uniq_array

  interface uniq
     module procedure :: f_i_uniq
     module procedure :: f_d_uniq
  end interface uniq


  interface nearless
     module procedure D_nearless, R_nearless, I_nearless
  end interface nearless

  interface assert_shape
     module procedure i_assert_shape_N1
     module procedure i_assert_shape_N2
     module procedure i_assert_shape_N3
     module procedure i_assert_shape_N4
     module procedure i_assert_shape_N5
     module procedure i_assert_shape_N6
     module procedure i_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure i_assert_shape_N8
#endif
     !
     module procedure d_assert_shape_N1
     module procedure d_assert_shape_N2
     module procedure d_assert_shape_N3
     module procedure d_assert_shape_N4
     module procedure d_assert_shape_N5
     module procedure d_assert_shape_N6
     module procedure d_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure d_assert_shape_N8
#endif
     !
     module procedure z_assert_shape_N1
     module procedure z_assert_shape_N2
     module procedure z_assert_shape_N3
     module procedure z_assert_shape_N4
     module procedure z_assert_shape_N5
     module procedure z_assert_shape_N6
     module procedure z_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
     module procedure z_assert_shape_N8
#endif
  end interface assert_shape


  interface reorder_array
     module procedure :: I_reshuffle
     module procedure :: D_reshuffle
     module procedure :: Z_reshuffle
     module procedure :: L_reshuffle
  end interface reorder_array

  interface reorder
     module procedure :: f_I_reshuffle
     module procedure :: f_D_reshuffle
     module procedure :: f_Z_reshuffle
     module procedure :: f_L_reshuffle
  end interface reorder

  interface sort_insertion
     module procedure :: sort_insertion_i
     module procedure :: sort_insertion_d
  end interface sort_insertion

  interface sort_quicksort
     module procedure :: sort_quicksort_i
     module procedure :: sort_quicksort_d
  end interface sort_quicksort

  interface sort
     module procedure :: sort_quicksort_i
     module procedure :: sort_quicksort_d
  end interface sort

  interface sort_array
     module procedure :: sort_quicksort_i
     module procedure :: sort_quicksort_d
  end interface sort_array


  interface sort_qsort
     module procedure :: sort_array_i
     module procedure :: sort_array_d
  end interface sort_qsort

  public :: sort
  public :: sort_array
  public :: sort_quicksort
  public :: sort_insertion
  public :: sort_qsort
  !
  public :: reorder_array
  public :: reorder
  !
  public :: assert_shape
  !
  public :: uniq_array
  public :: uniq
  !
  public :: uniinv
  public :: unista

contains



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: assert shape of a fiven array up to 7 dimensions (gnu gfortran max rank)
  !+-----------------------------------------------------------------------------+!
  subroutine i_assert_shape_N1(A,Ndim,routine,matname)
    integer,dimension(:),intent(in)          :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N1
  subroutine i_assert_shape_N2(A,Ndim,routine,matname)
    integer,dimension(:,:),intent(in)          :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N2
  subroutine i_assert_shape_N3(A,Ndim,routine,matname)
    integer,dimension(:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N3
  subroutine i_assert_shape_N4(A,Ndim,routine,matname)
    integer,dimension(:,:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N4
  subroutine i_assert_shape_N5(A,Ndim,routine,matname)
    integer,dimension(:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N5
  subroutine i_assert_shape_N6(A,Ndim,routine,matname)
    integer,dimension(:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N6
  subroutine i_assert_shape_N7(A,Ndim,routine,matname)
    integer,dimension(:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine i_assert_shape_N8(A,Ndim,routine,matname)
    integer,dimension(:,:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine i_assert_shape_N8
#endif
  !
  !
  !
  subroutine d_assert_shape_N1(A,Ndim,routine,matname)
    real(8),dimension(:),intent(in)            :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N1
  subroutine d_assert_shape_N2(A,Ndim,routine,matname)
    real(8),dimension(:,:),intent(in)          :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N2
  subroutine d_assert_shape_N3(A,Ndim,routine,matname)
    real(8),dimension(:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N3
  subroutine d_assert_shape_N4(A,Ndim,routine,matname)
    real(8),dimension(:,:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N4
  subroutine d_assert_shape_N5(A,Ndim,routine,matname)
    real(8),dimension(:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N5
  subroutine d_assert_shape_N6(A,Ndim,routine,matname)
    real(8),dimension(:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N6
  subroutine d_assert_shape_N7(A,Ndim,routine,matname)
    real(8),dimension(:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine d_assert_shape_N8(A,Ndim,routine,matname)
    real(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine d_assert_shape_N8
#endif
  !
  !
  !
  subroutine z_assert_shape_N1(A,Ndim,routine,matname)
    complex(8),dimension(:),intent(in)         :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N1
  subroutine z_assert_shape_N2(A,Ndim,routine,matname)
    complex(8),dimension(:,:),intent(in)          :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N2
  subroutine z_assert_shape_N3(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N3
  subroutine z_assert_shape_N4(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:,:),intent(in)        :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N4
  subroutine z_assert_shape_N5(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N5
  subroutine z_assert_shape_N6(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N6
  subroutine z_assert_shape_N7(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N7
#if defined __GFORTRAN__ &&  __GNUC__ > 8
  subroutine z_assert_shape_N8(A,Ndim,routine,matname)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: A
    integer,dimension(:),intent(in)            :: Ndim
    character(len=*),optional                  :: routine, matname
    if(any(shape(A) /= Ndim)) then
       if(present(routine).AND.present(matname))&
            write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
       stop "assert_shape error: wrong matrix shape"
    end if
  end subroutine z_assert_shape_N8
#endif




  ! SORTING 1D:
  !###################################################################
  !> sort using the insertion method (scale as N**2, small arrays)
  subroutine sort_insertion_i(a,indx_a)
    integer, dimension(:),intent(inout)      :: a
    integer,dimension(size(a)),intent(inout) :: indx_a
    integer                                  :: na
    real(8)                                  :: temp
    integer                                  :: i, j
    integer                                  :: idx_tmp
    !
    nA=size(a)
    do i=1,nA
       indx_a(i)=i
    enddo
    !
    do i = 2, nA
       j = i - 1
       temp = A(i)
       idx_tmp = indx_a(i)
       do
          if (j == 0) exit
          if (a(j) <= temp) exit
          A(j+1) = A(j)
          indx_a(j+1) = indx_a(j)
          j = j - 1
       end do
       a(j+1) = temp
       indx_a(j+1) = idx_tmp
    end do
  end subroutine sort_insertion_i
  !
  subroutine sort_insertion_d(a,indx_a)
    real(8), dimension(:),intent(inout)      :: a
    integer,dimension(size(a)),intent(inout) :: indx_a
    integer                                  :: na
    real(8)                                  :: temp
    integer                                  :: i, j
    integer                                  :: idx_tmp
    !
    nA=size(a)
    do i=1,nA
       indx_a(i)=i
    enddo
    !
    do i = 2, nA
       j = i - 1
       temp = A(i)
       idx_tmp = indx_a(i)
       do
          if (j == 0) exit
          if (a(j) <= temp) exit
          A(j+1) = A(j)
          indx_a(j+1) = indx_a(j)
          j = j - 1
       end do
       a(j+1) = temp
       indx_a(j+1) = idx_tmp
    end do
  end subroutine sort_insertion_d



  !> subroutine to sort using the quicksort algorithm
  subroutine sort_quicksort_i(a,indx)
    integer,dimension(:),intent(inout)                :: a
    integer,dimension(size(a)),intent(inout),optional :: indx
    integer,dimension(size(a))                        :: indx_a
    integer                                           :: i
    do i=1,size(a)
       indx_a(i)=i
    enddo
    !
    call quicksort_i(a,indx_a,size(a))
    if(present(indx))indx=indx_a
    !
  contains
    !
    recursive subroutine quicksort_i(list,order,nA)
      integer              :: nA
      integer              :: list(nA)
      integer              :: order(nA)
      integer              :: left, right, mid
      integer              :: pivot, temp
      integer              :: marker, idx_temp,i
      !
      if(na==1)return
      !
      mid = (1+nA)/2
      if (list(mid) >= list(1)) then
         if (list(mid) <= list(nA)) then
            pivot = list(mid)
         else if (list(nA) > list(1)) then
            pivot = list(nA)
         else
            pivot = list(1)
         end if
      else if (list(1) <= list(nA)) then
         pivot = list(1)
      else if (list(nA) > list(mid)) then
         pivot = list(nA)
      else
         pivot = list(mid)
      end if
      left  = 0
      right = nA + 1
      do while (left < right)
         right = right - 1
         do while (list(right) > pivot)
            right = right - 1
         end do
         left = left + 1
         do while (list(left) < pivot)
            left = left + 1
         end do
         if (left < right) then
            temp = list(left)
            list(left) = list(right)
            list(right) = temp
            !
            idx_temp     = order(left)
            order(left)  = order(right)
            order(right) = idx_temp
         end if
      end do
      !
      if (left == right) then
         marker = left + 1
      else
         marker = left
      end if
      !
      call quicksort_i(list(:marker-1),order(:marker-1),marker-1)
      call quicksort_i(list(marker:),order(marker:),nA-marker+1)
      !
    end subroutine quicksort_i
  end subroutine sort_quicksort_i

  !
  subroutine sort_quicksort_d(a,indx)
    real(8), dimension(:),intent(inout)               :: a
    integer,dimension(size(a)),intent(inout),optional :: indx
    integer,dimension(size(a))                        :: indx_a
    integer                                           :: i
    do i=1,size(a)
       indx_a(i)=i
    enddo
    !
    call quicksort_d(a,indx_a,size(a))
    if(present(indx))indx=indx_a
    !
  contains
    !
    recursive subroutine quicksort_d(list,order,nA)
      integer              :: nA
      real(8)              :: list(nA)
      integer              :: order(nA)
      integer              :: left, right, mid
      real(8)              :: pivot, temp
      integer              :: marker, idx_temp,i
      !
      if(na==1)return
      !
      mid = (1+nA)/2
      if (list(mid) >= list(1)) then
         if (list(mid) <= list(nA)) then
            pivot = list(mid)
         else if (list(nA) > list(1)) then
            pivot = list(nA)
         else
            pivot = list(1)
         end if
      else if (list(1) <= list(nA)) then
         pivot = list(1)
      else if (list(nA) > list(mid)) then
         pivot = list(nA)
      else
         pivot = list(mid)
      end if
      left  = 0
      right = nA + 1
      do while (left < right)
         right = right - 1
         do while (list(right) > pivot)
            right = right - 1
         end do
         left = left + 1
         do while (list(left) < pivot)
            left = left + 1
         end do
         if (left < right) then
            temp = list(left)
            list(left) = list(right)
            list(right) = temp
            !
            idx_temp     = order(left)
            order(left)  = order(right)
            order(right) = idx_temp
         end if
      end do
      !
      if (left == right) then
         marker = left + 1
      else
         marker = left
      end if
      !
      call quicksort_d(list(:marker-1),order(:marker-1),marker-1)
      call quicksort_d(list(marker:),order(marker:),nA-marker+1)
      !
    end subroutine quicksort_d
  end subroutine sort_quicksort_d





  !> subroutine to sort using the quicksort algorithm in a somehow slower implementation:
  subroutine sort_array_i(a,indx_a)
    integer, dimension(:),intent(inout)      :: a
    integer,dimension(size(a)),intent(inout) :: indx_a
    integer                                  :: na
    integer,dimension(size(a))               :: a_tmp
    integer                                  :: i
    nA=size(a)
    do i=1,nA
       indx_a(i)=i
    enddo
    call qsort_sort(a, indx_a, 1, nA)
    !
    do i=1,nA
       a_tmp(i) = a(indx_a(i))
    enddo
    a = a_tmp
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      implicit none
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left >= right ) return
      call swap_order( order, left, ran_uniform(left,right))
      last = left
      do i = left+1, right
         if ( array(order(i)) < array(order(left))  ) then
            last = last + 1
            call swap_order( order, last, i )
         endif
      enddo
      call swap_order( order, left, last )      
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !
    subroutine swap_order( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine swap_order
    !
    function ran_uniform(l,h) result(igrnd)
      integer,intent(in) :: l,h
      real(8)            :: u,r
      integer           :: igrnd
      call random_number(u)
      r=(h-l+1)*u+l
      igrnd=int(r)
    end function ran_uniform
  end subroutine sort_array_i


  subroutine sort_array_d(a,indx_a)
    real(8), dimension(:),intent(inout)      :: a
    integer,dimension(size(a)),intent(inout) :: indx_a
    integer                                  :: na
    real(8),dimension(size(a))               :: a_tmp
    integer                                  :: i
    nA=size(a)
    do i=1,nA
       indx_a(i)=i
    enddo
    call qsort_sort(a, indx_a, 1, nA)
    !
    do i=1,nA
       a_tmp(i) = a(indx_a(i))
    enddo
    a = a_tmp
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      implicit none
      real(8), dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left >= right ) return
      call swap_order( order, left, ran_uniform(left,right))
      last = left
      do i = left+1, right
         if ( array(order(i)) < array(order(left))  ) then
            last = last + 1
            call swap_order( order, last, i )
         endif
      enddo
      call swap_order( order, left, last )      
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !
    subroutine swap_order( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine swap_order
    !
    function ran_uniform(l,h) result(igrnd)
      integer,intent(in) :: l,h
      real(8)            :: u,r
      integer           :: igrnd
      call random_number(u)
      r=(h-l+1)*u+l
      igrnd=int(r)
    end function ran_uniform
  end subroutine sort_array_d















  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  subroutine i_uniq(AIN,AOUT,MASK)
    integer,dimension(:),intent(INOUT)                    :: AIN
    integer,dimension(:),allocatable,intent(OUT)          :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call i_unista(AIN,ndim,MASK)
    else
       call i_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end subroutine i_uniq
  subroutine d_uniq(AIN,AOUT,MASK)
    real(8),dimension(:),intent(INOUT)                    :: AIN
    real(8),dimension(:),allocatable,intent(OUT)          :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call d_unista(AIN,ndim,MASK)
    else
       call d_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end subroutine d_uniq

  function f_i_uniq(AIN,MASK) result(AOUT)
    integer,dimension(:),intent(INOUT)                    :: AIN
    integer,dimension(:),allocatable                      :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call i_unista(AIN,ndim,MASK)
    else
       call i_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end function f_i_uniq
  function f_d_uniq(AIN,MASK) result(AOUT)
    real(8),dimension(:),intent(INOUT)                    :: AIN
    real(8),dimension(:),allocatable                      :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call d_unista(AIN,ndim,MASK)
    else
       call d_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end function f_d_uniq




  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  subroutine I_reshuffle(Ain,Index)
    integer,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    integer,dimension(size(Ain)) :: Aout
    integer                        :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
    Ain = Aout
  end subroutine I_reshuffle
  subroutine D_reshuffle(Ain,Index)
    real(8),dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    real(8),dimension(size(Ain)) :: Aout
    integer                      :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
    Ain = Aout
  end subroutine D_reshuffle
  subroutine Z_reshuffle(Ain,Index)
    complex(8),dimension(:)         :: Ain
    integer,dimension(size(Ain))    :: Index
    complex(8),dimension(size(Ain)) :: Aout
    integer                         :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
    Ain = Aout
  end subroutine Z_reshuffle
  subroutine L_reshuffle(Ain,Index)
    logical,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    logical,dimension(size(Ain)) :: Aout
    integer                      :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
    Ain = Aout
  end subroutine L_reshuffle

  function f_I_reshuffle(Ain,Index) result(Aout)
    integer,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    integer,dimension(size(Ain)) :: Aout
    integer                        :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function f_I_reshuffle
  function f_D_reshuffle(Ain,Index) result(Aout)
    real(8),dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    real(8),dimension(size(Ain)) :: Aout
    integer                      :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function f_D_reshuffle
  function f_Z_reshuffle(Ain,Index) result(Aout)
    complex(8),dimension(:)         :: Ain
    integer,dimension(size(Ain))    :: Index
    complex(8),dimension(size(Ain)) :: Aout
    integer                         :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function f_Z_reshuffle
  function f_L_reshuffle(Ain,Index) result(Aout)
    logical,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    logical,dimension(size(Ain)) :: Aout
    integer                      :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function f_L_reshuffle

  ! subroutine reshuffle(array,order) 
  !   real(8),dimension(:)                    :: array
  !   integer,dimension(size(array))          :: order
  !   real(8),dimension(size(array))          :: dummy
  !   integer                                 :: i,Lk
  !   forall(i=1:size(array))dummy(i)=array(order(i))
  !   array=dummy
  ! end subroutine reshuffle



  subroutine d_uniinv (xdont, igoest)
    real (kind=8), dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out)      :: igoest
    real (kind=8) :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine d_uniinv
  !
  subroutine r_uniinv (xdont, igoest)
    real, dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out) :: igoest
    real    :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine r_uniinv
  subroutine i_uniinv (xdont, igoest)
    ! __________________________________________________________
    !   uniinv = merge-sort inverse ranking of an array, with removal of
    !   duplicate entries.
    !   the routine is similar to pure merge-sort ranking, but on
    !   the last pass, it sets indices in igoest to the rank
    !   of the value in the ordered set with duplicates removed.
    !   for performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    integer, dimension (:), intent (in)  :: xdont
    integer, dimension (:), intent (out) :: igoest
    ! __________________________________________________________
    integer :: xtst, xdona, xdonb
    !
    ! __________________________________________________________
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine i_uniinv




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  function d_nearless (xval) result (d_nl)
    !  nearest value less than given value
    real (kind=8), intent (in) :: xval
    real (kind=8) :: d_nl
    d_nl = nearest (xval, -1.0d0)
    return
  end function d_nearless
  function r_nearless (xval) result (r_nl)
    !  nearest value less than given value
    real, intent (in) :: xval
    real :: r_nl
    r_nl = nearest (xval, -1.0)
    return
  end function r_nearless
  function i_nearless (xval) result (i_nl)
    !  nearest value less than given value
    integer, intent (in) :: xval
    integer :: i_nl
    i_nl = xval - 1
    return
  end function i_nearless




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  subroutine d_unista (xdont, nuni, mask)
    real(kind=8), dimension (:), intent (inout) :: xdont
    integer, intent (out)                       :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false. !modify the mask to that next acces is false
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine d_unista

  subroutine r_unista (xdont, nuni, mask)
    real, dimension (:), intent (inout) :: xdont
    integer, intent (out)               :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine r_unista

  subroutine i_unista (xdont, nuni, mask)
    integer, dimension (:), intent (inout)  :: xdont
    integer, intent (out) :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine i_unista



end module SF_MISC
