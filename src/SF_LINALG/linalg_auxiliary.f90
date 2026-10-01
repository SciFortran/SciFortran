!-------------------------------------------------------------------------------------------
!PURPOSE: compute the determinant of a matrix using an LU factorization
!-------------------------------------------------------------------------------------------
function ddet(A) result(x)
  !This function returns the determinant of a real square matrix :f:var:`A`, computed from its LU factorization
  !(:f:func_inline:`dgetrf`) as the product of the diagonal elements of :math:`U`, with the sign changes due to the
  !row permutations. The complex version :f:func_inline:`zdet` uses :f:func_inline:`zgetrf`. Both are instances of the
  !generic interface :f:func_inline:`det`.
  !
  !The input matrix is copied and left untouched. The program stops if :f:var:`A` is not square, or if the LAPACK
  !routine returns an error code, which includes the case of an exactly singular matrix.
  !
  real(8), intent(in)  :: A(:, :)  ! square matrix [n,n], not modified
  real(8)              :: x        ! determinant of A
  integer              :: i
  integer              :: info, n
  integer, allocatable :: ipiv(:)
  real(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call dgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "dgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'det: dgetrf error'
  end if
  ! At now contains the LU of the factorization A = PLU
  ! as L has unit diagonal entries, the determinant can be computed
  ! from the product of U's diagonal entries. Additional sign changes
  ! stemming from the permutations P have to be taken into account as well.
  x = 1d0
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function ddet
!
function zdet(A) result(x)
  ! compute the determinant of a complex matrix using an LU factorization
  complex(8), intent(in)  :: A(:, :)
  complex(8)              :: x
  integer                 :: i
  integer                 :: info, n
  integer, allocatable    :: ipiv(:)
  complex(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call zgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "zgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'zdet error: zgetrf '
  end if
  ! for details on the computation, compare the comment in ddet().
  x = one
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function zdet





!-------------------------------------------------------------------------------------------
!PURPOSE:  construct real matrix from diagonal elements
!-------------------------------------------------------------------------------------------
pure function ddiag(x) result(A)
  !This function builds a real square matrix with the elements of the array :f:var:`x` on the main diagonal and zeros
  !elsewhere. The complex version is :f:func_inline:`zdiag`. Both are instances of the generic interface
  !:f:func_inline:`diag`, the inverse operation being :f:func_inline:`diagonal`.
  !
  real(8), intent(in)  :: x(:)    ! diagonal elements, size n
  real(8), allocatable :: A(:,:)  ! diagonal matrix, A(i,i) = x(i), [n,n]
  integer              :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = 0d0
  forall(i=1:n) A(i,i) = x(i)
end function ddiag
!-------------------------------------------------------------------------------------------
!PURPOSE:  construct complex matrix from diagonal elements
!-------------------------------------------------------------------------------------------!
pure function zdiag(x) result(A)
  complex(8), intent(in)  :: x(:)
  complex(8), allocatable :: A(:,:)
  integer                 :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = zero
  forall(i=1:n) A(i,i) = x(i)
end function zdiag





!-------------------------------------------------------------------------------------------
!PURPOSE:  return the diagonal of a matrix [real]
!-------------------------------------------------------------------------------------------
pure function d_diagonal(A) result(dd)
  !This function returns the main diagonal of a real matrix :f:var:`A`. The complex version is
  !:f:func_inline:`z_diagonal`. Both are instances of the generic interface :f:func_inline:`diagonal`, the inverse
  !operation being :f:func_inline:`diag`. The matrix is expected to be square (the size of the result is
  !:code:`size(A,1)`).
  !
  real(8),intent(in)           :: A(:,:)  ! square matrix [n,n]
  real(8),dimension(size(A,1)) :: dd      ! diagonal elements, dd(i) = A(i,i), size n
  integer                      :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function d_diagonal
!-------------------------------------------------------------------------------------------
!PURPOSE:  return the diagonal of a matrix [complex]
!-------------------------------------------------------------------------------------------
pure function z_diagonal(A) result(dd)
  complex(8),intent(in)           :: A(:,:)
  complex(8),dimension(size(A,1)) :: dd
  integer                         :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function z_diagonal




!-------------------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal [real]
!-------------------------------------------------------------------------------------------
pure function dtrace(A) result(t)
  !This function returns the trace of a real matrix :f:var:`A`, i.e. the sum of the elements along the main diagonal.
  !The complex version is :f:func_inline:`ztrace`. Both are instances of the generic interface :f:func_inline:`trace`.
  !For non-square matrices the sum runs up to :code:`min(m,n)`.
  !
  real(8), intent(in) :: A(:,:)  ! real matrix [m,n]
  real(8)             :: t       ! trace of A
  integer             :: i
  t = 0d0
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function dtrace
!-------------------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal [complex]
!-------------------------------------------------------------------------------------------
pure function ztrace(A) result(t)
  complex(8), intent(in) :: A(:,:)
  complex(8)             :: t
  integer                :: i
  t = zero
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function ztrace





!-------------------------------------------------------------------------------------------
!PURPOSE:  Returns the identity matrix of size n x n and type real.
!-------------------------------------------------------------------------------------------
pure function deye_matrix(n) result(A)
  !This function returns the real identity matrix of size :math:`n\times n`. It is an instance of the generic
  !interfaces :f:func_inline:`deye` and :f:func_inline:`eye`; the complex version is :f:func_inline:`zeye_matrix`.
  !
  integer, intent(in) :: n        ! size of the matrix
  real(8)             :: A(n, n)  ! identity matrix [n,n] or matrix element (scalar)
  integer             :: i
  A = 0d0
  do i = 1, n
     A(i,i) = 1d0
  end do
end function deye_matrix

pure function zeye_matrix(n) result(A)
  !This function returns the complex identity matrix of size :math:`n\times n`. It is an instance of the generic
  !interface :f:func_inline:`zeye`; the real version is :f:func_inline:`deye_matrix`.
  !
  integer, intent(in) :: n        ! size of the matrix
  complex(8)          :: A(n, n)  ! identity matrix [n,n] or matrix element (scalar)
  integer             :: i
  A = zero
  do i = 1, n
     A(i,i) = one
  end do
end function zeye_matrix

pure function deye_indices(i,j) result(a)
  !This function returns the real Kronecker delta :math:`\delta_{ij}`: 1 if :code:`i==j`, 0 otherwise. It is an
  !instance of the generic interfaces :f:func_inline:`deye` and :f:func_inline:`eye`; the complex version is
  !:f:func_inline:`zeye_indices`.
  !
  integer, intent(in) :: i  ! first index
  integer, intent(in) :: j  ! second index
  real(8)             :: a  ! delta_ij
  a = 0d0
  if(i==j)a=1d0
end function deye_indices

pure function zeye_indices(i,j) result(a)
  !This function returns the complex Kronecker delta :math:`\delta_{ij}`: 1 if :code:`i==j`, 0 otherwise. It is an
  !instance of the generic interface :f:func_inline:`zeye`; the real version is :f:func_inline:`deye_indices`.
  !
  integer, intent(in) :: i  ! first index
  integer, intent(in) :: j  ! second index
  complex(8)          :: a  ! delta_ij
  a = zero
  if(i==j)a=one
end function zeye_indices




!-------------------------------------------------------------------------------------------
!PURPOSE:  Returns an array of zeros of specified size from 1 to 7 dimension
!-------------------------------------------------------------------------------------------
pure function zzeros_1(n) result(A)
  !This function returns a rank-1 complex(8) array of shape :code:`[n]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n     ! extent of dimension 1
  complex(8)          :: A(n)  ! complex array, [n]
  A = zero
end function zzeros_1
!
pure function zzeros_2(n1,n2) result(A)
  !This function returns a rank-2 complex(8) array of shape :code:`[n1,n2]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1        ! extent of dimension 1
  integer, intent(in) :: n2        ! extent of dimension 2
  complex(8)          :: A(n1,n2)  ! complex array, [n1,n2]
  A = zero
end function zzeros_2
!
pure function zzeros_3(n1,n2,n3) result(A)
  !This function returns a rank-3 complex(8) array of shape :code:`[n1,n2,n3]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1           ! extent of dimension 1
  integer, intent(in) :: n2           ! extent of dimension 2
  integer, intent(in) :: n3           ! extent of dimension 3
  complex(8)          :: A(n1,n2,n3)  ! complex array, [n1,n2,n3]
  A = zero
end function zzeros_3
!
pure function zzeros_4(n1,n2,n3,n4) result(A)
  !This function returns a rank-4 complex(8) array of shape :code:`[n1,n2,n3,n4]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1              ! extent of dimension 1
  integer, intent(in) :: n2              ! extent of dimension 2
  integer, intent(in) :: n3              ! extent of dimension 3
  integer, intent(in) :: n4              ! extent of dimension 4
  complex(8)          :: A(n1,n2,n3,n4)  ! complex array, [n1,n2,n3,n4]
  A = zero
end function zzeros_4
!
pure function zzeros_5(n1,n2,n3,n4,n5) result(A)
  !This function returns a rank-5 complex(8) array of shape :code:`[n1,n2,n3,n4,n5]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                 ! extent of dimension 1
  integer, intent(in) :: n2                 ! extent of dimension 2
  integer, intent(in) :: n3                 ! extent of dimension 3
  integer, intent(in) :: n4                 ! extent of dimension 4
  integer, intent(in) :: n5                 ! extent of dimension 5
  complex(8)          :: A(n1,n2,n3,n4,n5)  ! complex array, [n1,n2,n3,n4,n5]
  A = zero
end function zzeros_5
!
pure function zzeros_6(n1,n2,n3,n4,n5,n6) result(A)
  !This function returns a rank-6 complex(8) array of shape :code:`[n1,n2,n3,n4,n5,n6]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                    ! extent of dimension 1
  integer, intent(in) :: n2                    ! extent of dimension 2
  integer, intent(in) :: n3                    ! extent of dimension 3
  integer, intent(in) :: n4                    ! extent of dimension 4
  integer, intent(in) :: n5                    ! extent of dimension 5
  integer, intent(in) :: n6                    ! extent of dimension 6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)  ! complex array, [n1,n2,n3,n4,n5,n6]
  A = zero
end function zzeros_6
!
pure function zzeros_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  !This function returns a rank-7 complex(8) array of shape :code:`[n1,n2,n3,n4,n5,n6,n7]` with all the elements set to 0.
  !It is an instance of the generic interface :f:func_inline:`zeros`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                       ! extent of dimension 1
  integer, intent(in) :: n2                       ! extent of dimension 2
  integer, intent(in) :: n3                       ! extent of dimension 3
  integer, intent(in) :: n4                       ! extent of dimension 4
  integer, intent(in) :: n5                       ! extent of dimension 5
  integer, intent(in) :: n6                       ! extent of dimension 6
  integer, intent(in) :: n7                       ! extent of dimension 7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)  ! complex array, [n1,n2,n3,n4,n5,n6,n7]
  A = zero
end function zzeros_7



pure function zones_1(n) result(A)
  !This function returns a rank-1 complex(8) array of shape :code:`[n]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n     ! extent of dimension 1
  complex(8)          :: A(n)  ! complex array, [n]
  A = one
end function zones_1
!
pure function zones_2(n1,n2) result(A)
  !This function returns a rank-2 complex(8) array of shape :code:`[n1,n2]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1        ! extent of dimension 1
  integer, intent(in) :: n2        ! extent of dimension 2
  complex(8)          :: A(n1,n2)  ! complex array, [n1,n2]
  A = one
end function zones_2
!
pure function zones_3(n1,n2,n3) result(A)
  !This function returns a rank-3 complex(8) array of shape :code:`[n1,n2,n3]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1           ! extent of dimension 1
  integer, intent(in) :: n2           ! extent of dimension 2
  integer, intent(in) :: n3           ! extent of dimension 3
  complex(8)          :: A(n1,n2,n3)  ! complex array, [n1,n2,n3]
  A = one
end function zones_3
!
pure function zones_4(n1,n2,n3,n4) result(A)
  !This function returns a rank-4 complex(8) array of shape :code:`[n1,n2,n3,n4]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1              ! extent of dimension 1
  integer, intent(in) :: n2              ! extent of dimension 2
  integer, intent(in) :: n3              ! extent of dimension 3
  integer, intent(in) :: n4              ! extent of dimension 4
  complex(8)          :: A(n1,n2,n3,n4)  ! complex array, [n1,n2,n3,n4]
  A = one
end function zones_4
!
pure function zones_5(n1,n2,n3,n4,n5) result(A)
  !This function returns a rank-5 complex(8) array of shape :code:`[n1,n2,n3,n4,n5]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                 ! extent of dimension 1
  integer, intent(in) :: n2                 ! extent of dimension 2
  integer, intent(in) :: n3                 ! extent of dimension 3
  integer, intent(in) :: n4                 ! extent of dimension 4
  integer, intent(in) :: n5                 ! extent of dimension 5
  complex(8)          :: A(n1,n2,n3,n4,n5)  ! complex array, [n1,n2,n3,n4,n5]
  A = one
end function zones_5
!
pure function zones_6(n1,n2,n3,n4,n5,n6) result(A)
  !This function returns a rank-6 complex(8) array of shape :code:`[n1,n2,n3,n4,n5,n6]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                    ! extent of dimension 1
  integer, intent(in) :: n2                    ! extent of dimension 2
  integer, intent(in) :: n3                    ! extent of dimension 3
  integer, intent(in) :: n4                    ! extent of dimension 4
  integer, intent(in) :: n5                    ! extent of dimension 5
  integer, intent(in) :: n6                    ! extent of dimension 6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)  ! complex array, [n1,n2,n3,n4,n5,n6]
  A = one
end function zones_6
!
pure function zones_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  !This function returns a rank-7 complex(8) array of shape :code:`[n1,n2,n3,n4,n5,n6,n7]` with all the elements set to 1.
  !It is an instance of the generic interface :f:func_inline:`ones`, which covers ranks from 1 to 7.
  !
  integer, intent(in) :: n1                       ! extent of dimension 1
  integer, intent(in) :: n2                       ! extent of dimension 2
  integer, intent(in) :: n3                       ! extent of dimension 3
  integer, intent(in) :: n4                       ! extent of dimension 4
  integer, intent(in) :: n5                       ! extent of dimension 5
  integer, intent(in) :: n6                       ! extent of dimension 6
  integer, intent(in) :: n7                       ! extent of dimension 7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)  ! complex array, [n1,n2,n3,n4,n5,n6,n7]
  A = one
end function zones_7





