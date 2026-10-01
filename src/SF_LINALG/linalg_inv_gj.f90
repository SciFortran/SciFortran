!+-----------------------------------------------------------------+
!PURPOSE  : Linear equation solution by Gauss-Jordan elimination
! a is an N x N input coefficient matrix. On output, a is replaced 
!by its matrix inverse.
!+-----------------------------------------------------------------+
subroutine Dinv_gj(a)
  !This subroutine inverts in place a real square matrix :f:var:`a` with the Gauss-Jordan elimination with full
  !pivoting. The complex version :f:func_inline:`Zinv_gj` works alike. Both are instances of the generic interface
  !:f:func_inline:`inv_gj`.
  !
  !On output :f:var:`a` is replaced by its inverse. The program stops if the matrix is singular.
  !
  real(8), dimension(:,:), intent(inout) :: a  ! in: square matrix [n,n]; out: its inverse
  integer, dimension(size(a,1))      :: ipiv,indxr,indxc
  !these arrays are used for bookkeeping on the pivoting.
  !integer                            :: nn
  logical, dimension(size(a,1))      :: lpiv
  real(8)                            :: pivinv
  real(8), dimension(size(a,1))      :: dumc
  integer, target                    :: irc(2)
  integer                            :: i,l,n
  integer, pointer                   :: irow,icol
  real(8)                            :: zero=0.d0, one=1.d0
  n=size(a,1)
  irow => irc(1)
  icol => irc(2)
  ipiv=0
  do i=1,n
     !main loop over columns to be reduced.
     lpiv = (ipiv == 0)
     !begin search for a pivot element.
     irc=maxloc(abs(a),outerand(lpiv,lpiv))
     ipiv(icol)=ipiv(icol)+1
     if (ipiv(icol) > 1) stop 'gaussj:singular matrix (1)'
     !we now have the pivot element, so we interchange
     !rows, if needed, to put the pivot element on the diagonal. the columns
     !are not physically interchanged, only relabeled:
     !indxc(i),the column of the ith pivot element, is the ith column that is
     !reduced, while indxr(i) is the row in which that pivot element was
     !originally located. if indxr(i) = indxc(i) there is an implied column
     !interchange. with this form of bookkeeping, the inverse matrix will be
     !scrambled by
     !columns.
     if (irow /= icol) call swap(a(irow,:),a(icol,:))
     indxr(i)=irow !we are now ready to divide the pivot row by the pivot element,
     !located at irow and icol.
     indxc(i)=icol
     if (a(icol,icol) == zero) stop 'gaussj:singular matrix (2)'
     pivinv=one/a(icol,icol)
     a(icol,icol)= one !cmplx(one,zero)
     a(icol,:)=a(icol,:)*pivinv
     dumc=a(:,icol)
     !next, we reduce the rows, except for the pivot one, of course.
     a(:,icol)     = zero !cmplx
     a(icol,icol)  = pivinv
     a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
     a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
  end do
  !it only remains to unscramble the solution in view of the column
  !interchanges.
  !we do this by interchanging pairs of columns in the reverse order that the
  !permutation
  !was built up.
  do l=n,1,-1
     call swap(a(:,indxr(l)),a(:,indxc(l)))
  end do
contains
  function outerand(a,b)
    implicit none
    logical, dimension(:), intent(in)   :: a,b
    logical, dimension(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)).and.spread(b,dim=1,ncopies=size(a))
  end function outerand
  function outerprod(a,b)
    real(8), dimension(:), intent(in) :: a,b
    real(8), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod
end subroutine Dinv_gj

subroutine Zinv_gj(a)
  complex(8), dimension(:,:), intent(inout) :: a
  integer, dimension(size(a,1))      :: ipiv,indxr,indxc
  !these arrays are used for bookkeeping on the pivoting.
  !integer                            :: nn
  logical, dimension(size(a,1))      :: lpiv
  complex(8)                         :: pivinv
  complex(8), dimension(size(a,1))   :: dumc
  integer, target                    :: irc(2)
  integer                            :: i,l,n
  integer, pointer                   :: irow,icol
  real(8)                            :: zero=0.d0, one=1.d0
  n=size(a,1)
  irow => irc(1)
  icol => irc(2)
  ipiv=0
  do i=1,n
     !main loop over columns to be reduced.
     lpiv = (ipiv == 0)
     !begin search for a pivot element.
     irc=maxloc(abs(a),outerand(lpiv,lpiv))
     ipiv(icol)=ipiv(icol)+1
     if (ipiv(icol) > 1) stop 'gaussj:singular matrix (1)'
     !we now have the pivot element, so we interchange
     !rows, if needed, to put the pivot element on the diagonal. the columns
     !are not physically interchanged, only relabeled:
     !indxc(i),the column of the ith pivot element, is the ith column that is
     !reduced, while indxr(i) is the row in which that pivot element was
     !originally located. if indxr(i) = indxc(i) there is an implied column
     !interchange. with this form of bookkeeping, the inverse matrix will be
     !scrambled by
     !columns.
     if (irow /= icol) call swap(a(irow,:),a(icol,:))
     indxr(i)=irow !we are now ready to divide the pivot row by the pivot element,
     !located at irow and icol.
     indxc(i)=icol
     if (a(icol,icol) == zero) stop 'gaussj:singular matrix (2)'
     pivinv=one/a(icol,icol)
     a(icol,icol)=cmplx(one,zero,8)
     a(icol,:)=a(icol,:)*pivinv
     dumc=a(:,icol)
     !next, we reduce the rows, except for the pivot one, of course.
     a(:,icol)     = cmplx(zero,zero,8)
     a(icol,icol)  = pivinv
     a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
     a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
  end do
  !it only remains to unscramble the solution in view of the column
  !interchanges.
  !we do this by interchanging pairs of columns in the reverse order that the
  !permutation
  !was built up.
  do l=n,1,-1
     call swap(a(:,indxr(l)),a(:,indxc(l)))
  end do
contains
  function outerand(a,b)
    implicit none
    logical, dimension(:), intent(in)   :: a,b
    logical, dimension(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)).and.spread(b,dim=1,ncopies=size(a))
  end function outerand
  function outerprod(a,b)
    complex(8), dimension(:), intent(in) :: a,b
    complex(8), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod
end subroutine Zinv_gj








!+-----------------------------------------------------------------+
!program  : swap
!+-----------------------------------------------------------------+
subroutine swap_i(a,b)
  integer, intent(inout) :: a,b
  integer :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_i
!-----------------------------
subroutine swap_r(a,b)
  real(8), intent(inout) :: a,b
  real(8) :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_r
!-----------------------------
subroutine swap_rv(a,b)
  real(8), dimension(:), intent(inout) :: a,b
  real(8), dimension(size(a)) :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_rv
!-----------------------------
subroutine swap_z(a,b)
  complex(8), intent(inout) :: a,b
  complex(8) :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_z
!-----------------------------
subroutine swap_zv(a,b)
  complex(8), dimension(:), intent(inout) :: a,b
  complex(8), dimension(size(a)) :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_zv
!-----------------------------
subroutine swap_zm(a,b)
  complex(8), dimension(:,:), intent(inout) :: a,b
  complex(8), dimension(size(a,1),size(a,2)) :: dum
  dum=a
  a=b
  b=dum
end subroutine swap_zm
