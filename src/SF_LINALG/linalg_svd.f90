!-------------------------------------------------------------------------------------------
!PURPOSE: compute the singular value decomposition A = U sigma Vtransp / U sigma V^H of 
! real/complex matrix A
!-------------------------------------------------------------------------------------------
subroutine dsvd(A, s, U, Vtransp)
  !This subroutine computes the singular value decomposition :math:`A = U \Sigma V^T` of a real :math:`m\times n`
  !matrix :f:var:`A`, with full :math:`U` and :math:`V^T`, using the LAPACK routine :f:func_inline:`dgesvd`. The
  !complex version :f:func_inline:`zsvd` uses :f:func_inline:`zgesvd` and returns :math:`V^H`. Both are instances of
  !the generic interface :f:func_inline:`svd`.
  !
  !The input matrix is copied and left untouched. :math:`\Sigma` is the :math:`m\times n` matrix with diagonal
  !elements :math:`\Sigma_{ii} = s_i`, the singular values :f:var:`s` being in descending order. Note that the
  !routine returns the transpose :f:var:`Vtransp` (conjugate transpose in the complex case) of :math:`V`, not :math:`V`.
  !
  !The program stops if :f:var:`U` does not have shape :code:`[m,m]`, if :f:var:`Vtransp` does not have shape
  !:code:`[n,n]`, or if the LAPACK routine returns an error code.
  !
  real(8), intent(in)  :: A(:,:)        ! real matrix [m,n], not modified
  real(8), intent(out) :: s(:)          ! singular values, descending, size min(m,n)
  real(8), intent(out) :: U(:,:)        ! left singular vectors (columns), [m,m]
  real(8), intent(out) :: Vtransp(:,:)  ! V^T (V^H in the complex case), [n,n]
  integer              :: info, lwork, m, n, ldu
  real(8), allocatable :: work(:), At(:,:)
  m = size(A(:,1))  ! = lda
  n = size(A(1,:))
  ldu = m
  allocate(At(m,n))
  At(:,:) = A(:,:)  ! use a temporary as dgesvd destroys its input
  call assert_shape(U, [m, m], "svd", "U")
  call assert_shape(Vtransp, [n, n], "svd", "Vtransp")
  ! query optimal lwork and allocate workspace:
  allocate(work(1))
  call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1, info)
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, lwork, info)
  if(info /= 0) then
     print *, "dgesvd returned info = ", info
     if(info < 0) then
        print *, "the ", -info, "-th argument had an illegal value"
     else
        print *, "DBDSQR did not converge, there are ", info
        print *, "superdiagonals of an intermediate bidiagonal form B"
        print *, "did not converge to zero. See the description of WORK"
        print *, "in DGESVD's man page for details."
     endif
     stop 'dsvd errpr: dgesvd'
  endif
end subroutine dsvd

subroutine zsvd(A, s, U, Vtransp)
  ! complex m x m matrix A
  ! U is m x min(m, n)
  ! Vtransp is n x n
  ! sigma is m x n with with sigma_ii = s_i
  ! note that this routine returns V^H, not V!
  complex(8), intent(in)  :: A(:,:)
  real(8), intent(out)    :: s(:)
  complex(8), intent(out) :: U(:,:), Vtransp(:,:)
  integer                 :: info, lwork, m, n, ldu, lrwork
  real(8), allocatable    :: rwork(:)
  complex(8), allocatable :: work(:), At(:,:)
  ! TODO: check shapes here and in other routines?
  m = size(A(:,1))  ! = lda
  n = size(A(1,:))
  ldu = m
  lrwork = 5*min(m,n)
  allocate(rwork(lrwork), At(m,n))
  At(:,:) = A(:,:)  ! use a temporary as zgesvd destroys its input
  call assert_shape(U, [m, m], "svd", "U")
  call assert_shape(Vtransp, [n, n], "svd", "Vtransp")
  ! query optimal lwork and allocate workspace:
  allocate(work(1))
  call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1,&
       rwork, info)
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, &
       lwork, rwork, info)
  if(info /= 0) then
     print *, "zgesvd returned info = ", info
     if(info < 0) then
        print *, "the ", -info, "-th argument had an illegal value"
     else
        print *, "ZBDSQR did not converge, there are ", info
        print *, "superdiagonals of an intermediate bidiagonal form B"
        print *, "did not converge to zero. See the description of WORK"
        print *, "in DGESVD's man page for details."
     endif
     stop 'zsvd error: zgesvd'
  endif
end subroutine zsvd
