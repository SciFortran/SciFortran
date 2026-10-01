
function dlstsq(A, b) result(x)
  !This function returns the least-squares solution :math:`x = \mathrm{argmin}_x \| A x - b \|_2` of the real linear
  !system :math:`A x = b`, with :math:`A` an :math:`m\times n` matrix, computed with the LAPACK routine
  !:f:func_inline:`dgelsy`, which uses a complete orthogonal factorization with column pivoting and is thus suitable
  !also for rank-deficient matrices (the effective rank is determined with :code:`rcond = 0`). The complex version
  !:f:func_inline:`zlstsq` uses :f:func_inline:`zgelsy`. Both are instances of the generic interface
  !:f:func_inline:`lstsq`.
  !
  !The input arrays are copied and left untouched. The solution is allocated with size :code:`n`. The program stops
  !if the LAPACK routine returns an error code.
  !
  real(8), intent(in)  :: A(:,:)  ! real matrix [m,n], not modified
  real(8), intent(in)  :: b(:)    ! right-hand side, size m, not modified
  real(8), allocatable :: x(:)    ! least-squares solution, size n
  integer              :: info, ldb, lwork, m, n, rank
  real(8)              :: rcond
  real(8), allocatable :: work(:), At(:,:), Bt(:,:)
  integer, allocatable :: jpvt(:)
  m = size(A(:,1)) ! = lda
  n = size(A(1,:))
  ldb = size(b)
  allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1))
  call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
       -1, info)  ! query optimal workspace size
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))  ! allocate with ideal size
  rcond = 0d0
  jpvt(:) = 0
  Bt(:,1) = b(:)  ! only one right-hand side
  At(:,:) = A(:,:)
  call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
       lwork, info)
  if(info /= 0) then
     print *, "dgelsy returned info = ", info
     print *, "the ", -info, "-th argument had an illegal value"
     stop 'dlstsq error: dgelsy'
  endif
  x(:) = Bt(1:n,1)
end function dlstsq

function zlstsq(A, b) result(x)
  complex(8), intent(in)  :: A(:,:), b(:)
  complex(8), allocatable :: x(:)
  integer                 :: info, ldb, lwork, m, n, rank
  real(8)                 :: rcond
  complex(8), allocatable :: At(:,:), Bt(:,:), work(:)
  real(8), allocatable    :: rwork(:)
  integer, allocatable    :: jpvt(:)
  m = size(A(:,1)) ! = lda
  n = size(A(1,:))
  ldb = size(b)
  allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1), rwork(2*n))
  call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
       -1, rwork, info)  ! query optimal workspace size
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))  ! allocate with ideal size
  rcond = 0d0
  jpvt(:) = 0
  Bt(:,1) = b(:)  ! only one right-hand side
  At(:,:) = A(:,:)
  call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
       lwork, rwork, info)
  if(info /= 0) then
     print *, "zgelsy returned info = ", info
     print *, "the ", -info, "-th argument had an illegal value"
     stop 'zlstsq error: zgelsy'
  endif
  x(:) = Bt(1:n,1)
end function zlstsq
