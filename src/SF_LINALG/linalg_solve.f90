subroutine Dsolve_1rhs(A,b,trans)
  !This subroutine solves the real linear system :math:`A x = b` with a single right-hand side using the LU
  !factorization, :f:func_inline:`dgetrf` + :f:func_inline:`dgetrs`. The complex version :f:func_inline:`Zsolve_1rhs`
  !uses :f:func_inline:`zgetrf` + :f:func_inline:`zgetrs`. Both are instances of the generic interface
  !:f:func_inline:`solve`.
  !
  !The solution :math:`x` overwrites :f:var:`b`. Note that :f:var:`A` is factorized in place by :f:func_inline:`dgetrf`:
  !despite the :code:`intent(in)` it contains the LU factors on exit, so a copy must be passed if the original matrix
  !is needed afterwards. If :f:var:`trans` is :code:`'T'` or :code:`'C'` the system :math:`A^T x = b` is solved instead.
  !The program stops if a LAPACK routine returns an error code, e.g. for a singular matrix.
  !
  real(8),dimension(:,:),intent(in)  :: A      ! square matrix [n,n]; overwritten by its LU factors
  real(8),dimension(:),intent(inout) :: b      ! in: right-hand side, size n; out: solution x
  real(8),dimension(:,:),allocatable :: b_
  character(len=1),optional          :: trans  ! optional: 'N' (default) solve A x = b, 'T'/'C' solve A^T x = b
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= 1
  allocate(ipvt(min(m,n)))
  call dgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  allocate(b_(ldb,nrhs))
  b_(:,1)=b
  call dgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  b=b_(:,1)
  deallocate(ipvt,b_)
end subroutine dsolve_1rhs

subroutine Zsolve_1rhs(A,b,trans)
  complex(8),dimension(:,:),intent(in)  :: A
  complex(8),dimension(:),intent(inout) :: b
  complex(8),dimension(:,:),allocatable :: b_
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= 1
  allocate(ipvt(min(m,n)))
  call zgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  allocate(b_(ldb,nrhs))
  b_(:,1)=b
  call zgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  b=b_(:,1)
  deallocate(ipvt,b_)
end subroutine Zsolve_1rhs

subroutine Dsolve_Mrhs(A,b,trans)
  !This subroutine solves the real linear system :math:`A X = B` with several right-hand sides, stored as the columns
  !of :f:var:`b`, using the LU factorization, :f:func_inline:`dgetrf` + :f:func_inline:`dgetrs`. The complex version
  !:f:func_inline:`Zsolve_Mrhs` uses :f:func_inline:`zgetrf` + :f:func_inline:`zgetrs`. Both are instances of the
  !generic interface :f:func_inline:`solve`.
  !
  !The solutions :math:`X` overwrite :f:var:`b`. Note that :f:var:`A` is factorized in place by :f:func_inline:`dgetrf`:
  !despite the :code:`intent(in)` it contains the LU factors on exit, so a copy must be passed if the original matrix
  !is needed afterwards. If :f:var:`trans` is :code:`'T'` or :code:`'C'` the system :math:`A^T X = B` is solved instead.
  !The program stops if a LAPACK routine returns an error code, e.g. for a singular matrix.
  !
  real(8),dimension(:,:),intent(in)    :: A    ! square matrix [n,n]; overwritten by its LU factors
  real(8),dimension(:,:),intent(inout) :: b    ! in: right-hand sides, [n,nrhs]; out: solutions X
  character(len=1),optional          :: trans  ! optional: 'N' (default) solve A X = B, 'T'/'C' solve A^T X = B
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B,1))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= size(B,2)
  allocate(ipvt(min(m,n)))
  call dgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  call dgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  deallocate(ipvt)
end subroutine dsolve_Mrhs

subroutine Zsolve_Mrhs(A,b,trans)
  complex(8),dimension(:,:),intent(in)    :: A
  complex(8),dimension(:,:),intent(inout) :: b
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B,1))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= size(B,2)   
  allocate(ipvt(n))
  call zgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  lda=n ; ldb=n ; nrhs=size(b,2)
  call zgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  deallocate(ipvt)
end subroutine Zsolve_Mrhs
