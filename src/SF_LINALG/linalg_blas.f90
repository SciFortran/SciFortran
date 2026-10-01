subroutine d_matmul(A,B,C,alfa,beta)
  !This subroutine computes the real matrix-matrix product :math:`C = \alpha A B + \beta C` with the BLAS routine
  !:f:func_inline:`dgemm`, :math:`A` being :math:`N\times K`, :math:`B` :math:`K\times M` and :math:`C` :math:`N\times M`.
  !The complex version is :f:func_inline:`z_matmul`. Both are instances of the generic interface
  !:f:func_inline:`mat_product`; see also the operator :code:`.x.` (:f:func_inline:`d_matmul_`).
  !
  !:f:var:`A` and :f:var:`B` are not modified. By default :math:`\alpha=1` and :math:`\beta=0`, so that :f:var:`C`
  !contains the plain product :math:`A B` on output. The program stops if :f:var:`B` or :f:var:`C` do not have the
  !shapes :code:`[K,M]` and :code:`[N,M]`.
  !
  real(8),dimension(:,:),intent(inout) :: A     ! matrix [N,K], not modified
  real(8),dimension(:,:),intent(inout) :: B     ! matrix [K,M], not modified
  real(8),dimension(:,:),intent(inout) :: C     ! in/out: C = alfa*A*B + beta*C, [N,M]
  real(8),optional                     :: alfa  ! optional scalar, default 1
  real(8),optional                     :: beta  ! optional scalar, default 0
  real(8)                              :: alfa_,beta_
  integer                              :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 ;if(present(alfa))alfa_=alfa
  beta_     = 0d0 ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call DGEMM('N', 'N', N, M, K, alfa_, A, N, B, K, beta_, C, M)
  !
  return
end subroutine d_matmul

subroutine z_matmul(A,B,C,alfa,beta)
  complex(8),dimension(:,:),intent(inout) :: A ![N,K]
  complex(8),dimension(:,:),intent(inout) :: B ![K,M]
  complex(8),dimension(:,:),intent(inout) :: C ![N,M]
  complex(8),optional                     :: alfa,beta
  complex(8)                              :: alfa_,beta_
  integer                                 :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0) ;if(present(alfa))alfa_=alfa
  beta_     = dcmplx(0d0,0d0) ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "z_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "z_matmul error: C has illegal shape"
  !
  call ZGEMM('N', 'N', N, M, K, alfa_, A, N, B, K, beta_, C, N)
  !
  return
end subroutine z_matmul




!############## OVERLOAD MATMUL OPERATOR --> .x. #################


function d_matmul_(A,B) result(C)
  !This function returns the real matrix-matrix product :math:`C = A B`, computed with the BLAS routine
  !:f:func_inline:`dgemm`, :math:`A` being :math:`N\times K` and :math:`B` :math:`K\times M`. It is the function behind
  !the operator :code:`.x.`, i.e. :code:`C = A .x. B`; the complex version is :f:func_inline:`z_matmul_`.
  !
  !The program stops if :f:var:`B` does not have the shape :code:`[K,M]`.
  !
  real(8),dimension(:,:),intent(in)      :: A  ! matrix [N,K]
  real(8),dimension(:,:),intent(in)      :: B  ! matrix [K,M]
  real(8),dimension(size(A,1),size(B,2)) :: C  ! product A*B, [N,M]
  integer                                :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call DGEMM('N', 'N', N, M, K, 1d0, A, N, B, K, 0d0, C, M)
  !
  return
end function d_matmul_

function z_matmul_(A,B) result(C)
  complex(8),dimension(:,:),intent(in)      :: A ![N,K]
  complex(8),dimension(:,:),intent(in)      :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  integer                                :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call ZGEMM('N', 'N', N, M, K, dcmplx(1d0,0d0), A, N, B, K, dcmplx(0d0,0d0), C, N)
  !
  return
end function z_matmul_




