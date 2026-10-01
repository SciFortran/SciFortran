!---------------------------------------------------------------------
!PURPOSE: Function to compute the tensor product (M1_kp_M2) of 
! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
! the number of rows and columns of the Matrix M1 and M2
!---------------------------------------------------------------------
function i_kronecker_product(A,B) result(AxB)
  !This function returns the Kronecker (tensor) product :math:`A \otimes B` of two matrices :f:var:`A` and :f:var:`B`,
  !i.e. the block matrix whose :math:`(i,j)` block is :math:`A_{ij} B`. The result has shape
  !:code:`[size(A,1)*size(B,1),size(A,2)*size(B,2)]`. This is the integer instance of the generic interfaces
  !:f:func_inline:`kron`, :f:func_inline:`kronecker_product` and of the operator :code:`.kx.`; the same product is
  !available for real (:f:func_inline:`d_kronecker_product`), real :math:`\otimes` complex
  !(:f:func_inline:`dc_kronecker_product`), complex :math:`\otimes` real (:f:func_inline:`cd_kronecker_product`) and
  !complex (:f:func_inline:`c_kronecker_product`) matrices.
  !
  integer,intent(in) :: A(:,:)                                        ! first factor, [rowA,colA]
  integer,intent(in) :: B(:,:)                                        ! second factor, [rowB,colB]
  integer            :: i,j
  integer            :: rowA,colA
  integer            :: rowB,colB
  integer            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))  ! Kronecker product, [rowA*rowB,colA*colB]
  AxB = 0
  rowA=size(A,1) ; colA=size(A,2)
  rowB=size(B,1) ; colB=size(B,2)
  forall(i=1:rowA,j=1:colA)
     AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  end forall
end function i_kronecker_product
!
function d_kronecker_product(A,B) result(AxB)
  real(8),intent(in) :: A(:,:), B(:,:)
  integer            :: i,j
  integer            :: rowA,colA
  integer            :: rowB,colB
  real(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
  AxB = 0d0
  rowA=size(A,1) ; colA=size(A,2)
  rowB=size(B,1) ; colB=size(B,2)
  forall(i=1:rowA,j=1:colA)
     AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  end forall
end function d_kronecker_product
!
function dc_kronecker_product(A,B) result(AxB)
  real(8),intent(in) :: A(:,:)
  complex(8),intent(in) :: B(:,:)
  integer            :: i,j
  integer            :: rowA,colA
  integer            :: rowB,colB
  complex(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
  AxB = zero
  rowA=size(A,1) ; colA=size(A,2)
  rowB=size(B,1) ; colB=size(B,2)
  forall(i=1:rowA,j=1:colA)
     AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  end forall
end function dc_kronecker_product
!
function cd_kronecker_product(A,B) result(AxB)
  complex(8),intent(in) :: A(:,:)
  real(8),intent(in)    :: B(:,:)
  integer               :: i,j
  integer               :: rowA,colA
  integer               :: rowB,colB
  complex(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
  AxB = zero
  rowA=size(A,1) ; colA=size(A,2)
  rowB=size(B,1) ; colB=size(B,2)
  forall(i=1:rowA,j=1:colA)
     AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  end forall
end function cd_kronecker_product
!
function c_kronecker_product(A,B) result(AxB)
  complex(8),intent(in) :: A(:,:), B(:,:)
  integer               :: i,j
  integer               :: rowA,colA
  integer               :: rowB,colB
  complex(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
  AxB = zero
  rowA=size(A,1) ; colA=size(A,2)
  rowB=size(B,1) ; colB=size(B,2)
  forall(i=1:rowA,j=1:colA)
     AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  end forall
end function c_kronecker_product






!+-----------------------------------------------------------------------------+!
!PURPOSE: Form a matrix A(:,:) from the outerproduct of two 1d arrays:
! A(i,j) = a_i*b_j
!+-----------------------------------------------------------------------------+!
function outerprod_d(a,b) result(outerprod)
  real(8), dimension(:), intent(in)   :: a,b
  real(8), dimension(size(a),size(b)) :: outerprod
  outerprod = spread(a,dim=2,ncopies=size(b)) * &
       spread(b,dim=1,ncopies=size(a))
end function outerprod_d
function outerprod_c(a,b) result(outerprod)
  complex(8), dimension(:), intent(in)   :: a,b
  complex(8), dimension(size(a),size(b)) :: outerprod
  outerprod = spread(a,dim=2,ncopies=size(b)) * &
       spread(b,dim=1,ncopies=size(a))
end function outerprod_c







!+-----------------------------------------------------------------------------+!
!PURPOSE:  cross or vector product for 2d and 3d vectors. 
!+-----------------------------------------------------------------------------+!
function cross_2d_d(a,b) result(c)
  real(8),dimension(2) :: a,b
  real(8)              :: c
  c = a(1)*b(2) - a(2)*b(1)
end function cross_2d_d
function cross_2d_c(a,b) result(c)
  complex(8),dimension(2) :: a,b
  complex(8)              :: c
  c = a(1)*b(2) - a(2)*b(1)
end function cross_2d_c
!
function cross_3d_d(a,b) result(c)
  real(8),dimension(3) :: a,b
  real(8),dimension(3) :: c
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_3d_d
function cross_3d_c(a,b) result(c)
  complex(8),dimension(3) :: a,b
  complex(8),dimension(3) :: c
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_3d_c





!+-----------------------------------------------------------------------------+!
!PURPOSE: evaluate the S3 product A.(BxC) for 3d vectors
!+-----------------------------------------------------------------------------+!
function s3_product_d(a,b,c) result(s3)
  real(8),dimension(3),intent(in) :: a,b,c
  real(8)                         :: s3
  s3 = dot_product(a,cross_product(b, c))
end function s3_product_d
function s3_product_c(a,b,c) result(s3)
  complex(8),dimension(3),intent(in) :: a,b,c
  real(8)                            :: s3
  s3 = dot_product(a,cross_product(b, c))
end function s3_product_c
