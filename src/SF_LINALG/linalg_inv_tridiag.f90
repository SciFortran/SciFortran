!+-----------------------------------------------------------------------------+!
!PURPOSE: return the N diagonal elements of the inverse matrix.
! b = sub-diagonal
! d = main-diagonal
! a = over-diagonal
!+-----------------------------------------------------------------------------+!
subroutine d_invert_tridiag_matrix(N,sub_,diag_,over_,Inv)
  !This subroutine returns the :code:`N` diagonal elements of the inverse of a real tridiagonal matrix, given by its
  !sub-, main and over-diagonals, without building the full inverse. The complex version is
  !:f:func_inline:`c_invert_tridiag_matrix`. Both are instances of the generic interface :f:func_inline:`inv_tridiag`.
  !
  !The diagonal of the inverse is obtained from a downward and an upward recursion,
  !:math:`d^L_i = d_i - s_{i-1} o_{i-1} / d^L_{i-1}`, :math:`d^R_i = d_i - o_i s_i / d^R_{i+1}`, as
  !:math:`(A^{-1})_{ii} = 1/(d^L_i + d^R_i - d_i)`, where :math:`d`, :math:`s`, :math:`o` are the main, sub- and
  !over-diagonal, :code:`sub_(i) = A(i+1,i)` and :code:`over_(i) = A(i,i+1)`. The program stops if a zero pivot is
  !found, i.e. if the matrix is ill-conditioned.
  !
  integer                         :: i,j
  integer                         :: N      ! size of the matrix
  real(8),dimension(N)            :: diag_  ! main diagonal, size N
  real(8),dimension(N-1)          :: sub_   ! sub-diagonal, A(i+1,i), size N-1
  real(8),dimension(N-1)          :: over_  ! over-diagonal, A(i,i+1), size N-1
  real(8),dimension(N)            :: Inv    ! diagonal elements of the inverse matrix, size N
  real(8)                         :: Foo,Cleft,Cright
  real(8),dimension(N)            :: dleft
  real(8),dimension(N)            :: dright
  !
  !DOWNWARD:
  dleft(1) = diag_(1)
  if(dleft(1)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  do i=2,N
     foo      = 1d0/dleft(i-1)
     cleft    = sub_(i-1)*foo
     dleft(i) = diag_(i) - cleft*over_(i-1)
     if(dleft(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  enddo
  !
  !UPWARD:
  dright(N) = diag_(N)
  if(dright(N)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  do i=N-1,1,-1
     foo       = 1d0/dright(i+1)
     cright    = over_(i)*foo
     dright(i) = diag_(i) - cright*sub_(i)
     if(dright(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  enddo
  !
  do i=1,N
     Foo    =  dleft(i) + dright(i) - diag_(i)
     Inv(i) = 1d0/Foo
  end do
end subroutine d_invert_tridiag_matrix

subroutine c_invert_tridiag_matrix(N,sub_,diag_,over_,Inv)
  integer                            :: i,j,N
  complex(8),dimension(N)            :: diag_
  complex(8),dimension(N-1)          :: sub_
  complex(8),dimension(N-1)          :: over_
  complex(8),dimension(N)            :: Inv
  complex(8)                         :: Foo,Cleft,Cright
  complex(8),dimension(N)            :: dleft
  complex(8),dimension(N)            :: dright
  !
  !DOWNWARD:
  dleft(1) = diag_(1)
  if(dleft(1)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  do i=2,N
     foo      = 1d0/dleft(i-1)
     cleft    = sub_(i-1)*foo
     dleft(i) = diag_(i) - cleft*over_(i-1) !over_(i-1)/dleft(i-1)*sub_(i)
     if(dleft(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  enddo
  !
  !UPWARD:
  dright(N) = diag_(N)
  if(dright(N)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  do i=N-1,1,-1
     foo       = 1d0/dright(i+1)
     cright    = over_(i)*foo
     dright(i) = diag_(i) - cright*sub_(i) !sub_(i+1)/dright(i+1)*over_(i)
     if(dright(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
  enddo
  !
  do i=1,N
     Foo    =  dleft(i) + dright(i) - diag_(i)
     Inv(i) = 1d0/Foo
  end do
end subroutine c_invert_tridiag_matrix





!+-----------------------------------------------------------------------------+!
!PURPOSE: return the Nb diagonal NxN blocks of the inverse matrix.
! b = sub-diagonal
! d = main-diagonal
! a = over-diagonal = sub-diagonal (symmetric matrix)
!+-----------------------------------------------------------------------------+!
subroutine d_invert_tridiag_block_matrix(Nb,N,sub_,diag_,over_,Ainv)
  !This subroutine returns the :code:`Nb` diagonal :math:`N\times N` blocks of the inverse of a real block tridiagonal
  !matrix, given by its sub-, main and over-diagonal blocks, without building the full inverse. The complex version is
  !:f:func_inline:`c_invert_tridiag_block_matrix`. Both are instances of the generic interface
  !:f:func_inline:`inv_tridiag`.
  !
  !It uses the block version of the downward and upward recursions of :f:func_inline:`d_invert_tridiag_matrix`,
  !:math:`(A^{-1})_{ii} = [D^L_i + D^R_i - D_i]^{-1}`, each block being inverted with :f:func_inline:`inv`. The blocks
  !are indexed as :code:`sub_(i,:,:) = A(i+1,i)` and :code:`over_(i,:,:) = A(i,i+1)`.
  !
  integer                              :: ib,i,j
  integer                              :: Nb     ! number of blocks
  integer                              :: N      ! size of each block
  real(8),dimension(Nb,N,N)            :: diag_  ! main diagonal blocks, [Nb,N,N]
  real(8),dimension(Nb-1,N,N)          :: sub_   ! sub-diagonal blocks, [Nb-1,N,N]
  real(8),dimension(Nb-1,N,N)          :: over_  ! over-diagonal blocks, [Nb-1,N,N]
  real(8),dimension(Nb,N,N)            :: Ainv   ! diagonal blocks of the inverse matrix, [Nb,N,N]
  real(8),dimension(N,N)               :: Foo,Cleft,Cright
  real(8),dimension(Nb,N,N)            :: Dleft
  real(8),dimension(Nb,N,N)            :: Dright
  !
  !DOWNWARD:
  dleft(1,:,:) = diag_(1,:,:)
  do i=2,Nb
     foo  = dleft(i-1,:,:) ; call inv(foo)
     cleft= matmul(sub_(i-1,:,:),foo)
     dleft(i,:,:) = diag_(i,:,:) - matmul(cleft,over_(i-1,:,:))
  enddo
  !
  !BACKWARD:
  dright(Nb,:,:) = diag_(Nb,:,:)
  do i=Nb-1,1,-1
     foo   = dright(i+1,:,:) ; call inv(foo)
     cright= matmul(over_(i,:,:),foo)
     dright(i,:,:) = diag_(i,:,:) - matmul(cright,sub_(i,:,:))
  enddo
  !
  do ib=1,Nb
     Ainv(ib,:,:)    =  dleft(ib,:,:) + dright(ib,:,:) - diag_(ib,:,:)
     call inv(Ainv(ib,:,:))
  end do
end subroutine d_invert_tridiag_block_matrix

subroutine c_invert_tridiag_block_matrix(Nb,N,sub_,diag_,over_,Ainv)
  integer                        :: ib,i,j,Nb,N
  complex(8),dimension(Nb,N,N)   :: diag_
  complex(8),dimension(Nb-1,N,N) :: sub_
  complex(8),dimension(Nb-1,N,N) :: over_
  complex(8),dimension(Nb,N,N)   :: Ainv
  complex(8),dimension(N,N)      :: Foo,Cleft,Cright
  complex(8),dimension(Nb,N,N)   :: Dleft
  complex(8),dimension(Nb,N,N)   :: Dright
  !
  !DOWNWARD:
  dleft(1,:,:) = diag_(1,:,:)
  do i=2,Nb
     foo  = dleft(i-1,:,:) ; call inv(foo)
     cleft= matmul(sub_(i-1,:,:),foo)
     dleft(i,:,:) = diag_(i,:,:) - matmul(cleft,over_(i-1,:,:))
  enddo
  !
  !BACKWARD:
  dright(Nb,:,:) = diag_(Nb,:,:)
  do i=Nb-1,1,-1
     foo   = dright(i+1,:,:) ; call inv(foo)
     cright= matmul(over_(i,:,:),foo)
     dright(i,:,:) = diag_(i,:,:) - matmul(cright,sub_(i,:,:))
  enddo
  !
  do ib=1,Nb
     Ainv(ib,:,:)    =  dleft(ib,:,:) + dright(ib,:,:) - diag_(ib,:,:)
     call inv(Ainv(ib,:,:))
  end do
end subroutine c_invert_tridiag_block_matrix


!
!
!

subroutine d_invert_tridiag_matrix_mat(Amat)
  !This subroutine works in place on a real tridiagonal matrix :f:var:`Amat`: the three diagonals are extracted with
  !:f:func_inline:`get_tridiag` and the diagonal elements of the inverse are computed as in
  !:f:func_inline:`d_invert_tridiag_matrix`. The complex version is :f:func_inline:`c_invert_tridiag_matrix_mat`. Both
  !are instances of the generic interface :f:func_inline:`inv_tridiag`.
  !
  !On output :f:var:`Amat` is a diagonal matrix, containing the diagonal elements of the inverse of the input matrix.
  !The tridiagonal structure of the input is not checked (see :f:func_inline:`check_tridiag`). The program stops if
  !:f:var:`Amat` is not square, or if a zero pivot is found.
  !
  real(8),dimension(:,:),intent(inout)         :: Amat  ! in: tridiagonal [N,N]; out: diag. matrix of inverse diagonal
  real(8),dimension(size(Amat,1))              :: diag_
  real(8),dimension(size(Amat,1)-1)            :: sub_
  real(8),dimension(size(Amat,1)-1)            :: over_
  real(8),dimension(size(Amat,1))              :: Inv
  integer                                      :: i,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"d_invert_tridiag_matrix_mat","Amat")
  call get_tridiag(Amat,sub_,diag_,over_)
  call inv_tridiag(N,sub_,diag_,over_,Inv)
  Amat= 0d0
  forall(i=1:N)Amat(i,i)=Inv(i)
end subroutine d_invert_tridiag_matrix_mat

subroutine c_invert_tridiag_matrix_mat(Amat)
  complex(8),dimension(:,:),intent(inout)         :: Amat
  complex(8),dimension(size(Amat,1))              :: diag_
  complex(8),dimension(size(Amat,1)-1)            :: sub_
  complex(8),dimension(size(Amat,1)-1)            :: over_
  complex(8),dimension(size(Amat,1))              :: Inv
  integer                                         :: i,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"d_invert_tridiag_matrix_mat","Amat")
  call get_tridiag(Amat,sub_,diag_,over_)
  call inv_tridiag(N,sub_,diag_,over_,Inv)
  Amat= dcmplx(0d0,0d0)
  forall(i=1:N)Amat(i,i)=Inv(i)
end subroutine c_invert_tridiag_matrix_mat

subroutine d_invert_tridiag_block_matrix_mat(Nb,N,Amat)
  !This subroutine works in place on a real block tridiagonal matrix :f:var:`Amat`, made of :code:`Nb` blocks of size
  !:math:`N\times N`: the three block diagonals are extracted with :f:func_inline:`get_tridiag` and the diagonal blocks
  !of the inverse are computed as in :f:func_inline:`d_invert_tridiag_block_matrix`. The complex version is
  !:f:func_inline:`c_invert_tridiag_block_matrix_mat`. Both are instances of the generic interface
  !:f:func_inline:`inv_tridiag`.
  !
  !On output :f:var:`Amat` is a block diagonal matrix, containing the diagonal blocks of the inverse of the input
  !matrix. The block tridiagonal structure of the input is not checked (see :f:func_inline:`check_tridiag`).
  !
  integer,intent(in)                         :: Nb    ! number of blocks
  integer,intent(in)                         :: N     ! size of each block
  real(8),dimension(Nb*N,Nb*N),intent(inout) :: Amat  ! in: block tridiagonal [Nb*N,Nb*N]; out: its inverse diag. blocks
  real(8),dimension(Nb-1,N,N)                :: sub_
  real(8),dimension(Nb,N,N)                  :: diag_
  real(8),dimension(Nb-1,N,N)                :: over_
  real(8),dimension(Nb,N,N)                  :: Inv
  integer                                    :: i,j,is,js,iblock
  call get_tridiag(Nb,N,Amat,sub_,diag_,over_)
  call inv_tridiag(Nb,N,sub_,diag_,over_,Inv)
  Amat=0d0
  do iblock=1,Nb
     do i=1,N
        do j=1,N
           is = i + (iblock-1)*N
           js = j + (iblock-1)*N
           Amat(is,js) = Inv(iblock,i,j)
        enddo
     enddo
  enddo
end subroutine d_invert_tridiag_block_matrix_mat

subroutine c_invert_tridiag_block_matrix_mat(Nb,N,Amat)
  integer,intent(in)                            :: Nb
  integer,intent(in)                            :: N
  complex(8),dimension(Nb*N,Nb*N),intent(inout) :: Amat
  complex(8),dimension(Nb-1,N,N)                :: sub_
  complex(8),dimension(Nb,N,N)                  :: diag_
  complex(8),dimension(Nb-1,N,N)                :: over_
  complex(8),dimension(Nb,N,N)                  :: Inv
  integer                                       :: i,j,is,js,iblock
  call get_tridiag(Nb,N,Amat,sub_,diag_,over_)
  call inv_tridiag(Nb,N,sub_,diag_,over_,Inv)
  Amat= dcmplx(0d0,0d0)
  do iblock=1,Nb
     do i=1,N
        do j=1,N
           is = i + (iblock-1)*N
           js = j + (iblock-1)*N
           Amat(is,js) = Inv(iblock,i,j)
        enddo
     enddo
  enddo
end subroutine c_invert_tridiag_block_matrix_mat
