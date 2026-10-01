!+-----------------------------------------------------------------------------+!
!PURPOSE: Get a the three (sub,main,over) diagonals of a Tridiagonal matrix
!+-----------------------------------------------------------------------------+!
subroutine d_get_tridiag(Amat,sub,diag,over)
  !This subroutine extracts the three diagonals of the real tridiagonal matrix :f:var:`Amat`: the sub-diagonal
  !:code:`sub(i) = A(i+1,i)`, the main diagonal :code:`diag(i) = A(i,i)` and, optionally, the over-diagonal
  !:code:`over(i) = A(i,i+1)`. The complex version is :f:func_inline:`c_get_tridiag`. Both are instances of the generic
  !interface :f:func_inline:`get_tridiag`, the inverse operation being :f:func_inline:`build_tridiag`.
  !
  !The input matrix is not modified and its tridiagonal structure is not checked (see
  !:f:func_inline:`check_tridiag`). The program stops if :f:var:`Amat` is not square.
  !
  real(8),dimension(:,:)                     :: Amat  ! tridiagonal matrix [N,N]
  real(8),dimension(size(Amat,1))            :: diag  ! main diagonal, size N
  real(8),dimension(size(Amat,1)-1)          :: sub   ! sub-diagonal, A(i+1,i), size N-1
  real(8),dimension(size(Amat,1)-1),optional :: over  ! optional out: over-diagonal, A(i,i+1), size N-1
  real(8),dimension(size(Amat,1)-1)          :: over_
  integer                                    :: i,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"d_get_tridiag","Amat")
  forall(i=1:N-1)
     sub(i)  = Amat(i+1,i)
     diag(i) = Amat(i,i)
     over_(i)= Amat(i,i+1)
  end forall
  diag(N) = Amat(N,N)
  if(present(over))over=over_
end subroutine d_get_tridiag
subroutine c_get_tridiag(Amat,sub,diag,over)
  complex(8),dimension(:,:)                     :: Amat
  complex(8),dimension(size(Amat,1))            :: diag
  complex(8),dimension(size(Amat,1)-1)          :: sub
  complex(8),dimension(size(Amat,1)-1),optional :: over
  complex(8),dimension(size(Amat,1)-1)          :: over_
  integer                                       :: i,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"d_get_tridiag","Amat")
  forall(i=1:N-1)
     sub(i)  = Amat(i+1,i)
     diag(i) = Amat(i,i)
     over_(i)= Amat(i,i+1)
  end forall
  diag(N) = Amat(N,N)
  if(present(over))over=over_
end subroutine c_get_tridiag
subroutine d_get_tridiag_block(Nblock,Nsize,Amat,sub,diag,over)
  !This subroutine extracts the three block diagonals of the real block tridiagonal matrix :f:var:`Amat`, made of
  !:code:`Nblock` blocks of size :math:`N_{size}\times N_{size}`: the sub-diagonal blocks :code:`sub(i,:,:) = A(i+1,i)`,
  !the main diagonal blocks :code:`diag(i,:,:) = A(i,i)` and, optionally, the over-diagonal blocks
  !:code:`over(i,:,:) = A(i,i+1)`. The complex version is :f:func_inline:`c_get_tridiag_block`. Both are instances of
  !the generic interface :f:func_inline:`get_tridiag`, the inverse operation being :f:func_inline:`build_tridiag`.
  !
  !The input matrix is not modified and its block tridiagonal structure is not checked.
  !
  integer                                          :: Nblock  ! number of blocks
  integer                                          :: Nsize   ! size of each block
  real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat    ! block tridiagonal matrix [Nblock*Nsize,Nblock*Nsize]
  real(8),dimension(Nblock-1,Nsize,Nsize)          :: sub     ! sub-diagonal blocks, [Nblock-1,Nsize,Nsize]
  real(8),dimension(Nblock,Nsize,Nsize)            :: diag    ! main diagonal blocks, [Nblock,Nsize,Nsize]
  real(8),dimension(Nblock-1,Nsize,Nsize),optional :: over    ! optional out: over-diagonal blocks, [Nblock-1,Nsize,Nsize]
  real(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
  integer                                          :: i,j,iblock,is,js
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Sub(iblock,i,j)   = Amat(Nsize+is,js)
           Diag(iblock,i,j)  = Amat(is,js)
           Over_(iblock,i,j) = Amat(is,Nsize+js)
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Diag(Nblock,i,j) = Amat(is,js)
     enddo
  enddo
  if(present(over))over=over_
end subroutine d_get_tridiag_block
subroutine c_get_tridiag_block(Nblock,Nsize,Amat,sub,diag,over)
  integer                                          :: Nblock
  integer                                          :: Nsize
  complex(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
  complex(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
  complex(8),dimension(Nblock,Nsize,Nsize)            :: diag
  complex(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
  complex(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
  integer                                          :: i,j,iblock,is,js
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Sub(iblock,i,j)   = Amat(Nsize+is,js)
           Diag(iblock,i,j)  = Amat(is,js)
           Over_(iblock,i,j) = Amat(is,Nsize+js)
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Diag(Nblock,i,j) = Amat(is,js)
     enddo
  enddo
  if(present(over))over=over_
end subroutine c_get_tridiag_block
