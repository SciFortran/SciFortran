!+-----------------------------------------------------------------------------+!
!PURPOSE: Build a tridiagonal Matrix Amat from the three (sub,main,over) diagonal.
! In this version the over-diagonal is optional
!+-----------------------------------------------------------------------------+!
function d_build_tridiag(sub,diag,over) result(Amat)
  !This function builds the real tridiagonal matrix :f:var:`Amat` from its three diagonals: the sub-diagonal
  !:code:`A(i+1,i) = sub(i)`, the main diagonal :code:`A(i,i) = diag(i)` and the over-diagonal
  !:code:`A(i,i+1) = over(i)`. If :f:var:`over` is not present the matrix is symmetric, :code:`over = sub`. The complex
  !version is :f:func_inline:`c_build_tridiag`. Both are instances of the generic interface
  !:f:func_inline:`build_tridiag`, the inverse operation being :f:func_inline:`get_tridiag`.
  !
  real(8),dimension(:)                     :: diag  ! main diagonal, size N
  real(8),dimension(size(diag)-1)          :: sub   ! sub-diagonal, size N-1
  real(8),dimension(size(diag)-1),optional :: over  ! optional: over-diagonal, size N-1 (default: over = sub)
  real(8),dimension(size(diag),size(diag)) :: Amat  ! tridiagonal matrix, [N,N]
  real(8),dimension(size(diag)-1)          :: over_
  integer                                  :: i,N
  over_=sub;if(present(over))over_=over
  N=size(diag)
  Amat=0d0
  forall(i=1:N-1)
     Amat(i+1,i) = sub(i)
     Amat(i,i)   = diag(i)
     Amat(i,i+1) = over_(i)
  end forall
  Amat(N,N)=diag(N)
end function d_build_tridiag
function c_build_tridiag(sub,diag,over) result(Amat)
  complex(8),dimension(:)                     :: diag
  complex(8),dimension(size(diag)-1)          :: sub
  complex(8),dimension(size(diag)-1),optional :: over
  complex(8),dimension(size(diag),size(diag)) :: Amat
  complex(8),dimension(size(diag)-1)          :: over_
  integer                                     :: i,N
  over_=sub;if(present(over))over_=over
  N=size(diag)
  Amat=dcmplx(0d0,0d0)
  forall(i=1:N-1)
     Amat(i+1,i) = sub(i)
     Amat(i,i)   = diag(i)
     Amat(i,i+1) = over_(i)
  end forall
  Amat(N,N)=diag(N)
end function c_build_tridiag
function d_build_tridiag_block(Nblock,Nsize,sub,diag,over) result(Amat)
  !This function builds the real block tridiagonal matrix :f:var:`Amat`, made of :code:`Nblock` blocks of size
  !:math:`N_{size}\times N_{size}`, from its three block diagonals: the sub-diagonal blocks
  !:code:`A(i+1,i) = sub(i,:,:)`, the main diagonal blocks :code:`A(i,i) = diag(i,:,:)` and the over-diagonal blocks
  !:code:`A(i,i+1) = over(i,:,:)`. If :f:var:`over` is not present :code:`over = sub`. The complex version is
  !:f:func_inline:`c_build_tridiag_block`. Both are instances of the generic interface
  !:f:func_inline:`build_tridiag`, the inverse operation being :f:func_inline:`get_tridiag`.
  !
  integer                                          :: Nblock  ! number of blocks
  integer                                          :: Nsize   ! size of each block
  real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat    ! block tridiagonal matrix, [Nblock*Nsize,Nblock*Nsize]
  real(8),dimension(Nblock-1,Nsize,Nsize)          :: sub     ! sub-diagonal blocks, [Nblock-1,Nsize,Nsize]
  real(8),dimension(Nblock,Nsize,Nsize)            :: diag    ! main diagonal blocks, [Nblock,Nsize,Nsize]
  real(8),dimension(Nblock-1,Nsize,Nsize),optional :: over    ! optional: over-diagonal blocks (default: over = sub)
  real(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
  integer                                          :: i,j,iblock,is,js
  over_=sub;if(present(over))over_=over
  !
  Amat=0d0
  !
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Amat(Nsize+is,js) = Sub(iblock,i,j)
           Amat(is,js)       = Diag(iblock,i,j)
           Amat(is,Nsize+js) = Over_(iblock,i,j)
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Amat(is,js)       = Diag(Nblock,i,j)
     enddo
  enddo
end function d_build_tridiag_block
function c_build_tridiag_block(Nblock,Nsize,sub,diag,over) result(Amat)
  integer                                          :: Nblock
  integer                                          :: Nsize
  complex(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
  complex(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
  complex(8),dimension(Nblock,Nsize,Nsize)            :: diag
  complex(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
  complex(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
  integer                                          :: i,j,iblock,is,js
  over_=sub;if(present(over))over_=over
  !
  Amat=0d0
  !
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Amat(Nsize+is,js) = Sub(iblock,i,j)
           Amat(is,js)       = Diag(iblock,i,j)
           Amat(is,Nsize+js) = Over_(iblock,i,j)
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Amat(is,js)       = Diag(Nblock,i,j)
     enddo
  enddo
end function c_build_tridiag_block
