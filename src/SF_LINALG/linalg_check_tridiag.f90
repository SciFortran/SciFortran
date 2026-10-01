!+-----------------------------------------------------------------------------+!
!PURPOSE:  comment
!+-----------------------------------------------------------------------------+!
function d_check_tridiag(Amat) result(Mcheck)
  !This function checks whether the real matrix :f:var:`Amat` is tridiagonal, i.e. whether all the elements outside
  !the main, sub- and over-diagonals are exactly zero. The complex version is :f:func_inline:`c_check_tridiag`. Both
  !are instances of the generic interface :f:func_inline:`check_tridiag`. The result is :code:`.true.` for a tridiagonal
  !matrix. The program stops if :f:var:`Amat` is not square.
  !
  real(8),dimension(:,:)                       :: Amat    ! square matrix [N,N] to be checked
  logical,dimension(size(Amat,1),size(Amat,2)) :: Lmat
  logical                                      :: Mcheck  ! .true. if Amat is tridiagonal
  integer                                      :: i,j,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"d_check_tridiag","Amat")
  Lmat=.true.
  forall(i=1:N-1)
     Lmat(i+1,i)=.false.
     Lmat(i,i)  =.false.
     Lmat(i,i+1)=.false.
  end forall
  Lmat(N,N)=.false.
  Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
end function d_check_tridiag
function c_check_tridiag(Amat) result(Mcheck)
  complex(8),dimension(:,:)                    :: Amat
  logical,dimension(size(Amat,1),size(Amat,2)) :: Lmat
  logical                                      :: Mcheck
  integer                                      :: i,N
  N=size(Amat,1)
  call assert_shape(Amat,[N,N],"c_check_tridiag","Amat")
  Lmat=.true.
  forall(i=1:N-1)
     Lmat(i+1,i)=.false.
     Lmat(i,i)  =.false.
     Lmat(i,i+1)=.false.
  end forall
  Lmat(N,N)=.false.
  Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
end function c_check_tridiag
function d_check_tridiag_block(Nblock,Nsize,Amat) result(Mcheck)
  !This function checks whether the real matrix :f:var:`Amat`, made of :code:`Nblock` blocks of size
  !:math:`N_{size}\times N_{size}`, is block tridiagonal, i.e. whether all the elements outside the main, sub- and
  !over-diagonal blocks are exactly zero. The complex version is :f:func_inline:`c_check_tridiag_block`. Both are
  !instances of the generic interface :f:func_inline:`check_tridiag`. The result is :code:`.true.` for a block
  !tridiagonal matrix.
  !
  integer                                          :: Nblock  ! number of blocks
  integer                                          :: Nsize   ! size of each block
  real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat    ! matrix [Nblock*Nsize,Nblock*Nsize] to be checked
  logical,dimension(Nblock*Nsize,Nblock*Nsize)     :: Lmat
  integer                                          :: i,j,iblock,is,js
  logical                                          :: Mcheck  ! .true. if Amat is block tridiagonal
  Lmat=.true.
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Lmat(Nsize+is,js) =.false.
           Lmat(is,js)       =.false.
           Lmat(is,Nsize+js) =.false.
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Lmat(is,js)=.false.
     enddo
  enddo
  Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
end function d_check_tridiag_block
function c_check_tridiag_block(Nblock,Nsize,Amat) result(Mcheck)
  integer                                         :: Nblock
  integer                                         :: Nsize
  complex(8),dimension(Nblock*Nsize,Nblock*Nsize) :: Amat
  logical,dimension(Nblock*Nsize,Nblock*Nsize)    :: Lmat
  integer                                         :: i,j,iblock,is,js
  logical                                         :: Mcheck
  Lmat=.true.
  do iblock=1,Nblock-1
     do i=1,Nsize
        do j=1,Nsize
           is = i + (iblock-1)*Nsize
           js = j + (iblock-1)*Nsize
           Lmat(Nsize+is,js) =.false.
           Lmat(is,js)       =.false.
           Lmat(is,Nsize+js) =.false.
        enddo
     enddo
  enddo
  do i=1,Nsize
     do j=1,Nsize
        is = i + (Nblock-1)*Nsize
        js = j + (Nblock-1)*Nsize
        Lmat(is,js)=.false.
     enddo
  enddo
  Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
end function c_check_tridiag_block
