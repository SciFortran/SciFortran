subroutine sreadA1_RR(pname,X,Y1)
  integer                             :: i,Np
  character(len=*)                    :: pname  !name of the input file
  real(8),dimension(:)                :: X
  real(8),dimension(size(X))          :: Y1
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  Np=size(X)
  do i=1,Np
     read(unit,*)X(i),Y1(i)
  enddo
  close(unit)
end subroutine sreadA1_RR

subroutine sreadA1_RC(pname,X,Y1)
  integer                       :: i,Np
  character(len=*)              :: pname  !name of the input file
  real(8),dimension(:)          :: X
  complex(8),dimension(size(X)) :: Y1
  real(8),dimension(size(X))    :: reY,imY
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  Np=size(X)
  do i=1,Np
     read(unit,*)X(i),imY(i),reY(i)
  enddo
  Y1=dcmplx(reY,imY)
  close(unit)
end subroutine sreadA1_RC



!----------------------------
!----------------------------
!----------------------------



subroutine sreadA2_RR(pname,X,Y1)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname  !name of the input file
  real(8),dimension(:,:)        :: Y1
  real(8),dimension(size(Y1,2)) :: X
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  !
  do i=1,Ny1
     do j=1,Ny2
        read(unit,*)X(j),Y1(i,j)
     enddo
  enddo
  !
  close(unit)
end subroutine sreadA2_RR

subroutine sreadA2_RC(pname,X,Y1)
  integer                                  :: i,j,Ny1,Ny2
  character(len=*)                         :: pname  !name of the input file
  complex(8),dimension(:,:)                :: Y1
  real(8),dimension(size(Y1,2))            :: X
  real(8),dimension(size(Y1,1),size(Y1,2)) :: reY,imY
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  do i=1,Ny1
     do j=1,Ny2
        read(unit,*)X(j),imY(i,j),reY(i,j)
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA2_RC


!----------------------------
!----------------------------
!----------------------------




subroutine sreadA3_RR(pname,X,Y1)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname  !name of the input file
  real(8),dimension(:,:,:)      :: Y1
  real(8),dimension(size(Y1,3)) :: X
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           read(unit,*)X(k),Y1(i,j,k)
        enddo
     enddo
  enddo
  close(unit)
end subroutine sreadA3_RR

subroutine sreadA3_RC(pname,X,Y1)
  integer                                             :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                    :: pname  !name of the input file
  complex(8),dimension(:,:,:)                         :: Y1
  real(8),dimension(size(Y1,3))                       :: X
  real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)) :: reY,imY
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           read(unit,*)X(k),imY(i,j,k),reY(i,j,k)
        enddo
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA3_RC


!----------------------------
!----------------------------
!----------------------------


subroutine sreadA4_RR(pname,X,Y1)
  integer                       :: Ny1,Ny2,Ny3,Ny4
  integer                       :: i1,i2,i3,i4
  character(len=*)              :: pname  !name of the input file
  real(8),dimension(:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,4)) :: X
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              read(unit,*)X(i4),Y1(i1,i2,i3,i4)
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine sreadA4_RR

subroutine sreadA4_RC(pname,X,Y1)
  integer                       :: Ny1,Ny2,Ny3,Ny4
  integer                       :: i1,i2,i3,i4
  character(len=*)              :: pname  !name of the input file
  complex(8),dimension(:,:,:,:) :: Y1
  real(8),dimension(size(Y1,4)) :: X
  real(8),dimension(&
       size(Y1,1),&
       size(Y1,2),&
       size(Y1,3),&
       size(Y1,4))              :: reY,imY
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              read(unit,*)X(i4),imY(i1,i2,i3,i4),reY(i1,i2,i3,i4)
           enddo
        enddo
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA4_RC


!----------------------------
!----------------------------
!----------------------------



subroutine sreadA5_RR(pname,X,Y1)
  integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                         :: i1,i2,i3,i4,i5
  character(len=*)                :: pname  !name of the input file
  real(8),dimension(:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,5))   :: X
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 read(unit,*)X(i5),Y1(i1,i2,i3,i4,i5)
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine sreadA5_RR

subroutine sreadA5_RC(pname,X,Y1)
  integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                         :: i1,i2,i3,i4,i5
  character(len=*)                :: pname  !name of the input file
  complex(8),dimension(:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,5))   :: X
  real(8),dimension(&
       size(Y1,1),&
       size(Y1,2),&
       size(Y1,3),&
       size(Y1,4),&
       size(Y1,5))                :: reY,imY
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 read(unit,*)X(i5),imY(i1,i2,i3,i4,i5),reY(i1,i2,i3,i4,i5)
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA5_RC


!----------------------------
!----------------------------
!----------------------------


subroutine sreadA6_RR(pname,X,Y1)
  integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                           :: i1,i2,i3,i4,i5,i6
  character(len=*)                  :: pname  !name of the input file
  real(8),dimension(:,:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,6))     :: X
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    read(unit,*)X(i6),Y1(i1,i2,i3,i4,i5,i6)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine sreadA6_RR

subroutine sreadA6_RC(pname,X,Y1)
  integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                           :: i1,i2,i3,i4,i5,i6
  character(len=*)                  :: pname  !name of the input file
  complex(8),dimension(:,:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,6))     :: X
  real(8),dimension(&
       size(Y1,1),&
       size(Y1,2),&
       size(Y1,3),&
       size(Y1,4),&
       size(Y1,5),&
       size(Y1,6))                  :: reY,imY
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    read(unit,*)X(i6),imY(i1,i2,i3,i4,i5,i6),reY(i1,i2,i3,i4,i5,i6)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA6_RC


!----------------------------
!----------------------------
!----------------------------



subroutine sreadA7_RR(pname,X,Y1)
  integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                    :: pname  !name of the input file
  real(8),dimension(:,:,:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,7))       :: X
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  Ny7=size(Y1,7)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    do i7=1,Ny7
                       read(unit,*)X(i7),Y1(i1,i2,i3,i4,i5,i6,i7)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine sreadA7_RR

subroutine sreadA7_RC(pname,X,Y1)
  integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                    :: pname  !name of the input file
  complex(8),dimension(:,:,:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,7))       :: X
  real(8),dimension(&
       size(Y1,1),&
       size(Y1,2),&
       size(Y1,3),&
       size(Y1,4),&
       size(Y1,5),&
       size(Y1,6),&
       size(Y1,7))                    :: reY,imY
  !
  call ioread_control(pname,control)
  open(free_unit(unit),file=reg(pname))
  !
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  Ny7=size(Y1,7)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    do i7=1,Ny7
                       read(unit,*)X(i7),imY(i1,i2,i3,i4,i5,i6,i7),reY(i1,i2,i3,i4,i5,i6,i7)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  Y1=dcmplx(reY,imY)
end subroutine sreadA7_RC
