subroutine splotA1_RR(pname,X,Y1,append)
  integer                       :: i,Np
  character(len=*)              :: pname   !name of the output file
  real(8),dimension(:)          :: X
  real(8),dimension(size(X))    :: Y1      !function values Y1(...,k) at X(k), real or complex, rank 1 to 7
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  Np=size(X)
  do i=1,Np
     write(unit,*)X(i),Y1(i)
  enddo
  close(unit)
end subroutine splotA1_RR
subroutine splotA1_RC(pname,X,Y1,append)
  integer                       :: i,j,Np
  character(len=*)              :: pname   !name of the output file
  real(8),dimension(:)          :: X
  complex(8),dimension(size(X)) :: Y1      !function values Y1(...,k) at X(k), real or complex, rank 1 to 7
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  Np=size(X)
  do i=1,Np
     write(unit,*)X(i),dimag(Y1(i)),dreal(Y1(i))
  enddo
  close(unit)
end subroutine splotA1_RC
!


!----------------------------
!----------------------------
!----------------------------


subroutine splotA2_RR(pname,X,Y1,append)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname   !name of the output file
  real(8),dimension(:,:)        :: Y1
  real(8),dimension(size(Y1,2)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  !
  do i=1,Ny1
     do j=1,Ny2
        write(unit,*)X(j),Y1(i,j)
     enddo
     write(unit,*)""
  enddo
  close(unit)
end subroutine splotA2_RR
!
subroutine splotA2_RC(pname,X,Y1,append)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname   !name of the output file
  complex(8),dimension(:,:)     :: Y1
  real(8),dimension(size(Y1,2)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  !
  do i=1,Ny1
     do j=1,Ny2
        write(unit,*)X(j),dimag(Y1(i,j)),dreal(Y1(i,j))
     enddo
     write(unit,*)""
  enddo
  close(unit)
end subroutine splotA2_RC

!----------------------------
!----------------------------
!----------------------------

subroutine splotA3_RR(pname,X,Y1,append)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname   !name of the output file
  real(8),dimension(:,:,:)      :: Y1
  real(8),dimension(size(Y1,3)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           write(unit,*)X(k),Y1(i,j,k)
        enddo
        write(unit,*)""
     enddo
  enddo
  close(unit)
end subroutine splotA3_RR
!
subroutine splotA3_RC(pname,X,Y1,append)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname   !name of the output file
  complex(8),dimension(:,:,:)   :: Y1
  real(8),dimension(size(Y1,3)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           write(unit,*)X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k))
        enddo
        write(unit,*)""
     enddo
  enddo
  close(unit)
end subroutine splotA3_RC

!----------------------------
!----------------------------
!----------------------------

subroutine splotA4_RR(pname,X,Y1,append)
  integer                       :: Ny1,Ny2,Ny3,Ny4
  integer                       :: i1,i2,i3,i4
  character(len=*)              :: pname   !name of the output file
  real(8),dimension(:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,4)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
              write(unit,*)X(i4),Y1(i1,i2,i3,i4)
           enddo
           write(unit,*)""
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA4_RR
!
subroutine splotA4_RC(pname,X,Y1,append)
  integer                       :: Ny1,Ny2,Ny3,Ny4
  integer                       :: i1,i2,i3,i4
  character(len=*)              :: pname   !name of the output file
  complex(8),dimension(:,:,:,:) :: Y1
  real(8),dimension(size(Y1,4)) :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional              :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                       :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
              write(unit,*)X(i4),dimag(Y1(i1,i2,i3,i4)),dreal(Y1(i1,i2,i3,i4))
           enddo
           write(unit,*)""
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA4_RC

!----------------------------
!----------------------------
!----------------------------



subroutine splotA5_RR(pname,X,Y1,append)
  integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                         :: i1,i2,i3,i4,i5
  character(len=*)                :: pname   !name of the output file
  real(8),dimension(:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,5))   :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                         :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
                 write(unit,*)X(i5),Y1(i1,i2,i3,i4,i5)
              enddo
              write(unit,*)""
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA5_RR
!
subroutine splotA5_RC(pname,X,Y1,append)
  integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                         :: i1,i2,i3,i4,i5
  character(len=*)                :: pname   !name of the output file
  complex(8),dimension(:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,5))   :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                         :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
                 write(unit,*)X(i5),dimag(Y1(i1,i2,i3,i4,i5)),dreal(Y1(i1,i2,i3,i4,i5))
              enddo
              write(unit,*)""
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA5_RC


!----------------------------
!----------------------------
!----------------------------


subroutine splotA6_RR(pname,X,Y1,append)
  integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                           :: i1,i2,i3,i4,i5,i6
  character(len=*)                  :: pname   !name of the output file
  real(8),dimension(:,:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,6))     :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                  :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                           :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
                    write(unit,*)X(i6),Y1(i1,i2,i3,i4,i5,i6)
                 enddo
                 write(unit,*)""
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA6_RR
!
subroutine splotA6_RC(pname,X,Y1,append)
  integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                           :: i1,i2,i3,i4,i5,i6
  character(len=*)                  :: pname   !name of the output file
  complex(8),dimension(:,:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,6))     :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                  :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                           :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
  !
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
                    write(unit,*)X(i6),dimag(Y1(i1,i2,i3,i4,i5,i6)),dreal(Y1(i1,i2,i3,i4,i5,i6))
                 enddo
                 write(unit,*)""
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA6_RC

!----------------------------
!----------------------------
!----------------------------



subroutine splotA7_RR(pname,X,Y1,append)
  integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                    :: pname   !name of the output file
  real(8),dimension(:,:,:,:,:,:,:)    :: Y1
  real(8),dimension(size(Y1,7))       :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                    :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                             :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
                       write(unit,*)X(i7),Y1(i1,i2,i3,i4,i5,i6,i7)
                    enddo
                    write(unit,*)""
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA7_RR
!
subroutine splotA7_RC(pname,X,Y1,append)
  integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                    :: pname   !name of the output file
  complex(8),dimension(:,:,:,:,:,:,:) :: Y1
  real(8),dimension(size(Y1,7))       :: X       !abscissa, one per element of the last dimension of Y1
  logical,optional                    :: append  !if T add to an existing file instead of overwriting it (default F)
  logical                             :: check,append_
  append_=.false.;if(present(append))append_=append
  if(append_)then
     inquire(file=reg(pname),exist=check)
     open(free_unit(unit),file=reg(pname),position="append")
     if(check)write(unit,*)
  else
     open(free_unit(unit),file=reg(pname))
  endif
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
                       write(unit,*)X(i7),dimag(Y1(i1,i2,i3,i4,i5,i6,i7)),dreal(Y1(i1,i2,i3,i4,i5,i6,i7))
                    enddo
                    write(unit,*)""
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
end subroutine splotA7_RC
