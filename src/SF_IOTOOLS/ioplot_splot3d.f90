subroutine d_splot3D(pname,X1,X2,Y,xmin,xmax,ymin,ymax,nosurface,wlines,nlines)
  integer                              :: i,j,Nx1,Nx2,count,Nl
  character(len=*)                     :: pname                !name of the data file, base name of the gnuplot scripts
  real(8),dimension(:)                 :: X1                   !grid points along x, first dimension of Y
  real(8),dimension(:)                 :: X2                   !grid points along y, second dimension of Y
  real(8),dimension(size(X1),size(X2)) :: Y                    !values Y(i,j)=F(X1(i),X2(j)); Y(i,j,m) for animate
  real(8),optional                     :: xmin,xmax,ymin,ymax  !plot ranges in x and y (default: extrema of X1, X2)
  integer,optional                     :: nlines               !step of the first index in the _withlines file (default 5)
  logical,optional                     :: wlines               !if present also write the _withlines file (any value)
  logical,optional                     :: nosurface            !if T do not write the gnuplot script of the surface (default F)
  logical                              :: nosurface_
  real(8)                              :: X1min,X1max
  real(8)                              :: X2min,X2max
  character(len=12)                    :: minx,miny,maxx,maxy
  character(len=:),allocatable         :: fname,dname
  fname=get_filename(reg(pname))
  dname=get_filepath(reg(pname))
  nosurface_=.false.;if(present(nosurface))nosurface_=nosurface
  Nx1=size(X1) ; Nx2=size(X2)
  Nl=5; if(present(nlines))Nl=nlines
  open(719,file=str(pname))
  if(present(wlines))open(720,file=str(pname)//"_withlines")
  do i=1,Nx1
     count=mod(i,Nl)
     do j=1,Nx2
        write(719,*)X1(i),X2(j),Y(i,j)
        if(present(wlines).AND.count==0)write(720,*)X1(i),X2(j),Y(i,j)
     enddo
     write(719,*)""
     if(present(wlines))write(720,*)""
     if(present(wlines))write(720,*)""
  enddo
  close(719)
  X1min=minval(X1);if(present(xmin))X1min=xmin
  X1max=maxval(X1);if(present(xmax))X1max=xmax
  X2min=minval(X2);if(present(ymin))X2min=ymin
  X2max=maxval(X2);if(present(ymax))X2max=ymax
  write(minx,"(f12.4)")X1min
  write(maxx,"(f12.4)")X1max
  write(miny,"(f12.4)")X2min
  write(maxy,"(f12.4)")X2max
  if(present(wlines))close(720)

  open(10,file=reg(pname)//"_map.gp")

  write(10,*)"set title '"//fname//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//minx//":"//maxx//"]"
  write(10,*)"set yrange ["//miny//":"//maxy//"]"
  write(10,*)"splot '"//trim(fname)//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out '"//fname//".png'"
  write(10,*)"#rep"
  close(10)
  !
  if(nosurface_)return
  !
  open(10,file=reg(pname)//"_surface.gp")
  write(10,*)"set title '"//fname//"'"
  write(10,*)"unset key"
  write(10,*)"unset grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot '"//fname//"' with pm3d"
  if(present(wlines))write(10,*)"rep '"//reg(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out '"//fname//".png'"
  write(10,*)"#rep"
  close(10)
end subroutine d_splot3D


subroutine c_splot3D(pname,X1,X2,Y,xmin,xmax,ymin,ymax,nosurface,wlines,nlines)
  integer                                 :: i,j,Nx1,Nx2,count,Nl
  character(len=*)                        :: pname                !name of the data file, base name of the gnuplot scripts
  real(8),dimension(:)                    :: X1                   !grid points along x, first dimension of Y
  real(8),dimension(:)                    :: X2                   !grid points along y, second dimension of Y
  complex(8),dimension(size(X1),size(X2)) :: Y                    !values Y(i,j)=F(X1(i),X2(j)); Y(i,j,m) for animate
  real(8),optional                        :: xmin,xmax,ymin,ymax  !plot ranges in x and y (default: extrema of X1, X2)
  integer,optional                        :: nlines               !step of the first index in the _withlines file (default 5)
  logical,optional                        :: wlines               !if present also write the _withlines file (any value)
  logical,optional                        :: nosurface            !if T do not write the gnuplot script of the surface (default F)
  logical                                 :: nosurface_
  real(8)                                 :: X1min,X1max
  real(8)                                 :: X2min,X2max
  character(len=12)                    :: minx,miny,maxx,maxy
  character(len=:),allocatable            :: fname,dname
  fname=get_filename(reg(pname))
  dname=get_filepath(reg(pname))
  nosurface_=.false.;if(present(nosurface))nosurface_=nosurface
  Nx1=size(X1) ; Nx2=size(X2)
  Nl=5; if(present(nlines))Nl=nlines
  open(619,file=dname//"re_"//fname)
  open(719,file=dname//"im_"//fname)
  if(present(wlines))then
     open(620,file=dname//"re_"//fname//"_withlines")
     open(720,file=dname//"im_"//fname//"_withlines")
  endif
  do i=1,Nx1
     count=mod(i,Nl)
     do j=1,Nx2
        write(619,*)X1(i),X2(j),dreal(Y(i,j))
        write(719,*)X1(i),X2(j),dimag(Y(i,j))
        if(present(wlines).AND.count==0)write(620,*)X1(i),X2(j),dreal(Y(i,j))
        if(present(wlines).AND.count==0)write(720,*)X1(i),X2(j),dimag(Y(i,j))
     enddo
     write(619,*)""
     write(719,*)""
     if(present(wlines))then
        write(620,*)""
        write(620,*)""
        write(720,*)""
        write(720,*)""
     endif
  enddo
  close(619)
  close(719)
  X1min=minval(X1);if(present(xmin))X1min=xmin
  X1max=maxval(X1);if(present(xmax))X1max=xmax
  X2min=minval(X2);if(present(ymin))X2min=ymin
  X2max=maxval(X2);if(present(ymax))X2max=ymax
  write(minx,"(f12.4)")X1min
  write(maxx,"(f12.4)")X1max
  write(miny,"(f12.4)")X2min
  write(maxy,"(f12.4)")X2max
  if(present(wlines))then
     close(620)
     close(720)
  endif


  !Re:
  open(10,file=pname//"_re_map.gp")
  write(10,*)"set title 'Re_"//fname//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//minx//":"//maxx//"]"
  write(10,*)"set yrange ["//miny//":"//maxy//"]"
  write(10,*)"splot 're_"//trim(fname)//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 're_"//fname//".png'"
  write(10,*)"#rep"
  close(10)
  !Im
  open(10,file=pname//"_im_map.gp")

  write(10,*)"set title 'Im_"//fname//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//minx//":"//maxx//"]"
  write(10,*)"set yrange ["//miny//":"//maxy//"]"
  write(10,*)"splot 'im_"//fname//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 'im_"//fname//".png'"
  write(10,*)"#rep"
  close(10)
  !
  if(nosurface_)return
  !Re
  open(10,file=pname//"_re_surface.gp")

  write(10,*)"set title 'Re_"//(fname)//"'"
  write(10,*)"unset key"
  write(10,*)"unset grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot 're_"//(fname)//"' with pm3d"
  if(present(wlines))write(10,*)"rep 're_"//str(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 're_"//(fname)//".png'"
  write(10,*)"#rep"
  close(10)
  !Im
  open(10,file=pname//"_im_surface.gp")

  write(10,*)"set title 'Im_"//(fname)//"'"
  write(10,*)"set nokey"
  write(10,*)"set grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot 'im_"//(fname)//"' with pm3d"
  if(present(wlines))write(10,*)"rep 'im_"//str(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 'im_"//(fname)//".png'"
  write(10,*)"#rep"
  close(10)
end subroutine c_splot3D



subroutine d_splot3d_animate(pname,X1,X2,Y,xmin,xmax,ymin,ymax)
  integer                  :: i,j,m,Nx1,Nx2,Nt
  character(len=*)         :: pname                !name of the data file, base name of the gnuplot scripts
  real(8),dimension(:)     :: X1                   !grid points along x, first dimension of Y
  real(8),dimension(:)     :: X2                   !grid points along y, second dimension of Y
  real(8),dimension(:,:,:) :: Y                    !values Y(i,j)=F(X1(i),X2(j)); Y(i,j,m) for animate
  real(8)                  :: X1min,X1max
  real(8)                  :: X2min,X2max
  real(8),optional         :: xmin,xmax,ymin,ymax  !plot ranges in x and y (default: extrema of X1, X2)
  real(8)                  :: Rmin(3),Rmax(3),Zmin,Zmax
  Character(len=256)       :: fname,dname
  fname=get_filename(reg(pname))
  dname=get_filepath(reg(pname))
  Nx1=size(X1) ; Nx2=size(X2)
  if(size(Y,1)/=Nx1) stop "Error Nx1"
  if(size(Y,2)/=Nx2) stop "Error Nx2"
  Nt=size(Y,3)
  write(*,*) "splot3d_animate "//str(fname)//" ("//str(Nx1)//","//str(Nx2)//","//str(Nt)//")"
  open(719,file=reg(pname))
  do m=1,Nt
     do i=1,Nx1
        do j=1,Nx2
           write(719,*)X1(i),X2(j),Y(i,j,m)
        enddo
        write(719,*)""
     enddo
  enddo
  close(719)
  X1min=minval(X1);if(present(xmin))X1min=xmin
  X1max=maxval(X1);if(present(xmax))X1max=xmax
  X2min=minval(X2);if(present(ymin))X2min=ymin
  X2max=maxval(X2);if(present(ymax))X2max=ymax
  Rmin=minval(Y);Rmax=maxval(Y)
  Zmin=minval(Rmin);Zmax=maxval(Rmax)
  open(10,file=reg(pname)//"_map.gp")
  ! write(10,*)"gnuplot -persist << EOF"
  write(10,*)"reset"
  write(10,*)"#set term gif animate"
  write(10,*)"#set output '"//trim(fname)//".gif'"

  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
  write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
  write(10,*)"n="//trim(txtfy(Nt))
  write(10,*)"do for [i=0:n-1]{"
  write(10,*)"set title sprintf('"//trim(fname)//"; i=%i',i+1)"
  write(10,*)"splot '"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
       "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
  write(10,*)"}"
  write(10,*)"#set output"
  ! write(10,"(A)")"EOF"
  close(10)
  ! call system("chmod +x "//reg(pname)//"_map.gp")
end subroutine d_splot3D_animate



subroutine c_splot3d_animate(pname,X1,X2,Y,xmin,xmax,ymin,ymax)
  integer                     :: i,j,m,Nx1,Nx2,Nt
  character(len=*)            :: pname                !name of the data file, base name of the gnuplot scripts
  real(8),dimension(:)        :: X1                   !grid points along x, first dimension of Y
  real(8),dimension(:)        :: X2                   !grid points along y, second dimension of Y
  complex(8),dimension(:,:,:) :: Y                    !values Y(i,j)=F(X1(i),X2(j)); Y(i,j,m) for animate
  real(8)                     :: X1min,X1max
  real(8)                     :: X2min,X2max
  real(8),optional            :: xmin,xmax,ymin,ymax  !plot ranges in x and y (default: extrema of X1, X2)
  real(8)                     :: Rmin(3),Rmax(3),Zmin,Zmax  
  character(len=256)          :: fname,dname
  fname=get_filename(reg(pname))
  dname=get_filepath(reg(pname))
  Nx1=size(X1) ; Nx2=size(X2)
  if(size(Y,1)/=Nx1) stop "Error Nx1"
  if(size(Y,2)/=Nx2) stop "Error Nx2"
  Nt=size(Y,3)
  write(*,*) "splot3d_animate "//str(fname)//" ("//str(Nx1)//","//str(Nx2)//","//str(Nt)//")"
  open(619,file=adjustl(trim(dname))//"re_"//adjustl(trim(fname)))
  open(719,file=adjustl(trim(dname))//"im_"//adjustl(trim(fname)))
  do m=1,Nt
     do i=1,Nx1
        do j=1,Nx2
           write(619,*)X1(i),X2(j),real(Y(i,j,m),8)
           write(719,*)X1(i),X2(j),aimag(Y(i,j,m))
        enddo
        write(619,*)""
        write(719,*)""
     enddo
  enddo
  close(619)
  close(719)
  X1min=minval(X1);if(present(xmin))X1min=xmin
  X1max=maxval(X1);if(present(xmax))X1max=xmax
  X2min=minval(X2);if(present(ymin))X2min=ymin
  X2max=maxval(X2);if(present(ymax))X2max=ymax
  Rmin=minval(dreal(Y));Rmax=maxval(dreal(Y))
  Zmin=minval(Rmin);Zmax=maxval(Rmax)
  open(10,file=reg(pname)//"_re_map.gp")
  ! write(10,*)"gnuplot -persist << EOF"
  write(10,*)"reset"
  write(10,*)"#set term gif animate"
  write(10,*)"#set output 're_"//trim(fname)//".gif'"

  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
  write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
  write(10,*)"n="//trim(txtfy(Nt))
  write(10,*)"do for [i=0:n-1]{"
  write(10,*)"set title sprintf('Re("//trim(fname)//"); i=%i',i+1)"
  write(10,*)"splot 're_"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
       "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
  write(10,*)"}"
  write(10,*)"#set output"
  ! write(10,"(A)")"EOF"
  close(10)

  Rmin=minval(dimag(Y));Rmax=maxval(dimag(Y))
  Zmin=minval(Rmin);Zmax=maxval(Rmax)
  open(10,file=reg(pname)//"_im_map.gp")
  ! write(10,*)"gnuplot -persist << EOF"
  write(10,*)"reset"
  write(10,*)"#set term gif animate"
  write(10,*)"#set output 'im_"//trim(fname)//".gif'"

  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
  write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
  write(10,*)"n="//trim(txtfy(Nt))
  write(10,*)"do for [i=0:n-1]{"
  write(10,*)"set title sprintf('Im("//trim(fname)//"); i=%i',i+1)"
  write(10,*)"splot 'im_"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
       "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
  write(10,*)"}"
  write(10,*)"#set output"
  ! write(10,"(A)")"EOF"
  close(10)
  ! call system("chmod +x "//reg(pname)//".gp")
end subroutine c_splot3d_animate


! subroutine splot3D__(pname,X1,X2,Y,wlines,nlines)
!   integer                              :: i,j,Nx1,count,Nl,Nk
!   character(len=*)                     :: pname
!   real(8),dimension(:)                 :: X1 !(0:Nt)
!   real(8),dimension(:,:)               :: X2 !(0:Nt,Lk)
!   real(8),dimension(size(X2,1),size(X2,2)) :: Y  !(0:Nt,Lk)
!   integer,optional                     :: nlines
!   logical,optional                     :: wlines
!   real(8)                              :: X1min,X1max
!   real(8)                              :: X2min,X2max
!   character(len=12)                   :: minx,miny,maxx,maxy
!   character(len=256)                   :: fname,dname
!   fname=get_filename(pname)
!   dname=get_filepath(pname)
!   Nx1=size(X1) ; Nk=size(X2,2)
!   Nl=5; if(present(nlines))Nl=nlines
!   open(719,file=adjustl(trim(pname)))
!   if(present(wlines))open(720,file=adjustl(trim(pname))//"_withlines")
!   do i=1,Nx1
!      count=mod(i,Nl)
!      do j=1,Nk
!         write(719,*)X1(i),X2(i,j),Y(i,j)
!         if(present(wlines).AND.count==0)write(720,*)X1(i),X2(i,j),Y(i,j)
!      enddo
!      write(719,*)""
!      if(present(wlines))write(720,*)""
!      if(present(wlines))write(720,*)""
!   enddo
!   close(719)
!   if(present(wlines))close(720)
!   X1min=minval(X1);X1max=maxval(X1)
!   X2min=minval(X2);X2max=maxval(X2)
!   write(minx,"(f12.4)")X1min
!   write(maxx,"(f12.4)")X1max
!   write(miny,"(f12.4)")X2min
!   write(maxy,"(f12.4)")X2max
!   call system("echo set term wxt 0 persist > "//adjustl(trim(pname))//".gp" )
!   call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set pm3d map >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set size square >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]>> "//adjustl(trim(pname))//".gp" )
!   call system("echo set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]>> "//adjustl(trim(pname))//".gp" )
!   call system("echo splot \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
!   !
!   call system("echo '#' >> "//adjustl(trim(pname))//".gp" )
!   call system("echo reset >> "//adjustl(trim(pname))//".gp" )
!   !
!   call system("echo set term wxt 1 persist >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set nokey >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set grid >> "//adjustl(trim(pname))//".gp" )
!   call system("echo set view 50,10,1,1 >> "//adjustl(trim(pname))//".gp" )
!   call system("echo splot \'"//trim(fname)//"\' with pm3d >> "//adjustl(trim(pname))//".gp" )
!   if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
!   call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
! end subroutine splot3D__
! !********************************************************************
! !********************************************************************
! !********************************************************************
