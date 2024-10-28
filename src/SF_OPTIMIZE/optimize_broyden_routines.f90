MODULE BROYDEN_FUNC_INTERFACE
  abstract interface 
     function broydn_func(x)
       real(8),dimension(:),intent(in) :: x
       real(8),dimension(size(x))      :: broydn_func
     end function broydn_func
  end interface
END MODULE BROYDEN_FUNC_INTERFACE


MODULE BROYDEN_ROUTINES
  USE BROYDEN_FUNC_INTERFACE
  implicit none
  integer, parameter             :: npar_arth=16,npar2_arth=8
  procedure(broydn_func),pointer :: funcv

contains

  subroutine fdjac(x,fvec,df)
    real(8), dimension(:), intent(in)    :: fvec
    real(8), dimension(:), intent(inout) :: x
    real(8), dimension(:,:), intent(out) :: df
    integer                              :: j,n
    real(8), parameter :: eps=1.0d-4
    real(8), dimension(size(x))          :: xsav,xph,h
    n=assert_eq4(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
    xsav=x
    h=eps*abs(xsav)
    where (h == 0.0) h=eps
    xph=xsav+h
    h=xph-xsav
    do j=1,n
       x(j)=xph(j)
       df(:,j)=(funcv(x)-fvec(:))/h(j)
       x(j)=xsav(j)
    end do
  end subroutine fdjac


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine lnsrch(xold,fold,g,p,x,f,stpmax,icheck,func)
    real(8), dimension(:), intent(in) :: xold,g
    real(8), dimension(:), intent(inout) :: p
    real(8), intent(in) :: fold,stpmax
    real(8), dimension(:), intent(out) :: x
    real(8), intent(out) :: f
    logical, intent(out) :: icheck
    interface
       function func(x)
         real(8)                          :: func
         real(8), dimension(:), intent(in):: x
       end function func
    end interface
    real(8), parameter :: alf=1.0e-4,tolx=epsilon(x)
    integer :: ndum
    real(8) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
    ndum=assert_eq4(size(g),size(p),size(x),size(xold),'lnsrch')
    icheck=.false.
    pabs=vabs(p(:))
    if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
    slope=dot_product(g,p)
    if (slope >= 0.0) then
       print*,'roundoff problem in lnsrch'
       stop
    endif
    alamin=tolx/maxval(abs(p(:))/max(abs(xold(:)),1.0d0))
    alam=1.0
    do
       x(:)=xold(:)+alam*p(:)
       f=func(x)
       if (alam < alamin) then
          x(:)=xold(:)
          icheck=.true.
          return
       else if (f <= fold+alf*alam*slope) then
          return
       else
          if (alam == 1.0) then
             tmplam=-slope/(2.0d0*(f-fold-slope))
          else
             rhs1=f-fold-alam*slope
             rhs2=f2-fold-alam2*slope
             a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
             b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
             if (a == 0.0) then
                tmplam=-slope/(2.0d0*b)
             else
                disc=b*b-3.0d0*a*slope
                if (disc < 0.0) then
                   tmplam=0.5d0*alam
                else if (b <= 0.0) then
                   tmplam=(-b+sqrt(disc))/(3.0d0*a)
                else
                   tmplam=-slope/(b+sqrt(disc))
                end if
             end if
             if (tmplam > 0.5d0*alam) tmplam=0.5d0*alam
          end if
       end if
       alam2=alam
       f2=f
       alam=max(tmplam,0.1d0*alam)
    end do
  end subroutine lnsrch


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine qrdcmp(a,c,d,sing)
    real(8), dimension(:,:), intent(inout) :: a
    real(8), dimension(:), intent(out) :: c,d
    logical, intent(out) :: sing
    integer :: k,n
    real(8) :: scale,sigma
    n=assert_eq4(size(a,1),size(a,2),size(c),size(d),'qrdcmp')
    sing=.false.
    do k=1,n-1
       scale=maxval(abs(a(k:n,k)))
       if (scale == 0.0) then
          sing=.true.
          c(k)=0.0
          d(k)=0.0
       else
          a(k:n,k)=a(k:n,k)/scale
          sigma=sign(vabs(a(k:n,k)),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),&
               matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
       end if
    end do
    d(n)=a(n,n)
    if (d(n) == 0.0) sing=.true.
  end subroutine qrdcmp


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine qrupdt(r,qt,u,v)
    real(8), dimension(:,:), intent(inout) :: r,qt
    real(8), dimension(:), intent(inout) :: u
    real(8), dimension(:), intent(in) :: v
    integer :: i,k,n
    n=assert_eqn((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),&
         size(v)/),'qrupdt')
    k=n+1-ifirstloc(u(n:1:-1) /= 0.0)
    if (k < 1) k=1
    do i=k-1,1,-1
       call rotate(r,qt,i,u(i),-u(i+1))
       u(i)=pythag(u(i),u(i+1))
    end do
    r(1,:)=r(1,:)+u(1)*v
    do i=1,k-1
       call rotate(r,qt,i,r(i,i),-r(i+1,i))
    end do
  contains
    subroutine rotate(r,qt,i,a,b)
      real(8), dimension(:,:), target, intent(inout) :: r,qt
      integer, intent(in) :: i
      real(8), intent(in) :: a,b
      real(8), dimension(size(r,1)) :: temp
      integer :: n
      real(8) :: c,fact,s
      n=assert_eq4(size(r,1),size(r,2),size(qt,1),size(qt,2),'rotate')
      if (a == 0.0) then
         c=0.0
         s=sign(1.0d0,b)
      else if (abs(a) > abs(b)) then
         fact=b/a
         c=sign(1.0d0/sqrt(1.0d0+fact**2),a)
         s=fact*c
      else
         fact=a/b
         s=sign(1.0d0/sqrt(1.0d0+fact**2),b)
         c=fact*s
      end if
      temp(i:n)=r(i,i:n)
      r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
      r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
      temp=qt(i,:)
      qt(i,:)=c*temp-s*qt(i+1,:)
      qt(i+1,:)=s*temp+c*qt(i+1,:)
    end subroutine rotate
    !    
    function pythag(a,b)
      real(8), intent(in) :: a,b
      real(8) :: pythag
      real(8) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if (absa > absb) then
         pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
         if (absb == 0.0) then
            pythag=0.0
         else
            pythag=absb*sqrt(1.0d0+(absa/absb)**2)
         end if
      end if
    end function pythag
  end subroutine qrupdt


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine rsolv(a,d,b)
    real(8), dimension(:,:), intent(in) :: a
    real(8), dimension(:), intent(in) :: d
    real(8), dimension(:), intent(inout) :: b
    integer :: i,n
    n=assert_eq4(size(a,1),size(a,2),size(b),size(d),'rsolv')
    b(n)=b(n)/d(n)
    do i=n-1,1,-1
       b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
    end do
  end subroutine rsolv


  !##################################################################
  !##################################################################
  !##################################################################



  function get_diag(mat) result(get_diag_dv)
    real(8),dimension(:,:),intent(in)      :: mat
    real(8), dimension(size(mat,1))        :: get_diag_dv
    integer                                :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  end function get_diag

  function assert_eq2(n1,n2,string)
    character(len=*), intent(in)           :: string
    integer, intent(in)                    :: n1,n2
    integer                                :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2

  function assert_eq4(n1,n2,n3,n4,string)
    character(len=*), intent(in)           :: string
    integer, intent(in)                    :: n1,n2,n3,n4
    integer                                :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq4'
    end if
  end function assert_eq4

  function assert_eqn(nn,string)
    character(len=*), intent(in) :: string
    integer, dimension(:), intent(in) :: nn
    integer :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eqn'
    end if
  end function assert_eqn

  function lower_triangle(j,k,extra)
    integer, intent(in)                    :: j,k
    integer, optional, intent(in)          :: extra
    logical, dimension(j,k)                :: lower_triangle
    integer                                :: n
    n=0;if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  contains
    function arth_i(first,increment,n)
      integer, intent(in)                  :: first,increment,n
      integer, dimension(n)                :: arth_i
      integer                              :: k,k2,temp
      if (n > 0) arth_i(1)=first
      if (n <= npar_arth) then
         do k=2,n
            arth_i(k)=arth_i(k-1)+increment
         end do
      else
         do k=2,npar2_arth
            arth_i(k)=arth_i(k-1)+increment
         end do
         temp=increment*npar2_arth
         k=npar2_arth
         do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
    end function arth_i
  end function lower_triangle

  function outerdiff(a,b) result(outerdiff_i)
    integer, dimension(:), intent(in) :: a,b
    integer, dimension(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  end function outerdiff

  function outerprod(a,b) result(outerprod_d)
    real(8), dimension(:), intent(in)      :: a,b
    real(8), dimension(size(a),size(b))    :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod

  subroutine put_diag(diagv,mat)
    real(8), dimension(:), intent(in)      :: diagv
    real(8), dimension(:,:), intent(inout) :: mat
    integer                                :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  end subroutine put_diag

  subroutine unit_matrix(mat)
    real(8), dimension(:,:), intent(out) :: mat
    integer :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0
    do i=1,n
       mat(i,i)=1.0
    end do
  end subroutine unit_matrix

  function vabs(v)
    real(8), dimension(:), intent(in) :: v
    real(8) :: vabs
    vabs=sqrt(dot_product(v,v))
  end function vabs

  function ifirstloc(mask)
    logical, dimension(:), intent(in) :: mask
    integer               :: ifirstloc
    integer, dimension(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  end function ifirstloc
END MODULE BROYDEN_ROUTINES



