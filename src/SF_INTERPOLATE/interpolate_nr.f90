module INTERPOLATE_NR
  implicit none
  private

  public :: locate
  public :: polint
  public :: polin2
contains



  function locate(xx,x)
    !This function performs a binary search in the table :f:var:`xx` and returns the index
    !:code:`j` such that :f:var:`x` lies between :code:`xx(j)` and :code:`xx(j+1)`. The table must
    !be monotonic, either ascending or descending. The function returns :code:`0` if :f:var:`x` lies
    !before :code:`xx(1)` and :code:`N` if it lies beyond :code:`xx(N)`, with :code:`N=size(xx)`, 
    !in the order of the table. The values :math:`x = xx(1)` and :math:`x = xx(N)` return :code:`1` 
    !and :code:`N-1` respectively. From Numerical Recipes.
    REAL(8), DIMENSION(:), INTENT(IN) :: xx  !monotonic table, ascending or descending
    REAL(8), INTENT(IN) :: x                 !value to locate in the table
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  END FUNCTION locate




  SUBROUTINE polint(xa,ya,x,y,dy)
    !This subroutine evaluates the polynomial of degree :code:`n-1` through the :code:`n` points
    !(:f:var:`xa`, :f:var:`ya`) at the point :f:var:`x`, using Neville's algorithm. It returns the
    !interpolated value :code:`y` and an estimate :code:`dy` of the error. The program stops if
    !the sizes of :f:var:`xa` and :f:var:`ya` differ, or if two points of :f:var:`xa` coincide.
    !From Numerical Recipes.
    REAL(8), DIMENSION(:), INTENT(IN) :: xa,ya  !data points: abscissas xa(i), all distinct, and ordinates ya(i)=f(xa(i))
    REAL(8), INTENT(IN)          :: x           !point where the polynomial is evaluated
    REAL(8), INTENT(OUT)         :: y,dy
    INTEGER                      :: m,n,ns
    REAL(8), DIMENSION(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.d0))then
          print*,'polint: calculation failure'
          stop
       endif
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  END SUBROUTINE polint


  subroutine polin2(x1a,x2a,ya,x1,x2,y,dy)
    !This subroutine evaluates the 2-dimensional polynomial interpolation of the function 
    !:f:var:`ya` tabulated on the grid :f:var:`x1a` :math:`\times` :f:var:`x2a`, at the point
    !(:f:var:`x1`, :f:var:`x2`). It applies :f:func_inline:`polint` along the second variable 
    !for each point of :f:var:`x1a`, and then along the first. It returns the interpolated value
    !:code:`y` and an estimate :code:`dy` of the error, obtained from the last interpolation.
    !The program stops if the shape of :f:var:`ya` is not :code:`(size(x1a),size(x2a))`.
    !From Numerical Recipes.
    real(8), dimension(:), intent(in)   :: x1a,x2a  !grid points along the first and the second variable
    real(8), dimension(:,:), intent(in) :: ya       !function values on the grid, ya(i,j)=f(x1a(i),x2a(j))
    real(8), intent(in)                 :: x1,x2    !coordinates of the point where the polynomial is evaluated
    real(8), intent(out)                :: y,dy
    integer                             :: j,m,ndum
    real(8), dimension(size(x1a))       :: ymtmp
    real(8), dimension(size(x2a))       :: yntmp
    m = size(x1a);if(m/=size(ya,1))stop "POLINT: wrong dimensions m"
    ndum=size(x2a);if(ndum/=size(ya,2))stop "POLINT: wrong dimensions ndum"
    do j=1,m
       yntmp=ya(j,:)
       call polint(x2a,yntmp,x2,ymtmp(j),dy)
    end do
    call polint(x1a,ymtmp,x1,y,dy)
  end subroutine polin2


  function iminloc(arr)
    real(8), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc

  function assert_eq(n1,n2,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq
    if (n1 == n2) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq'
    end if
  end function assert_eq

END MODULE INTERPOLATE_NR
