! This is an error function of a complex argument, which uses wpop(z).
elemental  function zerf( z )
  !This function returns the error function :math:`\mathrm{erf}(z)` of the complex argument :f:var:`z`, computed as
  !:math:`1 - w(iz)\,e^{-z^2}`, with :math:`w` the Faddeeva function. This is evaluated with a modified version of the algorithm
  !680 of G.P.M. Poppe and C.M.J. Wijers (ACM Trans. Math. Software 16, 1990), rewritten in Fortran 2008. The accuracy of
  !:math:`w(z)` is 14 significant digits in the first and second quadrants and 13 in the third and fourth ones, outside a
  !circular region of radius 0.126 around a zero of the function. The function is elemental.
  !
  complex(8), intent(in) :: z  !complex argument
  complex(8)             :: zerf
  zerf = one - wpop( cmplxj * z ) * exp( - z**2 )
end function zerf


! A modified version of algorithm 680, rewritten in Fortran 2008.
! G.P.M. Poppe, C.M.J. Wijers, More efficient computation of
! the complex error-function, ACM Trans. Math. Software 16:38-46, 1990.
!  and
! G.P.M. Poppe, C.M.J. Wijers, Algorithm 680, Evaluation of the
! complex error function, ACM Trans. Math. Software 16:47, 1990.
!
! Given a complex number z, this function computes
! the value of the Faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
! where erfc is the complex complementary error-function and i
! means sqrt(-1).  The accuracy of the algorithm for z in the 1st
! and 2nd quadrant is 14 significant digits; in the 3rd and 4th
! it is 13 significant digits outside a circular region with radius
! 0.126 around a zero of the function.
elemental function wpop( z )

  complex(8),intent(in) :: z
  complex(8)            :: wpop
  real(8),parameter     :: factor = 1.12837916709551257388d0    !  factor is 2/sqrt(pi)
  logical               :: a, b
  real(8)               :: xabs, yabs, x, y, qrho, xabsq, yquad, c, daux
  real(8)               :: h, h2, qlambda, rx, ry, sx, sy, tx, ty, u1, u2, v1, v2, w1, xaux
  real(8)               :: xquad, xsum, ysum, xi, yi, u, v
  integer               :: i, j, kapn, n, np1, nu
  !
  ! To avoid the complier uninitialised varning
  h2 = zero
  xi = dreal( z )
  yi = dimag( z )
  xabs = abs( xi )
  yabs = abs( yi )
  x = xabs / 6.3d0
  y = yabs / 4.4d0
  qrho = x**2 + y**2
  xabsq = xabs**2
  xquad = xabsq - yabs**2
  yquad = 2*xabs*yabs
  !
  a = qrho .lt. 0.085264d0
  !
  branch1: if ( a ) then

     ! If ( qrho .lt. 0.085264 ) then the Faddeeva-function is evaluated
     !  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
     !  n is the minimum number of terms needed to obtain the required
     !  accuracy

     qrho  = ( one - 0.85d0 * y ) * sqrt( qrho )
     n = nint( 6.0d0 + 72.0d0 * qrho )
     j = 2 * n + 1
     xsum = one / real( j, kind=8 )
     ysum = zero

     do i = n, 1, -1
        j = j - 2
        xaux = ( xsum * xquad - ysum * yquad) / real( i, kind=8 )
        ysum = ( xsum * yquad + ysum * xquad) / real( i, kind=8 )
        xsum = xaux + one / real( j, kind=8 )
     end do

     u1 = - factor * ( xsum * yabs + ysum * xabs ) + one
     v1 = factor * ( xsum * xabs - ysum * yabs )
     daux = exp( -xquad )
     u2 = daux * cos( yquad )
     v2 = -daux * sin( yquad )
     u = u1 * u2 - v1 * v2
     v = u1 * v2 + v1 * u2
  else

     bran2: if ( qrho .gt. one ) then

        ! If ( qrho .gt. 1) then w(z) is evaluated using the laplace
        ! continued fraction. nu is the minimum number of terms needed
        ! to obtain the required accuracy.

        h = zero
        kapn = 0
        qrho = sqrt( qrho )
        nu = int( 3.0d0 + (1442.0d0 / ( 26.0d0 * qrho &
             + 77.0d0 )))

     else

        ! If ( qrho .ge. 0.085264 .and. qrho .le. one ) then
        ! w(z) is evaluated by a truncated Taylor expansion,
        ! where the Laplace continued fraction is used to calculate
        ! the derivatives of w(z). KAPN is the minimum number of terms
        ! in the Taylor expansion needed to obtain the required accuracy.
        ! NU is the minimum number of terms of the continued fraction
        ! needed to calculate the derivatives with the required accuracy.

        qrho = ( one - y ) * sqrt( one - qrho )
        h = 1.88d0 * qrho
        h2 = two * h
        kapn = nint( 7.0d0 + 34.0d0 * qrho )
        nu   = nint( 16.0d0 + 26.0d0 * qrho )

     end if bran2

     b = h .gt. zero

     ! To avoid uninitialise compiler warning. qlambda is used
     ! only if (b), so can define to any value otherwise.
     qlambda = zero 
     if ( b ) qlambda = h2**kapn

     rx = zero
     ry = zero
     sx = zero
     sy = zero

     do n = nu, 0, -1
        np1 = n + 1
        tx = yabs + h + np1 * rx
        ty = xabs - np1 * ry
        c = half / (tx**2 + ty**2)
        rx = c * tx
        ry = c * ty
        if ( b .and. n .le. kapn ) then
           tx = qlambda + sx
           sx = rx*tx - ry*sy
           sy = ry*tx + rx*sy
           qlambda = qlambda / h2
        end if
     end do

     if ( h .eq. zero ) then
        u = factor * rx
        v = factor * ry
     else
        u = factor * sx
        v = factor * sy
     end if

     if ( yabs .eq. zero ) u = exp( -xabs**2 )

  end if branch1

  ! Evaluation of w(z) in the other quadrants

  if ( yi .lt. zero ) then

     if ( a ) then
        u2 = two * u2
        v2 = two * v2
     else
        xquad = -xquad
        w1 = two * exp( xquad )
        u2 = w1 * cos( yquad )
        v2 = -w1 * sin( yquad )
     end if

     u = u2 - u
     v = v2 - v
     if ( xi .gt. zero ) v = -v
  else
     if ( xi .lt. zero ) v = -v
  end if

  wpop = cmplx( u, v, kind=8 )

end function wpop















