  !*****************************************************************************80
  !
  !! AIRYA computes Airy functions and their derivatives.
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by 
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
  !    incorporate this routine into a user program that the copyright 
  !    is acknowledged.
  !
  !  Modified:
  !
  !    30 June 2012
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !       
  !  Parameters:
  !
  !    Input, real(8) X, the argument of the Airy function.
  !
  !    Output, real(8) AI, BI, AD, BD, the values of Ai(x), Bi(x),
  !    Ai'(x), Bi'(x).
  !
subroutine airya ( x, ai, bi, ad, bd )
  !AIRYA computes Airy functions and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang, Jianming Jin. FORTRAN90 version by John Burkardt.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument of the Airy function
  !@param ai: out: Ai(x)
  !@param bi: out: Bi(x)
  !@param ad: out: Ai'(x)
  !@param bd: out: Bi'(x)
  !
  implicit none

  real(8) :: ad
  real(8) :: ai
  real(8) :: bd
  real(8) :: bi
  real(8) :: c1
  real(8) :: c2
  real(8) :: pir
  real(8) :: sr3
  real(8) :: vi1
  real(8) :: vi2
  real(8) :: vj1
  real(8) :: vj2
  real(8) :: vk1
  real(8) :: vk2
  real(8) :: vy1
  real(8) :: vy2
  real(8) :: x
  real(8) :: xa
  real(8) :: xq
  real(8) :: z

  xa = abs ( x )
  pir = 0.318309886183891D+00
  c1 = 0.355028053887817D+00
  c2 = 0.258819403792807D+00
  sr3 = 1.732050807568877D+00
  z = xa ** 1.5D+00 / 1.5D+00
  xq = sqrt ( xa )

  call ajyik ( z, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )

  if ( x == 0.0D+00 ) then
     ai = c1
     bi = sr3 * c1
     ad = - c2
     bd = sr3 * c2
  else if ( 0.0D+00 < x ) then
     ai = pir * xq / sr3 * vk1
     bi = xq * ( pir * vk1 + 2.0D+00 / sr3 * vi1 )
     ad = - xa / sr3 * pir * vk2
     bd = xa * ( pir * vk2 + 2.0D+00 / sr3 * vi2 )
  else
     ai = 0.5D+00 * xq * ( vj1 - vy1 / sr3 )
     bi = - 0.5D+00 * xq * ( vj1 / sr3 + vy1 )
     ad = 0.5D+00 * xa * ( vj2 + vy2 / sr3 )
     bd = 0.5D+00 * xa * ( vj2 / sr3 - vy2 )
  end if

  return
end subroutine airya




  !*****************************************************************************80
  !
  !! AIRYB computes Airy functions and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 June 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, argument of Airy function.
  !
  !    Output, real(8) AI, Ai(x).
  !
  !    Output, real(8) BI, Bi(x).
  !
  !    Output, real(8) AD, Ai'(x).
  !
  !    Output, real(8) BD, Bi'(x).
  !
subroutine airyb ( x, ai, bi, ad, bd )
  !AIRYB computes Airy functions and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: argument of Airy function
  !@param ai: out: Ai(x)
  !@param bi: out: Bi(x)
  !@param ad: out: Ai'(x)
  !@param bd: out: Bi'(x)
  !
  implicit none

  real(8)    :: ad
  real(8)    :: ai
  real(8)    :: bd
  real(8)    :: bi
  real(8)    :: c1
  real(8)    :: c2
  real(8)    :: ck(41)
  real(8)    :: df
  real(8)    :: dg
  real(8)    :: dk(41)
  real(8)    :: eps
  real(8)    :: fx
  real(8)    :: gx
  integer(4) :: k
  integer(4) :: km
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: rp
  real(8)    :: sad
  real(8)    :: sai
  real(8)    :: sbd
  real(8)    :: sbi
  real(8)    :: sda
  real(8)    :: sdb
  real(8)    :: sr3
  real(8)    :: ssa
  real(8)    :: ssb
  real(8)    :: x
  real(8)    :: xa
  real(8)    :: xar
  real(8)    :: xcs
  real(8)    :: xe
  real(8)    :: xf
  real(8)    :: xm
  real(8)    :: xp1
  real(8)    :: xq
  real(8)    :: xr1
  real(8)    :: xr2
  real(8)    :: xss

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  c1 = 0.355028053887817D+00
  c2 = 0.258819403792807D+00
  sr3 = 1.732050807568877D+00
  xa = abs ( x )
  xq = sqrt ( xa )

  if ( x <= 0.0D+00 ) then
     xm = 8.0D+00
  else
     xm = 5.0D+00
  end if

  if ( x == 0.0D+00 ) then
     ai = c1
     bi = sr3 * c1
     ad = -c2
     bd = sr3 * c2
     return
  end if

  if ( xa <= xm ) then

     fx = 1.0D+00
     r = 1.0D+00
     do k = 1, 40
        r = r * x / ( 3.0D+00 * k ) * x / ( 3.0D+00 * k - 1.0D+00 ) * x
        fx = fx + r
        if ( abs ( r ) < abs ( fx ) * eps ) then
           exit
        end if
     end do

     gx = x
     r = x
     do k = 1, 40
        r = r * x / ( 3.0D+00 * k ) * x / ( 3.0D+00 * k + 1.0D+00 ) * x
        gx = gx + r
        if ( abs ( r ) < abs ( gx ) * eps ) then
           exit
        end if
     end do

     ai = c1 * fx - c2 * gx
     bi = sr3 * ( c1 * fx + c2 * gx )
     df = 0.5D+00 * x * x
     r = df
     do k = 1, 40
        r = r * x / ( 3.0D+00 * k ) * x / ( 3.0D+00 * k + 2.0D+00 ) * x
        df = df + r 
        if ( abs ( r ) < abs ( df ) * eps ) then
           exit
        end if
     end do

     dg = 1.0D+00
     r = 1.0D+00
     do k = 1, 40
        r = r * x / ( 3.0D+00 * k ) * x / ( 3.0D+00 * k - 2.0D+00 ) * x
        dg = dg + r
        if ( abs ( r ) < abs ( dg ) * eps ) then
           exit
        end if
     end do

     ad = c1 * df - c2 * dg
     bd = sr3 * ( c1 * df + c2 * dg )

  else

     xe = xa * xq / 1.5D+00
     xr1 = 1.0D+00 / xe
     xar = 1.0D+00 / xq
     xf = sqrt ( xar )
     rp = 0.5641895835477563D+00
     r = 1.0D+00
     do k = 1, 40
        r = r * ( 6.0D+00 * k - 1.0D+00 ) &
             / 216.0D+00 * ( 6.0D+00 * k - 3.0D+00 ) &
             / k * ( 6.0D+00 * k - 5.0D+00 ) / ( 2.0D+00 * k - 1.0D+00 )
        ck(k) = r
        dk(k) = - ( 6.0D+00 * k + 1.0D+00 ) / ( 6.0D+00 * k - 1.0D+00 ) * ck(k)
     end do

     km = int ( 24.5D+00 - xa )

     if ( xa < 6.0D+00 ) then
        km = 14
     end if

     if ( 15.0D+00 < xa ) then
        km = 10
     end if

     if ( 0.0D+00 < x ) then
        sai = 1.0D+00
        sad = 1.0D+00
        r = 1.0D+00
        do k = 1, km
           r = - r * xr1
           sai = sai + ck(k) * r
           sad = sad + dk(k) * r
        end do
        sbi = 1.0D+00
        sbd = 1.0D+00
        r = 1.0D+00
        do k = 1, km
           r = r * xr1
           sbi = sbi + ck(k) * r
           sbd = sbd + dk(k) * r
        end do
        xp1 = exp ( - xe )
        ai = 0.5D+00 * rp * xf * xp1 * sai
        bi = rp * xf / xp1 * sbi
        ad = -0.5D+00 * rp / xf * xp1 * sad
        bd = rp / xf / xp1 * sbd
     else
        xcs = cos ( xe + pi / 4.0D+00 )
        xss = sin ( xe + pi / 4.0D+00 )
        ssa = 1.0D+00
        sda = 1.0D+00
        r = 1.0D+00
        xr2 = 1.0D+00 / ( xe * xe )
        do k = 1, km
           r = - r * xr2
           ssa = ssa + ck(2*k) * r
           sda = sda + dk(2*k) * r
        end do
        ssb = ck(1) * xr1
        sdb = dk(1) * xr1
        r = xr1
        do k = 1, km
           r = - r * xr2
           ssb = ssb + ck(2*k+1) * r
           sdb = sdb + dk(2*k+1) * r
        end do
        ai = rp * xf * ( xss * ssa - xcs * ssb )
        bi = rp * xf * ( xcs * ssa + xss * ssb )
        ad = -rp / xf * ( xcs * sda + xss * sdb )
        bd =  rp / xf * ( xss * sda - xcs * sdb )
     end if

  end if

  return
end subroutine airyb




  !*****************************************************************************80
  !
  !! AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
  !
  !   Discussion:
  !
  !    Compute the first NT zeros of Airy functions Ai(x) and Ai'(x), 
  !    a and a', and the associated values of Ai(a') and Ai'(a); and 
  !    the first NT zeros of Airy functions Bi(x) and Bi'(x), b and
  !    b', and the associated values of Bi(b') and Bi'(b).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    14 March 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Input, integer(4) KF, the function code.
  !    1 for Ai(x) and Ai'(x);
  !    2 for Bi(x) and Bi'(x).
  !
  !    Output, real(8) XA(m), a, the m-th zero of Ai(x) or
  !    b, the m-th zero of Bi(x).
  !
  !    Output, real(8) XB(m), a', the m-th zero of Ai'(x) or
  !    b', the m-th zero of Bi'(x).
  !
  !    Output, real(8) XC(m), Ai(a') or Bi(b').
  !
  !    Output, real(8) XD(m), Ai'(a) or Bi'(b)
  !
subroutine airyzo ( nt, kf, xa, xb, xc, xd )
  !AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
  !
  !Compute the first NT zeros of Airy functions Ai(x) and Ai'(x),
  !a and a', and the associated values of Ai(a') and Ai'(a); and
  !the first NT zeros of Airy functions Bi(x) and Bi'(x), b and
  !b', and the associated values of Bi(b') and Bi'(b).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param nt: in: the number of zeros
  !@param kf: in: the function code. 1 for Ai(x) and Ai'(x); 2 for Bi(x) and Bi'(x)
  !@param xa: out: a, the m-th zero of Ai(x) or b, the m-th zero of Bi(x)
  !@param xb: out: a', the m-th zero of Ai'(x) or b', the m-th zero of Bi'(x)
  !@param xc: out: Ai(a') or Bi(b')
  !@param xd: out: Ai'(a) or Bi'(b)
  !
  implicit none

  integer(4) :: nt

  real(8)    :: ad
  real(8)    :: ai
  real(8)    :: bd
  real(8)    :: bi
  integer(4) :: i
  integer(4) :: kf
  real(8)    :: pi
  real(8)    :: rt
  real(8)    :: rt0
  real(8)    :: u
  real(8)    :: u1
  real(8)    :: x
  real(8)    :: xa(nt)
  real(8)    :: xb(nt)
  real(8)    :: xc(nt)
  real(8)    :: xd(nt)

  pi = 3.141592653589793D+00

  do i = 1, nt

     if (kf == 1) then
        u = 3.0D+00 * pi * ( 4.0D+00 * i - 1 ) / 8.0D+00
        u1 = 1.0D+00 / ( u * u )
        rt0 = - ( u * u ) ** ( 1.0 / 3.0 ) &
             * (((( -15.5902D+00 * u1 + 0.929844D+00 ) * u1 &
             - 0.138889D+00 ) * u1 + 0.10416667D+00 ) * u1 + 1.0D+00 )
     else if ( kf == 2 ) then
        if ( i == 1 ) then
           rt0 = -1.17371D+00
        else
           u = 3.0D+00 * pi * ( 4.0D+00 * i - 3.0D+00 ) / 8.0D+00
           u1 = 1.0D+00 / ( u * u )
           rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) &
                * (((( -15.5902D+00 * u1 + 0.929844D+00 ) * u1 &
                - 0.138889D+00 ) * u1 + 0.10416667D+00 ) * u1 + 1.0D+00 )
        end if
     end if

     do

        x = rt0
        call airyb ( x, ai, bi, ad, bd )

        if ( kf == 1 ) then
           rt = rt0 - ai / ad
        else
           rt = rt0 - bi / bd
        end if

        if ( abs ( ( rt - rt0 ) / rt ) <= 1.0D-09 ) then
           exit
        end if
        rt0 = rt

     end do

     xa(i) = rt
     if ( kf == 1 ) then
        xd(i) = ad
     else
        xd(i) = bd
     end if

  end do

  do i = 1, nt

     if ( kf == 1 ) then
        if ( i == 1 ) then
           rt0 = -1.01879D+00
        else
           u = 3.0D+00 * pi * ( 4.0D+00 * i - 3.0D+00 ) / 8.0D+00
           u1 = 1.0D+00 / ( u * u )
           rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) &
                * (((( 15.0168D+00 * u1 - 0.873954D+00 ) &
                * u1 + 0.121528D+00 ) * u1 - 0.145833D+00 ) * u1 + 1.0D+00 )
        end if
     else if ( kf == 2 ) then
        if ( i == 1 ) then
           rt0 = -2.29444D+00
        else
           u = 3.0D+00 * pi * ( 4.0D+00 * i - 1.0D+00 ) / 8.0D+00
           u1 = 1.0D+00 / ( u * u )
           rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) &
                * (((( 15.0168D+00 * u1 - 0.873954D+00 ) &
                * u1 + 0.121528D+00 ) * u1 - 0.145833D+00 ) * u1 + 1.0D+00 )
        end if
     end if

     do

        x = rt0
        call airyb ( x, ai, bi, ad, bd )

        if ( kf == 1 ) then
           rt = rt0 - ad / ( ai * x )
        else
           rt = rt0 - bd / ( bi * x )
        end if

        if ( abs ( ( rt - rt0 ) / rt ) <= 1.0D-09 ) then
           exit
        end if

        rt0 = rt

     end do

     xb(i) = rt
     if ( kf == 1 ) then
        xc(i) = ai
     else
        xc(i) = bi
     end if

  end do

  return
end subroutine airyzo




  !*****************************************************************************80
  !
  !! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
  !
  !  Discussion: 
  !
  !    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
  !    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.  X should not be zero.
  !
  !    Output, real(8) VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
  !    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
  !    K1/3(x), K2/3(x).
  !
subroutine ajyik ( x, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )
  !AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
  !
  !Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions
  !Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument. X should not be zero
  !@param vj1: out: J1/3(x)
  !@param vj2: out: J2/3(x)
  !@param vy1: out: Y1/3(x)
  !@param vy2: out: Y2/3(x)
  !@param vi1: out: I1/3(x)
  !@param vi2: out: I2/3(x)
  !@param vk1: out: K1/3(x)
  !@param vk2: out: K2/3(x)
  !
  implicit none

  real(8)    :: a0
  real(8)    :: b0
  real(8)    :: c0
  real(8)    :: ck
  real(8)    :: gn
  real(8)    :: gn1
  real(8)    :: gn2
  real(8)    :: gp1
  real(8)    :: gp2
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: l
  real(8)    :: pi
  real(8)    :: pv1
  real(8)    :: pv2
  real(8)    :: px
  real(8)    :: qx
  real(8)    :: r
  real(8)    :: rp
  real(8)    :: rp2
  real(8)    :: rq
  real(8)    :: sk
  real(8)    :: sum
  real(8)    :: uj1
  real(8)    :: uj2
  real(8)    :: uu0
  real(8)    :: vi1
  real(8)    :: vi2
  real(8)    :: vil
  real(8)    :: vj1
  real(8)    :: vj2
  real(8)    :: vjl
  real(8)    :: vk1
  real(8)    :: vk2
  real(8)    :: vl
  real(8)    :: vsl
  real(8)    :: vv
  real(8)    :: vv0
  real(8)    :: vy1
  real(8)    :: vy2
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xk

  if ( x == 0.0D+00 ) then
     vj1 = 0.0D+00
     vj2 = 0.0D+00
     vy1 = -1.0D+300
     vy2 = 1.0D+300
     vi1 = 0.0D+00
     vi2 = 0.0D+00
     vk1 = -1.0D+300
     vk2 = -1.0D+300
     return
  end if

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  gp1 = 0.892979511569249D+00
  gp2 = 0.902745292950934D+00
  gn1 = 1.3541179394264D+00
  gn2 = 2.678938534707747D+00
  vv0 = 0.444444444444444D+00
  uu0 = 1.1547005383793D+00
  x2 = x * x

  if ( x < 35.0D+00 ) then
     k0 = 12
  else if ( x < 50.0D+00 ) then
     k0 = 10
  else
     k0 = 8
  end if

  if ( x <= 12.0D+00 ) then

     do l = 1, 2
        vl = l / 3.0D+00
        vjl = 1.0D+00
        r = 1.0D+00
        do k = 1, 40
           r = -0.25D+00 * r * x2 / ( k * ( k + vl ) )
           vjl = vjl + r
           if ( abs ( r ) < 1.0D-15 ) then
              exit
           end if
        end do

        a0 = ( 0.5D+00 * x ) ** vl
        if ( l == 1 ) then
           vj1 = a0 / gp1 * vjl
        else
           vj2 = a0 / gp2 * vjl
        end if

     end do

  else

     do l = 1, 2

        vv = vv0 * l * l
        px = 1.0D+00
        rp = 1.0D+00

        do k = 1, k0
           rp = - 0.78125D-02 * rp &
                * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                / ( k * ( 2.0D+00 * k - 1.0D+00 ) * x2 )
           px = px + rp
        end do

        qx = 1.0D+00
        rq = 1.0D+00
        do k = 1, k0
           rq = - 0.78125D-02 * rq &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
                / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
           qx = qx + rq
        end do

        qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
        xk = x - ( 0.5D+00 * l / 3.0D+00 + 0.25D+00 ) * pi
        a0 = sqrt ( rp2 / x )
        ck = cos ( xk )
        sk = sin ( xk )
        if ( l == 1) then
           vj1 = a0 * ( px * ck - qx * sk )
           vy1 = a0 * ( px * sk + qx * ck )
        else
           vj2 = a0 * ( px * ck - qx * sk )
           vy2 = a0 * ( px * sk + qx * ck )
        end if

     end do

  end if

  if ( x <= 12.0D+00 ) then

     do l = 1, 2

        vl = l / 3.0D+00
        vjl = 1.0D+00
        r = 1.0D+00
        do k = 1, 40
           r = -0.25D+00 * r * x2 / ( k * ( k - vl ) )
           vjl = vjl + r
           if ( abs ( r ) < 1.0D-15 ) then
              exit
           end if
        end do

        b0 = ( 2.0D+00 / x ) ** vl
        if ( l == 1 ) then
           uj1 = b0 * vjl / gn1
        else
           uj2 = b0 * vjl / gn2
        end if

     end do

     pv1 = pi / 3.0D+00
     pv2 = pi / 1.5D+00
     vy1 = uu0 * ( vj1 * cos ( pv1 ) - uj1 )
     vy2 = uu0 * ( vj2 * cos ( pv2 ) - uj2 )

  end if

  if ( x <= 18.0D+00 ) then

     do l = 1, 2
        vl = l / 3.0D+00
        vil = 1.0D+00
        r = 1.0D+00
        do k = 1, 40
           r = 0.25D+00 * r * x2 / ( k * ( k + vl ) )
           vil = vil + r
           if ( abs ( r ) < 1.0D-15 ) then
              exit
           end if
        end do

        a0 = ( 0.5D+00 * x ) ** vl

        if ( l == 1 ) then
           vi1 = a0 / gp1 * vil
        else
           vi2 = a0 / gp2 * vil
        end if

     end do

  else

     c0 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )

     do l = 1, 2
        vv = vv0 * l * l
        vsl = 1.0D+00
        r = 1.0D+00
        do k = 1, k0
           r = - 0.125D+00 * r &
                * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
           vsl = vsl + r
        end do
        if ( l == 1 ) then
           vi1 = c0 * vsl
        else
           vi2 = c0 * vsl
        end if
     end do

  end if

  if ( x <= 9.0D+00 ) then

     do l = 1, 2
        vl = l / 3.0D+00
        if ( l == 1 ) then
           gn = gn1
        else
           gn = gn2
        end if
        a0 = ( 2.0D+00 / x ) ** vl / gn
        sum = 1.0D+00
        r = 1.0D+00
        do k = 1, 60
           r = 0.25D+00 * r * x2 / ( k * ( k - vl ) )
           sum = sum + r
           if ( abs ( r ) < 1.0D-15 ) then
              exit
           end if
        end do

        if ( l == 1 ) then
           vk1 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi1 )
        else
           vk2 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi2 )
        end if

     end do

  else

     c0 = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )

     do l = 1, 2
        vv = vv0 * l * l
        sum = 1.0D+00
        r = 1.0D+00
        do k = 1, k0
           r = 0.125D+00 * r * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
           sum = sum + r
        end do
        if ( l == 1 ) then
           vk1 = c0 * sum
        else
           vk2 = c0 * sum
        end if
     end do

  end if

  return
end subroutine ajyik




  !*****************************************************************************80
  !
  !! ASWFA: prolate and oblate spheroidal angular functions of the first kind.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    13 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter.
  !
  !    Input, integer(4) N, the mode parameter, with N = M, M+1, ...
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, real(8) X, the argument of the angular function.
  !    |X| < 1.0.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Output, real(8) S1F, S1D, the angular function of the first
  !    kind and its derivative.
  !
subroutine aswfa ( m, n, c, x, kd, cv, s1f, s1d )
  !ASWFA: prolate and oblate spheroidal angular functions of the first kind.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter
  !@param n: in: the mode parameter, with N = M, M+1, ...
  !@param c: in: the spheroidal parameter
  !@param x: in: the argument of the angular function. :math:`|X| < 1.0`
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param cv: in: the characteristic value
  !@param s1f: out: the angular function of the first kind
  !@param s1d: out: its derivative
  !
  implicit none

  real(8)    :: a0
  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: cv
  real(8)    :: d0
  real(8)    :: d1
  real(8)    :: df(200)
  real(8)    :: eps
  integer(4) :: ip
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm2
  real(8)    :: r
  real(8)    :: s1d
  real(8)    :: s1f
  real(8)    :: su1
  real(8)    :: su2
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x1

  eps = 1.0D-14
  x0 = x
  x = abs ( x )

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  nm = 10 + int ( ( n - m ) / 2 + c )
  nm2 = nm / 2 - 2 
  call sdmn ( m, n, c, cv, kd, df )
  call sckb ( m, n, c, df, ck )
  x1 = 1.0D+00 - x * x

  if ( m == 0 .and. x1 == 0.0D+00 ) then
     a0 = 1.0D+00
  else
     a0 = x1 ** ( 0.5D+00 * m )
  end if

  su1 = ck(1)
  do k = 1, nm2
     r = ck(k+1) * x1 ** k
     su1 = su1 + r
     if ( 10 <= k .and. abs ( r / su1 ) < eps ) then
        exit
     end if
  end do

  s1f = a0 * x ** ip * su1

  if ( x == 1.0D+00 ) then

     if ( m == 0 ) then
        s1d = ip * ck(1) - 2.0D+00 * ck(2)
     else if ( m == 1 ) then
        s1d = -1.0D+100
     else if ( m == 2 ) then
        s1d = -2.0D+00 * ck(1)
     else if ( 3 <= m ) then
        s1d = 0.0D+00
     end if

  else

     d0 = ip - m / x1 * x ** ( ip + 1.0D+00 )
     d1 = -2.0D+00 * a0 * x ** ( ip + 1.0D+00 )
     su2 = ck(2)
     do k = 2, nm2
        r = k * ck(k+1) * x1 ** ( k - 1.0D+00 )
        su2 = su2 + r
        if ( 10 <= k .and. abs ( r / su2 ) < eps ) then
           exit
        end if
     end do

     s1d = d0 * a0 * su1 + d1 * su2

  end if

  if ( x0 < 0.0D+00 ) then
     if ( ip == 0 ) then
        s1d = -s1d
     else if ( ip == 1 ) then
        s1f = -s1f
     end if
  end if

  x = x0

  return
end subroutine aswfa




  !*****************************************************************************80
  !
  !! ASWFB: prolate and oblate spheroidal angular functions of the first kind.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    20 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter, m = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M+1, M+2, ...
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, real(8) X, the argument, with |X| < 1.0.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Output, real(8) S1F, S1D, the angular function of the first
  !    kind and its derivative.
  !
subroutine aswfb ( m, n, c, x, kd, cv, s1f, s1d )
  !ASWFB: prolate and oblate spheroidal angular functions of the first kind.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter, m = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M+1, M+2, ...
  !@param c: in: the spheroidal parameter
  !@param x: in: the argument, with :math:`|X| < 1.0`
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param cv: in: the characteristic value
  !@param s1f: out: the angular function of the first kind
  !@param s1d: out: its derivative
  !
  implicit none

  real(8)    :: c
  real(8)    :: cv
  real(8)    :: df(200)
  real(8)    :: eps
  integer(4) :: ip
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: mk
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm2
  real(8)    :: pd(0:251)
  real(8)    :: pm(0:251)
  real(8)    :: s1d
  real(8)    :: s1f
  real(8)    :: su1
  real(8)    :: sw
  real(8)    :: x

  eps = 1.0D-14

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  nm = 25 + int ( ( n - m ) / 2 + c )
  nm2 = 2 * nm + m
  call sdmn ( m, n, c, cv, kd, df )
  call lpmns ( m, nm2, x, pm, pd )
  su1 = 0.0D+00
  do k = 1, nm
     mk = m + 2 * ( k - 1 ) + ip
     su1 = su1 + df(k) * pm(mk)
     if ( abs ( sw - su1 ) < abs ( su1 ) * eps ) then
        exit
     end if
     sw = su1
  end do

  s1f = ( -1.0D+00 ) ** m * su1

  su1 = 0.0D+00
  do k = 1, nm
     mk = m + 2 * ( k - 1 ) + ip
     su1 = su1 + df(k) * pd(mk)
     if ( abs ( sw - su1 ) < abs ( su1 ) * eps ) then
        exit
     end if
     sw = su1
  end do

  s1d = ( -1.0D+00 ) ** m * su1

  return
end subroutine aswfb




  !*****************************************************************************80
  !
  !! BERNOA computes the Bernoulli number Bn.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    11 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the index.
  !
  !    Output, real(8) BN, the value of the N-th Bernoulli number.
  !
subroutine bernoa ( n, bn )
  !BERNOA computes the Bernoulli number Bn.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the index
  !@param bn: out: the value of the N-th Bernoulli number
  !
  implicit none

  integer(4) :: n

  real(8)    :: bn(0:n)
  integer(4) :: j
  integer(4) :: k
  integer(4) :: m
  real(8)    :: r
  real(8)    :: s

  bn(0) = 1.0D+00
  bn(1) = -0.5D+00

  do m = 2, n
     s = - ( 1.0D+00 / ( m + 1.0D+00 ) - 0.5D+00 )
     do k = 2, m - 1
        r = 1.0D+00
        do j = 2, k
           r = r * ( j + m - k ) / j
        end do
        s = s - r * bn(k)
     end do
     bn(m) = s
  end do

  do m = 3, n, 2
     bn(m) = 0.0D+00
  end do

  return
end subroutine bernoa




  !*****************************************************************************80
  !
  !! BERNOB computes the Bernoulli number Bn.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    11 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the index.
  !
  !    Output, real(8) BN, the value of the N-th Bernoulli number.
  !
subroutine bernob ( n, bn )
  !BERNOB computes the Bernoulli number Bn.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the index
  !@param bn: out: the value of the N-th Bernoulli number
  !
  implicit none

  integer(4) :: n

  real(8)    :: bn(0:n)
  integer(4) :: k
  integer(4) :: m
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: s
  real(8)    :: tpi

  tpi = 6.283185307179586D+00
  bn(0) = 1.0D+00
  bn(1) = -0.5D+00
  bn(2) = 1.0D+00 / 6.0D+00
  r1 = ( 2.0D+00 / tpi )**2

  do m = 4, n, 2

     r1 = - r1 * ( m - 1 ) * m / ( tpi * tpi )
     r2 = 1.0D+00

     do k = 2, 10000
        s = ( 1.0D+00 / k ) ** m
        r2 = r2 + s
        if ( s < 1.0D-15 ) then
           exit
        end if
     end do

     bn(m) = r1 * r2

  end do

  return
end subroutine bernob




  !*****************************************************************************80
  !
  !! BETA computes the Beta function B(p,q).
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by 
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
  !    incorporate this routine into a user program that the copyright 
  !    is acknowledged.
  !
  !  Modified:
  !
  !    12 March 2012
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real(8) P, Q, the parameters.
  !    0 < P, 0 < Q.
  !
  !    Output, real(8) BT, the value of B(P,Q).
  !
subroutine betaf ( p, q, bt )
  !BETA computes the Beta function B(p,q).
  !
  !Original FORTRAN77 version by Shanjie Zhang, Jianming Jin. FORTRAN90 version by John Burkardt.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param p: in: the parameter P, 0 < P
  !@param q: in: the parameter Q, 0 < Q
  !@param bt: out: the value of B(P,Q)
  !
  implicit none

  real(8) :: bt
  real(8) :: gp
  real(8) :: gpq
  real(8) :: gq
  real(8) :: p
  real(8) :: ppq
  real(8) :: q

  call gammaf ( p, gp )
  call gammaf ( q, gq )
  ppq = p + q
  call gammaf ( ppq, gpq )
  bt = gp * gq / gpq

  return
end subroutine betaf




  !*****************************************************************************80
  !
  !! BJNDD computes Bessel functions Jn(x) and first and second derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    11 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BJ(N+1), DJ(N+1), FJ(N+1), the values of 
  !    Jn(x), Jn'(x) and Jn''(x) in the last entries.
  !
subroutine bjndd ( n, x, bj, dj, fj )
  !BJNDD computes Bessel functions Jn(x) and first and second derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param bj: out: Jn(x)
  !@param dj: out: Jn'(x)
  !@param fj: out: Jn''(x) in the last entries
  !
  implicit none

  integer(4) :: n

  real(8)    :: bj(n+1)
  real(8)    :: bs
  real(8)    :: dj(n+1)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: fj(n+1)
  integer(4) :: k
  integer(4) :: m
  integer(4) :: mt
  integer(4) :: nt
  real(8)    :: x

  do nt = 1, 900
     mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) &
          - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
     if ( 20 < mt ) then
        exit
     end if
  end do

  m = nt
  bs = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-35
  do k = m, 0, -1
     f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
     if ( k <= n ) then
        bj(k+1) = f
     end if
     if ( k == 2 * int ( k / 2 ) ) then
        bs = bs + 2.0D+00 * f
     end if
     f0 = f1
     f1 = f
  end do

  do k = 0, n
     bj(k+1) = bj(k+1) / ( bs - f )
  end do

  dj(1) = -bj(2)
  fj(1) = -1.0D+00 * bj(1) - dj(1) / x
  do k = 1, n
     dj(k+1) = bj(k) - k * bj(k+1) / x
     fj(k+1) = ( k * k / ( x * x ) - 1.0D+00 ) * bj(k+1) - dj(k+1) / x
  end do

  return
end subroutine bjndd




  !*****************************************************************************80
  !
  !! CBK computes coefficients for oblate radial functions with small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    20 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, real(8) QT, ?
  !
  !    Input, real(8) CK(*), ?
  !
  !    Output, real(8) BK(*), the coefficients.
  !
subroutine cbk ( m, n, c, cv, qt, ck, bk )
  !CBK computes coefficients for oblate radial functions with small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param cv: in: the characteristic value
  !@param qt: in: not documented in the original header
  !@param ck: in: not documented in the original header
  !@param bk: out: the coefficients
  !
  implicit none

  real(8)    :: bk(200)
  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: cv
  real(8)    :: eps
  integer(4) :: i
  integer(4) :: i1
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: m
  integer(4) :: n
  integer(4) :: n2
  integer(4) :: nm
  real(8)    :: qt
  real(8)    :: r1
  real(8)    :: s1
  real(8)    :: sw
  real(8)    :: t
  real(8)    :: u(200)
  real(8)    :: v(200)
  real(8)    :: w(200)

  eps = 1.0D-14
  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if
  nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
  u(1) = 0.0D+00
  n2 = nm - 2
  do j = 2, n2
     u(j) = c * c
  end do

  do j = 1, n2
     v(j) = ( 2.0D+00 * j - 1.0D+00 - ip ) &
          * ( 2.0D+00 * ( j - m ) - ip ) + m * ( m - 1.0D+00 ) - cv
  end do

  do j = 1, nm - 1
     w(j) = ( 2.0D+00 * j - ip ) * ( 2.0D+00 * j + 1.0D+00 - ip )
  end do

  if ( ip == 0 ) then

     do k = 0, n2 - 1

        s1 = 0.0D+00
        i1 = k - m + 1

        do i = i1, nm
           if ( 0 <= i ) then
              r1 = 1.0D+00
              do j = 1, k
                 r1 = r1 * ( i + m - j ) / j
              end do
              s1 = s1 + ck(i+1) * ( 2.0D+00 * i + m ) * r1
              if ( abs ( s1 - sw ) < abs ( s1 ) * eps ) then
                 exit
              end if
              sw = s1
           end if
        end do

        bk(k+1) = qt * s1

     end do

  else if ( ip == 1 ) then

     do k = 0, n2 - 1

        s1 = 0.0D+00
        i1 = k - m + 1

        do i = i1, nm

           if ( 0 <= i ) then

              r1 = 1.0D+00
              do j = 1, k
                 r1 = r1 * ( i + m - j ) / j
              end do

              if ( 0 < i ) then
                 s1 = s1 + ck(i) * ( 2.0D+00 * i + m - 1 ) * r1
              end if
              s1 = s1 - ck(i+1) * ( 2.0D+00 * i + m ) * r1
              if ( abs ( s1 - sw ) < abs ( s1 ) * eps ) then
                 exit
              end if
              sw = s1

           end if

        end do

        bk(k+1) = qt * s1

     end do

  end if

  w(1) = w(1) / v(1)
  bk(1) = bk(1) / v(1)
  do k = 2, n2
     t = v(k) - w(k-1) * u(k)
     w(k) = w(k) / t
     bk(k) = ( bk(k) - bk(k-1) * u(k) ) / t
  end do

  do k = n2 - 1, 1, -1
     bk(k) = bk(k) - w(k) * bk(k+1)
  end do

  return
end subroutine cbk




  !*****************************************************************************80
  !
  !! CCHG computes the confluent hypergeometric function.
  !
  !  Discussion:
  !
  !    This function computes the confluent hypergeometric function
  !    M(a,b,z) with real parameters a, b and complex argument z.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    26 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, parameter values.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CHG, the value of M(a,b,z).
  !
subroutine cchg ( a, b, z, chg )
  !CCHG computes the confluent hypergeometric function.
  !
  !This function computes the confluent hypergeometric function
  !M(a,b,z) with real parameters a, b and complex argument z.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param z: in: the argument
  !@param chg: out: the value of M(a,b,z)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: b
  real(8)    :: ba
  complex(8) :: cfac
  complex(8) :: chg
  complex(8) :: chg1
  complex(8) :: chg2
  complex(8) :: chw
  complex(8) :: ci
  complex(8) :: cr
  complex(8) :: cr1
  complex(8) :: cr2
  complex(8) :: crg
  complex(8) :: cs1
  complex(8) :: cs2
  complex(8) :: cy0
  complex(8) :: cy1
  real(8)    :: g1
  real(8)    :: g2
  real(8)    :: g3
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: la
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nl
  integer(4) :: ns
  real(8)    :: phi
  real(8)    :: pi
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: y
  complex(8) :: z
  complex(8) :: z0

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00 )
  a0 = a
  a1 = a
  z0 = z

  if ( b == 0.0D+00 .or. b == - int ( abs ( b ) ) ) then
     chg = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
  else if ( a == 0.0D+00 .or. z == 0.0D+00 ) then
     chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  else if ( a == -1.0D+00 ) then
     chg = 1.0D+00 - z / b
  else if ( a == b ) then
     chg = exp ( z )
  else if ( a - b == 1.0D+00 ) then
     chg = ( 1.0D+00 + z / b ) * exp ( z )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
     chg = ( exp ( z ) - 1.0D+00 ) / z
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
     m = int ( - a )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, m
        cr = cr * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * z
        chg = chg + cr
     end do
  else

     x0 = real ( z, kind = 8 )
     if ( x0 < 0.0D+00 ) then
        a = b - a
        a0 = a
        z = - z
     end if

     if ( a < 2.0D+00 ) then
        nl = 0
     else
        nl = 1
        la = int ( a )
        a = a - la - 1.0D+00
     end if

     do n = 0, nl

        if ( 2.0D+00 <= a0 ) then
           a = a + 1.0D+00
        end if

        if ( cdabs ( z ) < 20.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

           chg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           crg = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           do j = 1, 500
              crg = crg * ( a + j - 1.0D+00 ) / ( j * ( b + j - 1.0D+00 ) ) * z
              chg = chg + crg
              if ( abs ( ( chg - chw ) / chg ) < 1.0D-15 ) then
                 exit
              end if
              chw = chg
           end do

        else

           call gammaf ( a, g1 )
           call gammaf ( b, g2 )
           ba = b - a
           call gammaf ( ba, g3 )
           cs1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           cs2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           cr2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

           do i = 1, 8
              cr1 = - cr1 * (     a + i - 1.0D+00 ) * ( a - b + i ) / ( z * i )
              cr2 =   cr2 * ( b - a + i - 1.0D+00 ) * ( i - a ) / ( z * i )
              cs1 = cs1 + cr1
              cs2 = cs2 + cr2
           end do

           x = real ( z, kind = 8 )
           y = imag ( z )

           if ( x == 0.0D+00 .and. 0.0D+00 <= y ) then
              phi = 0.5D+00 * pi
           else if ( x == 0.0D+00 .and. y <= 0.0D+00 ) then
              phi = -0.5D+00 * pi
           else
              phi = atan ( y / x )
           end if

           if ( -1.5D+00 * pi < phi .and. phi <= -0.5 * pi ) then
              ns = -1
           else if ( -0.5D+00 * pi < phi .and. phi < 1.5D+00 * pi ) then
              ns = 1
           end if

           if ( y == 0.0D+00 ) then
              cfac = cos ( pi * a )
           else
              cfac = exp ( ns * ci * pi * a )
           end if

           chg1 = g2 / g3 * z ** ( - a ) * cfac * cs1
           chg2 = g2 / g1 * exp ( z ) * z ** ( a - b ) * cs2
           chg = chg1 + chg2

        end if

        if ( n == 0 ) then
           cy0 = chg
        else if ( n == 1 ) then
           cy1 = chg
        end if

     end do

     if ( 2.0D+00 <= a0 ) then
        do i = 1, la - 1
           chg = ( ( 2.0D+00 * a - b + z ) * cy1 + ( b - a ) * cy0 ) / a
           cy0 = cy1
           cy1 = chg
           a = a + 1.0D+00
        end do
     end if

     if ( x0 < 0.0D+00 ) then
        chg = chg * exp ( - z )
     end if

  end if

  a = a1
  z = z0

  return
end subroutine cchg




  !*****************************************************************************80
  !
  !! CERF computes the error function and derivative for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex(8), the argument.
  !
  !    Output, complex(8) CER, CDER, the values of erf(z) and erf'(z).
  !
subroutine cerf ( z, cer, cder )
  !CERF computes the error function and derivative for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: the argument
  !@param cer: out: erf(z)
  !@param cder: out: erf'(z)
  !
  implicit none

  complex(8) :: c0
  complex(8) :: cder
  complex(8) :: cer
  complex(8) :: cs
  real(8)    :: ei1
  real(8)    :: ei2
  real(8)    :: eps
  real(8)    :: er
  real(8)    :: er0
  real(8)    :: er1
  real(8)    :: er2
  real(8)    :: eri
  real(8)    :: err
  integer(4) :: k
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: ss
  real(8)    :: w
  real(8)    :: w1
  real(8)    :: w2
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: y
  complex(8) :: z

  eps = 1.0D-12
  pi = 3.141592653589793D+00
  x = real ( z, kind = 8 )
  y = imag ( z )
  x2 = x * x

  if ( x <= 3.5D+00 ) then

     er = 1.0D+00
     r = 1.0D+00
     do k = 1, 100
        r = r * x2 / ( k + 0.5D+00 )
        er = er + r
        if ( abs ( er - w ) <= eps * abs ( er ) ) then
           exit
        end if
        w = er
     end do

     c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
     er0 = c0 * er

  else

     er = 1.0D+00
     r = 1.0D+00
     do k = 1, 12
        r = - r * ( k - 0.5D+00 ) / x2
        er = er + r
     end do
     c0 = exp ( - x2 ) / ( x * sqrt ( pi ) )
     er0 = 1.0D+00 - c0 * er

  end if

  if ( y == 0.0D+00 ) then
     err = er0
     eri = 0.0D+00
  else
     cs = cos ( 2.0D+00 * x * y )
     ss = sin ( 2.0D+00 * x * y )
     er1 = exp ( - x2 ) * ( 1.0D+00 - cs ) / ( 2.0D+00 * pi * x )
     ei1 = exp ( - x2 ) * ss / ( 2.0D+00 * pi * x )
     er2 = 0.0D+00
     do n = 1, 100
        er2 = er2 + exp ( - 0.25D+00 * n * n ) &
             / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x &
             - 2.0D+00 * x * cosh ( n * y ) * cs &
             + n * sinh ( n * y ) * ss )
        if ( abs ( ( er2 - w1 ) / er2 ) < eps ) then
           exit
        end if
        w1 = er2
     end do

     c0 = 2.0D+00 * exp ( - x2 ) / pi
     err = er0 + er1 + c0 * er2
     ei2 = 0.0D+00
     do n = 1, 100
        ei2 = ei2 + exp ( - 0.25D+00 * n * n ) &
             / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x &
             * cosh ( n * y ) * ss + n * sinh ( n * y ) * cs )
        if ( abs ( ( ei2 - w2 ) / ei2 ) < eps ) then
           exit
        end if
        w2 = ei2
     end do

     eri = ei1 + c0 * ei2

  end if

  cer = cmplx ( err, eri, kind = 8 )
  cder = 2.0D+00 / sqrt ( pi ) * exp ( - z * z )

  return
end subroutine cerf




  !*****************************************************************************80
  !
  !! CERROR computes the error function for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CER, the function value.
  !
subroutine cerror ( z, cer )
  !CERROR computes the error function for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param cer: out: the function value
  !
  implicit none

  real(8)    :: a0
  complex(8) :: c0
  complex(8) :: cer
  complex(8) :: cl
  complex(8) :: cr
  complex(8) :: cs
  integer(4) :: k
  real(8)    :: pi
  complex(8) :: z
  complex(8) :: z1

  a0 = abs ( z )
  c0 = exp ( - z * z )
  pi = 3.141592653589793D+00
  z1 = z

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = - z
  end if

  if ( a0 <= 5.8D+00 ) then    

     cs = z1
     cr = z1
     do k = 1, 120
        cr = cr * z1 * z1 / ( k + 0.5D+00 )
        cs = cs + cr
        if ( abs ( cr / cs ) < 1.0D-15 ) then
           exit
        end if
     end do

     cer = 2.0D+00 * c0 * cs / sqrt ( pi )

  else

     cl = 1.0D+00 / z1              
     cr = cl
     do k = 1, 13
        cr = -cr * ( k - 0.5D+00 ) / ( z1 * z1 )
        cl = cl + cr
        if ( abs ( cr / cl ) < 1.0D-15 ) then
           exit
        end if
     end do

     cer = 1.0D+00 - c0 * cl / sqrt ( pi )

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     cer = -cer
  end if

  return
end subroutine cerror




  !*****************************************************************************80
  !
  !! CERZO evaluates the complex zeros of the error function.
  !
  !  Discussion:
  !
  !    The modified Newton method is used.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Output, complex(8) ZO(NT), the zeros.
  !
subroutine cerzo ( nt, zo )
  !CERZO evaluates the complex zeros of the error function.
  !
  !The modified Newton method is used.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param nt: in: the number of zeros
  !@param zo: out: the zeros
  !
  implicit none

  integer(4) :: nt

  integer(4) :: i
  integer(4) :: it
  integer(4) :: j
  integer(4) :: nr
  real(8)    :: pi
  real(8)    :: pu
  real(8)    :: pv
  real(8)    :: px
  real(8)    :: py
  real(8)    :: w
  real(8)    :: w0
  complex(8) :: z
  complex(8) :: zd
  complex(8) :: zf
  complex(8) :: zfd
  complex(8) :: zgd
  complex(8) :: zo(nt)
  complex(8) :: zp
  complex(8) :: zq
  complex(8) :: zw

  pi = 3.141592653589793D+00

  do nr = 1, nt

     pu = sqrt ( pi * ( 4.0D+00 * nr - 0.5D+00 ) )
     pv = pi * sqrt ( 2.0D+00 * nr - 0.25D+00 )
     px = 0.5D+00 * pu - 0.5D+00 * log ( pv ) / pu
     py = 0.5D+00 * pu + 0.5D+00 * log ( pv ) / pu
     z = cmplx ( px, py, kind = 8 )
     it = 0

     do

        it = it + 1
        call cerf ( z, zf, zd )
        zp = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do i = 1, nr - 1
           zp = zp * ( z - zo(i) )
        end do
        zfd = zf / zp

        zq = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        do i = 1, nr - 1
           zw = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           do j = 1, nr - 1
              if ( j /= i ) then
                 zw = zw * ( z - zo(j) )
              end if
           end do
           zq = zq + zw
        end do

        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = abs ( z )

        if ( 50 < it .or. abs ( ( w - w0 ) / w ) <= 1.0D-11 ) then
           exit
        end if

     end do

     zo(nr) = z

  end do

  return
end subroutine cerzo




  !*****************************************************************************80
  !
  !! CFC computes the complex Fresnel integral C(z) and C'(z).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    26 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) ZF, ZD, the values of C(z) and C'(z).
  !
subroutine cfc ( z, zf, zd )
  !CFC computes the complex Fresnel integral C(z) and C'(z).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param zf: out: C(z)
  !@param zd: out: C'(z)
  !
  implicit none

  complex(8) :: c
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cg
  complex(8) :: cr
  real(8)    :: eps
  integer(4) :: k
  integer(4) :: m
  real(8)    :: pi
  real(8)    :: w0
  real(8)    :: wa
  real(8)    :: wa0
  complex(8) :: z
  complex(8) :: z0
  complex(8) :: zd
  complex(8) :: zf
  complex(8) :: zp
  complex(8) :: zp2

  eps = 1.0D-14
  pi = 3.141592653589793D+00
  w0 = abs ( z )
  zp = 0.5D+00 * pi * z * z
  zp2 = zp * zp
  z0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( z .eq. z0 ) then

     c = z0

  else if ( w0 <= 2.5D+00 ) then

     cr = z
     c = cr
     do k = 1, 80
        cr = -0.5D+00 * cr * ( 4.0D+00 * k - 3.0D+00 ) &
             / k / ( 2.0D+00 * k - 1.0D+00 ) &
             / ( 4.0D+00 * k + 1.0D+00 ) * zp2
        c = c + cr
        wa = abs ( c )
        if ( abs ( ( wa - wa0 ) / wa ) < eps .and. 10 < k ) then
           exit
        end if
        wa0 = wa
     end do

  else if ( 2.5D+00 < w0 .and. w0 < 4.5D+00 ) then

     m = 85
     c = z0
     cf1 = z0
     cf0 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
     do k = m, 0, -1
        cf = ( 2.0D+00 * k + 3.0D+00 ) * cf0 / zp - cf1
        if ( k .eq. int ( k / 2 ) * 2 ) then
           c = c + cf
        end if
        cf1 = cf0
        cf0 = cf
     end do
     c = sqrt ( 2.0D+00 / ( pi * zp ) ) * sin ( zp ) / cf * c

  else

     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 20
        cr = - 0.25D+00 * cr * ( 4.0D+00 * k - 1.0D+00 ) &
             * ( 4.0D+00 * k - 3.0D+00 ) / zp2
        cf = cf + cr
     end do
     cr = 1.0D+00 / ( pi * z * z )
     cg = cr
     do k = 1, 12
        cr = - 0.25D+00 * cr * ( 4.0D+00 * k + 1.0D+00 ) &
             * ( 4.0D+00 * k - 1.0D+00 ) / zp2
        cg = cg + cr
     end do
     c = 0.5D+00 + ( cf * sin ( zp ) - cg * cos ( zp ) ) / ( pi * z )

  end if

  zf = c
  zd = cos ( 0.5D+00 * pi * z * z )

  return
end subroutine cfc




  !*****************************************************************************80
  !
  !! CFS computes the complex Fresnel integral S(z) and S'(z).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    24 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) ZF, ZD, the values of S(z) and S'(z).
  !
subroutine cfs ( z, zf, zd )
  !CFS computes the complex Fresnel integral S(z) and S'(z).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param zf: out: S(z)
  !@param zd: out: S'(z)
  !
  implicit none

  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cg
  complex(8) :: cr
  real(8)    :: eps
  integer(4) :: k
  integer(4) :: m
  real(8)    :: pi
  complex(8) :: s
  real(8)    :: w0
  real(8)    :: wb
  real(8)    :: wb0
  complex(8) :: z
  complex(8) :: z0
  complex(8) :: zd
  complex(8) :: zf
  complex(8) :: zp
  complex(8) :: zp2

  eps = 1.0D-14
  pi = 3.141592653589793D+00
  w0 = abs ( z )
  zp = 0.5D+00 * pi * z * z
  zp2 = zp * zp
  z0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( z == z0 ) then

     s = z0

  else if ( w0 <= 2.5D+00 ) then

     s = z * zp / 3.0D+00
     cr = s
     do k = 1, 80
        cr = -0.5D+00 * cr * ( 4.0D+00 * k - 1.0D+00 ) / k &
             / ( 2.0D+00 * k + 1.0D+00 ) &
             / ( 4.0D+00 * k + 3.0D+00 ) * zp2
        s = s + cr
        wb = abs ( s )
        if ( abs ( wb - wb0 ) < eps .and. 10 < k ) then
           exit
        end if
        wb0 = wb
     end do

  else if ( 2.5D+00 < w0 .and. w0 < 4.5D+00 ) then

     m = 85
     s = z0
     cf1 = z0
     cf0 = cmplx ( 1.0D-30, 0.0D+00, kind = 8  )
     do k = m, 0, -1
        cf = ( 2.0D+00 * k + 3.0D+00 ) * cf0 / zp - cf1
        if ( k /= int ( k / 2 ) * 2 ) then
           s = s + cf
        end if
        cf1 = cf0
        cf0 = cf
     end do
     s = sqrt ( 2.0D+00 / ( pi * zp ) ) * sin ( zp ) / cf * s

  else

     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8  )
     cf = cmplx ( 1.0D+00, 0.0D+00, kind = 8  )
     do k = 1, 20
        cr = -0.25D+00 * cr * ( 4.0D+00 * k - 1.0D+00 ) &
             * ( 4.0D+00 * k - 3.0D+00 ) / zp2
        cf = cf + cr
     end do
     cr = 1.0D+00 / ( pi * z * z )
     cg = cr
     do k = 1, 12
        cr = -0.25D+00 * cr * ( 4.0D+00 * k + 1.0D+00 ) &
             * ( 4.0D+00 * k - 1.0D+00 ) / zp2
        cg = cg + cr
     end do
     s = 0.5D+00 - ( cf * cos ( zp ) + cg * sin ( zp ) ) / ( pi * z )

  end if

  zf = s
  zd = sin ( 0.5D+00 * pi * z * z )

  return
end subroutine cfs




  !*****************************************************************************80
  !
  !! CGAMA computes the Gamma function for complex argument.
  !
  !  Discussion:
  !
  !    This procedcure computes the gamma function Γ(z) or ln[Γ(z)]
  !    for a complex argument
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    26 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, Y, the real and imaginary parts of 
  !    the argument Z.
  !
  !    Input, integer(4) KF, the function code.
  !    0 for ln[Γ(z)]
  !    1 for Γ(z)
  !
  !    Output, real(8) GR, GI, the real and imaginary parts of
  !    the selected function.
  !
subroutine cgama ( x, y, kf, gr, gi )
  !CGAMA computes the Gamma function for complex argument.
  !
  !This procedure computes the gamma function Gamma(z) or ln[Gamma(z)]
  !for a complex argument
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the real part of the argument Z
  !@param y: in: the imaginary part of the argument Z
  !@param kf: in: the function code. 0 for ln[Gamma(z)] 1 for Gamma(z)
  !@param gr: out: the real part of the selected function
  !@param gi: out: the imaginary part of the selected function
  !
  implicit none

  real(8), save, dimension ( 10 ) :: a = (/ &
       8.333333333333333D-02, -2.777777777777778D-03, &
       7.936507936507937D-04, -5.952380952380952D-04, &
       8.417508417508418D-04, -1.917526917526918D-03, &
       6.410256410256410D-03, -2.955065359477124D-02, &
       1.796443723688307D-01, -1.39243221690590D+00 /)
  real(8)    :: g0
  real(8)    :: gi
  real(8)    :: gi1
  real(8)    :: gr
  real(8)    :: gr1
  integer(4) :: j
  integer(4) :: k
  integer(4) :: kf
  integer(4) :: na
  real(8)    :: pi
  real(8)    :: si
  real(8)    :: sr
  real(8)    :: t
  real(8)    :: th
  real(8)    :: th1
  real(8)    :: th2
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x1
  real(8)    :: y
  real(8)    :: y1
  real(8)    :: z1
  real(8)    :: z2

  pi = 3.141592653589793D+00

  if ( y == 0.0D+00 .and. x == int ( x ) .and. x <= 0.0D+00 ) then
     gr = 1.0D+300
     gi = 0.0D+00
     return
  else if ( x < 0.0D+00 ) then
     x1 = x
     y1 = y
     x = -x
     y = -y
  end if

  x0 = x

  if ( x <= 7.0D+00 ) then
     na = int ( 7 - x )
     x0 = x + na
  end if

  z1 = sqrt ( x0 * x0 + y * y )
  th = atan ( y / x0 )
  gr = ( x0 - 0.5D+00 ) * log ( z1 ) - th * y - x0 &
       + 0.5D+00 * log ( 2.0D+00 * pi )
  gi = th * ( x0 - 0.5D+00 ) + y * log ( z1 ) - y

  do k = 1, 10
     t = z1 ** ( 1 - 2 * k )
     gr = gr + a(k) * t * cos ( ( 2.0D+00 * k - 1.0D+00 ) * th )
     gi = gi - a(k) * t * sin ( ( 2.0D+00 * k - 1.0D+00 ) * th )
  end do

  if ( x <= 7.0D+00 ) then
     gr1 = 0.0D+00
     gi1 = 0.0D+00
     do j = 0, na - 1
        gr1 = gr1 + 0.5D+00 * log ( ( x + j ) ** 2 + y * y )
        gi1 = gi1 + atan ( y / ( x + j ) )
     end do
     gr = gr - gr1
     gi = gi - gi1
  end if

  if ( x1 < 0.0D+00 ) then
     z1 = sqrt ( x * x + y * y )
     th1 = atan ( y / x )
     sr = - sin ( pi * x ) * cosh ( pi * y )
     si = - cos ( pi * x ) * sinh ( pi * y )
     z2 = sqrt ( sr * sr + si * si )
     th2 = atan ( si / sr )
     if ( sr < 0.0D+00 ) then
        th2 = pi + th2
     end if
     gr = log ( pi / ( z1 * z2 ) ) - gr
     gi = - th1 - th2 - gi
     x = x1
     y = y1
  end if

  if ( kf == 1 ) then
     g0 = exp ( gr )
     gr = g0 * cos ( gi )
     gi = g0 * sin ( gi )
  end if

  return
end subroutine cgama




  !*****************************************************************************80
  !
  !! CH12N computes Hankel functions of first and second kinds, complex argument.
  !
  !  Discussion:
  !
  !    Both the Hankel functions and their derivatives are computed.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    26 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) N, the order of the functions.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n), 
  !    the values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
  !
subroutine ch12n ( n, z, nm, chf1, chd1, chf2, chd2 )
  !CH12N computes Hankel functions of first and second kinds, complex argument.
  !
  !Both the Hankel functions and their derivatives are computed.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of the functions
  !@param z: in: the argument
  !@param nm: out: the highest order computed
  !@param chf1: out: Hn(1)(z)
  !@param chd1: out: Hn(1)'(z)
  !@param chf2: out: Hn(2)(z)
  !@param chd2: out: Hn(2)'(z)
  !
  implicit none

  integer(4) :: n

  complex(8) :: cbi(0:250)
  complex(8) :: cbj(0:250)
  complex(8) :: cbk(0:250)
  complex(8) :: cby(0:250)
  complex(8) :: cdi(0:250)
  complex(8) :: cdj(0:250)
  complex(8) :: cdk(0:250)
  complex(8) :: cdy(0:250)
  complex(8) :: chd1(0:n)
  complex(8) :: chd2(0:n)
  complex(8) :: cf1
  complex(8) :: cfac
  complex(8) :: chf1(0:n)
  complex(8) :: chf2(0:n)
  complex(8) :: ci
  integer(4) :: k
  integer(4) :: nm
  real(8)    :: pi
  complex(8) :: z
  complex(8) :: zi

  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  pi = 3.141592653589793D+00

  if ( imag ( z ) < 0.0D+00 ) then

     call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

     do k = 0, nm
        chf1(k) = cbj(k) + ci * cby(k)
        chd1(k) = cdj(k) + ci * cdy(k)
     end do

     zi = ci * z
     call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
     cfac = -2.0D+00 / ( pi * ci )

     do k = 0, nm
        chf2(k) = cfac * cbk(k)
        chd2(k) = cfac * ci * cdk(k)
        cfac = cfac * ci
     end do

  else if ( 0.0D+00 < imag ( z ) ) then

     zi = - ci * z
     call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
     cf1 = -ci
     cfac = 2.0D+00 / ( pi * ci )

     do k = 0, nm
        chf1(k) = cfac * cbk(k)
        chd1(k) = -cfac * ci * cdk(k)
        cfac = cfac * cf1
     end do

     call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

     do k = 0, nm
        chf2(k) = cbj(k) - ci * cby(k)
        chd2(k) = cdj(k) - ci * cdy(k)
     end do

  else

     call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

     do k = 0, nm
        chf1(k) = cbj(k) + ci * cby(k)
        chd1(k) = cdj(k) + ci * cdy(k)
        chf2(k) = cbj(k) - ci * cby(k)
        chd2(k) = cdj(k) - ci * cdy(k)
     end do

  end if

  return
end subroutine ch12n




  !*****************************************************************************80
  !
  !! CHGM computes the confluent hypergeometric function M(a,b,x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    27 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HG, the value of M(a,b,x).
  !
subroutine chgm ( a, b, x, hg )
  !CHGM computes the confluent hypergeometric function M(a,b,x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hg: out: the value of M(a,b,x)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: aa
  real(8)    :: b
  real(8)    :: hg
  real(8)    :: hg1
  real(8)    :: hg2
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: la
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nl
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: rg
  real(8)    :: sum1
  real(8)    :: sum2
  real(8)    :: ta
  real(8)    :: tb
  real(8)    :: tba
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: xg
  real(8)    :: y0
  real(8)    :: y1

  pi = 3.141592653589793D+00
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0D+00

  if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
     hg = 1.0D+300
  else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
     hg = 1.0D+00
  else if ( a == -1.0D+00 ) then
     hg = 1.0D+00 - x / b
  else if ( a == b ) then
     hg = exp ( x )
  else if ( a - b == 1.0D+00 ) then
     hg = ( 1.0D+00 + x / b ) * exp ( x )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
     hg = ( exp ( x ) - 1.0D+00 ) / x
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
     m = int ( - a )
     r = 1.0D+00
     hg = 1.0D+00
     do k = 1, m
        r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
        hg = hg + r
     end do
  end if

  if ( hg /= 0.0D+00 ) then
     return
  end if

  if ( x < 0.0D+00 ) then
     a = b - a
     a0 = a
     x = abs ( x )
  end if

  if ( a < 2.0D+00 ) then
     nl = 0
  end if

  if ( 2.0D+00 <= a ) then
     nl = 1
     la = int ( a )
     a = a - la - 1.0D+00
  end if

  do n = 0, nl

     if ( 2.0D+00 <= a0 ) then
        a = a + 1.0D+00
     end if

     if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

        hg = 1.0D+00
        rg = 1.0D+00
        do j = 1, 500
           rg = rg * ( a + j - 1.0D+00 ) &
                / ( j * ( b + j - 1.0D+00 ) ) * x
           hg = hg + rg
           if ( abs ( rg / hg ) < 1.0D-15 ) then
              exit
           end if
        end do

     else

        call gammaf ( a, ta )
        call gammaf ( b, tb )
        xg = b - a
        call gammaf ( xg, tba )
        sum1 = 1.0D+00
        sum2 = 1.0D+00
        r1 = 1.0D+00
        r2 = 1.0D+00
        do i = 1, 8
           r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
           r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
           sum1 = sum1 + r1
           sum2 = sum2 + r2
        end do
        hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
        hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
        hg = hg1 + hg2

     end if

     if ( n == 0 ) then
        y0 = hg
     else if ( n == 1 ) then
        y1 = hg
     end if

  end do

  if ( 2.0D+00 <= a0 ) then
     do i = 1, la - 1
        hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
        y0 = y1
        y1 = hg
        a = a + 1.0D+00
     end do
  end if

  if ( x0 < 0.0D+00 ) then
     hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end subroutine chgm




  !*****************************************************************************80
  !
  !! CHGU computes the confluent hypergeometric function U(a,b,x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    27 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HU, U(a,b,x).
  !
  !    Output, integer(4) MD, the method code.
  !
subroutine chgu ( a, b, x, hu, md )
  !CHGU computes the confluent hypergeometric function U(a,b,x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hu: out: U(a,b,x)
  !@param md: out: the method code
  !
  implicit none

  real(8)    :: a
  real(8)    :: a00
  real(8)    :: aa
  real(8)    :: b
  real(8)    :: b00
  logical    :: bl1
  logical    :: bl2
  logical    :: bl3
  logical    :: bn
  real(8)    :: hu
  real(8)    :: hu1
  integer(4) :: id
  integer(4) :: id1
  logical    :: il1
  logical    :: il2
  logical    :: il3
  integer(4) :: md
  real(8)    :: x

  aa = a - b + 1.0D+00
  il1 = a == int ( a ) .and. a <= 0.0D+00
  il2 = aa == int ( aa ) .and. aa <= 0.0D+00
  il3 = abs ( a * ( a - b + 1.0D+00 ) ) / x <= 2.0D+00
  bl1 = x <= 5.0D+00 .or. ( x <= 10.0D+00 .and. a <= 2.0D+00 )
  bl2 = ( 5.0D+00 < x .and. x <= 12.5D+00 ) .and. &
       ( 1.0D+00 <= a .and. a + 4.0D+00 <= b )
  bl3 = 12.5D+00 < x .and. 5.0D+00 <= a .and. a + 5.0D+00 <= b
  bn = b == int ( b ) .and. b .ne. 0.0D+00
  id1 = -100

  if ( b .ne. int ( b ) ) then
     call chgus ( a, b, x, hu, id1 )
     md = 1
     if ( 6 <= id1 ) then
        return
     end if
     hu1 = hu
  end if

  if ( il1 .or. il2 .or. il3 ) then
     call chgul ( a, b, x, hu, id )
     md = 2
     if ( 6 <= id ) then
        return
     end if
     if ( id < id1 ) then
        md = 1
        id = id1
        hu = hu1
     end if
  end if

  if ( 0.0D+00 <= a ) then
     if ( bn .and. ( bl1 .or. bl2 .or. bl3 ) ) then
        call chgubi ( a, b, x, hu, id )
        md = 3
     else
        call chguit ( a, b, x, hu, id )
        md = 4
     end if
  else
     if ( b <= a ) then
        a00 = a
        b00 = b
        a = a - b + 1.0D+00
        b = 2.0D+00 - b
        call chguit ( a, b, x, hu, id )
        hu = x ** ( 1.0D+00 - b00 ) * hu
        a = a00
        b = b00
        md = 4
     else if ( bn .and. ( .not. il1 ) ) then
        call chgubi ( a, b, x, hu, id )
        md = 3
     end if
  end if

  if ( id < 6 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'CHGU - Warning!'
     write ( *, '(a)' ) '  Accurate results were not obtained.'
  end if

  return
end subroutine chgu




  !*****************************************************************************80
  !
  !! CHGUBI: confluent hypergeometric function with integer argument B.
  !
  !  Discussion:
  !
  !    This procedure computes the confluent hypergeometric function
  !    U(a,b,x) with integer(4) b ( b = ±1,±2,... )
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HU, the value of U(a,b,x).
  !
  !    Output, integer(4) ID, the estimated number of significant
  !    digits.
  !
subroutine chgubi ( a, b, x, hu, id )
  !CHGUBI: confluent hypergeometric function with integer argument B.
  !
  !This procedure computes the confluent hypergeometric function
  !U(a,b,x) with integer(4) b ( b = +/-1,+/-2,... )
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hu: out: the value of U(a,b,x)
  !@param id: out: the estimated number of significant digits
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: a2
  real(8)    :: b
  real(8)    :: da1
  real(8)    :: da2
  real(8)    :: db1
  real(8)    :: db2
  real(8)    :: el
  real(8)    :: ga
  real(8)    :: ga1
  real(8)    :: h0
  real(8)    :: hm1
  real(8)    :: hm2
  real(8)    :: hm3
  real(8)    :: hmax
  real(8)    :: hmin
  real(8)    :: hu
  real(8)    :: hu1
  real(8)    :: hu2
  real(8)    :: hw
  integer(4) :: id
  integer(4) :: id1
  integer(4) :: id2
  integer(4) :: j
  integer(4) :: k
  integer(4) :: m
  integer(4) :: n
  real(8)    :: ps
  real(8)    :: r
  real(8)    :: rn
  real(8)    :: rn1
  real(8)    :: s0
  real(8)    :: s1
  real(8)    :: s2
  real(8)    :: sa
  real(8)    :: sb
  real(8)    :: ua
  real(8)    :: ub
  real(8)    :: x

  id = -100
  el = 0.5772156649015329D+00
  n = abs ( b - 1 )
  rn1 = 1.0D+00
  rn = 1.0D+00
  do j = 1, n
     rn = rn * j
     if ( j == n - 1 ) then
        rn1 = rn
     end if
  end do

  call psi ( a, ps )
  call gammaf ( a, ga )

  if ( 0.0D+00 < b ) then
     a0 = a
     a1 = a - n
     a2 = a1
     call gammaf ( a1, ga1 )
     ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga1 )
     ub = rn1 / ga * x ** ( - n )
  else
     a0 = a + n
     a1 = a0
     a2 = a
     call gammaf ( a1, ga1 )
     ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga ) * x ** n
     ub = rn1 / ga1
  end if

  hm1 = 1.0D+00
  r = 1.0D+00
  hmax = 0.0D+00
  hmin = 1.0D+300

  do k = 1, 150
     r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
     hm1 = hm1 + r
     hu1 = abs ( hm1 )
     hmax = max ( hmax, hu1 )
     hmin = min ( hmin, hu1 )
     if ( abs ( hm1 - h0 ) < abs ( hm1 ) * 1.0D-15 ) then
        exit
     end if
     h0 = hm1
  end do

  da1 = log10 ( hmax )
  if ( hmin /= 0.0D+00 ) then
     da2 = log10 ( hmin )
  end if
  id = 15 - abs ( da1 - da2 )
  hm1 = hm1 * log ( x )
  s0 = 0.0D+00
  do m = 1, n
     if ( 0.0D+00 <= b ) then
        s0 = s0 - 1.0D+00 / m
     else
        s0 = s0 + ( 1.0D+00 - a ) / ( m * ( a + m - 1.0D+00 ) )
     end if
  end do
  hm2 = ps + 2.0D+00 * el + s0
  r = 1.0D+00
  hmax = 0.0D+00
  hmin = 1.0D+300
  do k = 1, 150
     s1 = 0.0D+00
     s2 = 0.0D+00
     if ( 0.0D+00 < b ) then
        do m = 1, k
           s1 = s1 - ( m + 2.0D+00 * a - 2.0D+00 ) / ( m * ( m + a - 1.0D+00 ) )
        end do
        do m = 1, n
           s2 = s2 + 1.0D+00 / ( k + m )
        end do
     else
        do m = 1, k + n
           s1 = s1 + ( 1.0D+00 - a ) / ( m * ( m + a - 1.0D+00 ) )
        end do
        do m = 1, k
           s2 = s2 + 1.0D+00 / m
        end do
     end if
     hw = 2.0D+00 * el + ps + s1 - s2
     r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
     hm2 = hm2 + r * hw
     hu2 = abs ( hm2 )
     hmax = max ( hmax, hu2 )
     hmin = min ( hmin, hu2 )

     if ( abs ( ( hm2 - h0 ) / hm2 ) < 1.0D-15 ) then
        exit
     end if

     h0 = hm2

  end do

  db1 = log10 ( hmax )
  if ( hmin /= 0.0D+00 ) then
     db2 = log10 ( hmin )
  end if
  id1 = 15 - abs ( db1 - db2 )
  id = min ( id, id1 )

  if ( n == 0 ) then
     hm3 = 0.0D+00
  else
     hm3 = 1.0D+00
  end if

  r = 1.0D+00
  do k = 1, n - 1
     r = r * ( a2 + k - 1.0D+00 ) / ( ( k - n ) * k ) * x
     hm3 = hm3 + r
  end do

  sa = ua * ( hm1 + hm2 )
  sb = ub * hm3
  hu = sa + sb

  if ( sa /= 0.0D+00 ) then
     id1 = int ( log10 ( abs ( sa ) ) )
  end if

  if ( hu /= 0.0D+00 ) then
     id2 = int ( log10 ( abs ( hu ) ) )
  end if

  if ( sa * sb < 0.0D+00 ) then
     id = id - abs ( id1 - id2 )
  end if

  return
end subroutine chgubi




  !*****************************************************************************80
  !
  !! CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
  !
  !  Discussion:
  !
  !    This procedure computes the hypergeometric function U(a,b,x) by
  !    using Gaussian-Legendre integration (n = 60)
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, double precision A, B, parameters.
  !
  !    Input, double precision X, the argument.
  !
  !    Output, double precision HU, U(a,b,z).
  !
  !    Output, integer ID, the estimated number of significant digits.
  !
subroutine chguit ( a, b, x, hu, id )
  !CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
  !
  !This procedure computes the hypergeometric function U(a,b,x) by
  !using Gaussian-Legendre integration (n = 60)
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hu: out: U(a,b,z)
  !@param id: out: the estimated number of significant digits
  !
  implicit none

  double precision :: a
  double precision :: a1
  double precision :: b
  double precision :: b1
  double precision :: c
  double precision :: d
  double precision :: f1
  double precision :: f2
  double precision :: g
  double precision :: ga
  double precision :: hu
  double precision :: hu0
  double precision :: hu1
  double precision :: hu2
  integer          :: id
  integer          :: j
  integer          :: k
  integer          :: m
  double precision :: s
  double precision, save, dimension ( 30 ) :: t = (/ &
       0.259597723012478D-01, 0.778093339495366D-01, &
       0.129449135396945D+00, 0.180739964873425D+00, &
       0.231543551376029D+00, 0.281722937423262D+00, &
       0.331142848268448D+00, 0.379670056576798D+00, &
       0.427173741583078D+00, 0.473525841761707D+00, &
       0.518601400058570D+00, 0.562278900753945D+00, &
       0.604440597048510D+00, 0.644972828489477D+00, &
       0.683766327381356D+00, 0.720716513355730D+00, &
       0.755723775306586D+00, 0.788693739932264D+00, &
       0.819537526162146D+00, 0.848171984785930D+00, &
       0.874519922646898D+00, 0.898510310810046D+00, &
       0.920078476177628D+00, 0.939166276116423D+00, &
       0.955722255839996D+00, 0.969701788765053D+00, &
       0.981067201752598D+00, 0.989787895222222D+00, &
       0.995840525118838D+00, 0.999210123227436D+00 /)
  double precision :: t1
  double precision :: t2
  double precision :: t3
  double precision :: t4
  double precision, save, dimension ( 30 ) :: w = (/ &
       0.519078776312206D-01, 0.517679431749102D-01, &
       0.514884515009810D-01, 0.510701560698557D-01, &
       0.505141845325094D-01, 0.498220356905502D-01, &
       0.489955754557568D-01, 0.480370318199712D-01, &
       0.469489888489122D-01, 0.457343797161145D-01, &
       0.443964787957872D-01, 0.429388928359356D-01, &
       0.413655512355848D-01, 0.396806954523808D-01, &
       0.378888675692434D-01, 0.359948980510845D-01, &
       0.340038927249464D-01, 0.319212190192963D-01, &
       0.297524915007890D-01, 0.275035567499248D-01, &
       0.251804776215213D-01, 0.227895169439978D-01, &
       0.203371207294572D-01, 0.178299010142074D-01, &
       0.152746185967848D-01, 0.126781664768159D-01, &
       0.100475571822880D-01, 0.738993116334531D-02, &
       0.471272992695363D-02, 0.202681196887362D-02 /)
  double precision :: x

  id = 7
  a1 = a - 1.0D+00
  b1 = b - a - 1.0D+00
  c = 12.0D+00 / x

  do m = 10, 100, 5

     hu1 = 0.0D+00
     g = 0.5D+00 * c / m
     d = g
     do j = 1, m
        s = 0.0D+00
        do k = 1, 30
           t1 = d + g * t(k)
           t2 = d - g * t(k)
           f1 = exp ( - x * t1 ) * t1 ** a1 * ( 1.0D+00 + t1 ) ** b1
           f2 = exp ( - x * t2 ) * t2 ** a1 * ( 1.0D+00 + t2 ) ** b1
           s = s + w(k) * ( f1 + f2 )
        end do
        hu1 = hu1 + s * g
        d = d + 2.0D+00 * g
     end do

     if ( abs ( 1.0D+00 - hu0 / hu1 ) < 1.0D-07 ) then
        exit
     end if

     hu0 = hu1

  end do

  call gammaf ( a, ga )
  hu1 = hu1 / ga

  do m = 2, 10, 2
     hu2 = 0.0D+00
     g = 0.5D+00 / m
     d = g
     do j = 1, m
        s = 0.0D+00
        do k = 1, 30
           t1 = d + g * t(k)
           t2 = d - g * t(k)
           t3 = c / ( 1.0D+00 - t1 )
           t4 = c / ( 1.0D+00 - t2 ) 
           f1 = t3 * t3 / c * exp ( - x * t3 ) * t3 ** a1 * ( 1.0D+00 + t3 ) ** b1
           f2 = t4 * t4 / c * exp ( - x * t4 ) * t4 ** a1 * ( 1.0D+00 + t4 ) ** b1
           s = s + w(k) * ( f1 + f2 )
        end do
        hu2 = hu2 + s * g
        d = d + 2.0D+00 * g
     end do

     if ( abs ( 1.0D+00 - hu0 / hu2 ) < 1.0D-07 ) then
        exit
     end if

     hu0 = hu2

  end do

  call gammaf ( a, ga )
  hu2 = hu2 / ga
  hu = hu1 + hu2

  return
end subroutine chguit




  !*****************************************************************************80
  !
  !! CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HU, the value of U(a,b,x).
  !
  !    Output, integer(4) ID, the estimated number of 
  !    significant digits.
  !
subroutine chgul ( a, b, x, hu, id )
  !CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hu: out: the value of U(a,b,x)
  !@param id: out: the estimated number of significant digits
  !
  implicit none

  real(8)    :: a
  real(8)    :: aa
  real(8)    :: b
  real(8)    :: hu
  integer(4) :: id
  logical    :: il1
  logical    :: il2
  integer(4) :: k
  integer(4) :: nm
  real(8)    :: r
  real(8)    :: ra
  real(8)    :: r0
  real(8)    :: x

  id = -100
  aa = a - b + 1.0D+00
  il1 = ( a == int ( a ) ) .and. ( a <= 0.0D+00 )
  il2 = ( aa == int ( aa ) ) .and. ( aa <= 0.0D+00 )

  if ( il1 .or. il2 ) then

     if ( il1 ) then
        nm = abs ( a )
     end if

     if ( il2 ) then
        nm = abs ( aa )
     end if

     hu = 1.0D+00
     r = 1.0D+00
     do k = 1, nm
        r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
        hu = hu + r
     end do
     hu = x ** ( - a ) * hu
     id = 10

  else

     hu = 1.0D+00
     r = 1.0D+00
     do k = 1, 25
        r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
        ra = abs ( r )
        if ( ( 5 < k .and. r0 <= ra ) .or. ra < 1.0D-15 ) then
           exit
        end if
        r0 = ra
        hu = hu + r
     end do

     id = abs ( log10 ( ra ) )
     hu = x ** ( - a ) * hu

  end if

  return
end subroutine chgul




  !*****************************************************************************80
  !
  !! CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    27 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HU, U(a,b,x).
  !
  !    Output, integer(4) ID, the estimated number of 
  !    significant digits.
  !
subroutine chgus ( a, b, x, hu, id )
  !CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param hu: out: U(a,b,x)
  !@param id: out: the estimated number of significant digits
  !
  implicit none

  real(8)    :: a
  real(8)    :: b
  real(8)    :: d1
  real(8)    :: d2
  real(8)    :: ga
  real(8)    :: gab
  real(8)    :: gb
  real(8)    :: gb2
  real(8)    :: h0
  real(8)    :: hmax
  real(8)    :: hmin
  real(8)    :: hu
  real(8)    :: hu0
  real(8)    :: hua
  integer(4) :: id
  integer(4) :: j
  real(8)    :: pi
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: x
  real(8)    :: xg1
  real(8)    :: xg2

  id = -100
  pi = 3.141592653589793D+00
  call gammaf ( a, ga )
  call gammaf ( b, gb )
  xg1 = 1.0D+00 + a - b
  call gammaf ( xg1, gab )
  xg2 = 2.0D+00 - b
  call gammaf ( xg2, gb2 )
  hu0 = pi / sin ( pi * b )
  r1 = hu0 / ( gab * gb )
  r2 = hu0 * x ** ( 1.0D+00 - b ) / ( ga * gb2 )
  hu = r1 - r2
  hmax = 0.0D+00
  hmin = 1.0D+300
  do j = 1, 150
     r1 = r1 * ( a + j - 1.0D+00 ) / ( j * ( b + j - 1.0D+00 ) ) * x
     r2 = r2 * ( a - b + j ) / ( j * ( 1.0D+00 - b + j ) ) * x
     hu = hu + r1 - r2
     hua = abs ( hu )
     hmax = max ( hmax, hua )
     hmin = min ( hmin, hua )
     if ( abs ( hu - h0 ) < abs ( hu ) * 1.0D-15 ) then
        exit
     end if
     h0 = hu
  end do

  d1 = log10 ( hmax )
  if ( hmin /= 0.0D+00 ) then
     d2 = log10 ( hmin )
  end if
  id = 15 - abs ( d1 - d2 )

  return
end subroutine chgus




  !*****************************************************************************80
  !
  !! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
  !
  !  Discussion:
  !
  !    This procedure computes the modified Bessel functions I0(z), I1(z), 
  !    K0(z), K1(z), and their derivatives for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, 
  !    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), 
  !    and K1'(z).
  !
subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )
  !CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
  !
  !This procedure computes the modified Bessel functions I0(z), I1(z),
  !K0(z), K1(z), and their derivatives for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param cbi0: out: I0(z)
  !@param cdi0: out: I0'(z)
  !@param cbi1: out: I1(z)
  !@param cdi1: out: I1'(z)
  !@param cbk0: out: K0(z)
  !@param cdk0: out: K0'(z)
  !@param cbk1: out: K1(z)
  !@param cdk1: and K1'(z)
  !
  implicit none

  real(8), save, dimension ( 12 ) :: a = (/ &
       0.125D+00,           7.03125D-02,&
       7.32421875D-02,      1.1215209960938D-01,&
       2.2710800170898D-01, 5.7250142097473D-01,&
       1.7277275025845D+00, 6.0740420012735D+00,&
       2.4380529699556D+01, 1.1001714026925D+02,&
       5.5133589612202D+02, 3.0380905109224D+03 /)
  real(8)    :: a0
  real(8), save, dimension ( 10 ) :: a1 = (/ &
       0.125D+00,            0.2109375D+00, &
       1.0986328125D+00,     1.1775970458984D+01, &
       2.1461706161499D+002, 5.9511522710323D+03, &
       2.3347645606175D+05,  1.2312234987631D+07, &
       8.401390346421D+08,   7.2031420482627D+10 /)
  real(8), save, dimension ( 12 ) :: b = (/ &
       -0.375D+00,           -1.171875D-01, &
       -1.025390625D-01,     -1.4419555664063D-01, &
       -2.7757644653320D-01, -6.7659258842468D-01, &
       -1.9935317337513D+00, -6.8839142681099D+00, &
       -2.7248827311269D+01, -1.2159789187654D+02, &
       -6.0384407670507D+02, -3.3022722944809D+03 /)
  complex(8) :: ca
  complex(8) :: cb
  complex(8) :: cbi0
  complex(8) :: cbi1
  complex(8) :: cbk0
  complex(8) :: cbk1
  complex(8) :: cdi0
  complex(8) :: cdi1
  complex(8) :: cdk0
  complex(8) :: cdk1
  complex(8) :: ci
  complex(8) :: cr
  complex(8) :: cs
  complex(8) :: ct
  complex(8) :: cw
  integer(4) :: k
  integer(4) :: k0
  real(8)    :: pi
  real(8)    :: w0
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2
  complex(8) :: zr
  complex(8) :: zr2

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 .eq. 0.0D+00 ) then
     cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cbi1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cdi0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cdi1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     cbk0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cbk1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cdk0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cdk1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 <= 18.0D+00 ) then

     cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 50
        cr = 0.25D+00 * cr * z2 / ( k * k )
        cbi0 = cbi0 + cr
        if ( abs ( cr / cbi0 ) < 1.0D-15 ) then
           exit
        end if
     end do

     cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 50
        cr = 0.25D+00 * cr * z2 / ( k * ( k + 1 ) )
        cbi1 = cbi1 + cr
        if ( abs ( cr / cbi1 ) < 1.0D-15 ) then
           exit
        end if
     end do

     cbi1 = 0.5D+00 * z1 * cbi1

  else

     if ( a0 < 35.0D+00 ) then
        k0 = 12
     else if ( a0 < 50.0D+00 ) then
        k0 = 9
     else
        k0 = 7
     end if

     ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
     cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     zr = 1.0D+00 / z1
     do k = 1, k0
        cbi0 = cbi0 + a(k) * zr ** k
     end do
     cbi0 = ca * cbi0
     cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cbi1 = cbi1 + b(k) * zr ** k
     end do
     cbi1 = ca * cbi1

  end if

  if ( a0 <= 9.0D+00 ) then

     cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
     w0 = 0.0D+00
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 50
        w0 = w0 + 1.0D+00 / k
        cr = 0.25D+00 * cr / ( k * k ) * z2
        cs = cs + cr * ( w0 + ct )
        if ( abs ( ( cs - cw ) / cs ) < 1.0D-15 ) then
           exit
        end if
        cw = cs
     end do

     cbk0 = ct + cs

  else

     cb = 0.5D+00 / z1
     zr2 = 1.0D+00 / z2
     cbk0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 10
        cbk0 = cbk0 + a1(k) * zr2 ** k
     end do
     cbk0 = cb * cbk0 / cbi0

  end if

  cbk1 = ( 1.0D+00 / z1 - cbi1 * cbk0 ) / cbi0

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then

     if ( imag ( z ) < 0.0D+00 ) then
        cbk0 = cbk0 + ci * pi * cbi0
        cbk1 = - cbk1 + ci * pi * cbi1
     else
        cbk0 = cbk0 - ci * pi * cbi0
        cbk1 = - cbk1 - ci * pi * cbi1
     end if

     cbi1 = - cbi1

  end if

  cdi0 = cbi1
  cdi1 = cbi0 - 1.0D+00 / z * cbi1
  cdk0 = - cbk1
  cdk1 = - cbk0 - 1.0D+00 / z * cbk1

  return
end subroutine cik01




  !*****************************************************************************80
  !
  !! CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
  !
  !  Discussion:
  !
  !    This procedure computes modified Bessel functions Iv(z) and
  !    Kv(z) and their derivatives with a complex argument and a large order.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Iv(z) and Kv(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, real(8) CBIV, CDIV, CBKV, CDKV, the values of
  !    Iv(z), Iv'(z), Kv(z), Kv'(z).
  !
subroutine ciklv ( v, z, cbiv, cdiv, cbkv, cdkv )
  !CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
  !
  !This procedure computes modified Bessel functions Iv(z) and
  !Kv(z) and their derivatives with a complex argument and a large order.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Iv(z) and Kv(z)
  !@param z: in: the argument
  !@param cbiv: out: Iv(z)
  !@param cdiv: out: Iv'(z)
  !@param cbkv: out: Kv(z)
  !@param cdkv: out: Kv'(z)
  !
  implicit none

  real(8)    :: a(91)
  complex(8) :: cbiv
  complex(8) :: cbkv
  complex(8) :: cdiv
  complex(8) :: cdkv
  complex(8) :: ceta
  complex(8) :: cf(12)
  complex(8) :: cfi
  complex(8) :: cfk
  complex(8) :: csi
  complex(8) :: csk
  complex(8) :: ct
  complex(8) :: ct2
  complex(8) :: cws
  integer(4) :: i
  integer(4) :: k
  integer(4) :: km
  integer(4) :: l
  integer(4) :: l0
  integer(4) :: lf
  real(8)    :: pi
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vr
  complex(8) :: z

  pi = 3.141592653589793D+00
  km = 12
  call cjk ( km, a )

  do l = 1, 0, -1

     v0 = v - l
     cws = sqrt ( 1.0D+00 + ( z / v0 ) * ( z / v0 ) )
     ceta = cws + log ( z / v0 / ( 1.0D+00 + cws ) )
     ct = 1.0D+00 / cws
     ct2 = ct * ct
     do k = 1, km
        l0 = k * ( k + 1 ) / 2 + 1
        lf = l0 + k
        cf(k) = a(lf)
        do i = lf - 1, l0, -1
           cf(k) = cf(k) * ct2 + a(i)
        end do
        cf(k) = cf(k) * ct ** k
     end do
     vr = 1.0D+00 / v0
     csi = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, km
        csi = csi + cf(k) * vr ** k
     end do
     cbiv = sqrt ( ct / ( 2.0D+00 * pi * v0 ) ) * exp ( v0 * ceta ) * csi
     if ( l == 1 ) then
        cfi = cbiv
     end if
     csk = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, km
        csk = csk + ( - 1 ) ** k * cf(k) * vr ** k
     end do
     cbkv = sqrt ( pi * ct / ( 2.0D+00 * v0 ) ) * exp ( - v0 * ceta ) * csk

     if ( l == 1 ) then
        cfk = cbkv
     end if

  end do

  cdiv =   cfi - v / z * cbiv
  cdkv = - cfk - v / z * cbkv

  return
end subroutine ciklv




  !*****************************************************************************80
  !
  !! CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    30 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of In(z) and Kn(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CBI((0:N), CDI(0:N), CBK(0:N), CDK(0:N), 
  !    the values of In(z), In'(z), Kn(z), Kn'(z).
  !
subroutine cikna ( n, z, nm, cbi, cdi, cbk, cdk )
  !CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of In(z) and Kn(z)
  !@param z: in: the argument
  !@param nm: out: the highest order computed
  !@param cbi: In(z), n=0,...,N
  !@param cdi: In'(z), n=0,...,N
  !@param cbk: Kn(z), n=0,...,N
  !@param cdk: Kn'(z), n=0,...,N
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  complex(8) :: cbi(0:n)
  complex(8) :: cbi0
  complex(8) :: cbi1
  complex(8) :: cbk(0:n)
  complex(8) :: cbk0
  complex(8) :: cbk1
  complex(8) :: cdi(0:n)
  complex(8) :: cdi0
  complex(8) :: cdi1
  complex(8) :: cdk(0:n)
  complex(8) :: cdk0
  complex(8) :: cdk1
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: ckk
  complex(8) :: cs
  integer(4) :: k
  integer(4) :: m
  ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  complex(8) :: z

  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-100 ) then
     do k = 0, n
        cbi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cbk(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdk(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do
     cbi(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cdi(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     return
  end if

  call cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

  cbi(0) = cbi0
  cbi(1) = cbi1
  cbk(0) = cbk0
  cbk(1) = cbk1
  cdi(0) = cdi0
  cdi(1) = cdi1
  cdk(0) = cdk0
  cdk(1) = cdk1

  if ( n <= 1 ) then
     return
  end if

  m = msta1 ( a0, 200 )

  if ( m < n ) then
     nm = m
  else
     m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
  do k = m, 0, -1
     cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 + cf2
     if ( k <= nm ) then
        cbi(k) = cf
     end if
     cf2 = cf1
     cf1 = cf
  end do

  cs = cbi0 / cf
  do k = 0, nm
     cbi(k) = cs * cbi(k)
  end do

  do k = 2, nm
     if ( abs ( cbi(k-2) ) < abs ( cbi(k-1) ) ) then
        ckk = ( 1.0D+00 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
     else
        ckk = ( cbi(k) * cbk(k-2) + 2.0D+00 * ( k - 1.0D+00 ) &
             / ( z * z ) ) / cbi(k-2)
     end if
     cbk(k) = ckk
  end do

  do k = 2, nm
     cdi(k) =   cbi(k-1) - k / z * cbi(k)
     cdk(k) = - cbk(k-1) - k / z * cbk(k)
  end do

  return
end subroutine cikna




  !*****************************************************************************80
  !
  !! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
  !
  !  Discussion:
  !
  !    This procedure also evaluates the derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    30 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of In(z) and Kn(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CB((0:N), CDI(0:N), CBK(0:N), CDK(0:N), 
  !    the values of In(z), In'(z), Kn(z), Kn'(z).
  !
subroutine ciknb ( n, z, nm, cbi, cdi, cbk, cdk )
  !CIKNB computes complex modified Bessel functions In(z) and Kn(z).
  !
  !This procedure also evaluates the derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of In(z) and Kn(z)
  !@param z: in: the argument
  !@param nm: out: the highest order computed
  !@param cbi: In(z), n=0,...,N
  !@param cdi: In'(z), n=0,...,N
  !@param cbk: Kn(z), n=0,...,N
  !@param cdk: Kn'(z), n=0,...,N
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  complex(8) :: c
  complex(8) :: ca0
  complex(8) :: cbi(0:n)
  complex(8) :: cbkl
  complex(8) :: cbs
  complex(8) :: cdi(0:n)
  complex(8) :: cbk(0:n)
  complex(8) :: cdk(0:n)
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cg
  complex(8) :: cg0
  complex(8) :: cg1
  complex(8) :: ci
  complex(8) :: cr
  complex(8) :: cs0
  complex(8) :: csk0
  real(8)    :: el
  real(8)    :: fac
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: l
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: pi
  real(8)    :: vt
  complex(8) :: z
  complex(8) :: z1

  pi = 3.141592653589793D+00
  el = 0.57721566490153D+00
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-100 ) then
     do k = 0, n
        cbi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cbk(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdk(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do
     cbi(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cdi(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 ) 
     return
  end if

  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  else
     z1 = z
  end if

  if ( n == 0 ) then
     nm = 1
  end if

  m = msta1 ( a0, 200 )

  if ( m < nm ) then
     nm = m
  else
     m = msta2 ( a0, nm, 15 )
  end if

  cbs = 0.0D+00
  csk0 = 0.0D+00
  cf0 = 0.0D+00
  cf1 = 1.0D-100

  do k = m, 0, -1
     cf = 2.0D+00 * ( k + 1.0D+00 ) * cf1 / z1 + cf0
     if ( k <= nm ) then
        cbi(k) = cf
     end if
     if ( k /= 0 .and. k == 2 * int ( k / 2 ) ) then
        csk0 = csk0 + 4.0D+00 * cf / k
     end if
     cbs = cbs + 2.0D+00 * cf
     cf0 = cf1
     cf1 = cf
  end do

  cs0 = exp ( z1 ) / ( cbs - cf )

  do k = 0, nm
     cbi(k) = cs0 * cbi(k)
  end do

  if ( a0 <= 9.0D+00 ) then

     cbk(0) = - ( log ( 0.5D+00 * z1 ) + el ) * cbi(0) + cs0 * csk0
     cbk(1) = ( 1.0D+00 / z1 - cbi(1) * cbk(0) ) / cbi(0)

  else

     ca0 = sqrt ( pi / ( 2.0D+00 * z1 ) ) * exp ( -z1 )

     if ( a0 < 25.0D+00 ) then
        k0 = 16
     else if ( a0 < 80.0D+00 ) then
        k0 = 10
     else if ( a0 < 200.0D+00 ) then
        k0 = 8
     else
        k0 = 6
     end if

     do l = 0, 1
        cbkl = 1.0D+00
        vt = 4.0D+00 * l
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, k0
           cr = 0.125D+00 * cr &
                * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
           cbkl = cbkl + cr
        end do
        cbk(l) = ca0 * cbkl
     end do
  end if

  cg0 = cbk(0)
  cg1 = cbk(1)
  do k = 2, nm
     cg = 2.0D+00 * ( k - 1.0D+00 ) / z1 * cg1 + cg0
     cbk(k) = cg
     cg0 = cg1
     cg1 = cg
  end do

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     fac = 1.0D+00
     do k = 0, nm
        if ( imag ( z ) < 0.0D+00 ) then
           cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
        else
           cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
        end if
        cbi(k) = fac * cbi(k)
        fac = - fac
     end do
  end if

  cdi(0) = cbi(1)
  cdk(0) = -cbk(1)
  do k = 1, nm
     cdi(k) = cbi(k-1) - k / z * cbi(k)
     cdk(k) = - cbk(k-1) - k / z * cbk(k)
  end do

  return
end subroutine ciknb




  !*****************************************************************************80
  !
  !! CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
  !
  !  Discussion:
  !
  !    Compute the modified Bessel functions Iv(z), Kv(z)
  !    and their derivatives for an arbitrary order and
  !    complex argument
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:       
  !
  !    Input, real(8) V, the order of the functions.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
  !    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
  !
subroutine cikva ( v, z, vm, cbi, cdi, cbk, cdk )
  !CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
  !
  !Compute the modified Bessel functions Iv(z), Kv(z)
  !and their derivatives for an arbitrary order and
  !complex argument
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of the functions
  !@param z: in: the argument
  !@param vm: out: the highest order computed
  !@param cbi: out: In+v0(z)
  !@param cdi: out: In+v0'(z)
  !@param cbk: out: Kn+v0(z)
  !@param cdk: out: Kn+v0'(z)
  !
  implicit none

  real(8)    :: a0
  complex(8) :: ca
  complex(8) :: ca1
  complex(8) :: ca2
  complex(8) :: cb
  complex(8) :: cbi(0:*)
  complex(8) :: cbi0
  complex(8) :: cdi(0:*)
  complex(8) :: cbk(0:*)
  complex(8) :: cbk0
  complex(8) :: cbk1
  complex(8) :: cdk(0:*)
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: cg0
  complex(8) :: cg1
  complex(8) :: cgk
  complex(8) :: ci
  complex(8) :: ci0
  complex(8) :: cp
  complex(8) :: cr
  complex(8) :: cr1
  complex(8) :: cr2
  complex(8) :: cs
  complex(8) :: csu
  complex(8) :: ct
  complex(8) :: cvk
  real(8)    :: gan
  real(8)    :: gap
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: piv
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: v0n
  real(8)    :: v0p
  real(8)    :: vm
  real(8)    :: vt
  real(8)    :: w0
  real(8)    :: ws
  real(8)    :: ws0
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  piv = pi * v0
  vt = 4.0D+00 * v0 * v0

  if ( n == 0 ) then
     n = 1
  end if

  if ( a0 < 1.0D-100 ) then

     do k = 0, n
        cbi(k) = 0.0D+00
        cdi(k) = 0.0D+00
        cbk(k) = -1.0D+300
        cdk(k) = 1.0D+300
     end do

     if ( v0 == 0.0D+00 ) then
        cbi(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cdi(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     end if

     vm = v
     return

  end if

  if ( a0 < 35.0D+00 ) then
     k0 = 14
  else if ( a0 < 50.0D+00 ) then
     k0 = 10
  else
     k0 = 8
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 < 18.0D+00 ) then

     if ( v0 == 0.0D+00 ) then
        ca1 = cmplx (1.0D+00, 0.0D+00, kind = 8 )
     else
        v0p = 1.0D+00 + v0
        call gammaf ( v0p, gap )
        ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
     end if

     ci0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 50
        cr = 0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
        ci0 = ci0 + cr
        if ( abs ( cr ) < abs ( ci0 ) * 1.0D-15 ) then
           exit
        end if
     end do

     cbi0 = ci0 * ca1

  else

     ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
     cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cr = - 0.125D+00 * cr &
             * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
        cs = cs + cr
     end do
     cbi0 = ca * cs

  end if

  m = msta1 ( a0, 200 )

  if ( m < n ) then
     n = m
  else
     m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
  do k = m, 0, -1
     cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 + cf2
     if ( k <= n ) then
        cbi(k) = cf
     end if
     cf2 = cf1
     cf1 = cf
  end do

  cs = cbi0 / cf
  do k = 0, n
     cbi(k) = cs * cbi(k)
  end do

  if ( a0 <= 9.0D+00 ) then

     if ( v0 == 0.0D+00 ) then
        ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
        cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 50
           w0 = w0 + 1.0D+00 / k
           cr = 0.25D+00 * cr / ( k * k ) * z2
           cp = cr * ( w0 + ct )
           cs = cs + cp
           if ( 10 <= k .and. abs ( cp / cs ) < 1.0D-15 ) then
              exit
           end if
        end do

        cbk0 = ct + cs

     else

        v0n = 1.0D+00 - v0
        call gammaf ( v0n, gan )
        ca2 = 1.0D+00 / ( gan * ( 0.5D+00 * z1 ) ** v0 )
        ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
        csu = ca2 - ca1
        cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 50
           cr1 = 0.25D+00 * cr1 * z2 / ( k * ( k - v0 ) )
           cr2 = 0.25D+00 * cr2 * z2 / ( k * ( k + v0 ) )
           csu = csu + ca2 * cr1 - ca1 * cr2
           ws = abs ( csu )
           if ( 10 <= k .and. abs ( ws - ws0 ) / ws < 1.0D-15 ) then
              exit
           end if
           ws0 = ws
        end do

        cbk0 = 0.5D+00 * pi * csu / sin ( piv )

     end if

  else

     cb = exp ( - z1 ) * sqrt ( 0.5D+00 * pi / z1 )
     cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cr = 0.125D+00 * cr &
             * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
        cs = cs + cr
     end do
     cbk0 = cb * cs

  end if

  cbk1 = ( 1.0D+00 / z1 - cbi(1) * cbk0 ) / cbi(0)
  cbk(0) = cbk0
  cbk(1) = cbk1
  cg0 = cbk0
  cg1 = cbk1

  do k = 2, n
     cgk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / z1 * cg1 + cg0
     cbk(k) = cgk
     cg0 = cg1
     cg1 = cgk
  end do

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     do k = 0, n
        cvk = exp ( ( k + v0 ) * pi * ci )
        if ( imag ( z ) < 0.0D+00 ) then
           cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
           cbi(k) = cbi(k) / cvk
        else if ( 0.0D+00 < imag ( z ) ) then
           cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
           cbi(k) = cvk * cbi(k)
        end if
     end do
  end if

  cdi(0) = v0 / z * cbi(0) + cbi(1)
  cdk(0) = v0 / z * cbk(0) - cbk(1)
  do k = 1, n
     cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
     cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
  end do

  vm = n + v0

  return
end subroutine cikva




  !*****************************************************************************80
  !
  !! CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of the functions.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
  !    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
  !
subroutine cikvb ( v, z, vm, cbi, cdi, cbk, cdk )
  !CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of the functions
  !@param z: in: the argument
  !@param vm: out: the highest order computed
  !@param cbi: out: In+v0(z)
  !@param cdi: out: In+v0'(z)
  !@param cbk: out: Kn+v0(z)
  !@param cdk: out: Kn+v0'(z)
  !
  implicit none

  real(8)    :: a0
  complex(8) :: ca
  complex(8) :: ca1
  complex(8) :: ca2
  complex(8) :: cb
  complex(8) :: cbi(0:*)
  complex(8) :: cbi0
  complex(8) :: cdi(0:*)
  complex(8) :: cbk(0:*)
  complex(8) :: cbk0
  complex(8) :: cdk(0:*)
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: ci
  complex(8) :: ci0
  complex(8) :: ckk
  complex(8) :: cp
  complex(8) :: cr
  complex(8) :: cr1
  complex(8) :: cr2
  complex(8) :: cs
  complex(8) :: csu
  complex(8) :: ct
  complex(8) :: cvk
  real(8)    :: gan
  real(8)    :: gap
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: piv
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: v0n
  real(8)    :: v0p
  real(8)    :: vm
  real(8)    :: vt
  real(8)    :: w0
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2

  z1 = z
  z2 = z * z
  a0 = abs ( z )
  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  n = int ( v )
  v0 = v - n
  piv = pi * v0
  vt = 4.0D+00 * v0 * v0

  if ( n == 0 ) then
     n = 1
  end if

  if ( a0 < 1.0D-100 ) then
     do k = 0, n
        cbi(k) = 0.0D+00
        cdi(k) = 0.0D+00
        cbk(k) = -1.0D+300
        cdk(k) = 1.0D+300
     end do
     if ( v0 == 0.0D+00 ) then
        cbi(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cdi(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     end if
     vm = v
     return
  end if

  if ( a0 < 35.0D+00 ) then
     k0 = 14
  else if ( a0 < 50.0D+00 ) then
     k0 = 10
  else
     k0 = 8
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 < 18.0D+00 ) then

     if ( v0 == 0.0D+00 ) then
        ca1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     else
        v0p = 1.0D+00 + v0
        call gammaf ( v0p, gap )
        ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
     end if

     ci0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 50
        cr = 0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
        ci0 = ci0 + cr
        if ( abs ( cr / ci0 ) < 1.0D-15 ) then
           exit
        end if
     end do

     cbi0 = ci0 * ca1

  else

     ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
     cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cr = -0.125D+00 * cr &
             * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
        cs = cs + cr
     end do
     cbi0 = ca * cs

  end if

  m = msta1 ( a0, 200 )
  if ( m < n ) then
     n = m
  else
     m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
  do k = m, 0, -1
     cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 + cf2
     if ( k <= n ) then
        cbi(k) = cf
     end if
     cf2 = cf1
     cf1 = cf
  end do
  cs = cbi0 / cf

  do k = 0, n
     cbi(k) = cs * cbi(k)
  end do

  if ( a0 <= 9.0D+00 ) then

     if ( v0 == 0.0D+00 ) then

        ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
        cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 50
           w0 = w0 + 1.0D+00 / k
           cr = 0.25D+00 * cr / ( k * k ) * z2
           cp = cr * ( w0 + ct )
           cs = cs + cp
           if ( 10 <= k .and. abs ( cp / cs ) < 1.0D-15 ) then
              exit
           end if
        end do

        cbk0 = ct + cs

     else

        v0n = 1.0D+00 - v0
        call gammaf ( v0n, gan )
        ca2 = 1.0D+00 / ( gan * ( 0.5D+00 * z1 ) ** v0 )
        ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
        csu = ca2 - ca1
        cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 50
           cr1 = 0.25D+00 * cr1 * z2 / ( k * ( k - v0 ) )
           cr2 = 0.25D+00 * cr2 * z2 / ( k * ( k + v0 ) )
           cp = ca2 * cr1 - ca1 * cr2
           csu = csu + cp
           if ( 10 <= k .and. abs ( cp / csu ) < 1.0D-15 ) then
              exit
           end if
        end do

        cbk0 = 0.5D+00 * pi * csu / sin ( piv )

     end if

  else

     cb = exp ( -z1 ) * sqrt ( 0.5D+00 * pi / z1 )
     cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cr = 0.125D+00 * cr * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) &
             / ( k * z1 )
        cs = cs + cr
     end do

     cbk0 = cb * cs

  end if

  cbk(0) = cbk0

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     do k = 0, n
        cvk = exp ( ( k + v0 ) * pi * ci )
        if ( imag ( z ) < 0.0D+00 ) then
           cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
           cbi(k) = cbi(k) / cvk
        else if ( 0.0D+00 < imag ( z ) ) then
           cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
           cbi(k) = cvk * cbi(k)
        end if
     end do
  end if

  do k = 1, n
     ckk = ( 1.0D+00 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
     cbk(k) = ckk
  end do

  cdi(0) = v0 / z * cbi(0) + cbi(1)
  cdk(0) = v0 / z * cbk(0) - cbk(1)
  do k = 1, n
     cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
     cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
  end do

  vm = n + v0

  return
end subroutine cikvb




  !*****************************************************************************80
  !
  !! CISIA computes cosine Ci(x) and sine integrals Si(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument of Ci(x) and Si(x).
  !
  !    Output, real(8) CI, SI, the values of Ci(x) and Si(x).
  !
subroutine cisia ( x, ci, si )
  !CISIA computes cosine Ci(x) and sine integrals Si(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument of Ci(x) and Si(x)
  !@param ci: out: Ci(x)
  !@param si: out: Si(x)
  !
  implicit none

  real(8)    :: bj(101)
  real(8)    :: ci
  real(8)    :: el
  real(8)    :: eps
  integer(4) :: k
  integer(4) :: m
  real(8)    :: p2
  real(8)    :: si
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xa
  real(8)    :: xa0
  real(8)    :: xa1
  real(8)    :: xcs
  real(8)    :: xf
  real(8)    :: xg
  real(8)    :: xg1
  real(8)    :: xg2
  real(8)    :: xr
  real(8)    :: xs
  real(8)    :: xss

  p2 = 1.570796326794897D+00
  el = 0.5772156649015329D+00
  eps = 1.0D-15
  x2 = x * x

  if ( x == 0.0D+00 ) then

     ci = -1.0D+300
     si = 0.0D+00

  else if ( x <= 16.0D+00 ) then

     xr = -0.25D+00 * x2
     ci = el + log ( x ) + xr
     do k = 2, 40
        xr = -0.5D+00 * xr * ( k - 1 ) / ( k * k * ( 2 * k - 1 ) ) * x2
        ci = ci + xr
        if ( abs ( xr ) < abs ( ci ) * eps ) then
           exit
        end if
     end do

     xr = x
     si = x
     do k = 1, 40
        xr = -0.5D+00 * xr * ( 2 * k - 1 ) / k / ( 4 * k * k + 4 * k + 1 ) * x2
        si = si + xr
        if ( abs ( xr ) < abs ( si ) * eps ) then
           return
        end if
     end do

  else if ( x <= 32.0D+00 ) then

     m = int ( 47.2D+00 + 0.82D+00 * x )
     xa1 = 0.0D+00
     xa0 = 1.0D-100
     do k = m, 1, -1
        xa = 4.0D+00 * k * xa0 / x - xa1
        bj(k) = xa
        xa1 = xa0
        xa0 = xa
     end do
     xs = bj(1)
     do k = 3, m, 2
        xs = xs + 2.0D+00 * bj(k)
     end do
     bj(1) = bj(1) / xs
     do k = 2, m
        bj(k) = bj(k) / xs
     end do
     xr = 1.0D+00
     xg1 = bj(1)
     do k = 2, m
        xr = 0.25D+00 * xr * ( 2.0D+00 * k - 3.0D+00 ) **2 &
             / ( ( k - 1.0D+00 ) * ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) * x
        xg1 = xg1 + bj(k) * xr
     end do

     xr = 1.0D+00
     xg2 = bj(1)
     do k = 2, m
        xr = 0.25D+00 * xr * ( 2.0D+00 * k - 5.0D+00 )**2 &
             / ( ( k-1.0D+00 ) * ( 2.0D+00 * k - 3.0D+00 ) ** 2 ) * x
        xg2 = xg2 + bj(k) * xr
     end do

     xcs = cos ( x / 2.0D+00 )
     xss = sin ( x / 2.0D+00 )
     ci = el + log ( x ) - x * xss * xg1 + 2.0 * xcs * xg2 - 2.0 * xcs * xcs
     si = x * xcs * xg1 + 2.0 * xss * xg2 - sin ( x )

  else

     xr = 1.0D+00
     xf = 1.0D+00
     do k = 1, 9
        xr = -2.0D+00 * xr * k * ( 2 * k - 1 ) / x2
        xf = xf + xr
     end do
     xr = 1.0D+00 / x
     xg = xr
     do k = 1, 8
        xr = -2.0D+00 * xr * ( 2 * k + 1 ) * k / x2
        xg = xg + xr
     end do
     ci = xf * sin ( x ) / x - xg * cos ( x ) / x
     si = p2 - xf * cos ( x ) / x - xg * sin ( x ) / x

  end if

  return
end subroutine cisia




  !*****************************************************************************80
  !
  !! CISIB computes cosine and sine integrals.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    20 March 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument of Ci(x) and Si(x).
  !
  !    Output, real(8) CI, SI, the values of Ci(x) and Si(x).
  !
subroutine cisib ( x, ci, si )
  !CISIB computes cosine and sine integrals.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument of Ci(x) and Si(x)
  !@param ci: out: Ci(x)
  !@param si: out: Si(x)
  !
  implicit none

  real(8) :: ci
  real(8) :: fx
  real(8) :: gx
  real(8) :: si
  real(8) :: x
  real(8) :: x2

  x2 = x * x

  if ( x == 0.0D+00 ) then

     ci = -1.0D+300
     si = 0.0D+00

  else if ( x <= 1.0D+00 ) then

     ci = (((( -3.0D-08        * x2 &
          + 3.10D-06     ) * x2 &
          - 2.3148D-04   ) * x2 &
          + 1.041667D-02 ) * x2 &
          - 0.25D+00     ) * x2 + 0.577215665D+00 + log ( x )

     si = (((( 3.1D-07        * x2 &
          - 2.834D-05    ) * x2 &
          + 1.66667D-03  ) * x2 &
          - 5.555556D-02 ) * x2 + 1.0D+00 ) * x

  else

     fx = (((( x2              &
          + 38.027264D+00  ) * x2 &
          + 265.187033D+00 ) * x2 &
          + 335.67732D+00  ) * x2 &
          + 38.102495D+00  ) /    &
          (((( x2                 &
          + 40.021433D+00  ) * x2 &
          + 322.624911D+00 ) * x2 &
          + 570.23628D+00  ) * x2 &
          + 157.105423D+00 )

     gx = (((( x2               &
          + 42.242855D+00  ) * x2  &
          + 302.757865D+00 ) * x2  &
          + 352.018498D+00 ) * x2  &
          + 21.821899D+00 ) /      &
          (((( x2                  &
          + 48.196927D+00   ) * x2 &
          + 482.485984D+00  ) * x2 &
          + 1114.978885D+00 ) * x2 &
          + 449.690326D+00  ) / x

     ci = fx * sin ( x ) / x - gx * cos ( x ) / x

     si = 1.570796327D+00 - fx * cos ( x ) / x - gx * sin ( x ) / x

  end if

  return
end subroutine cisib




  !*****************************************************************************80
  !
  !! CJK: asymptotic expansion coefficients for Bessel functions of large order.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KM, the maximum value of K.
  !
  !    Output, real(8) A(L), the value of Cj(k) where j and k are 
  !    related to L by L = j+1+[k*(k+1)]/2; j,k = 0,1,...,Km.
  !
subroutine cjk ( km, a )
  !CJK: asymptotic expansion coefficients for Bessel functions of large order.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* A: the value of Cj(k) where j and k are related to L by L = j+1+[k*(k+1)]/2; j, k = 0, 1, ..., Km
  !
  !@param km: in: the maximum value of K
  !@param a: out: Cj(k), at L = j+1+k(k+1)/2, j,k = 0,1,...,Km
  !
  implicit none

  real(8)    :: a(*)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: g
  real(8)    :: g0
  integer(4) :: j
  integer(4) :: k
  integer(4) :: km
  integer(4) :: l1
  integer(4) :: l2
  integer(4) :: l3
  integer(4) :: l4

  a(1) = 1.0D+00
  f0 = 1.0D+00
  g0 = 1.0D+00
  do k = 0, km - 1
     l1 = ( k + 1 ) * ( k + 2 ) / 2 + 1
     l2 = ( k + 1 ) * ( k + 2 ) / 2 + k + 2
     f = ( 0.5D+00 * k + 0.125D+00 / ( k + 1 ) ) * f0
     g = - ( 1.5D+00 * k + 0.625D+00 &
          / ( 3.0D+00 * ( k + 1.0D+00 ) ) ) * g0
     a(l1) = f
     a(l2) = g
     f0 = f
     g0 = g
  end do

  do k = 1, km - 1
     do j = 1, k
        l3 = k * ( k + 1 ) / 2 + j + 1
        l4 = ( k + 1 ) * ( k + 2 ) / 2 + j + 1
        a(l4) = ( j + 0.5D+00 * k + 0.125D+00 &
             / ( 2.0D+00 * j + k + 1.0D+00 ) ) * a(l3) &
             - ( j + 0.5D+00 * k - 1.0D+00 + 0.625D+00 &
             / ( 2.0D+00 * j + k + 1.0D+00 ) ) * a(l3-1)
     end do
  end do

  return
end subroutine cjk




  !*****************************************************************************80
  !
  !! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, 
  !    CDY1, the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), 
  !    Y1(z), Y1'(z).
  !
subroutine cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )
  !CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param cbj0: out: J0(z)
  !@param cdj0: out: J0'(z)
  !@param cbj1: out: J1(z)
  !@param cdj1: out: J1'(z)
  !@param cby0: out: Y0(z)
  !@param cdy0: out: Y0'(z)
  !@param cby1: out: Y1(z)
  !@param cdy1: Y1'(z)
  !
  implicit none

  real(8), save, dimension ( 12 ) :: a = (/ &
       -0.703125D-01,0.112152099609375D+00, &
       -0.5725014209747314D+00,0.6074042001273483D+01, &
       -0.1100171402692467D+03,0.3038090510922384D+04, &
       -0.1188384262567832D+06,0.6252951493434797D+07, &
       -0.4259392165047669D+09,0.3646840080706556D+11, &
       -0.3833534661393944D+13,0.4854014686852901D+15 /)
  real(8)    :: a0
  real(8), save, dimension ( 12 ) :: a1 = (/ &
       0.1171875D+00,-0.144195556640625D+00, &
       0.6765925884246826D+00,-0.6883914268109947D+01, &
       0.1215978918765359D+03,-0.3302272294480852D+04, &
       0.1276412726461746D+06,-0.6656367718817688D+07, &
       0.4502786003050393D+09,-0.3833857520742790D+11, &
       0.4011838599133198D+13,-0.5060568503314727D+15 /)
  real(8), save, dimension ( 12 ) :: b = (/ &
       0.732421875D-01,-0.2271080017089844D+00, &
       0.1727727502584457D+01,-0.2438052969955606D+02, &
       0.5513358961220206D+03,-0.1825775547429318D+05, &
       0.8328593040162893D+06,-0.5006958953198893D+08, &
       0.3836255180230433D+10,-0.3649010818849833D+12, &
       0.4218971570284096D+14,-0.5827244631566907D+16 /)
  real(8), save, dimension ( 12 ) :: b1 = (/ &
       -0.1025390625D+00,0.2775764465332031D+00, &
       -0.1993531733751297D+01,0.2724882731126854D+02, &
       -0.6038440767050702D+03,0.1971837591223663D+05, &
       -0.8902978767070678D+06,0.5310411010968522D+08, &
       -0.4043620325107754D+10,0.3827011346598605D+12, &
       -0.4406481417852278D+14,0.6065091351222699D+16 /)
  complex(8) :: cbj0
  complex(8) :: cbj1
  complex(8) :: cby0
  complex(8) :: cby1
  complex(8) :: cdj0
  complex(8) :: cdj1
  complex(8) :: cdy0
  complex(8) :: cdy1
  complex(8) :: ci
  complex(8) :: cp
  complex(8) :: cp0
  complex(8) :: cp1
  complex(8) :: cq0
  complex(8) :: cq1
  complex(8) :: cr
  complex(8) :: cs
  complex(8) :: ct1
  complex(8) :: ct2
  complex(8) :: cu
  real(8)    :: el
  integer(4) :: k
  integer(4) :: k0
  real(8)    :: pi
  real(8)    :: rp2
  real(8)    :: w0
  real(8)    :: w1
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  rp2 = 2.0D+00 / pi
  ci = cmplx ( 0.0D+00, 1.0D+00 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 .eq. 0.0D+00 ) then
     cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cbj1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cdj0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cdj1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     cby0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cby1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cdy0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cdy1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 .le. 12.0D+00 ) then

     cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 40
        cr = -0.25D+00 * cr * z2 / ( k * k )
        cbj0 = cbj0 + cr
        if ( abs ( cr ) < abs ( cbj0 ) * 1.0D-15 ) then
           exit
        end if
     end do

     cbj1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 40
        cr = -0.25D+00 * cr * z2 / ( k * ( k + 1.0D+00 ) )
        cbj1 = cbj1 + cr
        if ( abs ( cr ) < abs ( cbj1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     cbj1 = 0.5D+00 * z1 * cbj1
     w0 = 0.0D+00
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 40
        w0 = w0 + 1.0D+00 / k
        cr = -0.25D+00 * cr / ( k * k ) * z2
        cp = cr * w0
        cs = cs + cp
        if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
           exit
        end if
     end do

     cby0 = rp2 * ( log ( z1 / 2.0D+00 ) + el ) * cbj0 - rp2 * cs
     w1 = 0.0D+00
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 40
        w1 = w1 + 1.0D+00 / k
        cr = -0.25D+00 * cr / ( k * ( k + 1 ) ) * z2
        cp = cr * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
        cs = cs + cp
        if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
           exit
        end if
     end do

     cby1 = rp2 * ( ( log ( z1 / 2.0D+00 ) + el ) * cbj1 &
          - 1.0D+00 / z1 - 0.25D+00 * z1 * cs )

  else

     if ( a0 < 35.0D+00 ) then
        k0 = 12
     else if ( a0 < 50.0D+00 ) then
        k0 = 10
     else
        k0 = 8
     end if

     ct1 = z1 - 0.25D+00 * pi
     cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
     end do
     cq0 = -0.125D+00 / z1
     do k = 1, k0
        cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
     end do
     cu = sqrt ( rp2 / z1 )
     cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
     cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
     ct2 = z1 - 0.75D+00 * pi
     cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
     end do
     cq1 = 0.375D+00 / z1
     do k = 1, k0
        cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
     end do
     cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
     cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     if ( imag ( z ) < 0.0D+00 ) then
        cby0 = cby0 - 2.0D+00 * ci * cbj0
        cby1 = - ( cby1 - 2.0D+00 * ci * cbj1 )
     else
        cby0 = cby0 + 2.0D+00 * ci * cbj0
        cby1 = - ( cby1 + 2.0D+00 * ci * cbj1 )
     end if
     cbj1 = -cbj1
  end if

  cdj0 = -cbj1
  cdj1 = cbj0 - 1.0D+00 / z * cbj1
  cdy0 = -cby1
  cdy1 = cby0 - 1.0D+00 / z * cby1

  return
end subroutine cjy01




  !*****************************************************************************80
  !
  !! CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Jv(z) and Yv(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CBJV, CDJV, CBYV, CDYV, the values of Jv(z), 
  !    Jv'(z), Yv(z), Yv'(z).
  !
subroutine cjylv ( v, z, cbjv, cdjv, cbyv, cdyv )
  !CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Jv(z) and Yv(z)
  !@param z: in: the argument
  !@param cbjv: out: Jv(z)
  !@param cdjv: out: Jv'(z)
  !@param cbyv: out: Yv(z)
  !@param cdyv: out: Yv'(z)
  !
  implicit none

  real(8)    :: a(91)
  complex(8) :: cbjv
  complex(8) :: cbyv
  complex(8) :: cdjv
  complex(8) :: cdyv
  complex(8) :: ceta
  complex(8) :: cf(12)
  complex(8) :: cfj
  complex(8) :: cfy
  complex(8) :: csj
  complex(8) :: csy
  complex(8) :: ct
  complex(8) :: ct2
  complex(8) :: cws
  integer(4) :: i
  integer(4) :: k
  integer(4) :: km
  integer(4) :: l
  integer(4) :: l0
  integer(4) :: lf
  real(8)    :: pi
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vr
  complex(8) :: z

  km = 12
  call cjk ( km, a )
  pi = 3.141592653589793D+00

  do l = 1, 0, -1

     v0 = v - l
     cws = sqrt ( 1.0D+00 - ( z / v0 ) * ( z / v0 ) )
     ceta = cws + log ( z / v0 / ( 1.0D+00 + cws ) )
     ct = 1.0D+00 / cws
     ct2 = ct * ct

     do k = 1, km
        l0 = k * ( k + 1 ) / 2 + 1
        lf = l0 + k
        cf(k) = a(lf)
        do i = lf - 1, l0, -1
           cf(k) = cf(k) * ct2 + a(i)
        end do
        cf(k) = cf(k) * ct ** k
     end do

     vr = 1.0D+00 / v0
     csj = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, km
        csj = csj + cf(k) * vr ** k
     end do
     cbjv = sqrt ( ct / ( 2.0D+00 * pi * v0 ) ) * exp ( v0 * ceta ) * csj
     if ( l == 1 ) then
        cfj = cbjv
     end if
     csy = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, km
        csy = csy + ( -1.0D+00 ) ** k * cf(k) * vr ** k
     end do
     cbyv = - sqrt ( 2.0D+00 * ct / ( pi * v0 ) ) * exp ( - v0 * ceta ) * csy
     if ( l == 1 ) then
        cfy = cbyv
     end if

  end do

  cdjv = - v / z * cbjv + cfj
  cdyv = - v / z * cbyv + cfy

  return
end subroutine cjylv




  !*****************************************************************************80
  !
  !! CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of Jn(z) and Yn(z).
  !
  !    Input, complex(8) Z, the argument of Jn(z) and Yn(z).
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8), CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
  !    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
  !
subroutine cjyna ( n, z, nm, cbj, cdj, cby, cdy )
  !CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of Jn(z) and Yn(z)
  !@param z: in: the argument of Jn(z) and Yn(z)
  !@param nm: out: the highest order computed
  !@param cbj: out: Jn(z), n=0,...,N
  !@param cdj: out: Jn'(z), n=0,...,N
  !@param cby: out: Yn(z), n=0,...,N
  !@param cdy: out: Yn'(z), n=0,...,N
  !
  implicit none
  integer(4) :: n
  real(8)    :: a0
  complex(8) :: cbj(0:n)
  complex(8) :: cbj0
  complex(8) :: cbj1
  complex(8) :: cby(0:n)
  complex(8) :: cby0
  complex(8) :: cby1
  complex(8) :: cdj(0:n)
  complex(8) :: cdj0
  complex(8) :: cdj1
  complex(8) :: cdy(0:n)
  complex(8) :: cdy0
  complex(8) :: cdy1
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: cg0
  complex(8) :: cg1
  complex(8) :: ch0
  complex(8) :: ch1
  complex(8) :: ch2
  complex(8) :: cj0
  complex(8) :: cj1
  complex(8) :: cjk
  complex(8) :: cp11
  complex(8) :: cp12
  complex(8) :: cp21
  complex(8) :: cp22
  complex(8) :: cs
  complex(8) :: cyk
  complex(8) :: cyl1
  complex(8) :: cyl2
  complex(8) :: cylk
  integer(4) :: k
  integer(4) :: lb
  integer(4) :: lb0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: pi
  real(8)    :: wa
  real(8)    :: ya0
  real(8)    :: ya1
  real(8)    :: yak
  complex(8) :: z

  pi = 3.141592653589793D+00
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-100 ) then
     do k = 0, n
        cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do
     cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     return
  end if

  call cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )
  cbj(0) = cbj0
  cbj(1) = cbj1
  cby(0) = cby0
  cby(1) = cby1
  cdj(0) = cdj0
  cdj(1) = cdj1
  cdy(0) = cdy0
  cdy(1) = cdy1

  if ( n <= 1 ) then
     return
  end if

  if ( n < int ( 0.25D+00 * a0 ) ) then

     cj0 = cbj0
     cj1 = cbj1
     do k = 2, n
        cjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cj1 - cj0
        cbj(k) = cjk
        cj0 = cj1
        cj1 = cjk
     end do

  else

     m = msta1 ( a0, 200 )

     if ( m < n ) then
        nm = m
     else
        m = msta2 ( a0, n, 15 )
     end if

     cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 ) 
     cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
     do k = m, 0, -1
        cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
        if ( k <= nm ) then
           cbj(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
     end do

     if ( abs ( cbj1 ) < abs ( cbj0 ) ) then
        cs = cbj0 / cf
     else
        cs = cbj1 / cf2
     end if

     do k = 0, nm
        cbj(k) = cs * cbj(k)
     end do

  end if

  do k = 2, nm
     cdj(k) = cbj(k-1) - k / z * cbj(k)
  end do
  ya0 = abs ( cby0 )
  lb = 0
  cg0 = cby0
  cg1 = cby1
  do k = 2, nm
     cyk = 2.0D+00 * ( k - 1.0D+00 ) / z * cg1 - cg0
     if ( abs ( cyk ) <= 1.0D+290 ) then         
        yak = abs ( cyk )
        ya1 = abs ( cg0 )
        if ( yak < ya0 .and. yak < ya1 ) then
           lb = k
        end if
        cby(k) = cyk
        cg0 = cg1
        cg1 = cyk
     end if
  end do

  if ( 4 < lb  .and. imag ( z ) /= 0.0D+00 ) then

     do

        if ( lb == lb0 ) then
           exit
        end if

        ch2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        ch1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        lb0 = lb
        do k = lb, 1, -1
           ch0 = 2.0D+00 * k / z * ch1 - ch2
           ch2 = ch1
           ch1 = ch0
        end do
        cp12 = ch0
        cp22 = ch2
        ch2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        ch1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = lb, 1, -1
           ch0 = 2.0D+00 * k / z * ch1 - ch2
           ch2 = ch1
           ch1 = ch0
        end do
        cp11 = ch0
        cp21 = ch2

        if ( lb == nm ) then
           cbj(lb+1) = 2.0D+00 * lb / z * cbj(lb) - cbj(lb-1)
        end if

        if ( abs ( cbj(1) ) < abs ( cbj(0) ) ) then
           cby(lb+1) = ( cbj(lb+1) * cby0 - 2.0D+00 * cp11 / ( pi * z ) ) / cbj(0)
           cby(lb) = ( cbj(lb) * cby0 + 2.0D+00 * cp12 / ( pi * z ) ) / cbj(0)
        else
           cby(lb+1) = ( cbj(lb+1) * cby1 - 2.0D+00 * cp21 / ( pi * z ) ) / cbj(1)
           cby(lb) = ( cbj(lb) * cby1 + 2.0D+00 * cp22 / ( pi * z ) )  / cbj(1)
        end if

        cyl2 = cby(lb+1)
        cyl1 = cby(lb)
        do k = lb - 1, 0, -1
           cylk = 2.0D+00 * ( k + 1.0D+00 ) / z * cyl1 - cyl2
           cby(k) = cylk
           cyl2 = cyl1
           cyl1 = cylk
        end do

        cyl1 = cby(lb)
        cyl2 = cby(lb+1)
        do k = lb + 1, nm - 1
           cylk = 2.0D+00 * k / z * cyl2 - cyl1
           cby(k+1) = cylk
           cyl1 = cyl2
           cyl2 = cylk
        end do

        do k = 2, nm
           wa = abs ( cby(k) )
           if ( wa < abs ( cby(k-1) ) ) then
              lb = k
           end if
        end do

     end do

  end if

  do k = 2, nm
     cdy(k) = cby(k-1) - k / z * cby(k)
  end do

  return
end subroutine cjyna




  !*****************************************************************************80
  !
  !! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) N, the order of Jn(z) and Yn(z).
  !
  !    Input, complex(8) Z, the argument of Jn(z) and Yn(z).
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N), 
  !    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
  !
subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )
  !CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of Jn(z) and Yn(z)
  !@param z: in: the argument of Jn(z) and Yn(z)
  !@param nm: out: the highest order computed
  !@param cbj: out: Jn(z)
  !@param cdj: out: Jn'(z)
  !@param cby: out: Yn(z)
  !@param cdy: out: Yn'(z)
  !
  implicit none

  integer(4) :: n

  real(8), save, dimension ( 4 ) :: a = (/ &
       -0.7031250000000000D-01, 0.1121520996093750D+00, &
       -0.5725014209747314D+00, 0.6074042001273483D+01 /)
  real(8)    :: a0
  real(8), save, dimension ( 4 ) :: a1 = (/ &
       0.1171875000000000D+00,-0.1441955566406250D+00, &
       0.6765925884246826D+00,-0.6883914268109947D+01 /)
  real(8), save, dimension ( 4 ) :: b = (/  &
       0.7324218750000000D-01,-0.2271080017089844D+00, &
       0.1727727502584457D+01,-0.2438052969955606D+02 /)
  real(8), save, dimension ( 4 ) :: b1 = (/ &
       -0.1025390625000000D+00,0.2775764465332031D+00, &
       -0.1993531733751297D+01,0.2724882731126854D+02 /)
  complex(8) :: cbj(0:n)
  complex(8) :: cbj0
  complex(8) :: cbj1
  complex(8) :: cbjk
  complex(8) :: cbs
  complex(8) :: cby(0:n)
  complex(8) :: cby0
  complex(8) :: cby1
  complex(8) :: cdj(0:n)
  complex(8) :: cdy(0:n)
  complex(8) :: ce
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: cp0
  complex(8) :: cp1
  complex(8) :: cq0
  complex(8) :: cq1
  complex(8) :: cs0
  complex(8) :: csu
  complex(8) :: csv
  complex(8) :: ct1
  complex(8) :: ct2
  complex(8) :: cu
  complex(8) :: cyy
  real(8)    :: el
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: pi
  real(8)    :: r2p
  real(8)    :: y0
  complex(8) :: z

  el = 0.5772156649015329D+00
  pi = 3.141592653589793D+00
  r2p = 0.63661977236758D+00
  y0 = abs ( imag ( z ) )
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-100 ) then
     do k = 0, n
        cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do
     cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     return
  end if

  if ( a0 <= 300.0D+00 .or. 80 < n ) then

     if ( n == 0 ) then
        nm = 1
     end if
     m = msta1 ( a0, 200 )
     if ( m < nm ) then
        nm = m
     else
        m = msta2 ( a0, nm, 15 )
     end if

     cbs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     csu = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     csv = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )

     do k = m, 0, -1
        cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
        if ( k <= nm ) then
           cbj(k) = cf
        end if
        if ( k == 2 * int ( k / 2 ) .and. k .ne. 0 ) then
           if ( y0 <= 1.0D+00 ) then
              cbs = cbs + 2.0D+00 * cf
           else
              cbs = cbs + ( -1.0D+00 ) ** ( k / 2 ) * 2.0D+00 * cf
           end if
           csu = csu + ( -1.0D+00 ) ** ( k / 2 ) * cf / k
        else if ( 1 < k ) then
           csv = csv + ( -1.0D+00 ) ** ( k / 2 ) * k / ( k * k - 1.0D+00 ) * cf
        end if
        cf2 = cf1
        cf1 = cf
     end do

     if ( y0 <= 1.0D+00 ) then
        cs0 = cbs + cf
     else
        cs0 = ( cbs + cf ) / cos ( z )
     end if

     do k = 0, nm
        cbj(k) = cbj(k) / cs0
     end do

     ce = log ( z / 2.0D+00 ) + el
     cby(0) = r2p * ( ce * cbj(0) - 4.0D+00 * csu / cs0 )
     cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0D+00 ) * cbj(1) &
          - 4.0D+00 * csv / cs0 )

  else

     ct1 = z - 0.25D+00 * pi
     cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 4
        cp0 = cp0 + a(k) * z ** ( - 2 * k )
     end do
     cq0 = -0.125D+00 / z
     do k = 1, 4
        cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
     end do
     cu = sqrt ( r2p / z )
     cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
     cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
     cbj(0) = cbj0
     cby(0) = cby0
     ct2 = z - 0.75D+00 * pi
     cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 4
        cp1 = cp1 + a1(k) * z ** ( - 2 * k )
     end do
     cq1 = 0.375D+00 / z
     do k = 1, 4
        cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
     end do
     cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
     cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )
     cbj(1) = cbj1
     cby(1) = cby1
     do k = 2, nm
        cbjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cbj1 - cbj0
        cbj(k) = cbjk
        cbj0 = cbj1
        cbj1 = cbjk
     end do
  end if

  cdj(0) = -cbj(1)
  do k = 1, nm
     cdj(k) = cbj(k-1) - k / z * cbj(k)
  end do

  if ( 1.0D+00 < abs ( cbj(0) ) ) then
     cby(1) = ( cbj(1) * cby(0) - 2.0D+00 / ( pi * z ) ) / cbj(0)
  end if

  do k = 2, nm
     if ( abs ( cbj(k-2) ) <= abs ( cbj(k-1) ) ) then
        cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
     else
        cyy = ( cbj(k) * cby(k-2) - 4.0D+00 * ( k - 1.0D+00 ) &
             / ( pi * z * z ) ) / cbj(k-2)
     end if
     cby(k) = cyy
  end do

  cdy(0) = -cby(1)
  do k = 1, nm
     cdy(k) = cby(k-1) - k / z * cby(k)
  end do

  return
end subroutine cjynb




  !*****************************************************************************80
  !
  !! CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Jv(z) and Yv(z).
  !
  !    Input, complex(8), the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
  !    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
  !
subroutine cjyva ( v, z, vm, cbj, cdj, cby, cdy )
  !CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Jv(z) and Yv(z)
  !@param z: the argument
  !@param vm: out: the highest order computed
  !@param cbj: out: Jn+v0(z)
  !@param cdj: out: Jn+v0'(z)
  !@param cby: out: Yn+v0(z)
  !@param cdy: out: Yn+v0'(z)
  !
  implicit none

  real(8)    :: a0
  complex(8) :: ca
  complex(8) :: ca0
  complex(8) :: cb
  complex(8) :: cbj(0:*)
  complex(8) :: cby(0:*)
  complex(8) :: cck
  complex(8) :: cdj(0:*)
  complex(8) :: cdy(0:*)
  complex(8) :: cec
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: cfac0
  complex(8) :: cfac1
  complex(8) :: cg0
  complex(8) :: cg1
  complex(8) :: ch0
  complex(8) :: ch1
  complex(8) :: ch2
  complex(8) :: ci
  complex(8) :: cju0
  complex(8) :: cju1
  complex(8) :: cjv0
  complex(8) :: cjv1
  complex(8) :: cjvl
  complex(8) :: cp11
  complex(8) :: cp12
  complex(8) :: cp21
  complex(8) :: cp22
  complex(8) :: cpz
  complex(8) :: cqz
  complex(8) :: cr
  complex(8) :: cr0
  complex(8) :: cr1
  complex(8) :: crp
  complex(8) :: crq
  complex(8) :: cs
  complex(8) :: cs0
  complex(8) :: cs1
  complex(8) :: csk
  complex(8) :: cyk
  complex(8) :: cyl1
  complex(8) :: cyl2
  complex(8) :: cylk
  complex(8) :: cyv0
  complex(8) :: cyv1
  real(8)    :: ga
  real(8)    :: gb
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: l
  integer(4) :: lb
  integer(4) :: lb0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: pv0
  real(8)    :: pv1
  real(8)    :: rp2
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vg
  real(8)    :: vl
  real(8)    :: vm
  real(8)    :: vv
  real(8)    :: w0
  real(8)    :: w1
  real(8)    :: wa
  real(8)    :: ya0
  real(8)    :: ya1
  real(8)    :: yak
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2
  complex(8) :: zk

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  pv0 = pi * v0
  pv1 = pi * ( 1.0D+00 + v0 )

  if ( a0 < 1.0D-100 ) then

     do k = 0, n
        cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do

     if ( v0 == 0.0D+00 ) then
        cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     else
        cdj(0) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end if

     vm = v                     
     return

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 <= 12.0D+00 ) then

     do l = 0, 1
        vl = v0 + l
        cjvl = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           cr = -0.25D+00 * cr * z2 / ( k * ( k + vl ) )
           cjvl = cjvl + cr
           if ( abs ( cr ) < abs ( cjvl ) * 1.0D-15 ) then
              exit
           end if
        end do

        vg = 1.0D+00 + vl
        call gammaf ( vg, ga )
        ca = ( 0.5D+00 * z1 ) ** vl / ga

        if ( l == 0 ) then
           cjv0 = cjvl * ca
        else
           cjv1 = cjvl * ca
        end if

     end do

  else

     if ( a0 < 35.0D+00 ) then
        k0 = 11
     else if ( a0 <50.0D+00 ) then
        k0 = 10
     else
        k0 = 8
     end if

     do j = 0, 1
        vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
        cpz = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        crp = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, k0
           crp = - 0.78125D-02 * crp &
                * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )  &
                / ( k * ( 2.0D+00 * k - 1.0D+00 ) * z2 )
           cpz = cpz + crp
        end do
        cqz = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        crq = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, k0
           crq = -0.78125D-02 * crq &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
                / ( k * ( 2.0D+00 * k + 1.0D+00 ) * z2 )
           cqz = cqz + crq
        end do
        cqz = 0.125D+00 * ( vv - 1.0D+00 ) * cqz / z1
        zk = z1 - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
        ca0 = sqrt ( rp2 / z1 )
        cck = cos ( zk )
        csk = sin ( zk )
        if ( j == 0 ) then
           cjv0 = ca0 * ( cpz * cck - cqz * csk )
           cyv0 = ca0 * ( cpz * csk + cqz * cck )
        else if ( j == 1 ) then
           cjv1 = ca0 * ( cpz * cck - cqz * csk )
           cyv1 = ca0 * ( cpz * csk + cqz * cck )
        end if
     end do

  end if

  if ( a0 <= 12.0D+00 ) then

     if ( v0 .ne. 0.0D+00 ) then

        do l = 0, 1
           vl = v0 + l
           cjvl = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           do k = 1, 40
              cr = -0.25D+00 * cr * z2 / ( k * ( k - vl ) )
              cjvl = cjvl + cr
              if ( abs ( cr ) < abs ( cjvl ) * 1.0D-15 ) then
                 exit
              end if
           end do

           vg = 1.0D+00 - vl
           call gammaf ( vg, gb )
           cb = ( 2.0D+00 / z1 ) ** vl / gb
           if ( l == 0 ) then
              cju0 = cjvl * cb
           else
              cju1 = cjvl * cb
           end if
        end do
        cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )
        cyv1 = ( cjv1 * cos ( pv1 ) - cju1 ) / sin ( pv1 )

     else

        cec = log ( z1 / 2.0D+00 ) + 0.5772156649015329D+00
        cs0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        w0 = 0.0D+00
        cr0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 30
           w0 = w0 + 1.0D+00 / k
           cr0 = -0.25D+00 * cr0 / ( k * k ) * z2
           cs0 = cs0 + cr0 * w0
        end do
        cyv0 = rp2 * ( cec * cjv0 - cs0 )
        cs1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        w1 = 0.0D+00
        cr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 30
           w1 = w1 + 1.0D+00 / k
           cr1 = -0.25D+00 * cr1 / ( k * ( k + 1 ) ) * z2
           cs1 = cs1 + cr1 * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
        end do
        cyv1 = rp2 * ( cec * cjv1 - 1.0D+00 / z1 - 0.25D+00 * z1 * cs1 )

     end if

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then

     cfac0 = exp ( pv0 * ci )
     cfac1 = exp ( pv1 * ci )

     if ( imag ( z ) < 0.0D+00 ) then
        cyv0 = cfac0 * cyv0 - 2.0D+00 * ci * cos ( pv0 ) * cjv0
        cyv1 = cfac1 * cyv1 - 2.0D+00 * ci * cos ( pv1 ) * cjv1
        cjv0 = cjv0 / cfac0
        cjv1 = cjv1 / cfac1
     else if ( 0.0D+00 < imag ( z ) ) then
        cyv0 = cyv0 / cfac0 + 2.0D+00 * ci * cos ( pv0 ) * cjv0
        cyv1 = cyv1 / cfac1 + 2.0D+00 * ci * cos ( pv1 ) * cjv1
        cjv0 = cfac0 * cjv0
        cjv1 = cfac1 * cjv1
     end if

  end if

  cbj(0) = cjv0
  cbj(1) = cjv1

  if ( 2 <= n .and. n <= int ( 0.25D+00 * a0 ) ) then

     cf0 = cjv0
     cf1 = cjv1
     do k = 2, n
        cf = 2.0D+00 * ( k + v0 - 1.0D+00 ) / z * cf1 - cf0
        cbj(k) = cf
        cf0 = cf1
        cf1 = cf
     end do

  else if ( 2 <= n ) then

     m = msta1 ( a0, 200 )
     if ( m < n ) then
        n = m
     else
        m = msta2 ( a0, n, 15 )
     end if
     cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
     do k = m, 0, -1
        cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z * cf1 - cf2
        if ( k <= n ) then
           cbj(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
     end do
     if ( abs ( cjv1 ) < abs ( cjv0 ) ) then
        cs = cjv0 / cf
     else
        cs = cjv1 / cf2
     end if

     do k = 0, n
        cbj(k) = cs * cbj(k)
     end do

  end if

  cdj(0) = v0 / z * cbj(0) - cbj(1)
  do k = 1, n
     cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
  end do

  cby(0) = cyv0
  cby(1) = cyv1
  ya0 = abs ( cyv0 )
  lb = 0
  cg0 = cyv0
  cg1 = cyv1
  do k = 2, n
     cyk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / z * cg1 - cg0
     if ( abs ( cyk ) <= 1.0D+290 ) then
        yak = abs ( cyk )
        ya1 = abs ( cg0 )
        if ( yak < ya0 .and. yak < ya1 ) then
           lb = k
        end if
        cby(k) = cyk
        cg0 = cg1
        cg1 = cyk
     end if
  end do

  if ( 4 < lb .and. imag ( z ) /= 0.0D+00 ) then

     do

        if ( lb == lb0 ) then
           exit
        end if

        ch2 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        ch1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        lb0 = lb
        do k = lb, 1, -1
           ch0 = 2.0D+00 * ( k + v0 ) / z * ch1 - ch2
           ch2 = ch1
           ch1 = ch0
        end do
        cp12 = ch0
        cp22 = ch2
        ch2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        ch1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = lb, 1, -1
           ch0 = 2.0D+00 * ( k + v0 ) / z * ch1 - ch2
           ch2 = ch1
           ch1 = ch0
        end do
        cp11 = ch0
        cp21 = ch2

        if ( lb == n ) then
           cbj(lb+1) = 2.0D+00 * ( lb + v0 ) / z * cbj(lb) - cbj(lb-1)
        end if

        if ( abs ( cbj(1) ) < abs ( cbj(0) ) ) then
           cby(lb+1) = ( cbj(lb+1) * cyv0 - 2.0D+00 * cp11 / ( pi * z ) ) &
                / cbj(0)
           cby(lb) = ( cbj(lb) * cyv0 + 2.0D+00 * cp12 / ( pi * z ) ) / cbj(0)
        else
           cby(lb+1) = ( cbj(lb+1) * cyv1 - 2.0D+00 * cp21 / ( pi * z ) ) &
                / cbj(1)
           cby(lb) = ( cbj(lb) * cyv1 + 2.0D+00 * cp22 / ( pi * z ) ) / cbj(1)
        end if

        cyl2 = cby(lb+1)
        cyl1 = cby(lb)
        do k = lb - 1, 0, -1
           cylk = 2.0D+00 * ( k + v0 + 1.0D+00 ) / z * cyl1 - cyl2
           cby(k) = cylk
           cyl2 = cyl1
           cyl1 = cylk
        end do

        cyl1 = cby(lb)
        cyl2 = cby(lb+1)
        do k = lb + 1, n - 1
           cylk = 2.0D+00 * ( k + v0 ) / z * cyl2 - cyl1
           cby(k+1) = cylk
           cyl1 = cyl2
           cyl2 = cylk
        end do

        do k = 2, n
           wa = abs ( cby(k) )
           if ( wa < abs ( cby(k-1) ) ) then
              lb = k
           end if
        end do

     end do

  end if

  cdy(0) = v0 / z * cby(0) - cby(1)
  do k = 1, n
     cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
  end do
  vm = n + v0

  return
end subroutine cjyva




  !*****************************************************************************80
  !
  !! CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Jv(z) and Yv(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
  !    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
  !
subroutine cjyvb ( v, z, vm, cbj, cdj, cby, cdy )
  !CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Jv(z) and Yv(z)
  !@param z: in: the argument
  !@param vm: out: the highest order computed
  !@param cbj: out: Jn+v0(z)
  !@param cdj: out: Jn+v0'(z)
  !@param cby: out: Yn+v0(z)
  !@param cdy: out: Yn+v0'(z)
  !
  implicit none

  real(8)    :: a0
  complex(8) :: ca
  complex(8) :: ca0
  complex(8) :: cb
  complex(8) :: cbj(0:*)
  complex(8) :: cby(0:*)
  complex(8) :: cck
  complex(8) :: cdj(0:*)
  complex(8) :: cdy(0:*)
  complex(8) :: cec
  complex(8) :: cf
  complex(8) :: cf1
  complex(8) :: cf2
  complex(8) :: cfac0
  complex(8) :: ci
  complex(8) :: cju0
  complex(8) :: cjv0
  complex(8) :: cjvn
  complex(8) :: cpz
  complex(8) :: cqz
  complex(8) :: cr
  complex(8) :: cr0
  complex(8) :: crp
  complex(8) :: crq
  complex(8) :: cs
  complex(8) :: cs0
  complex(8) :: csk
  complex(8) :: cyv0
  complex(8) :: cyy
  real(8)    :: ga
  real(8)    :: gb
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: pv0
  real(8)    :: rp2
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vg
  real(8)    :: vm
  real(8)    :: vv
  real(8)    :: w0
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2
  complex(8) :: zk

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  pv0 = pi * v0

  if ( a0 < 1.0D-100 ) then

     do k = 0, n
        cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do

     if ( v0 == 0.0D+00 ) then
        cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
     else
        cdj(0) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end if

     vm = v
     return

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
     z1 = -z
  end if

  if ( a0 <= 12.0D+00 ) then

     cjv0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 40
        cr = -0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
        cjv0 = cjv0 + cr
        if ( abs ( cr ) < abs ( cjv0 ) * 1.0D-15 ) then
           exit
        end if
     end do

     vg = 1.0D+00 + v0
     call gammaf ( vg, ga )
     ca = ( 0.5D+00 * z1 ) ** v0 / ga
     cjv0 = cjv0 * ca

  else

     if ( a0 < 35.0D+00 ) then
        k0 = 11
     else if ( a0 < 50.0D+00 ) then
        k0 = 10
     else
        k0 = 8
     end if

     vv = 4.0D+00 * v0 * v0
     cpz = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     crp = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        crp = -0.78125D-02 * crp &
             * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
             * ( vv - ( 4.0D+00 * k - 1.0D+00 ) **2 ) &
             / ( k * ( 2.0D+00 * k - 1.0D+00 ) * z2 )
        cpz = cpz + crp
     end do
     cqz = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     crq = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, k0
        crq = -0.78125D-02 * crq &
             * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
             * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
             / ( k * ( 2.0D+00 * k + 1.0D+00 ) * z2 )
        cqz = cqz + crq
     end do
     cqz = 0.125D+00 * ( vv - 1.0D+00 ) * cqz / z1
     zk = z1 - ( 0.5D+00 * v0 + 0.25D+00 ) * pi
     ca0 = sqrt ( rp2 / z1 )
     cck = cos ( zk )
     csk = sin ( zk )
     cjv0 = ca0 * ( cpz * cck - cqz * csk )
     cyv0 = ca0 * ( cpz * csk + cqz * cck )

  end if

  if ( a0 <= 12.0D+00 ) then

     if ( v0 .ne. 0.0D+00 ) then

        cjvn = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           cr = -0.25D+00 * cr * z2 / ( k * ( k - v0 ) )
           cjvn = cjvn + cr
           if ( abs ( cr ) < abs ( cjvn ) * 1.0D-15 ) then
              exit
           end if
        end do

        vg = 1.0D+00 - v0
        call gammaf ( vg, gb )
        cb = ( 2.0D+00 / z1 ) ** v0 / gb
        cju0 = cjvn * cb
        cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )

     else

        cec = log ( z1 / 2.0D+00 ) + 0.5772156649015329D+00
        cs0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        w0 = 0.0D+00
        cr0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 30
           w0 = w0 + 1.0D+00 / k
           cr0 = -0.25D+00 * cr0 / ( k * k ) * z2
           cs0 = cs0 + cr0 * w0
        end do
        cyv0 = rp2 * ( cec * cjv0 - cs0 )

     end if

  end if

  if ( n .eq. 0 ) then
     n = 1
  end if

  m = msta1 ( a0, 200 )
  if ( m < n ) then
     n = m
  else
     m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
  do k = m, 0, -1
     cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 - cf2
     if ( k <= n ) then
        cbj(k) = cf
     end if
     cf2 = cf1
     cf1 = cf
  end do

  cs = cjv0 / cf
  do k = 0, n
     cbj(k) = cs * cbj(k)
  end do

  if ( real ( z, kind = 8 ) < 0.0D+00) then

     cfac0 = exp ( pv0 * ci )
     if ( imag ( z ) < 0.0D+00 ) then
        cyv0 = cfac0 * cyv0 - 2.0D+00 * ci * cos ( pv0 ) * cjv0
     else if ( 0.0D+00 < imag ( z ) ) then
        cyv0 = cyv0 / cfac0 + 2.0D+00 * ci * cos ( pv0 ) * cjv0
     end if

     do k = 0, n
        if ( imag ( z ) < 0.0D+00) then
           cbj(k) = exp ( - pi * ( k + v0 ) * ci ) * cbj(k)
        else if ( 0.0D+00 < imag ( z ) ) then
           cbj(k) = exp ( pi * ( k + v0 ) * ci ) * cbj(k)
        end if
     end do

     z1 = z1

  end if

  cby(0) = cyv0
  do k = 1, n
     cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
     cby(k) = cyy
  end do

  cdj(0) = v0 / z * cbj(0) - cbj(1)
  do k = 1, n
     cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
  end do

  cdy(0) = v0 / z * cby(0) - cby(1)
  do k = 1, n
     cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
  end do

  vm = n + v0

  return
end subroutine cjyvb




  !*****************************************************************************80
  !
  !! CLPMN: associated Legendre functions and derivatives for complex argument.
  !
  !  Discussion:
  !
  !    Compute the associated Legendre functions Pmn(z)   
  !    and their derivatives Pmn'(z) for a complex argument
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) MM, the physical dimension of CPM and CPD.
  !
  !    Input, integer(4) M, N, the order and degree of Pmn(z).
  !
  !    Input, real(8) X, Y, the real and imaginary parts of 
  !    the argument Z.
  !
  !    Output, complex(8) CPM(0:MM,0:N), CPD(0:MM,0:N), the values of
  !    Pmn(z) and Pmn'(z).
  !
subroutine clpmn ( mm, m, n, x, y, cpm, cpd )
  !CLPMN: associated Legendre functions and derivatives for complex argument.
  !
  !Compute the associated Legendre functions Pmn(z)
  !and their derivatives Pmn'(z) for a complex argument
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param mm: in: the physical dimension of CPM and CPD
  !@param m: in: the order of Pmn(z)
  !@param n: in: the degree of Pmn(z)
  !@param x: in: the real part of the argument Z
  !@param y: in: the imaginary part of the argument Z
  !@param cpm: out: Pmn(z)
  !@param cpd: out: Pmn'(z)
  !
  implicit none

  integer(4) :: mm

  complex(8) :: cpd(0:mm,0:n)
  complex(8) :: cpm(0:mm,0:n)
  integer(4) :: i
  integer(4) :: j
  integer(4) :: ls
  integer(4) :: m
  integer(4) :: n
  real(8)    :: x
  real(8)    :: y
  complex(8) :: z
  complex(8) :: zq
  complex(8) :: zs

  z = cmplx ( x, y, kind = 8 )

  do i = 0, n
     do j = 0, m
        cpm(j,i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        cpd(j,i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     end do
  end do

  cpm(0,0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

  if ( abs ( x ) == 1.0D+00 .and. y == 0.0D+00 ) then

     do i = 1, n
        cpm(0,i) = x ** i
        cpd(0,i) = 0.5D+00 * i * ( i + 1 ) * x ** ( i + 1 )
     end do

     do j = 1, n
        do i = 1, m
           if ( i == 1 ) then
              cpd(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
           else if ( i == 2 ) then
              cpd(i,j) = -0.25D+00 &
                   * ( j + 2 ) * ( j + 1 ) * j * ( j - 1 ) * x ** ( j + 1 )
           end if
        end do
     end do

     return

  end if

  if ( 1.0D+00 < abs ( z ) ) then
     ls = -1
  else
     ls = 1
  end if

  zq = sqrt ( ls * ( 1.0D+00 - z * z ) )
  zs = ls * ( 1.0D+00 - z * z )
  do i = 1, m
     cpm(i,i) = -ls * ( 2.0D+00 * i - 1.0D+00 ) * zq * cpm(i-1,i-1)
  end do
  do i = 0, m
     cpm(i,i+1) = ( 2.0D+00 * i + 1.0D+00 ) * z * cpm(i,i)
  end do

  do i = 0, m
     do j = i + 2, n
        cpm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * z * cpm(i,j-1) &
             - ( i + j - 1.0D+00 ) * cpm(i,j-2) ) / ( j - i )
     end do
  end do

  cpd(0,0) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  do j = 1, n
     cpd(0,j) = ls * j * ( cpm(0,j-1) - z * cpm(0,j) ) / zs
  end do

  do i = 1, m
     do j = i, n
        cpd(i,j) = ls * i * z * cpm(i,j) / zs &
             + ( j + i ) * ( j - i + 1.0D+00 ) / zq * cpm(i-1,j)
     end do
  end do

  return
end subroutine clpmn




  !*****************************************************************************80
  !
  !! CLPN computes Legendre functions and derivatives for complex argument.
  !
  !  Discussion:
  !
  !    Compute Legendre polynomials Pn(z) and their derivatives Pn'(z) for 
  !    a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the degree.
  !
  !    Input, real(8) X, Y, the real and imaginary parts 
  !    of the argument.
  !
  !    Output, complex(8) CPN(0:N), CPD(0:N), the values of Pn(z)
  !    and Pn'(z).
  !
subroutine clpn ( n, x, y, cpn, cpd )
  !CLPN computes Legendre functions and derivatives for complex argument.
  !
  !Compute Legendre polynomials Pn(z) and their derivatives Pn'(z) for
  !a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the degree
  !@param x: in: the real part of the argument
  !@param y: in: the imaginary part of the argument
  !@param cpn: out: Pn(z)
  !@param cpd: out: Pn'(z)
  !
  implicit none

  integer(4) :: n

  complex(8) :: cp0
  complex(8) :: cp1
  complex(8) :: cpd(0:n)
  complex(8) :: cpf
  complex(8) :: cpn(0:n)
  integer(4) :: k
  real(8)    :: x
  real(8)    :: y
  complex(8) :: z

  z = cmplx ( x, y, kind = 8 )

  cpn(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  cpn(1) = z
  cpd(0) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  cpd(1) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

  cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  cp1 = z
  do k = 2, n
     cpf = ( 2.0D+00 * k - 1.0D+00 ) / k * z * cp1 - ( k - 1.0D+00 ) / k * cp0
     cpn(k) = cpf
     if ( abs ( x ) == 1.0D+00 .and. y == 0.0D+00 ) then
        cpd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
     else
        cpd(k) = k * ( cp1 - z * cpf ) / ( 1.0D+00 - z * z )
     end if
     cp0 = cp1
     cp1 = cpf
  end do

  return
end subroutine clpn




  !*****************************************************************************80
  !
  !! CLQMN: associated Legendre functions and derivatives for complex argument.
  !
  !  Discussion:
  !
  !    This procedure computes the associated Legendre functions of the second 
  !    kind, Qmn(z) and Qmn'(z), for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) MM, the physical dimension of CQM and CQD.
  !
  !    Input, integer(4) M, N, the order and degree of Qmn(z).
  !
  !    Input, real(8) X, Y, the real and imaginary parts of the 
  !    argument Z.
  !
  !    Output, complex(8) CQM(0:MM,0:N), CQD(0:MM,0:N), the values of
  !    Qmn(z) and Qmn'(z).
  !
subroutine clqmn ( mm, m, n, x, y, cqm, cqd )
  !CLQMN: associated Legendre functions and derivatives for complex argument.
  !
  !This procedure computes the associated Legendre functions of the second
  !kind, Qmn(z) and Qmn'(z), for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param mm: in: the physical dimension of CQM and CQD
  !@param m: in: the order of Qmn(z)
  !@param n: in: the degree of Qmn(z)
  !@param x: in: the real part of the argument Z
  !@param y: in: the imaginary part of the argument Z
  !@param cqm: out: Qmn(z)
  !@param cqd: out: Qmn'(z)
  !
  implicit none

  integer(4) :: mm
  integer(4) :: n

  complex(8) :: cq0
  complex(8) :: cq1
  complex(8) :: cq10
  complex(8) :: cqf
  complex(8) :: cqf0
  complex(8) :: cqf1
  complex(8) :: cqf2
  complex(8) :: cqm(0:mm,0:n)
  complex(8) :: cqd(0:mm,0:n)
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: km
  integer(4) :: ls
  integer(4) :: m
  real(8)    :: x
  real(8)    :: xc
  real(8)    :: y
  complex(8) :: z
  complex(8) :: zq
  complex(8) :: zs

  z = cmplx ( x, y, kind = 8 )

  if ( abs ( x ) == 1.0D+00 .and. y == 0.0D+00 ) then
     do i = 0, m
        do j = 0, n
           cqm(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
           cqd(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        end do
     end do
     return
  end if

  xc = abs ( z )

  if ( imag ( z ) == 0.0D+00 .or. xc < 1.0D+00 ) then
     ls = 1
  end if

  if ( 1.0D+00 < xc ) then
     ls = -1
  end if

  zq = sqrt ( ls * ( 1.0D+00 - z * z ) )
  zs = ls * ( 1.0D+00 - z * z )
  cq0 = 0.5D+00 * log ( ls * ( 1.0D+00 + z ) / ( 1.0D+00 - z ) )

  if ( xc < 1.0001D+00 ) then

     cqm(0,0) = cq0
     cqm(0,1) = z * cq0 - 1.0D+00
     cqm(1,0) = -1.0D+00 / zq
     cqm(1,1) = - zq * ( cq0 + z / ( 1.0D+00 - z * z ) )
     do i = 0, 1
        do j = 2, n
           cqm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * z * cqm(i,j-1) &
                - ( j + i - 1.0D+00 ) * cqm(i,j-2) ) / ( j - i )
        end do
     end do

     do j = 0, n
        do i = 2, m
           cqm(i,j) = -2.0D+00 * ( i - 1.0D+00 ) * z / zq * cqm(i-1,j) &
                - ls * ( j + i - 1.0D+00 ) * ( j - i + 2.0D+00 ) * cqm(i-2,j)
        end do
     end do

  else

     if ( 1.1D+00 < xc ) then
        km = 40 + m + n
     else
        km = ( 40 + m + n ) * int ( - 1.0D+00 - 1.8D+00 * log ( xc - 1.0D+00 ) )
     end if

     cqf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cqf1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = km, 0, -1
        cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1 &
             - ( k + 2.0D+00 ) * cqf2 ) / ( k + 1.0D+00 )
        if ( k <= n ) then
           cqm(0,k) = cqf0
        end if
        cqf2 = cqf1
        cqf1 = cqf0
     end do

     do k = 0, n
        cqm(0,k) = cq0 * cqm(0,k) / cqf0
     end do

     cqf2 = 0.0D+00
     cqf1 = 1.0D+00
     do k = km, 0, -1
        cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1 &
             - ( k + 1.0D+00 ) * cqf2 ) / ( k + 2.0D+00 )
        if ( k <= n ) then
           cqm(1,k) = cqf0
        end if
        cqf2 = cqf1
        cqf1 = cqf0
     end do

     cq10 = -1.0D+00 / zq
     do k = 0, n 
        cqm(1,k) = cq10 * cqm(1,k) / cqf0
     end do

     do j = 0, n
        cq0 = cqm(0,j)
        cq1 = cqm(1,j)
        do i = 0, m - 2
           cqf = -2.0D+00 * ( i + 1 ) * z / zq * cq1 &
                + ( j - i ) * ( j + i + 1.0D+00 ) * cq0
           cqm(i+2,j) = cqf
           cq0 = cq1
           cq1 = cqf
        end do
     end do

  end if

  cqd(0,0) = ls / zs
  do j = 1, n
     cqd(0,j) = ls * j * ( cqm(0,j-1) - z * cqm(0,j) ) / zs
  end do

  do j = 0, n
     do i = 1, m
        cqd(i,j) = ls * i * z / zs * cqm(i,j) &
             + ( i + j ) * ( j - i + 1.0D+00 ) / zq * cqm(i-1,j)
     end do
  end do

  return
end subroutine clqmn




  !*****************************************************************************80
  !
  !! CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the degree of Qn(z).
  !
  !    Input, real(8) X, Y, the real and imaginary parts of the 
  !    argument Z.
  !
  !    Output, complex(8) CQN(0:N), CQD(0:N), the values of Qn(z) 
  !    and Qn'(z.
  !
subroutine clqn ( n, x, y, cqn, cqd )
  !CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the degree of Qn(z)
  !@param x: in: the real part of the argument Z
  !@param y: in: the imaginary part of the argument Z
  !@param cqn: out: Qn(z)
  !@param cqd: out: Qn'(z
  !
  implicit none

  integer(4) :: n

  complex(8) :: cq0
  complex(8) :: cq1
  complex(8) :: cqf0
  complex(8) :: cqf1
  complex(8) :: cqf2
  complex(8) :: cqn(0:n)
  complex(8) :: cqd(0:n)
  integer(4) :: k
  integer(4) :: km
  integer(4) :: ls
  real(8)    :: x
  real(8)    :: y
  complex(8) :: z

  z = cmplx ( x, y, kind = 8 )

  if ( z == 1.0D+00 ) then
     do k = 0, n
        cqn(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cqd(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     end do
     return
  end if

  if ( 1.0D+00 < abs ( z ) ) then
     ls = -1
  else
     ls = +1
  end if

  cq0 = 0.5D+00 * log ( ls * ( 1.0D+00 + z ) / ( 1.0D+00 - z ) )
  cq1 = z * cq0 - 1.0D+00
  cqn(0) = cq0
  cqn(1) = cq1

  if ( abs ( z ) < 1.0001D+00 ) then

     cqf0 = cq0
     cqf1 = cq1
     do k = 2, n
        cqf2 = ( ( 2.0D+00 * k - 1.0D+00 ) * z * cqf1 &
             - ( k - 1.0D+00 ) * cqf0 ) / k
        cqn(k) = cqf2
        cqf0 = cqf1
        cqf1 = cqf2
     end do

  else

     if ( 1.1D+00 < abs ( z ) ) then
        km = 40 + n
     else
        km = ( 40 + n ) * int ( - 1.0D+00 &
             - 1.8D+00 * log ( abs ( z - 1.0D+00 ) ) )
     end if

     cqf2 = 0.0D+00
     cqf1 = 1.0D+00
     do k = km, 0, -1
        cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1 &
             - ( k + 2.0D+00 ) * cqf2 ) / ( k + 1.0D+00 )
        if ( k <= n ) then
           cqn(k) = cqf0
        end if
        cqf2 = cqf1
        cqf1 = cqf0
     end do
     do k = 0, n
        cqn(k) = cqn(k) * cq0 / cqf0
     end do
  end if

  cqd(0) = ( cqn(1) - z * cqn(0) ) / ( z * z - 1.0D+00 )
  do k = 1, n
     cqd(k) = ( k * z * cqn(k) - k * cqn(k-1) ) / ( z * z - 1.0D+00 )
  end do

  return
end subroutine clqn




  !*****************************************************************************80
  !
  !! COMELP computes complete elliptic integrals K(k) and E(k).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) HK, the modulus.  0 <= HK <= 1.
  !
  !    Output, real(8) CK, CE, the values of K(HK) and E(HK).
  !
subroutine comelp ( hk, ck, ce )
  !COMELP computes complete elliptic integrals K(k) and E(k).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param hk: in: the modulus. 0 <= HK <= 1
  !@param ck: out: K(HK)
  !@param ce: out: E(HK)
  !
  implicit none

  real(8) :: ae
  real(8) :: ak
  real(8) :: be
  real(8) :: bk
  real(8) :: ce
  real(8) :: ck
  real(8) :: hk
  real(8) :: pk

  pk = 1.0D+00 - hk * hk

  if ( hk == 1.0D+00 ) then

     ck = 1.0D+300
     ce = 1.0D+00

  else

     ak = ((( &
          0.01451196212D+00   * pk &
          + 0.03742563713D+00 ) * pk &
          + 0.03590092383D+00 ) * pk &
          + 0.09666344259D+00 ) * pk &
          + 1.38629436112D+00

     bk = ((( &
          0.00441787012D+00   * pk &
          + 0.03328355346D+00 ) * pk &
          + 0.06880248576D+00 ) * pk &
          + 0.12498593597D+00 ) * pk &
          + 0.5D+00

     ck = ak - bk * log ( pk )

     ae = ((( &
          0.01736506451D+00   * pk &
          + 0.04757383546D+00 ) * pk &
          + 0.0626060122D+00  ) * pk &
          + 0.44325141463D+00 ) * pk &
          + 1.0D+00

     be = ((( &
          0.00526449639D+00   * pk &
          + 0.04069697526D+00 ) * pk &
          + 0.09200180037D+00 ) * pk &
          + 0.2499836831D+00  ) * pk

     ce = ae - be * log ( pk )

  end if

  return
end subroutine comelp




  !*****************************************************************************80
  !
  !! CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CPB(0:N), CPD(0:N), the values of Dn(z) 
  !    and Dn'(z).
  !
subroutine cpbdn ( n, z, cpb, cpd )
  !CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param z: in: the argument
  !@param cpb: out: Dn(z)
  !@param cpd: out: Dn'(z)
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  complex(8) :: c0
  complex(8) :: ca0
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cfa
  complex(8) :: cfb
  complex(8) :: cpb(0:n)
  complex(8) :: cpd(0:n)
  complex(8) :: cs0
  integer(4) :: k
  integer(4) :: m
  integer(4) :: n0
  integer(4) :: n1
  integer(4) :: nm1
  real(8)    :: pi
  real(8)    :: x
  complex(8) :: z
  complex(8) :: z1

  pi = 3.141592653589793D+00
  x = real ( z, kind = 8 )
  a0 = abs ( z )
  c0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  ca0 = exp ( -0.25D+00 * z * z )

  if ( 0 <= n ) then

     cf0 = ca0
     cf1 = z * ca0
     cpb(0) = cf0
     cpb(1) = cf1
     do k = 2, n
        cf = z * cf1 - ( k - 1.0D+00 ) * cf0
        cpb(k) = cf
        cf0 = cf1
        cf1 = cf
     end do

  else

     n0 = -n

     if ( x <= 0.0D+00 .or. abs ( z ) .eq. 0.0D+00 ) then

        cf0 = ca0
        cpb(0) = cf0
        z1 = - z
        if ( a0 <= 7.0D+00 ) then
           call cpdsa ( -1, z1, cf1 )
        else
           call cpdla ( -1, z1, cf1 )
        end if
        cf1 = sqrt ( 2.0D+00 * pi ) / ca0 - cf1
        cpb(1) = cf1
        do k = 2, n0
           cf = ( - z * cf1 + cf0 ) / ( k - 1.0D+00 )
           cpb(k) = cf
           cf0 = cf1
           cf1 = cf
        end do

     else

        if ( a0 <= 3.0D+00 ) then

           call cpdsa ( -n0, z, cfa )
           cpb(n0) = cfa
           n1 = n0 + 1
           call cpdsa ( -n1, z, cfb )
           cpb(n1) = cfb
           nm1 = n0 - 1
           do k = nm1, 0, -1
              cf = z * cfa + ( k + 1.0D+00 ) * cfb
              cpb(k) = cf
              cfb = cfa
              cfa = cf
           end do

        else

           m = 100 + abs ( n )
           cfa = c0
           cfb = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
           do k = m, 0, -1
              cf = z * cfb + ( k + 1.0D+00 ) * cfa
              if ( k <= n0 ) then
                 cpb(k) = cf
              end if
              cfa = cfb
              cfb = cf
           end do
           cs0 = ca0 / cf
           do k = 0, n0
              cpb(k) = cs0 * cpb(k)
           end do

        end if

     end if

  end if

  cpd(0) = -0.5D+00 * z * cpb(0)

  if ( 0 <= n ) then
     do k = 1, n
        cpd(k) = -0.5D+00 * z * cpb(k) + k * cpb(k-1)
     end do
  else
     do k = 1, n0
        cpd(k) = 0.5D+00 * z * cpb(k) - cpb(k-1)
     end do
  end if

  return
end subroutine cpbdn




  !****************************************************************************80
  !
  !! CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer N, the order.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CDN, the function value.
  !
subroutine cpdla ( n, z, cdn )
  !CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param z: in: the argument
  !@param cdn: out: the function value
  !
  implicit none

  complex(8) :: cb0
  complex(8) :: cdn
  complex(8) :: cr
  integer(4) :: k
  integer(4) :: n
  complex(8) :: z

  cb0 = z ** n * exp ( -0.25D+00 * z * z )
  cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
  cdn = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

  do k = 1, 16

     cr = -0.5D+00 * cr * ( 2.0D+00 * k - n - 1.0D+00 ) &
          * ( 2.0D+00 * k - n - 2.0D+00 ) / ( k * z * z )

     cdn = cdn + cr

     if ( abs ( cr ) < abs ( cdn ) * 1.0D-12 ) then
        exit
     end if

  end do

  cdn = cb0 * cdn

  return
end subroutine cpdla




  !*****************************************************************************80
  !
  !! CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CDN, the value of DN(z).
  !
subroutine cpdsa ( n, z, cdn )
  !CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param z: in: the argument
  !@param cdn: out: the value of DN(z)
  !
  implicit none

  complex(8) :: ca0
  complex(8) :: cb0
  complex(8) :: cdn
  complex(8) :: cdw
  complex(8) :: cr
  real(8)    :: eps
  real(8)    :: g0
  real(8)    :: g1
  real(8)    :: ga0
  real(8)    :: gm
  integer(4) :: m
  integer(4) :: n
  real(8)    :: pd
  real(8)    :: pi
  real(8)    :: sq2
  real(8)    :: va0
  real(8)    :: vm
  real(8)    :: vt
  real(8)    :: xn
  complex(8) :: z

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  sq2 = sqrt ( 2.0D+00 )
  ca0 = exp ( - 0.25D+00 * z * z )
  va0 = 0.5D+00 * ( 1.0D+00 - n )

  if ( n == 0 ) then

     cdn = ca0

  else

     if ( abs ( z ) == 0.0D+00 ) then

        if ( va0 <= 0.0D+00 .and. va0 == int ( va0 ) ) then
           cdn = 0.0D+00
        else
           call gaih ( va0, ga0 )
           pd = sqrt ( pi ) / ( 2.0D+00 ** ( -0.5D+00 * n ) * ga0 )
           cdn = cmplx ( pd, 0.0D+00, kind = 8 )
        end if

     else

        xn = - n
        call gaih ( xn, g1 )
        cb0 = 2.0D+00 ** ( -0.5D+00 * n - 1.0D+00 ) * ca0 / g1
        vt = -0.5D+00 * n
        call gaih ( vt, g0 )
        cdn = cmplx ( g0, 0.0D+00, kind = 8 )
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

        do m = 1, 250
           vm = 0.5D+00 * ( m - n )
           call gaih ( vm, gm )
           cr = - cr * sq2 * z / m
           cdw = gm * cr
           cdn = cdn + cdw
           if ( abs ( cdw ) < abs ( cdn ) * eps ) then
              exit
           end if
        end do

        cdn = cb0 * cdn

     end if

  end if

  return
end subroutine cpdsa




  !*****************************************************************************80
  !
  !! CPSI computes the psi function for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, Y, the real and imaginary parts 
  !    of the argument.
  !
  !    Output, real(8) PSR, PSI, the real and imaginary parts
  !    of the function value.
  !
subroutine cpsi ( x, y, psr, psi )
  !CPSI computes the psi function for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the real part of the argument
  !@param y: in: the imaginary part of the argument
  !@param psr: out: the real part of the function value
  !@param psi: out: the imaginary part of the function value
  !
  implicit none

  real(8), save, dimension ( 8 ) :: a = (/ &
       -0.8333333333333D-01, 0.83333333333333333D-02, &
       -0.39682539682539683D-02, 0.41666666666666667D-02, &
       -0.75757575757575758D-02, 0.21092796092796093D-01, &
       -0.83333333333333333D-01, 0.4432598039215686D+00 /)
  real(8)    :: ct2
  integer(4) :: k
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: psi
  real(8)    :: psr
  real(8)    :: ri
  real(8)    :: rr
  real(8)    :: th
  real(8)    :: tm
  real(8)    :: tn
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x1
  real(8)    :: y
  real(8)    :: y1
  real(8)    :: z0
  real(8)    :: z2

  pi = 3.141592653589793D+00

  if ( y == 0.0D+00 .and. x == int ( x ) .and. x <= 0.0D+00 ) then

     psr = 1.0D+300
     psi = 0.0D+00

  else

     if ( x < 0.0D+00 ) then
        x1 = x
        y1 = y
        x = -x
        y = -y
     end if

     x0 = x

     if ( x < 8.0D+00 ) then
        n = 8 - int ( x )
        x0 = x + n
     end if

     if ( x0 == 0.0D+00 ) then
        if ( y /= 0.0D+00 ) then
           th = 0.5D+00 * pi
        else
           th = 0.0D+00
        end if
     else
        th = atan ( y / x0 )
     end if

     z2 = x0 * x0 + y * y
     z0 = sqrt ( z2 )
     psr = log ( z0 ) - 0.5D+00 * x0 / z2
     psi = th + 0.5D+00 * y / z2
     do k = 1, 8
        psr = psr + a(k) * z2 ** ( - k ) * cos ( 2.0D+00 * k * th )
        psi = psi - a(k) * z2 ** ( - k ) * sin ( 2.0D+00 * k * th )
     end do

     if ( x < 8.0D+00 ) then
        rr = 0.0D+00
        ri = 0.0D+00
        do k = 1, n
           rr = rr + ( x0 - k ) / ( ( x0 - k ) ** 2.0D+00 + y * y )
           ri = ri + y / ( ( x0 - k ) ** 2.0D+00 + y * y )
        end do
        psr = psr - rr
        psi = psi + ri
     end if

     if ( x1 < 0.0D+00 ) then
        tn = tan ( pi * x )
        tm = tanh ( pi * y )
        ct2 = tn * tn + tm * tm
        psr = psr + x / ( x * x + y * y ) + pi * ( tn - tn * tm * tm ) / ct2
        psi = psi - y / ( x * x + y * y ) - pi * tm * ( 1.0D+00 + tn * tn ) / ct2
        x = x1
        y = y1
     end if

  end if

  return
end subroutine cpsi




  !*****************************************************************************80
  !
  !! CSPHIK: complex modified spherical Bessel functions and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) N, the order of in(z) and kn(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CSI(0:N), CDI(0:N), CSK(0:N), CDK(0:N),
  !    the values of in(z), in'(z), kn(z), kn'(z).
  !
subroutine csphik ( n, z, nm, csi, cdi, csk, cdk )
  !CSPHIK: complex modified spherical Bessel functions and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of in(z) and kn(z)
  !@param z: in: the argument
  !@param nm: out: the highest order computed
  !@param csi: out: in(z)
  !@param cdi: out: in'(z)
  !@param csk: out: kn(z)
  !@param cdk: out: kn'(z)
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  complex(8) :: ccosh1
  complex(8) :: cdi(0:n)
  complex(8) :: cdk(0:n)
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: ci
  complex(8) :: cs
  complex(8) :: csi(0:n)
  complex(8) :: csi0
  complex(8) :: csi1
  complex(8) :: csinh1
  complex(8) :: csk(0:n)
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: pi
  complex(8) :: z

  pi = 3.141592653589793D+00
  a0 = abs ( z )    
  nm = n

  if ( a0 < 1.0D-60 ) then
     do k = 0, n
        csi(k) = 0.0D+00
        cdi(k) = 0.0D+00
        csk(k) = 1.0D+300
        cdk(k) = -1.0D+300
     end do
     csi(0) = 1.0D+00
     cdi(1) = 0.3333333333333333D+00
     return
  end if

  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  csinh1 = sin ( ci * z ) / ci
  ccosh1 = cos ( ci * z )
  csi0 = csinh1 / z
  csi1 = ( - csinh1 / z + ccosh1 ) / z
  csi(0) = csi0
  csi(1) = csi1

  if ( 2 <= n ) then

     m = msta1 ( a0, 200 )
     if ( m < n ) then
        nm = m
     else
        m = msta2 ( a0, n, 15 )
     end if

     cf0 = 0.0D+00
     cf1 = 1.0D+00-100
     do k = m, 0, -1
        cf = ( 2.0D+00 * k + 3.0D+00 ) * cf1 / z + cf0
        if ( k <= nm ) then
           csi(k) = cf
        end if
        cf0 = cf1
        cf1 = cf
     end do

     if ( abs ( csi0 ) <= abs ( csi1 ) ) then
        cs = csi1 / cf0
     else
        cs = csi0 / cf
     end if

     do k = 0, nm
        csi(k) = cs * csi(k)
     end do

  end if

  cdi(0) = csi(1)
  do k = 1, nm
     cdi(k) = csi(k-1) - ( k + 1.0D+00 ) * csi(k) / z
  end do

  csk(0) = 0.5D+00 * pi / z * exp ( - z )
  csk(1) = csk(0) * ( 1.0D+00 + 1.0D+00 / z )
  do k = 2, nm
     if ( abs ( csi(k-2) ) < abs ( csi(k-1) ) ) then
        csk(k) = ( 0.5D+00 * pi / ( z * z ) - csi(k) * csk(k-1) ) / csi(k-1)
     else
        csk(k) = ( csi(k) * csk(k-2) + ( k - 0.5D+00 ) * pi / z ** 3 ) / csi(k-2)
     end if
  end do

  cdk(0) = -csk(1)
  do k = 1, nm
     cdk(k) = - csk(k-1) - ( k + 1.0D+00 ) * csk(k) / z
  end do

  return
end subroutine csphik




  !*****************************************************************************80
  !
  !! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
  !
  !  Discussion:
  !
  !    This procedure computes spherical Bessel functions jn(z) and yn(z)
  !    and their derivatives for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) N, the order of jn(z) and yn(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, complex(8) CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
  !    the values of jn(z), jn'(z), yn(z), yn'(z).
  !
subroutine csphjy ( n, z, nm, csj, cdj, csy, cdy )
  !CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
  !
  !This procedure computes spherical Bessel functions jn(z) and yn(z)
  !and their derivatives for a complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of jn(z) and yn(z)
  !@param z: in: the argument
  !@param nm: out: the highest order computed
  !@param csj: jn(z), n=0,...,N
  !@param cdj: jn'(z), n=0,...,N
  !@param csy: yn(z), n=0,...,N
  !@param cdy: yn'(z), n=0,...,N
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  complex(8) :: csj(0:n)
  complex(8) :: cdj(0:n)
  complex(8) :: csy(0:n)
  complex(8) :: cdy(0:n)
  complex(8) :: cf
  complex(8) :: cf0
  complex(8) :: cf1
  complex(8) :: cs
  complex(8) :: csa
  complex(8) :: csb
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  complex(8) :: z

  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-60 ) then
     do k = 0, n
        csj(k) = 0.0D+00
        cdj(k) = 0.0D+00
        csy(k) = -1.0D+300
        cdy(k) = 1.0D+300
     end do
     csj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cdj(1) = cmplx ( 0.333333333333333D+00, 0.0D+00, kind = 8 )
     return
  end if

  csj(0) = sin ( z ) / z
  csj(1) = ( csj(0) - cos ( z ) ) / z

  if ( 2 <= n ) then
     csa = csj(0)
     csb = csj(1)
     m = msta1 ( a0, 200 )
     if ( m < n ) then
        nm = m
     else
        m = msta2 ( a0, n, 15 )
     end if
     cf0 = 0.0D+00
     cf1 = 1.0D+00-100
     do k = m, 0, -1
        cf = ( 2.0D+00 * k + 3.0D+00 ) * cf1 / z - cf0
        if ( k <= nm ) then
           csj(k) = cf
        end if
        cf0 = cf1
        cf1 = cf
     end do

     if ( abs ( csa ) <= abs ( csb ) ) then
        cs = csb / cf0
     else
        cs = csa / cf
     end if

     do k = 0, nm
        csj(k) = cs * csj(k)
     end do

  end if

  cdj(0) = ( cos ( z ) - sin ( z ) / z ) / z
  do k = 1, nm
     cdj(k) = csj(k-1) - ( k + 1.0D+00 ) * csj(k) / z
  end do
  csy(0) = - cos ( z ) / z
  csy(1) = ( csy(0) - sin ( z ) ) / z
  cdy(0) = ( sin ( z ) + cos ( z ) / z ) / z
  cdy(1) = ( 2.0D+00 * cdy(0) - cos ( z ) )  / z

  do k = 2, nm
     if ( abs ( csj(k-2) ) < abs ( csj(k-1) ) ) then 
        csy(k) = ( csj(k) * csy(k-1) - 1.0D+00 / ( z * z ) ) / csj(k-1)
     else
        csy(k) = ( csj(k) * csy(k-2) &
             - ( 2.0D+00 * k - 1.0D+00 ) / z ** 3 ) / csj(k-2)
     end if
  end do

  do k = 2, nm
     cdy(k) = csy(k-1) - ( k + 1.0D+00 ) * csy(k) / z
  end do

  return
end subroutine csphjy




  !*****************************************************************************80
  !
  !! CV0 computes the initial characteristic value of Mathieu functions.
  !
  !  Discussion:
  !
  !    This procedure computes the initial characteristic value of Mathieu 
  !    functions for m <= 12 or q <= 300 or q <= m*m.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !   03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code:
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the functions.
  !
  !    Input, real(8) Q, the parameter of the functions.
  !
  !    Output, real(8) A0, the characteristic value.
  !
subroutine cv0 ( kd, m, q, a0 )
  !CV0 computes the initial characteristic value of Mathieu functions.
  !
  !This procedure computes the initial characteristic value of Mathieu
  !functions for m <= 12 or q <= 300 or q <= m*m.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code: 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the functions
  !@param q: in: the parameter of the functions
  !@param a0: out: the characteristic value
  !
  implicit none

  real(8)    :: a0
  integer(4) :: kd
  integer(4) :: m
  real(8)    :: q
  real(8)    :: q2

  q2 = q * q

  if ( m == 0 ) then

     if ( q <= 1.0D+00 ) then

        a0 = ((( &
             0.0036392D+00   * q2 &
             - 0.0125868D+00 ) * q2 &
             + 0.0546875D+00 ) * q2 &
             - 0.5D+00 )       * q2

     else if ( q <= 10.0D+00 ) then

        a0 = (( &
             3.999267D-03   * q &
             - 9.638957D-02 ) * q &
             - 0.88297D+00 )  * q &
             + 0.5542818D+00 

     else

        call cvql ( kd, m, q, a0 )

     end if

  else if ( m == 1 ) then

     if ( q <= 1.0D+00 .and. kd == 2 ) then

        a0 = ((( &
             - 6.51D-04 * q &
             - 0.015625D+00 ) *  q &
             - 0.125D+00 ) * q &
             + 1.0D+00 ) * q &
             + 1.0D+00 

     else if ( q <= 1.0D+00 .and. kd == 3 ) then

        a0 = ((( &
             - 6.51D-04 * q &
             + 0.015625D+00 ) * q &
             - 0.125D+00 ) * q &
             - 1.0D+00 ) * q &
             + 1.0D+00

     else if ( q <= 10.0D+00 .and. kd == 2 ) then

        a0 = ((( &
             - 4.94603D-04 * q &
             + 1.92917D-02 ) * q &
             - 0.3089229D+00 ) * q &
             + 1.33372D+00 ) * q &
             + 0.811752D+00 

     else if ( q <= 10.0D+00 .and. kd == 3 ) then

        a0 = (( &
             1.971096D-03 * q &
             - 5.482465D-02 ) * q &
             - 1.152218D+00 ) * q &
             + 1.10427D+00 

     else

        call cvql ( kd, m, q, a0 )

     end if

  else if ( m == 2 ) then

     if ( q <= 1.0D+00 .and. kd == 1 ) then

        a0 = ((( &
             - 0.0036391D+00   * q2 &
             + 0.0125888D+00 ) * q2 &
             - 0.0551939D+00 ) * q2 &
             + 0.416667D+00 )  * q2 + 4.0D+00 

     else if ( q <= 1.0D+00 .and. kd == 4 ) then

        a0 = (  &
             0.0003617D+00 * q2  &
             - 0.0833333D+00 ) * q2 + 4.0D+00 

     else if ( q <= 15.0D+00 .and. kd == 1 ) then

        a0 = ((( &
             3.200972D-04    * q &
             - 8.667445D-03 )  * q &
             - 1.829032D-04 )  * q &
             + 0.9919999D+00 ) * q &
             + 3.3290504D+00 

     else if ( q <= 10.0D+00 .and. kd == 4 ) then

        a0 = (( &
             2.38446D-03 * q &
             - 0.08725329D+00 ) * q &
             - 4.732542D-03 ) * q &
             + 4.00909D+00 

     else

        call cvql ( kd, m, q, a0 )

     end if

  else if ( m == 3 ) then

     if ( q <= 1.0D+00 .and. kd == 2 ) then
        a0 = (( &
             6.348D-04 * q &
             + 0.015625D+00 ) * q &
             + 0.0625 ) * q2  &
             + 9.0D+00 
     else if ( q <= 1.0D+00 .and. kd == 3 ) then
        a0 = (( &
             6.348D-04 * q &
             - 0.015625D+00 ) * q &
             + 0.0625D+00 ) * q2 &
             + 9.0D+00 
     else if ( q <= 20.0D+00 .and. kd == 2 ) then
        a0 = ((( &
             3.035731D-04 * q &
             - 1.453021D-02 ) * q &
             + 0.19069602D+00 ) * q &
             - 0.1039356D+00 ) * q &
             + 8.9449274D+00 
     else if ( q <= 15.0D+00 .and. kd == 3 ) then
        a0 = (( &
             9.369364D-05 * q &
             - 0.03569325D+00 ) * q &
             + 0.2689874D+00 ) * q &
             + 8.771735D+00 
     else
        call cvql ( kd, m, q, a0 )
     end if

  else if ( m == 4 ) then

     if ( q <= 1.0D+00 .and. kd == 1 ) then
        a0 = (( &
             - 2.1D-06 * q2 &
             + 5.012D-04 ) * q2 &
             + 0.0333333 ) * q2 &
             + 16.0D+00
     else if ( q <= 1.0D+00 .and. kd == 4 ) then
        a0 = (( &
             3.7D-06 * q2 &
             - 3.669D-04 ) * q2 &
             + 0.0333333D+00 ) * q2 &
             + 16.0D+00
     else if ( q <= 25.0D+00 .and. kd == 1 ) then
        a0 = ((( &
             1.076676D-04 * q &
             - 7.9684875D-03 ) * q &
             + 0.17344854D+00 ) * q &
             - 0.5924058D+00 ) * q &
             + 16.620847D+00
     else if ( q <= 20.0D+00 .and. kd == 4 ) then
        a0 = (( &
             - 7.08719D-04 * q &
             + 3.8216144D-03 ) * q &
             + 0.1907493D+00 ) * q &
             + 15.744D+00
     else
        call cvql ( kd, m, q, a0 )
     end if

  else if ( m == 5 ) then

     if ( q <= 1.0D+00 .and. kd == 2 ) then
        a0 = (( &
             6.8D-6 * q &
             + 1.42D-05 ) * q2 &
             + 0.0208333D+00 ) * q2 &
             + 25.0D+00
     else if ( q <= 1.0D+00 .and. kd == 3 ) then
        a0 = (( &
             - 6.8D-06 * q &
             + 1.42D-05 ) * q2 &
             + 0.0208333D+00 ) * q2 &
             + 25.0D+00
     else if ( q <= 35.0D+00 .and. kd == 2 ) then
        a0 = ((( &
             2.238231D-05 * q &
             - 2.983416D-03 ) * q &
             + 0.10706975D+00 ) * q &
             - 0.600205D+00 ) * q &
             + 25.93515D+00
     else if ( q <= 25.0D+00 .and. kd == 3 ) then
        a0 = (( &
             - 7.425364D-04 * q &
             + 2.18225D-02 ) * q &
             + 4.16399D-02 ) * q &
             + 24.897D+00
     else
        call cvql ( kd, m, q, a0 )
     end if

  else if ( m == 6 ) then

     if ( q <= 1.0D+00 ) then
        a0 = ( 0.4D-06 * q2 + 0.0142857 ) * q2 + 36.0D+00
     else if ( q <= 40.0D+00 .and. kd == 1 ) then
        a0 = ((( &
             - 1.66846D-05 * q &
             + 4.80263D-04 ) * q &
             + 2.53998D-02 ) * q &
             - 0.181233D+00 ) * q  &
             + 36.423D+00
     else if ( q <= 35.0D+00 .and. kd == 4 ) then
        a0 = (( &
             - 4.57146D-04 * q &
             + 2.16609D-02 ) * q &
             - 2.349616D-02 ) * q &
             + 35.99251D+00
     else
        call cvql ( kd, m, q, a0 )
     end if

  else if ( m == 7 ) then

     if ( q <= 10.0D+00 ) then
        call cvqm ( m, q, a0 )
     else if ( q <= 50.0D+00 .and. kd == 2 ) then
        a0 = ((( &
             - 1.411114D-05 * q &
             + 9.730514D-04 ) * q &
             - 3.097887D-03 ) * q &
             + 3.533597D-02 ) * q &
             + 49.0547D+00
     else if ( q <= 40.0D+00 .and. kd == 3 ) then
        a0 = (( &
             - 3.043872D-04 * q &
             + 2.05511D-02 ) * q &
             - 9.16292D-02 ) * q &
             + 49.19035D+00
     else
        call cvql ( kd, m, q, a0 )
     end if

  else if ( 8 <= m ) then

     if ( q <= 3.0D+00 * m ) then
        call cvqm ( m, q, a0 )
     else if ( m * m .lt. q ) then
        call cvql ( kd, m, q, a0 )
     else if ( m == 8 .and. kd == 1 ) then
        a0 = ((( &
             8.634308D-06 * q &
             - 2.100289D-03 ) * q &
             + 0.169072D+00 ) * q &
             - 4.64336D+00 ) * q &
             + 109.4211D+00
     else if ( m == 8 .and. kd == 4 ) then
        a0 = (( &
             - 6.7842D-05 * q &
             + 2.2057D-03 ) * q &
             + 0.48296D+00 ) * q &
             + 56.59D+00
     else if ( m == 9 .and. kd == 2 ) then
        a0 = ((( &
             2.906435D-06 * q &
             - 1.019893D-03 ) * q &
             + 0.1101965D+00 ) * q &
             - 3.821851D+00 ) * q &
             + 127.6098D+00
     else if ( m == 9 .and. kd == 3 ) then
        a0 = (( &
             - 9.577289D-05 * q &
             + 0.01043839D+00 ) * q &
             + 0.06588934D+00 ) * q &
             + 78.0198D+00
     else if ( m == 10 .and. kd == 1 ) then
        a0 = ((( &
             5.44927D-07 * q &
             - 3.926119D-04 ) * q &
             + 0.0612099D+00 ) * q &
             - 2.600805D+00 ) * q &
             + 138.1923D+00
     else if ( m == 10 .and. kd == 4 ) then
        a0 = (( &
             - 7.660143D-05 * q &
             + 0.01132506D+00 ) * q &
             - 0.09746023D+00 ) * q &
             + 99.29494D+00
     else if ( m == 11 .and. kd == 2 ) then
        a0 = ((( &
             - 5.67615D-07 * q &
             + 7.152722D-06 ) * q &
             + 0.01920291D+00 ) * q &
             - 1.081583D+00 ) * q &
             + 140.88D+00
     else if ( m == 11 .and. kd == 3 ) then
        a0 = (( &
             - 6.310551D-05 * q &
             + 0.0119247D+00 ) * q &
             - 0.2681195D+00 ) * q &
             + 123.667D+00 
     else if ( m == 12 .and. kd == 1 ) then
        a0 = ((( &
             - 2.38351D-07 * q &
             - 2.90139D-05 ) * q &
             + 0.02023088D+00 ) * q &
             - 1.289D+00 ) * q &
             + 171.2723D+00
     else if ( m == 12 .and. kd == 4 ) then
        a0 = ((( &
             3.08902D-07 * q &
             - 1.577869D-04 ) * q &
             + 0.0247911D+00 ) * q &
             - 1.05454D+00 ) * q  &
             + 161.471D+00

     end if

  end if

  return
end subroutine cv0




  !*****************************************************************************80
  !
  !! CVA1 computes a sequence of characteristic values of Mathieu functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code.
  !    1, for cem(x,q)  ( m = 0,2,4,... )
  !    2, for cem(x,q)  ( m = 1,3,5,... )
  !    3, for sem(x,q)  ( m = 1,3,5,... )
  !    4, for sem(x,q)  ( m = 2,4,6,... )
  !
  !    Input, integer(4) M, the maximum order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Output, real(8) CV(*), characteristic values.
  !    For KD = 1, CV(1), CV(2), CV(3),..., correspond to
  !    the characteristic values of cem for m = 0,2,4,...
  !    For KD = 2, CV(1), CV(2), CV(3),..., correspond to
  !    the characteristic values of cem for m = 1,3,5,...
  !    For KD = 3, CV(1), CV(2), CV(3),..., correspond to
  !    the characteristic values of sem for m = 1,3,5,...
  !    For KD = 4, CV(1), CV(2), CV(3),..., correspond to
  !    the characteristic values of sem for m = 0,2,4,...
  !       
subroutine cva1 ( kd, m, q, cv )
  !CVA1 computes a sequence of characteristic values of Mathieu functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code. 1, for cem(x,q) ( m = 0,2,4,... ) 2, for cem(x,q) ( m = 1,3,5,... ) 3, for sem(x,q) ( m =
  !  1,3,5,... ) 4, for sem(x,q) ( m = 2,4,6,... )
  !* CV: characteristic values. For KD = 1, CV(1), CV(2), CV(3), ..., correspond to the characteristic values of cem for
  !  m = 0, 2, 4, ... For KD = 2, CV(1), CV(2), CV(3), ..., correspond to the characteristic values of cem for m = 1, 3,
  !  5, ... For KD = 3, CV(1), CV(2), CV(3), ..., correspond to the characteristic values of sem for m = 1, 3, 5, ...
  !  For KD = 4, CV(1), CV(2), CV(3), ..., correspond to the characteristic values of sem for m = 0, 2, 4, ...
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the maximum order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param cv: out: characteristic values: for KD=1,2,3,4 the cem (m even, odd) and sem (m odd, even) series
  !
  implicit none

  real(8)    :: cv(200)
  real(8)    :: d(500)
  real(8)    :: e(500)
  real(8)    :: eps
  real(8)    :: f(500)
  real(8)    :: g(200)
  real(8)    :: h(200)
  integer(4) :: i
  integer(4) :: ic
  integer(4) :: icm
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k1
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: nm
  integer(4) :: nm1
  real(8)    :: q
  real(8)    :: s
  real(8)    :: t
  real(8)    :: t1
  real(8)    :: x1
  real(8)    :: xa
  real(8)    :: xb

  eps = 1.0D-14

  if ( kd == 4 ) then
     icm = m / 2
  else
     icm = int ( m / 2 ) + 1
  end if

  if ( q == 0.0D+00 ) then

     if ( kd == 1 ) then
        do ic = 1, icm
           cv(ic) = 4.0D+00 * ( ic - 1.0D+00 ) ** 2
        end do
     else if ( kd /= 4 ) then
        do ic = 1, icm
           cv(ic) = ( 2.0D+00 * ic - 1.0D+00 ) ** 2
        end do
     else
        do ic = 1, icm
           cv(ic) = 4.0D+00 * ic * ic
        end do
     end if

  else

     nm = int ( 10D+00 + 1.5D+00 * m + 0.5D+00 * q )
     e(1) = 0.0D+00
     f(1) = 0.0D+00

     if ( kd == 1 ) then

        d(1) = 0.0D+00
        do i = 2, nm
           d(i) = 4.0D+00 * ( i - 1.0D+00 ) ** 2
           e(i) = q
           f(i) = q * q
        end do
        e(2) = sqrt ( 2.0D+00 ) * q
        f(2) = 2.0D+00 * q * q

     else if ( kd /= 4 ) then

        d(1) = 1.0D+00 + ( -1.0D+00 ) ** kd * q
        do i = 2, nm
           d(i) = ( 2.0D+00 * i - 1.0D+00 ) ** 2
           e(i) = q
           f(i) = q * q
        end do

     else

        d(1) = 4.0D+00
        do i = 2, nm
           d(i) = 4.0D+00 * i * i
           e(i) = q
           f(i) = q * q
        end do

     end if

     xa = d(nm) + abs ( e(nm) )
     xb = d(nm) - abs ( e(nm) )

     nm1 = nm - 1
     do i = 1, nm1
        t = abs ( e(i) ) + abs ( e(i+1) )
        t1 = d(i) + t
        xa = max ( xa, t1 )
        t1 = d(i) - t
        xb = min ( xb, t1 )
     end do

     do i = 1, icm
        g(i) = xa
        h(i) = xb
     end do

     do k = 1, icm

        do k1 = k, icm
           if ( g(k1) < g(k) ) then
              g(k) = g(k1)
              exit
           end if
        end do

        if ( k /= 1 .and. h(k) < h(k-1) ) then
           h(k) = h(k-1)
        end if

        do

           x1 = ( g(k) + h(k) ) /2.0D+00
           cv(k) = x1

           if ( abs ( ( g(k) - h(k) ) / x1 ) < eps ) then
              exit
           end if

           j = 0
           s = 1.0D+00
           do i = 1, nm
              if ( s == 0.0D+00 ) then
                 s = s + 1.0D-30
              end if
              t = f(i) / s
              s = d(i) - t - x1
              if ( s < 0.0D+00 ) then
                 j = j + 1
              end if
           end do

           if ( j < k ) then
              h(k) = x1
           else
              g(k) = x1
              if ( icm <= j ) then
                 g(icm) = x1
              else
                 h(j+1) = max ( h(j+1), x1 )
                 g(j) = min ( g(j), x1 )
              end if
           end if

        end do

        cv(k) = x1

     end do

  end if

  return
end subroutine cva1




  !*****************************************************************************80
  !
  !! CVA2 computes a specific characteristic value of Mathieu functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code:
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Output, real(8) A, the characteristic value.
  !
subroutine cva2 ( kd, m, q, a )
  !CVA2 computes a specific characteristic value of Mathieu functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code: 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param a: out: the characteristic value
  !
  implicit none

  real(8)    :: a
  real(8)    :: a1
  real(8)    :: a2
  real(8)    :: delta
  integer(4) :: i
  integer(4) :: iflag
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: ndiv
  integer(4) :: nn
  real(8)    :: q
  real(8)    :: q1
  real(8)    :: q2
  real(8)    :: qq

  if ( m <= 12 .or. q <= 3.0D+00 * m .or. m * m < q ) then

     call cv0 ( kd, m, q, a )

     if ( q /= 0.0D+00 ) then
        call refine ( kd, m, q, a, 1 )
     end if

  else

     ndiv = 10
     delta = ( m - 3.0D+00 ) * m / real ( ndiv, kind = 8 )

     if ( ( q - 3.0D+00 * m ) <= ( m * m - q ) ) then

        do

           nn = int ( ( q - 3.0D+00 * m ) / delta ) + 1
           delta = ( q - 3.0D+00 * m ) / nn
           q1 = 2.0D+00 * m
           call cvqm ( m, q1, a1 )
           q2 = 3.0D+00 * m
           call cvqm ( m, q2, a2 )
           qq = 3.0D+00 * m

           do i = 1, nn

              qq = qq + delta
              a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )

              if ( i == nn ) then
                 iflag = -1
              else
                 iflag = 1
              end if

              call refine ( kd, m, qq, a, iflag )
              q1 = q2
              q2 = qq
              a1 = a2
              a2 = a

           end do

           if ( iflag /= -10 ) then
              exit
           end if

           ndiv = ndiv * 2
           delta = ( m - 3.0D+00 ) * m / real ( ndiv, kind = 8 )

        end do

     else

        do

           nn = int ( ( m * m - q ) / delta ) + 1
           delta = ( m * m - q ) / nn
           q1 = m * ( m - 1.0D+00 )
           call cvql ( kd, m, q1, a1 )
           q2 = m * m
           call cvql ( kd, m, q2, a2 )
           qq = m * m

           do i = 1, nn

              qq = qq - delta
              a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )

              if ( i == nn ) then
                 iflag = -1
              else
                 iflag = 1
              end if

              call refine ( kd, m, qq, a, iflag )
              q1 = q2
              q2 = qq
              a1 = a2
              a2 = a

           end do

           if ( iflag /= -10 ) then
              exit
           end if

           ndiv = ndiv * 2
           delta = ( m - 3.0D+00 ) * m / real ( ndiv, kind = 8 )

        end do

     end if

  end if

  return
end subroutine cva2




  !*****************************************************************************80
  !
  !! CVF computes F for the characteristic equation of Mathieu functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code:
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Input, real(8) A, the characteristic value.
  !
  !    Input, integer(4) MJ, ?
  !
  !    Output, real(8) F, the value of the function for the
  !    characteristic equation.
  !
subroutine cvf ( kd, m, q, a, mj, f )
  !CVF computes F for the characteristic equation of Mathieu functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code: 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param a: in: the characteristic value
  !@param mj: in: not documented in the original header
  !@param f: out: the value of the function for the characteristic equation
  !
  implicit none

  real(8)    :: a
  real(8)    :: b
  real(8)    :: f
  integer(4) :: ic
  integer(4) :: j
  integer(4) :: j0
  integer(4) :: jf
  integer(4) :: kd
  integer(4) :: l
  integer(4) :: l0
  integer(4) :: m
  integer(4) :: mj
  real(8)    :: q
  real(8)    :: t0
  real(8)    :: t1
  real(8)    :: t2

  b = a
  ic = int ( m / 2 )
  l = 0
  l0 = 0
  j0 = 2
  jf = ic

  if ( kd == 1 ) then
     l0 = 2
     j0 = 3
  else if ( kd == 2 .or. kd == 3 ) then
     l = 1
  else if ( kd == 4 ) then
     jf = ic - 1
  end if

  t1 = 0.0D+00
  do j = mj, ic + 1, -1
     t1 = - q * q / ( ( 2.0D+00 * j + l ) ** 2 - b + t1 )
  end do

  if ( m <= 2 ) then

     t2 = 0.0D+00

     if ( kd == 1 ) then
        if ( m == 0 ) then
           t1 = t1 + t1
        else if ( m == 2 ) then
           t1 = - 2.0D+00 * q * q / ( 4.0D+00 - b + t1 ) - 4.0D+00
        end if
     else if ( kd == 2 ) then
        if ( m == 1 ) then
           t1 = t1 + q
        end if
     else if ( kd == 3 ) then
        if ( m == 1 ) then
           t1 = t1 - q
        end if
     end if

  else

     if ( kd == 1 ) then
        t0 = 4.0D+00 - b + 2.0D+00 * q * q / b
     else if ( kd == 2 ) then
        t0 = 1.0D+00 - b + q
     else if ( kd == 3 ) then
        t0 = 1.0D+00 - b - q
     else if ( kd == 4 ) then
        t0 = 4.0D+00 - b
     end if

     t2 = - q * q / t0
     do j = j0, jf
        t2 = - q * q / ( ( 2.0D+00 * j - l - l0 ) ** 2 - b + t2 )
     end do

  end if

  f = ( 2.0D+00 * ic + l ) ** 2 + t1 + t2 - b

  return
end subroutine cvf




  !*****************************************************************************80
  !
  !! CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    10 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code:
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter value.
  !
  !    Output, real(8) A0, the initial characteristic value.
  !
subroutine cvql ( kd, m, q, a0 )
  !CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code: 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter value
  !@param a0: out: the initial characteristic value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: c1
  real(8)    :: cv1
  real(8)    :: cv2
  real(8)    :: d1
  real(8)    :: d2
  real(8)    :: d3
  real(8)    :: d4
  integer(4) :: kd
  integer(4) :: m
  real(8)    :: p1
  real(8)    :: p2
  real(8)    :: q
  real(8)    :: w
  real(8)    :: w2
  real(8)    :: w3
  real(8)    :: w4
  real(8)    :: w6

  if ( kd == 1 .or. kd == 2 ) then
     w = 2.0D+00 * m + 1.0D+00
  else
     w = 2.0D+00 * m - 1.0D+00
  end if

  w2 = w * w
  w3 = w * w2
  w4 = w2 * w2
  w6 = w2 * w4
  d1 = 5.0D+00 + 34.0D+00 / w2 + 9.0D+00 / w4
  d2 = ( 33.0D+00 + 410.0D+00 / w2 + 405.0D+00 / w4 ) / w
  d3 = ( 63.0D+00 + 1260.0D+00 / w2 + 2943.0D+00 / w4 + 486.0D+00 / w6 ) / w2
  d4 = ( 527.0D+00 + 15617.0D+00 / w2 + 69001.0D+00 / w4 &
       + 41607.0D+00 / w6 ) / w3
  c1 = 128.0D+00
  p2 = q / w4
  p1 = sqrt ( p2 )
  cv1 = - 2.0D+00 * q + 2.0D+00 * w * sqrt ( q ) &
       - ( w2 + 1.0D+00 ) / 8.0D+00
  cv2 = ( w + 3.0D+00 / w ) + d1 / ( 32.0D+00 * p1 ) + d2 &
       / ( 8.0D+00 * c1 * p2 )
  cv2 = cv2 + d3 / ( 64.0D+00 * c1 * p1 * p2 ) + d4 &
       / ( 16.0D+00 * c1 * c1 * p2 * p2 )
  a0 = cv1 - cv2 / ( c1 * p1 )

  return
end subroutine cvql




  !*****************************************************************************80
  !
  !! CVQM computes the characteristic value of Mathieu functions for q <= m*m.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter value.
  !
  !    Output, real(8) A0, the initial characteristic value.
  !
subroutine cvqm ( m, q, a0 )
  !CVQM computes the characteristic value of Mathieu functions for q <= m*m.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter value
  !@param a0: out: the initial characteristic value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: hm1
  real(8)    :: hm3
  real(8)    :: hm5
  integer(4) :: m
  real(8)    :: q

  hm1 = 0.5D+00 * q / ( m * m - 1.0D+00 )
  hm3 = 0.25D+00 * hm1 ** 3 / ( m * m - 4.0D+00 )
  hm5 = hm1 * hm3 * q / ( ( m * m - 1.0D+00 ) * ( m * m - 9.0D+00 ) )
  a0 = m * m + q * ( hm1 + ( 5.0D+00 * m * m + 7.0D+00 ) * hm3 &
       + ( 9.0D+00 * m ** 4 + 58.0D+00 * m * m + 29.0D+00 ) * hm5 )

  return
end subroutine cvqm




  !*****************************************************************************80
  !
  !! CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer KF, the function choice.
  !    0 for ZF = Y0(z) and ZD = Y0'(z);
  !    1 for ZF = Y1(z) and ZD = Y1'(z);
  !    2 for ZF = Y1'(z) and ZD = Y1''(z).
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) ZF, ZD, the values of the requested function 
  !    and derivative.
  !
subroutine cy01 ( kf, z, zf, zd )
  !CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KF: the function choice. 0 for ZF = Y0(z) and ZD = Y0'(z); 1 for ZF = Y1(z) and ZD = Y1'(z); 2 for ZF = Y1'(z) and
  !  ZD = Y1''(z)
  !
  !@param kf: in: function choice: 0 for Y0, 1 for Y1, 2 for Y1' (ZF the function, ZD its derivative)
  !@param z: in: the argument
  !@param zf: out: the requested function
  !@param zd: out: derivative
  !
  implicit none

  real(8), save, dimension(12) :: a = (/ &
       -0.703125D-01, 0.112152099609375D+00, &
       -0.5725014209747314D+00, 0.6074042001273483D+01, &
       -0.1100171402692467D+03, 0.3038090510922384D+04, &
       -0.1188384262567832D+06, 0.6252951493434797D+07, &
       -0.4259392165047669D+09, 0.3646840080706556D+11, &
       -0.3833534661393944D+13, 0.4854014686852901D+15 /)
  real(8)    :: a0
  real(8), save, dimension(12) :: a1 = (/ &
       0.1171875D+00, -0.144195556640625D+00, &
       0.6765925884246826D+00, -0.6883914268109947D+01, &
       0.1215978918765359D+03, -0.3302272294480852D+04, &
       0.1276412726461746D+06, -0.6656367718817688D+07, &
       0.4502786003050393D+09, -0.3833857520742790D+11, &
       0.4011838599133198D+13, -0.5060568503314727D+15 /)
  real(8), save, dimension(12) :: b = (/ &
       0.732421875D-01, -0.2271080017089844D+00, &
       0.1727727502584457D+01, -0.2438052969955606D+02, &
       0.5513358961220206D+03, -0.1825775547429318D+05, &
       0.8328593040162893D+06, -0.5006958953198893D+08, &
       0.3836255180230433D+10, -0.3649010818849833D+12, &
       0.4218971570284096D+14, -0.5827244631566907D+16 /)
  real(8), save, dimension(12) :: b1 = (/ &
       -0.1025390625D+00, 0.2775764465332031D+00, &
       -0.1993531733751297D+01, 0.2724882731126854D+02, &
       -0.6038440767050702D+03, 0.1971837591223663D+05, &
       -0.8902978767070678D+06, 0.5310411010968522D+08, &
       -0.4043620325107754D+10, 0.3827011346598605D+12, &
       -0.4406481417852278D+14, 0.6065091351222699D+16 /)
  complex(8) :: cbj0
  complex(8) :: cbj1
  complex(8) :: cby0
  complex(8) :: cby1
  complex(8) :: cdy0
  complex(8) :: cdy1
  complex(8) :: ci
  complex(8) :: cp
  complex(8) :: cp0
  complex(8) :: cp1
  complex(8) :: cq0
  complex(8) :: cq1
  complex(8) :: cr
  complex(8) :: cs
  complex(8) :: ct1
  complex(8) :: ct2
  complex(8) :: cu
  real(8)    :: el
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: kf
  real(8)    :: pi
  real(8)    :: rp2
  real(8)    :: w0
  real(8)    :: w1
  complex(8) :: z
  complex(8) :: z1
  complex(8) :: z2
  complex(8) :: zd
  complex(8) :: zf

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  rp2 = 2.0D+00 / pi
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0D+00 ) then

     cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cbj1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     cby0 = cmplx ( -1.0D+30, 0.0D+00, kind = 8 )
     cby1 = cmplx ( -1.0D+30, 0.0D+00, kind = 8 )
     cdy0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
     cdy1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )

  else

     if ( real ( z, kind = 8 ) < 0.0D+00) then
        z1 = -z
     end if

     if ( a0 <= 12.0D+00 ) then

        cbj0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           cr = - 0.25D+00 * cr * z2 / ( k * k )
           cbj0 = cbj0 + cr
           if ( abs ( cr ) < abs ( cbj0 ) * 1.0D-15 ) then
              exit
           end if
        end do

        cbj1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           cr = -0.25D+00 * cr * z2 / ( k * ( k + 1.0D+00 ) )
           cbj1 = cbj1 + cr
           if ( abs ( cr ) < abs ( cbj1 ) * 1.0D-15 ) then
              exit
           end if
        end do

        cbj1 = 0.5D+00 * z1 * cbj1
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           w0 = w0 + 1.0D+00 / k
           cr = -0.25D+00 * cr / ( k * k ) * z2
           cp = cr * w0
           cs = cs + cp
           if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
              exit
           end if
        end do

        cby0 = rp2 * ( log ( z1 / 2.0D+00 ) + el ) * cbj0 - rp2 * cs
        w1 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        cs = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 40
           w1 = w1 + 1.0D+00 / k
           cr = - 0.25D+00 * cr / ( k * ( k + 1 ) ) * z2
           cp = cr * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
           cs = cs + cp
           if ( abs ( cp ) < abs ( cs ) * 1.0D-15 ) then
              exit
           end if
        end do

        cby1 = rp2 * ( ( log ( z1 / 2.0D+00 ) + el ) * cbj1 &
             - 1.0D+00 / z1 - 0.25D+00 * z1 * cs )

     else

        if ( a0 < 35.0D+00 ) then
           k0 = 12
        else if ( a0 < 50.0D+00 ) then
           k0 = 10
        else
           k0 = 8
        end if

        ct1 = z1 - 0.25D+00 * pi
        cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, k0
           cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
        end do
        cq0 = -0.125D+00 / z1
        do k = 1, k0
           cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
        end do
        cu = sqrt ( rp2 / z1 )
        cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
        cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
        ct2 = z1 - 0.75D+00 * pi
        cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, k0
           cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
        end do
        cq1 = 0.375D+00 / z1
        do k = 1, k0
           cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
        end do
        cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
        cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )

     end if

     if ( real ( z, kind = 8 ) < 0.0D+00 ) then

        if ( imag ( z ) < 0.0D+00 ) then
           cby0 = cby0 - 2.0D+00 * ci * cbj0
        else
           cby0 = cby0 + 2.0D+00 * ci * cbj0
        end if

        if ( imag ( z ) < 0.0D+00 ) then
           cby1 = - ( cby1 - 2.0D+00 * ci * cbj1 )
        else
           cby1 = - ( cby1 + 2.0D+00 * ci * cbj1 )
        end if
        cbj1 = - cbj1

     end if

     cdy0 = - cby1
     cdy1 = cby0 - 1.0D+00 / z * cby1

  end if

  if ( kf == 0 ) then
     zf = cby0
     zd = cdy0
  else if ( kf == 1 ) then
     zf = cby1
     zd = cdy1
  else if ( kf == 2 ) then
     zf = cdy1
     zd = - cdy1 / z - ( 1.0D+00 - 1.0D+00 / ( z * z ) ) * cby1
  end if

  return
end subroutine cy01




  !*****************************************************************************80
  !
  !! CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
  !
  !  Parameters:
  !
  !    Ths procedure computes the complex zeros of Y0(z), Y1(z) and Y1'(z), 
  !    and their associated values at the zeros using the modified Newton's 
  !    iteration method.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Input, integer(4) KF, the function choice.
  !    0 for Y0(z) and Y1(z0);
  !    1 for Y1(z) and Y0(z1);
  !    2 for Y1'(z) and Y1(z1').
  !
  !    Input, integer(4) KC, complex/real choice.
  !    0, for complex roots;
  !    1, for real roots.
  !
  !    Output, real(8) ZO(NT), ZV(NT), the zeros of Y0(z) or Y1(z) 
  !    or Y1'(z), and the value of Y0'(z) or Y1'(z) or Y1(z) at the L-th zero.
  !
subroutine cyzo ( nt, kf, kc, zo, zv )
  !CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KF: the function choice. 0 for Y0(z) and Y1(z0); 1 for Y1(z) and Y0(z1); 2 for Y1'(z) and Y1(z1')
  !* ZO, ZV: the zeros of Y0(z) or Y1(z) or Y1'(z), and the value of Y0'(z) or Y1'(z) or Y1(z) at the L-th zero
  !
  !@param nt: in: the number of zeros
  !@param kf: in: function choice: 0 for Y0(z) and Y1(z0), 1 for Y1(z) and Y0(z1), 2 for Y1'(z) and Y1(z1')
  !@param kc: in: complex/real choice. 0, for complex roots; 1, for real roots
  !@param zo: out: the zeros of Y0(z), Y1(z) or Y1'(z)
  !@param zv: out: the value of Y0'(z), Y1'(z) or Y1(z) at the L-th zero
  !
  implicit none

  integer(4) :: nt

  real(8)    :: h
  integer(4) :: i
  integer(4) :: it
  integer(4) :: j
  integer(4) :: kc
  integer(4) :: kf
  integer(4) :: nr
  real(8)    :: w
  real(8)    :: w0
  real(8)    :: x
  real(8)    :: y
  complex(8) :: z
  complex(8) :: zd
  complex(8) :: zero
  complex(8) :: zf
  complex(8) :: zfd
  complex(8) :: zgd
  complex(8) :: zo(nt)
  complex(8) :: zp
  complex(8) :: zq
  complex(8) :: zv(nt)
  complex(8) :: zw

  if ( kc == 0 ) then
     x = -2.4D+00
     y = 0.54D+00
     h = 3.14D+00
  else if ( kc == 1 ) then
     x = 0.89D+00
     y = 0.0D+00
     h = -3.14D+00
  end if

  if ( kf == 1 ) then
     x = -0.503D+00
  else if ( kf == 2 ) then
     x = 0.577D+00
  end if

  zero = cmplx ( x, y, kind = 8 )

  do nr = 1, nt

     if ( nr == 1 ) then
        z = zero
     else
        z = zo(nr-1) - h
     end if

     it = 0

     do

        it = it + 1
        call cy01 ( kf, z, zf, zd )

        zp = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do i = 1, nr - 1
           zp = zp * ( z - zo(i) )
        end do

        zfd = zf / zp

        zq = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        do i = 1, nr - 1
           zw = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           do j = 1, nr - 1
              if ( j /= i ) then
                 zw = zw * ( z - zo(j) )
              end if
           end do
           zq = zq + zw
        end do

        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = abs ( z )

        if ( 50 < it .or. abs ( ( w - w0 ) / w ) <= 1.0D-12 ) then
           exit
        end if

     end do

     zo(nr) = z

  end do

  do i = 1, nt
     z = zo(i)
     if ( kf == 0 .or. kf == 2 ) then
        call cy01 ( 1, z, zf, zd )
        zv(i) = zf
     else if ( kf == 1 ) then
        call cy01 ( 0, z, zf, zd )
        zv(i) = zf
     end if
  end do

  return
end subroutine cyzo




  !*****************************************************************************80
  !
  !! DVLA computes parabolic cylinder functions Dv(x) for large argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    06 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) VA, the order.
  !
  !    Output, real(8) PD, the function value.
  !
subroutine dvla ( va, x, pd )
  !DVLA computes parabolic cylinder functions Dv(x) for large argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param va: in: the order
  !@param x: in: the argument
  !@param pd: out: the function value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: ep
  real(8)    :: eps
  real(8)    :: gl
  integer(4) :: k
  real(8)    :: pd
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: va
  real(8)    :: vl
  real(8)    :: x
  real(8)    :: x1

  pi = 3.141592653589793D+00
  eps = 1.0D-12
  ep = exp ( -0.25D+00 * x * x )
  a0 = abs ( x ) ** va * ep
  r = 1.0D+00
  pd = 1.0D+00
  do k = 1, 16
     r = -0.5D+00 * r * ( 2.0D+00 * k - va - 1.0D+00 ) &
          * ( 2.0D+00 * k - va - 2.0D+00 ) / ( k * x * x )
     pd = pd + r
     if ( abs ( r / pd ) < eps ) then
        exit
     end if
  end do

  pd = a0 * pd

  if ( x < 0.0D+00 ) then
     x1 = - x
     call vvla ( va, x1, vl )
     call gammaf ( -va, gl )
     pd = pi * vl / gl + cos ( pi * va ) * pd
  end if

  return
end subroutine dvla




  !*****************************************************************************80
  !
  !! DVSA computes parabolic cylinder functions Dv(x) for small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) VA, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PD, the function value.
  !
subroutine dvsa ( va, x, pd )
  !DVSA computes parabolic cylinder functions Dv(x) for small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param va: in: the order
  !@param x: in: the argument
  !@param pd: out: the function value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: ep
  real(8)    :: eps
  real(8)    :: g0
  real(8)    :: g1
  real(8)    :: ga0
  real(8)    :: gm
  integer(4) :: m
  real(8)    :: pd
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: sq2
  real(8)    :: va
  real(8)    :: va0
  real(8)    :: vm
  real(8)    :: vt
  real(8)    :: x

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  sq2 = sqrt ( 2.0D+00 )
  ep = exp ( -0.25D+00 * x * x )
  va0 = 0.5D+00 * ( 1.0D+00 - va )

  if ( va == 0.0D+00 ) then

     pd = ep

  else

     if ( x == 0.0D+00 ) then
        if ( va0 <= 0.0D+00 .and. va0 == int ( va0 ) ) then
           pd = 0.0D+00
        else
           call gammaf ( va0, ga0 )
           pd = sqrt ( pi ) / ( 2.0D+00 ** ( -0.5D+00 * va ) * ga0 )
        end if

     else

        call gammaf ( -va, g1 )
        a0 = 2.0D+00 ** ( -0.5D+00 * va - 1.0D+00 ) * ep / g1
        vt = -0.5D+00 * va
        call gammaf ( vt, g0 )
        pd = g0
        r = 1.0D+00
        do m = 1, 250
           vm = 0.5D+00 * ( m - va )
           call gammaf ( vm, gm )
           r = -r * sq2 * x / m
           r1 = gm * r
           pd = pd + r1
           if ( abs ( r1 ) < abs ( pd ) * eps ) then
              exit
           end if
        end do

        pd = a0 * pd

     end if

  end if

  return
end subroutine dvsa




  !*****************************************************************************80
  !
  !! E1XA computes the exponential integral E1(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    06 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) E1, the function value.
  !
subroutine e1xa ( x, e1 )
  !E1XA computes the exponential integral E1(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param e1: out: the function value
  !
  implicit none

  real(8) :: e1
  real(8) :: es1
  real(8) :: es2
  real(8) :: x

  if ( x == 0.0D+00 ) then

     e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

     e1 = - log ( x ) + (((( &
          1.07857D-03 * x &
          - 9.76004D-03 ) * x &
          + 5.519968D-02 ) * x &
          - 0.24991055D+00 ) * x &
          + 0.99999193D+00 ) * x &
          - 0.57721566D+00

  else

     es1 = ((( x &
          + 8.5733287401D+00 ) * x &
          +18.059016973D+00  ) * x &
          + 8.6347608925D+00 ) * x &
          + 0.2677737343D+00

     es2 = ((( x &
          +  9.5733223454D+00 ) * x &
          + 25.6329561486D+00 ) * x &
          + 21.0996530827D+00 ) * x &
          +  3.9584969228D+00

     e1 = exp ( - x ) / x * es1 / es2

  end if

  return
end subroutine e1xa




  !*****************************************************************************80
  !
  !! E1XB computes the exponential integral E1(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    06 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) E1, the function value.
  !
subroutine e1xb ( x, e1 )
  !E1XB computes the exponential integral E1(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param e1: out: the function value
  !
  implicit none

  real(8)    :: e1
  real(8)    :: ga
  integer(4) :: k
  integer(4) :: m
  real(8)    :: r
  real(8)    :: t
  real(8)    :: t0
  real(8)    :: x

  if ( x == 0.0D+00 ) then

     e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

     e1 = 1.0D+00
     r = 1.0D+00

     do k = 1, 25
        r = -r * k * x / ( k + 1.0D+00 )**2
        e1 = e1 + r
        if ( abs ( r ) <= abs ( e1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     ga = 0.5772156649015328D+00
     e1 = - ga - log ( x ) + x * e1

  else

     m = 20 + int ( 80.0D+00 / x )
     t0 = 0.0D+00
     do k = m, 1, -1
        t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
     end do
     t = 1.0D+00 / ( x + t0 )
     e1 = exp ( -x ) * t

  end if

  return
end subroutine e1xb




  !*****************************************************************************80
  !
  !! E1Z computes the complex exponential integral E1(z).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) CE1, the function value.
  !
subroutine e1z ( z, ce1 )
  !E1Z computes the complex exponential integral E1(z).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param z: in: the argument
  !@param ce1: out: the function value
  !
  implicit none

  real(8)    :: a0
  complex(8) :: ce1
  complex(8) :: cr
  complex(8) :: ct
  complex(8) :: ct0
  real(8)    :: el
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: x
  complex(8) :: z

  pi = 3.141592653589793D+00
  el = 0.5772156649015328D+00
  x = real ( z, kind = 8 )
  a0 = abs ( z )

  if ( a0 == 0.0D+00 ) then
     ce1 = cmplx ( 1.0D+300, 0.0D+00, kind = 8 )
  else if ( a0 <= 10.0D+00 .or. &
       ( x < 0.0D+00 .and. a0 < 20.0D+00 ) ) then
     ce1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 150
        cr = - cr * k * z / ( k + 1.0D+00 )**2
        ce1 = ce1 + cr
        if ( abs ( cr ) <= abs ( ce1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     ce1 = - el - log ( z ) + z * ce1

  else

     ct0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     do k = 120, 1, -1
        ct0 = k / ( 1.0D+00 + k / ( z + ct0 ) )
     end do
     ct = 1.0D+00 / ( z + ct0 )

     ce1 = exp ( - z ) * ct
     if ( x <= 0.0D+00 .and. imag ( z ) == 0.0D+00 ) then
        ce1 = ce1 - pi * cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
     end if

  end if

  return
end subroutine e1z




  !*****************************************************************************80
  !
  !! EIX computes the exponential integral Ei(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    10 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) EI, the function value.
  !
subroutine eix ( x, ei )
  !EIX computes the exponential integral Ei(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ei: out: the function value
  !
  implicit none

  real(8)    :: ei
  real(8)    :: ga
  integer(4) :: k
  real(8)    :: r
  real(8)    :: x

  if ( x == 0.0D+00 ) then

     ei = -1.0D+300

  else if ( x <= 40.0D+00 ) then

     ei = 1.0D+00
     r = 1.0D+00
     do k = 1, 100
        r = r * k * x / ( k + 1.0D+00 )**2
        ei = ei + r
        if ( abs ( r / ei ) <= 1.0D-15 ) then
           exit
        end if
     end do

     ga = 0.5772156649015328D+00
     ei = ga + log ( x ) + x * ei

  else

     ei = 1.0D+00
     r = 1.0D+00
     do k = 1, 20
        r = r * k / x
        ei = ei + r
     end do
     ei = exp ( x ) / x * ei

  end if

  return
end subroutine eix




  !*****************************************************************************80
  !
  !! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    12 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) HK, the modulus, between 0 and 1.
  !
  !    Input, real(8) PHI, the argument in degrees.
  !
  !    Output, real(8) FE, EE, the values of F(k,phi) and E(k,phi).
  !
subroutine elit ( hk, phi, fe, ee )
  !ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param hk: in: the modulus, between 0 and 1
  !@param phi: in: the argument in degrees
  !@param fe: out: F(k,phi)
  !@param ee: out: E(k,phi)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: b
  real(8)    :: b0
  real(8)    :: c
  real(8)    :: ce
  real(8)    :: ck
  real(8)    :: d
  real(8)    :: d0
  real(8)    :: ee
  real(8)    :: fac
  real(8)    :: fe
  real(8)    :: g
  real(8)    :: hk
  integer(4) :: n
  real(8)    :: phi
  real(8)    :: pi
  real(8)    :: r

  g = 0.0D+00
  pi = 3.14159265358979D+00
  a0 = 1.0D+00
  b0 = sqrt ( 1.0D+00 - hk * hk )
  d0 = ( pi / 180.0D+00 ) * phi
  r = hk * hk

  if ( hk == 1.0D+00 .and. phi == 90.0D+00 ) then

     fe = 1.0D+300
     ee = 1.0D+00

  else if ( hk == 1.0D+00 ) then

     fe = log ( ( 1.0D+00 + sin ( d0 ) ) / cos ( d0 ) )
     ee = sin ( d0 )

  else

     fac = 1.0D+00
     do n = 1, 40
        a = ( a0 + b0 ) /2.0D+00
        b = sqrt ( a0 * b0 )
        c = ( a0 - b0 ) / 2.0D+00
        fac = 2.0D+00 * fac
        r = r + fac * c * c
        if ( phi /= 90.0D+00 ) then
           d = d0 + atan ( ( b0 / a0 ) * tan ( d0 ) )
           g = g + c * sin( d )
           d0 = d + pi * int ( d / pi + 0.5D+00 )
        end if
        a0 = a
        b0 = b
        if ( c < 1.0D-07 ) then
           exit
        end if
     end do

     ck = pi / ( 2.0D+00 * a )
     ce = pi * ( 2.0D+00 - r ) / ( 4.0D+00 * a )
     if ( phi == 90.0D+00 ) then
        fe = ck
        ee = ce
     else
        fe = d / ( fac * a )
        ee = fe * ce / ck + g
     end if

  end if

  return
end subroutine elit




  !*****************************************************************************80
  !
  !! ELIT3 computes the elliptic integral of the third kind.
  !
  !  Discussion:
  !
  !    Gauss-Legendre quadrature is used.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    14 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) PHI, the argument in degrees.
  !
  !    Input, real(8) HK, the modulus, between 0 and 1.
  !
  !    Input, real(8) C, the parameter, between 0 and 1.
  !
  !    Output, real(8) EL3, the value of the elliptic integral
  !    of the third kind.
  !
subroutine elit3 ( phi, hk, c, el3 )
  !ELIT3 computes the elliptic integral of the third kind.
  !
  !Gauss-Legendre quadrature is used.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param phi: in: the argument in degrees
  !@param hk: in: the modulus, between 0 and 1
  !@param c: in: the parameter, between 0 and 1
  !@param el3: out: the value of the elliptic integral of the third kind
  !
  implicit none

  real(8)    :: c
  real(8)    :: c0
  real(8)    :: c1
  real(8)    :: c2
  real(8)    :: el3
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: hk
  integer(4) :: i
  logical    :: lb1
  logical    :: lb2
  real(8)    :: phi
  real(8), dimension ( 10 ), save :: t = (/ &
       0.9931285991850949D+00, 0.9639719272779138D+00, &
       0.9122344282513259D+00, 0.8391169718222188D+00, &
       0.7463319064601508D+00, 0.6360536807265150D+00, &
       0.5108670019508271D+00, 0.3737060887154195D+00, &
       0.2277858511416451D+00, 0.7652652113349734D-01 /)
  real(8)    :: t1
  real(8)    :: t2
  real(8), dimension ( 10 ), save :: w = (/ &
       0.1761400713915212D-01, 0.4060142980038694D-01, &
       0.6267204833410907D-01, 0.8327674157670475D-01, &
       0.1019301198172404D+00, 0.1181945319615184D+00, &
       0.1316886384491766D+00, 0.1420961093183820D+00, &
       0.1491729864726037D+00, 0.1527533871307258D+00 /)

  lb1 = ( hk == 1.0D+00 ) .and. ( abs ( phi - 90.0D+00 ) <= 1.0D-08 )

  lb2 = c == 1.0D+00 .and. abs ( phi - 90.0D+00 ) <= 1.0D-08

  if ( lb1 .or. lb2 ) then
     el3 = 1.0D+300
     return
  end if

  c1 = 0.87266462599716D-02 * phi
  c2 = c1

  el3 = 0.0D+00
  do i = 1, 10
     c0 = c2 * t(i)
     t1 = c1 + c0
     t2 = c1 - c0
     f1 = 1.0D+00 / ( ( 1.0D+00 - c * sin(t1) * sin(t1) ) &
          * sqrt ( 1.0D+00 - hk * hk * sin ( t1 ) * sin ( t1 ) ) )
     f2 = 1.0D+00 / ( ( 1.0D+00 - c * sin ( t2 ) * sin ( t2 ) ) &
          * sqrt( 1.0D+00 - hk * hk * sin ( t2 ) * sin ( t2 ) ) )
     el3 = el3 + w(i) * ( f1 + f2 )
  end do

  el3 = c1 * el3

  return
end subroutine elit3




  !*****************************************************************************80
  !
  !! ENVJ is a utility function used by MSTA1 and MSTA2.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    14 March 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, ?
  !
  !    Input, real(8) X, ?
  !
  !    Output, real(8) ENVJ, ?
  !
function envj ( n, x )
  !ENVJ is a utility function used by MSTA1 and MSTA2.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param envj: out: the value of the function
  !
  implicit none

  real(8)    :: envj
  integer(4) :: n
  real(8)    :: x

  envj = 0.5D+00 * log10 ( 6.28D+00 * n ) - n * log10 ( 1.36D+00 * x / n )

  return
end function envj




  !*****************************************************************************80
  !
  !! ENXA computes the exponential integral En(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) EN(0:N), the function values.
  !
subroutine enxa ( n, x, en )
  !ENXA computes the exponential integral En(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param en: out: the function values
  !
  implicit none

  integer(4) :: n

  real(8)    :: e1
  real(8)    :: ek
  real(8)    :: en(0:n)
  integer(4) :: k
  real(8)    :: x

  en(0) = exp ( - x ) / x 
  call e1xb ( x, e1 )

  en(1) = e1
  do k = 2, n
     ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
     en(k) = ek
     e1 = ek
  end do

  return
end subroutine enxa




  !*****************************************************************************80
  !
  !! ENXB computes the exponential integral En(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    10 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) EN(0:N), the function values.
  !
subroutine enxb ( n, x, en )
  !ENXB computes the exponential integral En(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param en: out: the function values
  !
  implicit none

  integer(4) :: n

  real(8)    :: en(0:n)
  real(8)    :: ens
  integer(4) :: j
  integer(4) :: k
  integer(4) :: l
  integer(4) :: m
  real(8)    :: ps
  real(8)    :: r
  real(8)    :: rp
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: t
  real(8)    :: t0
  real(8)    :: x

  if ( x == 0.0D+00 ) then

     en(0) = 1.0D+300
     en(1) = 1.0D+300
     do k = 2, n
        en(k) = 1.0D+00 / ( k - 1.0D+00 )
     end do
     return

  else if ( x <= 1.0D+00 ) then

     en(0) = exp ( - x ) / x
     do l = 1, n
        rp = 1.0D+00
        do j = 1, l - 1
           rp = - rp * x / j
        end do
        ps = -0.5772156649015328D+00
        do m = 1, l - 1
           ps = ps + 1.0D+00 / m
        end do
        ens = rp * ( - log ( x ) + ps )
        s = 0.0D+00
        do m = 0, 20
           if ( m /= l - 1 ) then
              r = 1.0D+00
              do j = 1, m
                 r = - r * x / j
              end do
              s = s + r / ( m - l + 1.0D+00 )
              if ( abs ( s - s0 ) < abs ( s ) * 1.0D-15 ) then
                 exit
              end if
              s0 = s
           end if
        end do

        en(l) = ens - s

     end do

  else

     en(0) = exp ( - x ) / x
     m = 15 + int ( 100.0D+00 / x )
     do l = 1, n
        t0 = 0.0D+00
        do k = m, 1, -1
           t0 = ( l + k - 1.0D+00 ) / ( 1.0D+00 + k / ( x + t0 ) )
        end do
        t = 1.0D+00 / ( x + t0 )
        en(l) = exp ( - x ) * t
     end do

  end if

  return
end subroutine enxb




  !*****************************************************************************80
  !
  !! WERROR evaluates the error function.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) ERR, the function value.
  !
subroutine werror ( x, err )
  !WERROR evaluates the error function.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param err: out: the function value
  !
  implicit none

  real(8)    :: c0
  real(8)    :: eps
  real(8)    :: er
  real(8)    :: err
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: x
  real(8)    :: x2

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  x2 = x * x

  if ( abs ( x ) < 3.5D+00 ) then

     er = 1.0D+00
     r = 1.0D+00

     do k = 1, 50
        r = r * x2 / ( k + 0.5D+00 )
        er = er + r
        if ( abs ( r ) <= abs ( er ) * eps ) then
           exit
        end if
     end do

     c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
     err = c0 * er

  else

     er = 1.0D+00
     r = 1.0D+00
     do k = 1, 12
        r = - r * ( k - 0.5D+00 ) / x2
        er = er + r
     end do

     c0 = exp ( - x2 ) / ( abs ( x ) * sqrt ( pi ) )

     err = 1.0D+00 - c0 * er
     if ( x < 0.0D+00 ) then
        err = -err
     end if

  end if

  return
end subroutine werror




  !*****************************************************************************80
  !
  !! EULERA computes the Euler number En.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    10 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the index of the highest value to compute.
  !
  !    Output, real(8) EN(0:N), the Euler numbers up to the N-th value.
  !
subroutine eulera ( n, en )
  !EULERA computes the Euler number En.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the index of the highest value to compute
  !@param en: out: the Euler numbers up to the N-th value
  !
  implicit none

  integer(4) :: n

  real(8)    :: en(0:n)
  integer(4) :: j
  integer(4) :: k
  integer(4) :: m
  real(8)    :: r
  real(8)    :: s

  en(0) = 1.0D+00

  do m = 1, n / 2
     s = 1.0D+00
     do k = 1, m - 1
        r = 1.0D+00
        do j = 1, 2 * k
           r = r * ( 2.0D+00 * m - 2.0D+00 * k + j ) / j
        end do
        s = s + r * en(2*k)
     end do
     en(2*m) = -s
  end do

  return
end subroutine eulera




  !*****************************************************************************80
  !
  !! EULERB computes the Euler number En.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    09 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the index of the highest value to compute.
  !
  !    Output, real(8) EN(0:N), the Euler numbers up to the N-th value.
  !
subroutine eulerb ( n, en )
  !EULERB computes the Euler number En.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the index of the highest value to compute
  !@param en: out: the Euler numbers up to the N-th value
  !
  implicit none

  integer(4) :: n

  real(8)    :: en(0:n)
  real(8)    :: hpi
  real(8)    :: isgn
  integer(4) :: k
  integer(4) :: m
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: s

  hpi = 2.0D+00 / 3.141592653589793D+00
  en(0) = 1.0D+00
  en(2) = -1.0D+00
  r1 = -4.0D+00 * hpi ** 3

  do m = 4, n, 2
     r1 = - r1 * ( m - 1 ) * m * hpi * hpi
     r2 = 1.0D+00
     isgn = 1.0D+00
     do k = 3, 1000, 2
        isgn = - isgn
        s = ( 1.0D+00 / k ) ** ( m + 1 )
        r2 = r2 + isgn * s
        if ( s < 1.0D-15 ) then
           exit
        end if
     end do

     en(m) = r1 * r2

  end do

  return
end subroutine eulerb




  !*****************************************************************************80
  !
  !! FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code.
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the Mathieu function.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Input, real(8) A, the characteristic value of the Mathieu
  !    functions for given m and q.
  !
  !    Output, real(8) FC(*), the expansion coefficients of Mathieu
  !    functions ( k =  1,2,...,KM ).  FC(1),FC(2),FC(3),... correspond to
  !    A0,A2,A4,... for KD = 1 case, 
  !    A1,A3,A5,... for KD = 2 case,
  !    B1,B3,B5,... for KD = 3 case,
  !    B2,B4,B6,... for KD = 4 case.
  !
subroutine fcoef ( kd, m, q, a, fc )
  !FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code. 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !* FC: the expansion coefficients of Mathieu functions ( k = 1,2,...,KM ). FC(1), FC(2), FC(3), ... correspond to A0,
  !  A2, A4, ... for KD = 1 case, A1, A3, A5, ... for KD = 2 case, B1, B3, B5, ... for KD = 3 case, B2, B4, B6, ... for
  !  KD = 4 case
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the Mathieu function
  !@param q: in: the parameter of the Mathieu functions
  !@param a: in: the characteristic value of the Mathieu functions for given m and q
  !@param fc: out: expansion coefficients of the Mathieu functions, k=1,...,KM
  !
  implicit none

  real(8)    :: a
  real(8)    :: f
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: f3
  real(8)    :: fc(251)
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: kb
  integer(4) :: kd
  integer(4) :: km
  integer(4) :: l
  integer(4) :: m
  real(8)    :: q
  real(8)    :: qm
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: sp
  real(8)    :: ss
  real(8)    :: u
  real(8)    :: v

  if ( q <= 1.0D+00 ) then
     qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q &
          + 90.7D+00 * sqrt ( q ) * q
  else
     qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q &
          + 0.0037D+00 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5D+00 * m )

  if ( q == 0.0D+00 ) then

     do k = 1, km
        fc(k) = 0.0D+00
     end do

     if ( kd == 1 ) then
        fc((m+2)/2) = 1.0D+00
        if (m == 0 ) then
           fc(1) = 1.0D+00 / sqrt ( 2.0D+00 )
        end if
     else if ( kd == 4 ) then
        fc(m/2) = 1.0D+00
     else
        fc((m+1)/2) = 1.0D+00
     end if

     return

  end if

  kb = 0
  s = 0.0D+00
  f = 1.0D-100
  u = 0.0D+00
  fc(km) = 0.0D+00

  if ( kd == 1 ) then

     l = 0

     do k = km, 3, -1

        v = u
        u = f
        f = ( a - 4.0D+00 * k * k ) * u / q - v

        if ( abs ( f ) < abs ( fc(k+1) ) ) then

           kb = k
           fc(1) = 1.0D-100
           sp = 0.0D+00
           f3 = fc(k+1)
           fc(2) = a / q * fc(1)
           fc(3) = ( a - 4.0D+00 ) * fc(2) / q - 2.0D+00 * fc(1)
           u = fc(2)
           f1 = fc(3)

           do i = 3, kb
              v = u
              u = f1
              f1 = ( a - 4.0D+00 * ( i - 1.0D+00 ) ** 2 ) * u / q - v
              fc(i+1) = f1
              if ( i == kb ) then
                 f2 = f1
              else
                 sp = sp + f1 * f1
              end if
           end do

           sp = sp + 2.0D+00 * fc(1) ** 2 + fc(2) ** 2 + fc(3) ** 2
           ss = s + sp * ( f3 / f2 ) ** 2
           s0 = sqrt ( 1.0D+00 / ss )
           do j = 1, km
              if ( j <= kb + 1 ) then
                 fc(j) = s0 * fc(j) * f3 / f2
              else
                 fc(j) = s0 * fc(j)
              end if
           end do
           l = 1
           exit
        else
           fc(k) = f
           s = s + f * f
        end if

     end do

     if ( l == 0 ) then
        fc(2) = q * fc(3) / ( a - 4.0D+00 - 2.0D+00 * q * q / a )
        fc(1) = q / a * fc(2)
        s = s + 2.0D+00 * fc(1) ** 2 + fc(2) ** 2
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
           fc(k) = s0 * fc(k)
        end do
     end if

  else if ( kd == 2 .or. kd == 3 ) then

     l = 0

     do k = km, 3, -1

        v = u
        u = f
        f = ( a - ( 2.0D+00 * k - 1 ) ** 2 ) * u / q - v

        if ( abs ( fc(k) ) <= abs ( f ) ) then
           fc(k-1) = f
           s = s + f * f
        else
           kb = k
           f3 = fc(k)
           l = 1
           exit
        end if

     end do

     if ( l == 0 ) then

        fc(1) = q / ( a - 1.0D+00 - ( - 1 ) ** kd * q ) * fc(2)
        s = s + fc(1) * fc(1)
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
           fc(k) = s0 * fc(k)
        end do

     else

        fc(1) = 1.0D-100
        fc(2) = ( a - 1.0D+00 - ( - 1 ) ** kd * q ) / q * fc(1)
        sp = 0.0D+00
        u = fc(1)
        f1 = fc(2)
        do i = 2, kb - 1
           v = u
           u = f1
           f1 = ( a - ( 2.0D+00 * i - 1.0D+00 ) ** 2 ) * u / q - v
           if ( i /= kb - 1 ) then
              fc(i+1) = f1
              sp = sp + f1 * f1
           else
              f2 = f1
           end if
        end do

        sp = sp + fc(1) ** 2 + fc(2) ** 2
        ss = s + sp * ( f3 / f2 ) ** 2
        s0 = 1.0D+00 / sqrt ( ss )
        do j = 1, km
           if ( j < kb ) then
              fc(j) = s0 * fc(j) * f3 / f2
           else
              fc(j) = s0 * fc(j)
           end if
        end do

     end if

  else if ( kd == 4 ) then

     l = 0

     do k = km, 3, -1
        v = u
        u = f
        f = ( a - 4.0D+00 * k * k ) * u / q - v
        if ( abs ( fc(k) ) <= abs ( f ) ) then
           fc(k-1) = f
           s = s + f * f
        else
           kb = k
           f3 = fc(k)
           l = 1
           exit
        end if
     end do

     if ( l == 0 ) then

        fc(1) = q / ( a - 4.0D+00 ) * fc(2)
        s = s + fc(1) * fc(1)
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
           fc(k) = s0 * fc(k)
        end do

     else

        fc(1) = 1.0D-100
        fc(2) = ( a - 4.0D+00 ) / q * fc(1)
        sp = 0.0D+00
        u = fc(1)
        f1 = fc(2)

        do i = 2, kb - 1
           v = u
           u = f1
           f1 = ( a - 4.0D+00 * i * i ) * u / q - v
           if ( i /= kb - 1 ) then
              fc(i+1) = f1
              sp = sp + f1 * f1
           else
              f2 = f1
           end if
        end do

        sp = sp + fc(1) ** 2 + fc(2) ** 2
        ss = s + sp * ( f3 / f2 ) ** 2
        s0 = 1.0D+00 / sqrt ( ss )

        do j = 1, km
           if ( j < kb ) then
              fc(j) = s0 * fc(j) * f3 / f2
           else
              fc(j) = s0 * fc(j)
           end if
        end do

     end if

  end if

  if ( fc(1) < 0.0D+00 ) then
     do j = 1, km
        fc(j) = -fc(j)
     end do
  end if

  return
end subroutine fcoef




  !*****************************************************************************80
  !
  !! FCS computes Fresnel integrals C(x) and S(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    17 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) C, S, the function values.
  !
subroutine fcs ( x, c, s )
  !FCS computes Fresnel integrals C(x) and S(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param c: out: the Fresnel integral C(x)
  !@param s: out: the Fresnel integral S(x)
  !
  implicit none

  real(8)    :: c
  real(8)    :: eps
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: g
  integer(4) :: k
  integer(4) :: m
  real(8)    :: pi
  real(8)    :: px
  real(8)    :: q
  real(8)    :: r
  real(8)    :: s
  real(8)    :: su
  real(8)    :: t
  real(8)    :: t0
  real(8)    :: t2
  real(8)    :: x
  real(8)    :: xa

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  xa = abs ( x )
  px = pi * xa
  t = 0.5D+00 * px * xa
  t2 = t * t

  if ( xa == 0.0D+00 ) then

     c = 0.0D+00
     s = 0.0D+00

  else if ( xa < 2.5D+00 ) then

     r = xa
     c = r
     do k = 1, 50
        r = -0.5D+00 * r * ( 4.0D+00 * k - 3.0D+00 ) / k &
             / ( 2.0D+00 * k - 1.0D+00 ) / ( 4.0D+00 * k + 1.0D+00 ) * t2
        c = c + r
        if ( abs ( r ) < abs ( c ) * eps ) then
           exit
        end if
     end do

     s = xa * t / 3.0D+00
     r = s
     do k = 1, 50
        r = - 0.5D+00 * r * ( 4.0D+00 * k - 1.0D+00 ) / k &
             / ( 2.0D+00 * k + 1.0D+00 ) / ( 4.0D+00 * k + 3.0D+00 ) * t2
        s = s + r
        if ( abs ( r ) < abs ( s ) * eps ) then
           if ( x < 0.0D+00 ) then
              c = -c
              s = -s
           end if
           return
        end if
     end do

  else if ( xa < 4.5D+00 ) then

     m = int ( 42.0D+00 + 1.75D+00 * t )
     su = 0.0D+00
     c = 0.0D+00
     s = 0.0D+00
     f1 = 0.0D+00
     f0 = 1.0D-100

     do k = m, 0, -1
        f = ( 2.0D+00 * k + 3.0D+00 ) * f0 / t - f1
        if ( k == int ( k / 2 ) * 2 ) then
           c = c + f
        else
           s = s + f
        end if
        su = su + ( 2.0D+00 * k + 1.0D+00 ) * f * f
        f1 = f0
        f0 = f
     end do

     q = sqrt ( su )
     c = c * xa / q
     s = s * xa / q

  else

     r = 1.0D+00
     f = 1.0D+00
     do k = 1, 20
        r = -0.25D+00 * r * ( 4.0D+00 * k - 1.0D+00 ) &
             * ( 4.0D+00 * k - 3.0D+00 ) / t2
        f = f + r
     end do
     r = 1.0D+00 / ( px * xa )
     g = r
     do k = 1, 12
        r = -0.25D+00 * r * ( 4.0D+00 * k + 1.0D+00 ) &
             * ( 4.0D+00 * k - 1.0D+00 ) / t2
        g = g + r
     end do

     t0 = t - int ( t / ( 2.0D+00 * pi ) ) * 2.0D+00 * pi
     c = 0.5D+00 + ( f * sin ( t0 ) - g * cos ( t0 ) ) / px
     s = 0.5D+00 - ( f * cos ( t0 ) + g * sin ( t0 ) ) / px

  end if

  if ( x < 0.0D+00 ) then
     c = -c
     s = -s
  end if

  return
end subroutine fcs




  !*****************************************************************************80
  !
  !! FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    17 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KF, the function code.
  !    1 for C(z);
  !    2 for S(z)
  !
  !    Input, integer(4) NT, the total number of zeros desired.
  !
  !    Output, complex(8) Z0(NT), the zeros.
  !
subroutine fcszo ( kf, nt, zo )
  !FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param kf: in: the function code. 1 for C(z); 2 for S(z)
  !@param nt: in: the total number of zeros desired
  !@param zo: the zeros
  !
  implicit none

  integer(4) :: nt

  integer(4) :: i
  integer(4) :: it
  integer(4) :: j
  integer(4) :: kf
  integer(4) :: nr
  real(8)    :: pi
  real(8)    :: psq
  real(8)    :: px
  real(8)    :: py
  real(8)    :: w
  real(8)    :: w0
  complex(8) :: z
  complex(8) :: zd
  complex(8) :: zf
  complex(8) :: zfd
  complex(8) :: zgd
  complex(8) :: zo(nt)
  complex(8) :: zp
  complex(8) :: zq
  complex(8) :: zw

  pi = 3.141592653589793D+00

  do nr = 1, nt

     if ( kf == 1 ) then
        psq = sqrt ( 4.0D+00 * nr - 1.0D+00 )
     else
        psq = 2.0D+00 * sqrt ( real ( nr, kind = 8 ) )
     end if

     px = psq - log ( pi * psq ) / ( pi * pi * psq ** 3.0D+00 )
     py = log ( pi * psq ) / ( pi * psq )
     z = cmplx ( px, py )

     if ( kf == 2 ) then
        if ( nr == 2 ) then
           z = cmplx ( 2.8334D+00, 0.2443D+00 )
        else if ( nr == 3 ) then
           z = cmplx ( 3.4674D+00, 0.2185D+00 )
        else if ( nr == 4 ) then
           z = cmplx ( 4.0025D+00, 0.2008D+00 )
        end if
     end if

     it = 0

     do

        it = it + 1

        if ( kf == 1 ) then
           call cfc ( z, zf, zd )
        else
           call cfs ( z, zf, zd )
        end if

        zp = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do i = 1, nr - 1
           zp = zp * ( z - zo(i) )
        end do
        zfd = zf / zp
        zq = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

        do i = 1, nr - 1
           zw = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           do j = 1, nr - 1
              if ( j /= i ) then
                 zw = zw * ( z - zo(j) )
              end if
           end do
           zq = zq + zw
        end do

        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = cdabs ( z )

        if ( abs ( ( w - w0 ) / w ) <= 1.0D-12 ) then
           exit
        end if

        if ( 50 < it ) then
           exit
        end if

     end do

     zo(nr) = z

  end do

  return
end subroutine fcszo




  !*****************************************************************************80
  !
  !! FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    23 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KS, the sign code.
  !    0, to calculate F+(x) and K+(x);
  !    1, to calculate F_(x) and K_(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) FR, FI, FM, FA, the values of
  !    Re[F+/-(x)], Im[F+/-(x)], |F+/-(x)|, Arg[F+/-(x)]  (Degs.).
  !
  !    Output, real(8) GR, GI, GM, GA, the values of
  !    Re[K+/-(x)], Im[K+/-(x)], |K+/-(x)|, Arg[K+/-(x)]  (Degs.).
  !       
subroutine ffk ( ks, x, fr, fi, fm, fa, gr, gi, gm, ga )
  !FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param ks: in: the sign code. 0, to calculate F+(x) and K+(x); 1, to calculate F_(x) and K_(x)
  !@param x: in: the argument
  !@param fr: out: :math:`Re[F\pm(x)]`
  !@param fi: out: :math:`Im[F\pm(x)]`
  !@param fm: out: :math:`|F\pm(x)|`
  !@param fa: out: :math:`Arg[F\pm(x)]` (Degs.)
  !@param gr: out: :math:`Re[K\pm(x)]`
  !@param gi: out: :math:`Im[K\pm(x)]`
  !@param gm: out: :math:`|K\pm(x)|`
  !@param ga: out: :math:`Arg[K\pm(x)]` (Degs.)
  !
  implicit none

  real(8)    :: c1
  real(8)    :: cs
  real(8)    :: eps
  real(8)    :: fa
  real(8)    :: fi
  real(8)    :: fi0
  real(8)    :: fm
  real(8)    :: fr
  real(8)    :: ga
  real(8)    :: gi
  real(8)    :: gm
  real(8)    :: gr
  integer(4) :: k
  integer(4) :: ks
  integer(4) :: m
  real(8)    :: p2p
  real(8)    :: pi
  real(8)    :: pp2
  real(8)    :: s1
  real(8)    :: srd
  real(8)    :: ss
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: x4
  real(8)    :: xa
  real(8)    :: xc
  real(8)    :: xf
  real(8)    :: xf0
  real(8)    :: xf1
  real(8)    :: xg
  real(8)    :: xp
  real(8)    :: xq
  real(8)    :: xq2
  real(8)    :: xr
  real(8)    :: xs
  real(8)    :: xsu
  real(8)    :: xw

  srd = 57.29577951308233D+00
  eps = 1.0D-15
  pi = 3.141592653589793D+00
  pp2 = 1.2533141373155D+00
  p2p = 0.7978845608028654D+00
  xa = abs ( x )
  x2 = x * x
  x4 = x2 * x2

  if ( x == 0.0D+00 ) then

     fr = 0.5D+00 * sqrt ( 0.5D+00 * pi )
     fi = ( -1.0D+00 ) ** ks * fr
     fm = sqrt ( 0.25D+00 * pi )
     fa = ( -1.0D+00 ) ** ks * 45.0D+00
     gr = 0.5D+00
     gi = 0.0D+00
     gm = 0.5D+00
     ga = 0.0D+00

  else

     if ( xa <= 2.5D+00 ) then

        xr = p2p * xa
        c1 = xr
        do k = 1, 50
           xr = -0.5D+00 * xr * ( 4.0D+00 * k - 3.0D+00 ) / k &
                / ( 2.0D+00 * k - 1.0D+00 ) &
                / ( 4.0D+00 * k + 1.0D+00 ) * x4
           c1 = c1 + xr
           if ( abs ( xr / c1 ) < eps ) then
              exit
           end if
        end do

        s1 = p2p * xa * xa * xa / 3.0D+00
        xr = s1
        do k = 1, 50
           xr = -0.5D+00 * xr * ( 4.0D+00 * k - 1.0D+00 ) &
                / k / ( 2.0D+00 * k + 1.0D+00 ) &
                / ( 4.0D+00 * k + 3.0D+00 ) * x4
           s1 = s1 + xr
           if ( abs ( xr / s1 ) < eps ) then
              exit
           end if
        end do

     else if ( xa < 5.5D+00 ) then

        m = int ( 42.0D+00 + 1.75D+00 * x2 )
        xsu = 0.0D+00
        xc = 0.0D+00
        xs = 0.0D+00
        xf1 = 0.0D+00
        xf0 = 1.0D-100
        do k = m, 0, -1
           xf = ( 2.0D+00 * k + 3.0D+00 ) * xf0 / x2 - xf1
           if ( k == 2 * int ( k / 2 ) )  then
              xc = xc + xf
           else
              xs = xs + xf
           end if
           xsu = xsu + ( 2.0D+00 * k + 1.0D+00 ) * xf * xf
           xf1 = xf0
           xf0 = xf
        end do
        xq = sqrt ( xsu )
        xw = p2p * xa / xq
        c1 = xc * xw
        s1 = xs * xw

     else

        xr = 1.0D+00
        xf = 1.0D+00
        do k = 1, 12
           xr = -0.25D+00 * xr * ( 4.0D+00 * k - 1.0D+00 ) &
                * ( 4.0D+00 * k - 3.0D+00 ) / x4
           xf = xf + xr
        end do
        xr = 1.0D+00 / ( 2.0D+00 * xa * xa )
        xg = xr
        do k = 1, 12
           xr = -0.25D+00 * xr * ( 4.0D+00 * k + 1.0D+00 ) &
                * ( 4.0D+00 * k - 1.0D+00 ) / x4
           xg = xg + xr
        end do
        c1 = 0.5D+00 + ( xf * sin ( x2 ) - xg * cos ( x2 ) ) &
             / sqrt ( 2.0D+00 * pi ) / xa
        s1 = 0.5D+00 - ( xf * cos ( x2 ) + xg * sin ( x2 ) ) &
             / sqrt ( 2.0D+00 * pi ) / xa

     end if

     fr = pp2 * ( 0.5D+00 - c1 )
     fi0 = pp2 * ( 0.5D+00 - s1 )
     fi = ( -1.0D+00 ) ** ks * fi0
     fm = sqrt ( fr * fr + fi * fi )

     if ( 0.0D+00 <= fr ) then
        fa = srd * atan ( fi / fr )
     else if ( 0.0D+00 < fi ) then
        fa = srd * ( atan ( fi / fr ) + pi )
     else if ( fi < 0.0D+00 ) then
        fa = srd * ( atan ( fi / fr ) - pi )
     end if

     xp = x * x + pi / 4.0D+00
     cs = cos ( xp )
     ss = sin ( xp )
     xq2 = 1.0D+00 / sqrt ( pi )
     gr = xq2 * ( fr * cs + fi0 * ss )
     gi = ( -1.0D+00 ) ** ks * xq2 * ( fi0 * cs - fr * ss )
     gm = sqrt ( gr * gr + gi * gi )

     if ( 0.0D+00 <= gr ) then
        ga = srd * atan ( gi / gr )
     else if ( 0.0D+00 < gi ) then
        ga = srd * ( atan ( gi / gr ) + pi )
     else if ( gi < 0.0D+00 ) then
        ga = srd * ( atan ( gi / gr ) - pi )
     end if

     if ( x < 0.0D+00 ) then
        fr = pp2 - fr
        fi = ( -1.0D+00 ) ** ks * pp2 - fi
        fm = sqrt ( fr * fr + fi * fi )
        fa = srd * atan ( fi / fr )
        gr = cos ( x * x ) - gr
        gi = - ( -1.0D+00 ) ** ks * sin ( x * x ) - gi
        gm = sqrt ( gr * gr + gi * gi )
        ga = srd * atan ( gi / gr )
     end if

  end if

  return
end subroutine ffk




  !*****************************************************************************80
  !
  !! GAIH computes the GammaH function.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    09 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) GA, the function value.
  !
subroutine gaih ( x, ga )
  !GAIH computes the GammaH function.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ga: out: the function value
  !
  implicit none

  real(8)    :: ga
  integer(4) :: k
  integer(4) :: m
  integer(4) :: m1
  real(8)    :: pi
  real(8)    :: x

  pi = 3.141592653589793D+00

  if ( x == int ( x ) .and. 0.0 < x ) then
     ga = 1.0D+00
     m1 = int ( x - 1.0D+00 )
     do k = 2, m1
        ga = ga * k
     end do
  else if ( x + 0.5D+00 == int ( x + 0.5D+00) .and. 0.0D+00 < x ) then
     m = int ( x )
     ga = sqrt ( pi )
     do k = 1, m
        ga = 0.5D+00 * ga * ( 2.0D+00 * k - 1.0D+00 )
     end do
  end if

  return
end subroutine gaih




  !*****************************************************************************80
  !
  !! GAM0 computes the Gamma function for the LAMV function.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    09 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) GA, the function value.
  !   
subroutine gam0 ( x, ga )
  !GAM0 computes the Gamma function for the LAMV function.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ga: out: the function value
  !
  implicit none

  real(8), dimension ( 25 ) :: g = (/ &
       1.0D+00, &
       0.5772156649015329D+00, &
       -0.6558780715202538D+00, &
       -0.420026350340952D-01, &
       0.1665386113822915D+00, &
       -0.421977345555443D-01, &
       -0.96219715278770D-02, &
       0.72189432466630D-02, &
       -0.11651675918591D-02, &
       -0.2152416741149D-03, &
       0.1280502823882D-03, &
       -0.201348547807D-04, &
       -0.12504934821D-05, &
       0.11330272320D-05, &
       -0.2056338417D-06, &
       0.61160950D-08, &
       0.50020075D-08, &
       -0.11812746D-08, &
       0.1043427D-09, &
       0.77823D-11, &
       -0.36968D-11, &
       0.51D-12, &
       -0.206D-13, &
       -0.54D-14, &
       0.14D-14 /)
  real(8)    :: ga
  real(8)    :: gr
  integer(4) :: k
  real(8)    :: x

  gr = g(25)
  do k = 24, 1, -1
     gr = gr * x + g(k)
  end do

  ga = 1.0D+00 / ( gr * x )

  return
end subroutine gam0




  !*****************************************************************************80
  !
  !! GAMMA evaluates the Gamma function.
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by 
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
  !    incorporate this routine into a user program that the copyright 
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !    X must not be 0, or any negative integer.
  !
  !    Output, real(8) GA, the value of the Gamma function.
  !
subroutine gammaf ( x, ga )
  !GAMMA evaluates the Gamma function.
  !
  !Original FORTRAN77 version by Shanjie Zhang, Jianming Jin. FORTRAN90 version by John Burkardt.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument, which must not be 0, or any negative integer
  !@param ga: out: the value of the Gamma function
  !
  implicit none

  real(8), dimension ( 26 ) :: g = (/ &
       1.0D+00, &
       0.5772156649015329D+00, &
       -0.6558780715202538D+00, &
       -0.420026350340952D-01, &
       0.1665386113822915D+00, &
       -0.421977345555443D-01, &
       -0.96219715278770D-02, &
       0.72189432466630D-02, &
       -0.11651675918591D-02, &
       -0.2152416741149D-03, &
       0.1280502823882D-03, & 
       -0.201348547807D-04, &
       -0.12504934821D-05, &
       0.11330272320D-05, &
       -0.2056338417D-06, & 
       0.61160950D-08, &
       0.50020075D-08, &
       -0.11812746D-08, &
       0.1043427D-09, & 
       0.77823D-11, &
       -0.36968D-11, &
       0.51D-12, &
       -0.206D-13, &
       -0.54D-14, &
       0.14D-14, &
       0.1D-15 /)
  real(8)             :: ga
  real(8)             :: gr
  integer(4)          :: k
  integer(4)          :: m
  integer(4)          :: m1
  real(8), parameter  :: pi = 3.141592653589793D+00
  real(8)             :: r
  real(8)             :: x
  real(8)             :: z

  if ( x == aint ( x ) ) then

     if ( 0.0D+00 < x ) then
        ga = 1.0D+00
        m1 = int ( x ) - 1
        do k = 2, m1
           ga = ga * k
        end do
     else
        ga = 1.0D+300
     end if

  else

     if ( 1.0D+00 < abs ( x ) ) then
        z = abs ( x )
        m = int ( z )
        r = 1.0D+00
        do k = 1, m
           r = r * ( z - real ( k, kind = 8 ) )
        end do
        z = z - real ( m, kind = 8 )
     else
        z = x
     end if

     gr = g(26)
     do k = 25, 1, -1
        gr = gr * z + g(k)
     end do

     ga = 1.0D+00 / ( gr * z )

     if ( 1.0D+00 < abs ( x ) ) then
        ga = ga * r
        if ( x < 0.0D+00 ) then
           ga = - pi / ( x* ga * sin ( pi * x ) )
        end if
     end if

  end if

  return
end subroutine gammaf




  !*****************************************************************************80
  !
  !! GMN computes quantities for oblate radial functions with small argument.
  !
  !  Discussion:
  !
  !    This procedure computes Gmn(-ic,ix) and its derivative for oblate
  !    radial functions with a small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) BK(*), coefficients.
  !
  !    Output, real(8) GF, GD, the value of Gmn(-C,X) and Gmn'(-C,X).
  !
subroutine gmn ( m, n, c, x, bk, gf, gd )
  !GMN computes quantities for oblate radial functions with small argument.
  !
  !This procedure computes Gmn(-ic,ix) and its derivative for oblate
  !radial functions with a small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param bk: in: coefficients
  !@param gf: out: :math:`Gmn(-C,X)`
  !@param gd: out: :math:`Gmn'(-C,X)`
  !
  implicit none

  real(8)    :: bk(200)
  real(8)    :: c
  real(8)    :: eps
  real(8)    :: gd
  real(8)    :: gd0
  real(8)    :: gd1
  real(8)    :: gf
  real(8)    :: gf0
  real(8)    :: gw
  integer(4) :: ip
  integer(4) :: k
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  real(8)    :: x
  real(8)    :: xm

  eps = 1.0D-14

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
  xm = ( 1.0D+00 + x * x ) ** ( -0.5D+00 * m )
  gf0 = 0.0D+00
  do k = 1, nm
     gf0 = gf0 + bk(k) * x ** ( 2.0D+00 * k - 2.0D+00 )
     if ( abs ( ( gf0 - gw ) / gf0 ) < eps .and. 10 <= k ) then
        exit
     end if
     gw = gf0
  end do

  gf = xm * gf0 * x ** ( 1 - ip )

  gd1 = - m * x / ( 1.0D+00 + x * x ) * gf
  gd0 = 0.0D+00

  do k = 1, nm

     if ( ip == 0 ) then
        gd0 = gd0 + ( 2.0D+00 * k - 1.0D+00 ) * bk(k) &
             * x ** ( 2.0D+00 * k - 2.0D+00 )
     else
        gd0 = gd0 + 2.0D+00 * k * bk(k+1) * x ** ( 2.0D+00 * k - 1.0D+00 )
     end if

     if ( abs ( ( gd0 - gw ) / gd0 ) < eps .and. 10 <= k ) then
        exit
     end if

     gw = gd0

  end do

  gd = gd1 + xm * gd0

  return
end subroutine gmn




  !*****************************************************************************80
  !
  !! HERZO computes the zeros the Hermite polynomial Hn(x).
  !
  !  Discussion:
  !
  !    This procedure computes the zeros of Hermite polynomial Ln(x)
  !    in the interval [-1,+1], and the corresponding
  !    weighting coefficients for Gauss-Hermite integration.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of the polynomial.
  !
  !    Output, real(8) X(N), the zeros.
  !
  !    Output, real(8) W(N), the corresponding weights.
  !
subroutine herzo ( n, x, w )
  !HERZO computes the zeros the Hermite polynomial Hn(x).
  !
  !This procedure computes the zeros of Hermite polynomial Ln(x)
  !in the interval [-1,+1], and the corresponding
  !weighting coefficients for Gauss-Hermite integration.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of the polynomial
  !@param x: out: the zeros
  !@param w: out: the corresponding weights
  !
  implicit none

  integer(4) :: n

  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: fd
  real(8)    :: gd
  real(8)    :: hd
  real(8)    :: hf
  real(8)    :: hn
  integer(4) :: i
  integer(4) :: it
  integer(4) :: j
  integer(4) :: k
  integer(4) :: nr
  real(8)    :: p
  real(8)    :: q
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: w(n)
  real(8)    :: wp
  real(8)    :: x(n)
  real(8)    :: x0
  real(8)    :: z
  real(8)    :: z0
  real(8)    :: zl

  hn = 1.0D+00 / n
  zl = -1.1611D+00 + 1.46D+00 * sqrt ( real ( n, kind = 8 ) )

  do nr = 1, n / 2

     if ( nr == 1 ) then
        z = zl
     else
        z = z - hn * ( n / 2 + 1 - nr )
     end if

     it = 0

     do

        it = it + 1
        z0 = z
        f0 = 1.0D+00
        f1 = 2.0D+00 * z
        do k = 2, n
           hf = 2.0D+00 * z * f1 - 2.0D+00 * ( k - 1.0D+00 ) * f0
           hd = 2.0D+00 * k * f1
           f0 = f1
           f1 = hf
        end do

        p = 1.0D+00
        do i = 1, nr - 1
           p = p * ( z - x(i) )
        end do
        fd = hf / p

        q = 0.0D+00
        do i = 1, nr - 1
           wp = 1.0D+00
           do j = 1, nr - 1
              if ( j /= i ) then
                 wp = wp * ( z - x(j) )
              end if
           end do
           q = q + wp
        end do

        gd = ( hd - q * fd ) / p
        z = z - fd / gd

        if ( 40 < it .or. abs ( ( z - z0 ) / z ) <= 1.0D-15 ) then
           exit
        end if

     end do

     x(nr) = z
     x(n+1-nr) = -z
     r = 1.0D+00
     do k = 1, n
        r = 2.0D+00 * r * k
     end do
     w(nr) = 3.544907701811D+00 * r / ( hd * hd )
     w(n+1-nr) = w(nr)

  end do

  if ( n /= 2 * int ( n / 2 ) ) then
     r1 = 1.0D+00
     r2 = 1.0D+00
     do j = 1, n
        r1 = 2.0D+00 * r1 * j
        if ( ( n + 1 ) / 2 <= j ) then
           r2 = r2 * j
        end if
     end do
     w(n/2+1) = 0.88622692545276D+00 * r1 / ( r2 * r2 )
     x(n/2+1) = 0.0D+00
  end if

  return
end subroutine herzo




  !*****************************************************************************80
  !
  !! HYGFX evaluates the hypergeometric function F(A,B,C,X).
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by 
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
  !    incorporate this routine into a user program that the copyright 
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, C, X, the arguments of the function.
  !    C must not be equal to a nonpositive integer.
  !    X < 1.
  !
  !    Output, real HF, the value of the function.
  !
subroutine hygfx ( a, b, c, x, hf )
  !HYGFX evaluates the hypergeometric function F(A,B,C,X).
  !
  !Original FORTRAN77 version by Shanjie Zhang, Jianming Jin. FORTRAN90 version by John Burkardt.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: the argument A of the function
  !@param b: in: the argument B of the function
  !@param c: in: the argument C of the function, not a nonpositive integer
  !@param x: in: the argument :math:'X' of the function, :math:'X'  < 1
  !@param hf: out: the value of the function
  !
  implicit none

  real(8)             :: a
  real(8)             :: a0
  real(8)             :: aa
  real(8)             :: b
  real(8)             :: bb
  real(8)             :: c
  real(8)             :: c0
  real(8)             :: c1
  real(8), parameter  :: el = 0.5772156649015329D+00
  real(8)             :: eps
  real(8)             :: f0
  real(8)             :: f1
  real(8)             :: g0
  real(8)             :: g1
  real(8)             :: g2
  real(8)             :: g3
  real(8)             :: ga
  real(8)             :: gabc
  real(8)             :: gam
  real(8)             :: gb
  real(8)             :: gbm
  real(8)             :: gc
  real(8)             :: gca
  real(8)             :: gcab
  real(8)             :: gcb
  real(8)             :: gm
  real(8)             :: hf
  real(8)             :: hw
  integer(4)          :: j
  integer(4)          :: k
  logical             :: l0
  logical             :: l1
  logical             :: l2
  logical             :: l3
  logical             :: l4
  logical             :: l5
  integer(4)          :: m
  integer(4)          :: nm
  real(8)             :: pa
  real(8)             :: pb
  real(8), parameter  :: pi = 3.141592653589793D+00
  real(8)             :: r
  real(8)             :: r0
  real(8)             :: r1
  real(8)             :: rm
  real(8)             :: rp
  real(8)             :: sm
  real(8)             :: sp
  real(8)             :: sp0
  real(8)             :: x
  real(8)             :: x1

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 .or. l1 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFX - Fatal error!'
     write ( *, '(a)' ) '  The hypergeometric series is divergent.'
     return
  end if

  if ( 0.95D+00 < x ) then 
     eps = 1.0D-08
  else
     eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

     hf = 1.0D+00
     return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

     call gammaf ( c, gc )
     call gammaf ( c - a - b, gcab )
     call gammaf ( c - a, gca )
     call gammaf ( c - b, gcb )
     hf = gc * gcab /( gca *gcb )
     return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

     g0 = sqrt ( pi ) * 2.0D+00**( - a )
     call gammaf ( c, g1 )
     call gammaf ( 1.0D+00 + a / 2.0D+00 - b, g2 )
     call gammaf ( 0.5D+00 + 0.5D+00 * a, g3 )
     hf = g0 * g1 / ( g2 * g3 )
     return

  else if ( l2 .or. l3 ) then

     if ( l2 ) then
        nm = int ( abs ( a ) )
     end if

     if ( l3 ) then
        nm = int ( abs ( b ) )
     end if

     hf = 1.0D+00
     r = 1.0D+00

     do k = 1, nm
        r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x
        hf = hf + r
     end do

     return

  else if ( l4 .or. l5 ) then

     if ( l4 ) then
        nm = int ( abs ( c - a ) )
     end if

     if ( l5 ) then
        nm = int ( abs ( c - b ) )
     end if

     hf = 1.0D+00
     r  = 1.0D+00
     do k = 1, nm
        r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x
        hf = hf + r
     end do
     hf = ( 1.0D+00 - x )**( c - a - b ) * hf
     return

  end if

  aa = a
  bb = b
  x1 = x
  !
  !  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
  !
  if ( x < 0.0D+00 ) then
     x = x / ( x - 1.0D+00 )
     if ( a < c .and. b < a .and. 0.0D+00 < b ) then
        a = bb
        b = aa
     end if
     b = c - b
  end if

  if ( 0.75D+00 <= x ) then

     gm = 0.0D+00

     if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

        m = int ( c - a - b )
        call gammaf ( a, ga )
        call gammaf ( b, gb )
        call gammaf ( c, gc )
        call gammaf ( a + m, gam )
        call gammaf ( b + m, gbm )
        call psi ( a, pa )
        call psi ( b, pb )

        if ( m /= 0 ) then
           gm = 1.0D+00
        end if

        do j = 1, abs ( m ) - 1
           gm = gm * j
        end do

        rm = 1.0D+00
        do j = 1, abs ( m )
           rm = rm * j
        end do

        f0 = 1.0D+00
        r0 = 1.0D+00
        r1 = 1.0D+00
        sp0 = 0.0D+00
        sp = 0.0D+00

        if ( 0 <= m ) then

           c0 = gm * gc / ( gam * gbm )
           c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

           do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                   / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
           end do

           do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
                   + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
           end do

           f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
           hw = f1

           do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
                   + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                 sm = sm + ( 1.0D+00 - a ) &
                      / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
                      + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

              r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
                   / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                 exit
              end if

              hw = f1

           end do

           hf = f0 * c0 + f1 * c1

        else if ( m < 0 ) then

           m = - m
           c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
           c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

           do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
                   / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
           end do

           do k = 1, m
              sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
           end do

           f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )

           do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) &
                   / ( k * ( a + k - 1.0D+00 ) ) &
                   + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                 sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                   / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                 exit
              end if

              hw = f1

           end do

           hf = f0 * c0 + f1 * c1

        end if

     else

        call gammaf ( a, ga )
        call gammaf ( b, gb )
        call gammaf ( c, gc )
        call gammaf ( c - a, gca )
        call gammaf ( c - b, gcb )
        call gammaf ( c - a - b, gcab )
        call gammaf ( a + b - c, gabc )
        c0 = gc * gcab / ( gca * gcb )
        c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
        hf = 0.0D+00
        r0 = c0
        r1 = c1

        do k = 1, 250

           r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

           r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
                / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

           hf = hf + r0 + r1

           if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
              exit
           end if

           hw = hf

        end do

        hf = hf + c0 + c1

     end if

  else

     a0 = 1.0D+00

     if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

        a0 = ( 1.0D+00 - x )**( c - a - b )
        a = c - a
        b = c - b

     end if

     hf = 1.0D+00
     r = 1.0D+00

     do k = 1, 250

        r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x

        hf = hf + r

        if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
           exit
        end if

        hw = hf

     end do

     hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
     x = x1
     c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
     hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFX - Warning!'
     write ( *, '(a)' ) '  A large number of iterations were needed.'
     write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end subroutine hygfx




  !*****************************************************************************80
  !
  !! HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, C, parameters.
  !
  !    Input, complex(8) Z, the argument.
  !
  !    Output, complex(8) ZHF, the value of F(a,b,c,z).
  !
subroutine hygfz ( a, b, c, z, zhf )
  !HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param c: in: parameter C
  !@param z: in: the argument
  !@param zhf: out: the value of F(a,b,c,z)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: aa
  real(8)    :: b
  real(8)    :: bb
  real(8)    :: c
  real(8)    :: ca
  real(8)    :: cb
  real(8)    :: el
  real(8)    :: eps
  real(8)    :: g0
  real(8)    :: g1
  real(8)    :: g2
  real(8)    :: g3
  real(8)    :: ga
  real(8)    :: gab
  real(8)    :: gabc
  real(8)    :: gam
  real(8)    :: gb
  real(8)    :: gba
  real(8)    :: gbm
  real(8)    :: gc
  real(8)    :: gca
  real(8)    :: gcab
  real(8)    :: gcb
  real(8)    :: gcbk
  real(8)    :: gm
  integer(4) :: j
  integer(4) :: k
  logical    :: l0
  logical    :: l1
  logical    :: l2
  logical    :: l3
  logical    :: l4
  logical    :: l5
  logical    :: l6
  integer(4) :: m
  integer(4) :: mab
  integer(4) :: mcab
  integer(4) :: nca
  integer(4) :: ncb
  integer(4) :: nm
  real(8)    :: pa
  real(8)    :: pac
  real(8)    :: pb
  real(8)    :: pca
  real(8)    :: pi
  real(8)    :: rk1
  real(8)    :: rk2
  real(8)    :: rm
  real(8)    :: sj1
  real(8)    :: sj2
  real(8)    :: sm
  real(8)    :: sp
  real(8)    :: sp0
  real(8)    :: sq
  real(8)    :: t0
  real(8)    :: w0
  real(8)    :: ws
  real(8)    :: x
  real(8)    :: y
  complex(8) :: z
  complex(8) :: z00
  complex(8) :: z1
  complex(8) :: zc0
  complex(8) :: zc1
  complex(8) :: zf0
  complex(8) :: zf1
  complex(8) :: zhf
  complex(8) :: zp
  complex(8) :: zp0
  complex(8) :: zr
  complex(8) :: zr0
  complex(8) :: zr1
  complex(8) :: zw

  x = real ( z, kind = 8 )
  y = imag ( z )
  eps = 1.0D-15
  l0 = c == int ( c ) .and. c < 0.0D+00
  l1 = abs ( 1.0D+00 - x ) < eps .and. y == 0.0D+00 .and. &
       c - a - b <= 0.0D+00
  l2 = abs ( z + 1.0D+00 ) < eps .and. &
       abs ( c - a + b - 1.0D+00 ) < eps
  l3 = a == int ( a ) .and. a < 0.0D+00
  l4 = b == int ( b ) .and. b < 0.0D+00
  l5 = c - a == int ( c - a ) .and. c - a <= 0.0D+00
  l6 = c - b == int ( c - b ) .and. c - b <= 0.0D+00
  aa = a
  bb = b
  a0 = abs ( z )
  if ( 0.95D+00 < a0 ) then
     eps = 1.0D-08
  end if
  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( l0 .or. l1 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFZ - Fatal error!'
     write ( *, '(a)' ) '  The hypergeometric series is divergent.'
     stop
  end if

  if ( a0 == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

     zhf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

  else if ( z == 1.0D+00.and. 0.0D+00 < c - a - b ) then

     call gammaf ( c, gc )
     call gammaf ( c - a - b, gcab )
     call gammaf ( c - a, gca )
     call gammaf ( c - b, gcb )
     zhf = gc * gcab / ( gca * gcb )

  else if ( l2 ) then

     g0 = sqrt ( pi ) * 2.0D+00 ** ( - a )
     call gammaf ( c, g1 )
     call gammaf ( 1.0D+00 + a / 2.0D+00 - b, g2 )
     call gammaf ( 0.5D+00 + 0.5D+00 * a, g3 )
     zhf = g0 * g1 / ( g2 * g3 )

  else if ( l3 .or. l4 ) then

     if ( l3 ) then
        nm = int ( abs ( a ) )
     end if

     if ( l4 ) then
        nm = int ( abs ( b ) )
     end if

     zhf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     zr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, nm
        zr = zr * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * z
        zhf = zhf + zr
     end do

  else if ( l5 .or. l6 ) then

     if ( l5 ) then
        nm = int ( abs ( c - a ) )
     end if

     if ( l6 ) then
        nm = int ( abs ( c - b ) )
     end if

     zhf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     zr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, nm
        zr = zr * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * z
        zhf = zhf + zr
     end do
     zhf = ( 1.0D+00 - z ) ** ( c - a - b ) * zhf

  else if ( a0 <= 1.0D+00 ) then

     if ( x < 0.0D+00 ) then

        z1 = z / ( z - 1.0D+00 )
        if ( a < c .and. b < a .and. 0.0D+00 < b ) then  
           a = bb
           b = aa
        end if
        zc0 = 1.0D+00 / ( ( 1.0D+00 - z ) ** a )
        zhf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        zr0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        do k = 1, 500
           zr0 = zr0 * ( a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
                / ( k * ( c + k - 1.0D+00 ) ) * z1
           zhf = zhf + zr0
           if ( abs ( zhf - zw ) < abs ( zhf ) * eps ) then
              exit
           end if
           zw = zhf
        end do

        zhf = zc0 * zhf

     else if ( 0.90D+00 <= a0 ) then

        gm = 0.0D+00
        mcab = int ( c - a - b + eps * sign ( 1.0D+00, c - a - b ) )

        if ( abs ( c - a - b - mcab ) < eps ) then

           m = int ( c - a - b )
           call gammaf ( a, ga )
           call gammaf ( b, gb )
           call gammaf ( c, gc )
           call gammaf ( a + m, gam )
           call gammaf ( b + m, gbm ) 
           call psi ( a, pa )
           call psi ( b, pb )
           if ( m /= 0 ) then
              gm = 1.0D+00
           end if
           do j = 1, abs ( m ) - 1
              gm = gm * j
           end do
           rm = 1.0D+00
           do j = 1, abs ( m )
              rm = rm * j
           end do
           zf0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           zr0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           zr1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
           sp0 = 0.0D+00
           sp = 0.0D+00

           if ( 0 <= m ) then

              zc0 = gm * gc / ( gam * gbm )
              zc1 = - gc * ( z - 1.0D+00 ) ** m / ( ga * gb * rm )
              do k = 1, m - 1
                 zr0 = zr0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                      / ( k * ( k - m ) ) * ( 1.0D+00 - z )
                 zf0 = zf0 + zr0
              end do
              do k = 1, m
                 sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
                      + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / k
              end do
              zf1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - z )
              do k = 1, 500
                 sp = sp + ( 1.0D+00 - a ) &
                      / ( k * ( a + k - 1.0D+00 ) ) + ( 1.0D+00 - b ) &
                      / ( k * ( b + k - 1.0D+00 ) )
                 sm = 0.0D+00
                 do j = 1, m
                    sm = sm + ( 1.0D+00 - a ) / ( ( j + k ) &
                         * ( a + j + k - 1.0D+00 ) ) &
                         + 1.0D+00 / ( b + j + k - 1.0D+00 )
                 end do
                 zp = pa + pb + 2.0D+00 * el + sp + sm  + log ( 1.0D+00 - z )
                 zr1 = zr1 * ( a + m + k - 1.0D+00 ) &
                      * ( b + m + k - 1.0D+00 ) / ( k * ( m + k ) ) &
                      * ( 1.0D+00 - z )
                 zf1 = zf1 + zr1 * zp
                 if ( abs ( zf1 - zw ) < abs ( zf1 ) * eps ) then
                    exit
                 end if
                 zw = zf1
              end do

              zhf = zf0 * zc0 + zf1 * zc1

           else if ( m < 0 ) then

              m = - m
              zc0 = gm * gc / ( ga * gb * ( 1.0D+00 - z ) ** m )
              zc1 = - ( - 1.0D+00 ) ** m * gc / ( gam * gbm * rm )
              do k = 1, m - 1
                 zr0 = zr0 * ( a - m + k - 1.0D+00 ) &
                      * ( b - m + k - 1.0D+00 ) / ( k * ( k - m ) ) &
                      * ( 1.0D+00 - z )
                 zf0 = zf0 + zr0
              end do

              do k = 1, m
                 sp0 = sp0 + 1.0D+00 / k
              end do

              zf1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - z )

              do k = 1, 500
                 sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
                      + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
                 sm = 0.0D+00
                 do j = 1, m
                    sm = sm + 1.0D+00 / ( j + k )
                 end do
                 zp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - z )
                 zr1 = zr1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                      / ( k * ( m + k ) ) * ( 1.0D+00 - z )
                 zf1 = zf1 + zr1 * zp
                 if ( abs ( zf1 - zw ) < abs ( zf1 ) * eps ) then
                    exit
                 end if
                 zw = zf1

              end do

              zhf = zf0 * zc0 + zf1 * zc1

           end if

        else

           call gammaf ( a, ga )
           call gammaf ( b, gb )
           call gammaf ( c, gc )
           call gammaf ( c - a, gca )
           call gammaf ( c - b, gcb )
           call gammaf ( c - a - b, gcab )
           call gammaf ( a + b - c, gabc )
           zc0 = gc * gcab / ( gca * gcb )
           zc1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - z ) ** ( c - a - b )
           zhf = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
           zr0 = zc0
           zr1 = zc1
           do k = 1, 500
              zr0 = zr0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                   / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - z )
              zr1 = zr1 * ( c - a + k - 1.0D+00 ) &
                   * ( c - b + k - 1.0D+00 ) / ( k * ( c - a - b + k ) ) &
                   * ( 1.0D+00 - z )
              zhf = zhf + zr0 + zr1
              if ( abs ( zhf - zw ) < abs ( zhf ) * eps ) then
                 exit
              end if
              zw = zhf
           end do

           zhf = zhf + zc0 + zc1

        end if

     else

        z00 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

        if ( c - a < a .and. c - b < b ) then
           z00 = ( 1.0D+00 - z ) ** ( c - a - b )
           a = c - a
           b = c - b
        end if

        zhf = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
        zr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )

        do k = 1, 1500
           zr = zr * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( c + k - 1.0D+00 ) ) * z
           zhf = zhf + zr
           if ( abs ( zhf - zw ) <= abs ( zhf ) * eps ) then
              exit
           end if
           zw = zhf
        end do

        zhf = z00 * zhf

     end if

  else if ( 1.0D+00 < a0 ) then

     mab = int ( a - b + eps * sign ( 1.0D+00, a - b ) )

     if ( abs ( a - b - mab ) < eps .and. a0 <= 1.1D+00 ) then
        b = b + eps
     end if

     if ( eps < abs ( a - b - mab ) ) then

        call gammaf ( a, ga )
        call gammaf ( b, gb )
        call gammaf ( c, gc )
        call gammaf ( a - b, gab )
        call gammaf ( b - a, gba )
        call gammaf ( c - a, gca )
        call gammaf ( c - b, gcb )
        zc0 = gc * gba / ( gca * gb * ( - z ) ** a )
        zc1 = gc * gab / ( gcb * ga * ( - z ) ** b )
        zr0 = zc0
        zr1 = zc1
        zhf = cmplx ( 0.0D+00, 0.0D+00 )

        do k = 1, 500
           zr0 = zr0 * ( a + k - 1.0D+00 ) * ( a - c + k ) &
                / ( ( a - b + k ) * k * z )
           zr1 = zr1 * ( b + k - 1.0D+00 ) * ( b - c + k ) &
                / ( ( b - a + k ) * k * z )
           zhf = zhf + zr0 + zr1
           if ( abs ( ( zhf - zw ) / zhf ) <= eps ) then
              exit
           end if
           zw = zhf
        end do

        zhf = zhf + zc0 + zc1

     else

        if ( a - b < 0.0D+00 ) then
           a = bb
           b = aa
        end if

        ca = c - a
        cb = c - b
        nca = int ( ca + eps * sign ( 1.0D+00, ca ) )
        ncb = int ( cb + eps * sign ( 1.0D+00, cb ) )

        if ( abs ( ca - nca ) < eps .or. abs ( cb - ncb ) < eps ) then
           c = c + eps
        end if

        call gammaf ( a, ga )
        call gammaf ( c, gc )
        call gammaf ( c - b, gcb )
        call psi ( a, pa )
        call psi ( c - a, pca )
        call psi ( a - c, pac )
        mab = int ( a - b + eps )
        zc0 = gc / ( ga * ( - z ) ** b )
        call gammaf ( a - b, gm )
        zf0 = gm / gcb * zc0
        zr = zc0
        do k = 1, mab - 1
           zr = zr * ( b + k - 1.0D+00 ) / ( k * z )
           t0 = a - b - k
           call gammaf ( t0, g0 )
           call gammaf ( c - b - k, gcbk )
           zf0 = zf0 + zr * g0 / gcbk
        end do

        if ( mab == 0 ) then
           zf0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
        end if

        zc1 = gc / ( ga * gcb * ( - z ) ** a )
        sp = -2.0D+00 * el - pa - pca
        do j = 1, mab
           sp = sp + 1.0D+00 / j
        end do
        zp0 = sp + log ( - z )
        sq = 1.0D+00
        do j = 1, mab
           sq = sq * ( b + j - 1.0D+00 ) * ( b - c + j ) / j
        end do
        zf1 = ( sq * zp0 ) * zc1
        zr = zc1
        rk1 = 1.0D+00
        sj1 = 0.0D+00

        do k = 1, 10000
           zr = zr / z
           rk1 = rk1 * ( b + k - 1.0D+00 ) * ( b - c + k ) / ( k * k )
           rk2 = rk1
           do j = k + 1, k + mab
              rk2 = rk2 * ( b + j - 1.0D+00 ) * ( b - c + j ) / j
           end do
           sj1 = sj1 + ( a - 1.0D+00 ) / ( k * ( a + k - 1.0D+00 ) ) &
                + ( a - c - 1.0D+00 ) / ( k * ( a - c + k - 1.0D+00 ) )
           sj2 = sj1
           do j = k + 1, k + mab
              sj2 = sj2 + 1.0D+00 / j
           end do
           zp = -2.0D+00 * el - pa - pac + sj2 - 1.0D+00 / ( k + a - c ) &
                - pi / tan ( pi * ( k + a - c ) ) + log ( - z ) 
           zf1 = zf1 + rk2 * zr * zp
           ws = abs ( zf1 )
           if ( abs ( ( ws - w0 ) / ws ) < eps ) then
              exit
           end if
           w0 = ws
        end do

        zhf = zf0 + zf1

     end if

  end if

  a = aa
  b = bb
  if ( 150 < k ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFZ - Warning!'
     write ( *, '(a)' ) '  The solution returned may have low accuracy.'
  end if

  return
end subroutine hygfz




  !*****************************************************************************80
  !
  !! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
  !
  !  Discussion:
  !
  !    This procedure computes modified Bessel functions I0(x), I1(x),
  !    K0(x) and K1(x), and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
  !    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
  !
subroutine ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
  !IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
  !
  !This procedure computes modified Bessel functions I0(x), I1(x),
  !K0(x) and K1(x), and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param bi0: out: I0(x)
  !@param di0: out: I0'(x)
  !@param bi1: out: I1(x)
  !@param di1: out: I1'(x)
  !@param bk0: out: K0(x)
  !@param dk0: out: K0'(x)
  !@param bk1: out: K1(x)
  !@param dk1: out: K1'(x)
  !
  implicit none

  real(8), save, dimension ( 12 ) :: a = (/ &
       0.125D+00, 7.03125D-02, &
       7.32421875D-02, 1.1215209960938D-01, &
       2.2710800170898D-01, 5.7250142097473D-01, &
       1.7277275025845D+00, 6.0740420012735D+00, &
       2.4380529699556D+01, 1.1001714026925D+02, &
       5.5133589612202D+02, 3.0380905109224D+03 /)
  real(8), save, dimension ( 8 ) :: a1 = (/ &
       0.125D+00, 0.2109375D+00, &
       1.0986328125D+00, 1.1775970458984D+01, &
       2.1461706161499D+02, 5.9511522710323D+03, &
       2.3347645606175D+05, 1.2312234987631D+07 /)
  real(8), save, dimension ( 12 ) :: b = (/ &
       -0.375D+00, -1.171875D-01, &
       -1.025390625D-01, -1.4419555664063D-01, &
       -2.7757644653320D-01, -6.7659258842468D-01, &
       -1.9935317337513D+00, -6.8839142681099D+00, &
       -2.7248827311269D+01, -1.2159789187654D+02, &
       -6.0384407670507D+02, -3.3022722944809D+03 /)
  real(8)    :: bi0
  real(8)    :: bi1
  real(8)    :: bk0
  real(8)    :: bk1
  real(8)    :: ca
  real(8)    :: cb
  real(8)    :: ct
  real(8)    :: di0
  real(8)    :: di1
  real(8)    :: dk0
  real(8)    :: dk1
  real(8)    :: el
  integer(4) :: k
  integer(4) :: k0
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: w0
  real(8)    :: ww
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xr
  real(8)    :: xr2

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  x2 = x * x

  if ( x == 0.0D+00 ) then

     bi0 = 1.0D+00
     bi1 = 0.0D+00
     bk0 = 1.0D+300
     bk1 = 1.0D+300
     di0 = 0.0D+00
     di1 = 0.5D+00
     dk0 = -1.0D+300
     dk1 = -1.0D+300
     return

  else if ( x <= 18.0D+00 ) then

     bi0 = 1.0D+00
     r = 1.0D+00
     do k = 1, 50
        r = 0.25D+00 * r * x2 / ( k * k )
        bi0 = bi0 + r
        if ( abs ( r / bi0 ) < 1.0D-15 ) then
           exit
        end if
     end do

     bi1 = 1.0D+00
     r = 1.0D+00
     do k = 1, 50
        r = 0.25D+00 * r * x2 / ( k * ( k + 1 ) )
        bi1 = bi1 + r
        if ( abs ( r / bi1 ) < 1.0D-15 ) then
           exit
        end if
     end do

     bi1 = 0.5D+00 * x * bi1

  else

     if ( x < 35.0D+00 ) then
        k0 = 12
     else if ( x < 50.0D+00 ) then
        k0 = 9
     else
        k0 = 7
     end if

     ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
     bi0 = 1.0D+00
     xr = 1.0D+00 / x
     do k = 1, k0
        bi0 = bi0 + a(k) * xr ** k
     end do
     bi0 = ca * bi0
     bi1 = 1.0D+00
     do k = 1, k0
        bi1 = bi1 + b(k) * xr ** k
     end do
     bi1 = ca * bi1

  end if

  if ( x <= 9.0D+00 ) then

     ct = - ( log ( x / 2.0D+00 ) + el )
     bk0 = 0.0D+00
     w0 = 0.0D+00
     r = 1.0D+00
     do k = 1, 50
        w0 = w0 + 1.0D+00 / k
        r = 0.25D+00 * r / ( k * k ) * x2
        bk0 = bk0 + r * ( w0 + ct )
        if ( abs ( ( bk0 - ww ) / bk0 ) < 1.0D-15 ) then
           exit
        end if
        ww = bk0
     end do

     bk0 = bk0 + ct

  else

     cb = 0.5D+00 / x
     xr2 = 1.0D+00 / x2
     bk0 = 1.0D+00
     do k = 1, 8
        bk0 = bk0 + a1(k) * xr2 ** k
     end do
     bk0 = cb * bk0 / bi0

  end if

  bk1 = ( 1.0D+00 / x - bi1 * bk0 ) / bi0
  di0 = bi1
  di1 = bi0 - bi1 / x
  dk0 = - bk1
  dk1 = - bk0 - bk1 / x

  return
end subroutine ik01a




  !*****************************************************************************80
  !
  !! IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    17 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
  !    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
  !
subroutine ik01b ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
  !IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param bi0: out: I0(x)
  !@param di0: out: I0'(x)
  !@param bi1: out: I1(x)
  !@param di1: out: I1'(x)
  !@param bk0: out: K0(x)
  !@param dk0: out: K0'(x)
  !@param bk1: out: K1(x)
  !@param dk1: out: K1'(x)
  !
  implicit none

  real(8) :: bi0
  real(8) :: bi1
  real(8) :: bk0
  real(8) :: bk1
  real(8) :: di0
  real(8) :: di1
  real(8) :: dk0
  real(8) :: dk1
  real(8) :: t
  real(8) :: t2
  real(8) :: x

  if ( x == 0.0D+00 ) then

     bi0 = 1.0D+00
     bi1 = 0.0D+00
     bk0 = 1.0D+300
     bk1 = 1.0D+300
     di0 = 0.0D+00
     di1 = 0.5D+00
     dk0 = -1.0D+300
     dk1 = -1.0D+300
     return

  else if ( x <= 3.75D+00 ) then

     t = x / 3.75D+00
     t2 = t * t

     bi0 = ((((( &
          0.0045813D+00   * t2 &
          + 0.0360768D+00 ) * t2 &
          + 0.2659732D+00 ) * t2 &
          + 1.2067492D+00 ) * t2 &
          + 3.0899424D+00 ) * t2 &
          + 3.5156229D+00 ) * t2 &
          + 1.0D+00

     bi1 = x * (((((( &
          0.00032411D+00   * t2 &
          + 0.00301532D+00 ) * t2 &
          + 0.02658733D+00 ) * t2 &
          + 0.15084934D+00 ) * t2 &
          + 0.51498869D+00 ) * t2 &
          + 0.87890594D+00 ) * t2 &
          + 0.5D+00 )

  else

     t = 3.75D+00 / x

     bi0 = (((((((( &
          0.00392377D+00   * t &
          - 0.01647633D+00 ) * t &
          + 0.02635537D+00 ) * t &
          - 0.02057706D+00 ) * t &
          + 0.916281D-02 ) * t &
          - 0.157565D-02 ) * t &
          + 0.225319D-02 ) * t &
          + 0.01328592D+00 ) * t &
          + 0.39894228D+00 ) * exp ( x ) / sqrt ( x )

     bi1 = (((((((( &
          - 0.420059D-02     * t &
          + 0.01787654D+00 ) * t &
          - 0.02895312D+00 ) * t &
          + 0.02282967D+00 ) * t &
          - 0.01031555D+00 ) * t &
          + 0.163801D-02 ) * t &
          - 0.00362018D+00 ) * t &
          - 0.03988024D+00 ) * t &
          + 0.39894228D+00 ) * exp ( x ) / sqrt ( x )

  end if

  if ( x <= 2.0D+00 ) then

     t = x / 2.0D+00
     t2 = t * t

     bk0 = ((((( &
          0.0000074D+00   * t2 &
          + 0.0001075D+00 ) * t2 &
          + 0.00262698D+00 ) * t2 &
          + 0.0348859D+00 ) * t2 &
          + 0.23069756D+00 ) * t2 &
          + 0.4227842D+00 ) * t2 &
          - 0.57721566D+00 - bi0 * log ( t )

     bk1 = (((((( &
          - 0.00004686D+00   * t2 &
          - 0.00110404D+00 ) * t2 &
          - 0.01919402D+00 ) * t2 &
          - 0.18156897D+00 ) * t2 &
          - 0.67278579D+00 ) * t2 &
          + 0.15443144D+00 ) * t2 &
          + 1.0D+00 ) / x + bi1 * log ( t )

  else

     t = 2.0D+00 / x
     t2 = t * t

     bk0 = (((((( &
          0.00053208D+00   * t &
          - 0.0025154D+00 )  * t &
          + 0.00587872D+00 ) * t &
          - 0.01062446D+00 ) * t &
          + 0.02189568D+00 ) * t &
          - 0.07832358D+00 ) * t &
          + 1.25331414D+00 ) * exp ( - x ) / sqrt ( x )

     bk1 = (((((( &
          - 0.00068245D+00   * t &
          + 0.00325614D+00 ) * t &
          - 0.00780353D+00 ) * t &
          + 0.01504268D+00 ) * t &
          - 0.0365562D+00  ) * t & 
          + 0.23498619D+00 ) * t &
          + 1.25331414D+00 ) * exp ( - x ) / sqrt ( x )

  end if

  di0 = bi1
  di1 = bi0 - bi1 / x
  dk0 = -bk1
  dk1 = -bk0 - bk1 / x

  return
end subroutine ik01b




  !*****************************************************************************80
  !
  !! IKNA compute Bessel function In(x) and Kn(x), and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of In(x) and Kn(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) BI(0:N), DI(0:N), BK(0:N), DK(0:N),
  !    the values of In(x), In'(x), Kn(x), Kn'(x).
  !
subroutine ikna ( n, x, nm, bi, di, bk, dk )
  !IKNA compute Bessel function In(x) and Kn(x), and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of In(x) and Kn(x)
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param bi: out: In(x)
  !@param di: out: In'(x)
  !@param bk: out: Kn(x)
  !@param dk: out: Kn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: bi(0:n)
  real(8)    :: bi0
  real(8)    :: bi1
  real(8)    :: bk(0:n)
  real(8)    :: bk0
  real(8)    :: bk1
  real(8)    :: di(0:n)
  real(8)    :: di0
  real(8)    :: di1
  real(8)    :: dk(0:n)
  real(8)    :: dk0
  real(8)    :: dk1
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: g
  real(8)    :: g0
  real(8)    :: g1
  real(8)    :: h
  real(8)    :: h0
  real(8)    :: h1
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: s0
  real(8)    :: x

  nm = n

  if ( x <= 1.0D-100 ) then
     do k = 0, n
        bi(k) = 0.0D+00
        di(k) = 0.0D+00
        bk(k) = 1.0D+300
        dk(k) = -1.0D+300
     end do
     bi(0) = 1.0D+00
     di(1) = 0.5D+00
     return
  end if

  call ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
  bi(0) = bi0
  bi(1) = bi1
  bk(0) = bk0
  bk(1) = bk1
  di(0) = di0
  di(1) = di1
  dk(0) = dk0
  dk(1) = dk1

  if ( n <= 1 ) then
     return
  end if

  if ( 40.0D+00 < x .and. n < int ( 0.25D+00 * x ) ) then

     h0 = bi0
     h1 = bi1
     do k = 2, n
        h = -2.0D+00 * ( k - 1.0D+00 ) / x * h1 + h0
        bi(k) = h
        h0 = h1
        h1 = h
     end do

  else

     m = msta1 ( x, 200 )

     if ( m < n ) then
        nm = m
     else
        m = msta2 ( x, n, 15 )
     end if

     f0 = 0.0D+00
     f1 = 1.0D-100
     do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x + f0
        if ( k <= nm ) then
           bi(k) = f
        end if
        f0 = f1
        f1 = f
     end do
     s0 = bi0 / f
     do k = 0, nm
        bi(k) = s0 * bi(k)
     end do
  end if

  g0 = bk0
  g1 = bk1
  do k = 2, nm
     g = 2.0D+00 * ( k - 1.0D+00 ) / x * g1 + g0
     bk(k) = g
     g0 = g1
     g1 = g
  end do

  do k = 2, nm
     di(k) = bi(k-1) - k / x * bi(k)
     dk(k) = - bk(k-1) - k / x * bk(k)
  end do

  return
end subroutine ikna




  !*****************************************************************************80
  !
  !! IKNB compute Bessel function In(x) and Kn(x).
  !
  !  Discussion:
  !
  !    Compute modified Bessel functions In(x) and Kn(x),
  !    and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    17 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of In(x) and Kn(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) BI(0:N), DI(0:N), BK(0:N), DK(0:N),
  !    the values of In(x), In'(x), Kn(x), Kn'(x).
  !
subroutine iknb ( n, x, nm, bi, di, bk, dk )
  !IKNB compute Bessel function In(x) and Kn(x).
  !
  !Compute modified Bessel functions In(x) and Kn(x),
  !and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of In(x) and Kn(x)
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param bi: out: In(x)
  !@param di: out: In'(x)
  !@param bk: out: Kn(x)
  !@param dk: out: Kn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: a0
  real(8)    :: bi(0:n)
  real(8)    :: bk(0:n)
  real(8)    :: bkl
  real(8)    :: bs
  real(8)    :: di(0:n)
  real(8)    :: dk(0:n)
  real(8)    :: el
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: g
  real(8)    :: g0
  real(8)    :: g1
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: l
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: s0
  real(8)    :: sk0
  real(8)    :: vt
  real(8)    :: x

  pi = 3.141592653589793D+00
  el = 0.5772156649015329d0
  nm = n

  if ( x <= 1.0D-100 ) then
     do k = 0, n
        bi(k) = 0.0D+00
        di(k) = 0.0D+00
        bk(k) = 1.0D+300
        dk(k) = -1.0D+300
     end do
     bi(0) = 1.0D+00
     di(1) = 0.5D+00
     return
  end if

  if ( n == 0 ) then
     nm = 1
  end if

  m = msta1 ( x, 200 )
  if ( m < nm ) then
     nm = m
  else
     m = msta2 ( x, nm, 15 )
  end if

  bs = 0.0D+00
  sk0 = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-100
  do k = m, 0, -1
     f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 + f0
     if ( k <= nm ) then
        bi(k) = f
     end if
     if ( k /= 0 .and. k == 2 * int ( k / 2 ) ) then
        sk0 = sk0 + 4.0D+00 * f / k
     end if
     bs = bs + 2.0D+00 * f
     f0 = f1
     f1 = f
  end do

  s0 = exp ( x ) / ( bs - f )
  do k = 0, nm
     bi(k) = s0 * bi(k)
  end do

  if ( x <= 8.0D+00 ) then
     bk(0) = - ( log ( 0.5D+00 * x ) + el ) * bi(0) + s0 * sk0
     bk(1) = ( 1.0D+00 / x - bi(1) * bk(0) ) / bi(0)
  else
     a0 = sqrt ( pi / ( 2.0D+00 * x ) ) * exp ( - x ) 

     if ( x < 25.0D+00 ) then
        k0 = 16
     else if ( x < 80.0D+00 ) then
        k0 = 10
     else if ( x < 200.0D+00 ) then
        k0 = 8
     else
        k0 = 6
     end if

     do l = 0, 1
        bkl = 1.0D+00
        vt = 4.0D+00 * l
        r = 1.0D+00
        do k = 1, k0
           r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
           bkl = bkl + r
        end do
        bk(l) = a0 * bkl
     end do
  end if

  g0 = bk(0)
  g1 = bk(1)
  do k = 2, nm
     g = 2.0D+00 * ( k - 1.0D+00 ) / x * g1 + g0
     bk(k) = g
     g0 = g1
     g1 = g
  end do

  di(0) = bi(1)
  dk(0) = -bk(1)
  do k = 1, nm
     di(k) = bi(k-1) - k / x * bi(k)
     dk(k) = -bk(k-1) - k / x * bk(k)
  end do

  return
end subroutine iknb




  !*****************************************************************************80
  !
  !! IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    17 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Iv(x) and Kv(x).
  !    V = N + V0.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) BI(0:N), DI(0:N), BK(0:N), DK(0:N), the
  !    values of In+v0(x), In+v0'(x), Kn+v0(x), Kn+v0'(x).
  !
subroutine ikv ( v, x, vm, bi, di, bk, dk )
  !IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Iv(x) and Kv(x). V = N + V0
  !@param x: in: the argument
  !@param vm: out: the highest order computed
  !@param bi: out: In+v0(x)
  !@param di: out: In+v0'(x)
  !@param bk: out: Kn+v0(x)
  !@param dk: out: Kn+v0'(x)
  !
  implicit none

  real(8)    :: a1
  real(8)    :: a2
  real(8)    :: bi(0:*)
  real(8)    :: bi0
  real(8)    :: bk(0:*)
  real(8)    :: bk0
  real(8)    :: bk1
  real(8)    :: bk2
  real(8)    :: ca
  real(8)    :: cb
  real(8)    :: cs
  real(8)    :: ct
  real(8)    :: di(0:*)
  real(8)    :: dk(0:*)
  real(8)    :: f
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: gan
  real(8)    :: gap
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: piv
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: sum
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: v0n
  real(8)    :: v0p
  real(8)    :: vm
  real(8)    :: vt
  real(8)    :: w0
  real(8)    :: wa
  real(8)    :: ww
  real(8)    :: x
  real(8)    :: x2

  pi = 3.141592653589793D+00
  x2 = x * x
  n = int ( v )
  v0 = v - n
  if ( n == 0 ) then
     n = 1
  end if

  if ( x < 1.0D-100 ) then

     do k = 0, n
        bi(k) = 0.0D+00
        di(k) = 0.0D+00
        bk(k) = -1.0D+300
        dk(k) = 1.0D+300
     end do

     if ( v == 0.0D+00 ) then
        bi(0) = 1.0D+00
        di(1) = 0.5D+00
     end if

     vm = v
     return

  end if

  piv = pi * v0
  vt = 4.0D+00 * v0 * v0

  if ( v0 == 0.0D+00 ) then
     a1 = 1.0D+00
  else
     v0p = 1.0D+00 + v0
     call gammaf ( v0p, gap )
     a1 = ( 0.5D+00 * x ) ** v0 / gap
  end if

  if ( x < 35.0D+00 ) then
     k0 = 14
  else if ( x < 50.0D+00 ) then
     k0 = 10
  else
     k0 = 8
  end if

  if ( x <= 18.0D+00 ) then

     bi0 = 1.0D+00
     r = 1.0D+00
     do k = 1, 30
        r = 0.25D+00 * r * x2 / ( k * ( k + v0 ) )
        bi0 = bi0 + r
        if ( abs ( r / bi0 ) < 1.0D-15 ) then
           exit
        end if
     end do

     bi0 = bi0 * a1

  else

     ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
     sum = 1.0D+00
     r = 1.0D+00
     do k = 1, k0
        r = -0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
        sum = sum + r
     end do
     bi0 = ca * sum

  end if

  m = msta1 ( x, 200 )

  if ( m < n ) then
     n = m
  else
     m = msta2 ( x, n, 15 )
  end if

  f2 = 0.0D+00
  f1 = 1.0D-100
  do k = m, 0, -1
     f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 + f2
     if ( k <= n ) then
        bi(k) = f
     end if
     f2 = f1
     f1 = f
  end do

  cs = bi0 / f
  do k = 0, n
     bi(k) = cs * bi(k)
  end do

  di(0) = v0 / x * bi(0) + bi(1)
  do k = 1, n
     di(k) = - ( k + v0 ) / x * bi(k) + bi(k-1)
  end do

  if ( x <= 9.0D+00 ) then

     if ( v0 == 0.0D+00 ) then

        ct = - log ( 0.5D+00 * x ) - 0.5772156649015329D+00
        cs = 0.0D+00
        w0 = 0.0D+00
        r = 1.0D+00
        do k = 1, 50
           w0 = w0 + 1.0D+00 / k
           r = 0.25D+00 * r / ( k * k ) * x2
           cs = cs + r * ( w0 + ct )
           wa = abs ( cs )
           if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
              exit
           end if
           ww = wa
        end do

        bk0 = ct + cs

     else

        v0n = 1.0D+00 - v0
        call gammaf ( v0n, gan )
        a2 = 1.0D+00 / ( gan * ( 0.5D+00 * x ) ** v0 )
        a1 = ( 0.5D+00 * x ) ** v0 / gap
        sum = a2 - a1
        r1 = 1.0D+00
        r2 = 1.0D+00
        do k = 1, 120
           r1 = 0.25D+00 * r1 * x2 / ( k * ( k - v0 ) )
           r2 = 0.25D+00 * r2 * x2 / ( k * ( k + v0 ) )
           sum = sum + a2 * r1 - a1 * r2
           wa = abs ( sum )
           if ( abs ( ( wa - ww ) / wa ) < 1.0D-15 ) then
              exit
           end if
           ww = wa
        end do

        bk0 = 0.5D+00 * pi * sum / sin ( piv )

     end if

  else

     cb = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )
     sum = 1.0D+00
     r = 1.0D+00
     do k = 1, k0
        r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
        sum = sum + r
     end do
     bk0 = cb * sum

  end if

  bk1 = ( 1.0D+00 / x - bi(1) * bk0 ) / bi(0)
  bk(0) = bk0
  bk(1) = bk1
  do k = 2, n
     bk2 = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * bk1 + bk0
     bk(k) = bk2
     bk0 = bk1
     bk1 = bk2
  end do

  dk(0) = v0 / x * bk(0) - bk(1)
  do k = 1, n
     dk(k) = - ( k + v0 ) / x * bk(k) - bk(k-1)
  end do

  vm = n + v0

  return
end subroutine ikv




  !*****************************************************************************80
  !
  !! INCOB computes the incomplete beta function Ix(a,b).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, B, parameters.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BIX, the function value.
  !
subroutine incob ( a, b, x, bix )
  !INCOB computes the incomplete beta function Ix(a,b).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: parameter A
  !@param b: in: parameter B
  !@param x: in: the argument
  !@param bix: out: the function value
  !
  implicit none

  real(8)    :: a
  real(8)    :: b
  real(8)    :: bix
  real(8)    :: bt
  real(8)    :: dk(51)
  real(8)    :: fk(51)
  integer(4) :: k
  real(8)    :: s0
  real(8)    :: t1
  real(8)    :: t2
  real(8)    :: ta
  real(8)    :: tb
  real(8)    :: x

  s0 = ( a + 1.0D+00 ) / ( a + b + 2.0D+00 )
  call betaf ( a, b, bt )

  if ( x <= s0 ) then

     do k = 1, 20
        dk(2*k) = k * ( b - k ) * x / &
             ( a + 2.0D+00 * k - 1.0D+00 ) / ( a + 2.0D+00 * k )
     end do

     do k = 0, 20
        dk(2*k+1) = - ( a + k ) * ( a + b + k ) * x &
             / ( a + 2.0D+00 * k ) / ( a + 2.0D+00 * k + 1.0D+00 )
     end do

     t1 = 0.0D+00
     do k = 20, 1, -1
        t1 = dk(k) / ( 1.0D+00 + t1 )
     end do
     ta = 1.0D+00 / ( 1.0D+00 + t1 )
     bix = x ** a * ( 1.0D+00 - x ) ** b / ( a * bt ) * ta

  else

     do k = 1, 20
        fk(2*k) = k * ( a - k ) * ( 1.0D+00 - x ) &
             / ( b + 2.0D+00 * k - 1.0D+00 ) / ( b + 2.0D+00 * k )
     end do

     do k = 0,20
        fk(2*k+1) = - ( b + k ) * ( a + b + k ) * ( 1.0D+00 - x ) &
             / ( b + 2.0D+00 * k ) / ( b + 2.0D+00 * k + 1.0D+00 )
     end do

     t2 = 0.0D+00
     do k = 20, 1, -1
        t2 = fk(k) / ( 1.0D+00 + t2 )
     end do
     tb = 1.0D+00 / ( 1.0D+00 + t2 )
     bix = 1.0D+00 - x ** a * ( 1.0D+00 - x ) ** b / ( b * bt ) * tb

  end if

  return
end subroutine incob




  !*****************************************************************************80
  !
  !! INCOG computes the incomplete gamma function r(a,x), Γ(a,x), P(a,x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, the parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) GIN, GIM, GIP, the values of
  !    r(a,x), Γ(a,x), P(a,x).
  !
subroutine incog ( a, x, gin, gim, gip )
  !INCOG computes the incomplete gamma function r(a,x), Gamma(a,x), P(a,x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: the parameter
  !@param x: in: the argument
  !@param gin: out: r(a,x)
  !@param gim: out: Gamma(a,x)
  !@param gip: out: P(a,x)
  !
  implicit none

  real(8)    :: a
  real(8)    :: ga
  real(8)    :: gim
  real(8)    :: gin
  real(8)    :: gip
  integer(4) :: k
  real(8)    :: r
  real(8)    :: s
  real(8)    :: t0
  real(8)    :: x
  real(8)    :: xam

  xam = -  x + a * log ( x )

  if ( 700.0D+00 < xam .or. 170.0D+00 < a ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'INCOG - Fatal error!'
     write ( *, '(a)' ) '  A and/or X is too large!'
     stop
  end if

  if ( x == 0.0D+00 ) then

     gin = 0.0D+00
     call gammaf ( a, ga )
     gim = ga
     gip = 0.0D+00

  else if ( x <= 1.0D+00 + a ) then

     s = 1.0D+00 / a
     r = s
     do k = 1, 60
        r = r * x / ( a + k )
        s = s + r
        if ( abs ( r / s ) < 1.0D-15 ) then
           exit
        end if
     end do

     gin = exp ( xam ) * s
     call gammaf ( a, ga )
     gip = gin / ga
     gim = ga - gin

  else if ( 1.0D+00 + a < x ) then

     t0 = 0.0D+00
     do k = 60, 1, -1
        t0 = ( k - a ) / ( 1.0D+00 + k / ( x + t0 ) )
     end do
     gim = exp ( xam ) / ( x + t0 )
     call gammaf ( a, ga )
     gin = ga - gim
     gip = 1.0D+00 - gim / ga

  end if

  return
end subroutine incog




  !****************************************************************************80
  !
  !! ITAIRY computes the integrals of Airy functions.
  !
  !  Discussion:
  !
  !    Compute the integrals of Airy functions with respect to t,
  !    from 0 and x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    19 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) APT, BPT, ANT, BNT, the integrals, from 0 to x,
  !    of Ai(t), Bi(t), Ai(-t), and Bi(-t).
  !       
subroutine itairy ( x, apt, bpt, ant, bnt )
  !ITAIRY computes the integrals of Airy functions.
  !
  !Compute the integrals of Airy functions with respect to t,
  !from 0 and x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param apt: out: the integral from 0 to x of Ai(t)
  !@param bpt: out: the integral from 0 to x of Bi(t)
  !@param ant: out: the integral from 0 to x of Ai(-t)
  !@param bnt: out: the integral from 0 to x of Bi(-t)
  !
  implicit none

  real(8), save, dimension ( 16 ) :: a = (/ &
       0.569444444444444D+00, 0.891300154320988D+00, &
       0.226624344493027D+01, 0.798950124766861D+01, &
       0.360688546785343D+02, 0.198670292131169D+03, &
       0.129223456582211D+04, 0.969483869669600D+04, &
       0.824184704952483D+05, 0.783031092490225D+06, &
       0.822210493622814D+07, 0.945557399360556D+08, &
       0.118195595640730D+10, 0.159564653040121D+11, &
       0.231369166433050D+12, 0.358622522796969D+13 /)
  real(8)    :: ant
  real(8)    :: apt
  real(8)    :: bnt
  real(8)    :: bpt
  real(8)    :: c1
  real(8)    :: c2
  real(8)    :: eps
  real(8)    :: fx
  real(8)    :: gx
  integer(4) :: k
  integer(4) :: l
  real(8)    :: pi
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: q2
  real(8)    :: r
  real(8)    :: sr3
  real(8)    :: su1
  real(8)    :: su2
  real(8)    :: su3
  real(8)    :: su4
  real(8)    :: su5
  real(8)    :: su6
  real(8)    :: x
  real(8)    :: xe
  real(8)    :: xp6
  real(8)    :: xr1
  real(8)    :: xr2

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  c1 = 0.355028053887817D+00
  c2 = 0.258819403792807D+00
  sr3 = 1.732050807568877D+00

  if ( x == 0.0D+00 ) then

     apt = 0.0D+00
     bpt = 0.0D+00
     ant = 0.0D+00
     bnt = 0.0D+00

  else

     if ( abs ( x ) <= 9.25D+00 ) then

        do l = 0, 1

           x = ( -1.0D+00 ) ** l * x
           fx = x
           r = x

           do k = 1, 40
              r = r * ( 3.0D+00 * k - 2.0D+00 ) &
                   / ( 3.0D+00 * k + 1.0D+00 ) * x / ( 3.0D+00 * k ) &
                   * x / ( 3.0D+00 * k - 1.0D+00 ) * x 
              fx = fx + r
              if ( abs ( r ) < abs ( fx ) * eps ) then
                 exit
              end if
           end do

           gx = 0.5D+00 * x * x
           r = gx

           do k = 1, 40
              r = r * ( 3.0D+00 * k - 1.0D+00 ) &
                   / ( 3.0D+00 * k + 2.0D+00 ) * x / ( 3.0D+00 * k ) * x &
                   / ( 3.0D+00 * k + 1.0D+00 ) * x
              gx = gx + r
              if ( abs ( r ) < abs ( gx ) * eps ) then
                 exit
              end if
           end do

           ant = c1 * fx - c2 * gx
           bnt = sr3 * ( c1 * fx + c2 * gx )

           if ( l == 0 ) then
              apt = ant
              bpt = bnt
           else
              ant = -ant
              bnt = -bnt
              x = -x
           end if

        end do

     else

        q2 = 1.414213562373095D+00
        q0 = 0.3333333333333333D+00
        q1 = 0.6666666666666667D+00
        xe = x * sqrt ( x ) / 1.5D+00
        xp6 = 1.0D+00 / sqrt ( 6.0D+00 * pi * xe )
        su1 = 1.0D+00
        r = 1.0D+00
        xr1 = 1.0D+00 / xe
        do k = 1, 16
           r = - r * xr1
           su1 = su1 + a(k) * r
        end do
        su2 = 1.0D+00
        r = 1.0D+00
        do k = 1, 16
           r = r * xr1
           su2 = su2 + a(k) * r
        end do

        apt = q0 - exp ( - xe ) * xp6 * su1
        bpt = 2.0D+00 * exp ( xe ) * xp6 * su2
        su3 = 1.0D+00
        r = 1.0D+00
        xr2 = 1.0D+00 / ( xe * xe )
        do k = 1, 8
           r = - r * xr2
           su3 = su3 + a(2*k) * r
        end do
        su4 = a(1) * xr1
        r = xr1
        do k = 1, 7
           r = -r * xr2
           su4 = su4 + a(2*k+1) * r
        end do
        su5 = su3 + su4
        su6 = su3 - su4
        ant = q1 - q2 * xp6 * ( su5 * cos ( xe ) - su6 * sin ( xe ) )
        bnt = q2 * xp6 * ( su5 * sin ( xe ) + su6 * cos ( xe ) )

     end if

  end if

  return
end subroutine itairy




  !*****************************************************************************80
  !
  !! ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
  !
  !  Discussion:
  !
  !    This procedure integrates modified Bessel functions I0(t) and
  !    K0(t) with respect to t from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TI, TK, the integrals of I0(t) and K0(t)
  !    from 0 to X.
  !
subroutine itika ( x, ti, tk )
  !ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
  !
  !This procedure integrates modified Bessel functions I0(t) and
  !K0(t) with respect to t from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param ti: out: the integral of I0(t) from 0 to :math:`X`
  !@param tk: out: the integral of K0(t) from 0 to :math:`X`
  !
  implicit none

  real(8), save, dimension ( 10 ) :: a = (/ &
       0.625D+00,           1.0078125D+00, &
       2.5927734375D+00,    9.1868591308594D+00, &
       4.1567974090576D+01, 2.2919635891914D+02, &
       1.491504060477D+03,  1.1192354495579D+04, &
       9.515939374212D+04,  9.0412425769041D+05 /)
  real(8)    :: b1
  real(8)    :: b2
  real(8)    :: e0
  real(8)    :: el
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: rc1
  real(8)    :: rc2
  real(8)    :: rs
  real(8)    :: ti
  real(8)    :: tk
  real(8)    :: tw
  real(8)    :: x
  real(8)    :: x2

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( x == 0.0D+00 ) then

     ti = 0.0D+00
     tk = 0.0D+00
     return

  else if ( x < 20.0D+00 ) then

     x2 = x * x
     ti = 1.0D+00
     r = 1.0D+00
     do k = 1, 50
        r = 0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) &
             / ( k * k ) * x2
        ti = ti + r
        if ( abs ( r / ti ) < 1.0D-12 ) then
           exit
        end if
     end do

     ti = ti * x

  else

     ti = 1.0D+00
     r = 1.0D+00
     do k = 1, 10
        r = r / x
        ti = ti + a(k) * r
     end do
     rc1 = 1.0D+00 / sqrt ( 2.0D+00 * pi * x )
     ti = rc1 * exp ( x ) * ti

  end if

  if ( x < 12.0D+00 ) then

     e0 = el + log ( x / 2.0D+00 )
     b1 = 1.0D+00 - e0
     b2 = 0.0D+00
     rs = 0.0D+00
     r = 1.0D+00
     do k = 1, 50
        r = 0.25D+00 * r * ( 2 * k - 1.0D+00 ) &
             / ( 2 * k + 1.0D+00 ) / ( k * k ) * x2
        b1 = b1 + r * ( 1.0D+00 / ( 2 * k + 1 ) - e0 )
        rs = rs + 1.0D+00 / k
        b2 = b2 + r * rs
        tk = b1 + b2
        if ( abs ( ( tk - tw ) / tk ) < 1.0D-12 ) then
           exit
        end if
        tw = tk
     end do

     tk = tk * x

  else

     tk = 1.0D+00
     r = 1.0D+00
     do k = 1, 10
        r = -r / x
        tk = tk + a(k) * r
     end do
     rc2 = sqrt ( pi / ( 2.0D+00 * x ) )
     tk = pi / 2.0D+00 - rc2 * tk * exp ( - x )

  end if

  return
end subroutine itika




  !*****************************************************************************80
  !
  !! ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
  !
  !  Discussion:
  !
  !    This procedure integrates Bessel functions I0(t) and K0(t)
  !    with respect to t from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    24 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TI, TK, the integral of I0(t) and K0(t)
  !    from 0 to X.
  !
subroutine itikb ( x, ti, tk )
  !ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
  !
  !This procedure integrates Bessel functions I0(t) and K0(t)
  !with respect to t from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param ti: out: the integral of I0(t) from 0 to X
  !@param tk: out: the integral of K0(t) from 0 to X
  !
  implicit none

  real(8) :: pi
  real(8) :: t
  real(8) :: t1
  real(8) :: ti
  real(8) :: tk
  real(8) :: x

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

     ti = 0.0D+00

  else if ( x < 5.0D+00 ) then

     t1 = x / 5.0D+00
     t = t1 * t1
     ti = (((((((( &
          0.59434D-03 * t &
          + 0.4500642D-02 ) * t &
          + 0.044686921D+00 ) * t &
          + 0.300704878D+00 ) * t &
          + 1.471860153D+00 ) * t &
          + 4.844024624D+00 ) * t &
          + 9.765629849D+00 ) * t &
          +10.416666367D+00 ) * t &
          + 5.0D+00 ) * t1

  else if ( 5.0D+00 <= x .and. x <= 8.0D+00 ) then

     t = 5.0D+00 / x
     ti = ((( &
          - 0.015166D+00 * t &
          - 0.0202292D+00 ) * t &
          + 0.1294122D+00 ) * t &
          - 0.0302912D+00 ) * t &
          + 0.4161224D+00
     ti = ti * exp ( x ) / sqrt ( x )

  else

     t = 8.0D+00 / x
     ti = ((((( &
          - 0.0073995D+00 * t &
          + 0.017744D+00 ) * t &
          - 0.0114858D+00 ) * t &
          + 0.55956D-02 ) * t &
          + 0.59191D-02 ) * t &
          + 0.0311734D+00 ) * t &
          + 0.3989423D+00
     ti = ti * exp ( x ) / sqrt ( x )

  end if

  if ( x == 0.0D+00 ) then

     tk = 0.0D+00

  else if ( x <= 2.0D+00 ) then

     t1 = x / 2.0D+00
     t = t1 * t1
     tk = (((((( &
          0.116D-05        * t &
          + 0.2069D-04 )     * t &
          + 0.62664D-03 )    * t &
          + 0.01110118D+00 ) * t &
          + 0.11227902D+00 ) * t &
          + 0.50407836D+00 ) * t &
          + 0.84556868D+00 ) * t1
     tk = tk - log ( x / 2.0D+00 ) * ti

  else if ( 2.0D+00 < x .and. x <= 4.0D+00 ) then

     t = 2.0D+00 / x
     tk = ((( &
          0.0160395D+00   * t &
          - 0.0781715D+00 ) * t &
          + 0.185984D+00 )  * t &
          - 0.3584641D+00 ) * t &
          + 1.2494934D+00
     tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

  else if ( 4.0D+00 < x .and. x <= 7.0D+00 ) then

     t = 4.0D+00 / x
     tk = ((((( &
          0.37128D-02 * t &
          - 0.0158449D+00 ) * t &
          + 0.0320504D+00 ) * t &
          - 0.0481455D+00 ) * t &
          + 0.0787284D+00 ) * t &
          - 0.1958273D+00 ) * t &
          + 1.2533141D+00
     tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

  else

     t = 7.0D+00 / x
     tk = ((((( &
          0.33934D-03      * t &
          - 0.163271D-02 )   * t &
          + 0.417454D-02 )   * t &
          - 0.933944D-02 )   * t &
          + 0.02576646D+00 ) * t &
          - 0.11190289D+00 ) * t &
          + 1.25331414D+00
     tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

  end if

  return
end subroutine itikb




  !*****************************************************************************80
  !
  !! ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
  !
  !  Discussion:
  !
  !    This procedure integrates Bessel functions J0(t) and Y0(t) with
  !    respect to t from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TJ, TY, the integrals of J0(t) and Y0(t) 
  !    from 0 to x.
  !
subroutine itjya ( x, tj, ty )
  !ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
  !
  !This procedure integrates Bessel functions J0(t) and Y0(t) with
  !respect to t from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param tj: out: the integral of J0(t) from 0 to x
  !@param ty: out: the integral of Y0(t) from 0 to x
  !
  implicit none

  real(8)    :: a(18)
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: af
  real(8)    :: bf
  real(8)    :: bg
  real(8)    :: el
  real(8)    :: eps
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: r2
  real(8)    :: rc
  real(8)    :: rs
  real(8)    :: tj
  real(8)    :: ty
  real(8)    :: ty1
  real(8)    :: ty2
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xp

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  eps = 1.0D-12

  if ( x == 0.0D+00 ) then

     tj = 0.0D+00
     ty = 0.0D+00

  else if ( x <= 20.0D+00 ) then

     x2 = x * x
     tj = x
     r = x
     do k = 1, 60
        r = -0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) &
             / ( k * k ) * x2
        tj = tj + r
        if ( abs ( r ) < abs ( tj ) * eps ) then
           exit
        end if
     end do

     ty1 = ( el + log ( x / 2.0D+00 ) ) * tj
     rs = 0.0D+00
     ty2 = 1.0D+00
     r = 1.0D+00

     do k = 1, 60
        r = -0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) &
             / ( k * k ) * x2
        rs = rs + 1.0D+00 / k
        r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k + 1.0D+00 ) )
        ty2 = ty2 + r2
        if ( abs ( r2 ) < abs ( ty2 ) * eps ) then
           exit
        end if
     end do

     ty = ( ty1 - x * ty2 ) * 2.0D+00 / pi

  else

     a0 = 1.0D+00
     a1 = 5.0D+00 / 8.0D+00
     a(1) = a1

     do k = 1, 16
        af = ( ( 1.5D+00 * ( k + 0.5D+00 ) * ( k + 5.0D+00 / 6.0D+00 ) &
             * a1 - 0.5D+00 * ( k + 0.5D+00 ) * ( k + 0.5D+00 )  &
             * ( k - 0.5D+00 ) * a0 ) ) / ( k + 1.0D+00 )
        a(k+1) = af
        a0 = a1
        a1 = af
     end do

     bf = 1.0D+00
     r = 1.0D+00
     do k = 1, 8
        r = -r / ( x * x )
        bf = bf + a(2*k) * r
     end do
     bg = a(1) / x
     r = 1.0D+00 / x
     do k = 1, 8
        r = -r / ( x * x )
        bg = bg + a(2*k+1) * r
     end do
     xp = x + 0.25D+00 * pi
     rc = sqrt ( 2.0D+00 / ( pi * x ) )
     tj = 1.0D+00 - rc * ( bf * cos ( xp ) + bg * sin ( xp ) )
     ty = rc * ( bg * cos ( xp ) - bf * sin ( xp ) )

  end if

  return
end subroutine itjya




  !*****************************************************************************80
  !
  !! ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
  !
  !  Discussion:
  !
  !    This procedure integrates Bessel functions J0(t) and Y0(t)
  !    with respect to t from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TJ, TY, the integrals of J0(t) and Y0(t) 
  !    from 0 to x.
  !
subroutine itjyb ( x, tj, ty )
  !ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
  !
  !This procedure integrates Bessel functions J0(t) and Y0(t)
  !with respect to t from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param tj: out: the integral of J0(t) from 0 to x
  !@param ty: out: the integral of Y0(t) from 0 to x
  !
  implicit none

  real(8) :: f0
  real(8) :: g0
  real(8) :: pi
  real(8) :: t
  real(8) :: tj
  real(8) :: ty
  real(8) :: x
  real(8) :: x1
  real(8) :: xt

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

     tj = 0.0D+00
     ty = 0.0D+00

  else if ( x <= 4.0D+00 ) then

     x1 = x / 4.0D+00
     t = x1 * x1

     tj = ((((((( &
          - 0.133718D-03      * t &
          + 0.2362211D-02 )   * t &
          - 0.025791036D+00 ) * t &
          + 0.197492634D+00 ) * t &
          - 1.015860606D+00 ) * t &
          + 3.199997842D+00 ) * t &
          - 5.333333161D+00 ) * t &
          + 4.0D+00 ) * x1

     ty = (((((((( &
          0.13351D-04       * t &
          - 0.235002D-03 )    * t &
          + 0.3034322d-02 )   * t &
          - 0.029600855D+00 ) * t &
          + 0.203380298D+00 ) * t &
          - 0.904755062D+00 ) * t &
          + 2.287317974D+00 ) * t &
          - 2.567250468D+00 ) * t &
          + 1.076611469D+00 ) * x1

     ty = 2.0D+00 / pi * log ( x / 2.0D+00 ) * tj - ty

  else if ( x <= 8.0D+00 ) then

     xt = x - 0.25D+00 * pi
     t = 16.0D+00 / ( x * x )

     f0 = (((((( &
          0.1496119D-02     * t &
          - 0.739083D-02 )    * t &
          + 0.016236617D+00 ) * t &
          - 0.022007499D+00 ) * t &
          + 0.023644978D+00 ) * t &
          - 0.031280848D+00 ) * t &
          + 0.124611058D+00 ) * 4.0D+00 / x

     g0 = ((((( &
          0.1076103D-02     * t &
          - 0.5434851D-02 )   * t &
          + 0.01242264D+00 )  * t &
          - 0.018255209D+00 ) * t &
          + 0.023664841D+00 ) * t &
          - 0.049635633D+00 ) * t &
          + 0.79784879D+00

     tj = 1.0D+00 - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

     ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

  else

     t = 64.0D+00 / ( x * x )
     xt = x-0.25D+00 * pi

     f0 = ((((((( &
          - 0.268482D-04     * t &
          + 0.1270039D-03 )  * t &
          - 0.2755037D-03 )  * t &
          + 0.3992825D-03 )  * t &
          - 0.5366169D-03 )  * t &
          + 0.10089872D-02 ) * t &
          - 0.40403539D-02 ) * t &
          + 0.0623347304D+00 ) * 8.0D+00 / x

     g0 = (((((( &
          - 0.226238D-04        * t &
          + 0.1107299D-03 )     * t &
          - 0.2543955D-03 )     * t &
          + 0.4100676D-03 )     * t &
          - 0.6740148D-03 )     * t &
          + 0.17870944D-02 )    * t &
          - 0.01256424405D+00 ) * t &
          + 0.79788456D+00

     tj = 1.0D+00  - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

     ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

  end if

  return
end subroutine itjyb




  !*****************************************************************************80
  !
  !! ITSH0 integrates the Struve function H0(t) from 0 to x.
  !
  !  Discussion:
  !
  !    This procedure evaluates the integral of Struve function
  !    H0(t) with respect to t from 0 and x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    25 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TH0, the integral of H0(t) from 0 to x.
  !
subroutine itsh0 ( x, th0 )
  !ITSH0 integrates the Struve function H0(t) from 0 to x.
  !
  !This procedure evaluates the integral of Struve function
  !H0(t) with respect to t from 0 and x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param th0: out: the integral of H0(t) from 0 to x
  !
  implicit none

  real(8)    :: a(25)
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: af
  real(8)    :: bf
  real(8)    :: bg
  real(8)    :: el
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: rd
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: th0
  real(8)    :: ty
  real(8)    :: x
  real(8)    :: xp

  pi = 3.141592653589793D+00
  r = 1.0D+00            

  if ( x <= 30.0D+00 ) then

     s = 0.5D+00

     do k = 1, 100

        if ( k == 1 ) then
           rd = 0.5D+00
        else
           rd = 1.0D+00
        end if

        r = - r * rd * k / ( k + 1.0D+00 ) &
             * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
        s = s + r

        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if

     end do

     th0 = 2.0D+00 / pi * x * x * s

  else

     s = 1.0D+00
     do k = 1, 12
        r = - r * k / ( k + 1.0D+00 ) &
             * ( ( 2.0D+00 * k + 1.0D+00 ) / x ) ** 2
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     el = 0.57721566490153D+00
     s0 = s / ( pi * x * x ) + 2.0D+00 / pi &
          * ( log ( 2.0D+00 * x ) + el )
     a0 = 1.0D+00
     a1 = 5.0D+00 / 8.0D+00
     a(1) = a1
     do k = 1, 20
        af = ( ( 1.5D+00 * ( k + 0.5D+00 ) &
             * ( k + 5.0D+00 / 6.0D+00 ) * a1 - 0.5D+00 &
             * ( k + 0.5D+00 ) * ( k + 0.5D+00 ) &
             * ( k - 0.5D+00 ) * a0 ) ) / ( k + 1.0D+00 )
        a(k+1) = af
        a0 = a1
        a1 = af
     end do

     bf = 1.0D+00
     r = 1.0D+00
     do k = 1, 10
        r = - r / ( x * x )
        bf = bf + a(2*k) * r
     end do
     bg = a(1) / x
     r = 1.0D+00 / x
     do k = 1, 10
        r = - r / ( x * x ) 
        bg = bg + a(2*k+1) * r
     end do
     xp = x + 0.25D+00 * pi
     ty = sqrt ( 2.0D+00 / ( pi * x ) ) &
          * ( bg * cos ( xp ) - bf * sin ( xp ) )
     th0 = ty + s0

  end if

  return
end subroutine itsh0




  !*****************************************************************************80
  !
  !! ITSL0 integrates the Struve function L0(t) from 0 to x.
  !
  !  Discussion:
  !
  !    This procedure evaluates the integral of modified Struve function
  !    L0(t) with respect to t from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the upper limit of the integral.
  !
  !    Output, real(8) TL0, the integral of L0(t) from 0 to x.
  !
subroutine itsl0 ( x, tl0 )
  !ITSL0 integrates the Struve function L0(t) from 0 to x.
  !
  !This procedure evaluates the integral of modified Struve function
  !L0(t) with respect to t from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the upper limit of the integral
  !@param tl0: out: the integral of L0(t) from 0 to x
  !
  implicit none

  real(8)    :: a(18)
  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: af
  real(8)    :: el
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: rd
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: ti
  real(8)    :: tl0
  real(8)    :: x

  pi = 3.141592653589793D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then

     s = 0.5D+00
     do k = 1, 100

        if ( k == 1 ) then
           rd = 0.5D+00
        else
           rd = 1.0D+00
        end if
        r = r * rd * k / ( k + 1.0D+00 ) &
             * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     tl0 = 2.0D+00 / pi * x * x * s

  else

     s = 1.0D+00
     do k = 1, 10
        r = r * k / ( k + 1.0D+00 ) &
             * ( ( 2.0D+00 * k + 1.0D+00 ) / x ) ** 2
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     el = 0.57721566490153D+00
     s0 = - s / ( pi * x * x ) + 2.0D+00 / pi &
          * ( log ( 2.0D+00 * x ) + el )
     a0 = 1.0D+00
     a1 = 5.0D+00 / 8.0D+00
     a(1) = a1
     do k = 1, 10
        af = ( ( 1.5D+00 * ( k + 0.50D+00 ) &
             * ( k + 5.0D+00 / 6.0D+00 ) * a1 - 0.5D+00 &
             * ( k + 0.5D+00 ) ** 2 * ( k -0.5D+00 ) * a0 ) ) &
             / ( k + 1.0D+00 )
        a(k+1) = af
        a0 = a1
        a1 = af
     end do

     ti = 1.0D+00
     r = 1.0D+00
     do k = 1, 11
        r = r / x
        ti = ti + a(k) * r
     end do
     tl0 = ti / sqrt ( 2.0D+00 * pi * x ) * exp ( x ) + s0

  end if

  return
end subroutine itsl0




  !*****************************************************************************80
  !
  !! ITTH0 integrates H0(t)/t from x to oo.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    23 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the lower limit of the integral.
  !
  !    Output, real(8) TTH, the integral of H0(t)/t from x to oo.
  !
subroutine itth0 ( x, tth )
  !ITTH0 integrates H0(t)/t from x to oo.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the lower limit of the integral
  !@param tth: out: the integral of H0(t)/t from x to oo
  !
  implicit none

  real(8)    :: f0
  real(8)    :: g0
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: s
  real(8)    :: t
  real(8)    :: tth
  real(8)    :: tty
  real(8)    :: x
  real(8)    :: xt

  pi = 3.141592653589793D+00
  s = 1.0D+00
  r = 1.0D+00

  if ( x < 24.5D+00 ) then

     do k = 1, 60
        r = - r * x * x * ( 2.0D+00 * k - 1.0D+00 ) &
             / ( 2.0D+00 * k + 1.0D+00 ) ** 3
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     tth = pi / 2.0D+00 - 2.0D+00 / pi * x * s

  else

     do k = 1, 10
        r = - r * ( 2.0D+00 * k - 1.0D+00 ) ** 3 &
             / ( ( 2.0D+00 * k + 1.0D+00 ) * x * x )
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     tth = 2.0D+00 / ( pi * x ) * s
     t = 8.0D+00 / x
     xt = x + 0.25D+00 * pi
     f0 = ((((( &
          0.18118D-02 * t &
          - 0.91909D-02 ) * t &
          + 0.017033D+00 ) * t &
          - 0.9394D-03 ) * t &
          - 0.051445D+00 ) * t &
          - 0.11D-05 ) * t &
          + 0.7978846D+00
     g0 = ((((( &
          - 0.23731D-02 * t &
          + 0.59842D-02 ) * t &
          + 0.24437D-02 ) * t &
          - 0.0233178D+00 ) * t &
          + 0.595D-04 ) * t &
          + 0.1620695D+00 ) * t
     tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )
     tth = tth + tty

  end if

  return
end subroutine itth0




  !*****************************************************************************80
  !
  !! ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    23 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the integral limit.
  !
  !    Output, real(8) TTI, TTK, the integrals of [I0(t)-1]/t 
  !    from 0 to x, and of K0(t)/t from x to oo.
  !
subroutine ittika ( x, tti, ttk )
  !ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the integral limit
  !@param tti: out: the integral of [I0(t)-1]/t from 0 to x
  !@param ttk: out: the integral of K0(t)/t from x to oo
  !
  implicit none

  real(8)    :: b1
  real(8), save, dimension ( 8 ) :: c = (/ &
       1.625D+00, 4.1328125D+00, &
       1.45380859375D+01, 6.553353881835D+01, &
       3.6066157150269D+02, 2.3448727161884D+03, &
       1.7588273098916D+04, 1.4950639538279D+05 /)
  real(8)    :: e0
  real(8)    :: el
  integer(4) :: k
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: r2
  real(8)    :: rc
  real(8)    :: rs
  real(8)    :: tti
  real(8)    :: ttk
  real(8)    :: x

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( x == 0.0D+00 ) then
     tti = 0.0D+00
     ttk = 1.0D+300
     return
  end if

  if ( x < 40.0D+00 ) then
     tti = 1.0D+00
     r = 1.0D+00
     do k = 2, 50
        r = 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
        tti = tti + r
        if ( abs ( r / tti ) < 1.0D-12 ) then
           exit
        end if
     end do

     tti = tti * 0.125D+00 * x * x

  else

     tti = 1.0D+00
     r = 1.0D+00
     do k = 1, 8
        r = r / x
        tti = tti + c(k) * r
     end do
     rc = x * sqrt ( 2.0D+00 * pi * x )
     tti = tti * exp ( x ) / rc

  end if

  if ( x <= 12.0D+00 ) then

     e0 = ( 0.5D+00 * log ( x / 2.0D+00 ) + el ) &
          * log ( x / 2.0D+00 ) + pi * pi / 24.0D+00 + 0.5D+00 * el * el
     b1 = 1.5D+00 - ( el + log ( x / 2.0D+00 ) )
     rs = 1.0D+00
     r = 1.0D+00
     do k = 2, 50
        r = 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
        rs = rs + 1.0D+00 / k
        r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k ) &
             - ( el + log ( x / 2.0D+00 ) ) )
        b1 = b1 + r2
        if ( abs ( r2 / b1 ) < 1.0D-12 ) then
           exit
        end if
     end do

     ttk = e0 - 0.125D+00 * x * x * b1

  else

     ttk = 1.0D+00
     r = 1.0D+00
     do k = 1, 8
        r = - r / x
        ttk = ttk + c(k) * r
     end do
     rc = x * sqrt ( 2.0D+00 / pi * x )
     ttk = ttk * exp ( - x ) / rc

  end if

  return
end subroutine ittika




  !*****************************************************************************80
  !
  !! ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the integral limit.
  !
  !    Output, real(8) TTI, TTK, the integrals of
  !    [I0(t)-1]/t from 0 to x, and K0(t)/t from x to oo.
  !
subroutine ittikb ( x, tti, ttk )
  !ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the integral limit
  !@param tti: out: the integral of [I0(t)-1]/t from 0 to x
  !@param ttk: out: the integral of K0(t)/t from x to oo
  !
  implicit none

  real(8) :: e0
  real(8) :: el
  real(8) :: pi
  real(8) :: t
  real(8) :: t1
  real(8) :: tti
  real(8) :: ttk
  real(8) :: x
  real(8) :: x1

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( x == 0.0D+00 ) then

     tti = 0.0D+00

  else if ( x <= 5.0D+00 ) then

     x1 = x / 5.0D+00
     t = x1 * x1
     tti = ((((((( &
          0.1263D-03       * t &
          + 0.96442D-03 )    * t &
          + 0.968217D-02 )   * t &
          + 0.06615507D+00 ) * t &
          + 0.33116853D+00 ) * t &
          + 1.13027241D+00 ) * t &
          + 2.44140746D+00 ) * t &
          + 3.12499991D+00 ) * t

  else

     t = 5.0D+00 / x
     tti = ((((((((( &
          2.1945464D+00   * t &
          -  3.5195009D+00 ) * t &
          - 11.9094395D+00 ) * t &
          + 40.394734D+00  ) * t &
          - 48.0524115D+00 ) * t &
          + 28.1221478D+00 ) * t &
          -  8.6556013D+00 ) * t &
          +  1.4780044D+00 ) * t &
          -  0.0493843D+00 ) * t &
          +  0.1332055D+00 ) * t &
          +  0.3989314D+00
     tti = tti * exp ( x ) / ( sqrt ( x ) * x )

  end if

  if ( x == 0.0D+00 ) then

     ttk = 1.0D+300

  else if ( x <= 2.0D+00 ) then

     t1 = x / 2.0D+00
     t = t1 * t1
     ttk = ((((( &
          0.77D-06         * t &
          + 0.1544D-04 )     * t &
          + 0.48077D-03 )    * t &
          + 0.925821D-02 )   * t &
          + 0.10937537D+00 ) * t &
          + 0.74999993D+00 ) * t
     e0 = el + log ( x / 2.0D+00 )
     ttk = pi * pi / 24.0D+00 + e0 * ( 0.5D+00 * e0 + tti ) - ttk

  else if ( x <= 4.0D+00 ) then

     t = 2.0D+00 / x
     ttk = ((( &
          0.06084D+00    * t &
          - 0.280367D+00 ) * t &
          + 0.590944D+00 ) * t &
          - 0.850013D+00 ) * t &
          + 1.234684D+00
     ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

  else

     t = 4.0D+00 / x
     ttk = ((((( &
          0.02724D+00     * t &
          - 0.1110396D+00 ) * t &
          + 0.2060126D+00 ) * t &
          - 0.2621446D+00 ) * t &
          + 0.3219184D+00 ) * t &
          - 0.5091339D+00 ) * t &
          + 1.2533141D+00
     ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

  end if

  return
end subroutine ittikb




  !*****************************************************************************80
  !
  !! ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the integral limit.
  !
  !    Output, real(8) TTJ, TTY, the integrals of [1-J0(t)]/t 
  !    from 0 to x and of Y0(t)/t from x to oo.
  !
subroutine ittjya ( x, ttj, tty )
  !ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the integral limit
  !@param ttj: out: the integral of [1-J0(t)]/t from 0 to x
  !@param tty: out: the integral of Y0(t)/t from x to oo
  !
  implicit none

  real(8)    :: a0
  real(8)    :: b1
  real(8)    :: bj0
  real(8)    :: bj1
  real(8)    :: by0
  real(8)    :: by1
  real(8)    :: e0
  real(8)    :: el
  real(8)    :: g0
  real(8)    :: g1
  integer(4) :: k
  integer(4) :: l
  real(8)    :: pi
  real(8)    :: px
  real(8)    :: qx
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: rs
  real(8)    :: t
  real(8)    :: ttj
  real(8)    :: tty
  real(8)    :: vt
  real(8)    :: x
  real(8)    :: xk

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( x == 0.0D+00 ) then

     ttj = 0.0D+00
     tty = -1.0D+300

  else if ( x <= 20.0D+00 ) then

     ttj = 1.0D+00
     r = 1.0D+00
     do k = 2, 100
        r = - 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
        ttj = ttj + r
        if ( abs ( r ) < abs ( ttj ) * 1.0D-12 ) then
           exit
        end if
     end do

     ttj = ttj * 0.125D+00 * x * x
     e0 = 0.5D+00 * ( pi * pi / 6.0D+00 - el * el ) &
          - ( 0.5D+00 * log ( x / 2.0D+00 ) + el ) &
          * log ( x / 2.0D+00 )
     b1 = el + log ( x / 2.0D+00 ) - 1.5D+00
     rs = 1.0D+00
     r = -1.0D+00
     do k = 2, 100
        r = - 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
        rs = rs + 1.0D+00 / k
        r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k ) &
             - ( el + log ( x / 2.0D+00 ) ) ) 
        b1 = b1 + r2
        if ( abs ( r2 ) < abs ( b1 ) * 1.0D-12 ) then
           exit
        end if
     end do

     tty = 2.0D+00 / pi * ( e0 + 0.125D+00 * x * x * b1 )

  else

     a0 = sqrt ( 2.0D+00 / ( pi * x ) )

     do l = 0, 1

        vt = 4.0D+00 * l * l
        px = 1.0D+00
        r = 1.0D+00
        do k = 1, 14
           r = - 0.0078125D+00 * r &
                * ( vt - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
                / ( x * k ) * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                / ( ( 2.0D+00 * k - 1.0D+00 ) * x )
           px = px + r
           if ( abs ( r ) < abs ( px ) * 1.0D-12 ) then
              exit
           end if
        end do

        qx = 1.0D+00
        r = 1.0D+00
        do k = 1, 14
           r = -0.0078125D+00 * r &
                * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                / ( x * k ) * ( vt - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
                / ( 2.0D+00 * k + 1.0D+00 ) / x
           qx = qx + r
           if ( abs ( r ) < abs ( qx ) * 1.0D-12 ) then
              exit
           end if
        end do

        qx = 0.125D+00 * ( vt - 1.0D+00 ) / x * qx
        xk = x - ( 0.25D+00 + 0.5D+00 * l ) * pi
        bj1 = a0 * ( px * cos ( xk ) - qx * sin ( xk ) )
        by1 = a0 * ( px * sin ( xk ) + qx * cos ( xk ) )
        if ( l == 0 ) then
           bj0 = bj1
           by0 = by1
        end if

     end do

     t = 2.0D+00 / x
     g0 = 1.0D+00
     r0 = 1.0D+00
     do k = 1, 10
        r0 = - k * k * t * t *r0
        g0 = g0 + r0
     end do

     g1 = 1.0D+00
     r1 = 1.0D+00
     do k = 1, 10
        r1 = - k * ( k + 1.0D+00 ) * t * t * r1
        g1 = g1 + r1
     end do

     ttj = 2.0D+00 * g1 * bj0 / ( x * x ) - g0 * bj1 / x &
          + el + log ( x / 2.0D+00 )
     tty = 2.0D+00 * g1 * by0 / ( x * x ) - g0 * by1 / x

  end if

  return
end subroutine ittjya




  !*****************************************************************************80
  !
  !! ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the integral limit.
  !
  !    Output, real(8) TTJ, TTY, the integrals of [1-J0(t)]/t 
  !    from 0 to x and of Y0(t)/t from x to oo.
  !
subroutine ittjyb ( x, ttj, tty )
  !ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the integral limit
  !@param ttj: out: the integral of [1-J0(t)]/t from 0 to x
  !@param tty: out: the integral of Y0(t)/t from x to oo
  !
  implicit none

  real(8) :: e0
  real(8) :: el
  real(8) :: f0
  real(8) :: g0
  real(8) :: pi
  real(8) :: t
  real(8) :: t1
  real(8) :: ttj
  real(8) :: tty
  real(8) :: x
  real(8) :: x1
  real(8) :: xt

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00

  if ( x == 0.0D+00 ) then

     ttj = 0.0D+00
     tty = -1.0D+300

  else if ( x <= 4.0D+00 ) then

     x1 = x / 4.0D+00
     t = x1 * x1

     ttj = (((((( &
          0.35817D-04 * t &
          - 0.639765D-03 ) * t &
          + 0.7092535D-02 ) * t &
          - 0.055544803D+00 ) * t &
          + 0.296292677D+00 ) * t &
          - 0.999999326D+00 ) * t &
          + 1.999999936D+00 ) * t

     tty = ((((((( &
          - 0.3546D-05        * t &
          + 0.76217D-04 )     * t &
          - 0.1059499D-02 )   * t &
          + 0.010787555D+00 ) * t &
          - 0.07810271D+00 )  * t &
          + 0.377255736D+00 ) * t &
          - 1.114084491D+00 ) * t &
          + 1.909859297D+00 ) * t

     e0 = el + log ( x / 2.0D+00 )
     tty = pi / 6.0D+00 + e0 / pi * ( 2.0D+00 * ttj - e0 ) - tty

  else if ( x <= 8.0D+00 ) then

     xt = x + 0.25D+00 * pi
     t1 = 4.0D+00 / x
     t = t1 * t1

     f0 = ((((( &
          0.0145369D+00 * t &
          - 0.0666297D+00 ) * t &
          + 0.1341551D+00 ) * t &
          - 0.1647797D+00 ) * t &
          + 0.1608874D+00 ) * t &
          - 0.2021547D+00 ) * t &
          + 0.7977506D+00

     g0 = (((((( &
          0.0160672D+00   * t &
          - 0.0759339D+00 ) * t &
          + 0.1576116D+00 ) * t &
          - 0.1960154D+00 ) * t &
          + 0.1797457D+00 ) * t &
          - 0.1702778D+00 ) * t &
          + 0.3235819D+00 ) * t1

     ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) ) / ( sqrt ( x ) * x )
     ttj = ttj + el + log ( x / 2.0D+00 )
     tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )

  else

     t = 8.0D+00 / x
     xt = x + 0.25D+00 * pi

     f0 = ((((( &
          0.18118D-02    * t &
          - 0.91909D-02 )  * t &
          + 0.017033D+00 ) * t &
          - 0.9394D-03 )   * t &
          - 0.051445D+00 ) * t &
          - 0.11D-05 )     * t &
          + 0.7978846D+00

     g0 = ((((( &
          - 0.23731D-02     * t &
          + 0.59842D-02 )   * t &
          + 0.24437D-02 )   * t &
          - 0.0233178D+00 ) * t &
          + 0.595D-04 )     * t &
          + 0.1620695D+00 ) * t

     ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) )  &
          / ( sqrt ( x ) * x ) + el + log ( x / 2.0D+00 )
     tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) )  &
          / ( sqrt ( x ) * x )

  end if

  return
end subroutine ittjyb




  !*****************************************************************************80
  !
  !! JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
  !
  !  Discussion:
  !
  !    This procedure computes the zeros of Bessel functions Jn(x) and
  !    Jn'(x), and arrange them in the order of their magnitudes.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Output, integer(4) N(*), the  order of Jn(x) or Jn'(x) associated
  !    with the L-th zero.
  !
  !    Output, integer(4) M(*), the serial number of the zeros of Jn(x)
  !    or Jn'(x) associated with the L-th zero ( L is the serial number of all the
  !    zeros of Jn(x) and Jn'(x) ).
  !
  !    Output, character ( len = 4 ) P(L), 'TM' or 'TE', a code for designating 
  !    the zeros of Jn(x)  or Jn'(x).  In the waveguide applications, the zeros
  !    of Jn(x) correspond to TM modes and those of Jn'(x) correspond to TE modes.
  !
  !    Output, real(8) ZO(*), the zeros of Jn(x) and Jn'(x).
  !
subroutine jdzo ( nt, n, m, p, zo )
  !JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
  !
  !This procedure computes the zeros of Bessel functions Jn(x) and
  !Jn'(x), and arrange them in the order of their magnitudes.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* M: the serial number of the zeros of Jn(x) or Jn'(x) associated with the L-th zero ( L is the serial number of all
  !  the zeros of Jn(x) and Jn'(x) )
  !* P: 'TM' or 'TE', a code for designating the zeros of Jn(x) or Jn'(x). In the waveguide applications, the zeros of
  !  Jn(x) correspond to TM modes and those of Jn'(x) correspond to TE modes
  !
  !@param nt: in: the number of zeros
  !@param n: out: the order of Jn(x) or Jn'(x) associated with the L-th zero
  !@param m: out: serial number of the zeros of Jn(x) or Jn'(x) associated with the L-th zero
  !@param p: out: code, 'TM' for the zeros of Jn(x) or 'TE' for those of Jn'(x)
  !@param zo: out: the zeros of Jn(x) and Jn'(x)
  !
  implicit none

  real(8)               :: bj(101)
  real(8)               :: dj(101)
  real(8)               :: fj(101)
  integer(4)            :: i
  integer(4)            :: j
  integer(4)            :: k
  integer(4)            :: l
  integer(4)            :: l0
  integer(4)            :: l1
  integer(4)            :: l2
  integer(4)            :: m(1400)
  integer(4)            :: m1(70)
  integer(4)            :: mm
  integer(4)            :: n(1400)
  integer(4)            :: n1(70)
  integer(4)            :: nm
  integer(4)            :: nt
  character ( len = 4 ) :: p(1400)
  character ( len = 4 ) :: p1(70)
  real(8)               :: x
  real(8)               :: x0
  real(8)               :: x1
  real(8)               :: x2
  real(8)               :: xm
  real(8)               :: zo(1400)
  real(8)               :: zoc(70)

  if ( nt < 600 ) then
     xm = -1.0D+00 + 2.248485D+00 * real ( nt, kind = 8 ) ** 0.5D+00 &
          - 0.0159382D+00 * nt + 3.208775D-04 * real ( nt, kind = 8 ) ** 1.5D+00
     nm = int ( 14.5D+00 + 0.05875D+00 * nt )
     mm = int ( 0.02D+00 * nt ) + 6
  else
     xm = 5.0D+00 + 1.445389D+00 * ( real ( nt, kind = 8 ) ) ** 0.5D+00 &
          + 0.01889876D+00 * nt &
          - 2.147763D-04 * ( real ( nt, kind = 8 ) ) ** 1.5D+00
     nm = int ( 27.8D+00 + 0.0327D+00 * nt )
     mm = int ( 0.01088D+00 * nt ) + 10
  end if

  l0 = 0

  do i = 1,nm

     x1 = 0.407658D+00 + 0.4795504D+00 &
          * ( real ( i - 1, kind = 8 ) ) ** 0.5D+00 + 0.983618D+00 * ( i - 1 )
     x2 = 1.99535D+00 + 0.8333883 * ( real ( i - 1, kind = 8 ) ) ** 0.5D+00 &
          + 0.984584D+00 * ( i - 1 )
     l1 = 0

     do j = 1, mm

        if ( i == 1 .and. j == 1 ) then

           l1 = l1 + 1
           n1(l1) = i - 1
           m1(l1) = j
           if ( i == 1 ) then
              m1(l1) = j - 1
           end if
           p1(l1) = 'TE'
           zoc(l1) = x

           if ( i <= 15 ) then
              x1 = x + 3.057D+00 + 0.0122D+00 * ( i - 1 ) &
                   + ( 1.555D+00 + 0.41575D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
           else
              x1 = x + 2.918D+00 + 0.01924D+00 * ( i - 1 ) &
                   + ( 6.26D+00 + 0.13205D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
           end if

        else

           x = x1

           do

              call bjndd ( i, x, bj, dj, fj )
              x0 = x
              x = x - dj(i) / fj(i)

              if ( xm < x1 ) then
                 exit
              end if

              if ( abs ( x - x0 ) <= 1.0D-10 ) then
                 l1 = l1 + 1
                 n1(l1) = i - 1
                 m1(l1) = j
                 if ( i == 1 ) then
                    m1(l1) = j - 1
                 end if
                 p1(l1) = 'TE'
                 zoc(l1) = x

                 if ( i <= 15 ) then
                    x1 = x + 3.057D+00 + 0.0122D+00 * ( i - 1 ) &
                         + ( 1.555D+00 + 0.41575D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
                 else
                    x1 = x + 2.918D+00 + 0.01924D+00 * ( i - 1 ) &
                         + ( 6.26D+00 + 0.13205D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
                 end if
                 exit
              end if

           end do

        end if

        x = x2

        do

           call bjndd ( i, x, bj, dj, fj )
           x0 = x
           x = x - bj(i) / dj(i)

           if ( xm < x ) then
              exit
           end if

           if ( abs ( x - x0 ) <= 1.0D-10 ) then
              exit
           end if

        end do

        if ( x <= xm ) then

           l1 = l1 + 1
           n1(l1) = i - 1
           m1(l1) = j
           p1(l1) = 'TM'
           zoc(l1) = x
           if ( i <= 15 ) then
              x2 = x + 3.11D+00 + 0.0138D+00 * ( i - 1 ) &
                   + ( 0.04832D+00 + 0.2804D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
           else
              x2 = x + 3.001D+00 + 0.0105D+00 * ( i - 1 ) &
                   + ( 11.52D+00 + 0.48525D+00 * ( i - 1 ) ) / ( j + 3 ) ** 2
           end if

        end if

     end do

     l = l0 + l1
     l2 = l

     do

        if ( l0 == 0 ) then
           do k = 1, l
              zo(k) = zoc(k)
              n(k) = n1(k)
              m(k) = m1(k)
              p(k) = p1(k)
           end do
           l1 = 0
        else if ( l0 /= 0 ) then
           if ( zoc(l1) .le. zo(l0) ) then
              zo(l0+l1) = zo(l0)
              n(l0+l1) = n(l0)
              m(l0+l1) = m(l0)
              p(l0+l1) = p(l0)
              l0 = l0 - 1
           else
              zo(l0+l1) = zoc(l1)
              n(l0+l1) = n1(l1)
              m(l0+l1) = m1(l1)
              p(l0+l1) = p1(l1)
              l1 = l1 - 1
           end if
        end if

        if ( l1 == 0 ) then
           exit 
        end if

     end do

     l0 = l2

  end do

  return
end subroutine jdzo




  !*****************************************************************************80
  !
  !! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    08 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) U, the argument.
  !
  !    Input, real(8) HK, the modulus, between 0 and 1.
  !
  !    Output, real(8) ESN, ECN, EDN, EPH, the values of
  !    sn(u), cn(u), dn(u), and phi (in degrees).
  !
subroutine jelp ( u, hk, esn, ecn, edn, eph )
  !JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param u: in: the argument
  !@param hk: in: the modulus, between 0 and 1
  !@param esn: out: sn(u)
  !@param ecn: out: cn(u)
  !@param edn: out: dn(u)
  !@param eph: out: and phi (in degrees)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: b
  real(8)    :: b0
  real(8)    :: c
  real(8)    :: d
  real(8)    :: dn
  real(8)    :: ecn
  real(8)    :: edn
  real(8)    :: eph
  real(8)    :: esn
  real(8)    :: hk
  integer(4) :: j
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: r(40)
  real(8)    :: sa
  real(8)    :: t
  real(8)    :: u

  pi = 3.14159265358979D+00
  a0 = 1.0D+00
  b0 = sqrt ( 1.0D+00 - hk * hk )

  do n = 1, 40

     a = ( a0 + b0 ) / 2.0D+00
     b = sqrt ( a0 * b0 )
     c = ( a0 - b0 ) / 2.0D+00
     r(n) = c / a

     if ( c < 1.0D-07 ) then
        exit
     end if

     a0 = a
     b0 = b

  end do

  dn = 2.0D+00 ** n * a * u

  do j = n, 1, -1
     t = r(j) * sin ( dn )
     sa = atan ( t / sqrt ( abs ( 1.0D+00 - t * t )))
     d = 0.5D+00 * ( dn + sa )
     dn = d
  end do

  eph = d * 180.0D+00 / pi
  esn = sin ( d )
  ecn = cos ( d )
  edn = sqrt ( 1.0D+00 - hk * hk * esn * esn )

  return
end subroutine jelp




  !*****************************************************************************80
  !
  !! JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    01 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
  !    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
  !
subroutine jy01a ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
  !JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param bj0: out: J0(x)
  !@param dj0: out: J0'(x)
  !@param bj1: out: J1(x)
  !@param dj1: out: J1'(x)
  !@param by0: out: Y0(x)
  !@param dy0: out: Y0'(x)
  !@param by1: out: Y1(x)
  !@param dy1: out: Y1'(x)
  !
  implicit none

  real(8), save, dimension(12) :: a = (/ &
       -0.7031250000000000D-01, 0.1121520996093750D+00, &
       -0.5725014209747314D+00, 0.6074042001273483D+01, &
       -0.1100171402692467D+03, 0.3038090510922384D+04, &
       -0.1188384262567832D+06, 0.6252951493434797D+07, &
       -0.4259392165047669D+09, 0.3646840080706556D+11, &
       -0.3833534661393944D+13, 0.4854014686852901D+15 /)
  real(8), save, dimension(12) :: a1 = (/ &
       0.1171875000000000D+00, -0.1441955566406250D+00, &
       0.6765925884246826D+00, -0.6883914268109947D+01, &
       0.1215978918765359D+03, -0.3302272294480852D+04, &
       0.1276412726461746D+06, -0.6656367718817688D+07, &
       0.4502786003050393D+09, -0.3833857520742790D+11, &
       0.4011838599133198D+13, -0.5060568503314727D+15 /)
  real(8), save, dimension(12) :: b = (/ &
       0.7324218750000000D-01, -0.2271080017089844D+00, &
       0.1727727502584457D+01, -0.2438052969955606D+02, &
       0.5513358961220206D+03, -0.1825775547429318D+05, &
       0.8328593040162893D+06, -0.5006958953198893D+08, &
       0.3836255180230433D+10, -0.3649010818849833D+12, &
       0.4218971570284096D+14, -0.5827244631566907D+16 /)
  real(8), save, dimension(12) :: b1 = (/ &
       -0.1025390625000000D+00, 0.2775764465332031D+00, &
       -0.1993531733751297D+01, 0.2724882731126854D+02, &
       -0.6038440767050702D+03, 0.1971837591223663D+05, &
       -0.8902978767070678D+06, 0.5310411010968522D+08, &
       -0.4043620325107754D+10, 0.3827011346598605D+12, &
       -0.4406481417852278D+14, 0.6065091351222699D+16 /)
  real(8)    :: bj0
  real(8)    :: bj1
  real(8)    :: by0
  real(8)    :: by1
  real(8)    :: cs0
  real(8)    :: cs1
  real(8)    :: cu
  real(8)    :: dj0
  real(8)    :: dj1
  real(8)    :: dy0
  real(8)    :: dy1
  real(8)    :: ec
  integer(4) :: k
  integer(4) :: k0
  real(8)    :: p0
  real(8)    :: p1
  real(8)    :: pi
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: rp2
  real(8)    :: t1
  real(8)    :: t2
  real(8)    :: w0
  real(8)    :: w1
  real(8)    :: x
  real(8)    :: x2

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  x2 = x * x

  if ( x == 0.0D+00 ) then
     bj0 = 1.0D+00
     bj1 = 0.0D+00
     dj0 = 0.0D+00
     dj1 = 0.5D+00
     by0 = -1.0D+300
     by1 = -1.0D+300
     dy0 = 1.0D+300
     dy1 = 1.0D+300
     return
  end if

  if ( x <= 12.0D+00 ) then

     bj0 = 1.0D+00
     r = 1.0D+00
     do k = 1,30
        r = -0.25D+00 * r * x2 / ( k * k )
        bj0 = bj0 + r
        if ( abs ( r ) < abs ( bj0 ) * 1.0D-15 ) then
           exit
        end if
     end do

     bj1 = 1.0D+00
     r = 1.0D+00
     do k = 1, 30
        r = -0.25D+00 * r * x2 / ( k * ( k + 1.0D+00 ) )
        bj1 = bj1 + r
        if ( abs ( r ) < abs ( bj1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     bj1 = 0.5D+00 * x * bj1
     ec = log ( x / 2.0D+00 ) + 0.5772156649015329D+00
     cs0 = 0.0D+00
     w0 = 0.0D+00
     r0 = 1.0D+00
     do k = 1, 30
        w0 = w0 + 1.0D+00 / k
        r0 = -0.25D+00 * r0 / ( k * k ) * x2
        r = r0 * w0
        cs0 = cs0 + r
        if ( abs ( r ) < abs ( cs0 ) * 1.0D-15 ) then
           exit
        end if
     end do

     by0 = rp2 * ( ec * bj0 - cs0 )
     cs1 = 1.0D+00
     w1 = 0.0D+00
     r1 = 1.0D+00
     do k = 1, 30
        w1 = w1 + 1.0D+00 / k
        r1 = -0.25D+00 * r1 / ( k * ( k + 1 ) ) * x2
        r = r1 * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
        cs1 = cs1 + r
        if ( abs ( r ) < abs ( cs1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     by1 = rp2 * ( ec * bj1 - 1.0D+00 / x - 0.25D+00 * x * cs1 )

  else

     if ( x < 35.0D+00 ) then
        k0 = 12
     else if ( x < 50.0D+00 ) then
        k0 = 10
     else
        k0 = 8
     end if

     t1 = x - 0.25D+00 * pi
     p0 = 1.0D+00
     q0 = -0.125D+00 / x
     do k = 1, k0
        p0 = p0 + a(k) * x ** ( - 2 * k )
        q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
     end do
     cu = sqrt ( rp2 / x )
     bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
     by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
     t2 = x - 0.75D+00 * pi
     p1 = 1.0D+00
     q1 = 0.375D+00 / x
     do k = 1, k0
        p1 = p1 + a1(k) * x ** ( - 2 * k )
        q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
     end do
     cu = sqrt ( rp2 / x )
     bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
     by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )

  end if

  dj0 = - bj1
  dj1 = bj0 - bj1 / x
  dy0 = - by1
  dy1 = by0 - by1 / x

  return
end subroutine jy01a




  !*****************************************************************************80
  !
  !! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
  !    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
  !
subroutine jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
  !JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param bj0: out: J0(x)
  !@param dj0: out: J0'(x)
  !@param bj1: out: J1(x)
  !@param dj1: out: J1'(x)
  !@param by0: out: Y0(x)
  !@param dy0: out: Y0'(x)
  !@param by1: out: Y1(x)
  !@param dy1: out: Y1'(x)
  !
  implicit none

  real(8) :: a0
  real(8) :: bj0
  real(8) :: bj1
  real(8) :: by0
  real(8) :: by1
  real(8) :: dj0
  real(8) :: dj1
  real(8) :: dy0
  real(8) :: dy1
  real(8) :: p0
  real(8) :: p1
  real(8) :: pi
  real(8) :: q0
  real(8) :: q1
  real(8) :: t
  real(8) :: t2
  real(8) :: ta0
  real(8) :: ta1
  real(8) :: x

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

     bj0 = 1.0D+00
     bj1 = 0.0D+00
     dj0 = 0.0D+00
     dj1 = 0.5D+00
     by0 = -1.0D+300
     by1 = -1.0D+300
     dy0 = 1.0D+300
     dy1 = 1.0D+300
     return

  else if ( x <= 4.0D+00 ) then

     t = x / 4.0D+00
     t2 = t * t

     bj0 = (((((( &
          - 0.5014415D-03 * t2 &
          + 0.76771853D-02 ) * t2 &
          - 0.0709253492D+00 ) * t2 &
          + 0.4443584263D+00 ) * t2 &
          - 1.7777560599D+00 ) * t2 &
          + 3.9999973021D+00 ) * t2 &
          - 3.9999998721D+00 ) * t2 &
          + 1.0D+00

     bj1 = t * ((((((( &
          - 0.1289769D-03 * t2 &
          + 0.22069155D-02 ) * t2 &
          - 0.0236616773D+00 ) * t2 &
          + 0.1777582922D+00 ) * t2 &
          - 0.8888839649D+00 ) * t2 &
          + 2.6666660544D+00 ) * t2 &
          - 3.9999999710D+00 ) * t2 &
          + 1.9999999998D+00 )

     by0 = ((((((( &
          - 0.567433D-04 * t2 &
          + 0.859977D-03 ) * t2 &
          - 0.94855882D-02 ) * t2 &
          + 0.0772975809D+00 ) * t2 &
          - 0.4261737419D+00 ) * t2 &
          + 1.4216421221D+00 ) * t2 &
          - 2.3498519931D+00 ) * t2 &
          + 1.0766115157D+00 ) * t2 &
          + 0.3674669052D+00

     by0 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj0 + by0

     by1 = (((((((( &
          0.6535773D-03 * t2 &
          - 0.0108175626D+00 ) * t2 &
          + 0.107657606D+00 ) * t2 &
          - 0.7268945577D+00 ) * t2 &
          + 3.1261399273D+00 ) * t2 &
          - 7.3980241381D+00 ) * t2 &
          + 6.8529236342D+00 ) * t2 &
          + 0.3932562018D+00 ) * t2 &
          - 0.6366197726D+00 ) / x

     by1 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj1 + by1

  else

     t = 4.0D+00 / x
     t2 = t * t
     a0 = sqrt ( 2.0D+00 / ( pi * x ) )

     p0 = (((( &
          - 0.9285D-05 * t2 &
          + 0.43506D-04 ) * t2 &
          - 0.122226D-03 ) * t2 &
          + 0.434725D-03 ) * t2 &
          - 0.4394275D-02 ) * t2 &
          + 0.999999997D+00

     q0 = t * ((((( &
          0.8099D-05 * t2 &
          - 0.35614D-04 ) * t2 &
          + 0.85844D-04 ) * t2 &
          - 0.218024D-03 ) * t2 &
          + 0.1144106D-02 ) * t2 &
          - 0.031249995D+00 )

     ta0 = x - 0.25D+00 * pi
     bj0 = a0 * ( p0 * cos ( ta0 ) - q0 * sin ( ta0 ) )
     by0 = a0 * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )

     p1 = (((( &
          0.10632D-04 * t2 &
          - 0.50363D-04 ) * t2 &
          + 0.145575D-03 ) * t2 &
          - 0.559487D-03 ) * t2 &
          + 0.7323931D-02 ) * t2 &
          + 1.000000004D+00

     q1 = t * ((((( &
          - 0.9173D-05      * t2 &
          + 0.40658D-04 )   * t2 &
          - 0.99941D-04 )   * t2 &
          + 0.266891D-03 )  * t2 &
          - 0.1601836D-02 ) * t2 &
          + 0.093749994D+00 )

     ta1 = x - 0.75D+00 * pi
     bj1 = a0 * ( p1 * cos ( ta1 ) - q1 * sin ( ta1 ) )
     by1 = a0 * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )

  end if

  dj0 = - bj1
  dj1 = bj0 - bj1 / x
  dy0 = - by1
  dy1 = by0 - by1 / x

  return
end subroutine jy01b




  !*****************************************************************************80
  !
  !! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 April 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
  !    of Jn(x), Jn'(x), Yn(x), Yn'(x).
  !
subroutine jyna ( n, x, nm, bj, dj, by, dy )
  !JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param bj: out: Jn(x)
  !@param dj: out: Jn'(x)
  !@param by: out: Yn(x)
  !@param dy: out: Yn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: bj(0:n)
  real(8)    :: bj0
  real(8)    :: bj1
  real(8)    :: bjk
  real(8)    :: by(0:n)
  real(8)    :: by0
  real(8)    :: by1
  real(8)    :: cs
  real(8)    :: dj(0:n)
  real(8)    :: dj0
  real(8)    :: dj1
  real(8)    :: dy(0:n)
  real(8)    :: dy0
  real(8)    :: dy1
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: f2
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: x

  nm = n

  if ( x < 1.0D-100 ) then

     do k = 0, n
        bj(k) = 0.0D+00
        dj(k) = 0.0D+00
        by(k) = -1.0D+300
        dy(k) = 1.0D+300
     end do
     bj(0) = 1.0D+00
     dj(1) = 0.5D+00
     return

  end if

  call jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
  bj(0) = bj0
  bj(1) = bj1
  by(0) = by0
  by(1) = by1
  dj(0) = dj0
  dj(1) = dj1
  dy(0) = dy0
  dy(1) = dy1

  if ( n <= 1 ) then
     return
  end if

  if ( n < int ( 0.9D+00 * x) ) then

     do k = 2, n
        bjk = 2.0D+00 * ( k - 1.0D+00 ) / x * bj1 - bj0
        bj(k) = bjk
        bj0 = bj1
        bj1 = bjk
     end do

  else

     m = msta1 ( x, 200 )

     if ( m < n ) then
        nm = m
     else
        m = msta2 ( x, n, 15 )
     end if

     f2 = 0.0D+00
     f1 = 1.0D-100
     do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 - f2
        if ( k <= nm ) then
           bj(k) = f
        end if
        f2 = f1
        f1 = f
     end do

     if ( abs ( bj1 ) < abs ( bj0 ) ) then
        cs = bj0 / f
     else
        cs = bj1 / f2
     end if

     do k = 0, nm
        bj(k) = cs * bj(k)
     end do

  end if

  do k = 2, nm
     dj(k) = bj(k-1) - k / x * bj(k)
  end do

  f0 = by(0)
  f1 = by(1)
  do k = 2, nm
     f = 2.0D+00 * ( k - 1.0D+00 ) / x * f1 - f0
     by(k) = f
     f0 = f1
     f1 = f
  end do

  do k = 2, nm
     dy(k) = by(k-1) - k * by(k) / x
  end do

  return
end subroutine jyna




  !*****************************************************************************80
  !
  !! JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
  !    of Jn(x), Jn'(x), Yn(x), Yn'(x).
  !
subroutine jynb ( n, x, nm, bj, dj, by, dy )
  !JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param bj: out: Jn(x)
  !@param dj: out: Jn'(x)
  !@param by: out: Yn(x)
  !@param dy: out: Yn'(x)
  !
  implicit none

  integer(4) :: n

  real(8), save, dimension ( 4 ) :: a = (/ &
       -0.7031250000000000D-01, 0.1121520996093750D+00, &
       -0.5725014209747314D+00, 0.6074042001273483D+01 /)
  real(8), save, dimension ( 4 ) :: a1 = (/ &
       0.1171875000000000D+00, -0.1441955566406250D+00, &
       0.6765925884246826D+00, -0.6883914268109947D+01 /)
  real(8), save, dimension ( 4 ) :: b = (/ &
       0.7324218750000000D-01, -0.2271080017089844D+00, &
       0.1727727502584457D+01, -0.2438052969955606D+02 /)
  real(8), save, dimension ( 4 ) :: b1 = (/ &
       -0.1025390625000000D+00, 0.2775764465332031D+00, &
       -0.1993531733751297D+01, 0.2724882731126854D+02 /)
  real(8)    :: bj(0:n)
  real(8)    :: bj0
  real(8)    :: bj1
  real(8)    :: bjk
  real(8)    :: bs
  real(8)    :: by(0:n)
  real(8)    :: by0
  real(8)    :: by1
  real(8)    :: byk
  real(8)    :: cu
  real(8)    :: dj(0:n)
  real(8)    :: dy(0:n)
  real(8)    :: ec
  real(8)    :: f
  real(8)    :: f1
  real(8)    :: f2
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: p0
  real(8)    :: p1
  real(8)    :: pi
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: r2p
  real(8)    :: s0
  real(8)    :: su
  real(8)    :: sv
  real(8)    :: t1
  real(8)    :: t2
  real(8)    :: x

  pi = 3.141592653589793D+00
  r2p = 0.63661977236758D+00
  nm = n

  if ( x < 1.0D-100 ) then
     do k = 0, n
        bj(k) = 0.0D+00
        dj(k) = 0.0D+00
        by(k) = -1.0D+300
        dy(k) = 1.0D+300
     end do
     bj(0) = 1.0D+00
     dj(1) = 0.5D+00
     return
  end if

  if ( x <= 300.0D+00 .or. int ( 0.9D+00 * x ) < n ) then

     if ( n == 0 ) then
        nm = 1
     end if

     m = msta1 ( x, 200 )

     if ( m < nm ) then
        nm = m
     else
        m = msta2 ( x, nm, 15 )
     end if

     bs = 0.0D+00
     su = 0.0D+00
     sv = 0.0D+00
     f2 = 0.0D+00
     f1 = 1.0D-100

     do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 - f2
        if ( k <= nm ) then
           bj(k) = f
        end if
        if ( k == 2 * int ( k / 2 ) .and. k /= 0 ) then
           bs = bs + 2.0D+00 * f
           su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
        else if ( 1 < k ) then
           sv = sv + ( -1.0D+00 ) ** ( k / 2 ) * k / ( k * k - 1.0D+00 ) * f
        end if
        f2 = f1
        f1 = f
     end do

     s0 = bs + f
     do k = 0, nm
        bj(k) = bj(k) / s0
     end do

     ec = log ( x / 2.0D+00 ) + 0.5772156649015329D+00
     by0 = r2p * ( ec * bj(0) - 4.0D+00 * su / s0 )
     by(0) = by0
     by1 = r2p * ( ( ec - 1.0D+00 ) * bj(1) - bj(0) / x - 4.0D+00 * sv / s0 )
     by(1) = by1

  else

     t1 = x - 0.25D+00 * pi
     p0 = 1.0D+00
     q0 = -0.125D+00 / x
     do k = 1, 4
        p0 = p0 + a(k) * x ** ( - 2 * k )
        q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
     end do
     cu = sqrt ( r2p / x )
     bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
     by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
     bj(0) = bj0
     by(0) = by0
     t2 = x - 0.75D+00 * pi
     p1 = 1.0D+00
     q1 = 0.375D+00 / x
     do k = 1, 4
        p1 = p1 + a1(k) * x ** ( - 2 * k )
        q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
     end do
     bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
     by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )
     bj(1) = bj1
     by(1) = by1
     do k = 2, nm
        bjk = 2.0D+00 * ( k - 1.0D+00 ) / x * bj1 - bj0
        bj(k) = bjk
        bj0 = bj1
        bj1 = bjk
     end do
  end if

  dj(0) = -bj(1)
  do k = 1, nm
     dj(k) = bj(k-1) - k / x * bj(k)
  end do

  do k = 2, nm
     byk = 2.0D+00 * ( k - 1.0D+00 ) * by1 / x - by0
     by(k) = byk
     by0 = by1
     by1 = byk
  end do

  dy(0) = -by(1)
  do k = 1, nm
     dy(k) = by(k-1) - k * by(k) / x
  end do

  return
end subroutine jynb




  !*****************************************************************************80
  !
  !! JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BJN, DJN, FJN, BYN, DYN, FYN, the values of
  !    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
  !
subroutine jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
  !JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param bjn: out: Jn(x)
  !@param djn: out: Jn'(x)
  !@param fjn: out: Jn"(x)
  !@param byn: out: Yn(x)
  !@param dyn: out: Yn'(x)
  !@param fyn: out: Yn"(x)
  !
  implicit none

  real(8)    :: bj(102)
  real(8)    :: bjn
  real(8)    :: byn
  real(8)    :: bs
  real(8)    :: by(102)
  real(8)    :: djn
  real(8)    :: dyn
  real(8)    :: e0
  real(8)    :: ec
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: fjn
  real(8)    :: fyn
  integer(4) :: k
  integer(4) :: m
  integer(4) :: mt
  integer(4) :: n
  integer(4) :: nt
  real(8)    :: s1
  real(8)    :: su
  real(8)    :: x

  do nt = 1, 900
     mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) &
          - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
     if ( 20 < mt ) then
        exit
     end if
  end do

  m = nt
  bs = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-35
  su = 0.0D+00
  do k = m, 0, -1
     f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
     if ( k <= n + 1 ) then
        bj(k+1) = f
     end if
     if ( k == 2 * int ( k / 2 ) ) then
        bs = bs + 2.0D+00 * f
        if ( k /= 0 ) then
           su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
        end if
     end if
     f0 = f1
     f1 = f
  end do

  do k = 0, n + 1
     bj(k+1) = bj(k+1) / ( bs - f )
  end do

  bjn = bj(n+1)
  ec = 0.5772156649015329D+00
  e0 = 0.3183098861837907D+00
  s1 = 2.0D+00 * e0 * ( log ( x / 2.0D+00 ) + ec ) * bj(1)
  f0 = s1 - 8.0D+00 * e0 * su / ( bs - f )
  f1 = ( bj(2) * f0 - 2.0D+00 * e0 / x ) / bj(1)

  by(1) = f0
  by(2) = f1
  do k = 2, n + 1 
     f = 2.0D+00 * ( k - 1.0D+00 ) * f1 / x - f0
     by(k+1) = f
     f0 = f1
     f1 = f
  end do

  byn = by(n+1)
  djn = - bj(n+2) + n * bj(n+1) / x
  dyn = - by(n+2) + n * by(n+1) / x
  fjn = ( n * n / ( x * x ) - 1.0D+00 ) * bjn - djn / x
  fyn = ( n * n / ( x * x ) - 1.0D+00 ) * byn - dyn / x

  return
end subroutine jyndd




  !*****************************************************************************80
  !
  !! JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Jv(x) and Yv(x).
  !
  !    Input, real(8) X, the argument of Jv(x) and Yv(x).
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N),
  !    the values of Jn+v0(x), Jn+v0'(x), Yn+v0(x), Yn+v0'(x).
  !
subroutine jyv ( v, x, vm, bj, dj, by, dy )
  !JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Jv(x) and Yv(x)
  !@param x: in: the argument of Jv(x) and Yv(x)
  !@param vm: out: the highest order computed
  !@param bj: out: Jn+v0(x)
  !@param dj: out: Jn+v0'(x)
  !@param by: out: Yn+v0(x)
  !@param dy: out: Yn+v0'(x)
  !
  implicit none

  real(8)    :: a
  real(8)    :: a0
  real(8)    :: b
  real(8)    :: bj(0:*)
  real(8)    :: bju0
  real(8)    :: bju1
  real(8)    :: bjv0
  real(8)    :: bjv1
  real(8)    :: bjvl
  real(8)    :: by(0:*)
  real(8)    :: byv0
  real(8)    :: byv1
  real(8)    :: byvk
  real(8)    :: ck
  real(8)    :: cs
  real(8)    :: cs0
  real(8)    :: cs1
  real(8)    :: dj(0:*)
  real(8)    :: dy(0:*)
  real(8)    :: ec
  real(8)    :: el
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: ga
  real(8)    :: gb
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: l
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: pv0
  real(8)    :: pv1
  real(8)    :: px
  real(8)    :: qx
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: rp
  real(8)    :: rp2
  real(8)    :: rq
  real(8)    :: sk
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vg
  real(8)    :: vl
  real(8)    :: vm
  real(8)    :: vv
  real(8)    :: w0
  real(8)    :: w1
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xk

  el = 0.5772156649015329D+00
  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  x2 = x * x
  n = int ( v )
  v0 = v - n

  if ( x < 1.0D-100 ) then

     do k = 0, n
        bj(k) = 0.0D+00
        dj(k) = 0.0D+00
        by(k) = -1.0D+300
        dy(k) = 1.0D+300
     end do

     if ( v0 == 0.0D+00 ) then
        bj(0) = 1.0D+00
        dj(1) = 0.5D+00
     else
        dj(0) = 1.0D+300
     end if
     vm = v  
     return

  end if

  if ( x <= 12.0D+00 ) then

     do l = 0, 1
        vl = v0 + l
        bjvl = 1.0D+00
        r = 1.0D+00
        do k = 1, 40
           r = -0.25D+00 * r * x2 / ( k * ( k + vl ) )
           bjvl = bjvl + r
           if ( abs ( r ) < abs ( bjvl ) * 1.0D-15 ) then
              exit
           end if
        end do

        vg = 1.0D+00 + vl
        call gammaf ( vg, ga )
        a = ( 0.5D+00 * x ) ** vl / ga

        if ( l == 0 ) then
           bjv0 = bjvl * a
        else
           bjv1 = bjvl * a
        end if

     end do

  else

     if ( x < 35.0D+00 ) then
        k0 = 11
     else if ( x < 50.0D+00 ) then
        k0 = 10
     else
        k0 = 8
     end if

     do j = 0, 1

        vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
        px = 1.0D+00
        rp = 1.0D+00
        do k = 1, k0
           rp = -0.78125D-02 * rp &
                * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                / ( k * ( 2.0D+00 * k - 1.0D+00 ) * x2 )
           px = px + rp
        end do
        qx = 1.0D+00
        rq = 1.0D+00
        do k = 1, k0
           rq = -0.78125D-02 * rq &
                * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
                / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
           qx = qx + rq
        end do
        qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
        xk = x - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
        a0 = sqrt ( rp2 / x )
        ck = cos ( xk )
        sk = sin ( xk )
        if ( j == 0 ) then
           bjv0 = a0 * ( px * ck - qx * sk )
           byv0 = a0 * ( px * sk + qx * ck )
        else if ( j == 1 ) then
           bjv1 = a0 * ( px * ck - qx * sk )
           byv1 = a0 * ( px * sk + qx * ck )
        end if

     end do

  end if

  bj(0) = bjv0
  bj(1) = bjv1
  dj(0) = v0 / x * bj(0) - bj(1)
  dj(1) = - ( 1.0D+00 + v0 ) / x * bj(1) + bj(0)

  if ( 2 <= n .and. n <= int ( 0.9D+00 * x ) ) then
     f0 = bjv0
     f1 = bjv1
     do k = 2, n
        f = 2.0D+00 * ( k + v0 - 1.0D+00 ) / x * f1 - f0
        bj(k) = f
        f0 = f1
        f1 = f
     end do
  else if ( 2 <= n ) then
     m = msta1 ( x, 200 )
     if ( m < n ) then
        n = m
     else
        m = msta2 ( x, n, 15 )
     end if
     f2 = 0.0D+00
     f1 = 1.0D-100
     do k = m, 0, -1
        f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 - f2
        if ( k <= n ) then
           bj(k) = f
        end if
        f2 = f1
        f1 = f
     end do

     if ( abs ( bjv1 ) < abs ( bjv0 ) ) then
        cs = bjv0 / f
     else
        cs = bjv1 / f2
     end if
     do k = 0, n
        bj(k) = cs * bj(k)
     end do
  end if

  do k = 2, n
     dj(k) = - ( k + v0 ) / x * bj(k) + bj(k-1)
  end do

  if ( x <= 12.0D+00 ) then

     if ( v0 /= 0.0D+00 ) then

        do l = 0, 1

           vl = v0 + l
           bjvl = 1.0D+00
           r = 1.0D+00
           do k = 1, 40
              r = -0.25D+00 * r * x2 / ( k * ( k - vl ) )
              bjvl = bjvl + r
              if ( abs ( r ) < abs ( bjvl ) * 1.0D-15 ) then
                 exit
              end if
           end do

           vg = 1.0D+00 - vl
           call gammaf ( vg, gb )
           b = ( 2.0D+00 / x ) ** vl / gb

           if ( l == 0 ) then
              bju0 = bjvl * b
           else
              bju1 = bjvl * b
           end if

        end do

        pv0 = pi * v0
        pv1 = pi * ( 1.0D+00 + v0 )
        byv0 = ( bjv0 * cos ( pv0 ) - bju0 ) / sin ( pv0 )
        byv1 = ( bjv1 * cos ( pv1 ) - bju1 ) / sin ( pv1 )

     else

        ec = log ( x / 2.0D+00 ) + el
        cs0 = 0.0D+00
        w0 = 0.0D+00
        r0 = 1.0D+00
        do k = 1, 30
           w0 = w0 + 1.0D+00 / k
           r0 = -0.25D+00 * r0 / ( k * k ) * x2
           cs0 = cs0 + r0 * w0
        end do
        byv0 = rp2 * ( ec * bjv0 - cs0 )
        cs1 = 1.0D+00
        w1 = 0.0D+00
        r1 = 1.0D+00
        do k = 1, 30
           w1 = w1 + 1.0D+00 / k
           r1 = -0.25D+00 * r1 / ( k * ( k + 1 ) ) * x2
           cs1 = cs1 + r1 * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
        end do
        byv1 = rp2 * ( ec * bjv1 - 1.0D+00 / x - 0.25D+00 * x * cs1 )

     end if

  end if

  by(0) = byv0
  by(1) = byv1
  do k = 2, n
     byvk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * byv1 - byv0
     by(k) = byvk
     byv0 = byv1
     byv1 = byvk
  end do

  dy(0) = v0 / x * by(0) - by(1)
  do k = 1, n
     dy(k) = - ( k + v0 ) / x * by(k) + by(k-1)
  end do

  vm = n + v0

  return
end subroutine jyv




  !*****************************************************************************80
  !
  !! JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of the Bessel functions.
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Output, real(8) RJ0(NT), RJ1(NT), RY0(NT), RY1(NT), the zeros 
  !    of Jn(x), Jn'(x), Yn(x), Yn'(x).
  !
subroutine jyzo ( n, nt, rj0, rj1, ry0, ry1 )
  !JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of the Bessel functions
  !@param nt: in: the number of zeros
  !@param rj0: out: Jn(x)
  !@param rj1: out: Jn'(x)
  !@param ry0: out: Yn(x)
  !@param ry1: out: Yn'(x)
  !
  implicit none

  integer(4) :: nt

  real(8)    :: bjn
  real(8)    :: byn
  real(8)    :: djn
  real(8)    :: dyn
  real(8)    :: fjn
  real(8)    :: fyn
  integer(4) :: l
  integer(4) :: n
  real(8)    :: n_r8
  real(8)    :: rj0(nt)
  real(8)    :: rj1(nt)
  real(8)    :: ry0(nt)
  real(8)    :: ry1(nt)
  real(8)    :: x
  real(8)    :: x0

  n_r8 = real ( n, kind = 8 )

  if ( n <= 20 ) then
     x = 2.82141D+00 + 1.15859D+00 * n_r8 
  else
     x = n + 1.85576D+00 * n_r8 ** 0.33333D+00 &
          + 1.03315D+00 / n_r8 ** 0.33333D+00
  end if

  l = 0

  do

     x0 = x
     call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
     x = x - bjn / djn

     if ( 1.0D-09 < abs ( x - x0 ) ) then
        cycle
     end if

     l = l + 1
     rj0(l) = x
     x = x + 3.1416D+00 + ( 0.0972D+00 + 0.0679D+00 * n_r8 &
          - 0.000354D+00 * n_r8 ** 2 ) / l

     if ( nt <= l ) then
        exit
     end if

  end do

  if ( n <= 20 ) then
     x = 0.961587D+00 + 1.07703D+00 * n_r8 
  else
     x = n_r8 + 0.80861D+00 * n_r8 ** 0.33333D+00 &
          + 0.07249D+00 / n_r8 ** 0.33333D+00
  end if

  if ( n == 0 ) then
     x = 3.8317D+00
  end if

  l = 0

  do

     x0 = x
     call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
     x = x - djn / fjn
     if ( 1.0D-09 < abs ( x - x0 ) ) then
        cycle
     end if
     l = l + 1
     rj1(l) = x
     x = x + 3.1416D+00 + ( 0.4955D+00 + 0.0915D+00 * n_r8 &
          - 0.000435D+00 * n_r8 ** 2 ) / l

     if ( nt <= l ) then
        exit
     end if

  end do

  if ( n <= 20 ) then
     x = 1.19477D+00 + 1.08933D+00 * n_r8 
  else
     x = n_r8 + 0.93158D+00 * n_r8 ** 0.33333D+00 &
          + 0.26035D+00 / n_r8 ** 0.33333D+00
  end if

  l = 0

  do

     x0 = x
     call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
     x = x - byn / dyn

     if ( 1.0D-09 < abs ( x - x0 ) ) then
        cycle
     end if

     l = l + 1
     ry0(l) = x 
     x = x + 3.1416D+00 + ( 0.312D+00 + 0.0852D+00 * n_r8 &
          - 0.000403D+00 * n_r8 ** 2 ) / l

     if ( nt <= l ) then
        exit
     end if

  end do

  if ( n <= 20 ) then
     x = 2.67257D+00 + 1.16099D+00 * n_r8 
  else
     x = n_r8 + 1.8211D+00 * n_r8 ** 0.33333D+00 &
          + 0.94001D+00 / n_r8 ** 0.33333D+00
  end if

  l = 0

  do

     x0 = x
     call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
     x = x - dyn / fyn

     if ( 1.0D-09 < abs ( x - x0 ) ) then
        cycle
     end if

     l = l + 1
     ry1(l) = x
     x = x + 3.1416D+00 + ( 0.197D+00 + 0.0643D+00 * n_r8 &
          -0.000286D+00 * n_r8 ** 2 ) / l 

     if ( nt <= l ) then
        exit
     end if

  end do

  return
end subroutine jyzo




  !*****************************************************************************80
  !
  !! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
  !    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
  !
subroutine klvna ( x, ber, bei, ger, gei, der, dei, her, hei )
  !KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ber: out: :code:`ber(x)`
  !@param bei: out: :code:`bei(x)`
  !@param ger: out: :code:`ker(x)`
  !@param gei: out: :code:`kei(x)`
  !@param der: out: derivative of :code:`ber(x)`
  !@param dei: out: derivative of :code:`bei(x)`
  !@param her: out: derivative of :code:`ker(x)`
  !@param hei: out: derivative of :code:`kei(x)`
  !
  implicit none

  real(8)    :: bei
  real(8)    :: ber
  real(8)    :: cn0
  real(8)    :: cp0
  real(8)    :: cs
  real(8)    :: dei
  real(8)    :: der
  real(8)    :: el
  real(8)    :: eps
  real(8)    :: fac
  real(8)    :: gei
  real(8)    :: ger
  real(8)    :: gs
  real(8)    :: hei
  real(8)    :: her
  integer(4) :: k
  integer(4) :: km
  integer(4) :: m
  real(8)    :: pi
  real(8)    :: pn0
  real(8)    :: pn1
  real(8)    :: pp0
  real(8)    :: pp1
  real(8)    :: qn0
  real(8)    :: qn1
  real(8)    :: qp0
  real(8)    :: qp1
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: rc
  real(8)    :: rs
  real(8)    :: sn0
  real(8)    :: sp0
  real(8)    :: ss
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: x4
  real(8)    :: xc1
  real(8)    :: xc2
  real(8)    :: xd
  real(8)    :: xe1
  real(8)    :: xe2
  real(8)    :: xt

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  eps = 1.0D-15

  if ( x == 0.0D+00 ) then
     ber = 1.0D+00
     bei = 0.0D+00
     ger = 1.0D+300
     gei = -0.25D+00 * pi
     der = 0.0D+00
     dei = 0.0D+00
     her = -1.0D+300
     hei = 0.0D+00
     return
  end if

  x2 = 0.25D+00 * x * x
  x4 = x2 * x2

  if ( abs ( x ) < 10.0D+00 ) then

     ber = 1.0D+00
     r = 1.0D+00
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m - 1.0D+00 ) ** 2 * x4
        ber = ber + r
        if ( abs ( r ) < abs ( ber ) * eps ) then
           exit
        end if
     end do

     bei = x2
     r = x2
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
        bei = bei + r
        if ( abs ( r ) < abs ( bei ) * eps ) then
           exit
        end if
     end do

     ger = - ( log ( x / 2.0D+00 ) + el ) * ber + 0.25D+00 * pi * bei
     r = 1.0D+00
     gs = 0.0D+00
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m - 1.0D+00 ) ** 2 * x4
        gs = gs + 1.0D+00 / ( 2.0D+00 * m - 1.0D+00 ) + 1.0D+00 / ( 2.0D+00 * m )
        ger = ger + r * gs
        if ( abs ( r * gs ) < abs ( ger ) * eps ) then
           exit
        end if
     end do

     gei = x2 - ( log ( x / 2.0D+00 ) + el ) * bei - 0.25D+00 * pi * ber
     r = x2
     gs = 1.0D+00
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
        gs = gs + 1.0D+00 / ( 2.0D+00 * m ) + 1.0D+00 / ( 2.0D+00 * m + 1.0D+00 )
        gei = gei + r * gs
        if ( abs ( r * gs ) < abs ( gei ) * eps ) then
           exit
        end if
     end do

     der = -0.25D+00 * x * x2
     r = der
     do m = 1, 60
        r = -0.25D+00 * r / m / ( m + 1.0D+00 ) &
             / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
        der = der + r
        if ( abs ( r ) < abs ( der ) * eps ) then
           exit
        end if
     end do

     dei = 0.5D+00 * x
     r = dei
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m - 1.0D+00 ) &
             / ( 2.0D+00 * m + 1.0D+00 ) * x4
        dei = dei + r
        if ( abs ( r ) < abs ( dei ) * eps ) then
           exit
        end if
     end do

     r = -0.25D+00 * x * x2
     gs = 1.5D+00
     her = 1.5D+00 * r - ber / x &
          - ( log ( x / 2.0D+00 ) + el ) * der + 0.25D+00 * pi * dei
     do m = 1, 60
        r = -0.25D+00 * r / m / ( m + 1.0D+00 ) &
             / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
        gs = gs + 1.0D+00 / ( 2 * m + 1.0D+00 ) + 1.0D+00 &
             / ( 2 * m + 2.0D+00 )
        her = her + r * gs
        if ( abs ( r * gs ) < abs ( her ) * eps ) then
           exit
        end if
     end do

     r = 0.5D+00 * x
     gs = 1.0D+00
     hei = 0.5D+00 * x - bei / x &
          - ( log ( x / 2.0D+00 ) + el ) * dei - 0.25D+00 * pi * der
     do m = 1, 60
        r = -0.25D+00 * r / ( m * m ) / ( 2 * m - 1.0D+00 ) &
             / ( 2 * m + 1.0D+00 ) * x4
        gs = gs + 1.0D+00 / ( 2.0D+00 * m ) + 1.0D+00 &
             / ( 2 * m + 1.0D+00 )
        hei = hei + r * gs
        if ( abs ( r * gs ) < abs ( hei ) * eps ) then 
           return
        end if
     end do

  else

     pp0 = 1.0D+00
     pn0 = 1.0D+00
     qp0 = 0.0D+00
     qn0 = 0.0D+00
     r0 = 1.0D+00

     if ( abs ( x ) < 40.0D+00 ) then
        km = 18
     else
        km = 10
     end if

     fac = 1.0D+00
     do k = 1, km
        fac = -fac
        xt = 0.25D+00 * k * pi - int ( 0.125D+00 * k ) * 2.0D+00 * pi
        cs = cos ( xt )
        ss = sin ( xt )
        r0 = 0.125D+00 * r0 * ( 2.0D+00 * k - 1.0D+00 ) ** 2 / k / x
        rc = r0 * cs
        rs = r0 * ss
        pp0 = pp0 + rc
        pn0 = pn0 + fac * rc
        qp0 = qp0 + rs
        qn0 = qn0 + fac * rs
     end do

     xd = x / sqrt (2.0D+00 )
     xe1 = exp ( xd )
     xe2 = exp ( - xd )
     xc1 = 1.0D+00 / sqrt ( 2.0D+00 * pi * x )
     xc2 = sqrt ( 0.5D+00 * pi / x )
     cp0 = cos ( xd + 0.125D+00 * pi )
     cn0 = cos ( xd - 0.125D+00 * pi )
     sp0 = sin ( xd + 0.125D+00 * pi )
     sn0 = sin ( xd - 0.125D+00 * pi )
     ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
     gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
     ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
     bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
     pp1 = 1.0D+00
     pn1 = 1.0D+00
     qp1 = 0.0D+00
     qn1 = 0.0D+00
     r1 = 1.0D+00
     fac = 1.0D+00

     do k = 1, km
        fac = -fac
        xt = 0.25D+00 * k * pi - int ( 0.125D+00 * k ) * 2.0D+00 * pi
        cs = cos ( xt )
        ss = sin ( xt )
        r1 = 0.125D+00 * r1 &
             * ( 4.0D+00 - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / k / x
        rc = r1 * cs
        rs = r1 * ss
        pp1 = pp1 + fac * rc
        pn1 = pn1 + rc
        qp1 = qp1 + fac * rs
        qn1 = qn1 + rs
     end do

     her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
     hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
     der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
     dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

  end if

  return
end subroutine klvna




  !*****************************************************************************80
  !
  !! KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    03 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
  !    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
  !
subroutine klvnb ( x, ber, bei, ger, gei, der, dei, her, hei )
  !KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ber: out: :code:`ber(x)`
  !@param bei: out: :code:`bei(x)`
  !@param ger: out: :code:`ker(x)`
  !@param gei: out: :code:`kei(x)`
  !@param der: out: derivative of :code:`ber(x)`
  !@param dei: out: derivative of :code:`bei(x)`
  !@param her: out: derivative of :code:`ker(x)`
  !@param hei: out: derivative of :code:`kei(x)`
  !
  implicit none

  real(8)    :: bei
  real(8)    :: ber
  real(8)    :: csn
  real(8)    :: csp
  real(8)    :: dei
  real(8)    :: der
  real(8)    :: fxi
  real(8)    :: fxr
  real(8)    :: gei
  real(8)    :: ger
  real(8)    :: hei
  real(8)    :: her
  integer(4) :: l
  real(8)    :: pi
  real(8)    :: pni
  real(8)    :: pnr
  real(8)    :: ppi
  real(8)    :: ppr
  real(8)    :: ssn
  real(8)    :: ssp
  real(8)    :: t
  real(8)    :: t2
  real(8)    :: tni
  real(8)    :: tnr
  real(8)    :: tpi
  real(8)    :: tpr
  real(8)    :: u
  real(8)    :: v
  real(8)    :: x
  real(8)    :: yc1
  real(8)    :: yc2
  real(8)    :: yci
  real(8)    :: ye1
  real(8)    :: ye2
  real(8)    :: yei
  real(8)    :: yd

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

     ber = 1.0D+00
     bei = 0.0D+00
     ger = 1.0D+300
     gei = -0.25D+00 * pi
     der = 0.0D+00
     dei = 0.0D+00
     her = -1.0D+300
     hei = 0.0D+00

  else if ( x < 8.0D+00 ) then

     t = x / 8.0D+00
     t2 = t * t
     u = t2 * t2

     ber = (((((( &
          - 0.901D-05 * u &
          + 0.122552D-02 ) * u &
          - 0.08349609D+00 ) * u &
          + 2.64191397D+00 ) * u &
          - 32.36345652D+00 ) * u &
          + 113.77777774D+00 ) * u &
          - 64.0D+00 ) * u &
          + 1.0D+00

     bei = t * t * (((((( &
          0.11346D-03 * u &
          - 0.01103667D+00 ) * u &
          + 0.52185615D+00 ) * u &
          - 10.56765779D+00 ) * u &
          + 72.81777742D+00 ) * u &
          - 113.77777774D+00 ) * u &
          + 16.0D+00 )

     ger = (((((( &
          - 0.2458D-04 * u &
          + 0.309699D-02 ) * u &
          - 0.19636347D+00 ) * u &
          + 5.65539121D+00 ) * u &
          - 60.60977451D+00 ) * u &
          + 171.36272133D+00 ) * u &
          - 59.05819744D+00 ) * u &
          - 0.57721566D+00

     ger = ger - log ( 0.5D+00 * x ) * ber + 0.25D+00 * pi * bei

     gei = t2 * (((((( &
          0.29532D-03 * u &
          - 0.02695875D+00 ) * u &
          + 1.17509064D+00 ) * u &
          - 21.30060904D+00 ) * u &
          + 124.2356965D+00 ) * u &
          - 142.91827687D+00 ) * u &
          + 6.76454936D+00 )

     gei = gei - log ( 0.5D+00 * x ) * bei - 0.25D+00 * pi * ber

     der = x * t2 * (((((( &
          - 0.394D-05 * u &
          + 0.45957D-03 ) * u &
          - 0.02609253D+00 ) * u &
          + 0.66047849D+00 ) * u &
          - 6.0681481D+00 ) * u &
          + 14.22222222D+00 ) * u &
          - 4.0D+00 )

     dei = x * (((((( &
          0.4609D-04 * u &
          - 0.379386D-02 ) * u &
          + 0.14677204D+00 ) * u &
          - 2.31167514D+00 ) * u &
          + 11.37777772D+00 ) * u &
          - 10.66666666D+00 ) * u &
          + 0.5D+00 ) 

     her = x * t2 * (((((( &
          - 0.1075D-04 * u &
          + 0.116137D-02 ) * u &
          - 0.06136358D+00 ) * u &
          + 1.4138478D+00 ) * u &
          - 11.36433272D+00 ) * u &
          + 21.42034017D+00 ) * u &
          - 3.69113734D+00 )

     her = her - log ( 0.5D+00 * x ) * der - ber / x  &
          + 0.25D+00 * pi * dei

     hei = x * (((((( &
          0.11997D-03 * u &
          - 0.926707D-02 ) * u &
          + 0.33049424D+00 ) * u &
          - 4.65950823D+00 ) * u &
          + 19.41182758D+00 ) * u &
          - 13.39858846D+00 ) * u &
          + 0.21139217D+00 )

     hei = hei - log ( 0.5D+00 * x ) * dei - bei / x  &
          - 0.25D+00 * pi * der

  else

     t = 8.0D+00 / x

     do l = 1, 2

        v = ( -1.0D+00 ) ** l * t

        tpr = (((( &
             0.6D-06 * v &
             - 0.34D-05 ) * v &
             - 0.252D-04 ) * v &
             - 0.906D-04 ) * v * v &
             + 0.0110486D+00 ) * v

        tpi = (((( &
             0.19D-05 * v &
             + 0.51D-05 ) * v * v &
             - 0.901D-04 ) * v &
             - 0.9765D-03 ) * v &
             - 0.0110485D+00 ) * v &
             - 0.3926991D+00

        if ( l == 1 ) then
           tnr = tpr
           tni = tpi
        end if

     end do

     yd = x / sqrt ( 2.0D+00 )
     ye1 = exp ( yd + tpr )
     ye2 = exp ( - yd + tnr )
     yc1 = 1.0D+00 / sqrt ( 2.0D+00 * pi * x )
     yc2 = sqrt ( pi / ( 2.0D+00 * x ) )
     csp = cos ( yd + tpi )
     ssp = sin ( yd + tpi )
     csn = cos ( - yd + tni )
     ssn = sin ( - yd + tni )
     ger = yc2 * ye2 * csn
     gei = yc2 * ye2 * ssn
     fxr = yc1 * ye1 * csp
     fxi = yc1 * ye1 * ssp
     ber = fxr - gei / pi
     bei = fxi + ger / pi

     do l = 1, 2

        v = ( -1.0D+00 ) ** l * t

        ppr = ((((( &
             0.16D-05 * v &
             + 0.117D-04 ) * v &
             + 0.346D-04 ) * v &
             + 0.5D-06 ) * v &
             - 0.13813D-02 ) * v &
             - 0.0625001D+00 ) * v &
             + 0.7071068D+00

        ppi = ((((( &
             - 0.32D-05 * v &
             - 0.24D-05 ) * v &
             + 0.338D-04 ) * v &
             + 0.2452D-03 ) * v &
             + 0.13811D-02 ) * v &
             - 0.1D-06 ) * v &
             + 0.7071068D+00

        if ( l == 1 ) then
           pnr = ppr
           pni = ppi
        end if

     end do

     her =     gei * pni - ger * pnr
     hei = - ( gei * pnr + ger * pni )
     der = fxr * ppr - fxi * ppi - hei / pi
     dei = fxi * ppr + fxr * ppi + her / pi

  end if

  return
end subroutine klvnb




  !*****************************************************************************80
  !
  !! KLVNZO computes zeros of the Kelvin functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) NT, the number of zeros.
  !
  !    Input, integer(4) KD, the function code.
  !    1 for ber x, 
  !    2 for bei x,
  !    3 for ker x, 
  !    4 for kei x,
  !    5 for ber' x, 
  !    6 for bei' x,
  !    7 for ker' x, 
  !    8 for kei' x.
  !
  !    Output, real(8) ZO(NT), the zeros of the given Kelvin function.
  !
subroutine klvnzo ( nt, kd, zo )
  !KLVNZO computes zeros of the Kelvin functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Possible values of the function code :f:var:`kd`:
  !
  !* :code:`1`: :code:`ber(x)`
  !* :code:`2`: :code:`bei(x)`
  !* :code:`3`: :code:`ker(x)`
  !* :code:`4`: :code:`kei(x)`
  !* :code:`5`: derivative of :code:`ber(x)`
  !* :code:`6`: derivative of :code:`bei(x)`
  !* :code:`7`: derivative of :code:`ker(x)`
  !* :code:`8`: derivative of :code:`kei(x)`
  !
  !@param nt: in: the number of zeros
  !@param kd: in: function code
  !@param zo: out: the zeros of the given Kelvin function
  !
  implicit none

  integer(4) :: nt

  real(8)    :: bei
  real(8)    :: ber
  real(8)    :: ddi
  real(8)    :: ddr
  real(8)    :: dei
  real(8)    :: der
  real(8)    :: gdi
  real(8)    :: gdr
  real(8)    :: gei
  real(8)    :: ger
  real(8)    :: hei
  real(8)    :: her
  integer(4) :: kd
  integer(4) :: m
  real(8)    :: rt
  real(8)    :: rt0(8)
  real(8)    :: zo(nt)

  rt0(1) = 2.84891D+00
  rt0(2) = 5.02622D+00
  rt0(3) = 1.71854D+00
  rt0(4) = 3.91467D+00
  rt0(5) = 6.03871D+00
  rt0(6) = 3.77268D+00
  rt0(7) = 2.66584D+00
  rt0(8) = 4.93181D+00

  rt = rt0(kd)

  do m = 1, nt

     do

        call klvna ( rt, ber, bei, ger, gei, der, dei, her, hei )

        if ( kd == 1 ) then
           rt = rt - ber / der
        else if ( kd == 2 ) then
           rt = rt - bei / dei
        else if ( kd == 3 ) then
           rt = rt - ger / her
        else if ( kd == 4 ) then
           rt = rt - gei / hei
        else if ( kd == 5 ) then
           ddr = - bei - der / rt
           rt = rt - der / ddr
        else if ( kd == 6 ) then
           ddi = ber - dei / rt
           rt = rt - dei / ddi
        else if ( kd == 7 ) then
           gdr = - gei - her / rt
           rt = rt - her / gdr
        else
           gdi = ger - hei / rt
           rt = rt - hei / gdi
        end if

        if ( abs ( rt - rt0(kd) ) <= 5.0D-10 ) then
           exit
        end if

        rt0(kd) = rt

     end do

     zo(m) = rt
     rt = rt + 4.44D+00

  end do

  return
end subroutine klvnzo




  !*****************************************************************************80
  !
  !! KMN: expansion coefficients of prolate or oblate spheroidal functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    02 August 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Input, real(8) DF(*), the expansion coefficients.
  !
subroutine kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
  !KMN: expansion coefficients of prolate or oblate spheroidal functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param cv: in: the characteristic value
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param df: in: the expansion coefficients
  !@param dn: expansion coefficients
  !@param ck1: joining factor
  !@param ck2: joining factor
  !
  implicit none

  real(8)    :: c
  real(8)    :: ck1
  real(8)    :: ck2
  real(8)    :: cs
  real(8)    :: cv
  real(8)    :: df(200)
  real(8)    :: dn(200)
  real(8)    :: dnp
  real(8)    :: g0
  real(8)    :: gk0
  real(8)    :: gk1
  real(8)    :: gk2
  real(8)    :: gk3
  integer(4) :: i
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: l
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm1
  integer(4) :: nn
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: r3
  real(8)    :: r4
  real(8)    :: r5
  real(8)    :: rk(200)
  real(8)    :: sa0
  real(8)    :: sb0
  real(8)    :: su0
  real(8)    :: sw
  real(8)    :: t
  real(8)    :: tp(200)
  real(8)    :: u(200)
  real(8)    :: v(200)
  real(8)    :: w(200)

  nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
  nn = nm + m
  cs = c * c * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  do i = 1, nn + 3

     if ( ip == 0 ) then
        k = - 2 * ( i - 1 )
     else
        k = - ( 2 * i - 3 )
     end if

     gk0 = 2.0D+00 * m + k
     gk1 = ( m + k ) * ( m + k + 1.0D+00 )
     gk2 = 2.0D+00 * ( m + k ) - 1.0D+00
     gk3 = 2.0D+00 * ( m + k ) + 3.0D+00
     u(i) = gk0 * ( gk0 - 1.0D+00 ) * cs / ( gk2 * ( gk2 + 2.0D+00 ) )
     v(i) = gk1 - cv + ( 2.0D+00 * ( gk1 - m * m ) - 1.0D+00 ) * cs &
          / ( gk2 * gk3 )
     w(i) = ( k + 1.0D+00 ) * ( k + 2.0D+00 ) * cs / ( ( gk2 + 2.0D+00 ) * gk3 )

  end do

  do k = 1, m
     t = v(m+1)
     do l = 0, m - k - 1
        t = v(m-l) - w(m-l+1) * u(m-l) / t
     end do
     rk(k) = -u(k) / t
  end do

  r = 1.0D+00
  do k = 1, m
     r = r * rk(k)
     dn(k) = df(1) * r
  end do

  tp(nn) = v(nn+1)
  do k = nn - 1, m + 1,-1
     tp(k) = v(k+1) - w(k+2) * u(k+1) / tp(k+1)
     if ( m + 1 < k ) then
        rk(k) = -u(k) / tp(k)
     end if
  end do

  if ( m == 0 ) then
     dnp = df(1)
  else
     dnp = dn(m)
  end if

  dn(m+1) = ( - 1.0D+00 ) ** ip * dnp * cs &
       / ( ( 2.0D+00 * m - 1.0D+00 ) &
       * ( 2.0D+00 * m + 1.0D+00 - 4.0D+00 * ip ) * tp(m+1) )
  do k = m + 2, nn
     dn(k) = rk(k) * dn(k-1)
  end do

  r1 = 1.0D+00
  do j = 1, ( n + m + ip ) / 2
     r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
  end do
  nm1 = ( n - m ) / 2
  r = 1.0D+00
  do j = 1, 2 * m + ip
     r = r * j
  end do
  su0 = r * df(1)

  do k = 2, nm
     r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
          / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     su0 = su0 + r * df(k)
     if ( nm1 < k .and. &
          abs ( ( su0 - sw ) / su0 ) < 1.0D-14 ) then
        exit
     end if
     sw = su0
  end do

  if ( kd /= 1 ) then

     r2 = 1.0D+00
     do j = 1,m
        r2 = 2.0D+00 * c * r2 * j
     end do
     r3 = 1.0D+00
     do j = 1, ( n - m - ip ) / 2
        r3 = r3 * j
     end do
     sa0 = ( 2.0D+00 * ( m + ip ) + 1.0D+00 ) * r1 &
          / ( 2.0D+00 ** n * c ** ip * r2 * r3 * df(1) )
     ck1 = sa0 * su0

     if ( kd == -1 ) then
        return
     end if

  end if

  r4 = 1.0D+00
  do j = 1, ( n - m - ip ) / 2
     r4 = 4.0D+00 * r4 * j
  end do
  r5 = 1.0D+00
  do j = 1, m
     r5 = r5 * ( j + m ) / c
  end do

  if ( m == 0 ) then
     g0 = df(1)
  else
     g0 = dn(m)
  end if

  sb0 = ( ip + 1.0D+00 ) * c ** ( ip + 1 ) &
       / ( 2.0D+00 * ip * ( m - 2.0D+00 ) + 1.0D+00 ) &
       / ( 2.0D+00 * m - 1.0D+00 )

  ck2 = ( -1 ) ** ip * sb0 * r4 * r5 * g0 / r1 * su0

  return
end subroutine kmn




  !*****************************************************************************80
  !
  !! LAGZO computes zeros of the Laguerre polynomial, and integration weights.
  !
  !  Discussion:
  !
  !    This procedure computes the zeros of Laguerre polynomial Ln(x) in the 
  !    interval [0,∞], and the corresponding weighting coefficients for 
  !    Gauss-Laguerre integration.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of the Laguerre polynomial.
  !
  !    Output, real(8) X(N), the zeros of the Laguerre polynomial.
  !
  !    Output, real(8) W(N), the weighting coefficients.
  !
subroutine lagzo ( n, x, w )
  !LAGZO computes zeros of the Laguerre polynomial, and integration weights.
  !
  !This procedure computes the zeros of Laguerre polynomial Ln(x) in the
  !interval [0,infinity], and the corresponding weighting coefficients for
  !Gauss-Laguerre integration.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of the Laguerre polynomial
  !@param x: out: the zeros of the Laguerre polynomial
  !@param w: out: the weighting coefficients
  !
  implicit none

  integer(4) :: n

  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: fd
  real(8)    :: gd
  real(8)    :: hn
  integer(4) :: i
  integer(4) :: it
  integer(4) :: j
  integer(4) :: k
  integer(4) :: nr
  real(8)    :: p
  real(8)    :: pd
  real(8)    :: pf
  real(8)    :: q
  real(8)    :: w(n)
  real(8)    :: wp
  real(8)    :: x(n)
  real(8)    :: z
  real(8)    :: z0

  hn = 1.0D+00 / real ( n, kind = 8 )

  do nr = 1, n

     if ( nr == 1 ) then
        z = hn
     else
        z = x(nr-1) + hn * nr ** 1.27D+00
     end if

     it = 0

     do

        it = it + 1
        z0 = z
        p = 1.0D+00
        do i = 1, nr - 1
           p = p * ( z - x(i) )
        end do

        f0 = 1.0D+00
        f1 = 1.0D+00 - z
        do k = 2, n
           pf = (( 2.0D+00 * k - 1.0D+00 - z ) * f1 &
                - ( k - 1.0D+00 ) * f0 ) / k
           pd = k / z * ( pf - f1 )
           f0 = f1
           f1 = pf
        end do

        fd = pf / p

        q = 0.0D+00
        do i = 1, nr - 1
           wp = 1.0D+00
           do j = 1, nr - 1
              if ( j /= i ) then
                 wp = wp * ( z - x(j) )
              end if
           end do
           q = q + wp
        end do

        gd = ( pd - q * fd ) / p
        z = z - fd / gd

        if ( 40 < it .or. abs ( ( z - z0 ) / z ) <= 1.0D-15 ) then
           exit
        end if

     end do

     x(nr) = z
     w(nr) = 1.0D+00 / ( z * pd * pd )

  end do

  return
end subroutine lagzo




  !*****************************************************************************80
  !
  !! LAMN computes lambda functions and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    14 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) BL(0:N), DL(0:N), the
  !    value of the lambda function and its derivative of orders 0 through N.
  !
subroutine lamn ( n, x, nm, bl, dl )
  !LAMN computes lambda functions and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param bl: out: the lambda function
  !@param dl: out: its derivative of orders 0 through N
  !
  implicit none

  integer(4) :: n

  real(8)    :: bg
  real(8)    :: bk
  real(8)    :: bl(0:n)
  real(8)    :: bs
  real(8)    :: dl(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: i
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: uk
  real(8)    :: x
  real(8)    :: x2

  nm = n

  if ( abs ( x ) < 1.0D-100 ) then
     do k = 0, n
        bl(k) = 0.0D+00
        dl(k) = 0.0D+00
     end do
     bl(0) = 1.0D+00
     dl(1) = 0.5D+00
     return
  end if

  if ( x <= 12.0D+00 ) then

     x2 = x * x

     do k = 0, n
        bk = 1.0D+00
        r = 1.0D+00
        do i = 1, 50
           r = -0.25D+00 * r * x2 / ( i * ( i + k ) )
           bk = bk + r
           if ( abs ( r ) < abs ( bk ) * 1.0D-15 ) then
              exit
           end if
        end do

        bl(k) = bk
        if ( 1 <= k ) then
           dl(k-1) = - 0.5D+00 * x / k * bk
        end if

     end do

     uk = 1.0D+00
     r = 1.0D+00
     do i = 1, 50
        r = -0.25D+00 * r * x2 / ( i * ( i + n + 1.0D+00 ) )
        uk = uk + r
        if ( abs ( r ) < abs ( uk ) * 1.0D-15 ) then
           exit
        end if
     end do

     dl(n) = -0.5D+00 * x / ( n + 1.0D+00 ) * uk
     return

  end if

  if ( n == 0 ) then
     nm = 1
  end if

  m = msta1 ( x, 200 )

  if ( m < nm ) then
     nm = m
  else
     m = msta2 ( x, nm, 15 )
  end if

  bs = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-100
  do k = m, 0, -1
     f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
     if ( k <= nm ) then
        bl(k) = f
     end if
     if ( k == 2 * int ( k / 2 ) ) then
        bs = bs + 2.0D+00 * f
     end if
     f0 = f1
     f1 = f
  end do

  bg = bs - f
  do k = 0, nm
     bl(k) = bl(k) / bg
  end do

  r0 = 1.0D+00
  do k = 1, nm
     r0 = 2.0D+00 * r0 * k / x
     bl(k) = r0 * bl(k)
  end do

  dl(0) = -0.5D+00 * x * bl(1)
  do k = 1, nm
     dl(k) = 2.0D+00 * k / x * ( bl(k-1) - bl(k) )
  end do

  return
end subroutine lamn




  !*****************************************************************************80
  !
  !! LAMV computes lambda functions and derivatives of arbitrary order.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) VM, the highest order computed.
  !
  !    Output, real(8) VL(0:*), DL(0:*), the Lambda function and 
  !    derivative, of orders N+V0.
  !
subroutine lamv ( v, x, vm, vl, dl )
  !LAMV computes lambda functions and derivatives of arbitrary order.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order
  !@param x: in: the argument
  !@param vm: out: the highest order computed
  !@param vl: out: the Lambda function, of orders N+V0
  !@param dl: out: the derivative of the Lambda function, of orders N+V0
  !
  implicit none

  real(8)    :: v

  real(8)    :: a0
  real(8)    :: bjv0
  real(8)    :: bjv1
  real(8)    :: bk
  real(8)    :: ck
  real(8)    :: cs
  real(8)    :: dl(0:int(v))
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: fac
  real(8)    :: ga
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k0
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: px
  real(8)    :: qx
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: rc
  real(8)    :: rp
  real(8)    :: rp2
  real(8)    :: rq
  real(8)    :: sk
  real(8)    :: uk
  real(8)    :: v0
  real(8)    :: vk
  real(8)    :: vl(0:int(v))
  real(8)    :: vm
  real(8)    :: vv
  real(8)    :: x
  real(8)    :: x2
  real(8)    :: xk

  pi = 3.141592653589793D+00
  rp2 = 0.63661977236758D+00
  x = abs ( x )
  x2 = x * x
  n = int ( v )
  v0 = v - n
  vm = v

  if ( x <= 12.0D+00 ) then

     do k = 0, n

        vk = v0 + k
        bk = 1.0D+00
        r = 1.0D+00

        do i = 1, 50
           r = -0.25D+00 * r * x2 / ( i * ( i + vk ) )
           bk = bk + r
           if ( abs ( r ) < abs ( bk ) * 1.0D-15 ) then
              exit
           end if
        end do

        vl(k) = bk
        uk = 1.0D+00
        r = 1.0D+00
        do i = 1, 50
           r = -0.25D+00 * r * x2 / ( i * ( i + vk + 1.0D+00 ))
           uk = uk + r
           if ( abs ( r ) < abs ( uk ) * 1.0D-15 ) then
              exit
           end if
        end do

        dl(k) = - 0.5D+00 * x / ( vk + 1.0D+00 ) * uk

     end do

     return

  end if

  if ( x < 35.0D+00 ) then
     k0 = 11
  else if ( x < 50.0D+00 ) then
     k0 = 10
  else
     k0 = 8
  end if

  do j = 0, 1
     vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
     px = 1.0D+00
     rp = 1.0D+00
     do k = 1, k0
        rp = - 0.78125D-02 * rp * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
             * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
             / ( k * ( 2.0 * k - 1.0D+00 ) * x2 )
        px = px + rp
     end do
     qx = 1.0D+00
     rq = 1.0D+00
     do k = 1, k0
        rq = - 0.78125D-02 * rq * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
             * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
             / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
        qx = qx + rq
     end do
     qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
     xk = x - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
     a0 = sqrt ( rp2 / x )
     ck = cos ( xk )
     sk = sin ( xk )
     if ( j == 0 ) then
        bjv0 = a0 * ( px * ck - qx * sk )
     else
        bjv1 = a0 * ( px * ck - qx * sk )
     end if
  end do

  if ( v0 == 0.0D+00 ) then
     ga = 1.0D+00
  else
     call gam0 ( v0, ga )
     ga = v0 * ga
  end if

  fac = ( 2.0D+00 / x ) ** v0 * ga
  vl(0) = bjv0
  dl(0) = - bjv1 + v0 / x * bjv0
  vl(1) = bjv1
  dl(1) = bjv0 - ( 1.0D+00 + v0 ) / x * bjv1
  r0 = 2.0D+00 * ( 1.0D+00 + v0 ) / x

  if ( n <= 1 ) then
     vl(0) = fac * vl(0)
     dl(0) = fac * dl(0) - v0 / x * vl(0)
     vl(1) = fac * r0 * vl(1)
     dl(1) = fac * r0 * dl(1) - ( 1.0D+00 + v0 ) / x * vl(1)
     return
  end if

  if ( 2 <= n .and. n <= int ( 0.9D+00 * x ) ) then

     f0 = bjv0
     f1 = bjv1
     do k = 2, n
        f = 2.0D+00 * ( k + v0 - 1.0D+00 ) / x * f1 - f0
        f0 = f1
        f1 = f
        vl(k) = f
     end do

  else if ( 2 <= n ) then

     m = msta1 ( x, 200 )
     if ( m < n ) then
        n = m
     else
        m = msta2 ( x, n, 15 )
     end if
     f2 = 0.0D+00
     f1 = 1.0D-100
     do k = m, 0, -1
        f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 - f2
        if ( k <= n ) then
           vl(k) = f
        end if
        f2 = f1
        f1 = f
     end do

     if ( abs ( bjv0 ) <= abs ( bjv1 ) ) then
        cs = bjv1 / f2
     else
        cs = bjv0 / f
     end if

     do k = 0, n
        vl(k) = cs * vl(k)
     end do

  end if

  vl(0) = fac * vl(0)
  do j = 1, n
     rc = fac * r0
     vl(j) = rc * vl(j)
     dl(j-1) = - 0.5D+00 * x / ( j + v0 ) * vl(j)
     r0 = 2.0D+00 * ( j + v0 + 1 ) / x * r0
  end do
  dl(n) = 2.0D+00 * ( v0 + n ) * ( vl(n-1) - vl(n) ) / x
  vm = n + v0

  return
end subroutine lamv




  !*****************************************************************************80
  !
  !! LEGZO computes the zeros of Legendre polynomials, and integration weights.
  !
  !  Discussion:
  !
  !    This procedure computes the zeros of Legendre polynomial Pn(x) in the 
  !    interval [-1,1], and the corresponding weighting coefficients for 
  !    Gauss-Legendre integration.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    13 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of the polynomial.
  !
  !    Output, real(8) X(N), W(N), the zeros of the polynomial,
  !    and the corresponding weights.
  !
subroutine legzo ( n, x, w )
  !LEGZO computes the zeros of Legendre polynomials, and integration weights.
  !
  !This procedure computes the zeros of Legendre polynomial Pn(x) in the
  !interval [-1,1], and the corresponding weighting coefficients for
  !Gauss-Legendre integration.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of the polynomial
  !@param x: out: the polynomial
  !@param w: out: and the corresponding weights
  !
  implicit none

  integer(4) :: n

  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: fd
  real(8)    :: gd
  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: n0
  integer(4) :: nr
  real(8)    :: p
  real(8)    :: pd
  real(8)    :: pf
  real(8)    :: q
  real(8)    :: w(n)
  real(8)    :: wp
  real(8)    :: x(n)
  real(8)    :: z
  real(8)    :: z0

  n0 = ( n + 1 ) / 2

  do nr = 1, n0

     z = cos ( 3.1415926D+00 * ( nr - 0.25D+00 ) / n )

     do

        z0 = z
        p = 1.0D+00
        do i = 1, nr - 1
           p = p * ( z - x(i))
        end do
        f0 = 1.0D+00
        if ( nr == n0 .and. n /= 2 * int ( n / 2 ) ) then
           z = 0.0D+00
        end if
        f1 = z
        do k = 2, n
           pf = ( 2.0D+00 - 1.0D+00 / k ) * z * f1 &
                - ( 1.0D+00 - 1.0D+00 / k ) * f0
           pd = k * ( f1 - z * pf ) / ( 1.0D+00 - z * z )
           f0 = f1
           f1 = pf
        end do

        if ( z == 0.0D+00 ) then
           exit
        end if

        fd = pf / p
        q = 0.0D+00
        do i = 1, nr - 1
           wp = 1.0D+00
           do j = 1, nr - 1
              if ( j /= i ) then
                 wp = wp * ( z - x(j) )
              end if
           end do
           q = q + wp
        end do
        gd = ( pd - q * fd ) / p
        z = z - fd / gd

        if ( abs ( z - z0 ) < abs ( z ) * 1.0D-15 ) then
           exit
        end if

     end do

     x(nr) = z
     x(n+1-nr) = - z
     w(nr) = 2.0D+00 / ( ( 1.0D+00 - z * z ) * pd * pd )
     w(n+1-nr) = w(nr)

  end do

  return
end subroutine legzo




  !*****************************************************************************80
  !
  !! LGAMA computes the gamma function or its logarithm.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KF, the argument code.
  !    1, for gamma(x);
  !    2, for ln(gamma(x)).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) GL, the function value.
  !
subroutine lgama ( kf, x, gl )
  !LGAMA computes the gamma function or its logarithm.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param kf: in: the argument code. 1, for gamma(x); 2, for ln(gamma(x))
  !@param x: in: the argument
  !@param gl: out: the function value
  !
  implicit none

  real(8), save, dimension ( 10 ) :: a = (/ &
       8.333333333333333D-02, &
       -2.777777777777778D-03, &
       7.936507936507937D-04, &
       -5.952380952380952D-04, &
       8.417508417508418D-04, &
       -1.917526917526918D-03, &
       6.410256410256410D-03, &
       -2.955065359477124D-02, &
       1.796443723688307D-01, &
       -1.39243221690590D+00 /)
  real(8)    :: gl
  real(8)    :: gl0
  integer(4) :: k
  integer(4) :: kf
  integer(4) :: n
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x2
  real(8)    :: xp

  x0 = x

  if ( x == 1.0D+00 .or. x == 2.0D+00 ) then
     gl = 0.0D+00
     if ( kf == 1 ) then
        gl = 1.0D+00
     end if
     return
  else if ( x <= 7.0D+00 ) then
     n = int ( 7.0D+00 - x )
     x0 = x + n
  end if

  x2 = 1.0D+00 / ( x0 * x0 )
  xp = 6.283185307179586477D+00
  gl0 = a(10)

  do k = 9, 1, -1
     gl0 = gl0 * x2 + a(k)
  end do

  gl = gl0 / x0 + 0.5D+00 * log ( xp ) + ( x0 - 0.5D+00 ) * log ( x0 ) - x0

  if ( x <= 7.0D+00 ) then
     do k = 1, n
        gl = gl - log ( x0 - 1.0D+00 )
        x0 = x0 - 1.0D+00
     end do
  end if

  if ( kf == 1 ) then
     gl = exp ( gl )
  end if

  return
end subroutine lgama




  !*****************************************************************************80
  !
  !! LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    19 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) MM, the leading dimension of PM and PD.
  !
  !    Input, integer(4) M, the order of Pmn(x).
  !
  !    Input, integer(4) N, the degree of Pmn(x).
  !
  !    Input, real(8) X, the argument of Pmn(x).
  !
  !    Output, real(8) PM(0:MM,0:N), PD(0:MM,0:N), the
  !    values of Pmn(x) and Pmn'(x).
  !
subroutine lpmn ( mm, m, n, x, pm, pd )
  !LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param mm: in: the leading dimension of PM and PD
  !@param m: in: the order of Pmn(x)
  !@param n: in: the degree of Pmn(x)
  !@param x: in: the argument of Pmn(x)
  !@param pm: out: Pmn(x)
  !@param pd: out: Pmn'(x)
  !
  implicit none

  integer(4) :: mm
  integer(4) :: n

  integer(4) :: i
  integer(4) :: j
  integer(4) :: ls
  integer(4) :: m
  real(8)    :: pd(0:mm,0:n)
  real(8)    :: pm(0:mm,0:n)
  real(8)    :: x
  real(8)    :: xq
  real(8)    :: xs

  do i = 0, n
     do j = 0, m
        pm(j,i) = 0.0D+00
        pd(j,i) = 0.0D+00
     end do
  end do

  pm(0,0) = 1.0D+00

  if ( abs ( x ) == 1.0D+00 ) then

     do i = 1, n
        pm(0,i) = x ** i
        pd(0,i) = 0.5D+00 * i * ( i + 1.0D+00 ) * x ** ( i + 1 )
     end do

     do j = 1, n
        do i = 1, m
           if ( i == 1 ) then
              pd(i,j) = 1.0D+300
           else if ( i == 2 ) then
              pd(i,j) = -0.25D+00 * ( j + 2 ) * ( j + 1 ) * j &
                   * ( j - 1 ) * x ** ( j + 1 )
           end if
        end do
     end do

     return

  end if

  if ( 1.0D+00 < abs ( x ) ) then
     ls = -1
  else
     ls = +1
  end if

  xq = sqrt ( ls * ( 1.0D+00 - x * x ) )
  xs = ls * ( 1.0D+00 - x * x )
  do i = 1, m
     pm(i,i) = - ls * ( 2.0D+00 * i - 1.0D+00 ) * xq * pm(i-1,i-1)
  end do

  do i = 0, m
     pm(i,i+1) = ( 2.0D+00 * i + 1.0D+00 ) * x * pm(i,i)
  end do

  do i = 0, m
     do j = i + 2, n
        pm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * x * pm(i,j-1) - &
             ( i + j - 1.0D+00 ) * pm(i,j-2) ) / ( j - i )
     end do
  end do

  pd(0,0) = 0.0D+00
  do j = 1, n
     pd(0,j) = ls * j * ( pm(0,j-1) - x * pm(0,j) ) / xs
  end do

  do i = 1, m
     do j = i, n
        pd(i,j) = ls * i * x * pm(i,j) / xs + ( j + i ) &
             * ( j - i + 1.0D+00 ) / xq * pm(i-1,j)
     end do
  end do

  return
end subroutine lpmn




  !*****************************************************************************80
  !
  !! LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the order of Pmn(x).
  !
  !    Input, integer(4) N, the degree of Pmn(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PM(0:N), PD(0:N), the values and derivatives
  !    of the function from degree 0 to N.
  !
subroutine lpmns ( m, n, x, pm, pd )
  !LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the order of Pmn(x)
  !@param n: in: the degree of Pmn(x)
  !@param x: in: the argument
  !@param pm: out: the values of the function from degree 0 to N
  !@param pd: out: the derivatives of the function from degree 0 to N
  !
  implicit none

  integer(4) :: n

  integer(4) :: k
  integer(4) :: m
  real(8)    :: pm(0:n)
  real(8)    :: pm0
  real(8)    :: pm1
  real(8)    :: pm2
  real(8)    :: pmk
  real(8)    :: pd(0:n)
  real(8)    :: x
  real(8)    :: x0

  do k = 0, n
     pm(k) = 0.0D+00
     pd(k) = 0.0D+00
  end do

  if ( abs ( x ) == 1.0D+00 ) then

     do k = 0, n
        if ( m == 0 ) then
           pm(k) = 1.0D+00
           pd(k) = 0.5D+00 * k * ( k + 1.0D+00 )
           if ( x < 0.0D+00 ) then
              pm(k) = ( -1.0D+00 ) ** k * pm(k)
              pd(k) = ( -1.0D+00 ) ** ( k + 1 ) * pd(k)
           end if
        else if ( m == 1 ) then
           pd(k) = 1.0D+300
        else if ( m == 2 ) then
           pd(k) = -0.25D+00 * ( k + 2.0D+00 ) * ( k + 1.0D+00 ) &
                * k * ( k - 1.0D+00 )
           if ( x < 0.0D+00 ) then
              pd(k) = ( -1.0D+00 ) ** ( k + 1 ) * pd(k)
           end if
        end if
     end do
     return
  end if

  x0 = abs ( 1.0D+00 - x * x )
  pm0 = 1.0D+00
  pmk = pm0
  do k = 1, m
     pmk = ( 2.0D+00 * k - 1.0D+00 ) * sqrt ( x0 ) * pm0
     pm0 = pmk
  end do
  pm1 = ( 2.0D+00 * m + 1.0D+00 ) * x * pm0
  pm(m) = pmk
  pm(m+1) = pm1
  do k = m + 2, n
     pm2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * pm1 &
          - ( k + m - 1.0D+00 ) * pmk ) / ( k - m )
     pm(k) = pm2
     pmk = pm1
     pm1 = pm2
  end do

  pd(0) = ( ( 1.0D+00 - m ) * pm(1) - x * pm(0) ) &
       / ( x * x - 1.0D+00 )  
  do k = 1, n
     pd(k) = ( k * x * pm(k) - ( k + m ) * pm(k-1) ) &
          / ( x * x - 1.0D+00 )
  end do

  return
end subroutine lpmns




  !*****************************************************************************80
  !
  !! LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
  !
  !  Discussion:
  !
  !    Compute the associated Legendre function Pmv(x) with an integer order 
  !    and an arbitrary nonnegative degree v.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    19 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the degree of Pmv(x).
  !
  !    Input, integer(4) M, the order of Pmv(x).
  !
  !    Input, real(8) X, the argument of Pm(x).
  !
  !    Output, real(8) PMV, the value of Pm(x).
  !
subroutine lpmv ( v, m, x, pmv )
  !LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
  !
  !Compute the associated Legendre function Pmv(x) with an integer order
  !and an arbitrary nonnegative degree v.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the degree of Pmv(x)
  !@param m: in: the order of Pmv(x)
  !@param x: in: the argument of Pm(x)
  !@param pmv: out: the value of Pm(x)
  !
  implicit none

  real(8)    :: c0
  real(8)    :: el
  real(8)    :: eps
  integer(4) :: j
  integer(4) :: k
  integer(4) :: m
  integer(4) :: nv
  real(8)    :: pa
  real(8)    :: pi
  real(8)    :: pmv
  real(8)    :: pss
  real(8)    :: psv
  real(8)    :: pv0
  real(8)    :: qr
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: rg
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: s1
  real(8)    :: s2
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: vs
  real(8)    :: x
  real(8)    :: xq

  pi = 3.141592653589793D+00
  el = 0.5772156649015329D+00
  eps = 1.0D-14
  nv = int ( v )
  v0 = v - nv

  if ( x == -1.0D+00 .and. v /= nv ) then
     if ( m == 0 ) then
        pmv = -1.0D+300
     else
        pmv = 1.0D+300
     end if
     return
  end if

  c0 = 1.0D+00

  if ( m /= 0 ) then

     rg = v * ( v + m )
     do j = 1, m - 1 
        rg = rg * ( v * v - j * j )
     end do
     xq = sqrt ( 1.0D+00 - x * x )
     r0 = 1.0D+00
     do j = 1, m
        r0 = 0.5D+00 * r0 * xq / j
     end do
     c0 = r0 * rg

  end if

  if ( v0 == 0.0D+00 ) then

     pmv = 1.0D+00
     r = 1.0D+00
     do k = 1, nv - m
        r = 0.5D+00 * r * ( - nv + m + k - 1.0D+00 ) &
             * ( nv + m + k ) / ( k * ( k + m ) ) * ( 1.0D+00 + x )
        pmv = pmv + r
     end do
     pmv = ( -1.0D+00 ) ** nv * c0 * pmv

  else

     if ( -0.35D+00 <= x ) then

        pmv = 1.0D+00
        r = 1.0D+00
        do k = 1, 100
           r = 0.5D+00 * r * ( - v + m + k - 1.0D+00 ) &
                * ( v + m + k ) / ( k * ( m + k ) ) * ( 1.0D+00 - x )
           pmv = pmv + r
           if ( 12 < k .and. abs ( r / pmv ) < eps ) then
              exit
           end if
        end do

        pmv = ( -1.0D+00 ) ** m * c0 * pmv

     else

        vs = sin ( v * pi ) / pi
        pv0 = 0.0D+00

        if ( m /= 0 ) then

           qr = sqrt ( ( 1.0D+00 - x ) / ( 1.0D+00 + x ) )
           r2 = 1.0D+00
           do j = 1, m
              r2 = r2 * qr * j
           end do
           s0 = 1.0D+00
           r1 = 1.0D+00
           do k = 1, m - 1 
              r1 = 0.5D+00 * r1 * ( - v + k - 1 ) * ( v + k ) &
                   / ( k * ( k - m ) ) * ( 1.0D+00 + x )
              s0 = s0 + r1
           end do
           pv0 = - vs * r2 / m * s0

        end if

        call psi ( v, psv )
        pa = 2.0D+00 * ( psv + el ) + pi / tan ( pi * v ) &
             + 1.0D+00 / v

        s1 = 0.0D+00
        do j = 1, m
           s1 = s1 + ( j * j + v * v ) / ( j * ( j * j - v * v ) )
        end do

        pmv = pa + s1 - 1.0D+00 / ( m - v ) &
             + log ( 0.5D+00 * ( 1.0D+00 + x ) )
        r = 1.0D+00
        do k = 1, 100
           r = 0.5D+00 * r * ( - v + m + k - 1.0D+00 ) * ( v + m + k ) &
                / ( k * ( k + m ) ) * ( 1.0D+00 + x )
           s = 0.0D+00
           do j = 1, m
              s = s + ( ( k + j ) ** 2 + v * v ) &
                   / ( ( k + j ) * ( ( k + j ) ** 2 - v * v ) )
           end do
           s2 = 0.0D+00
           do j = 1, k
              s2 = s2 + 1.0D+00 / ( j * ( j * j - v * v ) )
           end do
           pss = pa + s + 2.0D+00 * v * v * s2 &
                - 1.0D+00 / ( m + k - v ) &
                + log ( 0.5D+00 * ( 1.0D+00 + x ) )
           r2 = pss * r
           pmv = pmv + r2
           if ( abs ( r2 / pmv ) < eps ) then
              exit
           end if
        end do

        pmv = pv0 + pmv * vs * c0

     end if

  end if

  return
end subroutine lpmv




  !*****************************************************************************80
  !
  !! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the maximum degree.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PN(0:N), PD(0:N), the values and derivatives
  !    of the polyomials of degrees 0 to N at X.
  !
subroutine lpn ( n, x, pn, pd )
  !LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the maximum degree
  !@param x: in: the argument
  !@param pn: out: the values of the polynomials of degrees 0 to N at X
  !@param pd: out: the derivatives of the polynomials of degrees 0 to N at X
  !
  implicit none

  integer(4) :: n

  integer(4) :: k
  real(8)    :: p0
  real(8)    :: p1
  real(8)    :: pd(0:n)
  real(8)    :: pf
  real(8)    :: pn(0:n)
  real(8)    :: x

  pn(0) = 1.0D+00
  pn(1) = x
  pd(0) = 0.0D+00
  pd(1) = 1.0D+00
  p0 = 1.0D+00
  p1 = x

  do k = 2, n

     pf = ( 2.0D+00 * k - 1.0D+00 ) / k * x * p1 &
          - ( k - 1.0D+00 ) / k * p0
     pn(k) = pf

     if ( abs ( x ) == 1.0D+00 ) then
        pd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
     else
        pd(k) = k * ( p1 - x * pf ) / ( 1.0D+00 - x * x )
     end if

     p0 = p1
     p1 = pf

  end do

  return
end subroutine lpn




  !*****************************************************************************80
  !
  !! LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
  !
  !  Discussion:
  !
  !    This routine computes Legendre polynomials Pn(x), Pn'(x)
  !    and the integral of Pn(t) from 0 to x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    13 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the maximum degree.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PN(0:N), PD(0:N), PL(0:N), the values, 
  !    derivatives and integrals of the polyomials of degrees 0 to N at X.
  !
subroutine lpni ( n, x, pn, pd, pl )
  !LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
  !
  !This routine computes Legendre polynomials Pn(x), Pn'(x)
  !and the integral of Pn(t) from 0 to x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the maximum degree
  !@param x: in: the argument
  !@param pn: out: the values of the polynomials of degrees 0 to N at X
  !@param pd: out: the derivatives of the polynomials of degrees 0 to N at X
  !@param pl: out: the integrals of the polynomials of degrees 0 to N at X
  !
  implicit none

  integer(4) :: n

  integer(4) :: j
  integer(4) :: k
  integer(4) :: n1
  real(8)    :: p0
  real(8)    :: p1
  real(8)    :: pd(0:n)
  real(8)    :: pf
  real(8)    :: pl(0:n)
  real(8)    :: pn(0:n)
  real(8)    :: r
  real(8)    :: x

  pn(0) = 1.0D+00
  pn(1) = x
  pd(0) = 0.0D+00
  pd(1) = 1.0D+00
  pl(0) = x
  pl(1) = 0.5D+00 * x * x
  p0 = 1.0D+00
  p1 = x

  do k = 2, n

     pf = ( 2.0D+00 * k - 1.0D+00 ) / k * x * p1 - ( k - 1.0D+00 ) / k * p0
     pn(k) = pf

     if ( abs ( x ) == 1.0D+00 ) then
        pd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
     else
        pd(k) = k * ( p1 - x * pf ) / ( 1.0D+00 - x * x )
     end if

     pl(k) = ( x * pn(k) - pn(k-1) ) / ( k + 1.0D+00 )
     p0 = p1
     p1 = pf

     if ( k /= 2 * int ( k / 2 ) ) then

        r = 1.0D+00 / ( k + 1.0D+00 )
        n1 = ( k - 1 ) / 2
        do j = 1, n1
           r = ( 0.5D+00 / j - 1.0D+00 ) * r
        end do
        pl(k) = pl(k) + r

     end if

  end do

  return
end subroutine lpni




  !*****************************************************************************80
  !
  !! LQMN computes associated Legendre functions Qmn(x) and derivatives.
  !
  !  Discussion:
  !
  !    This routine computes the associated Legendre functions of the
  !    second kind, Qmn(x) and Qmn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    13 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) MM, determines the leading dimension 
  !    of QM and QD.
  !
  !    Input, integer(4) M, the order of Qmn(x).
  !
  !    Input, integer(4) N, the degree of Qmn(x).
  !
  !    Output, real(8) QM(0:MM,0:N), QD(0:MM,0:N), contains the values
  !    of Qmn(x) and Qmn'(x).
  !
subroutine lqmn ( mm, m, n, x, qm, qd )
  !LQMN computes associated Legendre functions Qmn(x) and derivatives.
  !
  !This routine computes the associated Legendre functions of the
  !second kind, Qmn(x) and Qmn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param mm: in: determines the leading dimension of QM and QD
  !@param m: in: the order of Qmn(x)
  !@param n: in: the degree of Qmn(x)
  !@param x: the argument
  !@param qm: out: the values of Qmn(x)
  !@param qd: out: the values of Qmn'(x)
  !
  implicit none

  integer(4) :: mm
  integer(4) :: n

  integer(4) :: i
  integer(4) :: j
  integer(4) :: k
  integer(4) :: km
  integer(4) :: ls
  integer(4) :: m
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: q10
  real(8)    :: qd(0:mm,0:n)
  real(8)    :: qf
  real(8)    :: qf0
  real(8)    :: qf1
  real(8)    :: qf2
  real(8)    :: qm(0:mm,0:n)
  real(8)    :: x
  real(8)    :: xq
  real(8)    :: xs

  if ( abs ( x ) == 1.0D+00 ) then
     do i = 0, m
        do j = 0, n
           qm(i,j) = 1.0D+300
           qd(i,j) = 1.0D+300
        end do
     end do
     return
  end if

  if ( 1.0D+00 < abs ( x ) ) then
     ls = -1
  else
     ls = 1
  end if

  xs = ls * ( 1.0D+00 - x * x )
  xq = sqrt ( xs )
  q0 = 0.5D+00 * log ( abs ( ( x + 1.0D+00 ) / ( x - 1.0D+00 ) ) )

  if ( abs ( x ) < 1.0001D+00 ) then
     qm(0,0) = q0
     qm(0,1) = x * q0 - 1.0D+00
     qm(1,0) = -1.0D+00 / xq
     qm(1,1) = -xq * ( q0 + x / ( 1.0D+00 - x * x ) )
     do i = 0, 1
        do j = 2, n
           qm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * x * qm(i,j-1) &
                - ( j + i - 1.0D+00 ) * qm(i,j-2))/ ( j - i )
        end do
     end do

     do j = 0, n
        do i = 2, m
           qm(i,j) = -2.0D+00 * ( i - 1.0D+00 ) * x / xq * qm(i-1,j) &
                - ls * ( j + i - 1.0D+00 ) * ( j - i + 2.0D+00 ) * qm(i-2,j)
        end do
     end do

  else

     if ( 1.1D+00 < abs ( x ) ) then
        km = 40 + m + n
     else
        km = ( 40 + m + n ) &
             * int ( -1.0D+00 - 1.8D+00 * log ( x - 1.0D+00 ) )
     end if

     qf2 = 0.0D+00
     qf1 = 1.0D+00
     do k = km, 0, -1
        qf0 = ( ( 2 * k + 3.0D+00 ) * x * qf1 &
             - ( k + 2.0D+00 ) * qf2 ) / ( k + 1.0D+00 )
        if ( k <= n ) then
           qm(0,k) = qf0
        end if
        qf2 = qf1
        qf1 = qf0
     end do

     do k = 0, n
        qm(0,k) = q0 * qm(0,k) / qf0
     end do

     qf2 = 0.0D+00
     qf1 = 1.0D+00
     do k = km, 0, -1
        qf0 = ( ( 2 * k + 3.0D+00 ) * x * qf1 &
             - ( k + 1.0D+00 ) * qf2 ) / ( k + 2.0D+00 )
        if ( k <= n ) then
           qm(1,k) = qf0
        end if
        qf2 = qf1
        qf1 = qf0
     end do

     q10 = -1.0D+00 / xq
     do k = 0, n
        qm(1,k) = q10 * qm(1,k) / qf0
     end do

     do j = 0, n
        q0 = qm(0,j)
        q1 = qm(1,j)
        do i = 0, m - 2
           qf = -2.0D+00 * ( i + 1 ) * x / xq * q1 &
                + ( j - i ) * ( j + i + 1.0D+00 ) * q0
           qm(i+2,j) = qf
           q0 = q1
           q1 = qf
        end do
     end do

  end if

  qd(0,0) = ls / xs
  do j = 1, n
     qd(0,j) = ls * j * ( qm(0,j-1) - x * qm(0,j) ) / xs
  end do

  do j = 0, n
     do i = 1, m
        qd(i,j) = ls * i * x / xs * qm(i,j) &
             + ( i + j ) * ( j - i + 1.0D+00 ) / xq * qm(i-1,j)
     end do
  end do

  return
end subroutine lqmn




  !*****************************************************************************80
  !
  !! LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the order.
  !
  !    Input, integer(4) N, the degree.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) QM(0:N), QD(0:N), the values of Qmn(x) 
  !    and Qmn'(x).
  !
subroutine lqmns ( m, n, x, qm, qd )
  !LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the order
  !@param n: in: the degree
  !@param x: in: the argument
  !@param qm: out: Qmn(x)
  !@param qd: out: Qmn'(x)
  !
  implicit none

  integer(4) :: n

  integer(4) :: k
  integer(4) :: km
  integer(4) :: l
  integer(4) :: ls
  integer(4) :: m
  real(8)    :: q0
  real(8)    :: q00
  real(8)    :: q01
  real(8)    :: q0l
  real(8)    :: q10
  real(8)    :: q11
  real(8)    :: q1l
  real(8)    :: qd(0:n)
  real(8)    :: qf0
  real(8)    :: qf1
  real(8)    :: qf2
  real(8)    :: qg0
  real(8)    :: qg1
  real(8)    :: qh0
  real(8)    :: qh1
  real(8)    :: qh2
  real(8)    :: qm(0:n)
  real(8)    :: qm0
  real(8)    :: qm1
  real(8)    :: qmk
  real(8)    :: x
  real(8)    :: xq

  do k = 0, n
     qm(k) = 0.0D+00
     qd(k) = 0.0D+00
  end do

  if ( abs ( x ) == 1.0D+00 ) then
     do k = 0, n
        qm(k) = 1.0D+300
        qd(k) = 1.0D+300
     end do
     return
  end if

  if ( 1.0D+00 < abs ( x ) ) then
     ls = -1
  else
     ls = +1
  end if

  xq = sqrt ( ls * ( 1.0D+00 - x * x ) )
  q0 = 0.5D+00 * log ( abs ( ( x + 1.0D+00 ) / ( x - 1.0D+00 ) ) )
  q00 = q0
  q10 = -1.0D+00 / xq
  q01 = x * q0 - 1.0D+00
  q11 = - ls * xq * ( q0 + x / ( 1.0D+00 - x * x ) )
  qf0 = q00
  qf1 = q10
  do k = 2, m
     qm0 = -2.0D+00 * ( k - 1.0D+00 ) / xq * x * qf1 &
          - ls * ( k - 1.0D+00 ) * ( 2.0D+00 - k ) * qf0
     qf0 = qf1
     qf1 = qm0
  end do

  if ( m == 0 ) then
     qm0 = q00
  else if ( m == 1 ) then
     qm0 = q10
  end if

  qm(0) = qm0

  if ( abs ( x ) < 1.0001D+00 ) then

     if ( m == 0 .and. 0 < n ) then

        qf0 = q00
        qf1 = q01
        do k = 2, n
           qf2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qf1 &
                - ( k - 1.0D+00 ) * qf0 ) / k
           qm(k) = qf2
           qf0 = qf1
           qf1 = qf2
        end do

     end if
     qg0 = q01
     qg1 = q11
     do k = 2, m
        qm1 = - 2.0D+00 * ( k - 1.0D+00 ) / xq * x * qg1 &
             - ls * k * ( 3.0D+00 - k ) * qg0
        qg0 = qg1
        qg1 = qm1
     end do

     if ( m == 0 ) then
        qm1 = q01
     else if ( m == 1 ) then
        qm1 = q11
     end if
     qm(1) = qm1

     if ( m == 1 .and. 1 < n ) then

        qh0 = q10
        qh1 = q11
        do k = 2, n
           qh2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qh1 - k * qh0 ) &
                / ( k - 1.0D+00 )
           qm(k) = qh2
           qh0 = qh1
           qh1 = qh2
        end do

     else if ( 2 <= m ) then

        qg0 = q00
        qg1 = q01
        qh0 = q10
        qh1 = q11

        do l = 2, n
           q0l = ( ( 2.0D+00 * l - 1.0D+00 ) * x * qg1 &
                - ( l - 1.0D+00 ) * qg0 ) / l
           q1l = ( ( 2.0D+00 * l - 1.0D+00 ) * x * qh1 - l * qh0 ) &
                / ( l - 1.0D+00 )
           qf0 = q0l
           qf1 = q1l
           do k = 2, m
              qmk = - 2.0D+00 * ( k - 1.0D+00 ) / xq * x * qf1 &
                   - ls * ( k + l - 1.0D+00 ) * ( l + 2.0D+00 - k ) * qf0
              qf0 = qf1
              qf1 = qmk
           end do
           qm(l) = qmk
           qg0 = qg1
           qg1 = q0l
           qh0 = qh1
           qh1 = q1l
        end do

     end if

  else

     if ( 1.1D+00 < abs ( x ) ) then
        km = 40 + m + n
     else
        km = ( 40 + m + n ) * int ( - 1.0D+00 - 1.8D+00 * log ( x - 1.0D+00 ) )
     end if

     qf2 = 0.0D+00
     qf1 = 1.0D+00
     do k = km, 0, -1
        qf0 = ( ( 2.0D+00 * k + 3.0D+00 ) * x * qf1 &
             - ( k + 2.0D+00 - m ) * qf2 ) / ( k + m + 1.0D+00 )
        if ( k <= n ) then
           qm(k) = qf0
        end if
        qf2 = qf1
        qf1 = qf0
     end do

     do k = 0, n
        qm(k) = qm(k) * qm0 / qf0
     end do

  end if

  if ( abs ( x ) < 1.0D+00 ) then
     do k = 0, n
        qm(k) = ( -1 ) ** m * qm(k)
     end do
  end if

  qd(0) = ( ( 1.0D+00 - m ) * qm(1) - x * qm(0) )  / ( x * x - 1.0D+00 )
  do k = 1, n
     qd(k) = ( k * x * qm(k) - ( k + m ) * qm(k-1) ) / ( x * x - 1.0D+00 )
  end do

  return
end subroutine lqmns




  !*****************************************************************************80
  !
  !! LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    19 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the degree of Qn(x).
  !
  !    Input, real(8) X, the argument of Qn(x).
  !
  !    Output, real(8) QN(0:N), QD(0:N), the values of
  !    Qn(x) and Qn'(x).
  !
subroutine lqna ( n, x, qn, qd )
  !LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the degree of Qn(x)
  !@param x: in: the argument of Qn(x)
  !@param qn: out: Qn(x)
  !@param qd: out: Qn'(x)
  !
  implicit none

  integer(4) :: n

  integer(4) :: k
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: qd(0:n)
  real(8)    :: qf
  real(8)    :: qn(0:n)
  real(8)    :: x

  if ( abs ( x ) == 1.0D+00 ) then

     do k = 0, n
        qn(k) = 1.0D+300
        qd(k) = -1.0D+300
     end do

  else if ( abs ( x ) < 1.0D+00 ) then

     q0 = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )
     q1 = x * q0 - 1.0D+00
     qn(0) = q0
     qn(1) = q1
     qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
     qd(1) = qn(0) + x * qd(0)
     do k = 2, n
        qf = ( ( 2 * k - 1 ) * x * q1 - ( k - 1 ) * q0 ) / k
        qn(k) = qf
        qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0D+00 - x * x )
        q0 = q1
        q1 = qf
     end do

  end if

  return
end subroutine lqna




  !*****************************************************************************80
  !
  !! LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    19 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the degree of Qn(x).
  !
  !    Input, real(8) X, the argument of Qn(x).
  !
  !    Output, real(8) QN(0:N), QD(0:N), the values of
  !    Qn(x) and Qn'(x).
  !
subroutine lqnb ( n, x, qn, qd )
  !LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the degree of Qn(x)
  !@param x: in: the argument of Qn(x)
  !@param qn: out: Qn(x)
  !@param qd: out: Qn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: eps
  integer(4) :: j
  integer(4) :: k
  integer(4) :: l
  integer(4) :: nl
  real(8)    :: q0
  real(8)    :: q1
  real(8)    :: qc1
  real(8)    :: qc2
  real(8)    :: qd(0:n)
  real(8)    :: qf
  real(8)    :: qf0
  real(8)    :: qf1
  real(8)    :: qf2
  real(8)    :: qn(0:n)
  real(8)    :: qr
  real(8)    :: x
  real(8)    :: x2

  eps = 1.0D-14

  if ( abs ( x ) == 1.0D+00 ) then
     do k = 0, n
        qn(k) = 1.0D+300
        qd(k) = 1.0D+300
     end do
     return
  end if

  if ( x <= 1.021D+00 ) then

     x2 = abs ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )
     q0 = 0.5D+00 * log ( x2 )
     q1 = x * q0 - 1.0D+00
     qn(0) = q0
     qn(1) = q1
     qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
     qd(1) = qn(0) + x * qd(0)
     do k = 2, n
        qf = ( ( 2.0D+00 * k - 1.0D+00 ) * x * q1 &
             - ( k - 1.0D+00 ) * q0 ) / k
        qn(k) = qf
        qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0D+00 - x * x )
        q0 = q1
        q1 = qf
     end do

  else

     qc2 = 1.0D+00 / x
     do j = 1, n
        qc2 = qc2 * j / ( ( 2.0D+00 * j + 1.0D+00 ) * x )
        if ( j == n - 1 ) then
           qc1 = qc2
        end if
     end do

     do l = 0, 1

        nl = n + l
        qf = 1.0D+00
        qr = 1.0D+00
        do k = 1, 500
           qr = qr * ( 0.5D+00 * nl + k - 1.0D+00 ) &
                * ( 0.5D+00 * ( nl - 1 ) + k ) &
                / ( ( nl + k - 0.5D+00 ) * k * x * x )
           qf = qf + qr
           if ( abs ( qr / qf ) < eps ) then
              exit
           end if
        end do

        if ( l == 0 ) then
           qn(n-1) = qf * qc1
        else
           qn(n) = qf * qc2
        end if

     end do

     qf2 = qn(n)
     qf1 = qn(n-1)
     do k = n, 2, -1
        qf0 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qf1 - k * qf2 ) / ( k - 1.0D+00 )
        qn(k-2) = qf0
        qf2 = qf1
        qf1 = qf0
     end do

     qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
     do k = 1, n
        qd(k) = k * ( qn(k-1) - x * qn(k) ) / ( 1.0D+00 - x * x )
     end do

  end if

  return
end subroutine lqnb




  !*****************************************************************************80
  !
  !! MSTA1 determines a backward recurrence starting point for Jn(x).
  !
  !  Discussion:
  !
  !    This procedure determines the starting point for backward  
  !    recurrence such that the magnitude of    
  !    Jn(x) at that point is about 10^(-MP).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    08 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, integer(4) MP, the negative logarithm of the 
  !    desired magnitude.
  !
  !    Output, integer(4) MSTA1, the starting point.
  !
function msta1 ( x, mp )
  !MSTA1 determines a backward recurrence starting point for Jn(x).
  !
  !This procedure determines the starting point for backward
  !recurrence such that the magnitude of
  !Jn(x) at that point is about 10^(-MP).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param mp: in: the negative logarithm of the desired magnitude
  !@param msta1: out: the starting point
  !
  implicit none

  real(8)    :: a0
  ! real(8) envj
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: it
  integer(4) :: mp
  integer(4) :: msta1
  integer(4) :: n0
  integer(4) :: n1
  integer(4) :: nn
  real(8)    :: x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
     nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
     f = envj ( nn, a0 ) - mp
     if ( abs ( nn - n1 ) < 1 ) then
        exit
     end if
     n0 = n1
     f0 = f1
     n1 = nn
     f1 = f
  end do

  msta1 = nn

  return
end function msta1




  !*****************************************************************************80
  !
  !! MSTA2 determines a backward recurrence starting point for Jn(x).
  !
  !  Discussion:
  !
  !    This procedure determines the starting point for a backward
  !    recurrence such that all Jn(x) has MP significant digits.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    08 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument of Jn(x).
  !
  !    Input, integer(4) N, the order of Jn(x).
  !
  !    Input, integer(4) MP, the number of significant digits.
  !
  !    Output, integer(4) MSTA2, the starting point.
  !
function msta2 ( x, n, mp )
  !MSTA2 determines a backward recurrence starting point for Jn(x).
  !
  !This procedure determines the starting point for a backward
  !recurrence such that all Jn(x) has MP significant digits.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument of Jn(x)
  !@param n: in: the order of Jn(x)
  !@param mp: in: the number of significant digits
  !@param msta2: out: the starting point
  !
  implicit none

  real(8)    :: a0
  real(8)    :: ejn
  ! real(8) envj
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: hmp
  integer(4) :: it
  integer(4) :: mp
  integer(4) :: msta2
  integer(4) :: n
  integer(4) :: n0
  integer(4) :: n1
  integer(4) :: nn
  real(8)    :: obj
  real(8)    :: x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
     obj = mp
     n0 = int ( 1.1D+00 * a0 )
  else
     obj = hmp + ejn
     n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
     nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
     f = envj ( nn, a0 ) - obj
     if ( abs ( nn - n1 ) < 1 ) then
        exit
     end if
     n0 = n1
     f0 = f1
     n1 = nn
     f1 = f
  end do

  msta2 = nn + 10

  return
end function msta2




  !*****************************************************************************80
  !
  !! MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    20 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KF, the function code.
  !    1 for computing cem(x,q) and cem'(x,q)
  !    2 for computing sem(x,q) and sem'(x,q).
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Input, real(8) X, the argument of the Mathieu functions,
  !    in degrees.
  !
  !    Output, real(8) CSF, CSD, the values of cem(x,q) and cem'(x,q),
  !    or of sem(x,q) and sem'(x,q).
  !
subroutine mtu0 ( kf, m, q, x, csf, csd )
  !MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KF: the function code. 1 for computing cem(x,q) and cem'(x,q) 2 for computing sem(x,q) and sem'(x,q)
  !
  !@param kf: in: function code: 1 for cem(x,q) and cem'(x,q), 2 for sem(x,q) and sem'(x,q)
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param x: in: the argument of the Mathieu functions, in degrees
  !@param csf: out: cem(x,q) or sem(x,q)
  !@param csd: out: cem'(x,q) or sem'(x,q)
  !
  implicit none

  real(8)    :: a
  real(8)    :: csd
  real(8)    :: csf
  real(8)    :: eps
  real(8)    :: fg(251)
  integer(4) :: ic
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: kf
  integer(4) :: km
  integer(4) :: m
  real(8)    :: q
  real(8)    :: qm
  real(8)    :: rd
  real(8)    :: x
  real(8)    :: xr

  eps = 1.0D-14

  if ( kf == 1 ) then

     if ( m == 2 * int ( m / 2 ) ) then
        kd = 1
     else
        kd = 2
     end if

  else

     if ( m /= 2 * int ( m / 2 ) ) then
        kd = 3
     else
        kd = 4
     end if

  end if

  call cva2 ( kd, m, q, a )

  if ( q <= 1.0D+00 ) then
     qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q &
          + 90.7D+00 * sqrt ( q ) * q
  else
     qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q &
          + 0.0037D+00 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5D+00 * m )
  call fcoef ( kd, m, q, a, fg )
  ic = int ( m / 2 ) + 1
  rd = 1.74532925199433D-02
  xr = x * rd

  csf = 0.0D+00

  do k = 1, km

     if ( kd == 1 ) then
        csf = csf + fg(k) * cos ( ( 2.0D+00 * k - 2.0D+00 ) * xr )
     else if ( kd == 2 ) then
        csf = csf + fg(k) * cos ( ( 2.0D+00 * k - 1.0D+00 ) * xr )
     else if ( kd == 3 ) then
        csf = csf + fg(k) * sin ( ( 2.0D+00 * k - 1.0D+00 ) * xr )
     else if ( kd == 4 ) then
        csf = csf + fg(k) * sin ( 2.0D+00 * k * xr )
     end if

     if ( ic <= k .and. abs ( fg(k) ) < abs ( csf ) * eps ) then
        exit
     end if

  end do

  csd = 0.0D+00

  do k = 1, km

     if ( kd == 1 ) then
        csd = csd - ( 2 * k - 2 ) * fg(k) * sin ( ( 2 * k - 2 ) * xr )
     else if ( kd == 2 ) then
        csd = csd - ( 2 * k - 1 ) * fg(k) * sin ( ( 2 * k - 1 ) * xr )
     else if ( kd == 3 ) then
        csd = csd + ( 2 * k - 1 ) * fg(k) * cos ( ( 2 * k - 1 ) * xr )
     else if ( kd == 4 ) then
        csd = csd + 2.0D+00 * k * fg(k) * cos ( 2 * k * xr )
     end if

     if ( ic <= k .and. abs ( fg(k) ) < abs ( csd ) * eps ) then
        exit
     end if

  end do

  return
end subroutine mtu0




  !*****************************************************************************80
  !
  !! MTU12 computes modified Mathieu functions of the first and second kind.
  !
  !  Discussion:
  !
  !    This procedure computes modified Mathieu functions of the first and
  !    second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
  !    and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    31 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KF, the function code.
  !    1 for computing Mcm(x,q);
  !    2 for computing Msm(x,q).
  !
  !    Input, integer(4) KC, the function code.
  !    1, for computing the first kind
  !    2, for computing the second kind or Msm(2)(x,q) and Msm(2)'(x,q)
  !    3, for computing both the first and second kinds.
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Input, real(8) X, the argument of the Mathieu functions.
  !
  !    Output, real(8) F1R, D1R, F2R, D2R, the values of 
  !    Mcm(1)(x,q) or Msm(1)(x,q), Derivative of Mcm(1)(x,q) or Msm(1)(x,q),
  !    Mcm(2)(x,q) or Msm(2)(x,q), Derivative of Mcm(2)(x,q) or Msm(2)(x,q).
  !
subroutine mtu12 ( kf, kc, m, q, x, f1r, d1r, f2r, d2r )
  !MTU12 computes modified Mathieu functions of the first and second kind.
  !
  !This procedure computes modified Mathieu functions of the first and
  !second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
  !and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KC: the function code. 1, for computing the first kind 2, for computing the second kind or Msm(2)(x,q) and
  !  Msm(2)'(x,q) 3, for computing both the first and second kinds
  !* F1R, D1R, F2R, D2R: the values of Mcm(1)(x,q) or Msm(1)(x,q), Derivative of Mcm(1)(x,q) or Msm(1)(x,q), Mcm(2)(x,q)
  !  or Msm(2)(x,q), Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
  !
  !@param kf: in: the function code. 1 for computing Mcm(x,q); 2 for computing Msm(x,q)
  !@param kc: in: function code: 1 first kind, 2 second kind, 3 both
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param x: in: the argument of the Mathieu functions
  !@param f1r: out: Mcm(1)(x,q) or Msm(1)(x,q)
  !@param d1r: out: derivative of Mcm(1)(x,q) or Msm(1)(x,q)
  !@param f2r: out: Mcm(2)(x,q) or Msm(2)(x,q)
  !@param d2r: out: derivative of Mcm(2)(x,q) or Msm(2)(x,q)
  !
  implicit none

  real(8)    :: a
  real(8)    :: bj1(0:251)
  real(8)    :: bj2(0:251)
  real(8)    :: by1(0:251)
  real(8)    :: by2(0:251)
  real(8)    :: c1
  real(8)    :: c2
  real(8)    :: d1r
  real(8)    :: d2r
  real(8)    :: dj1(0:251)
  real(8)    :: dj2(0:251)
  real(8)    :: dy1(0:251)
  real(8)    :: dy2(0:251)
  real(8)    :: eps
  real(8)    :: f1r
  real(8)    :: f2r
  real(8)    :: fg(251)
  integer(4) :: ic
  integer(4) :: k
  integer(4) :: kc
  integer(4) :: kd
  integer(4) :: kf
  integer(4) :: km
  integer(4) :: m
  integer(4) :: nm
  real(8)    :: q
  real(8)    :: qm
  real(8)    :: u1
  real(8)    :: u2
  real(8)    :: w1
  real(8)    :: w2
  real(8)    :: x

  eps = 1.0D-14

  if ( kf == 1 ) then
     if ( m == 2 * int ( m / 2 ) ) then
        kd = 1
     else
        kd = 2
     end if
  else
     if ( m /= 2 * int ( m / 2 ) ) then
        kd = 3
     else
        kd = 4
     end if
  end if

  call cva2 ( kd, m, q, a )

  if ( q <= 1.0D+00 ) then
     qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q &
          + 90.7D+00 * sqrt ( q ) * q
  else
     qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q &
          + 0.0037D+00 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5D+00 * m )              
  call fcoef ( kd, m, q, a, fg )

  if ( kd == 4 ) then
     ic = m / 2
  else
     ic = int ( m / 2 ) + 1
  end if

  c1 = exp ( - x )
  c2 = exp ( x )
  u1 = sqrt ( q ) * c1
  u2 = sqrt ( q ) * c2

  call jynb ( km, u1, nm, bj1, dj1, by1, dy1 )
  call jynb ( km, u2, nm, bj2, dj2, by2, dy2 )

  if ( kc == 1 ) then

     f1r = 0.0D+00

     do k = 1, km

        if ( kd == 1 ) then
           f1r = f1r + ( - 1.0D+00 ) ** ( ic + k ) * fg(k) * bj1(k-1) * bj2(k-1)
        else if ( kd == 2 .or. kd == 3 ) then
           f1r = f1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) * ( bj1(k-1) * bj2(k) &
                + ( - 1.0D+00 ) ** kd * bj1(k) * bj2(k-1) )
        else
           f1r = f1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( bj1(k-1) * bj2(k+1) - bj1(k+1) * bj2(k-1) )
        end if

        if ( 5 <= k .and. abs ( f1r - w1 ) < abs ( f1r ) * eps ) then
           exit
        end if

        w1 = f1r

     end do

     f1r = f1r / fg(1)
     d1r = 0.0D+00
     do k = 1, km
        if ( kd == 1 ) then
           d1r = d1r + ( - 1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * bj1(k-1) * dj2(k-1) - c1 * dj1(k-1) * bj2(k-1) )
        else if ( kd == 2 .or. kd == 3 ) then
           d1r = d1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * ( bj1(k-1) * dj2(k) &
                + ( -1.0D+00 ) ** kd * bj1(k) * dj2(k-1) ) &
                - c1 * ( dj1(k-1) * bj2(k) &
                + ( -1.0D+00 ) ** kd * dj1(k) * bj2(k-1) ) )
        else
           d1r = d1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * ( bj1(k-1) * dj2(k+1) - bj1(k+1) * dj2(k-1) ) &
                - c1 * ( dj1(k-1) * bj2(k+1) - dj1(k+1) * bj2(k-1) ) )
        end if
        if ( 5 <= k .and. abs ( d1r - w2 ) < abs ( d1r ) * eps ) then
           exit
        end if
        w2 = d1r
     end do

     d1r = d1r * sqrt ( q ) / fg(1)

  else

     f2r = 0.0D+00

     do k = 1, km
        if ( kd == 1 ) then
           f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * bj1(k-1) * by2(k-1)
        else if ( kd == 2 .or. kd == 3 ) then
           f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) * ( bj1(k-1) * by2(k) &
                + ( -1.0D+00 ) ** kd * bj1(k) * by2(k-1) )
        else
           f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( bj1(k-1) * by2(k+1) - bj1(k+1) * by2(k-1) )
        end if
        if ( 5 <= k .and. abs ( f2r - w1 ) < abs ( f2r ) * eps ) then
           exit
        end if
        w1 = f2r
     end do

     f2r = f2r / fg(1)
     d2r = 0.0D+00

     do k = 1, km
        if ( kd == 1 ) then
           d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * bj1(k-1) * dy2(k-1) - c1 * dj1(k-1) * by2(k-1) )
        else if ( kd == 2 .or. kd == 3 ) then
           d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * ( bj1(k-1) * dy2(k) &
                + ( -1.0D+00 ) ** kd * bj1(k) * dy2(k-1) ) &
                - c1 * ( dj1(k-1) * by2(k) + ( -1.0D+00 ) ** kd &
                * dj1(k) * by2(k-1) ) )
        else
           d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) &
                * ( c2 * ( bj1(k-1) * dy2(k+1) - bj1(k+1) * dy2(k-1) ) &
                - c1 * ( dj1(k-1) * by2(k+1) - dj1(k+1) * by2(k-1) ) )
        end if

        if ( 5 <= k .and. abs ( d2r - w2 ) < abs ( d2r ) * eps ) then
           exit
        end if

        w2 = d2r

     end do

     d2r = d2r * sqrt ( q ) / fg(1)

  end if

  return
end subroutine mtu12




  !*****************************************************************************80
  !
  !! OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
  !
  !  Discussion:
  !
  !    This procedure computes orthogonal polynomials: Tn(x) or Un(x),
  !    or Ln(x) or Hn(x), and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    08 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KT, the function code:
  !    1 for Chebyshev polynomial Tn(x)
  !    2 for Chebyshev polynomial Un(x)
  !    3 for Laguerre polynomial Ln(x)
  !    4 for Hermite polynomial Hn(x)
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PL(0:N), DPL(0:N), the value and derivative of
  !    the polynomials of order 0 through N at X.
  !
subroutine othpl ( kf, n, x, pl, dpl )
  !OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
  !
  !This procedure computes orthogonal polynomials: Tn(x) or Un(x),
  !or Ln(x) or Hn(x), and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param kf: function code: 1 Chebyshev Tn(x), 2 Chebyshev Un(x), 3 Laguerre Ln(x), 4 Hermite Hn(x)
  !@param n: in: the order
  !@param x: in: the argument
  !@param pl: out: the value of the polynomials of order 0 through N at X
  !@param dpl: out: the derivative of the polynomials of order 0 through N at X
  !
  implicit none

  integer    :: n

  real(8)    :: a
  real(8)    :: b
  real(8)    :: c
  real(8)    :: dpl(0:n)
  real(8)    :: dy0
  real(8)    :: dy1
  real(8)    :: dyn
  integer(4) :: k
  integer(4) :: kf
  real(8)    :: pl(0:n)
  real(8)    :: x
  real(8)    :: y0
  real(8)    :: y1
  real(8)    :: yn

  a = 2.0D+00
  b = 0.0D+00
  c = 1.0D+00
  y0 = 1.0D+00
  y1 = 2.0D+00 * x
  dy0 = 0.0D+00
  dy1 = 2.0D+00
  pl(0) = 1.0D+00
  pl(1) = 2.0D+00 * x
  dpl(0) = 0.0D+00
  dpl(1) = 2.0D+00

  if ( kf == 1 ) then
     y1 = x
     dy1 = 1.0D+00
     pl(1) = x
     dpl(1) = 1.0D+00
  else if ( kf == 3 ) then
     y1 = 1.0D+00 - x
     dy1 = -1.0D+00
     pl(1) = 1.0D+00 - x
     dpl(1) = -1.0D+00
  end if

  do k = 2, n

     if ( kf == 3 ) then
        a = -1.0D+00 / k
        b = 2.0D+00 + a
        c = 1.0D+00 + a
     else if ( kf == 4 ) then
        c = 2.0D+00 * ( k - 1.0D+00 )
     end if

     yn = ( a * x + b ) * y1 - c * y0
     dyn = a * y1 + ( a * x + b ) * dy1 - c * dy0
     pl(k) = yn
     dpl(k) = dyn
     y0 = y1
     y1 = yn
     dy0 = dy1
     dy1 = dyn

  end do

  return
end subroutine othpl




  !*****************************************************************************80
  !
  !! PBDV computes parabolic cylinder functions Dv(x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) DV(0:*), DP(0:*), the values of
  !    Dn+v0(x), Dn+v0'(x).
  !
  !    Output, real(8) PDF, PDD, the values of Dv(x) and Dv'(x).
  !
subroutine pbdv ( v, x, dv, dp, pdf, pdd )
  !PBDV computes parabolic cylinder functions Dv(x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order
  !@param x: in: the argument
  !@param dv: out: Dn+v0(x)
  !@param dp: out: Dn+v0'(x)
  !@param pdf: out: Dv(x)
  !@param pdd: out: Dv'(x)
  !
  implicit none

  real(8)    :: dp(0:*)
  real(8)    :: dv(0:*)
  real(8)    :: ep
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: ja
  integer(4) :: k
  integer(4) :: l
  integer(4) :: m
  integer(4) :: na
  integer(4) :: nk
  integer(4) :: nv
  real(8)    :: pd
  real(8)    :: pd0
  real(8)    :: pd1
  real(8)    :: pdd
  real(8)    :: pdf
  real(8)    :: s0
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: v1
  real(8)    :: v2
  real(8)    :: vh
  real(8)    :: x
  real(8)    :: xa

  xa = abs ( x )
  vh = v
  v = v + sign ( 1.0D+00, v )
  nv = int ( v )
  v0 = v - nv
  na = abs ( nv )
  ep = exp ( -0.25D+00 * x * x )

  if ( 1 <= na ) then
     ja = 1
  end if

  if ( 0.0D+00 <= v ) then
     if ( v0 == 0.0D+00 ) then
        pd0 = ep
        pd1 = x * ep
     else
        do l = 0, ja
           v1 = v0 + l
           if ( xa <= 5.8D+00 ) then
              call dvsa ( v1, x, pd1 )
           else
              call dvla ( v1, x, pd1 )
           end if
           if ( l == 0 ) then
              pd0 = pd1
           end if
        end do
     end if

     dv(0) = pd0
     dv(1) = pd1
     do k = 2, na
        pdf = x * pd1 - ( k + v0 - 1.0D+00 ) * pd0
        dv(k) = pdf
        pd0 = pd1
        pd1 = pdf
     end do

  else

     if ( x <= 0.0D+00 ) then

        if ( xa <= 5.8D+00 )  then
           call dvsa ( v0, x, pd0 )
           v1 = v0 - 1.0D+00
           call dvsa ( v1, x, pd1 )
        else
           call dvla ( v0, x, pd0 )
           v1 = v0 - 1.0D+00
           call dvla ( v1, x, pd1 )
        end if

        dv(0) = pd0
        dv(1) = pd1
        do k = 2, na
           pd = ( - x * pd1 + pd0 ) / ( k - 1.0D+00 - v0 )
           dv(k) = pd
           pd0 = pd1
           pd1 = pd
        end do

     else if ( x <= 2.0D+00 ) then

        v2 = nv + v0
        if ( nv == 0 ) then
           v2 = v2 - 1.0D+00
        end if

        nk = int ( - v2 )
        call dvsa ( v2, x, f1 )
        v1 = v2 + 1.0D+00
        call dvsa ( v1, x, f0 )
        dv(nk) = f1
        dv(nk-1) = f0
        do k = nk - 2, 0, -1
           f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
           dv(k) = f
           f1 = f0
           f0 = f
        end do

     else

        if ( xa <= 5.8D+00 ) then
           call dvsa ( v0, x, pd0 )
        else
           call dvla ( v0, x, pd0 )
        end if

        dv(0) = pd0
        m = 100 + na
        f1 = 0.0D+00
        f0 = 1.0D-30
        do k = m, 0, -1
           f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
           if ( k <= na ) then
              dv(k) = f
           end if
           f1 = f0
           f0 = f
        end do
        s0 = pd0 / f
        do k = 0, na
           dv(k) = s0 * dv(k)
        end do

     end if

  end if

  do k = 0, na - 1
     v1 = abs ( v0 ) + k
     if ( 0.0D+00 <= v ) then
        dp(k) = 0.5D+00 * x * dv(k) - dv(k+1)
     else
        dp(k) = -0.5D+00 * x * dv(k) - v1 * dv(k+1)
     end if
  end do

  pdf = dv(na-1)
  pdd = dp(na-1)
  v = vh

  return
end subroutine pbdv




  !*****************************************************************************80
  !
  !! PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) VV(0:*), VP(0:*), the values of Vv(x), Vv'(x).
  !
  !    Output, real(8) PVF, PVD, the values of Vv(x) and Vv'(x).
  !
subroutine pbvv ( v, x, vv, vp, pvf, pvd )
  !PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order
  !@param x: in: the argument
  !@param vv: out: Vv(x)
  !@param vp: out: Vv'(x)
  !@param pvf: out: Vv(x)
  !@param pvd: out: Vv'(x)
  !
  implicit none

  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: ja
  integer(4) :: k
  integer(4) :: kv
  integer(4) :: l
  integer(4) :: m
  integer(4) :: na
  integer(4) :: nv
  real(8)    :: pi
  real(8)    :: pv0
  real(8)    :: pvd
  real(8)    :: pvf
  real(8)    :: q2p
  real(8)    :: qe
  real(8)    :: s0
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: v1
  real(8)    :: v2
  real(8)    :: vh
  real(8)    :: vp(0:*)
  real(8)    :: vv(0:*)
  real(8)    :: x
  real(8)    :: xa

  pi = 3.141592653589793D+00
  xa = abs ( x )
  vh = v
  v = v + sign ( 1.0D+00, v )
  nv = int ( v )
  v0 = v - nv
  na = abs ( nv )
  qe = exp ( 0.25D+00 * x * x )
  q2p = sqrt ( 2.0D+00 / pi )

  if ( 1 <= na ) then
     ja = 1
  end if

  if ( v <= 0.0D+00 ) then

     if ( v0 == 0.0D+00 ) then

        if ( xa <= 7.5D+00 ) then 
           call vvsa ( v0, x, pv0 )
        else
           call vvla ( v0, x, pv0 )
        end if

        f0 = q2p * qe
        f1 = x * f0
        vv(0) = pv0
        vv(1) = f0
        vv(2) = f1

     else

        do l = 0, ja
           v1 = v0 - l
           if ( xa <= 7.5D+00 ) then
              call vvsa ( v1, x, f1 )
           else
              call vvla ( v1, x, f1 )
           end if
           if ( l == 0 ) then
              f0 = f1
           end if
        end do

        vv(0) = f0
        vv(1) = f1

     end if

     if ( v0 == 0.0D+00 ) then
        kv = 3
     else
        kv = 2
     end if

     do k = kv, na
        f = x * f1 + ( k - v0 - 2.0D+00 ) * f0
        vv(k) = f
        f0 = f1
        f1 = f
     end do

  else

     if ( 0.0D+00 <= x .and. x <= 7.5D+00 ) then

        v2 = v
        if ( v2 < 1.0D+00 ) then
           v2 = v2 + 1.0D+00
        end if

        call vvsa ( v2, x, f1 )
        v1 = v2 - 1.0D+00
        kv = int ( v2 )
        call vvsa ( v1, x, f0 )
        vv(kv) = f1
        vv(kv-1) = f0
        do k = kv - 2, 0, - 1
           f = x * f0 - ( k + v0 + 2.0D+00 ) * f1
           if ( k <= na ) then
              vv(k) = f
           end if
           f1 = f0
           f0 = f
        end do

     else if ( 7.5D+00 < x ) then

        call vvla ( v0, x, pv0 )
        m = 100 + abs ( na )
        vv(1) = pv0
        f1 = 0.0D+00
        f0 = 1.0D-40
        do k = m, 0, -1
           f = x * f0 - ( k + v0 + 2.0D+00 ) * f1
           if ( k <= na ) then
              vv(k) = f
           end if
           f1 = f0
           f0 = f
        end do
        s0 = pv0 / f
        do k = 0, na
           vv(k) = s0 * vv(k)
        end do

     else

        if ( xa <= 7.5D+00 ) then
           call vvsa ( v0, x, f0 )
           v1 = v0 + 1.0D+00
           call vvsa ( v1, x, f1 )
        else
           call vvla ( v0, x, f0 )
           v1 = v0 + 1.0D+00
           call vvla ( v1, x, f1 )
        end if

        vv(0) = f0
        vv(1) = f1
        do k = 2, na
           f = ( x * f1 - f0 ) / ( k + v0 )
           vv(k) = f
           f0 = f1
           f1 = f
        end do

     end if

  end if

  do k = 0, na - 1
     v1 = v0 + k
     if ( 0.0D+00 <= v ) then
        vp(k) = 0.5D+00 * x * vv(k) - ( v1 + 1.0D+00 ) * vv(k+1)
     else
        vp(k) = - 0.5D+00 * x * vv(k) + vv(k+1)
     end if
  end do

  pvf = vv(na-1)
  pvd = vp(na-1)
  v = vh

  return
end subroutine pbvv




  !*****************************************************************************80
  !
  !! PBWA computes parabolic cylinder functions W(a,x) and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) A, the parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) W1F, W1D, W2F, W2D, the values of
  !    W(a,x), W'(a,x), W(a,-x), W'(a,-x).
  !
subroutine pbwa ( a, x, w1f, w1d, w2f, w2d )
  !PBWA computes parabolic cylinder functions W(a,x) and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param a: in: the parameter
  !@param x: in: the argument
  !@param w1f: out: W(a,x)
  !@param w1d: out: W'(a,x)
  !@param w2f: out: W(a,-x)
  !@param w2d: out: W'(a,-x)
  !
  implicit none

  real(8)    :: a
  real(8)    :: d(100)
  real(8)    :: d1
  real(8)    :: d2
  real(8)    :: dl
  real(8)    :: eps
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: g1
  real(8)    :: g2
  real(8)    :: h(100)
  real(8)    :: h0
  real(8)    :: h1
  real(8)    :: hl
  integer(4) :: k
  integer(4) :: l1
  integer(4) :: l2
  integer(4) :: m
  real(8)    :: p0
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: ugi
  real(8)    :: ugr
  real(8)    :: vgi
  real(8)    :: vgr
  real(8)    :: w1d
  real(8)    :: w1f
  real(8)    :: w2d
  real(8)    :: w2f
  real(8)    :: x
  real(8)    :: x1
  real(8)    :: x2
  real(8)    :: y1
  real(8)    :: y1d
  real(8)    :: y1f
  real(8)    :: y2d
  real(8)    :: y2f

  eps = 1.0D-15
  p0 = 0.59460355750136D+00

  if ( a == 0.0D+00 ) then
     g1 = 3.625609908222D+00
     g2 = 1.225416702465D+00
  else
     x1 = 0.25D+00
     y1 = 0.5D+00 * a
     call cgama ( x1, y1, 1, ugr, ugi )
     g1 = sqrt ( ugr * ugr + ugi * ugi )
     x2 = 0.75D+00
     call cgama ( x2, y1, 1, vgr, vgi )
     g2 = sqrt ( vgr * vgr + vgi * vgi )
  end if

  f1 = sqrt ( g1 / g2 )
  f2 = sqrt ( 2.0D+00 * g2 / g1 )
  h0 = 1.0D+00
  h1 = a
  h(1) = a
  do l1 = 4, 200, 2
     m = l1 / 2
     hl = a * h1 - 0.25D+00 * ( l1 - 2.0D+00 ) * ( l1 - 3.0D+00 ) * h0
     h(m) = hl
     h0 = h1
     h1 = hl
  end do
  y1f = 1.0D+00
  r = 1.0D+00
  do k = 1, 100
     r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k - 1.0D+00 ) )
     r1 = h(k) * r
     y1f = y1f + r1
     if ( abs ( r1 / y1f ) <= eps .and. 30 < k ) then
        exit
     end if
  end do

  y1d = a
  r = 1.0D+00
  do k = 1, 100
     r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k + 1.0D+00 ) )
     r1 = h(k+1) * r
     y1d = y1d + r1
     if ( abs ( r1 / y1d ) <= eps .and. 30 < k ) then
        exit
     end if
  end do

  y1d = x * y1d
  d1 = 1.0D+00
  d2 = a
  d(1) = 1.0D+00
  d(2) = a
  do l2 = 5, 160, 2
     m = ( l2 + 1 ) / 2
     dl = a * d2 - 0.25D+00 * ( l2 - 2.0D+00 ) * ( l2 - 3.0D+00 ) * d1
     d(m) = dl
     d1 = d2
     d2 = dl
  end do

  y2f = 1.0D+00
  r = 1.0D+00
  do k = 1, 100
     r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k + 1.0D+00 ) )
     r1 = d(k+1) * r
     y2f = y2f + r1
     if ( abs ( r1 / y2f ) <= eps .and. 30 < k ) then
        exit
     end if
  end do

  y2f = x * y2f
  y2d = 1.0D+00
  r = 1.0D+00
  do k = 1, 100
     r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k - 1.0D+00 ) )
     r1 = d(k+1) * r
     y2d = y2d + r1
     if ( abs ( r1 / y2d ) <= eps .and. 30 < k ) then
        exit
     end if
  end do

  w1f = p0 * ( f1 * y1f - f2 * y2f )
  w2f = p0 * ( f1 * y1f + f2 * y2f )
  w1d = p0 * ( f1 * y1d - f2 * y2d )
  w2d = p0 * ( f1 * y1d + f2 * y2d )

  return
end subroutine pbwa




  !*****************************************************************************80
  !
  !! PSI computes the PSI function.
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by 
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
  !    incorporate this routine into a user program that the copyright 
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) PS, the value of the PSI function.
  !
subroutine psi ( x, ps )
  !PSI computes the PSI function.
  !
  !Original FORTRAN77 by Shanjie Zhang, Jianming Jin. FORTRAN90 version by John Burkardt.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param ps: out: the value of the PSI function
  !
  implicit none

  real(8), parameter  :: a1 = -0.83333333333333333D-01
  real(8), parameter  :: a2 =  0.83333333333333333D-02
  real(8), parameter  :: a3 = -0.39682539682539683D-02
  real(8), parameter  :: a4 =  0.41666666666666667D-02
  real(8), parameter  :: a5 = -0.75757575757575758D-02
  real(8), parameter  :: a6 =  0.21092796092796093D-01
  real(8), parameter  :: a7 = -0.83333333333333333D-01
  real(8), parameter  :: a8 =  0.4432598039215686D+00
  real(8), parameter  :: el = 0.5772156649015329D+00
  integer(4)          :: k
  integer(4)          :: n
  real(8), parameter  :: pi = 3.141592653589793D+00
  real(8)             :: ps
  real(8)             :: s
  real(8)             :: x
  real(8)             :: x2
  real(8)             :: xa

  xa = abs ( x )
  s = 0.0D+00

  if ( x == aint ( x ) .and. x <= 0.0D+00 ) then

     ps = 1.0D+300
     return

  else if ( xa == aint ( xa ) ) then

     n = int ( xa )
     do k = 1, n - 1
        s = s + 1.0D+00 / real ( k, kind = 8 )
     end do

     ps = - el + s

  else if ( xa + 0.5D+00 == aint ( xa + 0.5D+00 ) ) then

     n = int ( xa - 0.5D+00 )

     do k = 1, n
        s = s + 1.0D+00 / real ( 2 * k - 1, kind = 8 )
     end do

     ps = - el + 2.0D+00 * s - 1.386294361119891D+00

  else

     if ( xa < 10.0D+00 ) then

        n = 10 - int ( xa )
        do k = 0, n - 1
           s = s + 1.0D+00 / ( xa + real ( k, kind = 8 ) )
        end do

        xa = xa + real ( n, kind = 8 )

     end if

     x2 = 1.0D+00 / ( xa * xa )

     ps = log ( xa ) - 0.5D+00 / xa + x2 * ((((((( &
          a8   &
          * x2 + a7 ) &
          * x2 + a6 ) &
          * x2 + a5 ) &
          * x2 + a4 ) &
          * x2 + a3 ) &
          * x2 + a2 ) &
          * x2 + a1 )

     ps = ps - s

  end if

  if ( x < 0.0D+00 ) then
     ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
  end if

  return
end subroutine psi




  !*****************************************************************************80
  !
  !! QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) CK(*), ?
  !
  !    Input, real(8) CK1, ?
  !
  !    Output, real(8) QS, ?
  !
  !    Output, real(8) QT, ?
  !
subroutine qstar ( m, n, c, ck, ck1, qs, qt )
  !QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param ck: in: not documented in the original header
  !@param ck1: in: not documented in the original header
  !@param qs: out: not documented in the original header
  !@param qt: out: not documented in the original header
  !
  implicit none

  real(8)    :: ap(200)
  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: ck1
  integer(4) :: i
  integer(4) :: ip
  integer(4) :: k
  integer(4) :: l
  integer(4) :: m
  integer(4) :: n
  real(8)    :: qs
  real(8)    :: qs0
  real(8)    :: qt
  real(8)    :: r
  real(8)    :: s
  real(8)    :: sk

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  r = 1.0D+00 / ck(1) ** 2
  ap(1) = r
  do i = 1, m
     s = 0.0D+00
     do l = 1, i
        sk = 0.0D+00
        do k = 0, l
           sk = sk + ck(k+1) * ck(l-k+1)
        end do
        s = s + sk * ap(i-l+1)
     end do
     ap(i+1) = -r * s
  end do

  qs0 = ap(m+1)     
  do l = 1, m
     r = 1.0D+00
     do k = 1, l
        r = r * ( 2.0D+00 * k + ip ) &
             * ( 2.0D+00 * k - 1.0D+00 + ip ) / ( 2.0D+00 * k ) ** 2
     end do
     qs0 = qs0 + ap(m-l+1) * r
  end do

  qs = ( -1.0D+00 ) ** ip * ck1 * ( ck1 * qs0 ) / c
  qt = - 2.0D+00 / ck1 * qs

  return
end subroutine qstar




  !*****************************************************************************80
  !
  !! RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of jn(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) RJ(0:N), the values of x jn(x).
  !
  !    Output, real(8) DJ(0:N), the values of [x jn(x)]'.
  !
subroutine rctj ( n, x, nm, rj, dj )
  !RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of jn(x)
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param rj: out: the values of x jn(x)
  !@param dj: out: the values of [x jn(x)]'
  !
  implicit none

  integer(4) :: n

  real(8)    :: cs
  real(8)    :: dj(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: rj(0:n)
  real(8)    :: rj0
  real(8)    :: rj1
  real(8)    :: x

  nm = n

  if ( abs ( x ) < 1.0D-100 ) then
     do k = 0, n
        rj(k) = 0.0D+00
        dj(k) = 0.0D+00
     end do
     dj(0) = 1.0D+00
     return
  end if

  rj(0) = sin ( x )
  rj(1) = rj(0) / x - cos ( x )
  rj0 = rj(0)
  rj1 = rj(1)

  if ( 2 <= n ) then

     m = msta1 ( x, 200 )

     if ( m < n ) then
        nm = m
     else
        m = msta2 ( x, n, 15 )
     end if

     f0 = 0.0D+00
     f1 = 1.0D-100
     do k = m, 0, -1
        f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
        if ( k <= nm ) then
           rj(k) = f
        end if
        f0 = f1
        f1 = f
     end do

     if ( abs ( rj1 ) < abs ( rj0 ) ) then
        cs = rj0 / f
     else
        cs = rj1 / f0
     end if

     do k = 0, nm
        rj(k) = cs * rj(k)
     end do

  end if

  dj(0) = cos ( x )
  do k = 1, nm
     dj(k) = - k * rj(k) / x + rj(k-1)
  end do

  return
end subroutine rctj




  !*****************************************************************************80
  !
  !! RCTY computes Riccati-Bessel function of the second kind, and derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of yn(x).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) RY(0:N), the values of x yn(x).
  !
  !    Output, real(8) DY(0:N), the values of [x yn(x)]'.
  !
subroutine rcty ( n, x, nm, ry, dy )
  !RCTY computes Riccati-Bessel function of the second kind, and derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of yn(x)
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param ry: out: the values of x yn(x)
  !@param dy: out: the values of [x yn(x)]'
  !
  implicit none

  integer(4) :: n

  real(8)    :: dy(0:n)
  integer(4) :: k
  integer(4) :: nm
  real(8)    :: rf0
  real(8)    :: rf1
  real(8)    :: rf2
  real(8)    :: ry(0:n)
  real(8)    :: x

  nm = n

  if ( x < 1.0D-60 ) then
     do k = 0, n
        ry(k) = -1.0D+300
        dy(k) = 1.0D+300
     end do
     ry(0) = -1.0D+00
     dy(0) = 0.0D+00
     return
  end if

  ry(0) = - cos ( x )
  ry(1) = ry(0) / x - sin ( x )
  rf0 = ry(0)
  rf1 = ry(1)
  do k = 2, n
     rf2 = ( 2.0D+00 * k - 1.0D+00 ) * rf1 / x - rf0
     if ( 1.0D+300 < abs ( rf2 ) ) then
        exit
     end if
     ry(k) = rf2
     rf0 = rf1
     rf1 = rf2
  end do

  nm = k - 1
  dy(0) = sin ( x )
  do k = 1, nm
     dy(k) = - k * ry(k) / x + ry(k-1)
  end do

  return
end subroutine rcty




  !*****************************************************************************80
  !
  !! REFINE refines an estimate of the characteristic value of Mathieu functions.
  !
  !  Discussion:
  !
  !    This procedure calculates the accurate characteristic value
  !    by the secant method.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    20 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) KD, the case code:
  !    1, for cem(x,q)  ( m = 0,2,4,...)
  !    2, for cem(x,q)  ( m = 1,3,5,...)
  !    3, for sem(x,q)  ( m = 1,3,5,...)
  !    4, for sem(x,q)  ( m = 2,4,6,...)
  !
  !    Input, integer(4) M, the order of the Mathieu functions.
  !
  !    Input, real(8) Q, the parameter of the Mathieu functions.
  !
  !    Input/output, real(8) A, the characteristic value, which
  !    should have been refined on output.
  !
subroutine refine ( kd, m, q, a, iflag )
  !REFINE refines an estimate of the characteristic value of Mathieu functions.
  !
  !This procedure calculates the accurate characteristic value
  !by the secant method.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KD: the case code: 1, for cem(x,q) ( m = 0,2,4,...) 2, for cem(x,q) ( m = 1,3,5,...) 3, for sem(x,q) ( m =
  !  1,3,5,...) 4, for sem(x,q) ( m = 2,4,6,...)
  !
  !@param kd: in: case code: 1 cem (m=0,2,4,..), 2 cem (m=1,3,5,..), 3 sem (m=1,3,5,..), 4 sem (m=2,4,6,..)
  !@param m: in: the order of the Mathieu functions
  !@param q: in: the parameter of the Mathieu functions
  !@param a: in/out: the characteristic value, which should have been refined on output
  !@param iflag: set to -10 on output if the refinement failed
  !
  implicit none

  real(8)    :: a
  real(8)    :: ca
  real(8)    :: delta
  real(8)    :: eps
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: it
  integer(4) :: iflag
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: mj
  real(8)    :: q
  real(8)    :: x
  real(8)    :: x0
  real(8)    :: x1

  eps = 1.0D-14
  mj = 10 + m
  ca = a
  delta = 0.0D+00
  x0 = a
  call cvf ( kd, m, q, x0, mj, f0 )
  x1 = 1.002D+00 * a
  call cvf ( kd, m, q, x1, mj, f1 )

  do

     do it = 1, 100
        mj = mj + 1
        x = x1 - ( x1 - x0 ) / ( 1.0D+00 - f0 / f1 )
        call cvf ( kd, m, q, x, mj, f )
        if ( abs ( 1.0D+00 - x1 / x ) < eps .or. f == 0.0D+00 ) then
           exit
        end if
        x0 = x1
        f0 = f1
        x1 = x
        f1 = f
     end do

     a = x

     if ( 0.05D+00 < delta ) then
        a = ca
        if ( iflag < 0 ) then
           iflag = -10
        end if
        return
     end if

     if ( abs ( ( a - ca ) / ca )  <= 0.05D+00 ) then
        exit
     end if

     x0 = ca
     delta = delta + 0.005D+00
     call cvf ( kd, m, q, x0, mj, f0 )
     x1 = ( 1.0D+00 + delta ) * ca
     call cvf ( kd, m, q, x1, mj, f1 )

  end do

  return
end subroutine refine




  !*****************************************************************************80
  !
  !! RMN1 computes prolate and oblate spheroidal functions of the first kind.
  !
  !  Discussion:
  !
  !    This procedure computes prolate and oblate spheroidal radial
  !    functions of the first kind for given m, n, c and x.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) DF(*), the expansion coefficients.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) R1F, R1D, the function and derivative.
  !
subroutine rmn1 ( m, n, c, x, df, kd, r1f, r1d )
  !RMN1 computes prolate and oblate spheroidal functions of the first kind.
  !
  !This procedure computes prolate and oblate spheroidal radial
  !functions of the first kind for given m, n, c and x.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param df: in: the expansion coefficients
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param r1f: out: the function value
  !@param r1d: out: the derivative value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: b0
  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: cx
  real(8)    :: df(200)
  real(8)    :: dj(0:251)
  real(8)    :: eps
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: l
  integer(4) :: lg
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm1
  integer(4) :: nm2
  integer(4) :: np
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r1
  real(8)    :: r1d
  real(8)    :: r1f
  real(8)    :: r2
  real(8)    :: r3
  real(8)    :: reg
  real(8)    :: sa0
  real(8)    :: sj(0:251)
  real(8)    :: suc
  real(8)    :: sud
  real(8)    :: sum
  real(8)    :: sw
  real(8)    :: sw1
  real(8)    :: x

  eps = 1.0D-14
  nm1 = int ( ( n - m ) / 2 )
  if ( n - m == 2 * nm1 ) then
     ip = 0
  else
     ip = 1
  end if
  nm = 25 + nm1 + int ( c )
  reg = 1.0D+00
  if ( 80 < m + nm ) then
     reg = 1.0D-200
  end if
  r0 = reg
  do j = 1, 2 * m + ip
     r0 = r0 * j
  end do
  r = r0    
  suc = r * df(1)
  do k = 2, nm
     r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
          / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     suc = suc + r * df(k)

     if ( nm1 < k .and. abs ( suc - sw ) < abs ( suc ) * eps ) then
        exit
     end if

     sw = suc

  end do

  if ( x == 0.0D+00 ) then

     call sckb ( m, n, c, df, ck )
     sum = 0.0D+00
     do j = 1, nm
        sum = sum + ck(j)
        if ( abs ( sum - sw1 ) < abs ( sum ) * eps ) then
           exit
        end if
        sw1 = sum
     end do

     r1 = 1.0D+00
     do j = 1, ( n + m + ip ) / 2
        r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
     end do

     r2 = 1.0D+00
     do j = 1, m
        r2 = 2.0D+00 * c * r2 * j
     end do

     r3 = 1.0D+00
     do j = 1, ( n - m - ip ) / 2
        r3 = r3 * j
     end do

     sa0 = ( 2.0D+00 * ( m + ip ) + 1.0D+00 ) * r1 &
          / ( 2.0D+00 ** n * c ** ip * r2 * r3 )

     if ( ip == 0 ) then
        r1f = sum / ( sa0 * suc ) * df(1) * reg
        r1d = 0.0D+00
     else if ( ip == 1 ) then
        r1f = 0.0D+00
        r1d = sum / ( sa0 * suc ) * df(1) * reg
     end if

     return

  end if

  cx = c * x
  nm2 = 2 * nm + m
  call sphj ( nm2, cx, nm2, sj, dj )
  a0 = ( 1.0D+00 - kd / ( x * x ) ) ** ( 0.5D+00 * m ) / suc  
  r1f = 0.0D+00
  do k = 1, nm
     l = 2 * k + m - n - 2 + ip
     if ( l == 4 * int ( l / 4 ) ) then
        lg = 1
     else
        lg = -1
     end if
     if ( k == 1 ) then
        r = r0
     else
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
             / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     end if
     np = m + 2 * k - 2 + ip
     r1f = r1f + lg * r * df(k) * sj(np)
     if ( nm1 < k .and. abs ( r1f - sw ) < abs ( r1f ) * eps ) then
        exit
     end if
     sw = r1f
  end do

  r1f = r1f * a0
  b0 = kd * m / x ** 3.0D+00 / ( 1.0D+00 - kd / ( x * x ) ) * r1f    
  sud = 0.0D+00

  do k = 1, nm

     l = 2 * k + m - n - 2 + ip

     if ( l == 4 * int ( l / 4 ) ) then
        lg = 1
     else
        lg = -1
     end if

     if ( k == 1 ) then
        r = r0
     else
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
             / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     end if

     np = m + 2 * k - 2 + ip
     sud = sud + lg * r * df(k) * dj(np)
     if ( nm1 < k .and. abs ( sud - sw ) < abs ( sud ) * eps ) then
        exit
     end if
     sw = sud
  end do

  r1d = b0 + a0 * c * sud

  return
end subroutine rmn1




  !*****************************************************************************80
  !
  !! RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
  !
  !  Discussion:
  !
  !    This procedure computes prolate and oblate spheroidal radial functions 
  !    of the second kind for given m, n, c and a large cx.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    30 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) DF(*), the expansion coefficients.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) R2F, R2D, the function and derivative values.
  !
subroutine rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
  !RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
  !
  !This procedure computes prolate and oblate spheroidal radial functions
  !of the second kind for given m, n, c and a large cx.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param df: in: the expansion coefficients
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param r2f: out: the function value
  !@param r2d: out: the derivative value
  !@param id: estimated accuracy, number of digits
  !
  implicit none

  real(8)    :: a0
  real(8)    :: b0
  real(8)    :: c
  real(8)    :: cx
  real(8)    :: df(200)
  real(8)    :: dy(0:251)
  real(8)    :: eps
  real(8)    :: eps1
  real(8)    :: eps2
  integer(4) :: id
  integer(4) :: id1
  integer(4) :: id2
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: l
  integer(4) :: lg
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm1
  integer(4) :: nm2
  integer(4) :: np
  real(8)    :: r
  real(8)    :: r0
  real(8)    :: r2d
  real(8)    :: r2f
  real(8)    :: reg
  real(8)    :: sw
  real(8)    :: suc
  real(8)    :: sud
  real(8)    :: sy(0:251)
  real(8)    :: x

  eps = 1.0D-14

  nm1 = int ( ( n - m ) / 2 )

  if ( n - m == 2 * nm1 ) then
     ip = 0
  else
     ip = 1
  end if
  nm = 25 + nm1 + int ( c )

  if ( 80 < m + nm ) then
     reg = 1.0D-200
  else
     reg = 1.0D+00
  end if
  nm2 = 2 * nm + m
  cx = c * x
  call sphy ( nm2, cx, nm2, sy, dy )
  r0 = reg
  do j = 1, 2 * m + ip
     r0 = r0 * j
  end do
  r = r0    
  suc = r * df(1)
  do k = 2, nm
     r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
          / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     suc = suc + r * df(k)
     if ( nm1 < k .and. abs ( suc - sw ) < abs ( suc ) * eps ) then
        exit
     end if
     sw = suc
  end do

  a0 = ( 1.0D+00 - kd / ( x * x ) ) ** ( 0.5D+00 * m ) / suc
  r2f = 0.0D+00
  do k = 1, nm
     l = 2 * k + m - n - 2 + ip
     if ( l == 4 * int ( l / 4 ) ) then
        lg = 1
     else
        lg = -1
     end if

     if ( k == 1 ) then
        r = r0
     else
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
             / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     end if

     np = m + 2 * k - 2 + ip
     r2f = r2f + lg * r * ( df(k) * sy(np) )
     eps1 = abs ( r2f - sw )
     if ( nm1 < k .and. eps1 < abs ( r2f ) * eps ) then
        exit
     end if
     sw = r2f
  end do

  id1 = int ( log10 ( eps1 / abs ( r2f ) + eps ) )
  r2f = r2f * a0

  if ( nm2 <= np ) then
     id = 10
     return
  end if

  b0 = kd * m / x ** 3.0D+00 / ( 1.0D+00 - kd / ( x * x ) ) * r2f                
  sud = 0.0D+00
  do k = 1, nm
     l = 2 * k + m - n - 2 + ip
     if ( l == 4 * int ( l / 4 ) ) then
        lg = 1
     else
        lg = -1
     end if
     if (k == 1) then
        r = r0
     else
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 ) &
             / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
     end if
     np = m + 2 * k - 2 + ip
     sud = sud + lg * r * ( df(k) * dy(np) )
     eps2 = abs ( sud - sw )
     if ( nm1 < k .and. eps2 < abs ( sud ) * eps ) then
        exit
     end if
     sw = sud
  end do

  r2d = b0 + a0 * c * sud
  id2 = int ( log10 ( eps2 / abs ( sud ) + eps ) )
  id = max ( id1, id2 )

  return
end subroutine rmn2l




  !*****************************************************************************80
  !
  !! RMN2SO: oblate radial functions of the second kind with small argument.
  !
  !  Discussion:
  !
  !    This procedure computes oblate radial functions of the second kind
  !    with a small argument, Rmn(-ic,ix) and Rmn'(-ic,ix).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    27 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, real(8) DF(*), the expansion coefficients.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) R2F, R2D, the values of Rmn(-ic,ix) 
  !    and Rmn'(-ic,ix).
  !
subroutine rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )
  !RMN2SO: oblate radial functions of the second kind with small argument.
  !
  !This procedure computes oblate radial functions of the second kind
  !with a small argument, Rmn(-ic,ix) and Rmn'(-ic,ix).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param cv: in: the characteristic value
  !@param df: in: the expansion coefficients
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param r2f: out: Rmn(-ic,ix)
  !@param r2d: out: Rmn'(-ic,ix)
  !
  implicit none

  real(8)    :: bk(200)
  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: ck1
  real(8)    :: ck2
  real(8)    :: cv
  real(8)    :: df(200)
  real(8)    :: dn(200)
  real(8)    :: eps
  real(8)    :: gd
  real(8)    :: gf
  real(8)    :: h0
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  real(8)    :: pi
  real(8)    :: qs
  real(8)    :: qt
  real(8)    :: r1d
  real(8)    :: r1f
  real(8)    :: r2d
  real(8)    :: r2f
  real(8)    :: sum
  real(8)    :: sw
  real(8)    :: x

  if ( abs ( df(1) ) <= 1.0D-280 ) then
     r2f = 1.0D+300
     r2d = 1.0D+300
     return
  end if

  eps = 1.0D-14
  pi = 3.141592653589793D+00
  nm = 25 + int ( ( n - m ) / 2 + c )
  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  call sckb ( m, n, c, df, ck )
  call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
  call qstar ( m, n, c, ck, ck1, qs, qt )
  call cbk ( m, n, c, cv, qt, ck, bk )

  if ( x == 0.0D+00 ) then

     sum = 0.0D+00
     do j = 1, nm
        sum = sum + ck(j)
        if ( abs ( sum - sw ) < abs ( sum ) * eps ) then
           exit
        end if
        sw = sum
     end do

     if ( ip == 0 ) then
        r1f = sum / ck1
        r2f = - 0.5D+00 * pi * qs * r1f
        r2d = qs * r1f + bk(1)
     else if ( ip == 1 ) then
        r1d = sum / ck1
        r2f = bk(1)
        r2d = -0.5D+00 * pi * qs * r1d
     end if

     return

  else

     call gmn ( m, n, c, x, bk, gf, gd )
     call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
     h0 = atan ( x ) - 0.5D+00 * pi
     r2f = qs * r1f * h0 + gf
     r2d = qs * ( r1d * h0 + r1f / ( 1.0D+00 + x * x ) ) + gd

  end if

  return
end subroutine rmn2so




  !*****************************************************************************80
  !
  !! RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, real(8) DF(*), the expansion coefficients.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) R2F, R2D, the values of the function and 
  !    its derivative.
  !
subroutine rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )
  !RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param cv: in: the characteristic value
  !@param df: in: the expansion coefficients
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param r2f: out: the values of the function
  !@param r2d: out: its derivative
  !
  implicit none

  real(8)    :: c
  real(8)    :: ck1
  real(8)    :: ck2
  real(8)    :: cv
  real(8)    :: df(200)
  real(8)    :: dn(200)
  real(8)    :: eps
  real(8)    :: ga
  real(8)    :: gb
  real(8)    :: gc
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: j1
  integer(4) :: j2
  integer(4) :: k
  integer(4) :: kd
  integer(4) :: ki
  integer(4) :: l1
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm1
  integer(4) :: nm2
  integer(4) :: nm3
  real(8)    :: pd(0:251)
  real(8)    :: pm(0:251)
  real(8)    :: qd(0:251)
  real(8)    :: qm(0:251)
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: r2d
  real(8)    :: r2f
  real(8)    :: r3
  real(8)    :: r4
  real(8)    :: sd
  real(8)    :: sd0
  real(8)    :: sd1
  real(8)    :: sd2
  real(8)    :: sdm
  real(8)    :: sf
  real(8)    :: spd1
  real(8)    :: spd2
  real(8)    :: spl
  real(8)    :: su0
  real(8)    :: su1
  real(8)    :: su2
  real(8)    :: sum
  real(8)    :: sw
  real(8)    :: x

  if ( abs ( df(1) ) < 1.0D-280 ) then
     r2f = 1.0D+300
     r2d = 1.0D+300
     return
  end if

  eps = 1.0D-14

  nm1 = int ( ( n - m ) / 2 )

  if ( n - m .eq. 2 * nm1 ) then
     ip = 0
  else
     ip = 1
  end if

  nm = 25 + nm1 + int ( c )
  nm2 = 2 * nm + m
  call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
  call lpmns ( m, nm2, x, pm, pd )
  call lqmns ( m, nm2, x, qm, qd )

  su0 = 0.0D+00
  do k = 1, nm
     j = 2 * k - 2 + m + ip
     su0 = su0 + df(k) * qm(j)
     if ( nm1 < k .and. abs ( su0 - sw ) < abs ( su0 ) * eps ) then
        exit                                                               
     end if
     sw = su0
  end do

  sd0 = 0.0D+00

  do k = 1, nm
     j = 2 * k - 2 + m + ip
     sd0 = sd0 + df(k) * qd(j)
     if ( nm1 < k .and. abs ( sd0 - sw ) < abs ( sd0 ) * eps ) then
        exit
     end if
     sw = sd0
  end do

  su1 = 0.0D+00
  sd1 = 0.0D+00
  do k = 1, m
     j = m - 2 * k + ip
     if ( j < 0 ) then
        j = - j - 1
     end if
     su1 = su1 + dn(k) * qm(j)
     sd1 = sd1 + dn(k) * qd(j)
  end do

  ga = ( ( x - 1.0D+00 ) / ( x + 1.0D+00 ) ) ** ( 0.5D+00 * m )

  do k = 1, m

     j = m - 2 * k + ip

     if ( 0 <= j ) then
        cycle
     end if

     if ( j < 0 ) then
        j = - j - 1
     end if
     r1 = 1.0D+00
     do j1 = 1, j
        r1 = ( m + j1 ) * r1
     end do
     r2 = 1.0D+00
     do j2 = 1, m - j - 2
        r2 = j2 * r2
     end do
     r3 = 1.0D+00
     sf = 1.0D+00
     do l1 = 1, j
        r3 = 0.5D+00 * r3 * ( - j + l1 - 1.0D+00 ) * ( j + l1 ) &
             / ( ( m + l1 ) * l1 ) * ( 1.0D+00 - x )
        sf = sf + r3
     end do

     if ( m - j <= 1 ) then
        gb = 1.0D+00
     else
        gb = ( m - j - 1.0D+00 ) * r2
     end if

     spl = r1 * ga * gb * sf
     su1 = su1 + ( -1 ) ** ( j + m ) * dn(k) * spl
     spd1 = m / ( x * x - 1.0D+00 ) * spl
     gc = 0.5D+00 * j * ( j + 1.0 ) / ( m + 1.0D+00 )
     sd = 1.0D+00
     r4 = 1.0D+00
     do l1 = 1, j - 1
        r4 = 0.5D+00 * r4 * ( - j + l1 ) * ( j + l1 + 1.0D+00 ) &
             / ( ( m + l1 + 1.0D+00 ) * l1 ) * ( 1.0D+00 - x )
        sd = sd + r4
     end do

     spd2 = r1 * ga * gb * gc * sd
     sd1 = sd1 + ( - 1 ) ** ( j + m ) * dn(k) * ( spd1 + spd2 )

  end do

  su2 = 0.0D+00
  ki = ( 2 * m + 1 + ip ) / 2
  nm3 = nm + ki
  do k = ki, nm3
     j = 2 * k - 1 - m - ip
     su2 = su2 + dn(k) * pm(j)
     if ( m < j .and. &
          abs ( su2 - sw ) < abs ( su2 ) * eps ) then
        exit
     end if
     sw = su2
  end do

  sd2 = 0.0D+00

  do k = ki, nm3
     j = 2 * k - 1 - m - ip
     sd2 = sd2 + dn(k) * pd(j)
     if ( m < j .and. &
          abs ( sd2 - sw ) < abs ( sd2 ) * eps ) then
        exit
     end if
     sw = sd2
  end do

  sum = su0 + su1 + su2
  sdm = sd0 + sd1 + sd2
  r2f = sum / ck2
  r2d = sdm / ck2

  return
end subroutine rmn2sp




  !*****************************************************************************80
  !
  !! RSWFO computes prolate spheroidal radial function of first and second kinds.
  !
  !  Discussion:
  !
  !    This procedure computes oblate radial functions of the first
  !    and second kinds, and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  ! 
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, integer(4) KF, the function code.
  !    1, for the first kind
  !    2, for the second kind
  !    3, for both the first and second kinds.
  !
  !    Output, real(8) R1F, the radial function of the first kind;
  !
  !    Output, real(8) R1D, the derivative of the radial function of
  !    the first kind;
  !
  !    Output, real(8) R2F, the radial function of the second kind;
  !
  !    Output, real(8) R2D, the derivative of the radial function of
  !    the second kind;
  !
subroutine rswfo ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )
  !RSWFO computes prolate spheroidal radial function of first and second kinds.
  !
  !This procedure computes oblate radial functions of the first
  !and second kinds, and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KF: the function code. 1, for the first kind 2, for the second kind 3, for both the first and second kinds
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument
  !@param cv: in: the characteristic value
  !@param kf: in: function code: 1 first kind, 2 second kind, 3 both
  !@param r1f: out: the radial function of the first kind;
  !@param r1d: out: the derivative of the radial function of the first kind;
  !@param r2f: out: the radial function of the second kind;
  !@param r2d: out: the derivative of the radial function of the second kind;
  !
  implicit none

  real(8)    :: c
  real(8)    :: cv
  real(8)    :: df(200)
  integer(4) :: id
  integer(4) :: kd
  integer(4) :: kf
  integer(4) :: m
  integer(4) :: n
  real(8)    :: r1d
  real(8)    :: r1f
  real(8)    :: r2d
  real(8)    :: r2f
  real(8)    :: x

  kd = -1
  call sdmn ( m, n, c, cv, kd, df )

  if ( kf /= 2 ) then
     call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
  end if

  if ( 1 < kf ) then
     id = 10
     if ( 1.0D-08 < x ) then
        call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
     end if
     if ( -1 < id ) then
        call rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )
     end if
  end if

  return
end subroutine rswfo




  !*****************************************************************************80
  !
  !! RSWFP computes prolate spheroidal radial function of first and second kinds.
  !
  !  Discussion:
  !
  !    This procedure computes prolate spheriodal radial functions of the
  !    first and second kinds, and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    07 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter;  M = 0, 1, 2, ...
  !
  !    Input, integer(4) N, mode parameter, N = M, M + 1, M + 2, ...
  !
  !    Input, real(8) C, spheroidal parameter.
  !
  !    Input, real(8) X, the argument of the radial function, 1 < X.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, integer(4) KF, the function code.
  !    1, for the first kind
  !    2, for the second kind
  !    3, for both the first and second kinds.
  !
  !    Output, real(8) R1F, the radial function of the first kind;
  !
  !    Output, real(8) R1D, the derivative of the radial function of
  !    the first kind;
  !
  !    Output, real(8) R2F, the radial function of the second kind;
  !
  !    Output, real(8) R2D, the derivative of the radial function of
  !    the second kind;
  !
subroutine rswfp ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )
  !RSWFP computes prolate spheroidal radial function of first and second kinds.
  !
  !This procedure computes prolate spheriodal radial functions of the
  !first and second kinds, and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* KF: the function code. 1, for the first kind 2, for the second kind 3, for both the first and second kinds
  !
  !@param m: in: the mode parameter; M = 0, 1, 2, ...
  !@param n: in: mode parameter, N = M, M + 1, M + 2, ...
  !@param c: in: spheroidal parameter
  !@param x: in: the argument of the radial function, 1 < X
  !@param cv: in: the characteristic value
  !@param kf: in: function code: 1 first kind, 2 second kind, 3 both
  !@param r1f: out: the radial function of the first kind;
  !@param r1d: out: the derivative of the radial function of the first kind;
  !@param r2f: out: the radial function of the second kind;
  !@param r2d: out: the derivative of the radial function of the second kind;
  !
  implicit none

  real(8)    :: c
  real(8)    :: cv
  real(8)    :: df(200)
  integer(4) :: id
  integer(4) :: kd
  integer(4) :: kf
  integer(4) :: m
  integer(4) :: n
  real(8)    :: r1d
  real(8)    :: r1f
  real(8)    :: r2d
  real(8)    :: r2f
  real(8)    :: x

  kd = 1
  call sdmn ( m, n, c, cv, kd, df )

  if ( kf /= 2 ) then
     call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
  end if

  if ( 1 < kf ) then
     call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
     if ( -8 < id ) then
        call rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )
     end if
  end if

  return
end subroutine rswfp




  !*****************************************************************************80
  !
  !! SCKA: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter.
  !
  !    Input, integer(4) N, the mode parameter.
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) CK(*), the expansion coefficients.
  !    CK(1), CK(2),... correspond to c0, c2,..., and so on.
  !       
subroutine scka ( m, n, c, cv, kd, ck )
  !SCKA: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter
  !@param n: in: the mode parameter
  !@param c: in: the spheroidal parameter
  !@param cv: in: the characteristic value
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param ck: out: the expansion coefficients. CK(1), CK(2), ... correspond to c0, c2, ..., and so on
  !
  implicit none

  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: cs
  real(8)    :: cv
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: fl
  real(8)    :: fs
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k1
  integer(4) :: kb
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: s0
  real(8)    :: su1
  real(8)    :: su2

  if ( c <= 1.0D-10 ) then
     c = 1.0D-10
  end if

  nm = 25 + int ( ( n - m ) / 2 + c )
  cs = c * c * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  fs = 1.0D+00
  f1 = 0.0D+00
  f0 = 1.0D-100
  kb = 0
  ck(nm+1) = 0.0D+00

  do k = nm, 1, -1

     f = ((( 2.0D+00 * k + m + ip ) &
          * ( 2.0D+00 * k + m + 1.0D+00 + ip ) - cv + cs ) * f0 &
          - 4.0D+00 * ( k + 1.0D+00 ) * ( k + m + 1.0D+00 ) * f1 ) / cs

     if ( abs ( ck(k+1) ) < abs ( f ) ) then

        ck(k) = f
        f1 = f0
        f0 = f

        if ( 1.0D+100 < abs ( f ) ) then
           do k1 = nm, k, -1
              ck(k1) = ck(k1) * 1.0D-100
           end do
           f1 = f1 * 1.0D-100
           f0 = f0 * 1.0D-100
        end if

     else

        kb = k
        fl = ck(k+1)
        f1 = 1.0D+00
        f2 = 0.25D+00 * ( ( m + ip ) * ( m + ip + 1.0D+00 ) &
             - cv + cs ) / ( m + 1.0D+00 ) * f1
        ck(1) = f1

        if ( kb == 1 ) then
           fs = f2
        else if (kb == 2 ) then
           ck(2) = f2
           fs = 0.125D+00 * ( ( ( m + ip + 2.0D+00 ) &
                * ( m + ip + 3.0D+00 ) - cv + cs ) * f2 &
                - cs * f1 ) / ( m + 2.0D+00 )
        else
           ck(2) = f2
           do j = 3, kb + 1
              f = 0.25D+00 * ( ( ( 2.0D+00 * j + m + ip - 4.0D+00 ) &
                   * ( 2.0D+00 * j + m + ip - 3.0D+00 ) - cv + cs ) * f2 &
                   - cs * f1 ) / ( ( j - 1.0D+00 ) * ( j + m - 1.0D+00 ) )
              if ( j <= kb ) then
                 ck(j) = f
              end if
              f1 = f2
              f2 = f
           end do
           fs = f
        end if

        exit

     end if

  end do

  su1 = 0.0D+00
  do k = 1, kb
     su1 = su1 + ck(k)
  end do

  su2 = 0.0D+00
  do k = kb + 1, nm
     su2 = su2 + ck(k)
  end do

  r1 = 1.0D+00
  do j = 1, ( n + m + ip ) / 2
     r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
  end do

  r2 = 1.0D+00
  do j = 1, ( n - m - ip ) / 2
     r2 = - r2 * j
  end do

  if ( kb == 0 ) then
     s0 = r1 / ( 2.0D+00 ** n * r2 * su2 )
  else
     s0 = r1 / ( 2.0D+00 ** n * r2 * ( fl / fs * su1 + su2 ) )
  end if

  do k = 1, kb
     ck(k) = fl / fs * s0 * ck(k)
  end do

  do k = kb + 1, nm
     ck(k) = s0 * ck(k)
  end do

  return
end subroutine scka




  !*****************************************************************************80
  !
  !! SCKB: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter.
  !
  !    Input, integer(4) N, the mode parameter.
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, real(8) DF(*), the expansion coefficients DK.
  !
  !    Output, real(8) CK(*), the expansion coefficients CK.
  !
subroutine sckb ( m, n, c, df, ck )
  !SCKB: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter
  !@param n: in: the mode parameter
  !@param c: in: the spheroidal parameter
  !@param df: in: the expansion coefficients DK
  !@param ck: out: the expansion coefficients CK
  !
  implicit none

  real(8)    :: c
  real(8)    :: ck(200)
  real(8)    :: d1
  real(8)    :: d2
  real(8)    :: d3
  real(8)    :: df(200)
  real(8)    :: fac
  integer(4) :: i
  integer(4) :: i1
  integer(4) :: i2
  integer(4) :: ip
  integer(4) :: k
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: reg
  real(8)    :: sum
  real(8)    :: sw

  c = max ( c, 1.0D-10 )

  nm = 25 + int ( 0.5D+00 * ( n - m ) + c )

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  if ( 80 < m + nm ) then
     reg = 1.0D-200
  else
     reg = 1.0D+00
  end if

  fac = - 0.5D+00 ** m

  do k = 0, nm - 1

     fac = - fac
     i1 = 2 * k + ip + 1
     r = reg
     do i = i1, i1 + 2 * m - 1
        r = r * i
     end do

     i2 = k + m + ip
     do i = i2, i2 + k - 1
        r = r * ( i + 0.5D+00 )
     end do

     sum = r * df(k+1)
     do i = k + 1, nm
        d1 = 2.0D+00 * i + ip
        d2 = 2.0D+00 * m + d1
        d3 = i + m + ip - 0.5D+00
        r = r * d2 * ( d2 - 1.0D+00 ) * i * ( d3 + k ) &
             / ( d1 * ( d1 - 1.0D+00 ) * ( i - k ) * d3 )
        sum = sum + r * df(i+1)
        if ( abs ( sw - sum ) < abs ( sum ) * 1.0D-14 ) then
           exit
        end if
        sw = sum
     end do

     r1 = reg
     do i = 2, m + k
        r1 = r1 * i
     end do

     ck(k+1) = fac * sum / r1

  end do

  return
end subroutine sckb




  !*****************************************************************************80
  !
  !! SDMN: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    29 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter.
  !
  !    Input, integer(4) N, the mode parameter.
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, real(8) CV, the characteristic value.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) DF(*), expansion coefficients;
  !    DF(1), DF(2), ... correspond to d0, d2, ... for even n-m and d1,
  !    d3, ... for odd n-m
  !
subroutine sdmn ( m, n, c, cv, kd, df )
  !SDMN: expansion coefficients for prolate and oblate spheroidal functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !Details of the arguments:
  !
  !* DF: expansion coefficients; DF(1), DF(2), ... correspond to d0, d2, ... for even n-m and d1, d3, ... for odd n-m
  !
  !@param m: in: the mode parameter
  !@param n: in: the mode parameter
  !@param c: in: the spheroidal parameter
  !@param cv: in: the characteristic value
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param df: out: expansion coefficients: d0, d2, ... for even n-m, d1, d3, ... for odd n-m
  !
  implicit none

  real(8)    :: a(200)
  real(8)    :: c
  real(8)    :: cs
  real(8)    :: cv
  real(8)    :: d(200)
  real(8)    :: d2k
  real(8)    :: df(200)
  real(8)    :: dk0
  real(8)    :: dk1
  real(8)    :: dk2
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  real(8)    :: f2
  real(8)    :: fl
  real(8)    :: fs
  real(8)    :: g(200)
  integer(4) :: i
  integer(4) :: ip
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k1
  integer(4) :: kb
  integer(4) :: kd
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  real(8)    :: r1
  real(8)    :: r3
  real(8)    :: r4
  real(8)    :: s0
  real(8)    :: su1
  real(8)    :: su2
  real(8)    :: sw

  nm = 25 + int ( 0.5D+00 * ( n - m ) + c )

  if ( c < 1.0D-10 ) then
     do i = 1, nm
        df(i) = 0D+00
     end do
     df((n-m)/2+1) = 1.0D+00
     return
  end if

  cs = c * c * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
     ip = 0
  else
     ip = 1
  end if

  do i = 1, nm + 2
     if ( ip == 0 ) then
        k = 2 * ( i - 1 )
     else
        k = 2 * i - 1
     end if
     dk0 = m + k
     dk1 = m + k + 1
     dk2 = 2 * ( m + k )
     d2k = 2 * m + k
     a(i) = ( d2k + 2.0D+00 ) * ( d2k + 1.0D+00 ) &
          / ( ( dk2 + 3.0D+00 ) * ( dk2 + 5.0D+00 ) ) * cs
     d(i) = dk0 * dk1 &
          + ( 2.0D+00 * dk0 * dk1 - 2.0D+00 * m * m - 1.0D+00 ) &
          / ( ( dk2 - 1.0D+00 ) * ( dk2 + 3.0D+00 ) ) * cs
     g(i) = k * ( k - 1.0D+00 ) / ( ( dk2 - 3.0D+00 ) &
          * ( dk2 - 1.0D+00 ) ) * cs
  end do

  fs = 1.0D+00
  f1 = 0.0D+00
  f0 = 1.0D-100
  kb = 0
  df(nm+1) = 0.0D+00

  do k = nm, 1, -1

     f = - ( ( d(k+1) - cv ) * f0 + a(k+1) * f1 ) / g(k+1)

     if ( abs ( df(k+1) ) < abs ( f ) ) then

        df(k) = f
        f1 = f0
        f0 = f
        if ( 1.0D+100 < abs ( f ) ) then
           do k1 = k, nm
              df(k1) = df(k1) * 1.0D-100
           end do
           f1 = f1 * 1.0D-100
           f0 = f0 * 1.0D-100
        end if

     else

        kb = k
        fl = df(k+1)
        f1 = 1.0D-100
        f2 = - ( d(1) - cv ) / a(1) * f1
        df(1) = f1

        if ( kb == 1 ) then

           fs = f2

        else if ( kb == 2 ) then

           df(2) = f2
           fs = - ( ( d(2) - cv ) * f2 + g(2) * f1 ) / a(2)

        else 

           df(2) = f2
           do j = 3, kb + 1
              f = - ( ( d(j-1) - cv ) * f2 + g(j-1) * f1 ) / a(j-1)
              if ( j <= kb ) then
                 df(j) = f
              end if
              if ( 1.0D+100 < abs ( f ) ) then
                 do k1 = 1, j
                    df(k1) = df(k1) * 1.0D-100
                 end do
                 f = f * 1.0D-100
                 f2 = f2 * 1.0D-100
              end if
              f1 = f2
              f2 = f
           end do
           fs = f

        end if

        exit

     end if

  end do

  su1 = 0.0D+00

  r1 = 1.0D+00
  do j = m + ip + 1, 2 * ( m + ip )
     r1 = r1 * j
  end do

  su1 = df(1) * r1
  do k = 2, kb
     r1 = - r1 * ( k + m + ip - 1.5D+00 ) / ( k - 1.0D+00 )
     su1 = su1 + r1 * df(k)
  end do

  su2 = 0.0D+00
  do k = kb + 1, nm
     if ( k /= 1 ) then
        r1 = - r1 * ( k + m + ip - 1.5D+00 ) / ( k - 1.0D+00 )
     end if
     su2 = su2 + r1 * df(k)
     if ( abs ( sw - su2 ) < abs ( su2 ) * 1.0D-14 ) then
        exit
     end if
     sw = su2
  end do

  r3 = 1.0D+00
  do j = 1, ( m + n + ip ) / 2
     r3 = r3 * ( j + 0.5D+00 * ( n + m + ip ) )
  end do

  r4 = 1.0D+00
  do j = 1, ( n - m - ip ) / 2
     r4 = -4.0D+00 * r4 * j
  end do

  s0 = r3 / ( fl * ( su1 / fs ) + su2 ) / r4
  do k = 1, kb
     df(k) = fl / fs * s0 * df(k)
  end do

  do k = kb + 1, nm
     df(k) = s0 * df(k)
  end do

  return
end subroutine sdmn




  !*****************************************************************************80
  !
  !! SEGV computes the characteristic values of spheroidal wave functions.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    28 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) M, the mode parameter.
  !
  !    Input, integer(4) N, the mode parameter.
  !
  !    Input, real(8) C, the spheroidal parameter.
  !
  !    Input, integer(4) KD, the function code.
  !    1, the prolate function.
  !    -1, the oblate function.
  !
  !    Output, real(8) CV, the characteristic value.
  !
  !    Output, real(8) EG(*), the characteristic value for 
  !    mode parameters m and n.  ( L = n - m + 1 )
  !
subroutine segv ( m, n, c, kd, cv, eg )
  !SEGV computes the characteristic values of spheroidal wave functions.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param m: in: the mode parameter
  !@param n: in: the mode parameter
  !@param c: in: the spheroidal parameter
  !@param kd: in: the function code. 1, the prolate function. -1, the oblate function
  !@param cv: out: the characteristic value
  !@param eg: out: the characteristic value for mode parameters m and n. ( L = n - m + 1 )
  !
  implicit none

  real(8)    :: a(300)
  real(8)    :: b(100)
  real(8)    :: c
  real(8)    :: cs
  real(8)    :: cv
  real(8)    :: cv0(100)
  real(8)    :: d(300)
  real(8)    :: d2k
  real(8)    :: dk0
  real(8)    :: dk1
  real(8)    :: dk2
  real(8)    :: e(300)
  real(8)    :: eg(200)
  real(8)    :: f(300)
  real(8)    :: g(300)
  real(8)    :: h(100)
  integer(4) :: i
  integer(4) :: icm
  integer(4) :: j
  integer(4) :: k
  integer(4) :: k1
  integer(4) :: kd
  integer(4) :: l
  integer(4) :: m
  integer(4) :: n
  integer(4) :: nm
  integer(4) :: nm1
  real(8)    :: s
  real(8)    :: t
  real(8)    :: t1
  real(8)    :: x1
  real(8)    :: xa
  real(8)    :: xb

  if ( c < 1.0D-10 ) then
     do i = 1, n
        eg(i) = ( i + m ) * ( i + m - 1.0D+00 )
     end do
     cv = eg(n-m+1)
     return
  end if

  icm = ( n - m + 2 ) / 2
  nm = 10 + int ( 0.5D+00 * ( n - m ) + c )
  cs = c * c * kd

  do l = 0, 1

     do i = 1, nm
        if ( l == 0 ) then
           k = 2 * ( i - 1 )
        else
           k = 2 * i - 1
        end if
        dk0 = m + k
        dk1 = m + k + 1
        dk2 = 2 * ( m + k )
        d2k = 2 * m + k
        a(i) = ( d2k + 2.0D+00 ) * ( d2k + 1.0D+00 ) &
             / ( ( dk2 + 3.0D+00 ) * ( dk2 + 5.0D+00 ) ) * cs
        d(i) = dk0 * dk1 + ( 2.0D+00 * dk0 * dk1 &
             - 2.0 * m * m - 1.0D+00 ) &
             / ( ( dk2 - 1.0D+00 ) * ( dk2 + 3.0D+00 ) ) * cs
        g(i) = k * ( k - 1.0D+00 ) / ( ( dk2 - 3.0D+00 ) &
             * ( dk2 - 1.0D+00 ) ) * cs
     end do

     do k = 2, nm
        e(k) = sqrt ( a(k-1) * g(k) )
        f(k) = e(k) * e(k)
     end do

     f(1) = 0.0D+00
     e(1) = 0.0D+00
     xa = d(nm) + abs ( e(nm) )
     xb = d(nm) - abs ( e(nm) )
     nm1 = nm - 1
     do i = 1, nm1
        t = abs ( e(i) ) + abs ( e(i+1) )
        t1 = d(i) + t
        if ( xa < t1 ) then
           xa = t1
        end if
        t1 = d(i) - t
        if ( t1 < xb ) then
           xb = t1
        end if
     end do

     do i = 1, icm
        b(i) = xa
        h(i) = xb
     end do

     do k = 1, icm

        do k1 = k, icm
           if ( b(k1) < b(k) ) then
              b(k) = b(k1)
              exit
           end if
        end do

        if ( k /= 1 .and. h(k) < h(k-1) ) then
           h(k) = h(k-1)
        end if

        do

           x1 = ( b(k) + h(k) ) /2.0D+00
           cv0(k) = x1

           if ( abs ( ( b(k) - h(k) ) / x1 ) < 1.0D-14 ) then
              exit
           end if

           j = 0
           s = 1.0D+00

           do i = 1, nm

              if ( s == 0.0D+00 ) then
                 s = s + 1.0D-30
              end if
              t = f(i) / s
              s = d(i) - t - x1
              if ( s < 0.0D+00 ) then
                 j = j + 1
              end if
           end do

           if ( j < k ) then

              h(k) = x1

           else

              b(k) = x1
              if ( icm <= j ) then
                 b(icm) = x1
              else
                 if ( h(j+1) < x1 ) then
                    h(j+1) = x1
                 end if
                 if ( x1 < b(j) ) then
                    b(j) = x1
                 end if
              end if

           end if

        end do

        cv0(k) = x1

        if ( l == 0 ) then
           eg(2*k-1) = cv0(k)
        else
           eg(2*k) = cv0(k)
        end if

     end do

  end do

  cv = eg(n-m+1)

  return
end subroutine segv




  !*****************************************************************************80
  !
  !! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    18 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order of In(X).
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) SI(0:N), DI(0:N), the values and derivatives
  !    of the function of orders 0 through N.
  !
subroutine sphi ( n, x, nm, si, di )
  !SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order of In(X)
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param si: out: the values of the function of orders 0 through N
  !@param di: out: the derivatives of the function of orders 0 through N
  !
  implicit none

  integer(4) :: n

  real(8)    :: cs
  real(8)    :: di(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: k
  integer(4) :: m
    ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: si(0:n)
  real(8)    :: si0
  real(8)    :: x

  nm = n

  if ( abs ( x ) < 1.0D-100 ) then
     do k = 0, n
        si(k) = 0.0D+00
        di(k) = 0.0D+00
     end do
     si(0) = 1.0D+00
     di(1) = 0.333333333333333D+00
     return
  end if

  si(0) = sinh ( x ) / x
  si(1) = -( sinh ( x ) / x - cosh ( x ) ) / x
  si0 = si(0)

  if ( 2 <= n ) then

     m = msta1 ( x, 200 )
     if ( m < n ) then
        nm = m
     else
        m = msta2 ( x, n, 15 )
     end if
     f0 = 0.0D+00
     f1 = 1.0D+00-100
     do k = m, 0, -1
        f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x + f0
        if ( k <= nm ) then
           si(k) = f
        end if
        f0 = f1
        f1 = f
     end do
     cs = si0 / f
     do k = 0, nm
        si(k) = cs * si(k)
     end do

  end if

  di(0) = si(1)
  do k = 1, nm
     di(k) = si(k-1) - ( k + 1.0D+00 ) / x * si(k)
  end do

  return
end subroutine sphi




  !*****************************************************************************80
  !
  !! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) SJ(0:N), the values of jn(x).
  !
  !    Output, real(8) DJ(0:N), the values of jn'(x).
  !
subroutine sphj ( n, x, nm, sj, dj )
  !SPHJ computes spherical Bessel functions jn(x) and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param sj: out: the values of jn(x)
  !@param dj: out: the values of jn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: cs
  real(8)    :: dj(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: k
  integer(4) :: m
  ! integer(4) msta1
  ! integer(4) msta2
  integer(4) :: nm
  real(8)    :: sa
  real(8)    :: sb
  real(8)    :: sj(0:n)
  real(8)    :: x

  nm = n

  if ( abs ( x ) <= 1.0D-100 ) then
     do k = 0, n
        sj(k) = 0.0D+00
        dj(k) = 0.0D+00
     end do
     sj(0) = 1.0D+00
     dj(1) = 0.3333333333333333D+00
     return
  end if

  sj(0) = sin ( x ) / x
  sj(1) = ( sj(0) - cos ( x ) ) / x

  if ( 2 <= n ) then

     sa = sj(0)
     sb = sj(1)
     m = msta1 ( x, 200 )
     if ( m < n ) then
        nm = m
     else
        m = msta2 ( x, n, 15 )
     end if

     f0 = 0.0D+00
     f1 = 1.0D+00-100
     do k = m, 0, -1
        f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
        if ( k <= nm ) then
           sj(k) = f
        end if
        f0 = f1
        f1 = f
     end do

     if ( abs ( sa ) <= abs ( sb ) ) then
        cs = sb / f0
     else
        cs = sa / f
     end if

     do k = 0, nm
        sj(k) = cs * sj(k)
     end do

  end if

  dj(0) = ( cos(x) - sin(x) / x ) / x
  do k = 1, nm
     dj(k) = sj(k-1) - ( k + 1.0D+00 ) * sj(k) / x
  end do

  return
end subroutine sphj




  !*****************************************************************************80
  !
  !! SPHK computes modified spherical Bessel functions kn(x) and derivatives.
  !
  !  Discussion:
  !
  !    This procedure computes modified spherical Bessel functions
  !    of the second kind, kn(x) and kn'(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) SK(0:N), DK(0:N), the values of kn(x) and kn'(x).
  !
subroutine sphk ( n, x, nm, sk, dk )
  !SPHK computes modified spherical Bessel functions kn(x) and derivatives.
  !
  !This procedure computes modified spherical Bessel functions
  !of the second kind, kn(x) and kn'(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param sk: out: kn(x)
  !@param dk: out: kn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: dk(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: k
  integer(4) :: nm
  real(8)    :: sk(0:n)
  real(8)    :: pi
  real(8)    :: x

  pi = 3.141592653589793D+00
  nm = n
  if ( x < 1.0D-60 ) then
     do k = 0,n
        sk(k) = 1.0D+300
        dk(k) = -1.0D+300
     end do
     return
  end if

  sk(0) = 0.5D+00 * pi / x * exp ( - x )
  sk(1) = sk(0) * ( 1.0D+00 + 1.0D+00 / x )
  f0 = sk(0)
  f1 = sk(1)
  do k = 2, n
     f = ( 2.0D+00 * k - 1.0D+00 ) * f1 / x + f0
     sk(k) = f
     if ( 1.0D+300 < abs ( f ) ) then
        exit
     end if
     f0 = f1
     f1 = f
  end do

  nm = k - 1

  dk(0) = -sk(1)
  do k = 1, nm
     dk(k) = -sk(k-1) - ( k + 1.0D+00 ) / x * sk(k)
  end do

  return
end subroutine sphk




  !*****************************************************************************80
  !
  !! SPHY computes spherical Bessel functions yn(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, integer(4) N, the order.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, integer(4) NM, the highest order computed.
  !
  !    Output, real(8) SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
  ! 
subroutine sphy ( n, x, nm, sy, dy )
  !SPHY computes spherical Bessel functions yn(x) and their derivatives.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param n: in: the order
  !@param x: in: the argument
  !@param nm: out: the highest order computed
  !@param sy: out: yn(x)
  !@param dy: out: yn'(x)
  !
  implicit none

  integer(4) :: n

  real(8)    :: dy(0:n)
  real(8)    :: f
  real(8)    :: f0
  real(8)    :: f1
  integer(4) :: k
  integer(4) :: nm
  real(8)    :: sy(0:n)
  real(8)    :: x

  nm = n

  if ( x < 1.0D-60 ) then
     do k = 0, n
        sy(k) = -1.0D+300
        dy(k) = 1.0D+300
     end do
     return
  end if

  sy(0) = - cos ( x ) / x
  sy(1) = ( sy(0) - sin ( x ) ) / x
  f0 = sy(0)
  f1 = sy(1)
  do k = 2, n
     f = ( 2.0D+00 * k - 1.0D+00 ) * f1 / x - f0
     sy(k) = f
     if ( 1.0D+300 <= abs ( f ) ) then
        exit
     end if
     f0 = f1
     f1 = f
  end do

  nm = k - 1
  dy(0) = ( sin ( x ) + cos ( x ) / x ) / x
  do k = 1, nm
     dy(k) = sy(k-1) - ( k + 1.0D+00 ) * sy(k) / x
  end do

  return
end subroutine sphy




  !*****************************************************************************80
  !
  !! STVH0 computes the Struve function H0(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) SH0, the value of H0(x).
  !
subroutine stvh0 ( x, sh0 )
  !STVH0 computes the Struve function H0(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param sh0: out: the value of H0(x)
  !
  implicit none

  real(8)    :: a0
  real(8)    :: by0
  integer(4) :: k
  integer(4) :: km
  real(8)    :: p0
  real(8)    :: pi
  real(8)    :: q0
  real(8)    :: r
  real(8)    :: s
  real(8)    :: sh0
  real(8)    :: t
  real(8)    :: t2
  real(8)    :: ta0
  real(8)    :: x

  pi = 3.141592653589793D+00
  s = 1.0D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then
     a0 = 2.0D+00 * x / pi
     do k = 1, 60
        r = - r * x / ( 2.0D+00 * k + 1.0D+00 ) * x &
             / ( 2.0D+00 * k + 1.0D+00 )
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     sh0 = a0 * s

  else

     if ( x < 50.0D+00 ) then
        km = int ( 0.5D+00 * ( x + 1.0D+00 ) )
     else
        km = 25
     end if

     do k = 1, km
        r = - r * ( ( 2.0D+00 * k - 1.0D+00 ) / x ) ** 2
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     t = 4.0D+00 / x
     t2 = t * t

     p0 = (((( &
          - 0.37043D-05     * t2 &
          + 0.173565D-04 )  * t2 &
          - 0.487613D-04 )  * t2 &
          + 0.17343D-03 )   * t2 &
          - 0.1753062D-02 ) * t2 &
          + 0.3989422793D+00

     q0 = t * ((((( &
          0.32312D-05     * t2 &
          - 0.142078D-04 )  * t2 &
          + 0.342468D-04 )  * t2 &
          - 0.869791D-04 )  * t2 &
          + 0.4564324D-03 ) * t2 &
          - 0.0124669441D+00 )

     ta0 = x - 0.25D+00 * pi
     by0 = 2.0D+00 / sqrt ( x ) &
          * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )
     sh0 = 2.0D+00 / ( pi * x ) * s + by0

  end if

  return
end subroutine stvh0




  !*****************************************************************************80
  !
  !! STVH1 computes the Struve function H1(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) SH1, the value of H1(x).
  !
subroutine stvh1 ( x, sh1 )
  !STVH1 computes the Struve function H1(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param sh1: out: the value of H1(x)
  !
  implicit none

  real(8)    :: a0
  real(8)    :: by1
  integer(4) :: k
  integer(4) :: km
  real(8)    :: p1
  real(8)    :: pi
  real(8)    :: q1
  real(8)    :: r
  real(8)    :: s
  real(8)    :: sh1
  real(8)    :: t
  real(8)    :: t2
  real(8)    :: ta1
  real(8)    :: x

  pi = 3.141592653589793D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then

     s = 0.0D+00
     a0 = - 2.0D+00 / pi
     do k = 1, 60
        r = - r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     sh1 = a0 * s

  else

     s = 1.0D+00

     if ( x <= 50.0D+00 ) then
        km = int ( 0.5D+00 * x )
     else
        km = 25
     end if

     do k = 1, km
        r = - r * ( 4.0D+00 * k * k - 1.0D+00 ) / ( x * x )
        s = s + r
        if ( abs ( r ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     t = 4.0D+00 / x
     t2 = t * t

     p1 = (((( &
          0.42414D-05      * t2 &
          - 0.20092d-04 )    * t2 &
          + 0.580759D-04 )   * t2 &
          - 0.223203D-03 )   * t2 &
          + 0.29218256D-02 ) * t2 &
          + 0.3989422819D+00

     q1 = t * ((((( &
          - 0.36594D-05     * t2 &
          + 0.1622D-04 )    * t2 &
          - 0.398708D-04 )  * t2 &
          + 0.1064741D-03 ) * t2 &
          - 0.63904D-03 )   * t2 &
          + 0.0374008364D+00 )

     ta1 = x - 0.75D+00 * pi
     by1 = 2.0D+00 / sqrt ( x ) * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )
     sh1 = 2.0D+00 / pi * ( 1.0D+00 + s / ( x * x ) ) + by1

  end if

  return
end subroutine stvh1




  !*****************************************************************************80
  !
  !! STVHV computes the Struve function Hv(x) with arbitrary order v.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    24 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of the function.
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) HV, the value of Hv(x).
  !
subroutine stvhv ( v, x, hv )
  !STVHV computes the Struve function Hv(x) with arbitrary order v.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of the function
  !@param x: in: the argument
  !@param hv: out: the value of Hv(x)
  !
  implicit none

  real(8)    :: bf
  real(8)    :: bf0
  real(8)    :: bf1
  real(8)    :: by0
  real(8)    :: by1
  real(8)    :: byv
  real(8)    :: ga
  real(8)    :: gb
  real(8)    :: hv
  integer(4) :: k
  integer(4) :: l
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: pu0
  real(8)    :: pu1
  real(8)    :: qu0
  real(8)    :: qu1
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: sa
  real(8)    :: sr
  real(8)    :: t0
  real(8)    :: t1
  real(8)    :: u
  real(8)    :: u0
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: va
  real(8)    :: vb
  real(8)    :: vt
  real(8)    :: x

  pi = 3.141592653589793D+00

  if ( x .eq. 0.0D+00 ) then
     if ( -1.0D+00 < v .or. int ( v ) - v == 0.5D+00 ) then
        hv = 0.0D+00
     else if ( v < -1.0D+00 ) then
        hv = ( -1.0D+00 ) ** ( int ( 0.5D+00 - v ) - 1 ) * 1.0D+300
     else if ( v == -1.0D+00 ) then
        hv = 2.0D+00 / pi
     end if
     return
  end if

  if ( x <= 20.0D+00 ) then

     v0 = v + 1.5D+00
     call gammaf ( v0, ga )
     s = 2.0D+00 / ( sqrt ( pi ) * ga )
     r1 = 1.0D+00

     do k = 1, 100
        va = k + 1.5D+00
        call gammaf ( va, ga )
        vb = v + k + 1.5D+00
        call gammaf ( vb, gb )
        r1 = -r1 * ( 0.5D+00 * x ) ** 2
        r2 = r1 / ( ga * gb )
        s = s + r2
        if ( abs ( r2 ) < abs ( s ) * 1.0D-12 ) then
           exit
        end if
     end do

     hv = ( 0.5D+00 * x ) ** ( v + 1.0D+00 ) * s

  else

     sa = ( 0.5D+00 * x ) ** ( v - 1.0D+00 ) / pi
     v0 = v + 0.5D+00
     call gammaf ( v0, ga )
     s = sqrt ( pi ) / ga
     r1 = 1.0D+00

     do k = 1, 12
        va = k + 0.5D+00
        call gammaf ( va, ga )
        vb = - k + v + 0.5D+00
        call gammaf ( vb, gb )
        r1 = r1 / ( 0.5D+00 * x ) ** 2
        s = s + r1 * ga / gb
     end do

     s0 = sa * s
     u = abs ( v )
     n = int ( u )
     u0 = u - n

     do l = 0, 1

        vt = 4.0D+00 * ( u0 + l ) ** 2
        r1 = 1.0D+00
        pu1 = 1.0D+00
        do k = 1, 12
           r1 = -0.0078125D+00 * r1 &
                * ( vt - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) &
                * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                / ( ( 2.0D+00 * k - 1.0D+00 ) * k * x * x )
           pu1 = pu1 + r1
        end do

        qu1 = 1.0D+00
        r2 = 1.0D+00
        do k = 1, 12
           r2 = -0.0078125D+00 * r2 &
                * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) &
                * ( vt - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) &
                / ( ( 2.0D+00 * k + 1.0D+00 ) * k * x * x )
           qu1 = qu1 + r2
        end do
        qu1 = 0.125D+00 * ( vt - 1.0D+00 ) / x * qu1

        if ( l == 0 ) then
           pu0 = pu1
           qu0 = qu1
        end if

     end do

     t0 = x - ( 0.5D+00 * u0 + 0.25D+00 ) * pi
     t1 = x - ( 0.5D+00 * u0 + 0.75D+00 ) * pi
     sr = sqrt ( 2.0D+00 / ( pi * x ) )
     by0 = sr * ( pu0 * sin ( t0 ) + qu0 * cos ( t0 ) )
     by1 = sr * ( pu1 * sin ( t1 ) + qu1 * cos ( t1 ) )
     bf0 = by0
     bf1 = by1
     do k = 2, n
        bf = 2.0D+00 * ( k - 1.0D+00 + u0 ) / x * bf1 - bf0
        bf0 = bf1
        bf1 = bf
     end do

     if ( n == 0 ) then
        byv = by0
     else if ( n == 1 ) then
        byv = by1
     else
        byv = bf
     end if
     hv = byv + s0
  end if

  return
end subroutine stvhv




  !*****************************************************************************80
  !
  !! STVL0 computes the modified Struve function L0(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    22 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) SL0, the function value.
  !
subroutine stvl0 ( x, sl0 )
  !STVL0 computes the modified Struve function L0(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param sl0: out: the function value
  !
  implicit none

  real(8)    :: a0
  real(8)    :: a1
  real(8)    :: bi0
  integer(4) :: k
  integer(4) :: km
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: s
  real(8)    :: sl0
  real(8)    :: x

  pi = 3.141592653589793D+00
  s = 1.0D+00
  r = 1.0D+00

  if ( x <= 20.0D+00 ) then

     a0 = 2.0D+00 * x / pi

     do k = 1, 60
        r = r * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     sl0 = a0 * s

  else

     if ( x < 50.0D+00 ) then
        km = int ( 0.5D+00 * ( x + 1.0D+00 ) )
     else
        km = 25
     end if

     do k = 1, km
        r = r * ( ( 2.0D+00 * k - 1.0D+00 ) / x ) ** 2
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     a1 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
     r = 1.0D+00
     bi0 = 1.0D+00
     do k = 1, 16
        r = 0.125D+00 * r * ( 2.0D+00 * k - 1.0D+00 ) ** 2 / ( k * x )
        bi0 = bi0 + r
        if ( abs ( r / bi0 ) < 1.0D-12 ) then
           exit
        end if
     end do

     bi0 = a1 * bi0
     sl0 = - 2.0D+00 / ( pi * x ) * s + bi0

  end if

  return
end subroutine stvl0




  !*****************************************************************************80
  !
  !! STVL1 computes the modified Struve function L1(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    05 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Output, real(8) SL1, the function value.
  !
subroutine stvl1 ( x, sl1 )
  !STVL1 computes the modified Struve function L1(x).
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param x: in: the argument
  !@param sl1: out: the function value
  !
  implicit none

  real(8)    :: a1
  real(8)    :: bi1
  integer(4) :: k
  integer(4) :: km
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: s
  real(8)    :: sl1
  real(8)    :: x

  pi = 3.141592653589793D+00
  r = 1.0D+00
  if ( x <= 20.0D+00 ) then
     s = 0.0D+00
     do k = 1, 60
        r = r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     sl1 = 2.0D+00 / pi * s

  else

     s = 1.0D+00
     km = int ( 0.50D+00 * x )
     km = min ( km, 25 )

     do k = 1, km
        r = r * ( 2.0D+00 * k + 3.0D+00 ) &
             * ( 2.0D+00 * k + 1.0D+00 ) / ( x * x )
        s = s + r
        if ( abs ( r / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     sl1 = 2.0D+00 / pi * ( -1.0D+00 + 1.0D+00 &
          / ( x * x ) + 3.0D+00 * s / x**4 )
     a1 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
     r = 1.0D+00
     bi1 = 1.0D+00
     do k = 1, 16
        r = -0.125D+00 * r &
             * ( 4.0D+00 - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
        bi1 = bi1 + r
        if ( abs ( r / bi1 ) < 1.0D-12 ) then
           exit
        end if
     end do

     sl1 = sl1 + a1 * bi1

  end if

  return
end subroutine stvl1




  !*****************************************************************************80
  !
  !! STVLV computes the modified Struve function Lv(x) with arbitary order.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    04 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) V, the order of Lv(x).
  !
  !    Input, real(8) X, the argument of Lv(x).
  !
  !    Output, real(8) SLV, the value of Lv(x).
  !
subroutine stvlv ( v, x, slv )
  !STVLV computes the modified Struve function Lv(x) with arbitary order.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param v: in: the order of Lv(x)
  !@param x: in: the argument of Lv(x)
  !@param slv: out: the value of Lv(x)
  !
  implicit none

  real(8)    :: bf
  real(8)    :: bf0
  real(8)    :: bf1
  real(8)    :: biv
  real(8)    :: biv0
  real(8)    :: ga
  real(8)    :: gb
  integer(4) :: k
  integer(4) :: l
  integer(4) :: n
  real(8)    :: pi
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: r2
  real(8)    :: s
  real(8)    :: s0
  real(8)    :: sa
  real(8)    :: slv
  real(8)    :: u
  real(8)    :: u0
  real(8)    :: v
  real(8)    :: v0
  real(8)    :: va
  real(8)    :: vb
  real(8)    :: vt
  real(8)    :: x

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

     if ( -1.0D+00 < v .or. int ( v ) - v == 0.5D+00 ) then
        slv = 0.0D+00
     else if ( v < -1.0D+00 ) then
        slv = ( -1 ) ** ( int ( 0.5D+00 - v ) - 1 ) * 1.0D+300
     else if ( v == -1.0D+00 ) then
        slv = 2.0D+00 / pi
     end if

  else if ( x <= 40.0D+00 ) then

     v0 = v + 1.5D+00
     call gammaf ( v0, ga )
     s = 2.0D+00 / ( sqrt ( pi ) * ga )
     r1 = 1.0D+00
     do k = 1, 100
        va = k + 1.5D+00
        call gammaf ( va, ga )
        vb = v + k + 1.5D+00
        call gammaf ( vb, gb )
        r1 = r1 * ( 0.5D+00 * x ) ** 2
        r2 = r1 / ( ga * gb )
        s = s + r2
        if ( abs ( r2 / s ) < 1.0D-12 ) then
           exit
        end if
     end do

     slv = ( 0.5D+00 * x ) ** ( v + 1.0D+00 ) * s

  else

     sa = -1.0D+00 / pi * ( 0.5D+00 * x ) ** ( v - 1.0D+00 )
     v0 = v + 0.5D+00
     call gammaf ( v0, ga )
     s = - sqrt ( pi ) / ga
     r1 = -1.0D+00
     do k = 1, 12
        va = k + 0.5D+00
        call gammaf ( va, ga )
        vb = - k + v + 0.5D+00
        call gammaf ( vb, gb )
        r1 = - r1 / ( 0.5D+00 * x ) ** 2
        s = s + r1 * ga / gb
     end do
     s0 = sa * s
     u = abs ( v )
     n = int ( u )
     u0 = u - n
     do l = 0, 1
        vt = u0 + l
        r = 1.0D+00
        biv = 1.0D+00
        do k = 1, 16
           r = -0.125D+00 * r * ( 4.0D+00 * vt * vt - &
                ( 2.0D+00 * k - 1.0D+00 )**2 ) / ( k * x )
           biv = biv + r
           if ( abs ( r / biv ) < 1.0D-12 ) then
              exit
           end if
        end do

        if ( l == 0 ) then
           biv0 = biv
        end if

     end do

     bf0 = biv0
     bf1 = biv
     do k = 2, n
        bf = - 2.0D+00 * ( k - 1.0D+00 + u0 ) / x * bf1 + bf0
        bf0 = bf1
        bf1 = bf
     end do

     if ( n == 0 ) then
        biv = biv0
     else if ( 1 < n ) then
        biv = bf
     end if

     slv = exp ( x ) / sqrt ( 2.0D+00 * pi * x ) * biv + s0

  end if

  return
end subroutine stvlv




  !*****************************************************************************80
  !
  !! VVLA computes parabolic cylinder function Vv(x) for large arguments.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    04 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) VA, the order nu.
  !
  !    Output, real(8) PV, the value of V(nu,x).
  !
subroutine vvla ( va, x, pv )
  !VVLA computes parabolic cylinder function Vv(x) for large arguments.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param va: in: the order nu
  !@param x: in: the argument
  !@param pv: out: the value of V(nu,x)
  !
  implicit none

  real(8)    :: a0
  real(8)    :: dsl
  real(8)    :: eps
  real(8)    :: gl
  integer(4) :: k
  real(8)    :: pdl
  real(8)    :: pi
  real(8)    :: pv
  real(8)    :: qe
  real(8)    :: r
  real(8)    :: va
  real(8)    :: x
  real(8)    :: x1

  pi = 3.141592653589793D+00
  eps = 1.0D-12
  qe = exp ( 0.25D+00 * x * x )
  a0 = abs ( x ) ** ( -va - 1.0D+00 ) * sqrt ( 2.0D+00 / pi ) * qe

  r = 1.0D+00
  pv = 1.0D+00
  do k = 1, 18
     r = 0.5D+00 * r * ( 2.0D+00 * k + va - 1.0D+00 ) &
          * ( 2.0D+00 * k + va ) / ( k * x * x )
     pv = pv + r
     if ( abs ( r / pv ) < eps ) then
        exit
     end if
  end do

  pv = a0 * pv

  if ( x < 0.0D+00 ) then
     x1 = -x
     call dvla ( va, x1, pdl )
     call gammaf ( -va, gl )
     dsl = sin ( pi * va ) * sin ( pi * va )
     pv = dsl * gl / pi * pdl - cos ( pi * va ) * pv
  end if

  return
end subroutine vvla




  !*****************************************************************************80
  !
  !! VVSA computes parabolic cylinder function V(nu,x) for small arguments.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    04 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real(8) X, the argument.
  !
  !    Input, real(8) VA, the order nu.
  !
  !    Output, real(8) PV, the value of V(nu,x).
  !
subroutine vvsa ( va, x, pv )
  !VVSA computes parabolic cylinder function V(nu,x) for small arguments.
  !
  !Original FORTRAN77 version by Shanjie Zhang and Jianming Jin.
  !Reference: Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996.
  !
  !@param va: in: the order nu
  !@param x: in: the argument
  !@param pv: out: the value of V(nu,x)
  !
  implicit none

  real(8)    :: a0
  real(8)    :: ep
  real(8)    :: eps
  real(8)    :: fac
  real(8)    :: g1
  real(8)    :: ga0
  real(8)    :: gm
  real(8)    :: gw
  integer(4) :: m
  real(8)    :: pi
  real(8)    :: pv
  real(8)    :: r
  real(8)    :: r1
  real(8)    :: sq2
  real(8)    :: sv
  real(8)    :: sv0
  real(8)    :: v1
  real(8)    :: va
  real(8)    :: va0
  real(8)    :: vb0
  real(8)    :: vm
  real(8)    :: x

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  ep = exp ( -0.25D+00 * x * x )
  va0 = 1.0D+00 + 0.5D+00 * va

  if ( x == 0.0D+00 ) then

     if ( ( va0 <= 0.0D+00 .and. va0 == int ( va0 ) ) .or. &
          va == 0.0D+00 ) then
        pv = 0.0D+00
     else
        vb0 = -0.5D+00 * va
        sv0 = sin ( va0 * pi )
        call gammaf ( va0, ga0 )
        pv = 2.0D+00 ** vb0 * sv0 / ga0
     end if

  else

     sq2 = sqrt ( 2.0D+00 )
     a0 = 2.0D+00 ** ( -0.5D+00 * va ) * ep / ( 2.0D+00 * pi )
     sv = sin ( - ( va + 0.5D+00 ) * pi )
     v1 = -0.5D+00 * va
     call gammaf ( v1, g1 )
     pv = ( sv + 1.0D+00 ) * g1
     r = 1.0D+00
     fac = 1.0D+00

     do m = 1, 250
        vm = 0.5D+00 * ( m - va )
        call gammaf ( vm, gm )
        r = r * sq2 * x / m
        fac = - fac
        gw = fac * sv + 1.0D+00
        r1 = gw * r * gm
        pv = pv + r1
        if ( abs ( r1 / pv ) < eps .and. gw /= 0.0D+00 ) then
           exit
        end if
     end do

     pv = a0 * pv

  end if

  return
end subroutine vvsa
