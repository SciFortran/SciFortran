!+-----------------------------------------------------------------+
!PURPOSE  : Build the BETHE Lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_lattice(dos,ome,Lk,D)
  !This subroutine builds the Bethe lattice structure of the problem: the energy grid :f:var:`ome` of :f:var:`Lk` equally spaced
  !points from :math:`-D` to :math:`D`, and the density of states times the energy step on that grid,
  !:code:`dos(i) = dens_bethe(ome(i),D)*de`, with :code:`de = 2D/(Lk-1)`. The density of states :code:`dos/de` is also written,
  !as two columns energy and density, in the file :code:`DOSbethe.lattice`, and a message is printed. Note that the energy step
  !is computed with :f:var:`D` itself, so :f:var:`D` has to be always passed, although it is declared optional.
  !
  real(8)          :: dos(Lk),ome(Lk)
  integer          :: Lk  !number of energy points
  real(8),optional :: D   !half bandwidth (pass it explicitly)
  integer          :: ie
  real(8)          :: de,e,D_
  complex(8)       :: gf,zeta
  D_=1.d0;if(present(D))D_=D
  de= 2.d0*D/dble(Lk-1)
  write(*,"(A,I8,A)")"Bethe Lattice with:",Lk," e-points"
  open(10,file="DOSbethe.lattice")
  do ie=1,Lk
     e=-D_ + dble(ie-1)*de
     dos(ie)=dens_bethe(e,D_)*de
     ome(ie)=e
     write(10,*)e,dos(ie)/de
  enddo
  close(10)
end subroutine bethe_lattice



!+-------------------------------------------------------------------+
!purpose  : calculate the non-interacting dos for BETHE lattice 
!+-------------------------------------------------------------------+
elemental function dens_bethe(x,D)
  !This function returns the density of states of the Bethe lattice with half bandwidth :f:var:`D` (default 1) at the energy
  !:f:var:`x`: :math:`\rho(x) = \frac{2}{\pi D}\sqrt{1-(x/D)^2}`. It is zero outside the band :math:`|x|>D`. The function is
  !elemental.
  !
  real(8),intent(in)          :: x  !energy
  real(8),intent(in),optional :: d  !half bandwidth (default 1)
  real(8)                     :: dens_bethe,d_
  complex(8)                  :: root,d2
  d_=1.d0;if(present(d))d_=d
  d2=dcmplx(d_,0.d0)
  root=dcmplx((1.d0-1.d0*((x/d_))**2),0.d0)
  root=sqrt(root)
  dens_bethe=(2.d0/(3.141592653589793238d0*d_))*root
end function dens_bethe




!+------------------------------------------------------------------+
!purpose  : get the hilber transfom of a given "zeta" with bethe dos
!+------------------------------------------------------------------+
elemental function gfbethe(w,zeta,d)
  !This function returns the Hilbert transform of the Bethe density of states with half bandwidth :f:var:`d` at the complex
  !energy :f:var:`zeta`, i.e. the local Green's function :math:`G(\zeta) = 2/(\zeta + s\sqrt{\zeta^2-d^2})`. The sign :math:`s`
  !selects the branch of the square root and it is fixed by the sign of :f:var:`w` times the imaginary part of the root, with
  !:f:var:`w` the real frequency. The function is elemental.
  !
  real(8),intent(in)    :: w     !real frequency, its sign selects the branch of the root
  real(8),intent(in)    :: d     !half bandwidth
  complex(8),intent(in) :: zeta  !complex energy
  complex(8)            :: gfbethe,sqroot
  real(8)               :: sq,sig
  sqroot=sqrt(zeta**2-d**2)
  sq=dimag(sqroot)
  sig=w*sq/abs(w*sq)
  gfbethe=2.d0/(zeta+sig*sqroot)
  return
end function gfbethe


!+------------------------------------------------------------------+
!purpose  : get the hilber transfom of a given "zeta" with bethe dos
!+------------------------------------------------------------------+
function gfbether(w,zeta,d)
  !This function returns the Hilbert transform of the Bethe density of states with half bandwidth :f:var:`d` at the complex
  !energy :f:var:`zeta`, i.e. the local Green's function :math:`G(\zeta) = 2/(\zeta + s\sqrt{\zeta^2-d^2})`, as
  !:f:func_inline:`gfbethe`, but the sign :math:`s` selecting the branch of the square root is the sign of the real part of
  !:f:var:`zeta`. The argument :f:var:`w` is not used. If the real part of :f:var:`zeta` is zero it is replaced by :code:`1d-8`,
  !so :f:var:`zeta` is modified.
  !
  real(8)               :: w     !not used
  real(8)               :: d     !half bandwidth
  complex(8)            :: zeta  !complex energy, modified if its real part is 0
  complex(8)            :: gfbether,sqroot
  real(8)               :: sig
  if(dreal(zeta)==0.d0)zeta=dcmplx(1.d-8,dimag(zeta))
  sqroot=sqrt(zeta**2-d**2)
  sig=dreal(zeta)/abs(dreal(zeta))
  gfbether=2.d0/(zeta+sig*sqroot)
end function gfbether




!+-----------------------------------------------------------------+
!purpose  : build the bethe lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_guess_g0(g0,d,beta,hloc) 
  !This subroutine returns in :f:var:`g0` a guess of the non-interacting local Green's function of the Bethe lattice with half
  !bandwidth :f:var:`d` at the inverse temperature :f:var:`beta`, on the first :code:`size(g0)` positive Matsubara frequencies
  !:math:`\omega_n = \pi(2n-1)/\beta`: :code:`g0(n) = gfbethe(`:math:`\omega_n`:code:`, i`:math:`\omega_n`:code:` - hloc, d)`,
  !where :f:var:`hloc` is the local energy level.
  !
  complex(8),dimension(:) :: g0
  real(8)                 :: d     !half bandwidth
  real(8)                 :: beta  !inverse temperature
  real(8)                 :: hloc  !local energy level
  integer                 :: i
  real(8)                 :: wm
  complex(8)              :: zeta
  do i=1,size(g0)
     wm    = 3.141592653589793238d0/beta*(2*i-1)
     zeta  = dcmplx(0.d0,wm) - hloc
     g0(i) = gfbethe(wm,zeta,d)
  enddo
end subroutine bethe_guess_g0
