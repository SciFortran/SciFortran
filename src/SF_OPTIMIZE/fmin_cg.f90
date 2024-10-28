!+-------------------------------------------------------------------+
!  PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
!     Adapted by FRPRM subroutine from NumRec (10.6)
!     Given a starting point P that is a vector of length N, 
!     the Fletcher-Reeves-Polak-Ribiere minimisation is performed 
!     n a functin FUNC,using its gradient as calculated by a 
!     routine DFUNC. The convergence tolerance on the function 
!     value is input as FTOL.  
!     Returned quantities are: 
!     - P (the location of the minimum), 
!     - ITER (the number of iterations that were performed), 
!     - FRET (the minimum value of the function). 
!     The routine LINMIN is called to perform line minimisations.
!     Minimisation routines: DFPMIN, D/LINMIN, MNBRAK, D/BRENT and D/F1DIM
!     come from Numerical Recipes.
!  NOTE: this routine makes use of abstract interface to communicate 
!     with routines contained elsewhere. an easier way would be to include
!     the routines inside each of the two following fmin_cg routines. 
!+-------------------------------------------------------------------+
subroutine fmin_cg_df(p,f,df,iter,fret,ftol,itmax,istop,iverbose,err)
  procedure(cgfit_func)                :: f
  procedure(cgfit_fjac)                :: df
  real(8), dimension(:), intent(inout) :: p
  integer, intent(out)                 :: iter
  real(8), intent(out)                 :: fret
  real(8),optional                     :: ftol
  real(8)                              :: ftol_,a,b
  integer, optional                    :: itmax,istop
  integer                              :: itmax_,istop_
  integer                              :: i,its
  real(8)                              :: dgg,fp,gam,gg,err_
  real(8), dimension(size(p))          :: g,h,xi,p_prev
  logical,optional                     :: iverbose
  logical                              :: iverbose_,converged
  real(8),dimension(2),optional        :: err
  !
  if(associated(func))nullify(func) ; func=>f
  if(associated(fjac))nullify(fjac) ; fjac=>df
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iverbose_)write(*,"(A,ES9.2)")"CG: ftol updated to:",ftol
  endif
  itmax_=500
  if(present(itmax))then
     itmax_=itmax
     if(iverbose_)write(*,"(A,I5)")"CG: itmax updated to:",itmax
  endif
  istop_=0
  if(present(istop))then
     istop_=istop
     if(iverbose_)write(*,"(A,I3)")"CG: istop update to:",istop
  endif
  !
  fp=func(p)
  xi=fjac(p)
  g=-xi
  h=g
  xi=h
  do its=1,itmax_
     iter = its
     p_prev = p
     call dlinmin(p,xi,fret,ftol_,itmax_)
     a = abs(fret-fp)/(1d0+abs(fp))
     b = dot_product(p-p_prev,p-p_prev)/(1d0+dot_product(p,p))
     if(present(err))err=[a,b]
     if(iverbose_)then
        write(*,*)"Iter,F_n =",iter,fp
        write(*,"(A10)")"    gradF:"
        do i=1,size(xi)
           write(*,*)xi(i)
        enddo
        write(*,*)"A,B,ftol =",a,b,ftol_          
     endif
     select case(istop_)
     case default
        converged = (a<ftol_) .AND. (b<ftol_)
        if( converged .AND. iverbose_ )print*,"Converged with (|F_n-F_n-1|/1+|F_n|), ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",a,b
     case(1)
        converged = a<ftol_
        if( converged .AND. iverbose_ )print*,"Converged with (|F_n-F_n-1|/1+|F_n|)< ftol_:",a
     case(2)
        converged = b<ftol_
        if( converged .AND. iverbose_ )print*,"Converged with ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",b
     end select
     if( converged )return
     if( iverbose_)write(*,*)""
     fp = fret
     xi = fjac(p)
     gg=dot_product(g,g)
     dgg=dot_product(xi+g,xi)  !polak-ribiere
     !dgg=dot_product(xi,xi)   !fletcher-reeves.
     if (gg == 0.d0) return
     gam=dgg/gg
     g=-xi
     h=g+gam*h
     xi=h
  end do
  if(iverbose_)write(*,*)"CG: MaxIter",itmax_," exceeded."
  nullify(func)
  nullify(fjac)
  return
end subroutine fmin_cg_df










!NUMERICAL EVALUATION OF THE GRADIENT:
subroutine fmin_cg_f(p,f,iter,fret,ftol,itmax,istop,deps,iverbose)
  procedure(cgfit_func)                :: f
  real(8), dimension(:), intent(inout) :: p
  integer, intent(out)                 :: iter
  real(8), intent(out)                 :: fret
  real(8),optional                     :: ftol,deps
  integer, optional                    :: itmax,istop
  real(8)                              :: ftol_,deps_,a,b
  integer                              :: itmax_,istop_
  integer                              :: i,its
  real(8)                              :: dgg,fp,gam,gg,err_
  real(8), dimension(size(p))          :: g,h,xi,p_prev
  logical,optional                     :: iverbose
  logical                              :: iverbose_,converged
  !
  !this is just to ensure that routine needing dfunc allocated
  !and properly definted will continue to work.
  if(associated(func))nullify(func) ; func=>f
  if(associated(fjac))nullify(fjac) ; fjac=>df
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  ftol_=1.d-5
  if(present(ftol))then
     ftol_=ftol
     if(iverbose_)write(*,"(A,ES9.2)")"CG: ftol updated to:",ftol
  endif
  itmax_=500
  if(present(itmax))then
     itmax_=itmax
     if(iverbose_)write(*,"(A,I5)")"CG: itmax updated to:",itmax
  endif
  istop_=0
  if(present(istop))then
     istop_=istop
     if(iverbose_)write(*,"(A,I3)")"CG: istop update to:",istop
  endif
  deps_=0d0
  if(present(deps))then
     deps_=deps
     if(iverbose_)write(*,"(A,ES9.2)")"CG: deps update to:",deps
  endif
  df_eps = deps_
  !
  fp=func(p)
  xi=fjac(p)
  g=-xi;
  h=g
  xi=h
  do its=1,itmax_
     iter=its
     p_prev = p
     call dlinmin(p,xi,fret,ftol_)
     a = abs(fret-fp)/(1d0+abs(fp))
     b = dot_product(p-p_prev,p-p_prev)/(1d0+dot_product(p,p))
     if(iverbose_)then
        write(*,*)"Iter,F_n =",iter,fp
        write(*,"(A10)")"    gradF:"
        do i=1,size(xi)
           write(*,*)xi(i)
        enddo
        write(*,*)"A,B,ftol =",a,b,ftol_
     endif
     select case(istop_)
     case default
        converged = (a<ftol_) .AND. (b<ftol_)
        if( converged )print*,"Converged with (|F_n-F_n-1|/1+|F_n|), ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",a,b
     case(1)
        converged = a<ftol_
        if( converged )print*,"Converged with (|F_n-F_n-1|/1+|F_n|)< ftol_:",a
     case(2)
        converged = b<ftol_
        if( converged )print*,"Converged with ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",b
     end select
     if( converged )return
     if( iverbose_)write(*,*)""
     fp=fret
     xi = fjac(p)
     gg=dot_product(g,g)
     dgg=dot_product(xi+g,xi)  !polak-ribiere
     if (gg == 0d0) return
     gam=dgg/gg
     g=-xi
     h=g+gam*h
     xi=h
  end do
  if(iverbose_)write(*,*)"CG: MaxIter",itmax_," exceeded."
  nullify(func)
  nullify(fjac)
  return
end subroutine fmin_cg_f
!
function df(p) 
  real(8),dimension(:)       :: p
  real(8),dimension(size(p)) :: df
  df=func_func_jacobian(func,size(p),p)
end function df
