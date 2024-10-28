subroutine broyden1(ff,x,icheck,maxits,tol,tol1,tolmin,stpmx,noexit)
  procedure(broydn_func)               :: ff
  real(8), dimension(:), intent(inout) :: x
  logical, optional                    :: noexit
  logical, optional                    :: icheck
  integer, optional                    :: maxits
  real(8), optional                    :: tol,tol1,tolmin,stpmx
  logical                              :: icheck_
  integer                              :: maxits_=200
  real(8)                              :: tolf_=1d-8,tolmin_=1d-7,stpmx_=100d0,tol1_
  real(8),parameter                    :: eps=epsilon(x),tolx=eps
  integer                              :: i,its,k,n
  real(8)                              :: f,fold,stpmax
  real(8), dimension(size(x)), target  :: fvec
  real(8), dimension(size(x))          :: c,d,fvcold,g,p,s,t,w,xold
  real(8), dimension(size(x),size(x))  :: qt,r
  logical                              :: restrt,sing,noexit_=.false. 
  if(present(noexit))   noexit_  = noexit
  !
  if(present(MAXITS))   MAXITS_  = MAXITS  
  !
  if(present(TOL))     TOLF_    = TOL
  !
  if(present(TOLMIN))   TOLMIN_  = TOLMIN  
  !
  if(present(STPMX))    STPMX_   = STPMX
  !
  tol1_ = 0.01d0*tolf_
  if(present(tol1))tol1_ = tol1
  !
  if(associated(funcv))nullify(funcv)
  funcv=>ff
  fmin_fvecp=>fvec
  !
  n=size(x)
  f=fmin_(x)
  ! if (maxval(abs(fvec(:))) < 0.01d0*TOLF_) then
  if (maxval(abs(fvec(:))) < Tol1_ ) then
     icheck_=.false.
     if(present(icheck))icheck=icheck_
     RETURN
  end if
  stpmax=STPMX_*max(vabs(x(:)),real(n,8))
  restrt=.true.
  do its=1,MAXITS_
     if (restrt) then
        call fdjac(x,fvec,r)
        call qrdcmp(r,c,d,sing)
        if (sing) then
           print*,'singular Jacobian in broydn'
           if(noexit_)then
              noexit=.false.
              return
           else
              open(534,file="BROYDEN1_ERROR.err");write(534,*)"";close(534)
              stop
           endif
        endif
        call unit_matrix(qt)
        do k=1,n-1
           if (c(k) /= 0.0) then
              qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
                   matmul(r(k:n,k),qt(k:n,:)))/c(k)
           end if
        end do
        where (lower_triangle(n,n)) r(:,:)=0.0
        call put_diag(d(:),r(:,:))
     else
        s(:)=x(:)-xold(:)
        do i=1,n
           t(i)=dot_product(r(i,i:n),s(i:n))
        end do
        w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
        where (abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) &
             w(:)=0.0
        if (any(w(:) /= 0.0)) then
           t(:)=matmul(qt(:,:),w(:))
           s(:)=s(:)/dot_product(s,s)
           call qrupdt(r,qt,t,s)
           d(:)=get_diag(r(:,:))
           if (any(d(:) == 0.0))then
              print*,'r singular in broydn'
              if(noexit_)then
                 noexit=.false.
                 return
              else
                 open(534,file="BROYDEN1_ERROR.err");write(534,*)"";close(534)
                 stop
              endif
           endif
        end if
     end if
     p(:)=-matmul(qt(:,:),fvec(:))
     do i=1,n
        g(i)=-dot_product(r(1:i,i),p(1:i))
     end do
     xold(:)=x(:)
     fvcold(:)=fvec(:)
     fold=f
     call rsolv(r,d,p)
     call lnsrch(xold,fold,g,p,x,f,stpmax,icheck_,fmin_)
     if (maxval(abs(fvec(:))) < TOLF_) then
        icheck_=.false.
        if(present(icheck))icheck=icheck_
        RETURN
     end if
     if (icheck_) then
        if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
             1.0d0)/max(f,0.5d0*n)) < TOLMIN_) RETURN
        restrt=.true.
     else
        restrt=.false.            
        if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
             1.0d0)) < TOLX) RETURN
     end if
  end do
  print*,'MAXITS exceeded in broydn'
  if(noexit_)then
     noexit=.false.
     if(present(icheck))icheck=icheck_
     return
  else
     if(present(icheck))icheck=icheck_
     open(534,file="BROYDEN1_ERROR.err");write(534,*)"";close(534)
     stop
  endif
end subroutine broyden1


function fmin_(x) result(fmin)
  real(8), dimension(:), intent(in) :: x
  real(8)                           :: fmin
  if (.not. associated(fmin_fvecp)) then
     print*,'fmin: problem with pointer for returned values'
     stop
  endif
  fmin_fvecp=funcv(x)
  fmin=0.5d0*dot_product(fmin_fvecp,fmin_fvecp)
end function fmin_

