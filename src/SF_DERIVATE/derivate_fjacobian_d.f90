!
!              N x N Jacobian (df_i/dx_j for i,j=1,...,N)
!-----------------------------------------------------------------------
subroutine fdjac_nn_func(funcv,x,fjac,ml,mu,epsfcn)
  interface 
     function funcv(x)
       real(8),dimension(:),intent(in) :: x
       real(8),dimension(size(x))      :: funcv !An external procedure which takes as input an array of :f:type:`real` :f:var:`x` 
                                                !and returns an array of :f:type:`real` 
     end function funcv
  end interface
  integer                              ::  n
  real(8),intent(in)                   ::  x(:) !An array containing the element in the function domain 
                                                !where the Jacobian matrix is to be calculated
  real(8)                              ::  x_(size(x))
  real(8)                              ::  fvec(size(x))
  real(8),dimension(size(x),size(x))   ::  fjac !The Jacobian matrix calculated at :f:var:`x`
  integer,optional                     ::  ml !Lower subdiagonal limit (number of lower subdiagonals to be calculated)
  integer,optional                     ::  mu !Upper subdiagonal limit (number of upper subdiagonals to be calculated)
  real(8),optional                     ::  epsfcn !Step size for the numerical calculation of the Jacobian (default :code:`epsilon(REAL)`)
  integer                              ::  ml_,mu_,msum
  real(8)                              ::  eps,eps_
  real(8)                              ::  epsmch
  real(8)                              ::  h,temp
  real(8)                              ::  wa1(size(x))
  real(8)                              ::  wa2(size(x))
  integer                              :: i,j,k
  n=size(x)
  x_ = x
  ml_ = n-1 ; if(present(ml))ml_=ml
  mu_ = n-1 ; if(present(mu))mu_=mu
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  msum = ml_ + mu_ + 1
  !  Evaluate the function
  fvec = funcv(x_)
  !  Computation of dense approximate jacobian.
  if(n <= msum)then
     do j=1,n
        temp = x_(j)
        h    = eps*abs(temp)
        if(h==0.d0) h = eps
        x_(j) = temp + h
        wa1  = funcv(x_)
        x_(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) )/h
     enddo
  else
     !  Computation of banded approximate jacobian.
     do k=1,msum
        do j=k,n,msum
           wa2(j) = x_(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
           x_(j) = wa2(j) + h
        end do
        wa1 = funcv(x_)
        do j=k,n,msum
           x_(j) = wa2(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
        enddo
        fjac(1:n,j)=0.d0
        do i=1,n
           if( (j-mu_<=i).AND.(i<=j+ml_) )then
              fjac(i,j) = ( wa1(i) - fvec(i) )/h
           end if
        end do
     end do
  end if
  return
end subroutine fdjac_nn_func

subroutine fdjac_nn_sub(funcv,x,fjac,ml,mu,epsfcn)
  interface 
     subroutine funcv(x,y)
       real(8),dimension(:),intent(in) :: x
       real(8),dimension(size(x))      :: y
     end subroutine funcv
  end interface
  integer          ::  n
  real(8),intent(in) ::  x(:)
  real(8)            ::  x_(size(x))
  real(8)          ::  fvec(size(x))
  real(8)          ::  fjac(size(x),size(x)) !The Jacobian matrix calculated at :f:var:`x`
  integer,optional ::  ml, mu    
  real(8),optional ::  epsfcn
  integer          ::  ml_,mu_,msum
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  real(8)          ::  wa1(size(x))
  real(8)          ::  wa2(size(x))
  integer          :: i,j,k
  n=size(x)
  x_ = x
  ml_ = n-1 ; if(present(ml))ml_=ml
  mu_ = n-1 ; if(present(mu))mu_=mu
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  msum = ml_ + mu_ + 1
  !  Evaluate the function
  call funcv(x_,fvec)
  !  Computation of dense approximate jacobian.
  if(n <= msum)then
     do j=1,n
        temp = x_(j)
        h    = eps*abs(temp)
        if(h==0.d0) h = eps
        x_(j) = temp + h
        call funcv(x_,wa1)
        x_(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) )/h
     enddo
  else
     !  Computation of banded approximate jacobian.
     do k=1,msum
        do j=k,n,msum
           wa2(j) = x_(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
           x_(j) = wa2(j) + h
        end do
        call funcv(x_,wa1)
        do j=k,n,msum
           x_(j) = wa2(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
        enddo
        fjac(1:n,j)=0.d0
        do i=1,n
           if( (j-mu_<=i).AND.(i<=j+ml_) )then
              fjac(i,j) = ( wa1(i) - fvec(i) )/h
           end if
        end do
     end do
  end if
  return
end subroutine fdjac_nn_sub

function f_jac_nn_func(funcv,x) result(df)
  interface
     function funcv(x)
       real(8), dimension(:),intent(in) :: x     !An array containing the element in the function domain where the Jacobian is to be calculated
       real(8), dimension(size(x))      :: funcv !The function of subroutine to be differentiated
     end function funcv
  end interface
  real(8),intent(in)                    ::  x(:) !An array containing the element in the function domain where the Jacobian is to be calculated
  real(8),dimension(size(x),size(x))    :: df    !The Jacobian matrix
  call fdjac_nn_func(funcv,x,df)
end function f_jac_nn_func

function f_jac_nn_sub(funcv,x) result(df)
  interface
     subroutine funcv(x,y)
       real(8), dimension(:),intent(in) :: x
       real(8), dimension(size(x))      :: y
     end subroutine funcv
  end interface
  real(8),dimension(:),intent(in)     :: x  !An array containing the element in the function domain where the Jacobian is to be calculated
  real(8),dimension(size(x),size(x))  :: df !The Jacobian matrix
  call fdjac_nn_sub(funcv,x,df)
end function f_jac_nn_sub






!
!              M x N Jacobian (df_i/dx_j for i=1,...,M;j=1,...,N)
!-----------------------------------------------------------------------
subroutine fdjac_mn_func(funcv,x,m,fjac,epsfcn)
  implicit none
  interface 
     function funcv(x,m)
       real(8),dimension(:),intent(in) :: x
       integer                         :: m !Dimension of the codomain of :f:func:`funcv`
       real(8),dimension(m)            :: funcv
     end function funcv
  end interface
  integer          ::  n
  integer          ::  m
  real(8),intent(in) ::  x(:)
  real(8)            ::  x_(size(x))
  real(8)          ::  fvec(m)
  real(8)          ::  fjac(m,size(x)) !The Jacobian matrix calculated at :f:var:`x`
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  real(8)          ::  wa1(m)
  real(8)          ::  wa2(m)
  integer          :: i,j,k
  n = size(x)
  x_ = x
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps    = sqrt(max(eps_,epsmch))
  fvec = funcv(x_,m)
  do j=1,n
     temp = x_(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x_(j) = temp + h
     wa1 = funcv(x_,m)
     x_(j) = temp
     fjac(1:m,j) = (wa1(1:m) - fvec(1:m))/h
  enddo
end subroutine fdjac_mn_func

subroutine fdjac_mn_sub(funcv,x,m,fjac,epsfcn)
  implicit none
  interface 
     subroutine funcv(x,m,y)
       implicit none
       integer                         :: m
       real(8),dimension(:),intent(in) :: x
       real(8),dimension(m)            :: y
     end subroutine funcv
  end interface
  integer          ::  n
  integer          ::  m
  real(8),intent(in) ::  x(:)
  real(8)            ::  x_(size(x))
  real(8)          ::  fvec(m)
  real(8)          ::  fjac(m,size(x)) !The Jacobian matrix calculated at :f:var:`x`
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  real(8)          ::  wa1(m)
  real(8)          ::  wa2(m)
  integer          :: i,j,k
  n=size(x)
  x_ = x
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps    = sqrt(max(eps_,epsmch))
  call funcv(x_,m,fvec)
  do j=1,n
     temp = x_(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x_(j) = temp + h
     call funcv(x_,m,wa1)
     x_(j) = temp
     fjac(1:m,j) = (wa1(1:m) - fvec(1:m))/h
  enddo
end subroutine fdjac_mn_sub

function f_jac_mn_func(funcv,x,m) result(df)
  interface 
     function funcv(x,m)
       real(8),dimension(:),intent(in) :: x 
       integer                         :: m !Dimension of the codomain of :f:func:`funcv` (optional)
       real(8),dimension(m)            :: funcv
     end function funcv
  end interface
  integer                           :: m  !Dimension of the codomain of :f:func:`funcv` (optional)
  integer                           :: n
  real(8),dimension(:),intent(in)   :: x  !An array containing the element in the function domain where the Jacobian is to be calculated
  real(8),dimension(m,size(x))      :: df !The Jacobian matrix
  call fdjac_mn_func(funcv,x,m,df)
end function f_jac_mn_func

function f_jac_mn_sub(funcv,x,m) result(df)
  interface
     subroutine funcv(x,m,y)
       implicit none
       integer                          :: m !Dimension of the codomain of :f:func:`funcv` (optional)
       real(8), dimension(:),intent(in) :: x
       real(8), dimension(m)            :: y
     end subroutine funcv
  end interface
  integer                           :: m  !Dimension of the codomain of :f:func:`funcv` (optional)
  real(8),dimension(:),intent(in)   :: x  !An array containing the element in the function domain where the Jacobian is to be calculated
  real(8),dimension(m,size(x))      :: df !The Jacobian matrix
  call fdjac_mn_sub(funcv,x,m,df)
end function f_jac_mn_sub





!
!              1 x N Jacobian (df_i/dx_j for i=1;j=1,...,N)
!-----------------------------------------------------------------------
subroutine fdjac_1n_func(funcv,x,fjac,epsfcn)
  implicit none
  interface 
     function funcv(x)
       implicit none
       real(8),dimension(:),intent(in) :: x
       real(8)                         :: funcv !An external procedure which takes as input an array of :f:type:`real` :f:var:`x` 
                                                !and returns a real number
     end function funcv
  end interface
  integer                         ::  n
  real(8),intent(in),dimension(:) ::  x !An array containing the element in the function domain where the gradient is to be calculated
  real(8)                         ::  x_(size(x))
  real(8)                         ::  fvec
  real(8),dimension(size(x))      ::  fjac !The gradient vector calculated at :f:var:`x`
  real(8),optional   ::  epsfcn !Step size for the numerical calculation of the Jacobian (default :code:`epsilon(REAL)`)
  real(8)            ::  eps,eps_
  real(8)            ::  epsmch
  real(8)            ::  h,temp
  real(8)            ::  wa1
  real(8)            ::  wa2
  integer            :: i,j,k
  n=size(x)
  x_ = x
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  !  Evaluate the function
  fvec = funcv(x_)
  do j=1,n
     temp = x_(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x_(j) = temp + h
     wa1  = funcv(x_)
     x_(j) = temp
     fjac(j) = (wa1 - fvec)/h
  enddo
end subroutine fdjac_1n_func

subroutine fdjac_1n_sub(funcv,x,fjac,epsfcn)
  interface 
     subroutine funcv(x,y)
       real(8),dimension(:),intent(in) :: x
       real(8)                         :: y
     end subroutine funcv
  end interface
  integer          ::  n
  real(8),intent(in) ::  x(:)
  real(8)            ::  x_(size(x))
  real(8)          ::  fvec
  real(8)          ::  fjac(size(x)) 
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  real(8)          ::  wa1
  real(8)          ::  wa2
  integer          :: i,j,k
  n=size(x)
  x_ = x
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  !  Evaluate the function
  call funcv(x_,fvec)
  !  Computation of dense approximate jacobian.
  do j=1,n
     temp = x_(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x_(j) = temp + h
     call funcv(x_,wa1)
     x_(j) = temp
     fjac(j) = (wa1-fvec)/h
  enddo
  return
end subroutine fdjac_1n_sub

function f_jac_1n_func(funcv,x) result(df)
  interface
     function funcv(x)
       real(8),dimension(:),intent(in) :: x
       real(8)                         :: funcv !An external procedure which takes as input an array of :f:type:`real` :f:var:`x`
                                                ! and returns a real number
     end function funcv
  end interface
  real(8), dimension(:), intent(in) :: x   !An array containing the element in the function domain where the gradient is to be calculated
  real(8), dimension(size(x))       :: df  !The gradient array
  call fdjac_1n_func(funcv,x,df)
end function f_jac_1n_func

function f_jac_1n_sub(funcv,x) result(df)
  interface
     subroutine funcv(x,y)
       real(8), dimension(:),intent(in) :: x
       real(8)                          :: y
     end subroutine funcv
  end interface
  real(8), dimension(:), intent(in) :: x   !An array containing the element in the function domain where the gradient is to be calculated
  real(8), dimension(size(x))       :: df
  call fdjac_1n_sub(funcv,x,df)
end function f_jac_1n_sub

