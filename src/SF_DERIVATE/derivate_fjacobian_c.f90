!
!              N x N Jacobian (df_i/dx_j for i,j=1,...,N)
!-----------------------------------------------------------------------
subroutine c_fdjac_nn_func(funcv,x,fjac,ml,mu,epsfcn)
  implicit none
  interface 
     function funcv(x)
       real(8),dimension(:),intent(in) :: x
       complex(8),dimension(size(x))   :: funcv !An external procedure which takes as input an array of :f:type:`real` :f:var:`x` 
                                                !and returns an array of :f:type:`complex`
     end function funcv
  end interface
  integer                                  ::  n
  real(8)                                  ::  x(:) !An array containing the element in the function domain 
                                                    !where the Jacobian matrix is to be calculated
  complex(8)                               ::  fvec(size(x))
  complex(8),dimension(size(x),size(x))    ::  fjac !The Jacobian matrix calculated at :f:var:`x`
  integer,optional                         ::  ml   !Lower subdiagonal limit (number of lower subdiagonals to be calculated)
  integer,optional                         ::  mu   !Upper subdiagonal limit (number of upper subdiagonals to be calculated)  
  real(8),optional                         ::  epsfcn !Step size for the numerical calculation of the Jacobian (default :code:`epsilon(REAL)`)
  integer                                  ::  ml_,mu_,msum
  real(8)                                  ::  eps,eps_
  real(8)                                  ::  epsmch
  real(8)                                  ::  h,temp
  complex(8)                               ::  wa1(size(x))
  complex(8)                               ::  wa2(size(x))
  integer                                  ::  i,j,k
  n=size(x)
  ml_ = n-1 ; if(present(ml))ml_=ml
  mu_ = n-1 ; if(present(mu))mu_=mu
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  msum = ml_ + mu_ + 1
  !  Evaluate the function
  fvec = funcv(x)
  !  Computation of dense approximate jacobian.
  if(n <= msum)then
     do j=1,n
        temp = x(j)
        h    = eps*abs(temp)
        if(h==0.d0) h = eps
        x(j) = temp + h
        wa1  = funcv(x)
        x(j) = temp
        fjac(1:n,j) = (wa1(1:n) - fvec(1:n))/h
     enddo
  else
     !  Computation of banded approximate jacobian.
     do k=1,msum
        do j=k,n,msum
           wa2(j) = x(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
           x(j) = wa2(j) + h
        end do
        wa1 = funcv(x)
        do j=k,n,msum
           x(j) = wa2(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
        enddo
        fjac(1:n,j)=dcmplx(0.d0,0.d0)
        do i=1,n
           if( (j-mu_<=i).AND.(i<=j+ml_) )then
              fjac(i,j) = ( wa1(i) - fvec(i) )/h
           end if
        end do
     end do
  end if
  return
end subroutine c_fdjac_nn_func

subroutine c_fdjac_nn_sub(funcv,x,fjac,ml,mu,epsfcn)
  implicit none
  interface 
     subroutine funcv(n,x,y)
       integer                         :: n
       real(8),dimension(n),intent(in) :: x
       complex(8),dimension(n)         :: y
     end subroutine funcv
  end interface
  integer          ::  n
  real(8)          ::  x(:)
  complex(8)       ::  fvec(size(x))
  complex(8)       ::  fjac(size(x),size(x))
  integer,optional ::  ml, mu    
  real(8),optional ::  epsfcn
  integer          ::  ml_,mu_,msum
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  complex(8)       ::  wa1(size(x))
  complex(8)       ::  wa2(size(x))
  integer          :: i,j,k
  n=size(x)
  ml_ = n-1 ; if(present(ml))ml_=ml
  mu_ = n-1 ; if(present(mu))mu_=mu
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  msum = ml_ + mu_ + 1
  !  Evaluate the function
  call funcv(n,x,fvec)
  !  Computation of dense approximate jacobian.
  if(n <= msum)then
     do j=1,n
        temp = x(j)
        h    = eps*abs(temp)
        if(h==0.d0) h = eps
        x(j) = temp + h
        call funcv(n,x,wa1)
        x(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) )/h
     enddo
  else
     !  Computation of banded approximate jacobian.
     do k=1,msum
        do j=k,n,msum
           wa2(j) = x(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
           x(j) = wa2(j) + h
        end do
        call funcv(n,x,wa1)
        do j=k,n,msum
           x(j) = wa2(j)
           h = eps*abs(wa2(j))
           if(h==0.d0)h = eps
        enddo
        fjac(1:n,j)=dcmplx(0.d0,0.d0)
        do i=1,n
           if( (j-mu_<=i).AND.(i<=j+ml_) )then
              fjac(i,j) = ( wa1(i) - fvec(i) )/h
           end if
        end do
     end do
  end if
  return
end subroutine c_fdjac_nn_sub

function c_f_jac_nn_func(funcv,n,x) result(df)
  interface
     function funcv(x)
       real(8), dimension(:),intent(in) :: x
       complex(8), dimension(size(x))   :: funcv
     end function funcv
  end interface
  integer                               :: n
  real(8), dimension(n), intent(inout)  :: x
  complex(8), dimension(n,n)            :: df
  call c_fdjac_nn_func(funcv,x,df)
end function c_f_jac_nn_func

function c_f_jac_nn_sub(funcv,n,x) result(df)
  interface
     subroutine funcv(n,x,y)
       integer                          :: n
       real(8), dimension(n),intent(in) :: x
       complex(8), dimension(n)         :: y
     end subroutine funcv
  end interface
  integer                              :: n
  real(8), dimension(n), intent(inout) :: x
  complex(8), dimension(n,n)           :: df
  call c_fdjac_nn_sub(funcv,x,df)
end function c_f_jac_nn_sub






!
!              M x N Jacobian (df_i/dx_j for i=1,...,M;j=1,...,N)
!-----------------------------------------------------------------------
subroutine c_fdjac_mn_func(funcv,n,x,m,fjac,epsfcn)
  implicit none
  interface 
     function funcv(n,x,m)
       integer                         :: m  !Dimension of the codomain of :f:func:`funcv`
       integer                         :: n  !Dimension of the domain :f:var:`x`
       real(8),dimension(n),intent(in) :: x
       complex(8),dimension(m)         :: funcv
     end function funcv
  end interface
  integer          ::  n
  integer          ::  m
  real(8)          ::  x(n)
  complex(8)       ::  fvec(m)
  complex(8)       ::  fjac(m,n)
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  complex(8)       ::  wa1(m)
  complex(8)       ::  wa2(m)
  integer          :: i,j,k
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps    = sqrt(max(eps_,epsmch))
  fvec = funcv(n,x,m)
  do j=1,n
     temp = x(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x(j) = temp + h
     wa1 = funcv(n,x,m)
     x(j) = temp
     fjac(1:m,j) = (wa1(1:m) - fvec(1:m))/h
  enddo
end subroutine c_fdjac_mn_func

subroutine c_fdjac_mn_sub(funcv,n,x,m,fjac,epsfcn)
  implicit none
  interface 
     subroutine funcv(n,x,m,y)
       integer                         :: n,m
       real(8),dimension(n),intent(in) :: x
       complex(8),dimension(m)         :: y
     end subroutine funcv
  end interface
  integer          ::  n
  integer          ::  m
  real(8)          ::  x(n)
  complex(8)       ::  fvec(m)
  complex(8)       ::  fjac(m,n)
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  complex(8)       ::  wa1(m)
  complex(8)       ::  wa2(m)
  integer          :: i,j,k
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps    = sqrt(max(eps_,epsmch))
  call funcv(n,x,m,fvec)
  do j=1,n
     temp = x(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x(j) = temp + h
     call funcv(n,x,m,wa1)
     x(j) = temp
     fjac(1:m,j) = (wa1(1:m) - fvec(1:m))/h
  enddo
end subroutine c_fdjac_mn_sub

function c_f_jac_mn_func(funcv,n,x,m) result(df)
  interface 
     function funcv(n,x,m)
       integer                         :: n,m
       real(8),dimension(n),intent(in) :: x
       complex(8),dimension(m)         :: funcv
     end function funcv
  end interface
  integer                               :: n,m
  real(8), dimension(n), intent(inout)  :: x
  complex(8), dimension(m,n)            :: df
  call c_fdjac_mn_func(funcv,n,x,m,df)    
end function c_f_jac_mn_func

function c_f_jac_mn_sub(funcv,n,x,m) result(df)
  interface
     subroutine funcv(n,x,m,y)
       implicit none
       integer                          :: n,m
       real(8), dimension(n),intent(in) :: x
       complex(8), dimension(m)         :: y
     end subroutine funcv
  end interface
  integer                               :: n,m
  real(8), dimension(n), intent(inout)  :: x
  complex(8), dimension(m,n)            :: df
  call c_fdjac_mn_sub(funcv,n,x,m,df)
end function c_f_jac_mn_sub





!
!              1 x N Jacobian (df_i/dx_j for i=1;j=1,...,N)
!-----------------------------------------------------------------------
subroutine c_fdjac_1n_func(funcv,x,fjac,epsfcn)
  implicit none
  interface 
     function funcv(x)
       real(8),dimension(:),intent(in) :: x
       complex(8)                      :: funcv   !An external procedure which takes as input an array of :f:type:`real` :f:var:`x` 
                                                  ! and returns a complex number
     end function funcv
  end interface
  integer                         ::  n
  real(8),dimension(:)            ::  x!An array containing the element in the function domain where the gradient is to be calculated
  complex(8)                      ::  fvec
  complex(8),dimension(size(x))   ::  fjac !The gradient vector calculated at :f:var:`x`
  real(8),optional                ::  epsfcn !Step size for the numerical calculation of the Jacobian (default :code:`epsilon(REAL)`)
  real(8)                         ::  eps,eps_
  real(8)                         ::  epsmch
  real(8)                         ::  h,temp
  complex(8)                      ::  wa1
  complex(8)                      ::  wa2
  integer                         :: i,j,k
  n=size(x)
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  !  Evaluate the function
  fvec = funcv(x)
  do j=1,n
     temp = x(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x(j) = temp + h
     wa1  = funcv(x)
     x(j) = temp
     fjac(j) = (wa1 - fvec)/h
  enddo
end subroutine c_fdjac_1n_func

subroutine c_fdjac_1n_sub(funcv,x,fjac,epsfcn)
  interface 
     subroutine funcv(n,x,y)
       integer                         :: n
       real(8),dimension(n),intent(in) :: x
       complex(8)                      :: y
     end subroutine funcv
  end interface
  integer          ::  n
  real(8)          ::  x(:)
  complex(8)       ::  fvec
  complex(8)       ::  fjac(size(x))
  real(8),optional ::  epsfcn
  real(8)          ::  eps,eps_
  real(8)          ::  epsmch
  real(8)          ::  h,temp
  complex(8)       ::  wa1
  complex(8)       ::  wa2
  integer          :: i,j,k
  n=size(x)
  eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
  epsmch = epsilon(epsmch)
  eps  = sqrt(max(eps_,epsmch))
  !  Evaluate the function
  call funcv(n,x,fvec)
  !  Computation of dense approximate jacobian.
  do j=1,n
     temp = x(j)
     h    = eps*abs(temp)
     if(h==0.d0) h = eps
     x(j) = temp + h
     call funcv(n,x,wa1)
     x(j) = temp
     fjac(j) = (wa1-fvec)/h
  enddo
  return
end subroutine c_fdjac_1n_sub

function c_f_jac_1n_func(funcv,n,x) result(df)
  interface
     function funcv(x)
       real(8),dimension(:),intent(in) :: x
       complex(8)                      :: funcv
     end function funcv
  end interface
  integer                               :: n
  real(8), dimension(n), intent(inout)  :: x
  complex(8), dimension(n)              :: df
  call c_fdjac_1n_func(funcv,x,df)
end function c_f_jac_1n_func

function c_f_jac_1n_sub(funcv,n,x) result(df)
  interface
     subroutine funcv(n,x,y)
       integer                          :: n
       real(8), dimension(n),intent(in) :: x
       complex(8)                       :: y
     end subroutine funcv
  end interface
  integer                               :: n
  real(8), dimension(n), intent(inout)  :: x
  complex(8), dimension(n)              :: df
  call c_fdjac_1n_sub(funcv,x,df)
end function c_f_jac_1n_sub

