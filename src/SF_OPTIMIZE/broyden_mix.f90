
!X =v[n]
!Fx=v[n+1]
subroutine d_broyden_mix(X,Fx,alpha,M,iter,w0)
  real(8),intent(inout),dimension(:)      :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x))   :: Fx     ! residual, x_out-x_in
  real(8),intent(in)                      :: alpha  ! mixing parameter
  integer,intent(in)                      :: M      ! number of previous iterations kept in memory, 0 is linear mixing
  integer,intent(in)                      :: iter   ! iteration number, 1 resets the memory
  real(8),optional                        :: w0     ! optional: weight of the initial iteration, default 0.01
  real(8)                                 :: wg0
  real(8),dimension(M)                    :: wg
  integer                                 :: N
  real(8)                                 :: Xsave(size(X))
  integer                                 :: iter_used,ipos,inext
  real(8),allocatable,dimension(:,:),save :: Df,Dv
  real(8),dimension(M,M)                  :: Beta
  real(8),dimension(M)                    :: Work
  real(8)                                 :: norm,GammaMix
  integer                                 :: i,j
  N   = size(X)
  wg0 = 0.01d0 ; if(present(w0))wg0=w0
  wg  = min(M,8)*1d0
  if(iter==1)then
     if(allocated(Df))deallocate(Df)
     if(allocated(Dv))deallocate(Dv)
     allocate(Df(M,N))
     allocate(Dv(M,N))
     Df   = 0d0
     Dv   = 0d0
  endif
  !
  Xsave = X
  !
  !linear mixing if M=0
  if(M==0)then
     X=X+alpha*Fx
     return
  endif
  !
  !Get broyden mixing pointer
  iter_used = min(iter-1,M)
  ipos      = iter-1-((iter-2)/M)*M
  !
  !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
  !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
  if(iter==1)then
     X=X+alpha*Fx               !(Fx-X)
     return
  else
     Df(ipos,:) = Fx - Df(ipos,:)
     Dv(ipos,:) = X  - Dv(ipos,:)
     norm = dot_product(Df(ipos,:),Df(ipos,:))
     norm = sqrt(norm)
     Df(ipos,:)=Df(ipos,:)/norm
     Dv(ipos,:)=Dv(ipos,:)/norm
  endif
  !
  !Build Beta
  !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
  beta=0d0
  do i=1,iter_used
     do j=i+1,iter_used
        beta(i,j) = wg(i)*wg(j)*dot_product(Df(j,:),Df(i,:))
     enddo
     beta(i,i) = wg0**2 + wg(i)**2
  enddo
  !
  call inv_sym(beta(1:iter_used,1:iter_used))
  !
  !Work(i) = Df(i)*f
  do i=1,iter_used
     work(i) = dot_product(Df(i,:),Fx)
  enddo
  !
  X=X+alpha*Fx
  !
  !v^{i+1} = v^i + alpha*f - sum_i g(i)*w(i)*u(i)
  !where:
  ! - g(i) = [sum_j work(j)*b(j,i)*w(j)]
  ! - u(i) = (alpha*DeltaF(i) + DeltaV(i))
  do i=1,iter_used
     GammaMix=0d0
     do j=1,iter_used
        GammaMix=GammaMix + beta(j,i)*wg(j)*Work(j)
     enddo
     !
     X = X - wg(i)*GammaMix*(alpha*Df(i,:) + Dv(i,:))
  enddo
  !
  !Store the next entries for DeltaF and DeltaV
  inext      = iter-((iter-1)/M)*M
  Df(inext,:)= Fx
  Dv(inext,:)= Xsave
end subroutine d_broyden_mix

subroutine c_broyden_mix(X,Fx,alpha,M,iter,w0)
  complex(8),intent(inout),dimension(:)      :: x
  complex(8),intent(in),dimension(size(x))   :: Fx
  real(8),intent(in)                         :: alpha
  integer,intent(in)                         :: M
  integer,intent(in)                         :: iter
  real(8),optional                           :: w0
  real(8)                                    :: wg0
  real(8),dimension(M)                       :: wg
  integer                                    :: N
  complex(8)                                 :: Xsave(size(X))
  integer                                    :: iter_used,ipos,inext
  complex(8),allocatable,dimension(:,:),save :: Df,Dv
  complex(8),dimension(M,M)                  :: Beta
  complex(8),dimension(M)                    :: Work
  complex(8)                                 :: GammaMix
  real(8)                                    :: norm
  integer                                    :: i,j
  N   = size(X)
  wg0 = 0.01d0 ; if(present(w0))wg0=w0
  wg  = min(M,8)*1d0
  if(iter==1)then
     if(allocated(Df))deallocate(Df)
     if(allocated(Dv))deallocate(Dv)
     allocate(Df(M,N))
     allocate(Dv(M,N))
     Df   = 0d0
     Dv   = 0d0
  endif
  !
  Xsave = X
  !
  !linear mixing if M=0
  if(M==0)then
     X=X+alpha*Fx
     return
  endif
  !
  !Get broyden mixing pointer
  iter_used = min(iter-1,M)
  ipos      = iter-1-((iter-2)/M)*M
  !
  !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
  !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
  if(iter==1)then
     X=X+alpha*Fx               !(Fx-X)
     return
  else
     Df(ipos,:) = Fx - Df(ipos,:)
     Dv(ipos,:) = X  - Dv(ipos,:)
     norm = dot_product(Df(ipos,:),Df(ipos,:))
     norm = sqrt(norm)
     Df(ipos,:)=Df(ipos,:)/norm
     Dv(ipos,:)=Dv(ipos,:)/norm
  endif
  !
  !Build Beta
  !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
  beta=0d0
  do i=1,iter_used
     do j=i+1,iter_used
        beta(i,j) = wg(i)*wg(j)*dot_product(Df(j,:),Df(i,:))
     enddo
     beta(i,i) = wg0**2 + wg(i)**2
  enddo
  !
  call inv_sym(beta(1:iter_used,1:iter_used))
  !
  !Work(i) = Df(i)*f
  do i=1,iter_used
     work(i) = dot_product(Df(i,:),Fx)
  enddo
  !
  X=X+alpha*Fx
  !
  !v^{i+1} = v^i + alpha*f - sum_i g(i)*w(i)*u(i)
  !where:
  ! - g(i) = [sum_j work(j)*b(j,i)*w(j)]
  ! - u(i) = (alpha*DeltaF(i) + DeltaV(i))
  do i=1,iter_used
     GammaMix=0d0
     do j=1,iter_used
        GammaMix=GammaMix + beta(j,i)*wg(j)*Work(j)
     enddo
     !
     X = X - wg(i)*GammaMix*(alpha*Df(i,:) + Dv(i,:))
  enddo
  !
  !Store the next entries for DeltaF and DeltaV
  inext      = iter-((iter-1)/M)*M
  Df(inext,:)= Fx
  Dv(inext,:)= Xsave
end subroutine c_broyden_mix
