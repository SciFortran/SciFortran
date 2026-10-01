subroutine d_adaptive_mix(x,Fx,alpha,iter)
   real(8),intent(inout),dimension(:)      :: x      ! in: x_in; out: mixed x
   real(8),intent(in),dimension(size(x))   :: Fx     ! residual, x_out-x_in
   real(8),intent(in)                      :: alpha  ! mixing parameter, initial and minimum speed
   integer,intent(in)                      :: iter   ! iteration number, 1 resets the memory
   integer                                 :: N
   real(8),parameter                       :: alpha_max = 1d0 
   real(8),allocatable,dimension(:),save   :: Fx_prev !prev iter
   real(8),allocatable,dimension(:),save   :: beta    !mixing β
   integer                                 :: j
   logical                                 :: not_oscillating
   !
   N = size(x)
   !
   if(iter==1)then !reset persistent variables
      if(allocated(Fx_prev))deallocate(Fx_prev)
      allocate(Fx_prev(N))
      Fx_prev = 0d0
      !
      if(allocated(beta))deallocate(beta)
      allocate(beta(N))
      beta = alpha
   endif
   !
   X = X + beta*Fx ! = (1-β)x_in + βx_out (β: effective speed)
   !
   if (iter > 1) then !update effective elemental mixing
      do j=1,N
         not_oscillating = Fx_prev(j)*Fx(j)>0 !same sign
         if(not_oscillating) then !increase mixing speed
            beta(j) = beta(j) + alpha !but with a cutoff
            if (beta(j) > alpha_max) beta(j) = alpha_max
         else
            beta(j) = alpha ! reset to input mixing α
         end if
      end do
   end if
   !
   Fx_prev = Fx !store Fx_prev for future calls
   !
end subroutine d_adaptive_mix

subroutine c_adaptive_mix(z,Fz,alpha,iter)
   complex(8),intent(inout),dimension(:)    :: z       !z_in
   complex(8),intent(in),dimension(size(z)) :: Fz      !z_out-z_in
   real(8),intent(in)                       :: alpha   !mixing α
   integer,intent(in)                       :: iter
   integer                                  :: N
   real(8),parameter                        :: alpha_max = 1d0 
   complex(8),allocatable,dimension(:),save :: Fz_prev !prev iter
   real(8),allocatable,dimension(:),save    :: beta    !mixing β
   integer                                  :: j
   logical                                  :: real_oscillates, imag_oscillates
   !
   N = size(z)
   !
   if(iter==1)then !reset persistent variables
      if(allocated(Fz_prev))deallocate(Fz_prev)
      allocate(Fz_prev(N))
      Fz_prev = zero
      !
      if(allocated(beta))deallocate(beta)
      allocate(beta(N))
      beta = alpha
   endif
   !
   Z = Z + beta*Fz ! = (1-β)z_in + βz_out (β: effective speed)            
   !
   if (iter > 1) then !update effective elemental mixing speed
      do j=1,N
         real_oscillates = dreal(Fz_prev(j))*dreal(Fz(j)) < 0
         imag_oscillates = dimag(Fz_prev(j))*dimag(Fz(j)) < 0
         if(real_oscillates.OR.imag_oscillates)then 
            beta(j) = alpha !reset to input mixing
         else !increase mixing speed by adding alpha
            beta(j) = beta(j) + alpha !but with a cutoff
            if (beta(j) > alpha_max) beta(j) = alpha_max
         end if
      end do
   end if
   !
   Fz_prev = Fz !store Fz_prev for future calls
   !
end subroutine c_adaptive_mix
