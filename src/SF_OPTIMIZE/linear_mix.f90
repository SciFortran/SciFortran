subroutine d_linear_mix_1(x,Fx,alpha)
  real(8),intent(inout),dimension(:)      :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x))   :: Fx     ! x_out-x_in
  real(8),intent(in)                      :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_1

subroutine d_linear_mix_2(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:)              :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2)) :: Fx     ! x_out-x_in
  real(8),intent(in)                                :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_2

subroutine d_linear_mix_3(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:,:)                      :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3)) :: Fx     ! x_out-x_in
  real(8),intent(in)                                          :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_3

subroutine d_linear_mix_4(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:,:,:)                              :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4)) :: Fx     ! x_out-x_in
  real(8),intent(in)                                                    :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_4

subroutine d_linear_mix_5(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:,:,:,:)                                      :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5)) :: Fx     ! x_out-x_in
  real(8),intent(in)                                                              :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_5

subroutine d_linear_mix_6(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:,:,:,:,:)                                               :: x      ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6)) :: Fx      ! x_out-x_in
  real(8),intent(in)                                                                         :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_6

subroutine d_linear_mix_7(x,Fx,alpha)
  real(8),intent(inout),dimension(:,:,:,:,:,:,:) :: x  ! in: x_in; out: mixed x
  real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7)) :: Fx     ! x_out-x_in
  real(8),intent(in)                                                                                  :: alpha  ! mixing parameter
  x = x + alpha*Fx
end subroutine d_linear_mix_7


subroutine c_linear_mix_1(x,Fx,alpha)
  complex(8),intent(inout),dimension(:)    :: x
  complex(8),intent(in),dimension(size(x)) :: Fx
  real(8),intent(in)                       :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_1

subroutine c_linear_mix_2(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:)              :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2)) :: Fx
  real(8),intent(in)                                   :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_2

subroutine c_linear_mix_3(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:,:)                      :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3)) :: Fx
  real(8),intent(in)                                             :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_3

subroutine c_linear_mix_4(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:,:,:)                              :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4)) :: Fx
  real(8),intent(in)                                                       :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_4

subroutine c_linear_mix_5(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:,:,:,:)                                      :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5)) :: Fx
  real(8),intent(in)                                                                 :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_5

subroutine c_linear_mix_6(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:,:,:,:,:)                                               :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6)) :: Fx
  real(8),intent(in)                                                                            :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_6

subroutine c_linear_mix_7(x,Fx,alpha)
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:)                                                      :: x
  complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7)) :: Fx
  real(8),intent(in)                                                                                     :: alpha
  x = x + alpha*Fx
end subroutine c_linear_mix_7
