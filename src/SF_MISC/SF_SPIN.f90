MODULE SF_SPIN
!SciFortan module for Pauli matrix algebra
  implicit none
  private

  complex(8),parameter :: zero =(0d0,0d0)
  complex(8),parameter :: xi   =(0d0,1d0)
  complex(8),parameter :: one  =(1d0,0d0)
  real(8),parameter    :: sqrt2=sqrt(2d0)
  real(8),parameter    :: sqrt3=sqrt(3d0)


  complex(8),dimension(2,2),parameter,public :: pauli_0=reshape([one,zero,zero,one],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_x=reshape([zero,one,one,zero],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_y=reshape([zero,xi,-xi,zero],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_z=reshape([one,zero,zero,-one],[2,2])
  !
  complex(8),dimension(2,2),parameter,public :: pauli_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_3=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_0=pauli_0
  complex(8),dimension(2,2),parameter,public :: pauli_tau_x=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_tau_y=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_z=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_tau_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_3=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_0=pauli_0
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_x=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_y=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_z=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_3=pauli_z
  !
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_plus =pauli_x+xi*pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_minus=pauli_x-xi*pauli_y
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_plus =pauli_x+xi*pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_minus=pauli_x-xi*pauli_y




  !SPIN 1
  complex(8),dimension(3,3),parameter,public :: spin1_0=reshape([&
       one ,zero,zero, &
       zero, one,zero, &
       zero,zero,one   &
       ], [3,3])
  complex(8),dimension(3,3),parameter,public :: spin1_x=reshape([&
       zero, one,zero,  &
       one ,zero, one,  &
       zero, one,zero   &
       ],[3,3])/sqrt2
  complex(8),dimension(3,3),parameter,public :: spin1_y=reshape([&
       zero, -xi ,zero,  &
       xi  ,zero , -xi,  &
       zero,xi   ,zero   &
       ], [3,3])/sqrt2
  complex(8),dimension(3,3),parameter,public :: spin1_z=reshape([&
       one ,zero,zero, &
       zero,zero,zero, &
       zero,zero,-one  &
       ],[3,3])
  complex(8),dimension(3,3),parameter,public :: spin1_plus  = spin1_x+xi*spin1_y
  complex(8),dimension(3,3),parameter,public :: spin1_minus = spin1_x-xi*spin1_y


  !SPIN 3/2
  complex(8),parameter :: two=(2d0,0d0)
  complex(8),parameter :: c2=(0d0,2d0)
  complex(8),parameter :: s3=(sqrt3,0d0)
  complex(8),parameter :: c3=(0d0,sqrt3)
  complex(8),parameter :: h12=(0.5d0,0d0)
  complex(8),parameter :: h32=(1.5d0,0d0)

  complex(8),dimension(4,4),parameter,public :: spin3half_0=reshape([&
       one,zero,zero,zero,  &
       zero,one,zero,zero,  &
       zero,zero,one,zero,  &
       zero,zero,zero,one   &
       ],[4,4])
  complex(8),dimension(4,4),parameter,public :: spin3half_x=reshape([&
       zero ,  s3 , zero , zero ,  &
       s3   , zero, two  , zero ,  &
       zero , two , zero , s3   ,  &
       zero , zero,  s3  , zero    &
       ],[4,4])/2d0
  complex(8),dimension(4,4),parameter,public :: spin3half_y=reshape([&
       zero , -c3 , zero , zero ,  &
       c3   , zero, -c2  , zero ,  &
       zero ,  c2 , zero , -c3  ,  &
       zero , zero, c3   , zero    &
       ],[4,4])/2d0
  complex(8),dimension(4,4),parameter,public :: spin3half_z=reshape([&
       h32,zero,zero,zero,  &
       zero,h12,zero,zero,  &
       zero,zero,-h12,zero,  &
       zero,zero,zero,-h32   &
       ],[4,4])
  complex(8),dimension(4,4),parameter,public :: spin3Half_plus  = spin3Half_x+xi*spin3Half_y
  complex(8),dimension(4,4),parameter,public :: spin3Half_minus = spin3Half_x-xi*spin3Half_y

END MODULE SF_SPIN





! !---------------------------------------------------------------------
! !PURPOSE: return the Kronecker's product of n Pauli's matrices. 
! ! The especification of the order of the matrices is given as input on the 
! ! vector vec_ord_pm, that has dimension npm.
! !---------------------------------------------------------------------
! function kronecker_product_pauli_recursive(vec_ord_pm) result(kron_prod_n_pauli_mat)
!   integer, intent(in)     :: vec_ord_pm(:)
!   complex(8)              :: kron_prod_n_pauli_mat(2**size(vec_ord_pm),2**size(vec_ord_pm))
!   integer                 :: d2
!   complex(8)              :: M2(2,2)
!   complex(8), allocatable :: M1(:,:), M1_kp_M2(:,:)
!   integer                 :: npm,d1, i
!   npm=size(vec_ord_pm)
!   d2=2
!   do i=1,npm-1
!      select case(vec_ord_pm(i+1))
!      case (0)
!         M2 = pauli_sigma_0
!      case (1)
!         M2 = pauli_sigma_1
!      case (2)
!         M2 = pauli_sigma_2
!      case (3)
!         M2 = pauli_sigma_3
!      end select
!      d1 = 2**i
!      if(i==1) then
!         allocate(M1(d1,d1))
!         select case(vec_ord_pm(i))
!         case (0) 
!            M1 = pauli_sigma_0
!         case (1) 
!            M1 = pauli_sigma_1
!         case (2) 
!            M1 = pauli_sigma_2
!         case (3) 
!            M1 = pauli_sigma_3
!         end select
!      endif
!      allocate(M1_kp_M2(d1*d2,d1*d2))
!      M1_kp_M2 = c_kronecker_product(M1,d1,d1,M2,d2,d2)  
!      deallocate(M1)
!      allocate(M1(1:d1*d2,1:d1*d2))
!      M1 = M1_kp_M2
!      deallocate(M1_kp_M2)
!   end do
!   kron_prod_n_pauli_mat = M1
!   deallocate(M1)
! end function kronecker_product_pauli_recursive

