MODULE SF_SPIN
!SciFortan module for Pauli matrix algebra
  implicit none
  private

  complex(8),parameter :: zero =(0d0,0d0)
  complex(8),parameter :: xi   =(0d0,1d0)
  complex(8),parameter :: one  =(1d0,0d0)
  real(8),parameter    :: sqrt2=sqrt(2d0)
  real(8),parameter    :: sqrt3=sqrt(3d0)


  complex(8),dimension(2,2),parameter,public :: pauli_0=&  !Pauli matrix 0 (2x2 identity) :math:`\sigma_0 = \begin{pmatrix}1&0\\0&1\end{pmatrix}`
                                                reshape([one,zero,zero,one],[2,2])   
  complex(8),dimension(2,2),parameter,public :: pauli_x=&  !Pauli matrix x :math:`\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}`
                                                reshape([zero,one,one,zero],[2,2])   
  complex(8),dimension(2,2),parameter,public :: pauli_y=&  !Pauli matrix y :math:`\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix}`
                                                reshape([zero,xi,-xi,zero],[2,2])    
  complex(8),dimension(2,2),parameter,public :: pauli_z=&  !Pauli matrix z :math:`\sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}`
                                                reshape([one,zero,zero,-one],[2,2])  
  
  complex(8),dimension(2,2),parameter,public :: pauli_1=&  !:math:`\sigma_1 = \sigma_x`
                                                pauli_x                            
  complex(8),dimension(2,2),parameter,public :: pauli_2=&  !:math:`\sigma_2 = \sigma_y`
                                                pauli_y                             
  complex(8),dimension(2,2),parameter,public :: pauli_3=&  !:math:`\sigma_3 = \sigma_z`
                                                pauli_z                              
  
  complex(8),dimension(2,2),parameter,public :: pauli_tau_0=& !:math:`\tau_0 = \sigma_0`
                                                pauli_0                          
  complex(8),dimension(2,2),parameter,public :: pauli_tau_x=& !:math:`\tau_x = \sigma_x`
                                                pauli_x                          
  complex(8),dimension(2,2),parameter,public :: pauli_tau_y=& !:math:`\tau_y = \sigma_y`
                                                pauli_y                          
  complex(8),dimension(2,2),parameter,public :: pauli_tau_z=& !:math:`\tau_z = \sigma_z`
                                                pauli_z                          
  
  complex(8),dimension(2,2),parameter,public :: pauli_tau_1=& !:math:`\tau_1 = \sigma_x`
                                                pauli_x                          
  complex(8),dimension(2,2),parameter,public :: pauli_tau_2=& !:math:`\tau_2 = \sigma_y`
                                                pauli_y                          
  complex(8),dimension(2,2),parameter,public :: pauli_tau_3=& !:math:`\tau_3 = \sigma_z`
                                                pauli_z                          
  
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_0=& !:math:`\sigma_0 = \begin{pmatrix}1&0\\0&1\end{pmatrix}`
                                                pauli_0                        
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_x=& !:math:`\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}`
                                                pauli_x                        
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_y=& !:math:`\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix}`
                                                pauli_y                        
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_z=& !:math:`\sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}`
                                                pauli_z                        
  
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_1=& !:math:`\sigma_1 = \sigma_x`
                                                pauli_x                        
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_2=& !:math:`\sigma_2 = \sigma_y`
                                                pauli_y                        
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_3=& !:math:`\sigma_3 = \sigma_z`
                                                pauli_z                        
  
  
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_plus =&    !2 x raising operator :math:`\sigma_+ = \sigma_x + i \sigma_y = \begin{pmatrix}0&2\\0&0\end{pmatrix}`
                                                pauli_x+xi*pauli_y         
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_minus=&    !2 x lowering operator :math:`\sigma_- = \sigma_x - i \sigma_y = \begin{pmatrix}0&0\\2&0\end{pmatrix}`
                                                pauli_x-xi*pauli_y         
  
  complex(8),dimension(2,2),parameter,public :: pauli_tau_plus =&      !2 x raising operator :math:`\tau_+ = \tau_x + i \tau_y = \begin{pmatrix}0&2\\0&0\end{pmatrix}`
                                                pauli_x+xi*pauli_y         
  complex(8),dimension(2,2),parameter,public :: pauli_tau_minus=&      !2 x lowering operator :math:`\tau_- = \tau_x - i \tau_y = \begin{pmatrix}0&0\\2&0\end{pmatrix}`
                                                pauli_x-xi*pauli_y         




  !SPIN 1
  complex(8),dimension(3,3),parameter,public :: spin1_0=reshape([&        !3x3 identity matrix :math:`\mathbb{1}_3 = \begin{pmatrix}1&0&0\\0&1&0\\0&0&1\end{pmatrix}`
       one ,zero,zero, &
       zero, one,zero, &
       zero,zero,one   &
       ], [3,3])
  complex(8),dimension(3,3),parameter,public :: spin1_x=reshape([&        !spin 1 matrix :math:`S_x`, basis m=1,0,-1 :math:`S_x = \frac{1}{\sqrt{2}}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix}`
       zero, one,zero,  &
       one ,zero, one,  &
       zero, one,zero   &
       ],[3,3])/sqrt2 
  complex(8),dimension(3,3),parameter,public :: spin1_y=reshape([&        !spin 1 matrix :math:`S_y`, basis m=1,0,-1 :math:`S_y = \frac{1}{\sqrt{2}}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix}`
       zero, -xi ,zero,  &
       xi  ,zero , -xi,  &
       zero,xi   ,zero   &
       ], [3,3])/sqrt2 
  complex(8),dimension(3,3),parameter,public :: spin1_z=reshape([&        !spin 1 matrix :math:`S_z`, basis m=1,0,-1 :math:`S_z = \begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix}`
       one ,zero,zero, &
       zero,zero,zero, &
       zero,zero,-one  &
       ],[3,3]) 
  complex(8),dimension(3,3),parameter,public :: spin1_plus  = &           !spin 1 raising operator :math:`S_+ = S_x + i S_y = \sqrt{2}\begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix}`
                                                spin1_x+xi*spin1_y             
  
  complex(8),dimension(3,3),parameter,public :: spin1_minus = &           !spin 1 lowering operator :math:`S_- = S_x - i S_y = \sqrt{2}\begin{pmatrix}0&0&0\\1&0&0\\0&1&0\end{pmatrix}`
                                                spin1_x-xi*spin1_y             


  !SPIN 3/2
  complex(8),parameter :: two=(2d0,0d0)
  complex(8),parameter :: c2=(0d0,2d0)
  complex(8),parameter :: s3=(sqrt3,0d0)
  complex(8),parameter :: c3=(0d0,sqrt3)
  complex(8),parameter :: h12=(0.5d0,0d0)
  complex(8),parameter :: h32=(1.5d0,0d0)

  complex(8),dimension(4,4),parameter,public :: spin3half_0=reshape([&                      !4x4 identity matrix :math:`\mathbb{1}_4 = \begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&1&0\\0&0&0&1\end{pmatrix}`
       one,zero,zero,zero,  &
       zero,one,zero,zero,  &
       zero,zero,one,zero,  &
       zero,zero,zero,one   &
       ],[4,4])            
  complex(8),dimension(4,4),parameter,public :: spin3half_x=reshape([&                      !spin 3/2 matrix :math:`S_x`, basis m=3/2,...,-3/2 :math:`S_x = \frac{1}{2}\begin{pmatrix}0&\sqrt{3}&0&0\\\sqrt{3}&0&2&0\\0&2&0&\sqrt{3}\\0&0&\sqrt{3}&0\end{pmatrix}`
       zero ,  s3 , zero , zero ,  &
       s3   , zero, two  , zero ,  &
       zero , two , zero , s3   ,  &
       zero , zero,  s3  , zero    &
       ],[4,4])/2d0        
  complex(8),dimension(4,4),parameter,public :: spin3half_y=reshape([&                      !spin 3/2 matrix :math:`S_y`, basis m=3/2,...,-3/2 :math:`S_y = \frac{1}{2}\begin{pmatrix}0&-i\sqrt{3}&0&0\\i\sqrt{3}&0&-2i&0\\0&2i&0&-i\sqrt{3}\\0&0&i\sqrt{3}&0\end{pmatrix}`
       zero , -c3 , zero , zero ,  &
       c3   , zero, -c2  , zero ,  &
       zero ,  c2 , zero , -c3  ,  &
       zero , zero, c3   , zero    &
       ],[4,4])/2d0        
  complex(8),dimension(4,4),parameter,public :: spin3half_z=reshape([&                      !spin 3/2 matrix :math:`S_z`, basis m=3/2,...,-3/2 :math:`S_z = \frac{1}{2}\begin{pmatrix}3&0&0&0\\0&1&0&0\\0&0&-1&0\\0&0&0&-3\end{pmatrix}`
       h32,zero,zero,zero,  &
       zero,h12,zero,zero,  &
       zero,zero,-h12,zero,  &
       zero,zero,zero,-h32   &
       ],[4,4])            
  complex(8),dimension(4,4),parameter,public :: spin3Half_plus  = &        !spin 3/2 raising operator :math:`S_+ = S_x + i S_y = \begin{pmatrix}0&\sqrt{3}&0&0\\0&0&2&0\\0&0&0&\sqrt{3}\\0&0&0&0\end{pmatrix}`
                                                spin3Half_x+xi*spin3Half_y 
  complex(8),dimension(4,4),parameter,public :: spin3Half_minus = &        !spin 3/2 lowering operator :math:`S_- = S_x - i S_y = \begin{pmatrix}0&0&0&0\\\sqrt{3}&0&0&0\\0&2&0&0\\0&0&\sqrt{3}&0\end{pmatrix}`
                                                spin3Half_x-xi*spin3Half_y 

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
