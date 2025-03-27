MODULE SF_SPARSE_ARRAY_CSC
  USE SF_LINALG, only: deye,zeye
  USE SF_SPARSE_COMMON
  !MPI STILL TO BE DONE
#ifdef _MPI
  USE MPI
#endif
  implicit none

  private
  !CSC COLS
  type sparse_col_csc
     integer                                   :: size
     integer,dimension(:),allocatable          :: rows     
  end type sparse_col_csc
  !
  type, extends(sparse_col_csc) :: sparse_dcol_csc
     real(8),dimension(:),allocatable          :: vals
  end type sparse_dcol_csc
  !
  type, extends(sparse_col_csc) :: sparse_zcol_csc
     complex(8),dimension(:),allocatable       :: vals
  end type sparse_zcol_csc

  type, extends(sparse_matrix) :: sparse_dmatrix_csc
     type(sparse_dcol_csc),dimension(:),allocatable :: col
   contains
     procedure,pass :: init      => init_dmatrix_csc
     procedure,pass :: free      => free_dmatrix_csc
     procedure,pass :: load      => load_dmatrix_csc
     procedure,pass :: dump      => dump_dmatrix_csc
     procedure,pass :: insert    => insert_delement_csc
     procedure,pass :: get       => get_delement_csc
     procedure,pass :: as_matrix => as_dmatrix_csc
     procedure,pass :: dgr       => dgr_dmatrix_csc
     procedure,pass :: transpose => transpose_dmatrix_csc
  end type sparse_dmatrix_csc
  
  type, extends(sparse_matrix) :: sparse_zmatrix_csc
     type(sparse_zcol_csc),dimension(:),allocatable :: col
   contains
     procedure,pass :: init      => init_zmatrix_csc
     procedure,pass :: free      => free_zmatrix_csc
     procedure,pass :: load      => load_zmatrix_csc
     procedure,pass :: dump      => dump_zmatrix_csc
     procedure,pass :: insert    => insert_zelement_csc
     procedure,pass :: get       => get_zelement_csc
     procedure,pass :: as_matrix => as_zmatrix_csc
     procedure,pass :: dgr       => dgr_zmatrix_csc
     procedure,pass :: transpose => transpose_zmatrix_csc
     procedure,pass :: conjg     => conjg_zmatrix_csc
  end type sparse_zmatrix_csc

  
  
  interface sparse_matrix
     module procedure :: construct_dmatrix_csc
     module procedure :: construct_zmatrix_csc
  end interface sparse_matrix
  !
  interface as_sparse
     module procedure :: construct_dmatrix_csc
     module procedure :: construct_zmatrix_csc
  end interface as_sparse
  !
  interface sparse
     module procedure :: construct_dmatrix_csc
     module procedure :: construct_zmatrix_csc
  end interface sparse
  
  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
     module procedure :: dmatrix_equal_scalar_csc
     module procedure :: dmatrix_equal_dmatrix_csc
     module procedure :: zmatrix_equal_scalar_csc
     module procedure :: zmatrix_equal_zmatrix_csc
  end interface assignment(=)

  !ADDITION
  interface operator(+)
     module procedure :: plus_dmatrix_csc
     module procedure :: plus_zmatrix_csc
  end interface operator(+)

  !SUBTRACTION
  interface operator(-)
     module procedure :: minus_dmatrix_csc
     module procedure :: minus_zmatrix_csc
  end interface operator(-)
  
  !SCALAR PRODUCT
  interface operator(*)
     module procedure :: left_product_dmatrix_i_csc
     module procedure :: left_product_dmatrix_d_csc
     module procedure :: left_product_zmatrix_i_csc
     module procedure :: left_product_zmatrix_d_csc
     module procedure :: left_product_zmatrix_z_csc
     !
     module procedure :: right_product_dmatrix_i_csc
     module procedure :: right_product_dmatrix_d_csc
     module procedure :: right_product_zmatrix_i_csc
     module procedure :: right_product_zmatrix_d_csc
     module procedure :: right_product_zmatrix_z_csc
  end interface operator(*)

  !SCALAR DIVISION
  interface operator(/)
     module procedure :: right_division_dmatrix_i_csc
     module procedure :: right_division_dmatrix_d_csc
     module procedure :: right_division_zmatrix_i_csc
     module procedure :: right_division_zmatrix_d_csc
     module procedure :: right_division_zmatrix_z_csc
  end interface operator(/)

  !MATRIX PRODUCT
  interface matmul
     module procedure :: left_matmul_dmatrix_darray_csc
     module procedure :: left_matmul_zmatrix_zarray_csc 
  end interface matmul
  

  !KRONECKER PRODUCT
  interface operator(.x.)
     module procedure :: kron_dmatrix_csc
     module procedure :: kron_zmatrix_csc
  end interface operator(.x.)

  interface sp_kron
     module procedure :: restricted_kron_dmatrix_csc
     module procedure :: restricted_kron_zmatrix_csc
  end interface sp_kron

  interface transpose
     module procedure :: transpose_dmatrix_csc
     module procedure :: transpose_zmatrix_csc
  end interface transpose

  interface hconjg
     module procedure :: dgr_dmatrix_csc
     module procedure :: dgr_zmatrix_csc
  end interface hconjg
  
  public :: sparse_dmatrix_csc, sparse_zmatrix_csc  
  public :: as_sparse
  public :: sparse
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.x.)
  public :: sp_kron
  public :: transpose
  public :: hconjg
  public :: matmul


  
contains
  
  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !          init:           empty matrix of dimenson (N,N1)
  !          constructor:    dump a dense matrix into the sparse
  !+------------------------------------------------------------------+
  !REAL
  elemental subroutine init_dmatrix_csc(sparse,N,N1)
    class(sparse_dmatrix_csc),intent(inout) :: sparse
    integer,intent(in)                    :: N
    integer,intent(in),optional           :: N1
    integer                               :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)call sparse%free()    
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%col(sparse%Ncol))
    do i=1,sparse%Ncol
       sparse%col(i)%size=0
       allocate(sparse%col(i)%vals(0)) !empty array
       allocate(sparse%col(i)%rows(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine init_dmatrix_csc
  !
  function construct_dmatrix_csc(matrix) result(self)
    real(8),dimension(:,:),intent(in)    :: matrix
    type(sparse_dmatrix_csc)        :: self
    call self%load(matrix)
  end function construct_dmatrix_csc
  !
  subroutine load_dmatrix_csc(sparse,matrix)
    class(sparse_dmatrix_csc),intent(inout) :: sparse
    real(8),dimension(:,:),intent(in)  :: matrix
    integer                            :: i,j,Ndim1,Ndim2
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2)
    do j=1,Ndim2
       do i=1,Ndim1
          if(matrix(i,j)/=0d0)call insert_delement_csc(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine load_dmatrix_csc
  !
  !COMPLEX
  elemental subroutine init_zmatrix_csc(sparse,N,N1)
    class(sparse_zmatrix_csc),intent(inout) :: sparse
    integer,intent(in)                      :: N
    integer,intent(in),optional             :: N1
    integer                                 :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)call sparse%free()
    !
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%col(sparse%Ncol))
    do i=1,sparse%Ncol
       sparse%col(i)%size=0
       allocate(sparse%col(i)%vals(0)) !empty array
       allocate(sparse%col(i)%rows(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine init_zmatrix_csc
  !
  function construct_zmatrix_csc(matrix) result(self)
    complex(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_zmatrix_csc)              :: self
    call self%load(matrix)
  end function construct_zmatrix_csc
  !
  subroutine load_zmatrix_csc(sparse,matrix)
    class(sparse_zmatrix_csc),intent(inout) :: sparse
    complex(8),dimension(:,:),intent(in)  :: matrix
    integer                            :: i,j,Ndim1,Ndim2
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2)
    do j=1,Ndim2
       do i=1,Ndim1
          if(matrix(i,j)/=0d0)call insert_zelement_csc(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine load_zmatrix_csc

  !+------------------------------------------------------------------+
  !PURPOSE: frees a sparse matrix
  !+------------------------------------------------------------------+
  !REAL
  elemental subroutine free_dmatrix_csc(sparse)    
    class(sparse_dmatrix_csc),intent(inout) :: sparse
    integer                            :: i
    !
    do i=1,sparse%Ncol
       sparse%col(i)%Size  = 0
       if(allocated(sparse%col(i)%vals))deallocate(sparse%col(i)%vals)
       if(allocated(sparse%col(i)%rows))deallocate(sparse%col(i)%rows)
    enddo
    if(allocated(sparse%col))deallocate(sparse%col)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine free_dmatrix_csc
  !
  !COMPLEX
  elemental subroutine free_zmatrix_csc(sparse)    
    class(sparse_zmatrix_csc),intent(inout) :: sparse
    integer                                 :: i
    !
    do i=1,sparse%Ncol
       sparse%col(i)%Size  = 0
       if(allocated(sparse%col(i)%vals))deallocate(sparse%col(i)%vals)
       if(allocated(sparse%col(i)%rows))deallocate(sparse%col(i)%rows)
    enddo
    if(allocated(sparse%col))deallocate(sparse%col)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine free_zmatrix_csc



  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  !REAL
  subroutine insert_delement_csc(sparse,value,i,j)
    class(sparse_dmatrix_csc),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    integer                               :: row,pos
    logical                               :: iadd
    !
    row = i
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%col(j)%rows == row))then         !
       pos = binary_search(sparse%col(j)%rows,row) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%col(j)%vals(pos)=sparse%col(j)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it
       call append(sparse%col(j)%vals,value)
       call append(sparse%col(j)%rows,row)
       sparse%col(j)%Size = sparse%col(j)%Size + 1
       sparse%Nelements = sparse%Nelements+1
    endif
    !
    if(sparse%col(j)%Size > sparse%Nrow)stop "insert_delement_csc ERROR: row%Size > sparse%Ncol"
  end subroutine insert_delement_csc
  !
  !COMPLEX
  subroutine insert_zelement_csc(sparse,value,i,j)
    class(sparse_zmatrix_csc),intent(inout) :: sparse
    complex(8),intent(in)                 :: value
    integer,intent(in)                    :: i,j
    integer                               :: row,pos
    logical                               :: iadd
    !
    !
    row = i
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%col(j)%rows == row))then         !
       pos = binary_search(sparse%col(j)%rows,row) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%col(j)%vals(pos)=sparse%col(j)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call append(sparse%col(j)%vals,value)
       call append(sparse%col(j)%rows,row)
       sparse%col(j)%Size = sparse%col(j)%Size + 1
       sparse%Nelements = sparse%Nelements+1
    endif
    !
    if(sparse%col(j)%Size > sparse%Nrow)stop "insert_zelement_csc ERROR: row%Size > sparse%Ncol"
    !
  end subroutine insert_zelement_csc

    !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  function get_delement_csc(sparse,i,j) result(value)
    class(sparse_dmatrix_csc),intent(inout) :: sparse    
    integer,intent(in)                 :: i,j
    real(8)                            :: value
    integer                            :: pos
    value=0d0
    do pos=1,sparse%col(j)%size
       if(i==sparse%col(j)%rows(pos))value=sparse%col(j)%vals(pos)
    enddo
  end function get_delement_csc

  function get_zelement_csc(sparse,i,j) result(value)
    class(sparse_zmatrix_csc),intent(inout) :: sparse    
    integer,intent(in)                 :: i,j
    complex(8)                         :: value
    integer                            :: pos
    value=0d0
    do pos=1,sparse%col(j)%size
       if(i==sparse%col(j)%rows(pos))value=sparse%col(j)%vals(pos)
    enddo
  end function get_zelement_csc



  

  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array summing elements
  !+------------------------------------------------------------------+
  !REAL
  subroutine dump_dmatrix_csc(sparse,matrix)
    class(sparse_dmatrix_csc),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    !
    !
    Ndim1=size(matrix,1);  Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_dmatrix: dimensions error"
    !
    do j=1,Ndim2
       do i=1,sparse%col(j)%Size
          matrix(sparse%col(j)%rows(i),j) = matrix(sparse%col(j)%rows(i),j) + sparse%col(j)%vals(i)
       enddo
    enddo
  end subroutine dump_dmatrix_csc
  !
  !COMPLEX
  subroutine dump_zmatrix_csc(sparse,matrix)
    class(sparse_zmatrix_csc),intent(in)      :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    integer                                 :: i,j,Ndim1,Ndim2
    !
    !
    Ndim1=size(matrix,1);   Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_zmatrix: dimensions error"
    !
    matrix=0.d0
    do j=1,Ndim2
       do i=1,sparse%col(j)%Size
          matrix(sparse%col(j)%rows(i),j) = matrix(sparse%col(j)%rows(i),j) + sparse%col(j)%vals(i)
       enddo
    enddo
  end subroutine dump_zmatrix_csc


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array replacing elements
  !+------------------------------------------------------------------+
  !REAL
  function as_dmatrix_csc(sparse) result(matrix)
    class(sparse_dmatrix_csc),intent(in)            :: sparse
    real(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
    integer                                    :: i,j
    matrix = 0d0
    do j=1,sparse%Ncol
       do i=1,sparse%col(j)%Size
          matrix(sparse%col(j)%rows(i),j) = sparse%col(j)%vals(i)
       enddo
    enddo
  end function as_dmatrix_csc
  
  !COMPLEX
  function as_zmatrix_csc(sparse) result(matrix)
    class(sparse_zmatrix_csc),intent(in)            :: sparse
    complex(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
    integer                                    :: i,j
    matrix = 0d0
    do j=1,sparse%Ncol
       do i=1,sparse%col(j)%Size
          matrix(sparse%col(j)%rows(i),j) = sparse%col(j)%vals(i)
       enddo
    enddo
  end function as_zmatrix_csc


  !##################################################################
  !               SPARSE MATRIX KRON PRODUCT 
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
  !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  !+------------------------------------------------------------------+
  !REAL
  function kron_dmatrix_csc(A,B) result(AxB)
    type(sparse_dmatrix_csc), intent(in) :: A,B
    type(sparse_dmatrix_csc)             :: AxB
    integer                         :: i,j,k,l,jrow,krow
    integer                         :: indx_row,indx_col
    real(8)                         :: value
    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    do indx_col =1,A%Ncol*B%Ncol
       k=mod(indx_col,B%Ncol); if(k==0)k=B%Ncol
       j= (indx_col-1)/B%Ncol+1
       do jrow=1,A%col(j)%size
          i=A%col(j)%rows(jrow)
          do krow=1,B%col(k)%size
             l = B%col(k)%rows(krow)
             indx_row = l + (i-1)*B%Nrow
             value = A%col(j)%vals(jrow)*B%col(k)%vals(krow)
             !
             call append(AxB%col(indx_col)%vals,value)
             call append(AxB%col(indx_col)%rows,indx_row)
             AxB%col(indx_col)%Size = AxB%col(indx_col)%Size + 1
             AxB%Nelements = AxB%Nelements+1
          end do
       end do
    end do
  end function kron_dmatrix_csc
  !
  function restricted_kron_dmatrix_csc(A,B,states) result(AxB)
    type(sparse_dmatrix_csc), intent(in)  :: A,B
    integer,dimension(:),intent(in)  :: states
    type(sparse_dmatrix_csc)      :: AxB
    integer                          :: i,irow,j,k,krow,l,istate,jstate
    integer                          :: indx_row,indx_col
    real(8)                          :: value
    integer,dimension(:),allocatable :: inv_states
    call AxB%free()
    call AxB%init(size(states),size(states))
    allocate(inv_states(A%Ncol*B%Ncol))
    do i=1,size(states)
       inv_states(states(i)) = i
    enddo
    do istate = 1,size(states)
       indx_col=states(istate)
       k = mod(indx_col,B%Ncol);if(k==0)k=B%Ncol
       i = (indx_col-1)/B%Ncol+1
       do irow=1,A%col(i)%size
          j = A%col(i)%rows(irow)
          do krow=1,B%col(k)%size
             l = B%col(k)%rows(krow)
             indx_row = l + (j-1)*B%Nrow
             jstate   = inv_states(indx_row)
             value    = A%col(i)%vals(irow)*B%col(k)%vals(krow)
             !
             call append(AxB%col(istate)%vals,value)
             call append(AxB%col(istate)%rows,jstate)
             AxB%col(istate)%Size = AxB%col(istate)%Size + 1
             AxB%Nelements = AxB%Nelements + 1
             !
          enddo
       enddo
    enddo
  end function restricted_kron_dmatrix_csc
  !
  !COMPLEX
  function kron_zmatrix_csc(A,B) result(AxB)
    type(sparse_zmatrix_csc), intent(in) :: A,B
    type(sparse_zmatrix_csc)             :: AxB
    integer                         :: i,jrow,j,k,krow,l
    integer                         :: indx_row,indx_col
    complex(8)                      :: value
    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    do indx_col =1,A%Ncol*B%Ncol
       k=mod(indx_col,B%Ncol); if(k==0)k=B%Ncol
       j= (indx_col-1)/B%Ncol+1
       do jrow=1,A%col(j)%size
          i=A%col(j)%rows(jrow)
          do krow=1,B%col(k)%size
             l = B%col(k)%rows(krow)
             indx_row = l + (i-1)*B%Nrow
             value = A%col(j)%vals(jrow)*B%col(k)%vals(krow)
             !
             call append(AxB%col(indx_col)%vals,value)
             call append(AxB%col(indx_col)%rows,indx_row)
             AxB%col(indx_col)%Size = AxB%col(indx_col)%Size + 1
             AxB%Nelements = AxB%Nelements+1
          end do
       end do
    end do
  end function kron_zmatrix_csc
  !
  function restricted_kron_zmatrix_csc(A,B,states) result(AxB)
    type(sparse_zmatrix_csc), intent(in)  :: A,B
    integer,dimension(:),intent(in)  :: states
    type(sparse_zmatrix_csc)      :: AxB
    integer                          :: i,irow,j,k,krow,l,istate,jstate
    integer                          :: indx_row,indx_col
    complex(8)                          :: value
    integer,dimension(:),allocatable :: inv_states
    call AxB%free()
    call AxB%init(size(states),size(states))
    allocate(inv_states(A%Ncol*B%Ncol))
    do i=1,size(states)
       inv_states(states(i)) = i
    enddo
    do istate = 1,size(states)
       indx_col=states(istate)
       k = mod(indx_col,B%Ncol);if(k==0)k=B%Ncol
       i = (indx_col-1)/B%Ncol+1
       do irow=1,A%col(i)%size
          j = A%col(i)%rows(irow)
          do krow=1,B%col(k)%size
             l = B%col(k)%rows(krow)
             indx_row = l + (j-1)*B%Nrow
             jstate   = inv_states(indx_row)
             value    = A%col(i)%vals(irow)*B%col(k)%vals(krow)
             !
             call append(AxB%col(istate)%vals,value)
             call append(AxB%col(istate)%rows,jstate)
             AxB%col(istate)%Size = AxB%col(istate)%Size + 1
             AxB%Nelements = AxB%Nelements + 1
             !
          enddo
       enddo
    enddo
  end function restricted_kron_zmatrix_csc
  !
  

  
  !##################################################################
  !               SPARSE MATRIX BASIC ALGEBRA 
  !##################################################################
  !REAL
  function dgr_dmatrix_csc(a) result(c)
    class(sparse_dmatrix_csc), intent(in) :: a
    type(sparse_dmatrix_csc)              :: c
    integer                          :: row
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)
          call c%insert(val,j,row)
       enddo
    enddo
  end function dgr_dmatrix_csc
  !
  function transpose_dmatrix_csc(a) result(c)
    class(sparse_dmatrix_csc), intent(in) :: a
    type(sparse_dmatrix_csc)              :: c
    integer                          :: row
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)
          call c%insert(val,j,row)
       enddo
    enddo
  end function transpose_dmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine dmatrix_equal_dmatrix_csc(a,b)
    type(sparse_dmatrix_csc),intent(inout) :: a
    type(sparse_dmatrix_csc),intent(in)    :: b
    integer                           :: row
    real(8)                           :: val
    integer                           :: i,j    
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = b%col(j)%vals(i)
          call a%insert(val,row,j)
       enddo
    enddo
  end subroutine dmatrix_equal_dmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine dmatrix_equal_scalar_csc(a,c)
    type(sparse_dmatrix_csc),intent(inout) :: a
    real(8),intent(in)                :: c
    integer                           :: i,j    
    ! if(.not.a%status)stop "dmatrix_equal_scalar error: a is not allocated"
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          a%col(j)%vals(i) = c
       enddo
    enddo
  end subroutine dmatrix_equal_scalar_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function plus_dmatrix_csc(a,b) result(c)
    type(sparse_dmatrix_csc), intent(in) :: a !First addend
    type(sparse_dmatrix_csc), intent(in) :: b !Second addend
    type(sparse_dmatrix_csc)             :: c !Result
    integer                         :: row
    real(8)                         :: val
    integer                         :: i,j
    ! if(.not.a%status)stop "plus_dmatrix error: a is not allocated"
    ! if(.not.b%status)stop "plus_dmatrix error: b is not allocated"
    if(a%Nrow/=b%Nrow)stop "plus_dmatrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_dmatrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = b%col(j)%vals(i)
          call c%insert(val,row,j)
       enddo
    enddo
  end function plus_dmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function minus_dmatrix_csc(a,b) result(c)
    type(sparse_dmatrix_csc), intent(in) :: a !Minuend
    type(sparse_dmatrix_csc), intent(in) :: b !Subtrahend
    type(sparse_dmatrix_csc)             :: c !Result
    integer                         :: row
    real(8)                         :: val
    integer                         :: i,j    
    if(a%Nrow/=b%Nrow)stop "minus_dmatrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "minus_dmatrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = -b%col(j)%vals(i)
          call c%insert(val,row,j)
       enddo
    enddo
  end function minus_dmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function left_product_dmatrix_i_csc(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                        :: row
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function left_product_dmatrix_i_csc
  !
  function left_product_dmatrix_d_csc(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                             :: row
    real(8)                             :: val
    integer                             :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function left_product_dmatrix_d_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_product_dmatrix_i_csc(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                        :: row
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_product_dmatrix_i_csc
  !
  function right_product_dmatrix_d_csc(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                        :: row
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_product_dmatrix_d_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_division_dmatrix_i_csc(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                        :: row
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)/C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_division_dmatrix_i_csc
  !
  function right_division_dmatrix_d_csc(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    type(sparse_dmatrix_csc)            :: B
    integer                        :: row
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)/C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_division_dmatrix_d_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  !function  deye_csc(ndim) result(self)
  !   type(sparse_dmatrix_csc) :: self
  !   integer,intent(in)       :: ndim
  !   integer                  :: i
  !   call self%init(ndim,ndim)
  !   do i=1,ndim
  !      call self%insert(1.d0,i,i)
  !   end do
  ! end func
  !tion deye_csc

  !COMPLEX
  function dgr_zmatrix_csc(a) result(c)
    class(sparse_zmatrix_csc), intent(in) :: a
    type(sparse_zmatrix_csc)              :: c
    integer                               :: row
    complex(8)                            :: val
    integer                               :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = conjg(a%col(j)%vals(i))
          call c%insert(val,j,row)
       enddo
    enddo
  end function dgr_zmatrix_csc
  !
  function transpose_zmatrix_csc(a) result(c)
    class(sparse_zmatrix_csc), intent(in) :: a
    type(sparse_zmatrix_csc)              :: c
    integer                               :: row
    complex(8)                            :: val
    integer                               :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)
          call c%insert(val,j,row)
       enddo
    enddo
  end function transpose_zmatrix_csc
  !
  function conjg_zmatrix_csc(a) result(c)
    class(sparse_zmatrix_csc), intent(in) :: a
    type(sparse_zmatrix_csc)              :: c
    integer                          :: row
    complex(8)                       :: val
    integer                          :: i,j    
    call c%init(a%Nrow,a%Ncol)      !tranpose
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = conjg(a%col(j)%vals(i))
          call c%insert(val,row,j)
       enddo
    enddo
  end function conjg_zmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine zmatrix_equal_zmatrix_csc(a,b)
    type(sparse_zmatrix_csc),intent(inout) :: a
    type(sparse_zmatrix_csc),intent(in)    :: b
    integer                           :: row
    complex(8)                        :: val
    integer                           :: i,j    
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = b%col(j)%vals(i)
          call a%insert(val,row,j)
       enddo
    enddo
  end subroutine zmatrix_equal_zmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine zmatrix_equal_scalar_csc(a,c)
    type(sparse_zmatrix_csc),intent(inout) :: a
    complex(8),intent(in)                :: c
    integer                           :: i,j    
    ! if(.not.a%status)stop "zmatrix_equal_scalar error: a is not allocated"
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          a%col(j)%vals(i) = c
       enddo
    enddo
  end subroutine zmatrix_equal_scalar_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function plus_zmatrix_csc(a,b) result(c)
    type(sparse_zmatrix_csc), intent(in) :: a,b
    type(sparse_zmatrix_csc)             :: c
    integer                         :: row
    complex(8)                      :: val
    integer                         :: i,j
    ! if(.not.a%status)stop "plus_zmatrix error: a is not allocated"
    ! if(.not.b%status)stop "plus_zmatrix error: b is not allocated"
    if(a%Nrow/=b%Nrow)stop "plus_zmatrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_zmatrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = b%col(j)%vals(i)
          call c%insert(val,row,j)
       enddo
    enddo
  end function plus_zmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function minus_zmatrix_csc(a,b) result(c)
    type(sparse_zmatrix_csc), intent(in) :: a,b
    type(sparse_zmatrix_csc)             :: c
    integer                         :: row
    complex(8)                      :: val
    integer                         :: i,j    
    if(a%Nrow/=b%Nrow)stop "plus_zmatrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_zmatrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do j=1,b%Ncol
       do i=1,b%col(j)%size
          row = b%col(j)%rows(i)
          val = -b%col(j)%vals(i)
          call c%insert(val,row,j)
       enddo
    enddo
  end function minus_zmatrix_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function left_product_zmatrix_i_csc(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function left_product_zmatrix_i_csc
  !
  function left_product_zmatrix_d_csc(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function left_product_zmatrix_d_csc
  !
  function left_product_zmatrix_z_csc(C,A) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function left_product_zmatrix_z_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_product_zmatrix_i_csc(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_product_zmatrix_i_csc
  !
  function right_product_zmatrix_d_csc(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_product_zmatrix_d_csc
  !
  function right_product_zmatrix_z_csc(A,C) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)*C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_product_zmatrix_z_csc
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_division_zmatrix_i_csc(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)/C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_division_zmatrix_i_csc
  !
  function right_division_zmatrix_d_csc(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)/C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_division_zmatrix_d_csc
  !
  function right_division_zmatrix_z_csc(A,C) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    type(sparse_zmatrix_csc)            :: B
    integer                        :: row
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do j=1,a%Ncol
       do i=1,a%col(j)%size
          row = a%col(j)%rows(i)
          val = a%col(j)%vals(i)/C
          call b%insert(val,row,j)
       enddo
    enddo
  end function right_division_zmatrix_z_csc

    function left_matmul_dmatrix_darray_csc(C,A) result(B)
    real(8),intent(in),dimension(:)     :: C
    type(sparse_dmatrix_csc),intent(in) :: A
    real(8),dimension(size(C))          :: B
    integer                             :: i,j
    B=0.d0
    if(A%Nrow/=size(C)) stop "size(sparse,1)!=size(array) in left_matmul_dmatrix_darray_csc" 
    do j=1,A%Ncol
       do i=1,A%col(j)%size
          B(j)=B(j)+A%col(j)%vals(i)*C(A%col(j)%rows(i))
       enddo
    end do
  end function left_matmul_dmatrix_darray_csc
  !
  function left_matmul_zmatrix_zarray_csc(C,A) result(B)
    complex(8),intent(in),dimension(:)  :: C
    type(sparse_zmatrix_csc),intent(in) :: A
    complex(8),dimension(size(C))       :: B
    integer                             :: i,j
    B=0.d0
    if(A%Ncol/=size(C)) stop "size(sparse,2)!=size(array) in sparse_matmul_zmatrix_zarray_csc" 
    do j=1,A%Ncol
       do i=1,A%col(j)%size
          B(j)=B(j)+A%col(j)%vals(i)*C(A%col(j)%rows(i))
       enddo
    end do
  end function left_matmul_zmatrix_zarray_csc

  
  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  ! function zeye_csc(ndim) result(self)
  !   type(sparse_zmatrix_csc) :: self
  !   integer,intent(in)       :: ndim
  !   integer                  :: i
  !   call self%init(ndim,ndim)
  !   do i=1,ndim
  !      call self%insert((1.d0,0.d0),i,i)
  !   end do
  ! end function zeye_csc
  



END MODULE SF_SPARSE_ARRAY_CSC
