MODULE SF_SPARSE_ARRAY_CSR
  USE SF_LINALG, only: zeye, deye
  USE SF_SPARSE_COMMON
  !MPI STILL TO BE DONE
#ifdef _MPI
  USE MPI
#endif
  implicit none

  private

  !CSR ROWS
  type sparse_row_csr
     integer                                   :: size
     integer,dimension(:),allocatable          :: cols     
  end type sparse_row_csr
  !
  type, extends(sparse_row_csr) :: sparse_drow_csr
     real(8),dimension(:),allocatable          :: vals
  end type sparse_drow_csr
  !
  type, extends(sparse_row_csr) :: sparse_zrow_csr
     complex(8),dimension(:),allocatable       :: vals
  end type sparse_zrow_csr

  type, extends(sparse_matrix) :: sparse_dmatrix_csr
     type(sparse_drow_csr),dimension(:),allocatable :: row
   contains
     procedure,pass :: init      => init_dmatrix_csr
     procedure,pass :: free      => free_dmatrix_csr
     procedure,pass :: load      => load_dmatrix_csr
     procedure,pass :: dump      => dump_dmatrix_csr
     procedure,pass :: insert    => insert_delement_csr
     procedure,pass :: get       => get_delement_csr
     procedure,pass :: as_matrix => as_dmatrix_csr
     procedure,pass :: dgr       => dgr_dmatrix_csr
     procedure,pass :: transpose => transpose_dmatrix_csr
  end type sparse_dmatrix_csr
  
  type, extends(sparse_matrix) :: sparse_zmatrix_csr
     type(sparse_zrow_csr),dimension(:),allocatable :: row
   contains
     procedure,pass :: init      => init_zmatrix_csr
     procedure,pass :: free      => free_zmatrix_csr
     procedure,pass :: load      => load_zmatrix_csr
     procedure,pass :: dump      => dump_zmatrix_csr
     procedure,pass :: insert    => insert_zelement_csr
     procedure,pass :: get       => get_zelement_csr
     procedure,pass :: as_matrix => as_zmatrix_csr
     procedure,pass :: dgr       => dgr_zmatrix_csr
     procedure,pass :: transpose => transpose_zmatrix_csr
     procedure,pass :: conjg     => conjg_zmatrix_csr
  end type sparse_zmatrix_csr

  
  interface sparse_matrix
     module procedure :: construct_dmatrix_csr
     module procedure :: construct_zmatrix_csr
  end interface sparse_matrix
  !
  interface as_sparse
     module procedure :: construct_dmatrix_csr
     module procedure :: construct_zmatrix_csr
  end interface as_sparse
  !
  interface sparse
     module procedure :: construct_dmatrix_csr
     module procedure :: construct_zmatrix_csr
  end interface sparse
  
  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
     module procedure :: dmatrix_equal_scalar_csr
     module procedure :: dmatrix_equal_dmatrix_csr
     module procedure :: zmatrix_equal_scalar_csr
     module procedure :: zmatrix_equal_zmatrix_csr
  end interface assignment(=)

  !ADDITION
  interface operator(+)
     module procedure :: plus_dmatrix_csr
     module procedure :: plus_zmatrix_csr
  end interface operator(+)

  !SUBTRACTION
  interface operator(-)
     module procedure :: minus_dmatrix_csr
     module procedure :: minus_zmatrix_csr
  end interface operator(-)
  
  !SCALAR PRODUCT
  interface operator(*)
     module procedure :: left_product_dmatrix_i_csr
     module procedure :: left_product_dmatrix_d_csr
     module procedure :: left_product_zmatrix_i_csr
     module procedure :: left_product_zmatrix_d_csr
     module procedure :: left_product_zmatrix_z_csr
     !
     module procedure :: right_product_dmatrix_i_csr
     module procedure :: right_product_dmatrix_d_csr
     module procedure :: right_product_zmatrix_i_csr
     module procedure :: right_product_zmatrix_d_csr
     module procedure :: right_product_zmatrix_z_csr
  end interface operator(*)

  !SCALAR DIVISION
  interface operator(/)
     module procedure :: right_division_dmatrix_i_csr
     module procedure :: right_division_dmatrix_d_csr
     module procedure :: right_division_zmatrix_i_csr
     module procedure :: right_division_zmatrix_d_csr
     module procedure :: right_division_zmatrix_z_csr
  end interface operator(/)

  !MATRIX PRODUCT
  interface matmul
     module procedure :: right_matmul_dmatrix_darray_csr
     module procedure :: right_matmul_zmatrix_zarray_csr 
  end interface matmul
  

  !KRONECKER PRODUCT
  interface operator(.x.)
     module procedure :: kron_dmatrix_csr
     module procedure :: kron_zmatrix_csr
  end interface operator(.x.)

  interface kron
     module procedure :: kron_dmatrix_csr
     module procedure :: kron_zmatrix_csr
     module procedure :: restricted_kron_dmatrix_csr
     module procedure :: restricted_kron_zmatrix_csr
  end interface kron

  interface transpose
     module procedure :: transpose_dmatrix_csr
     module procedure :: transpose_zmatrix_csr
  end interface transpose

  interface hconjg
     module procedure :: dgr_dmatrix_csr
     module procedure :: dgr_zmatrix_csr
  end interface hconjg

  
  public :: sparse_dmatrix_csr, sparse_zmatrix_csr  
  public :: as_sparse
  public :: sparse
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.x.)
  public :: kron
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
  elemental subroutine init_dmatrix_csr(sparse,N,N1)
    class(sparse_dmatrix_csr),intent(inout) :: sparse
    integer,intent(in)                      :: N
    integer,intent(in),optional             :: N1
    integer                                 :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)  call sparse%free() 
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine init_dmatrix_csr
  !
  function construct_dmatrix_csr(matrix) result(self)
    real(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_dmatrix_csr)              :: self
    call self%load(matrix)
  end function construct_dmatrix_csr
  !
  subroutine load_dmatrix_csr(sparse,matrix)
    class(sparse_dmatrix_csr),intent(inout) :: sparse
    real(8),dimension(:,:),intent(in)  :: matrix
    integer                            :: i,j,Ndim1,Ndim2
    !
    call sparse%free() 
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2)
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0d0)call insert_delement_csr(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine load_dmatrix_csr
  !
  !COMPLEX
  elemental subroutine init_zmatrix_csr(sparse,N,N1)
    class(sparse_zmatrix_csr),intent(inout) :: sparse
    integer,intent(in)                      :: N
    integer,intent(in),optional             :: N1
    integer                                 :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)  call sparse%free() 
    !
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine init_zmatrix_csr
  !
  function construct_zmatrix_csr(matrix) result(self)
    complex(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_zmatrix_csr)              :: self
    call self%load(matrix)
  end function construct_zmatrix_csr
  !
  subroutine load_zmatrix_csr(sparse,matrix)
    class(sparse_zmatrix_csr),intent(inout) :: sparse
    complex(8),dimension(:,:),intent(in)  :: matrix
    integer                            :: i,j,Ndim1,Ndim2
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2)
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0d0)call insert_zelement_csr(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine load_zmatrix_csr

  !+------------------------------------------------------------------+
  !PURPOSE: frees a sparse matrix
  !+------------------------------------------------------------------+
  !REAL
  elemental subroutine free_dmatrix_csr(sparse)    
    class(sparse_dmatrix_csr),intent(inout) :: sparse
    integer                            :: i
    !
    do i=1,sparse%Nrow
       sparse%row(i)%Size  = 0
       if(allocated(sparse%row(i)%vals))deallocate(sparse%row(i)%vals)
       if(allocated(sparse%row(i)%cols))deallocate(sparse%row(i)%cols)
    enddo
    if(allocated(sparse%row))deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine free_dmatrix_csr
  !
  !COMPLEX
  elemental subroutine free_zmatrix_csr(sparse)    
    class(sparse_zmatrix_csr),intent(inout) :: sparse
    integer                            :: i
    !
    do i=1,sparse%Nrow
       sparse%row(i)%Size  = 0
       if(allocated(sparse%row(i)%vals))deallocate(sparse%row(i)%vals)
       if(allocated(sparse%row(i)%cols))deallocate(sparse%row(i)%cols)
    enddo
    if(allocated(sparse%row))deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine free_zmatrix_csr


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  !REAL
  subroutine insert_delement_csr(sparse,value,i,j)
    class(sparse_dmatrix_csr),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    integer                               :: column,pos
    logical                               :: iadd
    !
    !
    column = j
    iadd = .false.                          !check if column already exist
    if(any(sparse%row(i)%cols == column))then         !
       pos = binary_search(sparse%row(i)%cols,column) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%row(i)%vals(pos)=sparse%row(i)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it
       call append(sparse%row(i)%vals,value)
       call append(sparse%row(i)%cols,column)
       sparse%row(i)%Size = sparse%row(i)%Size + 1
       sparse%Nelements = sparse%Nelements+1
    endif
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "insert_delement_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine insert_delement_csr
  !
  !COMPLEX
  subroutine insert_zelement_csr(sparse,value,i,j)
    class(sparse_zmatrix_csr),intent(inout) :: sparse
    complex(8),intent(in)                 :: value
    integer,intent(in)                    :: i,j
    !
    integer                               :: column,pos
    logical                               :: iadd
    !
    !
    column = j
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%row(i)%cols == column))then         !
       pos = binary_search(sparse%row(i)%cols,column) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%row(i)%vals(pos)=sparse%row(i)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call append(sparse%row(i)%vals,value)
       call append(sparse%row(i)%cols,column)
       sparse%row(i)%Size = sparse%row(i)%Size + 1
       sparse%Nelements = sparse%Nelements+1
    endif
    !
    if(sparse%row(i)%Size > sparse%Ncol)stop "insert_zelement_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine insert_zelement_csr
  
  !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  function get_delement_csr(sparse,i,j) result(value)
    class(sparse_dmatrix_csr),intent(inout) :: sparse    
    integer,intent(in)                 :: i,j
    real(8)                            :: value
    integer                            :: pos
    value=0d0
    do pos=1,sparse%row(i)%size
       if(j==sparse%row(i)%cols(pos))value=sparse%row(i)%vals(pos)
    enddo
  end function get_delement_csr

  function get_zelement_csr(sparse,i,j) result(value)
    class(sparse_zmatrix_csr),intent(inout) :: sparse    
    integer,intent(in)                 :: i,j
    complex(8)                         :: value
    integer                            :: pos
    value=0d0
    do pos=1,sparse%row(i)%size
       if(j==sparse%row(i)%cols(pos))value=sparse%row(i)%vals(pos)
    enddo
  end function get_zelement_csr


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  !REAL
  subroutine dump_dmatrix_csr(sparse,matrix)
    class(sparse_dmatrix_csr),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    !
    !
    Ndim1=size(matrix,1);  Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine dump_dmatrix_csr
  !
  !COMPLEX
  subroutine dump_zmatrix_csr(sparse,matrix)
    class(sparse_zmatrix_csr),intent(in)      :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    integer                                 :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1);   Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine dump_zmatrix_csr

  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  !REAL
  function as_dmatrix_csr(sparse) result(matrix)
    class(sparse_dmatrix_csr),intent(in)            :: sparse
    real(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
    integer                                    :: i,j
    matrix = 0d0
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = sparse%row(i)%vals(j)
       enddo
    enddo
  end function as_dmatrix_csr
  
  !COMPLEX
  function as_zmatrix_csr(sparse) result(matrix)
    class(sparse_zmatrix_csr),intent(in)            :: sparse
    complex(8),dimension(sparse%Nrow,sparse%Ncol) :: matrix
    integer                                    :: i,j
    matrix = 0d0
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = sparse%row(i)%vals(j)
       enddo
    enddo
  end function as_zmatrix_csr
  
  !##################################################################
  !               SPARSE MATRIX KRON PRODUCT 
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
  !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  !+------------------------------------------------------------------+
  !REAL
  function kron_dmatrix_csr(A,B) result(AxB)
    type(sparse_dmatrix_csr), intent(in) :: A,B
    type(sparse_dmatrix_csr)             :: AxB
    integer                         :: i,icol,j,k,kcol,l
    integer                         :: indx_row,indx_col
    real(8)                         :: value
!    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    do indx_row = 1,A%Nrow*B%Nrow
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(indx_row)%vals,value)
             call append(AxB%row(indx_row)%cols,indx_col)
             AxB%row(indx_row)%Size = AxB%row(indx_row)%Size + 1
             AxB%Nelements = AxB%Nelements+1
             !
          enddo
       enddo
    enddo
  end function kron_dmatrix_csr
  !
  function restricted_kron_dmatrix_csr(A,B,states) result(AxB)
    type(sparse_dmatrix_csr), intent(in)  :: A,B
    integer,dimension(:),intent(in)  :: states
    type(sparse_dmatrix_csr)      :: AxB
    integer                          :: i,icol,j,k,kcol,l,istate,jstate
    integer                          :: indx_row,indx_col
    real(8)                          :: value
    integer,dimension(:),allocatable :: inv_states
!    call AxB%free()
    call AxB%init(size(states),size(states))
    allocate(inv_states(A%Ncol*B%Ncol))
    inv_states=0
    do i=1,size(states)
       inv_states(states(i)) = i
    enddo
    do istate = 1,size(states)
       indx_row=states(istate)
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             jstate   = inv_states(indx_col)
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(istate)%vals,value)
             call append(AxB%row(istate)%cols,jstate)
             AxB%row(istate)%Size = AxB%row(istate)%Size + 1
             AxB%Nelements = AxB%Nelements + 1
             !
          enddo
       enddo
    enddo
  end function restricted_kron_dmatrix_csr
  !
  !COMPLEX
  function kron_zmatrix_csr(A,B) result(AxB)
    type(sparse_zmatrix_csr), intent(in) :: A,B
    type(sparse_zmatrix_csr)             :: AxB
    integer                         :: i,icol,j,k,kcol,l
    integer                         :: indx_row,indx_col
    complex(8)                      :: value
!    call AxB%free()
    call AxB%init(a%Nrow*b%Nrow,a%Ncol*b%Ncol)
    do indx_row = 1,A%Nrow*B%Nrow
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(indx_row)%vals,value)
             call append(AxB%row(indx_row)%cols,indx_col)
             AxB%row(indx_row)%Size = AxB%row(indx_row)%Size + 1
             AxB%Nelements = AxB%Nelements+1
             !
          enddo
       enddo
    enddo
  end function kron_zmatrix_csr
  !
  function restricted_kron_zmatrix_csr(A,B,states) result(AxB)
    type(sparse_zmatrix_csr), intent(in)  :: A,B
    integer,dimension(:),intent(in)  :: states
    type(sparse_zmatrix_csr)      :: AxB
    integer                          :: i,icol,j,k,kcol,l,istate,jstate
    integer                          :: indx_row,indx_col
    complex(8)                          :: value
    integer,dimension(:),allocatable :: inv_states
!    call AxB%free()
    call AxB%init(size(states),size(states))
    allocate(inv_states(A%Ncol*B%Ncol))
    inv_states=0
    do i=1,size(states)
       inv_states(states(i)) = i
    enddo
    do istate = 1,size(states)
       indx_row=states(istate)
       k = mod(indx_row,B%Nrow);if(k==0)k=B%Nrow
       i = (indx_row-1)/B%Nrow+1
       do icol=1,A%row(i)%size
          j = A%row(i)%cols(icol)
          do kcol=1,B%row(k)%size
             l = B%row(k)%cols(kcol)
             indx_col = l + (j-1)*B%Ncol
             jstate   = inv_states(indx_col)
             value    = A%row(i)%vals(icol)*B%row(k)%vals(kcol)
             !
             call append(AxB%row(istate)%vals,value)
             call append(AxB%row(istate)%cols,jstate)
             AxB%row(istate)%Size = AxB%row(istate)%Size + 1
             AxB%Nelements = AxB%Nelements + 1
             !
          enddo
       enddo
    enddo
  end function restricted_kron_zmatrix_csr
  !
  

  
  !##################################################################
  !               SPARSE MATRIX BASIC ALGEBRA 
  !##################################################################
  !REAL
  function dgr_dmatrix_csr(a) result(c)
    class(sparse_dmatrix_csr), intent(in) :: a
    type(sparse_dmatrix_csr)              :: c
    integer                          :: col
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)
          call c%insert(val,col,i)
       enddo
    enddo
  end function dgr_dmatrix_csr
  !
  function transpose_dmatrix_csr(a) result(c)
    class(sparse_dmatrix_csr), intent(in) :: a
    type(sparse_dmatrix_csr)              :: c
    integer                          :: col
    real(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)
          call c%insert(val,col,i)
       enddo
    enddo
  end function transpose_dmatrix_csr
  !
  
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine dmatrix_equal_dmatrix_csr(a,b)
    type(sparse_dmatrix_csr),intent(inout) :: a
    type(sparse_dmatrix_csr),intent(in)    :: b
    integer                           :: col
    real(8)                           :: val
    integer                           :: i,j    
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call a%insert(val,i,col)
       enddo
    enddo
  end subroutine dmatrix_equal_dmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine dmatrix_equal_scalar_csr(a,c)
    type(sparse_dmatrix_csr),intent(inout) :: a
    real(8),intent(in)                :: c
    integer                           :: i,j    
    ! if(.not.a%status)stop "matrix_equal_scalar error: a is not allocated"
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          a%row(i)%vals(j) = c
       enddo
    enddo
  end subroutine dmatrix_equal_scalar_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function plus_dmatrix_csr(a,b) result(c)
    type(sparse_dmatrix_csr), intent(in) :: a !First addendum
    type(sparse_dmatrix_csr), intent(in) :: b !Second addendum
    type(sparse_dmatrix_csr)             :: c !Result
    integer                         :: col
    real(8)                         :: val
    integer                         :: i,j
    ! if(.not.a%status)stop "plus_matrix error: a is not allocated"
    ! if(.not.b%status)stop "plus_matrix error: b is not allocated"
    if(a%Nrow/=b%Nrow)stop "plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function plus_dmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function minus_dmatrix_csr(a,b) result(c)
    type(sparse_dmatrix_csr), intent(in) :: a !Minuend
    type(sparse_dmatrix_csr), intent(in) :: b !Subtrahend
    type(sparse_dmatrix_csr)             :: c !Result
    integer                         :: col
    real(8)                         :: val
    integer                         :: i,j    
    if(a%Nrow/=b%Nrow)stop "plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val =-b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function minus_dmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function left_product_dmatrix_i_csr(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function left_product_dmatrix_i_csr
  !
  function left_product_dmatrix_d_csr(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function left_product_dmatrix_d_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_product_dmatrix_i_csr(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_product_dmatrix_i_csr
  !
  function right_product_dmatrix_d_csr(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_product_dmatrix_d_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_division_dmatrix_i_csr(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_division_dmatrix_i_csr
  !
  function right_division_dmatrix_d_csr(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B
    integer                        :: col
    real(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_division_dmatrix_d_csr

  function right_matmul_dmatrix_darray_csr(A,C) result(B)
    real(8),intent(in),dimension(:)     :: C
    type(sparse_dmatrix_csr),intent(in) :: A
    real(8),dimension(size(C))          :: B
    integer                             :: i,j
    B=0.d0
    if(A%Ncol/=size(C)) stop "size(sparse,2)!=size(array) in right_matmul_dmatrix_darray_csr" 
    do i=1,A%Nrow
       do j=1,A%row(i)%size
          B(i)=B(i)+A%row(i)%vals(j)*C(A%row(i)%cols(j))
       enddo
    end do
  end function right_matmul_dmatrix_darray_csr
  !
  function right_matmul_zmatrix_zarray_csr(A,C) result(B)
    complex(8),intent(in),dimension(:)  :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    complex(8),dimension(size(C))       :: B
    integer                             :: i,j
    B=0.d0
    if(A%Ncol/=size(C)) stop "size(sparse,2)!=size(array) in right_matmul_zmatrix_zarray_csr" 
    do i=1,A%Nrow
       do j=1,A%row(i)%size
          B(i)=B(i)+A%row(i)%vals(j)*C(A%row(i)%cols(j))
       enddo
    end do
  end function right_matmul_zmatrix_zarray_csr

  
  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  ! function deye_csr(ndim) result(self)
  !   type(sparse_dmatrix_csr) :: self
  !   integer,intent(in)       :: ndim
  !   integer                  :: i
  !   call self%init(ndim,ndim)
  !   do i=1,ndim
  !      call self%insert(1.d0,i,i)
  !   end do
  ! end function deye_csr

  !COMPLEX
  function dgr_zmatrix_csr(a) result(c)
    class(sparse_zmatrix_csr), intent(in) :: a
    type(sparse_zmatrix_csr)              :: c
    integer                          :: col
    complex(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = conjg(a%row(i)%vals(j))
          call c%insert(val,col,i)
       enddo
    enddo
  end function dgr_zmatrix_csr
  !
  function transpose_zmatrix_csr(a) result(c)
    class(sparse_zmatrix_csr), intent(in) :: a
    type(sparse_zmatrix_csr)              :: c
    integer                          :: col
    complex(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Ncol,a%Nrow)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)
          call c%insert(val,col,i)
       enddo
    enddo
  end function transpose_zmatrix_csr
  !
  function conjg_zmatrix_csr(a) result(c)
    class(sparse_zmatrix_csr), intent(in) :: a
    type(sparse_zmatrix_csr)              :: c
    integer                          :: col
    complex(8)                          :: val
    integer                          :: i,j    
    call c%init(a%Nrow,a%Ncol)      !tranpose
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = conjg(a%row(i)%vals(j))
          call c%insert(val,i,col)
       enddo
    enddo
  end function conjg_zmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+
  subroutine zmatrix_equal_zmatrix_csr(a,b)
    type(sparse_zmatrix_csr),intent(inout) :: a
    type(sparse_zmatrix_csr),intent(in)    :: b
    integer                           :: col
    complex(8)                        :: val
    integer                           :: i,j    
    call a%free()
    call a%init(b%Nrow,b%Ncol)
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call a%insert(val,i,col)
       enddo
    enddo
  end subroutine zmatrix_equal_zmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine zmatrix_equal_scalar_csr(a,c)
    type(sparse_zmatrix_csr),intent(inout) :: a
    complex(8),intent(in)                :: c
    integer                           :: i,j    
    ! if(.not.a%status)stop "matrix_equal_scalar error: a is not allocated"
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          a%row(i)%vals(j) = c
       enddo
    enddo
  end subroutine zmatrix_equal_scalar_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function plus_zmatrix_csr(a,b) result(c)
    type(sparse_zmatrix_csr), intent(in) :: a,b
    type(sparse_zmatrix_csr)             :: c
    integer                         :: col
    complex(8)                         :: val
    integer                         :: i,j
    ! if(.not.a%status)stop "plus_matrix error: a is not allocated"
    ! if(.not.b%status)stop "plus_matrix error: b is not allocated"
    if(a%Nrow/=b%Nrow)stop "plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val = b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function plus_zmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function minus_zmatrix_csr(a,b) result(c)
    type(sparse_zmatrix_csr), intent(in) :: a,b
    type(sparse_zmatrix_csr)             :: c
    integer                         :: col
    complex(8)                         :: val
    integer                         :: i,j    
    if(a%Nrow/=b%Nrow)stop "plus_matrix error: a.Nrow != b.Nrow"
    if(a%Ncol/=b%Ncol)stop "plus_matrix error: a.Ncol != b.Ncol"
    c=a                         !copy a into c
    do i=1,b%Nrow
       do j=1,b%row(i)%size
          col = b%row(i)%cols(j)
          val =-b%row(i)%vals(j)
          call c%insert(val,i,col)
       enddo
    enddo
  end function minus_zmatrix_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function left_product_zmatrix_i_csr(C,A) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function left_product_zmatrix_i_csr
  !
  function left_product_zmatrix_d_csr(C,A) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function left_product_zmatrix_d_csr
  !
  function left_product_zmatrix_z_csr(C,A) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function left_product_zmatrix_z_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_product_zmatrix_i_csr(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_product_zmatrix_i_csr
  !
  function right_product_zmatrix_d_csr(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_product_zmatrix_d_csr
  !
  function right_product_zmatrix_z_csr(A,C) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)*C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_product_zmatrix_z_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function right_division_zmatrix_i_csr(A,C) result(B)
    integer,intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                     :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_division_zmatrix_i_csr
  !
  function right_division_zmatrix_d_csr(A,C) result(B)
    real(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_division_zmatrix_d_csr
  !
  function right_division_zmatrix_z_csr(A,C) result(B)
    complex(8),intent(in)             :: C
    type(sparse_zmatrix_csr),intent(in) :: A
    type(sparse_zmatrix_csr)            :: B
    integer                        :: col
    complex(8)                        :: val
    integer                        :: i,j   
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          val = a%row(i)%vals(j)/C
          call b%insert(val,i,col)
       enddo
    enddo
  end function right_division_zmatrix_z_csr
  !
  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  ! !+------------------------------------------------------------------+
  ! function zeye_csr(ndim) result(self)
  !   type(sparse_zmatrix_csr) :: self
  !   integer,intent(in)       :: ndim
  !   integer                  :: i
  !   call self%init(ndim,ndim)
  !   do i=1,ndim
  !      call self%insert((1.d0,0.d0),i,i)
  !   end do
  ! end function zeye_csr
    
END MODULE SF_SPARSE_ARRAY_CSR
