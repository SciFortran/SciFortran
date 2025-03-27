MODULE SF_SPARSE_ARRAY_ALGEBRA
  USE SF_SPARSE_COMMON
  USE SF_SPARSE_ARRAY_CSC
  USE SF_SPARSE_ARRAY_CSR
  implicit none

  interface matmul
     module procedure :: dmatmul_csr_csr
     module procedure :: zmatmul_csr_csr
     !
     module procedure :: dmatmul_csc_csc
     module procedure :: zmatmul_csc_csc
     !
     module procedure :: dmatmul_csc_csr_2csr
     module procedure :: zmatmul_csc_csr_2csr
     !
     module procedure :: dmatmul_csc_csr_2csc
     module procedure :: zmatmul_csc_csr_2csc
  end interface matmul

  public :: matmul

  
contains

  function dmatmul_csr_csr(A,B) result(AxB)
    type(sparse_dmatrix_csr), intent(in) :: A,B
    type(sparse_dmatrix_csr)             :: AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: irow,j,jcol,k,kcol
    real(8)                              :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in dmatmul_csr_csr"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do irow=1,Na(1)
       do j=1,A%row(irow)%Size
          jcol=A%row(irow)%cols(j)
          aval=A%row(irow)%vals(j)
          do k=1,B%row(jcol)%Size
             kcol=B%row(jcol)%cols(k)
             bval=B%row(jcol)%vals(k)
             AxB%insert(aval*bval,irow,kcol)
          end do
       end do
    end do
  end function dmatmul_csr_csr

  
  function zmatmul_csr_csr(A,B) result(AxB)
    type(sparse_zmatrix_csr), intent(in) :: A,B
    type(sparse_zmatrix_csr)             :: AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: irow,j,jcol,k,kcol
    complex(8)                           :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csr_csr"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do irow=1,Na(1)
       do j=1,A%row(irow)%Size
          jcol=A%row(irow)%cols(j)
          aval=A%row(irow)%vals(j)
          do k=1,B%row(jcol)%Size
             kcol=B%row(jcol)%cols(k)
             bval=B%row(jcol)%vals(k)
             AxB%insert(aval*bval,irow,kcol)
          end do
       end do
    end do
  end function zmatmul_csr_csr

  
  function dmatmul_csc_csc(A,B) result(AxB)
    type(sparse_dmatrix_csc), intent(in) :: A,B
    type(sparse_dmatrix_csc)             :: AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    real(8)                              :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in dmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do icol=1,Nb(2)
       do j=1,B%col(icol)%Size
          jrow=B%col(icol)%rows(j)
          bval=B%col(icol)%vals(j)
          do k=1,A%col(jrow)%Size
             krow=A%col(jrow)%rows(k)
             aval=A%col(jrow)%vals(k)
             AxB%insert(aval*bval,krow,icol)
          end do
       end do
    end do
  end function dmatmul_csc_csc

  function zmatmul_csc_csc(A,B) result(AxB)
    type(sparse_zmatrix_csc), intent(in) :: A,B
    type(sparse_zmatrix_csc)             :: AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    complex(8)                              :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do icol=1,Nb(2)
       do j=1,B%col(icol)%Size
          jrow=B%col(icol)%rows(j)
          bval=B%col(icol)%vals(j)
          do k=1,A%col(jrow)%Size
             krow=A%col(jrow)%rows(k)
             aval=A%col(jrow)%vals(k)
             AxB%insert(aval*bval,krow,icol)
          end do
       end do
    end do
  end function zmatmul_csc_csc

  
  
  function dmatmul_csc_csr_2csc(A,B) result(AxB)
    type(sparse_dmatrix_csc), intent(in) :: A,AxB
    type(sparse_dmatrix_csr)             :: B
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    real(8)                              :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do i=1,Na(2)
       do j=1,A%col(i)%Size
          jrow=A%col(i)%rows(j)
          aval=A%col(i)%vals(j)
          do k=1,B%row(i)%Size
             kcol=B%row(i)%cols(k)
             bval=B%row(i)%vals(k)
             AxB%insert(aval*bval,jrow,kcol)
          end do
       end do
    end do
  end function zmatmul_csc_csr_2csc

  
  function zmatmul_csc_csr_2csc(A,B) result(AxB)
    type(sparse_zmatrix_csc), intent(in) :: A,AxB
    type(sparse_zmatrix_csr)             :: B
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    complex(8)                           :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do i=1,Na(2)
       do j=1,A%col(i)%Size
          jrow=A%col(i)%rows(j)
          aval=A%col(i)%vals(j)
          do k=1,B%row(i)%Size
             kcol=B%row(i)%cols(k)
             bval=B%row(i)%vals(k)
             AxB%insert(aval*bval,jrow,kcol)
          end do
       end do
    end do
  end function zmatmul_csc_csr_2csc
  
  function dmatmul_csc_csr_2csr(A,B) result(AxB)
    type(sparse_dmatrix_csc), intent(in) :: A
    type(sparse_dmatrix_csr)             :: B,AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    real(8)                              :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do i=1,Na(2)
       do j=1,A%col(i)%Size
          jrow=A%col(i)%rows(j)
          aval=A%col(i)%vals(j)
          do k=1,B%row(i)%Size
             kcol=B%row(i)%cols(k)
             bval=B%row(i)%vals(k)
             AxB%insert(aval*bval,jrow,kcol)
          end do
       end do
    end do
  end function dmatmul_csc_csr_2csr

  
  function zmatmul_csc_csr_2csr(A,B) result(AxB)
    type(sparse_zmatrix_csc), intent(in) :: A
    type(sparse_zmatrix_csr)             :: B,AxB
    integer                              :: Na(2),Nb(2)
    integer                              :: icol,j,jrow,k,krow
    complex(8)                           :: aval,bval
    
    Na = A%shape(); Nb = B%shape()
    if(Na(2)/=Nb(1))stop "Matrix not matching dimension in zmatmul_csc_csc"
    call AxB%free()
    call AxB%init(Na(1),Nb(2))
    do i=1,Na(2)
       do j=1,A%col(i)%Size
          jrow=A%col(i)%rows(j)
          aval=A%col(i)%vals(j)
          do k=1,B%row(i)%Size
             kcol=B%row(i)%cols(k)
             bval=B%row(i)%vals(k)
             AxB%insert(aval*bval,jrow,kcol)
          end do
       end do
    end do
  end function zmatmul_csc_csr_2csr
  
END MODULE SF_SPARSE_ARRAY_ALGEBRA
