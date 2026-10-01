!-------------------------------------------------------------------------------------------
!PURPOSE  : Invert a general triangular m*n matrix using LAPACK library 
! M is destroyed and replaces by its inverse M^-1
! M on output is the square matrix n*n
!-------------------------------------------------------------------------------------------
subroutine Dinv_triang(A,uplo,diag)
  !This subroutine inverts in place a real triangular matrix :f:var:`A` using the LAPACK routine
  !:f:func_inline:`dtrtri`. The complex version :f:func_inline:`Zinv_triang` uses :f:func_inline:`ztrtri`. Both are
  !instances of the generic interface :f:func_inline:`inv_triang`.
  !
  !On output the triangle selected by :f:var:`uplo` is replaced by the triangle of the inverse; the other triangle is
  !not referenced. If :code:`diag='U'` the matrix is assumed unit triangular and its diagonal is not referenced.
  !The program stops if the LAPACK routine returns an error code, e.g. for a singular matrix.
  !
  real(8),dimension(:,:)           :: A     ! in: triangular matrix [n,n]; out: its inverse (same triangle)
  character(len=1),optional        :: uplo  ! optional: 'U' (default) upper, 'L' lower triangular
  character(len=1),optional        :: diag  ! optional: 'N' (default) non-unit, 'U' unit triangular
  character(len=1)                 :: uplo_
  character(len=1)                 :: diag_
  integer                          :: n,lda,info
  uplo_="U";if(present(uplo))uplo_=uplo
  diag_="N";if(present(diag))diag_=diag !not a unit triangular matrix
  lda = max(1,size(A,1))
  n   = size(A,2)
  call dtrtri(uplo_,diag_,n,A,lda,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertTRIANG: dtrtri"
end subroutine Dinv_triang
!
subroutine Zinv_triang(A,uplo,diag)
  complex(8),dimension(:,:)           :: A
  character(len=1),optional        :: uplo,diag
  character(len=1)                 :: uplo_
  character(len=1)                 :: diag_
  integer                          :: ndim1,ndim2
  integer                          :: n,lda,info
  uplo_="U";if(present(uplo))uplo_=uplo
  diag_="N";if(present(diag))diag_=diag !not a unit triangular matrix
  lda = max(1,size(A,1))
  n   = size(A,2)
  call ztrtri(uplo_,diag_,n,A,lda,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertTRIANG: ztrtri"
end subroutine Zinv_triang
