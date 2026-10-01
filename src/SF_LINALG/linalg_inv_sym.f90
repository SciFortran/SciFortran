!-------------------------------------------------------------------------------------------
!PURPOSE  : Invert a symmetric n*n matrix using LAPACK library 
! M is destroyed and replaces by its inverse M^-1
!-------------------------------------------------------------------------------------------
subroutine Dinv_sym(A,uplo)
  !This subroutine inverts in place a real symmetric matrix :f:var:`A` using the Bunch-Kaufman factorization,
  !:f:func_inline:`dsytrf` + :f:func_inline:`dsytri`. The complex version :f:func_inline:`Zinv_sym` (complex symmetric,
  !NOT hermitian, matrices; see :f:func_inline:`zinv_her`) uses :f:func_inline:`zsytrf` + :f:func_inline:`zsytri`. Both
  !are instances of the generic interface :f:func_inline:`inv_sym`.
  !
  !Only the triangle selected by :f:var:`uplo` is referenced on input. On output :f:var:`A` contains the full inverse:
  !the triangle computed by LAPACK is mirrored onto the other one. The program stops if any of the LAPACK routines
  !returns an error code.
  !
  real(8),dimension(:,:)           :: A     ! in: symmetric matrix; out: its inverse (full matrix)
  character(len=*),optional        :: uplo  ! optional: 'U' (default) or 'L', triangle of A referenced
  character(len=1)                 :: uplo_
  integer                          :: n,lda,info,lwork,i,j
  integer,dimension(:),allocatable :: ipvt
  real(8),dimension(:),allocatable :: work
  real(8),dimension(1)             :: lwork_guess
  uplo_="U";if(present(uplo))uplo_=uplo
  lda  = max(1,size(A,1))
  n    = size(A,2)
  allocate(ipvt(n))
  call dsytrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 1st call dsytrf"
  lwork=lwork_guess(1);allocate(work(lwork))
  call dsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 2nd call dsytrf"
  call dsytri(uplo_,n,A,lda,ipvt,work,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: dsytri"
  deallocate(ipvt,work)
  if(uplo_=="U")then
     forall(i=1:size(A,1),j=1:size(A,2),i>j)A(i,j)=A(j,i)
  elseif(uplo_=="L")then
     forall(i=1:size(A,1),j=1:size(A,2),i<j)A(i,j)=A(j,i)
  endif
  !
end subroutine Dinv_SYM
!-----------------------------
subroutine Zinv_sym(A,uplo)
  complex(8),dimension(:,:)           :: A
  character(len=*),optional           :: uplo
  character(len=1)                    :: uplo_
  integer                             :: n,lda,info,lwork,i,j
  integer,dimension(:),allocatable    :: ipvt
  complex(8),dimension(:),allocatable :: work
  complex(8),dimension(1)             :: lwork_guess
  uplo_="U";if(present(uplo))uplo_=uplo
  lda  = max(1,size(A,1))
  n    = size(A,2)
  allocate(ipvt(n))
  call zsytrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 1st call zsytrf"
  lwork=lwork_guess(1);allocate(work(lwork))
  call zsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 2nd call zsytrf"
  call zsytri(uplo_,n,A,lda,ipvt,work,info)
  if(info/=0)stop "Error MATRIX/D_mat_invertSYM: zsytri"
  deallocate(ipvt,work)
  if(uplo_=="U")then
     forall(i=1:size(A,1),j=1:size(A,2),i>j)A(i,j)=A(j,i)
  elseif(uplo_=="L")then
     forall(i=1:size(A,1),j=1:size(A,2),i<j)A(i,j)=A(j,i)
  endif
end subroutine Zinv_sym
