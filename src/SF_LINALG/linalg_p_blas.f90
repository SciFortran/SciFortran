subroutine p_d_matmul(A,B,C,Nblock,alfa,beta)
  !This subroutine is the distributed-memory (ScaLAPACK) counterpart of :f:func_inline:`d_matmul`: it computes the real
  !matrix-matrix product :math:`C = \alpha A B + \beta C` with the routine :f:func_inline:`PDGEMM`. The complex version
  !is :f:func_inline:`p_z_matmul`. Both are instances of the generic interface :f:func_inline:`p_mat_product`,
  !available only when compiling with :code:`_SCALAPACK`.
  !
  !The matrices :f:var:`A` and :f:var:`B` (full copies are expected on every process) are distributed on the BLACS
  !process grid of the module (:code:`p_context`, :code:`p_Nx` x :code:`p_Ny`) in a block-cyclic layout with square
  !blocks of size :f:var:`Nblock`, and the product is gathered back into :f:var:`C`. By default :math:`\alpha=1` and
  !:math:`\beta=0`. Timings are written by the master process to the file :code:`p_dgemm.info`. The program stops if
  !:f:var:`B` or :f:var:`C` do not have the shapes :code:`[K,M]` and :code:`[N,M]`.
  !
  real(8),dimension(:,:),intent(inout) :: A       ! matrix [N,K]
  real(8),dimension(:,:),intent(inout) :: B       ! matrix [K,M]
  real(8),dimension(:,:),intent(inout) :: C       ! out: product alfa*A*B, [N,M]
  integer                              :: Nblock  ! block size of the block-cyclic distribution
  real(8),optional                     :: alfa    ! optional scalar, default 1
  real(8),optional                     :: beta    ! optional scalar, default 0
  real(8)                              :: alfa_,beta_
  integer                              :: Nb
  integer                              :: N,K,M
  integer                              :: Nrows,Kcols,Krows,Mcols
  !
  integer                              :: i,j,lda,info
  integer,external                     :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                 :: descAloc,descBloc,descCloc,descC
  real(8),dimension(:,:),allocatable   :: Aloc,Bloc,Cloc
  integer                              :: rankX,rankY
  integer                              :: sendR,sendC,Nipiv
  integer                              :: Nrow,Ncol
  integer                              :: myi,myj,unit
  real(8)                              :: t_stop,t_start
  logical                              :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 ;if(present(alfa))alfa_=alfa
  beta_     = 0d0 ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo(p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 100
  !
  Nb = Nblock
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_dgemm.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Nrows x Kcols=",Nrows,"x",Kcols
     write(unit,"(A20,I8,A5,I8)")"Krows x Mcols=",Krows,"x",Mcols
  endif
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols))
  allocate(Bloc(Krows,Mcols))
  allocate(Cloc(Nrows,Mcols))
  call Distribute_BLACS(A,Aloc,descAloc,unit)
  call Distribute_BLACS(B,Bloc,descBloc,unit)
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A,B:",t_stop-t_start
  !
  if(master)call cpu_time(t_start)
  call PDGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time matmul C_loc:",t_stop-t_start
  !
  C=0d0
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  call Gather_BLACS(Cloc,C,descC,unit)
  !
  if(master)close(unit)
100 continue
  return
  !
end subroutine p_d_matmul

subroutine p_z_matmul(A,B,C,Nblock,alfa,beta)
  complex(8),dimension(:,:),intent(in)                    :: A ![N,K]
  complex(8),dimension(:,:),intent(in)                    :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)),intent(inout) :: C ![N,M]
  integer                                                 :: Nblock
  real(8),optional                                        :: alfa,beta
  real(8)                                                 :: alfa_,beta_
  integer                                                 :: Nb
  integer                                                 :: N,K,M
  integer                                                 :: Nrows,Kcols,Krows,Mcols
  !
  integer                                                 :: i,j,lda,info
  integer,external                                        :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                                    :: descAloc,descBloc,descCloc,descC
  complex(8),dimension(:,:),allocatable                   :: Aloc,Bloc,Cloc
  integer                                                 :: rankX,rankY
  integer                                                 :: sendR,sendC,Nipiv
  integer                                                 :: Nrow,Ncol
  integer                                                 :: myi,myj,unit
  real(8)                                                 :: t_stop,t_start
  logical                                                 :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0) ;if(present(alfa))alfa_=alfa
  beta_     = dcmplx(0d0,0d0) ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  !
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 1010
  !
  Nb = Nblock
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_dgemm.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Nrows x Kcols=",Nrows,"x",Kcols
     write(unit,"(A20,I8,A5,I8)")"Krows x Mcols=",Krows,"x",Mcols
  endif
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=zero
  allocate(Bloc(Krows,Mcols));Bloc=zero
  allocate(Cloc(Nrows,Mcols));Cloc=zero
  call Distribute_BLACS(A,Aloc,descAloc,unit)
  call Distribute_BLACS(B,Bloc,descBloc,unit)
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A,B:",t_stop-t_start
  !
  if(master)call cpu_time(t_start)
  call PZGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time matmul C_loc:",t_stop-t_start
  !
  C=dcmplx(0d0,0d0)
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  call Gather_BLACS(Cloc,C,descC,unit)
  !
  if(master)close(unit)
1010 continue
  return
  !
end subroutine p_z_matmul






!###### OVERLOAD PARALLEL MATMUL OPERATOR --> .px. #########

function p_d_matmul_f(A,B) result(C)
  !This function returns the real matrix-matrix product :math:`C = A B`, computed in parallel with the ScaLAPACK
  !routine :f:func_inline:`PDGEMM`, :math:`A` being :math:`N\times K` and :math:`B` :math:`K\times M`. It is the
  !function behind the operator :code:`.Px.`, i.e. :code:`C = A .Px. B`; the complex version is
  !:f:func_inline:`p_z_matmul_f`. Available only when compiling with :code:`_SCALAPACK`.
  !
  !The block size of the block-cyclic distribution is chosen automatically, as the largest power of 2 smaller than
  !:code:`min(N,K,M)`, up to 64. The product is gathered back into :f:var:`C`. The program stops if :f:var:`B` does not
  !have the shape :code:`[K,M]`.
  !
  real(8),dimension(:,:),intent(in)      :: A  ! matrix [N,K]
  real(8),dimension(:,:),intent(in)      :: B  ! matrix [K,M]
  real(8),dimension(size(A,1),size(B,2)) :: C  ! product A*B, [N,M]
  real(8)                                :: alfa_,beta_
  integer                                :: Nb
  integer                                :: N,K,M
  integer                                :: Nrows,Kcols,Krows,Mcols
  !
  integer                                :: i,j,lda,info
  integer,external                       :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                   :: descAloc,descBloc,descCloc,descC
  real(8),dimension(:,:),allocatable     :: Aloc,Bloc,Cloc
  integer                                :: rankX,rankY
  integer                                :: sendR,sendC,Nipiv
  integer                                :: Nrow,Ncol
  integer                                :: myi,myj,unit
  real(8)                                :: t_stop,t_start
  logical                                :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 
  beta_     = 0d0 
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 201
  !
  I  = min(N,K,M)
  Nb = min(shift_dw(I)/2,64)
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=0d0
  allocate(Bloc(Krows,Mcols));Bloc=0d0
  allocate(Cloc(Nrows,Mcols));Cloc=0d0
  call Distribute_BLACS(A,Aloc,descAloc)
  call Distribute_BLACS(B,Bloc,descBloc)
  !
  call PDGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  C=0d0
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  call Gather_BLACS(Cloc,C,descC)
201 continue
  return
  !
end function p_d_matmul_f

function p_z_matmul_f(A,B) result(C)
  complex(8),dimension(:,:),intent(in)      :: A ![N,K]
  complex(8),dimension(:,:),intent(in)      :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  real(8)                                   :: alfa_,beta_
  integer                                   :: Nb
  integer                                   :: N,K,M
  integer                                   :: Nrows,Kcols,Krows,Mcols
  !
  integer                                   :: i,j,lda,info
  integer,external                          :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                      :: descAloc,descBloc,descCloc,descC
  complex(8),dimension(:,:),allocatable     :: Aloc,Bloc,Cloc
  integer                                   :: rankX,rankY
  integer                                   :: sendR,sendC,Nipiv
  integer                                   :: Nrow,Ncol
  integer                                   :: myi,myj,unit
  real(8)                                   :: t_stop,t_start
  logical                                   :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0)
  beta_     = dcmplx(0d0,0d0)
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 101
  !
  I  = min(N,K,M)
  Nb = min(shift_dw(I)/2,64)
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=zero
  allocate(Bloc(Krows,Mcols));Bloc=zero
  allocate(Cloc(Nrows,Mcols));Cloc=zero
  call Distribute_BLACS(A,Aloc,descAloc)
  call Distribute_BLACS(B,Bloc,descBloc)
  !
  call PZGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  C=dcmplx(0d0,0d0)
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, Nrows, info)
  call Gather_BLACS(Cloc,C,descC)
101 continue
  return
  !
end function p_z_matmul_f
