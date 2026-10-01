subroutine p_Dinv(A,Nblock)
  !This subroutine is the distributed-memory (ScaLAPACK) counterpart of :f:func_inline:`dinv`: it inverts in place a
  !general real square matrix :f:var:`A` with the LU factorization, :f:func_inline:`PDGETRF` + :f:func_inline:`PDGETRI`.
  !The complex version :f:func_inline:`p_Zinv` uses :f:func_inline:`PZGETRF` + :f:func_inline:`PZGETRI`. Both are
  !instances of the generic interface :f:func_inline:`p_inv`, available only when compiling with :code:`_SCALAPACK`.
  !
  !The matrix :f:var:`A` (a full copy is expected on every process) is distributed on the BLACS process grid of the
  !module (:code:`p_context`, :code:`p_Nx` x :code:`p_Ny`) in a block-cyclic layout with square blocks of size
  !:f:var:`Nblock`; the inverse is gathered back into :f:var:`A`. Timings are written by the master process to the
  !file :code:`p_inv.info`. The program stops if :f:var:`A` is not square, or if a ScaLAPACK routine fails.
  !
  real(8),dimension(:,:),intent(inout)       :: A       ! in: square matrix to be inverted; out: its inverse
  integer                                    :: Nblock  ! block size of the block-cyclic distribution
  integer                                    :: Nb
  integer                                    :: Ns
  integer                                    :: Qrows,Qcols
  integer                                    :: i,j,lda,info
  !
  real(8)                                    :: guess_lwork(1)
  integer                                    :: guess_liwork(1)
  real(8),dimension(:),allocatable           :: work
  integer,dimension(:),allocatable           :: Iwork 
  integer                                    :: lwork
  integer                                    :: liwork
  integer, allocatable                       :: ipiv(:)
  !
  integer,external                           :: numroc,indxG2L,indxG2P,indxL2G
  real(8),external                           :: dlamch,Pilaenvx
  !
  real(8),dimension(:,:),allocatable         :: A_loc,Z_loc
  integer                                    :: rankX,rankY
  integer                                    :: sendR,sendC,Nipiv
  integer                                    :: Nrow,Ncol
  integer                                    :: myi,myj,unit
  integer,dimension(9)                       :: descA,descAloc,descZloc
  real(8)                                    :: t_stop,t_start
  logical                                    :: master
  !
  !
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  !
  if(rankX<0.AND.rankY<0)goto 100
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_inv.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Qrows x Qcols=",Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols))
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  !< Distribute A
  call Distribute_BLACS(A,A_loc,descAloc,unit)
  !
  if(master)call cpu_time(t_start)
  allocate(Ipiv(Qrows*Nb))
  call PDGETRF(Ns, Ns, A_loc, 1, 1, descAloc, IPIV, INFO)
  if(info /= 0) then
     print*, "PDGETRF returned info =", info
     stop
  end if
  !
  !
  call PDGETRI( Ns, A_loc, 1, 1, descAloc, IPIV, guess_lWORK, -1, guess_LIWORK, -1, INFO )
  lwork = guess_lwork(1)
  liwork= guess_liwork(1)
  allocate(work(lwork))
  allocate(iwork(liwork))
  call PDGETRI( Ns, A_loc, 1, 1, descAloc, IPIV, Work, Lwork, Iwork, LIwork, INFO )
  if(info /= 0) then
     print*, "PDGETRI ERROR. returned info =", info
     stop
  end if
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Invert :",t_stop-t_start
  !
  A=0d0
  call Gather_BLACS(A_loc,A,descA,unit)
  !
  if(master)close(unit)
  !
100 continue
  return
  !
end subroutine p_Dinv





subroutine p_Zinv(A,Nblock)
  complex(8),dimension(:,:),intent(inout) :: A
  integer                                 :: Nblock
  integer                                 :: Nb
  integer                                 :: Ns
  integer                                 :: Qrows,Qcols
  integer                                 :: i,j,lda,info
  !
  complex(8)                              :: guess_lwork(1)
  integer                                 :: guess_liwork(1)
  complex(8),dimension(:),allocatable     :: work
  integer,dimension(:),allocatable        :: Iwork 
  integer                                 :: lwork
  integer                                 :: liwork
  integer, allocatable                    :: ipiv(:)
  !
  integer,external                        :: numroc,indxG2L,indxG2P,indxL2G
  !
  complex(8),dimension(:,:),allocatable   :: A_loc,Z_loc
  integer                                 :: rankX,rankY
  integer                                 :: sendR,sendC,Nipiv
  integer                                 :: Nrow,Ncol
  integer                                 :: myi,myj,unit
  integer,dimension(9)                    :: descA,descAloc,descZloc
  real(8)                                 :: t_stop,t_start
  logical                                 :: master
  !
  !
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  !
  if(rankX<0.AND.rankY<0)goto 200
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_inv.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Qrows x Qcols=",Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols))
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  !< Distribute A
  call Distribute_BLACS(A,A_loc,descAloc,unit)
  !
  !< Allocate distributed eigenvector matrix
  allocate(Z_loc(Qrows,Qcols));Z_loc=0d0
  call descinit( descZloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)call cpu_time(t_start)
  allocate(Ipiv(Qrows*Nb))
  call PZGETRF(Ns, Ns, A_loc, 1, 1, descAloc, IPIV, INFO)
  if(info /= 0) then
     print*, "PZGETRF returned info =", info
     stop
  end if
  !
  !
  call PZGETRI( Ns, A_loc, 1, 1, descAloc, IPIV, guess_lWORK, -1, guess_LIWORK, -1, INFO )
  lwork = guess_lwork(1)
  liwork= guess_liwork(1)
  allocate(work(lwork))
  allocate(iwork(liwork))
  call PZGETRI( Ns, A_loc, 1, 1, descAloc, IPIV, Work, Lwork, Iwork, LIwork, INFO )
  if(info /= 0) then
     print*, "PZGETRI returned info =", info
     stop
  end if
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time inv A_loc:",t_stop-t_start
  !
  A=0d0
  call Gather_BLACS(A_loc,A,descA,unit)
  if(master)close(unit)
200 continue
  return
  !
end subroutine p_Zinv
