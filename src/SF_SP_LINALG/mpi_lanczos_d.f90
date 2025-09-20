!---------------------------------------------------------------------
!Purpose: use plain lanczos to get the groundstate energy
!---------------------------------------------------------------------
subroutine mpi_lanczos_eigh_d(MpiComm,MatVec,Egs,Vect,Nitermax,iverbose,threshold,ncheck,vrandom,NumOp)
  integer                              :: MpiComm
  interface 
     subroutine MatVec(Nloc,vin,vout)
       integer :: Nloc
       real(8) :: vin(Nloc)
       real(8) :: vout(Nloc)
     end subroutine MatVec
  end interface
  real(8)                              :: egs
  real(8),dimension(:)                 :: vect !Nloc
  integer                              :: Nitermax
  !
  integer                              :: Nloc
  real(8),optional                     :: threshold
  integer,optional                     :: ncheck,NumOp
  logical,optional                     :: iverbose
  logical,optional                     :: vrandom
  !
  real(8),dimension(size(vect))        :: vin,vout
  integer                              :: iter,nlanc,nhxv
  real(8),dimension(Nitermax+1)        :: alanc,blanc
  real(8),dimension(Nitermax,Nitermax) :: Z
  real(8),dimension(Nitermax)          :: diag,subdiag,esave
  real(8)                              :: a_,b_,norm,diff,norm_tmp
  integer                              :: i,ierr
  logical                              :: vran=.true.
  !
  logical                              :: mpi_master
  !
  if(MpiComm==Mpi_Comm_Null)return
  !
  mpi_master=get_master_MPI(MpiComm)
  !
  Nloc = size(vect)
  !
  if(present(iverbose))verb=iverbose
  if(present(threshold))threshold_=threshold
  if(present(ncheck))ncheck_=ncheck
  if(present(vrandom))vran=vrandom
  !
  norm_tmp=dot_product(vect,vect); norm=0d0
  call AllReduce_MPI(MpiComm,norm_tmp,norm)
  !
  if(norm==0d0)then
     if(vran)then
        ! call random_seed(size=nrandom)
        ! if(allocated(seed_random))deallocate(seed_random)
        ! allocate(seed_random(nrandom))
        ! seed_random=1234567
        ! call random_seed(put=seed_random)
        ! deallocate(seed_random)
        ! call random_number(vect)
        call mt_random(vect)
        if(verb.AND.mpi_master)write(*,*)"MPI_LANCZOS_EIGH: random initial vector generated:"
     else
        vect = 1d0
        if(verb.AND.mpi_master)write(*,*)"MPI_LANCZOS_EIGH: unitary initial vector generated:"
     endif
     norm_tmp=dot_product(vect,vect); norm=0d0
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     vect=vect/sqrt(norm)
  endif
  !
  !============= LANCZOS LOOP =====================
  !
  Vin  = Vect                   !save input vector for Eigenvector calculation:
  Vout = 0d0
  alanc= 0d0
  blanc= 0d0
  nlanc= 0
  nhxv =0
  !
  lanc_loop: do iter=1,Nitermax
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,a_,b_)
     nhxv=nhxv+1
     if(abs(b_)<threshold_)exit lanc_loop
     !
     nlanc=nlanc+1
     !
     alanc(iter) = a_ ; blanc(iter+1) = b_
     !
     diag(1:Nlanc)    = alanc(1:Nlanc)
     subdiag(2:Nlanc) = blanc(2:Nlanc)
     call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
     !
     if(nlanc >= Ncheck_)then
        esave(nlanc-(Ncheck_-1))=diag(1)
        if(nlanc >= (Ncheck_+1))then
           diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
           if(verb.AND.mpi_master)write(*,*)'Iter, E0, deltaE = ',iter,diag(1),diff
           if(abs(diff).le.threshold_)exit lanc_loop
        endif
     endif
  enddo lanc_loop
  if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  !
  !============== END LANCZOS LOOP ======================
  !
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
  !
  !Get the Eigenvalues:
  egs = diag(1)
  !
  !Get the Eigenvector:
  Vin = Vect
  vout= 0d0
  vect= 0d0
  do iter=1,nlanc
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,alanc(iter),blanc(iter))
     nhxv=nhxv+1
     vect = vect + vin*Z(iter,1)
  end do
  norm_tmp=dot_product(vect,vect); norm=0d0
  call Allreduce_MPI(MpiComm,norm_tmp,norm)
  vect=vect/sqrt(norm)
  if(verb)then
     call MatVec(Nloc,vect,vout)
     if(mpi_master)write(*,*)"|H*v-E*v|=",sum(abs(vout-egs*vect))/Nloc
  endif
  Nitermax=Nlanc
  if(present(NumOp))NumOp=nhxv
end subroutine mpi_lanczos_eigh_d






!---------------------------------------------------------------------
!Purpose: use simple Lanczos to tri-diagonalize a matrix H (defined 
! in the H*v function).
!---------------------------------------------------------------------
subroutine mpi_lanczos_tridiag_d(MpiComm,MatVec,vin,alanc,blanc,threshold)
  integer                                      :: MpiComm
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  real(8),dimension(:),intent(inout)           :: vin !Nloc
  real(8),dimension(size(vin))                 :: vout,vtmp
  real(8),dimension(:),intent(inout)           :: alanc
  real(8),dimension(size(alanc)),intent(inout) :: blanc
  integer                                      :: Nitermax,Nloc,i
  integer                                      :: iter
  real(8)                                      :: a_,b_
  real(8),optional                             :: threshold
  !
  logical                                      :: mpi_master
  !
  if(MpiComm==Mpi_Comm_Null)return
  !
  mpi_master=get_master_MPI(MpiComm)
  !
  Nloc = size(vin)
  !
  if(present(threshold))threshold_=threshold
  !
  vtmp = vin
  Nitermax = size(alanc)
  a_=0d0
  b_=0d0
  vout=0d0
  do iter=1,Nitermax
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vtmp,vout,a_,b_)
     alanc(iter)=a_
     if(abs(b_)<threshold_)exit
     if(iter<nitermax)blanc(iter+1)=b_
  enddo
  if(iter==nitermax.AND.mpi_master)write(*,"(A)")"MPI_LANCZOS_TRIDIAG_D: reach Nitermax"
end subroutine mpi_lanczos_tridiag_d




!---------------------------------------------------------------------
!Purpose: plain homebrew lanczos iteration (no orthogonalization)
!note: the a,b variables are real, even in the complex matrix case
!to understand why check out the Gollub-Van Loan textbook.
!a it is easy: hermiticity->diag\in\RRR
!b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
!the absolute value. A lemma shows that the phase can be chosen 
!identically zero
!MPI VERSION
!---------------------------------------------------------------------
subroutine mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,alfa,beta)
  integer                                    :: MpiComm
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  real(8),dimension(:),intent(inout)         :: vin !Nloc
  real(8),dimension(size(vin)),intent(inout) :: vout
  real(8),dimension(size(vin))               :: tmp
  real(8),intent(inout)                      :: alfa,beta
  real(8)                                    :: atmp,btmp
  integer                                    :: iter,nloc
  real(8)                                    :: norm,norm_tmp
  !
  logical                                    :: mpi_master
  !
  if(MpiComm==Mpi_Comm_Null)return
  !
  nloc=size(vin)
  !
  mpi_master=get_master_MPI(MpiComm)
  !
  if(iter==1)then
     norm_tmp=dot_product(vin,vin)
     norm = 0d0 ; call AllReduce_MPI(MpiComm,norm_tmp,norm)
     if(mpi_master.AND.norm==0d0)stop "MPI_LANCZOS_ITERATION_D: norm = 0!!"
     vin=vin/sqrt(norm)
  else
     tmp = vin
     vin = vout/beta
     vout= -beta*tmp
  endif
  call MatVec(nloc,vin,tmp)
  vout = vout + tmp
  atmp = dot_product(vin,vout) ; alfa = 0d0; call AllReduce_MPI(MpiComm,atmp,alfa)
  vout = vout - alfa*vin
  btmp = dot_product(vout,vout); beta = 0d0; call AllReduce_MPI(MpiComm,btmp,beta)
  beta = sqrt(beta)
end subroutine mpi_lanczos_iteration_d
