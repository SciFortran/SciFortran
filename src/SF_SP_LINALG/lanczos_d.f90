!---------------------------------------------------------------------
!Purpose: use plain lanczos to get the groundstate energy
!---------------------------------------------------------------------
subroutine lanczos_eigh_d(MatVec,Egs,Vect,Nitermax,iverbose,threshold,ncheck,vrandom,NumOp)
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  real(8)                              :: egs
  real(8),dimension(:)                 :: vect
  integer                              :: Nitermax
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
  real(8)                              :: a_,b_
  real(8)                              :: norm,diff
  integer                              :: i,ierr
  logical                              :: vran=.true.
  !
  if(present(iverbose))verb=iverbose
  if(present(threshold))threshold_=threshold
  if(present(ncheck))ncheck_=ncheck
  if(present(vrandom))vran=vrandom
  !
  norm=dot_product(vect,vect)
  if(norm==0d0)then
     if(vran)then
        call mt_random(vect)
        if(verb)write(*,*)"LANCZOS_EIGH_D: random initial vector generated:"
     else
        vect = 1d0
        if(verb)write(*,*)"LANCZOS_EIGH_D: unitary initial vector generated:"
     endif
     vect=vect/sqrt(dot_product(vect,vect))
  endif
  !
  !============= LANCZOS LOOP =====================
  !
  vin = vect
  vout= 0.d0
  alanc=0.d0
  blanc=0.d0
  nlanc=0
  nhxv =0
  !
  lanc_loop: do iter=1,Nitermax
     call lanczos_iteration_d(MatVec,iter,vin,vout,a_,b_)
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
           if(verb)write(*,*)'iter, E0, deltaE = ',iter,diag(1),diff
           if(abs(diff).le.threshold_)exit lanc_loop
        endif
     endif
     !
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
  vin =vect
  vout=0d0
  vect=0d0
  do iter=1,Nlanc
     call lanczos_iteration_d(MatVec,iter,vin,vout,alanc(iter),blanc(iter))
     nhxv=nhxv+1
     vect = vect + vin*Z(iter,1)
  end do
  norm=sqrt(dot_product(vect,vect))
  vect=vect/norm
  !
  if(verb)then
     call MatVec(size(vect),vect,vout)
     write(*,*)"|H*v-E*v|=",sum(abs(vout-egs*vect))/size(vect)
  endif
  Nitermax=Nlanc
  if(present(NumOp))NumOp=nhxv
end subroutine lanczos_eigh_d













!---------------------------------------------------------------------
!Purpose: use simple Lanczos to tri-diagonalize a matrix H (defined 
! in the H*v function).
!---------------------------------------------------------------------
subroutine lanczos_tridiag_d(MatVec,vin,alanc,blanc,threshold)
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  real(8),dimension(:),intent(inout)           :: vin
  real(8),dimension(size(vin))                 :: vout
  real(8),dimension(:),intent(inout)           :: alanc
  real(8),dimension(size(alanc)),intent(inout) :: blanc
  integer                                      :: Nitermax
  integer                                      :: iter
  real(8)                                      :: a_,b_
  real(8),optional                             :: threshold
  !
  if(present(threshold))threshold_=threshold
  !
  Nitermax = size(alanc)
  a_=0d0
  b_=0d0
  vout=0d0
  do iter=1,Nitermax
     call lanczos_iteration_d(MatVec,iter,vin,vout,a_,b_)
     alanc(iter)=a_
     if(abs(b_)<threshold_)exit
     if(iter<nitermax)blanc(iter+1)=b_
  enddo
  if(iter==nitermax)write(*,"(A)")"LANCZOS_TRIDIAG_D: reach Nitermax"
end subroutine lanczos_tridiag_d








!---------------------------------------------------------------------
!Purpose: plain homebrew lanczos iteration (no orthogonalization)
!note: the a,b variables are real, even in the complex matrix case
!to understand why check out the Gollub-Van Loan textbook.
!a it is easy: hermiticity->diag\in\RRR
!b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
!the absolute value. A lemma shows that the phase can be chosen 
!identically zero
!---------------------------------------------------------------------
subroutine lanczos_iteration_d(MatVec,iter,vin,vout,alfa,beta)
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  integer                                    :: iter
  real(8),dimension(:),intent(inout)         :: vin
  real(8),dimension(size(vin)),intent(inout) :: vout
  real(8),dimension(size(vin))               :: tmp
  real(8),intent(inout)                      :: alfa,beta
  integer                                    :: nloc
  real(8)                                    :: norm
  !
  nloc=size(vin)
  !
  if(iter==1)then
     norm=sqrt(dot_product(vin,vin))
     if(norm==0.d0)stop "LANCZOS_ITERATION_D: norm =0!!"
     vin    = vin/norm
  else
     tmp = vin
     vin = vout/beta
     vout= -beta*tmp
  endif
  call MatVec(nloc,vin,tmp)
  vout = vout + tmp
  alfa = dot_product(vin,vout)
  vout = vout - alfa*vin
  beta = sqrt(dot_product(vout,vout))
end subroutine lanczos_iteration_d
