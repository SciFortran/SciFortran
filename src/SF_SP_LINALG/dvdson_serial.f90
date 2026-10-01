! *   on entry                                                            
! *   -------                                                             
! *   OP          User supplied routine with calling sequence OP(N,M,B,C).
! *               B and C are N x M matrices and C stores the result AxB. 
! *               It should be declared external in the main program.     
! *   N           Order of the matrix.                                    
! *   LIM         The upper limit on the dimension of the expanding basis.
! *               NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
! *               only for LIM=NUME=N. The choice of LIM depends on the   
! *               available workspace (see below). If the space is        
! *               available it is preferable to have a large LIM, but not 
! *               larger than NUME$+$40.                                  
! *   DIAG        Array of size N with the diagonal elements of the       
! *               matrix A.                                               
! *   ILOW        The index of the lowest eigepair to be computed. If     
! *               (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs     
! *               to be computed should be contained in array ISELEC.     
! *               (Modified on exit).                                     
! *   IHIGH       The index of the highest eigenpair to be computed.      
! *               Considered ONLY when ILOW is in the range               
! *               (0.LT.ILOW.LE.N). (Modified on exit).                   
! *   ISELEC      Array of size LIM holding the user specified indices    
! *               for the eigenpairs to be computed. Considered only when 
! *               (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from   
! *               the first position until a non positive integer is met. 
! *                  Example: if N=500, ILOW=0, and ISELEC(1)=495,        
! *                  ISELEC(2)=497, ISELEC(3)=-1, the program will find   
! *                  2 of the highest eigenpairs, pairs 495 and 497.      
! *               Any order of indices is acceptable (Modified on exit).  
! *   NIV         Number of Initial Vector estimates provided by the user.
! *               If NIV is in the range:  (NUME).LE.(NIV).LE.(LIM),      
! *               the first NIV columns of size N of WORK should contain  
! *               the estimates (see below). In all other cases of NIV,   
! *               the program generates initial estimates.                
! *   MBLOCK      Number of vectors to be targeted in each iteration.     
! *               1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.     
! *               Large block size reduces the number of iterations       
! *               (matrix acceses) but increases the matrix-vector        
! *               multiplies. It should be used when the matrix accese    
! *               is expensive (disc, recomputed or distributed).         
! *   CRITE       Convergence threshold for eigenvalues.                  
! *               If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted 
! *               eigenvalues, convergence is signaled.                   
! *   CRITC       Convergence threshold for the coefficients of the last  
! *               added basis vector(s). If all of those corresponding to 
! *               unconverged eigenpairs are less than CRITC convergence  
! *               is signaled.                                            
! *   CRITR       Convergence threshold for residual vector norms. If     
! *               all the residual norms ||Ax_i-l_ix_i|| of the targeted  
! *               x_i are less than CRITR convergence is signaled.        
! *               If ANY of the criteria are satisfied the algorithm stops
! *   ORTHO       The threshold over which loss of orthogonality is       
! *               assumed. Usually ORTHO.LE.CRITR*10 but the process can  
! *               be skipped by setting ORTHO to a large number(eg,1.D+3).
! *   MAXITER     Upper bound on the number of iterations of the          
! *               algorithm. When MAXITER is exceeded the algorithm stops.
! *               A typical MAXITER can be MAX(200,NUME*40), but it can   
! *               be increased as needed.                                 
! *   WORK        Real array of size IWRSZ. Used for both input and output
! *               If NIV is in ((NUME).LE.(NIV).LE.(LIM)), on input, WORK 
! *               must have the NIV initial estimates. These NIV N-element
! *               vectors start from WORK(1) and continue one after the   
! *               other. They must form an orthonormal basis.             
! *   IWRSZ       The size of the real workspace. It must be at least as  
! *               large as:                                               
! *                                                                       
! *                       2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME        
! *                                                                       
! *   IWORK       Integer work array of size IIWSZ. Used as scrath array  
! *               for indices and for use in the LAPACK routines.         
! *   IIWSZ       The size of the integer workspace. It must be at least  
! *               as large as:                                            
! *                                    6*LIM + NUME                       
! *                                                                       
! *               If LIM or NUME needs to be increased, the space should  
! *               also be increased accordingly. For given IWRSZ and      
! *               IIWSZ one can calculate how big a problem one can       
! *               solve (LIM,NUME).                                       
! *                                                                       
! *   on exit                                                             
! *   -------                                                             
! *   WORK(1)     The first NUME*N locations contain the approximations to
! *               the NUME extreme eigenvectors. If the lowest eigenpairs 
! *               are required, (HIEND=false), eigenvectors appear in     
! *               ascending order, otherwise (HIEND=false), they appear in
! *               descending order. If only some are requested, the order 
! *               is the above one for all the NUME extreme eigenvectors, 
! *               but convergence has been reached only for the selected  
! *               ones. The rest are the current approximations to the    
! *               non-selected eigenvectors.                              
! *   WORK(NUME*N+1)                                                      
! *               The next NUME locations contain the approximations to   
! *               the NUME extreme eigenvalues, corresponding to the above
! *               NUME eigenvectors. The same ordering and convergence    
! *               status applies here as well.                            
! *   WORK(NUME*N+NUME+1)                                                 
! *               The next NUME locations contain the corresponding values
! *               of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of 
! *               the last step of the algorithm.                         
! *   WORK(NUME*N+NUME+NUME+1)                                            
! *               The next NUME locations contain the corresponding       
! *               residual norms of the NUME above eigenvectors, of the   
! *               last step.                                              
! *   HIEND       Logical. If .true. on exit the highest eigenpairs are   
! *               found in descending order. Otherwise, the lowest        
! *               eigenpairs are arranged in ascending order.             
! *   NLOOPS      The number of iterations it took to reach convergence.  
! *               This is also the number of matrix references.           
! *   NMV         The number of Matrix-vector(M-V) multiplies. Each matrix
! *               reference can have up to size(block) M-V multiplies.    
! *   IERR        An integer denoting the completions status:             
! *               IERR = 0        denotes normal completion.              
! *               IERR = -k       denotes error in DSPEVX (k eigenpairs   
! *                               not converged)

subroutine dvdson_eigh_d(MatVec,eval,evec,Nblock,Nitermax,Tol)
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(Nloc,vin,vout)            !vout=H*vin on vectors of size Nloc
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  !Arguments
  real(8),intent(inout)            :: eval(:)    !lowest eigenvalues in ascending order: [Neigen]
  real(8),intent(inout)            :: evec(:,:)  !in: 1st column = diagonal of H (computed if 0); out: eigenvectors [Ns,Neigen]
  integer,optional                 :: Nblock     !number of eigenpairs targeted in each iteration (default Neigen)
  integer,optional                 :: Nitermax   !maximum number of iterations (default max(200,40*Neigen))
  real(8),optional                 :: Tol        !convergence threshold on the eigenvalues (default 1d-15)
  !Dimensions:
  integer                          :: Ns
  integer                          :: Neigen
  integer                          :: Nume
  integer                          :: Lim
  integer                          :: Ilow,Ihigh
  integer                          :: Niv
  integer                          :: Ncv,Mblock
  integer                          :: MaxIter
  integer                          :: Iwrsz
  integer                          :: Iiwsz
  !Arrays:
  integer,dimension(:),allocatable :: Iselec
  real(8),dimension(:),allocatable :: Work
  integer,dimension(:),allocatable :: Iwork
  !Control Vars:
  integer                          :: i,j,ie,ival,ivec
  real(8)                          :: Crite
  real(8)                          :: Critc
  real(8)                          :: Critr
  real(8)                          :: Ortho
  logical                          :: Hiend
  integer                          :: Nloops
  integer                          :: Nmv
  integer                          :: Ierr
  !
  real(8),dimension(:),allocatable :: Adiag

  Ncv = Neigen ; if(present(Nblock))Ncv = Nblock

  Ns     = size(evec,1)
  Neigen = size(eval)
  if(any(shape(Evec)/=[Ns,Neigen]))stop "dvdson_eigh_d ERROR: Evec wrong shape"
  !
  !
  Nume   = Neigen !Distance of the eigenpair index ik farthest from the extreme
  !               !we look for the first, i.e. lowest, Neigen Eigenpairs here.
  Lim    = min(Nume + 40,Ns)
  Ilow   = 1
  Ihigh  = Neigen
  Niv    = 0                    !to be used for initial vectors
  Mblock = max(1,min(Ncv,Neigen))
  Crite  = 1d-15    ; if(present(tol))Crite=tol
  Critc  = 1d-15
  Critr  = 1d-15
  Ortho  = Critr
  MaxIter = max(200,Nume*40) ; if(present(Nitermax))MaxIter = Nitermax
  Iwrsz  = 2*Ns*Lim + Lim*Lim + (Nume+10)*Lim + Nume
  Iiwsz  = 6*Lim + Nume
  !
  allocate(Iselec(lim))
  allocate(Work(Iwrsz))
  allocate(Iwork(Iiwsz))
  Iselec = 0
  Work   = 0d0
  Iwork  = 0

  allocate(Adiag(Ns))
  if(sum(abs(Evec(:,1)))==0d0)then
     call MatVec_to_diagonal(Adiag)
  else
     Adiag = Evec(:,1)
  endif
  !
  call DVDSON(block_MatVec, Ns, Lim, Adiag,&
       Ilow, Ihigh, Iselec, NIV, Mblock,&
       Crite, Critc, Critr, ortho, Maxiter,&
       Work,Iwrsz,Iwork,Iiwsz,HIEND,NLOOPS,NMV,IERR)

  !POST PROCESSING:
  if(ierr/=0)then
     if(ierr<0)then
        write(*,*)"IERR = ",ierr
        stop "sp_dvdson_eigh: IERR denotes error in DSPEVX (k eigenpairs not converged)"
     else
        write(*,*)"IERR = ",ierr
        stop "sp_dvdson_eigh: IERR denotes some error, see documentation"
        ! if (INT( MOD(IERR,  2)/1) == 0 )stop "sp_dvdson_eigh: IERR signals N < LIM"
        ! If (INT( MOD(IERR,  4)/2) == 0 ) stop "sp_dvdson_eigh: IERR signals LIM < 1"
        ! If (INT( MOD(IERR,  8)/4) == 0  ) stop "sp_dvdson_eigh: IERR signals ISELEC(1)<1, and no range specified"
        ! If (INT( MOD(IERR, 16)/8) == 0  ) stop "sp_dvdson_eigh: IERR signals IHIGH > N (in range or ISELEC)"
        ! If (INT( MOD(IERR, 32)/16) == 0 ) stop "sp_dvdson_eigh: IERR signals IHIGH < ILOW (Invalid range)"       
        ! If (INT( MOD(IERR, 64)/32) == 0 ) stop "sp_dvdson_eigh: IERR signals NEIG >= LIM (Too many wanted)"      
        ! If (INT( MOD(IERR,128)/64) == 0 ) stop "sp_dvdson_eigh: IERR signals Probable duplication in ISELEC"     
        ! If (INT( MOD(IERR,256)/128) == 0) stop "sp_dvdson_eigh: IERR signals NUME >= LIM (max eigen very far)"   
        ! If (INT( MOD(IERR,512)/256) == 0) stop "sp_dvdson_eigh: IERR signals MBLOCK is out of bounds"   
        ! If (INT( MOD(IERR,1024)/512) == 0) stop "sp_dvdson_eigh: IERR signals IWRSZ or IIWSZ is not enough"
        ! If (INT( MOD(IERR,2048)/1024) == 0) stop "sp_dvdson_eigh: IERR signals Orthogonalization Failed"     
        ! If (INT( MOD(IERR,4096)/2048) == 0) stop "sp_dvdson_eigh: IERR signals NLOOPS > MAXITER"        
     endif
  endif

  !Retrieve the Evals,Evecs results from Work array
  do ie=1,Nume
     ivec = 1  + (ie-1)*Ns
     ival = ie + Nume*Ns 
     Eval(ie) = Work(ival)
     Evec(:,ie) = Work(ivec:ivec+Ns-1)
  enddo


contains

  subroutine block_MatVec(N,Ne,B,C)
    integer :: N,Ne
    real(8) :: B(N,Ne)
    real(8) :: C(N,Ne)
    integer :: ie
    !
    do ie=1,Ne
       call MatVec(N,B(:,ie),C(:,ie))
    enddo
    !
  end subroutine block_MatVec

  subroutine MatVec_to_diagonal(diag)
    real(8) :: diag(Ns)
    integer :: i
    real(8) :: vin(Ns),vout(Ns)
    do i=1,Ns
       vin    = 0d0
       vin(i) = 1d0
       call MatVec(Ns,Vin,Vout)
       diag(i)= Vout(i)
    enddo
  end subroutine MatVec_to_diagonal


end subroutine DVDSON_EIGH_D
