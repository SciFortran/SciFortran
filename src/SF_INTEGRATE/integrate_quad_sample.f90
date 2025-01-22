subroutine quad_sample(fsample,a,b,&
     epsabs,epsrel,&
     key,&
     singular_endpoint,&
     singular_points,&
     cpole,&
     alfa,beta,&
     omega,&
     weight_func,&
     verbose,&
     Ninterp,&
     result)
  real(8),dimension(:)             :: fsample
  real(8)                          :: a
  real(8)                          :: b
  !optional variables
  real(8),optional                 :: epsabs
  real(8),optional                 :: epsrel
  integer,optional                 :: key               !order of GK rule in QAG
  logical,optional                 :: singular_endpoint !T if singular_endpoint exists (QAGS)
  real(8),dimension(:),optional    :: singular_points   !location of singular points in QAGP
  real(8),optional                 :: cpole             !location of the pole in QAWC f(x)/x-cpole
  real(8),optional                 :: alfa,beta         !parameters for QAWS
  real(8),optional                 :: omega             !frequency of the weight funcion in QAWO
  integer,optional                 :: weight_func       !if(QAWF)then
  !                                                  !  weight_func=1,2 -> cos(omega*x),sin(omega*x)
  !                                                  !if(QAWS)then
  !                                                  !  weight_func = 1  (x-a)**alfa*(b-x)**beta
  !                                                  !  weight_func = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
  !                                                  !  weight_func = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
  !                                                  !  weight_func = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
  logical,optional                 :: verbose
  integer,optional                 :: Ninterp
  real(8)                          :: result
  !actual default variables
  real(8)                          :: epsabs_
  real(8)                          :: epsrel_
  logical                          :: verbose_
  real(8)                          :: abserr
  real(8),dimension(size(fsample)) :: xsample
  integer                          :: Neval
  integer                          :: Ier,i
  integer                          :: Nsingular_points
  integer                          :: Ninterp_,Lin
  character(len=20)                :: routine_name
  !
  type finter_type
     real(8),allocatable :: X(:)
     real(8),allocatable :: F(:)
     integer             :: Imin=0,Imax=0,N=0
     logical             :: status=.false.
  end type finter_type
  type(finter_type)                :: Finterp
  !
  !
  epsabs_ = 1d-12 ; if(present(epsabs))epsabs_=epsabs
  epsrel_ = 1d-6  ; if(present(epsrel))epsrel_=epsrel
  verbose_=.false.; if(present(verbose))verbose_=verbose
  Ninterp_=3      ; if(present(Ninterp))Ninterp_=Ninterp
  !
  !Build the interpolating function to be integrated:
  Lin = size(fsample)
  xsample = sf_integrate_linspace(a,b,Lin)
  call init_finter(Finterp,xsample,fsample,Ninterp_)
  !
  if(present(key).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: key & alfa,beta"
  if(present(key).AND.present(cpole))stop "ERROR in quad: key & cpole"
  if(present(key).AND.present(singular_points))stop "ERROR in quad: key & singular_points"
  if(present(key).AND.present(singular_endpoint))stop "ERROR in quad: key & singular_endpoint"
  if(present(key).AND.present(omega))stop "ERROR in quad: key & omega"
  if(present(singular_endpoint).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_endpoint & alfa,beta"
  if(present(singular_endpoint).AND.present(cpole))stop "ERROR in quad: singular_endpoint & cpole"
  if(present(singular_endpoint).AND.present(singular_points))stop "ERROR in quad: singular_endpoint & singular_points"
  if(present(singular_endpoint).AND.present(omega))stop "ERROR in quad: singular_endpoint & omega"
  if(present(cpole).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: cpole & alfa,beta"
  if(present(cpole).AND.present(singular_points))stop "ERROR in quad: cpole & singular_points"
  if(present(cpole).AND.present(omega))stop "ERROR in quad: cpole & omega"
  if(present(omega).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: omega & alfa,beta"
  if(present(omega).AND.present(singular_points))stop "ERROR in quad: omega & singular_points"
  if(present(singular_points).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_points & alfa,beta"
  !
  !QNG:  no need
  !QAGS: need singular_endpoint=T
  !QAG : need key=1,6
  !QAGP: need singular_points
  !QAWC: need cpole
  !QAWO: need omega + weight_func=1,2
  !QAWS: need alfa + beta + weight_func=1,2,3,4
  routine_name='QNG'
  if(present(singular_endpoint))routine_name='QAGS'
  if(present(key))routine_name='QAG'
  if(present(singular_points))routine_name='QAGP'
  if(present(cpole))routine_name='QAWC'
  if(present(omega))routine_name='QAWO'
  if(present(alfa).AND.present(beta))routine_name='QAWS'
  !
  if(verbose_)write(*,"(A,A)")"QUAD: selected procedure =", routine_name
  !
  select case(routine_name)
  case ('QNG')
     call QNG(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAGS')
     call QAGS(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'(Singular) endpoints (a*,b*)                =', a,b
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAG')
     if(key<1.OR.key>6)stop "ERROR in quad: use QAG, key !=[1:6]"
     call QAG(func,a,b,epsabs_,epsrel_,key,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,I8)')    'Approximation or                            =', key
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAGP')
     Nsingular_points=size(singular_points)
     call QAGP(func,a,b,Nsingular_points,singular_points,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,I8)')'N_singular_points                           =', Nsingular_points
        do i=1,Nsingular_points
           write(*,"(A,I6,F14.6)")'I_singular_point                            =',i,singular_points(i)
        enddo
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAWC')
     call QAWC(func,a,b,cpole,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,F14.6)')'C_pole                                      =', cpole
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAWO')
     if(weight_func<1.OR.weight_func>2)stop "ERROR in quad: use QAWO, weight_func !=[1,2]"
     if(verbose_)then
        if(weight_func==1)then
           write(*,'(A)')'W(x) function is cos(omega*x)'
        elseif(weight_func==2)then
           write(*,'(A)')'W(x) function is sin(omega*x)'
        endif
     endif
     call QAWO(func,a,b,omega,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then 
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,F14.6)') 'Omega                                       =', omega
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  case ('QAWS')
     if(weight_func<1.OR.weight_func>4)stop "ERROR in quad: use QAWS, weight_func !=[1,..,4]"
     if(verbose_)then
        if(weight_func==1)then
           write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta         '
        elseif(weight_func==2)then
           write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)'
        elseif(weight_func==3)then
           write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(b-x) '
        elseif(weight_func==4)then
           write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)'
        endif
     endif
     call QAWS(func,a,b,alfa,beta,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,2F14.6)')'Parameters (alfa,beta)                      =', alfa,beta
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier/=0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           stop
        endif
     endif
     !
     !
  end select
  !
  call delete_finter(finterp)
  !
  !
contains
  !
  !
  subroutine init_finter(func,Xin,Fin,N)
    type(finter_type) :: func
    real(8)           :: xin(:)
    real(8)           :: fin(size(xin))
    integer           :: N,Lin
    if(func%status)deallocate(func%x,func%f)
    Lin=size(xin)
    allocate(func%x(Lin),func%f(Lin))
    func%X    = Xin
    func%F    = Fin
    func%Imax = Lin
    func%Imin = 1
    func%N    = N
    func%status=.true.
  end subroutine init_finter
  !
  subroutine delete_finter(func)
    type(finter_type) :: func
    if(allocated(func%x))deallocate(func%x)
    if(allocated(func%f))deallocate(func%f)
    func%imax=0
    func%imin=0
    func%N   =0
    func%status=.false.
  end subroutine delete_finter
  !
  function func(x)
    real(8)           :: x
    real(8)           :: func
    real(8)           :: dy
    integer           :: j,k,k0,k1
    integer           :: n
    N=finterp%N    !order of polynomial interpolation
    func=0d0
    j=locate(finterp%X(finterp%Imin:finterp%Imax),x)
    !
    k=max(j-(N-1)/2,1)
    k0=k
    if(k0 < finterp%Imin)k0=finterp%Imin
    k1=k+N+1
    if(k1 > finterp%Imax)then
       k1=finterp%Imax
       k0=k1-N-1
    endif
    call polint(finterp%X(k0:k1),finterp%F(k0:k1),x,func,dy)
  end function func
end subroutine quad_sample
