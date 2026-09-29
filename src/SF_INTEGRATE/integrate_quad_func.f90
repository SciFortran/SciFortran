subroutine quad_func(func,a,b,epsabs,epsrel,&
     key,&
     inf,&
     singular_endpoint,&
     singular_points,&
     cpole,&
     alfa,beta,&
     omega,&
     weight_func,&
     verbose,&
     strict,&
     result)
  interface
     function func(x)   !integrand f(x)
       real(8) :: x
       real(8) :: func
     end function func
  end interface
  !optional variables
  real(8)                       :: a                 !lower integration limit (finite limit of the range if b is absent)
  real(8),optional              :: b                 !upper integration limit (absent: semi-infinite range, needs inf or omega)
  real(8),optional              :: epsabs            !absolute accuracy requested (default 1d-12)
  real(8),optional              :: epsrel            !relative accuracy requested (default 1d-6)
  integer,optional              :: key               !Gauss-Kronrod rule of QAG: 1..6 = 15,21,31,41,51,61 points
  integer,optional              :: inf               !range of QAGI: -1 (-inf:a), 1 (a:inf), 2 (-inf:inf)
  logical,optional              :: singular_endpoint !T if singular_endpoint exists (QAGS)
  real(8),dimension(:),optional :: singular_points   !location of singular points in QAGP
  real(8),optional              :: cpole             !location of the pole in QAWC f(x)/x-cpole
  real(8),optional              :: alfa,beta         !exponents of the QAWS weight (x-a)**alfa*(b-x)**beta, both > -1
  real(8),optional              :: omega             !frequency of the weight function cos/sin(omega*x) in QAWO,QAWF
  integer,optional              :: weight_func       !weight function index
  logical,optional              :: verbose           !if T print the selected routine and the integration report (default F)
  logical,optional              :: strict            !T (default) stop on QUADPACK error code IER (only if verbose=F)
  real(8)                       :: result
  !actual default variables
  real(8)                       :: epsabs_
  real(8)                       :: epsrel_
  logical                       :: verbose_
  logical                       :: strict_
  real(8)                       :: abserr
  integer                       :: Neval
  integer                       :: Ier,i
  integer                       :: Nsingular_points
  character(len=20)             :: routine_name
  !
  !
  !
  epsabs_ = 1d-12 ; if(present(epsabs))epsabs_=epsabs
  epsrel_ = 1d-6  ; if(present(epsrel))epsrel_=epsrel
  verbose_=.false.; if(present(verbose))verbose_=verbose
  strict_ =.true. ; if(present(strict))strict_=strict
  !
  if(present(key).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: key & alfa,beta"
  if(present(key).AND.present(cpole))stop "ERROR in quad: key & cpole"
  if(present(key).AND.present(singular_points))stop "ERROR in quad: key & singular_points"
  if(present(key).AND.present(singular_endpoint))stop "ERROR in quad: key & singular_endpoint"
  if(present(key).AND.present(omega))stop "ERROR in quad: key & omega"
  if(present(key).AND.present(inf))stop "ERROR in quad: key & inf"
  if(present(singular_endpoint).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_endpoint & alfa,beta"
  if(present(singular_endpoint).AND.present(cpole))stop "ERROR in quad: singular_endpoint & cpole"
  if(present(singular_endpoint).AND.present(singular_points))stop "ERROR in quad: singular_endpoint & singular_points"
  if(present(singular_endpoint).AND.present(omega))stop "ERROR in quad: singular_endpoint & omega"
  if(present(singular_endpoint).AND.present(inf))stop "ERROR in quad: singular_endpoint & inf"
  if(present(cpole).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: cpole & alfa,beta"
  if(present(cpole).AND.present(singular_points))stop "ERROR in quad: cpole & singular_points"
  if(present(cpole).AND.present(omega))stop "ERROR in quad: cpole & omega"
  if(present(cpole).AND.present(inf))stop "ERROR in quad: cpole & inf"
  if(present(omega).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: omega & alfa,beta"
  if(present(omega).AND.present(singular_points))stop "ERROR in quad: omega & singular_points"
  if(present(omega).AND.present(inf))stop "ERROR in quad: omega & inf"
  if(present(singular_points).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_points & alfa,beta"
  if(present(singular_points).AND.present(inf))stop "ERROR in quad: singular_points & inf"
  if(present(inf).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: inf & alfa,beta"
  !
  if(present(b))then
     !QNG>QAGS:  no need
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
  else
     !QAGI: need inf
     !QAWF: need omega + weight_func=1,2
     if(present(inf))then
        routine_name='QAGI'
     elseif(present(omega))then
        routine_name='QAWF'
     else
        stop 'ERROR in quad: b not present but neither inf (QAGI) nor omega (QAWF) are given. stop.'
     endif
  endif
  !
  if(verbose_)write(*,"(A,A)")"QUAD: selected procedure =", routine_name
  !
  result=0d0
  !
  select case(routine_name)
  case ('QNG')
     ! call QNG(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
     ! call QAG(func,a,b,epsabs_,epsrel_,1,result,abserr,neval,ier)
     call QAGS(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier>0)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a
        write(*,'(A,F14.6)') 'Omega                                       =', omega
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
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
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
        endif
     endif
     !
     !
  case ('QAGI')
     if(inf<-1.OR.inf>2)stop "ERROR in quad: use QAGI, inf !=[-1,1,2]"
     call QAGI(func,a,inf,epsabs_,epsrel_,result,abserr,neval,ier)
     if(verbose_)then
        if(inf==-1)then
           write(*,'(A,F14.6)')'Integration interval is (-infty:a)',a
        elseif(inf==1)then
           write(*,'(A,F14.6)')'Integration interval is (a:infty) ',a
        elseif(inf==2)then
           write(*,'(A)')'Integration interval is (-infty:infty)'
        endif
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
        endif
     endif
     !
     !
  case ('QAWF')
     if(weight_func<1.OR.weight_func>2)stop "ERROR in quad: use QAWF, weight_func !=[1,2]"
     if(verbose_)then
        if(weight_func==1)then
           write(*,'(A)')'W(x) function is cos(omega*x)'
        elseif(weight_func==2)then
           write(*,'(A)')'W(x) function is sin(omega*x)'
        endif
     endif
     call QAWF(func,a,omega,weight_func,epsabs_,result,abserr,neval,ier)
     if(verbose_)then
        write(*,'(A,2F14.6)')'Endpoints (a,infty)                         =', a
        write(*,'(A,F14.6)') 'Omega                                       =', omega
        write(*,'(A,F14.6)') 'Estimated integral                          =', result
        write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
        write(*,'(A,I8)')    'Error return code IER                       =', ier
     else
        if(ier>2)then
           write(*,'(A,I8)') 'Error return code IER =', ier
           if(strict_)stop
        endif
     endif
     !
     !
  end select
end subroutine quad_func
