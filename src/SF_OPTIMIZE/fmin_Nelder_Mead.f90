
!+-------------------------------------------------------------------+
!PURPOSE  : 
! NELMIN minimizes a function using the Nelder-Mead algorithm.
!    This routine seeks the minimum value of a user-specified function.
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
! !    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
! !    is estimated to minimize the function.
!
!  OPTIONAL:
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!+-------------------------------------------------------------------+
subroutine fmin(fn,start,&
     lambda,tol,conv_check,max_fun_calls,fun_calls,num_restart,ierr)
  !This subroutine minimizes a function of several variables with the Nelder-Mead downhill simplex algorithm (:code:`NELMIN`),
  !a simplex function minimisation procedure due to Nelder and Mead (1965), as implemented by O'Neill (Applied Statistics 20,
  !338, 1971). The convergence check, carried out every :f:var:`conv_check` iterations, is on the variance of the function
  !values at the vertices of the simplex; if the simplex is found not to be at a minimum it is restarted. The minimum value of
  !the function is not returned.
  !
  !@param fn: function to be minimized, real(8) function fn(x) with real(8) x(:)
  !@param start: in: starting point; out: coordinates of the minimum
  !@param lambda: optional: size and shape of the initial simplex, one step per variable, default 1
  !@param tol: optional: terminating limit for the variance of the function values, > 0, default 1d-8
  !@param conv_check: optional: the convergence is checked every conv_check iterations, >= 1, default 10
  !@param max_fun_calls: optional: maximum number of function evaluations, default 500
  !@param fun_calls: optional out: number of function evaluations used
  !@param num_restart: optional out: number of restarts
  !@param ierr: optional out: 0 ok, 1 illegal tol or conv_check, 2 max_fun_calls exceeded
  !
  interface
     function fn(x)
       real(8),dimension(:) :: x
       real(8)              :: fn
     end function fn
  end interface
  real(8)            :: start(:)
  real(8),optional   :: lambda(size(start)) !--> step
  real(8),optional   :: tol                 !--> reqmin 
  integer,optional   :: conv_check          !--> konvge
  integer,optional   :: max_fun_calls       !--> kcount
  integer,optional   :: fun_calls           !--> icount
  integer,optional   :: num_restart         !--> numres
  integer,optional   :: ierr                !--> ifault
  !
  real(8)            :: step(size(start))
  real(8)            :: reqmin
  integer            :: konvge   
  integer            :: kcount
  integer            :: icount
  integer            :: numres    
  integer            :: ifault
  !
  real(8)            :: xmin(size(start))
  integer            :: n
  real(8), parameter :: ccoeff = 0.5D+00
  real(8)            :: del
  real(8), parameter :: ecoeff = 2.0D+00
  real(8), parameter :: eps = 0.001D+00
  integer            :: i
  integer            :: ihi
  integer            :: ilo
  integer            :: j
  integer            :: jcount
  integer            :: l
  real(8)            :: p(size(start),size(start)+1)
  real(8)            :: p2star(size(start))
  real(8)            :: pbar(size(start))
  real(8)            :: pstar(size(start))
  real(8),parameter  :: rcoeff = 1.0D+00
  real(8)            :: rq
  real(8)            :: x
  real(8)            :: y(size(start)+1)
  real(8)            :: y2star
  real(8)            :: ylo
  real(8)            :: ynewlo
  real(8)            :: ystar
  real(8)            :: z
  !
  n = size(start)
  !
  step=1d0;if(present(lambda))step=lambda
  reqmin=1d-8;if(present(tol))reqmin=tol
  konvge=10;if(present(conv_check))konvge=conv_check
  kcount=500;if(present(max_fun_calls))kcount=max_fun_calls    
  !
  !  Check the input parameters.
  !
  if ( reqmin <= 0.0D+00 ) then
     ifault = 1
     return
  end if
  if ( n < 1 ) then
     ifault = 1
     return
  end if
  if ( konvge < 1 ) then
     ifault = 1
     return
  end if
  !
  !  Initialization.
  !
  icount = 0
  numres = 0
  jcount = konvge
  del = 1.0D+00
  rq = reqmin * real ( n, kind = 8 )
  !
  !  Initial or restarted loop.
  !
  do
     p(1:n,n+1) = start(1:n)
     y(n+1) = fn ( start )
     icount = icount + 1
     !
     !  Define the initial simplex.
     !
     do j = 1, n
        x = start(j)
        start(j) = start(j) + step(j) * del
        p(1:n,j) = start(1:n)
        y(j) = fn ( start )
        icount = icount + 1
        start(j) = x
     end do
     !
     !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
     !  the vertex of the simplex to be replaced.
     !
     ilo = minloc ( y(1:n+1), 1 )
     ylo = y(ilo)
     !
     !  Inner loop.
     !
     do while ( icount < kcount )
        !
        !  YNEWLO is, of course, the HIGHEST value???
        !
        ihi = maxloc ( y(1:n+1), 1 )
        ynewlo = y(ihi)
        !
        !  Calculate PBAR, the centroid of the simplex vertices
        !  excepting the vertex with Y value YNEWLO.
        !
        do i = 1, n
           pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
        end do
        !
        !  Reflection through the centroid.
        !
        pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
        ystar = fn ( pstar )
        icount = icount + 1
        !
        !  Successful reflection, so extension.
        !
        if ( ystar < ylo ) then
           p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
           y2star = fn ( p2star )
           icount = icount + 1
           !
           !  Retain extension or contraction.
           !
           if ( ystar < y2star ) then
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
           else
              p(1:n,ihi) = p2star(1:n)
              y(ihi) = y2star
           end if
           !
           !  No extension.
           !
        else
           l = 0
           do i = 1, n + 1
              if ( ystar < y(i) ) then
                 l = l + 1
              end if
           end do
           if ( 1 < l ) then
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
              !
              !  Contraction on the Y(IHI) side of the centroid.
              !
           else if ( l == 0 ) then
              p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
              !
              !  Contract the whole simplex.
              !
              if ( y(ihi) < y2star ) then
                 do j = 1, n + 1
                    p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                    xmin(1:n) = p(1:n,j)
                    y(j) = fn ( xmin )
                    icount = icount + 1
                 end do
                 ilo = minloc ( y(1:n+1), 1 )
                 ylo = y(ilo)
                 cycle
                 !
                 !  Retain contraction.
                 !
              else
                 p(1:n,ihi) = p2star(1:n)
                 y(ihi) = y2star
              end if
              !
              !  Contraction on the reflection side of the centroid.
              !
           else if ( l == 1 ) then
              p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
              !
              !  Retain reflection?
              !
              if ( y2star <= ystar ) then
                 p(1:n,ihi) = p2star(1:n)
                 y(ihi) = y2star
              else
                 p(1:n,ihi) = pstar(1:n)
                 y(ihi) = ystar
              end if
           end if
        end if
        !
        !  Check if YLO improved.
        !
        if ( y(ihi) < ylo ) then
           ylo = y(ihi)
           ilo = ihi
        end if
        jcount = jcount - 1
        if ( 0 < jcount ) then
           cycle
        end if
        !
        !  Check to see if minimum reached.
        !
        if ( icount <= kcount ) then
           jcount = konvge
           x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
           z = sum ( ( y(1:n+1) - x )**2 )
           if ( z <= rq ) then
              exit
           end if
        end if
     end do
     !
     !  Factorial tests to check that YNEWLO is a local minimum.
     !
     xmin(1:n) = p(1:n,ilo)
     ynewlo = y(ilo)
     if ( kcount < icount ) then
        ifault = 2
        exit
     end if
     ifault = 0
     do i = 1, n
        del = step(i) * eps
        xmin(i) = xmin(i) + del
        z = fn ( xmin )
        icount = icount + 1
        if ( z < ynewlo ) then
           ifault = 2
           exit
        end if
        xmin(i) = xmin(i) - del - del
        z = fn ( xmin )
        icount = icount + 1
        if ( z < ynewlo ) then
           ifault = 2
           exit
        end if
        xmin(i) = xmin(i) + del
     end do
     if ( ifault == 0 ) then
        exit
     end if
     !
     !  Restart the procedure.
     !
     start(1:n) = xmin(1:n)
     del = eps
     numres = numres + 1
  end do
  if(present(fun_calls))fun_calls=icount
  if(present(num_restart))num_restart=numres
  if(present(ierr))ierr=ifault
  start = xmin
  return
end subroutine fmin
