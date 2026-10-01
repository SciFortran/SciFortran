  !+-------------------------------------------------------------------+
  !PURPOSE  : Use non-linear least squares to fit a function, f, to data.
  !+-------------------------------------------------------------------+
  !LMDIF INTERFACE
  subroutine curvefit_lmdif_func(model_func,a,xdata,ydata,tol,info)
    interface
       function model_func(x,a)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: model_func
       end function model_func
    end interface
    real(8),dimension(:)           :: a      ! in: starting guess of the parameters; out: fitted parameters, size n
    real(8),dimension(:)           :: xdata  ! abscissas of the data, size m
    real(8),dimension(size(xdata)) :: ydata  ! ordinates of the data, size m
    integer                        :: m
    real(8),optional               :: tol    ! optional: tolerance, default 1d-15
    integer,optional               :: info   ! optional out: status of the MINPACK routine
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    !
    n=size(a)
    m=size(xdata)
    !
    call lmdif1(curvefit_lmdif_func2sub,m,n,a,fvec,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmdif_func2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      integer ::  iflag
      fvec = model_func(xdata,a) - ydata
      if(iflag<0)stop "CURVEFIT_LMDIF_func2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmdif_func2sub
  end subroutine curvefit_lmdif_func

  subroutine curvefit_lmdif_sub(model_func,a,xdata,ydata,tol,info)
    interface
       subroutine model_func(x,a,f)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: f
       end subroutine model_func
    end interface
    real(8),dimension(:)           :: a      ! in: starting guess of the parameters; out: fitted parameters, size n
    real(8),dimension(:)           :: xdata  ! abscissas of the data, size m
    real(8),dimension(size(xdata)) :: ydata  ! ordinates of the data, size m
    integer                        :: m
    real(8),optional               :: tol    ! optional: tolerance, default 1d-15
    integer,optional               :: info   ! optional out: status of the MINPACK routine
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    !
    n=size(a)
    m=size(xdata)
    !
    call lmdif1(curvefit_lmdif_sub2sub,m,n,a,fvec,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmdif_sub2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m),fvec_(m)
      integer ::  iflag
      call model_func(xdata,a,fvec_)
      fvec = fvec_ - ydata
      if(iflag<0)stop "CURVEFIT_LMDIF_sub2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmdif_sub2sub
  end subroutine curvefit_lmdif_sub


  !LMDER INTERFACE:
  subroutine curvefit_lmder_func(model_func,model_dfunc,a,xdata,ydata,tol,info)
    interface
       function model_func(x,a)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: model_func
       end function model_func
       !
       function model_dfunc(x,a)
         real(8),dimension(:)               :: x
         real(8),dimension(:)               :: a
         real(8),dimension(size(x),size(a)) :: model_dfunc
       end function model_dfunc
    end interface
    real(8),dimension(:)                   :: a      ! in: starting guess of the parameters; out: fitted parameters, size n
    real(8),dimension(:)                   :: xdata  ! abscissas of the data, size m
    real(8),dimension(size(xdata))         :: ydata  ! ordinates of the data, size m
    integer                                :: m
    real(8),optional                       :: tol    ! optional: tolerance, default 1d-15
    integer,optional                       :: info   ! optional out: status of the MINPACK routine
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    !
    call lmder1(curvefit_lmder1_func2sub,m,n,a,fvec,fjac,m,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmder1_func2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         fvec = model_func(xdata,a) - ydata
      elseif(iflag==2)then
         fjac = model_dfunc(xdata,a)
      endif
      if(iflag<0)stop "CURVEFIT_LMDER1_func2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmder1_func2sub
  end subroutine curvefit_lmder_func

  subroutine curvefit_lmder_sub(model_func,model_dfunc,a,xdata,ydata,tol,info)
    interface
       subroutine model_func(x,a,f)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: f
       end subroutine model_func
       !
       subroutine model_dfunc(x,a,df)
         real(8),dimension(:)               :: x
         real(8),dimension(:)               :: a
         real(8),dimension(size(x),size(a)) :: df
       end subroutine model_dfunc
    end interface
    real(8),dimension(:)                   :: a      ! in: starting guess of the parameters; out: fitted parameters, size n
    real(8),dimension(:)                   :: xdata  ! abscissas of the data, size m
    real(8),dimension(size(xdata))         :: ydata  ! ordinates of the data, size m
    integer                                :: m
    real(8),optional                       :: tol    ! optional: tolerance, default 1d-15
    integer,optional                       :: info   ! optional out: status of the MINPACK routine
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmder1(curvefit_lmder1_sub2sub,m,n,a,fvec,fjac,m,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine curvefit_lmder1_sub2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m),fvec_(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         call model_func(xdata,a,fvec_)
         fvec = fvec_ - ydata
      elseif(iflag==2)then
         call model_dfunc(xdata,a,fjac)
      endif
      if(iflag<0)stop "CURVEFIT_LMDER1_sub2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmder1_sub2sub
  end subroutine curvefit_lmder_sub
