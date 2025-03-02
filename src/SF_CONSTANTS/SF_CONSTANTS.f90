module SF_CONSTANTS
!SciFortran module for physical and mathematical constants
  implicit none

  !COMMONLY USED PARAMETERS
  complex(8),parameter,public :: zero=(0.d0,0.d0)   !Complex number :math:`0+i0`
  complex(8),parameter,public :: xi=(0.d0,1.d0)     !Complex number :math:`0+i1`
  complex(8),parameter,public :: one=(1.d0,0.d0)    !Complex number :math:`1+i0`
  real(8),parameter,public    :: sqrt2 = 1.41421356237309504880169d0  !Real number :math:`\sqrt(2)`
  real(8),parameter,public    :: sqrt3 = 1.73205080756887729352745d0  !Real number :math:`\sqrt(3)`
  real(8),parameter,public    :: sqrt6 = 2.44948974278317809819728d0  !Real number :math:`\sqrt(6)`
  real(8),parameter,public    :: pi    = 3.14159265358979323846264338327950288419716939937510d0  !Real number :math:`\pi`
  real(8),parameter,public    :: pi2   = 6.28318530717959d0 !Real number :math:`2\pi`
  real(8),parameter,public    :: gamma_euler = 0.57721566490153286060d0  !Real number Euler-Mascheroni constant :math:`\gamma`
  real(8),parameter,public    :: euler= 2.7182818284590452353602874713526624977572470936999595749669676277240766303535d0 !Real Euler number :math:`e`
  integer,parameter,public    :: max_int  = huge(1) !Largest integer number
  real(8),parameter,public    :: max_real = huge(1.d0) !Largest real number
  real(8),parameter,public    :: epsilonr=epsilon(1.d0) !Smallest double precision real number
  real(8),parameter,public    :: epsilonq=1.d-30   ! Smallest quad precision real number
  integer,parameter,public    :: dbl=8             ! "double" precision
  integer,parameter,public    :: dp=8              ! "double" precision, alternative naming
  integer,parameter,public    :: ddp=16            ! "quad"   precision
  integer,parameter,public    :: sp = kind(1.0)    ! "single" precision



  !PHYSICAL CONSTANTS (expressed in the SI unit)
  real(8),parameter,public ::                                 Avogadro_constant=  0.602214129000D+24 !Number of particles per mole :math:`N_{A} [mol^{-1}]`
  real(8),parameter,public ::                                     Bohr_magneton=  0.927400968000D-23 !Natural unit of the electron magnetic moment 
                                                                                                     !:math:`\mu_{B} = \frac{e\hbar}{2m_{e}} [J\cdot T^{-1}]`
  real(8),parameter,public ::                             Bohr_magneton_in_eVoT=  0.578838180660D-04 !Natural unit of the electron magnetic moment  
                                                                                                     !:math:`\mu_{B} = \frac{e\hbar}{2m_{e}} [eV\cdot T^{-1}]`
  real(8),parameter,public ::                             Bohr_magneton_in_HzoT=  0.139962455500D+11 !Natural unit of the electron magnetic moment 
                                                                                                     !:math:`\mu_{B} = \frac{e\hbar}{2m_{e}} [Hz\cdot T^{-1}]`
  real(8),parameter,public ::         Bohr_magneton_in_inverse_meters_per_tesla=     46.6864498D0000 !Natural unit of the electron magnetic moment  
                                                                                                     !:math:`\mu_{B} = \frac{e\hbar}{2m_{e}} [m^{-1}\cdot T^{-1}]`
  real(8),parameter,public ::                              Bohr_magneton_in_KoT=      0.67171388D000 !Natural unit of the electron magnetic moment  
                                                                                                     !:math:`\mu_{B} = \frac{e\hbar}{2m_{e}} [K \cdot T^{-1}]`
  real(8),parameter,public ::                                       Bohr_radius=  0.529177210920D-10 !Approximate nucleus-electron distance in an Hydrogen atom 
                                                                                                     !:math:`a_{0} = \frac{4\pi\epsilon_{0}\hbar^{2}}{e^{2}m_{e}} [m]`
  real(8),parameter,public ::                                Boltzmann_constant=  0.138064880000D-22 !Proportionality constant between energy and temperature of an ideal gas 
                                                                                                     !:math:`k_{B} [J \cdot K^{-1}]`
  real(8),parameter,public ::                        Boltzmann_constant_in_eVoK=  0.861733240000D-04 !Proportionality constant between energy and temperature of an ideal gas 
                                                                                                     !:math:`k_{B} [eV \cdot K^{-1}]`
  real(8),parameter,public ::                        Boltzmann_constant_in_HzoK=  0.208366180000D+11 !Proportionality constant between energy and temperature of an ideal gas 
                                                                                                     !:math:`k_{B} [Hz \cdot K^{-1}]`
  real(8),parameter,public ::   Boltzmann_constant_in_inverse_meters_per_kelvin=     69.503476D00000 !Proportionality constant between energy and temperature of an ideal gas 
                                                                                                     !:math:`k_{B} [m^{-1} \cdot K^{-1}]`
  real(8),parameter,public ::                                Compton_wavelength=  0.242631023890D-11 !Wavelength of a photon having the same rest energy as the electron 
                                                                                                     !:math:`\lambda_{c} = \frac{h}{m_{e}c} [m]`
  real(8),parameter,public ::                      Compton_wavelength_over_2_pi=  0.386159268000D-12 !Normalized Compton wavelength :math:`\lambda_{c} = \frac{\hbar}{m_{e}c} [m]`
  real(8),parameter,public ::                                 electric_constant=  0.885418781700D-11 !Vacuum permittivity :math:`\epsilon_{0} [F\cdot m^{-1}]`
  real(8),parameter,public ::                  electron_charge_to_mass_quotient= -0.175882008800D+12 !Ratio between electronic charge and mass :math:`\frac{q_{e}}{m_{e}} [C \cdot kg^{-1}]`
  real(8),parameter,public ::                                 electron_g_factor= -0.200231930436D+01 !Electron g-factor :math:`g`  (Dimensionless)
  real(8),parameter,public ::                           electron_gyromag__ratio=  0.176085970800D+12 !Electron gyromagnetic ratio 
                                                                                                     !:math:`\gamma_{e} = \frac{g_{e}\mu_{B}}{\hbar}[rad\cdot s^{-1}\cdot T^{-1}]`
  real(8),parameter,public ::                 electron_gyromag__ratio_over_2_pi=  0.280249526600D+05 !Normalized electron gyromagnetic ratio :math:`\gamma_{e} = 
                                                                                                     !\frac{g_{e}\mu_{B}}{2\pi\hbar}[rad\cdot s^{-1}\cdot T^{-1}]`
  real(8),parameter,public ::                                electron_mag__mom_= -0.928476430000D-23 !Electron magnetic moment :math:`\mu_{e} [J\cdot T^{-1}]`
  real(8),parameter,public ::         electron_mag__mom__to_Bohr_magneton_ratio= -0.100115965218D+01 !Ratio between electron magnetic moment and Bohr magneton 
                                                                                                     !:math:`\frac{\mu_{e}}{\mu_{e}}=\frac{g}{2}`
  real(8),parameter,public ::                                     electron_mass=  0.910938291000D-30 !Mass of the electron :math:`m_{e}\;[kg]`
  real(8),parameter,public ::                   electron_mass_energy_equivalent=  0.818710506000D-13 !Electron rest energy :math:`[J]`
  real(8),parameter,public ::            electron_mass_energy_equivalent_in_MeV=      0.510998928D00 !Electron rest energy :math:`[MeV]`
  real(8),parameter,public ::                                     electron_volt=  0.160217656500D-18 !Equivalent of :math:`1\;eV\;[J]`
  real(8),parameter,public ::       electron_volt_atomic_mass_unit_relationship=  0.107354415000D-08 !Atomic mass unit equivalence of :math:`1\;eV`
  real(8),parameter,public ::                electron_volt_hartree_relationship=     0.03674932379D0 !Equivalent of :math:`1\;eV\;[Ha]`
  real(8),parameter,public ::                  electron_volt_hertz_relationship=  0.241798934800D+15 !Equivalent of :math:`1\;eV\;[Hz]`
  real(8),parameter,public ::          electron_volt_inverse_meter_relationship=  0.806554429000D+06 !Equivalent of :math:`1\;eV\;[m^-1]`
  real(8),parameter,public ::                  electron_volt_joule_relationship=  0.160217656500D-18 !Equivalent of :math:`1\;eV\;[J]`
  real(8),parameter,public ::                 electron_volt_kelvin_relationship=  0.116045190000D+05 !Equivalent of :math:`1\;eV\;[K]`
  real(8),parameter,public ::               electron_volt_kilogram_relationship=  0.178266184500D-35 !Equivalent of :math:`1\;eV\;[kg]`
  real(8),parameter,public ::                                 elementary_charge=  0.160217656500D-18 !Charge unit (charge of the proton) :math:`[C]`
  real(8),parameter,public ::                          elementary_charge_over_h=  0.241798934800D+15 !Charge unit (charge of the proton) divided by :math:`h\;[m^{-2} kg^{-1} s^{2} A]`
  real(8),parameter,public ::                                  Faraday_constant=  0.964853328900D+05 !Quotient between total charge and number of charge carriers (2019 SI units) 
                                                                                                     !:math:`\frac{q}{n}=eN_{A}\;[C\;mol^{-1}]`
  real(8),parameter,public ::Faraday_constant_for_conventional_electric_current=  0.964853251000D+05 !Quotient between total charge and 
                                                                                                     !number of charge carriers (1990 "conventional" values) 
                                                                                                     !:math:`\frac{q}{n}=eN_{A}\;[C\;mol^{-1}]`
  real(8),parameter,public ::                           fine_structure_constant=  0.729735256980D-02 !Fundamental constant characterizing the strength of 
                                                                                                     !the electromagnetic interaction between elementary 
                                                                                                     !charged particles :math:`\alpha=\frac{e^2}{2\epsilon_{0}hc}`
  real(8),parameter,public ::                                Josephson_constant=  0.483597870000D+15 !Inverse of the magnetic flux quantum :math:`K_{J}=\frac{2e}{h}\;[Hz\;V^{-1}]`
  real(8),parameter,public ::                  joule_electron_volt_relationship=  0.624150934000D+19 !Equivalent of :math:`1J` in :math:`[eV]`
  real(8),parameter,public ::                          joule_hertz_relationship=  0.150919031100D+34 !Equivalent of :math:`1J` in :math:`[Hz]`, from :math:`f=\frac{E}{h}`
  real(8),parameter,public ::                  joule_inverse_meter_relationship=  0.503411701000D+25 !Equivalent of :math:`1J` in :math:`[m^{-1}]`, from :math:`\frac{1}{\lambda}=\frac{E}{hc}`
  real(8),parameter,public ::                         joule_kelvin_relationship=  0.724297160000D+23 !Equivalent of :math:`1J` in :math:`[K]`, from :math:`T=\frac{E}{k_{B}}`
  real(8),parameter,public ::                       joule_kilogram_relationship=  0.111265005600D-16 !Equivalent of :math:`1J` in :math:`[kg]`, from :math:`m=\frac{E}{c^{2}}`
  real(8),parameter,public ::              kelvin_atomic_mass_unit_relationship=  0.925108680000D-13 !Equivalent of :math:`1K` in :math:`[amu]`
  real(8),parameter,public ::                 kelvin_electron_volt_relationship=  0.861733240000D-04 !Equivalent of :math:`1K` in :math:`[eV]`
  real(8),parameter,public ::                       kelvin_hartree_relationship=  0.316681140000D-05 !Equivalent of :math:`1K` in :math:`[Ha]`
  real(8),parameter,public ::                         kelvin_hertz_relationship=  0.208366180000D+11 !Equivalent of :math:`1K` in :math:`[Hz]`
  real(8),parameter,public ::                 kelvin_inverse_meter_relationship=     69.503476D00000 !Equivalent of :math:`1K` in :math:`[m^{-1}]`
  real(8),parameter,public ::                         kelvin_joule_relationship=  0.138064880000D-22 !Equivalent of :math:`1K` in :math:`[J]` (Boltzmann constant)
  real(8),parameter,public ::                      kelvin_kilogram_relationship=  0.153617900000D-39 !Equivalent of :math:`1K` in :math:`[kg]`
  real(8),parameter,public ::            kilogram_atomic_mass_unit_relationship=  0.602214129000D+27 !Equivalent of :math:`1kg` in :math:`[amu]`
  real(8),parameter,public ::               kilogram_electron_volt_relationship=  0.560958885000D+36 !Equivalent of :math:`1kg` in :math:`[eV]`
  real(8),parameter,public ::                     kilogram_hartree_relationship=  0.206148596800D+35 !Equivalent of :math:`1kg` in :math:`[Ha]`
  real(8),parameter,public ::                       kilogram_hertz_relationship=  0.135639260800D+50 !Equivalent of :math:`1kg` in :math:`[Hz]`
  real(8),parameter,public ::               kilogram_inverse_meter_relationship=  0.452443873000D+42 !Equivalent of :math:`1kg` in :math:`[m^{-1}]`
  real(8),parameter,public ::                       kilogram_joule_relationship=  0.898755178700D+17 !Equivalent of :math:`1kg` in :math:`[J]`
  real(8),parameter,public ::                      kilogram_kelvin_relationship=  0.650965820000D+40 !Equivalent of :math:`1kg` in :math:`[K]`
  real(8),parameter,public ::                      lattice_parameter_of_silicon=  0.543102050400D-09 !Lattice parameter of diamond cubic silicon lattice :math:`[m]`
  real(8),parameter,public ::                            natural_unit_of_action=  0.105457172600D-33 !Value of :math:`\hbar = \frac{h}{2\pi}` :math:`[J \cdot s]`
  real(8),parameter,public ::                    natural_unit_of_action_in_eV_s=  0.658211928000D-15 !Value of :math:`\hbar = \frac{h}{2\pi}` :math:`[eV \cdot s]`
  real(8),parameter,public ::                            natural_unit_of_energy=  0.818710506000D-13 !Rest energy of the electron :math:`[J]`
  real(8),parameter,public ::                     natural_unit_of_energy_in_MeV=      0.510998928D00 !Rest energy of the electron :math:`[MeV]`
  real(8),parameter,public ::                            natural_unit_of_length=  0.386159268000D-12 !Reduce Compton wavelength of the electron :math:`\lambda_{C}=\frac{\hbar}{m_{e}c}\;[m]`
  real(8),parameter,public ::                              natural_unit_of_mass=  0.910938291000D-30 !Mass of the electron :math:`m_{e}\;[kg]`
  real(8),parameter,public ::                            natural_unit_of_mom_um=  0.273092429000D-21 !Natural momentum unit :math:`m_{e}c\;[kg \cdot m \cdot s^{-1}]`
  real(8),parameter,public ::                   natural_unit_of_mom_um_in_MeVoc=      0.510998928D0  !Natural momentum unit in MeV over c :math:`m_{e}c\;[MeV\cdot c^{-1}]`
  real(8),parameter,public ::                              natural_unit_of_time=  0.128808866833D-20 !Reduced Compton time of the electron :math:`\t_{C}=\frac{\hbar}{m_{e}c^{2}}\;[s]`
  real(8),parameter,public ::                          natural_unit_of_velocity=  0.299792458000D+09 !Speed of light :math:`c\;[m\cdot s^{-1}]`
  real(8),parameter,public ::                 Newtonian_constant_of_gravitation=  0.667384000000D-10 !Gravitational constant :math:`G\;[m^{3}\cdot kg^{-1} \cdot s^{-2}]`
  real(8),parameter,public ::                                   Planck_constant=  0.662606957000D-33 !Planck constant :math:`h\;[J \cdot s]`
  real(8),parameter,public ::                           Planck_constant_in_eV_s=  0.413566751600D-14 !Planck constant :math:`h\;[eV \cdot s]`
  real(8),parameter,public ::                         Planck_constant_over_2_pi=  0.105457172600D-33 !Reduced Planck constant :math:`h\;[J \cdot s]`
  real(8),parameter,public ::                                  Rydberg_constant=  0.109737315685D+08 !Highest wavenumber of a photon emitted by a Hydrogen atom 
                                                                                                     !:math:`R_{\infty}=\frac{m_{e}e^{4}}{8\epsilon_{0}^{2}h^{3}c}\;[m^{-1}]`
  real(8),parameter,public ::                    Rydberg_constant_times_c_in_Hz=  0.328984196036D+16 !Rydberg constant times :math:`c` : 
                                                                                                     !:math:`R_{\infty}c=\frac{m_{e}e^{4}}{8\epsilon_{0}^{2}h^{3}}\;[Hz]`
  real(8),parameter,public ::                   Rydberg_constant_times_hc_in_eV=     13.60569253d000 !Rydberg constant times :math:`hc` : 
                                                                                                     !:math:`R_{\infty}c=\frac{m_{e}e^{4}}{8\epsilon_{0}^{2}h^{2}}\;[eV]`
  real(8),parameter,public ::                    Rydberg_constant_times_hc_in_J=  0.217987217100D-17 !Rydberg constant times :math:`hc` : 
                                                                                                     !:math:`R_{\infty}c=\frac{m_{e}e^{4}}{8\epsilon_{0}^{2}h^{2}}\;[J]`
  real(8),parameter,public ::                          speed_of_light_in_vacuum=  0.299792458000D+09 !Speed of light :math:`c\;[m\cdot s^{-1}]`
  real(8),parameter,public ::                  standard_acceleration_of_gravity=      9.80665D000000 !Acceleration of gravity on Earth :math:`g\;[m\cdot s^{-2}]`
  real(8),parameter,public ::                         Stefan_Boltzmann_constant=  0.567037300000D-07 !Total energy radiated per unit surface 
                                                                                                     !of a black body :math:`\sigma\;[W \cdot m^{2} \cdot K^{-4}]`



  public :: timestamp

  public :: stop_error
  
  interface isinfty
  !
  !Evaluates an integer, real or complex number and returns
  !a logical expression. If :code:`.true.` the input is infinite.
  !
     module procedure :: i_isinfty
     module procedure :: d_isinfty
     module procedure :: z_isinfty
  end interface isinfty
  public :: isinfty  

  interface isnan
  !
  !Evaluates an integer, real or complex quantity and returns
  !a logical expression. If :code:`.true.` the input is :code:`NaN`.
  !
     module procedure :: i_isnan
     module procedure :: d_isnan
     module procedure :: z_isnan
  end interface isnan
  public :: isnan


  interface wait
  !
  !Waits for a given time interval (integer, real or double) expressed
  !in milliseconds before continuing execution.
  !
     module procedure :: i_wait
     module procedure :: r_wait
     module procedure :: d_wait
  end interface wait
  public :: wait


contains


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: test if a given number if infinity
  !+-----------------------------------------------------------------------------+!
  elemental function i_isinfty(a) result(bool)
    integer,intent(in) :: a      !Number to evaluate 
    logical            :: bool   !Result
    bool = (a-1 == a)
  end function i_isinfty
  elemental function d_isinfty(a) result(bool)
    real(8),intent(in) :: a
    logical            :: bool
    bool = (a-1 == a)
  end function d_isinfty
  elemental function z_isinfty(a) result(bool)
    complex(8),intent(in) :: a
    logical               :: bool
    bool = (a-1 == a)
  end function z_isinfty


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: test if a given number is actually NaN
  !+-----------------------------------------------------------------------------+!
  elemental function i_isnan(a) result(bool)
    integer,intent(in)    :: a !Number to evaluate
    logical               :: bool !Result
    bool = (a /= a) .OR. (a-1 == a)
  end function i_isnan
  elemental function d_isnan(a) result(bool)
    real(8),intent(in)    :: a
    logical               :: bool
    bool = (a /= a) .OR. (a-1 == a)
  end function d_isnan
  elemental function z_isnan(a) result(bool)
    complex(8),intent(in) :: a
    logical               :: bool
    bool = (a /= a) .OR. (a-1 == a)
  end function z_isnan


  !+-------------------------------------------------------------------+
  !PURPOSE  : prints the current YMDHMS date as a time stamp.
  !+-------------------------------------------------------------------+
  subroutine timestamp(unit)
  !
  !Prints the date and time to :code:`unit`. 
  !The date is in the format :code:`Timestamp: [day][month name][year][h]:[m]:[s].[ms]`
  !
    integer,optional        :: unit  !Output unit. Default :code:`6`
    integer                 :: unit_
    integer(4),dimension(8) :: data
    unit_=6;if(present(unit))unit_=unit
    call date_and_time(values=data)
    call print_date(data,unit_)
  end subroutine timestamp



  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual date
  !+-------------------------------------------------------------------+
  subroutine print_date(dummy,unit)
    integer(4),dimension(8) :: dummy
    integer                 :: unit
    integer(4)                          :: year
    integer(4)                          :: mese
    integer(4)                          :: day
    integer(4)                          :: h
    integer(4)                          :: m
    integer(4)                          :: s
    integer(4)                          :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(unit,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(unit,*)""
  end subroutine print_date



  subroutine stop_error(msg)
    !
    !Aborts the program with nonzero exit code. 
    !The statement :code:`STOP "msg"` will return exit code :code:`0` when compiled using
    !gfortran. :f:func:`stop_error` uses the statement :code:`stop 1` which returns an exit code
    !:code:`1` and a print statement to print the message.
    !
    character(len=*) :: msg ! Message to print on stderr
    write(0,*) msg
    stop 1
  end subroutine stop_error


  subroutine i_wait(time)
    integer              :: time ! desired sleep interval [ms]
    integer,dimension(8) :: t    ! arguments for date_and_time
    integer              :: s1,s2,ms1,ms2 ! start and end times [ms]
    ! Get start time:
    call date_and_time(values=t)
    ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
    !
    do ! check time:
       call date_and_time(values=t)
       ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
       if(ms2-ms1>=time)exit
    enddo
    return
  end subroutine i_wait

  subroutine r_wait(time)
    real                 :: time ! desired sleep interval [ms]
    integer,dimension(8) :: t    ! arguments for date_and_time
    integer              :: s1,s2,ms1,ms2 ! start and end times [ms]
    ! Get start time:
    call date_and_time(values=t)
    ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
    !
    do ! check time:
       call date_and_time(values=t)
       ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
       if(ms2-ms1>=time)exit
    enddo
    return
  end subroutine r_wait

  subroutine d_wait(time)
    real(8)              :: time ! desired sleep interval [ms]
    integer,dimension(8) :: t    ! arguments for date_and_time
    integer              :: s1,s2,ms1,ms2 ! start and end times [ms]
    ! Get start time:
    call date_and_time(values=t)
    ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
    !
    do ! check time:
       call date_and_time(values=t)
       ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
       if(ms2-ms1>=time)exit
    enddo
    return
  end subroutine d_wait

END MODULE SF_CONSTANTS


