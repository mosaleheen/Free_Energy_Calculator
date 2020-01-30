module constants
  implicit none

  public
  save


    ! Intrinsic datatypes
    integer, parameter  :: dp = selected_real_kind (15, 307)
    integer, parameter  :: ip = selected_real_kind (15)
    integer, parameter  :: sl = 1000
    
    ! Number constants
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one  = 1.0_dp 

    ! Physical constants and conversion factors
    ! Values obtained from Peter J. Mohr et al.,Rev. Mod. Phys., 2016    

    real(dp), parameter :: Pi           =  3.1415926536_dp
    real(dp), parameter :: kB           =  8.6173303e-5_dp                                             	! boltzmann constant in eV/K
    real(dp), parameter :: planck       =  4.135667662e-15_dp                                         	! planck constant in eV.s
    real(dp), parameter :: light_m      =  2.99792458e8_dp                                             	! velocity of light in m/s
    real(dp), parameter :: Na    	=  6.022140857e23_dp						! avogradro's number in (1/mol) 
    real(dp), parameter :: B2A		=  0.52917721067_dp						! 1 Bohr     = X Angstrom
    real(dp), parameter :: A2B		=  one / B2A							! 1 Angstrom = X Bohr
    real(dp), parameter :: Ha2eV	=  27.21138602_dp						! 1 Hartree  = X eV
    real(dp), parameter :: eV2Ha	=  one / Ha2eV							! 1 eV       = X Hartree
    real(dp), parameter :: kBha		=  kB * eV2Ha  							! boltzmann constant in Ha/K	
    real(dp), parameter :: planckha     =  planck * eV2ha  						! planck constant Ha.s
    real(dp), parameter :: HaB2eVA	=  B2A * eV2Ha							! Hartree/Bohr to eV/Angstrom conversion
    real(dp), parameter :: eVA2HaB	=  one / HaB2eVA						! eV/Angstrom to Hartree/Bohr conversion
    real(dp), parameter :: J2eV		=  6.241509126e18_dp						! 1 Joules   = X eV  
    real(dp), parameter :: eV2J         =  one / J2eV							! 1 eV       = X Joules
    real(dp), parameter :: light_cm     =  2.99792458e10_dp                                             ! velocity of light in cm/s
    real(dp), parameter :: planck_SI    =  planck * eV2J                                                ! planck constant in J.s
    real(dp), parameter :: kB_SI        =  kB * eV2J                                                    ! boltzmann constant in J/K

end module constants    
  
