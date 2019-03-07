module FESpecificProcedures
use constants,                   only: dp, sl, zero, one, Pi, kB, planck, light_cm, Na, kB_SI
use FEGeneralProcedures
implicit none

! Selective parameters
  real(dp),     parameter :: BartoPa = 100000                 
  real(dp),     parameter :: SPA = 1.5816077992244E+19

contains


!*******************************************************************************************!
!                      Procedure for dealing with surface species                           !
!*******************************************************************************************!

   Subroutine surface_species (filename, s, j_temp, temp, zpctot, Qvib_temp)
   implicit none
   character(sl), intent(in)  :: filename
   integer,       intent(in)  :: s
   integer,       intent(in)  :: j_temp
   real(dp),      intent(in)  :: temp(j_temp)
   real(dp),      intent(out) :: zpctot
   real(dp),      intent(out) :: Qvib_temp(j_temp)
   character(sl)              :: Fname, string
   integer                    :: i
   real(dp),      allocatable :: zpe(:), zpc(:), qvib(:)
   integer                    :: alloc_err
   integer                    :: nless
   real(dp)                   :: freqtot, scf
   
     zpctot = zero     
     allocate (zpe(s))
     allocate (zpc(s))
     allocate (qvib(s))

       do i = 1, s
          write (unit = string, fmt = '(I1)') i
          Fname = trim (filename) // trim (string)
          call countline (Fname, nless)
          call deal_frequency_surface (Fname, i, j_temp, s, temp, nless, scf, freqtot, qvib)
          zpe(i) = 0.50 * planck * light_cm * freqtot
          zpc(i) = scf + zpe(i)
          zpctot = zpctot + zpc(i)
       end do
     Qvib_temp(j_temp) = PRODUCT(qvib)
     deallocate (zpe,  stat = alloc_err)
     deallocate (zpc,  stat = alloc_err)
     deallocate (qvib, stat = alloc_err)
           
   End Subroutine surface_species

!*******************************************************************************************!
!                      Procedure for dealing with gaseous species                           !
!*******************************************************************************************!

   Subroutine gas_species (filename, s, j_temp, temp, zpctot, Qgas, mass)
   implicit none
   character(sl), intent(in)  :: filename
   integer,       intent(in)  :: s
   integer,       intent(in)  :: j_temp
   real(dp),      intent(in)  :: temp(j_temp)
   real(dp),      intent(out) :: zpctot
   real(dp),      intent(out) :: Qgas(j_temp)
   real(dp),      intent(out) :: mass
   character(sl)              :: Fname, string
   integer                    :: i
   real(dp),      allocatable :: zpe(:), zpc(:), qtot(:)
   integer                    :: ierror, alloc_err
   real(dp)                   :: scf
   real(dp)                   :: filler, diff_filler

   zpctot = zero
   allocate (zpe(s))
   allocate (zpc(s))
   allocate (qtot(s))

     do i = 1, s   
        write (unit = string, fmt = '(I1)') i
        Fname = trim (filename) // trim (string)
        OPEN(unit = 10, file = Fname, status = 'old', action = 'read',&
             iostat = ierror)
             read(10, *, iostat = ierror) scf
             read(10, *, iostat = ierror) zpe(i)
             read(10, *, iostat = ierror) mass
               do 
                  read(10, *, iostat = ierror) filler, qtot(i)
                  if (ierror /= 0) exit
                      diff_filler = filler - temp(j_temp)
                      if ( abs(diff_filler) .gt. 0.001 ) cycle
                      if ( abs(diff_filler) .lt. 0.001 ) exit
               end do 
             zpc(i) = scf + zpe(i)
             zpctot = zpctot + zpc(i)
      end do   
    Qgas(j_temp) = PRODUCT(qtot)
    CLOSE(10)
    deallocate (zpe,  stat = alloc_err)
    deallocate (zpc,  stat = alloc_err)
    deallocate (qtot, stat = alloc_err)

   End Subroutine

!*******************************************************************************************!
!             Procedure for calculating and writing surface thermodynamics                  !
!*******************************************************************************************! 

   Subroutine calculate_and_write_thermo(zpcrtot, zpcptot, j_temp, temp, qvibr, qvibp, delE)
   implicit none
   real(dp),    intent(in)   :: zpcrtot, zpcptot
   integer,     intent(in)   :: j_temp
   real(dp),    intent(in)   :: temp(j_temp)
   real(dp),    intent(in)   :: qvibr(j_temp), qvibp(j_temp)
   real(dp),    intent(out)  :: delE
   real(dp)                  :: delG
   real(dp)                  :: qelec, Keq

        delE     = zpcptot - zpcrtot
        qelec    = exp ( - delE / kB / temp(j_temp) )
        Keq      = qelec * qvibp(j_temp) / qvibr(j_temp)
        delG     = -kB * temp(j_temp) * log (Keq)
        
        write(20, 600) temp(j_temp), qelec, qvibr(j_temp), qvibp(j_temp)
        600 format (2X, F6.2, 3X, 3 (ES20.13,3X))
        write(30,610) temp(j_temp), Keq, delG
        610 format (2X, F6.2, 3X, 2 (ES20.13,3X))         
              
   End Subroutine calculate_and_write_thermo

!*******************************************************************************************!
!             Procedure for calculating and writing surface kinetics                        !
!*******************************************************************************************! 

   Subroutine calculate_and_write_kin(zpcrtot, zpcptot, zpctstot, j_temp, temp, qvibr, qvibp, qvibts, delE, delEF, delER)
   implicit none
   real(dp),    intent(in)   :: zpcrtot, zpcptot, zpctstot
   integer,     intent(in)   :: j_temp
   real(dp),    intent(in)   :: temp(j_temp)
   real(dp),    intent(in)   :: qvibr(j_temp), qvibp(j_temp), qvibts(j_temp)
   real(dp),    intent(out)  :: delE, delEF, delER
   real(dp)                  :: delG, delGF
   real(dp)                  :: qelec, qelecF, qelecR, Keq, kf, kr

        delE     =  zpcptot - zpcrtot
        delEF    =  zpctstot - zpcrtot
        delER    =  zpctstot - zpcptot
        qelec    =  exp ( - delE / kB / temp(j_temp) )
        qelecF   =  exp ( - delEF / kB / temp(j_temp) )
        qelecR   =  exp ( - delER / kB / temp(j_temp) )
        Keq      =  qelec * qvibp(j_temp) / qvibr(j_temp)
        kf       =  kB * temp(j_temp) * qelecF * qvibts(j_temp) / qvibr(j_temp) / planck
        kr       =  kB * temp(j_temp) * qelecR * qvibts(j_temp) / qvibp(j_temp) / planck
        delG     = -kB * temp(j_temp) * log (Keq)
        delGF    = -kB * temp(j_temp) * log ( kf * planck / ( kB * temp(j_temp) ) )

        write(20, 600) temp(j_temp), qelec, qelecF, qelecR, qvibr(j_temp), qvibp(j_temp), qvibts(j_temp)
        600 format (2X, F6.2, 3X, 6 (ES20.13,3X))
        write(30,610) temp(j_temp), Keq, kf, kr, delG, delGF
        610 format (2X, F6.2, 3X, 5 (ES20.13,3X))

   End Subroutine calculate_and_write_kin

!*******************************************************************************************!
!             Procedure for calculating and writing adsorption processes                    !
!*******************************************************************************************! 

   Subroutine calculate_and_write_ads(zpcgtot, zpcrtot, zpcTRtot, zpcptot, mass, j_temp, temp, Qgas, qvibr, qvibp, delE)
   implicit none

   real(dp),    intent(in)   :: zpcgtot, zpcrtot, zpcptot
   real(dp),    intent(out)  :: zpcTRtot
   real(dp),    intent(in)   :: mass
   integer,     intent(in)   :: j_temp
   real(dp),    intent(in)   :: temp(j_temp)
   real(dp),    intent(in)   :: Qgas(j_temp), qvibr(j_temp), qvibp(j_temp)
   real(dp),    intent(out)  :: delE
   real(dp)                  :: qelec, QRtot
   real(dp)                  :: Keq, kf, kr
   real(dp)                  :: delG, delGF

        zpcTRtot =  zpcgtot + zpcrtot
        delE     =  zpcptot - zpcTRtot
        qelec    =  exp ( - delE / kB / temp(j_temp) )
        QRtot    =  Qgas(j_temp) * qvibr(j_temp)
        Keq      =  qelec * ( qvibp(j_temp) / QRtot )
        kf       =  ( BartoPa * sqrt(Na) ) / ( SPA * sqrt ( 2.0 * Pi * mass * kB_SI * temp(j_temp) ) )
        kr       =  kf / Keq
        delG     = -kB * temp(j_temp) * log (Keq)
        delGF    = -kB * temp(j_temp) * log ( ( kf * planck ) / ( kB * temp(j_temp) ) ) 

        write(20, 600) temp(j_temp), qelec, QRtot, qvibp(j_temp)
        600 format (2X, F6.2, 3X, 3 (ES20.13,3X))
        write(30,610) temp(j_temp), Keq, kf, kr, delG, delGF
        610 format (2X, F6.2, 3X, 5 (ES20.13,3X))

   End Subroutine calculate_and_write_ads



!*******************************************************************************************!
!             Procedure for calculating and writing desorption processes                    !
!*******************************************************************************************! 

   Subroutine calculate_and_write_des(zpcgtot, zpcrtot, zpcTPtot, zpcptot, mass, j_temp, temp, Qgas, qvibr, qvibp, delE)
   implicit none

   real(dp),    intent(in)   :: zpcgtot, zpcrtot, zpcptot
   real(dp),    intent(out)  :: zpcTPtot
   real(dp),    intent(in)   :: mass
   integer,     intent(in)   :: j_temp
   real(dp),    intent(in)   :: temp(j_temp)
   real(dp),    intent(in)   :: Qgas(j_temp), qvibr(j_temp), qvibp(j_temp)
   real(dp),    intent(out)  :: delE
   real(dp)                  :: qelec, QPtot
   real(dp)                  :: Keq, kf, kr
   real(dp)                  :: delG, delGF

        zpcTPtot =  zpcgtot + zpcptot
        delE     =  zpcTPtot - zpcrtot
        qelec    =  exp ( - delE / kB / temp(j_temp) )
        QPtot    =  Qgas(j_temp) * qvibp(j_temp)
        Keq      =  qelec * QPtot / qvibr(j_temp)
        kr       =  ( BartoPa * sqrt(Na) ) / ( SPA * sqrt ( 2.0 * Pi * mass * kB_SI * temp(j_temp) ) )
        kf       =  kr * Keq
        delG     = -kB * temp(j_temp) * log (Keq)
        delGF    = -kB * temp(j_temp) * log ( kf * planck / ( kB * temp(j_temp) ) )

        write(20, 600) temp(j_temp), qelec, qvibr(j_temp), QPtot
        600 format (2X, F6.2, 3X, 3 (ES20.13,3X))
        write(30,610) temp(j_temp), Keq, kf, kr, delG, delGF
        610 format (2X, F6.2, 3X, 5 (ES20.13,3X))

   End Subroutine calculate_and_write_des

!*******************************************************************************************!
!             Procedure for calculating and writing gaseous processes                       !
!*******************************************************************************************! 

   Subroutine calculate_and_write_gas(zpcrtot, zpcptot, j_temp, temp, qgasr, qgasp, delE)
   implicit none

   real(dp),    intent(in)   :: zpcrtot, zpcptot
   integer,     intent(in)   :: j_temp
   real(dp),    intent(in)   :: temp(j_temp)
   real(dp),    intent(in)   :: qgasr(j_temp), qgasp(j_temp)
   real(dp),    intent(out)  :: delE
   real(dp)                  :: qelec
   real(dp)                  :: Keq
   real(dp)                  :: delG

        delE     =  zpcptot - zpcrtot
        qelec    =  exp ( - delE / kB / temp(j_temp) )
        Keq      =  qelec * qgasp(j_temp) / qgasr(j_temp)
        delG     = -kB * temp(j_temp) * log (Keq)

        write(20, 600) temp(j_temp), qelec, qgasr(j_temp), qgasp(j_temp)
        600 format (2X, F6.2, 3X, 3 (ES20.13,3X))
        write(30,610) temp(j_temp), Keq, delG
        610 format (2X, F6.2, 3X, 2 (ES20.13,3X))

   End Subroutine calculate_and_write_gas

end module FESpecificProcedures
