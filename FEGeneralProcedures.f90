module FEGeneralProcedures
use constants,                only : dp, sl, zero, one, kB, planck, light_cm
implicit none
! Selective format specifiers
  character(sl), parameter :: F001 = '(A)'

contains


!*****************************************************************************************!
!         Procedure for counting lines (One less than total number of lines)              !
!*****************************************************************************************!

   Subroutine countline(filename,nless)
   implicit none
   character(sl), intent(in)   :: filename
   integer,       intent(out)  :: nless
   integer                     :: nlines, ierror

     nlines = 0
     OPEN(unit = 10, file = filename, status = 'old', action = 'read',&
          iostat = ierror)
          do
              read(10, *, iostat = ierror)
              if (ierror /= 0) exit
              nlines = nlines + 1
          end do
          nless = nlines - 1
     CLOSE(unit = 10)
   End Subroutine countline


!******************************************************************************************!
!                              Procedure for counting files                                !
!******************************************************************************************!

   Subroutine countfiles(filename, nfile)
   implicit none
   character(sl), intent(in)   :: filename
   integer,       intent(out)  :: nfile
   logical                     :: fexist
   character(sl)               :: Fname, string

     nfile = 1
          do 
              write(unit = string, fmt = '(I1)') nfile
              Fname = trim(filename) // trim(string)
              inquire (file = Fname, exist = fexist)
              if (.not. fexist) then
              nfile = nfile - 1
              exit
              end if
              nfile = nfile + 1
          end do
   End Subroutine countfiles

!*******************************************************************************************!
!       Procedure for reading frequencies and calculating total frequency and               !
!                            vibrational partition functions                                !
!*******************************************************************************************! 

   Subroutine deal_frequency_surface(filename, i, j_temp, s, temp, nless, scf, freqtot, qvib)
   implicit none

   integer,      intent(in)    :: i, j_temp
   integer,      intent(in)    :: s
   character(sl),intent(in)    :: filename
   integer,      intent(in)    :: nless
   real(dp),     intent(in)    :: temp(j_temp)
   real(dp),     intent(out)   :: scf
   real(dp),     intent(out)   :: freqtot
   real(dp),     intent(inout) :: qvib(s)
   real(dp),     allocatable   :: freq(:)
   integer                     :: ierror, alloc_err
   real(dp)                    :: vib
   real(dp)                    :: vibmulti

      freqtot  = zero
      vibmulti = one
      allocate (freq(nless))

      OPEN(unit = 10, file = filename, status = 'old', action = 'read',&
           iostat = ierror)
           read(10, *, iostat = ierror) scf
             do
                if (nless == 0) exit
                read(10, *, iostat = ierror) freq(nless)
                if (ierror /= 0) exit
                if (freq(nless) == 0) exit
                if (freq(nless) < 100) then
                 freq(nless) = 100
                end if
                freqtot   = freqtot + freq(nless)
                vib       = one / ( one - exp ( -planck * light_cm * freq(nless) / kB / temp(j_temp) ) )
                vibmulti  = vib * vibmulti
             end do
       qvib(i) = vibmulti

       CLOSE(10)
       deallocate (freq, stat = alloc_err)
   End Subroutine deal_frequency_surface


!***************************************************************************************************!
!       Procedure for opening and formatting files to write outputs of surface thermodynamics       !
!***************************************************************************************************!  

   Subroutine fileforwrite_thermo()
   implicit none
   integer            :: ierror

      OPEN (unit = 20, file = 'MKM_Partition_Func', access = 'sequential',&
            status = 'replace', action = 'write', position ='append',&
            iostat = ierror)
                write (20,900)
                900 format (3X,"Temp",5X,"Elec.Part.Function",5X,"Reac.Vib.Part.Fun.",&
                            5X,"Prod.Vib.Part.Fun.")
                write (20,910)
                910 format (2X,"======",4X,"===================",4X,"===================",&
                            4X,"===================")

      OPEN (unit = 30, file = 'MKM_Rate_Constants', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                 write (30,920)
                 920 format (3X,"Temp",4X,"Equilibrium Constant",4X,"Free Energy of Rxn.")
                 write (30,930)
                 930 format (2X,"=====",4X,"===================",4X,"===================")

      OPEN (unit = 40, file = 'MKM_Reaction_Energy', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                 write (40,940)
                 940 format (2X,"ZP Corr.Tot.Reactant",3X,"ZP Corr.Tot.Product",3X,&
                                 "Reaction Energy")
                 write (40,950)
                 950 format (2X,"====================",3X,"===================",3X,&
                                 "==================")
   End subroutine fileforwrite_thermo


!*********************************************************************************************!
!       Procedure for opening and formatting files to write outputs of surface kinetics       !
!*********************************************************************************************!       

   Subroutine fileforwrite_kinetics()
   implicit none
   integer              :: ierror
     
      OPEN (unit = 20, file = 'MKM_Partition_Func', access = 'sequential',&
            status = 'replace', action = 'write', position ='append',&
            iostat = ierror)
                 write (20,900)
                 900 format (3X,"Temp",5X,"Elec.Part.Function",5X,"Forward Elec.Part.",&
                             5X,"Reverse Elec.Part.",5X,"Reac.Vib.Part.Fun.",&
                             5X,"Prod.Vib.Part.Fun.",5X,"Trans.Vib.Part.Fun.")
                 write (20,910)
                 910 format (2X,"======",4X,"===================",4X,"===================",&
                             4X,"===================",4X,"===================",&
                             4X,"===================",4X,"===================")

      OPEN (unit = 30, file = 'MKM_Rate_Constants', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                  write (30,920)
                  920 format (3X,"Temp",4X,"Equilibrium Constant",4X,"Forward Rate.Cons.",&
                              4X,"Reverse Rate.Cons.",4X,"Free Energy of Rxn.",&
                              4X,"Forward Free Energy Barrier")
                  write (30,930)
                  930 format (2X,"======",4X,"===================",4X,"===================",&
                              4X,"===================",4X,"==================",&
                              4X,"===========================")

      OPEN (unit = 40, file = 'MKM_Reaction_Energy', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                   write (40,940)
                   940 format (2X,"ZP Corr.Tot.Reactant",3X,"ZP Corr.Tot.Product",3X,&
                               "ZP Corr.Tot.TransState",3X,"Reaction Energy",3X,&
                               "Forward Ac. Barrier",3X,"Reverse Ac. Barrier")
                   write (40,950)
                   950 format (2X,"====================",3X,"===================",3X,&
                               "======================",3X,"================",3X,&
                               "====================",3X,"==================")
   End subroutine fileforwrite_kinetics

!************************************************************************************************!
!       Procedure for opening and formatting files to write outputs of adsorption processes      !
!************************************************************************************************!       

   Subroutine fileforwrite_adsdes()
   implicit none
   integer              :: ierror

      OPEN (unit = 20, file = 'MKM_Partition_Func', access = 'sequential',&
            status = 'replace', action = 'write', position ='append',&
            iostat = ierror)
                 write (20,900)
                 900 format (3X,"Temp",5X,"Elec.Part.Function",5X,"Reac.Tot.Part.Fun.",&
                             5X,"Prod.Tot.Part.Fun.")
                 write (20,910)
                 910 format (2X,"======",4X,"===================",4X,"===================",&
                             4X,"===================")

      OPEN (unit = 30, file = 'MKM_Rate_Constants', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                  write (30,920)
                  920 format (3X,"Temp",4X,"Equilibrium Constant",4X,"Forward Rate.Cons.",&
                              4X,"Reverse Rate.Cons.",4X,"Free Energy of Rxn.",&
                              4X,"Forward Free Energy Barrier")
                  write (30,930)
                  930 format (2X,"======",4X,"===================",4X,"===================",&
                              4X,"===================",4X,"==================",&
                              4X,"===========================")

      OPEN (unit = 40, file = 'MKM_Reaction_Energy', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                   write (40,940)
                   940 format (2X,"ZP Corr.Tot.Reactant",3X,"ZP Corr.Tot.Product",3X,&
                               "Reaction Energy")
                   write (40,950)
                   950 format (2X,"====================",3X,"===================",3X,&
                               "================")
   End subroutine fileforwrite_adsdes

!*******************************************************************************************!
!   Procedure for opening and formatting files to write outputs of gaseous thermodynamics   !
!*******************************************************************************************!  

   Subroutine fileforwrite_gas()
   implicit none
   integer            :: ierror

      OPEN (unit = 20, file = 'MKM_Partition_Func', access = 'sequential',&
            status = 'replace', action = 'write', position ='append',&
            iostat = ierror)
                write (20,900)
                900 format (3X,"Temp",5X,"Elec.Part.Function",5X,"Reac.Tot.Part.Fun.",&
                            5X,"Prod.Tot.Part.Fun.")
                write (20,910)
                910 format (2X,"======",4X,"===================",4X,"===================",&
                            4X,"===================")

      OPEN (unit = 30, file = 'MKM_Rate_Constants', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                 write (30,920)
                 920 format (3X,"Temp",4X,"Equilibrium Constant",4X,"Free Energy of Rxn.")
                 write (30,930)
                 930 format (2X,"=====",4X,"===================",4X,"===================")

      OPEN (unit = 40, file = 'MKM_Reaction_Energy', access = 'sequential',&
            status = 'replace', action = 'write', position = 'append',&
            iostat = ierror)
                 write (40,940)
                 940 format (2X,"ZP Corr.Tot.Reactant",3X,"ZP Corr.Tot.Product",3X,&
                                 "Reaction Energy")
                 write (40,950)
                 950 format (2X,"====================",3X,"===================",3X,&
                                "==================")
   End subroutine fileforwrite_gas

!***********************************************************************!
!                       Procedure for help message                      !
!***********************************************************************!

   Subroutine print_help()
   implicit none
   character(sl)       :: msgs(1:100)
   integer              :: ii, nn

   ! Setup the help message
   nn = 0
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = '******************             USAGE Policy            *****************************'                
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = 'Free-energy-calculator -surf/ads/des/gas -thermo/kin(if the first argument is "surf"'
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = ' See InputFiles directory for input file format and README file for applicability   '
   nn = nn + 1; msgs(nn) = '                                                                                    '
   
   ! print help and terminate the program
   do ii = 1, nn
      write (*,F001) trim(msgs(ii))
   end do
   stop
   End Subroutine print_help

end module FEGeneralProcedures
