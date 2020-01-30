!********************************************************************************************************************************!
!												                                 !
!		  	   Purpose: Calculate free energy differences and rate constants for elementary processes                !
!                           Author: Mohammad Shamsus Saleheen, Department of Chemical Engineering,USC                            !  
!	       	              Date: 06.18.16							                                 !
!   	              Modification:							                                         !
!          Reasons of Modification:                                                                                              !
!                            Notes: Please pay attention to the customization of this code, i.e. temperarure range, input        !
!                                   file names for reactant states, transition states, and product states, file names for        !
!                                   when gas phase species are involved, and the format of the input files (both for siurface    !
!                                   and gaseous species                                                                          ! 
!                                                               	      					                 !
!********************************************************************************************************************************!
program FEMain
use constants,                  only : kB, planck, light_cm, zero, one 
use FEGeneralProcedures
use FESpecificProcedures

implicit none
  integer                       :: iag, num_args
  character(sl), allocatable    :: args(:)
  integer                       :: nor, nots, nop
  integer                       :: nogr, nogp
  character(sl)                 :: FileR = 'reac-', FileP = 'prod-', FileTS = 'ts-'
  character(sl)                 :: FileGR = 'Greac-', FileGP = 'Gprod-'
  real(dp)                      :: Tmin  = 273.0, Tint = 25.0
  integer                       :: Tpoints = 29
  real(dp), allocatable         :: temp(:)
  integer                       :: j_temp
  real(dp)                      :: zpcRtot, zpcPtot, zpcTStot
  real(dp)                      :: zpcGtot
  real(dp)                      :: zpcTRtot, zpcTPtot
  real(dp), allocatable         :: QvibR(:), QvibP(:), QvibTS(:)
  real(dp), allocatable         :: Qgas(:)
  real(dp), allocatable         :: QgasR(:), QgasP(:)
  real(dp)                      :: mass
  real(dp)                      :: delE, delEF, delER                                                   ! zero point corrected reaction energy, forward, and reverse activation barrier
  integer                       :: alloc_err
 
  allocate (temp(Tpoints))
  allocate (QvibR(Tpoints))
  allocate (QvibP(Tpoints))
  allocate (QvibTS(Tpoints))
  allocate (Qgas(Tpoints))
  allocate (QgasR(Tpoints))
  allocate (QgasP(Tpoints))
  
  do j_temp = 1 , Tpoints
     temp(j_temp) = Tmin
     Tmin = Tmin + Tint
  end do
  
  num_args = command_argument_count()
  allocate (args(num_args))  
  do iag = 1, num_args
     call get_command_argument(iag, args(iag)) 
  end do  
    

  if       (args(1) == '-help')  then
            call  print_help()
  else if  (args(1) == '-surf')  then
        
        if (args(2) == '-thermo') then
               call countfiles(FileR, nor)
               call countfiles(FileP, nop)
               call fileforwrite_thermo()
               do j_temp = 1, Tpoints                               
                  call surface_species (FileR, nor, j_temp, temp, zpcRtot, QvibR)
                  call surface_species (FileP, nop, j_temp, temp, zpcPtot, QvibP)
                  call calculate_and_write_thermo(zpcRtot, zpcPtot, j_temp, temp, QvibR, QvibP, delE)
               end do
               write (40, 100) zpcRtot, zpcPtot, delE
               100 format (2X, 3 (F18.12,4X))
               close(20)
               close(30)
               close(40)

        else if  (args(2) == '-kin') then
               call countfiles(FileR, nor)
               call countfiles(FileP, nop)
               call countfiles(FileTS, nots)
               call fileforwrite_kinetics()
               do j_temp = 1, Tpoints
                  call surface_species (FileR, nor, j_temp, temp, zpcRtot, QvibR)
                  call surface_species (FileP, nop, j_temp, temp, zpcPtot, QvibP)
                  call surface_species (FileTS, nots, j_temp, temp, zpcTStot, QvibTS)
                  call calculate_and_write_kin(zpcRtot, zpcPtot, zpcTStot, j_temp, temp, QvibR, QvibP, QvibTS, delE, delEF, delER)
               end do
               write (40, 150) zpcRtot, zpcPtot, zpcTStot, delE, delEF, delER
               150 format (2X, 6 (F18.12,4X))
               close(20)
               close(30)
               close(40)
        else 
               write (*,*) 'Need a valid second argument: Check -help'
               stop    
        end if
  else if (args(1) == '-ads') then
               call countfiles(FileGR, nogr)
               call countfiles(FileR, nor)
               call countfiles(FileP, nop)
               call fileforwrite_adsdes()
               do j_temp = 1, Tpoints
                  call gas_species (FileGR, nogr, j_temp, temp, zpcGtot, Qgas, mass)
                  call surface_species (FileR, nor, j_temp, temp, zpcRtot, QvibR)
                  call surface_species (FileP, nop, j_temp, temp, zpcPtot, QvibP)
                  call calculate_and_write_ads(zpcGtot, zpcRtot, zpcTRtot, zpcPtot, mass, j_temp, temp, Qgas, QvibR, QvibP, delE)
               end do
               write (40, 200) zpcTRtot, zpcPtot, delE
               200 format (2X, 3 (F18.12,4X))
               close(20)
               close(30)
               close(40)
  else if (args(1) == '-des') then
               call countfiles(FileR, nor)
               call countfiles(FileP, nop)
               call countfiles(FileGP, nogp)
               call fileforwrite_adsdes()
               do j_temp = 1, Tpoints
                  call surface_species (FileR, nor, j_temp, temp, zpcRtot, QvibR)
                  call surface_species (FileP, nop, j_temp, temp, zpcPtot, QvibP)
                  call gas_species (FileGP, nogp, j_temp, temp, zpcGtot, Qgas, mass)
                  call calculate_and_write_des(zpcGtot, zpcRtot, zpcTPtot, zpcPtot, mass, j_temp, temp, Qgas, QvibR, QvibP, delE)
               end do
               write (40, 250) zpcRtot, zpcTPtot, delE
               250 format (2X, 3 (F18.12,4X))
               close(20)
               close(30)
               close(40)
  else if (args(1) == '-gas') then
               call countfiles(FileGR, nogr)
               call countfiles(FileGP, nogp)
               call fileforwrite_gas()
               do j_temp = 1, Tpoints
                  call gas_species (FileGR, nogr, j_temp, temp, zpcRtot, QgasR, mass)
                  call gas_species (FileGP, nogp, j_temp, temp, zpcPtot, QgasP, mass)
                  call calculate_and_write_gas(zpcRtot, zpcPtot, j_temp, temp, QgasR, QgasP, delE)
               end do
               write (40, 260) zpcRtot, zpcPtot, delE
               260 format (2X, 3 (F18.12,4X))
               close(20)
               close(30)
               close(40)
  else
       write (*,*) 'Answer not recognized, type -help for usage policy'
       stop  
  end if

  deallocate (temp,   stat = alloc_err)
  deallocate (QvibR,  stat = alloc_err)
  deallocate (QvibP,  stat = alloc_err)
  deallocate (QvibTS, stat = alloc_err)
  deallocate (Qgas,   stat = alloc_err)
  deallocate (QgasR,  stat = alloc_err)
  deallocate (QgasP,  stat = alloc_err)

end program FEMain
