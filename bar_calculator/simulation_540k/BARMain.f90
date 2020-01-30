!********************************************************************************************************************************!      
!												                                 !
!		  	   Purpose: Use Bennett Acceptance Ratio to calculate free energy difference between two states		 !
!                           Author: Mohammad Shamsus Saleheen, Department of Chemical Engineering,USC                            !  
!	       	              Date: 06.18.17							                                 !
!   	              Modification:							                                         !
!          Reasons of Modification:                                                                                              !
!                           Notes: False position method was used to solve the non-linear equation, the implicit function        !
!                                   adopted here is defined in Shirts, M. R. et. al.; "An introduction to best practices         !
!                                   in free energy calculations", equation 4                                                     !
!                                                               	      					                 !
!********************************************************************************************************************************!


program BARMain
use constants
use BARSubroutines
implicit none

        !**********************************************************************************************************************************
        !      						Member Elements                                                                   !
        !                                           **********************                                                                !
        !       'Ffile/RFile'   --- forward and reverse workfunction input files                                                          !       
        !       'nF/nR'	        --- number of forward and reverse workfunction values (not thermal energy weighted), the input files      !
        !                           contain 2 header lines, hence, nF,nR = nline -2. 							  !				
        !       'uF/uR'	        --- forward and reverse workfunctions (1D array), units in hartree                           		  !
        !       'M'             --- log ratio number of forward and reverse workfunctions                                                 !
        !       'zF/zR'	        --- exponential averaging results using forward and reverse workfunctions                                 !   		
        !       'ub/lb'         --- upper and lower bound for deltaF                                                                      !
        !       'Fu/Fl'         --- bar function value using upper and lower bounds                                                       !
        !       'deltaF'        --- free energy difference using BAR method                                                       	  !	
        !       'Fsum/Rsum'	--- sum of the forward and reverse workfunction exponential terms, as defined in Shirts, M. R.	          !
        !		            et. al.; "An introduction to best practices in free energy calculations", equation 4                  !
        !       'boundstat'     --- boundstat returns a value depending on where the bound of the solution lies, that's how it            !             
        !                           writes down whether the solution lies within first, second, or third bound in the output file         !
        !                                                                                                                                 !
        !					       Member Procedures                                                                  !
        !                                           **********************                                                                !
        !       'countline'     --- procedure to count number of lines in a file                                                          !
        !       'readarray'     --- reads the values of forward and reverse work functions                                                !
        !       'zwanzig'       --- procedure to employ zwanzig equation                                                                  !
        !       'BARequation'   --- procedure to solve the BAR implicit equation                                                          ! 
        !       'BARSolver'     --- procedure to solve BAR equation using false position method                                           !
        !                                                                                                                                 !
        !**********************************************************************************************************************************



! Input and output files

    character(20) :: Ffile = 'DeltaFF'
    character(20) :: Rfile = 'DeltaFR'

! Declaration of variables

    integer		  :: nline, alloc_err, ierror                   
    integer		  :: nF, nR
    real(dp), allocatable :: uF(:)
    real(dp), allocatable :: uR(:)
    real(dp)		  :: M
    real(dp)		  :: zF, zR
    real(dp)		  :: ub, lb
    real(dp)              :: Fu, Fl
    real(dp)		  :: deltaF
    integer               :: boundstat

! Format specifier
    character(*), parameter :: M001 = '(3X, "deltaF: ", 3X, F20.15)'
	
! Execution section

    call countline (Ffile, nline)
    nF = nline - 2
    call countline (Rfile, nline)
    nR = nline - 2
	
    allocate (uF(nF))
    allocate (uR(nR))
	
    call readarray (Ffile, nF, uF)                                      ! Read the values of forward work function
    call readarray (Rfile, nR, uR)	                                ! Read the values of reverse work function

    call zwanzig (uF, nF, zF)                                           ! Perform exponential averaging for forward work functions 
    call zwanzig (uR, nR, zR)                                           ! Perform exponential averaging for reverse work functions

    ub =   MAX (zF, zR)                                                 ! Upper bound of the solution
    lb =   MIN (zF, zR)                                                 ! Lower bound of the solution

    call BARequation (uF, uR, nF, nR, ub, Fu)
    call BARequation (uF, uR, nF, nR, lb, Fl)

    if  (Fu * Fl < zero) then
        boundstat = 1
        call BARSolver (uF, uR, nF, nR, ub, lb, boundstat, deltaF)
        write(*, M001) deltaF
        write(*,*) 'Solution is within first bound'
    else
        ub =   ABS (MAX (zF, zR))                                       ! If the first bound doesn't work, widening the bounds
        lb = - ABS (MAX (zF, zR))
        call BARequation (uF, uR, nF, nR, ub, Fu)
        call BARequation (uF, uR, nF, nR, lb, Fl)

        if  (Fu * Fl < zero) then
                boundstat = 2   
                call BARSolver (uF, uR, nF, nR, ub, lb, boundstat, deltaF)
                write(*, M001)  deltaF
                write(*,*) 'Solution is within second bound'
        else
                ub =   1000*(ABS (MAX (zF, zR)))                         ! If the second bound doesn't work, widening the bounds further
                lb = - 1000*(ABS (MAX (zF, zR)))
                call BARequation (uF, uR, nF, nR, ub, Fu)
                call BARequation (uF, uR, nF, nR, lb, Fl)

            if  (Fu * Fl < zero) then  
                        boundstat = 3 
                        call BARSolver (uF, uR, nF, nR, ub, lb, boundstat, deltaF)
                        write(*, M001)  deltaF
                        write(*,*)'Solution is within third bound'
            else 
            write(*,*)'Bounds are still not good enough! Increase it further more.'
            end if   

        end if

     end if




    deallocate (uF, stat = alloc_err)           
    deallocate (uR, stat = alloc_err)  

end program BARMain
