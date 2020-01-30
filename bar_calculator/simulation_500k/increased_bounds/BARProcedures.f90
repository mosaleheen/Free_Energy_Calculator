module BARSubroutines
use constants
implicit none

! Selective parameters
    real(dp), parameter :: T     	= 500.00									! K
    real(dp), parameter :: beta  	= one / (kBha * T)								! beta in units of hartree
    real(dp), parameter :: tol   	= 1.0e-15									! tolerance for BAR free energy solver
    integer,  parameter :: iter  	= 1000										! number of iterations for BAR solver
    character(*), parameter :: F001     = '(4X, I3, 6X, F22.15, 7X, F22.15)'
	
contains

!****************************************************************************************
!	                   Procedure for counting lines of a file 			!
!****************************************************************************************
    
    subroutine countline (filename, nline)			        				
    implicit none

	character(20), intent(in) 	:: filename
	integer,       intent(out)	:: nline
	integer			        :: ierror
	
	nline=0
     	OPEN(UNIT = 10, file = filename, status = 'old',&
             action = 'read', iostat = ierror)
         	 do
           		read(10, *, iostat = ierror)
           		if (ierror /= 0) exit
		   	nline = nline + 1
         	 end do
     	CLOSE(UNIT = 10)
    
    end subroutine

!****************************************************************************************
!			 Procedure for reading workfunction array		        !
!****************************************************************************************

        
    subroutine readarray (filename, Narr, wf)
    implicit none
	character(20), intent(in) 	:: filename
    	integer,       intent(in) 	:: Narr	
    	real(dp),      intent(out) 	:: wf(Narr)
        integer			        :: i, ierror
        
        OPEN(Unit = 10, file = filename, status = 'old', &
             action = 'read', iostat = ierror)
             read(10, * , iostat = ierror)
             read(10, * , iostat = ierror)	
             do i = 1, Narr
            	read(10, *, iostat = ierror) wf(i)
             if  (ierror /= 0) exit
			end do
    	CLOSE(10)  
    
    end subroutine


!****************************************************************************************
!			    Procedure to employ Zwanzig equation 			!
!****************************************************************************************


    
    subroutine zwanzig	(wf, Narr, Expavg)									
    implicit none
    	integer,  intent(in)		:: Narr	
    	real(dp), intent(in)		:: wf(Narr)
    	real(dp), intent(out)		:: Expavg
    	integer				:: i
	real(dp)			:: xsum, xavg
    
    	xsum = 0.0
        do i = 1, Narr
           xsum = xsum + EXP(-beta*wf(i))
        end do  
			
	xavg 	= xsum / REAL(Narr)
	Expavg	= -LOG (xavg) / beta

    end subroutine


!****************************************************************************************
!		            Procedure to employ BAR equation 				!
!******************************************************************************************************************   
!                                                                                                                 !
! Notes: Something very important to understand here is the sign used for reverse workfunction. The sign is       !
! opposite from Shirts et.al. because a forward work function gives you U(2)-U(1), while reverse work function    !
! gives you U(1)-U(2); if you want to find the optimum of workfunction that minimizes variance using this         !
! procedure, you have to use both forward and reverse as U(2)-U(1) or U(1)-U(2). In this case, I switched         !
! the sign of the reverse one.                                                                                    !
!                                                                                                                 !
!******************************************************************************************************************



    subroutine BARequation (uF, uR, nF, nR, deltaA, Fzero)			
    implicit none
    	integer,  intent(in)		:: nF, nR
        real(dp), intent(in)		:: uF(nF), uR(nR)
        real(dp), intent(in)		:: deltaA
        real(dp), intent(out)		:: Fzero
        real(dp)			:: M
        real(dp)			:: Fsum, Rsum
        integer				:: i
        				 
	M       = log (REAL(nF) / REAL(nR)) 
	Fsum	= 0.0
        do i = 1, nF
           Fsum = Fsum + one / ( 1 + exp ( M + beta * uF(i) - beta * deltaA  ) ) 
	end do
            
	Rsum 	= 0.0
        do i = 1, nR
           Rsum = Rsum + one / ( 1 + exp ( - ( M - beta * uR(i) - beta * deltaA  ) ) ) 
	end do            			
	Fzero = ( Fsum / REAL(nF) ) - ( Rsum / REAL (nR) )
            
    end subroutine



!****************************************************************************************
!		 Procedure to solve BAR using False position method			!
!********************************************************************************************************* 
!                                                                                                        !
! Notes: Depending on the 'boundstat' value, it first writes down the place                              !
! where the solution lies, then it employs false position method to find the                             !
! solution within that bound.                                                                            !
! For false position method,                                                                             !
! Xupper - Xlower     Xupper - Xnew                                                                      !
!----------------- = ---------------, which gives you Xnew, also defined as deltaAr in this procedure    !
! Fupper - Flower     Fupper - 0                                                                         !
!                                                                                                        !
!*********************************************************************************************************

    subroutine BARSolver (uF, uR, nF, nR, ub, lb, boundstat, deltaF)				 
    implicit none
	integer,  intent(in)	:: nF, nR
    	real(dp), intent(in)	:: uF(nF), uR(nR)
    	real(dp), intent(in)	:: ub, lb
        integer,  intent(in)    :: boundstat
	real(dp), intent(out)   :: deltaF
	real(dp)		:: deltaAu, deltaAl, deltaAr 
        real(dp)                :: Fu, Fl, Fr
    	integer			:: i, ierror
        character(20)           :: Ofile = 'BARdelta'
        
        OPEN(Unit = 20, file = Ofile, access = 'sequential', &
             status = 'replace', action = 'write', position = 'append', &
             iostat = ierror)

             if ( boundstat == 1 ) then
             write(20,900) 
             900 format  (2X, "Solution lies within first bound")
             elseif ( boundstat == 2 ) then
             write(20,910)
             910 format  (2X, "Solution lies within second bound")
             else 
             write(20,920)
             920 format  (2X, "Solution lies within second bound")
             end if
 
             write(20, 1000)
             1000 format (2X, 'Iteration', 9X, 'Function Value',&
                          8X, 'Perturbation Energy(Hartree)')
             write(20,1010)
             1010 format (2X, '=========', 5X, '=====================',&
                          5X, '============================')
         deltaAu = ub
         deltaAl = lb
	 do i = 1, iter 
           call BARequation (uF, uR, nF, nR, deltaAu, Fu)
           call BARequation (uF, uR, nF, nR, deltaAl, Fl)
           deltaAr = deltaAu - Fu*((deltaAu - deltaAl)/(Fu - Fl))
           call BARequation (uF, uR, nF, nR, deltaAr, Fr)
           
           if   (abs (zero - Fr) < tol)	 then
	        deltaF = deltaAr
                write(20,F001) i, Fr, deltaF
                exit
           elseif   (Fr * Fu < zero) 	 then    
            	deltaAl = deltaAr
                write(20,F001) i, Fr, deltaAl
           elseif   (Fr * Fl < zero)	 then	
                deltaAu = deltaAr
                write(20,F001) i, Fr, deltaAu
           else
                deltaF = deltaAr
                write(20,1020)
                1020 format  (3X, "Convergence couldn't be reached within specified interation limit")
                write(*,1030) 
                1030 format  (3X, "Value of deltaF from last iteration:")
                write(20,F001) i, Fr, deltaF         
           end if
           
          end do
         CLOSE(20)            
         
   end subroutine

 
end module BARSubroutines
