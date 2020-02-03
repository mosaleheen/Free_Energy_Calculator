module confidenceroutines
   implicit none
   integer, parameter    :: dp = selected_real_kind(15,307)

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
!			 Procedure for reading Grxn and Gactivetion		        !
!****************************************************************************************

        
   subroutine readvalues (filename, Narr, wi, wf)
   implicit none
   character(30), intent(in) 	:: filename
   integer,       intent(in) 	:: Narr
   real(dp),	  intent(out)   :: wi(Narr)	
   real(dp),      intent(out) 	:: wf(Narr)
   integer				        :: i, ierror
        
   OPEN(Unit = 10, file = filename, status = 'old', &
        action = 'read', iostat = ierror)
        do i = 1, Narr
           read(10, *, iostat = ierror) wi(i), wf(i)
           if  (ierror /= 0) exit
	    end do
   CLOSE(10)      
   end subroutine


!****************************************************************************************
!			 Procedure for calculating mean and variance		        !
!****************************************************************************************
   subroutine calculatestat ( Narr, Gval, mean, variance)
   implicit none
   integer,  intent(in)  :: Narr
   real(dp), intent(in)  :: Gval(Narr)
   real(dp), intent(out) :: mean
   real(dp), intent(out) :: variance
   integer				 :: i
   real(dp)				 :: Gsum, Gsumsq

   Gsum = 0.0
   do i = 1, Narr
     Gsum = Gsum + Gval(i)
   end do
   mean = Gsum / REAL(Narr)
   Gsumsq = 0.0
   do i = 1, Narr
     Gsumsq = Gsumsq + ( Gval(i) - mean ) ** 2
   end do
   variance = Gsumsq / REAL ( Narr - 1.0 )              
   end subroutine

!****************************************************************************************
!			 Procedure for calculating bounds for mean and variance		!
!****************************************************************************************  
   subroutine calculatebounds ( Narr, cmean, cvarl, cvaru, variance, CIMean, Ubvar, Lbvar, Ubmean, Lbmean)
   implicit none
   integer,  intent(in)  :: Narr
   real(dp), intent(in)  :: cmean, cvarl, cvaru, variance
   real(dp), intent(out) :: CIMean, Ubvar, Lbvar, Ubmean, Lbmean
   
   CIMean = cmean * sqrt( Variance ) / sqrt( REAL(Narr) )
   Ubvar  = ( REAL(Narr) - 1.0 ) * Variance / cvarl
   Lbvar  = ( REAL(Narr) - 1.0 ) * Variance / cvaru
   Ubmean = cmean * sqrt( Ubvar ) / sqrt( REAL(Narr) )
   Lbmean = cmean * sqrt( Lbvar ) / sqrt( REAL(Narr) )

   end subroutine


!****************************************************************************************
!			 Procedure for opening output file		                !
!****************************************************************************************  

  subroutine fileforwrite ( Ofile )
  implicit none
  character(30), intent(in) :: Ofile
  integer                   :: ierror

  OPEN (UNIT = 20, file = Ofile, access= 'sequential', &
    	status = 'replace', action = 'write', position = 'append', &
        iostat = ierror)
  write (20,500)
 		 500 format(16X, "Mean", 11X, "Variance", 6X, "Confidence Mean", &
         7X, "UB Variance", 7X, "LB Variance", 7X, &
         "UB Mean", 10X, "LB mean")
  write (20,510)
 		 510 format(16X, "====", 11X, "========", 6X, "===============", &
         7X, "===========", 7X, "===========", 7X, &
         "=======", 10X, "=======")
  end subroutine
   
end module confidenceroutines
