program confidencemain
  use confidenceroutines
  implicit none



  integer               :: npoints
  integer               :: alloc_err, status = 0
  real(dp), allocatable :: Grxn(:), Gact(:)
  real(dp)	        :: cmean, cvaru, cvarl
  real(dp)		:: Mean, Variance
  real(dp)		:: CIMean
  real(dp)              :: Ubmean, Lbmean
  real(dp)		:: Ubvar, Lbvar
  character(30)		:: Ifile = 'InputCONF'
  character(30)		:: Ofile = 'OutputCONF'

  call countline(Ifile, npoints)
  allocate (Grxn(npoints))
  allocate (Gact(npoints))
  if ( npoints == 3 ) then
      cmean = 4.30
      cvarl = 0.05
      cvaru = 7.38
  elseif ( npoints == 5 ) then
      cmean = 2.78
      cvarl = 0.48
      cvaru = 11.14
  elseif ( npoints == 10 ) then
      cmean = 2.26
      cvarl = 2.70
      cvaru = 19.02
  else           
     write (*,*) "This script is written only for 3, 5, or 10 datapoints,", &
     	         " To use it for other number of data points, check the", &
                 " reference listed in header"
  	 call exit (status)
  end if  
   
  call readvalues ( Ifile, npoints, Grxn, Gact )  
  
  call calculatestat ( npoints, Grxn, Mean, Variance)
  call calculatebounds ( npoints, cmean, cvarl, cvaru, Variance, CIMean, Ubvar, Lbvar, Ubmean, Lbmean )
  call fileforwrite ( Ofile )
  write (20, 100) 'Grxn', Mean, Variance, CIMean, Ubvar, Lbvar, Ubmean, Lbmean
  100 format (1X, A5, 2X, F15.8, 2X, F15.8, 2X, F15.8, 5X, F15.8, 3X, F15.8, 2X, F15.8, 2X, F15.8)
  
  call calculatestat ( npoints, Gact, Mean, Variance)
  call calculatebounds ( npoints, cmean, cvarl, cvaru, Variance, CIMean, Ubvar, Lbvar, Ubmean, Lbmean )
  write (20, 110) 'Gact', Mean, Variance, CIMean, Ubvar, Lbvar, Ubmean, Lbmean
  110 format (1X, A5, 2X, F15.8, 2X, F15.8, 2X, F15.8, 5X, F15.8, 3X, F15.8, 2X, F15.8, 2X, F15.8)
  
  deallocate (Grxn, stat = alloc_err)
  deallocate (Gact, stat = alloc_err)
  CLOSE(20)

end program   
  
