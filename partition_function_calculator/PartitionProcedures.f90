module PartitionProcedures
use constants
implicit none

! Selective format specifiers
  character(sl), parameter :: F001 = '(A)'

contains



!*******************************************************************************************!
!             Procedure for calculating eigenvalues using Jacobi's algorithm                !
!*******************************************************************************************! 

   Subroutine eigenval(A,nnn,E,niter)
   implicit none

   integer(dp), intent(in)      :: nnn
   real(dp),    intent(inout)   :: A(nnn,nnn)
   real(dp),    intent(out)     :: E(nnn)
   integer(dp), intent(out)     :: niter
   real(dp)                     :: ADM(nnn,nnn),S1(nnn,nnn),S1T(nnn,nnn),D1(nnn,nnn),D(nnn,nnn)
   integer(dp)                  :: m(2)
   real(dp)                     :: Theta, t, c, s, Amax
   integer(dp)                  :: i, j

   niter=0

   do
      niter = niter + 1
      do i = 1, nnn
         do j = 1, nnn
            if ( i == j ) then
               ADM(i,j) = 0
            else
               ADM(i,j) = abs(A(i,j))
            end if
         end do
      end do
   
      Amax = maxval(ADM)

      do i = 1, nnn
         do j = 1, nnn
            if (abs(ADM(i,j)) == Amax) then
                m(1) = j
                m(2) = i
            end if
         end do
      end do 
      if (Amax < 0.0001) exit
      Theta = (A(m(1),m(1)) - A(m(2),m(2))) /2 / A(m(2),m(1))
      if (Theta < 0) then
          t = -1 / (abs(Theta) + (Theta * Theta + 1)**0.5)
      else
          t =  1 / (abs(Theta) + (Theta * Theta + 1)**0.5)
      end if
      c = 1 / (t**2 + 1)**0.5
      s = c * t

      do i = 1, nnn
         do j =1, nnn
            if (i==j) then
               S1(i,j) = 1
            else
               S1(i,j) = 0
            end if
         end do
      end do
      S1(m(1),m(1)) =  c
      S1(m(2),m(2)) =  c
      S1(m(2),m(1)) =  s
      S1(m(1),m(2)) = -s
      S1T = transpose(S1)
      D1  = matmul(S1T,A)
      D   = matmul(D1,S1)
      do i = 1, nnn
         do j = 1, nnn
            if (abs(D(i,j)) < 0.0001) then
                D(i,j) = 0
            else
                D(i,j) = D(i,j)
            end if
         end do
      end do
      A = D
      end do
      do i =1, nnn
         do j = 1, nnn
            if ( i == j ) then
                E(i) = A(i,j)
            end if
         end do
      end do

   End Subroutine eigenval

!************************************************************************************************!
!            Procedure for opening and formatting files to write outputs                         !
!************************************************************************************************!

   Subroutine fileforwrite()
   implicit none
   integer      :: ierror
   
   OPEN ( unit = 20, file = 'PartitonFunctions', access = 'sequential',&
          status = 'replace', action = 'write', position ='append',&
          iostat = ierror)

          write (20, 900)
          900 format(7X,"T",15X,"qrot",16X,"qtrans",16X,"qvib",16X,"ln(qrot)",16X,"ln(qtrans)",16X,"ln(qvib)",16X,"lnq")
          write (20, 910)
          910 format (3X,'===========================================================', &
                        '===========================================================', &
                        '===========================================================')
   OPEN ( unit = 30, file = 'Gspecies', access = 'sequential',&
          status = 'replace', action = 'write', position ='append',&
          iostat = ierror)

   End Subroutine fileforwrite

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
   nn = nn + 1; msgs(nn) = '**************************        USAGE Policy           ***************************'                
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = '                  PartitionCalculator -ndiff -linearity -symmetry                   '
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = '            -ndiff =  different types of atoms present in the system                '
   nn = nn + 1; msgs(nn) = '                                                                                    '
   nn = nn + 1; msgs(nn) = ' See InputFiles directory for input file format and README file for applicability   '
   nn = nn + 1; msgs(nn) = '                                                                                    '
   
   ! print help and terminate the program
   do ii = 1, nn
      write (*,F001) trim(msgs(ii))
   end do
   stop
   End Subroutine print_help

end module PartitionProcedures
