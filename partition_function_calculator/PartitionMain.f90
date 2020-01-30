!********************************************************************************************************************************!      
!                                                                                                                                !
!     Purpose: Calculate rotational, translational and vibrational partition functions from optimized CONTCAR and wave numbers   !
!                           Author: Osman Mamun, Department of Chemical Engineering,USC                                          !
!                             Date: 07.15.15                                                                                     !
!                     Modification:                                                                                              !
!          Reasons of Modification:                                                                                              !
!                            Notes:                                                                                              !
!                                                                                                                                !
!********************************************************************************************************************************!

Program partcalc
use constants,          only : dp, sl, Pi, kB, kB_SI, planck, planck_SI, light_cm  
use PartitionProcedures

implicit none

   character(len=1)      :: name(3)
   real(dp), allocatable :: T(:)
   real(dp), allocatable :: x(:), y(:), z(:)
   real(dp), allocatable :: qrot(:), qvib(:), qtrans(:), lnqrot(:), lnqtrans(:), f(:), lnqvib(:), lnq(:), qtot(:)
   real(dp)              :: multiplier, IXX, IYY, IZZ, IXY, IYZ, IZX, xmass, ymass, zmass, qvi, ZPE
   real(dp)              :: Imatrix(3,3), A(3,3)
   real(dp)              :: n(3), m(3), E(3)
   real(dp)              :: msum, mxsum, mysum, mzsum, mu, sigma, d, P
   integer(dp)           :: nTemp, ntotal, i, niter, nnn, Linearity, ndiff, nvib, j
   real(dp)              :: LT = 273, HT = 973, IT = 25

   P   = 101325


    write(*,*) 'How many different types of atoms do you have in your system? (Maximum is 3)'
    read (*,*) ndiff  
    write(*,*) 'What is the linearity of the molecule? 1=linear 2=non-linear'
    read (*,*) Linearity
    nTemp=((HT-LT)/IT)+1
   
    allocate (T(nTemp))
    allocate (qrot(nTemp))
    allocate (lnqrot(nTemp))
    allocate (qvib(nTemp))
    allocate (lnqvib(nTemp))
    allocate (qtrans(nTemp))
    allocate (lnqtrans(nTemp))
    allocate (lnq(nTemp))
    allocate (qtot(nTemp))

    do i = 1, nTemp
       T(i) = LT + (i-1) * IT
    end do
    write(*,*) 'What is the symmetry number of the molecule involved?'
    read(*,*) sigma

    Open(unit = 10,file = 'CONTCAR',status = 'old')
    read(10,*)
    read(10,*) multiplier
    read(10,*)
    read(10,*)
    read(10,*)    
    if (ndiff == 1) then
          read(10,*) name(1)
          read(10,*) n(1)
          n(2) = 0
          n(3) = 0
          m(2) = 0
          m(3) = 0
          ntotal = n(1)
    else if (ndiff == 2) then
          read(10,*) name(1), name(2)
          read(10,*) n(1), n(2)
          ntotal = n(1) + n(2)
          n(3) = 0
          m(3) = 0
    else
          read(10,*) name(1), name(2), name(3)
          read(10,*) n(1), n(2), n(3)
          ntotal = n(1) + n(2) + n(3)
    end if
    allocate (x(ntotal))
    allocate (y(ntotal))
    allocate (z(ntotal))
    read(10,*)

    do i = 1, ntotal
      read(10,*) x(i), y(i), z(i)
      x(i) = x(i) * multiplier
      y(i) = y(i) * multiplier
      z(i) = z(i) * multiplier
    end do
    close(10)

    if (Linearity == 1) then
        nvib = 3 * ntotal - 5
    else 
        nvib = 3 * ntotal - 6
    end if

    allocate (f(nvib))
    ZPE = 0
    Open(unit = 20, file = 'freq.dat', status = 'old')

    do i = 1, nvib
      read(20,*) f(i)
      ZPE = ZPE + 0.5 * planck * light_cm * f(i)
    end do
    close(20)

    do i = 1, nTemp
      qvi  = 1
      do j = 1, nvib
        qvi = qvi * (1 / (1 - exp( -planck * f(j) * light_cm / kB / T(i) ) ) )
      end do
      qvib(i)   = qvi
      lnqvib(i) = log(qvib(i))
    end do


    do i = 1, ndiff
       if (name(i) == 'C') then
           m(i)=12.0107
       else if (name(i) == 'O') then
           m(i)=15.9994
       else if (name(i) == 'H') then
           m(i)=1.00794
       else 
           write(*,*) "Not eligible to use this code, please see the README file"
           stop
       end if
    end do
    mxsum = 0
    mysum = 0
    mzsum = 0
    do i = 1, ntotal
       if (i <= n(1)) then
           mxsum = mxsum + m(1) * x(i)
           mysum = mysum + m(1) * y(i)
           mzsum = mzsum + m(1) * z(i)
       else if (i > n(1) .and. i <= n(1) + n(2)) then
           mxsum = mxsum + m(2) * x(i)
           mysum = mysum + m(2) * y(i)
           mzsum = mzsum + m(2) * z(i)
       else
           mxsum = mxsum + m(3) * x(i)
           mysum = mysum + m(3) * y(i)
           mzsum = mzsum + m(3) * z(i)
       end if
    end do
    msum  = m(1) * n(1) + m(2) * n(2) + m(3) * n(3)
    xmass = mxsum / msum
    ymass = mysum / msum
    zmass = mzsum / msum
    if (Linearity == 2) then
        IXX = 0
        IYY = 0
        IZZ = 0
        IXY = 0
        IYZ = 0
        IZX = 0
        do i = 1, ntotal
           if (i <= n(1)) then
               IXX = IXX + m(1) * ((y(i)-ymass)**2+(z(i)-zmass)**2)
               IYY = IYY + m(1) * ((z(i)-zmass)**2+(x(i)-xmass)**2)
               IZZ = IZZ + m(1) * ((x(i)-xmass)**2+(y(i)-ymass)**2)
               IXY = IXY + m(1) * (x(i)-xmass)*(y(i)-ymass)
               IYZ = IYZ + m(1) * (y(i)-ymass)*(z(i)-zmass)
               IZX = IZX + m(1) * (z(i)-zmass)*(x(i)-xmass)
           else if (i > n(1) .and. i <= n(1) + n(2)) then
               IXX = IXX + m(2) * ((y(i)-ymass)**2+(z(i)-zmass)**2)
               IYY = IYY + m(2) * ((z(i)-zmass)**2+(x(i)-xmass)**2)
               IZZ = IZZ + m(2) * ((x(i)-xmass)**2+(y(i)-ymass)**2)
               IXY = IXY + m(2) * (x(i)-xmass)*(y(i)-ymass)
               IYZ = IYZ + m(2) * (y(i)-ymass)*(z(i)-zmass)
               IZX = IZX + m(2) * (z(i)-zmass)*(x(i)-xmass)
           else
               IXX = IXX + m(3) * ((y(i)-ymass)**2+(z(i)-zmass)**2)
               IYY = IYY + m(3) * ((z(i)-zmass)**2+(x(i)-xmass)**2)
               IZZ = IZZ + m(3) * ((x(i)-xmass)**2+(y(i)-ymass)**2)
               IXY = IXY + m(3) * (x(i)-xmass)*(y(i)-ymass)
               IYZ = IYZ + m(3) * (y(i)-ymass)*(z(i)-zmass)
               IZX = IZX + m(3) * (z(i)-zmass)*(x(i)-xmass)
           end if
        end do
        Imatrix(1,1) = IXX
        Imatrix(2,2) = IYY
        Imatrix(3,3) = IZZ
        Imatrix(1,2) = IXY
        Imatrix(1,3) = IZX
        Imatrix(2,3) = IYZ
        Imatrix(2,1) = IXY
        Imatrix(3,2) = IYZ
        Imatrix(3,1) = IZX  
        A   = Imatrix
        nnn = 3
        call Eigenval(A,nnn,E,niter)
        do i = 1, ntemp
          qrot(i) = (Pi*Pi*Pi*Pi**0.5/sigma)*((8*kB_SI*T(i)*1e-23/Na)**1.5)*((E(1)*E(2)*E(3))**0.5)/planck_SI/planck_SI/planck_SI
        end do
    else 
        if (ndiff == 1) then
            do i = 1, ntemp
               d = ((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
               qrot(i) = (1/sigma)*(8*Pi*Pi*T(i)*kB_SI/Na/planck_SI/planck_SI)*(m(1)/2)*d*1e-23
            end do
        else 
            mu = (m(1) * m(2)) / (m(1)+m(2))
            do i = 1, ntemp
               d = ((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
               qrot(i) = (1/sigma)*(8*Pi*Pi*T(i)*kB_SI/Na/planck_SI/planck_SI)*(mu/4)*d*1e-23
            end do 
        end if
    end if
   
    call fileforwrite()
    write(30, 610) ZPE
    610 format (F17.15)
    write(30,620) msum/1000
    620 format (F7.5)

    do i = 1, nTemp
        lnqrot(i)   = log(qrot(i))
        qtrans(i)   = ((2*Pi*msum*kB_SI*T(i)/planck_SI/planck_SI/1000/Na)**1.5)*(kB_SI*T(i)/P)
        lnqtrans(i) = log(qtrans(i))
        lnq(i)      = lnqtrans(i) + lnqrot(i) + lnqvib(i)
        qtot(i)     = exp(lnq(i))

        write(20,500) T(i),qrot(i),qtrans(i),qvib(i),lnqrot(i),lnqtrans(i),lnqvib(i),lnq(i)
              500 format (3X,F8.2,3X,7(ES20.13,3X))
        write(30,630) T(i), qtot(i)
              630 format (F6.2, 3X, ES20.13)
    end do

    write(20,510)
          510 format (3X,'===========================================================', &
                         '===========================================================', &
                         '===========================================================')
    write(20,*) 'ZPE=',ZPE,'eV'

    close(20) 
    close(30)

end program partcalc
