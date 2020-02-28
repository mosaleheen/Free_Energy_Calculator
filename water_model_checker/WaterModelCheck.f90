!********************************************************************************************************************************!      
!                                                                                                                                !
!                          Purpose: Check the water model from a CONFIG file format                                              !
!                           Author: Mohammad Shamsus Saleheen, Department of Chemical Engineering,USC                            !
!                             Date: 04.11.15                                                                                     !
!                     Modification: Need to account for the periodicity of the box                                                                                             !
!          Reasons of Modification:                                                                                              !
!                            Notes:                                                                                              !
!                                                                                                                                !
!********************************************************************************************************************************!


program Watercheck
implicit none
integer, parameter     :: dp = selected_real_kind(15, 307)
integer                :: ierror, j, totalwater, oxygenid, hydrogen1id, hydrogen2id
character(5),parameter :: oxygen ='Ow'
character(5)           :: atom
real(dp)               :: xo, yo, zo, xh1, yh1, zh1, xh2, yh2, zh2
real(dp)               :: oh1, oh2, h1h2
character(30)          :: filename

write(*,*)"Enter you coordinate filename:"
read(*,*)filename
OPEN(UNIT=20,file='CONFIG_Water_Model',access='sequential',&
     status='replace',action='write',position='append',&
     iostat=ierror)
     write(20,100)
     100 format(1X,'Atom Id ',3X,' OH1-distance ',5X,' OH2-distance',5X,&
               ' H1H2-distance ')
     write(20,110)
     110 format(1X,'========',3X,'==============',5X,'==============',5X,&
               '================')
OPEN(UNIT=10,file=filename,status='old',action='read',&
     iostat=ierror)
     do j=1,5
     read(10,*,iostat=ierror)
     end do
     totalwater=0
     do 
     read(10,*,iostat=ierror)atom,oxygenid
     read(10,*,iostat=ierror)xo,yo,zo
        if (ierror /= 0) exit
          if ( atom .eq. oxygen) then
               read(10,*,iostat=ierror)atom,hydrogen1id
               read(10,*,iostat=ierror)xh1,yh1,zh1
               read(10,*,iostat=ierror)atom,hydrogen2id
               read(10,*,iostat=ierror)xh2,yh2,zh2
               oh1=sqrt((xo-xh1)**2+(yo-yh1)**2+(zo-zh1)**2)
               oh2=sqrt((xo-xh2)**2+(yo-yh2)**2+(zo-zh2)**2)
               h1h2=sqrt((xh1-xh2)**2+(yh1-yh2)**2+(zh1-zh2)**2)
               totalwater=totalwater+1
               write(20,120)oxygenid,oh1,oh2,h1h2
               120 format(1X,I9,3X,F15.12,5X,F15.12,5X,F15.12)
          end if
      end do
      write(*,*)"Total number of water molecules found: ",totalwater
         
CLOSE(10)
CLOSE(20)        

end program Watercheck
