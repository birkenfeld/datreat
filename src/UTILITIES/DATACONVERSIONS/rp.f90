! read SANS DATA: 
! use gfortran rp.f90
! a.out filename <conc> <collen> <numor> > datfile
program rp
implicit none
double precision :: q, iq, err, dx, conc, l2
integer          :: ios
integer          :: numor
CHARACTER(len=255) :: cmd, lin, fil

     CALL get_command_argument(1,fil)
     open(10, file=trim (fil))
     CALL get_command_argument(2,cmd)
     read(cmd,*) conc
     CALL get_command_argument(3,cmd)
     read(cmd,*) l2
     CALL get_command_argument(4,cmd)
     read(cmd,*) numor
     write(*,'(a)')trim(fil)
     write(*,'(a,a,i9)')fil(1:8),'  "I(q)"  vs q  ', numor
     write(*,'(a,f12.6)') "conc ", conc
     write(*,'(a,f12.6)') "l2 ", l2
     write(*,*)"values"
     do
       read(10,'(a)',iostat=ios)lin
       if(ios .ne.0 ) exit
       read(lin,*,iostat=ios) q, iq, err, dx
       if(ios .ne. 0) cycle
       if(iq  == 0d0) cycle
       write(*,'(3f14.7)') q, iq, err
     enddo 
     write(*,'(a)')"#nxt" 
     close(10)
end program rp 
