! read wims rheo data example 
! use gfortran rp_g.f90
! a.out filename <conc> <collen> <numor> > datfile
program rp
implicit none
double precision, dimension(1000) :: Freq=0, reG=0, imG=0, tan_delta=0, Temp=0, Eta=0
double precision :: tave, ttolerance = 3d0
integer          :: ios, i, anfang, ende, n
integer          :: numor
CHARACTER(len=255) :: cmd, lin, fil

     

     CALL get_command_argument(1,fil)
     open(10, file=trim (fil))
!     CALL get_command_argument(2,cmd)
!     read(cmd,*) conc
!     CALL get_command_argument(3,cmd)
!     read(cmd,*) l2
     CALL get_command_argument(2,cmd)
     read(cmd,*) numor

  
dr:  do i=1,size(Freq)
       read(10,'(a)',iostat=ios)lin
!write(*,*)i,lin(1:80)
       if(ios .ne.0 ) exit
       read(lin,*,iostat=ios) Freq(i), reG(i), imG(i), tan_delta(i), Temp(i), Eta(i) 
!       write(*,*) Freq(i), reG(i), imG(i), tan_delta(i), Temp(i), Eta(i) 
       if(ios .ne. 0) cycle
     enddo dr
     close(10)

   
     do i=1,size(Freq)
       if(Freq(i) .ne. 0d0) exit
     enddo
     anfang = i
dw:  do
       tave = temp(anfang)
       if(tave == 0d0) exit dw
da:       do i=anfang, size(temp)
write(*,*)"T:",i,tave,temp(i)
            if(abs(tave-temp(i)) > ttolerance) exit da
          enddo da
       ende = i
       n    = ende-anfang+1
       tave = sum(temp(anfang:ende))/n
      
       numor = numor + 1
       write(*,'(a)')trim(fil)
       write(*,'(a,a,i9)')fil(1:8),'  reG  vs q  ', numor
       write(*,'(a,f12.6)') "temp ",tave 
       write(*,'(a,f12.6)') "reG  ",1d0 
       write(*,'(a,f12.6)') "imG  ",0d0 
       write(*,*)"values"    
       do i=anfang, ende
         write(*,'(3e15.7)') Freq(i), reG(i), 0.02*abs(reG(i))
       enddo

       write(*,'(a)')"#nxt" 
      
       numor = numor + 1
       write(*,'(a)')trim(fil)
       write(*,'(a,a,i9)')fil(1:8),'  reG  vs q  ', numor
       write(*,'(a,f12.6)') "temp ",tave 
       write(*,'(a,f12.6)') "reG  ",0d0 
       write(*,'(a,f12.6)') "imG  ",1d0 
       write(*,*)"values"    
       do i=anfang, ende
         write(*,'(3e15.7)') Freq(i), imG(i), 0.02*abs(reG(i))
       enddo

       write(*,'(a)')"#nxt" 

       if(ende >= size(Freq)) exit dw
       anfang = ende+1

     enddo dw


end program rp 
