program tstread
 implicit none
 character(len=256) :: line
 double precision :: x(10)
 integer          :: ios
 integer          :: i, j, n

 do i=1,10
   line =  " "
   read(5,'(a)',iostat=ios,end=99) line
   write(*,*)"ios=",ios
   x = 0
   read(line,* ,iostat=ios) x
   write(*,*)"ios=",ios
   write(*,*)"line=",trim(line)
   write(*,'(a,10f10.4)')"x   =",x

 enddo
99 continue
end program tstread 
