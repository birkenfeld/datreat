program testfor
implicit none
real :: x
integer          :: i
character(len=100) :: c

do i=-10,10
  x = 5d0**i
  write(6,'(i8,3g18.8)')i,x,Huge(x),-Huge(x)
enddo

c="abc"//char(10)//" def"
write(6,'(a)') trim(c)

end program testfor
