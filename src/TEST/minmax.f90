program minmax 
implicit none
double precision :: y(10,10),ymi,yma
integer          :: isel(10)
integer          :: i, j

do i=1,10
 do j=1,10
    y(i,j) = i+10*j
 enddo
enddo

isel = [3,4,5,0,0,0,0,0,0,0]
ymi  = minval(y(:,isel(1:3)))
yma  = maxval(y(:,isel(1:3)))

write(*,*)ymi,yma

end program minmax
