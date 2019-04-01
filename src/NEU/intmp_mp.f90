program intmp_mp
use integration_mp

implicit none

integer, parameter :: N = 100
double precision   :: a, b, z(N), par(N)
integer            :: i

call integral_setmethod_mp(5)
call setmaxiteration_mp(1000)
call integral_setaccuracy_mp(1d-8)

do i=1,N
  par(i) = 0.1d0*i
enddo

!$OMP PARALLEL DO PRIVATE(i,a,b) SHARED(z) SCHEDULE(dynamic)
do i=1,N
  a    = 0d0
  b    = 0.1d0*i
  z(i) = integral_mp(f,i,a,b)
enddo
!$OMP END PARALLEL DO

write(*,'(i8,f18.9)')(i,z(i),i=1,N)
write(*,*)
write(*,'(f18.9)') sum(z)

contains

function f(x,i) result(y)
 implicit none
 double precision, intent(in)  :: x
 integer         , intent(in)  :: i
 double precision              :: y
 integer                       :: ii
 
 y = 0
 do ii=1,10000
   y = y+1d0/sqrt(x*ii) 
 enddo 

 y = y * par(i)

end function f

end program intmp_mp
