program testit
use reptation
implicit none

double precision :: N          =  12790
double precision :: Ne         =  205.228
double precision :: lseg       =  4.05
double precision :: wl4        =  61785.41
double precision :: Re         =  48.97
double precision :: alpha0     =  0.0651165
double precision :: talphamax  =  8.5
double precision :: talphawd   =  0.6717

double precision :: q          =  0.096
double precision :: t
double precision :: sq(2)

integer :: i


open(10,file='rptest.dtr')
write(10,'(a)')"rptest.dtr"
write(10,'(a)')"rptest  sqt vs t  1"
write(10,'(a,f12.6,i6)')"q ", q, 0
write(10,*)
do i=0,150
  t  = 0.001d0 * exp(0.1d0*i)
  sq =  reptation_sqt(q,t, N, lseg, Ne, Re, wl4, alpha0, talphamax,talphawd) 
  write(*,*) t, sq, sq(2)/sq(1)
  write(10,'(3f14.6)') t, sq(2)/sq(1), 0.001
enddo
write(10,*)
write(10,'(a)')"#nxt"
close(10)


end program testit
