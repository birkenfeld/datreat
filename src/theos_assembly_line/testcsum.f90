program testcsum
implicit none
double precision, parameter :: Pi=4*atan(1d0)
double precision            :: cn, cm, csum
double precision            :: an, am
integer                     :: n, nn, mm, p

write(6,*)"Enter n, nn > :"
read(5,*) n, nn

do mm=1,n

csum = 0
do p=1,n-1
  cn   = cos(p*Pi*nn/dfloat(n))
  cm   = cos(p*Pi*mm/dfloat(n))
  csum = csum + (cn-cm)**2 / p**2
enddo

csum = csum * 2* n /Pi**2

write(6,*) nn, mm, csum

enddo


end program testcsum

