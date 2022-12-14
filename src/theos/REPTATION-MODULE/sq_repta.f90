! Compile and link with:
! gfortran -o sq_repta reptation_module.f90 sq_repta.f90 -fopenmp -O2 -static-libgfortran
!
program sq_repta
use reptation
implicit none

integer, parameter :: nt = 200    ! number of timepoints in a log scale

double precision :: N          =  12790
double precision :: Ne         =  205.228
double precision :: lseg       =  4.05
double precision :: wl4        =  61785.41
double precision :: Re         =  48.97
double precision :: alpha0     =  0.0651165
double precision :: talphamax  =  8.5
double precision :: talphawd   =  0.6717

double precision :: q          =  0.096
double precision :: t,  tmin, tmax, tadvance
double precision :: sq(2)

integer :: i, length, stat
character(len=128) :: val
  

call get_command_argument(1,val,length, stat);  
if(index(val,"help")>0) then
 write(*,'(a)') "s(q,t) modelling for repating polymer chains, call:" 
 write(*,'(a)') "sq_repta q t N lseg Ne Re Wl4 alpha0 tmax tw" 
 write(*,'(a)') "above numerical values are to be given for:" 
 write(*,'(a)') "q            =  momentum transfer (q-value)" 
 write(*,'(a)') "tmin         =  (Fourier)-time start" 
 write(*,'(a)') "tmax         =  (Fourier)-time end" 
 write(*,'(a)') "N            =  total number of polymer chain segments " 
 write(*,'(a)') "lseg         =  segment length (monomer) " 
 write(*,'(a)') "Ne           =  number of segments pre entanglement blob" 
 write(*,'(a)') "Wl4          =  Rouse rate" 
 write(*,'(a)') "alpha0       =  Non-Gaussianity corr. amplitued" 
 write(*,'(a)') "talphamax    =  Non-Gaussianity alpha(t) peak position" 
 write(*,'(a)') "talphawd     =  Non-Gaussianity alpha(t) peak width (log)" 
 stop
else
                                                read(val,*,iostat=stat) q
endif
call get_command_argument(2,val,length, stat);  read(val,*,iostat=stat) tmin
call get_command_argument(3,val,length, stat);  read(val,*,iostat=stat) tmax
call get_command_argument(4,val,length, stat);  read(val,*,iostat=stat) N
call get_command_argument(5,val,length, stat);  read(val,*,iostat=stat) lseg
call get_command_argument(6,val,length, stat);  read(val,*,iostat=stat) Ne
call get_command_argument(7,val,length, stat);  read(val,*,iostat=stat) Re
call get_command_argument(8,val,length, stat);  read(val,*,iostat=stat) Wl4
call get_command_argument(9,val,length, stat);  read(val,*,iostat=stat) alpha0
call get_command_argument(10,val,length, stat);  read(val,*,iostat=stat) talphamax
call get_command_argument(11,val,length,stat);  read(val,*,iostat=stat) talphawd

if(stat .ne.0) stop "wrong argument list"

tadvance = exp(log(tmax/tmin)/nt)

do i=0,nt
 t = tmin * tadvance**i
 sq =  reptation_sqt(q,t, N, lseg, Ne, Re, wl4, alpha0, talphamax,talphawd) 
 write(*,'(2f18.6)') t, sq(2)/sq(1)
enddo

end program sq_repta
