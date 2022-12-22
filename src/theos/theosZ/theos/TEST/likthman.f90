program likthman_tst
implicit none    
double precision :: tx, sl, q, d, W, l, alpha, n
integer          :: i
write(*,*)
q = 0.1d0
d = 50d0
l = 4d0
W = 300d0
n = 500
alpha = 1.5d0
do i=1,200
  tx = i*5
  sl = likhtman_clf(tx,q,d,W,n,l,alpha) 
  write(*,'(i8,2f12.6)')i, tx, sl
enddo

contains
 
function likhtman_clf(t,q,d,W,n,l,alpha) result(sumlik)
  implicit none 
  double precision, intent(in ) :: t     ! time
  double precision, intent(in ) :: q     ! momentum transfer
  double precision, intent(in ) :: d     ! tube diameter
  double precision, intent(in ) :: W     ! Rouse rate
  double precision, intent(in ) :: n     ! Rouse rate
  double precision, intent(in ) :: l     ! segment length 
  double precision, intent(in ) :: alpha ! CLF exponent (default approx 1.5)

  double precision, parameter :: Pi = 3.141592654d0
  double precision :: sumlik
  double precision :: sumlik0, st
  double precision :: mue
  double precision :: taue
  double precision :: tau0
  double precision :: num
  double precision :: ne


       ne   =  (d/l)**2
       tau0 =  36.d0/(W*((q*l)**4))
       taue =  ne**2/((pi**2)*W)

       mue  =  q**2*n*l**2/12.d0
       num  =  n/ne
       st   =  min(0.5d0*((alpha/num)*(t/taue)**0.25d0), 0.5d0)
       sumlik  = (n/(2*mue**2))*(2*mue+exp(-2*mue)+2.d0-4*mue*st-4 * exp(-2*mue*st)+exp(-4d0*mue*st))
       sumlik0 = (n/(2*mue**2))*(2*mue+exp(-2*mue)-1.d0)
       sumlik  = sumlik/sumlik0

end function likhtman_clf

end program likthman_tst
 
