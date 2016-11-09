program tpgf
  implicit none
  double precision :: q         = 0.1d0 
  double precision :: tau       = 100d-3
  double precision :: diffusion = 3d0
  double precision :: lz        = 100d0

  double precision :: pfg1d, sqtx

  integer :: i
  do i=0,1000
     tau = i
     sqtx = pfg1d(q,tau,diffusion,lz)
     write(6,*) q, tau, sqtx
  enddo

end program tpgf
 

double precision function pfg1d(q,tau,diffusion,lz)

! ref. Hall P-L. and Ross D.K. , Molecular Physics , 1978, p 1549-1554
  implicit none

  double precision, parameter    :: Pi = 4*atan(1d0)
  double precision, intent(in)               :: q           ! q-value
  double precision, intent(in)               :: tau         ! t
  double precision, intent(in)               :: diffusion   ! diffusion
  double precision, intent(in)               :: lz          ! channel length
  
  double precision               :: eps   = 1d-8
  integer                        :: maxit = 1000
  double precision               :: erracc, adapint
 
 
  pfg1d = 0.5d0 * adapint(sqt_ker,0d0, Pi, eps, maxit, erracc) 



contains

double precision function a_n(n,qzl)
! --------------------------------
! Eqns (11)+(12) of ref
  implicit none
  integer         , intent(in)  :: n
  double precision, intent(in)  :: qzl

  double precision, parameter   :: Pi = 4*atan(1d0)
  double precision, parameter   :: eps = 1d-9

  double precision              :: denom

  if(n==0) then
    denom = qzl**2
    if(denom < eps) then
      a_n = 1d0
    else
      a_n = 2/denom * (1d0-cos(qzl))
    endif 
    return
  endif

  denom = (qzl**2 - (n*Pi)**2)**2
  if(denom < eps) then
    a_n = 0.5d0
  else
    a_n = (2*qzl)**2/denom * (1d0 - cos(qzl)*(-1d0)**n)
  endif

end function a_n


double precision function lorenzian(dqq, omega)
   implicit none
   double precision,  intent(in)  :: dqq
   double precision,  intent(in)  :: omega
   double precision, parameter    :: Pi = 4*atan(1d0)


   lorenzian   = dqq/Pi/(dqq**2+omega**2)

end function lorenzian

double precision function ftlor(dqq, t)
   implicit none
   double precision,  intent(in)  :: dqq
   double precision,  intent(in)  :: t

   ftlor       = exp(-abs(dqq*t))
    
end function ftlor


double precision function sqt(diff, qz, l, t)
   implicit none
   double precision, intent(in) :: diff
   double precision, intent(in) :: qz
   double precision, intent(in) :: l
   double precision, intent(in) :: t

   double precision, parameter    :: Pi = 4*atan(1d0)
   double precision :: a, eps = 1d-8
   integer :: i, imax = 200

   sqt  = 0 
   imax = nint(sqrt(2d0 * sqrt((qz*l)**2 * eps**3)) / (Pi*eps)) + 1
   do i=0,imax
     a = a_n(i,qz*l)
     sqt = sqt + a * ftlor(diff*(i*Pi/l)**2,t)
   enddo

end function sqt


double precision function sqt_ker(theta)
   implicit none 
   double precision, intent(in) ::  theta
   
   sqt_ker = sqt(diffusion, q * cos(theta), lz, tau) * sin(theta)


end function sqt_ker



end function pfg1d
