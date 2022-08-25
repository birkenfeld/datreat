program dt

double precision :: t,q,Dlong,tstar,beta, sqt, diffusion_beta

q = 0.1
Dlong = 10
tstar = 10
beta  = 0.5

do i=0,500
t   = i*0.1
sqt = diffusion_beta(t,q,Dlong,tstar,beta)
enddo


end program dt


function diffusion_beta(t,q,Dlong,tstar,beta) result(sqt)
!-----------------------------------------------------------------------------
!
! diffusion scattering function under Gaussian approximation
! but with subdiffusive start
! 
! the physical parameters are the classic long time diffusion constant Dlong
!                         and a transition time tstar
!                         below tstar the diffusion is sub (super?) diffusive 
!                         with <r**2>(t) ~ t**beta
!
! the transition is smooth (stetig 1 x differenzierbar) at tstar
! this implies that the long time diffusion is linear with an offset, which
! implicitly is controlle by the value of tstar and the smoothness condition
!
! Michael Monkenbusch, JCNS-1, June 2015
!-----------------------------------------------------------------------------
implicit none
 

double precision, intent(in)  :: t       ! time in the time units used in the diffusion constant units
double precision, intent(in)  :: q       ! "q", momentum transfer in the length units used for Dlong
double precision, intent(in)  :: Dlong   ! long time diffusion constant in units consistent wit t and q
double precision, intent(in)  :: tstar   ! transition time between short and long time diffusion (units as t)
double precision, intent(in)  :: beta    ! streching exponent for the sublinear initial part of diffusion

double precision              :: sqt     ! return value = intermediate scattering factor of the diffusion

double precision              :: r2


if(t < tstar) then
   r2  = 6*t**beta*tstar**(-beta+1)*Dlong/beta
else
   r2  = 6*Dlong*(beta*t-beta*tstar+tstar)/beta
endif

sqt = exp(-q*q*r2/6d0) 

! testing
! write(6,'(3F18.9)') t, r2, sqt
! testing

end function diffusion_beta 
