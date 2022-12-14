program testviscosity
implicit none
double precision               :: T, eta
double precision               :: toluene_viscosity
 double precision, parameter   :: factor_for_deuterated = 1.05d0

integer :: i

  write(6,*)" T/K           eta / mPa s     eta(d) / mPa s"    

  do i=240,400,1
    T   = i
    eta = toluene_viscosity(T) * 1d3
    write(6,'(f12.2,3x,2f12.5)') T, eta, eta*factor_for_deuterated
  enddo


end program testviscosity

function toluene_viscosity( temperature ) result(eta)
! Ref: Santos FJV et al., J. Phys. Chem. Ref. Data, Vol. 35, No. 1, 20

   implicit none
   double precision, intent(in)  :: temperature
   double precision              :: eta



   double precision              :: eta_star
   double precision              :: T_star
   double precision              :: lneta_star

   double precision, parameter   :: Tref   = 298.15d0
   double precision, parameter   :: etaref = 554.2d-6
   double precision, parameter   :: A      = -5.2203d0
   double precision, parameter   :: B      =  8.964d0
   double precision, parameter   :: C      = -5.834d0
   double precision, parameter   :: D      =  2.098d0

  

   T_star     = temperature / Tref
   
   lneta_star = A + B/T_star  + C/T_star**2  + D/T_star**3

   eta        = etaref * exp(lneta_star)
   
end function toluene_viscosity 
