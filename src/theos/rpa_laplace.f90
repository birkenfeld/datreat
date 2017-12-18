
!----------------------------------------------------------------------------------------------------------------
module rpa_laplace
!----------------------------------------------------------------------------------------------------------------

integer,          parameter :: XPREC = 8   !! TBD hier kind DBL einsetzen
real(kind=XPREC), parameter :: Pi = 4d0*atan(1d0)

real(kind=XPREC)            :: xS00
real(kind=XPREC)            :: xa1
real(kind=XPREC)            :: xa2
real(kind=XPREC)            :: xr1
real(kind=XPREC)            :: xr2
real(kind=XPREC)            :: xr3
real(kind=XPREC)            :: xScc00
real(kind=XPREC)            :: xb1
real(kind=XPREC)            :: xb2
real(kind=XPREC)            :: xg1
real(kind=XPREC)            :: xg2
real(kind=XPREC)            :: xg3
real(kind=XPREC)            :: xphi

real(kind=XPREC)            :: xil     = 0.01d0    ! distance of path for inv-laplace integration
real(kind=XPREC)            :: epap    = 1d-8      ! Apodisationfactor for numerical integrand
real(kind=XPREC)            :: epsrpa  = 1d-5      ! accuracy parameter



real   (kind=XPREC) :: limit_scale_factor=150d0   ! determines limit of integration
real   (kind=XPREC) :: rlow = 1d0/10000d0         ! lowest expected detectable rate




real(kind=XPREC)            :: t_param
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following variables are use by the 2D (ie. 3 component version)
  integer, private     :: is_sel      ! selector 1: S11, 2: S12, 3: S12, 4: S22 
  real   (kind=XPREC)  :: phi1        ! volume fraction of polymer component 1
  real   (kind=XPREC)  :: phi2        ! volume fraction of polymer component 2
  real   (kind=XPREC)  :: Scc00       ! unperturbed structure factor S(Q) of "matrix" polymers
  real   (kind=XPREC)  :: S0011       ! unperturbed structure factor S(Q) of polymer 1
  real   (kind=XPREC)  :: S0022       ! unperturbed structure factor S(Q) of polymer 2
  integer              :: nexpcc      ! number of exp-functions to describe background
  integer              :: nexp1       ! number of exp-functions to describe component1
  integer              :: nexp2       ! number of exp-functions to describe component2
  integer, parameter   :: mex = 32    ! maximum number of exponentials (typ 3-6 !)
  real   (kind=XPREC)  :: aexp_cc(mex)! amplitude coeffs for laplace-exp representation of "matrix"
  real   (kind=XPREC)  :: rexp_cc(mex)! rate      coeffs for laplace-exp representation of "matrix"
  real   (kind=XPREC)  :: aexp_s1(mex)! amplitude coeffs for laplace-exp representation of polymer 1
  real   (kind=XPREC)  :: rexp_s1(mex)! rate      coeffs for laplace-exp representation of polymer 1
  real   (kind=XPREC)  :: aexp_s2(mex)! amplitude coeffs for laplace-exp representation of polymer 2
  real   (kind=XPREC)  :: rexp_s2(mex)! rate      coeffs for laplace-exp representation of polymer 2
 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

 function Fn(s,n,a,r)
!--------------------
  implicit none
  complex(kind=XPREC)               :: Fn    ! Laplace transform of the n-exp timefunctio
  complex(kind=XPREC) , intent(in)  :: s     ! Laplace 'frequency' variable         
  integer             , intent(in)  :: n     ! number of single simple exponentials represented
  real   (kind=XPREC) , intent(in)  :: a(n)  ! amplitude coefficients  
  real   (kind=XPREC) , intent(in)  :: r(n)  ! rate coefficients

  integer :: i
  
  Fn = (0d0,0d0)
  do i = 1, n
    Fn = Fn +  a(i) / ( s + r(i) ) 
  enddo

! eventuell geht mittlerweile auch:  Fn = sum( a/(s+r) )
             
 end function Fn



  function Ss_kernel_2D(u) result(v)
 !----------------------------------
  implicit none 
 
  real   (kind=XPREC)              :: v           ! real part of integrand
  real   (kind=XPREC), intent(in)  :: u           ! integration path variable
 
  complex(kind=XPREC)  :: s           ! Laplace variable
  complex(kind=XPREC)  :: Fc          ! relaxation function "matrix" polymer (cc-component) 
  complex(kind=XPREC)  :: Fs1         ! relaxation function of polymer component 1
  complex(kind=XPREC)  :: Fs2         ! relaxation function of polymer component 2
  complex(kind=XPREC)  :: Ss11        ! S_11(q,s)
  complex(kind=XPREC)  :: Ss22        ! S_22(q,s)
  complex(kind=XPREC)  :: Ss12        ! S_12(q,s)
  
  complex(kind=XPREC)  :: Ssxx

! Reference Maple: RpaLaplace_corr_2D_3.mw  (23)-(27)

!! the following variables shoul be global for the module --> move to module section
!!  integer              :: is_sel      ! selector 1: S11, 2: S12, 3: S12, 4: S22 
!!  complex(kind=XPREC)  :: Ssxx        ! return variable
!!  real   (kind=XPREC)  :: phi1        ! volume fraction of polymer component 1
!!  real   (kind=XPREC)  :: phi2        ! volume fraction of polymer component 2
!!  real   (kind=XPREC)  :: Scc00       ! unperturbed structure factor S(Q) of "matrix" polymers
!!  real   (kind=XPREC)  :: S0011       ! unperturbed structure factor S(Q) of polymer 1
!!  real   (kind=XPREC)  :: S0022       ! unperturbed structure factor S(Q) of polymer 2
!!  integer              :: nexpcc      ! number of exp-functions to describe background
!!  integer              :: nexp1       ! number of exp-functions to describe component1
!!  integer              :: nexp2       ! number of exp-functions to describe component2
!!  integer, parameter   :: mex = 32    ! maximum number of exponentials (typ 3-6 !)
!!  real   (kind=XPREC)  :: aexp_cc(mex)! amplitude coeffs for laplace-exp representation of "matrix"
!!  real   (kind=XPREC)  :: rexp_cc(mex)! rate      coeffs for laplace-exp representation of "matrix"
!!  real   (kind=XPREC)  :: aexp_s1(mex)! amplitude coeffs for laplace-exp representation of polymer 1
!!  real   (kind=XPREC)  :: rexp_s1(mex)! rate      coeffs for laplace-exp representation of polymer 1
!!  real   (kind=XPREC)  :: aexp_s2(mex)! amplitude coeffs for laplace-exp representation of polymer 2
!!  real   (kind=XPREC)  :: rexp_s2(mex)! rate      coeffs for laplace-exp representation of polymer 2

!! note expternal parameter: t_param for the time in the resulting S(Q,t)
!! and the apodisation epsilon parameter epap
!! xil the parallel shift of the path
     
 

  s     = cmplx(xil,u)

 
  Fc   =  Fn (s, nexpcc , aexp_cc, rexp_cc )
  Fs1  =  Fn (s, nexp1  , aexp_s1, rexp_s1 )
  Fs2  =  Fn (s, nexp2  , aexp_s2, rexp_s2 )


  select case(is_sel) 

  case(1)
   Ss11 = &
     phi1*S0011*(-Fs1*(-1+phi1+phi2)**2*(Fc*s-1)*Scc00**2+(-1+phi1+phi2)*(((-2+(Fc+Fs2)*s)* &
     S0022*phi2+Fc*s*phi1*S0011)*Fs1-Fc*phi1*S0011)*Scc00-phi2*((phi2*S0022*(Fs2*s-1)+Fs2* &
     S0011*phi1*s)*Fs1-Fs2*phi1*S0011)*S0022)/((-1+phi1+phi2)*Scc00*(Fc*s-1)-Fs1*S0011*phi1*s+ &
     (-Fs2*S0022*s+S0022)*phi2+phi1*S0011)/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
   Ssxx = Ss11
   

   case(2:3) 
    Ss12 = &
     phi1*S0011*phi2*S0022*(-(-1+phi1+phi2)*(Fs1*Fs2*s+Fc-Fs1-Fs2)*Scc00+Fs2*S0011*(Fs1*s-1)* &
     phi1+Fs1*phi2*S0022*(Fs2*s-1))/((-1+phi1+phi2)*Scc00*(Fc*s-1)+(-Fs1*S0011*s+S0011)*phi1- &
     phi2*S0022*(Fs2*s-1))/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
    Ssxx = Ss12
    

   case(4)
    Ss22 = &
     (-Fs2*(-1+phi1+phi2)**2*(Fc*s-1)*Scc00**2+(((-2+(Fc+Fs1)*s)*S0011*phi1+phi2*S0022*Fc*s)* &
     Fs2-Fc*phi2*S0022)*(-1+phi1+phi2)*Scc00-phi1*S0011*((S0011*(Fs1*s-1)*phi1+Fs1*S0022*phi2*s)* &
     Fs2-Fs1*phi2*S0022))*phi2*S0022/((-1+phi1+phi2)*Scc00*(Fc*s-1)-Fs2*S0022*phi2*s+ &
     (-Fs1*S0011*s+S0011)*phi1+phi2*S0022)/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
    Ssxx = Ss22
   
   
   case default
    write(6,*)"Error in Ss_kernel, invalid selection (outside 1:4), is_sel=", is_sel
    stop   ! this indicates a programming error, thus we stop here -> fix the bug !

  end select

    Ssxx = Ssxx * exp(t_param * s) * exp(-epap*u*u)
    Ssxx = Ssxx / (2 * Pi ) 

    v = Realpart(Ssxx)
 

end function Ss_kernel_2D


 function St_rpa(t,i1,i2) result(yr)
!-----------------------------------
  use integration
  implicit none

  real   (kind=XPREC), intent(in) :: t 
  integer            , intent(in) :: i1
  integer            , intent(in) :: i2
  real   (kind=XPREC)             :: yr


  integer, parameter  :: idx(2,2) = reshape([1,2,3,4],[2,2])
  real   (kind=XPREC) :: limit_scale
  real   (kind=XPREC) :: rmin
  real   (kind=XPREC) :: a, b

  is_sel  = idx(i1,i2)

  call integral_setmethod(4)
  call integral_setaccuracy(1d-10)
  call intprint(0)
!
 
  rmin = max(  rlow,                            &
               min( minval(rexp_s1(1:nexp1)),   &
                    minval(rexp_s1(1:nexp2)),   &
                    minval(rexp_cc(1:nexpcc)) ))

  limit_scale = limit_scale_factor / rmin

  t_param = t

  b  =  max(limit_scale/t, limit_scale)
  a  = -b
! yr =   integral( Ss_kernel_2d ,a  ,b,10000,1d-8)
  yr = 2*integral( Ss_kernel_2d ,0d0,b,500000,epsrpa)

end function St_rpa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine PolyRootCoeffs2d(coeffs, S0011, a1, a2, r1, r2, r3,  &
                                    S0022, Scc00, b1, b2, g1, g2, g3, phi1, phi2)
  implicit none
   complex(kind=XPREC), intent(out) :: coeffs(0:6) ! Coefficients for the root polynome that results
                                                   ! from the inversel Laplace Transform
                                                   ! needed to determine the zero location that enter
                                                   ! the exp sum that yields the results of the invlaplace tr
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

 
   real   (kind=XPREC), intent(in) :: S0011    ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: S0022    ! S(q) value of the same as matrix but maybe with different label                 
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi1     ! volume fraction of sample polymer
   real   (kind=XPREC), intent(in) :: phi2     ! volume fraction of labelled matrix polymer

   complex(kind=XPREC)  :: C0, C1, C2, C3, C4, C5, C6

   real   (kind=XPREC) :: a3, b3
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)



C0=-r3*((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)*g3*r2*g2*r1*g1



C1=(((((-(-1+phi1+phi2)*(b3-1)*Scc00-phi1*S0011+phi2*S0022*(b3-1))*g2-((-1+phi1+phi2)* &
(b2-1)*Scc00+phi1*S0011-phi2*S0022*(b2-1))*g3)*g1-((b1-1)*(-1+phi1+phi2)*Scc00+phi1* &
S0011-phi2*S0022*(b1-1))*g3*g2)*r3+g1*((-1+phi1+phi2)*Scc00+S0011*(a3-1)*phi1-phi2*S0022)* &
g3*g2)*r2+((-1+phi1+phi2)*Scc00+S0011*(a2-1)*phi1-phi2*S0022)*g1*r3*g3*g2)*r1+g1*r3* &
((-1+phi1+phi2)*Scc00+S0011*(a1-1)*phi1-phi2*S0022)*g3*r2*g2


C2=((((-(b2+b3-1)*(-1+phi1+phi2)*Scc00-phi1*S0011+phi2*S0022*(b2+b3-1))*g1+(-(b1+b3-1)* &
(-1+phi1+phi2)*Scc00-phi1*S0011+phi2*S0022*(b1+b3-1))*g2-((b1+b2-1)*(-1+phi1+phi2)*Scc00+ &
phi1*S0011-S0022*(b1+b2-1)*phi2)*g3)*r3+((-(-1+phi1+phi2)*(b3-1)*Scc00+S0011*(a3-1)*phi1+ &
phi2*S0022*(b3-1))*g2+(-(-1+phi1+phi2)*(b2-1)*Scc00+S0011*(a3-1)*phi1+phi2*S0022*(b2-1))* &
g3)*g1+(-(b1-1)*(-1+phi1+phi2)*Scc00+S0011*(a3-1)*phi1+phi2*S0022*(b1-1))*g3*g2)*r2+(((- &
(-1+phi1+phi2)*(b3-1)*Scc00+S0011*(a2-1)*phi1+phi2*S0022*(b3-1))*g2+(-(-1+phi1+phi2)*(b2-1)* &
Scc00+S0011*(a2-1)*phi1+phi2*S0022*(b2-1))*g3)*g1+(-(b1-1)*(-1+phi1+phi2)*Scc00+S0011*(a2-1)* &
phi1+phi2*S0022*(b1-1))*g3*g2)*r3+g1*((-1+phi1+phi2)*Scc00+S0011*(a2+a3-1)*phi1-phi2*S0022)* &
g3*g2)*r1+((((-(-1+phi1+phi2)*(b3-1)*Scc00+S0011*(a1-1)*phi1+phi2*S0022*(b3-1))*g2+(-(-1+phi1+ &
phi2)*(b2-1)*Scc00+S0011*(a1-1)*phi1+phi2*S0022*(b2-1))*g3)*g1+(-(b1-1)*(-1+phi1+phi2)*Scc00+ &
S0011*(a1-1)*phi1+phi2*S0022*(b1-1))*g3*g2)*r3+g1*((-1+phi1+phi2)*Scc00+S0011*(a1+a3-1)*phi1- &
phi2*S0022)*g3*g2)*r2+g1*((-1+phi1+phi2)*Scc00+phi1*S0011*(a1+a2-1)-phi2*S0022)*r3*g3*g2



C3=-(-1+phi1+phi2)*((((b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r2+((b2+b3-1)* &
g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r3+((b3-1)*g2+(b2-1)*g3)*g1+g3*g2*(b1-1))*r1+(((b2+b3-1)*g1+ &
(b1+b3-1)*g2+(b1+b2-1)*g3)*r3+((b3-1)*g2+(b2-1)*g3)*g1+g3*g2*(b1-1))*r2+(((b3-1)*g2+(b2-1) &
*g3)*g1+g3*g2*(b1-1))*r3-g3*g1*g2)*Scc00+(((-phi1*S0011+S0022*(b1+b2+b3-1)*phi2)*r3+(S0011* &
(a3-1)*phi1+phi2*S0022*(b2+b3-1))*g1+(S0011*(a3-1)*phi1+phi2*S0022*(b1+b3-1))*g2+(S0011* &
(a3-1)*phi1+S0022*(b1+b2-1)*phi2)*g3)*r2+((S0011*(a2-1)*phi1+phi2*S0022*(b2+b3-1))*g1+ &
(S0011*(a2-1)*phi1+phi2*S0022*(b1+b3-1))*g2+(S0011*(a2-1)*phi1+S0022*(b1+b2-1)*phi2)*g3)* &
r3+((S0011*(a2+a3-1)*phi1+phi2*S0022*(b3-1))*g2+(S0011*(a2+a3-1)*phi1+phi2*S0022*(b2-1))* &
g3)*g1+(S0011*(a2+a3-1)*phi1+phi2*S0022*(b1-1))*g3*g2)*r1+(((S0011*(a1-1)*phi1+phi2*S0022* &
(b2+b3-1))*g1+(S0011*(a1-1)*phi1+phi2*S0022*(b1+b3-1))*g2+(S0011*(a1-1)*phi1+S0022*(b1+b2-1)* &
phi2)*g3)*r3+((S0011*(a1+a3-1)*phi1+phi2*S0022*(b3-1))*g2+(S0011*(a1+a3-1)*phi1+phi2*S0022* &
(b2-1))*g3)*g1+(S0011*(a1+a3-1)*phi1+phi2*S0022*(b1-1))*g3*g2)*r2+(((phi1*S0011*(a1+a2-1)+ &
phi2*S0022*(b3-1))*g2+(phi1*S0011*(a1+a2-1)+phi2*S0022*(b2-1))*g3)*g1+(phi1*S0011*(a1+a2-1)+ &
phi2*S0022*(b1-1))*g3*g2)*r3+(phi1*S0011*(a1+a2+a3-1)-phi2*S0022)*g1*g3*g2


C4=-(((b1+b2+b3-1)*r2+(b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r1+((b1+b2+b3-1) &
*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r2+((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r3+ &
((b3-1)*g2+(b2-1)*g3)*g1+g3*g2*(b1-1))*(-1+phi1+phi2)*Scc00+(((a3-1)*r2+r3*(a2-1)+(a2+a3-1)* &
(g3+g1+g2))*r1+((a1-1)*r3+(a1+a3-1)*(g3+g1+g2))*r2+(a1+a2-1)*(g3+g1+g2)*r3+(a1+a2+a3-1)* &
((g3+g2)*g1+g3*g2))*S0011*phi1+(((b1+b2+b3-1)*r2+(b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+ &
(b1+b2-1)*g3)*r1+((b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r2+((b2+b3-1)*g1+ &
(b1+b3-1)*g2+(b1+b2-1)*g3)*r3+((b3-1)*g2+(b2-1)*g3)*g1+g3*g2*(b1-1))*S0022*phi2



C5=-(-1+phi1+phi2)*((r1+r2+r3+g2+g3)*b1+(r1+r2+r3+g1+g3)*b2+(r1+r2+r3+g1+g2)* &
b3-r1-r2-r3-g1-g2-g3)*Scc00+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)+(a1+a2+a3-1)* &
(g3+g1+g2))*S0011*phi1+phi2*S0022*((r1+r2+r3+g2+g3)*b1+(r1+r2+r3+g1+g3)*b2+ &
(r1+r2+r3+g1+g2)*b3-r1-r2-r3-g1-g2-g3)


C6=-(b1+b2+b3-1)*(-1+phi1+phi2)*Scc00+S0022*(b1+b2+b3-1)*phi2+phi1*S0011*(a1+a2+a3-1)


coeffs=[C0,C1,C2,C3,C4,C5,C6]

end subroutine PolyRootCoeffs2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function InvLaplace2d11(t, S0011, a1, a2, r1, r2, r3, S0022, Scc00, b1, b2, g1, g2, g3, phi1, phi2 ) result (St)
!---------------------------------------------------------------------------------------------------------------

  implicit none
   complex(kind=XPREC)  :: St                  ! S(Q,t)  (Q only implicit)
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   real   (kind=XPREC), intent(in) :: t 
   real   (kind=XPREC), intent(in) :: S0011    ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: S0022    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi1     ! volume fraction of sample polymer
   real   (kind=XPREC), intent(in) :: phi2     ! volume fraction of sample polymer
 
  integer, parameter   :: maxdeg = 6
  integer              :: degree = maxdeg
  complex(kind=XPREC)  :: coeffs(0:maxdeg), pzero(1:maxdeg), Z 
  real   (kind=XPREC)  :: Prefak
  integer              :: i

  real   (kind=XPREC), parameter :: epsilon = 1d-8
  real   (kind=XPREC) :: a3, b3
   
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)

! write(6,*)'enter invaplace..'

  call PolyRootCoeffs2d(coeffs, S0011, a1,a2,r1,r2,r3, S0022, Scc00, b1,b2,g1,g2,g3, phi1,phi2)

!! if the above conditions for a3 and b3 are enforceed (i.e. a1+a2+a3=1, ...)
!! the the highest coefficients are identical to zero, which is not well received
!! by the roots search, therefore we check here for this and in case reduce the
!! degree of the polynomial
dl:  do i=degree,1,-1
       if(abs(coeffs(i)) > epsilon * sum(abs(coeffs(0:degree)))) then
         degree = i
         exit dl
       endif
     enddo dl

  call cmplx_roots_gen(pzero, coeffs, degree, .true. , .false. )

! write(6,*)'roots: ',pzero

  call check_degeneracy(degree, pzero)



prefak=-phi1*S0011*((-1+phi1+phi2)*Scc00-phi2*S0022)/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)


  St = 0

  do i = 1, degree
   Z = pzero(i)
   St= St + &
exp(Z*t)*(((b1+b2+b3-1)*(a1+a2+a3)*(-1+phi1+phi2)*Scc00-phi1*S0011*(b1+b2+b3)*(a1+a2+a3-1)                 & 
-S0022*(b1+b2+b3-1)*(a1+a2+a3)*phi2)*Z**5+(((b2+b3-1)*(a1+a2+a3)*g1+(b1+b3-1)*(a1+a2+a3)*g2+(b1+b2-1)*            &
(a1+a2+a3)*g3+(b1+b2+b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*(-1+phi1+phi2)*Scc00-((b2+b3)*(a1+a2+a3-1)*        &
g1+(b1+b3)*(a1+a2+a3-1)*g2+(b1+b2)*(a1+a2+a3-1)*g3+(b1+b2+b3)*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)))*          &
phi1*S0011-((b2+b3-1)*(a1+a2+a3)*g1+(b1+b3-1)*(a1+a2+a3)*g2+(b1+b2-1)*(a1+a2+a3)*g3+(b1+b2+b3-1)*((a2+a3)*        &
r1+(a1+a3)*r2+r3*(a1+a2)))*S0022*phi2)*Z**4+((-1+phi1+phi2)*(((b3-1)*(a1+a2+a3)*g2+(b2-1)*(a1+a2+a3)*             &
g3+(b2+b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g1+((b1-1)*(a1+a2+a3)*g3+(b1+b3-1)*((a2+a3)*r1+(a1+a3)*r2+       &
r3*(a1+a2)))*g2+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b2-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b2+b3-1))*        &
Scc00-phi1*((b3*(a1+a2+a3-1)*g2+b2*(a1+a2+a3-1)*g3+(b2+b3)*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)))*g1+          &
(b1*(a1+a2+a3-1)*g3+(b1+b3)*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)))*g2+(b1+b2)*((a2+a3-1)*r1+(a1+a3-1)*         &
r2+r3*(a1+a2-1))*g3+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*(b1+b2+b3))*S0011-S0022*(((b3-1)*(a1+a2+a3)*          &
g2+(b2-1)*(a1+a2+a3)*g3+(b2+b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g1+((b1-1)*(a1+a2+a3)*g3+(b1+b3-1)*         &
((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g2+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b2-1)*g3+((a2*r3+a3*r2)*r1+         &
a1*r3*r2)*(b1+b2+b3-1))*phi2)*Z**3+(((((-a1-a2-a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b3-1))*g2+              &
((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2-1)*g3+(b2+b3-1)*((a2*r3+a3*r2)*r1+a1*r3*r2))*g1+((b1-1)*((a2+a3)*          &
r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b3-1))*g2+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*            &
(b1+b2-1))*(-1+phi1+phi2)*Scc00-((b3*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)          &
*r2+r3*(a1+a2-1))*b2*g3+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*(b2+b3))*g1+(b1*((a2+a3-1)*r1+(a1+a3-1)*          &
r2+r3*(a1+a2-1))*g3+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*(b1+b3))*g2+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*          &
(a1-1))*(b1+b2)*g3-r3*r2*r1*(b1+b2+b3))                                                                           &
*phi1*S0011-((((-a1-a2-a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b3-1))*g2+((a2+a3)*r1+(a1+a3)*r2+r3*            &
(a1+a2))*(b2-1)*g3+(b2+b3-1)*((a2*r3+a3*r2)*r1+a1*r3*r2))*g1+((b1-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))           &
*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b3-1))*g2+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b1+b2-1))*S0022*phi2)*              &
Z**2+((((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b3-1))*g2+((a2*r3+a3*r2)            &
*r1+a1*r3*r2)*g3*(b2-1))*g1+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b1-1)*g2)*(-1+phi1+phi2)*Scc00+((-b3*(((a3-1)         &
*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g2-(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*b2*g3+r3*r2*r1*(b2+b3))*g1+            &
(-b1*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3+r3*r2*r1*(b1+b3))*g2+g3*r1*r2*r3*(b1+b2))*phi1*S0011-             &
(((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b3-1))*g2+((a2*r3+a3*r2)*r1+a1*r3*r2)     &
*g3*(b2-1))*g1+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b1-1)*g2)*S0022*phi2)*Z-((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*g2          &
*(-1+phi1+phi2)*g1*Scc00+r3*r2*phi1*((b2*g3+b3*g2)*g1+g3*b1*g2)*S0011*r1+((a2*r3+a3*r2)*r1+a1*r3*r2)*S0022*g3     &
*g2*g1*phi2)/((6*(b1+b2+b3-1)*(-1+phi1+phi2)*Scc00-6*phi1*S0011*(a1+a2+a3-1)-6*S0022*(b1+b2+b3-1)*phi2)*          &
Z**5+(5*((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*(r1+r2+r3))*(-1+phi1+phi2)*Scc00-5*                  &
((a1+a2+a3-1)*g1+(a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*phi1*S0011-5*            &
((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*(r1+r2+r3))*S0022*phi2)*Z**4+(4*(((b3-1)*g2+                 &
(b2-1)*g3+(b2+b3-1)*(r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*(r1+r2+r3))*g2+(b1+b2-1)*                                 &
(r1+r2+r3)*g3+((r2+r3)*r1+r2*r3)*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00-4*(((a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+          &
(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g1+((a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+           &
((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*phi1*S0011-4*S0022*            &
(((b3-1)*g2+(b2-1)*g3+(b2+b3-1)*(r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*(r1+r2+r3))*g2+(b1+b2-1)*(r1+r2+r3)           &
*g3+((r2+r3)*r1+r2*r3)*(b1+b2+b3-1))*phi2)*Z**3+(3*(((-g3+(b3-1)*r1+(b3-1)*r2+r3*(b3-1))*g2+(b2-1)*               &
(r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((b1-1)*(r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)*r1+r2*r3))*g2+            &
((r2+r3)*r1+r2*r3)*(b1+b2-1)*g3+r3*r2*r1*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00-3*phi1*((((a1+a2+a3-1)*g3+            &
(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*             &
(a2-1))*r1+r2*r3*(a1-1))*g1+(((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*          &
(a1-1))*g2+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3-r3*r2*r1)*S0011-3*(((-g3+(b3-1)*r1+(b3-1)*r2+r3*            &
(b3-1))*g2+(b2-1)*(r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((b1-1)*(r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)         &
*r1+r2*r3))*g2+((r2+r3)*r1+r2*r3)*(b1+b2-1)*g3+r3*r2*r1*(b1+b2+b3-1))*S0022*phi2)*Z**2+(2*((((-r1-r2-r3)          &
*g3+((r2+r3)*r1+r2*r3)*(b3-1))*g2+((r2+r3)*r1+r2*r3)*(b2-1)*g3+r3*r2*r1*(b2+b3-1))*g1+(((r2+r3)*r1+r2*r3)*        &
(b1-1)*g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*(-1+phi1+phi2)*Scc00+2*(((((-a2-a3+1)*r1+(-a1-a3+1)       &
*r2-r3*(a1+a2-1))*g3+((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g2+(((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*         &
g3+r3*r2*r1)*g1+((((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g2+g3*r1*r2*r3)*phi1*S0011-2*S0022*        &
((((-r1-r2-r3)*g3+((r2+r3)*r1+r2*r3)*(b3-1))*g2+((r2+r3)*r1+r2*r3)*(b2-1)*g3+r3*r2*r1*(b2+b3-1))*g1+              &
(((r2+r3)*r1+r2*r3)*(b1-1)*g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*phi2)*Z+(((((-r2-r3)*r1-              &
r2*r3)*g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*(b2-1))*g1+g3*r1*r2*r3*g2*(b1-1))*(-1+phi1+phi2)*Scc00+phi1*            &
((((((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g2+g3*r1*r2*r3)*g1+g3*r1*r2*r3*g2)*S0011-S0022*          &
(((((-r2-r3)*r1-r2*r3)*g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*(b2-1))*g1+g3*r1*r2*r3*g2*(b1-1))*phi2)

 enddo

   St = St * Prefak 

end function InvLaplace2d11


function InvLaplace2d12(t, S0011, a1, a2, r1, r2, r3, S0022, Scc00, b1, b2, g1, g2, g3, phi1, phi2 ) result (St)
!---------------------------------------------------------------------------------------------------------------

  implicit none
   complex(kind=XPREC)  :: St                  ! S(Q,t)  (Q only implicit)
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   real   (kind=XPREC), intent(in) :: t 
   real   (kind=XPREC), intent(in) :: S0011    ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: S0022    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi1     ! volume fraction of sample polymer
   real   (kind=XPREC), intent(in) :: phi2     ! volume fraction of sample polymer
 
  integer, parameter   :: maxdeg = 6
  integer              :: degree = maxdeg
  complex(kind=XPREC)  :: coeffs(0:maxdeg), pzero(1:maxdeg), Z 
  real   (kind=XPREC)  :: Prefak
  integer              :: i

  real   (kind=XPREC), parameter :: epsilon = 1d-8
  real   (kind=XPREC) :: a3, b3
   
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)

! write(6,*)'enter invaplace..'

  call PolyRootCoeffs2d(coeffs, S0011, a1,a2,r1,r2,r3, S0022, Scc00, b1,b2,g1,g2,g3, phi1,phi2)

!! if the above conditions for a3 and b3 are enforceed (i.e. a1+a2+a3=1, ...)
!! the the highest coefficients are identical to zero, which is not well received
!! by the roots search, therefore we check here for this and in case reduce the
!! degree of the polynomial
dl:  do i=degree,1,-1
       if(abs(coeffs(i)) > epsilon * sum(abs(coeffs(0:degree)))) then
         degree = i
         exit dl
       endif
     enddo dl

  call cmplx_roots_gen(pzero, coeffs, degree, .true. , .false. )

! write(6,*)'roots: ',pzero

  call check_degeneracy(degree, pzero)



Prefak=-phi1*S0011*phi2*S0022/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)


  St = 0

  do i = 1, degree
   Z = pzero(i)
   St= St + &
      (((b1+b2+b3-1)*(a1+a2+a3)*(-1+phi1+phi2)*Scc00-phi1*S0011*(b1+b2+b3)*(a1+a2+a3-1)&
        -S0022*(b1+b2+b3-1)*(a1+a2+a3)*phi2)*Z**5+(((b2+b3-1)*(a1+a2+a3)*g1+(b1+b3-1)*&
        (a1+a2+a3)*g2+(b1+b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b2+b3-1))&
       *(-1+phi1+phi2)*Scc00-phi1*S0011*((b2+b3)*(a1+a2+a3-1)*g1+(b1+b3)*(a1+a2+a3-1)*g2+(b1+b2)&
       *(a1+a2+a3-1)*g3+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b2+b3))-S0022*((b2+b3-1)*&
       (a1+a2+a3)*g1+(b1+b3-1)*(a1+a2+a3)*g2+(b1+b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*&
       (a1+a2))*(b1+b2+b3-1))*phi2)*Z**4+((((b3-1)*(a1+a2+a3)*g2+(b2-1)*(a1+a2+a3)*g3+&
       ((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2+b3-1))*g1+((b1-1)*(a1+a2+a3)*g3+((a2+a3)*r1+&
       (a1+a3)*r2+r3*(a1+a2))*(b1+b3-1))*g2+(b1+b2-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+&
       ((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00-((b3*(a1+a2+a3-1)*g2+&
       b2*(a1+a2+a3-1)*g3+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b2+b3))*g1+(b1*(a1+a2+a3-1)&
       *g3+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b3))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*&
       (a1+a2-1))*(b1+b2)*g3+(b1+b2+b3)*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1)))*phi1*S0011-&
       S0022*(((b3-1)*(a1+a2+a3)*g2+(b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*&
       (b2+b3-1))*g1+((b1-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b3-1))*g2+&
       (b1+b2-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*&
       (b1+b2+b3-1))*phi2)*Z**3+(((((-a1-a2-a3)*g3+(b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))&
       *g2+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b2+b3-1))*g1+&
       (((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b3-1))*g2+&
       (b1+b2-1)*((a2*r3+a3*r2)*r1+a1*r3*r2)*g3)*(-1+phi1+phi2)*Scc00-phi1*((((a2+a3-1)*r1+&
       (a1+a3-1)*r2+r3*(a1+a2-1))*b3*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*b2*g3+(b2+b3)*&
       (((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1)))*g1+(((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*b1*&
       g3+(b1+b3)*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1)))*g2+(b1+b2)*(((a3-1)*r2+r3*(a2-1))*r1+&
       r2*r3*(a1-1))*g3-r3*r2*r1*(b1+b2+b3))*S0011-S0022*((((-a1-a2-a3)*g3+(b3-1)*((a2+a3)*r1+&
       (a1+a3)*r2+r3*(a1+a2)))*g2+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2-1)*g3+((a2*r3+a3*r2)*&
       r1+a1*r3*r2)*(b2+b3-1))*g1+(((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1-1)*g3+((a2*r3+a3*r2)*&
       r1+a1*r3*r2)*(b1+b3-1))*g2+(b1+b2-1)*((a2*r3+a3*r2)*r1+a1*r3*r2)*g3)*phi2)*Z**2+&
       ((((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+(b3-1)*((a2*r3+a3*r2)*r1+a1*r3*r2))*g2+&
       ((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b2-1))*g1+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*g2*(b1-1))*&
       (-1+phi1+phi2)*Scc00+((-b3*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g2-b2*&
       (((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3+r3*r2*r1*(b2+b3))*g1+(-b1*(((a3-1)*&
       r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3+r3*r2*r1*(b1+b3))*g2+g3*r1*r2*r3*(b1+b2))*phi1*&
       S0011-S0022*(((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+(b3-1)*((a2*r3+a3*r2)*r1+a1*r3*r2))&
       *g2+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b2-1))*g1+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*g2*(b1-1))*phi2)&
       *Z-((a2*r3+a3*r2)*r1+a1*r3*r2)*g1* &
       g3*g2*(-1+phi1+phi2)*Scc00+r3*phi1*((b2*g3+b3*g2)*g1+g3*b1*g2)*r2*r1*S0011+S0022*&
       ((a2*r3+a3*r2)*r1+a1*r3*r2)*g1*g3*g2*phi2)*exp(Z*t)/((6*(b1+b2+b3-1)*&
       (-1+phi1+phi2)*Scc00-6*phi1*S0011*(a1+a2+a3-1)-6*S0022*(b1+b2+b3-1)*phi2)*Z**5+&
       (5*((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*(r1+r2+r3))*(-1+phi1+phi2)*Scc00-&
       5*((a1+a2+a3-1)*g1+(a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*&
       phi1*S0011-5*((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*(r1+r2+r3))*S0022*phi2)&
       *Z**4+(4*(((b3-1)*g2+(b2-1)*g3+(b2+b3-1)*(r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*(r1+r2+r3))&
       *g2+(b1+b2-1)*(r1+r2+r3)*g3+(b1+b2+b3-1)*((r2+r3)*r1+r2*r3))*(-1+phi1+phi2)*Scc00-&
       4*(((a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g1+((a1+a2+a3-1)&
       *g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+&
       ((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*phi1*S0011-4*S0022*(((b3-1)*g2+(b2-1)*g3+(b2+b3-1)*&
       (r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*(r1+r2+r3))*g2+(b1+b2-1)*(r1+r2+r3)*g3+(b1+b2+b3-1)*&
       ((r2+r3)*r1+r2*r3))*phi2)*Z**3+(3*(((-g3+(b3-1)*r1+(b3-1)*r2+r3*(b3-1))*g2+(b2-1)*&
       (r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((b1-1)*(r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)*&
       r1+r2*r3))*g2+(b1+b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00&
       -3*((((a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*&
       r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g1+(((a2+a3-1)*r1+(a1+a3-1)*r2+&
       r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g2+(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*&
       (a1-1))*g3-r3*r2*r1)*phi1*S0011-3*S0022*(((-g3+(b3-1)*r1+(b3-1)*r2+r3*(b3-1))*g2+(b2-1)*&
       (r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((b1-1)*&
       (r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)*r1+r2*r3))*g2+(b1+b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*&
       (b1+b2+b3-1))*phi2)*Z**2+(2*(-1+phi1+phi2)*((((-r1-r2-r3)*g3+(b3-1)*((r2+r3)*r1+r2*r3))&
       *g2+(b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b2+b3-1))*g1+((b1-1)*((r2+r3)*r1+r2*r3)*&
       g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*Scc00+2*phi1*(((((-a2-a3+1)*r1+&
       (-a1-a3+1)*r2-r3*(a1+a2-1))*g3+((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g2+(((-a3+1)*&
       r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g1+((((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*&
       (a1-1))*g3+r3*r2*r1)*g2+g3*r1*r2*r3)*S0011-2*S0022*((((-r1-r2-r3)*g3+(b3-1)*&
       ((r2+r3)*r1+r2*r3))*g2+(b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b2+b3-1))*g1+((b1-1)*&
       ((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*phi2)*Z+&
       (((((-r2-r3)*r1-r2*r3)*g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*(b2-1))*g1+g3*r1*r2*r3*g2*&
       (b1-1))*(-1+phi1+phi2)*Scc00+phi1*((((((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*&
       g3+r3*r2*r1)*g2+g3*r1*r2*r3)*g1+g3*r1*r2*r3*g2)*S0011-S0022*(((((-r2-r3)*r1-r2*r3)&
       *g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*(b2-1))*g1+g3*r1*r2*r3*g2*(b1-1))*phi2)
       
  enddo

   St = St * Prefak 

end function InvLaplace2d12


function InvLaplace2d22(t, S0011, a1, a2, r1, r2, r3, S0022, Scc00, b1, b2, g1, g2, g3, phi1, phi2 ) result (St)
!---------------------------------------------------------------------------------------------------------------

  implicit none
   complex(kind=XPREC)  :: St                  ! S(Q,t)  (Q only implicit)
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   real   (kind=XPREC), intent(in) :: t 
   real   (kind=XPREC), intent(in) :: S0011    ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: S0022    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi1     ! volume fraction of sample polymer
   real   (kind=XPREC), intent(in) :: phi2     ! volume fraction of sample polymer
 
  integer, parameter   :: maxdeg = 6
  integer              :: degree = maxdeg
  complex(kind=XPREC)  :: coeffs(0:maxdeg), pzero(1:maxdeg), Z 
  real   (kind=XPREC)  :: Prefak
  integer              :: i

  real   (kind=XPREC), parameter :: epsilon = 1d-8
  real   (kind=XPREC) :: a3, b3, aplus
   
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)

! write(6,*)'enter invaplace..'

  call PolyRootCoeffs2d(coeffs, S0011, a1,a2,r1,r2,r3, S0022, Scc00, b1,b2,g1,g2,g3, phi1,phi2)

!! if the above conditions for a3 and b3 are enforceed (i.e. a1+a2+a3=1, ...)
!! the the highest coefficients are identical to zero, which is not well received
!! by the roots search, therefore we check here for this and in case reduce the
!! degree of the polynomial
dl:  do i=degree,1,-1
       if(abs(coeffs(i)) > epsilon * sum(abs(coeffs(0:degree)))) then
         degree = i
         exit dl
       endif
     enddo dl

  call cmplx_roots_gen(pzero, coeffs, degree, .true. , .false. )

! write(6,*)'roots: ',pzero

  call check_degeneracy(degree, pzero)



  prefak = -phi1*S0022**2*phi2**2*S0011/(((-1+phi1+phi2)*Scc00-phi2*S0022)* &
           ((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022))

  St = 0

  do i = 1, degree
   Z = pzero(i)
   St= St + &
       (((b1+b2+b3-1)*(a1+a2+a3)*(-1+phi1+phi2)*Scc00-phi1*S0011*(b1+b2+b3)*(a1+a2+a3-1)-S0022 &
     *(b1+b2+b3-1)*(a1+a2+a3)*phi2)*Z**5+(((b2+b3-1)*(a1+a2+a3)*g1+(b1+b3-1)*(a1+a2+a3)*g2+&
     (b1+b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00&
     -phi1*S0011*((b2+b3)*(a1+a2+a3-1)*g1+(b1+b3)*(a1+a2+a3-1)*g2+(b1+b2)*(a1+a2+a3-1)*g3+((a2+a3-1)&
     *r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b2+b3))-S0022 &
     *((b2+b3-1)*(a1+a2+a3)*g1+(b1+b3-1)*(a1+a2+a3)*g2+(b1+b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*&
     r2+r3*(a1+a2))*(b1+b2+b3-1))*phi2)*Z**4+((((b3-1)*(a1+a2+a3)*g2+(b2-1)*(a1+a2+a3)*g3+&
     ((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2+b3-1))*g1+((b1-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2&
     +r3*(a1+a2))*(b1+b3-1))*g2+(b1+b2-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+&
     a1*r3*r2)*(b1+b2+b3-1))*(-1+phi1+phi2)*Scc00-((b3*(a1+a2+a3-1)*g2+b2*(a1+a2+a3-1)*g3+((a2+a3-1)&
     *r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b2+b3))*g1+(b1*(a1+a2+a3-1)*g3+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*&
     (a1+a2-1))*(b1+b3))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b2)*g3+(b1+b2+b3)*(((a3-1)&
     *r2+r3*(a2-1))*r1+r2*r3*(a1-1)))*phi1*S0011-S0022 &
     *(((b3-1)*(a1+a2+a3)*g2+(b2-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b2+b3-1))*g1+&
     ((b1-1)*(a1+a2+a3)*g3+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*(b1+b3-1))*g2+(b1+b2-1)*((a2+a3)*r1+&
     (a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b2+b3-1))*phi2)*Z**3+(((((-a1-&
     a2-a3)*g3+(b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g2+((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*&
     (b2-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b2+b3-1))*g1+(((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*&
     (b1-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b3-1))*g2+(b1+b2-1)*((a2*r3+a3*r2)*r1+a1*r3*r2)&
     *g3)*(-1+phi1+phi2)*Scc00-phi1*((((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*b3*g2+((a2+a3-1)&
     *r1+(a1+a3-1)*r2+r3*(a1+a2-1))*b2*g3+(b2+b3)*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1)))*g1+&
     (((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*b1*g3+(b1+b3)*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*&
     (a1-1)))*g2+(b1+b2)*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3-r3*r2*r1*(b1+b2+b3))*S0011-S0022 &
     *((((-a1-a2-a3)*g3+(b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g2+((a2+a3)*r1+(a1+a3)*r2+&
     r3*(a1+a2))*(b2-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b2+b3-1))*g1+(((a2+a3)*r1+(a1+a3)*r2+r3*&
     (a1+a2))*(b1-1)*g3+((a2*r3+a3*r2)*r1+a1*r3*r2)*(b1+b3-1))*g2+(b1+b2-1)*((a2*r3+a3*r2)*&
     r1+a1*r3*r2)*g3)*phi2)*Z**2+((((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+(b3-1)*&
     ((a2*r3+a3*r2)*r1+a1*r3*r2))*g2+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b2-1))*g1+((a2*r3+a3*r2)&
     *r1+a1*r3*r2)*g3*g2*(b1-1))*(-1+phi1+phi2)*Scc00+((-b3*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*&
     (a1-1))*g2-b2*(((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3+r3*r2*r1*(b2+b3))*g1+(-b1*(((a3-1)&
     *r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3+r3*r2*r1*(b1+b3))*g2+g3*r1*r2*r3*(b1+b2))*phi1*S0011-S0022 &
     *(((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))*g3+(b3-1)*((a2*r3+a3*r2)*r1+a1*r3*r2))*g2+&
     ((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*(b2-1))*g1+((a2*r3+a3*r2)*r1+a1*r3*r2)*g3*g2*(b1-1))*phi2)*&
     Z-((a2*r3+a3*r2)*r1+a1*r3*r2)*g1*g3*g2*(-1+phi1+phi2)*Scc00+r3*phi1*((b2*g3+b3*g2)*&
     g1+g3*b1*g2)*r2*r1*S0011+S0022 &
     *((a2*r3+a3*r2)*r1+a1*r3*r2)*g1*g3*g2*phi2)*exp(Z*t)/((6*(b1+b2+b3-1)*&
     (-1+phi1+phi2)*Scc00-6*phi1*S0011*(a1+a2+a3-1)-6*S0022 &
     *(b1+b2+b3-1)*phi2)*Z**5+(5*((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*&
     (r1+r2+r3))*(-1+phi1+phi2)*Scc00-5*((a1+a2+a3-1)*g1+(a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)&
     *r1+(a1+a3-1)*r2+r3*(a1+a2-1))*phi1*S0011-5*((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)&
     *(r1+r2+r3))*S0022 &
     *phi2)*Z**4+(4*(((b3-1)*g2+(b2-1)*g3+(b2+b3-1)*(r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*&
     (r1+r2+r3))*g2+(b1+b2-1)*(r1+r2+r3)*g3+(b1+b2+b3-1)*((r2+r3)*r1+r2*r3))*(-1+phi1+phi2)*&
     Scc00-4*(((a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g1+&
     ((a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+&
     r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*phi1*S0011-4*S0022 &
     *(((b3-1)*g2+(b2-1)*g3+(b2+b3-1)*(r1+r2+r3))*g1+((b1-1)*g3+(b1+b3-1)*(r1+r2+r3))*g2+(b1+b2-1)&
     *(r1+r2+r3)*g3+(b1+b2+b3-1)*((r2+r3)*r1+r2*r3))*phi2)*Z**3+(3*(((-g3+(b3-1)*r1+(b3-1)*&
     r2+r3*(b3-1))*g2+(b2-1)*(r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((b1-1)*(r1+r2+r3)*g3+&
     (b1+b3-1)*((r2+r3)*r1+r2*r3))*g2+(b1+b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b1+b2+b3-1))*&
     (-1+phi1+phi2)*Scc00-3*((((a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+&
     ((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g1+&
     (((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g2+&
     (((a3-1)*r2+r3*(a2-1))*r1+r2*r3*(a1-1))*g3-r3*r2*r1)*phi1*S0011-3*S0022 &
     *(((-g3+(b3-1)*r1+(b3-1)*r2+r3*(b3-1))*g2+(b2-1)*(r1+r2+r3)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))&
     *g1+((b1-1)*(r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)*r1+r2*r3))*g2+(b1+b2-1)*((r2+r3)*r1+r2*r3)*&
     g3+r3*r2*r1*(b1+b2+b3-1))*phi2)*Z**2+(2*(-1+phi1+phi2)*((((-r1-r2-r3)*g3+(b3-1)*&
     ((r2+r3)*r1+r2*r3))*g2+(b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b2+b3-1))*g1+((b1-1)*&
     ((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*Scc00+2*phi1*&
     (((((-a2-a3+1)*r1+(-a1-a3+1)*r2-r3*(a1+a2-1))*g3+((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*&
     (a1-1))*g2+(((-a3+1)*r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g1+((((-a3+1)*&
     r2-r3*(a2-1))*r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g2+g3*r1*r2*r3)*S0011-2*S0022 &
     *((((-r1-r2-r3)*g3+(b3-1)*((r2+r3)*r1+r2*r3))*g2+(b2-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*&
     r1*(b2+b3-1))*g1+((b1-1)*((r2+r3)*r1+r2*r3)*g3+r3*r2*r1*(b1+b3-1))*g2+g3*r1*r2*r3*&
     (b1+b2-1))*phi2)*Z+(((((-r2-r3)*r1-r2*r3)*g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*&
     (b2-1))*g1+g3*r1*r2*r3*g2*(b1-1))*(-1+phi1+phi2)*Scc00+phi1*((((((-a3+1)*r2-r3*(a2-1))&
     *r1-r2*r3*(a1-1))*g3+r3*r2*r1)*g2+g3*r1*r2*r3)*g1+g3*r1*r2*r3*g2)*S0011-S0022 &
     *(((((-r2-r3)*r1-r2*r3)*g3+r3*r2*r1*(b3-1))*g2+g3*r1*r2*r3*(b2-1))*g1+g3*r1*r2*r3*g2*(b1-1))*phi2)
     
  enddo
   
   aplus =-((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)*&
           (b3*exp(-g3*t)+b1*exp(-g1*t)+b2*exp(-g2*t))*Scc00*(-1+phi1+phi2)
   St = St + aplus/(S0011*phi1*S0022*phi2)
   St = St * Prefak 

end function InvLaplace2d22





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





 function Ss4(s, S00, a1, a2, r1, r2, r3, Scc00, b1, b2, g1, g2, g3, phi) result(Ss)
!---------===-----------------------------------------------------------------------
   implicit none
   complex(kind=XPREC)             :: Ss       ! Result is the Laplace transform of the RPA expression
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   complex(kind=XPREC), intent(in) :: s        ! Laplace "frequency" variable
   real   (kind=XPREC), intent(in) :: S00      ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi      ! volume fraction of sample polymer


   complex(kind=XPREC)             :: Fs, Fc   ! Laplace transform of time functions (later possibly external fkts)


   real   (kind=XPREC) :: a3, b3
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)



   Fs = a1/(s+r1) + a2/(s+r2) + a3/(s+r3)
   Fc = b1/(s+g1) + b2/(s+b2) + b3/(s+g3)

   Ss = 1d0/(1d0/((1d0/Fs-s)*phi*S00*(-1d0/((1d0-phi)*Scc00*Fc*s-(1d0-phi)*Scc00) - &
        1d0/(1d0-phi)/Scc00)/s+1d0)*(1d0/Fs-s)*(phi*S00/(1d0-phi)/Scc00+1d0)+s)  /  &
        (1d0/phi/S00+1d0/(1d0-phi)/Scc00)


  end function Ss4


  function Ss_kernel_re(u) result(v)
    implicit none
    real(kind=XPREC)             :: v
    real(kind=XPREC), intent(in) :: u

    complex(kind=XPREC)          :: s, ssval

    s     = cmplx(xil,u)

    ssval = Ss4(s, xS00, xa1, xa2, xr1, xr2, xr3, xScc00, xb1, xb2, xg1, xg2, xg3, xphi)
    ssval = ssval * exp(t_param * s) * exp(-epap*u*u)
    ssval = ssval / (2 * Pi ) 

    v = Realpart(ssval)

  end function Ss_kernel_re



  function Ss_kernel_im(u) result(v)
    implicit none
    real(kind=XPREC)             :: v
    real(kind=XPREC), intent(in) :: u

    complex(kind=XPREC)          :: s, ssval

    s     = cmplx(xil,u)

    ssval = Ss4(s, xS00, xa1, xa2, xr1, xr2, xr3, xScc00, xb1, xb2, xg1, xg2, xg3, xphi) 
    ssval = ssval * exp(t_param * s) * exp(-epap*u*u)
    ssval = ssval / (2 * Pi ) 

    v = Imagpart(ssval)

  end function Ss_kernel_im









 function Ss_test4(s, S00, a1, a2, r1, r2, r3, Scc00, b1, b2, g1, g2, g3, phi) result(Ss)
!---------========-----------------------------------------------------------------------
   implicit none
   complex(kind=XPREC)             :: Ss       ! Result is the Laplace transform of the RPA expression
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   complex(kind=XPREC), intent(in) :: s        ! Laplace "frequency" variable
   real   (kind=XPREC), intent(in) :: S00      ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi      ! volume fraction of sample polymer

   real   (kind=XPREC)             :: a3, b3

   complex(kind=XPREC)  ::  t4  
   complex(kind=XPREC)  ::  t5 
   complex(kind=XPREC)  ::  t7 
   complex(kind=XPREC)  ::  t8 
   complex(kind=XPREC)  ::  t10 
   complex(kind=XPREC)  ::  t15 
   complex(kind=XPREC)  ::  t16 
   complex(kind=XPREC)  ::  t17 
   complex(kind=XPREC)  ::  t19 
   complex(kind=XPREC)  ::  t22 
   complex(kind=XPREC)  ::  t25 
   complex(kind=XPREC)  ::  t28 
   complex(kind=XPREC)  ::  t30 
   complex(kind=XPREC)  ::  t33 
   complex(kind=XPREC)  ::  t34 
   complex(kind=XPREC)  ::  t36 
   complex(kind=XPREC)  ::  t38 
   complex(kind=XPREC)  ::  t41 
   complex(kind=XPREC)  ::  t44 
   complex(kind=XPREC)  ::  t47 
   complex(kind=XPREC)  ::  t48 
   complex(kind=XPREC)  ::  t49 
   complex(kind=XPREC)  ::  t50 
   complex(kind=XPREC)  ::  t52 
   complex(kind=XPREC)  ::  t53 
   complex(kind=XPREC)  ::  t62 
   complex(kind=XPREC)  ::  t65 
   complex(kind=XPREC)  ::  t69 
   complex(kind=XPREC)  ::  t71 
   complex(kind=XPREC)  ::  t75 
   complex(kind=XPREC)  ::  t78 
   complex(kind=XPREC)  ::  t83 
   complex(kind=XPREC)  ::  t85 
   complex(kind=XPREC)  ::  t86 
   complex(kind=XPREC)  ::  t103
   complex(kind=XPREC)  ::  t106
   complex(kind=XPREC)  ::  t108
   complex(kind=XPREC)  ::  t110
   complex(kind=XPREC)  ::  t111
   complex(kind=XPREC)  ::  t122
   complex(kind=XPREC)  ::  t125
   complex(kind=XPREC)  ::  t129
   complex(kind=XPREC)  ::  t130
   complex(kind=XPREC)  ::  t136
   complex(kind=XPREC)  ::  t139
   complex(kind=XPREC)  ::  t156
   complex(kind=XPREC)  ::  t174
   complex(kind=XPREC)  ::  t178
   complex(kind=XPREC)  ::  t181
   complex(kind=XPREC)  ::  t182
   complex(kind=XPREC)  ::  t183
   complex(kind=XPREC)  ::  t200
   complex(kind=XPREC)  ::  t221
   complex(kind=XPREC)  ::  t248
   complex(kind=XPREC)  ::  t249
   complex(kind=XPREC)  ::  t253
   complex(kind=XPREC)  ::  t254
   complex(kind=XPREC)  ::  t268
   complex(kind=XPREC)  ::  t272
   complex(kind=XPREC)  ::  t276
   complex(kind=XPREC)  ::  t277
   complex(kind=XPREC)  ::  t279
   complex(kind=XPREC)  ::  t280
   complex(kind=XPREC)  ::  t286
   complex(kind=XPREC)  ::  t287
   complex(kind=XPREC)  ::  t308
   complex(kind=XPREC)  ::  t309
   complex(kind=XPREC)  ::  t315
   complex(kind=XPREC)  ::  t318
   complex(kind=XPREC)  ::  t320
   complex(kind=XPREC)  ::  t328
   complex(kind=XPREC)  ::  t335
   complex(kind=XPREC)  ::  t346
   complex(kind=XPREC)  ::  t350
   complex(kind=XPREC)  ::  t352
   complex(kind=XPREC)  ::  t359
   complex(kind=XPREC)  ::  t362
   complex(kind=XPREC)  ::  t371
   complex(kind=XPREC)  ::  t375
   complex(kind=XPREC)  ::  t379
   complex(kind=XPREC)  ::  t380
   complex(kind=XPREC)  ::  t396
   complex(kind=XPREC)  ::  t400
   complex(kind=XPREC)  ::  t406
   complex(kind=XPREC)  ::  t433
   complex(kind=XPREC)  ::  t443


!! -- start code (created by maple Python(Ss_test4,optimize);)

    a3 = 1 -(a1+a2)
    b3 = 1 -(b1+b2)


    t4 = b1 + b2 + b3 - 1
    t5 = a1 + a2 + a3
    t7 = Scc00 * t5 * t4
    t8 = b1 + b2 + b3
    t10 = a1 + a2 + a3 - 1
    t15 = s ** 2
    t16 = t15 ** 2
    t17 = t16 * s
    t19 = b2 + b3 - 1
    t22 = b1 + b3 - 1
    t25 = b1 + b2 - 1
    t28 = a2 + a3
    t30 = a1 + a3
    t33 = r3 * (a1 + a2)
    t34 = r1 * t28 + r2 * t30 + t33
    t36 = -g1 * t19 * t5 - g2 * t5 * t22 - g3 * t5 * t25 - t34 * t4
    t38 = b2 + b3
    t41 = b1 + b3
    t44 = b1 + b2
    t47 = a2 + a3 - 1
    t48 = r1 * t47
    t49 = a1 + a3 - 1
    t50 = r2 * t49
    t52 = r3 * (a1 + a2 - 1)
    t53 = t48 + t50 + t52
    t62 = -1 + b3
    t65 = -1 + b2
    t69 = -g2 * t62 * t5 - g3 * t65 * t5 - t34 * t19
    t71 = -1 + b1
    t75 = -g3 * t5 * t71 - t34 * t22
    t78 = g3 * t34 * t25
    t83 = r2 * r3
    t85 = (a2 * r3 + a3 * r2) * r1 + a1 * t83
    t86 = t4 * t85
    t103 = -1 + a3
    t106 = r3 * (a2 - 1)
    t108 = r1 * (r2 * t103 + t106)
    t110 = (-1 + a1) * t83
    t111 = t108 + t110
    t122 = t15 * s
    t125 = t34 * t62
    t129 = g3 * t34 * t65
    t130 = t19 * t85
    t136 = -g3 * t34 * t71 - t22 * t85
    t139 = g3 * t25 * t85
    t156 = r1 * r2
    t174 = t85 * t62
    t178 = g3 * t65 * t85
    t181 = t85 * g2
    t182 = g3 * t71
    t183 = t182 * t181
    t200 = g3 * r1
    t221 = g1 * g3 * Scc00 * t181
    t248 = r1 + r2 + r3
    t249 = t248 * t4
    t253 = g2 * t10
    t254 = g3 * t10
    t268 = t19 * t248
    t272 = t248 * t22
    t276 = g3 * t248 * t25
    t277 = r2 + r3
    t279 = r1 * t277 + t83
    t280 = t279 * t4
    t286 = g2 * (t254 + t48 + t50 + t52)
    t287 = g3 * t53
    t308 = g3 * t65 * t248
    t309 = t279 * t19
    t315 = -g3 * t248 * t71 - t279 * t22
    t318 = g3 * t279 * t25
    t320 = t4 * r3 * t156
    t328 = r3 * t156
    t335 = t62 * r3
    t346 = t279 * t62
    t350 = g3 * t65 * t279
    t352 = t19 * r3 * t156
    t359 = -g3 * t279 * t71 - t22 * r3 * t156
    t362 = t25 * t83 * t200
    t371 = r1 * (-r2 * t103 - t106)
    t375 = g3 * (t371 - t110)
    t379 = g2 * (t375 + t328)
    t380 = t83 * t200
    t396 = t335 * t156
    t400 = t65 * t83 * t200
    t406 = t71 * g2 * r3 * r2 * t200
    t433 = (S00 - Scc00) * phi + Scc00

    t443=-1/t433/(t16*t15*(phi*(t10*S00-Scc00*t4)+Scc00*t4)+t17*(phi*(Scc00*(-g1*t19-g2*t22-g3*t25-t249)&
        +S00*(g1*t10+t253+t254+t48+t50+t52))+(g1*t19+g2*t22+g3*t25+t249)*Scc00)+t16*(phi*(Scc00*(g1*(-g2&
        *t62-g3*t65-t268)+g2*(-g3*t71-t272)-t276-t280)+(g1*(t253+t254+t48+t50+t52)+t286+t287+t108+t110)*&
        S00)+(g1*(g2*t62+g3*t65+t268)+g2*(t182+t272)+t276+t280)*Scc00)+t122*(phi*(Scc00*(g1*(g2*(-r1*t62&
        -r2*t62-t62*r3+g3)-t308-t309)+g2*t315-t318-t320)+(g1*(t286+t287+t108+t110)+g2*(t287+t108+t110)+&
        g3*t111-t328)*S00)+(g1*(g2*(r1*t62+r2*t62-g3+t335)+t308+t309)-g2*t315+t318+t320)*Scc00)+t15*(phi&
        *(Scc00*(g1*(g2*(g3*t248-t346)-t350-t352)+g2*t359-t362)-(g1*(g2*(g3*(-r1*t47-r2*t49-t52)+t371-&
        t110)+t375+t328)+t379+t380)*S00)+(g1*(g2*(-g3*t248+t346)+t350+t352)-g2*t359+t362)*Scc00)+s*(phi*&
        (Scc00*(g1*(g2*(g3*t279-t396)-t400)-t406)-(g1*(t379+t380)+g2*t83*t200)*S00)+(g1*(g2*(g3*(-r1*t277&
        -t83)+t396)+t400)+t406)*Scc00)-r1*t433*g3*r3*g2*r2*g1)*(t17*(phi*(t10*t8*S00-t7)+t7)+t16*(phi*&
        (Scc00*t36+(g1*t10*t38+g2*t10*t41+g3*t10*t44+t8*t53)*S00)-Scc00*t36)+t122*(phi*(Scc00*(g1*t69+g2&
        *t75-t78-t86)+(g1*(g3*t10*b2+g2*t10*b3+t53*t38)+g2*(g3*t10*b1+t41*t53)+g3*t44*t53+t8*t111)*S00)+&
       (-g1*t69-g2*t75+t78+t86)*Scc00)+t15*(phi*(Scc00*(g1*(g2*(g3*t5-t125)-t129-t130)+g2*t136-t139)+(g1*&
       (g3*t53*b2+g2*t53*b3+t111*t38)+g2*(g3*t53*b1+t41*t111)+g3*t44*t111-t8*r3*t156)*S00)+(g1*(g2*(-g3*&
        t5+t125)+t129+t130)-g2*t136+t139)*Scc00)+s*(phi*(Scc00*(g1*(g2*(g3*t34-t174)-t178)-t183)-(g1*(-g3&
        *t111*b2-g2*t111*b3+t38*r3*t156)+g2*(-g3*t111*b1+t41*r3*t156)+t44*t83*t200)*S00)+(g1*(g2*(g3*(-r1&
        *t28-r2*t30-t33)+t174)+t178)+t183)*Scc00)+phi*(t221-S00*r2*((b2*g3+b3*g2)*g1+g3*b1*g2)*r3*r1)-&
        t221)*(-1+phi)*Scc00*S00*phi

    
    Ss = t443
    
end function Ss_test4





subroutine PolyRootCoeffs4(coeffs, S00, a1, a2, r1, r2, r3, Scc00, b1, b2, g1, g2, g3, phi)
  implicit none
   complex(kind=XPREC), intent(out) :: coeffs(0:6) ! Coefficients for the root polynome that results
                                                   ! from the inversel Laplace Transform
                                                   ! needed to determine the zero location that enter
                                                   ! the exp sum that yields the results of the invlaplace tr
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

 
   real   (kind=XPREC), intent(in) :: S00      ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi      ! volume fraction of sample polymer


   real   (kind=XPREC) :: a3, b3
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)

   coeffs    = 0

   coeffs(0) = -S00*g1*g2*g3*phi*r1*r2*r3+Scc00*g1*g2*g3*phi*r1*r2*r3-Scc00*g1*g2*g3*r1*r2*r3

   coeffs(1) =(((((-(-1+b3)*(-1+phi)*Scc00-phi*S00)*g2-g3*((-1+b2)*(-1+phi)*Scc00+phi*S00))*g1-&
             g3*((-1+b1)*(-1+phi)*Scc00+phi*S00)*g2)*r3+g1*g3*g2*((-1+phi)*Scc00+S00*phi*(-1+a3)))&
             *r2+g1*g3*r3*((-1+phi)*Scc00+S00*phi*(a2-1))*g2)*r1+g1*g3*r3*g2*((-1+phi)*Scc00+S00*phi*(-1+a1))*r2
 
   coeffs(2) =(((-((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*(-1+phi)*Scc00-S00*phi*(g1+g2+g3))*r3-&
              (((-1+b3)*g2+(-1+b2)*g3)*g1+g2*g3*(-1+b1))*(-1+phi)*Scc00+phi*(-1+a3)*((g2+g3)*g1+g2*g3)&
              *S00)*r2+(-(((-1+b3)*g2+(-1+b2)*g3)*g1+g2*g3*(-1+b1))*(-1+phi)*Scc00+phi*(a2-1)*((g2+g3)&
              *g1+g2*g3)*S00)*r3+g1*g3*g2*((-1+phi)*Scc00+S00*phi*(a2+a3-1)))*r1+((-(((-1+b3)*g2+(-1+b2)&
              *g3)*g1+g2*g3*(-1+b1))*(-1+phi)*Scc00+phi*((g2+g3)*g1+g2*g3)*S00*(-1+a1))*r3+g1*g3*g2*&
              ((-1+phi)*Scc00+S00*phi*(a1+a3-1)))*r2+g1*g3*r3*g2*((-1+phi)*Scc00+S00*phi*(a1+a2-1))

   coeffs(3) =-((((b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r2+((b2+b3-1)*g1+(b1+b3-1)&
              *g2+(b1+b2-1)*g3)*r3+((-1+b3)*g2+(-1+b2)*g3)*g1+g2*g3*(-1+b1))*r1+(((b2+b3-1)*g1+(b1+b3-1)&
              *g2+(b1+b2-1)*g3)*r3+((-1+b3)*g2+(-1+b2)*g3)*g1+g2*g3*(-1+b1))*r2+(((-1+b3)*g2+(-1+b2)*g3)&
              *g1+g2*g3*(-1+b1))*r3-g1*g3*g2)*(-1+phi)*Scc00+phi*S00*(((-r3+(-1+a3)*g1+(-1+a3)*g2+(-1+a3)&
              *g3)*r2+(g1+g2+g3)*(a2-1)*r3+(a2+a3-1)*((g2+g3)*g1+g2*g3))*r1+((g1+g2+g3)*(-1+a1)*r3+((g2+g3)&
              *g1+g2*g3)*(a1+a3-1))*r2+(a1+a2-1)*((g2+g3)*g1+g2*g3)*r3+g1*g2*g3*(a1+a2+a3-1))

   coeffs(4) =-(((b1+b2+b3-1)*r2+(b1+b2+b3-1)*r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r1+((b1+b2+b3-1)&
              *r3+(b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r2+((b2+b3-1)*g1+(b1+b3-1)*g2+(b1+b2-1)*g3)*r3+&
              ((-1+b3)*g2+(-1+b2)*g3)*g1+g2*g3*(-1+b1))*(-1+phi)*Scc00+phi*S00*(((-1+a3)*r2+r3*(a2-1)+&
              (g1+g2+g3)*(a2+a3-1))*r1+((-1+a1)*r3+(g1+g2+g3)*(a1+a3-1))*r2+(g1+g2+g3)*(a1+a2-1)*r3+&
              (a1+a2+a3-1)*((g2+g3)*g1+g2*g3))

   coeffs(5) =-((b1+b2+b3-1)*r1+(b1+b2+b3-1)*r2+(b1+b2+b3-1)*r3+(g2+g3)*b1+(g1+g3)*b2+(g1+g2)*b3-g1-g2-g3)&
              *(-1+phi)*Scc00+phi*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)+(g1+g2+g3)*(a1+a2+a3-1))*S00

   coeffs(6) =-(-1 + phi) * (b1 + b2 + b3 - 1) * Scc00 + S00 * phi * (a1 + a2 + a3 - 1)

end subroutine PolyRootCoeffs4




function InvLaplace4a(t, S00, a1, a2, r1, r2, r3, Scc00, b1, b2, g1, g2, g3, phi ) result (St)
!---------------------------------------------------------------------------------------------

  implicit none
   complex(kind=XPREC)  :: St                  ! S(Q,t)  (Q only implicit)
                                               ! for a 2-component system for any give q and with
                                               ! time functions of the undisturbed S modelled by 3 exp's
                                               ! see: Akcasu + Tombakoglu, Macromolecules 23, (1990) 607
                                               ! Maple WS_ RpaLaplace_corr5.mw

   real   (kind=XPREC), intent(in) :: t 
   real   (kind=XPREC), intent(in) :: S00      ! S(q) value of the "sample" component                  
   real   (kind=XPREC), intent(in) :: a1,a2    ! Amplitudes of the exp time fkt components (a3=1-a1-a2) inferred
   real   (kind=XPREC), intent(in) :: r1,r2,r3 ! correstonding rates for the time-fkts  (for sample)
   real   (kind=XPREC), intent(in) :: Scc00    ! S(q) value of the "matrix" component                  
   real   (kind=XPREC), intent(in) :: b1,b2    ! Amplitudes of the exp time fkt components (b3=1-b1-b2) inferred
   real   (kind=XPREC), intent(in) :: g1,g2,g3 ! correstonding rates for the time-fkts  (for matrix)
   real   (kind=XPREC), intent(in) :: phi      ! volume fraction of sample polymer
 
  integer, parameter   :: maxdeg = 6
  integer              :: degree = maxdeg
  complex(kind=XPREC)  :: coeffs(0:maxdeg), pzero(1:maxdeg), Z 
  real   (kind=XPREC)  :: Prefak
  integer              :: i

  real   (kind=XPREC), parameter :: epsilon = 1d-8
  real   (kind=XPREC) :: a3, b3
   
   
   a3 = 1 -(a1+a2)
   b3 = 1 -(b1+b2)

! write(6,*)'enter invaplace..'

  call PolyRootCoeffs4(coeffs, S00, a1, a2, r1, r2, r3, Scc00, b1, b2, g1, g2, g3, phi)

!! if the above conditions for a3 and b3 are enforceed (i.e. a1+a2+a3=1, ...)
!! the the highest coefficients are identical to zero, which is not well received
!! by the roots search, therefore we check here for this and in case reduce the
!! degree of the polynomial
dl:  do i=degree,1,-1
       if(abs(coeffs(i)) > epsilon * sum(abs(coeffs(0:degree)))) then
         degree = i
         exit dl
       endif
     enddo dl

  call cmplx_roots_gen(pzero, coeffs, degree, .true. , .false. )

! write(6,*)'roots: ',pzero

  call check_degeneracy(degree, pzero)


  Prefak =  1 / (-phi * S00 + Scc00 * phi - Scc00)

! write(6,*)'Prefak=', Prefak

  St = 0

  do i = 1, degree
   Z = pzero(i)
   St= St + &
      (((b1+b2+b3-1)*(a1+a2+a3)*(-1+phi)*Scc00-S00*phi*(b1+b2+b3)*(a1+a2+a3-1))*Z**5+(((a1+a2+a3)*(b2+b3-1)*g1+&
      (b1+b3-1)*(a1+a2+a3)*g2+(b1+b2-1)*(a1+a2+a3)*g3+(b1+b2+b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*(-1+phi)&
      *Scc00-S00*((b2+b3)*(a1+a2+a3-1)*g1+(b1+b3)*(a1+a2+a3-1)*g2+(b1+b2)*(a1+a2+a3-1)*g3+((a2+a3-1)*r1+(a1+a3-1)&
      *r2+r3*(a1+a2-1))*(b1+b2+b3))*phi)*Z**4+((-1+phi)*(((a1+a2+a3)*(-1+b3)*g2+(a1+a2+a3)*(-1+b2)*g3+(b2+b3-1)*&
      ((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))*g1+((-1+b1)*(a1+a2+a3)*g3+(b1+b3-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2)))&
      *g2+(b1+b2-1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+r2*r3*a1)*(b1+b2+b3-1))*Scc00-S00*phi&
      *((b3*(a1+a2+a3-1)*g2+b2*(a1+a2+a3-1)*g3+(b2+b3)*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1)))*g1+(b1*(a1+a2+a3-1)&
      *g3+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b3))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*(b1+b2)*g3+&
      (((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*(b1+b2+b3)))*Z**3+(((((-a1-a2-a3)*g3+(-1+b3)*((a2+a3)*r1+(a1+a3)*&
      r2+r3*(a1+a2)))*g2+(-1+b2)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+r2*r3*a1)*(b2+b3-1))*g1+&
      ((-1+b1)*((a2+a3)*r1+(a1+a3)*r2+r3*(a1+a2))*g3+((a2*r3+a3*r2)*r1+r2*r3*a1)*(b1+b3-1))*g2+((a2*r3+a3*r2)*r1+&
      r2*r3*a1)*(b1+b2-1)*g3)*(-1+phi)*Scc00-S00*((b3*((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+b2*((a2+a3-1)*&
      r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+(b2+b3)*(((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1)))*g1+(b1*((a2+a3-1)*r1+&
      (a1+a3-1)*r2+r3*(a1+a2-1))*g3+(((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*(b1+b3))*g2+(((-1+a3)*r2+r3*(a2-1))&
      *r1+r2*r3*(-1+a1))*(b1+b2)*g3-r1*r2*r3*(b1+b2+b3))*phi)*Z**2+((-1+phi)*(((((-a2-a3)*r1+(-a1-a3)*r2-r3*(a1+a2))&
      *g3+(-1+b3)*((a2*r3+a3*r2)*r1+r2*r3*a1))*g2+((a2*r3+a3*r2)*r1+r2*r3*a1)*(-1+b2)*g3)*g1+g2*((a2*r3+a3*r2)*r1&
      +r2*r3*a1)*(-1+b1)*g3)*Scc00+S00*((-b3*(((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*g2-b2*(((-1+a3)*r2+r3*(a2-1))&
      *r1+r2*r3*(-1+a1))*g3+r1*r2*r3*(b2+b3))*g1+(-b1*(((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*g3+r1*r2*r3*(b1+b3))&
      *g2+g3*r1*r2*r3*(b1+b2))*phi)*Z-g2*(-1+phi)*g1*((a2*r3+a3*r2)*r1+r2*r3*a1)*g3*Scc00+S00*((b2*g3+b3*g2)*g1+g3*&
      b1*g2)*r2*phi*r3*r1)*exp(Z*t)/((6*(b1+b2+b3-1)*(-1+phi)*Scc00-6*S00*phi*(a1+a2+a3-1))*Z**5+(5*((b2+b3-1)*g1+&
      (b1+b3-1)*g2+(b1+b2-1)*g3+(b1+b2+b3-1)*(r1+r2+r3))*(-1+phi)*Scc00-5*((a1+a2+a3-1)*g1+(a1+a2+a3-1)*g2+&
      (a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*S00*phi)*Z**4+(4*(((-1+b3)*g2+(-1+b2)*g3+(r1+r2+r3)&
      *(b2+b3-1))*g1+((-1+b1)*g3+(b1+b3-1)*(r1+r2+r3))*g2+(b1+b2-1)*(r1+r2+r3)*g3+(b1+b2+b3-1)*((r2+r3)*r1+r2*r3))&
      *(-1+phi)*Scc00-4*S00*(((a1+a2+a3-1)*g2+(a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g1+&
      ((a1+a2+a3-1)*g3+(a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))&
      *g3+((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*phi)*Z**3+(3*(((-g3+(-1+b3)*r1+(-1+b3)*r2+r3*(-1+b3))*g2+&
      (r1+r2+r3)*(-1+b2)*g3+(b2+b3-1)*((r2+r3)*r1+r2*r3))*g1+((-1+b1)*(r1+r2+r3)*g3+(b1+b3-1)*((r2+r3)*r1+r2*r3))&
      *g2+(b1+b2-1)*((r2+r3)*r1+r2*r3)*g3+r1*r2*r3*(b1+b2+b3-1))*(-1+phi)*Scc00-3*S00*((((a1+a2+a3-1)*g3+(a2+a3-1)&
      *r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g2+((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((-1+a3)*r2+r3*(a2-1))&
      *r1+r2*r3*(-1+a1))*g1+(((a2+a3-1)*r1+(a1+a3-1)*r2+r3*(a1+a2-1))*g3+((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))&
      *g2+(((-1+a3)*r2+r3*(a2-1))*r1+r2*r3*(-1+a1))*g3-r1*r2*r3)*phi)*Z**2+(2*(-1+phi)*((((-r1-r2-r3)*g3+(-1+b3)&
      *((r2+r3)*r1+r2*r3))*g2+((r2+r3)*r1+r2*r3)*(-1+b2)*g3+r1*r2*r3*(b2+b3-1))*g1+((-1+b1)*((r2+r3)*r1+r2*r3)*g3&
      +r1*r2*r3*(b1+b3-1))*g2+g3*r1*r2*r3*(b1+b2-1))*Scc00+2*S00*(((((-a2-a3+1)*r1+(-a1-a3+1)*r2-r3*(a1+a2-1))*g3+&
      ((1-a3)*r2-r3*(a2-1))*r1-r2*r3*(-1+a1))*g2+(((1-a3)*r2-r3*(a2-1))*r1-r2*r3*(-1+a1))*g3+r1*r2*r3)*g1+((((1-a3)&
      *r2-r3*(a2-1))*r1-r2*r3*(-1+a1))*g3+r1*r2*r3)*g2+g3*r1*r2*r3)*phi)*Z+(((((-r2-r3)*r1-r2*r3)*g3+r1*r2*r3*&
      (-1+b3))*g2+g3*r1*r2*r3*(-1+b2))*g1+g3*r1*r2*r3*g2*(-1+b1))*(-1+phi)*Scc00+S00*phi*((((((1-a3)*r2-r3*(a2-1))&
      *r1-r2*r3*(-1+a1))*g3+r1*r2*r3)*g2+g3*r1*r2*r3)*g1+g3*r1*r2*r3*g2))
   enddo

   St = St * Prefak * (phi**2-phi)*Scc00*S00
                     !! dieser Faktor fehlte noch! Noch andere??



end function InvLaplace4a
 

subroutine check_degeneracy(degree, roots)
  implicit none
  integer, intent(in)                :: degree
  complex(kind=XPREC), intent(inout) :: roots(degree)

  real(kind=XPREC),parameter    :: epsilon = 1_XPREC-6
  integer                       :: i,j
  complex(kind=XPREC)           :: dr

  do i=1,degree-1
    do j=i+1,degree
      dr = roots(i) - roots(j)
      if(abs(dr) < epsilon) then
        write(6,*)"suspected degeneracy of roots:"
        write(6,*) i, roots(i)
        write(6,*) j, roots(j)
        roots(j) = roots(j) + epsilon  !! tentativ replacement for Invlaplace
      endif
    enddo
  enddo

end subroutine check_degeneracy



!---------------- for the time being include the root finder here -----------------------------------------------
! ---> own module


!
!   Licensed under the Apache License, Version 2.0 (the "License").
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!
!-------------------------------------------------------------------!
!
! The authors also make this file available under the terms of
! GNU Lesser General Public License version 2 or any later version.
! (text of the LGPL licence 2 in NOTICE file)
!
!-------------------------------------------------------------------!
!
! A custom in the scientific comunity is (regardless of the licence
! you chose to use or distribute this software under)
! that if this code was important in the scientific process or
! for the results of your scientific work, we kindly ask you for the
! appropriate citation of the Paper (Skowron & Gould 2012), and
! we would be greatful if you pass the information about
! the proper citation to anyone whom you redistribute this software to.
!
!-------------------------------------------------------------------!
!
!    No    Subroutine
!
!     1   cmplx_roots_gen               - general polynomial solver, works for random degree, not as fast or robust as cmplx_roots_5
!     2   cmplx_roots_5                 - complex roots finding algorithm taylored for 5th order polynomial (with failsafes for polishing)
!     3   sort_5_points_by_separation   - sorting of an array of 5 points, 1st most isolated, 4th and 5th - closest
!     4   sort_5_points_by_separation_i - sorting same as above, returns array of indicies rather than sorted array
!     5   find_2_closest_from_5         - finds closest pair of 5 points
!     6   cmplx_laguerre                - Laguerre's method with simplified Adams' stopping criterion 
!     7   cmplx_newton_spec             - Newton's method with stopping criterion calculated every 10 steps
!     8   cmplx_laguerre2newton         - three regime method: Laguerre's, Second-order General method and Newton's
!     9   solve_quadratic_eq            - quadratic equation solver
!    10   solve_cubic_eq                - cubic equation solver based on Lagrange's method
!    11   divide_poly_1                 - division of the polynomial by (x-p)
!
! fortran 90 code
!
! Paper:  Skowron & Gould 2012
!         "General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses"
!
! for a full text see:
!     http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
!     or http://arxiv.org/find/astro-ph
!     or http://www.adsabs.harvard.edu/abstract_service.html
! see also file NOTICE and LICENSE
!
! ver. 2012.03.03 initial
! ver. 2014.03.12 bug fix
! ver. 2016.01.21 bug fix
! ver. 2016.04.28 bug fix
!
!-------------------------------------------------------------------!
! _1_                  CMPLX_ROOTS_GEN                              !
!-------------------------------------------------------------------!
subroutine cmplx_roots_gen(roots, poly, degree, polish_roots_after, use_roots_as_starting_points)
  ! This subroutine finds roots of a complex polynomial. 
  ! It is general, however less fast or robust than cmplx_roots_5
  ! which contains failsafe checks in the polishing stage, but is
  ! designed only for 5th order polynomials.
  ! It uses a new dynamic root finding algorithm (see the Paper).
  !
  ! It can use Laguerre's method (subroutine cmplx_laguerre)
  ! or Laguerre-.gt.SG-.gt.Newton method (subroutine 
  ! cmplx_laguerre2newton - this is default choice) to find 
  ! roots. It divides polynomial one by one by found roots. At the 
  ! end it finds last root from Viete's formula for quadratic 
  ! equation. Finally, it polishes all found roots using a full
  ! polynomial and Newton's or Laguerre's method (default is
  ! Laguerre's - subroutine cmplx_laguerre). 
  ! You can change default choices by commenting out and uncommenting
  ! certain lines in the code below.
  !
  ! Note:
  ! - we solve for the last root with Viete's formula rather 
  !   than doing full Laguerre step (which is time consuming
  !   and unnecessary)
  ! - we do not introduce any preference to real roots
  ! - in Laguerre implementation we omit unneccesarry calculation of
  !   absolute values of denominator
  ! - we do not sort roots. If you need to sort 
  !   roots - we have provided sorting subroutine called:
  !   sort_5_points_by_separation, which sorts points from most 
  !   isolated to most close. Algorithm in this routine can be 
  !   easily used for number of points different than 5.
  !
  implicit none
  ! roots  - array which will hold all roots that had been found.
  !          If the flag 'use_roots_as_starting_points' is set to 
  !          .true., then instead of point (0,0) we use value from
  !          this array as starting point for cmplx_laguerre
  ! poly -   is an array of polynomial cooefs, length = degree+1, 
  !          poly(1) is a constant term:
  !               1              2             3
  !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  ! degree - degree of the polynomial and size of 'roots' array
  ! polish_roots_after - after all roots have been found by dividing
  !          original polynomial by each root found,
  !          you can opt in to polish all roots using full  
  !          polynomial
  ! use_roots_as_starting_points - usually we start Laguerre's 
  !          method from point (0,0), but you can decide to use the 
  !          values of 'roots' array as starting point for each new
  !          root that is searched for. This is useful if you have
  !          very rough idea where some of the roots can be. 
  !
  integer, parameter :: RK = 8
  integer, intent(in) :: degree
  complex(kind=RK), dimension(degree+1), intent(in) :: poly ! coeffs of the polynomial
  complex(kind=RK), dimension(degree), intent(inout) :: roots
  logical, intent(in) :: polish_roots_after, use_roots_as_starting_points

  complex(kind=RK), dimension(degree+1) :: poly2
  complex(kind=RK), parameter :: zero = cmplx(0d0,0d0,RK)
  integer :: i, n, iter
  logical :: success
  complex(kind=RK) :: coef, prev

  poly2=poly
  
  ! initialize starting points
  if(.not.use_roots_as_starting_points) roots=zero

  ! skip small degree polynomials from doing Laguerre's method
  if(degree <= 1)then
    if(degree==1) roots(1)=-poly(1)/poly(2)
    return
  endif


  do n=degree, 3, -1

    ! find root with Laguerre's method
    !call cmplx_laguerre(poly2, n, roots(n), iter, success) 
    ! or
    ! find root with (Laguerre's method -.gt. SG method -.gt. Newton's method)
    call cmplx_laguerre2newton(poly2, n, roots(n), iter, success, 2)
    if(.not.success) then
      roots(n)=zero
      call cmplx_laguerre(poly2, n, roots(n), iter, success)
    endif

    ! divide the polynomial by this root
    coef=poly2(n+1)
    do i=n,1,-1
      prev=poly2(i)
      poly2(i)=coef
      coef=prev+roots(n)*coef
    enddo
    ! variable coef now holds a remainder - should be close to 0

  enddo

  ! find all but last root with Laguerre's method
  !call cmplx_laguerre(poly2, 2, roots(2), iter, success)
  ! or
  call cmplx_laguerre2newton(poly2, 2, roots(2), iter, success, 2)
  if(.not.success) then
    call solve_quadratic_eq(roots(2),roots(1),poly2)
  else
    ! calculate last root from Viete's formula
    roots(1)=-(roots(2)+poly2(2)/poly2(3))
  endif


  if(polish_roots_after)then
    do n=1, degree ! polish roots one-by-one with a full polynomial 
      call cmplx_laguerre(poly, degree, roots(n), iter, success) 
      !call cmplx_newton_spec(poly, degree, roots(n), iter, success)
    enddo
  endif  
 
  return
end


!-------------------------------------------------------------------!
! _2_                     CMPLX_ROOTS_5                             !
!-------------------------------------------------------------------!
subroutine cmplx_roots_5(roots, first_3_roots_order_changed, poly, polish_only)
  implicit none
  ! Subroutine finds or polishes roots of a complex polynomial 
  ! (degree=5)
  ! This routine is especially tailored for solving binary lens 
  ! equation in form of 5th order polynomial. 
  !
  ! Use of this routine, in comparission to 'cmplx_roots_gen' can yield
  ! consideribly faster code, because it makes polishing of the roots 
  ! (that come in as a guess from previous solutions) secure by
  ! implementing additional checks on the result of polishing. 
  ! If those checks are not satisfied then routine reverts to the 
  ! robust algorithm. These checks are designed to work for 5th order 
  ! polynomial originated from binary lens equation.
  !
  ! Usage:
  !
  ! polish_only == false - I do not know the roots, routine should  
  !                find them from scratch. At the end it
  !                sorts roots from the most distant to closest.
  !                Two last roots are the closest (in no particular
  !                order).
  ! polish_only = true - I do know the roots pretty well, for e&e
  !                I have changed the coefficiens of the polynomial 
  !                only a bit, so the two closest roots are 
  !                most likely still the closest ones.
  !                If the output flag 'first_3_roots_order_changed'
  !                is returned as 'false', then first 3 returned roots
  !                are in the same order as initialy given to the 
  !                routine. The last two roots are the closest ones, 
  !                but in no specific order (!).
  !                If 'first_3_roots_order_changed' is 'true' then
  !                it means that all roots had been resorted.
  !                Two last roots are the closest ones. First is most 
  !                isolated one.
  !
  !
  ! If you do not know the position of the roots just use flag
  ! polish_only=.false. In this case routine will find the roots by
  ! itself.
  
  ! Returns all five roots in the 'roots' array.
  !
  ! poly  - is an array of polynomial cooefs, length = degree+1 
  !       poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + poly(4) x^3 + ...
  ! roots - roots of the polynomial ('out' and optionally 'in')
  !
  !
  ! Jan Skowron 2011
  !
  integer, parameter :: degree=5
  integer, parameter :: RK=8  ! kind for real and complex variables
  complex(kind=RK), dimension(degree), intent(inout) :: roots
  logical, intent(out) :: first_3_roots_order_changed
  complex(kind=RK), dimension(degree+1), intent(in)  :: poly
  logical, intent(in) :: polish_only
  !------------------------------------------------------------------
  !
  complex(kind=RK) :: remainder, roots_robust(degree)
  real(kind=RK) :: d2min
  integer :: iter, loops, go_to_robust, m, root4, root5, i, i2
  complex(kind=RK), dimension(degree+1) :: poly2
  !integer, dimension(degree) :: sorted_ind
  complex(kind=RK), parameter :: zero=cmplx(0d0,0d0,RK)
  logical :: succ


  !---------------------------------------
  roots_robust = roots
  
  go_to_robust=0
  if(.not.polish_only) then
    ! initialize roots
    roots=zero
    go_to_robust=1
  endif
  first_3_roots_order_changed=.false.

  do loops=1,3

    ! ROBUST
    ! (we do not know the roots)
    if(go_to_robust.gt.0) then

      if(go_to_robust.gt.2)then  ! something is wrong
        roots=roots_robust    ! return not-polished roots, because polishing creates errors
        return
      endif

      poly2=poly ! copy coeffs
      do m=degree,4,-1 ! find the roots one-by-one (until 3 are left to be found)
        call cmplx_laguerre2newton(poly2, m, roots(m), iter, succ, 2)
        if(.not.succ)then
          roots(m)=zero
          call cmplx_laguerre(poly2, m, roots(m), iter, succ)
        endif
        ! divide polynomial by this root
        call divide_poly_1(poly2, remainder, roots(m), poly2, m)
      enddo
      ! find last 3 roots with cubic euqation solver (Lagrange's method)
      call solve_cubic_eq(roots(1),roots(2),roots(3),poly2)
      ! all roots found
      
      ! sort roots - first will be most isolated, last two will be the closest
      call sort_5_points_by_separation(roots) 
      ! copy roots in case something will go wrong during polishing
      roots_robust=roots

      ! set flag, that roots have been resorted
      first_3_roots_order_changed=.true.
    endif  ! go_to_robust.gt.0

    ! POLISH 
    ! (we know the roots approximately, and we guess that last two are closest)
    !---------------------
      poly2=poly ! copy coeffs

      do m=1,degree-2
      !do m=1,degree                      ! POWN - polish only with Newton (option)

        ! polish roots with full polynomial
        call cmplx_newton_spec(poly2, degree, roots(m), iter, succ)

        if(.not.succ)then
          ! go back to robust
          go_to_robust=go_to_robust+1
          roots=zero
          exit
        endif
      enddo ! m=1,degree-2


      if(succ) then

        ! comment out division and quadratic if you (POWN) polish with Newton only
        do m=1,degree-2
          call divide_poly_1(poly2, remainder, roots(m), poly2, degree-m+1)
        enddo
        ! last two roots are found with quadratic equation solver 
        ! (this is faster and more robust, although little less accurate)
        call solve_quadratic_eq(roots(degree-1), roots(degree  ),poly2)
        ! all roots found and polished

        ! TEST ORDER
        ! test closest roots if they are the same pair as given to polish
        call find_2_closest_from_5(root4,root5, d2min, roots)

        ! check if the closest roots are not too close, this could happen
        ! when using polishing with Newton only, when two roots erroneously 
        ! colapsed to the same root. This check is not needed for polishing
        ! 3 roots by Newton and using quadratic for the remaining two.
        ! If the real roots are so close indeed (very low probability), this will just 
        ! take more time and the unpolished result be returned at the end
        ! but algorithm will work, and will return accurate enough result
        !if(d2min.lt.1d-18) then             ! POWN - polish only with Newton 
        !  go_to_robust=go_to_robust+1    ! POWN - polish only with Newton
        !else                             ! POWN - polish only with Newton

        if((root4.lt.degree-1).or.(root5.lt.degree-1)) then
          ! after polishing some of the 3 far roots become one of the 2 closest ones
          ! go back to robust
          if(go_to_robust.gt.0)then
            ! if came from robust 
            ! copy two most isolated roots as starting points for new robust
            do i=1,degree-3
              roots(degree-i+1)=roots_robust(i)
            enddo 
          else
            ! came from users initial guess
            ! copy some 2 roots (except the closest ones)
            i2=degree
            do i=1,degree
              if((i/=root4).and.(i/=root5))then
                roots(i2)=roots(i)
                i2=i2-1
              endif
              if(i2.le.3) exit ! do not copy those that will be done by cubic in robust
            enddo
          endif
          go_to_robust=go_to_robust+1
        else
          ! root4 and root5 comes from the initial closest pair
          ! most common case
          return
        endif

        !endif                            ! POWN - polish only with Newton
      endif !  
    !---------------------
  enddo ! loops

  return
end


!-------------------------------------------------------------------!
! _3_                SORT_5_POINTS_BY_SEPARATION                    !
!-------------------------------------------------------------------!
subroutine sort_5_points_by_separation(points)
  ! Sort array of five points 
  ! Most isolated point will become the first point in the array
  ! The closest points will be the last two points in the array
  !
  ! Algorithm works well for all dimensions. We put n=5 as 
  ! a hardcoded value just for optimization purposes.
  implicit none
  integer, parameter :: RK=8
  integer, parameter :: n=5 !  works for different n as well, but is faster for n as constant (optimization)
  complex(kind=RK), dimension(n), intent(inout) :: points

  integer, dimension(n) :: sorted_points
  complex(kind=RK), dimension(n) :: savepoints
  integer :: i

  call sort_5_points_by_separation_i(sorted_points, points)
  savepoints=points
  do i=1,n
    points(i)=savepoints(sorted_points(i))
  enddo
  return
end

!-------------------------------------------------------------------!
! _4_              SORT_5_POINTS_BY_SEPARATION_I                    !
!-------------------------------------------------------------------!
subroutine sort_5_points_by_separation_i(sorted_points, points)  
  ! Return index array that sorts array of five points 
  ! Index of the most isolated point will appear on the firts place 
  ! of the output array.
  ! The indices of the closest 2 points will be at the last two 
  ! places in the 'sorted_points' array
  !
  ! Algorithm works well for all dimensions. We put n=5 as 
  ! a hardcoded value just for optimization purposes.
  implicit none
  integer, parameter :: RK=8
  integer, parameter :: n=5 !  works for different n as well, but is faster for n as constant (optimization)
  integer, dimension(n), intent(out) :: sorted_points
  complex(kind=RK), dimension(n), intent(in) :: points
   
  real(kind=RK) :: dmin, d1, d2, d
  real(kind=RK), dimension(n,n) :: distances2
  integer :: ki, kj, ind2, put
  real(kind=RK), dimension(n) :: neigh1st, neigh2nd
  complex(kind=RK) :: p

  distances2=1d100
  dmin=1d100

  do kj=1, n              
    do ki=1, kj-1    
      p=points(ki)-points(kj)
      d=real(conjg(p)*p)
      distances2(ki,kj)=d
      distances2(kj,ki)=d
    enddo
  enddo

  ! find neighbours  
  neigh1st=1d100
  neigh2nd=1d100
  do kj=1, n
    do ki=1, n
      d=distances2(kj,ki)
      if(d.lt.neigh2nd(kj))then
        if(d.lt.neigh1st(kj))then
          neigh2nd(kj)=neigh1st(kj)
          neigh1st(kj)=d
        else
          neigh2nd(kj)=d
        endif
      endif
    enddo
  enddo    

  ! initialize sorted_points
  do ki=1,n
    sorted_points(ki)=ki
  enddo   
 
  ! sort the rest 1..n-2
  do kj=2,n
    d1=neigh1st(kj)
    d2=neigh2nd(kj)
    put=1
    do ki=kj-1,1,-1
      ind2=sorted_points(ki)
      d=neigh1st(ind2)
      if(d.ge.d1) then
        if(d==d1)then
          if(neigh2nd(ind2).gt.d2)then
            put=ki+1
            exit
          endif
        else
          put=ki+1
          exit
        endif
      endif
      sorted_points(ki+1)=sorted_points(ki)
    enddo
    sorted_points(put)=kj
  enddo
    
  return
end


!-------------------------------------------------------------------!
! _5_                     FIND_2_CLOSEST_FROM_5                     !
!-------------------------------------------------------------------!
subroutine find_2_closest_from_5(i1,i2, d2min, points) 
  ! Returns indices of the two closest points out of array of 5
  implicit none
  integer, parameter :: RK=8
  integer, parameter :: n=5 ! will work for other n too, but it is faster with n as constant
  integer, intent(out) :: i1, i2
  !real(kind=RK), dimension(n,n) :: distances2
  complex(kind=RK), dimension(n), intent(in) :: points
  real(kind=RK), intent(out) :: d2min  ! square of minimal distance
   
  real(kind=RK) :: d2min1, d2
  integer :: i,j
  complex(kind=RK) :: p

  d2min1=1d100
  do j=1,n
    !distances2(j,j)=0d0
    do i=1,j-1
      p=points(i)-points(j)
      d2=real(conjg(p)*p)
      !distances2(i,j)=d2
      !distances2(j,i)=d2
      if(d2.le.d2min1)then
        i1=i
        i2=j
        d2min1=d2
      endif
    enddo
  enddo
  
  d2min=d2min1
  
end



!-------------------------------------------------------------------!
! _6_                     CMPLX_LAGUERRE                            !
!-------------------------------------------------------------------!
recursive subroutine cmplx_laguerre(poly, degree, root, iter, success)
  implicit none
  ! Subroutine finds one root of a complex polynomial using 
  ! Laguerre's method. In every loop it calculates simplified 
  ! Adams' stopping criterion for the value of the polynomial.
  !
  ! Uses 'root' value as a starting point (!!!!!)
  ! Remember to initialize 'root' to some initial guess or to 
  ! point (0,0) if you have no prior knowledge.
  !
  ! poly - is an array of polynomial cooefs
  !        length = degree+1, poly(1) is constant 
  !               1              2             3
  !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  ! degree - a degree of the polynomial
  ! root - input: guess for the value of a root
  !        output: a root of the polynomial
  ! iter - number of iterations performed (the number of polynomial
  !        evaluations and stopping criterion evaluation)
  ! success - is false if routine reaches maximum number of iterations
  !
  ! For a summary of the method go to: 
  ! http://en.wikipedia.org/wiki/Laguerre's_method
  !
  integer, parameter :: RK=8  ! kind for real and complex variables
  integer, parameter :: MAX_ITERS=200   ! Laguerre is used as a failsafe
  ! constants needed to break cycles in the scheme
  integer, parameter :: FRAC_JUMP_EVERY=10
  integer, parameter :: FRAC_JUMP_LEN=10
  real(kind=RK), dimension(FRAC_JUMP_LEN), parameter :: FRAC_JUMPS=(/0.64109297d0, &
                                        0.91577881d0, 0.25921289d0,  0.50487203d0, &
                                        0.08177045d0, 0.13653241d0,  0.306162d0  , &
                                        0.37794326d0, 0.04618805d0,  0.75132137d0/) ! some random numbers
  real(kind=RK), parameter :: pi = 3.141592653589793d0
  real(kind=RK) :: faq ! jump length
  real(kind=RK), parameter :: FRAC_ERR = 2.0d-15  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

  integer, intent(in) :: degree
  complex(kind=RK), dimension(degree+1), intent(in)  :: poly
  integer, intent(out) :: iter
  complex(kind=RK), intent(inout) :: root
  logical, intent(out) :: success


  complex(kind=RK) :: p         ! value of polynomial
  complex(kind=RK) :: dp        ! value of 1st derivative 
  complex(kind=RK) :: d2p_half  ! value of 2nd derivative
  integer :: i, k
  logical :: good_to_go
  !complex(kind=RK) :: G, H, G2
  complex(kind=RK) :: denom, denom_sqrt, dx, newroot
  real(kind=RK) :: ek, absroot, abs2p
  complex(kind=RK) :: fac_netwon, fac_extra, F_half, c_one_nth
  real(kind=RK) :: one_nth, n_1_nth, two_n_div_n_1
  complex(kind=RK), parameter :: c_one=cmplx(1d0,0d0,RK)
  complex(kind=RK), parameter :: zero=cmplx(0d0,0d0,RK)
  real(kind=RK) :: stopping_crit2

  !---------------------------------------
 
  iter=0
  success=.true.

  ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
  if(.false.)then ! change false--.gt.true if you would like to use caution about having first coefficient == 0
    if(degree.lt.0) then
      write(*,*) 'Error: cmplx_laguerre: degree.lt.0'
      return
    endif
    if(poly(degree+1)==zero) then
      if(degree==0) return
      call cmplx_laguerre(poly, degree-1, root, iter, success)
      return
    endif
    if(degree.le.1)then
      if(degree==0) then  ! we know from previous check than poly(1) not equal zero
        success=.false.
        write(*,*) 'Warning: cmplx_laguerre: degree=0 and poly(1)/=0, no roots'
        return
      else
        root=-poly(1)/poly(2)
        return
      endif
    endif
  endif
  !  end EXTREME failsafe
    
  good_to_go=.false.
  one_nth=1d0/degree
  n_1_nth=(degree-1d0)*one_nth
  two_n_div_n_1=2d0/n_1_nth
  c_one_nth=cmplx(one_nth,0d0,RK)


  do i=1,MAX_ITERS
    ! prepare stoping criterion
    ek=abs(poly(degree+1))
    absroot=abs(root)
    ! calculate value of polynomial and its first two derivatives
    p  =poly(degree+1)
    dp =zero
    d2p_half=zero
    do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
      d2p_half=dp + d2p_half*root
      dp =p + dp*root
      p  =poly(k)+p*root    ! b_k
      ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
      ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
      ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
      ! Eq 8.
      ek=absroot*ek+abs(p)
    enddo
    iter=iter+1
    
    abs2p=real(conjg(p)*p)
    if(abs2p==0d0) return
    stopping_crit2=(FRAC_ERR*ek)**2
    if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
      ! do additional iteration if we are less than 10x from stopping criterion
      if(abs2p.lt.0.01d0*stopping_crit2) then
        return ! return immediately, because we are at very good place
      else
        good_to_go=.true. ! do one iteration more
      endif
    else
      good_to_go=.false.  ! reset if we are outside the zone of the root
    endif
  
    faq=1d0
    denom=zero
    if(dp/=zero)then
      fac_netwon=p/dp
      fac_extra=d2p_half/dp
      F_half=fac_netwon*fac_extra
      
      denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

      !G=dp/p  ! gradient of ln(p)
      !G2=G*G
      !H=G2-2d0*d2p_half/p  ! second derivative of ln(p)
      !denom_sqrt=sqrt( (degree-1)*(degree*H-G2) )
    
      ! NEXT LINE PROBABLY CAN BE COMMENTED OUT  
      if(real(denom_sqrt).ge.0d0)then
        ! real part of a square root is positive for probably all compilers. You can 
        ! test this on your compiler and if so, you can omit this check
        denom=c_one_nth+n_1_nth*denom_sqrt
      else
        denom=c_one_nth-n_1_nth*denom_sqrt
      endif
    endif
    if(denom==zero)then !test if demoninators are .gt. 0.0 not to divide by zero
      dx=(absroot+1d0)*exp(cmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,RK)) ! make some random jump
    else
      dx=fac_netwon/denom
      !dx=degree/denom
    endif
    
    newroot=root-dx
    if(newroot==root) return ! nothing changes -.gt. return
    if(good_to_go)then       ! this was jump already after stopping criterion was met
      root=newroot
      return
    endif

    if(mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
      faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
      newroot=root-faq*dx ! do jump of some semi-random length (0.lt.faq.lt.1)
    endif
    root=newroot
  enddo
  success=.false.
  ! too many iterations here  
end


!-------------------------------------------------------------------!
! _7_                     CMPLX_NEWTON_SPEC                         !
!-------------------------------------------------------------------!
recursive subroutine cmplx_newton_spec(poly, degree, root, iter, success)
  implicit none
  ! Subroutine finds one root of a complex polynomial using 
  ! Newton's method. It calculates simplified Adams' stopping 
  ! criterion for the value of the polynomial once per 10 iterations (!),
  ! after initial iteration. This is done to speed up calculations
  ! when polishing roots that are known preety well, and stopping
  ! criterion does significantly change in their neighborhood.
  !
  ! Uses 'root' value as a starting point (!!!!!)
  ! Remember to initialize 'root' to some initial guess.
  ! Do not initilize 'root' to point (0,0) if the polynomial 
  ! coefficients are strictly real, because it will make going 
  ! to imaginary roots impossible.
  !
  ! poly - is an array of polynomial cooefs
  !        length = degree+1, poly(1) is constant 
  !               1              2             3
  !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  ! degree - a degree of the polynomial
  ! root - input: guess for the value of a root
  !        output: a root of the polynomial
  ! iter - number of iterations performed (the number of polynomial
  !        evaluations)
  ! success - is false if routine reaches maximum number of iterations
  !
  ! For a summary of the method go to: 
  ! http://en.wikipedia.org/wiki/Newton's_method
  !
  integer, parameter :: RK=8  ! kind for real and complex variables
  integer, parameter :: MAX_ITERS=50
  ! constants needed to break cycles in the scheme
  integer, parameter :: FRAC_JUMP_EVERY=10
  integer, parameter :: FRAC_JUMP_LEN=10
  real(kind=RK), dimension(FRAC_JUMP_LEN), parameter :: FRAC_JUMPS=(/0.64109297d0, &
                                        0.91577881d0, 0.25921289d0,  0.50487203d0, &
                                        0.08177045d0, 0.13653241d0,  0.306162d0  , &
                                        0.37794326d0, 0.04618805d0,  0.75132137d0/) ! some random numbers
  real(kind=8), parameter :: pi = 3.141592653589793d0
  real(kind=RK) :: faq ! jump length
  real(kind=RK), parameter :: FRAC_ERR = 2.0d-15  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

  integer, intent(in) :: degree
  complex(kind=RK), dimension(degree+1), intent(in)  :: poly
  integer, intent(out) :: iter
  complex(kind=RK), intent(inout) :: root
  logical, intent(out) :: success


  complex(kind=RK) :: p    ! value of polynomial
  complex(kind=RK) :: dp   ! value of 1st derivative 
  integer :: i, k
  logical :: good_to_go
  complex(kind=RK) :: dx, newroot
  real(kind=RK) :: ek, absroot, abs2p
  complex(kind=RK), parameter :: zero=cmplx(0d0,0d0,RK)
  real(kind=RK) :: stopping_crit2

  !---------------------------------------


  iter=0
  success=.true.

  ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
  if(.false.)then ! change false--.gt.true if you would like to use caution about having first coefficient == 0
    if(degree.lt.0) then
      write(*,*) 'Error: cmplx_newton_spec: degree.lt.0'
      return
    endif
    if(poly(degree+1)==zero) then
      if(degree==0) return
      call cmplx_newton_spec(poly, degree-1, root, iter, success)
      return
    endif
    if(degree.le.1)then
      if(degree==0) then  ! we know from previous check than poly(1) not equal zero
        success=.false.
        write(*,*) 'Warning: cmplx_newton_spec: degree=0 and poly(1)/=0, no roots'
        return
      else
        root=-poly(1)/poly(2)
        return
      endif
    endif
  endif
  !  end EXTREME failsafe

  good_to_go=.false.

  stopping_crit2 = 0d0  ! value not importat, will be initialized anyway on the first loop (because mod(1,10)==1)
  do i=1,MAX_ITERS
    faq=1d0

    ! prepare stoping criterion
    ! calculate value of polynomial and its first two derivatives
    p  =poly(degree+1)
    dp =zero

    if(mod(i,10)==1) then ! calculate stopping criterion every tenth iteration
      ek=abs(poly(degree+1))
      absroot=abs(root)
      do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
        dp =p + dp*root
        p  =poly(k)+p*root    ! b_k
        ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
        ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
        ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
        ! Eq 8.
        ek=absroot*ek+abs(p)
      enddo
      stopping_crit2=(FRAC_ERR*ek)**2
    else               ! calculate just the value and derivative
      do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
        dp =p + dp*root
        p  =poly(k)+p*root    ! b_k
      enddo
    endif
    iter=iter+1

    
    abs2p=real(conjg(p)*p)
    if(abs2p==0d0) return

    if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
      if(dp==zero) return ! if we have problem with zero, but we are close to the root, just accept
      ! do additional iteration if we are less than 10x from stopping criterion
      if(abs2p.lt.0.01d0*stopping_crit2) then
        return ! return immediately, because we are at very good place
      else
        good_to_go=.true. ! do one iteration more
      endif
    else
      good_to_go=.false. ! reset if we are outside the zone of the root
    endif


    if(dp==zero)then
      ! problem with zero
      dx=(abs(root)+1d0)*exp(cmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,RK)) ! make some random jump
    else
      dx=p/dp  ! Newton method, see http://en.wikipedia.org/wiki/Newton's_method
    endif

  

    newroot=root-dx
    if(newroot==root) return ! nothing changes -.gt. return
    if(good_to_go)then       ! this was jump already after stopping criterion was met
      root=newroot
      return
    endif

    if(mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
      faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
      newroot=root-faq*dx ! do jump of some semi-random length (0.lt.faq.lt.1)
    endif
    root=newroot
  enddo
  success=.false.
  return
  ! too many iterations here  
end



!-------------------------------------------------------------------!
! _8_                     CMPLX_LAGUERRE2NEWTON                     !
!-------------------------------------------------------------------!
recursive subroutine cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode)
  implicit none
  ! Subroutine finds one root of a complex polynomial using 
  ! Laguerre's method, Second-order General method and Newton's
  ! method - depending on the value of function F, which is a 
  ! combination of second derivative, first derivative and
  ! value of polynomial [F=-(p"*p)/(p'p')].
  ! 
  ! Subroutine has 3 modes of operation. It starts with mode=2
  ! which is the Laguerre's method, and continues until F
  ! becames F.lt.0.50, at which point, it switches to mode=1, 
  ! i.e., SG method (see paper). While in the first two
  ! modes, routine calculates stopping criterion once per every
  ! iteration. Switch to the last mode, Newton's method, (mode=0) 
  ! happens when becomes F.lt.0.05. In this mode, routine calculates
  ! stopping criterion only once, at the beginning, under an 
  ! assumption that we are already very close to the root.
  ! If there are more than 10 iterations in Newton's mode, 
  ! it means that in fact we were far from the root, and
  ! routine goes back to Laguerre's method (mode=2).
  !
  ! Uses 'root' value as a starting point (!!!!!)
  ! Remember to initialize 'root' to some initial guess or to 
  ! point (0,0) if you have no prior knowledge.
  !
  ! poly - is an array of polynomial cooefs
  !        length = degree+1, poly(1) is constant 
  !               1              2             3
  !          poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  ! degree - a degree of the polynomial
  ! root - input: guess for the value of a root
  !        output: a root of the polynomial
  ! iter - number of iterations performed (the number of polynomial
  !        evaluations and stopping criterion evaluation)
  ! success - is false if routine reaches maximum number of iterations
  ! starting_mode - this should be by default = 2. However if you  
  !                 choose to start with SG method put 1 instead. 
  !                 Zero will cause the routine to 
  !                 start with Newton for first 10 iterations, and
  !                 then go back to mode 2.
  !                 
  !
  ! For a summary of the method see the paper: Skowron & Gould (2012)
  !
  integer, parameter :: RK=8  ! kind for real and complex variables
  integer, parameter :: MAX_ITERS=50
  ! constants needed to break cycles in the scheme
  integer, parameter :: FRAC_JUMP_EVERY=10
  integer, parameter :: FRAC_JUMP_LEN=10
  real(kind=RK), dimension(FRAC_JUMP_LEN), parameter :: FRAC_JUMPS=(/0.64109297d0, &
                                        0.91577881d0, 0.25921289d0,  0.50487203d0, &
                                        0.08177045d0, 0.13653241d0,  0.306162d0  , &
                                        0.37794326d0, 0.04618805d0,  0.75132137d0/) ! some random numbers
  real(kind=RK), parameter :: pi = 3.141592653589793d0
  real(kind=RK) :: faq ! jump length
  real(kind=RK), parameter :: FRAC_ERR = 2.0d-15  ! fractional error for kind=8 (see. Adams 1967 Eqs 9 and 10)

  integer, intent(in) :: degree
  complex(kind=RK), dimension(degree+1), intent(in)  :: poly
  integer, intent(out) :: iter
  complex(kind=RK), intent(inout) :: root
  integer, intent(in) :: starting_mode
  logical, intent(out) :: success


  complex(kind=RK) :: p         ! value of polynomial
  complex(kind=RK) :: dp        ! value of 1st derivative 
  complex(kind=RK) :: d2p_half  ! value of 2nd derivative
  integer :: i, j, k
  logical :: good_to_go
  !complex(kind=RK) :: G, H, G2
  complex(kind=RK) :: denom, denom_sqrt, dx, newroot
  real(kind=RK) :: ek, absroot, abs2p, abs2_F_half
  complex(kind=RK) :: fac_netwon, fac_extra, F_half, c_one_nth
  real(kind=RK) :: one_nth, n_1_nth, two_n_div_n_1
  integer :: mode
  complex(kind=RK), parameter :: c_one=cmplx(1d0,0d0,RK)
  complex(kind=RK), parameter :: zero=cmplx(0d0,0d0,RK)
  real(kind=RK) :: stopping_crit2
 
  iter=0
  success=.true.
  stopping_crit2 = 0d0  !  value not importat, will be initialized anyway on the first loop (because mod(1,10)==1)

  ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
  if(.false.)then ! change false--.gt.true if you would like to use caution about having first coefficient == 0
    if(degree.lt.0) then
      write(*,*) 'Error: cmplx_laguerre2newton: degree.lt.0'
      return
    endif
    if(poly(degree+1)==zero) then
      if(degree==0) return 
      call cmplx_laguerre2newton(poly, degree-1, root, iter, success, starting_mode)
      return
    endif
    if(degree.le.1)then
      if(degree==0) then  ! we know from previous check than poly(1) not equal zero
        success=.false.
        write(*,*) 'Warning: cmplx_laguerre2newton: degree=0 and poly(1)/=0, no roots'
        return
      else
        root=-poly(1)/poly(2)
        return
      endif  
    endif
  endif
  !  end EXTREME failsafe

  j=1
  good_to_go=.false.

  mode=starting_mode  ! mode=2 full laguerre, mode=1 SG, mode=0 newton

  do ! infinite loop, just to be able to come back from newton, if more than 10 iteration there

  !------------------------------------------------------------- mode 2
  if(mode.ge.2) then  ! LAGUERRE'S METHOD
    one_nth=1d0/degree
    n_1_nth=(degree-1d0)*one_nth
    two_n_div_n_1=2d0/n_1_nth
    c_one_nth=cmplx(one_nth,0d0,RK)
  
    do i=1,MAX_ITERS  !
      faq=1d0

      ! prepare stoping criterion
      ek=abs(poly(degree+1))
      absroot=abs(root)
      ! calculate value of polynomial and its first two derivatives
      p  =poly(degree+1)
      dp =zero
      d2p_half=zero
      do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
        d2p_half=dp + d2p_half*root
        dp =p + dp*root
        p  =poly(k)+p*root    ! b_k
        ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
        ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
        ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
        ! Eq 8.
        ek=absroot*ek+abs(p)
      enddo
      abs2p=real(conjg(p)*p) !abs(p)
      iter=iter+1
      if(abs2p==0d0) return
      
      stopping_crit2=(FRAC_ERR*ek)**2
      if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
        ! do additional iteration if we are less than 10x from stopping criterion
        if(abs2p.lt.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
          return ! return immediately, because we are at very good place
        else
          good_to_go=.true. ! do one iteration more
        endif
      else
        good_to_go=.false. ! reset if we are outside the zone of the root
      endif
    
      denom=zero
      if(dp/=zero)then 
        fac_netwon=p/dp
        fac_extra=d2p_half/dp
        F_half=fac_netwon*fac_extra

        abs2_F_half=real(conjg(F_half)*F_half)
        if(abs2_F_half.le.0.0625d0)then     ! F.lt.0.50, F/2.lt.0.25
          ! go to SG method
          if(abs2_F_half.le.0.000625d0)then ! F.lt.0.05, F/2.lt.0.025
            mode=0 ! go to Newton's
          else
            mode=1 ! go to SG
          endif
        endif
        
        
        denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

        ! NEXT LINE PROBABLY CAN BE COMMENTED OUT 
        if(real(denom_sqrt).ge.0d0)then
          ! real part of a square root is positive for probably all compilers. You can 
          ! test this on your compiler and if so, you can omit this check
          denom=c_one_nth+n_1_nth*denom_sqrt
        else
          denom=c_one_nth-n_1_nth*denom_sqrt
        endif
      endif
      if(denom==zero)then !test if demoninators are .gt. 0.0 not to divide by zero
        dx=(abs(root)+1d0)*exp(cmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,RK)) ! make some random jump
      else
        dx=fac_netwon/denom
      endif
 
      newroot=root-dx
      if(newroot==root) return ! nothing changes -.gt. return
      if(good_to_go)then       ! this was jump already after stopping criterion was met
        root=newroot
        return
      endif

      if(mode/=2) then 
        root=newroot
        j=i+1    ! remember iteration index
        exit     ! go to Newton's or SG
      endif

      if(mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
        faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
        newroot=root-faq*dx ! do jump of some semi-random length (0.lt.faq.lt.1)
      endif
      root=newroot
    enddo ! do mode 2

    if(i.ge.MAX_ITERS) then
      success=.false.
      return
    endif

  endif ! if mode 2


  !------------------------------------------------------------- mode 1
  if(mode==1) then  ! SECOND-ORDER GENERAL METHOD (SG)

    do i=j,MAX_ITERS  !
      faq=1d0

      ! calculate value of polynomial and its first two derivatives
      p  =poly(degree+1)
      dp =zero
      d2p_half=zero
      if(mod(i-j,10)==0)then
        ! prepare stoping criterion
        ek=abs(poly(degree+1))
        absroot=abs(root)
        do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
          d2p_half=dp + d2p_half*root
          dp =p + dp*root
          p  =poly(k)+p*root    ! b_k
          ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
          ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
          ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
          ! Eq 8.  
          ek=absroot*ek+abs(p)
        enddo
        stopping_crit2=(FRAC_ERR*ek)**2
      else
        do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
          d2p_half=dp + d2p_half*root
          dp =p + dp*root
          p  =poly(k)+p*root    ! b_k
        enddo
      endif


      abs2p=real(conjg(p)*p) !abs(p)**2
      iter=iter+1
      if(abs2p==0d0) return


      if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967) 
        if(dp==zero) return
        ! do additional iteration if we are less than 10x from stopping criterion
        if(abs2p.lt.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
          return ! return immediately, because we are at very good place
        else
          good_to_go=.true. ! do one iteration more
        endif
      else
        good_to_go=.false. ! reset if we are outside the zone of the root
      endif
      
      if(dp==zero)then !test if demoninators are .gt. 0.0 not to divide by zero
        dx=(abs(root)+1d0)*exp(cmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,RK)) ! make some random jump
      else
        fac_netwon=p/dp
        fac_extra=d2p_half/dp
        F_half=fac_netwon*fac_extra

        abs2_F_half=real(conjg(F_half)*F_half)
        if(abs2_F_half.le.0.000625d0)then ! F.lt.0.05, F/2.lt.0.025
          mode=0 ! set Newton's, go there after jump
        endif
      
        dx=fac_netwon*(c_one+F_half)  ! SG
      endif
      
      newroot=root-dx
      if(newroot==root) return ! nothing changes -.gt. return
      if(good_to_go)then       ! this was jump already after stopping criterion was met
        root=newroot
        return
      endif

      if(mode/=1) then 
        root=newroot
        j=i+1    ! remember iteration number
        exit     ! go to Newton's
      endif

      if(mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
        faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
        newroot=root-faq*dx ! do jump of some semi-random length (0.lt.faq.lt.1)
      endif
      root=newroot



    enddo ! do mode 1

    if(i.ge.MAX_ITERS) then
      success=.false.
      return
    endif
 
  endif ! if mode 1


  !------------------------------------------------------------- mode 0
  if(mode==0) then  ! NEWTON'S METHOD

    do i=j,j+10  ! do only 10 iterations the most, then go back to full Laguerre's
      faq=1d0    

      
      ! calculate value of polynomial and its first two derivatives
      p  =poly(degree+1)
      dp =zero
      if(i==j)then ! calculate stopping crit only once at the begining
        ! prepare stoping criterion
        ek=abs(poly(degree+1))
        absroot=abs(root)
        do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
          dp =p + dp*root
          p  =poly(k)+p*root    ! b_k
          ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
          ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
          ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
          ! Eq 8.  
          ek=absroot*ek+abs(p)
        enddo
        stopping_crit2=(FRAC_ERR*ek)**2
      else        !
        do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
          dp =p + dp*root
          p  =poly(k)+p*root    ! b_k
        enddo
      endif
      abs2p=real(conjg(p)*p) !abs(p)**2
      iter=iter+1
      if(abs2p==0d0) return


      if(abs2p.lt.stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
        if(dp==zero) return 
        ! do additional iteration if we are less than 10x from stopping criterion
        if(abs2p.lt.0.01d0*stopping_crit2) then ! ten times better than stopping criterion
          return ! return immediately, because we are at very good place
        else
          good_to_go=.true. ! do one iteration more
        endif
      else
        good_to_go=.false. ! reset if we are outside the zone of the root
      endif
    
      if(dp==zero)then ! test if demoninators are .gt. 0.0 not to divide by zero
        dx=(abs(root)+1d0)*exp(cmplx(0d0,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,RK)) ! make some random jump
      else
        dx=p/dp
      endif
      
      newroot=root-dx
      if(newroot==root) return ! nothing changes -.gt. return
      if(good_to_go)then
        root=newroot
        return
      endif

      ! this loop is done only 10 times. So skip this check
      !if(mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
      !  faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
      !  newroot=root-faq*dx ! do jump of some semi-random length (0.lt.faq.lt.1)
      !endif
      root=newroot

    enddo ! do mode 0 10 times

    if(iter.ge.MAX_ITERS) then
      ! too many iterations here  
      success=.false.
      return
    endif
    mode=2 ! go back to Laguerre's. This happens when we were unable to converge in 10 iterations with Newton's 

  endif ! if mode 0

  enddo ! end of infinite loop

  !------------------------------------------------------------- 
  success=.false.
end


!-------------------------------------------------------------------!
! _9_                     SOLVE_QUADRATIC_EQ                        !
!-------------------------------------------------------------------!
subroutine solve_quadratic_eq(x0,x1,poly)
  ! Quadratic equation solver for complex polynomial (degree=2)
  implicit none
  complex(kind=8), intent(out) :: x0, x1
  complex(kind=8), dimension(*), intent(in) :: poly ! coeffs of the polynomial
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 
  complex(kind=8) :: a, b, c, b2, delta

  complex(kind=8) :: val, x
  integer :: i

  a=poly(3)
  b=poly(2)
  c=poly(1)
  ! quadratic equation: a z^2 + b z + c = 0  

  b2=b*b
  delta=sqrt(b2-4d0*(a*c))
  if( real(conjg(b)*delta).ge.0d0 )then  ! scallar product to decide the sign yielding bigger magnitude
    x0=-0.5d0*(b+delta)
  else
    x0=-0.5d0*(b-delta)
  endif
  if(x0==cmplx(0d0,0d0,8))then
    x1=cmplx(0d0,0d0,8)
  else ! Viete's formula
    x1=c/x0
    x0=x0/a
  endif


  if(.false.)then  ! print the results

    x=x0
    val=poly(3)
    do i=2,1,-1
      val=val*x+poly(i)
    enddo
    write(*,'(2f19.15,a3,2f19.15)') x,' -.gt.',val

    x=x1
    val=poly(3)
    do i=2,1,-1
      val=val*x+poly(i)
    enddo
    write(*,'(2f19.15,a3,2f19.15)') x,' -.gt.',val

  endif

end

!-------------------------------------------------------------------!
! _10_                    SOLVE_CUBIC_EQ                            !
!-------------------------------------------------------------------!
subroutine solve_cubic_eq(x0,x1,x2,poly)
  ! Cubic equation solver for complex polynomial (degree=3)
  ! http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
  implicit none
  complex(kind=8), intent(out) :: x0, x1, x2
  complex(kind=8), dimension(*), intent(in) :: poly ! coeffs of the polynomial
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3             4
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + poly(4) x^3
  complex(kind=8), parameter :: zeta =cmplx(-0.5d0, 0.8660254037844386d0, 8)  ! sqrt3(1)
  complex(kind=8), parameter :: zeta2=cmplx(-0.5d0,-0.8660254037844386d0, 8)  ! sqrt3(1)**2
  real(kind=8), parameter    :: third=0.3333333333333333d0                    ! 1/3
  complex(kind=8) :: s0, s1, s2
  complex(kind=8) :: E1 ! x0+x1+x2
  complex(kind=8) :: E2 ! x0x1+x1x2+x2x0
  complex(kind=8) :: E3 ! x0x1x2
  complex(kind=8) :: A, B, a_1, E12
  complex(kind=8) :: delta, A2

  complex(kind=8) :: val, x
  integer :: i

  a_1=poly(4)**(-1)
  E1=-poly(3)*a_1
  E2=poly(2)*a_1
  E3=-poly(1)*a_1
  
  s0=E1
  E12=E1*E1
  A=2d0*E1*E12-9d0*E1*E2+27d0*E3  ! =  s1^3 + s2^3
  B=E12-3d0*E2                    !  = s1 s2
  ! quadratic equation: z^2-Az+B^3=0  where roots are equal to s1^3 and s2^3
  A2=A*A
  delta=sqrt(A2-4d0*(B*B*B))
  if( real(conjg(A)*delta).ge.0d0 )then ! scallar product to decide the sign yielding bigger magnitude
    s1=(0.5d0*(A+delta))**third
  else
    s1=(0.5d0*(A-delta))**third
  endif
  if(s1==cmplx(0d0,0d0,8))then
    s2=cmplx(0d0,0d0,8)
  else
    s2=B/s1
  endif


  x0=third*(s0+s1+s2)
  x1=third*(s0+s1*zeta2+s2*zeta )
  x2=third*(s0+s1*zeta +s2*zeta2)


  if(.false.)then  ! print the results

    x=x0
    val=poly(4)
    do i=3,1,-1
      val=val*x+poly(i)
    enddo
    write(*,'(2f19.15,a3,2f19.15)') x,' -.gt.',val

    x=x1
    val=poly(4)
    do i=3,1,-1
      val=val*x+poly(i)
    enddo
    write(*,'(2f19.15,a3,2f19.15)') x,' -.gt.',val

    x=x2
    val=poly(4)
    do i=3,1,-1
      val=val*x+poly(i)
    enddo
    write(*,'(2f19.15,a3,2f19.15)') x,' -.gt.',val


  endif
end


!-------------------------------------------------------------------!
! _11_                    DIVIDE_POLY_1                             !
!-------------------------------------------------------------------!
subroutine divide_poly_1(polyout, remainder, p, polyin, degree)
  ! Subroutine will divide polynomial 'polyin' by (x-p)
  ! results will be returned in polynomial 'polyout' of degree-1
  ! The remainder of the division will be returned in 'remainder'
  !
  ! You can provide same array as 'polyin' and 'polyout' - this
  ! routine will work fine, though it will not set to zero the 
  ! unused, highest coefficient in the output array. You just have
  ! remember the proper degree of a polynomial.
  implicit none
  integer, intent(in) :: degree
  complex(kind=8), dimension(degree), intent(out) :: polyout
  complex(kind=8), intent(out) :: remainder
  complex(kind=8), intent(in) :: p
  complex(kind=8), dimension(degree+1), intent(in) :: polyin ! coeffs of the polynomial
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  integer :: i
  complex(kind=8) :: coef, prev

  coef=polyin(degree+1)
  polyout=polyin(1:degree)
  do i=degree,1,-1
    prev=polyout(i)
    polyout(i)=coef
    coef=prev+p*coef
  enddo
  remainder=coef
  return
end


!-------------------------------------------------------------------!
! _12_                    EVAL_POLY                                 !
!-------------------------------------------------------------------!
complex(kind=8) function eval_poly(x, poly, degree, errk)
  ! Evaluation of the complex polynomial 'poly' of a given degree 
  ! at the point 'x'. This routine calculates also the simplified
  ! Adams' (1967) stopping criterion. ('errk' should be multiplied 
  ! by 2d-15 for double precision, real*8, arithmetic)
  implicit none
  complex(kind=8), intent(in) :: x
  integer, intent(in) :: degree
  real(kind=8), intent(out) :: errk ! 
  complex(kind=8), dimension(degree+1), intent(in) :: poly ! coeffs of the polynomial
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  integer :: i
  complex(kind=8) :: val
  real(kind=8) :: absx

  ! prepare stoping criterion
  errk=abs(poly(degree+1))
  val=poly(degree+1)
  absx=abs(x)
  do i=degree,1,-1  ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
    val=val*x+poly(i)
    ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
    ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
    ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
    ! Eq 8.
    errk=errk*absx+abs(val)
  enddo
  eval_poly=val

  ! if(abs(val).lt.2d-15*errk) return  ! (simplified a little Eq. 10 of Adams 1967)


  return
end

!-------------------------------------------------------------------!
! _13_                    MULTIPLY_POLY_1                           !
!-------------------------------------------------------------------!
subroutine multiply_poly_1(polyout, p, polyin, degree)
  ! Subroutine will multiply polynomial 'polyin' by (x-p)
  ! results will be returned in polynomial 'polyout' of degree+1
  !
  ! You can provide same array as 'polyin' and 'polyout' - this
  ! routine will work fine.
  implicit none
  integer, intent(in) :: degree  ! OLD degree, new will be +1
  complex(kind=8), dimension(degree+2), intent(out) :: polyout
  complex(kind=8), intent(in) :: p
  complex(kind=8), dimension(degree+1), intent(in) :: polyin ! coeffs of the polynomial
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  integer :: i

  polyout(1:degree+1)=polyin(1:degree+1) ! copy

  polyout(degree+2)=polyout(degree+1)
  do i=degree+1,2,-1
    polyout(i)=polyout(i-1)-polyout(i)*p
  enddo
  polyout(1)=-polyout(1)*p
  return
end


!-------------------------------------------------------------------!
! _14_                    CREATE_POLY_FROM_ROOTS                    !
!-------------------------------------------------------------------!
subroutine create_poly_from_roots(poly, degree, a, roots)
  ! Routine will build polynomial from a set of points given in 
  ! the array 'roots'. These points will be zeros of the resulting 
  ! polynomial.
  !
  ! poly - is an array of polynomial cooefs, length = degree+1, poly(1) is constant 
  !             1              2             3
  !        poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
  ! degree - is and integer denoting size of the 'roots' array
  ! a - gives the leading coefficient of the resutling polynomial
  ! roots - input array of points, size=degree
  !
  ! This subroutine works, but it is not optimal - it will work  
  ! up to a polynomial of degree~50, if you would like to have 
  ! more robust routine, up to a degree ~ 2000, you should 
  ! split your factors into a binary tree, and then multiply 
  ! leafs level by level from down up with subroutine like: 
  ! multiply_poly for arbitraty polynomial multiplications
  ! not multiply_poly_1.
  !
  implicit none
  integer, intent(in) :: degree
  complex(kind=8), dimension(degree+1), intent(out) :: poly
  complex(kind=8), intent(in) :: a
  complex(kind=8), dimension(degree), intent(in) :: roots
  !
  integer :: i
  !
  poly=cmplx(0d0,0d0,8)
  poly(1)=a  ! leading coeff of the polynomial
  do i=1,degree
    call multiply_poly_1(poly, roots(i), poly, i-1)
  enddo
end




!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------
end module rpa_laplace
!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------





!! program testint
!! 
!!  use integration
!!  implicit none
!! 
!!  double precision, external  :: f, ff
!! 
!!  double precision :: a=0d0, b=1d0
!!  double precision :: y, yy
!! 
!!   call integral_setmethod(5)
!!   call integral_setaccuracy(1d-11)
!!   call intprint(0)
!! 
!!  
!!   y = integral(f,a,b)
!!  
!!   write(6,*)
!!   write(6,*)' integral = ', y, '+-', integral_getaccuracy()
!!  
!!   y = integral(f,a,b,5)
!!  
!!   write(6,*)
!!   write(6,*)' integral = ', y, '+-', integral_getaccuracy()
!!  
!!  
!!   y = integral(f,a,b,1000,1d-8)
!!  
!!   write(6,*)
!!   write(6,*)' integral = ', y, '+-', integral_getaccuracy()
!!  
!!   yy = 0
!!   y = integral(f,a,b,1000,1d-8,yy)
!!  
!!   write(6,*)
!!   write(6,*)' integral = ', y, '+-', integral_getaccuracy(),yy
!!  
!!    yy = 0
!!    y = integral(ff,a,b)
!!   
!!    write(6,*)
!!    write(6,*)' integral = ', y, '+-', integral_getaccuracy()
!!  
!!    yy = 0
!!    y = integral(ff,a,b,300,1d-15,yy)
!!   
!!    write(6,*)
!!    write(6,*)' integral = ', y, '+-', integral_getaccuracy(),yy
!!  
!! 
!!    call integral_setmininterval(1d-25)
!! 
!!    yy = 0
!!    y = integral(ff,a,b,10000,1d-20,yy)
!!   
!!    write(6,*)
!!    write(6,*)' integral = ', y, '+-', integral_getaccuracy(),yy
!!  
!! 
!! end program testint
!! 
!!  
!! 
!! 
!! 
!! double precision function f(x)
!!  double precision :: x
!!  f = 1/sqrt(sqrt(x))
!! end function f
!! 
!! double precision function ff(t)
!!   use integration
!!   double precision :: t
!!   double precision, external :: f
!!   ff = integral(f,0d0,t)
!! end function ff
!! 
!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------
!! 
!! 
!! 
!! program testing
!! 
!! use rpa_laplace
!! use integration
!! 
!! implicit none
!! 
!! complex(kind=XPREC) :: s, S1, S2, S11, S11i
!! real   (kind=XPREC) :: t, yr, yi, range, yyr, yyi 
!! real   (kind=XPREC) :: a = -3000d0
!! real   (kind=XPREC) :: b =  3000d0
!! 
!! 
!! 
!! ! real   (kind=XPREC) :: limit_scale_factor=100d0  ! determines limit of integration
!! real   (kind=XPREC) :: limit_scale
!! 
!! ! real   (kind=XPREC) :: rlow = 1d0/2000d0         ! lowest expected detectable rate
!! real   (kind=XPREC) :: rmin
!! 
!! 
!! integer :: i, j, ni = 10
!! 
!! 
!!  xS00   = 20
!!  xa1    = 0.3
!!  xa2    = 0.4
!!  xr1    = 0.02
!!  xr2    = 1
!!  xr3    = 12
!!  xScc00 = 1
!!  xb1    = 0.35
!!  xb1    = 0.35
!!  xg1    = 2
!!  xg2    = 5
!!  xg3    = 7
!!  xphi   = 0.3
!! 
!!  is_sel   = 1        ! selector 1: S11, 2: S12, 3: S12, 4: S22 
!!  phi1     = 0.3      ! volume fraction of polymer component 1
!!  phi2     = 0.1      ! volume fraction of polymer component 2
!!  Scc00    = xScc00   ! unperturbed structure factor S(Q) of "matrix" polymers
!!  S0011    = xS00     ! unperturbed structure factor S(Q) of polymer 1
!!  S0022    = xScc00*2   ! unperturbed structure factor S(Q) of polymer 2
!!  nexpcc   = 3        ! number of exp-functions to describe background
!!  nexp1    = 3        ! number of exp-functions to describe component1
!!  nexp2    = 3        ! number of exp-functions to describe component2
!!  aexp_cc(1:3) =[xb1,xb2,1-xb1-xb2] ! amplitude coeffs for laplace-exp representation of "matrix"
!!  rexp_cc(1:3) =[xg1,xg2,xg3]       ! rate      coeffs for laplace-exp representation of "matrix"
!!  aexp_s1(1:3) =[xa1,xa2,1-xa1-xa2] ! amplitude coeffs for laplace-exp representation of polymer 1
!!  rexp_s1(1:3) =[xr1,xr2,xr3]       ! rate      coeffs for laplace-exp representation of polymer 1
!!  aexp_s2      = aexp_cc            ! amplitude coeffs for laplace-exp representation of polymer 2
!!  rexp_s2      = rexp_cc            ! rate      coeffs for laplace-exp representation of polymer 2
!!  
!! 
!!  xil         =  0.0001d0    ! distance of path for inv-laplace integration
!!  epap        =  1d-9
!! 
!!   call integral_setmethod(4)
!!   call integral_setaccuracy(1d-10)
!!   call intprint(0)
!! !
!!   
!! 
!! 
!! !!  do i=1,100
!! !! 
!! !!  s = cmplx(0.03d0*i,1d0)
!! !! 
!! !! 
!! !!  S1 = Ss4     (s, xS00, xa1, xa2, xr1, xr2, xr3, xScc00, xb1, xb2, xg1, xg2, xg3, xphi)
!! !!  S2 = Ss_test4(s, xS00, xa1, xa2, xr1, xr2, xr3, xScc00, xb1, xb2, xg1, xg2, xg3, xphi)
!! !!   
!! !!  write(6,*) i
!! !!  write(6,*) S1
!! !!  write(6,*) S2
!! !! 
!! !!  enddo
!!   
!!   rmin = max(rlow, min( minval(rexp_s1(1:nexp1)), minval(rexp_s1(1:nexp2)),  minval(rexp_cc(1:nexpcc)) ))
!! 
!! 
!!   limit_scale = limit_scale_factor / rmin
!! 
!!   write(6,*)"rmin...",rmin,limit_scale
!! 
!! 
!!  do i=1,200
!!   t = 0.001d0 * i
!!   S1   = InvLaplace4a(t, xS00, xa1, xa2, xr1, xr2, xr3, xScc00, xb1, xb2, xg1, xg2, xg3, xphi )
!!   S11  = InvLaplace2d(t, S0011, xa1, xa2, xr1, xr2, xr3, S0022, Scc00, xb1, xb2, xg1, xg2, xg3, phi1, phi2 )
!! 
!!   t_param = t
!! 
!!   b =  max(limit_scale/t, limit_scale)
!!   a = -b
!! 
!!   yr = integral( Ss_kernel_Re ,a,b,10000,1d-8)
!! 
!!   yi = 2*integral( Ss_kernel_Re ,0d0,b,10000,1d-8)
!! 
!! 
!! !  yyr = integral( Ss_kernel_2d ,a,b,10000,1d-8)
!! 
!!   yyr = st_rpa( t, 1,1)
!! 
!!   yyi = 2*integral( Ss_kernel_2d ,0d0,b,10000,1d-8)
!! 
!!    
!! 
!! !!  do j=1,ni
!! !!    range=3*(b-a)/ni
!! !!    yi = yi + integral( Ss_kernel_Re ,3*a+(j-1)*range,3*a+j*range,10000,1d-9)
!! !!  enddo
!! 
!! !  write(6,*)"a:",i,S1, yr, yi
!!   write(6,*)"b:",i,S11, yyr, yyi
!!   
!! 
!!  enddo
!! 
!! end program testing
!! 
