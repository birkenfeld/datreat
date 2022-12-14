program testlocrep
implicit none

double precision :: q = 0.1d0
double precision :: t = 1d0
double precision :: a = 5d0
double precision :: W = 10d0
double precision :: n = 1d4
double precision :: ne = 60d0
double precision :: b = 1d0
double precision :: t0 = 1d0
double precision :: tw = 1d0
double precision :: tww 
double precision :: ww0 = 0.02d0
integer          :: modelr = 0


double precision :: y1, y2


integer :: i, j


! do j=1,10

! a = 3*j

tww = tw
do i=1,100
  t = 1d-8*1.2d0**i
  y1 = local_reptationdr(q,  t, a, W, n, ne, b,t0,tw,ww0) 
  y2 = local_reptationdr2(q, t, a, W, n, ne, b, modelr) 

  write(*,*) a, t, y1 , y2, y1-y2
  
enddo

! enddo




contains
 
  function local_reptationdr(q, t, a, W, n, ne, b,t0,tw,ww0) result(val)
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: a        ! step (entanglement?) length
    double precision, intent(in)   :: W        ! (Rouse?) rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments/entanglement
    double precision, intent(in)   :: b        ! b-fluctuation intensity
    double precision, intent(in)   :: t0       ! w cross over time
    double precision, intent(in)   :: tw       ! w cross-over width
    double precision, intent(in)   :: ww0      ! transition depth
    double precision               :: wt
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)

    double precision :: T1, T2, T3
    wt = w*(1d0 - (1d0-ww0)/(1d0+exp((t-t0)/tw)))

    T1 = 72d0 * (exp(-q**2 * n * a**2 / 6d0) + q**2 * n * a**2 / 6d0 - 1d0) / q**4 / a**4
    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*wt*t))-1)*2/3*ne*sqrt(wt*t/Pi)+ &
        ne*(n/3+a**2*q**2*wt*t/9)*exp(a**4*q**4*wt*t/36)* &
         (erfc(a**2*q**2*sqrt(wt*t)/6)-erfc(n/(2*sqrt(wt*t))+a**2*q**2*sqrt(wt*t)/6)))
    T3 = 72d0/((q**4*a**4))*(exp(-(q**2*n*a**2)/6d0)+(q**2*n*a**2)/6d0-1d0)+b*ne*n/3d0

    val = (T1+T2)/T3


  end function local_reptationdr







  function local_reptationdr2(q, t, a, W0, n, ne, b, modelr) result(val)   !! with series expansion
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: a        ! step (entanglement?) length
    double precision, intent(in)   :: W0        ! (Rouse?) rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments/entanglement
    double precision, intent(in)   :: b        ! b-fluctuation intensity
    integer         , intent(in)   :: modelr   ! full deGennes 0, or only T2 (locrep) 1
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)

    double precision :: T1, T2, T3, ec1, ec2, dec
    double precision :: x, xa
    double precision, parameter :: xlimser = 15d0
    double precision :: W

    W = w0*(1 - (1-ww0)/(1+exp((t-t0)/tww)))

    T1 = 72d0 * (exp(-q**2 * n * a**2 / 6d0) + q**2 * n * a**2 / 6d0 - 1d0) / q**4 / a**4
!    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
!        ne*(n/3+a**2*q**2*w*t/9)*exp(a**4*q**4*w*t/36)* &
!         (erfc(a**2*q**2*sqrt(w*t)/6)-erfc(n/(2*sqrt(w*t))+a**2*q**2*sqrt(w*t)/6)))
    T3 = 72d0/((q**4*a**4))*(exp(-(q**2*n*a**2)/6d0)+(q**2*n*a**2)/6d0-1d0)+b*ne*n/3d0
!
!    val = (T1+T2)/T3

    x = sqrt(a**4*q**4*w*t/36)
    xa= n/(2*sqrt(w*t))
    
    if( abs(x) < xlimser ) then
      dec = exp(x**2) *(erfc(x)-erfc(x+xa))
    else
      dec =  1d0/(sqrt(Pi)*x)-1d0/(2*sqrt(Pi)*x**3)+3d0/(4*sqrt(Pi)*x**5) &
          - (x**4-xa*x**3+(xa**2-1d0/2d0)*x**2+(-xa**3+3d0/2d0*xa)*x+xa**4-3*xa**2+3d0/4d0) &
           * exp(-xa**2)*exp(-2*xa*x)/(sqrt(Pi)*x**5)
    endif


    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
        ne*(n/3+a**2*q**2*w*t/9) &
        * dec )

    if(modelr == 1) then
      val = T2
    else
      val = (T1+T2)/T3
    endif


  end function local_reptationdr2



end program testlocrep
 
