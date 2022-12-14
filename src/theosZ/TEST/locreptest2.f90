PROGRAM locreptest
 use  ADAPTIVE_INTEGRAL
 implicit none
 
 double precision :: q, t, A, B, L
 double precision :: nu
 double precision :: y, y2, y3, y4, y5
 double precision :: t0, tinf
 double precision :: eps, erracc
 integer          :: maxit
 integer          :: M
 integer          :: nt
 integer          :: i


 double precision :: local_reptation, local_reptation_lr
! double precision, external :: lr_kernel,lrring_kernel


 eps    = 1d-9
 maxit  = 1000

 nu     = 0.5d0
 M      = 100

 q      = 3d0
 A      = 1d0
 B      = 10d0
 L      = 10d0
 nt     = 30
 t0     = 0.001d0
 tinf   = 10000d0

 

 do i=0,nt
   t = t0 * exp( (log(tinf/t0)/nt) * i )
   y = local_reptation(q, t, B, L) 
!   y2= adapint (lr_kernel, 0d0, L, eps, maxit, erracc) 
!   y3= adapint (lrring_kernel, 0d0, L, eps, maxit, erracc)    !! still mission periodic
   y4= local_reptation_lr(q, t, B, L, nu, M, 0) 
   y5= local_reptation_lr(q, t, B, L, nu, M, 1) 

   write(*,*) t, y,y4,y5
 enddo

 

END PROGRAM locreptest


function local_reptation_lr(q, t, B, L, nu, M, mode) result(val)
 use  ADAPTIVE_INTEGRAL

 implicit none
 double precision, intent(in)  :: q
 double precision, intent(in)  :: t
 double precision, intent(in)  :: B
 double precision, intent(in)  :: L
 double precision, intent(in)  :: nu
 integer         , intent(in)  :: M
 integer         , intent(in)  :: mode   ! 0=std, 1=ring

 double precision :: val
 double precision :: A

 double precision :: eps 
 double precision :: erracc
 integer          :: maxit

 A     =   1d0
 eps   =   1d-9
 maxit = 1000

 if(mode == 0) then
   val = adapint (lr_kernel    , 0d0, L, eps, maxit, erracc) 
 else
   val = adapint (lrring_kernel, 0d0, L, eps, maxit, erracc) 
 endif
  


CONTAINS

  function lr_kernel(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)

!   val = (A + B*exp(-u**2/(4d0*t))/(2d0*sqrt(Pi*t)))*exp(-q**2*abs(u)**(2*nu)/6d0)
    val = (A + B* diff_1d(u)                        )*exp(-q**2*abs(u)**(2*nu)/6d0)
    val = 2*val*(L-u)/L

  end function lr_kernel



  function lrring_kernel(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)

    val = (A + B*periodic_diff(u))*exp(-q**2*(abs(u)*(1d0-abs(u)/L))**(2*nu)/6d0)
    val = 2*val*(L-u)/L

  end function lrring_kernel


  function periodic_diff(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)
    integer                      :: k

    val = 0d0
    do k=-M,M
      val = val + cos(2*Pi*u*k/L) * exp(-k*k*t*2d0/5d0)
    enddo 
    val = val / (L)

  end function periodic_diff


  function diff_1d(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)
    integer                      :: k

    val = exp(-u**2/(4d0*t))/(2d0*sqrt(Pi*t))

  end function diff_1d

end function local_reptation_lr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function local_reptation(q, t, B, L) result(val)
    implicit none
    double precision, intent(in)   :: q, t, B, L
    double precision, parameter    :: A = 1d0
    double precision               :: val
    double precision, parameter    :: sqp = sqrt(4*atan(1d0))

    double precision               :: ec1, ec2, dec

    ec1 = erfc((q ** 2 * t + 0.3D1 * L) * t ** (-0.1D1 / 0.2D1) / 0.6D1)
    ec2 = erfc(sqrt(t) * q ** 2 / 0.6D1)

    if( isnan(ec1) ) ec1 = 1d0 - erf((q ** 2 * t + 0.3D1 * L) * t ** (-0.1D1 / 0.2D1) / 0.6D1)
    if( isnan(ec2) ) ec2 = 1d0 - erf(sqrt(t) * q ** 2 / 0.6D1)

    dec =   (-ec1 + ec2 )

    if( isnan(dec) ) dec = 0

      val= 0.72D2 * (sqrt(t) * q ** 4 * &
           exp((-0.2D1 * L * q ** 2 * t - 0.3D1 * L ** 2) &
           / t / 0.12D2) / 0.36D2 + sqrt(0.3141592653589793D1) &
           * (q ** 2 * t / 0.3D1 + L) * q ** 4 * &
           exp(t * q ** 4 / 0.36D2) * &
          ( dec ) / 0.72D2 - sqrt(t) * q ** 4 / 0.36D2) &
           * B / q ** 4 * 0.3141592653589793D1 ** (-0.1D1 / 0.2D1) / L + &
           0.72D2 * (A * &
           exp(-q ** 2 * L / 0.6D1)  &
         * sqrt(0.3141592653589793D1) &
         + (A * L * q ** 2 / 0.6D1 - A) * sqrt(0.3141592653589793D1)) * &
         0.3141592653589793D1 ** (-0.1D1 / 0.2D1) / q ** 4 / L

  end function local_reptation 
