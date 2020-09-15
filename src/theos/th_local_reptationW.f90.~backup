 FUNCTION th_locrepw(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrepw
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor                                                                       
     double precision :: b          ! fluctuation intensity (relative)                                                
     double precision :: a          ! length scale                                                                    
     double precision :: tau        ! timescale                                                                       
     double precision :: lz         ! total length                                                                    
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: local_reptation_lw, sqt, sqt0
     double precision :: teps = 1d-6
 
     double precision :: nu
     double precision :: W, wfac
     double precision :: Ne, lseg   !! ??

     integer          :: m, mode
     

     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locrepw'
       nparx =        10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrepw = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " but with finite summation of integrals and lenthscale," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters" //cr//parspace//&
                                " tau from W and Ne: tau = (Ne*Pi)**2/W; Ne = (a/lseg)**2 "
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'b       '  ! fluctuation intensity (relative)                                                
        parnam ( 3) = 'a       '  ! length scale                                                                    
        parnam ( 4) = 'w       '  ! Rouserate                                                                    
        parnam ( 5) = 'lz      '  ! total length                                                                    
        parnam ( 6) = 'teps    '  ! near zero for t=0 evaluation 
        parnam ( 7) = 'nu      '  ! default 0.5
        parnam ( 8) = 'm       '  ! sharpness of transition 
        parnam ( 9) = 'mode    '  ! 0=std, 1=ring, 2= ring with 1d-diff, 3= std analytic
        parnam (10) = 'wfac    '  ! W Modification factor 
      
                                                                         
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fluctuation intensity (relative)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 5,idesc) = "total length" !//cr//parspace//&
        th_param_desc( 6,idesc) = "teps near zero for t=0 evaluation sor S(Q,t=0) evaluation" !//cr//parspace//&
        th_param_desc( 7,idesc) = "chain expansion exponent defaut=0.5" !//cr//parspace//&
        th_param_desc( 8,idesc) = "diffusion mode summation number (only for mode=1)" !//cr//parspace//&
        th_param_desc( 9,idesc) = "mode = 0: std locrep expr via numerical integration (allows nu .ne. 1/2)," //cr//parspace//&
                                  "mode = 1: ring variant," //cr//parspace//&
                                  "mode = 2: ring with nonperiodic diff" //cr//parspace// &
                                  "mode = 3: analytical form of locrep (original, assumes nu=0.5)"
        th_param_desc(10,idesc) = "Factor to modify W, W ==> wfac*W" !//cr//parspace//& 
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
        th_file_param(  2,idesc) = "lseg     > segment length used to compute Ne and from that tau(W)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(1,idesc)   = "a2sqt = a**2/sqrt(tau) "
! 
        th_locrepw = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =         pa( 1)
      b        =   abs ( pa( 2) )
      a        =   abs ( pa( 3) )
      W        =   abs ( pa( 4) )
      lz       =   abs ( pa( 5) )
      teps     =         pa( 6)
      nu       =   abs ( pa(7)  )
      M        =   nint( pa(8)  )
      mode     =   nint( pa(9)  )
      wfac     =   abs ( pa(10) )
      if (teps <= 0d0) teps = 1d-6
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: lseg-value
      xh =     3.5d0
      call parget('lseg    ',xh,iadda,ier)
      if(ier .ne. 0) Write(*,*)"WARNING: locrepw: parameter lseg in data is missing!! ", &
                                           "         assuming lseg = ", lseg
      lseg     = xh

      if(wfac <=0d0) wfac = 1d0

      Ne       = (a/lseg)**2
      tau      = (Ne**2 * pi**2)  / (wfac * W) 

      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     sqt0 =  local_reptation_lw( q*a, teps,  b, lz, nu, M, mode)
     sqt  =  local_reptation_lw( q*a, t/tau, b, lz, nu, M, mode)
     th_locrepw = ampli * sqt/sqt0


!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda,ier)
                   
 end function th_locrepw




function local_reptation_lw(q, tx, B, L, nu, M, mode) result(val)
 use  ADAPTIVE_INTEGRAL

 implicit none
 double precision, intent(in)  :: q
 double precision, intent(in)  :: tx
 double precision, intent(in)  :: B
 double precision, intent(in)  :: L
 double precision, intent(in)  :: nu
 integer         , intent(in)  :: M
 integer         , intent(in)  :: mode   ! 0=std, 1=ring

 double precision :: val
 double precision :: A, La, Lb, t
 double precision, parameter :: tzero = 1d-6

 double precision :: eps 
 double precision :: erracc
 integer          :: maxit


 t     = max(tzero,tx)
 A     =   1d0
 eps   =   1d-8
 maxit = 10000


 La    = 10d0 * sqrt(t)
 Lb    = L-2*La
 if(Lb <= 0d0) then
  La = L/2
  Lb = 0
 endif


 if(mode == 0) then
   val =       adapint (lr_kernel    , 0d0  , La,    eps, maxit, erracc) 
   val = val + adapint (lr_kernel    , La   , L-La,  eps, maxit, erracc) 
   val = val + adapint (lr_kernel    , L-La , L,     eps, maxit, erracc) 
 else if(mode == 1) then
   val =       adapint (lrring_kernel    , 0d0  , La,    eps, maxit, erracc) 
   val = val + adapint (lrring_kernel    , La   , L-La,  eps, maxit, erracc) 
   val = val + adapint (lrring_kernel    , L-La , L,     eps, maxit, erracc) 
 else if(mode == 2) then
   val =       adapint (lrringd_kernel    , 0d0  , La,    eps, maxit, erracc) 
   val = val + adapint (lrringd_kernel    , La   , L-La,  eps, maxit, erracc) 
   val = val + adapint (lrringd_kernel    , L-La , L,     eps, maxit, erracc) 
 else if(mode == 3) then
   if( abs(nu-0.5d0) > epsilon(1.0) ) Write(*,*)"WARNING: locrepw mode = 3",mode," ignores nu= ",nu
   val = local_reptation(q, t, B, L)
 else
   write(*,*)"WARNING: locrepw mode = ",mode," NOT KNOWN! "
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


  function lrringd_kernel(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)

    val = (A + B*diff_1d(u))*exp(-q**2*(abs(u)*(1d0-abs(u)/L))**(2*nu)/6d0)
    val = 2*val*(L-u)/L

  end function lrringd_kernel


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

end function local_reptation_lw


