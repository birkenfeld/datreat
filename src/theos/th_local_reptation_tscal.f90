 FUNCTION th_locrept(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrept
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
     double precision :: tscale     ! total length                                                                    
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: local_reptation, sqt, sqt0
     double precision, parameter   :: teps = 1d-6
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locrept'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrept = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " but with finite summation of integrals and lenthscale," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters"
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'b       '  ! fluctuation intensity (relative)                                                
        parnam ( 3) = 'a       '  ! length scale                                                                    
        parnam ( 4) = 'tau     '  ! timescale                                                                       
        parnam ( 5) = 'lz      '  ! total length                                                                    
        parnam ( 6) = 'tscale  '  ! time streching factor                                                          
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fluctuation intensity (relative)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 5,idesc) = "total length" !//cr//parspace//&
        th_param_desc( 6,idesc) = "time streching factor" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(1,idesc)   = "a2sqt = a**2/sqrt(tau) "
! 
        th_locrept = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      b        =   abs( pa( 2) )
      a        =   abs( pa( 3) )
      tau      =   abs( pa( 4) )
      lz       =   abs( pa( 5) )
      tscale   =   abs( pa( 6) )
      if(tscale == 0d0) tscale = 1d0
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      t        = x * tscale
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     sqt0 =  local_reptation( q*a, teps,  b, lz)
     sqt  =  local_reptation( q*a, t/tau, b, lz)
     th_locrept = ampli * sqt/sqt0

!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda,ier)
                   
 end function th_locrept

! subroutines and functions entered here are private to this theory and share its variables 
 
  function local_reptation_t(q, t, B, L) result(val)
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

  end function local_reptation_t
