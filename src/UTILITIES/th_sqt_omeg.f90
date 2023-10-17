 FUNCTION th_sqt_omeg(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  combined s(q,omega) and s(q,t) data from an s(q,t) model
! 
      use theory_description 
      implicit none 
      real    :: th_sqt_omeg
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
     double precision :: intensit   ! prefactor                                                                       
     double precision :: tau0       ! KWW time-constant, prefactor in front of q-dependence                           
     double precision :: beta       ! streched exp, prefactor in front of-q-dependence                                
     double precision :: epsilon    ! accuracy parameter for FT-integrations (DO NOT FIT)                             
     double precision :: omega0     ! omega scale zero shift                                                          
     double precision :: u_sqr      ! <u**2> value for Debye-Waller-Factor                                            
     double precision :: j0         ! jump length (if applicable)                                                     
     double precision :: beta0      ! beta offset                                                                     
     double precision :: qexp_t0    ! q-exponent for tau0                                                             
     double precision :: qexp_bet   ! beta-exponent                                                                   
     double precision :: bkgr       ! constant background                                                             
     double precision :: n_mg       ! fraction of hydrogen in side chains                                             
! the recin parameter representation 
     double precision :: q          ! q-value    default value                                                        
     double precision :: nse        ! explict switch between sqt and sqomega (if nse=1 ==> sqt)                       
     double precision :: bk1level   ! resolution description bgr-levwl                                                
     double precision :: bk1slope   ! resolution description slope of bgr                                             
     double precision :: ga1inten   ! resolution description 1. gaussian amplitude                                    
     double precision :: ga1width   ! resolution description 1. gaussian width                                        
     double precision :: ga1cente   ! resolution description 1. gaussian center                                       
     double precision :: ga2inten   ! resolution description 2. gaussian amplitude                                    
     double precision :: ga2width   ! resolution description 2. gaussian width                                        
     double precision :: ga2cente   ! resolution description 2. gaussian center                                       
     double precision :: ga3inten   ! resolution description 3. gaussian amplitude                                    
     double precision :: ga3width   ! resolution description 3. gaussian width                                        
     double precision :: ga3cente   ! resolution description 3. gaussian center                                       
     double precision :: ga4inten   ! resolution description 4. gaussian amplitude                                    
     double precision :: ga4width   ! resolution description 4. gaussian width                                        
     double precision :: ga4cente   ! resolution description 4. gaussian center                                       
     double precision :: ga5inten   ! resolution description 5. gaussian amplitude                                    
     double precision :: ga5width   ! resolution description 5. gaussian width                                        
     double precision :: ga5cente   ! resolution description 5. gaussian center                                       
     double precision :: ga6inten   ! resolution description 6. gaussian amplitude                                    
     double precision :: ga6width   ! resolution description 6. gaussian width                                        
     double precision :: ga6cente   ! resolution description 6. gaussian center                                       
     double precision :: ga7inten   ! resolution description 7. gaussian amplitude                                    
     double precision :: ga7width   ! resolution description 7. gaussian width                                        
     double precision :: ga7cente   ! resolution description 7. gaussian center                                       
     double precision :: ga8inten   ! resolution description 8. gaussian amplitude                                    
     double precision :: ga8width   ! resolution description 8. gaussian width                                        
     double precision :: ga8cente   ! resolution description 8. gaussian center                                       
! the reout parameter representation 
     double precision :: tau0_q     ! scaled: str_tau0 * (1+jlen*qz**qexpt0)/(qz**qexpt0)                             
     double precision :: beta_q     ! scaled: str_beta * (qz**qexpbeta + beta0)                                       
 
     double precision :: th
 
     integer            :: maxit = 1000
     double precision   :: rsum, res
     double precision   :: erraccu
     integer            :: i
     double precision   :: t
     double precision   :: omega
     double precision   :: o0
     double precision   :: a0
     double precision   :: rsum
     double precision   :: str_tau0, str_beta
     double precision   :: str_delta
     double precision   :: a,b
     double precision   :: res
     double precision   :: gampli(8)
     double precision   :: gwidth(8)
     double precision   :: gcenter(8)

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sqt_omeg'
       nparx =       12
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sqt_omeg = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " combined s(q,omega) and s(q,t) data from an s(q,t) model"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'intensit'  ! prefactor                                                                       
        parnam ( 2) = 'tau0    '  ! KWW time-constant, prefactor in front of q-dependence                           
        parnam ( 3) = 'beta    '  ! streched exp, prefactor in front of-q-dependence                                
        parnam ( 4) = 'epsilon '  ! accuracy parameter for FT-integrations (DO NOT FIT)                             
        parnam ( 5) = 'omega0  '  ! omega scale zero shift                                                          
        parnam ( 6) = 'u_sqr   '  ! <u**2> value for Debye-Waller-Factor                                            
        parnam ( 7) = 'j0      '  ! jump length (if applicable)                                                     
        parnam ( 8) = 'beta0   '  ! beta offset                                                                     
        parnam ( 9) = 'qexp_t0 '  ! q-exponent for tau0                                                             
        parnam (10) = 'qexp_bet'  ! beta-exponent                                                                   
        parnam (11) = 'bkgr    '  ! constant background                                                             
        parnam (12) = 'n_mg    '  ! fraction of hydrogen in side chains                                             
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "KWW time-constant, prefactor in front of q-dependence" !//cr//parspace//&
        th_param_desc( 3,idesc) = "streched exp, prefactor in front of-q-dependence" !//cr//parspace//&
        th_param_desc( 4,idesc) = "accuracy parameter for FT-integrations (DO NOT FIT)" !//cr//parspace//&
        th_param_desc( 5,idesc) = "omega scale zero shift" !//cr//parspace//&
        th_param_desc( 6,idesc) = "<u**2> value for Debye-Waller-Factor" !//cr//parspace//&
        th_param_desc( 7,idesc) = "jump length (if applicable)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "beta offset" !//cr//parspace//&
        th_param_desc( 9,idesc) = "q-exponent for tau0" !//cr//parspace//&
        th_param_desc(10,idesc) = "beta-exponent" !//cr//parspace//&
        th_param_desc(11,idesc) = "constant background" !//cr//parspace//&
        th_param_desc(12,idesc) = "fraction of hydrogen in side chains" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
        th_file_param(  2,idesc) = "nse      > explict switch between sqt and sqomega (if nse=1 ==> sqt)"
        th_file_param(  3,idesc) = "bk1level > resolution description bgr-levwl"
        th_file_param(  4,idesc) = "bk1slope > resolution description slope of bgr"
        th_file_param(  5,idesc) = "ga1inten > resolution description 1. gaussian amplitude"
        th_file_param(  6,idesc) = "ga1width > resolution description 1. gaussian width"
        th_file_param(  7,idesc) = "ga1cente > resolution description 1. gaussian center"
        th_file_param(  8,idesc) = "ga2inten > resolution description 2. gaussian amplitude"
        th_file_param(  9,idesc) = "ga2width > resolution description 2. gaussian width"
        th_file_param( 10,idesc) = "ga2cente > resolution description 2. gaussian center"
        th_file_param( 11,idesc) = "ga3inten > resolution description 3. gaussian amplitude"
        th_file_param( 12,idesc) = "ga3width > resolution description 3. gaussian width"
        th_file_param( 13,idesc) = "ga3cente > resolution description 3. gaussian center"
        th_file_param( 14,idesc) = "ga4inten > resolution description 4. gaussian amplitude"
        th_file_param( 15,idesc) = "ga4width > resolution description 4. gaussian width"
        th_file_param( 16,idesc) = "ga4cente > resolution description 4. gaussian center"
        th_file_param( 17,idesc) = "ga5inten > resolution description 5. gaussian amplitude"
        th_file_param( 18,idesc) = "ga5width > resolution description 5. gaussian width"
        th_file_param( 19,idesc) = "ga5cente > resolution description 5. gaussian center"
        th_file_param( 20,idesc) = "ga6inten > resolution description 6. gaussian amplitude"
        th_file_param( 21,idesc) = "ga6width > resolution description 6. gaussian width"
        th_file_param( 22,idesc) = "ga6cente > resolution description 6. gaussian center"
        th_file_param( 23,idesc) = "ga7inten > resolution description 7. gaussian amplitude"
        th_file_param( 24,idesc) = "ga7width > resolution description 7. gaussian width"
        th_file_param( 25,idesc) = "ga7cente > resolution description 7. gaussian center"
        th_file_param( 26,idesc) = "ga8inten > resolution description 8. gaussian amplitude"
        th_file_param( 27,idesc) = "ga8width > resolution description 8. gaussian width"
        th_file_param( 28,idesc) = "ga8cente > resolution description 8. gaussian center"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "tau0_q   > scaled: str_tau0 * (1+jlen*qz**qexpt0)/(qz**qexpt0)"
        th_out_param(  2,idesc) = "beta_q   > scaled: str_beta * (qz**qexpbeta + beta0)"
! 
        th_sqt_omeg = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      intensit =      pa( 1)
      tau0     =      pa( 2)
      beta     =      pa( 3)
      epsilon  =      pa( 4)
      omega0   =      pa( 5)
      u_sqr    =      pa( 6)
      j0       =      pa( 7)
      beta0    =      pa( 8)
      qexp_t0  =      pa( 9)
      qexp_bet =      pa(10)
      bkgr     =      pa(11)
      n_mg     =      pa(12)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: explict switch between sqt and sqomega (if nse=1 ==> sqt)
      xh =      0
      call parget('nse     ',xh,iadda,ier)
      nse      = xh
! >>> extract: resolution description bgr-levwl
      xh =      0
      call parget('bk1level',xh,iadda,ier)
      bk1level = xh
! >>> extract: resolution description slope of bgr
      xh =      0
      call parget('bk1slope',xh,iadda,ier)
      bk1slope = xh
! >>> extract: resolution description 1. gaussian amplitude
      xh =      0
      call parget('ga1inten',xh,iadda,ier)
      ga1inten = xh
! >>> extract: resolution description 1. gaussian width
      xh =      0
      call parget('ga1width',xh,iadda,ier)
      ga1width = xh
! >>> extract: resolution description 1. gaussian center
      xh =      0
      call parget('ga1cente',xh,iadda,ier)
      ga1cente = xh
! >>> extract: resolution description 2. gaussian amplitude
      xh =      0
      call parget('ga2inten',xh,iadda,ier)
      ga2inten = xh
! >>> extract: resolution description 2. gaussian width
      xh =      0
      call parget('ga2width',xh,iadda,ier)
      ga2width = xh
! >>> extract: resolution description 2. gaussian center
      xh =      0
      call parget('ga2cente',xh,iadda,ier)
      ga2cente = xh
! >>> extract: resolution description 3. gaussian amplitude
      xh =      0
      call parget('ga3inten',xh,iadda,ier)
      ga3inten = xh
! >>> extract: resolution description 3. gaussian width
      xh =      0
      call parget('ga3width',xh,iadda,ier)
      ga3width = xh
! >>> extract: resolution description 3. gaussian center
      xh =      0
      call parget('ga3cente',xh,iadda,ier)
      ga3cente = xh
! >>> extract: resolution description 4. gaussian amplitude
      xh =      0
      call parget('ga4inten',xh,iadda,ier)
      ga4inten = xh
! >>> extract: resolution description 4. gaussian width
      xh =      0
      call parget('ga4width',xh,iadda,ier)
      ga4width = xh
! >>> extract: resolution description 4. gaussian center
      xh =      0
      call parget('ga4cente',xh,iadda,ier)
      ga4cente = xh
! >>> extract: resolution description 5. gaussian amplitude
      xh =      0
      call parget('ga5inten',xh,iadda,ier)
      ga5inten = xh
! >>> extract: resolution description 5. gaussian width
      xh =      0
      call parget('ga5width',xh,iadda,ier)
      ga5width = xh
! >>> extract: resolution description 5. gaussian center
      xh =      0
      call parget('ga5cente',xh,iadda,ier)
      ga5cente = xh
! >>> extract: resolution description 6. gaussian amplitude
      xh =      0
      call parget('ga6inten',xh,iadda,ier)
      ga6inten = xh
! >>> extract: resolution description 6. gaussian width
      xh =      0
      call parget('ga6width',xh,iadda,ier)
      ga6width = xh
! >>> extract: resolution description 6. gaussian center
      xh =      0
      call parget('ga6cente',xh,iadda,ier)
      ga6cente = xh
! >>> extract: resolution description 7. gaussian amplitude
      xh =      0
      call parget('ga7inten',xh,iadda,ier)
      ga7inten = xh
! >>> extract: resolution description 7. gaussian width
      xh =      0
      call parget('ga7width',xh,iadda,ier)
      ga7width = xh
! >>> extract: resolution description 7. gaussian center
      xh =      0
      call parget('ga7cente',xh,iadda,ier)
      ga7cente = xh
! >>> extract: resolution description 8. gaussian amplitude
      xh =      0
      call parget('ga8inten',xh,iadda,ier)
      ga8inten = xh
! >>> extract: resolution description 8. gaussian width
      xh =      0
      call parget('ga8width',xh,iadda,ier)
      ga8width = xh
! >>> extract: resolution description 8. gaussian center
      xh =      0
      call parget('ga8cente',xh,iadda,ier)
      ga8cente = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     str_tau0 = str_tau0 * (1+j0*q**qexp_t0)/(q**qexp_t0)
     str_beta = str_beta * (q**qexp_bet + beta0)

     if( nse == 1 .or. q < 0.2d0) then
        gampli(1)   = ga1inten
        gwidth(1)   = ga1width
        gcenter(1)  = ga1cente

        gampli(2)   = ga2inten
        gwidth(2)   = ga2width
        gcenter(2)  = ga2cente

        gampli(3)   = ga3inten
        gwidth(3)   = ga3width
        gcenter(3)  = ga3cente

        gampli(4)   = ga4inten
        gwidth(4)   = ga4width
        gcenter(4)  = ga4cente

        gampli(5)   = ga5inten
        gwidth(5)   = ga5width
        gcenter(5)  = ga5cente

        gampli(6)   = ga6inten
        gwidth(6)   = ga6width
        gcenter(6)  = ga6cente

        gampli(7)   = ga7inten
        gwidth(7)   = ga7width
        gcenter(7)  = ga7cente

        gampli(8)   = ga8inten
        gwidth(8)   = ga8width
        gcenter(8)  = ga8cente

        o0         = x - omega0

        rsum = 0
drs:    do i=1,size(gampli)
         if(gampli(i) == 0d0) cycle
         omega     = o0 - gcenter(i)
         str_delta = abs(gwidth(i))
         a         = 0d0
         b         = 9.0d0/str_delta
         res       = adapint(fqt_kernel,a,b,epsilon,maxit,erraccu)*2
         rsum      = rsum + gampli(i)*res/(2*Pi)*sqrt(Pi)
        enddo drs

        dwf        = exp(-u_sqr*q*q/3.0d0)
        eisf       = (1.0/3.0)*(1.0+2.0*(sin(q*1.78)/(q*1.78)))

        th = intensit * dwf * (1-nmg+(nmg*eisf))*rsum + bkgr &
                      + bgr_level + bgr_slope*o0

    else
        t  = x
        th = intensit * fqt(t)
    endif
     th_sqt_omeg = th
 
! ---- writing computed parameters to the record >>>  
      call parset('tau0_q  ',sngl(tau0_q),iadda,ier)
      call parset('beta_q  ',sngl(beta_q),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

       function fqt_kernel(t)
!      -----------------------
       implicit none
       double precision, intent(in) :: t
       double precision             :: strex_kernel2

       fqt_kernel = fqt(t) * exp(-1d0/4d0*(str_delta*t)**2) * cos(t*Omega)*str_delta

       end function fqt_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! die eigentliche Zeitfunktion                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function fqt(t)
!      ---------------
       implicit none
       double precision, intent(in) :: t
       double precision             :: fqt

       fqt =  exp(-(t/str_tau0)**str_beta)

       end function fqt

 end function th_sqt_omeg
