 FUNCTION th_codif(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  General diffusion expression assuming a nonlinear rr rr = ((a*(t/tref)**alpha + b*(t/tref)**beta)/r0**2)**gamma th   = int(w(D) * exp(-q^(2*delta)/6 *rr )
! 
      use theory_description 
      implicit none 
      real    :: th_codif
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
     double precision :: tref       ! reference time                                                                  
     double precision :: r0         ! reference length                                                                
     double precision :: a          ! first displacement parameter                                                    
     double precision :: alpha      ! exponent for a part                                                             
     double precision :: b          ! second displacement parameter                                                   
     double precision :: beta       ! exponent for b-part                                                             
     double precision :: gamma      ! exponent for all rr                                                             
     double precision :: delta      ! extra exponent for q   Default 1                                                
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
 
     double precision   :: t, rr
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'codif'
       nparx =        9
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_codif = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " General diffusion expression assuming a nonlinear rr " // &
                                "rr = ((a*(t/tref)**alpha + b*(t/tref)**beta)/r0**2)**gamma "//&
                                "th   = ampli* exp(-q^(2*delta)/6 *rr )"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'tref    '  ! reference time                                                                  
        parnam ( 3) = 'r0      '  ! reference length                                                                
        parnam ( 4) = 'a       '  ! first displacement parameter                                                    
        parnam ( 5) = 'alpha   '  ! exponent for a part                                                             
        parnam ( 6) = 'b       '  ! second displacement parameter                                                   
        parnam ( 7) = 'beta    '  ! exponent for b-part                                                             
        parnam ( 8) = 'gamma   '  ! exponent for all rr                                                             
        parnam ( 9) = 'delta   '  ! extra exponent for q   Default 1                                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "reference time" !//cr//parspace//&
        th_param_desc( 3,idesc) = "reference length" !//cr//parspace//&
        th_param_desc( 4,idesc) = "first displacement parameter" !//cr//parspace//&
        th_param_desc( 5,idesc) = "exponent for a part" !//cr//parspace//&
        th_param_desc( 6,idesc) = "second displacement parameter" !//cr//parspace//&
        th_param_desc( 7,idesc) = "exponent for b-part" !//cr//parspace//&
        th_param_desc( 8,idesc) = "exponent for all rr" !//cr//parspace//&
        th_param_desc( 9,idesc) = "extra exponent for q   Default 1" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_codif = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      tref     =      pa( 2)
      r0       =      pa( 3)
      a        =      pa( 4)
      alpha    =      pa( 5)
      b        =      pa( 6)
      beta     =      pa( 7)
      gamma    =      pa( 8)
      delta    =      pa( 9)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =  0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     rr = ((a*(t/tref)**alpha + b*(t/tref)**beta)/r0**2)**gamma
     th_codif = ampli* exp(-q**(2*delta)/6d0 *rr )
 
! ---- writing computed parameters to the record >>>  
 
 end function th_codif
