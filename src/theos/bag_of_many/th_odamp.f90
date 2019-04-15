      FUNCTION th_odamp (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!
!     
!     Combination of streched exponential with a oscillationg function
!     a building blocks for any combination of exp and oscillations        
!     y=amplitu*exp([x/(tau0)]^beta) * cos(freq*x+phase)
!     units: nothing special  
!     Reference : unknown author in unknown Journal

		implicit none          ! to make it easier to understand                                                          
		real x                                     ! xvalue
		real    pa (20)                            ! parameter values of theory
		CHARACTER(8) thnam                         ! theory    name in datreat
		CHARACTER(8) parnam (20)                   ! parameter namen of theory
		integer      npar                          ! maximum number of theory parameters in main program
		integer      ini                           ! flag initialisation or calculation
 		integer, intent(inout) :: nopar            ! number of Parameter in data
		character*80, intent(inout) :: napar(mbuf) ! name of  parameters(n)
		real, intent(inout) :: params(mbuf)        ! value of parameters(n)
		integer :: mbuf                            ! maximum length of data arrays in calling main program   
  
		! define local parametrs here
		real amp, tau0, beta, pi, xoffs, coffs, freq, phase, omega 
		real th_odamp
		integer nparx                             ! number of theory parameters
        integer ier                               ! error variable 		                                                        
        DATA pi/3.141592653589793/
!                                                                       
! ----- initialisation of theory-----                                            
      IF (ini.eq.0) then 
         thnam = 'odamp  '      ! name of theory, max 8*char
         nparx = 7            ! number of theory parameters
         IF (npar.lt.nparx) then 
            WRITE (6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)") thnam, nparx, npar 
            th_odamp = 0 
            RETURN 
         ENDIF 
!        --------------> set the number of parameters and add a short comment here
         npar = nparx              ! to return number of parameters
         parnam (1) = 'amplitu'    ! amplitude   should be 1
         parnam (2) = 'tau'        ! relaxation time 
         parnam (3) = 'beta'       ! stretched exponent
         parnam (4) = 'freq'       ! exponent of q
         parnam (5) = 'phase'      ! phase of cosine
         parnam (6) = 'exoffset'   ! exponent x offset
         parnam (7) = 'cosyoffs'   ! cos-term y-offset
!                                                                       
         th_odamp = 0 
         RETURN 
      ENDIF ! of initialisation
			
      !     get the parameters from the data 
      !     searched by the name of the parameters and assigned to variables
      !     old was call to  parget from data

      amp   = pa (1) 
      tau0  = pa (2)  
      beta  = pa (3) 
      freq  = pa (4)
      phase = pa (5) 
      xoffs = pa (6)
      coffs = pa (7)
 
                                                                        
  ! ----NOW calculate theory here -----
      omega = freq * 2*pi

      th_odamp = amp * exp ( -((x-xoffs)/tau0)**beta ) * (coffs + cos( omega * x - phase ))/(1d0+coffs) 
			                                                               
      RETURN 
      END FUNCTION th_odamp        
