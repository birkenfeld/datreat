      FUNCTION strexpo (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!     also a nice axample for programming a theory                                                                 
! -------> strexpo  <--------    stretched exponential function                                       
!                                                                       
!      
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
		real amp, tau0, bet, qexp, q, tau, pi, arg
		real strexpo
		integer nparx                             ! number of theory parameters
    integer ier                               ! error variable 		                                                        
    DATA pi/3.141592653589793/ 
!                                                                       
! ----- initialisation of theory-----                                            
      IF (ini.eq.0) then 
         thnam = 'strexpo'      ! name of theory, max 8*char
         nparx = 4              ! number of theory parameters
         IF (npar.lt.nparx) then 
            WRITE (6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)") thnam, nparx, npar 
            strexpo = 0 
            RETURN 
         ENDIF 
!        --------------> set the number of parameters and add a short comment here
         npar = nparx              ! to return number of parameters
         parnam (1) = 'amplitu'    ! amplitude   1
         parnam (2) = 'tau'         
         parnam (3) = 'beta' 
         parnam (4) = 'qexp' 
!                                                                       
         strexpo = 0 
         RETURN 
      ENDIF ! of initialisation
			
      !     get the parameters from the data 
      !     searched by the name of the parameters and assigned to variables
      !     old was call to  parget from data
		call getpar('q       ',q,nopar ,params,napar,mbuf, ier)
		if (ier.gt.0) write(*,*) 'q not found in data parameters   ier=',ier      ! only test it
			                                              
		! get the parameters from the Theory    (better names)                                
      amp  = pa (1) 
      tau0 = abs (pa (2) ) 
      bet  = pa (3) 
      qexp = pa (4) 
                                                                        
! ----NOW calculate theory here -----  
      q = abs (q) 
      tau = tau0 * q**qexp 
                                                                        
      IF (bet.lt.0.0) bet = 0.0 
      IF (bet.gt.1.0) bet = 1.0 
      arg = abs (x / tau) **bet 
      IF (arg.lt. - 30.0) arg = - 30.0 
      IF (arg.gt.30.0) arg = 30.0 
!    FINALY the result of this funktion  is ============                                                                        
      strexpo = amp * exp ( - arg) 
      !         for demonstration of calculation of a new parameter accessible by ?? +<name> and writen with save together with all data        
		!call setpar ('test    ',q ,nopar ,params,napar,mbuf, ier) 
		!if (ier.gt.0) write(*,*) 'could not set value ',q,'   ier=',ier
			
			                                                               
      RETURN 
      END FUNCTION strexpo        