      FUNCTION th40 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================               
!                                                                       
! -------> expq <--------                                               
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
   			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			                                                                     
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'expq' 
         nparx = 4 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th40 = 0.0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitu' 
         parnam (2) = 'tau0' 
         parnam (3) = 'q0' 
         parnam (4) = 'q_exp' 
!                                                                       
         th40 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      CALL getpar ('q       ', qp,nopar ,params,napar,mbuf, ier) 
	
      t = x 
      amp = pa (1) 
      tau0 = pa (2) 
      q0 = pa (3) 
      q_exp = pa (4) 
                                                                        
                                                                        
      IF (q_exp.lt. - 30.0) q_exp = - 30.0 
      IF (q_exp.gt.30.0) q_exp = 30.0 
      tau = ( (abs (qp / q0) ) **q_exp) * abs (tau0) 
                                                                        
      arg = - t / tau 
                                                                        
      IF (arg.lt. - 50.0) arg = - 50.0 
                                                                        
                                                                        
      th40 = amp * exp (arg) 
!                                                                       
      RETURN 
      END FUNCTION th40                             
