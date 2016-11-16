     FUNCTION zimm (x, pa, thnam, parnam, npar,  ini, nopar ,params,napar,mbuf) 
!                                                                       
!  zimm 
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_zimm, a 
      REAL(8) SQ_zimmT 
      REAL(8) epsilon, diff 
      REAL qget, tget 
                                                                        
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'zimm' 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            zimm = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'eta_solv' 
         parnam (3) = 'epsilon ' 
         parnam (4) = 'temp    ' 
         parnam (5) = 'com_diff' 
                                                                        
!                                                                       
         zimm = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau = x 
      a = pa (1) 
      eta = abs (pa (2) ) 
      epsilon = abs (pa (3) ) 
      temp = pa (4) 
                              ! in cm**2/sec                            
      diff = abs (pa (5) ) 
                                                                        
                                       ! in A**2/ns                     
      diff = diff * 1d-9 / 1d-16 
                                                                        
      IF (epsilon.eq.0.0d0) epsilon = 1.0d-3 
                                                                        
      qget = 0.01 
      CALL getpar ('q       ', qget, nopar ,params,napar,mbuf, ier) 
      qz = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier)  
         temp = tget 
      ENDIF 
                                                                        
! --- include center of mass diffusion ---                              
      a = a * dexp ( - qz * qz * diff * tau) 
                                                                        
      IF (pa (3) .lt.0) then 
         WRITE (6, * ) 'Full zimm computation from scratch!' 
         zimm = a * SQ_zimm (tau, qz, temp, eta, epsilon) 
      ELSE 
         zimm = a * SQ_zimmT (tau, qz, temp, eta, epsilon) 
      ENDIF 
                                                                        
!                                                                       
      RETURN 
      END FUNCTION zimm            
