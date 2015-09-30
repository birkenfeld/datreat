      FUNCTION th32 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rouse <--------                                              
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) SQ_rouseT 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff 
      REAL qget, tget 
                                                                        
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'rouse' 
         nparx = 7 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th32 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'xi_frict' 
         parnam (3) = 'b_segmnt' 
         parnam (4) = 'epsilon ' 
      parnam (5)  = 'temp    ' 
         parnam (6) = 'com_diff' 
         parnam (7) = 'q_width' 
                                                                        
!                                                                       
         th32 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau = x 
      a0 = pa (1) 
      xi = abs (pa (2) ) 
      b = abs (pa (3) ) 
      epsilon = abs (pa (4) ) 
      temp = pa (5) 
                              ! in cm**2/sec                            
      diff = abs (pa (6) ) 
      q_width = pa (7) 
                                                                        
                                       ! in A**2/ns                     
      diff = diff * 1d-9 / 1d-16 
                                                                        
      IF (epsilon.eq.0.0d0) epsilon = 1.0d-3 
                                                                        
      qget = 0.01 
      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier)  
      qz = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier) 
         temp = tget 
      ENDIF 
                                                                        
                                                                        
      th32 = 0 
      sum = 0 
      sumnorm = 0 
      nqw = 15 
      dqw = 4 * q_width / nqw 
      IF (q_width.eq.0) then 
         nqw = 0 
         q_width = 1.0d0 
      ENDIF 
                                                                        
      DO i = - nqw, nqw 
      qzz = qz + i * dqw 
      IF (qzz.gt.0) then 
         fn = dexp ( - (i * dqw / q_width) **2) 
         sumnorm = sumnorm + fn 
                                                                        
! --- include center of mass diffusion ---                              
         a = fn * a0 * dexp ( - qzz * qzz * diff * tau) 
                                                                        
         IF (pa (4) .lt.0) then 
            WRITE (6, * ) 'Full rouse computation from scratch!' 
            sum = sum + a * SQ_rouse (tau, qzz, temp, xi, b, epsilon) 
         ELSE 
            sum = sum + a * SQ_rouseT (tau, qzz, temp, xi, b, epsilon) 
         ENDIF 
      ENDIF 
      enddo 
                                                                        
      th32 = sum / sumnorm 
!                                                                       
      RETURN 
      END FUNCTION th32   