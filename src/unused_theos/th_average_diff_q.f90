      FUNCTION diffqav (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
! diffqav 
! average over exp(-q**2 t) assuming Intensity goes like 1/q**2
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
         thnam = 'diffqav' 
         nparx = 3 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th26 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplit' 
         parnam (2) = 'd' 
         parnam (3) = 'qwid/2' 
                                                                        
                                                                        
         th26 = 0.0 
                                                                        
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q = 0.0 
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)  
      IF (q.eq.0) write (6, * ) 'ERROR: q not found' 
      IF (q.le.abs (pa (3) ) ) q = abs (pa (3) ) + 0.001 
                                                                        
      t = x 
      D = abs (pa (2) ) 
      q1 = q - abs (pa (3) ) 
      q2 = q + abs (pa (3) ) 
                                                                        
                                                                        
      IF (q2.le.q1) q2 = q1 + 0.0001 
      IF (t.le.1e-6) t = 1e-6 
                                                                        
      IF (d * q * q * t.gt.50.0) d = 50.0 / (q * q * t) 
                                                                        
      t1 = q2**2 
      t7 = sqrt (D * t) 
      t11 = sqrt (0.3141593E1) 
      t17 = q1**2 
      t35 = (q1 * exp ( - D * t1 * t) * t7 + q1 * D * t * t11 * erf (t7 &
      * q2) * q2 - q2 * exp ( - D * t17 * t) * t7 - q2 * D * t * t11 *  &
      erf (t7 * q1) * q1) / t7 / ( - q2 + q1)                           
                                                                        
      diffqav = t35 * pa (1) 
!                                                                       
      RETURN 
      END FUNCTION diffqav    
