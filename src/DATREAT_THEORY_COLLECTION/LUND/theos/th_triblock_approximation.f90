                                                      
                                                        
                                                                        
      FUNCTION th5 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> triblock1 approximation  <------                             
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			                                                                        
      DATA zpi / 6.283185 / 
                                                                        
      DOUBLEPRECISION t, q, rate, dispe, uq, a1, a2 
      REAL qp 
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'tribapr' 
         nparx = 3 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th5 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'uquadrat' 
         parnam (2) = 'rate' 
         parnam (3) = 'dispers' 
!                                                                       
         th5 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      t = x 
      uq = pa (1) 
      rate = pa (2) 
      dispe = pa (3) 
                                                                        
! extract parameter !                                                   
      CALL getpar ('q       ', qp,nopar ,params,napar,mbuf, ier)  
      IF (ier.ne.0) then 
         WRITE (6, * ) 'th5-brush: q-parameter not found! set q=0.1' 
         qp = 0.1 
      ENDIF 
      q = qp 
                                                                        
      IF (dispe.gt.10.0d0) then 
         WRITE (6, * ) 'warning: dispe=', dispe, ' readjusted' 
         dispe = 10 
      ENDIF 
      IF (dispe.lt. - 10.0d0) then 
         WRITE (6, * ) 'warning: dispe=', dispe, ' readjusted' 
         dispe = - 10 
      ENDIF 
                                                                        
      a1 = - (q * q * uq * uq) 
      a2 = - abs (rate) * (q**dispe) * t 
      IF (a1.lt. - 50d0) a1 = - 50.0d0 
      IF (a2.lt. - 50d0) a2 = - 50.0d0 
                                                                        
      th5 = exp (a1) + (1 - exp (a1) ) * exp (a2) 
!                                                                       
      RETURN 
      END FUNCTION th5                              
!               