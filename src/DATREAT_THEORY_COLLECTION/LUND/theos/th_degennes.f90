      FUNCTION th38 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ========================================                          
!                                                                       
!                          -------> degennes <--------                  
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION Pi 
      Parameter (Pi = 3.141592654d0) 
                                                                        
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL th38, x, pa, qq, zpi, xh, vol_frac 
      INTEGER ini, npar, nparx 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
                                                                        
      DOUBLEPRECISION w, d, q, n, ne, degennes, t, a, l 
                                                                        
                                                                        
                                                                        
      INTEGER  ier, iot 
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'degennes ' 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th38 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitud' 
      parnam (2)  = 'w       ' 
      parnam (3)  = 'dtube   ' 
      parnam (4)  = 'n       ' 
      parnam (5)  = 'l       ' 
!                                                                       
         th38 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
                          ! amplitude (should be 1)                     
      a = pa (1) 
                          ! Rouse rate                                  
      w = pa (2) 
                          ! Tube diameter                               
      d = pa (3) 
                          ! no of segments of polymer                   
      n = pa (4) 
                          ! segment length                              
      l = pa (5) 
                                                                        
      CALL getpar ('q       ', xh,nopar ,params,napar,mbuf, ier)  
      IF (ier.eq.0) then 
         q = xh 
      ELSE 
         q = 1.0 
      ENDIF 
                                                                        
      t = x 
                                                                        
      th38 = degennes (t, q, d, W, n, l, ne) 
                                                                        
      th38 = th38 * a 
                                                                        
      CALL setpar ('ne      ', sngl (ne) ,nopar ,params,napar,mbuf, ier)  
                                                                        
      RETURN 
      END FUNCTION th38                             
                                                                        
                                                                        
      DOUBLEPRECISION function degennes (t, q, d, W, n, l, ne) 
!      ----------------------------------------------------             
!                                                                       
! siehe P.Schleger et. al. PRL 81, p.124 (1998)                         
!                                                                       
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION t, q, d, W 
      DOUBLEPRECISION derfc, pi 
      DOUBLEPRECISION tau0, taud, t0, td 
      DOUBLEPRECISION n, ne, l, y, sum, eqd 
      INTEGER m, i, j 
                                                                        
      PARAMETER (pi = 3.141592654d0) 
                                                                        
!      l    =   d/sqrt(ne)                                              
      ne = (d / l) **2 
      tau0 = 36 / (W * ( (q * l) **4) ) 
      taud = 3 * (n**3) * (l**2) / ( (pi**2) * W * (d**2) ) 
                                                                        
      td = t / taud 
      t0 = t / tau0 
      eqd = exp ( - ( (q * d / 6) **2) ) 
                                                                        
      m = 10 / sqrt (td+0.001D0) + 2 
      sum = 0 
      DO i = 1, m, 2 
      sum = sum + exp ( - i * i * td) / (i * i) 
      enddo 
      sum = sum * (8 / (pi**2) ) * eqd 
                                                                        
                                                                        
      degennes = (1 - eqd) * exp (t0) * derfc (sqrt (t0) ) + sum 
                                                                        
      RETURN 
      END FUNCTION degennes                         
