                                                                       
      FUNCTION th4 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3)
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
 
      REAL(8) m1, m2, m3, m4, m5, m6, m7, m8, m9 
      REAL(8) c1, c2, c3, c4, c5, c6, c7, c8, c9 
      REAL(8) ampl, r, y 
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'approxi' 
         nparx = 15 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th4 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'ampl' 
         parnam (2) = 'm1' 
         parnam (3) = 'm2' 
         parnam (4) = 'm3' 
         parnam (5) = 'm4' 
         parnam (6) = 'm5' 
         parnam (7) = 'c1' 
         parnam (8) = 'c2' 
         parnam (9) = 'c3' 
         parnam (10) = 'c4' 
         parnam (11) = 'c5' 
         parnam (12) = 'm6' 
         parnam (13) = 'c6' 
         parnam (14) = 'm7' 
         parnam (15) = 'c7' 
!                                                                       
         th4 = 0 
         RETURN 
      ENDIF 
!                                                                       
      r = exp (x) 
      ampl = pa (1) 
      m1 = pa (2) 
      m2 = pa (3) 
      m3 = pa (4) 
      m4 = pa (5) 
      m5 = pa (6) 
      c1 = pa (7) 
      c2 = pa (8) 
      c3 = pa (9) 
      c4 = pa (10) 
      c5 = pa (11) 
      m6 = pa (12) 
      c6 = pa (13) 
      m7 = pa (14) 
      c7 = pa (15) 
                                                                        
! ---- calculate theory here -----                                      
      y = ampl * (r**m1 + c5 * r**m3 + c4 * r**m5) / (c1 + c2 * r**m2 + &
      c3 * r**m4 + c6 * r**m6 + c7 * r**m7)                             
      th4 = dlog (y) 
!      -----------------------------------------------                  
!                                                                       
      RETURN 
      END FUNCTION th4                              
!                  