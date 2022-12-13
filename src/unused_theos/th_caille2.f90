      FUNCTION caillexx (x, pa, thnam, parnam, npar, ini)
!     ===================================================               
!                                                                       
!     ===================================================               
!                                                                       
! -------> caille2<------                                               
! wege in q entlang schraeger schnitte durch die qz-qortho-ebene        
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (40), qq (3) 
      DATA zpi / 6.283185 / 
      REAL(8) qortho, qz, a, alpha, d, q0, dwf, peakf, scale 
      REAL(8) eps, eps2, xmuswit, strechf, rellim, splitf4, rres 
      REAL(8) caill4a, caille3, caille4, caille2, caillr3, caills3 
      REAL(8) caille 
                                                                        
      COMMON / thiadd / iadda 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'caille2' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            caillexx = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'scale' 
         parnam (2) = 'qz:qx' 
         parnam (3) = 'chi' 
         parnam (4) = 'q0' 
         parnam (5) = 'a' 
         parnam (6) = 'ln_alpha' 
         parnam (7) = 'd' 
         parnam (8) = '4a:3r:3s' 
         parnam (9) = 'nmax3' 
         parnam (10) = 'dwf' 
         parnam (11) = 'peakf' 
         parnam (12) = 'eps' 
         parnam (13) = 'eps2' 
         parnam (14) = 'xmuswit' 
         parnam (15) = 'strechf' 
         parnam (16) = 'rellim' 
         parnam (17) = 'splitf4' 
         parnam (18) = 'e1switch' 
         parnam (19) = 'polpow' 
         parnam (20) = 'rres' 
!                                                                       
         caillexx = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
      scale = pa (1) 
      qz_qx = pa (2) 
      chi = pa (3) 
      q0 = pa (4) 
      a = pa (5) 
      alpha = pa (6) 
      d = pa (7) 
      i4a3r3s = pa (8) + 0.1 
      nmax3 = pa (9) + 0.1 
      dwf = pa (10) 
      peakf = pa (11) 
      eps = pa (12) 
      eps2 = pa (13) 
      xmuswit = pa (14) 
      strechf = pa (15) 
      rellim = pa (16) 
      splitf4 = pa (17) 
      ie1switch = pa (18) + 0.1 
      ipolpow = pa (19) + 0.1 
      rres = pa (20) 
                                                                        
      IF (alpha.gt.69.0) alpha = 69 
      IF (alpha.lt. - 69.0) alpha = - 69 
                                                                        
      alpha = exp (alpha) 
                                                                        
                                                                        
                                                                        
      CALL setcai (eps, eps2, xmuswit, strechf, rellim, splitf4, rres,  &
      itch, ipolpow)                                                    
                                                                        
! --- look for qz or qortho parameters in the files ---                 
      ier = - 1 
      CALL parget ('qz:qx   ', qz_qx, iadda, ier) 
! --- and look for a, alpha and q0                                      
      ier = - 1 
      aa = a 
      CALL parget ('a       ', aa, iadda, ier) 
      a = aa 
      ier = - 1 
      aa = alpha 
      CALL parget ('alpha   ', aa, iadda, ier) 
      alpha = aa 
      ier = - 1 
      aa = q0 
      CALL parget ('q0      ', aa, iadda, ier) 
      q0 = aa 
                                                                        
      CALL qcalc (x, qz_qx, chi, qz_s, qortho_s) 
      qz = qz_s 
      qortho = qortho_s 
                                                                        
      GOTO (10, 20, 30), i4a3r3s 
   10 CONTINUE 
      caille = caill4a (qortho, qz, q0, d, a, alpha) 
      GOTO 100 
   20 CONTINUE 
      caille = caillr3 (qortho, qz, q0, d, a, alpha, nmax3) 
      GOTO 100 
   30 CONTINUE 
      caille = caills3 (qortho, qz, q0, d, a, alpha, peakf, dwf) 
  100 CONTINUE 
                                                                        
                                                                        
                                                                        
      caillexx = caille * scale
!                                                                       
      RETURN 
      END FUNCTION caillexx
