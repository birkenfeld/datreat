      FUNCTION th34 (x, pa, thnam, parnam, npar, ini) 
!     ===================================================               
!                                                                       
!     ===================================================               
!                                                                       
! -------> caille <------                                               
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
         thnam = 'caille' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th34 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'scale' 
         parnam (2) = 'qz:ortho' 
         parnam (3) = 'iz:ortho' 
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
         th34 = 0 
         RETURN 
      ENDIF 
                                                                        
      iot = iout () 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
      scale = pa (1) 
      qz_ortho = pa (2) 
      iz_ortho = pa (3) + 0.1 
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
      CALL parget ('qz      ', qz_ortho, iadda, ier) 
      IF (ier.eq.0) iz_ortho = 0 
!     --> qortho overrides qz                                           
      ier = - 1 
      CALL parget ('qortho  ', qz_ortho, iadda, ier) 
      IF (ier.eq.0) iz_ortho = 1 
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
      ier = - 1 
      CALL parget ('q_rock  ', q_rock, iadda, ier) 
      IF (ier.eq.0) then 
         i_rock = 1 
      ELSE 
         i_rock = 0 
      ENDIF 
                                                                        
                                                                        
      IF (i_rock.eq.0) then 
         IF (iz_ortho.eq.0) then 
            qz = qz_ortho 
            qortho = x 
      IF (iot.gt.0) write (6,  * ) 'qortho =', x, ' at qz     =', qz_ort&
     &ho                                                                
         ELSE 
            qz = x 
            qortho = qz_ortho 
      IF (iot.gt.0) write (6,  * ) 'qz     =', x, ' at q_ortho=', qz_ort&
     &ho                                                                
         ENDIF 
                                                                        
         GOTO (10, 20, 30), i4a3r3s 
   10    CONTINUE 
         caille = caill4a (qortho, qz, q0, d, a, alpha) 
         GOTO 100 
   20    CONTINUE 
         caille = caillr3 (qortho, qz, q0, d, a, alpha, nmax3) 
         GOTO 100 
   30    CONTINUE 
         caille = caills3 (qortho, qz, q0, d, a, alpha, peakf, dwf) 
  100    CONTINUE 
                                                                        
         IF (iot.gt.0) write (6, * ) 'caille = ', caille, ' scale= ',   &
         scale                                                          
         th34 = caille * scale 
                                                                        
      ELSE 
                                                                        
!  this is a rocking curve, x is psi in degrees                         
         sr = sin (x * 3.14159 / 180.0) 
         cr = cos (x * 3.14159 / 180.0) 
         qz = q_rock * cr 
         qortho = q_rock * sr 
         IF (iot.gt.0) write (6, * ) 'Rocking curve x=', x, qz, qortho 
                                                                        
         ier = - 1 
         xnsum = 1.0 
      CALL parget ('nsum    ', xnsum, iadda, ier) 
         nsum = xnsum + 0.1 
                                                                        
         ier = - 1 
         sig_qr = 1.0 
      CALL parget ('sig_qr  ', sig_qr, iadda, ier) 
                                                                        
                                                                        
!! here we introduce some special resolution convolution for the        
!! rocking curve                                                        
         caille = 0.0 
         fsum = 0.0 
                                                                        
         IF (iot.gt.0) write (6, * ) 'nsum=', nsum, ' sig_qr=', sig_qr 
                                                                        
         DO isum = - nsum / 2, nsum / 2 
                                                                        
         IF (nsum.gt.1) then 
            dqr = sig_qr * 2 / (nsum / 2) 
         ELSE 
            dqr = 1.0 
         ENDIF 
         q_r = q_rock + isum * dqr 
                                                                        
         IF (q_r.le.0.0) goto 1111 
                                                                        
         qz = q_r * cr 
         qortho = q_r * sr 
                                                                        
         qz = abs (qz) 
         qortho = abs (qortho) 
                                                                        
                                                                        
         f = exp ( - (isum * dqr / sig_qr) **2) 
                                                                        
         GOTO (110, 120, 130), i4a3r3s 
  110    CONTINUE 
         caille = caille+f * caill4a (qortho, qz, q0, d, a, alpha) 
         GOTO 1100 
  120    CONTINUE 
         caille = caille+f * caillr3 (qortho, qz, q0, d, a, alpha,      &
         nmax3)                                                         
         GOTO 1100 
  130    CONTINUE 
         caille = caille+f * caills3 (qortho, qz, q0, d, a, alpha,      &
         peakf, dwf)                                                    
 1100    CONTINUE 
                                                                        
         fsum = fsum + f 
                                                                        
         IF (iot.gt.1) then 
            WRITE (6, * ) 'isum=', isum, ' qr=', q_r 
            WRITE (6, * ) 'f=', f, ' caille=', caille 
         ENDIF 
                                                                        
 1111    CONTINUE 
                                                                        
              ! isum                                                    
         enddo 
                                                                        
         IF (iot.gt.0) write (6, * ) 'fsum = ', fsum 
         IF (fsum.gt.0.0) caille = caille / fsum 
         IF (iot.gt.0) write (6, * ) 'caille = ', caille, ' scale= ',   &
         scale                                                          
                                                                        
         th34 = caille * scale 
                                                                        
              !< rocking curve !                                        
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION th34                  