      FUNCTION th23 (x, pa, thnam, parnam, npar, ini) 
!     ===================================================               
!                                                                       
!     ===================================================               
!                                                                       
! -------> modified caille : caillm <---------                          
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (40), qq (3) 
      DATA zpi / 6.283185 / 
      REAL(8) qortho, qz, a, alpha, d, q0, dwf, peakf, scale 
      REAL(8) eps, eps2, xmuswit, strechf, rellim, splitf4, rres 
      REAL(8) caill4a, caille3, caille4, caille2, caillr3, caills3 
      REAL(8) caille 
      REAL(8) kappa_kt, b_si, temperatur, kb, d_si, Pi 
      REAL(8) temp, len_unit, eta 
                                                                        
      PARAMETER (Pi = 3.141592654d0) 
                                     ! Boltzmann k in Joule/K           
      PARAMETER (kb = 1.380662d-23) 
                                     ! m --> Angstroem                  
      PARAMETER (len_unit = 1.0d-10) 
                                                                        
      COMMON / thiadd / iadda 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'caillm' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th23 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'scale' 
         parnam (2) = 'qz:ortho' 
         parnam (3) = 'iz:ortho' 
         parnam (4) = 'q0' 
         parnam (5) = 'kappa_kt' 
         parnam (6) = 'b_si' 
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
         th23 = 0 
         RETURN 
      ENDIF 
                                                                        
      iot = iout () 
!                                                                       
! ---- calculate theory here -----                                      
      temperatur = 300.0d0 
      temp = 300.0 
                                                                        
      scale = pa (1) 
      qz_ortho = pa (2) 
      iz_ortho = pa (3) + 0.1 
      q0 = pa (4) 
      kappa_kt = pa (5) 
      b_si = pa (6) 
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
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
! --- look for qz or qortho parameters in the files ---                 
      ier = - 1 
      CALL parget ('qz      ', qz_ortho, iadda, ier) 
      IF (ier.eq.0) iz_ortho = 0 
!     --> qortho overrides qz                                           
      ier = - 1 
      CALL parget ('qortho  ', qz_ortho, iadda, ier) 
      IF (ier.eq.0) iz_ortho = 1 
! --- and look for a, alpha and q0                                      
                                                                        
      temp = temperatur 
      CALL parget ('temp    ', temp, iadda, ier) 
      temperatur = temp 
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
                                                                        
                                                                        
      b_si = abs (b_si) 
      kappa_kt = abs (kappa_kt) 
      d_si = 1d-10 * d 
                                                                        
      alpha = kappa_kt * kb * temperatur / (d_si * b_si) 
                                                                        
      a = kb * temperatur / (4 * Pi * b_si * sqrt (alpha) ) / (         &
      len_unit**2)                                                      
      alpha = alpha / (len_unit**2) 
                                                                        
                                                                        
      eta = 0.5d0 * kb * temperatur * Pi / (b_si * d_si * d_si * sqrt ( &
      kappa_kt * kb * temperatur / (d_si * b_si) ) )                    
                                                                        
      CALL parset ('eta     ', sngl (eta) , iadda) 
      CALL parset ('a       ', sngl (a) , iadda) 
      CALL parset ('ln_alpha ', sngl (dlog (alpha) ) , iadda) 
                                                                        
      IF (iot.gt.0) then 
      WRITE (6, '(a,e17.8)') 'b_si     = ', b_si 
         WRITE (6, '(a,e17.8)') 'kappa_kt = ', kappa_kt 
      WRITE (6, '(a,e17.8)') 'd_si     = ', d_si 
      WRITE (6, '(a,e17.8)') 'a        = ', a 
         WRITE (6, '(a,e17.8)') 'ln(alpha)= ', dlog (alpha) 
      WRITE (6, '(a,e17.8)') 'alpha    = ', alpha 
      WRITE (6, '(a,e17.8)') 'temp     = ', temperatur 
      WRITE (6, '(a,e17.8)') 'eta      = ', eta 
      ENDIF 
                                                                        
                                                                        
                                                                        
      CALL setcai (eps, eps2, xmuswit, strechf, rellim, splitf4, rres,  &
      itch, ipolpow)                                                    
                                                                        
                                                                        
                                                                        
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
         caille = cail2r3 (qortho, qz, q0, d, a, alpha, nmax3) 
!!         caille = caills3(qortho,qz,q0,d,a,alpha,peakf,dwf)           
  100    CONTINUE 
!         write(6,'(7e15.6)')qortho,qz,q0,d,a,alpha,caille              
                                                                        
         IF (iot.gt.0) write (6, * ) 'caille = ', caille, ' scale= ',   &
         scale                                                          
         th23 = caille * scale 
                                                                        
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
                                                                        
         th23 = caille * scale 
                                                                        
              !< rocking curve !                                        
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION th23                             
