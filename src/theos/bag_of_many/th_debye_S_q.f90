      FUNCTION debye2 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!  Debye function 
!                                                                       

                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION Pi, Navog 
      Parameter (Pi = 3.141592654d0) 
      Parameter (Navog = 6.022045d23) 
                                                                        
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL debye2, x, pa, qq, zpi, xh, vol_frac 
      INTEGER ini, npar, nparx 
      DIMENSION pa (20), qq (3) 
		integer :: mbuf
		integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 /
                                                                        
      DOUBLEPRECISION v, nphi 
                                                                        
      DOUBLEPRECISION ampli, q, rg, qrg, gamma, sq, sqi 
      COMMON / cdebik / gamma, qrg 
                                                                        
      DOUBLEPRECISION mw, bpoly, drho, conc, rhopoly, fac 
      REAL bsolv 
                                                                        
                                                                        
      DOUBLEPRECISION dbik, adapint, erra 
      EXTERNAL dbik 
                                                                        
                                                                        
      INTEGER  ier,  iot 
!                                                                       
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
      thnam = 'debye2  ' 
         nparx = 8 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            debye2 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
      parnam (1)  = 'intensit' 
      parnam (2)  = 'rg      ' 
      parnam (3)  = 'gamma   ' 
      parnam (4)  = 'molwght ' 
      parnam (5)  = 'density ' 
      parnam (6)  = 'bpolym  ' 
      parnam (7)  = 'v       ' 
         parnam (8) = 'volfrac ' 
!                                                                       
         debye2 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      ampli = abs (pa (1) ) 
      rg = abs (pa (2) ) 
      gamma = abs (pa (3) ) 
      mw = abs (pa (4) ) 
      rhopoly = abs (pa (5) ) 
      bpoly = pa (6) 
      v = pa (7) 
      conc = abs (pa (8) ) 
                                                                        
                                                                        
                                                                        
      CALL getpar ('bsolv   ', bsolv, nopar ,params,napar,mbuf, ier) 
      drho = bpoly - bsolv 
                                                                        
      fac = conc * mw / rhopoly / Navog 
                                                                        
      q = x 
      qrg = q * q * rg * rg 
                                                                        
      IF (gamma.eq.0.0d0) then 
         sq = (2 / (qrg**2) ) * (qrg - 1 + exp ( - qrg) ) 
      ELSE 
         sq = adapint (dbik, 0.0d0, 1.0d0, 1d-8, 1000, erra) 
      ENDIF 
                                                                        
! Zimm formula for interacting chains                                   
      sqi = 1.0d0 / (fac * sq) + v * conc * conc 
      IF (abs (sqi) .lt.1d-20) sqi = 1d-20 
                                                                        
      debye2 = ampli * drho**2 / sqi 
                                                                        
      RETURN 
      END FUNCTION debye2                             
                                                                        
                                                                        
      DOUBLEPRECISION function dbik (u) 
!       ---------------------------------                               
!                                                                        
      DOUBLEPRECISION u 
      DOUBLEPRECISION qrg, gamma 
      COMMON / cdebik / gamma, qrg 
                                                                        
      dbik = 2 * (1 - u) * exp ( - (u**gamma) * qrg) 
                                                                        
      RETURN 
      END FUNCTION dbik                             
