      FUNCTION debyerpa (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!  RPA with Debye function 
!                                                                       

                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION Pi, Navog 
      Parameter (Pi = 3.141592654d0) 
      Parameter (Navog = 6.022045d23) 
                                                                        
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL debyerpa, x, pa, qq, zpi, xh, vol_frac 
      INTEGER ini, npar, nparx 
      DIMENSION pa (20), qq (3) 
		integer :: mbuf
		integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 /
                                                                        
      DOUBLEPRECISION nphi, chi 
                                                                        
      DOUBLEPRECISION ampli, q, rg, qrg, gamma, sq, sqi 
      COMMON / cdebik / gamma, qrg 
                                                                        
      DOUBLEPRECISION mw, mn, bpoly, drho, conc, rhopoly, fac
      REAL bsolv 
                                                                        
                                                                        
      DOUBLEPRECISION dbik1, adapint, erra, z 
      EXTERNAL dbik1 
                                                                        
                                                                        
      INTEGER  ier,  iot 
!                                                                       
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'debyerpa  ' 
         nparx = 9 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            debyerpa = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
      parnam (1)  = 'intensit' 
      parnam (2)  = 'rg      ' 
      parnam (3)  = 'gamma   ' 
      parnam (4)  = 'molwght '
      parnam (5)  = 'monwght '
      parnam (6)  = 'density ' 
      parnam (7)  = 'bpolym  ' 
      parnam (8)  = 'chi     ' 
      parnam (9)  = 'volfrac '     
!                                                                       
      debyerpa = 0 
      RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      ampli = abs (pa (1) ) 
      rg = abs (pa (2) ) 
      gamma = abs (pa (3) ) 
      mw = abs (pa (4) ) 
      mn = abs( pa(5) )
      rhopoly = abs (pa (6) ) 
      bpoly = pa (7) 
      chi = pa (8)  
      conc = abs (pa (9) ) 
                                                                       
        
      z = mw / mn    ! Nr. monomers                                                                 
                                                                        
      CALL getpar ('bsolv   ', bsolv, nopar ,params,napar,mbuf, ier) 
      drho = bpoly - bsolv 
                                                                        
      fac = mn / rhopoly / Navog 
                                                                        
      q = x 
      qrg = q * q * rg * rg 
                                                                        
      IF (gamma.eq.0.0d0) then 
         sq = (2 / (qrg**2) ) * (qrg - 1 + exp ( - qrg) ) 
      ELSE 
         sq = adapint (dbik1, 0.0d0, 1.0d0, 1d-8, 1000, erra) 
      ENDIF 
                                                                        
! Zimm formula for interacting chains                                   
!      sqi = 1.0d0 / (fac1 * sq) + 1.0d0 / (fac2 * sq) 
      sqi = 1.0d0 / (z* conc * sq) + 1.0d0 / (z* (1-conc) * sq) -2d0 * chi
      IF (abs (sqi) .lt.1d-20) sqi = 1d-20 
                                                                        
      debyerpa = ampli * drho**2 * fac / sqi 
                                                                        
      RETURN 
      END FUNCTION debyerpa                             
                                                                        
                                                                        
      DOUBLEPRECISION function dbik1 (u) 
!       ---------------------------------                               
!                                                                        
      DOUBLEPRECISION u 
      DOUBLEPRECISION qrg, gamma 
      COMMON / cdebik / gamma, qrg 
                                                                        
      dbik1 = 2 * (1 - u) * exp ( - (u**gamma) * qrg) 
                                                                        
      RETURN 
      END FUNCTION dbik1                             
