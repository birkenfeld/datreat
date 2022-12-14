      FUNCTION saxabs (q, R, db, n) 
!     =============================                                     
! ---> (KWS)-Streuintensitaet eines Schichtsystems auf runden Plaettchen
!      ---------------------------------------------------------------  
!      q          = Impulsuebertrag                                     
!      R          = Plaettchenradius ( R neg. ==> kleiner Integrationsb.
!      db(2,1..n) = Schichtdicken, Streulaengendichte                   
!      n          = Anzahl Schichten                                    
!      ---------------------------------------------------------------  
                                                                        
      IMPLICIT real (a - h, o - z) 
      DIMENSION db (2, n) 
      PARAMETER (maxl = 20) 
      COMMON / cuj0ui / qc, rc, dbc (2, maxl), nc 
      EXTERNAL fsfo 
      DATA pi / 3.141592654 /, qmin / 1d-30 / 
! ---                                                                   
! to be sure that Table in uj0uf is initialized                         
! otherwise there is a recursive call of qdag !                         
      xyz = 0.0 
      dummy = uj0uf (xyz) 
! --- Data-transfer to integrand fsfo via common:                       
      nc = n 
      rc = abs (r) 
      qc = abs (q) 
      IF (qc.lt.qmin) qc = qmin 
      DO i = 1, n 
      dbc (1, i) = db (1, i) 
      dbc (2, i) = db (2, i) 
      enddo 
! --- Integration ---                                                   
      errabs = 1e-5 
      errrel = 1e-4 
      irule = 2 
      a = 0.0 
      b = pi * 0.5 
      IF (r.lt.0.0) then 
!                  ----> nur relevante Zone integrieren !               
         as = (4 * pi) / (qc * rc) 
         IF (as.lt.1.0) b = asin (as) 
      ENDIF 
      CALL qdag (fsfo, a, b, errabs, errrel, irule, result, errest) 
      saxabs = result 
!     ---------------                                                   
                                                                        
      RETURN 
      END FUNCTION saxabs                           
                                                                        
      FUNCTION fsfo (u) 
!     ================                                                  
! --- Integrand function for angular averaging ---                      
      IMPLICIT real (a - h, o - z) 
      PARAMETER (maxl = 20) 
      COMMON / cuj0ui / qc, rc, dbc (2, maxl), nc 
! ---                                                                   
      v = sin (u) 
      qr = qc * rc * v 
      apara = uj0uf (qr) * (rc**4) 
      q = qc * cos (u) 
      aorth = orthin (q, dbc, nc) 
                                                                        
      fsfo = apara * aorth * v 
!      -------------------------                                        
      RETURN 
      END FUNCTION fsfo                             
                                                                        
                                                                        
      FUNCTION uj0uf (qr) 
!      ==================                                               
! ----------------------------------------------------------------------
! --- berechnet nach Tabelle :                                         C
!     (1/x**2) Integral(0..x)(u*J0(u)*du)                              C
!     ------------------------------------------------------------     C
! --->|das entspricht der "in-plane" Streuung eines Zylinders    |     C
!     |die Streuintensitaet ist proportional: R**2 * uj0uf( R*Q )|     C
!     ------------------------------------------------------------     C
!     x = qr                                                           C
! ----------------------------------------------------------------------
!                                                                       
      IMPLICIT real (a - h, o - z) 
      PARAMETER (min = - 5, max = 1000) 
      DIMENSION tab (min:max) 
      SAVE tab, dx 
      LOGICAL init / .false. / 
      SAVE init 
      EXTERNAL fu 
      DATA pi / 3.141592654 / 
! ---                                                                   
      IF (.not.init) then 
!             berechne Tabelle von :                                    
!             (1/x**2) Integral(0..x)(u*J0(u)*du)                       
!             -----------------------------------                       
         Write (6, * ) 'Start Uj0u-table generation no-elems=', max 
         dx = 0.3 
         xstart = 0.000001 
                                                                        
         errabs = 1e-4 
         errrel = 1e-3 
         irule = 4 
                                                                        
         a = 0 
         b = xstart 
         s = 0 
         DO i = 0, max 
         CALL qdag (fu, a, b, errabs, errrel, irule, result, errest) 
         s = s + result 
         tab (i) = 2 * pi * s / (b * b) 
         a = b 
         b = b + dx 
         IF (mod (i, 100) .eq.0) Write (6, * ) i, a, tab (i), errest 
         enddo 
!--                                                                     
         ntab = max 
         DO i = 1, - min 
         tab ( - i) = tab (i) 
         enddo 
         init = .true. 
         WRITE (6, * ) 'Uj0uf: Table has been initialized' 
      ENDIF 
! -----------------------------------------------------------------     
! Interpolation                                                         
! -----------------------------------------------------------------     
      n = qr / dx 
      IF (n.gt.ntab - 3) then 
         uj0uf = 0.d0 
         WRITE (6, * ) 'Uj0u: Table-Limits exceeded ! ' 
         RETURN 
      ENDIF 
                                                                        
      p = qr / dx - n 
      pp = p * p 
      uj0uf = (pp - 1.d0) * p * (p - 2.d0) * tab (n - 2) / 24.d0 -      &
      (p - 1.d0) * p * (pp - 4.d0) * tab (n - 1) / 6.d0 + (pp - 1.d0)   &
      * (pp - 4.d0) * tab (n) / 4.d0 - (p + 1.d0) * p * (pp - 4.d0)     &
      * tab (n + 1) / 6.d0 + (pp - 1.d0) * p * (p + 2.d0) * tab (n + 2) &
      / 24.d0                                                           
      uj0uf = uj0uf**2 
!      ---------------> Intensity !                                     
                                                                        
      RETURN 
      END FUNCTION uj0uf                            
!                                                                       
      FUNCTION fu (u) 
!      -------------- ( Integrand fuer Tabellengenerator UJ0UF )        
      IMPLICIT real (a - h, o - z) 
      fu = u * Bsj0 (u) 
      RETURN 
      END FUNCTION fu                               
!                                                                       
                                                                        
      FUNCTION orthin (q, db, n) 
!      ===========================                                      
! ---> 1-D Sequenz von Schichtdicken und Streukontrast --> Intensitaet  
!      ---------------------------------------------------------------  
!      q          = Impulsuebertrag                                     
!      db(2,1..n) = Schichtdicken, Streulaengendichte                   
!      n          = Anzahl Schichten                                    
!      ---------------------------------------------------------------  
                                                                        
      IMPLICIT real (a - h, o - z) 
      DIMENSION db (2, n) 
      COMPLEX csum, cexp, ci, conjg, cabs 
      DATA ci / (0.0, 1.0) /, pi / 3.141592654 /, qmin / 1d-30 / 
!      ------                                                           
      qq = abs (q) 
      IF (qq.lt.qmin) qq = qmin 
      csum = (0.0, 0.0) 
      x = 0.0 
      DO i = 1, n 
      d = db (1, i) 
      b = db (2, i) 
      csum = csum + b * (cexp (ci * qq * (x + d) ) - cexp (ci * qq * x) &
      ) / qq                                                            
      x = x + d 
      enddo 
                                                                        
      orthin = csum * conjg (csum) 
                                                                        
      RETURN 
      END FUNCTION orthin                           
