      FUNCTION thdatr1 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     =================================================                 
!                                                                       
! simple gauss 
!                                                                       
!                                                                       
                                                                       
      PARAMETER (pi = 3.141592654) 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) q, rave, drho, z, phi, den, rp, bq, pq, sq, r3, alpha,    &
      sqs, scale                                                        
      REAL(8) vol, cma, gu, go, epsilon, erroraccup, erroraccub 
                                                        !! aix.sp extchk
      REAL(8) pschulz, pschj1, betaj1, adapint, peryev 
                                                        !! aix          
      REAL(8) fiq1sch, fiq2sch 
      COMMON / cschulz / z, rave, drho, q 
      EXTERNAL fiq1sch 
      EXTERNAL fiq2sch 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'j1ave' 
!        ----------------> give here 8 character name of theory         
!                          sherical formfactor with averaging           
!                          by schulz-distribution and                   
!                          s(q) according to percus & yevik             
! see. m.kotlarchyk and h.s.chen, j.chem.phys 79 (1983) 2461            
         nparx = 9 
!        ----------------> number of required parameters                
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            thdatr1 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'drho*a*a' 
         parnam (2) = 'rave/a' 
         parnam (3) = 'z' 
         parnam (4) = 'den*a**3' 
         parnam (5) = 'rp/r3' 
         parnam (6) = 'scale' 
         parnam (7) = 'v/cm**3' 
         parnam (8) = 'select' 
         parnam (9) = 'rmax/int' 
!                                                                       
         thdatr1 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
!                                                                       
      isel = pa (8) + 0.5 
      rmax = pa (9) 
! --- selektion der funktion: isel=0 : j1ave --> streuung               
!                             isel=1 : s(q)  (percus yevik)             
!                             isel=2 : s'(q) (percus yevik)             
!                             isel=3 : den*v*p(q)                       
!                             isel=4 : beta(q)                          
!                             isel=5 : pschulz(r)                       
!                             isel=6 : log(j1ave)                       
                                                                        
      cma = 1d8 
!               ---> converison a into cm                               
      vol = pa (7) 
      scale = pa (6) 
      q = x * scale 
      drho = pa (1) 
      rave = pa (2) / scale 
      z = pa (3) 
      den = pa (4) 
      alpha = pa (5) 
      r3 = ( (z + 2) * (z + 3) * rave**3 / (z + 1) **2) ** (1.d0 / 3.d0) 
      rp = alpha * r3 
                                                                        
      IF (rmax.le.0.0) then 
         pq = pschj1 (q, drho, rave, z) 
         bq = betaj1 (q, rave, z) 
      ELSE 
!        ---> numerical integration                                     
         maxit = 50 
         epsilon = 1d-6 
         gu = 0.d0 
         go = rmax 
         pq = adapint (fiq2sch, gu, go, epsilon, maxit, erroraccup) 
         bq = adapint (fiq1sch, gu, go, epsilon, maxit, erroraccub) 
         bq = bq * bq / pq 
      ENDIF 
                                                                        
      sq = peryev (q, rp, den, - 1d-7) 
      sqs = 1 + bq * (sq - 1) 
                                                                        
      thdatr1 = vol * cma * den * pq * sqs * scale**3 
!        -------------------------                                      
      IF (isel.eq.6) thdatr1 = alog (thdatr1) 
      IF (isel.eq.1) thdatr1 = sq 
      IF (isel.eq.2) thdatr1 = sqs 
      IF (isel.eq.3) thdatr1 = vol * cma * den * pq * scale**3 
      IF (isel.eq.4) thdatr1 = bq 
      IF (isel.eq.5) thdatr1 = pschulz (q, rave, z) 
!                                                                       
      RETURN 
      END FUNCTION thdatr1                          
!                                                                       
!                                                                       
!                                                                       
                            