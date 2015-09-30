      FUNCTION th36 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ========================================                          
!                                                                       
! -------> srouse <--------                                             
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) R, W, Sq0, Sq, Sqt 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff, dr, scos, SRouse 
      REAL qget, tget 
      INTEGER n 
                                                                        
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'srouse' 
         nparx = 8 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th36 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'xi_frict' 
         parnam (3) = 'n_segmnt' 
                                  ! : re**2 = N * b**2                  
      parnam (4)  = 're      ' 
      parnam (5)  = 'temp    ' 
         parnam (6) = 'com_diff' 
         parnam (7) = 'q_width ' 
                                  ! : projected stretch                 
      parnam (8)  = 'scos    ' 
                                                                        
!                                                                       
         th36 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau = x 
      a0 = pa (1) 
      xi = abs (pa (2) ) 
      N = nint (pa (3) ) 
      R = abs (pa (4) ) 
      temp = pa (5) 
                              ! in cm**2/sec                            
      diff = abs (pa (6) ) 
      q_width = pa (7) 
      scos = pa (8) 
                                                                        
                                       ! in A**2/ns                     
      diff = diff * 1d-9 / 1d-16 
                                                                        
      IF (epsilon.eq.0.0d0) epsilon = 1.0d-3 
                                                                        
      qget = 0.01 
      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier)  
      qz = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier)  
         temp = tget 
      ENDIF 
                                                                        
                                                                        
      th36 = 0 
      sum = 0 
      sumnorm = 0 
      nqw = 15 
      dqw = 4 * q_width / nqw 
      IF (q_width.eq.0) then 
         nqw = 0 
         q_width = 1.0d0 
      ENDIF 
                                                                        
      DO i = - nqw, nqw 
      qzz = qz + i * dqw 
      IF (qzz.gt.0) then 
         fn = dexp ( - (i * dqw / q_width) **2) 
                                                                        
! --- include center of mass diffusion ---                              
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)                         
         dr = diff 
                                                                        
         CALL Srouse1 (qzz, tau, temp, dr, xi, N, R, W, scos, Sq, Sqt) 
         sum = sum + fn * Sqt 
         sumnorm = sumnorm + Sq * fn 
                                                                        
      ENDIF 
      enddo 
                                                                        
      IF (sumnorm.gt.0.0d0) then 
         th36 = sum / sumnorm 
      ELSE 
         th36 = 0 
      ENDIF 
                                                                        
      CALL setpar ('w       ', sngl (W) , nopar ,params,napar,mbuf, ier) 
                                        ! in cm**2/s                    
      dr = dr / (1d-9 / 1d-16) 
      CALL setpar ('diff    ', sngl (dr) ,nopar ,params,napar,mbuf, ier) 
!                                                                       
      RETURN 
      END FUNCTION th36                             
                                                                        
                                                                        
                                                                        
                                                                        
      SUBROUTINE Srouse1 (q, t, temp, Dr, xi, N, R, W, scos, Sq, Sqt) 
!      =============================================------==            
!                                                                       
! Rouse expression for a chain of finite length:                        
! Input parameters:                                                     
!    q     ----> momentum transfer in A**-1                             
!    t     ----> time in nano-sec                                       
!    temp  ----> temperature in K                                       
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- 
!    xi    ----> friction coefficient (Si: Ns/m=kg/s, hier: kg/ns)      
!    N     ----> number of chain segments                               
!    R     ----> end-to-end distance of the polymer molecule            
!    scso  ----> projected stretch                                      
! Output parameters:                                                    
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2               
!    Sq    <--- S(Q)                                                    
!    Sqt   <--- S(Q,t)                                                  
! ------------------------------------------------------------          
!                                                                       
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION kb, pi 
      PARAMETER (kb = 1.380662d-23) 
      PARAMETER (pi = 3.141592654d0) 
                                                                        
      DOUBLEPRECISION q, t, temp, Dr, xi, R, W, Sq, Sqt 
      INTEGER N, nn, mm, p 
                                                                        
      DOUBLEPRECISION l, tau_p, kbt, Sq0, arg1, arg2 
      DOUBLEPRECISION lcos, scos, sig, t0, srouse_stretch 
                                                                        
      IF (N.le.0) then 
         W = 999 
         Sq = 999 
         Sqt = 999 
         WRITE (6, * ) 'Error Number of chain segments is <= 0!', N 
         RETURN 
      ENDIF 
                                                                        
! ---- determine the segment length l ----                              
      l = sqrt (R**2 / N) 
                                                                        
! ---- and the Rousefactor ----                                         
                                ! in Joule = kg*m**2/s**2               
      kbt = temp * kb 
                                ! in         kg*A**2/ns**2              
      kbt = kbt * 100 
                                ! in 1/ns                               
      W = 3 * kbt / (xi * (l**2) ) 
                                                                        
! ---- set the diffusion constant if input is zero --- !                
      IF (Dr.eq.0.0d0) then 
         Dr = kbt / (N * xi) 
      ENDIF 
                                                                        
                                                                        
      lcos = scos * l 
      sig = l 
      t0 = 1d-6 
                                                                        
                                                                        
      Sq = srouse_stretch (q, t0, lcos, sig, W, N) 
      Sqt = srouse_stretch (q, t, lcos, sig, W, N) 
                                                                        
      RETURN 
      END SUBROUTINE Srouse1                        
                                                                        
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function srouse_stretch (q, t, lcos, sig, W, N) 
!      ----------------------------------------------------------       
! deGennes form for streched chain (finite sum instead of inf!)         
! deGennes, Physics 3,pp 37-45, (1967)                                  
!                                                                       
!      q    = momentum transfer                                         
!      t    = time                                                      
!      lcos = projected strech l=<a_n>                                  
!      sig  = sigma in degennes paper --> l0                            
!      W    = rate                                                      
!      N    = chain length (integer)                                    
                                                                        
                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION q, t, lcos, sig, W 
      INTEGER N 
                                                                        
      DOUBLEPRECISION sum, f, s, r, g_rs, Pi, arg 
      INTEGER i, j 
      PARAMETER (Pi = 3.141592654) 
                                                                        
                                                                        
      f = ( - q * q * sig * sig / 6.0d0) 
      r = Sqrt (4 * W * t / Pi) 
                                                                        
      sum = N * exp (f * r * g_rs (0.0d0) ) 
                                                                        
      DO i = 1, N 
      arg = (f * (i + r * g_rs (i * i / (4 * W * t) ) ) ) 
      IF (arg.lt. - 50.0d0) arg = - 50.0d0 
      sum = sum + 2 * (N - i) * exp (arg) * cos (i * q * lcos) 
      enddo 
                                                                        
      srouse_stretch = sum / (N * N) 
                                                                        
      RETURN 
      END FUNCTION srouse_stretch                   
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function g_rs (u) 
!      ---------------------------------                                
!                                                                       
! This is the integral g(u) from deGennes, Physics 3,pp 37-45, (1967)   
!                                                                       
! $\int_1^{\infty} d \tau \, \tau^{-2} \exp(-\tau^2 u) $                
!                                                                       
! $=e^{-u} + \sqrt{\pi u} [ \erf(\sqrt{u})-1 ] $                        
!                                                                       
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION u, derf, Pi, s15adf 
      PARAMETER (Pi = 3.141592654) 
      INTEGER ifail 
                                                                        
      IF (u.gt.40.0d0) then 
         g_rs = 0.0d0 
         RETURN 
      ENDIF 
                                                                        
!                                                                       
!      g_rs = exp(-u)+sqrt(Pi*u)*(derf(sqrt(u))-1.0d0)                  
!                                                                       
      ifail = 0 
                                                          ! NAG subrouti
      g_rs = exp ( - u) + sqrt (Pi * u) * ( - s15adf (sqrt (u), ifail) ) 
                                                                        
      RETURN 
      END FUNCTION g_rs                             
