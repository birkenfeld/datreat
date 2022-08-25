      FUNCTION th_degennesA (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!  degennes 
!                                                                       
                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION Pi 
      Parameter (Pi = 3.141592654d0) 
                                                                        
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL th_degennesA, x, pa, qq, zpi, xh, vol_frac 
      INTEGER ini, npar, nparx 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
                                                                        
      DOUBLEPRECISION w, d, q, n, ne, degennesA, t, a, l
      double precision :: dp
         
                                                                        
                                                                        
                                                                        
      INTEGER  ier, iot 
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'degennea' 
         nparx = 6 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_degennesA = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
      parnam (1) = 'amplitud' 
      parnam (2)  = 'w       ' 
      parnam (3)  = 'dtube   ' 
      parnam (4)  = 'n       ' 
      parnam (5)  = 'l       '
      parnam (6)  = 'dpow    '  
!                                                                       
         th_degennesA = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
                          ! amplitude (should be 1)                     
      a = pa (1) 
                          ! Rouse rate                                  
      w = abs(pa (2)) 
                          ! Tube diameter                               
      d = abs(pa (3)) 
                          ! no of segments of polymer                   
      n = abs(pa (4)) 
                          ! segment length                              
      l = abs(pa (5)) 

      dp = abs(pa(6))      ! Non Gaussian power (dp = 1 = std)

      if(dp == 0d0) dp = 1d0
                                                                        
      CALL getpar ('q       ', xh,nopar ,params,napar,mbuf, ier)  
      IF (ier.eq.0) then 
         q = xh 
      ELSE 
         q = 1.0 
      ENDIF 
                                                                        
      t = x 
                                                                        
      th_degennesA = degennesA (t, q, d, W, n, l, ne, dp) 
                                                                        
      th_degennesA = th_degennesA * a 
                                                                        
      CALL setpar ('ne      ', sngl (ne) ,nopar ,params,napar,mbuf, ier)  
                                                                        
      RETURN 
      END FUNCTION th_degennesA                             
                                                                        
                                                                        
      DOUBLEPRECISION function degennesA (t, q, d, W, n, l, ne, dpow) 
!      ----------------------------------------------------             
!                                                                       
! siehe P.Schleger et. al. PRL 81, p.124 (1998)                         
!                                                                       
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION t, q, d, W , dpow
      DOUBLEPRECISION derfc, pi 
      DOUBLEPRECISION tau0, taud, t0, td 
      DOUBLEPRECISION n, ne, l, y, sum, eqd , eqd0
      double precision :: adapint, err, qx, dx
      INTEGER m, i, j 
                                                                        
      PARAMETER (pi = 3.141592654d0) 
                                                                        
!      l    =   d/sqrt(ne)                                              
      ne = (d / l) **2 
      tau0 = 36 / (W * ( (q * l) **4) ) 
      taud = 3 * (n**3) * (l**2) / ( (pi**2) * W * (d**2) ) 
                                                                        
      td = t / taud 
      t0 = t / tau0 
!      eqd = exp ( - ( (q * d / 6) **2) )
      qx = 0
      dx = d
      eqd0 = adapint(tube_kernel, 0d0, 5*d, 1d-8, 1000, err)
      qx = q
      eqd = adapint(tube_kernel, 0d0, 5*d, 1d-8, 1000, err)/eqd0      



!      write(6,*)'td,t0,eqd ',td, t0, eqd
                                                                        
!      m = 10 / sqrt (td+0.001D0) + 2 
      m = 30 / sqrt (td+0.001D0) + 2 
      sum = 0 
      DO i = 1, m, 2 
      sum = sum + exp ( - i * i * td) / (i * i) 
      enddo 
      sum = sum * (8 / (pi**2) ) * eqd
!     write(6,*)'m, sum ',m,sum
                                                                        
                                                                        
      degennesA = (1 - eqd) * exp (t0) * erfc (sqrt (t0) ) + sum 

!      write(6,*)'degennes ',degennes, erfc(sqrt(t0))
      
      RETURN 

      contains
        double precision function tube_kernel(r)
          implicit none
          double precision, intent(in) :: r
 
          tube_kernel = cos(qx*r) * exp(- ( ((3*r/dx)**2 ) ** dpow ) )         

        end function tube_kernel
 


      END FUNCTION degennesA                         
