      FUNCTION th20 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> sofq <--------                                               
!                                                                       
!                                                                       
      PARAMETER (pi = 3.141592654) 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) qqr (1), sq, eta, scl, gamma, r 
      REAL(8) d, den 
                                                        !! aix.sp extchk
      REAL(8) pschulz, pschj1, betaj1, adapint, peryev 
                                                        !! aix          
      REAL(8) q 
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'sofq' 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th20 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'volfrac' 
         parnam (2) = 'scalelen' 
         parnam (3) = 'gamma' 
         parnam (4) = 'r' 
         parnam (5) = 'mode' 
!                                                                       
         th20 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      eta = pa (1) 
      scl = pa (2) 
      gamma = pa (3) 
      r = pa (4) 
      mode = pa (5) + 0.1 
                                                                        
      IF (mode.ge.1.and.mode.le.5) then 
         qqr (1) = x * r * 2.0 
         CALL sofq (qqr, sq, 1, eta, scl, gamma, r, mode, ierr) 
         IF (ierr.lt.0) write (6, * ) 'sofq: ierr=', ierr 
      ELSE 
         IF (mode.eq.6) then 
            q = x 
            d = 2 * r 
            den = 6.0d0 * eta / (pi * d * d * d) 
            sq = peryev (q, r, den, - 1d-7) 
         ELSE 
            WRITE (6, * ) 'sofq: mode=', mode, ' is out of range' 
            sq = 1.0 
         ENDIF 
      ENDIF 
      th20 = sq 
!                                                                       
      RETURN 
      END FUNCTION th20          