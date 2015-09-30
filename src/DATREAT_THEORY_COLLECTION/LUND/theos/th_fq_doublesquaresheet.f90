      FUNCTION th35 (x, pa, thnam, parnam, npar,  ini, nopar ,params,napar,mbuf)
!     ===================================================               
!                                                                       
! -------> formfactor squared of a double square-well                   
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (40), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'doubshee' 
         nparx = 6 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th35 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'd_outer' 
         parnam (2) = 'd_inner' 
         parnam (3) = 'b_outer' 
         parnam (4) = 'b_inner' 
         parnam (5) = 'qz:qx' 
         parnam (6) = 'chi' 
!                                                                       
                                                                        
         th35 = 0 
         RETURN 
      ENDIF 
!                                                                       
      qz_qx = pa (5) 
      chi = pa (6) 
      ier = - 1 
      CALL getpar ('qz:qx   ', qz_qx,nopar ,params,napar,mbuf, ier) 
      CALL qcalc (x, qz_qx, chi, qz, qortho) 
                                                                        
                                                                        
! ---- calculate theory here -----                                      
      dd = 0.5 * (pa (1) + pa (2) ) 
!       th35 = (2/qz) * (2*pa(3)*cos(dd*qz)*sin(pa(1)*qz*0.5) +         
!     *                   pa(4)*sin(pa(2)*qz*0.5) )                     
                                                                        
      hv = (2 / qz) * pa (4) * sin (0.5 * pa (2) * qz) + 2 * pa (3)     &
      * ( (sin (0.5 * (2 * pa (1) + pa (2) ) * qz) ) / qz - sin (0.5 *  &
      pa (2) * qz) / qz)                                                
                                                                        
      th35 = hv * hv 
                                                                        
!                                                                       
      RETURN 
      END FUNCTION th35                           