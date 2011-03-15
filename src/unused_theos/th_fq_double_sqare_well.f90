     FUNCTION doubsheq (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> formfactor squared of a double square-well                   
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'doubsheq' 
         nparx = 4 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            doubsheq = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'd_outer' 
         parnam (2) = 'd_inner' 
         parnam (3) = 'b_outer' 
         parnam (4) = 'b_inner' 
!                                                                       
         doubsheq = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      dd = 0.5 * (pa (1) + pa (2) ) 
      doubsheq = (2 / x) * (2 * pa (3) * cos (dd * x) * sin (pa (1) * x *   &
      0.5) + pa (4) * sin (pa (2) * x * 0.5) )                          
      doubsheq = doubsheq * doubsheq
                                                                        
      RETURN 
      END FUNCTION doubsheq
