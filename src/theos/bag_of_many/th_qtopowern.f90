                                                                       
      FUNCTION th10 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> q**n <--------                                               
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
                                     !!<-------!!                       
      REAL n 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'q**n' 
         nparx = 3 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th10 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'an' 
         parnam (2) = 'n' 
         parnam (3) = 'offset' 
!                                                                       
         th10 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      n = pa (2) 
      offset = pa (3) 
      th10 = pa (1) * (abs (x - offset) ) **n 
!                                                                       
      RETURN 
      END FUNCTION th10       