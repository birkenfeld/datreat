      FUNCTION sqRDG (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!     RDG formfactor of a sphere 
!                                                                       
!                                                                       
      PARAMETER (pi = 3.141592654) 
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
         thnam = 'sqRDG' 
         nparx = 3 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            sqRDG = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'density' 
         parnam (2) = 'r' 
         parnam (3) = 'drho' 
!                                                                       
         sqRDG = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      r = pa (2) 
      drho = pa (3) 
      z = abs (x * r) 
      aj1 = (sin (z) - z * cos (z) ) / (z**3) 
      fj1qr = 4 * pi / 3 * r * r * r * drho * 3 * aj1 
      sqRDG = fj1qr**2 * pa (1) 
!                                                                       
      RETURN 
      END FUNCTION sqRDG          
