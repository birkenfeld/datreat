                                                                 
      FUNCTION gauspolynom (x, pa, thnam, parnam, npar,  ini, nopar ,params,napar,mbuf)
!     ===================================================               
!                                                                       
! -------> gaupol <--------                                             
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
         thnam = 'gaupol' 
         nparx = 10 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            gauspolynom = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'center' 
         parnam (3) = 'a1' 
         parnam (4) = 'a2' 
         parnam (5) = 'a3' 
         parnam (6) = 'a4' 
         parnam (7) = 'a5' 
         parnam (8) = 'a6' 
         parnam (9) = 'a7' 
         parnam (10) = 'a8' 
!                                                                       
         gauspolynom = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      y = x - pa (2) 
      arg = pa (10) * y 
      arg = (arg + pa (9) ) * y 
      arg = (arg + pa (8) ) * y 
      arg = (arg + pa (7) ) * y 
      arg = (arg + pa (6) ) * y 
      arg = (arg + pa (5) ) * y 
      arg = (arg + pa (4) ) * y 
      arg = (arg + pa (3) ) * y 
      IF (arg.gt.50.0) arg = 50.0 
      gauspolynom = pa (1) * exp ( - arg)
!                                                                       
      RETURN 
      END FUNCTION gauspolynom