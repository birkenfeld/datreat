      FUNCTION th27 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)  
!     ===================================================               
!                                                                       
! -------> cumulant <--------                                           
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			REAL mu, ampli, sum, ff 
      DIMENSION mu (20) 
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'cumulant' 
         nparx = 7 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th27 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitu' 
         parnam (2) = 'mu1' 
         parnam (3) = 'mu2' 
         parnam (4) = 'mu3' 
         parnam (5) = 'mu4' 
         parnam (6) = 'mu5' 
         parnam (7) = 'mu6' 
                                                                        
!                                                                       
         th27 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      ampli = pa (1) 
      DO i = 2, nparx 
      mu (i + 1 - 2) = pa (i) 
      enddo 
                                                                        
      ff = 1.0 
      sum = 0.0 
      DO i = 1, 6 
      ff = ff * i 
      sum = sum + mu (i) * (x**i) / ff 
      enddo 
      IF (sum.gt.20.0) sum = 20.0 
      IF (sum.lt. - 20.0) sum = - 20.0 
                                                                        
      th27 = ampli * exp ( - sum) 
!                                                                       
      RETURN 
      END FUNCTION th27        