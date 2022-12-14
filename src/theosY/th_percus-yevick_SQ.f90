                                                                    
      FUNCTION th7 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!  test fuer percus-yevick s(q) 
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
                                                        !! aix.sp extchk
      REAL(8) pschulz, pschj1, betaj1, adapint
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
		    DATA zpi / 6.283185 /

      double precision :: q, rr, den, eps, peryev

                          
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'peryev' 
         nparx = 4 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th7 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intens' 
         parnam (2) = 'radius' 
         parnam (3) = 'density' 
         parnam (4) = 'epsilon' 
!                                                                       
         th7 = 0 
         RETURN 
      ENDIF 

      q   = x
      rr  = abs(pa(2))
      den = abs(pa(3))
      eps = abs(pa(4))
    
      if( eps == 0 ) eps = 1d-11
!                                                                       
! ---- calculate theory here -----                                      
!                                                                       
      th7 = pa (1) * peryev (q, rr, den, eps)
                                                     
                                                                        
      RETURN 


      END FUNCTION th7      
