      FUNCTION subdiff (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
! diffqav 
! average over exp(-q**2 t) assuming Intensity goes like 1/q**2
!
!                                                    
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 /

      double precision :: uu, tau, beta


!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'subdiff' 
         nparx = 4
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            subdiff = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplit' 
         parnam (2) = 'uu' 
         parnam (3) = 'tau' 
         parnam (4) = 'beta' 
                                                                        
                                                                        
         subdiff = 0.0 
                                                                        
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q = 0.0 
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)  
      IF (q.eq.0) write (6, * ) 'ERROR: q not found' 

      uu   = abs(pa(2))
      tau  = abs(pa(3))
      beta = abs(pa(4))
     


      subdiff = pa(1) * exp(-uu*q*q*(x/tau)**beta) 
!                                                                       
      RETURN 
      END FUNCTION subdiff
