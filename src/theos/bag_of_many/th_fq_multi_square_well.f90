      FUNCTION fqsheet (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================               
!  formfactor squared of a multiple square well
!                                                                       
!     |-dn-|....|-d2-|--d1--|---d0---|--d1--|-d2-|....|-dn-|
!       bn        b2    b1      b0      b1    b2        bn
!                                                                      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      DOUBLEPRECISION d (10), b (10), c (10), sum, q 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'fqsheet' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            fqsheet = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'd0' 
         parnam (2) = 'd1' 
         parnam (3) = 'd2' 
         parnam (4) = 'd3' 
         parnam (5) = 'd4' 
         parnam (6) = 'd5' 
         parnam (7) = 'd6' 
         parnam (8) = 'd7' 
         parnam (9) = 'd8' 
         parnam (10) = 'd9' 
         parnam (11) = 'b0' 
         parnam (12) = 'b1' 
         parnam (13) = 'b2' 
         parnam (14) = 'b3' 
         parnam (15) = 'b4' 
         parnam (16) = 'b5' 
         parnam (17) = 'b6' 
         parnam (18) = 'b7' 
         parnam (19) = 'b8' 
         parnam (20) = 'b9' 
!                                                                       
         fqsheet = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q = x 
      IF (abs (q) .lt.1d-8) q = 1d-8 
                                                                        
      DO i = 1, 10 
      d (i) = pa (i) 
      b (i) = pa (i + 10) 
      enddo 
                                                                        
      b (1) = b (1) / 2 
      c (1) = 0.0 
      DO i = 2, 10 
      c (i) = c (i - 1) + (d (i - 1) + d (i) ) / 2 
      enddo 
                                                                        
      sum = 0.0 
      DO i = 1, 10 
      sum = sum + b (i) * 2 * cos (c (i) * q) * sin (q * d (i) / 2) 
      enddo 
                                                                        
      fqsheet = (2 / q) * sum 
                                                                        
      fqsheet = fqsheet * fqsheet 
                                                                        
      RETURN 
      END FUNCTION fqsheet         
