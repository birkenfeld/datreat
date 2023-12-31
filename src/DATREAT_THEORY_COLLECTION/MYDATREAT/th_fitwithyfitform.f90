      function th3 (x, pa, thnam, parnam, npar,idum, ini) 
!     ===================================================               
!                                                                       
! -------> evaluate from function-buffer <-----                         
      COMMON / outlev / iout, ibild, ierrs, inka, iibuf, xxxx, yyyy, ptxf (20)                                                         
      CHARACTER(132) xformel, yformel, yfitform 
      COMMON / formul / xformel, yformel, yfitform 
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      REAL(8) val8y 
			integer idum
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
				thnam = 'eval  ' 
				nparx = 10 
				IF (npar.lt.nparx) then 
						WRITE (6, 1) thnam, nparx, npar 
		1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)                              
						th3 = 0 
						RETURN 
					ENDIF 
				npar = nparx 
	!        --------------> set the number of parameters                   
				parnam (1)  = 'p(1)  ' 
				parnam (2)  = 'p(2)  ' 
				parnam (3)  = 'p(3)  ' 
				parnam (4)  = 'p(4)  ' 
				parnam (5)  = 'p(5)  ' 
				parnam (6)  = 'p(6)  ' 
				parnam (7)  = 'p(7)  ' 
				parnam (8)  = 'p(8)  ' 
				parnam (9)  = 'p(9)  ' 
				parnam (10) = 'p(10) ' 
	!                                                                       
				th3 = 0 
				RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      DO i = 1, nparx 
      	ptxf (i) = pa (i) 
      enddo 
      xxxx = x 
      iout = iout - 5 
      CALL evaluate (yfitform, val8y, iery) 
      iout = iout + 5 
      th3 = val8y 
!                                                                       
      RETURN 
      END FUNCTION th3                              
