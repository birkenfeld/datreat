      FUNCTION thlorenz (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================               
!     lorenzian 
!        amp= pa (1)
!        XX = (x - pa (3) )
!        g =  pa (2)                                                           
!        thlorenz = amp/pi *(g/2)**2/(XX**2+(g/2)**2)
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
		real        amp,XX,g 
		DATA pi / 3.14159265 /
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'lorenz'
         nparx = 3
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           thlorenz = 0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'width' 
         parnam (3) = 'center' 
!                                                                       
         thlorenz = 0
         RETURN 
      ENDIF 
!                                                                       
! define your parameters here
		amp= pa (1)
		XX = (x - pa (3) )
		g =  pa (2)
! ---- calculate theory here -----
  		thlorenz = amp/pi *(g/2)**2/(XX**2+(g/2)**2)
    
!                                                                       
      RETURN 
      END FUNCTION thlorenz