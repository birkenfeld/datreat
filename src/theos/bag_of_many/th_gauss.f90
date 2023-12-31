      FUNCTION thgauss (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!      gaussian
!      Y= intensit*exp((x-center)/width)                                                                
      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'gauss' 
         nparx = 3
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           thgauss = 0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'width' 
         parnam (3) = 'center' 
!                                                                       
         thgauss = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      
      arg = ( (x - pa (3) ) / (pa (2)) ) **2
      IF (arg.gt.50.0) arg = 50.0 
      thgauss = pa (1) * exp ( - arg) 
!                                                                       
      RETURN 
      END FUNCTION thgauss