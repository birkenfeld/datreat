      FUNCTION thgauss (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!      gaussian
!      Y= intensit*exp((x-center)/width) 
      use theory_description                                                               
      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n

                integer :: idesc
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

        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "simple Gaussian function:                               "//cr//&
                                 " th = intensity * exp[-{(x-center)/withd}**2]           "
            


         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' ; th_param_desc(1,idesc) = "prefactor "
         parnam (2) = 'width'    ; th_param_desc(2,idesc) = "width = sqrt(2)*sigma "
         parnam (3) = 'center'   ; th_param_desc(3,idesc) = "center of Gaussian "
!
         th_citation(idesc) = " ... den alten Carl Friedrich Gauss"
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
