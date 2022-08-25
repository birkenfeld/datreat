      FUNCTION thgauss2 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!      gaussian
!      Y= intensit*exp((x-center)/width)                                                                

       double precision :: xx, y, d, width, gw

       double precision, parameter :: Pi=3.141592654d0
      
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
         thnam = 'gauss2' 
         nparx = 4
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           thgauss2 = 0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'width' 
         parnam (3) = 'center'
         parnam (4) = 'xbinwid'    ! Breite der x-bins  
!                                                                       
         thgauss2 = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----
      xx    = x-pa(3)                                      
      width = pa(2)
      d     = pa(4)

      arg = ( (x - pa (3) ) / (pa (2)) ) **2
      IF (arg.gt.50.0)then
         arg = 50.0 
         thgauss2 = pa (1) * exp ( - arg)
      else          
         gw = -(0.5d0*sqrt(Pi)*width/d)*(erf(0.5d0*(2*xx-d)/width)-erf(0.5d0*(2*xx+d)/width))
         thgauss2 = pa(1) * gw
      endif
                                                             
      RETURN 
      END FUNCTION thgauss2
