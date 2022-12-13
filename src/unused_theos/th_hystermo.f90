      FUNCTION hystermo (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> hystereseis-modelling  <------                                            
!                                                                       
        !
      implicit none  
      real  :: hystermo  
      CHARACTER(8) thnam, parnam (20)
      real :: x, pa, qq, zpi
      DIMENSION pa (20), qq (3) 
      integer :: mbuf, nparx, ini, npar
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) ani, ri, ra, centerx, centery, centerz, xlen, b (3),      &
      r (3), xa (3), xe (3)                                             
      REAL(8) epsilon, rsp, rssep, shift, cur
      integer :: isym , isel, maxita  
      COMMON / tceps / epsilon, maxita 
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'hystermo' 
         nparx = 10 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            hystermo = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'centerx' 
         parnam (2) = 'length' 
         parnam (3) = 'shift' 
         parnam (4) = 'rssep' 
         parnam (5) = 'ri' 
         parnam (6) = 'ra' 
         parnam (7) = 'cur' 
         parnam (8) = 'symm' 
         parnam (9) = 'maxit' 
         parnam (10) = 'epsilon' 
!                                                                       
         hystermo = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      centerx = pa (1) 
      xlen    = pa (2)
      shift   = pa (3)
      rssep   = pa (4)
      ri      = pa (5) 
      ra      = pa (6) 
      cur     = pa (7)
      isym    = nint(pa(8))
      isel   =  1
      maxita = nint(pa (9)) 
      epsilon = pa (10) 
      r (1) = x 
      r (2) = 0.0 
      r (3) = 0.0
      
      xa (1) = centerx - xlen * 0.5d0 + 0.5d0*shift
      xa (2) = 0 
      xa (3) = 0 
      xe (1) = centerx + xlen * 0.5d0 + 0.5d0*shift
      xe (2) = 0 
      xe (3) = 0
      rsp    = ri + rssep * (ra-ri)
      ani    = cur 
      CALL thick_coil_field (r, xa, xe, ri, rsp, ani, b) 
      hystermo = b (isel)
      xa(1) = xa(1) - shift
      xe(1) = xe(1) - shift
      CALL thick_coil_field (r, xa, xe, rsp, ra, ani, b) 
      hystermo = hystermo - b (isel)

      if( isym == 1) then
        xa (1) = centerx - xlen * 0.5d0 + 0.5d0*shift
        xa (2) = 0 
        xa (3) = 0 
        xe (1) = centerx + xlen * 0.5d0 + 0.5d0*shift
        xe (2) = 0 
        xe (3) = 0
        rsp    = ri + rssep * (ra-ri)
        ani    = cur 
        CALL thick_coil_field (r, -xa, -xe, ri, rsp, ani, b) 
        hystermo = hystermo + (-b (isel))
        xa(1) = xa(1) - shift
        xe(1) = xe(1) - shift
        CALL thick_coil_field (r, -xa, -xe, rsp, ra, ani, b) 
        hystermo = hystermo - (-b (isel))
     endif   

      
!           ----         ---                                            
! ---> diese routine + zubehoer befindet sich im modul thiceld fortran -
                                                                        
!                                                                       
      RETURN 
      END FUNCTION hystermo
