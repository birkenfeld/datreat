      FUNCTION thick_co2 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> thickcoil <------                                            
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
                        integer :: mbuf
                        integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
                        real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) :: ani, ri, ra, centerx, centery, centerz, xlen, b(3),  r(3), xa(3), xe(3) 
      REAL(8) :: pathdir(3), pathstart(3)                                            
      REAL(8) epsilon 
      COMMON / tceps / epsilon, maxita 
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'thickco2' 
         nparx = 17
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            thick_co2 = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'centerx' 
         parnam (2) = 'centery' 
         parnam (3) = 'centerz' 
         parnam (4) = 'length' 
         parnam (5) = 'ri' 
         parnam (6) = 'ra' 
         parnam (7) = 'nw' 
         parnam (8) = 'cur' 
         parnam (9) = 'componen' 
         parnam (10) = 'maxit' 
         parnam (11) = 'epsilon'
         parnam (12) = 'dirx' 
         parnam (13) = 'diry' 
         parnam (14) = 'dirz' 
         parnam (15) = 'startx' 
         parnam (16) = 'starty' 
         parnam (17) = 'startz' 
!                                                                       
         thick_co2 = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      centerx       = pa (1) 
      centery       = pa (2) 
      centerz       = pa (3) 
      xlen          = pa (4) 
      ri            = pa (5) 
      ra            = pa (6) 
      ani           = pa (7) 
      cur           = pa (8) 
      isel          = pa (9) + 0.1 
      maxita        = pa (10) + 0.1 
      epsilon       = pa (11) 
      pathdir(1)    = pa (12)
      pathdir(2)    = pa (13)
      pathdir(3)    = pa (14)
      pathstart(1)  = pa (15)
      pathstart(2)  = pa (16)
      pathstart(3)  = pa (17)

      pathdir    = pathdir / sqrt(dot_product(pathdir,pathdir))

      IF (isel.lt.1.or.isel.gt.3) isel = 1 
!      r (1) = x 
!      r (2) = 0.0 
!      r (3) = 0.0

      r      = pathstart + pathdir * x
 
      xa (1) = centerx - xlen * 0.5d0 
      xa (2) = centery 
      xa (3) = centerz 
      xe (1) = centerx + xlen * 0.5d0 
      xe (2) = centery 
      xe (3) = centerz 
      ani = ani * cur 
      CALL thick_coil_field (r, xa, xe, ri, ra, ani, b) 
!           ----         ---                                            
! ---> diese routine + zubehoer befindet sich im modul thiceld fortran -
                                                                        
      thick_co2 = b (isel)
!                                                                       
      RETURN 
      END FUNCTION thick_co2
