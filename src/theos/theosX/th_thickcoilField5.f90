      FUNCTION thick_co5 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> thickcoil <------                                            
!                                                                       
      implicit none
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
                        integer :: mbuf
                        integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
                        real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) :: ani, ri, ra, centerx, centery, centerz, xlen, b(3),  r(3), xa(3), xe(3) 
      REAL(8) :: pathdir(3), pathstart(3,2), orirad, oridev 
      integer :: ityp, isel, ier, nparx, maxita, npar, ini                                     
      REAL(8) epsilon 
      real :: qget, pa, qq, thick_co5, x
      double precision, parameter :: pi=4d0*atan(1d0)
      double precision :: bval, cur

       integer iadda
       common/thiadd/iadda

      COMMON / tceps / epsilon, maxita 
                                                                        
    
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'thickco5' 
         nparx = 19
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            thick_co5 = 0
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
         parnam (9) = 'componen'   ! ignore componen here 
         parnam (10) = 'maxit' 
         parnam (11) = 'epsilon'
         parnam (12) = 'startx1' 
         parnam (13) = 'starty1' 
         parnam (14) = 'startz1' 
         parnam (15) = 'startx2' 
         parnam (16) = 'starty2' 
         parnam (17) = 'startz2' 
         parnam (18) = 'oridev'     ! Abweichung in degree
         parnam (19) = 'orirad'     ! Abweichung in degree  von der radialen Orientierung
!                                                                       
         thick_co5 = 0
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
      isel          = nint(pa (9))
      maxita        = nint(pa (10)) 
      epsilon       = pa (11) 
      pathstart(1,1)    = pa (12)
      pathstart(2,1)    = pa (13)
      pathstart(3,1)    = pa (14)
      pathstart(1,2)    = pa (15)
      pathstart(2,2)    = pa (16)
      pathstart(3,2)    = pa (17)
      oridev           = pa (18)
      orirad           = pa (19)


      pathdir = [1d0, 0d0, 0d0]


      IF (isel.lt.1.or.isel.gt.3) isel = 1 
!      r (1) = x 
!      r (2) = 0.0 
!      r (3) = 0.0

       qget = 1
       call        parget('pathtyp ',qget,iadda,ier)
       ityp = nint(qget)

       qget = 1
       call        parget('bcompon ',qget,iadda,ier)
       isel = nint(qget)


 !!!! if parameter has pathtyp=2 the offaxis paths is assumed
 !!!! if bcompión  parameter is 2 radial field      

      r      = pathstart(:,ityp) + pathdir * x
 
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

      select case (isel) 
      case(1)  
         bval = b(1)
      case(2)
         bval = b(2) / cos(pi/180d0*orirad) + sin(pi/180d0*oridev) * b(1)
      end select

                                                                        
      thick_co5 = bval
!                                                                       
      RETURN 
      END FUNCTION thick_co5
