      FUNCTION thick_co6 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
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
      real    :: qget
      double precision :: asymmrl,asymmcm, f1, f2, f3, bpart(3), x1(3), x2(3)
      double precision :: oridev, orirad
      COMMON / tceps / epsilon, maxita 
      double precision, parameter :: pi=4d0*atan(1d0)
      double precision :: bval, cur
      integer          :: show

                                                                  
      DATA zpi / 6.283185 / 

       integer iadda
       common/thiadd/iadda

!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'thickco6' 
         nparx = 21
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            thick_co6 = 0
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
         parnam (18) = 'asymmrl' ! reight left
         parnam (19) = 'asymmom' ! middle center
         parnam (20) = 'oridev'     ! Abweichung in degree
         parnam (21) = 'orirad'     ! Abweichung in degree  von der radialen Orientierung
!
!                                                                       
         thick_co6 = 0
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
      pathdir(1)    = pa (12)
      pathdir(2)    = pa (13)
      pathdir(3)    = pa (14)
      pathstart(1)  = pa (15)
      pathstart(2)  = pa (16)
      pathstart(3)  = pa (17)
      asymmrl       = pa (18)
      asymmcm       = pa (19)
      oridev           = pa (20)
      orirad           = pa (21)



        qget = pathstart(1)
        call        parget('startx    ',qget,iadda,ier)
        pathstart(1) = qget

        qget = pathstart(2)
        call        parget('starty    ',qget,iadda,ier)
        pathstart(2) = qget

        qget = pathstart(3)
        call        parget('startz    ',qget,iadda,ier)
        pathstart(3) = qget

        qget = isel
        call        parget('compo     ',qget,iadda,ier)
        isel = nint(qget) 

        qget = 0
        call        parget('show     ',qget,iadda,ier)
        show = nint(qget) 






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

      b = 0

      f1 = 1d0/3d0 -asymmrl -asymmcm/2
      f2 = 1d0/3d0          +asymmcm
      f3 = 1d0/3d0 +asymmrl -asymmcm/2

        

     

! part 1
      x1  = xa
      x2  = xa + (xe-xa) / 3d0
      CALL thick_coil_field (r, x1, x2, ri, ra, ani * f1, bpart) 
      if(abs(x)<1e-2 .and. show==1) then
        write(*,'("A  ",f12.6,5x,5f18.9)')x,x1(1), x2(1), ri, ra, ani * f1 * 1d-4
      endif
      b   = b + bpart
! part 2
      x1  = x2
      x2  = x1 +  (xe-xa) / 3d0
      CALL thick_coil_field (r, x1(1), x2(1), ri, ra, ani * f2, bpart) 
      if(abs(x)<1e-2 .and. show == 1) then
        write(*,'("B  ",f12.6,5x,5f18.9)')x,x1(1), x2(1), ri, ra, ani * f2 * 1d-4
      endif
      b   = b + bpart
! part 3
      x1  = x2
      x2  = x1 +  (xe-xa) / 3d0
      CALL thick_coil_field (r, x1, x2, ri, ra, ani * f3, bpart) 
      if(abs(x)<1e-2 .and. show==1) then
        write(*,'("C  ",f12.6,5x,5f18.9)')x,x1(1), x2(1), ri, ra, ani * f3 * 1d-4
      endif
      b   = b + bpart
!           ----         ---                                            
! ---> diese routine + zubehoer befindet sich im modul thiceld fortran -

                    

      select case (isel) 
      case(1)  
         bval = b(1)
      case(2)
         bval = b(2) / cos(pi/180d0*orirad) + sin(pi/180d0*oridev) * b(1)
      end select


                                           
      thick_co6 = bval
!                                                                       
      RETURN 
      END FUNCTION thick_co6
