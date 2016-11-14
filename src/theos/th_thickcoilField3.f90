      FUNCTION thick_co3 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf)
!     ===================================================
!
! -------> thickcoil <------
!
!
      use theory_description
      CHARACTER(8) thnam, parnam (20)
      DIMENSION pa (20), qq (3)
                        integer :: mbuf
                        integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
                        real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) :: ani, ri, ra, centerx, centery, centerz, xlen, b(3),  r(3), xa(3), xe(3)
      REAL(8) :: pathdir(3), pathstart(3)
      REAL(8) epsilon
      double precision :: asymm, bpart(3), x1(3), x2(3)
      COMMON / tceps / epsilon, maxita

      DATA zpi / 6.283185 /
!
! ----- initialisation -----
      IF (ini.eq.0) then
         thnam = 'thickco3'
         nparx = 18
         IF (npar.lt.nparx) then
            WRITE (6, 1) thnam, nparx, npar
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)
            thick_co3 = 0
            RETURN
         ENDIF
         npar = nparx
! >>>>> describe theory with   18 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "magnetic field of a cylindrical coil arround the x-axis   "//cr//&
                                 "computed values correspond to the selected field component in Tesla "!
 !
        th_citation(idesc) = "CITATIONS OF LIT HERE"
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
         parnam (18) = 'asymm'
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " x-coordinate of coil center units=m       " !//cr//parspace//&
        th_param_desc( 2,idesc) = " y-coordinate of coil center               " !//cr//parspace//&
        th_param_desc( 3,idesc) = " z-coordinate of coil center               " !//cr//parspace//&
        th_param_desc( 4,idesc) = " length of cylindrical coil axis=x         " !//cr//parspace//&
        th_param_desc( 5,idesc) = " inner rdius of cylindrical coil units=m   " !//cr//parspace//&
        th_param_desc( 6,idesc) = " outer radius of cylindrical coil          " !//cr//parspace//&
        th_param_desc( 7,idesc) = " effective number of turns (may be fractional) " !//cr//parspace//&
        th_param_desc( 8,idesc) = " current units=Amperes                     " !//cr//parspace//&
        th_param_desc( 9,idesc) = " select 1=x, 2=y, 3=z component of field DO NOT FIT " !//cr//parspace//&
        th_param_desc(10,idesc) = " max. iterations of adaptive integration DO NOT FIT " !//cr//parspace//&
        th_param_desc(11,idesc) = " requested integration accuracy          DO NOT FIT " !//cr//parspace//&
        th_param_desc(12,idesc) = " x-component of direction of path parametrized by x " !//cr//parspace//&
        th_param_desc(13,idesc) = " y-component of direction of path parametrized by x " !//cr//parspace//&
        th_param_desc(14,idesc) = " z-component of direction of path parametrized by x " !//cr//parspace//&
        th_param_desc(15,idesc) = " x-component of start of path              " !//cr//parspace//&
        th_param_desc(16,idesc) = " y-component of start of path              " !//cr//parspace//&
        th_param_desc(17,idesc) = " z-component of start of path              " !//cr//parspace//&
        th_param_desc(18,idesc) = " left right asymmetry of winding density   " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
!
!
         thick_co3 = 0
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
      asymm         = pa (18)

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

! part 1
      x1  = xa
      x2  = xa + 0.5d0 * (xe-xa)
      CALL thick_coil_field (r, x1, x2, ri, ra, ani * (1d0-asymm)*0.5d0, bpart)
      b   = b + bpart
! part 2
      x1  = x2
      x2  = xe
      CALL thick_coil_field (r, x1, x2, ri, ra, ani * (1d0+asymm)*0.5d0, bpart)
      b   = b + bpart
!           ----         ---
! ---> diese routine + zubehoer befindet sich im modul thiceld fortran -

      thick_co3 = b (isel)
!
      RETURN
      END FUNCTION thick_co3
