      FUNCTION debye2 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!  Debye function
!


      use theory_description
      IMPLICIT none

      DOUBLE PRECISION :: Pi, Navog
      Parameter (Pi = 3.141592654d0)
      Parameter (Navog = 6.022045d23)


      CHARACTER(8) thnam, parnam (20)
      REAL debye2, x, pa, qq, zpi, xh, vol_frac
      INTEGER ini, npar, nparx
      DIMENSION pa (20), qq (3)
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character(len=80), intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
      DATA zpi / 6.283185 /

      DOUBLE PRECISION :: v, nphi

      DOUBLE PRECISION :: ampli, q, rg, qrg, gam_exp, sq, sqi
      COMMON / cdebik / gam_exp, qrg

      DOUBLE PRECISION :: mw, bpoly, drho, conc, rhopoly, fac
      REAL :: bsolv = 0d0


      DOUBLE PRECISION :: dbik, adapint, erra
      EXTERNAL dbik


      INTEGER  ier,  iot
!
!
! ----- initialisation -----
      IF (ini.eq.0) then
      thnam = 'debye2  '
         nparx = 8
         IF (npar.lt.nparx) then
            WRITE (6, 1) thnam, nparx, npar
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)
            debye2 = 0
            RETURN
         ENDIF
         npar = nparx
! >>>>> describe theory with    8 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "Debye polymer structure factor with molecular parameters and a2 correction"//cr//&
                                 "for non-Gaussian chain statistics the parameter gam_exp may deviate from 1  "//cr//&
                                 "for gam_exp=0:  sq = (2 / (qrg**2) ) * (qrg - 1 + exp ( - qrg) )            "//cr//&
                                 "              with qrg = q**2 * Rg**2                                     "//cr//&
                                 "              with fac = conc*(1-conc) * mw / rhopoly / Navog             "//cr//&
                                 "a2 corr is: sqi = intensit * drho**2 / (fac * sq) + v * conc * conc       "//cr//&
                                 "if gam_exp not =0 instead of the Debye form ad direct integration is used   "//cr//&
                                 "   sq = int(0..1){2 * (1 - u) * exp ( - (u**gam_exp) * qrg)}                "//cr//&
                                 "        (compare also benoitmf theory) " 
 !
        th_citation(idesc) = "CITATIONS OF LIT HERE"
!        --------------> set the number of parameters
      parnam (1)  = 'intensit'
      parnam (2)  = 'rg      '
      parnam (3)  = 'gamma   '
      parnam (4)  = 'molwght '
      parnam (5)  = 'density '
      parnam (6)  = 'bpolym  '
      parnam (7)  = 'v       '
      parnam (8) = 'volfrac '
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " intensity (amplitide) prefactor  " !//cr//parspace//&
        th_param_desc( 2,idesc) = " Radius of gyration (units as 1/q) " !//cr//parspace//&
        th_param_desc( 3,idesc) = " gam_exp = 2*chain expansion coefficient " //cr//parspace//&
                                  " if gam_exp = 0, simple Debye is assumed "
        th_param_desc( 4,idesc) = " Molecular weight of polymer chain " !//cr//parspace//&
        th_param_desc( 5,idesc) = " Density of polymer " !//cr//parspace//&
        th_param_desc( 6,idesc) = " scattering length density " !//cr//parspace//&
        th_param_desc( 7,idesc) = " effective volume entering the a2 expression " !//cr//parspace//&
        th_param_desc( 8,idesc) = " concentration in terms of volume fractio " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
       th_file_param(1,idesc) = "bsolv   solvent scattering length density" 
! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
!
!
         debye2 = 0
         RETURN
      ENDIF
!
! ---- calculate theory here -----
      ampli   =     (pa (1) )
      rg      = abs (pa (2) )
      gam_exp   = abs (pa (3) )
      mw      = abs (pa (4) )
      rhopoly = abs (pa (5) )
      bpoly   = pa (6)
      v       = pa (7)
      conc    = abs (pa (8) )



      CALL getpar ('bsolv   ', bsolv, nopar ,params,napar,mbuf, ier)
      drho = bpoly - bsolv

      fac = conc*(1d0-conc) * mw / rhopoly / Navog

      q = x
      qrg = q * q * rg * rg

      IF (gam_exp.eq.0.0d0) then
         sq = (2 / (qrg**2) ) * (qrg - 1 + exp ( - qrg) )
      ELSE
         sq = adapint (dbik, 0.0d0, 1.0d0, 1d-7, 20000, erra)
      ENDIF

! Zimm formula for interacting chains
      sqi = 1.0d0 / (fac * sq) + v * conc * conc
      IF (abs (sqi) .lt. epsilon(sqi)) sqi = epsilon(sqi)

      debye2 = ampli * drho**2 / sqi

      RETURN
      END FUNCTION debye2


      DOUBLEPRECISION function dbik (u)
!       ---------------------------------
!
      DOUBLEPRECISION u
      DOUBLEPRECISION qrg, gam_exp
      COMMON / cdebik / gam_exp, qrg

      dbik = 2 * (1 - u) * exp ( - (u**gam_exp) * qrg)

      RETURN
      END FUNCTION dbik
