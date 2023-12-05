      FUNCTION th_diffusion (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================
!
!
      use theory_description
      implicit none
      real, intent(in)             :: x
      real, intent(in)             :: pa (20)
      character(8), intent(inout)  :: thnam
      character(8), intent(inout)  :: parnam(20)
      integer, intent(in)          :: mbuf
      integer, intent(inout)       :: npar
      integer, intent(in)          :: ini
      integer, intent(inout)       :: nopar                 ! Anzahl der Parameter data
      character(80), intent(inout) :: napar(mbuf)           ! name des parameters n
      real, intent(inout)          :: params(mbuf)          ! value des parameters n

      real                         :: th_diffusion


      integer          :: ier, nparx, show_rr
      double precision :: tau, diff, tc, nu_subdiff, a_cross, r02, rr, ew_zero_lim, beta_q
      real             :: q

!
! ----- initialisation -----
      IF (ini.eq.0) then
         thnam = 'diffusio'
         nparx = 8
         IF (npar.lt.nparx) then
            WRITE (6, 1) thnam, nparx, npar
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)
            th_diffusion = 0
            RETURN
         ENDIF
         npar = nparx
! >>>>> describe theory with    8 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "General diffusion Ansatz                           "//cr//&
                                 "rr = ( (r02*(tau/tc)**nu_subdiff)**a_cross + (6*diff*tau)**a_cross)**(1d0/a_cross)"//cr//&
                                 "th = A0 * exp( -(q*q*rr/6d0)**beta_q )"//cr//&
                                 "   except if ( show_rr == 1 ) then"    //cr//&
                                 "th = rr                           "    
 !
        th_citation(idesc) = "CITATIONS OF LIT HERE"
!        --------------> set the number of parameters
         parnam (1) = 'intensit'
         parnam (2) = 'd0'
         parnam (3) = 'r02'
         parnam (4) = 'tref'
         parnam (5) = 'alpha'
         parnam (6) = 'a_cross'
         parnam (7) = 'beta_q'
         parnam (8) = 'show_rr'
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " prefactor A0  (not if show_rr)          " !//cr//parspace//&
        th_param_desc( 2,idesc) = " long time asymptotic diffusion constant " !//cr//parspace//&
        th_param_desc( 3,idesc) = " reference mean squared displacemengt for short times (at tref) " !//cr//parspace//&
        th_param_desc( 4,idesc) = " reference time for short time diffusion " !//cr//parspace//&
        th_param_desc( 5,idesc) = " = nu_subdiff  exponent " !//cr//parspace//&
        th_param_desc( 6,idesc) = " transition exponent between short and long time diffusuion " !//cr//parspace//&
        th_param_desc( 7,idesc) = " heterogeneity exponent  " !//cr//parspace//&
        th_param_desc( 8,idesc) = " switch to extract rr directly (if=1) " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
       th_file_param(1,idesc) = "q = q-value " 
! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
!

         th_diffusion = 0.0

         RETURN
      ENDIF
!
! ---- calculate theory here -----
      q = 0.0
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)
      IF (q.eq.0) write (6, * ) 'ERROR: q not found'

      tau         = x
      diff        = abs( pa (2) )    ! in A**2/ns
      r02         = abs( pa (3) )
      tc          = abs( pa (4) )
      nu_subdiff  = abs( pa (5) )
      a_cross     = abs( pa (6) )
      beta_q      = abs( pa (7) )
      show_rr     = nint(pa (8) )

      rr = ( (r02*(tau/tc)**nu_subdiff)**a_cross + (6*diff*tau)**a_cross)**(1d0/a_cross)
      if(show_rr == 1 ) then
         th_diffusion = rr
      else
         th_diffusion = pa(1) * exp( -(q*q*rr/6d0)**beta_q )
      endif

!
      RETURN
      END FUNCTION th_diffusion
