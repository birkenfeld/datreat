      FUNCTION th_diffusion (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
! 
!                              
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
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'd0' 
         parnam (3) = 'r02' 
         parnam (4) = 'tref' 
         parnam (5) = 'alpha' 
         parnam (6) = 'a_cross' 
         parnam (7) = 'beta_q'
         parnam (8) = 'show_rr'
                                                                        
                                                                        
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
