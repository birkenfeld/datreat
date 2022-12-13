      FUNCTION strexpsqq (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf)
!
!     strexpsqq 
!     stretched exponential function
!     y=amplitu*exp([x/(tau0 * q**qexp)]^beta)
!     units: nothing special
!     Reference : unknown author in unknown Journal
!     author r.biehl
!     tested
      use theory_description
		implicit none          ! to make it easier to understand
		real x                                     ! xvalue
		real    pa (20)                            ! parameter values of theory
		CHARACTER(8) thnam                         ! theory    name in datreat
		CHARACTER(8) parnam (20)                   ! parameter namen of theory
		integer      npar                          ! maximum number of theory parameters in main program
		integer      ini                           ! flag initialisation or calculation
 		integer, intent(inout) :: nopar            ! number of Parameter in data
		character*80, intent(inout) :: napar(mbuf) ! name of  parameters(n)
		real, intent(inout) :: params(mbuf)        ! value of parameters(n)
		integer :: mbuf                            ! maximum length of data arrays in calling main program

		! define local parametrs here
		real amp, tau0, bet, qexp, q, tau, pi, arg
		real strexpsqq
                real :: tauave
		integer nparx                             ! number of theory parameters
        integer ier                               ! error variable 
        double precision :: sq, rr, qq
        double precision :: qqr(1), scl, gam, eta 
        integer          :: mode, ierr

        integer          :: iadda 
        common/thiadd/iadda
		
        DATA pi/3.141592653589793/
!
! ----- initialisation of theory-----
      IF (ini.eq.0) then
         thnam = 'strexsqq'     ! name of theory, max 8*char
         nparx = 9              ! number of theory parameters
! >>>>> describe theory with    4 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "streched exponential                                  "//cr//&
                                 "   a0 * exp(-(t/tau)**beta)                           "//cr//&
                                 "        with tau = tau0*(q**qexp)/s(q)                "//cr//&
                                 "   S(Q) model paraeters are:                          "//cr//&
                                 "   eta = volume fraction                              "//cr//&      
                                 "   SCL=characteristic length of the interaction       "//cr//&      
                                 "           screening length for mode=1,2              "//cr//&     
                                 "           the width of the potential well for        "//cr//&     
                                 "                           mode=3                     "//cr//&     
                                 "   GAMMA=strength of the interaction                  "//cr//&     
                                 "           mode=1,2 GEK=gamma*exp(-r/scl) contact pot."//cr//&     
                                 "           mode=3 depth of the potential well         "//cr//&     
                                 "           in units of kT, + ==> repulsive            "//cr//&     
                                 "                           - ==> attractive           "//cr//&     
                                 "   R=outer radius of the spheres                      "//cr//&     
                                 "   MODE=   1       H.-P. form                         "//cr//&     
                                 "           2       Shieu-Chen form                    "//cr//&     
                                 "           3       Sharma-Sharma form                 "//cr//&     
                                 "           4       Critical diverg.                   "//cr//&     
                                 "           5       none (sq=1.)                       "    
                                 
 !
        th_citation(idesc) = " "
 
         IF (npar.lt.nparx) then
            WRITE (6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)") thnam, nparx, npar
            strexpsqq = 0
            RETURN
         ENDIF
!        --------------> set the number of parameters and add a short comment here
         npar = nparx              ! to return number of parameters
         parnam (1) = 'amplitu'    ! amplitude   should be 1
         parnam (2) = 'tau'        ! relaxation time
         parnam (3) = 'beta'       ! stretched exponent
         parnam (4) = 'qexp'       ! exponent of q
         parnam (5) = 'eta'        ! volumen fraction
         parnam (6) = 'scl'        ! characteristic length
         parnam (7) = 'gamma'      ! strebgth of interactio
         parnam (8) = 'r'          ! radius
         parnam (9) = 'mode'       ! selection of S(Q) type
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " prefactor a0" !//cr//parspace//&
        th_param_desc( 2,idesc) = " tau0 in time units of the x-axis (abs is applied) " !//cr//parspace//&
        th_param_desc( 3,idesc) = " streching exponent (abs is applied) " !//cr//parspace//&
        th_param_desc( 4,idesc) = " q-exponent (tau(q)=tau0*(q**qexp)   " !//cr//parspace//&
        th_param_desc( 5,idesc) = " volumen fraction                    " !//cr//parspace//&
        th_param_desc( 6,idesc) = " characteristic length               " !//cr//parspace//&
        th_param_desc( 7,idesc) = " strebgth of interactio              " !//cr//parspace//&
        th_param_desc( 8,idesc) = " radius                              " !//cr//parspace//&
        th_param_desc( 9,idesc) = " selection of S(Q) type              " !//cr//parspace//

! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
       th_file_param(1,idesc) = "q > q-value for computation of tau0, units in q0 : tau0=tau(q0)  " 
 ! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
       th_out_param(1,idesc)  = "tauave = average tau = tau * Gamma(1/beta)/beta "  
!
!
         strexpsqq = 0
         RETURN
      ENDIF ! of initialisation
			
      !     get the parameters from the data
      !     searched by the name of the parameters and assigned to variables
      !     old was call to  parget from data
		call getpar('q       ',q,nopar ,params,napar,mbuf, ier)
		if (ier.gt.0) write(*,*) 'q not found in data parameters   ier=',ier      ! only test it
			
		! get the parameters from the Theory    (better names)
      amp  = pa (1)
      tau0 = abs (pa (2) )
      bet  = abs (pa (3) )
      qexp = pa (4)
      eta  = abs(pa(5))
      scl  = abs(pa(6))
      gam  = abs(pa(7))
      rr   = abs(pa(8))
      mode = nint(pa(9))
   
! ----NOW calculate theory here -----
      q = abs (q)
      qq = q
      qqr(1) = qq * rr * 2
      CALL sofq (qqr, sq, 1, eta, scl, gam, rr, mode, ierr) 
      if(ierr .ne. 0 ) write(6,*)"warning: sofd determination failed"
 !     tau = (tau0 * q**qexp) * peryev( qq, rr, den, eps ) 
      tau = (tau0 * q**qexp) * sq 

!      IF (bet.lt.0.0) bet = 0.0
!      IF (bet.gt.1.0) bet = 1.0

      arg = abs (x / tau) **bet
!      IF (arg.lt. - 30.0) arg = - 30.0
!      IF (arg.gt.30.0) arg = 30.0
!    FINALY the result of this funktion  is ============
      strexpsqq = amp * exp ( - arg)
!    for demonstration of calculation of a new parameter accessible by ?? +<name> and writen with save together with all data
      tauave = tau0 * Gamma(1.0/bet)/bet 
      call parset('tauave  ',tauave,iadda,ier)		
			
      RETURN
      END FUNCTION strexpsqq
