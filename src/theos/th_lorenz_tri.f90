      FUNCTION th_lortri (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
! Lorenzian (diffusiontype) convoluted with triangular shape with FWHM = delta 
! 
!
!          
      use    PhysicalConstantsPlus  
                                          
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 /

      double precision :: D, omega, dE, e1, e0, de0, dee
      double precision :: lambda0, lambda1
      double precision :: delta_lambda
      real             :: q, temp, alam, dlam
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'lortri' 
         nparx = 4 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_lortri = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplit' 
         parnam (2) = 'd'       ! in units derived form q and x units
         parnam (3) = 'fwhm'    ! resolution widthd 
         parnam (4) = 'offset'  ! omega offset                                                               
                                                                        
         th_lortri = 0.0 
                                                                        
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q = 0.0 
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)  
      IF (q.eq.0) write (6, * ) 'ERROR: q not found' 

      CALL getpar ('angle   ', angle_2tht ,nopar ,params,napar,mbuf, ier)  
      if(ier.ne.0) then
        write(6,*)'ERROR: Parameter angle is missing in parameterlist of record!'
        th_lortri = 0
        return
      endif

      CALL getpar ('lambda  ', alam ,nopar ,params,napar,mbuf, ier)
      if(ier.ne.0) then
        write(6,*)'ERROR: Parameter lambda is missing in parameterlist of record!'
        th_lortri = 0
        return
      endif

      CALL getpar ('dlambda ', dlam ,nopar ,params,napar,mbuf, ier)
      if(ier.ne.0) then
        write(6,*)'ERROR: Parameter dlambda (channel dist) is missing in parameterlist of record!'
        th_lortri = 0
        return
      endif

      CALL getpar ('temp    ', temp ,nopar ,params,napar,mbuf, ier)
      if(ier.ne.0) then
        write(6,*)'ERROR: Parameter temp is missing in parameterlist of record!'
        th_lortri = 0
        return
      endif
!
! and assume x is lambda
!
      alam = alam - pa(4)

      lambda0        = x    * 1d-10                      ! das ist lambda_final
      lambda1        = alam * 1d-10                      ! das ist lambda_in
      delta_lambda   = dlam * 1d-10
      e0      = NeutronEnergy_fromLambda(lambda0)
      e1      = NeutronEnergy_fromLambda(lambda1)
      dE      = e1-e0                                    ! hquer*omega = E_in - E_final
      de0     = NeutronEnergy_fromLambda(lambda1-delta_lambda/2)-NeutronEnergy_fromLambda(lambda1+delta_lambda/2)
      dee     = NeutronEnergy_fromLambda(lambda0-delta_lambda/2)-NeutronEnergy_fromLambda(lambda0+delta_lambda/2)
      dbf     = exp(0.5d0*dE/(temp*Boltzmannkonstante))

!      omega = x - pa(4)

      omega = dE*2*Pi/Planckkonstante / 1d9




      D     = pa (2)**2        ! square of diffusion constatnt ! 
      delta = abs (pa (3) )

      q = ((2*Pi/alam)**2+(2*Pi/x)**2 - 2*COS(angle_2tht*Pi/180d0)*(2*Pi/alam)*(2*Pi/x))  ! q-squared !!

! convolution of a lorenzian (simple diffusion, Gamma proto D*q**2) 
! with a simple triangular resolution function with fwhm = delta
   
      th_lortri = (-q*sqrt(D)*log(D*q**2+delta**2+2*omega*delta+omega**2)-   &
                    q*sqrt(D)*log(D*q**2+delta**2-2*omega*delta+omega**2)+   & 
                    2*q*sqrt(D)*log(D*q**2+omega**2)+2*delta*atan((delta+omega)/sqrt(D)/q)+ &
                    2*omega*atan((delta+omega)/sqrt(D)/q)+2*delta*atan((-omega+delta)/q/sqrt(D))- &
                    2*omega*atan((-omega+delta)/q/sqrt(D))-4*omega*atan(omega/q/sqrt(D)))  &
                    /0.3141592653589793D1/delta**2/2d0 &
                 * pa(1)

!          
!    Detailed balance factor, kinematic factor and channelwidth  if available
!          
!
! and assume x is lambda
!
      lambda0        = x    * 1d-10                      ! das ist lambda_final
      lambda1        = alam * 1d-10                      ! das ist lambda_in
      delta_lambda   = dlam * 1d-10
      e0      = NeutronEnergy_fromLambda(lambda0)
      e1      = NeutronEnergy_fromLambda(lambda1)
      dE      = e1-e0                                    ! hquer*omega = E_in - E_final
      de0     = NeutronEnergy_fromLambda(lambda1-delta_lambda/2)-NeutronEnergy_fromLambda(lambda1+delta_lambda/2)
      dee     = NeutronEnergy_fromLambda(lambda0-delta_lambda/2)-NeutronEnergy_fromLambda(lambda0+delta_lambda/2)
      dbf     = exp(0.5d0*dE/(temp*Boltzmannkonstante))

      th_lortri = th_lortri*dbf * abs(dee/de0) * lambda1/lambda0     
!                                                ==> (kf/ki)=(lam_i/lam_f)
!                                    ====> channel width transformation
!                           ===> detailed balance factor
!
!     envetually missing: detector efficiency correction
!     maybe we do it elsewhere
!


!      write(6,'(11e14.7)')x,temp,dbf, abs(dee/de0), lambda1, delta_lambda, e0, e1, de, de0, dee
              
                                                    
      RETURN 
      END FUNCTION th_lortri 
   
