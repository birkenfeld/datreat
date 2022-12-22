      FUNCTION m_exp_q (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!   multi exp decay exp_q 
!!----------------------------------------------------------            
!! 2-Exponentials with Q-dependend Amplitude-Ratio and Decay            
!! Amplitude = p1                                                       
!! Aratio    = p2 + p3*Q**p4+p5*Q**p6+p7*Q**p8                          
!! Tau1      = p9 + p10*Q**p11 + p12*Q**p13                             
!! Beta1     = p14                                                      
!! Tau2      = p15+ p16*Q**p17 + p18*Q**p19                             
!! Beta2     = p20                                                      
!!----------------------------------------------------------            
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'exp_q' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            m_exp_q = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitud' 
         parnam (2) = 'aratioc0' 
         parnam (3) = 'aratioc1' 
         parnam (4) = 'aratioe1' 
         parnam (5) = 'aratioc2' 
         parnam (6) = 'aratioe2' 
         parnam (7) = 'aratioc3' 
         parnam (8) = 'aratioe3' 
      parnam (9)  = 'tau1c0  ' 
      parnam (10)  = 'tau1c1  ' 
      parnam (11)  = 'tau1e1  ' 
      parnam (12)  = 'tau1c2  ' 
      parnam (13)  = 'tau1e2  ' 
      parnam (14)  = 'beta1   ' 
      parnam (15)  = 'tau1c0  ' 
      parnam (16)  = 'tau2c1  ' 
      parnam (17)  = 'tau2e1  ' 
      parnam (18)  = 'tau2c2  ' 
      parnam (19)  = 'tau2e2  ' 
      parnam (20)  = 'beta2   ' 
!                                                                       
         m_exp_q = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
      amplitud = pa (1) 
      aratioc0 = pa (2) 
      aratioc1 = pa (3) 
      aratioe1 = pa (4) 
      aratioc2 = pa (5) 
      aratioe2 = pa (6) 
      aratioc3 = pa (7) 
      aratioe3 = pa (8) 
      tau1c0 = pa (9) 
      tau1c1 = pa (10) 
      tau1e1 = pa (11) 
      tau1c2 = pa (12) 
      tau1e2 = pa (13) 
      beta1 = pa (14) 
      tau2c0 = pa (15) 
      tau2c1 = pa (16) 
      tau2e1 = pa (17) 
      tau2c2 = pa (18) 
      tau2e2 = pa (19) 
      beta2 = pa (20) 
                                                                        
      q = 1.0 
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)  
      q = abs (q) 
      t = abs (x) 
                                                                        
                                                                        
      aratio = aratioc0 + aratioc1 * q**aratioe1 + aratioc2 * q**       &
      aratioe2 + aratioc3 * q**aratioe3                                 
                                                                        
      a1 = 1.0 / (1.0 + aratio) 
      a2 = 1.0 - a1 
                                                                        
      tau1 = tau1c0 + tau1c1 * q**tau1e1 + tau1c2 * q**tau1e2 
      tau2 = tau2c0 + tau2c1 * q**tau2e1 + tau2c2 * q**tau2e2 
                                                                        
      IF (tau1.gt.0.0d0) then 
         arg1 = - (t / tau1) **beta1 
      ELSE 
         arg1 = - 30.0 
      ENDIF 
                                                                        
      IF (tau2.gt.0.0d0) then 
         arg2 = - (t / tau2) **beta2 
      ELSE 
         arg2 = - 30.0 
      ENDIF 
                                                                        
      IF (arg1.lt. - 30.0) arg1 = - 30.0 
      IF (arg2.lt. - 30.0) arg2 = - 30.0 
                                                                        
                                                                        
      m_exp_q = amplitud * (a1 * exp (arg1) + a2 * exp (arg2) ) 
!                                                                       
      RETURN 
      END FUNCTION m_exp_q      
