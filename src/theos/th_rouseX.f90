      FUNCTION rouseX (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rousex <-------- with wavelength spread ---                                              
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) SQ_rouseT 
      REAL(8) a0, sum, sumnorm, qzz, fn 
      REAL(8) epsilon, diff 
      REAL qget, tget , vget
                                                                        

      double precision :: lam_min, lam_max, lam0, lam, dlam, tau0
      double precision :: muelam
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'rousex' 
         nparx = 10
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            rousex = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'xi_frict' 
         parnam (3) = 'b_segmnt' 
         parnam (4) = 'epsilon ' 
         parnam (5) = 'temp    ' 
         parnam (6) = 'com_diff' 
         parnam (7) = 'lam_min' 
         parnam (8) = 'lam_max' 
         parnam (9) = 'n_lam' 
         parnam (10)= 'muelam'   ! exponent for intensity variation with lambda i.e. Phi~lam**-4,SQ~lam**2 => -2  
                                                                        
!                                                                       
         rousex = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau = x 
      a0      = pa (1) 
      xi      = abs (pa (2) ) 
      b       = abs (pa (3) ) 
      epsilon = abs (pa (4) ) 
      temp    = pa (5) 
                              ! in cm**2/sec                            
      diff    = abs (pa (6) ) 
      lam_min = pa (7) 
      lam_max = pa (8)
      nqw     = NINT(pa(9))
      muelam  = pa(10) 
                                                                        
                                       ! in A**2/ns                     
      diff = diff * 1d-9 / 1d-16 
                                                                        
      IF (epsilon.eq.0.0d0) epsilon = 1.0d-3 
                                                                        
      qget = 0.01 
      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier)  
      qz = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier) 
         temp = tget 
      ENDIF
 
      lam0 = (lam_max+lam_min)/2
      if(lam_min.eq.0d0) then
!!! if lambda-min and lambda-max are present in the parameter list we take them from that source
        vget    = 1
        CALL getpar ('lam_min ', vget ,nopar ,params,napar,mbuf, ier)
        lam_min = vget
        vget    = 1        
        CALL getpar ('lam_max ', vget ,nopar, params,napar,mbuf, ier)  
        lam_max = vget
        vget    = 1      
        CALL getpar ('lambda  ', vget ,nopar ,params,napar,mbuf, ier)
        lam0    = vget        
      endif
                                                                  
                                                                         
      tau0 = tau
      dlam = (lam_max-lam_min)/(2*nqw) 


      IF (nqw.eq.0) then 
        dlam =0d0
      else 
        dlam = (lam_max-lam_min)/(2*nqw)
      ENDIF 
                        

      rousex = 0 
      sum = 0 
      sumnorm = 0 
                                                
      DO i = -nqw, nqw
 
        lam = lam0+i*dlam

        qzz = qz *lam0/lam
        tau = tau0*(lam/lam0)**3

!        write(6,*)i,lam,qzz,tau


!
         fn = (lam/lam0)**muelam   ! here: possible intensity variation due lambda spectrum times S(Q)

         sumnorm = sumnorm + fn 
! --- include center of mass diffusion ---                              
         a = fn * a0 * dexp ( - qzz * qzz * diff * tau) 
                                                                        
         IF (pa (4) .lt.0) then 
            WRITE (6, * ) 'Full rouse computation from scratch!' 
            sum = sum + a * SQ_rouse (tau, qzz, temp, xi, b, epsilon) 
         ELSE 
            sum = sum + a * SQ_rouseT (tau, qzz, temp, xi, b, epsilon) 
         ENDIF 


      enddo 
                 

      tau = tau0
                                                       
      rousex = sum / sumnorm 
!                                                                       
      RETURN 
      END FUNCTION rousex   
