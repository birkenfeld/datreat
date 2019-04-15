      FUNCTION roused (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!     rouse 
!                                                                       
!       
      implicit none    
      real :: roused                                                            
      CHARACTER(8) thnam, parnam (20) 
      real :: x, pa, qq
      integer :: ier, nparx, npar, ini
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) SQ_rouseT 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL qget, tget 
      REAL kbolz 
    double precision              :: Dcm                ! centre-of-mass diffusion               
    double precision, parameter   :: a1=0.5871d0, a2=0.3282d0, a3=0.08475d0   !! approximation parameters for ROUSE F(x)
    double precision, parameter   :: x1=4.112d0,  x2=1.780d0,  x3=0.5282d0     
    double precision, parameter   :: b1=0.7484d0, b2=0.8973d0, b3=1.0d0
    double precision              :: Fx
    double precision              :: wl4, diff, nsegment, l0, xg, t , S_coh_Polymer, Q
      Parameter (kbolz = 1.380662e-23) 
                                                                        
                                                                        
     double precision, parameter ::  zpi = 6.283185d0 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'rouse_ap'       ! approximation with stretched exponentials 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            roused = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'wl4' 
         parnam (3) = 'nsegment'
         parnam (4) = 'l0'
         parnam (5) = 'diff'                                                                        
!                                                                       
         roused = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau      = x 
      a0       = pa (1) 
      wl4      = abs (pa (2) )
      nsegment = pa(3)
      l0       = pa(4)
 
      diff = abs (pa (5) ) 

      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier) 
      Q = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier)  
         temp = tget 
      ENDIF 
                                                                        
      
      if(diff .le. 0d0) then                                                                  
        diff = Wl4/(3d0*l0**2)/nsegment
      endif

!      diff = diff * 1d-9 / 1d-16 

      CALL setpar ('diff    ', sngl (diff) ,nopar ,params,napar,mbuf, ier) 
        
      t    = tau
      Dcm  = diff

      Xg   = Wl4 * (Q**4) * t / 36d0                                                                

      if( t > 0d0 ) then
      Fx = a1*exp(-(Xg/x1)**b1) +  a2*exp(-(Xg/x2)**b2) + a3*exp(-(Xg/x3)**b3)  
        S_coh_Polymer = Fx * exp(-Q*Q*Dcm*t)
      else
       S_coh_Polymer = 1d0
      endif

      roused = a0*S_coh_Polymer
                                                                  
      RETURN 
      END FUNCTION roused                             
