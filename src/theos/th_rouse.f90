      FUNCTION rouse (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rouse <--------                                              
!            
      use theory_description                                                           
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) SQ_rouseT 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff 
      REAL qget, tget 
                                                                        
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'rouse' 
         nparx = 7 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            rouse = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters    

        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "Rouse single chain S(Q,t) deGennes integral expression   "//cr//&
                                 "  with extra center-of-mass diffusion exp(-D_cm*d**2*t) as factor"//cr//&
                                 "  with optional averaging over q-width"
 !
        th_citation(idesc) = "Doi and Edwards Book and refernces therein"
           

               
         parnam (1) = 'intensit' 
         parnam (2) = 'xi_frict' 
         parnam (3) = 'b_segmnt' 
         parnam (4) = 'epsilon ' 
         parnam (5) = 'temp    ' 
         parnam (6) = 'com_diff' 
         parnam (7) = 'q_width' 
!
         th_param_desc(1,idesc) = "prefactor "  
         th_param_desc(2,idesc) = "friction coefficient:  Wl4   = 36 *   (kb * temp)*b*b/(12*xi)"  
         th_param_desc(3,idesc) = "segment length b in A"  
         th_param_desc(4,idesc) = "integration accuracy parameter. DO NOT FIT. " //cr//parspace//& 
                                  "if epsilon < 0, numerical integration with eps=1e-8 is enforced"
         th_param_desc(5,idesc) = "temperature overrides parameter in data records"//cr//parspace//&
                                  "set to 0 if  extraction from parameters in records is desired"  
         th_param_desc(6,idesc) = "center-of-mass diffusion in units cm**2/sec (assuming A, ns units for q,t)"  
         th_param_desc(7,idesc) = "if non-zero a 15 points averaging over the q-range is perfromed."  

         th_file_param(:,idesc) = " "
         th_file_param(1,idesc) = " q     =  q value of the S(q,t), units A**-1"
         th_file_param(2,idesc) = " temp  =  T in K "

         th_out_param(:,idesc)  = " "
         th_out_param(1,idesc)  = " wl4 = Rouse rate in A**4/ns"
                                                                         
!                                                                       
         rouse = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau = x 
      a0 = pa (1) 
      xi = abs (pa (2) ) 
      b = abs (pa (3) ) 
      epsilon = abs (pa (4) ) 
      temp = pa (5) 
                              ! in cm**2/sec                            
      diff = abs (pa (6) ) 
      q_width = pa (7) 
                                                                        
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
                                                                        
                                                                        
      rouse = 0 
      sum = 0 
      sumnorm = 0 
      nqw = 15 
      dqw = 4 * q_width / nqw 
      IF (q_width.eq.0) then 
         nqw = 0 
         q_width = 1.0d0 
      ENDIF 
                                                                        
      DO i = - nqw, nqw 
      qzz = qz + i * dqw 
      IF (qzz.gt.0) then 
         fn = dexp ( - (i * dqw / q_width) **2) 
         sumnorm = sumnorm + fn 
                                                                        
! --- include center of mass diffusion ---                              
         a = fn * a0 * dexp ( - qzz * qzz * diff * tau) 
                                                                        
         IF (pa (4) .lt.0) then 
            WRITE (6, * ) 'Full rouse computation from scratch!' 
            sum = sum + a * SQ_rouse (tau, qzz, temp, xi, b, 1d-8) 
         ELSE 
            sum = sum + a * SQ_rouseT (tau, qzz, temp, xi, b, epsilon) 
         ENDIF 
      ENDIF 
      enddo 
                                                                        
      rouse = sum / sumnorm 
!                                                                       
      RETURN 
      END FUNCTION rouse   
