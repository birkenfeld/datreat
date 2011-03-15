      FUNCTION th_diff123(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ==============================================================================               
!                                                                       
! -------> diff in 1 2 or3 D <--------                                            
!                                                                       
!   latex(u);
!   \sin \left( \theta \right) 
!   latex(y1d);
!   {e^{- \left( \cos \left( \theta \right)  \right) ^{2}{q}^{2}t}}
!   latex(y);
!   {e^{-{q}^{2} \left( \sin \left( \theta \right)  \right) ^{2}t}}
!   
!                     
!  latex(d2d);
!  {\frac {{\it dawson} \left( q\sqrt {t} \right) }{q\sqrt {t}}}
!  latex(d1d);
!  1/2\,{\frac {{\it erf} \left( q\sqrt {t} \right) \sqrt {\pi }}{q\sqrt 
!  {t}}}
!  latex(int(u*y1d,theta=0..Pi/2));
!  1/2\,{\frac {{\it erf} \left( q\sqrt {t} \right) \sqrt {\pi }}{q\sqrt 
!  {t}}}
!  latex(int(u*y,theta=0..Pi/2));
!                              
!  1/2\,{\frac {\sqrt {\pi }{e^{-{q}^{2}t}}{\it erfi} \left( q\sqrt {t}
!   \right) }{q\sqrt {t}}}
!  

                                             
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			DATA zpi / 6.283185 / 

      double precision ddaws,  derf     ! Slatec
      double precision s15aff           ! Nag Dawson Integral
      double precision s15aef           ! Nag erf
      double precision d, xx, ff

      double precision pi      
      parameter (pi=3.141592654d0)
                                                                 
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'diff123' 
         nparx = 3 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_diff123 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplit' 
         parnam (2) = 'd' 
         parnam (3) = 'dim' 
                                                                        
                                                                        
         th_diff123 = 0.0 
                                                                        
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q = 0.0 
      CALL getpar ('q       ', q,nopar ,params,napar,mbuf, ier)  
      IF (q.eq.0) write (6, * ) 'ERROR: q not found' 
     
      

                                                                  
      t = x 
      D = abs (pa (2) ) 
      
      xx = sqrt(q*q*D*t)     

      if(xx.lt.1d-6) xx = 1d-6  

      idim = NINT(pa(3))
      if(idim.lt.1) idim = 3
      if(idim.gt.3) idim = 3

      ifail = 0
      if(idim.eq.1) then
!       ff = (sqrt(Pi)/2)*derf(xx)/xx  
        ff = (sqrt(Pi)/2)*s15aef(xx,ifail)/xx  
      endif
      if(idim.eq.2) then
!       ff = dawson(xx)/xx  
        ff = s15aff(xx,ifail)/xx  
      endif                                      
      if(idim.eq.3) then
        ff = exp(-xx*xx)  
      endif
                                        
                                                                        
                                                                        
      th_diff123 = pa (1) * ff 
!                                                                       
      RETURN 
      END FUNCTION th_diff123       
