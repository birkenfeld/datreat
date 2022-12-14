      FUNCTION ashelli (r) 
!      -------------------                                              
! -- integrand function for shell-scattering ---                        
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (mcashi = 10) 
      COMMON / cashi / qcashi, r0cashi, w0cashi, r1cashi, w1cashi,      &
      rhocashi, acashi (mcashi), ecashi (mcashi), ncashi                
                                                                        
      sum = 0.0d0 
      DO i = 1, ncashi 
      IF (acashi (i) .ne.0.d0) sum = sum + acashi (i) * r**ecashi (i) 
      enddo 
                                                                        
! -- aeussere kante mit fermi-funktion abschneiden ..                   
      sum = sum / (1.0d0 + dexp ( (r - r0cashi) / w0cashi) ) 
! -- innere kante der schale mit fermifunktion anschalten...            
! -- und innere streulaengendichte dabei abschalten....                 
      fi = 1.0d0 / (1.0d0 + dexp ( (r1cashi - r) / w1cashi) ) 
      sum = sum * fi + (1.d0 - fi) * rhocashi 
                                                                        
      qr = qcashi * r 
      IF (qr.lt.1.0d-6) then 
!ccc     write(6,*)'ashelli: qr=',qr,' replaced by qr=1e-6'             
         qr = 1.0d-6 
      ENDIF 
                                                                        
      ashelli = sum * r * r * sin (qr) / qr 
!                      !!!                                              
!--- vergessen! eingefuegt 12.11.92 16:34                               
      RETURN 
      END FUNCTION ashelli                          
