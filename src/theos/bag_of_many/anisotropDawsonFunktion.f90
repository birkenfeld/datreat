      FUNCTION anidawson (q, R1, R2) 
!      ---------------------------                                      
! Anisotrope Dawson-Funktion mit Skalierung auf 1 bei q=0               
!                                                                       
                                                                        
      DOUBLEPRECISION pi 
      PARAMETER (pi = 3.141592654d0) 
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION q, R1, R2, phi, anidawson 
      DOUBLEPRECISION xx, DDAWS, adwsker 
      DOUBLEPRECISION qq, rr1, rr2 
      DOUBLEPRECISION adapint, erraccu, errest, gu, go, fn 
                                                                        
      COMMON / cadws / qq, rr1, rr2 
      EXTERNAL adwsker 
                                                                        
      maxit = 1000 
      errest = 1d-6 
      gu = 0.0d0 
      go = pi 
      fn = 1.0d0 / (go - gu) 
                                                                        
      qq = q 
      rr1 = r1 
      rr2 = r2 
                                                                        
      anidawson = fn * adapint (adwsker, gu, go, errest, maxit, erraccu) 
                                                                        
      RETURN 
      END FUNCTION anidawson                        
                                                                        
                                                                        
      FUNCTION adwsker (phi) 
!      ---------------------                                            
!  Kernel for integration for anisotrop Dawson function !!              
! -------------------------------------------------------               
                                                                        
      IMPLICIT none 
      DOUBLEPRECISION q, R1, R2, phi, xx, ddaws, adwsker 
      COMMON / cadws / q, R1, R2 
                                                                        
      xx = sqrt ( (q * R1 * cos (phi) ) **2 + (q * R2 * sin (phi) ) **2)&
      / 2                                                               
                                 ! imslroutine                          
      adwsker = ddaws (xx) / xx 
                                                                        
      RETURN 
      END FUNCTION adwsker                          
