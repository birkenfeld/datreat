!*** schulz verteilung ***                                              
      FUNCTION pschulz (r, rave, z) 
!      ----------------------------                                     
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      pschulz = 0.0 
      z1 = z + 1 
      IF (z1.lt.0.0) return 
      IF (rave.lt.0.0) return 
!      pschulz = (z1/rave)**z1 * r**z * exp(-z1*r/rave) / gamma(z1)     
      arg1 = dlog (z1 / rave) 
      arg2 = dlog (r) 
      arg3 = - z1 * r / rave 
      arg4 = - dlngam (z1) 
      pschulz = dexp (z1 * arg1 + z * arg2 + arg3 + arg4) 
      RETURN 
      END FUNCTION pschulz                          
!**                                                                     
!** hilfsfunktionen: averaging ueber kugelformfaktoren **               
!** vg. m.kotlarchyk,s.h.chen,j.chem.phys 79 (1983) 2461                
!**                                                                     
      FUNCTION g1schul (q, rave, z) 
!      ----------------------------                                     
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      z1 = z + 1 
      z2 = z + 2 
      z3 = z + 3 
      a = z1 / (q * rave) 
      a4 = 4 + a * a 
      at = atan (2 / a) 
      g1schul = a** ( - z1) - a4** ( - z1 / 2) * cos (z1 * at) + z2 *   &
      z1 * (a** ( - z3) + a4** ( - z3 / 2) * cos (z3 * at) ) - 2 * z1 * &
      a4** ( - z2 / 2) * sin (z2 * at)                                  
      RETURN 
      END FUNCTION g1schul                          
                                                                        
      FUNCTION g2schul (q, rave, z) 
!      ----------------------------                                     
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      z1 = z + 1 
      z2 = z + 2 
      a = z1 / (q * rave) 
      a1 = 1 + a * a 
      at = atan (1 / a) 
      g2schul = sin (z1 * at) - z1 / sqrt (a1) * cos (z2 * at) 
      RETURN 
      END FUNCTION g2schul                          
!**                                                                     
      FUNCTION pschj1 (q, drho, rave, z) 
!      ------------------------------                                   
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
      z1 = z + 1 
      z7 = z + 7 
      a = z1 / (q * rave) 
!      pschj1 = 8*(pi*drho)**2*(rave/z1)**6*a**z7*g1schul(q,rave,z)     
      arg1 = dlog ( (rave / z1) ) 
      arg2 = dlog (a) 
      arg3 = dlog (g1schul (q, rave, z) ) 
      arg4 = 6 * arg1 + z7 * arg2 + arg3 
      pschj1 = 8 * (pi * drho) **2 * dexp (arg4) 
                                                                        
      RETURN 
      END FUNCTION pschj1                           
!**                                                                     
      FUNCTION betaj1 (q, rave, z) 
!      ------------------------------                                   
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      z1 = z + 1 
      a = z1 / (q * rave) 
      g1 = g1schul (q, rave, z) 
      g2 = g2schul (q, rave, z) 
                                                                        
!      betaj1 = 2*a**z1*(1+a*a)**(-z1)*g2*g2/g1                         
      arg1 = dlog (a) 
      arg2 = dlog (1 + a * a) 
      arg3 = z1 * (arg1 - arg2) 
!      betaj1 = 2*a**z1*(1+a*a)**(-z1)*g2*g2/g1                         
      betaj1 = 2 * dexp (arg3) * g2 * g2 / g1 
      RETURN 
      END FUNCTION betaj1                           
!**                                                                     
      FUNCTION fj1qr (q, r, drho) 
!      ------------------------                                         
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      x = q * r 
      aj1 = (sin (x) - x * cos (x) ) / (x**3) 
      fj1qr = 4 * pi / 3 * r * r * r * drho * 3 * aj1 
      RETURN 
      END FUNCTION fj1qr                            
!**                                                                     
!** integrandfunktionen fuer numerische integration                     
!**                                                                     
      FUNCTION fiq1sch (r) 
!      -------------------                                              
!      f(q,r)*pschulz(r)                                                
                                                                        
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      COMMON / cschulz / zzc, ravec, drho, q 
                                                                        
      fiq1sch = fj1qr (q, r, drho) * pschulz (r, ravec, zzc) 
                                                                        
      RETURN 
      END FUNCTION fiq1sch                          
!**                                                                     
      FUNCTION fiq2sch (r) 
!      -------------------                                              
!      f(q,r)**2*pschulz(r)                                             
                                                                        
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (pi = 3.141592654) 
                                                                        
      COMMON / cschulz / zzc, ravec, drho, q 
                                                                        
      fiq2sch = fj1qr (q, r, drho) **2 * pschulz (r, ravec, zzc) 
                                                                        
      RETURN 
      END FUNCTION fiq2sch                          
