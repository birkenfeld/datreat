      SUBROUTINE qcalc (q, qz_qx, chi, qz, qortho) 
!    ------------------------------------------                         
!  -Berechnet qz und qortho fuer schraege schnitte aus 2d daten         
!  -fuer  caille theorie                                                
                                                                        
!--  Berechnung der schraegen Schnitte--------------                    
      xmx = sqrt (1.0 + qz_qx**2) 
      qz0 = q / xmx 
      qortho0 = qz0 * qz_qx 
!--  Kippung des Probenhalters entspricht Streuebene-------------------- 
      chirad = chi * 3.1415926535 / 180.0 
      qz = cos (chirad) * qz0 
      qortho = sqrt (qortho0**2 + (sin (chirad) * qz0) **2) 
                                                                        
      RETURN 
      END SUBROUTINE qcalc 