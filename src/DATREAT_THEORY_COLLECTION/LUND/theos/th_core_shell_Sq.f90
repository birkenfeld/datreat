      FUNCTION th39 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)  
!     ===================================================               
!                                                                       
!      shell3                                                           
!      ------                                                           
!                                                                       
! -------> core-shell-shell-modell <--------                            
!      auessere Schale kann parabolisch oder starlike                   
!      gemacht werden                                                   
!                                                                       
!                                                                       
      IMPLICIT none 
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL pa (20) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			                                                                        
      REAL(8) epsilon, errac, eps 
      INTEGER maxit 
      PARAMETER (epsilon = 1.0d-10, maxit = 100) 
                                                                        
      REAL(8) pi 
      PARAMETER (pi = 3.141592654d0) 
      REAL(8) Navog 
      PARAMETER (Navog = 6.022045d23) 
                                                                        
                                                                        
      REAL(8) adapint, rmina, rmaxa 
                                                                        
      REAL(8) amplitu, mcore, mbrush, rhocore, rhobrsh 
      REAL(8) bcore, bbrush, naggr, bkuhn, sigcore, sigbrush 
      REAL(8) rcore, rbrush, alpha, qa 
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION phi_01, r_end1, sigma_1 
      COMMON / cshell1 / phi_01, r_end1, sigma_1 
                                                                        
      DOUBLEPRECISION phi_02, r_end2, sigma_2 
      COMMON / cshell2 / phi_02, r_end2, sigma_2 
                                                                        
      DOUBLEPRECISION phi_03, r_end3, sigma_3 
      COMMON / cshell3 / phi_03, r_end3, sigma_3 
      DOUBLEPRECISION fstar, dparabol, xnustar 
      COMMON / chs3a / fstar, dparabol, xnustar 
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
                                              !!                        
      DOUBLEPRECISION q 
                           ! to communicate q !!                        
      COMMON / cqval / q 
                                              !!                        
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                              !!                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
      DOUBLEPRECISION V1, V2, V3, R1, R2, R3, a1, a2, a3 
      DOUBLEPRECISION swellc, swells2, swells3 
                                                                        
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION d1plus, part23, d2plus, d3plus, sigma1, sigma2,   &
      sigma3, partstar, gamma, lparabol, pi0_blob, xi_blob              
                                                                        
      REAL(8) Vc, Vb, Rc, Rm, Rm0 
      REAL(8) afactor, dbb, dbc, rhosolv 
      REAL(8) fc, fb, fbcont, fcs, fbs, fbrc 
      REAL(8) pi0, xi0, xi_est, phi, mu 
      REAL(8) reff, xi0_est, phi_est, blob3 
                                                                        
      REAL rhosolv4, phi_nr, x, th39 
      LOGICAL compute_phi 
                                                                        
      INTEGER  inew, ini, nparx, npar, ier, iphi 
                                                                        
      DOUBLEPRECISION Rc0, VcS, fshell2, z1t, Vs2, Rc2, Vs2S, fshell3 
      DOUBLEPRECISION Vs3, Rc3, Vs3S, f_brush, f_core 
                                                                        
      DOUBLEPRECISION r, erracc 
                                                                        
      DOUBLEPRECISION phi1, phi2, phi3, f1, f2, f3, delta 
      DOUBLEPRECISION vol1_ker, vol2_ker, vol3_ker 
      DOUBLEPRECISION sca1_ker, sca2_ker, sca3_ker 
                                                                        
      INTEGER ichange 
      SAVE V1, V2, V3 
                                                                        
                                                                        
      EXTERNAL vol1_ker, vol2_ker, vol3_ker 
      EXTERNAL sca1_ker, sca2_ker, sca3_ker 
                                                                        
                                                                        
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'shell3_f' 
         nparx = 20 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th39 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
                                  ! concentration factor                
         parnam (1) = 'amplitu' 
                                  ! molecular weight core   (g/mol)     
         parnam (2) = 'mcore' 
                                  ! molecular weight brush  (g/mol)     
         parnam (3) = 'mbrush' 
                                  ! density core matter     (g/cm**3)   
         parnam (4) = 'rhocore' 
                                  ! density brush matter    (g/cm**3)   
         parnam (5) = 'rhobrsh' 
                                  ! scattering length density core mat. 
         parnam (6) = 'bcore' 
                                  ! scattering length density brush mat.
         parnam (7) = 'bbrush' 
                                  ! aggregation number real to be fitted
         parnam (8) = 'naggr' 
                                  ! extra radius of shell1=core (compare
         parnam (9) = 'd1plus' 
                                  ! relative distribution of shell amoun
         parnam (10) = 'part23' 
                                  ! extra radius of shell2 (compared to 
         parnam (11) = 'd2plus' 
                                  ! extra radius of shell3 (compared to 
         parnam (12) = 'd3plus' 
                                  ! core   smearing                     
         parnam (13) = 'sigma1' 
                                  ! shell2 smearing                     
         parnam (14) = 'sigma2' 
                                  ! shell3 smearing                     
         parnam (15) = 'sigma3' 
                                  ! relative distribution of parbolic:st
         parnam (16) = 'partstar' 
                                  ! star-like exponent would be 4/3, con
         parnam (17) = 'gamma' 
                                  ! length=thickness of parabolic brush 
         parnam (18) = 'lparabol' 
                                  ! scattering length density correction
         parnam (19) = 'f_brush' 
                                  ! scattering length density correction
         parnam (20) = 'f_core' 
!                                                                       
         th39 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here length unit is cm -----                    
!                                                                       
      xunit = 1d-8 
                           ! q in cm**-1                                
      q = x / xunit 
                                                                        
      amplitu = pa (1) 
      mcore = abs (pa (2) ) 
      mbrush = abs (pa (3) ) 
      rhocore = pa (4) 
      rhobrsh = pa (5) 
      bcore = pa (6) 
      bbrush = pa (7) 
      naggr = abs (pa (8) ) 
      d1plus = abs (pa (9) ) * xunit 
      part23 = abs (pa (10) ) 
      d2plus = abs (pa (11) ) * xunit 
      d3plus = abs (pa (12) ) * xunit 
      sigma1 = abs (pa (13) ) * xunit 
      sigma2 = abs (pa (14) ) * xunit 
      sigma3 = abs (pa (15) ) * xunit 
      partstar = abs (pa (16) ) 
                                ! star-like exponent                    
      gamma = pa (17) 
      lparabol = abs (pa (18) ) * xunit 
      f_brush = pa (19) 
      f_core = pa (20) 
                                                                        
! hole Streulaengendichte des Loesungsmittels aus der Parameterliste    
      CALL getpar ('bsolv   ', rhosolv4, nopar ,params,napar,mbuf, ier) 
      rhosolv = rhosolv4 
                                                                        
!                                                                       
! falls vorhanden Umschalten auf Phi-Plot                               
!                                                                       
      CALL getpar ('phi_nr  ', phi_nr, nopar ,params,napar,mbuf, ier) 
      IF (ier.eq.0) then 
         compute_phi = .true. 
         iphi = NINT (phi_nr) 
      ELSE 
         compute_phi = .false. 
      ENDIF 
                                                                        
!                                                                       
! bestimme die Radien aus den anderen Parametern                        
!                                                                       
! 1. Kompaktvolumina von core und brush                                 
!                                                                       
      Vc = naggr * mcore / (rhocore * Navog) 
      Vb = naggr * mbrush / (rhobrsh * Navog) 
      CALL setpar ('vc      ', sngl (Vc * 1d24) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('vb      ', sngl (Vb * 1d24) , nopar ,params,napar,mbuf, ier) 
! 2. bestimme core-Kompaktradius                                        
      Rc0 = (3d0 * Vc / (4d0 * Pi) ) ** (0.33333333d0) 
! 3. damit den core-Schwellfaktor                                       
      R1 = Rc0 + d1plus 
      VcS = (4d0 * Pi / 3.0d0) * R1**3 
      swellc = VcS / Vc 
! 5. der relative Anteil der 2.ten Schale am Kompaktvolumen der Schale  
!!     fshell2 =  z1t(part23)+1d-5 ! to avoid zerodiv         !! z1t is 
                                          ! to avoid zerodiv  !! easier 
      fshell2 = part23 / (1 + part23) + 1d-5 
      Vs2 = Vb * fshell2 
      Rc2 = (3d0 * (VcS + Vs2) / (4d0 * Pi) ) ** (0.33333333d0) 
      R2 = Rc2 + d2plus 
      Vs2S = (4d0 * Pi / 3.0d0) * R2**3 
      swells2 = (Vs2S - VcS) / Vs2 
! 6. und die entsprechenden Daten fuer die 3te Schale                   
      fshell3 = 1.0d0 - fshell2 
      Vs3 = Vb * fshell3 
      Rc3 = (3d0 * (Vs2S + Vs3) / (4d0 * Pi) ) ** (0.33333333d0) 
      R3 = Rc3 + d3plus 
      Vs3S = (4d0 * Pi / 3.0d0) * R3**3 
      swells3 = (Vs3S - Vs2S) / Vs3 
!                                                                       
      CALL setpar ('r1      ', sngl (R1 / xunit) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('r2      ', sngl (R2 / xunit) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('r3      ', sngl (R3 / xunit) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('fshell2 ', sngl (fshell2) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('fshell3 ', sngl (fshell3) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('swell1  ', sngl (swellc) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('swell2  ', sngl (swells2) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('swell3  ', sngl (swells3) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('vs2     ', sngl (Vs2 * 1d24) , nopar ,params,napar,mbuf, ier) 
      CALL setpar ('vs3     ', sngl (Vs3 * 1d24) , nopar ,params,napar,mbuf, ier) 
                                                                        
! 7. Umrechnung des relativen Anteils star/parolic                      
      fstar = z1t (partstar) 
      CALL setpar ('fstar   ', sngl (fstar) , nopar ,params,napar,mbuf, ier) 
! 8. Uebertragung der Daten auf die Phi-common-Bloecke                  
      r_end1 = R1 
      r_end2 = R2 
      r_end3 = R3 
      sigma_1 = sigma1 
      sigma_2 = sigma2 
      sigma_3 = sigma3 
                                                                        
!      fstar      = fstar                                               
      dparabol = lparabol 
      xnustar = gamma 
                                                                        
      phi_01 = 1.0d0 / swellc 
      phi_02 = 1.0d0 / swells2 
      phi_03 = 1.0d0 / swells3 
                                                                        
      alpha = 1.0d0 
!                                                                       
! 9. Sicherstellen der Massenerhaltung                                  
!                                                                       
      eps = epsilon * 4 * Pi * (R2**3) / 3.0d0 
                                                                        
      IF (ichange (pa, nparx) .ne.0) then 
         WRITE (6, * ) 'new volumes' 
         Rm = r_end1 + 5 * sigma_1 
         V1 = adapint (vol1_ker, 0.0d0, Rm, eps, maxit, erracc) 
         Rm = r_end2 + 5 * sigma_2 
         V2 = adapint (vol2_ker, 0.0d0, Rm, eps, maxit, erracc) 
         Rm = r_end3 + 5 * sigma_3 
         V3 = adapint (vol3_ker, 0.0d0, Rm, eps, maxit, erracc) 
                                                                        
      ENDIF 
                                                                        
      CALL setpar ('v1i     ', sngl (V1 * 1d24) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('v2i     ', sngl (V2 * 1d24) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('v3i     ', sngl (V3 * 1d24) , nopar ,params,napar,mbuf, ier) 
                                                                        
! erstens: erzwinge Kompaktvolumenvehaeltnis 2te/3te Schale             
      IF (V2.gt.0.0d0) then 
         delta = part23 * V3 / V2 
      ELSE 
         delta = 0.0d0 
      ENDIF 
! zweitens: erzwinge Kompaktvolumen der gesamten Schale                 
      alpha = Vb / (V3 * (1 + part23) ) 
! reskaliere die Schalen, so dass die Volumina konserviert sind         
      a1 = Vc / V1 
      a2 = delta * alpha 
      a3 = alpha 
                                                                        
      CALL setpar ('alpha   ', sngl (alpha) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('delta   ', sngl (delta) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('vc_chk  ', sngl (a1 * V1 * 1d24) ,nopar ,params,napar,mbuf, ier) 
      CALL setpar ('vb_chk  ', sngl ( (a2 * V2 + a3 * V3)  * 1d24) , nopar ,params,napar,mbuf, ier)                       
                                                                        
                                                                        
      dbc = bcore-rhosolv * (1.0d0 - f_core) 
      dbb = bbrush - rhosolv * (1.0d0 - f_brush) 
                                                                        
                                                                        
                                         ! --> Streuung                 
      IF (.not.compute_phi) then 
!                                                                       
! ---------------- die Formfaktoren --------------------------------    
!                                                                       
         eps = epsilon * 4 * Pi * (R2**3) / 3.0d0 
                                                                        
         Rm = r_end1 + 5 * sigma_1 
         f1 = a1 * adapint (sca1_ker, 0.0d0, Rm, eps, maxit, erracc) 
         Rm = r_end2 + 5 * sigma_2 
         f2 = a2 * adapint (sca2_ker, 0.0d0, Rm, eps, maxit, erracc) 
         Rm = r_end3 + 5 * sigma_3 
         f3 = a3 * adapint (sca3_ker, 0.0d0, Rm, eps, maxit, erracc) 
                                                                        
                                                                        
! und Kontraste anbrigen                                                
                                                                        
         afactor = amplitu / (Vc + Vb) 
                                                                        
         th39 = afactor * (dbc * f1 + dbb * (f2 + f3) ) **2 
                                                                        
!! add the blob-scattering  (xi=xi0*phi**(1/(mu-2)))                    
!       mu=0.6d0                                                        
!       phi=1.0d0                                                       
!                                                                       
!       th39 = th39 +                                                   
!     *        amplitu * (bbrush-rhosolv)**2*                           
!     *        blob3(q,xi_blob,phi,Pi0_blob,mu)                         
                                                                        
                                                                        
                                         ! --> ansonsten Phi(r)         
      ELSE 
         r = x * xunit 
                                         ! --> Phi Schale 1             
         IF (iphi.eq.1) then 
            th39 = a1 * phi1 (r) 
                                         ! --> Phi Schale 2             
         ELSEIF (iphi.eq.2) then 
            th39 = a2 * phi2 (r) 
                                         ! --> Phi Schale 3             
         ELSEIF (iphi.eq.3) then 
            th39 = a3 * phi3 (r) 
                                         ! --> Phi Schale 1+2+3         
         ELSEIF (iphi.eq.4) then 
            th39 = a1 * phi1 (r) + a2 * phi2 (r) + a3 * phi3 (r) 
                                         ! --> radiale Massenverteilung 
         ELSEIF (iphi.eq.5) then 
            th39 = (a1 * phi1 (r) + a2 * phi2 (r) + a3 * phi3 (r) )     &
            * ( (r / xunit) **2)                                        
                                         ! --> Streulaengendichteprofil 
         ELSEIF (iphi.eq.6) then 
            th39 = a1 * phi1 (r) * dbc + (a2 * phi2 (r) + a3 * phi3 (r) &
            ) * dbb                                                     
         ELSE 
            th39 = 0 
         ENDIF 
                                                                        
      ENDIF 
                                                                        
      RETURN 
      END FUNCTION th39                             
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function fermi (x, sigma) 
!      ----------------------------------------                         
!                                                                       
! ---> fermi Funkktion 1/(1+exp(x/sigma)  ROBUST                        
!                                                                       
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION x, sigma, arg 
!                                                                       
      IF (sigma.gt.0.0d0) then 
         arg = x / sigma 
         IF (arg.lt. - 100d0) then 
            fermi = 1.0d0 
            RETURN 
         ENDIF 
         IF (arg.gt.100d0) then 
            fermi = 0.0d0 
            RETURN 
         ENDIF 
         fermi = 1.0d0 / (1.0d0 + exp (arg) ) 
      ELSE 
         fermi = 0.0d0 
         RETURN 
      ENDIF 
                                                                        
      RETURN 
      END FUNCTION fermi                            
                                                                        
                                                                        
      DOUBLEPRECISION function z1t (x) 
!      --------------------------------                                 
!                                                                       
! ---> Transformation, die beliebige x auf [0,1] abbildet               
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION Pi, x 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      z1t = 2.0d0 * atan (abs (x) ) / pi 
                                                                        
      RETURN 
      END FUNCTION z1t                              
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!  Modellierung 1.Schale : hier einfache Fermi-Fktn               !!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function phi1 (r) 
!     ---------------------------------                                 
!                                                                       
! ---> Volumefraction 1.Schale                                          
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, fermi 
      DOUBLEPRECISION phi_01, r_end1, sigma_1 
      COMMON / cshell1 / phi_01, r_end1, sigma_1 
                                                                        
      phi1 = phi_01 * fermi (r - r_end1, sigma_1) 
                                                                        
      RETURN 
      END FUNCTION phi1                             
                                                                        
                                                                        
      DOUBLEPRECISION function vol1_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Volumenintegration 1.Schale                           
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi1, Pi 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      vol1_ker = (4 * Pi) * phi1 (r) * r * r 
                                                                        
      RETURN 
      END FUNCTION vol1_ker                         
                                                                        
                                                                        
      DOUBLEPRECISION function sca1_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Streuintegration 1.Schale                             
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi1, Pi, qq 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      DOUBLEPRECISION q 
      COMMON / cqval / q 
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                                                        
                                                                        
      IF (q.lt.1d-7 / xunit) then 
         qq = 1d-7 / xunit 
      ELSE 
         qq = q 
      ENDIF 
                                                                        
      sca1_ker = (4 * Pi) * phi1 (r) * sin (qq * r) * r / qq 
                                                                        
      RETURN 
      END FUNCTION sca1_ker                         
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!  Modellierung 2.Schale : hier einfache Fermi-Fktn               !!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function phi2 (r) 
!     ---------------------------------                                 
!                                                                       
! ---> Volumefraction 2.Schale                                          
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, fermi 
      DOUBLEPRECISION phi_01, r_end1, sigma_1 
      COMMON / cshell1 / phi_01, r_end1, sigma_1 
                                                                        
      DOUBLEPRECISION phi_02, r_end2, sigma_2 
      COMMON / cshell2 / phi_02, r_end2, sigma_2 
                                                                        
      phi2 = phi_02 * (1 - fermi (r - r_end1, sigma_1) ) * fermi (r -   &
      r_end2, sigma_2)                                                  
                                                                        
      RETURN 
      END FUNCTION phi2                             
                                                                        
                                                                        
      DOUBLEPRECISION function vol2_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Volumenintegration 2.Schale                           
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi2, Pi 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      vol2_ker = (4 * Pi) * phi2 (r) * r * r 
                                                                        
      RETURN 
      END FUNCTION vol2_ker                         
                                                                        
                                                                        
      DOUBLEPRECISION function sca2_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Streuintegration 2.Schale                             
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi2, Pi, qq 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      DOUBLEPRECISION q 
      COMMON / cqval / q 
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                                                        
                                                                        
      IF (q.lt.1d-7 / xunit) then 
         qq = 1d-7 / xunit 
      ELSE 
         qq = q 
      ENDIF 
                                                                        
      sca2_ker = (4 * Pi) * phi2 (r) * sin (qq * r) * r / qq 
                                                                        
      RETURN 
      END FUNCTION sca2_ker                         
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!  Modellierung 3.Schale :  Fermi-Fktn * Parabol oder star        !!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                                                                        
                                                                        
                                                                        
      DOUBLEPRECISION function phi3 (r) 
!     ---------------------------------                                 
!                                                                       
! ---> Volumefraction 3.Schale                                          
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, fermi, parabol_brush 
      DOUBLEPRECISION phi_01, r_end1, sigma_1 
      COMMON / cshell1 / phi_01, r_end1, sigma_1 
                                                                        
      DOUBLEPRECISION phi_02, r_end2, sigma_2 
      COMMON / cshell2 / phi_02, r_end2, sigma_2 
                                                                        
      DOUBLEPRECISION phi_03, r_end3, sigma_3 
      COMMON / cshell3 / phi_03, r_end3, sigma_3 
      DOUBLEPRECISION fstar, dparabol, xnustar 
      COMMON / chs3a / fstar, dparabol, xnustar 
                                                                        
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                                                        
                                                                        
      IF (r.gt.r_end2) then 
                                                                        
         phi3 = phi_03 * (1 - fermi (r - r_end2, sigma_2) ) * fermi (r -&
         r_end3, sigma_3) * ( (1.d0 - fstar) * parabol_brush (r -       &
         r_end2, dparabol) + fstar / ( (r / xunit + 1d-4) **xnustar) )  
                                                                        
      ELSE 
         phi3 = phi_03 * (1 - fermi (r - r_end2, sigma_2) ) * fermi (r -&
         r_end3, sigma_3) * ( (1.d0 - fstar) * parabol_brush (r_end2 -  &
         r_end2, dparabol) + fstar / ( (r_end2 / xunit + 1d-4) **       &
         xnustar) )                                                     
      ENDIF 
                                                                        
      RETURN 
      END FUNCTION phi3                             
                                                                        
      DOUBLEPRECISION function parabol_brush (x, l) 
!     --------------------------------------------                      
                                                                        
      IMPLICIT none 
      DOUBLEPRECISION x, l, yl 
                                                                        
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                                                        
      IF (l.lt.1d-4 * xunit) then 
         yl = 1d-4 * xunit 
      ELSE 
         yl = l 
      ENDIF 
                                                                        
      IF (x.lt.0.0d0) then 
         parabol_brush = 1.0d0 
         RETURN 
      ENDIF 
                                                                        
      IF (x.gt.yl) then 
         parabol_brush = 0.0d0 
         RETURN 
      ENDIF 
                                                                        
      parabol_brush = 1.0d0 - (x / yl) **2 
                                                                        
      RETURN 
      END FUNCTION parabol_brush                    
                                                                        
                                                                        
      DOUBLEPRECISION function vol3_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Volumenintegration 3.Schale                           
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi3, Pi 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      vol3_ker = (4 * Pi) * phi3 (r) * r * r 
                                                                        
      RETURN 
      END FUNCTION vol3_ker                         
                                                                        
                                                                        
      DOUBLEPRECISION function sca3_ker (r) 
!     -------------------------------------                             
!                                                                       
! --> Kernel fuer Streuintegration 2.Schale                             
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION r, phi3, Pi, qq 
      PARAMETER (Pi = 3.141592654d0) 
                                                                        
      DOUBLEPRECISION q 
      COMMON / cqval / q 
                                              !!                        
      DOUBLEPRECISION xunit 
                                              !!                        
      COMMON / c_unit1 / xunit 
                                                                        
      IF (q.lt.1d-7 / xunit) then 
         qq = 1d-7 / xunit 
      ELSE 
         qq = q 
      ENDIF 
                                                                        
      sca3_ker = (4 * Pi) * phi3 (r) * sin (qq * r) * r / qq 
                                                                        
      RETURN 
      END FUNCTION sca3_ker                         
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      REAL(8) function blob3 (q, xi0, phi, Pi0, mu) 
!      =====================================                            
!                                                                       
! Blobformel nach Dozier mit konsistenter Skalierung bzgl.              
! Scaling laws (Daoud et al, Macromolecules 8, p804, 1975               
!                                                                       
! so dass I(0) --> Kompressiblitaetsstreuung (nach Anbringen des Kontras
! alles auf Volumenbruch Phi des Polymers bezogen.                      
!                                                                       
! siehe auch Maple Worksheet ~/SEMIDILUTE_SOLUTIONS/Schneiders_blob2.mws
!                                                                       
                                                                        
                                                                        
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION q, xi0, phi, Pi0, mu, X 
                                                                        
! q    --> Q                                                            
! xi0  --> referenz Blob-Groesse sodass: xi = xi0* phi**(-nu/3*nu-1), mi
! phi  --> Volumenanteil                                                
! Pi0  --> OsmoticPressure/kT ("phi=1") --> dient als eff. Intesitaetssk
! mue  --> Exponent mue = 1/nu-1    (Flory excl. Vol. mue = 2/3)        
                                                                        
      X = q * xi0 
                                                                        
      blob3 = - Phi** (2 * (mu - 1) / (mu - 2) ) * (mu - 2) * sin (mu * &
      atan (q * xi0 * Phi** (1 / (mu - 2) ) ) ) / Pi0 / q / xi0 *       &
      (xi0**2 * Phi** (2 / (mu - 2) ) * q**2 + 1) ** ( - mu / 2)        &
      / mu / 3                                                          
                                                                        
                                                                        
!      write(6,'(A,6E13.6)')'blob3: ',q,xi0,phi,Pi0,mu,blob3            
                                                                        
      RETURN 
      END FUNCTION blob3                            
                                                                        
                                                                        
                                                                        
      INTEGER function ichange (p, n) 
!      -----------------------------                                    
      INTEGER n, j, i, max 
      PARAMETER (max = 200) 
      REAL p (n), p0 (max) 
      SAVE p0 
                                                                        
      IF (n.gt.max) then 
         ichange = max 
         RETURN 
      ENDIF 
                                                                        
      j = 0 
      DO i = 1, n 
      IF (p (i) .ne.p0 (i) ) j = j + 1 
      p0 (i) = p (i) 
      enddo 
                                                                        
      ichange = j 
      RETURN 
      END FUNCTION ichange                          
