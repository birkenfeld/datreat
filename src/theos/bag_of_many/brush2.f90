!*ds                                                                    
!*ed                                                                    
      SUBROUTINE thbru2  (q, tau, r0, sigmar0, rho_innen, ratio, allfak, &
      k_max, ll_max, iflcorr, ninter, nqres, dqres, aiqt, aiqt0, sq,    &
      ier)                                                              
!                                                                       
! **** version mit wellenlaengenverbreiterung!                          
! **** input parameter:                                                 
!      q ......... skalierter impulsuebertrag (--> schichtdicke=1)      
!      tau ....... skalierte zeit                                       
!      r0 ........ skalierte dicke der festen innenkugel                
!      sigmar0 ... verteilungsbreite (gauss) der innenkugel             
!      rho_innen . relative streulaengendichte der innenkugel           
!      ratio ..... verhaeltnis s(q) zu inelastish                       
!      allfak .... elam(n) --> elam(n)+allfak*l^2                       
!      k_max ..... hoechster eigenwert der radialen loesung             
!      ll_max .... hoechstes beruecksichtigtes l                        
!      iflcorr ... flag: benutzung von korrekturen fuer l<>0            
!      ninter .... anzahl stuetzstellen fuer eigenfunktionen            
!      nqres ..... anzahl stuetzstellen der wellenlaengenmittelung      
!      dqres ..... fwhm (relativ!) der wellenlaengenverteilung          
! **** output parameter:                                                
!      aiqt ...... i(q,tau)                                             
!      aiqt0...... i(q,tau=0)                                           
!      sq ........ i(q,tau=+infinity)                                   
!      ier ....... fehlerindikator                                      
!==================================================================     
!                                                                       
!                                                                       
      IMPLICIT real (8)(a - h, o - z) 
      CHARACTER(8) filna, chrnxt, chrval 
      LOGICAL found, newsol, newsoq, newwoq, newfl, newcof 
                                                                        
      PARAMETER (pi = 3.141592654d0) 
                                                                        
      PARAMETER (mint = 10) 
      PARAMETER (mint2 = ( (mint + 1) * mint) / 2) 
      DIMENSION cint (mint2), aint (mint + 1), bint (mint + 1) 
                                                                        
      PARAMETER (nmax = 100, kdom = 10, llmax = 15, kmax = 10 + kdom) 
      DIMENSION x (nmax, kmax), y (nmax, kmax), elams ( - 2:kdom),      &
      n (kmax)                                                          
      COMMON / c2report / x, y, elams, n, i, ireport 
                                                                        
      DIMENSION ipermu (nmax) 
                                                                        
                                                                        
! ******************************************************************    
! * speicher fuer die wellenlaengenmittelung                       *    
! ******************************************************************    
      PARAMETER (mqres = 20) 
      COMMON / c2qres / qarray (mqres), warray (mqres) 
                                                                        
! ******************************************************************    
! * speicher fuer korrekturfaktoren fuer hoehere l-werte           *    
! ******************************************************************    
      DIMENSION fl (0:llmax, 0:kmax) 
      DIMENSION fe (0:llmax, 0:kmax) 
      COMMON / c2flfe / fl, fe 
                                                                        
! ******************************************************************    
! * save-commons                                                   *    
! ******************************************************************    
      COMMON / i2save / kmold1, ninterold, iflcold, kmold3, llmold3,    &
      nqrold                                                            
      COMMON / r2save / r0old1, qold2, r0old2, sigmaold, rhoold, soq,   &
      qold3, told4, aiqti, aiqt0i, qold10, dqold                        
      COMMON / a2save / afakv (0:llmax, 0:kmax), tau_i (0:llmax, 0:kmax) 
                                                                        
! ******************************************************************    
! * datenfelder fuer nag-sturm-lioville dgl-solver                 *    
! ******************************************************************    
      PARAMETER (mbrkp = 20) 
      DIMENSION xpoint (mbrkp) 
      DIMENSION hmax (2, mbrkp) 
! ******************************************************************    
! * common kommunikation mit externen funktionen des dgl-solvers   *    
! ******************************************************************    
      COMMON / c2oefpa / eta 
      COMMON / c2monit / infomo 
! ******************************************************************    
! * external functions fuer             dgl-solver                 *    
! ******************************************************************    
      EXTERNAL c2oeffn 
      EXTERNAL b2dyval 
      EXTERNAL m2onit 
      EXTERNAL r2eport 
                                                                        
                                                                        
! ******************************************************************    
! * sektion 1 initialisieren der dgl-loesungen                     *    
! ******************************************************************    
!                                                                       
! ueberspringen dieser sektion falls die parameter nicht                
! geaendert wurden.                                                     
!                                                                       
      newsol = .false. 
                                                                        
      IF (k_max.gt.kmold1.or.ninterp.ne.ninterold.or.r0.ne.r0old1) then 
                                                                        
      WRITE (6,  * ) 'k_max, km_old     ', k_max, kmold1 
      WRITE (6,  * ) 'ninterp .....     ', ninterp, ninterold 
      WRITE (6,  * ) 'r0 ..........     ', r0, r0old1 
                                                                        
         kmold1 = k_max 
         ninterold = ninterp 
         r0old1 = r0 
         newsol = .true. 
                                                                        
! ******************************************************************    
! * initialisationen   fuer             dgl-solver                 *    
! ******************************************************************    
         tol = 1.0d-4 
         elams ( - 2) = 10.0d0 
         elams ( - 1) = 18.0d0 
                                                                        
         DO i = 1, 2 
         DO j = 1, mbrkp 
         hmax (i, j) = 0.0d0 
         enddo 
         enddo 
                                                                        
! --- setup xpoint ---                                                  
         xpoint (1) = r0 
         xpoint (2) = r0 
         xpoint (3) = r0 + 0.10d0 
         xpoint (4) = r0 + 0.50d0 
         xpoint (5) = r0 + 0.90d0 
         xpoint (6) = r0 + 0.95d0 
         xpoint (7) = r0 + 0.99d0 
         xpoint (8) = r0 + 1.00d0 
         xpoint (9) = r0 + 1.00d0 
         nbrkp = 9 
         match = 4 
                                                                        
         maxit = 0 
         maxfun = 0 
         infomo = iout () 
         ireport = 0 
         eta = 1.0d0 
         nint = 4 
         ninterp = ninter 
         IF (ninterp.le.10) then 
            ninterp = 50 
            WRITE (6, * ) 'thbrush: set ninter to ', ninterp 
         ENDIF 
         c_exp = - 4.0d0 / 3.0d0 
! --- check der dimensionen ---                                         
         IF (k_max.gt.kdom) then 
            WRITE (6, * ) 'thbrush: k_max is out of range:', k_max, '>',&
            kdom                                                        
            ier = 1 
            RETURN 
         ENDIF 
         IF (ll_max.gt.llmax) then 
      WRITE (6,  * ) 'thbrush: ll_max is out of range:', ll_max, '>', ll&
     &max                                                               
            ier = 2 
            RETURN 
         ENDIF 
! ******************************************************************    
                                                                        
                                                                        
!      ===============================================                  
!      =  brush-dgl fuer radiale geometrie (star)    =                  
!      =  grad e(r) div u                            =                  
!      =  u = loesung / r**2 --> report              =                  
!      ===============================================                  
                                                                        
         DO kelam = 0, k_max 
                                                                        
         n (1) = 0 
         delam = elams (kelam - 1) - elams (kelam - 2) 
         elam = elams (kelam - 1) + delam * 2.0d0 
         CALL d02kef (xpoint, nbrkp, match, c2oeffn, b2dyval, kelam,    &
         tol, elam, delam, hmax, maxit, maxfun, m2onit, r2eport, ifail) 
         IF (ifail.ne.0) write (6, * ) 'd02kef-->ifail = ', ifail 
         IF (iout () .ge.0) then 
            WRITE (6, 6003) kelam, elam, delam 
 6003 FORMAT   (' lambda(',i3,') = ',f14.8,' +- ',f14.8,' ..!') 
         ENDIF 
!                                                                       
!cc      write(6,*)'..... sorting (1)                -->(2)'            
         n (2) = n (1) 
         ipath = 2 
         DO i = 1, n (1) 
         ipermu (i) = i 
         enddo 
         CALL dsvrgp (n (1), x (1, 1), x (1, 2), ipermu) 
         CALL dpermu (n (1), y (1, 1), ipermu, ipath, y (1, 2) ) 
!                                                                       
!        write(6,*)'..... remove double x-values (2) -->(3)'            
         nn = 1 
         x (1, 3) = 0.0d0 
         y (1, 3) = 0.0d0 
         DO i = 1, n (2) 
         IF (x (nn, 3) .ne.x (i, 2) ) then 
            nn = nn + 1 
            x (nn, 3) = x (i, 2) 
            y (nn, 3) = y (i, 2) 
         ENDIF 
         enddo 
         n (3) = nn 
!                                                                       
!        write(6,*)'..... interpolate (3)            -->(4)'            
!        write(6,*)'..... (4) * c(z)                 -->(5)'            
         IF (nint.gt.mint) nint = mint 
         IF (ninterp.gt.nmax) ninterp = nmax 
         nint1 = nint + 1 
         nint2 = (nint * nint1) / 2 
         DO i = 1, nmax 
         y (i, 4) = 0.d0 
!cc       ys(i,4)  = 0.d0                                               
         enddo 
         n (4) = ninterp 
         dx = 1.0d0 / (ninterp - 1) 
         is = 1 
!cc      xinteg = 0.0d0                                                 
         DO i = 1, ninterp 
         xh = dx * (i - 1) + r0 
         x (i, 4) = xh 
         x (i, 5) = xh 
!cc        write(6,*)i,xh                                               
         IF (xh.lt.x (1, 3) ) then 
            j = 1 
            GOTO 1010 
         ENDIF 
         DO j = is, n (3) - 1 
         IF ( (x (j, 3) .le.xh) .and. (x (j + 1, 3) .gt.xh) ) then 
!cc            write(6,'(1x,2i3,3f16.6)')j,is,x(j,3),xh,x(j+1,3)        
            GOTO 1010 
         ENDIF 
         enddo 
 1010    CONTINUE 
         is = j 
         i0 = is - nint / 2 
         IF (i0.lt.1) i0 = 1 
         IF (i0.gt.n (3) - nint - 1) i0 = n (3) - nint - 1 
         DO j = 0, nint 
         aint (j + 1) = x (i0 + j, 3) 
         bint (j + 1) = y (i0 + j, 3) 
         enddo 
         CALL e01aaf (aint, bint, cint, nint1, nint2, nint, xh) 
         y (i, 4) = cint (nint2) 
         IF (xh.gt.0.d0) then 
            y (i, 5) = y (i, 4) / (xh** (4.0d0 / 3.0d0) ) 
! --- multipliziere mit fermi-abschneidefunktion -----                  
!cc           arg    = (rfermi - xh) / sfermi                           
!cc           y(i,5) = y(i,5) / ( 1.0d0 + dexp(-arg))                   
!cc           xinteg = xinteg + y(i,4)**2                               
         ELSE 
            y (i, 5) = 0.0d0 
         ENDIF 
         enddo 
         n (5) = n (4) 
!cc      xinteg = xinteg * dx                                           
!cc      write(6,*)'....> integral( u(z)  **2 dz ) = ',xinteg           
!cc     write(6,*)'===================================================='
!                                                                       
! --- speichern der werte -----------                                   
         elams (kelam) = elam 
         kelam_u = 10 + kelam 
         n (kelam_u) = n (5) 
         DO i = 1, n (5) 
         x (i, kelam_u) = x (i, 5) 
         y (i, kelam_u) = y (i, 5) 
         enddo 
         enddo 
      ENDIF 
!    -----<--- init-if <---                                             
                                                                        
                                                                        
                                                                        
! ******************************************************************    
! * berechnung neuer q-gewichte falls erforderlich                 *    
! ******************************************************************    
      IF (nqres.gt.mqres) nqres = mqres 
      newwoq = .false. 
                                                                        
      IF (q.ne.qold10.or.dqres.ne.dqold.or.nqres.ne.nqrold) then 
                                                                        
         newwoq = .true. 
         qold10 = q 
         dqold = dqres 
         nqrold = nqres 
                                                                        
         dqabs = q * dqres 
         dqh = 2 * dqabs / (nqres + 1) 
                                                                        
         wnorm = 0.0d0 
         DO i = 1, nqres 
         qact = q - dqabs + i * dqh 
         IF (qact.le.q) then 
            wact = i * dqh 
         ELSE 
            wact = dqabs + dqabs - i * dqh 
         ENDIF 
         wnorm = wnorm + wact 
         qarray (i) = qact 
         warray (i) = wact 
         enddo 
         DO i = 1, nqres 
         warray (i) = warray (i) / wnorm 
         enddo 
                                                                        
         IF (iout () .gt.0) then 
            WRITE (6, * ) 'dqres=', dqres 
            WRITE (6, * ) 'dqabs=', dqabs 
      WRITE (6,  * ) 'dqh  =', dqh 
            DO i = 1, nqres 
      WRITE (6,  * ) i, ': q=', qarray (i) , '  w=', warray (i) 
            enddo 
         ENDIF 
                                                                        
      ENDIF 
                                                                        
! ******************************************************************    
! * berechnung von s(q) falls erforderlich                         *    
! ******************************************************************    
      newsoq = .false. 
                                                                        
      IF (newwoq.or.q.ne.qold2.or.r0.ne.r0old2.or.sigmar0.ne.sigmaold.or&
     &.rho_innen.ne.rhoold) then                                        
                                                                        
         qold2 = q 
         r0old2 = r0 
         sigmaold = sigmar0 
         rhoold = rho_innen 
         newsoq = .true. 
                                                                        
         c_exp = - 4.0d0 / 3.0d0 
         soq = 0.d0 
         DO i = 1, nqres 
         soq = soq + warray (i) * s2_of_q (qarray (i), r0, sigmar0,     &
         rho_innen, c_exp)                                              
         enddo 
                                                                        
         IF (iout () .ge.0) then 
            WRITE (6, * ) 'thbrush new s(q=', q, ')=', soq 
         ENDIF 
                                                                        
      ENDIF 
                                                                        
                                                                        
                                                                        
! ******************************************************************    
! * berechnung der koeffizienten fuer i(q,t) inelastisch           *    
! ******************************************************************    
      newfl = .false. 
! --- a. muss fl auf 1 initialisiert werden ?? ---                      
      IF (fl (0, 0) .ne.1.0d0.and.iflcorr.eq.0) then 
         iflcold = iflcorr 
         newfl = .true. 
         WRITE (6, * ) 'init fl, fe to 1 !' 
         DO i = 0, kmax 
         DO j = 0, llmax 
         fl (i, j) = 1.0d0 
         fe (i, j) = 1.0d0 
         enddo 
         enddo 
      ENDIF 
! --- b. muss fl eingelsesn werden? ---------                           
      IF (iflcorr.ne.iflcold) then 
         iflcold = iflcorr 
         newfl = .true. 
                                                                        
         DO i = 0, kmax 
         DO j = 0, llmax 
         fl (i, j) = 1.0d0 
         fe (i, j) = 1.0d0 
         enddo 
         enddo 
                                                                        
         filna = 'fldat' 
         OPEN (66, file = filna) 
         READ (66, * ) kmi 
         lli = 0 
 6601    CONTINUE 
         READ (66, *, end = 6602) (fl (lli, j), j = 0, kmi) 
         READ (66, *, end = 6602) (fe (lli, j), j = 0, kmi) 
         lli = lli + 1 
         IF (lli.gt.llmax) goto 6602 
         GOTO 6601 
 6602    CONTINUE 
         CLOSE (66) 
!c       write(6,*)'kmi=',kmi,'  lli=',lli,'  llmax=',llmax             
         DO j = 0, kmi 
         DO l = lli, llmax 
         fl (l, j) = fl (lli - 1, j) 
         fe (l, j) = fe (lli - 1, j) 
         enddo 
         enddo 
 6604    CONTINUE 
         DO l = 0, llmax 
         WRITE (6, 6603) l, (fl (l, j), j = 0, 5), (fe (l, j), j = 0, 5) 
 6603 FORMAT    (' l=',i2,1x,6f7.2,' ....'/                             &
     &           '   ',2x,1x,6f7.2)                                     
         enddo 
      ENDIF 
!      ------ of section b. of fl-inits ---                             
                                                                        
! --------> so jetzt kommen die exp(-tau/tau_i) koeffizientenarrays --- 
      newcof = .false. 
                                                                        
      IF (k_max.ne.kmold3.or.ll_max.ne.llmold3.or.newsol.or.newwoq.or.ne&
     &wfl.or.q.ne.qold3) then                                           
                                                                        
         kmold3 = k_max 
         llmold3 = ll_max 
         qold3 = q 
         newcof = .true. 
                                                                        
         DO j = 0, k_max 
         DO l = 0, ll_max 
         fcorrl = fl (0, j) / fl (l, j) 
         fcorre = fe (0, j) / fe (l, j) 
         dx = x (2, 10 + j) - x (1, 10 + j) 
         elam = elams (j) + allfak * l * l 
         afak0 = elams (0) / elam 
         tau_i (l, j) = elams (0) * fcorrl / elam 
         afakv (l, j) = 0.0d0 
         DO i = 1, nqres 
         ffq = b2sj_int (qarray (i), l, dx, y (1, 10 + j), n (10 + j),  &
         r0)                                                            
         sffq = (2 * l + 1) * ffq * ffq 
         afak = afak0 * sffq * 4.0d0 * pi * fcorre * fcorrl 
         afakv (l, j) = afakv (l, j) + warray (i) * afak 
         enddo 
         enddo 
         enddo 
         IF (iout () .ge.0) then 
            WRITE (6, * ) 'thbrush new afak and tau- coeffs' 
         ENDIF 
                                                                        
      ENDIF 
!      -----> if init sektion                                           
                                                                        
! ********************************************************************* 
!  berechnungssektion                                                   
! ********************************************************************* 
      IF (newcof.or.tau.ne.told4) then 
         told4 = tau 
         aiqti = 0.0d0 
         aiqt0i = 0.0d0 
         DO j = 0, k_max 
         DO l = 0, ll_max 
         etau = dexp ( - tau / tau_i (l, j) ) 
         aiqti = aiqti + afakv (l, j) * etau 
         aiqt0i = aiqt0i + afakv (l, j) 
         enddo 
         enddo 
         IF (iout () .gt.0) then 
            WRITE (6, * ) 'thbrush(t=', tau, ')-->', aiqti, aiqt0i 
         ENDIF 
      ENDIF 
!                                                                       
! -------------- final results ---------------------------              
      ier = 0 
      sq = soq 
      aiqt = aiqti + ratio * soq 
      aiqt0 = aiqt0i + ratio * soq 
                                                                        
      RETURN 
      END SUBROUTINE thbru2                         
                                                                        
                                                                        
                                                                        
! ******************************************************************    
! * external sbrs    nag-sturm-lioville dgl-solver                 *    
! ******************************************************************    
      SUBROUTINE c2oeffn (p, q, dqdl, x, elam, jint) 
!     =======================================                           
      IMPLICIT real (8)(a - h, o - z) 
      COMMON / c2oefpa / eta 
                                                                        
! ----- die funktion q(x,lambda) ----                                   
      qq = eta / (x**4) 
      q = elam * qq 
! ----- die ableitung von dq/dlambda ----                               
      dqdl = qq 
! ----- die funktion p(x) -----                                         
      ENTRY pfp2 (x, p) 
      p = 1.d0 / (x**5) 
                                                                        
      RETURN 
      END SUBROUTINE c2oeffn                        
!                                                                       
!                                                                       
      SUBROUTINE b2dyval (xl, xr, elam, yl, yr) 
!     =================================== boundary values               
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION yl (3), yr (3) 
!                                                                       
      yl (1) = 0.d0 
      yl (2) = 1.0d0 
!                                                                       
      yr (1) = 1.0d0 
      yr (2) = 0.0d0 
                                                                        
      RETURN 
      END SUBROUTINE b2dyval                        
!                                                                       
!                                                                       
      SUBROUTINE m2onit (maxit, iflag, elam, finfo) 
!     ========================================                          
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION finfo (15) 
      COMMON / c2monit / infomo 
                                                                        
      IF (infomo.le.0) return 
                                                                        
      IF (iflag.lt.0) then 
         WRITE (6, * ) 'monit iflag=', iflag, ' error exit will follow' 
         GOTO 900 
      ENDIF 
                                                                        
      IF (iflag.eq.1) write (6,  * ) 'monit iflag=', iflag, ' trying to &
     &bracket lambda =', elam                                           
                                                                        
      IF (iflag.eq.2) write (6,  * ) 'monit iflag=', iflag, ' converging&
     & to lambda =', elam                                               
                                                                        
      IF (infomo.le.1) return 
  900 CONTINUE 
!     finfo dump                                                        
      DO i = 1, 15 
      WRITE (6, * ) 'monit finfo(', i, ') = ', finfo (i) 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE m2onit                         
                                                                        
      SUBROUTINE r2eport (xp, vp, jint) 
!     ===============================                                   
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION vp (3) 
                                                                        
      PARAMETER (nmax = 100, kdom = 10, llmax = 15, kmax = 10 + kdom) 
      DIMENSION x (nmax, kmax), y (nmax, kmax), elams ( - 2:kdom),      &
      n (kmax)                                                          
      COMMON / c2report / x, y, elams, n, i, ireport 
                                                                        
      yp = dexp (vp (3) * 0.5d0) * dsin (vp (2) * 0.5d0) / dsqrt (vp (1)&
      )                                                                 
      CALL pfp2 (xp, px) 
      yd = dsqrt (vp (1) ) * dexp (vp (3) * 0.5d0) * dcos (vp (2)       &
      * 0.5d0) / px                                                     
                                                                        
      IF (ireport.gt.0) write (6, 6001) n (1), xp, yp, yd 
 6001 FORMAT(' report(',i3,'): ',3f16.7) 
                                                                        
      IF (n (1) .ge.nmax) then 
         WRITE (6, * ) 'report no of x points is exhausted max=', nmax 
         RETURN 
      ENDIF 
                                                                        
      n (1) = n (1) + 1 
      x (n (1), 1) = xp 
      y (n (1), 1) = yp / (xp**2) 
!ccc  ys(n(1),1) = yd                                                   
                                                                        
      RETURN 
      END SUBROUTINE r2eport                        
                                                                        
!*ds                                                                    
!*ed                                                                    
!                                                                       
! ----- spherical besselfunction of order m>=0 -------------------------
!                                                                       
      FUNCTION b2sjn (m, z) 
!      ==================                                               
!                                                                       
      IMPLICIT real (8)(a - h, o - z) 
!                                                                       
      DATA eps / 1d-6 / 
!                                                                       
! --- look for bad m value ---                                          
      IF (m.lt.0) then 
      WRITE (6,  * ) ' error: b2sjn called with neg m =', m, ' b2sjn = 0&
     &'                                                                 
         b2sjn = 0 
         RETURN 
      ENDIF 
!                                                                       
      IF (z.eq.0) then 
         b2sjn = 1.d0 
         IF (m.gt.0) b2sjn = 0 
         RETURN 
      ENDIF 
! --- prepare factorial for limiting formula ---                        
      f = 1d0 
      DO 101 i = 1, 2 * m + 1, 2 
         f = f * i 
  101 END DO 
!                                                                       
! ---- treat case z small ----                                          
      zzt = (z**m) / f 
      IF (dabs (zzt) .lt.eps) then 
!        -----> use limiting formula                                    
         b2sjn = (z**m / f) * (1d0 - z**2 / (2 * (2 * m + 3) ) + z**4 / &
         (8 * (2 * m + 3) * (2 * m + 5) ) )                             
!        write(6,*)' asymptotic formula'                                
         RETURN 
      ENDIF 
!                                                                       
! --- treat cases m=0 and m=1 ---                                       
      IF (m.lt.2) then 
         IF (m.eq.0) then 
            b2sjn = dsin (z) / z 
            RETURN 
         ENDIF 
         IF (m.eq.1) then 
            b2sjn = (dsin (z) / z - dcos (z) ) / z 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
! --- treat cases m >= 2 ---                                            
!                                                                       
      fnm1 = 1.d0 / z 
      fn = fnm1**2 
      fmn = 0 
      fmnp1 = fnm1 
      maxm = m - 1 
!                                                                       
      DO 1 n = 1, maxm 
         fnp1 = (2 * n + 1) * fn / z - fnm1 
         fmnm1 = (1 - 2 * n) * fmn / z - fmnp1 
         fnm1 = fn 
         fn = fnp1 
         fmnp1 = fmn 
         fmn = fmnm1 
    1 END DO 
      fmnm1 = (1 - 2 * m) * fmn / z - fmnp1 
      b2sjn = fn * dsin (z) - ( - 1) **m * fmnm1 * dcos (z) 
!                                                                       
      RETURN 
      END FUNCTION b2sjn                            
                                                                        
                                                                        
!*ds                                                                    
!*ed                                                                    
      FUNCTION b2sj_int (q, l, dr, y, ny, r0) 
!     ------------------------------------------------------------      
!     berechnet :                                                       
!                                                                       
!          integral(r0..r0+1) j_l(qr)*(2*r*y(r)+r*r*y'(r))*dr           
!                                                                       
!     ------------------------------------------------------------      
                                                                        
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION y (ny) 
                                                                        
!cc   if(q.gt.30.d0) then                                               
!cc     b2sj_int = 0.d0                                                 
!cc     return                                                          
!cc   endif                                                             
                                                                        
      sum = 0.d0 
      DO i = 1, ny 
      r = r0 + (i - 1) * dr 
      IF (i.lt.ny.and.i.gt.1) then 
         deriv = (y (i + 1) - y (i - 1) ) / (2.0d0 * dr) 
      ELSE 
         IF (i.eq.1) deriv = y (1) / dr 
         IF (i.eq.ny) deriv = - y (ny) / dr 
      ENDIF 
      sum = sum + b2sjn (l, q * r) * (2.0d0 * r * y (i) + r * r * deriv) 
      enddo 
      b2sj_int = sum * dr 
                                                                        
      RETURN 
      END FUNCTION b2sj_int                         
                                                                        
!*ds                                                                    
!*ed                                                                    
                     
                                                                        
                 
                                                                        
                                                                        
      FUNCTION s2_of_q (q, r0, sigmar0, rho_innen, exp_aussen) 
!      ==================================================               
      IMPLICIT real (8)(a - h, o - z) 
      PARAMETER (epsilon = 1.0d-5) 
      PARAMETER (maxiter = 100) 
!ccc   parameter(ninteg  =25)                                           
      PARAMETER (ninteg = 2) 
      PARAMETER (asigma = 3.d0) 
      PARAMETER (qmin = 1.0d-6) 
      PARAMETER (rmin0 = 0.1d0) 
      PARAMETER (pi = 3.141592654d0) 
      COMMON / c2sint / expo, qq 
      EXTERNAL s2int 
                                                                        
      expo = exp_aussen 
      IF (q.lt.qmin) then 
         qq = qmin 
      ELSE 
         qq = q 
      ENDIF 
                                                                        
      rmin = r0 - asigma * sigmar0 
      rmax = r0 + asigma * sigmar0 
      IF (rmin.lt.rmin0 * r0) rmin = rmin0 * r0 
      dr = (rmax - rmin) / ninteg 
      weight = 0.0d0 
      soq = 0.0d0 
                                                                        
      DO i = 0, ninteg 
      r0h = rmin + i * dr 
      r1 = r0h + 1.0d0 
      qr0 = r0h * qq 
      sext = a2dapint (s2int, r0h, r1, epsilon, maxiter, erroraccu) 
      sinn = rho_innen * (dsin (qr0) - qr0 * dcos (qr0) ) / (qq**3) 
      wexp = dexp ( - ( ( (r0h - r0) / sigmar0) **2) ) 
      weight = weight + wexp 
!ccc     write(6,*)wexp,weight, soq                                     
      soq = soq + wexp * (sext + sinn) **2 
      enddo 
!ccc   write(6,*)'weight = ',weight                                     
!ccc   write(6,*)'soq    = ',soq                                        
      soq = soq / weight 
! -- zusammen:                                                          
      s2_of_q = soq * (4 * pi) **2 
                                                                        
      RETURN 
      END FUNCTION s2_of_q                          
                                                                        
      FUNCTION s2int (r) 
!      ================ integrandfunktion fuer s2_of_q                  
      IMPLICIT real (8)(a - h, o - z) 
      COMMON / c2sint / expo, qq 
                                                                        
      qr = qq * r 
      s2int = r** (expo + 2.0d0) * dsin (qr) / qr 
                                                                        
      RETURN 
      END FUNCTION s2int                            
