      SUBROUTINE SOFQ (QQ, SQ, NPT, eta, SCL, GAMMA, R, MODE, IERR) 
!**********************************************************             
!     Subroutine to calculate S(q)                                      
!     Input parameters:                                                 
!               QQ      array of dimension NPT with the values q*r*2.   
!                       where to calculate S(q)                         
!               NPT     number of points                                
!               eta=volume fraction                                     
!               SCL=characteristic length of the interaction            
!                       screening length for mode=1,2                   
!                       the width of the potential well for             
!                                       mode=3                          
!               GAMMA=strength of the interaction                       
!                       mode=1,2 GEK=gamma*exp(-r/scl) contact pot.     
!                       mode=3 depth of the potential well              
!                       in units of kT, + ==> repulsive                 
!                                       - ==> attractive                
!               R=outer radius of the spheres                           
!               MODE=   1       H.-P. form                              
!                       2       Shieu-Chen form                         
!                       3       Sharma-Sharma form                      
!                       4       Critical diverg.                        
!                       5       none (sq=1.)                            
!       Works in double precision                                       
!***********************************************************************
                                                                        
      IMPLICIT REAL (8)(A - H, O - Z) 
      DIMENSION QQ (NPT), SQ (NPT) 
!     COMMON /SQ/ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1      
                                                                        
      IF (scl.ne.0.0) then 
         AK = 2. * R / SCL 
!         AK=R/SCL                                                      
      ELSE 
         AK = 1e20 
      ENDIF 
      GOTO (1, 2, 3, 4, 5) MODE 
      STOP ' Invalid MODE' 
    1 GEK = GAMMA * EXP ( - AK) 
      CALL SQHPA (QQ, SQ, NPT, IERR, eta, ak, gek) 
      RETURN 
    2 GEK = GAMMA * EXP ( - AK) 
      CALL SQMQ (QQ, SQ, NPT, IERR, eta, ak, gek) 
      RETURN 
    3 GEK = GAMMA 
      CALL SQSW (QQ, SQ, NPT, IERR, eta, ak, gek) 
      RETURN 
    4 GEK = GAMMA 
      CALL SQCR (QQ, SQ, NPT, IERR, eta, ak, gek) 
      RETURN 
    5 DO i = 1, npt 
      sq (i) = 1.0 
      enddo 
      ierr = 0 
      RETURN 
      END SUBROUTINE SOFQ                           
!   ################################################################################                                                                        
      SUBROUTINE SQMQ (QQ, SQ, NPT, IERR, eta, ak, gek) 
                                                                        
!       Modified version of the Hayter-Penfold program                  
!       bye E.Sheu      (Exxon Res. & Eng.)                             
                                                                        
      IMPLICIT REAL (8)(A - H, O - Z) 
      DIMENSION RR (4), COEF (5), W (5), RTR (4), RTI (4) 
      DIMENSION QQ (NPT), SQ (NPT) 
      REAL(8) K 
                                                                        
!     Notation is the same as in the H-P program                        
                                                                        
      DATA PI / 3.141592654 / 
                                                                        
!     Calculation of the coefficients of the Beta equation              
                                                                        
      IERR = 0 
      Z = AK 
      EZ = EXP (Z) 
      EK = - GEK 
      Z2 = Z * Z 
      ETA2 = ETA * ETA 
      ETA01 = 1.0 / (1.0 - ETA) 
      ETA02 = ETA01 * ETA01 
      ETA12 = ETA * ETA02 
      ETA22 = ETA * ETA12 
      ZE = EXP ( - Z) 
      TMP1 = 18. * ETA22 / Z2 * (2.0 - Z - ZE * (2.0 + Z) ) 
      TMP2 = ZE * (1. + Z) 
      TMP3 = 1. + 2. * ETA - 6. * ETA / Z 
      TMP4 = 1.5 * ETA + (1. - 4. * ETA) / Z 
      SS = 12. * ETA * (1. - Z2 / 2. - TMP2) / Z2 
      TT = 12. * ETA * (1. - Z - ZE) / Z 
      APY = ETA02 + ETA12 * 2. 
      BPY = - 1.5 * ETA12 
      XX = 6.0 * ETA * (Z * ZE-6. * ETA * ETA01 / Z2 * (2. - 2. * Z -   &
      ZE * (2. - Z2) ) - TMP1)                                          
      YY = Z - 6. * ETA * ETA01 / Z2 * (2. - Z2 - 2. * TMP2) - TMP1 
      DD = Z - APY * SS - BPY * TT 
      EE = - 6. * ETA + 12. * ETA12 / Z * (TMP3 * SS - TMP4 * TT) 
      FF = - 6. * ETA * (1. - ZE) * (1. - ZE) + 12. * ETA12 / Z *       &
      ( (TMP3 * (1. - TMP2) + 3. * ETA * Z * ZE) * SS - (TMP4 * (1. -   &
      TMP2) - (1. - 4. * ETA) / 2 * Z * ZE) * TT)                       
      NR = 1 
                                                                        
!     Coefficients of beta equation                                     
                                                                        
      COEF (1) = EK * EK 
      COEF (2) = - EK * YY 
      COEF (3) = 12. * ETA * EK 
      COEF (4) = - XX 
      COEF (5) = 36. * ETA2 
                                                                        
!     Find the root for beta                                            
                                                                        
      CALL POLRT (COEF, W, 4, RTR, RTI, IERR) 
                                                                        
!     Select the real and physical root                                 
      IF (IERR.NE.0) THEN 
         PRINT * , ' Error in POLRT', IERR 
         RETURN 
      ENDIF 
                                                                        
      DO 2 LL = 1, 4 
         IF (RTI (LL) .EQ.0.0) GOTO 1 
         GOTO 2 
    1    RR (NR) = RTR (LL) 
         NR = NR + 1 
    2 END DO 
      PRINT *, nr 
      DO itmp = 1, nr 
      PRINT *, itmp, abs (rr (itmp) ) 
      enddo 
      IF (NR.EQ.1) GOTO 4 
      NN = NR - 1 
      RT = RR (1) 
      ART = ABS (RR (1) ) 
      DO 3 LN = 2, NN 
         RS = ABS (RR (LN) ) 
         IF (RS.GT.ART) GOTO 3 
         ART = RS 
         RT = RR (LN) 
    3 END DO 
      GOTO 6 
    4 PRINT * , ' No real root was found !' 
      IERR = 5 
      RETURN 
    6 BETA = RT 
!     print *,BETA,ff                                                   
                                                                        
!     Calculate the coefficients of q(r)                                
                                                                        
      D = ( ( - EK + BETA * DD) * ZE+EE * BETA * BETA) / (FF * BETA *   &
      BETA)                                                             
      A = APY + 12. * BETA * EZ * ETA12 / Z * (TMP3 * (D-1. - D * TMP2) &
      + 3. * ETA * D * Z * ZE)                                          
      B = BPY - 12. * BETA * EZ * ETA12 / Z * (TMP4 * (D-1. - D * TMP2) &
      - (1. - 4. * ETA) * .5 * D * Z * ZE)                              
      Q0 = 1. - 12. * ETA * ( - A / 3. - .5 * B + D * BETA * (1. + 1. / &
      Z) - BETA / Z * EXP (Z) * (D-1.) )                                
                                                                        
!     Calculate the structure factor                                    
                                                                        
      DO 7 I = 1, NPT 
         K = QQ (I) 
                                        ! take this as min. q           
         IF (K.LT..005) K = .005 
         F2 = K * K 
         F3 = F2 * K 
         SK = DSIN (K) 
         CK = DCOS (K) 
         BK = BETA * D / (Z * Z + F2) 
         QR = 1. - 12. * ETA * (A * CK / F2 - A * SK / F3 + B * CK / F2 &
         - B / F2 + BETA * D * SK / K + BK * (Z * CK - K * SK - Z * EZ) &
         + BK * EZ * Z / D)                                             
         QM = - 12. * ETA * (A * (CK - 1.) / F3 + A * SK / F2 - A * .5 /&
         K + B * SK / F2 - B / K + BETA * D / K * (1. - CK) + BK *      &
         (Z * SK + K * CK - K * EZ) + K * EZ * BK / D)                  
         SQ (I) = 1. / (QR * QR + QM * QM) 
    7 END DO 
      RETURN 
      END SUBROUTINE SQMQ                           
                                                                        
                                                                        
      SUBROUTINE POLRT (XCOF, COF, M, ROOTR, ROOTI, IER) 
      DIMENSION XCOF ( * ), COF ( * ), ROOTR ( * ), ROOTI ( * ) 
      DOUBLEPRECISION XO, YO, X, Y, XPR, YPR, UX, UY, V, YT, XT, U, XT2,&
      YT2, SUMSQ, DX, DY, TEMP, ALPHA                                   
      DOUBLEPRECISION XCOF, COF, ROOTR, ROOTI 
                                                                        
      IFIT = 0 
      N = M 
      IER = 0 
      IF (XCOF (N + 1) ) 10, 25, 10 
   10 IF (N) 15, 15, 32 
   15 IER = 1 
   20 RETURN 
   25 IER = 4 
      GOTO 20 
   30 IER = 2 
      GOTO 20 
   32 IF (N - 36) 35, 35, 30 
   35 NX = N 
      NXX = N + 1 
      N2 = 1 
      KJ1 = N + 1 
      DO 40 L = 1, KJ1 
         MT = KJ1 - L + 1 
         COF (MT) = XCOF (L) 
   40 END DO 
   45 XO = .00500101 
      YO = .01000101 
      IN = 0 
   50 X = XO 
      XO = - 10.0 * YO 
      YO = - 10. * X 
      X = XO 
      Y = YO 
      IN = IN + 1 
      GOTO 59 
   55 IFIT = 1 
      XPR = X 
      YPR = Y 
   59 ICT = 0 
   60 UX = 0.0 
      UY = 0.0 
      V = 0.0 
      YT = 0.0 
      XT = 1.0 
      U = COF (N + 1) 
      IF (U) 65, 130, 65 
   65 DO 70 I = 1, N 
         L = N - I + 1 
         TEMP = COF (L) 
         XT2 = X * XT - Y * YT 
         YT2 = X * YT + Y * XT 
         U = U + TEMP * XT2 
         V = V + TEMP * YT2 
         FI = I 
         UX = UX + FI * XT * TEMP 
         UY = UY - FI * YT * TEMP 
         XT = XT2 
         YT = YT2 
   70 END DO 
      SUMSQ = UX * UX + UY * UY 
      IF (SUMSQ) 75, 110, 75 
   75 DX = (V * UY - U * UX) / SUMSQ 
      X = X + DX 
      DY = - (U * UY + V * UX) / SUMSQ 
      Y = Y + DY 
   78 IF (DABS (DY) + DABS (DX) - 1.0D-5) 100, 80, 80 
   80 ICT = ICT + 1 
      IF (ICT - 500) 60, 85, 85 
   85 IF (IFIT) 100, 90, 100 
   90 IF (IN - 5) 50, 95, 95 
   95 IER = 3 
      GOTO 20 
  100 DO 105 L = 1, NXX 
         MT = KJ1 - L + 1 
         TEMP = XCOF (MT) 
         XCOF (MT) = COF (L) 
         COF (L) = TEMP 
  105 END DO 
      ITEMP = N 
      N = NX 
      NX = ITEMP 
      IF (IFIT) 120, 55, 120 
  110 IF (IFIT) 115, 50, 115 
  115 X = XPR 
      Y = YPR 
  120 IFIT = 0 
  122 IF (DABS (Y) - 1.0D-4 * DABS (X) ) 135, 125, 125 
  125 ALPHA = X + X 
      SUMSQ = X * X + Y * Y 
      N = N - 2 
      GOTO 140 
  130 X = 0.0 
      NX = NX - 1 
      NXX = NXX - 1 
  135 Y = 0.0 
      SUMSQ = 0.0 
      ALPHA = X 
      N = N - 1 
  140 COF (2) = COF (2) + ALPHA * COF (1) 
  145 DO 150 L = 2, N 
         COF (L + 1) = COF (L + 1) + ALPHA * COF (L) - SUMSQ * COF (L - &
         1)                                                             
  150 END DO 
  155 ROOTI (N2) = Y 
      ROOTR (N2) = X 
      N2 = N2 + 1 
      IF (SUMSQ) 160, 165, 160 
  160 Y = - Y 
      SUMSQ = 0.0 
      GOTO 155 
  165 IF (N) 20, 20, 45 
      END SUBROUTINE POLRT                          
                                                                        
      SUBROUTINE SQHPA (QQ, SQ, NPT, IERR, eta, ak, gek) 
!     ROUTINE TO CALCULATE S(Q*SIG) FOR A SCREENED COULOMB              
!     POTENTIAL BETWEEN FINITE PARTICLES OF DIAMETER 'SIG'              
!     AT ANY VOLUME FRACTION. THIS ROUTINE IS MUCH MORE POWER-          
!     FUL THAN "SQHP" AND SHOULD BE USED TO REPLACE THE LATTER          
!     IN EXISTING PROGRAMS. NOTE THAT THE COMMON AREA IS                
!     CHANGED; IN PARTICULAR, THE POTENTIAL IS PASSED                   
!     DIRECTLY AS 'GEK' = GAMMA*EXP(-K) IN THE PRESENT ROUTINE.         
!                                                                       
!     JOHN B.HAYTER (I.L.L.) 19-AUG-81                                  
!                                                                       
!     CALLING SEQUENCE:                                                 
!                                                                       
!      CALL SQHPA(QQ,SQ,NPT,IERR)                                       
!                                                                       
!     QQ: ARRAY OF DIMENSION NPT CONTAINING THE VALUES                  
!      OF Q*SIG AT WHICH S(Q*SIG) WILL BE CALCULATED.                   
!     SQ: ARRAY OF DIMENSION NPT INTO WHICH VALUES OF                   
!      S(Q*SIG) WILL BE RETURNED.                                       
!     NPT: NUMBER OF VALUES OF Q*SIG.                                   
!                                                                       
!     IERR > 0: NORMAL EXIT; IERR=NUMBER OF ITERATIONS.                 
!      -1: NEWTON ITERATION NON-CONVERGENT IN "SQCOEF"                  
!      -2: NEWTON ITERATION NON-CONVERGENT IN "SQFUN".                  
!      -3: CANNOT RESCALE TO G(1+) > 0.                                 
!                                                                       
!     ALL OTHER PARAMETERS ARE TRANSMITTED THROUGH A SINGLE             
!     NAMED COMMON AREA:                                                
!                                                                       
!     REAL*8 A,B,C,F                                                    
!     COMMON /   / ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1    
!                                                                       
!     ON ENTRY:                                                         
!                                                                       
!     ETA: VOLUME FRACTION                                              
!     GEK: THE CONTACT POTENTIAL GAMMA*EXP(-K)                          
!     AK: THE DIMENSIONLESS SCREENING CONSTANT                          
!     AK = KAPPA*SIG WHERE KAPPA IS THE INVERSE SCREENING               
!      LENGTH AND SIG IS THE PARTICLE DIAMETER.                         
!                                                                       
!     ON EXIT:                                                          
!                                                                       
!     GAMK IS THE COUPLING: 2*GAMMA*S*EXP(-K/S), S=ETA**(1/3).          
!     SETA, SGEK AND SAK ARE THE RESCALED INPUT PARAMETERS.             
!     SCAL IS THE RESCALING FACTOR: (ETA/SETA)**(1/3).                  
!     G1=G(1+), THE CONTACT VALUE OF G(R/SIG).                          
!     A,B,C,F,U,V ARE THE CONSTANTS APPEARING IN THE ANALYTIC           
!     SOLUTION OF THE MSA (HAYTER-PENFOLD; MOL.PHYS. 42: 109 (1981))    
!                                                                       
!     NOTES:                                                            
!                                                                       
!     (A) AFTER THE FIRST CALL TO SQHPA, S(Q*SIG) MAY BE EVALUATED      
!     AT OTHER Q*SIG VALUES BY REDEFINING THE ARRAY QQ AND CALLING      
!     "SQHCAL" DIRECTLY FROM THE MAIN PROGRAM.                          
!                                                                       
!     (B) THE RESULTING S(Q*SIG) MAY BE TRANSFORMED TO G(R/SIG)         
!     USING THE ROUTINE "TROGS".                                        
!                                                                       
!     (C) NO ERROR CHECKING OF INPUT PARAMETERS IS PERFORMED;           
!     IT IS THE RESPONSIBILITY OF THE CALLING PROGRAM TO VERIFY         
!     VALIDITY.                                                         
!     SUBROUTINES REQUIRED BY SQHPA:                                    
!                                                                       
!     (1) SQCOEF RESCALES THE PROBLEM AND CALULATES THE                 
!      APPROPRIATE COEFFICIENTS FOR "SQHCAL".                           
!                                                                       
!     (2) SQFUN CALCULATES VARIOUS VALUES FOR "SQCOEF".                 
!                                                                       
!     (3) SQHCAL CALCULATES H-P S(Q*SIG) GIVEN A,B,C,F.                 
!                                                                       
!                                                                       
!                                                                       
      IMPLICIT REAL (8)(A - H, O - Z) 
      DIMENSION QQ ( * ), SQ ( * ) 
!                                                                       
!     COMMON /SQ/ ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1     
!                                                                       
      IERR = 0 
!                                                                       
!     FIRST CALCULATE COUPLING.                                         
!                                                                       
      S = ETA**0.3333333 
      GAMK = 2.0 * S * GEK * EXP (AK - AK / S) 
!                                                                       
!     CALCULATE COEFFICIENTS, CHECK ALL IS WELL                         
!                                                                       
   10 IERR = 0 
      CALL SQCOEF (IERR, ETA, GEK, AK, A, B, C, F, U, V, GAMK, SETA,    &
      SGEK, SAK, SCAL, G1)                                              
                                                                        
      IF (IERR.GE.0) GOTO 20 
      WRITE ( *, 9) IERR 
    9 FORMAT (' ERROR IN SQCOEF: IERR=',I5) 
      RETURN 
!                                                                       
!     CALCULATE S(Q*SIG)                                                
!                                                                       
   20 CALL SQHCAL (QQ, SQ, NPT, ETA, GEK, AK, A, B, C, F, U, V, GAMK,   &
      SETA, SGEK, SAK, SCAL, G1)                                        
!                                                                       
      RETURN 
      END SUBROUTINE SQHPA                          
      SUBROUTINE SQCOEF (IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK, SETA,&
      SGEK, SAK, SCAL, G1)                                              
!     CALCULATES RESCALED VOLUME FRACTION AND CORRESPONDING             
!     COEFFICIENTS FOR "SQHPA".                                         
!                                                                       
!     JOHN B.HAYTER (I.L.L.) 14-SEP-81                                  
!                                                                       
!     ON EXIT:                                                          
!                                                                       
!     SETA IS THE RESCALED VOLUME FRACTION.                             
!     SGEK IS THE RESCALED CONTACT POTENTIAL.                           
!     SAK IS THE RESCALED SCREENING CONSTANT.                           
!     A,B,C,F,U,V ARE THE MSA COEFFICIENTS.                             
!     G1=G(1+) IS THE CONTACT VALUE OF G(R/SIG);                        
!     FOR THE GILLAN CONDITION, THE DIFFERENCE FROM                     
!     ZERO INDICATES THE COMPUTATIONAL ACCURACY.                        
!                                                                       
!     IR > 0: NORMAL EXIT, IR IS THE NUMBER OF ITERATIONS.              
!      < 0: FAILED TO CONVERGE.                                         
!                                                                       
!                                                                       
      IMPLICIT REAL (8)(A - H, O - Z) 
      REAL(8) A, B, C, F 
!     COMMON /SQ/ ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1     
!                                                                       
      DATA ITM / 40 /, ACC / 5.0D-6 / 
      IG = 1 
      IF (AK.LT. (1.0 + 8.0 * ETA) ) GOTO 10 
      IG = 0 
      CALL SQFUN (G1, ETA, 1, IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK, &
      SETA, SGEK, SAK, SCAL, G1)                                        
      IF (IR.LT.0.OR.G1.GE.0.0) RETURN 
   10 SETA = DMIN1 (ETA, 2.D-1) 
      IF (IG.EQ.1.AND.GAMK.LT.0.15) GOTO 50 
      J = 0 
   20 J = J + 1 
      IF (J.GT.ITM) GOTO 30 
      IF (SETA.LE.0.0) SETA = ETA / J 
      IF (SETA.GT.0.6) SETA = 0.35 / J 
      E1 = SETA 
      CALL SQFUN (F1, E1, 2, IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK,  &
      SETA, SGEK, SAK, SCAL, G1)                                        
      E2 = SETA * 1.01 
      CALL SQFUN (F2, E2, 2, IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK,  &
      SETA, SGEK, SAK, SCAL, G1)                                        
      E2 = E1 - (E2 - E1) * F1 / (F2 - F1) 
      SETA = E2 
      DEL = ABS ( (E2 - E1) / E1) 
      IF (DEL.GT.ACC) GOTO 20 
      CALL SQFUN (G1, E2, 4, IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK,  &
      SETA, SGEK, SAK, SCAL, G1)                                        
      IR = J 
      IF (IG.EQ.1) GOTO 40 
      RETURN 
   30 IR = - 1 
      RETURN 
   40 IF (SETA.GE.ETA) RETURN 
   50 CALL SQFUN (G1, ETA, 3, IR, ETA, GEK, AK, A, B, C, F, U, V, GAMK, &
      SETA, SGEK, SAK, SCAL, G1)                                        
      IF (IR.LT.0.OR.G1.GE.0.0) RETURN 
      IR = - 3 
      RETURN 
!                                                                       
      END SUBROUTINE SQCOEF                         
      SUBROUTINE SQFUN (FVAL, EVAR, IX, IR, RETA, RGEK, RAK, A, B, C, F,&
      U, V, GAMK, SETA, SGEK, SAK, SCAL, dummy)                         
                                                                 !! aix 
!     CALCULATES VARIOUS COEFFICIENTS AND FUNCTION                      
!     VALUES FOR "SQCOEF" (USED BY "SQHPA").                            
!                                                                       
!     ** THIS ROUTINE WORKS LOCALLY IN DOUBLE PRECISION **              
!                                                                       
!     JOHN B.HAYTER (I.L.L.) 23-MAR-82                                  
!                                                                       
!     IX=1: SOLVE FOR LARGE K, RETURN G(1+).                            
!      2: RETURN FUNCTION TO SOLVE FOR ETA(GILLAN).                     
!      3: ASSUME NEAR GILLAN, SOLVE, RETURN G(1+).                      
!      4: RETURN G(1+) FOR ETA=ETA(GILLAN).                             
!                                                                       
!                                                                       
      IMPLICIT REAL (8)(A - H, O - Z) 
!                                                                       
!     COMMON /SQ/RETA,RGEK,RAK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL      
!                                                                       
!     CALCULATE CONSTANTS; NOTATION IS HAYTER-PENFOLD (1981).           
!                                                                       
      DATA ACC / 1.0D-6 /, ITM / 40 / 
      ETA = EVAR 
!     WRITE(*,9)ETA                                                     
!9    FORMAT(' ETA= ',F12.7)                                            
      ETA2 = ETA * ETA 
      ETA3 = ETA * ETA2 
      E12 = 12.0D0 * ETA 
      E24 = E12 + E12 
      SCAL = (RETA / EVAR) **0.33333333 
      SAK = RAK / SCAL 
      IBIG = 0 
      IF (SAK.GT.15.0.AND.IX.EQ.1) IBIG = 1 
      SGEK = RGEK * SCAL * EXP (RAK - SAK) 
      GEK = SGEK 
      AK = SAK 
      AK2 = AK * AK 
      AK1 = 1.0D0 + AK 
      DAK2 = 1.0D0 / AK2 
      DAK4 = DAK2 * DAK2 
      D = 1.0D0 - ETA 
      D2 = D * D 
      DAK = D / AK 
      DD2 = 1.0D0 / D2 
      DD4 = DD2 * DD2 
      DD45 = DD4 * 2.0D-1 
      ETA3D = 3.0D0 * ETA 
      ETA6D = ETA3D+ETA3D 
      ETA32 = ETA3 + ETA3 
      ETA2D = ETA + 2.0D0 
      ETA2D2 = ETA2D * ETA2D 
      ETA21 = 2.0D0 * ETA + 1.0D0 
      ETA22 = ETA21 * ETA21 
!                                                                       
!     ALPHA(I)                                                          
!                                                                       
      AL1 = - ETA21 * DAK 
      AL2 = (14.0D0 * ETA2 - 4.0D0 * ETA - 1.0D0) * DAK2 
      AL3 = 36.0D0 * ETA2 * DAK4 
!                                                                       
!     BETA(I)                                                           
!                                                                       
      BE1 = - (ETA2 + 7.0D0 * ETA + 1.0D0) * DAK 
      BE2 = 9.0D0 * ETA * (ETA2 + 4.0D0 * ETA - 2.0D0) * DAK2 
      BE3 = 12.0D0 * ETA * (2.0D0 * ETA2 + 8.0D0 * ETA - 1.0D0) * DAK4 
!                                                                       
!     NU(I)                                                             
!                                                                       
      VU1 = - (ETA3 + 3.0D0 * ETA2 + 45.0D0 * ETA + 5.0D0) * DAK 
      VU2 = (ETA32 + 3.0D0 * ETA2 + 42.0D0 * ETA - 2.0D1) * DAK2 
      VU3 = (ETA32 + 3.0D1 * ETA - 5.0D0) * DAK4 
      VU4 = VU1 + E24 * AK * VU3 
      VU5 = ETA6D * (VU2 + 4.0D0 * VU3) 
!                                                                       
!     PHI(I)                                                            
!                                                                       
      PH1 = ETA6D / AK 
      PH2 = D-E12 * DAK2 
!                                                                       
!     TAU(I)                                                            
!                                                                       
      TA1 = (ETA + 5.0D0) / (5.0D0 * AK) 
      TA2 = ETA2D * DAK2 
      TA3 = - E12 * GEK * (TA1 + TA2) 
      TA4 = ETA3D * AK2 * (TA1 * TA1 - TA2 * TA2) 
      TA5 = ETA3D * (ETA + 8.0D0) * 1.0D-1 - 2.0D0 * ETA22 * DAK2 
!                                                                       
!     DOUBLE PRECISION SINH(K), COSH(K)                                 
!                                                                       
      EX1 = DEXP (AK) 
      EX2 = 0.0D0 
      IF (SAK.LT.20.0) EX2 = DEXP ( - AK) 
      SK = 0.5D0 * (EX1 - EX2) 
      CK = 0.5D0 * (EX1 + EX2) 
      CKMA = CK - 1.0D0 - AK * SK 
      SKMA = SK - AK * CK 
!                                                                       
!     A(I)                                                              
!                                                                       
      A1 = (E24 * GEK * (AL1 + AL2 + AK1 * AL3) - ETA22) * DD4 
      IF (IBIG.NE.0) GOTO 10 
      A2 = E24 * (AL3 * SKMA + AL2 * SK - AL1 * CK) * DD4 
      A3 = E24 * (ETA22 * DAK2 - 0.5D0 * D2 + AL3 * CKMA - AL1 * SK +   &
      AL2 * CK) * DD4                                                   
!                                                                       
!     B(I)                                                              
!                                                                       
   10 B1 = (1.5D0 * ETA * ETA2D2 - E12 * GEK * (BE1 + BE2 + AK1 * BE3) )&
      * DD4                                                             
      IF (IBIG.NE.0) GOTO 20 
      B2 = E12 * ( - BE3 * SKMA - BE2 * SK + BE1 * CK) * DD4 
      B3 = E12 * (0.5D0 * D2 * ETA2D-ETA3D * ETA2D2 * DAK2 - BE3 * CKMA &
      + BE1 * SK - BE2 * CK) * DD4                                      
!                                                                       
!     V(I)                                                              
!                                                                       
   20 V1 = (ETA21 * (ETA2 - 2.0D0 * ETA + 1.0D1) * 2.5D-1 - GEK *       &
      (VU4 + VU5) ) * DD45                                              
      IF (IBIG.NE.0) GOTO 30 
      V2 = (VU4 * CK - VU5 * SK) * DD45 
      V3 = ( (ETA3 - 6.0D0 * ETA2 + 5.0D0) * D-ETA6D * (2.0D0 * ETA3 -  &
      3.0D0 * ETA2 + 18.0D0 * ETA + 1.0D1) * DAK2 + E24 * VU3 + VU4 *   &
      SK - VU5 * CK) * DD45                                             
!                                                                       
   30 PP1 = PH1 * PH1 
      PP2 = PH2 * PH2 
      PP = PP1 + PP2 
      P1P2 = PH1 * PH2 * 2.0D0 
!                                                                       
!     P(I)                                                              
!                                                                       
      P1 = (GEK * (PP1 + PP2 - P1P2) - 0.5D0 * ETA2D) * DD2 
      IF (IBIG.NE.0) GOTO 40 
      P2 = (PP * SK + P1P2 * CK) * DD2 
      P3 = (PP * CK + P1P2 * SK + PP1 - PP2) * DD2 
!                                                                       
!     T(I)                                                              
!                                                                       
   40 T1 = TA3 + TA4 * A1 + TA5 * B1 
      IF (IBIG.NE.0) GOTO 140 
      T2 = TA4 * A2 + TA5 * B2 + E12 * (TA1 * CK - TA2 * SK) 
      T3 = TA4 * A3 + TA5 * B3 + E12 * (TA1 * SK - TA2 * (CK - 1.0D0) ) &
      - 4.0D-1 * ETA * (ETA + 1.0D1) - 1.0D0                            
!                                                                       
!     MU(I)                                                             
!                                                                       
      UM1 = T2 * A2 - E12 * V2 * V2 
      UM2 = T1 * A2 + T2 * A1 - E24 * V1 * V2 
      UM3 = T2 * A3 + T3 * A2 - E24 * V2 * V3 
      UM4 = T1 * A1 - E12 * V1 * V1 
      UM5 = T1 * A3 + T3 * A1 - E24 * V1 * V3 
      UM6 = T3 * A3 - E12 * V3 * V3 
!                                                                       
!     GILLAN CONDITION?                                                 
!                                                                       
!     YES - G(X=1+)=0                                                   
!                                                                       
!     COEFFICIENTS AND FUNCTION VALUE.                                  
!                                                                       
      IF (IX.NE.0.AND.1.NE.0) GOTO 70 
      CA = AK2 * P1 + 2.0D0 * (B3 * P1 - B1 * P3) 
      CA = - CA / (AK2 * P2 + 2.0D0 * (B3 * P2 - B2 * P3) ) 
      FA = - (P1 + P2 * CA) / P3 
      IF (IX.EQ.2) FVAL = UM1 * CA * CA + (UM2 + UM3 * FA) * CA + UM4 + &
      UM5 * FA + UM6 * FA * FA                                          
      IF (IX.EQ.4) FVAL = - (P1 + P2 * CA + P3 * FA) 
   50 F = FA 
      C = CA 
      B = B1 + B2 * CA + B3 * FA 
      A = A1 + A2 * CA + A3 * FA 
      V = (V1 + V2 * CA + V3 * FA) / A 
   60 G24 = E24 * GEK * EX1 
      U = (AK * AK2 * CA - G24) / (AK2 * G24) 
      RETURN 
!                                                                       
!     NO - CALCULATE REMAINING COEFFICIENTS.                            
!                                                                       
!     LAMBDA(I)                                                         
!                                                                       
   70 AL1 = E12 * P2 * P2 
      AL2 = E24 * P1 * P2 - B2 - B2 
      AL3 = E24 * P2 * P3 
      AL4 = E12 * P1 * P1 - B1 - B1 
      AL5 = E24 * P1 * P3 - B3 - B3 - AK2 
      AL6 = E12 * P3 * P3 
!                                                                       
!     OMEGA(I,J)                                                        
!                                                                       
      W16 = UM1 * AL6 - AL1 * UM6 
      W15 = UM1 * AL5 - AL1 * UM5 
      W14 = UM1 * AL4 - AL1 * UM4 
      W13 = UM1 * AL3 - AL1 * UM3 
      W12 = UM1 * AL2 - AL1 * UM2 
      W26 = UM2 * AL6 - AL2 * UM6 
      W25 = UM2 * AL5 - AL2 * UM5 
      W24 = UM2 * AL4 - AL2 * UM4 
      W36 = UM3 * AL6 - AL3 * UM6 
      W35 = UM3 * AL5 - AL3 * UM5 
      W34 = UM3 * AL4 - AL3 * UM4 
      W32 = UM3 * AL2 - AL3 * UM2 
      W46 = UM4 * AL6 - AL4 * UM6 
      W56 = UM5 * AL6 - AL5 * UM6 
      W3526 = W35 + W26 
      W3425 = W34 + W25 
!                                                                       
!     QUARTIC COEFFICIENTS W(I)                                         
!                                                                       
      W4 = W16 * W16 - W13 * W36 
      W3 = 2.0D0 * W16 * W15 - W13 * W3526 - W12 * W36 
      W2 = W15 * W15 + 2.0D0 * W16 * W14 - W13 * W3425 - W12 * W3526 
      W1 = 2.0D0 * W15 * W14 - W13 * W24 - W12 * W3425 
      W0 = W14 * W14 - W12 * W24 
!                                                                       
!     ESTIMATE THE STARTING VALUE OF F                                  
!                                                                       
      IF (IX.EQ.1) GOTO 90 
!                                                                       
!     ASSUME NOT TOO FAR FROM GILLAN CONDITION.                         
!     IF BOTH GEK AND AK ARE SMALL, USE P-W ESTIMATE.                   
!                                                                       
      G1 = 0.5D0 * ETA2D * DD2 * DEXP ( - GEK) 
      IF (SGEK.GT.2.0.OR.SAK.GT.1.0) GOTO 80 
      E24G = E24 * GEK * DEXP (AK) 
      PWK = DSQRT (E24G) 
      QPW = (1.0D0 - DSQRT (1.0D0 + 2.0D0 * D2 * D * PWK / ETA22) )     &
      * ETA21 / D                                                       
      G1 = - QPW * QPW / E24 + 0.5D0 * ETA2D * DD2 
   80 PG = P1 + G1 
      CA = AK2 * PG + 2.0D0 * (B3 * PG - B1 * P3) + E12 * G1 * G1 * P3 
      CA = - CA / (AK2 * P2 + 2.0D0 * (B3 * P2 - B2 * P3) ) 
      FAP = - (PG + P2 * CA) / P3 
      GOTO 100 
!                                                                       
!     LARGE K.                                                          
!                                                                       
   90 FAP = (W14 - W34 - W46) / (W12 - W15 + W35 - W26 + W56 - W32) 
!                                                                       
!     AND REFINE IT ACCORDING TO NEWTON                                 
!                                                                       
  100 I = 0 
  110 I = I + 1 
      IF (I.GT.ITM) GOTO 120 
      FA = FAP 
      FUN = W0 + (W1 + (W2 + (W3 + W4 * FA) * FA) * FA) * FA 
      FUND = W1 + (2.0D0 * W2 + (3.0D0 * W3 + 4.0D0 * W4 * FA) * FA)    &
      * FA                                                              
      FAP = FA - FUN / FUND 
      IF (DABS ( (FAP - FA) / FA) .GT.ACC) GOTO 110 
      IR = IR + I 
      GOTO 130 
!                                                                       
!     FAILED TO CONVERGE IN ITM ITERATIONS                              
!                                                                       
  120 IR = - 2 
      RETURN 
!                                                                       
  130 FA = FAP 
      CA = - (W16 * FA * FA + W15 * FA + W14) / (W13 * FA + W12) 
      G1 = - (P1 + P2 * CA + P3 * FA) 
      FVAL = G1 
      IF (ABS (FVAL) .LT.1.0E-3) FVAL = 0.0 
      SETA = EVAR 
      GOTO 50 
!                                                                       
!     VERY LARGE SCREENING.                                             
!                                                                       
  140 V3 = ( (ETA3 - 6.0D0 * ETA2 + 5.0D0) * D-ETA6D * (2.0D0 * ETA3 -  &
      3.0D0 * ETA2 + 18.0D0 * ETA + 1.0D1) * DAK2 + E24 * VU3) * DD45   
      T3 = TA4 * A3 + TA5 * B3 + E12 * TA2 - 4.0D-1 * ETA * (ETA +      &
      1.0D1) - 1.0D0                                                    
      P3 = (PP1 - PP2) * DD2 
      B3 = E12 * (0.5D0 * D2 * ETA2D-ETA3D * ETA2D2 * DAK2 + BE3)       &
      * DD4                                                             
      A3 = E24 * (ETA22 * DAK2 - 0.5D0 * D2 - AL3) * DD4 
      UM6 = T3 * A3 - E12 * V3 * V3 
      UM5 = T1 * A3 + A1 * T3 - E24 * V1 * V3 
      UM4 = T1 * A1 - E12 * V1 * V1 
      AL6 = E12 * P3 * P3 
      AL5 = E24 * P1 * P3 - B3 - B3 - AK2 
      AL4 = E12 * P1 * P1 - B1 - B1 
      W56 = UM5 * AL6 - AL5 * UM6 
      W46 = UM4 * AL6 - AL4 * UM6 
      FA = - W46 / W56 
      CA = - FA 
      F = FA 
      C = CA 
      B = B1 + B3 * FA 
      A = A1 + A3 * FA 
      V = V1 + V3 * FA 
      G1 = - (P1 + P3 * FA) 
      FVAL = G1 
      IF (ABS (FVAL) .LT.1.0E-3) FVAL = 0.0 
      SETA = EVAR 
      GOTO 60 
!                                                                       
      END SUBROUTINE SQFUN                          
      SUBROUTINE SQHCAL (QQ, SQ, NPT, ETAZ, GEKZ, AKZ, A, B, C, F, U, V,&
      GAMK, SETA, SGEK, SAK, SCAL, dummy)                               
                                                                 !! aix 
!     CALCULATES S(Q) FOR "SQHPA"                                       
!                                                                       
!     ** THIS ROUTINE WORKS LOCALLY IN DOUBLE PRECISION **              
!                                                                       
!     JOHN B.HAYTER (I.L.L.) 19-AUG-81                                  
!                                                                       
!                                                                       
      IMPLICIT REAL (8)(A - H, O - Z) 
!     REAL*4 QQ(1),SQ(1),ETAZ,GEKZ,AKZ,U,V,GAMK,SETA,SGEK,SAK,SCAL      
      DIMENSION QQ ( * ), SQ ( * ) 
!                                                                       
!     COMMON /SQ/ ETAZ,GEKZ,AKZ,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL     
!                                                                       
      ETA = SETA 
      AK = SAK 
      GEK = SGEK 
      E24 = 24.0D0 * ETA 
      X1 = DEXP (AK) 
      X2 = 0.0D0 
      IF (SAK.LT.20.0) X2 = DEXP ( - AK) 
      CK = 0.5D0 * (X1 + X2) 
      SK = 0.5D0 * (X1 - X2) 
      AK2 = AK * AK 
!                                                                       
      DO 20 I = 1, NPT 
         IF (QQ (I) .LE.0.0) GOTO 10 
         QK = QQ (I) / SCAL 
         Q2K = QK * QK 
         QK2 = 1.0D0 / Q2K 
         QK3 = QK2 / QK 
         QQK = 1.0D0 / (QK * (Q2K + AK2) ) 
         SINK = DSIN (QK) 
         COSK = DCOS (QK) 
         ASINK = AK * SINK 
         QCOSK = QK * COSK 
         AQK = (A * (SINK - QCOSK) + B * ( (2.0D0 * QK2 - 1.0D0)        &
         * QCOSK + 2.0D0 * SINK - 2.0D0 / QK) + 0.5D0 * ETA * A *       &
         (24.0D0 * QK3 + 4.0D0 * (1.0D0 - 6.0D0 * QK2) * SINK - (1.0D0 -&
         12.0D0 * QK2 + 24.0D0 * QK2 * QK2) * QCOSK) ) * QK3 + C *      &
         (CK * ASINK - SK * QCOSK) * QQK + F * (SK * ASINK - QK *       &
         (CK * COSK - 1.0D0) ) * QQK + F * (COSK - 1.0D0) * QK2 - GEK * &
         (ASINK + QCOSK) * QQK                                          
         SQ (I) = 1.0D0 / (1.0D0 - E24 * AQK) 
         GOTO 20 
   10    SQ (I) = - 1.0D0 / A 
   20 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE SQHCAL                         
                                                                        
                                                                        
      SUBROUTINE SQSW (QQ, SQ, NPT, IERR, eta, ak, gek) 
!******************************************************                 
!     Subroutine to calculate S(q) for a square well potential          
!       after:  R.V.Sharma,K.C.Sharma Physica 89A 213 (1977)            
!       variables as in the Hayter-Penfold SQHPA and passed             
!       through the COMMON /SQ/                                         
!******************************************************                 
      IMPLICIT REAL (8)(A - H, O - Z) 
      DIMENSION QQ (NPT), SQ (NPT) 
!     COMMON /SQ/ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1      
      FLAMB = 1. + 1. / AK 
      ETA4 = (1. - ETA) **4 
      ALPHA = ( (1. + 2. * ETA) **2 + ETA**3 * (ETA - 4.) ) / ETA4 
      BETA = - 1. / 3. * ETA * (18. + 20. * ETA - 12 * ETA**2 + ETA**4) &
      / ETA4                                                            
      GAMMA = .5 * ETA * ALPHA 
                                                                        
      DO I = 1, NPT 
      QD = QQ (I) 
                                    !say the smallest QD                
      IF (QD.LT.0.01) QD = .01 
      QDL = QD * FLAMB 
      SIQD = DSIN (QD) 
      COQD = DCOS (QD) 
      SIQDL = DSIN (QDL) 
      COQDL = DCOS (QDL) 
      S1 = ALPHA * (SIQD-QD * COQD) / QD**3 
      S2 = BETA * (2. * QD * SIQD- (QD**2 - 2.) * COQD-2.) / QD**4 
      S3 = GAMMA * ( (4. * QD**3 - 24. * QD) * SIQD- (QD**4 - 12. * QD**&
      2 + 24.) * COQD+24.) / QD**6                                      
                                                    ! the '-' is in GEK 
      S4 = GEK * (SIQDL - QDL * COQDL + QD * COQD-SIQD) / QD**3 
      CK = - (S1 + S2 + S3 + S4) * 24. * ETA 
      IF (CK.NE.1.0) THEN 
         SQ (I) = 1. / (1. - CK) 
         IERR = 0 
      ELSE 
         SQ (I) = 1. 
         IERR = 1 
         PRINT * , ' S(q) singular !' 
      ENDIF 
      ENDDO 
      RETURN 
      END SUBROUTINE SQSW                           
                                                                        
                                                                        
                                                                        
      SUBROUTINE SQCR (QQ, SQ, NPT, IERR, eta, ak, gek) 
!******************************************************                 
!     Subroutine to calculate S(q) for a critical system                
!                                                                       
!       variables as in the Hayter-Penfold SQHPA and passed             
!       through the COMMON /SQ/                                         
!******************************************************                 
      IMPLICIT REAL (8)(A - H, O - Z) 
      DIMENSION QQ (NPT), SQ (NPT) 
!     COMMON /SQ/ETA,GEK,AK,A,B,C,F,U,V,GAMK,SETA,SGEK,SAK,SCAL,G1      
      ETA4 = (1. - ETA) **4 
      ALPHA = ( (1. + 2. * ETA) **2 + ETA**3 * (ETA - 4.) ) / ETA4 
      BETA = - 1. / 3. * ETA * (18. + 20. * ETA - 12 * ETA**2 + ETA**4) &
      / ETA4                                                            
      GAMMA = .5 * ETA * ALPHA 
      DO I = 1, NPT 
      QD = QQ (I) 
      critic = - gek / (1. + (qd / ak) **2) 
                                        !say the smallest QD for HS term
      IF (QD.LT.0.003) qd = .003 
      SIQD = DSIN (QD) 
      COQD = DCOS (QD) 
      CK = ALPHA * QD**3 * (SIQD-QD * COQD) + BETA * QD**2 * (2. * QD * &
      SIQD- (QD**2 - 2.) * COQD-2.) + GAMMA * ( (4. * QD**3 - 24. * QD) &
      * SIQD- (QD**4 - 12. * QD**2 + 24.) * COQD+24.)                   
      CK = - CK * 24. * ETA / QD**6 
      IF (CK.NE.1.0) THEN 
         SQ (I) = 1. / (1. - CK) + critic 
         IERR = 0 
      ELSE 
         SQ (I) = 1. 
         IERR = 1 
         PRINT * , ' S(q) singular !' 
      ENDIF 
      ENDDO 
      RETURN 
      END SUBROUTINE SQCR                           
