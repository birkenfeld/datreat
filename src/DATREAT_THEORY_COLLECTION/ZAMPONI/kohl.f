C                                                                       TH111170
      FUNCTION TH15(QQ,OM,TEMP,PA,THNAM,PARNAM,NPAR,INI)                TH111180
C     ==================================================                TH111190
C                                                                       TH111200
C -------> KOHLRAUSCH <-----
C
C     KOHLRAUSCH-RELAXATION WITH ARBITRARY BETA FROM [0.1,1.0]:
C        PHI(T)=F0*EXP(-(T/TAU)**BETA)
C     A AND F0 ARE RELATED BY A=2*S*F0 WHERE S IS THE STATIC
C     STRUCTURE FACTOR (DEFINITION COMPATIBLE TO TH18). THE
C     SCATTERING FUNCTION IS AVERAGED OVER THE CHANNEL WIDTH
C     IN ORDER TO PROTECT AGAINST INTENSITY LOSS AROUND OMEGA=0
C     (LIKE IN THLOR).
C     PARAMETER 5 DETERMINES WHETHER TAU IS THE KOHLRAUSCH-
C     TAU AS DEFINED ABOVE (PA(5)=0.) OR THE AVERAGE RELAXATION
C     TIME DEFINED BY <TAU>=(1/BETA)*GAMMA(1/BETA)*TAU (PA(5)<>0.).
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
      IMPLICIT NONE
!                                                                       !
      INCLUDE  'natur_konstanten_defs.f'
!                                                                       !
      DOUBLE PRECISION om, qq, temp, th15
      DOUBLE PRECISION pa
      DOUBLE PRECISION skohlA, skohl1
!                                                                       !
      INTEGER*4        ini, npar, nparx
!                                                                       !
      CHARACTER*8      thnam, parnam(10)
!                                                                       !
      DIMENSION        pa(10)
!                                                                       !
      EXTERNAL         diro 
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
C
C ----- INITIALISATION -----                                            TH111270
       IF(INI.EQ.0) THEN                                                TH111280
         THNAM = 'KOHL '                                                TH111290
C        ----------------> GIVE YOUR 8 CHARACTER NAME OF THEORY         TH111300
         NPARX = 7                                                      TH111310
C        ----------------> NUMBER OF REQUIRED PARAMETERS                TH111320
         IF(NPAR.LT.NPARX) THEN                                         TH111330
           WRITE(6,1)THNAM,NPARX,NPAR                                   TH111340
1          FORMAT(' THEORY: ',A8,' NO OF PARAMETRS=',I8,                TH111350
     *      ' EXCEEDS CURRENT MAX. = ',I8)                              TH111360
           th15 = 0.0d0
           RETURN                                                       TH111370
         ENDIF                                                          TH111380
         NPAR = NPARX                                                   TH111390
C        --------------> SET THE NUMBER OF PARAMETERS                   TH111400
         PARNAM(1) = 'A'                                                TH111410
         PARNAM(2) = 'TAU0'                                             TH111410
         PARNAM(3) = 'BETA'                                             TH111410
         PARNAM(4) = 'OMOFF'                                            TH111410
         PARNAM(5) = 'TAUTYP'                                           TH111410
         PARNAM(6) = 'TAUQEXP'                                          TH111410
         PARNAM(7) = 'UQUADRAT'
         TH15      = 0.0D0
         RETURN                                                         TH111440
       ENDIF                                                            TH111450
C                                                                       TH111460
C ---- CALCULATE THEORY HERE -----                                      TH111470
C                                                                       TH111480
       TH15 = SKOHL1(QQ,OM,TEMP,PA)                                     TH111490
C      ============================                                     TH111500
C                                                                       TH111510
       TH15 = TH15 * EXP(OM*DBF/TEMP)                                   TH111550
       TH15 = TH15 * EXP(-qq*PA(7)) *2.0d-9
C                                                                       TH111560
       RETURN                                                           TH111570
       END                                                              TH111580
 
 
 
 
      DOUBLE PRECISION FUNCTION SKOHL1(QQ,OM,TEMP,PA)
C     ===============================================
C
C     SCATTERING FROM KOHLRAUSCH RELAXATION
C     (INTERFACE TOFSYS-->SKOHL)
C
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
      IMPLICIT NONE
      INCLUDE  'natur_konstanten_defs.f'
!                                                                       !
      DOUBLE PRECISION omold, omold2, betaold, taufac, qq
      DOUBLE PRECISION a, beta, dgamma, omr, om, temp, tau
      DOUBLE PRECISION pa
      DOUBLE PRECISION skohl
      DOUBLE PRECISION dom
!                                                                       !
      DIMENSION        pa(10)
!                                                                       !
      DATA OMOLD/1.D30/,OMOLD2/1.D30/,BETAOLD/1.D0/,TAUFAC/1.D0/
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
 
      BETA=ABS(PA(3))
      IF (BETA.NE.BETAOLD) THEN
        IF (BETA.GT.0.02D0) THEN
          TAUFAC=BETA/DGAMMA(1.D0/BETA)
        ELSE
          TAUFAC=1.D30
        ENDIF
        BETAOLD=BETA
      END IF
      OMR=(OM-PA(4))*1519.27D0
C     ----------------------> OMEGA IN 1/NSEC
      A=PA(1)
 
c------------------------------------------------------------------
c     a = a*(-0.17D0*QQ+1.04D0)
c
c------------------------------------------------------------------
 
      IF (PA(5).EQ.0.0D0) THEN
 
cccccccccccccccccccccccccccMARKERccccccccccccccccccccccccccccccccccccc
c       TAU=PA(2)/(QQ**(2.0D0/beta))
        TAU=PA(2) * dsqrt(qq)**(-pa(6)/beta)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      ELSEif(PA(5).lt.2) THEN
        TAU=TAUFAC*PA(2)
      else
        tau = pa(2)/ 
     *   ((1.0001d0-sin(dsqrt(qq)*pa(6))/(dsqrt(qq)*pa(6)))**(1/beta))
      END IF
C     ----------------------> TAU IN NSEC
      BETA=PA(3)
 
      IF (OMR.EQ.OMOLD) OMOLD=OMOLD2
 
      IF (OMR.GT.OMOLD) THEN
        dom = OMR-OMOLD
        SKOHL1=SKOHL(A,BETA,TAU,OMOLD,dom)
      ELSE
        dom = 1.d-3*ABS(OMR)
        SKOHL1=SKOHL(A,BETA,TAU,OMR,dom)
      END IF
 
      OMOLD2=OMOLD
      OMOLD=OMR
 
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION SKOHL(A,BETA,TAU,OMEGA,DOM)
 
C     Calculation of the average value of the scattering
C     function resulting from a Kohlrausch-type relaxation
C     in a frequency interval [OMEGA,OMEGA+DOM]
C         A      Prefactor (=2*S*F0, i.e. integral intens. =A/2)
C         BETA   Kohlrausch-exponent (0.1...1)
C         TAU    Kohlrausch-time
C         OMEGA  Lower limit of frequency interval
C         DOM    width of frequency interval (>1E-8*OMEGA)
 
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
      IMPLICIT NONE
!                                                                       !
      DOUBLE PRECISION a, beta, tau,omega, dom
      DOUBLE PRECISION ftkint, o10, o20, o1, o2
      DOUBLE PRECISION dprod
!                                                                       !
      EXTERNAL dprod
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
 
      IF (BETA.LT.0.1D0.OR.BETA.GT.1.0D0) THEN
        WRITE(*,*)
     2    'SKOHL: (Fatal) beta<=0.1 or >1, calculation aborted!'
        SKOHL=1.0
        RETURN
      END IF
      IF (DOM.LT.1D-8*ABS(OMEGA))
     2  WRITE(*,*) 'SKOHL: (Warning) DOM too small, accuracy loss!'
      IF (TAU.LE.0.0D0) THEN
        WRITE(*,*)
     2    'SKOHL: (Fatal) tau<=0, calculation aborted!'
        SKOHL=1.0D0
        RETURN
      END IF
 
      O10 = OMEGA*TAU
      O20 = O10 + DOM*TAU
      IF (O10.GT.0.D0) THEN
        O1=O10
        O2=O20
      ELSE
        IF (O20.LT.0.D0) THEN
          O1=-O20
          O2=-O10
        ELSE
          O1=O20
          O2=-O10
          GOTO 1
        END IF
      END IF
      SKOHL=0.159154943D0*A/DOM*(FTKINT(BETA,O2)-FTKINT(BETA,O1))
      RETURN
 
1     SKOHL=0.159154943D0*A/DOM*(FTKINT(BETA,O1)+FTKINT(BETA,O2))
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION FTKINT(BETA,OT)
 
C     Calculation of the 'masterfunction' of the integral from
C     0 to omega*tau (OT) of the scattering function.
C     The first call to this function causes an initialisation
C     which implies reading of two-dimensional B-spline
C     coefficients from the FILE FTKSPC. The function is
C     evaluated by interpolation using those coefficients.
 
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
      IMPLICIT NONE
!                                                                       !
      INCLUDE  'tofsys_parameter_defs.f'
      INCLUDE  'tofsys_comdefs.f'
      INCLUDE  'theory_comdefs.f'
!                                                                       !
      INTEGER*4 korder, nb, no, no1, i, init, nb1, ko1
      PARAMETER (korder = 5, nb = 40, no = 40)
!                                                                       !
      DOUBLE PRECISION beta
      DOUBLE PRECISION BSCOEF,XKNOT,YKNOT,BETA1,OT,BLOT,DBS2VL,BSC
      DOUBLE PRECISION dgamma
!                                                                       !
      DIMENSION XKNOT(NB+KORDER),YKNOT(NO+KORDER),
     *          BSCOEF(NB,NO),BSC(NB*NO)
!                                                                       !
      COMMON /SPLPAR/ BSC,XKNOT,YKNOT,INIT
!                                                                       !
      EQUIVALENCE(BSCOEF(1,1),BSC(1))
!                                                                       !
      EXTERNAL DBS2VL
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
 
      IF (INIT.NE.221134719) THEN
 
C       Initialisation code, to be run only once:
 
        WRITE(*,*) 'FTKINT: (Info) Reading Coefficient file.'
 
        OPEN(UNIT=67,FILE=kohl_data,STATUS='OLD')
        READ(67,*) KO1,NB1,NO1
        IF (KO1.NE.KORDER.OR.NB1.NE.NB.OR.NO1.NE.NO) GOTO 3
 
        DO I=1,NB+KORDER
          READ(67,*) XKNOT(I)
        END DO
 
        DO I=1,NO+KORDER
          READ(67,*) YKNOT(I)
        END DO
 
        DO I=1,NB*NO
          READ(67,*) BSC(I)
        END DO
 
        INIT=221134719
 
      END IF
 

      ftkint = 0.0d0

      IF (OT.EQ.0.D0) GOTO 4
 
      BETA1=BETA
 
      BLOT=BETA*DLOG10(OT)
      IF (BLOT.LT.-1.5D0) GOTO 1
      IF (BLOT.GT.2.8D0) GOTO 2
 
      FTKINT=DEXP(DBS2VL(BETA1,BLOT,KORDER,KORDER,
     2            XKNOT,YKNOT,NB,NO,BSCOEF))
      RETURN
 
1     FTKINT=OT/BETA1*(DGAMMA(1.D0/BETA1)
     2       -DGAMMA(3.D0/BETA1)/6.D0*OT*OT)
      RETURN
 
2     FTKINT=1.57079632679489656-DGAMMA(BETA1)*DSIN(
     2       1.57079632679489656*BETA1)*OT**(-BETA1)
      RETURN
 
3     WRITE(*,*) 'FTKINT: (Fatal) Coefficient file incompatible!'
4     FTKINT=0.D0
      RETURN
      END
