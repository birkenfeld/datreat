   program tests
!  double precision :: x,y(3),de1,ddaws,derf
!
!   x    = 1d0
!   y(1) = derf(x)
!   y(2) = ddaws(x)
!   y(3) = de1(x)
!
!   write(6,*)x,y
!
!   x    = 2d0
!   y(1) = derf(x)
!   y(2) = ddaws(x)
!   y(3) = de1(x)
!
!   write(6,*)x,y
!
!   double precision :: bsjn, z, dz, y(-2:10), r1, r2
!   double precision :: radial_bsjn2_integral
!   integer          :: m,i,j
!
!   dz =0.2d0
!   do i=0,100
!     z = i*dz
!     do m=-2,10
!       y(m)=bsjn(m,z)      
!     enddo
!     write(6,'(f10.6,2x,13f10.6)')z,y
!   enddo
!
!
!   write(6,*)'-------------------------------'
!
!   dz =0.2d0
!   r1 =0.3d0
!   r2 =1.2d0
!   do i=0,100
!     z = i*dz
!     do m=0,10
!       y(m)=radial_bsjn2_integral(m,z,r1,r2)
!     enddo
!     write(6,'(f10.6,2x,11f10.6)')z,y(0:10)
!   enddo
!
    double precision :: q,t,d,a,zx,y,DINSPH_t 
    integer          :: lx
1   write(6,*) 'Q,t,D,A,LX,ZX > '
    read(5,*,end=999) Q,t,D,A,LX,ZX


!    call VDCOEF(ZX,LX) 


    y = DINSPH_t(Q,t,D,A,LX,ZX)     
    write(6,*) y
    goto 1
999 continue
   return
   end program tests
  

      FUNCTION BSJN(M,X)                                               
!      ==================                                               
!                                                                       
! ----- SPHERICAL BESSELFUNCTION OF ORDER M>=0 -------------------------
!                                                                       
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
        DOUBLE PRECISION :: bsjn, x, eps, f, fmn, fmnm1, fmnp1, fn, fnm1, fnp1, z, zzt
!                                                                       
        INTEGER          :: i, m, maxm, n
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
       DATA EPS/1.0D-6/                                        
!                          
!
       z = x
       if(z.eq.0d0) z = 1d-8
!                                                   
! --- LOOK FOR BAD M VALUE ---                                          
       IF(M.LT.-2) THEN                                                  
        WRITE(6,*)' ERROR: BSJN CALLED M =',M,' < -2 --> BSJN = 0'      
        BSJN = 0.0D0
        RETURN   
       else if(m.eq.-2) then
         BSJN = -cos(z)/(z*z)-sin(z)/z
         return
       else if(m.eq.-1) then          
         BSJN =  cos(z)/z
         return                                                      
       ENDIF                                                            
!                                                                       
       IF(Z.EQ.0.0D0) THEN 
         BSJN = 1.00D0
         IF(M.GT.0) BSJN = 0.0D0
         RETURN                                                         
       ENDIF                                                            
! --- PREPARE FACTORIAL FOR LIMITING FORMULA ---                        
       F = 1.00D0
       DO 101 I=1,2*M+1,2                                               
         F = F * I                                                      
101    CONTINUE                                                         
!                                                                       
! ---- TREAT CASE Z SMALL ----                                          
       ZZT = (Z**M) / F                                                 
       IF(ABS(ZZT).LT.EPS) THEN                                         
!        -----> USE LIMITING FORMULA                                    
           BSJN = (Z**M / F ) * (1.0D0-Z**2/(2*(2*M+3)) +     &
                                 Z**4/(8*(2*M+3)*(2*M+5)) )             
!        WRITE(6,*)' ASYMPTOTIC FORMULA'                                
         RETURN                                                         
       ENDIF                                                            
!                                                                       
! --- TREAT CASES M=0 AND M=1 ---                                       
       IF(M.LT.2) THEN                                                  
        IF(M.EQ.0) THEN                                                 
          BSJN = SIN(Z) / Z                                             
          RETURN                                                        
        ENDIF                                                           
        IF(M.EQ.1) THEN                                                 
          BSJN = ( SIN(Z) / Z  - COS(Z) ) / Z                           
          RETURN                                                        
        ENDIF                                                           
      ENDIF                                                             
!                                                                       
! --- TREAT CASES M >= 2 ---                                            
!                                                                       
      FNM1 = 1.00D0/Z   
      FN   = FNM1**2                                                    
      FMN  = 0.0D0
      FMNP1= FNM1                                                       
      MAXM = M - 1                                                      
!                                                                       
      DO 1 N=1,MAXM                                                     
        FNP1   = (2*N+1) * FN/Z - FNM1                                  
        FMNM1  = (1-2*N) * FMN/Z- FMNP1                                 
        FNM1   = FN                                                     
        FN     = FNP1                                                   
        FMNP1  = FMN                                                    
        FMN    = FMNM1                                                  
1     CONTINUE                                                          
      FMNM1 = (1-2*M) * FMN/Z - FMNP1                                   
      BSJN  = FN*SIN(Z) - (-1)**M * FMNM1*COS(Z)                        
!                                                                       
      RETURN                                                            
      END  function bsjn                                                



      function radial_bsjn2_integral(n,q,r1,r2)
!     -----------------------------------------
!
!     Integral over a shell from r1..r2 of bsjn(n,q*r)**2  as is needed in the rotational diffusion
!     modelling for a rotating protein in incoherent approximation
!     The integral is normalized to the shell volume
!
!     i.e. integral_[r1..r2](bsjn(x,q*r)**2 *r**2 * dr) / integral_[r1..r2](r**2 *dr) 
!
      implicit none
      double precision              :: radial_bsjn2_integral
      double precision, intent(in)  :: q, r1, r2
      integer         , intent(in)  :: n
    
      double precision              :: bsjn

      radial_bsjn2_integral = -(0.5d0)*R1**3*bsjn(n  , R1*q)**2                 &
                              +(0.5d0)*R1**3*bsjn(n-1, R1*q)  *bsjn(n+1, R1*q)  &
                              +(0.5d0)*R2**3*bsjn(n  , R2*q)**2                 &
                              -(0.5d0)*R2**3*bsjn(n-1, R2*q)  *bsjn(n+1, R2*q) 


!     normalizing to volume
      
      radial_bsjn2_integral =  radial_bsjn2_integral*3/(r2**3-r1**3)

      return
      end function  radial_bsjn2_integral





!                                                                       
        DOUBLE PRECISION FUNCTION DISPH(QQ,OM,TEMP,PA) 
!       ==============================================
!                                                                       
! ---- ADAPTION FOR DIFFUSION IN SPHERE (DINSPH) -----                  
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!        INCLUDE  'natur_konstanten_defs.f'
!                                                                       
        DOUBLE PRECISION a, d, o, q, zx
        DOUBLE PRECISION dinsph
!                                                                       
        DOUBLE PRECISION om, qq, temp
        DOUBLE PRECISION pa
!                                                                       
        INTEGER*4        lx
!                                                                       
        DIMENSION        PA(10)  
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
!                                                                       
         O = OM * 1.51925D12                                             
!        O = OM * om2rad_sec                      
!       -------------------> OMEGA IN SEC**-1                           
        Q   = SQRT(QQ)                                                  
        D   = PA(2) * 1.0D16                                            
        A   = PA(3)                                                     
        LX  = NINT(PA(4))
        ZX  = PA(5)                                                     
!                                                                       
        DISPH = DINSPH(Q,O,D,A,LX,ZX) * PA(1)                           
!                                                                       
        RETURN                                                          
        END                                                             
!                                                                       
!                                                                       
!                                                                       
! ---- SCATTERING LAW OF A DIFFUSION INSIDE A SPHERE ----               
!      ACCORDING TO F.VOLINO & A.J. DIANOUX MOL.PHYS (1980),VOL41,271   
!                                                                       
       FUNCTION DINSPH(Q,OM,D,A,LX,ZX)                                  
!      -------------------------------                                  
!                                                                       
!   Q   = MOMENTUM TRANSFER                                             
!   OM  = FREQUENCY (ENERGY TRANSFER)                                   
!   D   = DIFFUSION CONSTANT                                            
!   A   = RADIUS OF SPHERE                                              
!   LX  = MAX L TO BE USED                                              
!   ZX  = MAX VALUE OF ZERO TO BE USED                                  
!                                                                       
!   ZEROS(I) = VALUE OF I-TH ZERO = XNL                                 
!   LOFZ (I) = L-VALUE FOR THIS ZERO                                    
!   NOFZ (I) = NO OF THIS ZERO                                          
!   NZERO    = TOTAL NO OF ZEROS                                        
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!        INCLUDE  'th1_vdzero_comdef.f'

        DOUBLE PRECISION        zeros
        INTEGER*4               lofz, nofz, nzero, maxz
        PARAMETER               (maxz = 5000)
        DIMENSION               zeros(maxz), lofz(maxz), nofz(maxz)
!                                                                           !
        COMMON /vdzero/         zeros, lofz, nofz, nzero

        DOUBLE PRECISION        zmx
        INTEGER*4               lmx
!                                                                           !
        COMMON /dinbuf/         zmx, lmx


!        INCLUDE  'th1_dinbuf_comdef.f'
!        INCLUDE 'natur_konstanten_defs.f'
!                                                                       
        DOUBLE PRECISION a, aq, d, da, delta, dinsph, dom, om, oml, &
                         q, qa, sum, x, xx, zx
        DOUBLE PRECISION bsjn
!                                                                       
        INTEGER*4        i, l, lx, n
!                                                                       
        DATA             OML/1.0D50/
        double precision, parameter :: pi=3.141592654d0

!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
       DELTA = 1.0D-6
!      -------------> LIMIT FOR AVOIDING SINGULARITY                    
! ---- OML IS NEEDED TO IDENTIFY THE ZERO POSITION & GET THE FACTOR FOR 
!      THE DELTA FUNCTION                                               
       IF(OML.GT.OM) OML = OM                                           
!                    --------> RESETS OML IF A NEW DETECTOR IS TAKEN    
!                                                                       
! ---- GENERATE ZEROS IF NECESSARY ----                                 
!                                                                       
       IF(LX.NE.LMX.OR.DABS(ZX-ZMX).GT.0.2D0.OR.NZERO.EQ.0) THEN        
         LMX = LX                                                       
         ZMX = ZX                                                       
         CALL VDCOEF(ZMX,LMX)                                           
         WRITE(6,*)' .....',NZERO,' NEW VDCOEFFS GENERATED'             
       ENDIF                                                            
!                                                                       
! ---- DO THE SUMMATION ----                                            
!                                                                       
       QA = Q*A                                                         
       DA = D / A**2                                                    
       SUM= 0.0D0
!                                                                       
       DO 10 I=1,NZERO                                                  
        N = NOFZ(I)                                                     
        L = LOFZ(I) 
        if(l.eq.0 .and.n.eq.0) cycle                                               
        X = ZEROS(I)                                                    
        XX= X**2                                                        
        IF(DABS(QA-ZEROS(I)).GT.DELTA) THEN                             
          AQ = ( QA * BSJN(L+1,QA) - L * BSJN(L,QA) ) / (QA**2-XX)
          AQ = AQ**2 * 6.0D0 * XX / ( XX - L*(L+1) )
        ELSE                                                            
          AQ = 1.5D0 * BSJN(L,X)**2 * (XX - L*(L+1) ) / XX 
        ENDIF                                                           
!                                                                       
        SUM = SUM + (2*L+1)*AQ *    XX*DA / ( (XX*DA)**2 + OM**2 )       ! hier Zeitfunktion einsetzen  
!                                                                       
10     CONTINUE                                                         
       SUM = SUM / PI                                                   
!                                                                       
! ---- ADD THE DELTA FUNCTION ----                                      
!                                                                       
       IF(OM*OML.LT.0.0D0.OR.OM.EQ.0.0D0) THEN 
        DOM = DABS(OM-OML)                                              
        SUM = SUM + (3.0D0*BSJN(1,QA)/QA)**2 / DOM 
       ENDIF                                                            
       DINSPH = SUM                                                     
!      ------------                                                     
!                                                                       
       RETURN                                                           
       END                                                              


!                                                                       
!                                                                       
!                                                                       
! ---- SCATTERING LAW OF A DIFFUSION INSIDE A SPHERE ----    
!      time domain: intermediate scattering function           
!      ACCORDING TO F.VOLINO & A.J. DIANOUX MOL.PHYS (1980),VOL41,271   
!                                                                       
       FUNCTION DINSPH_t(Q,t,D,A,LX,ZX)                                  
!      --------------------------------                                  
!                                                                       
!   Q   = MOMENTUM TRANSFER                                             
!   OM  = FREQUENCY (ENERGY TRANSFER)                                   
!   D   = DIFFUSION CONSTANT                                            
!   A   = RADIUS OF SPHERE                                              
!   LX  = MAX L TO BE USED                                              
!   ZX  = MAX VALUE OF ZERO TO BE USED                                  
!                                                                       
!   ZEROS(I) = VALUE OF I-TH ZERO = XNL                                 
!   LOFZ (I) = L-VALUE FOR THIS ZERO                                    
!   NOFZ (I) = NO OF THIS ZERO                                          
!   NZERO    = TOTAL NO OF ZEROS                                        
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!        INCLUDE  'th1_vdzero_comdef.f'

        DOUBLE PRECISION        zeros
        INTEGER*4               lofz, nofz, nzero, maxz
        PARAMETER               (maxz = 5000)
        DIMENSION               zeros(maxz), lofz(maxz), nofz(maxz)
!                                                                           !
        COMMON /vdzero/         zeros, lofz, nofz, nzero

        DOUBLE PRECISION        zmx
        INTEGER*4               lmx
!                                                                           !
        COMMON /dinbuf/         zmx, lmx


!        INCLUDE  'th1_dinbuf_comdef.f'
!        INCLUDE 'natur_konstanten_defs.f'
!                                                                       
        DOUBLE PRECISION a, aq, d, da, delta, dinsph_t, dom, om, oml, &
                         q, qa, sum, x, xx, zx, t
        DOUBLE PRECISION bsjn

        double precision, parameter :: pi=3.141592654d0
!                                                                       
        INTEGER*4        i, l, lx, n
!                                                                       
!        DATA             OML/1.0D50/
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
       DELTA = 1.0D-6
!      -------------> LIMIT FOR AVOIDING SINGULARITY                    
! ---- OML IS NEEDED TO IDENTIFY THE ZERO POSITION & GET THE FACTOR FOR 
!      THE DELTA FUNCTION                                               
!       IF(OML.GT.OM) OML = OM                                           
!                    --------> RESETS OML IF A NEW DETECTOR IS TAKEN    
!                                                                       
! ---- GENERATE ZEROS IF NECESSARY ----                                 
!                                                                       
       IF(LX.NE.LMX.OR.DABS(ZX-ZMX).GT.0.2D0.OR.NZERO.EQ.0) THEN        
         LMX = LX                                                       
         ZMX = ZX                                                       
         CALL VDCOEF(ZMX,LMX)                                           
         WRITE(6,*)' .....',NZERO,' NEW VDCOEFFS GENERATED'             
       ENDIF                                                            
!                                                                       
! ---- DO THE SUMMATION ----                                            
!                                                                       
       QA = Q*A                                                         
       DA = D / A**2                                                    
       SUM= 0.0D0
!                                                                       
       DO 10 I=1,NZERO                                                  
        N = NOFZ(I)                                                     
        L = LOFZ(I)  
        if(l.eq.0 .and.n.eq.0) cycle                                                   
        X = ZEROS(I)                                                    
        XX= X**2                                                        
        IF(DABS(QA-ZEROS(I)).GT.DELTA) THEN                             
          AQ = ( QA * BSJN(L+1,QA) - L * BSJN(L,QA) ) / (QA**2-XX)
          AQ = AQ**2 * 6.0D0 * XX / ( XX - L*(L+1) )
        ELSE                                                            
          AQ = 1.5D0 * BSJN(L,X)**2 * (XX - L*(L+1) ) / XX 
        ENDIF                                                           
!                                                                       
!        SUM = SUM + (2*L+1)*AQ *    XX*DA / ( (XX*DA)**2 + OM**2 )       ! hier Zeitfunktion einsetzen  
        SUM = SUM + (2*L+1)*AQ *    Pi*exp(-XX*DA*abs(t))
!
!   write(6,*)'1: ',i,n,l,qa,aq,xx,sum                                                                       
10     CONTINUE                                                         
       SUM = SUM / PI  
!   write(6,*)'1a: ',sum                                                 
!                                                                       
! ---- ADD THE DELTA FUNCTION ----                                      
!                                                                       
!       IF(OM*OML.LT.0.0D0.OR.OM.EQ.0.0D0) THEN 
!        DOM = DABS(OM-OML)                                              
!        SUM = SUM + (3.0D0*BSJN(1,QA)/QA)**2 / DOM 
!       ENDIF                                                            
!
        SUM = SUM + (3.0D0*BSJN(1,QA)/QA)**2 
!    write(6,*)'2: ',qa,aq,xx,sum   
       DINSPH_t = SUM                           
!      --------------  
!   write(6,*)'2a: ',sum                                  
       RETURN                                
       END                                  



       DOUBLE PRECISION FUNCTION BLJ(Z,L)  
!      ==================================
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
        DOUBLE PRECISION z 
        DOUBLE PRECISION bsjn
!                                                                       
        INTEGER*4       l
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
       IF(L.EQ.0) THEN                                                  
         BLJ = BSJN(1,Z)  
         RETURN                                                         
       ENDIF                                                            
!                                                                       
       BLJ = L*BSJN(L,Z) - Z * BSJN(L+1,Z)    
!      -------------------------------------
       RETURN                                                           
       END                                                              
!                                                                       
!                                                                       
!                                                                       

!                                                                       
! ---- DETERMINE COEFFICIENTS XNL FOR VOLINO DIANOUX THEORY FOR         
!      THE SCATTERING OF DIFFUSION IN A SPHERICAL POTENTIAL WELL ----   
!                                                                       
       SUBROUTINE VDCOEF(ZMX,LMX)                                       
!      --------------------------                                       
!                                                                       
! ---- ZMX = UPPER LIMIT FOR ZEROES ----                                
!      LMX = UPPER LIMIT FOR L                                          
!                                                                       
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!       INCLUDE  'outlev_comdef.f'
!       INCLUDE  'th1_vdzero_comdef.f'
        DOUBLE PRECISION        zeros
        INTEGER*4               lofz, nofz, nzero, maxz
        PARAMETER               (maxz = 5000)
        DIMENSION               zeros(maxz), lofz(maxz), nofz(maxz)
!                                                                           !
        COMMON /vdzero/         zeros, lofz, nofz, nzero
!                                                                       
        DOUBLE PRECISION dz, dz0, eps, f, f0, f1, z0, zmx
        DOUBLE PRECISION blj
!                                                                       
        INTEGER*4        i, j, l, lmx, maxit, ncnt, nl
        integer  :: ierrs
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
!                                                                       
! ---- SET PARAMETER-VALUES ---                                         
       MAXIT = 50                                                       
!      -----------> MAX NO. OF NEWTON ITERATIONS                        
       DZ    = 0.001D0
!      -----------> SETP FOR NEWTON DIFFERENTIATION                     
       EPS   = 1.0D-9  
!      -----------> REQUIRED ACCURACY FOR NEWTON ITERATION              
       DZ0   = 0.2D0
!      -----------> INCREMENT FOR EXHAUSTIVE ZERO SEARCH                
!                                                                       
       WRITE(6,*)' VDCOEF GENERATING COEFFICIENTS FOR ',LMX,ZMX,' ...'  
!                                                                       
       NL = ZMX / DZ0                                                   
       NZERO = 0                                                        
!                                                                       
       DO 5 L = 0,LMX                                                   
         NCNT = 0                                                       
         DO 10 J = 1,NL                                                 
          Z0 = (J-1) * DZ0                                              
          F  = BLJ(Z0,L) * BLJ(Z0-DZ0,L)                                
          IF(F.LT.0) THEN                                               
!                                                                       
! ---- NEWTON-ITERATION ----                                            
            DO 100 I=1,MAXIT                                            
             F0 = BLJ(Z0,L)                                             
             F1 = BLJ(Z0+DZ,L)                                          
             Z0 = Z0 - F0 * DZ / (F1-F0)                                
!            WRITE(6,*)I,F0,Z0                                          
             IF(DABS(F0).LT.EPS) GOTO 101                               
100         CONTINUE                                                    
            WRITE(6,*)' WARNING BAD CONVERGENCE , RESIDUAL=',F0         
            IERRS = 55                                                  
!                                                                       
101         CONTINUE                                                    
            NZERO = NZERO + 1                                           
            IF(NZERO.GT.MAXZ) THEN                                      
              WRITE(6,*)' NO OF ZEROS EXCEEDS MAX DIMENSION = ',MAXZ    
              IERRS = 777                                               
              NZERO = MAXZ                                              
            ENDIF                                                       
!                                                                       
            ZEROS(NZERO) = Z0                                           
            LOFZ(NZERO)  = L                                            
            NOFZ(NZERO)  = NCNT                                         
!                                                                       
!       WRITE(6,*)L,NCNT, Z0                          
            NCNT = NCNT + 1                                             
          ENDIF                                                         
10       CONTINUE                                                       
5      CONTINUE                                                         
!                                                                       
       RETURN                                                           
!                                                                       
       END                                                              

