                                                      
                                                        
                                                                        
      FUNCTION th_undul_inc (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> triblock undulation incoherent approximation see EPL (Montes et al)  <------                             
!                                                                       
!
      
!     IMPLICIT NONE                                                                 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
			                                                                        
      DATA zpi / 6.283185 / 
                                                                        
       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

      DOUBLEPRECISION t, q, gam, xi, l0, eta, f, v1, v2, tau1, tau2
      double precision DERF, DDAWS, DE1
      INTEGER IFAIL
      REAL qp 
      REAL S13AAF, S15AEF
!      EXTERNAL  
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'undulinc' 
         nparx = 5
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th5 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intens'     ! Scaling Prefactor
         parnam (2) = 'gam'        ! Surface tension  
         parnam (3) = 'xi'         ! correlation length
         parnam (4) = 'l0'         ! minimal length
         parnam (5) = 'eta'        ! effective viscosity  
!                                                                       
         th_undul_inc = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      t     = x*1e-9 
      a0    = pa (1)
      gam   = pa (2) 
      xi    = pa (3) 
      l0    = pa (4) 
      eta   = pa (5)
      
      IFAIL = 0 
                                                                        
! extract parameter !                                                   
      CALL getpar ('q       ', qp,nopar ,params,napar,mbuf, ier)  
      IF (ier.ne.0) then 
         WRITE (6, * ) 'undulinc: q-parameter not found! set q=0.1' 
         qp = 0.1 
      ENDIF 
      q = qp*1e10 

      CALL getpar ('temp    ', temp,nopar ,params,napar,mbuf, ier)  
      IF (ier.ne.0) then 
         WRITE (6, * ) 'undulinc: temp-parameter not found! set temp=400' 
         temp=393 
      ENDIF 
      tau1 = (4*eta/gam)/(2*pi/xi)
      tau2 = (4*eta/gam)/(2*pi/l0)
      v1 = t/tau1        !  t/(2*pi*gam/xi/4/eta)
      v2 = t/tau2          !t/(2*pi*gam/l0/4/eta)
!      write(6,*)'1:',  tau1, tau2       


!      f = kb*temp/(4*Pi*gam)
!      f = f*(log(xi/l0)-DEI(-v1)+DEI(-v2))


!      f = kb*temp/(4*pi*gam) !! Falsch, KORR
       f = 2* kb*temp/(2*pi*gam) !! 6.8.09 und 14.6

!      f = f*(log(xi/l0) -S13AAF(v1,IFAIL) +S13AAF(v2,IFAIL))                              
      f = f*(log(xi/l0) -DE1(v1) +DE1(v2))                              
      ifail = 0
!      write(6,*)'2: ',S13AAF(v1,IFAIL), S13AAF(v2,IFAIL), sqrt(f)
      ifail = 0
                           

      f = f/2    ! da war noch ein Fehler im EPL-Paper


!      th_undul_inc = a0*0.5d0*sqrt(Pi)*S15AEF(q*sqrt(f),IFAIL)/(q*sqrt(f))
      th_undul_inc = a0*0.5d0*sqrt(Pi)*DERF(q*sqrt(f))/(q*sqrt(f))
!      write(6,*)'3: ',th_undul_inc,q*sqrt(f)                                   
!                                    
      RETURN 
      END FUNCTION th_undul_inc
!               
