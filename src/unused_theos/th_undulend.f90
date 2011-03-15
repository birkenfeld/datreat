                                                      
                                                        
                                                                        
      FUNCTION undulend2 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> triblock undulation incoherent approximation 
!          end modifying factor
!          +                          
!                                                                       
!
      
      IMPLICIT NONE
      
      real undulend2, pa, zpi, qq, x
      integer npar, ini, nparx
                                                                 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
			      

      integer ier
      double precision a0, sq , qf2, sigma
      real   temp

      double precision s13aaf, s15aff, s15aef

      double precision qexp, qoffs, const, frac
                       

                                           
      DATA zpi / 6.283185 / 
                                                                        
       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

      DOUBLE PRECISION t, q, gam, xi, l0, eta, f, v1, v2, tau1, tau2
      double precision uuz, uux, wl4
      INTEGER IFAIL, iadda
       common/thiadd/iadda
 
      REAL qp 
!      REAL S13AAF, S15AEF
!      EXTERNAL  
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'undulen2' 
         nparx = 7
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            undulend2 = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intens'     ! Scaling Prefactor
         parnam (2) = 'gam'        ! Surface tension     [SI]
         parnam (3) = 'xi'         ! correlation length  [SI]
         parnam (4) = 'l0'         ! minimal length      [SI]
         parnam (5) = 'eta'        ! effective viscosity [SI]  
         parnam (6) = 'sigma'      ! influence length in fraction of chainlength
         parnam (7) = 'wl4'        ! in A*4/ns  !!!!!!   [A,ns]
!                                                                       
         undulend2  = 0
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
      
      sigma = pa (6)
      wl4   = pa (7)


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

! Undulationsbeitrag:

      tau1 = (4*eta/gam)/(2*pi/xi)
      tau2 = (4*eta/gam)/(2*pi/l0)
      v1 = t/tau1        !  t/(2*pi*gam/xi/4/eta)
      v2 = t/tau2          !t/(2*pi*gam/l0/4/eta)
!      write(6,*)'1:',  tau1, tau2       


!      f = kb*temp/(4*Pi*gam)
!      f = f*(log(xi/l0)-DEI(-v1)+DEI(-v2))

!      f = kb*temp/(4*pi*gam) !! Falsch, KORR
       f =  2 * kb*temp/(2*pi*gam) !! 6.8.09+14.6.

      f = f*(log(xi/l0) -S13AAF(v1,IFAIL) +S13AAF(v2,IFAIL))                              
      ifail = 0
!      write(6,*)'2: ',S13AAF(v1,IFAIL), S13AAF(v2,IFAIL), sqrt(f)
      ifail = 0


      uuz = f      ! --> u**2


! Rousebeitrag (in eine Richtung daher *1/3) aber 2D in der Ebene

      uux = 2.0d0/3.0d0*sqrt(x*wl4/Pi)*1d-20     !! --> damit ist uux in m**2


!!      qf2 = (uux+uuz)*q*q/2 --> inkonsistent da nur 1 der 2D richtungen uux mus 2 mal kommen ...
      qf2 = (2*uux+uuz)*q*q/2


      ifail = 0
   

      undulend2 = -sigma*S13AAF(qf2,IFAIL)+sigma*S13AAF(qf2*exp(-1.0d0/sigma),IFAIL)
      undulend2 = a0* undulend2

 
!                                    
      RETURN 
      END FUNCTION undulend2
!               
