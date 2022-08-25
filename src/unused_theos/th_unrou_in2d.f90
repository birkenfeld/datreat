                                                      
                                                        
                                                                        
      FUNCTION th_unrou_inc3(x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> triblock undulation incoherent approximation see EPL (Montes et al)  <------    
!          + orthogonal 2D Rouse Segment motion   ... Factor 1/2 u**2 check older Version                         
!                                                                       
!
      
      IMPLICIT NONE
      
      real th_unrou_inc3, pa, zpi, qq, x
      integer npar, ini, nparx
                                                                 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
			      

      integer ier
      double precision a0, a01, sq 
      real   temp

      double precision s13aaf, s15aff, s15aef

      double precision qexp, qoffs, const, frac
                       

                                           
      DATA zpi / 6.283185 / 
                                                                        
       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

      DOUBLE PRECISION t, q, gam, xi, l0, eta, f, v1, v2, tau1, tau2
      double precision uuz, uux, wl4, Dcm
      INTEGER IFAIL, iadda
       common/thiadd/iadda
 
      REAL qp 
!      REAL S13AAF, S15AEF
!      EXTERNAL  
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'unrouin3' 
         nparx = 10
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_unrou_inc3 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intens'     ! Scaling Prefactor
         parnam (2) = 'gam'        ! Surface tension     [SI]
         parnam (3) = 'xi'         ! correlation length  [SI]
         parnam (4) = 'l0'         ! minimal length      [SI]
         parnam (5) = 'eta'        ! effective viscosity [SI]  
         parnam (6) = 'wl4'        ! in A*4/ns  !!!!!!   [A,ns]
         parnam (7) = 'const'      ! Anteil Konstante (Multiplikator)
         parnam (8) = 'qexp'       ! Qexp des Konstante Anteils
         parnam (9) = 'qoffs'      ! qoffset
         parnam (10) =  'Dcm'      ! Centre of mass diff of whole molecule (DB)
!
         th_unrou_inc3 = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      t     = x*1e-9 
      a01   = pa (1)
      a0    = 1
      gam   = pa (2) 
      xi    = pa (3) 
      l0    = pa (4) 
      eta   = pa (5)
      
      wl4   = pa (6)

      const = pa(7)
      qexp  = pa(8)
      qoffs = pa(9)
      Dcm  =  pa(10)*1d-9 / 1d-16 ! in A2/ns

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
       f = 2* kb*temp/(2*pi*gam) !! 6.8.09 und 14.6
  
      f = f*(log(xi/l0) -S13AAF(v1,IFAIL) +S13AAF(v2,IFAIL))                              
      ifail = 0
!      write(6,*)'2: ',S13AAF(v1,IFAIL), S13AAF(v2,IFAIL), sqrt(f)
      ifail = 0


      uuz = f      ! --> u**2


! Rousebeitrag (in eine Richtung daher *1/3) aber 2D in der Ebene

      uux = 2.0d0/3.0d0*sqrt(x*wl4/Pi)*1d-20 + 4*Dcm*x*1d-20    !! --> damit ist uux in m**2


      write(6,*)uuz,uux

! Winkelmittelung

      ifail = 0

      if(abs((uuz-uux)*q*q).lt.1d-6) then
        th_unrou_inc3 = a0*exp(-q*q*uuz/2.0d0)
        goto 200
      endif

      if(uuz.gt.uux) then
        sq = q*sqrt(uuz/2-uux/2)
        write(6,*)'sq=',sq
        th_unrou_inc3 = a0*S15AEF(sq,ifail)*sqrt(Pi)*exp(-q*q*uux/2d0)/(2*sq)
!                         Erf
        goto 200
      endif

      if(uuz.lt.uux) then
        sq = q*sqrt(-uuz/2+uux/2)
        th_unrou_inc3 = a0* S15AFF(sq,ifail)*2*exp(sq**2-q*q*uux/2d0)/(2*sq)
!                          Dawson, wegen erf(i*x)=2*i*exp(x**2)*Dawson(x)/sqrt(Pi)
        goto 200
      endif

200   continue


      frac = (q+qoffs)**qexp/(const+(q+qoffs)**qexp)
                           
      th_unrou_inc3 = frac*th_unrou_inc3 +(1-frac)

      th_unrou_inc3 = th_unrou_inc3 * a01

      call  parset('frac    ',sngl(frac),iadda,ier) 
 
!                                    
      RETURN 
      END FUNCTION th_unrou_inc3
!               
