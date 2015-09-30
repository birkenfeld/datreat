c*ds
c*ds
      function th17(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> fatkim <--------
! 
! Fatkullin & Kimmich Expression versus Doi-Edwards Q**0.25 regime (rms--> Sinc(Q,t))
! short time Rouse rms (->gaussian->)                                  --> Sinc(Q,t)
! soft crossover function
! ability to modify by extra term (deGennes plateau): exp(-Q**2*d**2/36)
! width d = a*rw  (a  =parameter: step length of random walk of contorted tube)
!                  rw =modification factor: 0=off (original formulation))
! other parameters:     an  = normalization factor should always be 1
!                       wl4 = Rouse W*l^4 (in A ns units)
!                       a   = step length (= tube diameter)
!                       nu  = exponent modificator in Fatkullin-Kimmich expr. , should be 1
!                       fde = 1 selects Doi-Edwards expression for Q**1/4
!                       ffk = 1 selects Fatkullin Kimmich
!                       tr  = smooth--->sharp transition ( about 10)
c
c
       
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
       
       double precision aint
       double precision a, Wl4, Q, t
       double precision soft_crossover, rms_tr, AFK, RSINC, DESINC
       double precision spar, addexp, rw
       double precision dq, dt, y1,y2,y3,y4,rms1,rms2,rms3,rms4
       double precision t0, tf
       double precision deGfac, y23, f2, f3
 
       common/thiadd/iadda
!
! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'fatkim'
         nparx = 8
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th17  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'an'
         parnam(2) = 'wl4'
         parnam(3) = 'a'
         parnam(4) = 'nu'
         parnam(5) = 'tr'
         parnam(6) = 'rw'   ! extra tube width
         parnam(7) = 'fde'  ! Doi-Edwards rms-faktor
         parnam(8) = 'ffk'  ! Faktkullin Kimmich-faktor
c
         th17 = 0
         return
       endif
!
! ---- calculate theory here -----
! Unit System: A, ns

       t = x

       aint      = pa(1)
       wl4       = pa(2)  ! stdvalue Pe(509K)=70000
       a         = pa(3)  ! stdvalue Pe      =46
       addexp    = pa(4)  ! stdvalue 1
       spar      = pa(5)  ! stdvalue 10    sharpness of transition 
       rw        = pa(6)  ! stdvalue rw tube diameter factor
       fde       = pa(7)  ! selects Doi& Edwards
       ffk       = pa(8)  ! selects Fatkullin&Kimmich
 
       Q         = dparam('q       ')
       if(Q.eq.0.0d0) then
        write(6,*)'fatkim: parameter q is missing!'
        Q = 0.1  
       endif                         

       deGfac = exp(-Q**2*(a*rw)**2/36.0d0)

       y1   = RSINC(t,Q,Wl4)
       y2   = DESINC(t,Q,a,Wl4)
       y3   = AFK(t,Q,a,Wl4,addexp)

       y23  = (fde*y2+ffk*y3)*deGfac 
       
       write(*,*) t,(y23-y1)**2

       y4   = soft_crossover(y1,y23,spar)
    
       th17 = y4

! this would yield alternative versions...
!       rms1 = rms_tr(y1,Q)
!       rms2 = rms_tr(y2,Q)
!       rms3 = rms_tr(y3,Q)
!       rms4 = rms_tr(y4,Q)


!       th17 = rms4
!       th17 = soft_crossover(rms1,rms3,-spar)


       th17 = aint * th17
 
!
       return
       end


!         implicit none
!
!         double precision a, Wl4, Q, t
!         double precision soft_crossover, rms_tr, AFK, RSINC, DESINC
!         double precision spar, addexp
!         double precision dq, dt, y1,y2,y3,y4,rms1,rms2,rms3,rms4
!         double precision t0, tf
!         integer i
!
!         spar = 10.0d0
!
!
!         Wl4    = 7.0d4
!         a      = 46.0d0
!         addexp = 1.0d0
!
!         Q   = 0.15d0
!         tf  = 1.1
!         t   = 0.01d0
!
!         open(10,file='fk.dgli',status='unknown')
!         do i=0,100
!           t    = t*tf
!           y1   = RSINC(t,Q,Wl4)
!           y2   = DESINC(t,Q,a,Wl4)
!           y3   = AFK(t,Q,a,Wl4,addexp)
!           y4   = soft_crossover(y1,y3,spar)
!           rms1 = rms_tr(y1,Q)
!           rms2 = rms_tr(y2,Q)
!           rms3 = rms_tr(y3,Q)
!           rms4 = rms_tr(y4,Q)
!
!           write(6,'(f10.4,3x,4f8.4,5x,4f10.2)')t,y1,y2,y3,y4,
!     *              rms1,rms2,rms3,rms4
!           write(10,'(f10.4,3x,4f8.4,5x,4f10.2)')t,y1,y2,y3,y4,
!     *              rms1,rms2,rms3,rms4
!
!         enddo
!         close(10)
!
!         end
!



         function soft_crossover(fa,fb,spar)
!        ----------------------------------> soft crossover function
!                                            Schaerfe durch spar
         implicit none
         double precision soft_crossover, fa,fb,spar
         integer ifail
         double precision S15AEF         ! erf(x,ifail)
         double precision h
 
         ifail = 1
         h = (1.0d0 + S15AEF(spar*(fb-fa),ifail))/2.0d0
         soft_crossover = fb*h+fa*(1.0d0-h)

         return
         end



         function rms_tr(y,Q)
!        -------------------- Transformation zu mean squared displacements

         implicit none
         double precision rms_tr, y, Q

         rms_tr = -6*log(y)/Q**2

         return
         end




! Fatkullin-Kimmich-Form von Sinc im local Reptation Regime !
         function AFK(t,Q,a,Wl4,addexp)

         implicit none
         double precision Pi
         parameter (Pi=3.141592654d0)

         double precision Wl4,t,PhiRouse,SS,a,AFK,Q,addexp
         integer          ifail 
         double precision S15ADF     ! Nag's erfc(x,ifail)
         double precision arg1, arg2

         PhiRouse = 2*sqrt(Wl4*t/Pi)
         SS       = PhiRouse/3.0d0
         SS       = SS**addexp
         arg1     = Q**4*a**2*SS/72.0d0
         arg2     = Q**2*a*sqrt(SS)/6.0d0/sqrt(2.0d0)
         ifail    = 1
         AFK      = exp(arg1) * S15ADF(arg2,ifail)

         return
         end


! Rouse-Form von Sinc !
         function RSINC(t,Q,Wl4)

         implicit none
         double precision Pi
         parameter (Pi=3.141592654d0)

         double precision Wl4,t,PhiRouse,RSINC,Q
         double precision arg1

         PhiRouse = 2*sqrt(Wl4*t/Pi)
         arg1     = -Q**2 * PhiRouse/6.0d0
         RSINC    = exp(arg1) 

         return
         end

! DOI-EDWARDS for Sinc im local reptation regime!
         function DESINC(t,Q,a,Wl4)

         implicit none
         double precision Pi
         parameter (Pi=3.141592654d0)

         double precision a,Wl4,t,PhiRouse,DESINC,Q
         double precision arg1, PhiDE

         PhiRouse = 2*sqrt(Wl4*t/Pi)
         PhiDE    = a * sqrt(PhiRouse/3.0d0)
         arg1     = -Q**2 * PhiDE/6.0d0
         DESINC   = exp(arg1) 

         return
         end

