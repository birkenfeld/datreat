      function th17(x,pa,thnam,parnam,npar,ini)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ========================================
!
!         -------> degennes <--------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       geaendert nach likhtman - mz !!!


         implicit none

         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)


         CHARACTER*8     thnam,parnam(20)
         REAL            th17, x, pa, qq, zpi, xh, vol_frac,sel
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision w,d,q,n,ne,likhtman,t,a,l,taud,tau0,td2
         double precision Hs,lab
!test
         double precision adw,dw
 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'likhtman '
!         nparx =  8
!test
         nparx =  9
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th17  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitud'          
         parnam(2) = 'w       '
         parnam(3) = 'dtube   '
         parnam(4) = 'n       '
         parnam(5) = 'l       '
         parnam(6) = 'td2     '
         parnam(7) = 'label   '
         parnam(8) = 'sel     '
!         parnam(9) = 'adw     '
         parnam(9) = 'cnu     '
c
         th17 = 0
         return
       endif
c
c ---- calculate theory here -----
        a       = pa(1)   ! amplitude (should be 1)
        w       = pa(2)   ! Rouse rate
        d       = pa(3)   ! Tube diameter
        n       = pa(4)   ! no of segments of polymer
        l       = pa(5)   ! segment length 
        td2     = pa(6)   ! umfunktioniert zur unterscheidung theorie
        lab     = pa(7)   ! label length (0-0.5)
        sel     = pa(8)   ! select label modification
        adw     = pa(9)   ! debye-waller

        call parget('q       ',xh ,iadda,ier)
        if(ier.eq.0) then
          q       = xh
        else
          q       = 1.0  
        endif
 
        t = x 
 
        th17 = likhtman(t,q,d,W,n,l,ne,taud,tau0,td2,lab,sel,adw)

        th17 = th17 * a

        call parset('ne      ',sngl(ne),iadda)
        call parset('taud    ',sngl(taud),iadda)
        call parset('tau0    ',sngl(tau0),iadda)

        return
        end 


       double precision function likhtman(t,q,d,W,n,l,ne,taud,tau0,
     *     td2,lab,sel,adw)
!      -------------------------------------------------------------
!
! siehe P.Schleger et. al. PRL 81, p.124 (1998)
! creep-Term korrigiert nach Doi und Edwards !!
!
       implicit none
       double precision t,q,d,W
       double precision derfc, pi,arg1
       double precision tau0, taud, t0, td,td2
       double precision n, ne, l, y, sum, eqd,deb,argd,npl

       double precision mue,alphap,sumt0,alp
       double precision diffe
c mz
       double precision taue, num, st, sumlik, sumlik0
       double precision x,tx
       double precision Hs,lab
       real sel

!test
       double precision adw,dw
       double precision cr,cnu
c mz

       integer m,i,j

       parameter (pi=3.141592654d0)

       ne   =   (d/l)**2
       tau0 =  36.d0/(W*((q*l)**4))
       taud = 3.d0*(n**3)*(l**2)/( (pi**2)*W*(d**2)) 
c mz
       taue = ne**2/((pi**2)*W)


C setze taud als Fitparamter oder nicht (td2=0)
c td2 umfunktioniert

c       if(td2.eq.0.d0) then
        td = t/taud
c       else
c        td = t/abs(td2)
c       endif

       t0 = t/tau0
       eqd = exp(-((q*d/6.d0)**2))

       m=n
   
       sum = 0.d0
       sumt0=0.d0

C ********
C        npl=10.d0
        npl=n
C *******

        mue=q**2*npl*l**2/12.d0
c mz
        num = n/ne
        st = 0.5d0*((1.5d0/num)*(t/taue)**0.25)
!test
!        cnu=adw
!        st = 0.5d0*(cnu*(1.5d0/num)*(t/taue)**0.25)

! test CR
!        cnu=adw
!        cr=1.8/num*(cnu*t/taue)**0.25
!        st = 0.5d0*((1.5d0/num)*(t/taue)**0.25)*cr
!
        if(sel.eq.0) then ! full chain
           sumlik =(n/(2*mue**2))*(2.*mue+dexp(-2*mue)+2.d0-4.*mue*st-4*
     +          dexp(-2*mue*st)+dexp(-4*mue*st))
           sumlik0 = (n/(2*mue**2))*(2*mue+dexp(-2*mue)-1.d0)
           sumlik = sumlik/sumlik0
        elseif(sel.eq.1) then !CENTER LABEL
c  probleme bei grossen n!
c  ausserdem Hs(-1+lab+st)=0 fuer lab und st <=0.5!
C           sumlik = (2.d0+2.*mue-4.*st*mue+dexp(-2*mue*(1-2*lab))-4.*
C     #          dexp(-2*(st-lab)*mue)+dexp(-4*(st-lab)*mue)+Hs(-st+lab)*
C     #          (4.*dexp(-2*(st-lab)*mue)-dexp(-4*(st-lab)*mue)+4.*
C     #          (st-lab)*mue-3.d0)+Hs(-1+lab+st)*(dexp(-4*(-1+lab+st)*
C     #          mue)+3.*dexp(-2*(st-lab)*mue)-2.*dexp(-2*mue*(2*st-1))-
C     #          1.d0+4.*mue*(-1+2.*lab)-2.*dexp(2*mue*(-1+2*lab))+
C     #          dexp(2*(st-lab)*mue)))/(2.*mue**2)
c  nur umsortiert!
           sumlik = (2.d0+2.*mue-4.*st*mue+dexp(-2*mue*(1-2*lab))-
     #          (1-Hs(-st+lab))*4.*dexp(-2*(st-lab)*mue)+
     #          (1-Hs(-st+lab))*dexp(-4*(st-lab)*mue)+
     #          Hs(-st+lab)*(4.*(st-lab)*mue-3.d0))/(2.*mue**2)

           sumlik0 = (-1.d0+2.*mue*(1-2.*lab)+dexp(-2*mue*(1-2*lab)))
     #          /(2.*mue**2)
           sumlik = sumlik/sumlik0

        elseif(sel.eq.2) then  !END LABEL
           sumlik=(1.d0+dexp(-2*mue)+dexp(-4*mue*st)-4.*dexp(-2*mue*st)
     #          -2.*dexp(2*mue*(lab-1))+dexp(2*mue*(2*lab-1))
     #          +2.*dexp(-2*mue*lab)+4.*mue*(lab-st)+Hs(st-lab)*(-1.d0
     #          +4.*dexp(-2*mue*st)-2.*dexp(2*mue*(-2*st+lab))
     #          +dexp(4*mue*(lab-st))-2.*dexp(-2*mue*lab)
     #          +4.*mue*(st-lab)))/(2.*mue**2)
           
           sumlik0   = (2.*dexp(-2*mue*lab)+4.*mue*lab-2.d0+dexp(-2*mue)
     #         -2*dexp(2*mue*(lab-1))+dexp(2*mue*(2*lab-1)))/(2.*mue**2)
           
          sumlik = sumlik/sumlik0
        endif
c     mz        

      if(t0.gt.50.d0) then
      x = t0 
      tx = 1/sqrt(0.3141592653589793D1)*sqrt(1/x)-1/sqrt(0.3141592653589
     #793D1)*sqrt(1/x)**3/2+3.D0/4.D0/sqrt(0.3141592653589793D1)*sqrt(1/
     #x)**5-15.D0/8.D0/sqrt(0.3141592653589793D1)*sqrt(1/x)**7
      else
         tx = exp(t0)*derfc(sqrt(t0))
      endif

      
c mz
      if(td2.eq.0) then
         likhtman = (1.d0-eqd)*tx + eqd*sumlik
!test debye-waller
!         likhtman = dw*(1.d0-eqd)*tx + eqd*sumlik
      elseif (td2.eq.1) then
         likhtman = eqd*sumlik
      else
         likhtman = (1.d0-eqd)*tx
c         likhtman = (1.d0-eqd)*tx+eqd
!test
!         dw=dexp(-adw*q**2*t**0.25)
!         likhtman = dw*((1.d0-eqd)*tx+eqd)
      endif
c      write(*,*) 's(t) ',st
c mz

       return
       end


c ============================================================================== 
      double precision function Hs(x)
 
      implicit none
 
      double precision x
 
      if (x.le.0.) then
         Hs = 0.
      else
         Hs = 1.
      endif
 
      return
      end
 
c ==============================================================================
