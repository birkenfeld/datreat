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
         REAL            th17, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision w,d,q,n,ne,likhtman,t,a,l,taud,tau0,td2


 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'likhtman '
         nparx =  6
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
        td2     = pa(6)

        call parget('q       ',xh ,iadda,ier)
        if(ier.eq.0) then
          q       = xh
        else
          q       = 1.0  
        endif
 
        t = x 
 
        th17 = likhtman(t,q,d,W,n,l,ne,taud,tau0,td2)

        th17 = th17 * a

        call parset('ne      ',sngl(ne),iadda)
        call parset('taud    ',sngl(taud),iadda)
        call parset('tau0    ',sngl(tau0),iadda)

        return
        end 


       double precision function likhtman(t,q,d,W,n,l,ne,taud,tau0,td2)
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
c mz

       integer m,i,j

       parameter (pi=3.141592654d0)

       ne   =   (d/l)**2
       tau0 =  36.d0/(W*((q*l)**4))
       taud = 3.d0*(n**3)*(l**2)/( (pi**2)*W*(d**2)) 
c mz
       taue = ne**2/((pi**2)*W)

C setze taud als Fitparamter oder nicht (td2=0)


       if(td2.eq.0.d0) then
        td = t/taud
       else
        td = t/abs(td2)
       endif

       t0 = t/tau0
       eqd = exp(-((q*d/6.d0)**2))
  

C       m   = 10/sqrt(t/(abs(taud)+0.001D0)+0.001d0)+2
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

        sumlik = (n/(2*mue**2))*(2*mue+dexp(-2*mue)+2.d0-4*mue*st-4*
     +       dexp(-2*mue*st)+dexp(-4*mue*st))

        sumlik0 = (n/(2*mue**2))*(2*mue+dexp(-2*mue)-1.d0)

        sumlik = sumlik/sumlik0

c mz
        
c        do i=1,m

cC Loese Gleichung alphap*tan(alphap)-mue=0 mit Newton'schen Naeherungsverfahren
cC setze Startwert fuer alphap nahe (i-0.5)*pi
c          alphap=(i-0.5d0)*pi-0.1d-3

c          do j=1,150
          
c          diffe=alphap*dtan(alphap)-mue


c          if(abs(diffe).lt.0.1d-1) goto 122

c          alphap=alphap-(alphap*dtan(alphap)-mue)/
c     _       (dtan(alphap)+alphap*(1.d0+(dtan(alphap))**2))
         
c          enddo

c 122    CONTINUE
c        alp=alphap

c         arg1=-td*4.d0*alp**2/pi**2
c         if(arg1.lt.-300.d0) arg1=-300.d0
c         if(arg1.gt.300.d0) arg1=300.d0

c         sum = sum + 1.d0/(alp**2*(mue**2+alp**2+mue))*
c     _               (dsin(alp))**2*
c     _               dexp(arg1)


c         sumt0 = sumt0 + 1.d0/(alp**2*(mue**2+alp**2+mue))*
c     _                   (dsin(alp))**2


c        enddo


cC       sum = sum/sumt0 *eqd

c       sum = sum/sumt0

 
C alternative Berechnung des Roehrenformfaktors   
C       argd=q**2*n*13.05/6.d0
C       deb=n*2.d0/argd**2*(dexp(-argd)-1.d0+argd)
C       eqd=dexp(-(n*d**2/3.d0/13.05/deb))


      if(t0.gt.50.d0) then
      x = t0 
      tx = 1/sqrt(0.3141592653589793D1)*sqrt(1/x)-1/sqrt(0.3141592653589
     #793D1)*sqrt(1/x)**3/2+3.D0/4.D0/sqrt(0.3141592653589793D1)*sqrt(1/
     #x)**5-15.D0/8.D0/sqrt(0.3141592653589793D1)*sqrt(1/x)**7
      else
      tx = exp(t0)*derfc(sqrt(t0))
      endif


c mz
       likhtman = (1.d0-eqd)*tx + eqd*sumlik
c mz

       return
       end












