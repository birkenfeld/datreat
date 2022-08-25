c*ds
c*ed
c *********** thdatr1 fortran ******************************************
c                                                                      *
c ---> version for general data fitting !!!!
c                                                                      *
c***********************************************************************
c*                                                                     *
c*    provides theories for triax program system                      *
c*    please use the convention explained by the comments in the       *
c*    sample theories                                                  *
c*                                                                     *
c*    s(q,ome) is expected to be in units of sec                       *
c*             the scaling factor may be interpreted as                *
c*             n*sigma-inc in cm**2                                    *
c*             --> omit the 1/hquer factor of the usual def. of s(q,o) *
c*                                                                     *
c*    on input q (via qq(3) ) is given in miller-indizes               *
c*    om       is given in thz      (neg. = energy gain of neutron)    *
c*                                                                     *
c*                                                                     *
c***********************************************************************
c
c*ds
c*ds
      function thdatr1(x,pa,thnam,parnam,npar,ini)
c     =================================================
c
c -------> simple gauss :  g <-------
c          only energy dependence
c
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c -----> general construction scheme for theory functions <----
c        ------------------------------------------------
c
c ---- names must be thdatr1 , th2 .... thn
c      max. no is given by parameter mth in monito
c      to increase max no. change all mth='s & extend all lists
c      beginning with thdatr1( ...
c
c ---- inputs:
c      x    = x-value
c      pa(*)= parameters
c      ini  = initialisation flag (0 init , 1 calculate)
c ---- ouputs: (only for init=0)
c      thnam= name of theory
c      parnam=names of parameters
c      npar = no. of req. parameters (must be max. no. on input)
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
       parameter(pi=3.141592654)
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8 q,rave,drho,z,phi,den,rp,bq,pq,sq,r3,alpha,sqs,scale
       real*8 vol, cma, gu, go, epsilon, erroraccup, erroraccub
       real*8 pschulz, pschj1, betaj1, adapint, peryev  !! aix.sp extchk !!
       real*8 fiq1sch, fiq2sch                          !! aix
       common/cschulz/z  ,rave ,drho,q
       external fiq1sch
       external fiq2sch
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'j1ave'
c        ----------------> give here 8 character name of theory
c                          sherical formfactor with averaging
c                          by schulz-distribution and
c                          s(q) according to percus & yevik
c see. m.kotlarchyk and h.s.chen, j.chem.phys 79 (1983) 2461
         nparx = 9
c        ----------------> number of required parameters
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           thdatr1  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'drho*a*a'
         parnam(2) = 'rave/a'
         parnam(3) = 'z'
         parnam(4) = 'den*a**3'
         parnam(5) = 'rp/r3'
         parnam(6) = 'scale'
         parnam(7) = 'v/cm**3'
         parnam(8) = 'select'
         parnam(9) = 'rmax/int'
c
         thdatr1 = 0
         return
       endif
c
c ---- calculate theory here -----
c
         isel = pa(8)+0.5
         rmax = pa(9)
c --- selektion der funktion: isel=0 : j1ave --> streuung
c                             isel=1 : s(q)  (percus yevik)
c                             isel=2 : s'(q) (percus yevik)
c                             isel=3 : den*v*p(q)
c                             isel=4 : beta(q)
c                             isel=5 : pschulz(r)
c                             isel=6 : log(j1ave)
 
         cma  = 1d8
c               ---> converison a into cm
         vol  = pa(7)
         scale= pa(6)
         q    = x*scale
         drho = pa(1)
         rave = pa(2)/scale
         z    = pa(3)
         den  = pa(4)
         alpha= pa(5)
         r3   = ( (z+2)*(z+3)*rave**3/(z+1)**2 )**(1.d0/3.d0)
         rp   = alpha*r3
 
         if(rmax.le.0.0) then
           pq = pschj1(q,drho,rave,z)
           bq = betaj1(q,rave,z)
         else
c        ---> numerical integration
           maxit   = 50
           epsilon = 1d-6
           gu      = 0.d0
           go      = rmax
           pq  = adapint(fiq2sch,gu,go,epsilon,maxit,erroraccup)
           bq  = adapint(fiq1sch,gu,go,epsilon,maxit,erroraccub)
           bq  = bq *bq /pq
         endif
 
         sq = peryev(q,rp,den,-1d-7)
         sqs= 1+bq*(sq-1)
 
         thdatr1 = vol * cma * den * pq * sqs *scale**3
c        -------------------------
         if(isel.eq.6) thdatr1 = alog(thdatr1)
         if(isel.eq.1) thdatr1 = sq
         if(isel.eq.2) thdatr1 = sqs
         if(isel.eq.3) thdatr1 = vol * cma * den * pq * scale**3
         if(isel.eq.4) thdatr1 = bq
         if(isel.eq.5) thdatr1 = pschulz(q, rave, z)
c
       return
       end
c
c
c
c*ds
c*ds
      function th2(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> bkgr  <--------
c
c ------ linear background ---------------------------------------------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'bkgr    '
         nparx = 2
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th2  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'level '
         parnam(2) = 'slope '
c
         th2 = 0
         return
       endif
c
c ---- calculate theory here -----
c
       th2 = pa(1) + pa(2) * x
c      ------------------------
c
       return
       end
c
c
c*ds
c*ds
      function th3(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> evaluate from function-buffer <-----
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
       character*132 xformel,yformel,yfitform
       common/formul/xformel,yformel,yfitform
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8 val8y
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'eval  '
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th3  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'p(1)  '
         parnam(2) = 'p(2)  '
         parnam(3) = 'p(3)  '
         parnam(4) = 'p(4)  '
         parnam(5) = 'p(5)  '
         parnam(6) = 'p(6)  '
         parnam(7) = 'p(7)  '
         parnam(8) = 'p(8)  '
         parnam(9) = 'p(9)  '
         parnam(10)= 'p(10) '
c
         th3 = 0
         return
       endif
c
c ---- calculate theory here -----
       do i=1,nparx
         ptxf(i) = pa(i)
       enddo
       xxxx = x
       iout = iout - 5
       call evaluate(yfitform,val8y,iery)
       iout = iout + 5
       th3 = val8y
c
       return
       end
c
c
c
c*ds
c*ds
!      function th4(x,pa,thnam,parnam,npar,ini)
!c     ===================================================
!c
!c -------> 1/(a0+a2*q^2+a4*q^4) <-------
!c
!c
!       character*8 thnam,parnam(20)
!       dimension pa(20),qq(3)
!       data zpi/6.283185/
!c
!c ----- initialisation -----
!       if(ini.eq.0) then
!         thnam = 'emul'
!         nparx = 3
!         if(npar.lt.nparx) then
!           write(6,1)thnam,nparx,npar
!1          format(' theory: ',a8,' no of parametrs=',i8,
!     *      ' exceeds current max. = ',i8)
!           th4  = 0
!           return
!         endif
!         npar = nparx
!c        --------------> set the number of parameters
!         parnam(1) = 'a0'
!         parnam(2) = 'a2'
!         parnam(3) = 'a4'
!c
!         th4 = 0
!         return
!       endif
!c
!c ---- calculate theory here -----
!       xqq = x*x
!       th4 = 1.0 / (pa(1)+pa(2)*xqq+pa(3)*xqq*xqq)
!c      -------------------------------------------
!c
!       return
!       end
!


      function th4(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8 m1,m2,m3,m4,m5,m6,m7,m8,m9
       real*8 c1,c2,c3,c4,c5,c6,c7,c8,c9
       real*8 ampl, r, y

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'approxi'
         nparx = 15
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th4  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'ampl'
         parnam(2) = 'm1'
         parnam(3) = 'm2'
         parnam(4) = 'm3'
         parnam(5) = 'm4'
         parnam(6) = 'm5'
         parnam(7) = 'c1'
         parnam(8) = 'c2'
         parnam(9) = 'c3'
         parnam(10)= 'c4'
         parnam(11)= 'c5'
         parnam(12)= 'm6'
         parnam(13)= 'c6'
         parnam(14)= 'm7'
         parnam(15)= 'c7'
c
         th4 = 0
         return
       endif
c
        r    = exp(x)
        ampl = pa(1 )
        m1   = pa(2 )
        m2   = pa(3 )
        m3   = pa(4 )
        m4   = pa(5 )
        m5   = pa(6 )
        c1   = pa(7 )
        c2   = pa(8 )
        c3   = pa(9 )
        c4   = pa(10)
        c5   = pa(11)
        m6   = pa(12)
        c6   = pa(13)
        m7   = pa(14)
        c7   = pa(15)
        
c ---- calculate theory here -----
       y = ampl*(r**m1+c5*r**m3+c4*r**m5)/
     *          (c1+c2*r**m2+c3*r**m4+c6*r**m6+c7*r**m7)
       th4 = dlog(y)
c      -----------------------------------------------
c
       return
       end
c
c
c*ds
c*ds
!      function th5(x,pa,thnam,parnam,npar,ini)
!c     ===================================================
!c
!c -------> saxabs <------
!c
!c
!       character*8 thnam,parnam(20)
!       dimension pa(20),qq(3)
!       data zpi/6.283185/
!c
!c ----- initialisation -----
!       if(ini.eq.0) then
!         thnam = 'saxabs'
!         nparx = 19
!         if(npar.lt.nparx) then
!           write(6,1)thnam,nparx,npar
!1          format(' theory: ',a8,' no of parametrs=',i8,
!     *      ' exceeds current max. = ',i8)
!           th5  = 0
!           return
!         endif
!         npar = nparx
!c        --------------> set the number of parameters
!         parnam(1) = 'i(q=0)'
!         parnam(2) = 'r0'
!         parnam(3) = 'n'
!         parnam(4) = 'd1'
!         parnam(5) = 'b1'
!         parnam(6) = 'd2'
!         parnam(7) = 'b2'
!         parnam(8) = 'd3'
!         parnam(9) = 'b3'
!         parnam(10)= 'd4'
!         parnam(11)= 'b4'
!         parnam(12)= 'd5'
!         parnam(13)= 'b5'
!         parnam(14)= 'd6'
!         parnam(15)= 'b6'
!         parnam(16)= 'd7'
!         parnam(17)= 'b7'
!         parnam(18)= 'd8'
!         parnam(19)= 'b8'
!c
!         th5 = 0
!         return
!       endif
!c
!c ---- calculate theory here -----
!       r0   = pa(2)
!       n    = pa(3) + 0.1
!       if(n.gt.8) then
!         n = 8
!         write(6,*)'th5(saxabs): n-restricted to 8 !'
!       endif
!       th5 = pa(1) * saxabs(x, r0, pa(4), n)
!c
!       return
!       end

      function th5(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> triblock1 approximation  <------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)

       data zpi/6.283185/

       double precision t, q, rate, dispe, uq, a1, a2
       real*4 qp
       common/thiadd/iadda

c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'tribapr'
         nparx = 3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th5  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'uquadrat'
         parnam(2) = 'rate'
         parnam(3) = 'dispers'
c
         th5 = 0
         return
       endif
c
c ---- calculate theory here -----
       t     =  x
       uq    = pa(1)
       rate  = pa(2)
       dispe = pa(3)

! extract parameter ! 
       call        parget('q       ',qp,iadda,ier)
       if(ier.ne.0) then
         write(6,*)'th5-brush: q-parameter not found! set q=0.1'
         qp = 0.1
       endif
       q = qp   

       if(dispe.gt.10.0d0) then
        write(6,*)'warning: dispe=',dispe,' readjusted'
        dispe = 10
       endif
       if(dispe.lt.-10.0d0) then
        write(6,*)'warning: dispe=',dispe,' readjusted'
         dispe = -10
       endif

       a1 = -(q*q*uq*uq)
       a2 = -abs(rate)*(q**dispe)*t
       if(a1.lt.-50d0) a1 = -50.0d0
       if(a2.lt.-50d0) a2 = -50.0d0

       th5 = exp(a1)+(1-exp(a1))*exp(a2)
c
       return
       end
c
c*ds
c*ds
!tillNOV00!      function th6(x,pa,thnam,parnam,npar,ini)
!tillNOV00!c     ===================================================
!tillNOV00!c
!tillNOV00!c -------> echo <--------
!tillNOV00!c
!tillNOV00!c
!tillNOV00!       character*8 thnam,parnam(20)
!tillNOV00!       dimension pa(20),qq(3)
!tillNOV00!       data zpi/6.283185/
!tillNOV00!       data gam_larmor/183.0333e6/
!tillNOV00!       data amass_n/0.16749286e-26/
!tillNOV00!       data hplanck/0.66260755e-33/
!tillNOV00!c
!tillNOV00!c ----- initialisation -----
!tillNOV00!       if(ini.eq.0) then
!tillNOV00!         thnam = 'echo '
!tillNOV00!         nparx = 11
!tillNOV00!         if(npar.lt.nparx) then
!tillNOV00!           write(6,1)thnam,nparx,npar
!tillNOV00!1          format(' theory: ',a8,' no of parametrs=',i8,
!tillNOV00!     *      ' exceeds current max. = ',i8)
!tillNOV00!           th6  = 0
!tillNOV00!           return
!tillNOV00!         endif
!tillNOV00!         npar = nparx
!tillNOV00!c        --------------> set the number of parameters
!tillNOV00!         parnam(1)  = 'amplitu'
!tillNOV00!         parnam(2)  = 'average'
!tillNOV00!         parnam(3)  = 'derotat'
!tillNOV00!         parnam(4)  = 'lambda0'
!tillNOV00!         parnam(5)  = 'lamfwhm'
!tillNOV00!         parnam(6)  = 'ishift '
!tillNOV00!         parnam(7)  = 'nturns '
!tillNOV00!         parnam(8)  = 'rcoil  '
!tillNOV00!         parnam(9)  = 'lcoil  '
!tillNOV00!         parnam(10) = 'acoil  '
!tillNOV00!         parnam(11) = 'bcoil  '
!tillNOV00!      
!tillNOV00!c
!tillNOV00!         th6 = 0
!tillNOV00!         return
!tillNOV00!       endif
!tillNOV00!c
!tillNOV00!c ---- calculate theory here -----
!tillNOV00!       xs  = x-pa(6)
!tillNOV00!       alam = pa(4)*1e-10
!tillNOV00!       dlam = pa(5)*1e-10*0.5/sqrt(alog(2.0))
!tillNOV00!       xtur = pa(7)
!tillNOV00!       rcoil= pa(8)
!tillNOV00!       coill= pa(9)
!tillNOV00!       acoil= pa(10)
!tillNOV00!       bcoil= pa(11) 
!tillNOV00!       omeg = gam_larmor*amass_n/hplanck* (2*zpi*1e-7)
!tillNOV00!       omeg = omeg * xtur * (blinteg(acoil,bcoil,coill,rcoil)-
!tillNOV00!!     *                      blinteg(-100.0,acoil,coill,rcoil))/2        
!tillNOV00!     *                      blinteg(-0.4535,acoil,coill,rcoil))/2        
!tillNOV00!
!tillNOV00!       write(6,*)'omeg: ',omeg,omeg*alam,omeg*alam*57.29578 
!tillNOV00!       th6 = pa(2) + pa(1) * 
!tillNOV00!     *       exp(-(xs*omeg*dlam/2)**2) * cos(omeg*alam*(xs+pa(3)))
!tillNOV00!c
!tillNOV00!       return
!tillNOV00!       end
!tillNOV00!c*ds
c*ds
      function th6(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> echo <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       common/thiadd/iadda
       data zpi/6.283185/
       data gam_larmor/183.0333e6/
       data amass_n/0.16749286e-26/
       data hplanck/0.66260755e-33/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'echo '
         nparx = 12
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th6  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1)  = 'amplitu'
         parnam(2)  = 'average'
         parnam(3)  = 'derotat'
         parnam(4)  = 'lambda0'
         parnam(5)  = 'lamfwhm'
         parnam(6)  = 'ishift '
         parnam(7)  = 'nturns '
         parnam(8)  = 'rcoil  '
         parnam(9)  = 'lcoil  '
         parnam(10) = 'acoil  '
         parnam(11) = 'bcoil  '
         parnam(12) = 'a2coil  '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Zugrunde liegende Geometrie                                                      !!
!!                                                                                  !!
!!                                                                                  !!
!!  a2coil                 acoil                   bcoil                            !!
!!                                 lcoil                                            !!
!!  |----------------------|-------====------------|                                !!
!!                                 ^                                                !!
!!  pi/2                   pi      0-pkt           pi/2                             !!
!!                                                                                  !!
!!                          <acoil><====bcoil=====>                                 !!
!!  <=====acoil2==================>                                                 !!
!!                                 <lcoil>                                          !!
!!                                                                                  !!
!!                                                                                  !!
!!  rcoil: effektiver Spulenradius                                                  !!
!!                                                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
c
         th6 = 0
         return
       endif
c
c ---- calculate theory here -----
       xs  = x-pa(6)
       alam = pa(4)*1e-10
       dlam = pa(5)*1e-10*0.5/sqrt(alog(2.0))
       xtur = pa(7)
       rcoil= pa(8)
       coill= pa(9)
       acoil= pa(10)
       bcoil= pa(11)
       acoil2=pa(12) 
       omeg = gam_larmor*amass_n/hplanck* (2*zpi*1e-7)
       omeg = omeg * xtur * (blinteg(acoil,bcoil,coill,rcoil)-
     *                       blinteg(acoil2,acoil,coill,rcoil))/2        


!      write(6,*)'omeg: ',omeg,omeg*alam,omeg*alam*57.29578 
       call parset('pha_sens ',sngl(omeg*alam*57.29578d0),iadda)

       th6 = pa(2) + pa(1) * 
     *       exp(-(xs*omeg*dlam/2)**2) * cos(omeg*alam*(xs+pa(3)))
c
       return
       end


       function blinteg(a,b,xl,r)

         blinteg = 
     *           sqrt( (r/xl)**2 + ((xl-a)/xl)**2 ) 
     *         - sqrt( (r/xl)**2 + ((xl-b)/xl)**2 )   
     *         + sqrt( (r/xl)**2 + (b/xl)**2 )   
     *         - sqrt( (r/xl)**2 + (a/xl)**2 )   
       return
       end 
c
c
c*ds
c*ds
      function th7(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> test fuer percus-yevick s(q) --------------------------------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8 pschulz, pschj1, betaj1, adapint, peryev  !! aix.sp extchk !!
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'peryev'
         nparx = 4
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th7  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intens'
         parnam(2) = 'radius'
         parnam(3) = 'density'
         parnam(4) = 'epsilon'
c
         th7 = 0
         return
       endif
c
c ---- calculate theory here -----
c
       th7 = pa(1)*peryev(dble(x),dble(pa(2)),dble(pa(3)),dble(pa(4)))
 
       return
       end
c
c
c*ds
c*ds
c ---> th8 extern dazulinken !
c ---> th8 extern dazulinken !
 
      function th10(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> q**n <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real n                        !!<-------!!
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'q**n'

c mz: erweitert um bk-add

         nparx = 4

         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th10  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'an'
         parnam(2) = 'n'
         parnam(3) = 'offset'
         parnam(4) = 'bk    '
c
         th10 = 0
         return
       endif
c
c ---- calculate theory here -----
       n = pa(2) 
       offset = pa(3)
       th10 = pa(1) * (abs(x-offset))**n + pa(4)
c
       return
       end
c*ds
c*ds
!!
!!!!!! <- th9 dazulinken <-
!!
c*ds
c*ds
      function th12(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> scattering of spheres (j1) -------------------------
c
c
       parameter(pi=3.141592654)
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'j1qr'
         nparx = 3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th12  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'density'
         parnam(2) = 'r'
         parnam(3) = 'drho'
c
         th12 = 0
         return
       endif
c
c ---- calculate theory here -----
       r   = pa(2)
       drho= pa(3)
       z   = abs(x*r)
       aj1 = (sin(z)-z*cos(z))/(z**3)
       fj1qr= 4*pi/3*r*r*r*drho*3*aj1
       th12 = fj1qr**2 * pa(1)
c
       return
       end
c*ds
c*ds
      function th13(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> gaussian <----
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'gauss'
         nparx = 3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th13  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'width'
         parnam(3) = 'center'
c
         th13 = 0
         return
       endif
c
c ---- calculate theory here -----
       arg  = ((x-pa(3))/pa(2))**2
       if(arg .gt. 50.0 ) arg = 50.0
       th13 = pa(1) * exp( - arg )
c
       return
       end
!! th16 (ashelli)
!!
!! th18
!! th19
!! th20
c*ds
c*ds
      function th20(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> sofq <--------
c
c
       parameter(pi=3.141592654)
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8 qqr(1),sq,eta,scl,gamma,r
       real*8 d, den
       real*8 pschulz, pschj1, betaj1, adapint, peryev  !! aix.sp extchk !!
       real*8 q                                         !! aix
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'sofq'
         nparx = 5
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th20  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'volfrac'
         parnam(2) = 'scalelen'
         parnam(3) = 'gamma'
         parnam(4) = 'r'
         parnam(5) = 'mode'
c
         th20 = 0
         return
       endif
c
c ---- calculate theory here -----
         eta       = pa(1)
         scl       = pa(2)
         gamma     = pa(3)
         r         = pa(4)
         mode      = pa(5) + 0.1
 
         if(mode.ge.1 .and. mode. le. 5) then
           qqr(1)    = x * r * 2.0
           call sofq(qqr,sq,1,eta,scl,gamma,r,mode,ierr)
           if(ierr.lt.0) write(6,*)'sofq: ierr=',ierr
         else
           if(mode.eq.6) then
             q  = x
             d  = 2*r
             den= 6.0d0 * eta / (pi*d*d*d)
             sq = peryev(q,r ,den,-1d-7)
           else
             write(6,*)'sofq: mode=',mode,' is out of range'
             sq = 1.0
           endif
         endif
       th20 = sq
c
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     th23 aus externem Prog. per link                            !!                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function th25(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> formfactor squared of a multiple square well
c
!!     |-dn-|....|-d2-|--d1--|---d0---|--d1--|-d2-|....|-dn-|   
!!       bn        b2    b1      b0      b1    b2        bn    
!!
!!
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       double precision d(10),b(10),c(10),sum,q
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'fqsheet'
         nparx = 20
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th25  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'd0'
         parnam(2) = 'd1'
         parnam(3) = 'd2'
         parnam(4) = 'd3'
         parnam(5) = 'd4'
         parnam(6) = 'd5'
         parnam(7) = 'd6'
         parnam(8) = 'd7'
         parnam(9) = 'd8'
         parnam(10)= 'd9'
         parnam(11)= 'b0'
         parnam(12)= 'b1'
         parnam(13)= 'b2'
         parnam(14)= 'b3'
         parnam(15)= 'b4'
         parnam(16)= 'b5'
         parnam(17)= 'b6'
         parnam(18)= 'b7'
         parnam(19)= 'b8'
         parnam(20)= 'b9'
c
         th25 = 0
         return
       endif
c
c ---- calculate theory here -----
       q = x
       if(abs(q).lt.1d-8) q = 1d-8
     
       do i=1,10
        d(i) = pa(i)
        b(i) = pa(i+10)
       enddo

       b(1) = b(1)/2
       c(1) = 0.0
       do i=2,10
         c(i) = c(i-1)+(d(i-1)+d(i))/2
       enddo

       sum = 0.0
       do i=1,10
        sum = sum + b(i)*2*cos(c(i)*q)*sin(q*d(i)/2)
       enddo

       th25 = (2/q)*sum

       th25 = th25*th25

       return
       end
c*ds
c*ds
      function th24(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> formfactor squared of a double square-well
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'doubsheq'
         nparx = 4
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th24  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'd_outer'
         parnam(2) = 'd_inner'
         parnam(3) = 'b_outer'
         parnam(4) = 'b_inner'
c
         th24 = 0
         return
       endif
c
c ---- calculate theory here -----
       dd   = 0.5 * ( pa(1) + pa(2) )
       th24 = (2/x) * (2*pa(3)*cos(dd*x)*sin(pa(1)*x*0.5) +
     *                   pa(4)*sin(pa(2)*x*0.5) )
       th24 = th24*th24

       return
       end
c*ds
c*ds
      function th26(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> diffqav <--------
c
c average over exp(-Ð*q**2 t) assuming Intensity goes like 1/q**2
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       common/thiadd/iadda
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'diffqav'
         nparx = 3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th26  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplit'
         parnam(2) = 'd'
         parnam(3) = 'qwid/2'


         th26 = 0.0

         return
       endif
c
c ---- calculate theory here -----
        q = 0.0
        call parget('q       ',q,iadda,ier)
        if(q.eq.0) write(6,*)'ERROR: q not found'
        if(q.le.abs(pa(3))) q = abs(pa(3))+0.001

        t = x
        D = abs(pa(2))
        q1= q - abs(pa(3))
        q2= q + abs(pa(3))

   
        if(q2.le.q1) q2 = q1+0.0001
        if(t.le.1e-6)t  = 1e-6
        
        if(d*q*q*t.gt.50.0) d=50.0/(q*q*t)

        t1 = q2**2
        t7 = sqrt(D*t)
        t11 = sqrt(0.3141593E1)
        t17 = q1**2
        t35 = (q1*exp(-D*t1*t)*t7+q1*D*t*t11*erf(t7*q2)*q2
     *       -q2*exp(-D*t17*t)*t7-q2*D*t*t11*erf(t7*q1)*q1)
     *       /t7/(-q2+q1)

         th26 = t35*pa(1)
c
       return
       end
c*ds
c*ds
      function th27(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> cumulant <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real      mu, ampli, sum, ff
       dimension mu(20)
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'cumulant'
         nparx = 7
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th27  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1)  = 'amplitu'
         parnam(2)  = 'mu1'
         parnam(3)  = 'mu2'
         parnam(4)  = 'mu3'
         parnam(5)  = 'mu4'
         parnam(6)  = 'mu5'
         parnam(7)  = 'mu6'

c
         th27 = 0
         return
       endif
c
c ---- calculate theory here -----
       ampli = pa(1)
       do i=2,nparx
         mu(i+1-2) = pa(i)
       enddo

       ff   = 1.0
       sum  = 0.0
       do i=1,6
        ff = ff*i
        sum = sum + mu(i)*(x**i)/ff
       enddo
       if(sum.gt. 20.0) sum =  20.0
       if(sum.lt.-20.0) sum = -20.0

       th27 = ampli * exp(-sum)
c
       return
       end
c*ds
c*ds
      function th28(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> strexpo  <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)

       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'strexpo'
         nparx = 4
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th28  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitu'
         parnam(2) = 'tau'
         parnam(3) = 'beta'
         parnam(4) = 'qexp'
c
         th28 = 0
         return
       endif
c
c ---- calculate theory here -----
       amp  = pa(1)
       tau0 = abs(pa(2))
       bet  = pa(3)
       qexp = pa(4)

       q = 1.0
       call  parget('q       ',q,iadda,ier)
       q = abs(q)
 
       tau = tau0 * q**qexp
        
       if(bet.lt.  0.0) bet =  0.0
       if(bet.gt.  1.0) bet =  1.0
       arg = abs(x/tau)**bet
       if(arg.lt.-30.0) arg = -30.0
       if(arg.gt. 30.0) arg =  30.0  
      
       th28 = amp * exp(-arg)
c
       return
       end
c*ds
c*ds
      function th29(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> exp_q <--------
!!----------------------------------------------------------
!! 2-Exponentials with Q-dependend Amplitude-Ratio and Decay
!! Amplitude = p1
!! Aratio    = p2 + p3*Q**p4+p5*Q**p6+p7*Q**p8
!! Tau1      = p9 + p10*Q**p11 + p12*Q**p13
!! Beta1     = p14
!! Tau2      = p15+ p16*Q**p17 + p18*Q**p19
!! Beta2     = p20
!!----------------------------------------------------------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
       common/thiadd/iadda

c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'exp_q'
         nparx = 20
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th29  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitud'
         parnam(2) = 'aratioc0'
         parnam(3) = 'aratioc1'
         parnam(4) = 'aratioe1'
         parnam(5) = 'aratioc2'
         parnam(6) = 'aratioe2'
         parnam(7) = 'aratioc3'
         parnam(8) = 'aratioe3'
         parnam(9) = 'tau1c0  '
         parnam(10)= 'tau1c1  '
         parnam(11)= 'tau1e1  '
         parnam(12)= 'tau1c2  '
         parnam(13)= 'tau1e2  '
         parnam(14)= 'beta1   '
         parnam(15) ='tau1c0  '
         parnam(16)= 'tau2c1  '
         parnam(17)= 'tau2e1  '
         parnam(18)= 'tau2c2  '
         parnam(19)= 'tau2e2  '
         parnam(20)= 'beta2   '
c
         th29 = 0
         return
       endif
c
c ---- calculate theory here -----

        amplitud  = pa(1) 
        aratioc0  = pa(2) 
        aratioc1  = pa(3) 
        aratioe1  = pa(4) 
        aratioc2  = pa(5) 
        aratioe2  = pa(6) 
        aratioc3  = pa(7) 
        aratioe3  = pa(8) 
        tau1c0    = pa(9) 
        tau1c1    = pa(10)
        tau1e1    = pa(11)
        tau1c2    = pa(12)
        tau1e2    = pa(13)
        beta1     = pa(14)
        tau2c0    = pa(15)
        tau2c1    = pa(16)
        tau2e1    = pa(17)
        tau2c2    = pa(18)
        tau2e2    = pa(19)
        beta2     = pa(20)

        q = 1.0
        call        parget('q       ',q,iadda,ier)
        q = abs(q)
        t = abs(x)


        aratio = aratioc0 + aratioc1*q**aratioe1 + aratioc2*q**aratioe2 
     *                    + aratioc3*q**aratioe3

        a1     = 1.0/(1.0+aratio)
        a2     = 1.0 - a1

        tau1   = tau1c0   + tau1c1*q**tau1e1     + tau1c2*q**tau1e2
        tau2   = tau2c0   + tau2c1*q**tau2e1     + tau2c2*q**tau2e2

        if (tau1.gt.0.0d0) then 
          arg1   = -(t/tau1)**beta1
        else
          arg1   = -30.0
        endif

        if (tau2.gt.0.0d0) then 
          arg2   = -(t/tau2)**beta2
        else
          arg2   = -30.0
        endif

        if(arg1.lt.-30.0) arg1 = -30.0
        if(arg2.lt.-30.0) arg2 = -30.0

        
        th29 = amplitud*(a1*exp(arg1)+a2*exp(arg2))
c
       return
       end
c*ds
c*ds
      function th30(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> gaupol <--------
c
c
        character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'gaupol'
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th30  = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'center'
         parnam(3) = 'a1'
         parnam(4) = 'a2'
         parnam(5) = 'a3'
         parnam(6) = 'a4'
         parnam(7) = 'a5'
         parnam(8) = 'a6'
         parnam(9) = 'a7'
         parnam(10)= 'a8'
c
         th30 = 0
         return
       endif
c
c ---- calculate theory here -----
       y = x-pa(2)
       arg  =      pa(10)*y
       arg =  (arg+pa(9))*y
       arg =  (arg+pa(8))*y
       arg =  (arg+pa(7))*y
       arg =  (arg+pa(6))*y
       arg =  (arg+pa(5))*y
       arg =  (arg+pa(4))*y
       arg =  (arg+pa(3))*y
       if(arg .gt. 50.0 ) arg = 50.0
       th30 = pa(1) * exp( - arg )
c
       return
       end
c
c*ds
c*ds
      function th31(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> zimm <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_zimm, a
       real*8                            SQ_zimmT
       real*8    epsilon, diff
       real*4    qget, tget

 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'zimm'
         nparx = 5
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th31 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'eta_solv'
         parnam(3) = 'epsilon '
         parnam(4) = 'temp    '
         parnam(5) = 'com_diff'

c
         th31 = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a        = pa(1)
       eta      = abs(pa(2))
       epsilon  = abs(pa(3))
       temp     = pa(4)
       diff     = abs(pa(5))  ! in cm**2/sec

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) write(6,*)'Warning q not found' 
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif

c --- include center of mass diffusion ---
       a = a * dexp(-qz*qz*diff*tau)
       
       if(pa(3).lt.0) then
         write(6,*)'Full zimm computation from scratch!'
         th31 = a * SQ_zimm(tau,qz,temp,eta,epsilon)
       else
         th31 = a * SQ_zimmT(tau,qz,temp,eta,epsilon)
       endif

c
       return
       end
c*ds
c*ds
      function th32(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> rouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8                            SQ_rouseT
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff
       real*4    qget, tget

 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'rouse'
         nparx = 7
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th32 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'xi_frict'
         parnam(3) = 'b_segmnt'
         parnam(4) = 'epsilon '
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width'

c
         th32 = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       xi       = abs(pa(2))
       b        = abs(pa(3))
       epsilon  = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) write(6,*)'Warning q not found' 
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif


      th32    = 0
      sum     = 0
      sumnorm = 0
      nqw     = 15
      dqw     = 4*q_width/nqw
      if(q_width.eq.0) then
        nqw = 0
        q_width = 1.0d0
      endif 

      do i=-nqw,nqw
       qzz = qz + i*dqw
       if(qzz.gt.0) then 
         fn   = dexp( -(i*dqw/q_width)**2)
         sumnorm = sumnorm + fn

c --- include center of mass diffusion ---
         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
       
         if(pa(4).lt.0) then
           write(6,*)'Full rouse computation from scratch!'
           sum = sum + a * SQ_rouse(tau,qzz,temp,xi,b,epsilon)
         else
           sum = sum + a * SQ_rouseT(tau,qzz,temp,xi,b,epsilon)
         endif
        endif
       enddo

       th32 = sum/sumnorm
c
       return
       end


      function th33(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> rouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8                            SQ_rouseT
       real*8    epsilon, diff
       real*4    qget, tget

 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'rouse61'
         nparx = 6
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th33 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'xi_frict'
         parnam(3) = 'b_segmnt'
         parnam(4) = 'epsilon '
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'

c
         th33 = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a        = pa(1)
       xi       = abs(pa(2))
       b        = abs(pa(3))
       epsilon  = abs(pa(4))
       temp     = pa(5)
       diff     = pa(6)  ! in cm**2/sec

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) write(6,*)'Warning q not found' 
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif

c --- include center of mass diffusion ---
       a = a * dexp(-qz*qz*diff*tau)
       
       th33 = a * rouse61(tau, qz, xi, b, temp, epsilon )
       write(6,*)'rouse61(',tau,qz,')=',th33
c
       return
       end
c*ds
c*ds
c*ds
c*ds
      function th34(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c     ===================================================
c
c -------> caille <------
c
c
       character*8 thnam,parnam(20)
       dimension pa(40),qq(3)
       data zpi/6.283185/
       real*8 qortho, qz, a, alpha, d, q0, dwf, peakf, scale
       real*8 eps, eps2, xmuswit, strechf, rellim, splitf4,rres
       real*8 caill4a,caille3,caille4,caille2,caillr3,caills3
       real*8 caille
 
       common/thiadd/iadda
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'caille'
         nparx = 20
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th34 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'scale'
         parnam(2) = 'qz:ortho'
         parnam(3) = 'iz:ortho'
         parnam(4) = 'q0'
         parnam(5) = 'a'
         parnam(6) = 'ln_alpha'
         parnam(7) = 'd'
         parnam(8) = '4a:3r:3s'
         parnam(9) = 'nmax3'
         parnam(10)= 'dwf'
         parnam(11)= 'peakf'
         parnam(12)= 'eps'
         parnam(13)= 'eps2'
         parnam(14)= 'xmuswit'
         parnam(15)= 'strechf'
         parnam(16)= 'rellim'
         parnam(17)= 'splitf4'
         parnam(18)= 'e1switch'
         parnam(19)= 'polpow'
         parnam(20)= 'rres'
c
         th34 = 0
         return
       endif

       iot = iout()
c
c ---- calculate theory here -----
 
        scale      = pa( 1)
        qz_ortho   = pa( 2)
        iz_ortho   = pa( 3)+0.1
        q0         = pa( 4)
        a          = pa( 5)
        alpha      = pa( 6)
        d          = pa( 7)
        i4a3r3s    = pa( 8)+0.1
        nmax3      = pa( 9)+0.1
        dwf        = pa(10)
        peakf      = pa(11)
        eps        = pa(12)
        eps2       = pa(13)
        xmuswit    = pa(14)
        strechf    = pa(15)
        rellim     = pa(16)
        splitf4    = pa(17)
       ie1switch   = pa(18)+0.1
       ipolpow     = pa(19)+0.1
        rres       = pa(20)
 
        if (alpha.gt.69.0) alpha = 69
        if (alpha.lt.-69.0) alpha = -69

        alpha = exp(alpha)

        call setcai(eps,eps2,xmuswit,strechf,rellim,splitf4,
     *             rres,itch,ipolpow)
 
c --- look for qz or qortho parameters in the files ---
        ier = -1
        call parget('qz      ',qz_ortho,iadda,ier)
        if(ier.eq.0) iz_ortho = 0
c     --> qortho overrides qz
        ier = -1
        call parget('qortho  ',qz_ortho,iadda,ier)
        if(ier.eq.0) iz_ortho = 1
c --- and look for a, alpha and q0
        ier = -1
        aa = a
        call parget('a       ',aa      ,iadda,ier)
        a  = aa
        ier = -1
        aa  = alpha
        call parget('alpha   ',aa      ,iadda,ier)
        alpha = aa
        ier = -1
        aa  = q0
        call parget('q0      ',aa      ,iadda,ier)
        q0  = aa
        ier = -1
        call parget('q_rock  ',q_rock  ,iadda,ier)
        if(ier.eq.0) then 
          i_rock = 1 
        else
          i_rock = 0
        endif

 
        if(i_rock.eq.0) then
         if(iz_ortho.eq.0) then
           qz      = qz_ortho
           qortho  = x
           if(iot.gt.0)write(6,*)'qortho =',x,' at qz     =',qz_ortho
         else
           qz      = x
           qortho  = qz_ortho
           if(iot.gt.0)write(6,*)'qz     =',x,' at q_ortho=',qz_ortho
         endif

        goto (10,20,30),i4a3r3s
10      continue
         caille = caill4a(qortho,qz,q0,d,a,alpha)
         goto 100
20      continue
         caille = caillr3(qortho,qz,q0,d,a,alpha,nmax3)
         goto 100
30      continue
         caille = caills3(qortho,qz,q0,d,a,alpha,peakf,dwf)
100     continue
 
         if(iot.gt.0)write(6,*)'caille = ',caille,' scale= ', scale
         th34 = caille * scale

        else

c  this is a rocking curve, x is psi in degrees
           sr = sin(x*3.14159/180.0)
           cr = cos(x*3.14159/180.0)
           qz = q_rock*cr
           qortho = q_rock*sr
           if(iot.gt.0)write(6,*)'Rocking curve x=',x,qz,qortho      
 
           ier = -1
           xnsum = 1.0
           call parget('nsum    ',xnsum  ,iadda,ier)
           nsum = xnsum+0.1

           ier = -1
           sig_qr = 1.0
           call parget('sig_qr  ',sig_qr  ,iadda,ier)
         

!! here we introduce some special resolution convolution for the
!! rocking curve
           caille = 0.0        
           fsum   = 0.0
        
           if(iot.gt.0)write(6,*)'nsum=',nsum,' sig_qr=',sig_qr               

           do isum = -nsum/2,nsum/2
            
            if(nsum.gt.1) then
              dqr = sig_qr*2/(nsum/2)
            else
              dqr = 1.0
            endif
            q_r = q_rock+isum*dqr
            
            if(q_r.le.0.0) goto 1111 
               
            qz     = q_r*cr
            qortho = q_r*sr

            qz     = abs(qz)
            qortho = abs(qortho)
           

            f = exp(-(isum*dqr/sig_qr)**2)        
 
           goto (110,120,130),i4a3r3s
110        continue
           caille = caille+f*caill4a(qortho,qz,q0,d,a,alpha)
           goto 1100
120        continue
           caille = caille+f*caillr3(qortho,qz,q0,d,a,alpha,nmax3)
           goto 1100
130        continue
           caille = caille+f*caills3(qortho,qz,q0,d,a,alpha,peakf,dwf)
1100       continue
 
           fsum = fsum + f
           
           if(iot.gt.1) then
             write(6,*)'isum=',isum,' qr=',q_r
             write(6,*)'f=',f,' caille=',caille 
           endif

 1111      continue

        enddo ! isum

        if(iot.gt.0) write(6,*)'fsum = ',fsum
        if(fsum.gt.0.0) caille = caille/fsum
        if(iot.gt.0) write(6,*)'caille = ',caille,' scale= ', scale
 
        th34 = caille * scale

       endif  !< rocking curve !
c
       return
       end
c*ds
c*ds
      function th35(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> formfactor squared of a double square-well
c
c
       character*8 thnam,parnam(20)
       dimension pa(40),qq(3)
       common/thiadd/iadda
       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'doubshee'
         nparx = 6
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th35 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'd_outer'
         parnam(2) = 'd_inner'
         parnam(3) = 'b_outer'
         parnam(4) = 'b_inner'
         parnam(5) = 'qz:qx'
         parnam(6) = 'chi'
c

         th35 = 0
         return
       endif
c
        qz_qx = pa(5)
        chi   = pa(6)
        ier = -1
        call parget('qz:qx   ',qz_qx   ,iadda,ier)
        call qcalc(x,qz_qx,chi,qz,qortho)
          

c ---- calculate theory here -----
       dd   = 0.5 * ( pa(1) + pa(2) )
c       th35 = (2/qz) * (2*pa(3)*cos(dd*qz)*sin(pa(1)*qz*0.5) +
c     *                   pa(4)*sin(pa(2)*qz*0.5) )
     
        hv = (2/qz)*pa(4)*sin(0.5*pa(2)*qz) +
     _  2*pa(3)*((sin(0.5*(2*pa(1)+pa(2))*qz))/qz -
     _  sin(0.5*pa(2)*qz)/qz)
    
        th35=hv*hv

c
       return
       end
c
!!       th36 (alt rectab durch externes th36 ersetzt z.B pepep_nr2_th36.f )

!!       th37 (alt caille, ersetzt 21 plabu3, neu hier pepep3 )
!!
!!       th38 (alt dreiexp, ersetzt durck pepep_4, extern)
c*ds
c*ds

c*ds
c*ds
      function th40(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> expq <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       common/thiadd/iadda


       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'expq'
         nparx = 4
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th40 = 0.0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitu'
         parnam(2) = 'tau0'
         parnam(3) = 'q0'
         parnam(4) = 'q_exp'
c
         th40 = 0
         return
       endif
c
c ---- calculate theory here -----
       call        parget('q       ',qp,iadda,ier)
       t      = x
       amp    = pa(1)
       tau0   = pa(2)
       q0     = pa(3)
       q_exp  = pa(4)


       if(q_exp .lt. -30.0) q_exp = -30.0
       if(q_exp .gt.  30.0) q_exp =  30.0
       tau    = ((abs(qp/q0))**q_exp)*abs(tau0)

       arg    = -t/tau

       if(arg .lt. -50.0) arg = -50.0
       

       th40 = amp * exp(arg)
c
       return
       end
