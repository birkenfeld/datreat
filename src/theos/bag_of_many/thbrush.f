c*ds
c*ed
       subroutine thbrush( q, tau, r0, sigmar0, rho_innen, ratio,
     *                     k_max, ll_max, iflcorr, ninter,
     *                     aiqt, aiqt0, sq, ier)
c
c **** input parameter:
c      q ......... skalierter impulsuebertrag (--> schichtdicke=1)
c      tau ....... skalierte zeit
c      r0 ........ skalierte dicke der festen innenkugel
c      sigmar0 ... verteilungsbreite (gauss) der innenkugel
c      rho_innen . relative streulaengendichte der innenkugel
c      ratio ..... verhaeltnis s(q) zu inelastish
c      k_max ..... hoechster eigenwert der radialen loesung
c      ll_max .... hoechstes beruecksichtigtes l
c      iflcorr ... flag: benutzung von korrekturen fuer l<>0
c      ninter .... anzahl stuetzstellen fuer eigenfunktionen
c **** output parameter:
c      aiqt ...... i(q,tau)
c      aiqt0...... i(q,tau=0)
c      sq ........ i(q,tau=+infinity)
c      ier ....... fehlerindikator
c==================================================================
c
c
       implicit     real*8 (a-h,o-z)
       character*8  filna, chrnxt, chrval
       logical      found, newsol, newsoq, newfl, newcof
 
       parameter    (pi=3.141592654d0)
 
       parameter    (mint=10)
       parameter    (mint2=((mint+1)*mint)/2)
       dimension    cint(mint2),aint(mint+1),bint(mint+1)
 
       parameter (nmax=300,kdom=50,llmax=30,kmax=10+kdom)
       dimension    x(nmax,kmax),y(nmax,kmax),elams(-2:kdom),n(kmax)
       common /creport/ x,y,elams,n,i,ireport
 
       dimension    ipermu(nmax)
       integer      kflag
 
c ******************************************************************
c * speicher fuer korrekturfaktoren fuer hoehere l-werte           *
c ******************************************************************
       dimension fl(0:llmax,0:kmax)
       dimension fe(0:llmax,0:kmax)
       common/cflfe/ fl, fe
 
c ******************************************************************
c * save-commons                                                   *
c ******************************************************************
       common /isave/kmold1,ninterold,iflcold,
     *               kmold3,llmold3
       common /rsave/r0old1,qold2,r0old2,sigmaold,rhoold,soq,
     *               qold3,told4,aiqti,aiqt0i
       common /asave/afakv(0:llmax,0:kmax),tau_i(0:llmax,0:kmax)
 
c ******************************************************************
c * datenfelder fuer nag-sturm-lioville dgl-solver                 *
c ******************************************************************
       parameter(mbrkp=20)
       dimension xpoint(mbrkp)
       dimension hmax(2,mbrkp)
c ******************************************************************
c * common kommunikation mit externen funktionen des dgl-solvers   *
c ******************************************************************
       common /coefpa/ eta
       common /cmonit/ infomo
c ******************************************************************
c * external functions fuer             dgl-solver                 *
c ******************************************************************
       external coeffn
       external bdyval
       external monit
       external report
 
 
c ******************************************************************
c * sektion 1 initialisieren der dgl-loesungen                     *
c ******************************************************************
c
c ueberspringen dieser sektion falls die parameter nicht
c geaendert wurden.
c
       newsol = .false.
 
       if( k_max    .gt.   kmold1     .or.
     *     ninterp  .ne.   ninterold  .or.
     *     r0       .ne.   r0old1         ) then
 
        write(6,*)'k_max, km_old     ',   k_max,   kmold1
        write(6,*)'ninterp .....     ',   ninterp, ninterold
        write(6,*)'r0 ..........     ',   r0, r0old1
 
         kmold1    = k_max
         ninterold = ninterp
         r0old1    = r0
         newsol    = .true.
 
c ******************************************************************
c * initialisationen   fuer             dgl-solver                 *
c ******************************************************************
       tol       = 1.0d-4
       elams(-2) = 10.0d0
       elams(-1) = 18.0d0
 
       do i=1,2
        do j=1,mbrkp
         hmax(i,j) = 0.0d0
        enddo
       enddo
 
c --- setup xpoint ---
         xpoint(1 ) = r0
         xpoint(2 ) = r0
         xpoint(3 ) = r0+0.10d0
         xpoint(4 ) = r0+0.50d0
         xpoint(5 ) = r0+0.90d0
         xpoint(6 ) = r0+0.95d0
         xpoint(7 ) = r0+0.99d0
         xpoint(8 ) = r0+1.00d0
         xpoint(9 ) = r0+1.00d0
         nbrkp      = 9
         match      = 4
 
       maxit     = 0
       maxfun    = 0
       infomo    = iout()
       ireport   = 0
       eta       = 1.0d0
       nint      = 4
       ninterp   = ninter
       if(ninterp.le.10) then
          ninterp = 50
          write(6,*)'thbrush: set ninter to ',ninterp
       endif
       c_exp     = -4.0d0/3.0d0
c --- check der dimensionen ---
       if(k_max.gt.kdom) then
          write(6,*)'thbrush: k_max is out of range:',k_max,'>',kdom
          ier = 1
          return
       endif
       if(ll_max.gt.llmax) then
          write(6,*)'thbrush: ll_max is out of range:',ll_max,'>',llmax
          ier = 2
          return
       endif
c ******************************************************************
 
 
c      ===============================================
c      =  brush-dgl fuer radiale geometrie (star)    =
c      =  grad e(r) div u                            =
c      =  u = loesung / r**2 --> report              =
c      ===============================================
 
       do kelam=0,k_max
 
         n(1) = 0
         delam = elams(kelam-1) - elams(kelam-2)
         elam  = elams(kelam-1) + delam*2.0d0
         call d02kef( xpoint,nbrkp,match,coeffn,bdyval,kelam,tol,
     *                elam,delam,hmax,maxit,maxfun,
     *                monit, report, ifail)
        if(ifail.ne.0)write(6,*)'d02kef-->ifail = ',ifail
        if(iout().ge.0) then
         write(6,6003) kelam,elam,delam
6003     format(' lambda(',i3,') = ',f14.8,' +- ',f14.8,' ..!')
        endif
c
ccc      write(6,*)'..... sorting (1)                -->(2)'
         n(2)  = n(1)
         ipath = 2
         do i=1,n(1)
            ipermu(i) = i
         enddo
 !        call dsvrgp(n(1),x(1,1),x(1,2),ipermu)
 !        call dpermu(n(1),y(1,1),ipermu,ipath,y(1,2))

!! new !!
         KFLAG = 1
         call DPSORT (X(1,1), N(1), IPERMU, KFLAG, IER)    ! replaces dsvrgp
         x(1:N(1),2) = x(ipermu(1:N(1)),1)                 ! modern fortran can do the permutation directly
         y(1:N(1),2) = y(ipermu(1:N(1)),1)                 ! so dpermu is no longer needed


c
c        write(6,*)'..... remove double x-values (2) -->(3)'
         nn = 1
         x(1,3) = 0.0d0
         y(1,3) = 0.0d0
         do i=1,n(2)
           if(x(nn,3).ne.x(i,2)) then
              nn = nn + 1
              x(nn,3) = x(i,2)
              y(nn,3) = y(i,2)
           endif
         enddo
         n(3) = nn
c
c        write(6,*)'..... interpolate (3)            -->(4)'
c        write(6,*)'..... (4) * c(z)                 -->(5)'
         if(nint.gt.mint)     nint = mint
         if(ninterp.gt.nmax)  ninterp = nmax
         nint1 = nint + 1
         nint2 = (nint*nint1)/2
         do i=1,nmax
          y(i,4)   = 0.d0
ccc       ys(i,4)  = 0.d0
         enddo
         n(4) = ninterp
         dx   = 1.0d0/(ninterp-1)
         is   = 1
ccc      xinteg = 0.0d0
         do i=1,ninterp
           xh     = dx*(i-1) + r0
           x(i,4) = xh
           x(i,5) = xh
ccc        write(6,*)i,xh
           if(xh.lt.x(1,3)) then
              j = 1
              goto 1010
           endif
           do j=is,n(3)-1
             if((x(j,3).le.xh) .and. (x(j+1,3).gt.xh) ) then
ccc            write(6,'(1x,2i3,3f16.6)')j,is,x(j,3),xh,x(j+1,3)
               goto 1010
             endif
           enddo
1010       continue
           is = j
           i0 = is-nint/2
           if(i0 .lt. 1)       i0 = 1
           if(i0 .gt. n(3)-nint-1)  i0 = n(3)-nint-1
           do j=0,nint
             aint(j+1) = x(i0+j,3)
             bint(j+1) = y(i0+j,3)
           enddo
           call e01aaf(aint,bint,cint,nint1,nint2,nint,xh)
           y(i,4) = cint(nint2)
           if(xh.gt.0.d0) then
              y(i,5) = y(i,4) / (xh**(4.0d0/3.0d0))
c --- multipliziere mit fermi-abschneidefunktion -----
ccc           arg    = (rfermi - xh) / sfermi
ccc           y(i,5) = y(i,5) / ( 1.0d0 + dexp(-arg))
ccc           xinteg = xinteg + y(i,4)**2
           else
              y(i,5) = 0.0d0
           endif
         enddo
         n(5) = n(4)
ccc      xinteg = xinteg * dx
ccc      write(6,*)'....> integral( u(z)  **2 dz ) = ',xinteg
ccc     write(6,*)'===================================================='
c
c --- speichern der werte -----------
         elams(kelam) = elam
         kelam_u = 10 + kelam
         n(kelam_u) = n(5)
         do i=1,n(5)
          x(i,kelam_u) = x(i,5)
          y(i,kelam_u) = y(i,5)
         enddo
       enddo
      endif
c    -----<--- init-if <---
 
 
 
c ******************************************************************
c * berechnung von s(q) falls erforderlich                         *
c ******************************************************************
       newsoq = .false.
 
       if( q         .ne.  qold2    .or.
     *     r0        .ne.  r0old2   .or.
     *     sigmar0   .ne.  sigmaold .or.
     *     rho_innen .ne.  rhoold    ) then
 
           qold2     = q
           r0old2    = r0
           sigmaold  = sigmar0
           rhoold    = rho_innen
           newsoq    = .true.
 
           c_exp     = -4.0d0/3.0d0
           soq    = s_of_q(q,r0,sigmar0,rho_innen,c_exp)
 
           if(iout().ge.0) then
             write(6,*)'thbrush new s(q=',q,')=',soq
           endif
 
       endif
 
 
 
c ******************************************************************
c * berechnung der koeffizienten fuer i(q,t) inelastisch           *
c ******************************************************************
       newfl = .false.
c --- a. muss fl auf 1 initialisiert werden ?? ---
       if( fl(0,0).ne.1.0d0 .and. iflcorr.eq. 0) then
         iflcold = iflcorr
         newfl   = .true.
         write(6,*)'init fl, fe to 1 !'
            do i=0,kmax
              do j=0,llmax
                fl(i,j) = 1.0d0
                fe(i,j) = 1.0d0
              enddo
            enddo
       endif
c --- b. muss fl eingelsesn werden? ---------
       if(iflcorr.ne.iflcold) then
           iflcold = iflcorr
           newfl   = .true.
 
         do i=0,kmax
           do j=0,llmax
             fl(i,j) = 1.0d0
             fe(i,j) = 1.0d0
           enddo
         enddo
 
         filna = 'fldat'
         open(66,file=filna)
         read(66,*) kmi
         lli = 0
6601     continue
         read(66,*,end=6602)(fl(lli,j),j=0,kmi)
         read(66,*,end=6602)(fe(lli,j),j=0,kmi)
         lli = lli+1
         if(lli.gt.llmax) goto 6602
         goto 6601
6602     continue
         close(66)
cc       write(6,*)'kmi=',kmi,'  lli=',lli,'  llmax=',llmax
         do j=0,kmi
           do l=lli,llmax
              fl(l,j) = fl(lli-1,j)
              fe(l,j) = fe(lli-1,j)
           enddo
         enddo
6604     continue
         do l=0,llmax
          write(6,6603)l,(fl(l,j),j=0,5),(fe(l,j),j=0,5)
6603      format(' l=',i2,1x,6f7.2,' ....'/
     *           '   ',2x,1x,6f7.2)
         enddo
       endif
c      ------ of section b. of fl-inits ---
 
c --------> so jetzt kommen die exp(-tau/tau_i) koeffizientenarrays ---
       newcof = .false.
 
       if( k_max     .ne. kmold3    .or.
     *     ll_max    .ne. llmold3   .or.
     *     newsol                   .or.
     *     newfl                    .or.
     *     q         .ne. qold3         ) then
 
           kmold3   = k_max
           llmold3  = ll_max
           qold3    = q
           newcof   = .true.
 
            do j=0,k_max
                do l = 0,ll_max
                  fcorrl = fl(0,j)/fl(l,j)
                  fcorre = fe(0,j)/fe(l,j)
                  dx     = x(2,10+j)-x(1,10+j)
                  ffq    = bsj_int(q,l,dx,y(1,10+j),n(10+j),r0)
cccccc            sffq   = sffq + (2*l+1)*ffq*ffq
                  sffq   =        (2*l+1)*ffq*ffq
cccccc          enddo
                  afak   = elams(0) / elams(j)
                  afak   = afak*sffq*4.0d0*pi * fcorre * fcorrl
                  afakv(l,j) = afak
                  tau_i(l,j) = elams(0) * fcorrl / elams(j)
                enddo
            enddo
            if(iout().ge.0) then
              write(6,*)'thbrush new afak and tau- coeffs'
            endif
 
       endif
c      -----> if init sektion
 
c *********************************************************************
c  berechnungssektion
c *********************************************************************
       if( newcof             .or.
     *     tau    .ne. told4      ) then
           told4  = tau
           aiqti  = 0.0d0
           aiqt0i = 0.0d0
           do j=0,k_max
            do l=0,ll_max
               etau   = dexp(-tau / tau_i(l,j) )
               aiqti  = aiqti   + afakv(l,j) * etau
               aiqt0i = aiqt0i  + afakv(l,j)
            enddo
           enddo
           if(iout().gt.0) then
             write(6,*)'thbrush(t=',tau,')-->',aiqti,aiqt0i
           endif
       endif
c
c -------------- final results ---------------------------
       ier     = 0
       sq      = soq
       aiqt    = aiqti  + ratio * soq
       aiqt0   = aiqt0i + ratio * soq
 
       return
       end
 
 
 
c ******************************************************************
c * external sbrs    nag-sturm-lioville dgl-solver                 *
c ******************************************************************
      subroutine coeffn(p,q,dqdl,x,elam,jint)
c     =======================================
      implicit real*8 (a-h,o-z)
      common /coefpa/ eta
 
c ----- die funktion q(x,lambda) ----
      qq= eta / (x**4)
      q = elam * qq
c ----- die ableitung von dq/dlambda ----
      dqdl = qq
c ----- die funktion p(x) -----
      entry pfp(x,p)
      p = 1.d0 / (x**5)
 
      return
      end
c
c
      subroutine bdyval(xl,xr,elam,yl,yr)
c     =================================== boundary values
      implicit real*8 (a-h,o-z)
      dimension yl(3), yr(3)
c
      yl(1) = 0.d0
      yl(2) = 1.0d0
c
      yr(1) = 1.0d0
      yr(2) = 0.0d0
 
      return
      end
c
c
      subroutine monit(maxit,iflag,elam,finfo)
c     ========================================
      implicit real*8 (a-h,o-z)
      dimension finfo(15)
      common/cmonit/infomo
 
      if(infomo.le.0) return
 
      if(iflag.lt.0) then
       write(6,*)'monit iflag=',iflag,' error exit will follow'
       goto 900
      endif
 
      if(iflag.eq.1)
     * write(6,*)'monit iflag=',iflag,' trying to bracket lambda =',
     *            elam
 
      if(iflag.eq.2)
     * write(6,*)'monit iflag=',iflag,' converging to lambda =',
     *            elam
 
      if(infomo.le.1) return
900   continue
c     finfo dump
      do i=1,15
        write(6,*)'monit finfo(',i,') = ',finfo(i)
      enddo
 
      return
      end
 
      subroutine report( xp,vp,jint )
c     ===============================
      implicit real*8 (a-h,o-z)
      dimension vp(3)
 
      parameter  (nmax=300,kdom=50,llmax=30,kmax=10+kdom)
      dimension    x(nmax,kmax),y(nmax,kmax),elams(-2:kdom),n(kmax)
      common /creport/ x,y,elams,n,i,ireport
 
      yp = dexp(vp(3)*0.5d0)*dsin(vp(2)*0.5d0)/dsqrt(vp(1))
      call pfp(xp,px)
      yd = dsqrt(vp(1))*dexp(vp(3)*0.5d0)*dcos(vp(2)*0.5d0)/px
 
      if(ireport.gt.0) write(6,6001)n(1),xp,yp,yd
6001  format(' report(',i3,'): ',3f16.7)
 
      if(n(1).ge.nmax) then
        write(6,*)'report no of x points is exhausted max=',nmax
        return
      endif
 
      n(1) = n(1)+1
      x(n(1),1)  = xp
      y(n(1),1)  = yp / (xp**2)
cccc  ys(n(1),1) = yd
 
      return
      end
 
!! c*ds                                                                    th121140
!! c*ed                                                                    th121150
!! c                                                                       th121160
!! c ----- spherical besselfunction of order m>=0 -------------------------th121170
!! c                                                                       th121180
!!        function bsjn(m,z)                                               th121190
!! c      ==================                                               th121200
!! c                                                                       th121210
!!        implicit real*8 (a-h,o-z)                                        th121220
!! c                                                                       th121230
!!        data eps/1d-6/                                                   th121240
!! c                                                                       th121250
!! c --- look for bad m value ---                                          th121260
!!        if(m.lt.0) then                                                  th121270
!!         write(6,*)' error: bsjn called with neg m =',m,' bsjn = 0'      th121280
!!         bsjn = 0                                                        th121290
!!         return                                                          th121300
!!        endif                                                            th121310
!! c                                                                       th121320
!!        if(z.eq.0) then                                                  th121330
!!          bsjn = 1.d0                                                    th121340
!!          if(m.gt.0) bsjn = 0                                            th121350
!!          return                                                         th121360
!!        endif                                                            th121370
!! c --- prepare factorial for limiting formula ---                        th121380
!!        f = 1d0                                                          th121390
!!        do 101 i=1,2*m+1,2                                               th121400
!!          f = f * i                                                      th121410
!! 101    continue                                                         th121420
!! c                                                                       th121430
!! c ---- treat case z small ----                                          th121440
!!        zzt = (z**m) / f                                                 th121450
!!        if(dabs(zzt).lt.eps) then                                        th121460
!! c        -----> use limiting formula                                    th121470
!!            bsjn = (z**m / f ) * (1d0-z**2/(2*(2*m+3)) +                 th121480
!!      *                           z**4/(8*(2*m+3)*(2*m+5)) )             th121490
!! c        write(6,*)' asymptotic formula'                                th121500
!!          return                                                         th121510
!!        endif                                                            th121520
!! c                                                                       th121530
!! c --- treat cases m=0 and m=1 ---                                       th121540
!!        if(m.lt.2) then                                                  th121550
!!         if(m.eq.0) then                                                 th121560
!!           bsjn = dsin(z) / z                                            th121570
!!           return                                                        th121580
!!         endif                                                           th121590
!!         if(m.eq.1) then                                                 th121600
!!           bsjn = ( dsin(z) / z  - dcos(z) ) / z                         th121610
!!           return                                                        th121620
!!         endif                                                           th121630
!!       endif                                                             th121640
!! c                                                                       th121650
!! c --- treat cases m >= 2 ---                                            th121660
!! c                                                                       th121670
!!       fnm1 = 1.d0/z                                                     th121680
!!       fn   = fnm1**2                                                    th121690
!!       fmn  = 0                                                          th121700
!!       fmnp1= fnm1                                                       th121710
!!       maxm = m - 1                                                      th121720
!! c                                                                       th121730
!!       do 1 n=1,maxm                                                     th121740
!!         fnp1   = (2*n+1) * fn/z - fnm1                                  th121750
!!         fmnm1  = (1-2*n) * fmn/z- fmnp1                                 th121760
!!         fnm1   = fn                                                     th121770
!!         fn     = fnp1                                                   th121780
!!         fmnp1  = fmn                                                    th121790
!!         fmn    = fmnm1                                                  th121800
!! 1     continue                                                          th121810
!!       fmnm1 = (1-2*m) * fmn/z - fmnp1                                   th121820
!!       bsjn  = fn*dsin(z) - (-1)**m * fmnm1*dcos(z)                      th121830
!! c                                                                       th121840
!!       return                                                            th121850
!!       end                                                               th121860
!!  
!!  
!! c*ds                                                                    th121140
c*ed                                                                    th121150
      function bsj_int(q,l,dr,y,ny,r0)
c     ------------------------------------------------------------
c     berechnet :
c
c          integral(r0..r0+1) j_l(qr)*(2*r*y(r)+r*r*y'(r))*dr
c
c     ------------------------------------------------------------
 
      implicit real*8 (a-h,o-z)
      dimension y(ny)
 
ccc   if(q.gt.30.d0) then
ccc     bsj_int = 0.d0
ccc     return
ccc   endif
 
      sum = 0.d0
      do i=1,ny
         r = r0 + (i-1)*dr
         if(i.lt.ny .and. i.gt.1) then
           deriv = (y(i+1)-y(i-1))/(2.0d0*dr)
         else
           if(i.eq.1)   deriv = y(1)/dr
           if(i.eq.ny)  deriv = -y(ny)/dr
         endif
         sum = sum + bsjn(l,q*r)*(2.0d0*r*y(i)+r*r*deriv)
      enddo
      bsj_int = sum*dr
 
      return
      end
 
 
       function s_of_q(q,r0,sigmar0,rho_innen,exp_aussen)
c      ==================================================
       implicit real*8 (a-h,o-z)
       parameter(epsilon =1.0d-5)
       parameter(maxiter =100)
cccc   parameter(ninteg  =25)
       parameter(ninteg  =2 )
       parameter(asigma  =3.d0)
       parameter(qmin    =1.0d-6)
       parameter(rmin0   =0.1d0)
       parameter(pi      =3.141592654d0)
       common /csint/expo,qq
       external sint
 
       expo  = exp_aussen
       if(q.lt.qmin) then
          qq = qmin
       else
          qq = q
       endif
 
       rmin   = r0-asigma*sigmar0
       rmax   = r0+asigma*sigmar0
       if(rmin.lt.rmin0*r0) rmin = rmin0*r0
       dr     = (rmax-rmin)/ninteg
       weight = 0.0d0
       soq    = 0.0d0
 
       do i=0,ninteg
         r0h    = rmin + i*dr
         r1     = r0h+1.0d0
         qr0    = r0h * qq
         sext   = adapint(sint,r0h,r1,epsilon,maxiter,erroraccu)
         sinn   = rho_innen * (dsin(qr0)-qr0*dcos(qr0))/(qq**3)
         wexp   = dexp(-( ((r0h-r0)/sigmar0 )**2 )  )
         weight = weight + wexp
cccc     write(6,*)wexp,weight, soq
         soq    = soq + wexp * (sext+sinn)**2
       enddo
cccc   write(6,*)'weight = ',weight
cccc   write(6,*)'soq    = ',soq
       soq = soq / weight
c -- zusammen:
       s_of_q = soq * (4*pi)**2
 
       return
       end
 
       function sint(r)
c      ================ integrandfunktion fuer s_of_q
       implicit real*8 (a-h,o-z)
       common /csint/expo,qq
 
       qr   = qq*r
       sint = r**(expo+2.0d0) * dsin(qr)/qr
 
       return
       end


!! Older Version used Library calls that must be replaced
!! 
!! ------------------------------------------------------------------------
!! IMSL Name:  PERMU/DPERMU (Single/Double precision version)
!!  
!! Revised:    March 19, 1991
!!  
!! Purpose:    Rearrange the elements of an array as specified by a
!!             permutation.
!!  
!! Usage:      CALL PERMU (N, X, IPERMU, IPATH, XPERMU)
!!  
!! Arguments:
!!    N      - Length of the arrays X and XPERMU.  (Input)
!!    X      - Real vector of length N containing the array to be
!!             permuted.  (Input)
!!    IPERMU - Integer vector of length N containing a permutation
!!             IPERMU(1), ..., IPERMU(N) of the integers 1, ..., N.
!!             (Input)
!!    IPATH  - Integer flag.  (Input)
!!             IPATH = 1 means IPERMU represents a forward permutation,
!!             i.e., X(IPERMU(I)) is moved to XPERMU(I).  IPATH = 2
!!             means IPERMU represents a backward permutation, i.e.,
!!             X(I) is moved to XPERMU(IPERMU(I)).
!!    XPERMU - Real vector of length N containing the array X
!!             permuted.  (Output)
!!             If X is not needed, X and XPERMU can share the same
!!             storage locations.
!!  
!! Keywords:   Utilities; Forward permutation; Backward permutation
!!  
!! GAMS:       N8
!!  
!! Chapters:   MATH/LIBRARY Utilities
!!             STAT/LIBRARY Utilities
!!  
!! Page No.:   MATH/LIBRARY User's Manual page 1274
!!             STAT/LIBRARY User's Manual page 1387
!!  
!! !------------------------------------------------
!! 
!! ------------------------------------------------------------------------
!! IMSL Name:  SVRGP/DSVRGP (Single/Double precision version)
!!  
!! Revised:    May 30, 1991
!!  
!! Purpose:    Sort a real array by algebraically increasing value and
!!             return the permutation that rearranges the array.
!!  
!! Usage:      CALL SVRGP (N, RA, RB, IPERM)
!!  
!! Arguments:
!!    N      - Number of elements in the array to be sorted.  (Input)
!!    RA     - Vector of length N containing the array to be sorted.
!!             (Input)
!!    RB     - Vector of length N containing the sorted array.
!!             (Output)
!!             If RA is not needed, RA and RB can share the same
!!             storage locations.
!!    IPERM  - Vector of length N.  (Input/Output)
!!             On input, IPERM should be initialized to the values 1,
!!             2, ..., N.  On output, IPERM contains a record of
!!             permutations made on the vector RA.
!!  
!! Remark:
!!    For wider applicability, integers (1,2,...,N) that are to be
!!    associated with RA(I) for I = 1, 2, ..., N may be entered into
!!    IPERM(I) in any order.  Note that these integers must be unique.
!!  
!! Keyword:    Utilities
!!  
!! GAMS:       N6a1b; N6a2b
!!  
!! Chapters:   MATH/LIBRARY Utilities
!!             STAT/LIBRARY Utilities
!!  
!! Page No.:   MATH/LIBRARY User's Manual page 1280
!!             STAT/LIBRARY User's Manual page 1399
!!  
!! ----------------------------------------------------------------
!! NAG: E01AAF interpolates at a given point x from a table of function values yi 
!! evaluated at equidistant or non-equidistant points xi, for i=1,2,í°ƒ,n+1,  
!! using Aitken's technique of successive linear interpolations.2  Specification
!! SUBROUTINE E01AAF (A, B, C, N1, N2, N, X)
!! INTEGERN1, N2, N
!! double precisionA(N1), B(N1), C(N2), X3  
!! Description
!! E01AAF interpolates at a given point x from a table of values xi and yi, 
!! for  i=1,2,í°ƒ,n+1 using Aitken's method.  
!! The intermediate values of linear interpolations are stored to enable 
!! an estimate of the accuracy of the results to be made.4  ReferencesFrÅˆberg C E (1970)  Introduction to Numerical Analysis Addiso 5  Parameters
!! 1: A(N1) double precision arrayInput/OutputOn entry: Ai must contain the x-component of the ith data point, xi, for i=1,2,í°ƒ,n+1.
!! On exit: Ai contains the value xi-x, for i=1,2,í°ƒ,n+1.
!! 2: B(N1) double precision arrayInput/OutputOn entry: Bi must contain the y-component (function value) of the ith data point, yi, for i=1,2,í°ƒ,n+1.
!! On exit: the contents of B are unspecified.3: C(N2) double precision array OutputOn exit: 
!! C1,í°ƒ,Cn contain the first set of linear interpolations,Cn+1,í°ƒ,C2Å◊n-1 contain the second set of linear interpolations CnÅ◊n+1/2 contains the interpolated function value at the point x.
!! 4:   N1 INTEGERInputOn entry: 
!! 
!! the value n+1 where n is the number of intervals; that is, N1 is the number of data points.5:  N2 INTEGERInput On entry: 
!! 
!! the value nÅ◊n+1/2 where n is the number of intervals.6: N  INTEGERInputOn entry: the number of intervals which are to be used in interpolating the value at x; that is, there are n+1 data points xi,yi.7: X double precisionInputOn entry: the point x at which the interpolation is required.6  Error Indicators and Warnings
!! None.
!! 


!! SLATEC

C      SUBROUTINE DPSORT (DX, N, IPERM, KFLAG, IER)
C***BEGIN PROLOGUE  DPSORT
C***PURPOSE  Return the permutation vector generated by sorting a given
C            array and, optionally, rearrange the elements of the array.
C            The array may be sorted in increasing or decreasing order.
C            A slightly modified quicksort algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A1B, N6A2B
C***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
C***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
C***AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DPSORT returns the permutation vector IPERM generated by sorting
C   the array DX and, optionally, rearranges the values in DX.  DX may
C   be sorted in increasing or decreasing order.  A slightly modified
C   quicksort algorithm is used.
C
C   IPERM is such that DX(IPERM(I)) is the Ith value in the
C   rearrangement of DX.  IPERM may be applied to another array by
C   calling IPPERM, SPPERM, DPPERM or HPPERM.
C
C   The main difference between DPSORT and its active sorting equivalent
C   DSORT is that the data are referenced indirectly rather than
C   directly.  Therefore, DPSORT should require approximately twice as
C   long to execute as DSORT.  However, DPSORT is more general.
C
C   Description of Parameters
C      DX - input/output -- double precision array of values to be
C           sorted.  If ABS(KFLAG) = 2, then the values in DX will be
C           rearranged on output; otherwise, they are unchanged.
C      N  - input -- number of values in array DX to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the value in the original order of the
C              DX array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting DX in increasing order and sort DX also.
C            =  1  means return the permutation vector resulting from
C                  sorting DX in increasing order and do not sort DX.
C            = -1  means return the permutation vector resulting from
C                  sorting DX in decreasing order and do not sort DX.
C            = -2  means return the permutation vector resulting from
C                  sorting DX in decreasing order and sort DX also.
C      IER - output -- error indicator:
C          =  0  if no error,
C          =  1  if N is zero or negative,
C          =  2  if KFLAG is not 2, 1, -1, or -2.
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified by John A. Wisniewski to use the Singleton
C           quicksort algorithm.
C   870423  Modified by Gregory S. Rhoads for passive sorting with the
C           option for the rearrangement of the original data.
C   890619  Double precision version of SPSORT created by D. W. Lozier.
C   890620  Algorithm for rearranging the data vector corrected by R.
C           Boisvert.
C   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
C   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
C   920507  Modified by M. McClain to revise prologue text.
C   920818  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
C***END PROLOGUE  DPSORT
