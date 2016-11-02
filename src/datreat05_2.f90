! INCLUDE "commons.h"

! **********************************************************************
! *    treating any data                                               *
! *    imsl version 10 fitprogram test                                 *
! *                                                                    *
! *    datreat version from 18.5.88                                    *
! *    m. monkenbusch, iff, kfa-juelich                                *
! *    neu 3.11.93: multiply flag fuer theorien                        *
! *    neu 28.1.99: thc -n x1 <x1> x2 <x2> on log scale                *
! *    neu 09.3.99: fun zimm, rouse etc + error                        *
! *                 plot noparplo , achsen lange strings               *
! *                 plot error bars bei log ok                         *
! *                 purge loescht tatsaechlich (auch errors)           *
! *                 sel add                                            *
! *    neu 04.02.00:seterr                                             *
! *                 generiere Fehler aus Daten (Fitweighting)          *
! *                 ueber yformel                            (RCS:1.3) *
! *                                                                    *
! *    neu 23.02.00:fit map ..                                         *
! *    neu 19.05.00:dparam('name   ') fuer th-programmierung verfuegbar*
! *    neu         :msave             fuer Mehrfachsave                *
! *                                                                    *
! *                                                                    *
! *    BSS-Version datreat05_1                                         *
! *    1. sel with new options: next, all, range                       *
! *                                                                    *
! **********************************************************************
!
  program datreat
      use dimensions
      use cincom
      use cincoc
      use cmargs
      use xoutxx
      use cdata
      use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
      use constants
      use PhysicalConstantsPlus
      use unift
      implicit none

      integer iadda
      common /thiadd/iadda


       integer isigint, ifound
       real yyy, zpar, ywr, yw2, yy, ysum, ynno, yval, ymm2, ymaxf,ymax
       real yw1, yminf, ym1, yfamp, yerr, yer, yampf, y3, yctab, y11, y1
       real y00, y2, xz, y, y0, xx, xval, xsum, yinterp, xr, xnc, xmon2
       real xmonitor, xmon1, xn, xmax, xk00, fsplin, xm, ximaxf, xi, ximaxq
       real xcut, x1raster, xk0, x2raster, x, val, tt, wl4, thval, x01, x02
       real theta2, thick, temp_z, trans, t1, sump, sumb, sx, slope, smpar
       real t2, sy, sigma, selpartol, r0, selparval, ri, result, qend, qstart
       real pkave, pfkave, pfave, qz, qrs, qmax, pxxd, pkqave, qortho
       real pfqave, p3qave, p3ave, parval_x, gunten, fy, fwid, laenge, goben
       real fia, famp, facto2, eta_s, f1, errmx, facto1, errrec, errmax
       real erraec, echo, dxx, fx, smirro, fai, esum, errest, errret, df
       real delqv, ddj, delqh, dj, draster, d0, d1, d, d2, blim, bklen, a2
       real a, alam, aa, akb_z, ai0, ampli, amp, apar, alim, an0, bfave, bkgr
       real bkqave, bkave, dq, dt, dx, bpar, bfkave, aver, bqkave, detdis, rrv
       integer nw, nx, numr,  num2, num1, nra, npp, npoint, nraster
       integer npk, npay, npl, npax, np, nold, nnpar, nnn, nnumi, nn4, nn3
       integer nnw, nn1, nn, nn2, nfits, newadd, nneu, ncp, ncpl, ncut, nech
       integer nbuff, newnum, nb, n1, n2, n3, n0, nfft, n, msel, mpk, modeqc
       integer maxint, m, lyy, lxx, ll, kk, mcut, k1, k2, kz, k, jpar, jl, l
       integer j, j2, ival, kan1, ithc, iy, j1, iss, ispc, isl, ira, iseed, ir
       integer intval, ito, isel, ip, imx, imaxf, iminf, ikz, ilas, inew, ii
       integer ifunx, ifuny, ik, iery, ies, idimux, ier, iprs, ifrom, ierr
       integer ierx, iequal, ides, icut, id, ib, ileng, ic, iauto, iaddp 
       integer iaddf, iad2, iadd, iad1, ia1, ia2, ia3, ia4, ia5, ia6, i, ia
       integer itcal, ias, igrand

       real*8 getval, dble

       character :: fpnam*8,xpnam*8,ypnam*8
       character fsname*1024
       character*1 csel
       real*8 val8x,val8y

       complex ca(mwert), cb(mwert), cmplx, conjg
       dimension alim(2), blim(2)
!
       dimension nnumi(minc)
       logical :: withmo = .true.,fonpla = .false.,found

       parameter( mcut=200 )
       dimension rrv(0:mcut)

! ---> transfer der addresse der gerade bearbeiteten kurve nach thXX
!      damit koennen dann ggf. weitere parameter ueber parget gewonnen
!      werden


       double precision qortho_0, qortho_1,qz_0, qz_1
       integer          nz, northo, iortho, iz, ith
       character*40     filnam
       character*8      filext, chrval
       character*8      selparna

                                      ! dir-command --> comment-length
       integer          len_comm

       character*420    pastring
       character*80     fostring
       character*8      dirpan(20)
       real*4           dirpav(20)

       real          :: y_scaling
       real          :: diffscal, diffrouse



       character(len=80) :: cbuffer


!
       external fdes
       external f
!
!
       real :: errabs=1.0,errel=1.e-2
       real :: xcopy1, xcopy2

!
! lokale Hilfsvariablen fuer arit2 (geht noch besser) !
!
       integer :: nout
       integer, parameter :: mwork=c_MWERT
       real :: yout1(mwork), yerout1(mwork), yout2(mwork), yerout2(mwork)



       real    :: x1int, x2int, sum, sumer
       integer :: i1int, i2int

!
       double precision :: Ommax = 1000d0
       integer          :: nomega = 100

       real             :: He3_pressure_bar = 4.75  ! (IN5 value in 2013)
       real             :: Tube_diameter_m  = 0.025 ! (IN5 1 inch tubes)
       real             :: lenfrac          = 1.0
       real             :: detector_sensitivity, detsens, lambdaA


       double precision :: unift_range_expand     = 1d0
       double precision :: unift_resolution_limit = 0.1d0
       

         real  :: lambda0, lambda1, delta_lambda, temp,angle_2tht 
!?         real  :: e0, e1, dE, de0, dee, dbb, omega, om_cm, qrc, qa, q
         real  :: e0, e1, dE, de0, dee, dbb, omega, om_cm, qa, q
         real  :: kinetic_factor, channel_factor,channel_width,channel_width0
         real  :: sume, sumq
  
       integer :: irstp, recstep = 1   


       real    ::   lower_range 
       real    ::   upper_range
       logical ::   range_is_y  
 

! ---- initialisations ----
! ---- error-set ----------
       call erset (0,1,0)

                        !! activates Crtl-C signal handler
       call sigset(-1)
                        !! Clears Signals
       call sig_Reset()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        write(6,*)
                        write(6,*)'======================================================='
                        write(6,*)'=   datreat12_2     Version: mm-develop 2.4b          ='
                        write(6,*)'=   -----------     --------                          ='
                        write(6,*)'=   Author: M.Monkenbusch  R. Biehl, O.Holderer, JCNS ='
                        write(6,*)'======================================================='
                        prompt = "#mm-develop 2.4b -> " 
                        write(6,*)
                        write(6,*)
                        write(6,*)
                        write(6,*)'=================================================================='
                        write(6,*)'=  NEW NEW NEW NEW NEW                                           ='
                        write(6,*)'=  sel : advanced functions, type "sel help" to learn more       ='
                        write(6,*)'=  average : combine data,   type "average help"                 ='
                        write(6,*)'=  copy : copy selected, if x1 <x1> x2 <x2> are given            ='
                        write(6,*)'=         copy only that range                                   ='
                        write(6,*)'=  sequence: replace x-values by sequence nubers 1..n            ='
                        write(6,*)'=  swapxy:   exchange x and y values, discard errors             ='
                        write(6,*)'=  numorchg: change numor of selected record                     ='
                        write(6,*)'=  arit:     back to old one                                     ='
                        write(6,*)'=  arit2:    note needs sorted data...                           ='
                        write(6,*)'=  get_th:   activates theory of msave files:: gth <filenam>     ='
                        write(6,*)'=  msave:    unlimited filename length                           ='
                        write(6,*)'=  .... up to 80 theories, 40 th-params, 300 uservars            ='
                        write(6,*)'=  fit: maxfn , ngood  parameters revitalized                    ='
                        write(6,*)'=  .... more commands with help option, try <cmd> help           ='
                        write(6,*)'=  .... improved output formats in plots and lastth              ='
                        write(6,*)'=  .... observes ERROR (and stops makros in expression eval)     ='
                        write(6,*)'=  HINT: use sys to issue system commands including -,+ ...      ='
                        write(6,*)'=        do this always in makros                                ='
                        write(6,*)'=  In parameter list quotes are respected to protect strings     ='
                        write(6,*)'=  (internal incom: we only assume evaluation if 1st char is     ='
                        write(6,*)'=  .(+-1..9)                                                     =' 
                        write(6,*)'=  with fit: parameter parwght  , go, help, couple               =' 
                        write(6,*)'=================================================================='
                        write(6,*)' Pi = ',pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       pi   = 4 * atan(1.0)

       facto1 = 1
       facto2 = 2

! ---- noise amplitude ----
       ampli   = 1.0
       iseed   = 12347897
! ---- nse-echoform ----
       nech    = 64
       nra     = 10
       amp     = 10
       ncp     = 1
       mpk     = 3
       ic      = 1
       j1echo  = 100.e0
       j2      = 100.e0
       j0delta = 0.e0
       cdelta  = 1.e-6
       alam0   = 8.e0
       dalam   = 0.5e0
       tau     = 1.e-8
       dj      = 0
       ddj     = 1.358237/(4*alam0)
       famp    = 20.0
       fwid    = 2.0

       erraec  = 1.e-4
       errrec  = 1.e-4
       maxint  = 2000

       sigma   = 1d0

! ---- ecc -------------
      n = 200
      ri= 5.0
      r0= 18.9
      d0= 0.2
      a = 1.2e-3
      ncut = 7
!!!!!!!!!!
       msel = 0
       selpartol = 0.001
       selparval = 0
!!!!!!!!!!
!-------- cailtab -------------
       northo   = 11
       nz       = 11
       qortho_0 = 0.0001d0
       qortho_1 = 0.1001d0
       qz_0     = 0.01d0
       qz_1     = 0.15d0
       filext   = 'last'
!
!-------- fun Rouse, Zimm Scaling -----
       wl4      = 1.0
       eta_s    = 1e-3
                                     ! fuer A**-3ns
       akb_z    = 1e21*1.380662e-23
!

!--------------------------------------
       ! comment-length in dir (display)
       len_comm = 30

       title='Datreat Plot'

       call init_theories(thenam,thparn,nthpar,thrapar,thramin,thramax,mth,mtpar)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Anfang der Hauptschleife
 2000  continue
       call sig_Reset()
       if(ierrs.ne.0) then
         write(6,*)' Final: error return code:',ierrs
         ioldc = 0
!        ----------> if error forget rest of commandline (reslin)
         if(inka.ne.5) then
            close(inka)
            inka = 5
!           --------> go back to input from keyboard
         endif
       endif
       ierrs = 0
!       mask_err = .false.
!

       iibuf= isels(1)
       iadda= isels(1)
       mask_err = .false.   ! expression evaluation error detection is activated
                            ! since evaluation is done in an early stage of line interpretation
                            ! this error should be masked if comment lines etc. are read
                            ! e.g. in input in order to avoid erroneous error messages
                            ! when constructions that resemble expression are contained.

       call incom(comand)
!      ------------------
       iout = iot
!
!
       if(comand.eq.'?       '.or.comand.eq.'help    ') then
!                    -                       ----
!         call system('firefox '//datreat_path//'/doc/DatreatManual6.html &')
         call command_list()
         write(*,*)"Please check datreat../doc/DatreatManual6.html for more help"
         goto 2000
       endif
!
!
       if(comand.eq.'in      '.or.comand.eq.'input   ') then
!                    --                      -----
         mask_err = .true.
         call input
         nsel = 1
         isels(1) = nbuf
         ifits(1) = 0
         mask_err = .false.
         goto 2000
       endif
!
       if(comand.eq.'out_gli '.or.comand.eq.'gli     ') then
!                    --                      -----
         call out_gli
         goto 2000
       endif
!
       if(comand.eq.'inscn   ') then
!                    -----
         call inscn
         nsel = 1
         isels(1) = nbuf
         ifits(1) = 0
         goto 2000
       endif
!
!
       if(comand.eq.'cs      '.or.comand.eq.'clearsel') then
!                    --                      --------
! 2020    continue
!   <------ dieser sprung ist unelegant und nicht systemkonform
!           in der naechsten version subroutine cs verwenden !!
         nsel = 0
         do 333 i=1,msel
          isels(i) = 0
          ifits(i) = 0
  333    continue
         write(6,*)'selections are removed ....'
         goto 2000
       endif
!
!
       if(comand.eq.'z       '.or.comand.eq.'zero    ') then
!                    -                       ----
         nbuf = 0
         write(6,*)'nbuf has been resetted ....'
         nsel = 0
         do  i=1,msel
          isels(i) = 0
          ifits(i) = 0
         enddo
         write(6,*)'selections are removed ....'
       endif
!
!
       if(comand.eq.'noise   ') then
!                    -----
         if(nsel.eq.0) then
           write(6,*)'no item selected'
           goto 2000
         endif
         iadd = isels(1)
         nbuf = nbuf + 1
         if(nbuf.gt.mbuf) nbuf = mbuf
         ides = nbuf
         call txfera(iadd,ides)
         nn   = nwert(iadd)
         ampli = getval('a       ',dble(ampli),inew)
         iseed = intval('iseed   ',iseed,inew)
!cccc    call rnopt(6)
         if(inew.ne.0) then
            call rnset(iseed)
            write(6,*)'new seed=',iseed
         endif
         if(found('errors  ')) then
           do i=1,nn
             yerror(i,ides) = sqrt(ampli*ywerte(i,ides))/ampli
           enddo
         else
           do i=1,nn
            ywerte(i,ides) = igrand(ywerte(i,iadd)*ampli)
            yerror(i,ides) = sqrt(ywerte(i,ides))
           enddo
         endif
         isels(1) = ides
         nsel     = 1
         goto 2000
       endif
!
!
       if(comand.eq.'gen.res ') then
!                    -------
! ---> create a file in in11 spin-echo result format
         if(nsel.eq.0) then
           write(6,*)'no item selected'
           goto 2000
         endif
         iadd = isels(1)
         nn   = nwert(iadd)
         do i=1,8
           if(name(iadd)(i:i).eq.' ') name(iadd)(i:i)='x'
         enddo
         open(33,file=name(iadd)//'.res',status='UNKNOWN')
         ynno = ywerte(1,iadd)
         call parget('theta   ',theta2,iadd,ier)
         call parget('alam    ',alam  ,iadd,ier)
         call parget('bkgr    ',bkgr  ,iadd,ier)
         alam = alam*1.0e10
         write(33,'(a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')                &
     &   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',               &
     &   bkgr
         write(6,*)'----------------------------------'
         write(6 ,'(a,a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')              &
     &   ' file=',                                                      &
     &   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',               &
     &   bkgr
         write(6,'(1x,a)')coment(iadd)
         do i=1,nn
          write(33,*)xwerte(i,iadd),ywerte(i,iadd)/ynno,                &
     &                              yerror(i,iadd)/ynno
         enddo
         close(33)
         goto 2000
       endif
!
!

       if(comand.eq.'echocurv') then
!                    ---------
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= echocurv j <j> ......                                                       '
           write(6,*)'=    computation of simulated phase scan of NSE for a simple Lorenzian spect. '
           write(6,*)'=    ft by integration assuming gaussian distribution functions               '
           write(6,*)'=    parameters:                                                              '
           write(6,*)'=      j   <val>        : field integral                                      '
           write(6,*)'=      dj  <val>        : asymmetry (phase0)                                  '
           write(6,*)'=      j0delta <val>    : inhomogeneity offset (external fields..)            '
           write(6,*)'=      ddj <val>        : phase scan step (field integral)                    '
           write(6,*)'=      n <val>          : number of points in scan                            '
           write(6,*)'=      cdelta <val>     : relative inhomogeneity                              '
           write(6,*)'=      lambda0 <val>    : wavelength                                          '
           write(6,*)'=      dlambda <val>    : wavelength width                                    '
           write(6,*)'=      tau <val>        : relaxation time --> spectrum                        '
           write(6,*)'=      errabs <> errel <> maxint  : integration parameters                    '
           write(6,*)'=============================================================================='
        endif

         if(nbuf.lt.mbuf) then
            ia1  =  nbuf+1
            nbuf =  ia1
         else
            call errsig(100,'buffer full!$')
            goto 2000
         endif

         j1echo = getval('j       ',dble(j1echo)    ,inew)
         dj     = getval('dj      ',dble(dj )       ,inew)
         j0delta= getval('j0delta ',dble(j0delta)   ,inew)
         ddj    = getval('ddj     ',dble(ddj)       ,inew)
         nech   = intval('n       ',nech  ,inew)
         cdelta = getval('cdelta  ',dble(cdelta)    ,inew)
         alam0  = getval('lambda0 ',dble(alam0)     ,inew)
         dalam  = getval('dlambda ',dble(dalam)     ,inew)
         tau    = getval('tau     ',dble(tau  )     ,inew)
         erraec = getval('errabs  ',dble(erraec)    ,inew)
         errrec = getval('errrel  ',dble(errrec)    ,inew)
         maxint = intval('maxint  ',maxint,inew)

         call     parset('j       ',j1echo     ,ia1)
         call     parset('dj      ',dj         ,ia1)
         call     parset('j0delta ',j0delta    ,ia1)
         call     parset('ddj     ',ddj        ,ia1)
         call     parset('cdelta  ',cdelta     ,ia1)
         call     parset('lambda0 ',alam0      ,ia1)
         call     parset('dlambda ',dalam      ,ia1)
         call     parset('tau     ',tau        ,ia1)
         call     parset('errabs  ',erraec     ,ia1)
         call     parset('errrel  ',errrec     ,ia1)

         alim(1) = alam0 - 3*dalam
         blim(1) = alam0 + 3*dalam

         alim(2) = -5.d0/tau
         blim(2) = -alim(2)

         alim(1)= getval('a1      ',dble(alim(1)),inew)
         alim(2)= getval('a2      ',dble(alim(2)),inew)
         blim(1)= getval('b1      ',dble(blim(1)),inew)
         blim(2)= getval('b2      ',dble(blim(2)),inew)


         write(6,*)      'j       ',j1echo
         write(6,*)      'dj      ',dj
         write(6,*)      'ddj     ',ddj
         write(6,*)      'j0delta ',j0delta
         write(6,*)      'cdelta  ',cdelta
         write(6,*)      'lambda0 ',alam0
         write(6,*)      'dlambda ',dalam
         write(6,*)      'tau     ',tau
         write(6,*)      'errabs  ',erraec
         write(6,*)      'errrel  ',errrec
         write(6,*)      'maxint  ',maxint
         write(6,*)      'a1      ',alim(1)
         write(6,*)      'a2      ',alim(2)
         write(6,*)      'b1      ',blim(1)
         write(6,*)      'b2      ',blim(2)
         j2echo = j1echo + dj

         call  qand(f,2,alim,blim,erraec,errrec,maxint,                 &
     &              result, errest )
         write(6,*)'-------------------------------------------'
         write(6,*)'result =',result
         write(6,*)'errest =',errest
         write(6,*)'-------------------------------------------'
         n0 =     nech/2
         n0     = intval('n0      ',n0    ,inew)
                                                       !! aix.sp
         call     parset('n0      ',(real(n0)) ,ia1)

         nwert(ia1) = nech
         do i=1,nech
           j2echo = j1echo + (i-n0) * ddj + dj
           xwerte(i,ia1) = j2echo-j1echo
! --- imsl mehrdimensionale integration ---
           call  qand(f,2,alim,blim,erraec,errrec,maxint,               &
     &                result, errest )
           write(6,*)'result =',result,'     errest=',errest
           ywerte(i,ia1) = result
           yerror(i,ia1) = 0
         enddo
         xname(ia1) = 'dj/gm'
         yname(ia1) = 'counts'
         name(ia1)  = 'echo'
         coment(ia1)= 'echo'
         numor(ia1) = intval('numor   ',1,inew)
         write(6,*)'ok'
         goto 2000
       endif


       if(comand.eq.'fft     ') then
!                    -------
         ia1    = isels(1)
         call parget ('n0      ',an0     ,ia1 ,ier)
         if(ier.eq.0) then
           n0 = an0+0.1
         else
           n0 = nwert(ia1)/2
         endif
         mpk = intval('mpeak   ',mpk,inew)
         famp= getval('famp    ',dble(famp),inew)
         fwid= getval('fwid    ',dble(fwid),inew)
         write(6,*)'n0 ....... =',n0
         write(6,*)'mpeak..... =',mpk
         write(6,*)'famp ..... =',famp
         write(6,*)'fwid ..... =',fwid
         if(nbuf.le.mbuf-5) then
           ia2 = nbuf+1
           ia3 = nbuf+2
           ia4 = nbuf+3
           ia5 = nbuf+4
           ia6 = nbuf+5
           nbuf= ia6
         else
           call errsig(100,'buffer full!$')
           goto 2000
         endif

         write(6,*)'fft    ',ia1,' to  ',ia2,' phase -->',ia3
         write(6,*)'filter ',ia4,' filtered: ',ia5
         call txfera(ia1,ia2)
         call txfera(ia1,ia3)
         call txfera(ia1,ia4)
         call txfera(ia1,ia5)
         call txfera(ia1,ia6)
         n1 = nwert(ia1)
         dt = xwerte(2,ia1)-xwerte(1,ia1)
         do i=1,n1
          ca(i)    =  cmplx( ywerte(i,ia1)*dt, 0.0 )
         enddo
         call  fftcf( n1, ca, cb )
         tt = xwerte(n1,ia1) - xwerte(1,ia1)
         df = 6.283185307/tt
         do i=1,n1
          xwerte(i,ia2) = (i-1)*df
          xwerte(i,ia3) = xwerte(i,ia2)
          ywerte(i,ia2) = sqrt(real(cb(i)*conjg(cb(i))))
          ywerte(i,ia3) = atan2(real(cb(i)),aimag(cb(i)))
         enddo
! --- und filtern ---
         do i=1,n1/2
           ywerte(i,ia4) = 2*atan(famp*exp(-((i-1-n1/4)/fwid)**2))/pi
           xwerte(i,ia4) = xwerte(i,ia2)
           cb(i) = cb(i) * ywerte(i,ia4)
           xwerte(n1+1-i,ia4) = xwerte(n1+1-i,ia2)
           if(i.ge.2) then
             cb(n1+2-i)= cb(n1+2-i) * ywerte(i,ia4)
             ywerte(n1+2-i,ia4) = ywerte(i,ia4)
           endif
         enddo
         ywerte(n1/2+1,ia4) = 0
         cb(n1/2+1) = (0.0,0.0)
         call  fftcb( n1, cb, ca )
         ymaxf = 0
         yminf = 1.e30
         do i=1,n1
          ywerte(i,ia5) = real(ca(i))/(dt*n1)
          xwerte(i,ia5) = xwerte(i,ia1)
          ywerte(i,ia6) = aimag(ca(i))/(dt*n1)
          xwerte(i,ia6) = xwerte(i,ia1)
          if(ymaxf.lt.ywerte(i,ia5)) then
             ymaxf = ywerte(i,ia5)
             imaxf = i
             ymm2  = ywerte(i,ia5)-ywerte(i+2,ia5)
             ym1   = ywerte(i+1,ia5)
          endif
          if(yminf.gt.ywerte(i,ia5)) then
             yminf = ywerte(i,ia5)
             iminf = i
          endif
         enddo
         yampf = sqrt( (ymm2/2.0)**2 + ym1**2)

         name(ia2) = 'fftabs '
         name(ia3) = 'fftphas'
         name(ia4) = 'filter '
         name(ia5) = 'filtrd.r'
         name(ia6) = 'filtrd.i'
         xname(ia2)= '_1/dj*gm'
         xname(ia3)= '_1/dj*gm'
         xname(ia4)= '_1/dj*gm'
         xname(ia5)= '   dj/gm'
         xname(ia6)= '   dj/gm'
         numor(ia2)= numor(ia1)+1000
         numor(ia3)= numor(ia1)+2000
         numor(ia4)= numor(ia1)+3000
         numor(ia5)= numor(ia1)+4000
         numor(ia6)= numor(ia1)+5000
         nwert(ia4)= n1
         nwert(ia5)= n1
         nwert(ia6)= n1
! --- speziell fuer echo-auswertung ----
         nnn = nwert(ia2)
         nn1 = 3
         nn2 = nnn/4-mpk
         nn3 = nnn/4+mpk
         nn4 = nnn/2
         sumb= 0
         nb  = 0
         sump= 0
         npk = 0
         do i=nn1,nn2
          sumb = sumb + ywerte(i,ia2)**2
          nb   = nb+1
         enddo
         do i=nn3,nn4
          sumb = sumb + ywerte(i,ia2)**2
          nb   = nb+1
         enddo
         sumb = sumb/nb
         do i=nn2,nn3
           yyy  = ywerte(i,ia2)**2 - sumb
           if(yyy.ge.0.0) then
              sump = sump + sqrt(yyy)
           else
              sump = sump - sqrt(-yyy)
           endif
           npk = npk+1
         enddo
         sumb  = sqrt(sumb)
         sump  = sump*df/pi
         write(6,*)'peak ---> ',sump
         write(6,*)'noise --> ',sumb
         call parset('peak    ',sump,ia2)
         call parset('mpk     ',float(mpk),ia2)
         call parset('noisef  ',sumb,ia2)
! --- echoauswertung konventionelle methode ---
         y1 = ywerte(n0-2,ia1)
         y2 = ywerte(n0-1,ia1)
         y3 = ywerte(n0  ,ia1)
         aver = (y1+y3)*0.5
         echo = sqrt( ((y3-y1)*0.5)**2 + (y2-aver)**2 )
         write(6,*)'three point echo = ',echo
         call parset('echo3   ',echo,ia2)
         write(6,*)'filter maximum   = ',ymaxf,'   (min=',yminf,')'
         write(6,*)'filter amplitude = ',yampf,'   (min=',yminf,')'

         goto 2000
       endif


       if(comand.eq.'fftrand ') then
!                    -------
         ia1    = isels(1)
         call parget ('n0      ',an0     ,ia1 ,ier)
         if(ier.eq.0) then
           n0 = an0+0.1
         else
           n0 = nwert(ia1)/2
         endif
         if(nbuf.le.mbuf-5) then
           ia2 = nbuf+1
           ia3 = nbuf+2
           ia4 = nbuf+3
           ia5 = nbuf+4
           ia6 = nbuf+5
           nbuf= ia6
         else
           call errsig(100,'buffer full!$')
           goto 2000
         endif
         mpk = intval('mpeak   ',mpk,inew)
         amp = getval('a       ',dble(amp),inew)
         nra = intval('nrand   ',nra,inew)
         famp= getval('famp    ',dble(famp),inew)
         fwid= getval('fwid    ',dble(fwid),inew)

         write(6,*)'fft  ',ia1,' to  ',ia2,' phase -->',ia3
         write(6,*)'filter ',ia4,' filtered: ',ia5
         write(6,*)'(a)mplitude=',amp
         write(6,*)'nrand      =',nra
         write(6,*)'mpeak      =',mpk
         write(6,*)'n0         =',n0
         write(6,*)'famp       =',famp
         write(6,*)'fwid       =',fwid
         call txfera(ia1,ia2)
         call txfera(ia1,ia3)
         call txfera(ia1,ia4)
         call txfera(ia1,ia5)
         call txfera(ia1,ia6)
         n1 = nwert(ia1)
         dt = xwerte(2,ia1)-xwerte(1,ia1)
! --- speziell fuer echo-auswertung ----
         nnn = nwert(ia2)
         nn1 = 3
         nn2 = nnn/4-mpk
         nn3 = nnn/4+mpk
         nn4 = nnn/2
         do i=1,n1
          ywerte(i,ia3) = 0
          yerror(i,ia3) = 0
          ywerte(i,ia6) = 0
          yerror(i,ia6) = 0
         enddo
! ---> fuer den filter
         do i=1,n1/2
           ywerte(i,ia4) = 2*atan(famp*exp(-((i-1-n1/4)/fwid)**2))/pi
           if(i.ge.2) then
             ywerte(n1+2-i,ia4) = ywerte(i,ia4)
           endif
         enddo
         ywerte(n1/2+1,ia4) = 0
! <---
         bkave = 0
         bkqave= 0
         pkave = 0
         pkqave= 0
         p3ave = 0
         p3qave= 0
         pfave = 0
         pfqave= 0
         ximaxf = 0
         ximaxq = 0

         do ira=1,nra
           do i=1,n1
             ywr = igrand(ywerte(i,ia1)*amp)*dt
             ca(i)    =  cmplx( ywr, 0.0 )
           enddo
           call  fftcf( n1, ca, cb )
           tt = xwerte(n1,ia1) - xwerte(1,ia1)
           df = 6.283185307/tt
           do i=1,n1
            ywerte(i,ia2) = sqrt(real(cb(i)*conjg(cb(i))))
            ywerte(i,ia3) = ywerte(i,ia3) + ywerte(i,ia2)
            yerror(i,ia3) = yerror(i,ia3) + ywerte(i,ia2)**2
           enddo
! --- und filtern ---
         do i=1,n1
           cb(i) = cb(i) * ywerte(i,ia4)
         enddo
         call  fftcb( n1, cb, ca )
         ymaxf = 0
         yminf = 1.e30
         do i=1,n1
          ywerte(i,ia5) = ca(i)/(dt*n1)
          ywerte(i,ia6) = ywerte(i,ia6) + ywerte(i,ia5)
          yerror(i,ia6) = yerror(i,ia6) + ywerte(i,ia5)**2
          if(ymaxf.lt.ywerte(i,ia5)) then
             ymaxf = ywerte(i,ia5)
             imaxf = i
          endif
          if(yminf.gt.ywerte(i,ia5)) then
             yminf = ywerte(i,ia5)
             iminf = i
          endif
         enddo
         imaxf = intval('imaxf   ',imaxf,inew)
         ximaxf = ximaxf + imaxf
         ximaxq = ximaxq + imaxf**2
         yfamp =                                                        &
     &         ywerte(imaxf,ia5)**2 + ywerte(imaxf+1,ia5)**2            &
     &    +    ywerte(imaxf,ia5)**2 + ywerte(imaxf-1,ia5)**2
         yfamp = sqrt( yfamp/2.0)
!        ------> filtern: bis hier
! ---- echoauswertung ----
           sumb= 0
           nb  = 0
           sump= 0
           npk = 0
           do i=nn1,nn2
            sumb = sumb + ywerte(i,ia2)
            nb   = nb+1
           enddo
           do i=nn3,nn4
            sumb = sumb + ywerte(i,ia2)
            nb   = nb+1
           enddo
           sumb = sumb/nb
           do i=nn2,nn3
             yyy  = ywerte(i,ia2)**2 - sumb
             if(yyy.ge.0.0) then
                sump = sump + sqrt(yyy)
             else
                sump = sump - sqrt(-yyy)
             endif
             npk = npk+1
           enddo
           sump  = sump*df/pi
! --- 3-pt-konventionell ---
           y1 = igrand(ywerte(n0-2,ia1)*amp*n1/3.0)
           y2 = igrand(ywerte(n0-1,ia1)*amp*n1/3.0)
           y3 = igrand(ywerte(n0,ia1)*amp*n1/3.0 )
           aver = (y1+y3)*0.5
           echo = sqrt( ((y3-y1)*0.5)**2 + (y2-aver)**2 )
           if(iout.gt.0)write(6,*)ira,'  peak:',sump,'  noise:',sumb,   &
     &                            ' peakf:',yfamp
!
           bkave = bkave + sumb
           bkqave= bkqave+ sumb**2
           pkave = pkave + sump
           pkqave= pkqave+ sump**2
           p3ave = p3ave + echo
           p3qave= p3qave+ echo**2
           pfave = pfave + yfamp
           pfqave= pfqave+ yfamp**2
!
         enddo
         bkave = bkave/nra
         pkave = pkave/nra
         p3ave = p3ave/nra
         pfave = pfave/nra
         bkqave= bkqave/nra
         pkqave= pkqave/nra
         p3qave= p3qave/nra
         pfqave= pfqave/nra
         ximaxf = ximaxf/nra
         ximaxq = ximaxq/nra
         bkqave= sqrt(bkqave - bkave**2)
         pkqave= sqrt(pkqave - pkave**2)
         p3qave= sqrt(p3qave - p3ave**2)*3.0/n1
         p3ave = p3ave *3.0/n1
         pfqave= sqrt(pfqave - pfave**2)
         ximaxq= sqrt(ximaxq-ximaxf**2)
         write(6,*)'noise average=',bkave,'  deviation=',bkqave
         write(6,*)'peak  average=',pkave,'  deviation=',pkqave
         write(6,*)'peak3 average=',p3ave,'  deviation=',p3qave
         write(6,*)'peakf average=',pfave,'  deviation=',pfqave
         write(6,*)'peakf posit. =',ximaxf,'  deviation=',ximaxq
         call parset('peak    ',pkave,ia3)
         call parset('peak3   ',p3ave,ia3)
         call parset('peakdev ',pkqave,ia3)
         call parset('peak3dev',p3qave,ia3)
         call parset('mpk     ',float(mpk),ia3)
         call parset('noisef  ',bkave,ia3)
         call parset('noisedev',bqkave,ia3)
         call parset('peakf   ',pfave,ia5)
         call parset('peakfdev',pfkave,ia5)
         call parset('famp    ',pfave,ia5)
         call parset('fwid    ',bfkave,ia5)
         call parset('famp    ',bfave,ia4)
         call parset('fwid    ',pfkave,ia4)
         call parset('peakf   ',pfave,ia6)
         call parset('peakfdev',pfkave,ia6)
         call parset('famp    ',pfave,ia6)
         call parset('fwid    ',pfkave,ia6)


         do i=1,n1
            xwerte(i,ia2) = (i-1)*df
            xwerte(i,ia3) = xwerte(i,ia2)
            xwerte(i,ia4) = xwerte(i,ia2)
            xwerte(i,ia5) = xwerte(i,ia1)
            xwerte(i,ia6) = xwerte(i,ia1)
            ywerte(i,ia3) = ywerte(i,ia3)/nra
            ywerte(i,ia6) = ywerte(i,ia6)/nra
            yerror(i,ia3) = yerror(i,ia3)/nra
            yerror(i,ia6) = yerror(i,ia6)/nra
            yerror(i,ia3) = sqrt(yerror(i,ia3)-ywerte(i,ia3)**2)
            yerror(i,ia6) = sqrt(yerror(i,ia6)-ywerte(i,ia6)**2)
         enddo

         name(ia2) = 'fftabs '
         name(ia3) = 'fftave '
         name(ia4) = 'filter '
         name(ia5) = 'filtered'
         name(ia6) = 'filtaver'
         xname(ia2)= '_1/dj*gm'
         xname(ia3)= '_1/dj*gm'
         xname(ia4)= '_1/dj*gm'
         xname(ia5)= '   dj/gm'
         xname(ia6)= '   dj/gm'
         numor(ia2)= numor(ia1)+1000
         numor(ia3)= numor(ia1)+2000
         numor(ia4)= numor(ia1)+3000
         numor(ia5)= numor(ia1)+4000
         numor(ia6)= numor(ia1)+5000
         nwert(ia4)= n1
         nwert(ia5)= n1
         nwert(ia6)= n1
         write(6,*)'peak ---> ',sump
         write(6,*)'noise --> ',sumb
         call parset('peak    ',sump,ia2)
         call parset('mpk     ',float(mpk),ia2)
         call parset('noisef  ',sumb,ia2)

         goto 2000
       endif

       if(comand.eq.'uni_ft   ') then
!                    -------
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= uni_ft [store_at <n>] [resnorm]  [rexpand <val>]  [reslim <val>]            '
           write(6,*)'=      performs fourier transform of the first selected records               '
           write(6,*)'=      the second selected record (may) hold the resolution                   '
           write(6,*)'=   options:                                                                  '
           write(6,*)'=            resnorm : normalisation to resolution intensity                  '
           write(6,*)'=            rexpand : multiplies (instead of adding)                         '
           write(6,*)'=            reslim:   limit for resolution at deconvolution                  '
           write(6,*)'=   spectrum x-axis unit given by valid names: micro-eV, meV, GHz, omega      '
           write(6,*)'=   intermediate results (Ft(data), Ft(res), ..) are stored on etra records   '
           write(6,*)'=============================================================================='
           goto 2000
        endif

         if(nsel.le.0) then
           write(6,*)'error: no data records selected: no action !'
           write(6,*)
           write(6,*)'select the data record and the resolution record '
           write(6,*)'using the sel command; sequence of selcted items matters! '
           write(6,*)'if only one record is selected, only simple FT without '
           write(6,*)'deconvolution is performed'
           ierrs = 1
           goto 2000
         endif 
  
         if(nsel.gt.2) then
           write(6,*)'error: too many data records selected: no action !'
           write(6,*)
           write(6,*)'select only the data record and the resolution record '
           write(6,*)' ... for the rest its not clear what to do... '
           write(6,*)' ... in order to process several q-values reiterate '
           write(6,*)' ... directly or by using a makro '
           ierrs = 2
           goto 2000
         endif   

         ia1 = nbuf+1
         ia1 = intval('store_at',ia1,inew)

         if(found('resnorm  ') ) then
            write(6,*)"Option:   resnorm    ==> use resolution for intensity normalisation also"
         else
            write(6,*)"NO option: resnorm  specified ==> use resolution for shape deconvolution ONLY"
         endif

         unift_range_expand     = getval('rexpand ',unift_range_expand,ier)
         unift_resolution_limit = getval('reslim  ',unift_resolution_limit,ier)

         
         write(6,*)'uni_ft: selected records and store data beginning from recordnr.:',ia1
         if (nsel.eq.1) then 
           isels(2) = 0
         endif
         call uni_ft(isels(1), isels(2), ia1, found('resnorm  '), unift_range_expand, unift_resolution_limit )

         goto 2000
       endif

!
       if(comand.eq.'arit    ') then ! alte version von arit 
!                    ----
         newnum = 777777
         idimux = 0

         num1 = 0
         num2 = 0
         do 4711 i=1,inames
          j = inapa(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'div     ') idimux = 1
          if(vname(i).eq.'mult    ') idimux = 2
          if(vname(i).eq.'to      ') newnum = Nint(rpar(j))
          if(vname(i).eq.'sc      ') then
            num1 = Nint(rpar(j))
            num2 = Nint(rpar(j+1))
          endif
          if(vname(i).eq.'factor1 '.or.vname(i).eq.'f1      ')          &
     &                                                 facto1 = rpar(j)
          if(vname(i).eq.'factor2 '.or.vname(i).eq.'f2      ')          &
     &                                                 facto2 = rpar(j)
 4711    continue
! --- figure out the adresses ---
         if(nsel.eq.2) then
           iad1 = isels(1)
           iad2 = isels(2)
         else
           iad1 = 0
           iad2 = 0
         endif
da1:     do i=1,nbuf
          if(numor(i) == 0 ) cycle da1
          if(numor(i).eq.num1) iad1 = i
          if(numor(i).eq.num2) iad2 = i
         enddo  da1
         if(iad1.eq.0) then
           write(6,*)'file :',num1,' not found'
           goto 2000
         endif
         if(iad2.eq.0) then
           write(6,*)'file :',num2,' not found'
           goto 2000
         endif
! --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4715 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
 4715    continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.mbuf) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
! ---- normalize to monitor-values ----
         if(withmo) then
           write(6,*)'normalizing to monitor-values ....'
           write(6,*)'use option nonorm to switch off, norm to swit. on'
! ---- extract monitorvalues ----
            xmon1 = 0
           do 4717 i=1,nopar(iad1)
             if(napar(i,iad1).eq.'monitor  ') xmon1 = params(i,iad1)
 4717      continue
           xmon2 = 0
           do 4719 i=1,nopar(iad2)
             if(napar(i,iad2).eq.'monitor  ') xmon2 = params(i,iad2)
 4719      continue
           if(xmon1.eq.0) then
             write(6,*)'monitor of first file is lacking take 1'
             xmon1 = 1
           endif
           if(xmon2.eq.0) then
             write(6,*)'monitor of 2nd   file is lacking take 1'
             xmon2 = 1
           endif
         else
           xmon1 =1
           xmon2 =1
         endif
!        ------> of if(withmo ..
! --- start the computation ---
         n1 = nwert(iad1)
         n2 = nwert(iad2)
         n3 = 0
         do 4721 i=1,n1
          x1 = xwerte(i,iad1)
!         -----> look for appropriate values in the seconf file
          j1 = 0
          j2 = 0
          d1 = 1.e30
          d2 = 1.e30
          do 4723 j=1,n2
           x2 = xwerte(j,iad2)
           if(abs(x1-x2).lt.abs(d1).and.x1.ge.x2) then
            j1 = j
            d1 = x1-x2
           endif
           if(abs(x1-x2).lt.abs(d2).and.x1.lt.x2) then
            j2 = j
            d2 = x1-x2
           endif
 4723     continue
          if(j1.eq.0.or.j2.eq.0) goto 4721
           n3 = n3+1
           xwerte(n3,newadd) = x1
           xx = x1
           y0 = ywerte(i ,iad1)
           y00= y0
           y1 = ywerte(j1,iad2)
           y2 = ywerte(j2,iad2)
           x1 = xwerte(j1,iad2)
           x2 = xwerte(j2,iad2)
           y11= y2 + (y1-y2)/(x1-x2)*(xx-x2)
           if(idimux.eq.0 ) y  = facto1 * y0/xmon1 +                    &
     &          facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2
           if(idimux.eq.1 ) y  = (facto1 * y0/xmon1) /                  &
     &        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           if(idimux.eq.2 ) y  = (facto1 * y0/xmon1) *                  &
     &        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           ywerte(n3,newadd) = y
           yyy=y
! --- and the errors ! ----
           y0 = yerror(i ,iad1)
           y1 = yerror(j1,iad2)
           y2 = yerror(j2,iad2)
           x1 = xwerte(j1,iad2)
           x2 = xwerte(j2,iad2)
           if(idimux.eq.0) y = (facto1 * y0/xmon1)**2 +                 &
     &          (facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)**2
           if(idimux.eq.1) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           if(idimux.eq.2) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           yerror(n3,newadd) = sqrt(y)
 4721    continue
         call txfpar(iad1,newadd)
         nwert(newadd) = n3
         numor(newadd) = newnum
         xname(newadd) = xname(iad1)
         yname(newadd) = yname(iad1)
          name(newadd) =  name(iad1)
         coment(newadd)= coment(iad1)(1:40)//coment(iad2)(1:40)
         nnpar=nopar(iad1)
       if(nnpar.le.mpar-2) then
         nnpar=nnpar+2
         nopar(newadd) = nnpar
         params(1,newadd) = facto1
         params(2,newadd) = facto2
         napar(1,newadd)  = name(iad1)
         napar(2,newadd)  = name(iad2)
         np=3
       endif
         do 4728 i=np,nnpar
           params(i,newadd)= params(i+1-np,iad1)
           napar(i,newadd) = napar (i+1-np,iad1)
 4728    continue
         write(6,*)n3,' points resulting from arit are stored as sc',   &
     &                  newnum
         if(idimux.eq.0)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' + ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.1)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' / ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.2)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' * ',facto2,' * ',   &
     &             name(iad2)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
!
!
       if(comand.eq.'arit2    ') then  ! deals better with interpolation at boundaries of range (vorher arit2)!
!                    ----
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= arit f1 <factor1> f2 <factor2> [to <numor>]  [options]                      '
           write(6,*)'=      adds the two selected records with factors 1,2 and stores as numor     '
           write(6,*)'=      for nonmatchig x-values interpolation is applied                       '
           write(6,*)'=   options:                                                                  '
           write(6,*)'=            div     : divides (instead of adding)                            '
           write(6,*)'=            mult    : multiplies (instead of adding)                         '
           write(6,*)'=            sc <numor1> <numor2>  : selection by numor match instead of sel  '
           write(6,*)'=   CHECK FUNCTIONALITY!!!!!!                                                 '
           write(6,*)'=============================================================================='
        endif

         newnum = 777777
         idimux = 0

         num1 = 0
         num2 = 0
         do i=1,inames
          j = inapa(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'div     ') idimux = 1
          if(vname(i).eq.'mult    ') idimux = 2
          if(vname(i).eq.'to      ') newnum = Nint(rpar(j))
          if(vname(i).eq.'sc      ') then
            num1 = Nint(rpar(j))
            num2 = Nint(rpar(j+1))
          endif
          if(vname(i).eq.'factor1 '.or.vname(i).eq.'f1      ')          &
     &                                                 facto1 = rpar(j)
          if(vname(i).eq.'factor2 '.or.vname(i).eq.'f2      ')          &
     &                                                 facto2 = rpar(j)
        enddo
! --- figuer out the adresses ---
         if(nsel.eq.2) then
           iad1 = isels(1)
           iad2 = isels(2)
         else
           iad1 = 0
           iad2 = 0
         endif
da12:    do i=1,nbuf
          if(numor(i) == 0 ) cycle da12
          if(numor(i).eq.num1) iad1 = i
          if(numor(i).eq.num2) iad2 = i
         enddo da12
         if(iad1.eq.0) then
           write(6,*)'file :',num1,' not found'
           goto 2000
         endif
         if(iad2.eq.0) then
           write(6,*)'file :',num2,' not found'
           goto 2000
         endif
! --- if result is already there , replace it ---
         newadd = nbuf+1
         do i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
         enddo
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.mbuf) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
! ---- normalize to monitor-values ----
         if(withmo) then
           write(6,*)'normalizing to monitor-values ....'
           write(6,*)'use option nonorm to switch off, norm to swit. on'
! ---- extract monitorvalues ----
            xmon1 = 0
           do i=1,nopar(iad1)
             if(napar(i,iad1).eq.'monitor  ') xmon1 = params(i,iad1)
           enddo
           xmon2 = 0
           do i=1,nopar(iad2)
             if(napar(i,iad2).eq.'monitor  ') xmon2 = params(i,iad2)
           enddo
           if(xmon1.eq.0) then
             write(6,*)'monitor of first file is lacking take 1'
             xmon1 = 1
           endif
           if(xmon2.eq.0) then
             write(6,*)'monitor of 2nd   file is lacking take 1'
             xmon2 = 1
           endif
         else
           xmon1 =1
           xmon2 =1
         endif
!        ------> of if(withmo ..
! --- start the computation ---
         n1 = nwert(iad1)
         n2 = nwert(iad2)
         n3 = n1
         if(n3.gt.mwork) then
           write(6,*)'arit2 workspace insufficient'
           goto 2000
         endif
         call get_pair_of_points(xwerte(1,iad1),ywerte(1,iad1),yerror(1,iad1),n1, &
    &                            xwerte(1,iad2),ywerte(1,iad2),yerror(1,iad2),n3, &
    &                            xwerte(1,newadd), yout1, yerout1, yout2, yerout2, nout)

         nwert(newadd) = nout
         do i=1,nout
          x1 = xwerte(i,newadd)
          y0 = yout1(i)
          y1 = yout2(i)
          y00 = y0
          y11 = y1
           if(idimux.eq.0 ) y  = facto1 * y0/xmon1 +                    &
     &          facto2 * y1 / xmon2
           if(idimux.eq.1 ) y  = (facto1 * y0/xmon1) /                  &
     &        (  facto2 *y1 / xmon2)
           if(idimux.eq.2 ) y  = (facto1 * y0/xmon1) *                  &
     &        (  facto2 * y1 / xmon2)
           ywerte(i,newadd) = y
           yyy=y
! --- and the errors ! ----
           y0 = yerout1(i)
           y1 = yerout2(i)
           if(idimux.eq.0) y = (facto1 * y0/xmon1)**2 +                 &
     &          (facto2 * y1 / xmon2)**2
           if(idimux.eq.1) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * y1 / xmon2)/y11)
           if(idimux.eq.2) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * y1 / xmon2)/y11)
           yerror(i,newadd) = sqrt(y)
         enddo
         call txfpar(iad1,newadd)
         nwert(newadd) = n3
         numor(newadd) = newnum
         xname(newadd) = xname(iad1)
         yname(newadd) = yname(iad1)
          name(newadd) =  name(iad1)
         coment(newadd)= coment(iad1)(1:40)//coment(iad2)(1:40)
         nnpar=nopar(iad1)
       if(nnpar.le.mpar-2) then
         nnpar=nnpar+2
         nopar(newadd) = nnpar
         params(1,newadd) = facto1
         params(2,newadd) = facto2
         napar(1,newadd)  = name(iad1)
         napar(2,newadd)  = name(iad2)
         np=3
       endif
         do  i=np,nnpar
           params(i,newadd)= params(i+1-np,iad1)
           napar(i,newadd) = napar (i+1-np,iad1)
         enddo
         write(6,*)n3,' points resulting from arit are stored as sc',   &
     &                  newnum
         if(idimux.eq.0)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' + ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.1)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' / ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.2)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' * ',facto2,' * ',   &
     &             name(iad2)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
!
!
       if(comand.eq.'average ') then
!                    -------
         call average_data()

         goto 2000
       endif
!
       if(comand.eq.'addval  ') then
!                    ------
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= addval <x> <y> <err>                                                        '
           write(6,*)'=    appends one data point to selected record                                '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
        endif

         if(nsel.eq.1) then
           iad1 = isels(1)
         else
           write(6,*)'select one and only one record'
           ierrs= 1
           goto 2000
         endif
! --- start the computation ---
         
         nwert(iad1) = nwert(iad1)+1
         xwerte( nwert(iad1), iad1 ) = rpar(1)
         ywerte( nwert(iad1), iad1 ) = rpar(2)
         yerror( nwert(iad1), iad1 ) = rpar(3)

         goto 2000
       endif
!
       if(comand.eq.'integrat') then
!                    --------           approximates integral from binned data
!                                       by simple summation of data times bin-width
!                                       there may be better interpolation based schemes
!                                       but this is clear an simple                        
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= integrat x1 <x1> x2 <x2>                                                    '
           write(6,*)'=        approximates integral from binned data over the given range          '
           write(6,*)'=        by simple summation of data times bin-width                          '
           write(6,*)'=        there may be better interpolation based schemes                      '
           write(6,*)'=        but this is clear an simple                                          '
            write(6,*)'=============================================================================='
           goto 2000
        endif
         if(nsel<1) then
           call errsig(200,"ERROR: integrat slect at least one record$")
           goto 2000
         endif
! --- start the computation --- these are the integration limits
 
         x1int        = getval('x1      ',dble(x1int),inew)
         x2int        = getval('x2      ',dble(x2int),inew)
 
         do j=1,nsel
            iad1 = isels(j)         
            i1int = 1
            i2int = nwert(iad1)
            
            do i=1,nwert(iad1)-1
             if(xwerte(i,iad1).le.x1int .and.xwerte(i+1,iad1).gt.x1int) i1int=i
             if(xwerte(i,iad1).le.x2int .and.xwerte(i+1,iad1).gt.x2int) i2int=i
            enddo
       
            
            sum   = (xwerte(i1int,iad1)-x1int)*ywerte(i1int,iad1)
            sumer = ((xwerte(i1int,iad1)-x1int)*yerror(i1int,iad1))**2       
       
            sum   = sum + (x2int-xwerte(i2int,iad1))*ywerte(i2int,iad1)
            sumer = sumer + ((x2int-xwerte(i2int,iad1))*yerror(i2int,iad1))**2
            do i=i1int+1,i2int-1
              sum   = sum   + (xwerte(i+1,iad1)-xwerte(i,iad1))*ywerte(i+1,iad1)
              sumer = sumer + ((xwerte(i+1,iad1)-xwerte(i,iad1))*yerror(i+1,iad1))**2
            enddo
           
            sumer = sqrt(sumer)
       
            write(6,*)"================= Record: ",iad1,"  ===================="
            write(6,'(a,i9,1x,e13.6,a,e13.6,a)')'integral from: ',i1int,x1int,' (',xwerte(i1int,iad1),' )'
            write(6,'(a,i9,1x,e13.6,a,e13.6,a)')'           to: ',i2int,x2int,' (',xwerte(i2int,iad1),' )'
          
            write(6,*                  )'         ===>: ',sum, '+-', sumer
!            xwerte( nwert(iad1), iad1 ) = rpar(1)
!            ywerte( nwert(iad1), iad1 ) = rpar(2)
!            yerror( nwert(iad1), iad1 ) = rpar(2)
       
            call parset ('integral',sum        ,iad1 ) 
            call parset ('x1integ ',x1int      ,iad1 ) 
            call parset ('x2integ ',x2int      ,iad1 ) 
         enddo 
         goto 2000
       endif
!

!
!
       if(comand.eq.'combine ') then
!                    -------
         newnum = 777777
         do 4811 i=1,inames
          j = inapa(i)
          k = inpar(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'raster  ') then
           x1raster = rpar(j)
           x2raster = rpar(j+1)
           nraster  = Nint(rpar(j+2))
          endif
          if(vname(i).eq.'to      ') newnum = Nint(rpar(j))
          if(vname(i).eq.'sc      ') then
             do 48181 l=1,k
               nnumi(l) = Nint(rpar(l-1+j))
48181        continue
             call search(nnumi,k)
          endif
 4811    continue
! --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4815 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
 4815    continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.mbuf) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
! --- create the raster ----------------------------------------------
         if(nraster.gt.mwert) nraster = mwert
         if(nraster.lt.2    ) nraster = 2
         draster = (x2raster-x1raster)/(nraster-1)
         do  i=1,nraster
           xwerte(i,newadd) = x1raster + (i-1)*draster
           ywerte(i,newadd) = 0.d0
           yerror(i,newadd) = 0.d0
         enddo
         write(6,*)'raster generated from x1=',x1raster
         write(6,*)'                 to   x2=',x2raster
         write(6,*)'with a number of points of n = ',nraster
         write(6,*)'the newaddress is at ......... ',newadd
! --- loop over the selectionlist ----------------------------
         nwert(newadd) = nraster
         do isl=1,nsel
           iadd = isels(isl)
           xmonitor = 1.d0
           call parget ('monitor ',xmonitor,iadd,ier)
           n1 = nwert(newadd)
           n2 = nwert(iadd)
           do 4821 i=1,n1
            x1 = xwerte(i,newadd)
!         -----> look for appropriate values in the second file
            j1 = 0
            j2 = 0
            d1 = 1.e30
            d2 = 1.e30
            do 4823 j=1,n2
             x2 = xwerte(j,iadd)
             if(abs(x1-x2).lt.abs(d1).and.x1.ge.x2) then
              j1 = j
              d1 = x1-x2
             endif
             if(abs(x1-x2).lt.abs(d2).and.x1.lt.x2) then
              j2 = j
              d2 = x1-x2
             endif
 4823       continue
            if(j1.eq.0.or.j2.eq.0) goto 4821
            xx = x1
            y0 = ywerte(i ,newadd)
            y1 = ywerte(j1,iadd)
            y2 = ywerte(j2,iadd)
            x1 = xwerte(j1,iadd)
            x2 = xwerte(j2,iadd)
            y  = y0 + (y2 + (y1-y2)/(x1-x2)*(xx-x2))
            ywerte(i,newadd) = y
            yerror(i,newadd) = yerror(i,newadd) + xmonitor
!           ----- temorary misuse of the error-filed to collect
!                 monitor-values
 4821      continue
         enddo
! ---- normalize the data and create errors --------
         do i=1,n1
           if(yerror(i,newadd).gt.0.0) then
             yer     = sqrt(ywerte(i,newadd))/ yerror(i,newadd)
             ywerte(i,newadd) = ywerte(i,newadd) / yerror(i,newadd)
             yerror(i,newadd) = yer
           endif
         enddo
! --- transfer parameters ---
         call txfpar(isels(1),newadd)
         nwert(newadd) = n1
         numor(newadd) = newnum
         xname(newadd) = xname(isels(1))
         yname(newadd) = yname(isels(1))
         name(newadd)  =  name(isels(1))
         coment(newadd)= coment(isels(1))(1:40)//coment(isels(2))(1:40)
         nnpar=nopar(isels(1))
         call  parset ('monitor ',1.0,newadd)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
!
!
!
!
       if(comand.eq.'m       '.or.comand.eq.'mirror  ') then
!                    -                       ------
         xk0= rpar(1)
         k1 = rpar(2)+0.0001
         k2 = rpar(3)+0.0001
         if(ipars.le.0) xk0 = 32
         if(ipars.le.1) k1  = xk0 - 20
         if(ipars.le.1) k2  = xk0 - 10
         if(nsel.eq.0) then
           write(6,*)'no scan selected ; plot one to select it'
           goto 2000
         endif
         ia = isels(1)
         if(numor(ia).ge.numpls) then
           write(6,*) numor(ia),'=numor > ',numpls,'; warning no origina&
     &                data'
         endif
         if(xname(ia).ne.'kanal   '.and.xname(ia).ne.'kanal   ') then
           write(6,*)'warning x-values are not kanal !!!'
         endif
         nw = nwert(ia)
! --- determine the center ---
         xk00 = xk0
         call detxk0(ywerte(1,ia),ywerte(1,nbuf+1),xk0,k1,k2,nw)
         if(xk0.le.0)then
            xk0 = xk00
            call detxk1(ywerte(1,ia),ywerte(1,nbuf+1),xk0,k1,k2,nw)
         endif
         write(6,*)'center channel is = ',xk0
         write(6,*)'generate the mirror-image ...'
         f1 = smirro(ywerte(1,ia),ywerte(1,nbuf+1),xk0,1,nw,nw)
         do 1702 i=1,nw
           xwerte(i,nbuf+1) = xwerte(i,ia)
 1702    continue
         nbuf       = nbuf + 1
         name(nbuf) = name(ia)
         xname(nbuf)= xname(ia)
         yname(nbuf)= 'sym-'//yname(ia)(1:4)
         numor(nbuf)=mod(numor(ia),numpls)+numpls
         coment(nbuf)=coment(ia)
         nwert(nbuf)= nw
         write(6,*)'ok'
         isels(2) = ia
         isels(1) = nbuf
         ifits(1)=0
         ifits(2)=0
         nsel     = 2
         npp=nopar(ia)
         if(npp.ge.mpar) then
           write(6,*)'xk0 cannot be stored on params because too many p'
           goto 2000
         endif
!
         call parset ('xk0     ',xk0,ia)
         nopar(nbuf)= npp+1
         do 2007 i=1,npp+1
          params(i,nbuf) = params(i,ia)
          napar (i,nbuf) = napar (i,ia)
 2007    continue
         goto 2000
       endif
!
!
!
       if(comand.eq.'ercorrc '.or.comand.eq.'ecc     ') then
!                    -------                 ---
          if(nbuf.eq.mbuf) then
            write(6,*)' no space for an additional item'
            goto 2000
          endif
         nbuf = nbuf+1
         nsel = 1
         isels(1)= nbuf
         nfits        = 0
         numor(nbuf)  = 100*nbuf
         xname(nbuf)  = 'x'
         yname(nbuf)  = 'err'
         name(nbuf)   = 'cc_err'
         coment(nbuf)= 'cc_err'
         n            = intval('n       ',n,inew)
         aa           = getval('a       ',dble(aa),inew)
         ri           = getval('ri      ',dble(ri),inew)
         r0           = getval('r0      ',dble(r0),inew)
         d0           = getval('d0      ',dble(d0),inew)
         ncut         = intval('ncut    ',ncut,inew)

         if(n.gt.mwert) then
           write(6,*) ' n is too large'
           goto 2000
         endif

         nwert(nbuf)  = n

         if(ncut.gt.mcut) then
           write(6,*)' ncut is too large max=',mcut
           goto 2000
         endif

         write(6,*)' a ........ = ',a
         write(6,*)' ri ....... = ',ri
         write(6,*)' r0 ....... = ',r0
         write(6,*)' d0 ....... = ',d0
         write(6,*)' ncut ..... = ',ncut
         write(6,*)' n ........ = ',n


         call parset ('a       ',a  ,nbuf)
         call parset ('ri      ',ri ,nbuf)
         call parset ('r0      ',r0 ,nbuf)
         call parset ('d0      ',d0 ,nbuf)
         xnc = ncut
         call parset ('ncut    ',xnc ,nbuf)

         a2 = a/2.0
         x1 = 0.0
         xm = r0*sqrt(xnc)
         dx = xm  / n
         do i=1,ncut+1
           xi = i
           rrv(i) = r0 *sqrt( xi )
         enddo
         rrv(0) = ri
         do i=0,n
           xr = x1 + i*dx
           j  = 0
           y  = 0.0
           do j=0,ncut
            if( xr .lt. rrv(j)) goto 20022
           enddo
20022      continue
           if(j.gt.0) then
             xn = a2*(rrv(j)**2-rrv(j-1)**2)+d0*alog(rrv(j)/rrv(j-1))
             xz = a2*(xr**2    -rrv(j-1)**2)+d0*alog(xr/rrv(j-1))
             y  =  xz/xn + j-1
           endif
20032      continue

           ywerte(i,nbuf) = y
           yerror(i,nbuf) = 0
           xwerte(i,nbuf) = xr

         enddo

         goto 2000
       endif
!
!
       if(comand.eq.'sym     ') then
!                    ---->  symmetrize
! ---- build now the symmetric average on the right side ---
         if(nsel.lt.2) then
           write(6,*)'invalid isel-list, start with mirro'
           goto 2000
         endif
         ib = isels(1)
         ia = isels(2)
         nn1 = numor(ia)
         nn2 = numor(ib)
         if(nn1+numpls.ne.nn2) then
            write(6,*)'no pair present, use mirro'
            goto 2000
         endif
! --- look for xk0 ---
         do 2077 i=1,nopar(ib)
            if(napar(i,ib).eq.'xk0     ') then
              xk0 = params(i,ib)
              npp = i
            endif
 2077    continue
            if(vname(1).eq.'xk0     ') xk0=rpar(1)
         write(6,*)'xk0 = ',xk0
         nw = nwert(ia)
         kan1 = xk0 + 1
         do 2010 i=kan1,nw
           yw1 = ywerte(i,ia)
           yw2 = ywerte(i,ib)
           yy  = ( yw1 + yw2 ) / 2
           if(yw1.eq.0) yy = yw2
           if(yw2.eq.0) yy = yw1
           ywerte(i+1-kan1,ib) = yy
 2010    continue
         nwert(ib) = nw-kan1+1
         xk0            = xk0 - kan1 + 1
         params(npp,ib) = xk0
         numor(ib)=mod(numor(ib),numpls)+2*numpls
         write(6,*)'ok symm-ds has ',nwert(ib),' points'
         goto 2000
       endif
!
!
       if(comand.eq.'kz      ') then
!                    -->    kanalzusammenfassung
! ---- build now the symmetric average on the right side ---
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= kz <nk>                                                                     '
           write(6,*)'=    contraction of spectra by adding <nk> adjacent channele                  '
           write(6,*)'=    channel width paramezters etc. are update in the parameterlist           '
           write(6,*)'=    such that spectra fitting stays with the proper _xwidth                  '
           write(6,*)'=    The contratction is on place (i.e.)                                      '
           write(6,*)'=    All selected items are treated                                           '
           write(6,*)'=============================================================================='
        endif

         kz = NINT(rpar(1))
         if(kz.lt.2) goto 2000
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          ikz=kz
          if(ikz.gt.n)ikz=n
          j  = 0
          m  = n/ikz
          ir = n-ikz*m
          do l=1,m
            xsum=0.0
            ysum=0.0
            esum=0.0
            do k=1,ikz
              ik=(l-1)*ikz+k
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ikz
            ywerte(j,ia) = ysum/ikz
            yerror(j,ia) = sqrt(esum)/ikz
          enddo
         if(ir.gt.0) then
            xsum=0.0
            ysum=0.0
            esum=0.0
            do ik=m*ikz+1,n
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ir
            ywerte(j,ia) = ysum/ir
            yerror(j,ia) = sqrt(esum)/ir
         endif
         nwert(ia)=j
         call parget ('_rebin  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * ikz
         call parset ('_rebin  ',xx,ia)

         call parget ('_xwidth ',xx,ia,ier)
         if(ier.ne.0) then
           xx = xwerte(m/2+1,ia)-xwerte(m/2+1,ia)
         else
           xx = xx * ikz
         endif
         call parset ('_xwidth ',xx,ia)

         write(6,*)'channels contracted for ',ia,' new no. of ch=',j,   &
     &             ' (',ir,')'
         enddo
         goto 2000
       endif



       if(comand.eq.'rebin   ') then
!                    ----->    kanalzusammenfassung mit Kopie auf
!                              neuen Platz
! 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rebin <nk>                                                                  '
           write(6,*)'=    contraction of spectra by adding <nk> adjacent channele                  '
           write(6,*)'=    channel width paramezters etc. are update in the parameterlist           '
           write(6,*)'=    such that spectra fitting stays with the proper _xwidth                  '
           write(6,*)'=    The contracted spectra are stored on new records.                        '
           write(6,*)'=    All selected items are treated   (similar kz)                            '
           write(6,*)'=============================================================================='
        endif
         kz = NINT(rpar(1))
         if(kz.lt.2) goto 2000
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

!     >  copy first  
         if(nsel+nbuf.gt.mbuf) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          ikz=kz
          if(ikz.gt.n)ikz=n
          j  = 0
          m  = n/ikz
          ir = n-ikz*m
          do l=1,m
            xsum=0.0
            ysum=0.0
            esum=0.0
            do k=1,ikz
              ik=(l-1)*ikz+k
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ikz
            ywerte(j,ia) = ysum/ikz
            yerror(j,ia) = sqrt(esum)/ikz
          enddo
         if(ir.gt.0) then
            xsum=0.0
            ysum=0.0
            esum=0.0
            do ik=m*ikz+1,n
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ir
            ywerte(j,ia) = ysum/ir
            yerror(j,ia) = sqrt(esum)/ir
         endif
         nwert(ia)=j
         numor(ia)=numor(ia)+100000
         call parget ('_rebin  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * ikz
         call parset ('_rebin  ',xx,ia)

         call parget ('_xwidth ',xx,ia,ier)
         if(ier.ne.0) then
           xx = xwerte(m/2+1,ia)-xwerte(m/2+1,ia)
         else
           xx = xx * ikz
         endif
         call parset ('_xwidth ',xx,ia)

         

         write(6,'(a,i4,a,i4)')'channels of[',ia,'] contracted by ',ikz 
         enddo
         goto 2000
       endif
!


       if(comand.eq.'rerange ') then
!                    ----->    neuer x-range durch selektion
!                              neuen Platz
!           
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rerange <min> <max>     [y]                                                 '
           write(6,*)'=    copy of selected records such that                                       '
           write(6,*)'=    only values betwen min < x|y < max  are retained , y: yrange             '
           write(6,*)'=    without option x-range is pertained                                      '
           write(6,*)'=============================================================================='
        endif


         lower_range = rpar(1)
         upper_range = rpar(2)
         range_is_y  = found('y       ')

!     >  copy first  
         if(nsel+nbuf.gt.mbuf) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          j  = 0
drer1:    do ik=1,n
             if(range_is_y) then
               if(ywerte(ik,ia) < lower_range .or. ywerte(ik,ia) > upper_range) cycle drer1
             else
               if(xwerte(ik,ia) < lower_range .or. xwerte(ik,ia) > upper_range) cycle drer1
             endif
             j = j+1
             xwerte(j,ia) = xwerte(ik,ia)
             ywerte(j,ia) = ywerte(ik,ia)
             yerror(j,ia) = yerror(ik,ia)
          enddo drer1

          nwert(ia)=j
          numor(ia)=numor(ia)+10000
          if(range_is_y) write(6,'(a)',advance='no') "Y-"
          write(6,'(a,i4,a,2e14.6)')'range of[',ia,'] limited to ',lower_range, upper_range 
         enddo
         goto 2000
       endif
!




       if(comand.eq.'swapxy   ') then
!                    ----->    vertauschen von x und y
!                              
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= swapxy                                                                      '
           write(6,*)'=    role of x amd y-values are exchanged                                     '
           write(6,*)'=    NOTE: y-errors are lost                                                  '
           write(6,*)'=    operation is on place no new copies!                                     '
           write(6,*)'=============================================================================='
        endif
! 
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do l=1,n
              xsum = xwerte(l,ia)
              ysum = ywerte(l,ia)
              xwerte(l,ia) = ysum
              ywerte(l,ia) = xsum             
              yerror(l,ia) = 0
           enddo
           cbuffer   = xname(ia)
           xname(ia) = yname(ia)
           yname(ia) = cbuffer
           write(6,'(a,i5,a)')" x and y values of ",ia," have been swapped, errors set to zero!" 
         enddo
         goto 2000
       endif
!

       if(comand.eq.'sequence   ') then
!                    ----->  sequence  
!                                       
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= sequence                                                                    '
           write(6,*)'=    replace x-values by sequence numbers 1..n                                '
           write(6,*)'=    operation is on place, previous x-values are lost                        '
           write(6,*)'=============================================================================='
        endif

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do l=1,n
             xwerte(l,ia) = l
          enddo
           xname(ia) = "#"
           write(6,'(a,i5,a)')" x values of ",ia," have been replace by sequence numbers" 
         enddo
         goto 2000
       endif
!
!
       if(comand.eq.'copy   ') then
!                    ----->  copy 
!                                       
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= copy [x1 <min> x2 <max>                                                     '
           write(6,*)'=    copy of selected records to new record buffers                           '
           write(6,*)'=    if x1, x2 are given only x-values within the interval are copied         '
           write(6,*)'=============================================================================='
        endif

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         xcopy1 = getval('x1      ',dble(xcopy1),ier)
         xcopy2 = getval('x2      ',dble(xcopy2),ier)

         do i=1,nsel
          ia = isels(i)
          if(nbuf < mbuf) then
             nbuf = nbuf + 1
          else
             Write(6,*)"ERROR: too many records !"
             goto 2000
          endif

           call DataCopy(ia,nbuf)
           write(6,'(a,i5,a,i5)')" copy ", ia," to ",nbuf 
           isels(i) = nbuf 

           if(found('x1      ') .and. found('x2      ')) then
             Write(6,'(a,2e14.6)')" ...restricted to the interval: ",xcopy1, xcopy2
             nwert(nbuf) = 0
             do j=1,nwert(ia)
               if(xwerte(j,ia) >= xcopy1 .and. xwerte(j,ia) <= xcopy2) then
                  nwert(nbuf) = nwert(nbuf) + 1
                  xwerte(nwert(nbuf) ,ia)  = xwerte(j,ia)
                  ywerte(nwert(nbuf) ,ia)  = ywerte(j,ia)
                  yerror(nwert(nbuf) ,ia)  = yerror(j,ia)
               endif
             enddo
           endif

         enddo
         goto 2000
       endif
!
       

       if(comand.eq.'scale   ') then
!                    ----->    Skalierung der Intesnitaet          

        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= scale <factor>                                                              '
           write(6,*)'=    scales selected spectra with a factor, bookkeeping in parameter section  '
           write(6,*)'=    also for gaussian resolution parameters and errors                       '
           write(6,*)'=============================================================================='
        endif
                              
         if(ipars.ne.1) then
           write(6,*)'ERROR: need just one value as parameter..'
           ierrs = 1
           goto 2000
         endif

         y_scaling =  rpar(1)

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

!     >  copy first  
         if(nsel+nbuf.gt.mbuf) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do ik=1,n
              ywerte(ik,ia) = ywerte(ik,ia) *y_scaling
              yerror(ik,ia) = yerror(ik,ia) *y_scaling
          enddo
         numor(ia)=numor(ia)+200000

         call par_gaint_scale(y_scaling,ia)


         call parget ('_scale  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * y_scaling
         call parset ('_scale  ',xx,ia)

         write(6,'(a,i4,a,e13.6)')'values of[',ia,'] scaled by ',y_scaling 
         enddo
         goto 2000
       endif
!
!
!
!
       if(comand.eq.'gaiscale') then
!                    -------->    Skalierung der Aufloesungsparameter
!                              
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= gaiscale <factor>                                                           '
           write(6,*)'=    scaling of the gaussian resolution parameters associated                 '
           write(6,*)'=    with the selected spectra                                                '
           write(6,*)'=============================================================================='
        endif

         if(ipars.ne.1) then
           write(6,*)'ERROR: need just one value as parameter..'
           ierrs = 1
           goto 2000
         endif

         y_scaling =  rpar(1)

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif
 
         do i=1,nsel
           call par_gaint_scale(y_scaling,isels(i))
           write(6,'(a,i4,a,e13.6)')'ga#inten values of[',isels(i),'] scaled by ',y_scaling 
         enddo
         goto 2000
       endif
!




!
       if(comand.eq.'interpol') then
!                    ---->  interpolate
! ---- build now the symmetric average on the right side ---
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= interpol <n>                                                                '
           write(6,*)'=    creates a clone of the firts (only) item of the selection list           '
           write(6,*)'=    by interpolating y-values such that the x-range is covered by <n>        '
           write(6,*)'=    equidistant points.                                                      '
           write(6,*)'=    Errors are NOT treated (yet)                                             '
           write(6,*)'=============================================================================='
        endif

         ia = isels(1)
         ib = nbuf+1
         if(ib.gt.mbuf) ib=mbuf
         nbuf = ib
         call txfpar(ia,ib)
         nn = rpar(1)
         if(nn.le.1) nn=mwert
         dxx= (xwerte(nwert(ia),ia)-xwerte(1,ia))/(nn-1)
         do i=1,nn
           xx = xwerte(1,ia)+(i-1)*dxx
           xwerte(i,ib) = xx
           ywerte(i,ib) = yinterp(xx)
         enddo
         nwert(ib) = nn
         isels(1)  = ib
         nsel      = 1
         goto 2000
       endif
!
!
       if(comand.eq.'addsels ') then
!                    ---->  summation over selected records : will become obsolete due to name
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= addsels                                                                     '
           write(6,*)'=    creates record containing the sum of all selected (compatible)           '
           write(6,*)'=    records.      (similar recsum)                                           '
           write(6,*)'=============================================================================='
        endif
         call sumseldat
         goto 2000
       endif
!
!
       if(comand.eq.'recsum  ') then
!                    ---->  summation over selected records: the new one
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= recsum from <n1> to <n2>                                                    '
           write(6,*)'=    creates record containing the sum of all records from a continuous       '
           write(6,*)'=    address range <n1> .. <n2>      (similar addsels)                        '
           write(6,*)'=============================================================================='
        endif
         if(found('from    ')) then
           ia = intval('from    ',isels(1),inew)
           ib = intval('to      ',isels(nsel),inew)
           nsel = ib-ia+1
           do i=1,nsel
             isels(i) = ia-1+i
           enddo
           write(6,'(a,i4,a,i4)')'sum records from ',ia,' ... ',ib
         endif
         call sumseldat
         goto 2000
       endif
!
!
       if(comand.eq.'parextra') then
!                    ---> parabola extrapolation to q--> 0

        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= parextra                                                                    '
           write(6,*)'=    adds a point with x = 0 by parabola extrapolation of data starting with  '
           write(6,*)'=    finite x-values (e.g. extrapolation of SANS data to q=0)                 '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
        endif

         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
         ia    = isels(1)
         slope = (ywerte(3,ia)-ywerte(1,ia))/(xwerte(3,ia)-xwerte(1,ia))
         apar  = slope / (2.0*xwerte(2,ia))
         bpar  = ywerte(2,ia) - slope*0.5*xwerte(2,ia)
         dx    = xwerte(2,ia)-xwerte(1,ia)
         npl   = xwerte(1,ia) / dx + 1
         write(6,*)'parabola: a*x^2 + b'
         write(6,*)'a = ',apar
         write(6,*)'b = ',bpar
         write(6,*)'dx= ',dx
         write(6,*)'n+= ',npl
         nold = nwert(ia)
         if(nwert(ia)+npl.gt.mwert) then
            write(6,*)'too many points'
            ierrs = 1
            goto 2000
         endif
         nwert(ia) = nwert(ia) + npl
         yerr      = yerror(1,ia)
         do i=nwert(ia),npl+1,-1
            ywerte(i,ia) = ywerte(i-npl,ia)
            xwerte(i,ia) = xwerte(i-npl,ia)
            yerror(i,ia) = yerror(i-npl,ia)
         enddo
         do i=1,npl
           xwerte(i,ia) = (i-1)*dx
           ywerte(i,ia) = apar*xwerte(i,ia)**2 + bpar
           yerror(i,ia) = yerr
         enddo
         write(6,*)'done'
         goto 2000
        endif
!
!
       if(comand.eq.'hiqextra') then
!                    ---> extrapolation for high q .
!                         data == a*q^z
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= hiqextra <xcut>                                                             '
           write(6,*)'=    adds points from hiq x (q) extrapolation by a*q^z  up to x=xcut          '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
        endif
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
!        ---------------
         xcut  = rpar(1)
!        ---------------
         ia    = isels(1)
         do i=1,nwert(ia)
            if(xwerte(i,ia).gt.xcut) goto 77771
         enddo
77771    continue
         icut = i-1
         slope = (ywerte(icut,ia)-ywerte(icut-2,ia))/                   &
     &           (xwerte(icut,ia)-xwerte(icut-2,ia))
         zpar  = xwerte(icut-1,ia)*slope / ywerte(icut-1,ia)
         apar  = ywerte(icut-1,ia)*xwerte(icut-1,ia)**(-zpar)
         dx    = xwerte(icut,ia)-xwerte(icut-1,ia)
         do i=icut,mwert
           xx = xwerte(icut,ia)+(i-icut)*dx
           ywerte(i,ia) = apar*xx**zpar
           yerror(i,ia) = yerror(icut,ia)
           xwerte(i,ia) = xx
         enddo
         write(6,*)'high-q extrapolation :  a*q^z'
         write(6,*)'a =  ',apar
         write(6,*)'z =  ',zpar
         write(6,*)'dx=  ',dx
         write(6,*)'ncut=',icut
         nwert(ia) = mwert
         write(6,*)'done'
         goto 2000
        endif
!
       if(comand.eq.'qc      '.or.comand.eq.'q-conv  ') then
!                    ---> channels to q-values
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia = isels(1)
        modeqc=0
        if(xname(ia)(1:1).eq.'k'.or.xname(ia)(1:1).eq.'c')modeqc=1
        if(xname(ia).eq.'rad*1.e-6')modeqc=2
        if(modeqc.eq.0) then
           write(6,*)'data do not have the proper x-values:',xname(ia)
           ierrs=100
           goto 2000
         endif
         do 2001 i=1,nopar(ia)
           if(napar(i,ia).eq.'xk0      ') xk0   = params(i, ia)
           if(napar(i,ia).eq.'lambda   ') alam0 = params(i, ia)
           if(napar(i,ia).eq.'detdis   ') detdis= params(i, ia)
           if(napar(i,ia).eq.'bklen    ') bklen = params(i, ia)
 2001    continue
           if(nbuf.lt.mbuf)then
           nbuf=nbuf+1
           ib=nbuf
        endif
         do 2002 i=1,nwert(ia)
           if(modeqc.eq.1) then
             xwerte(i,ib)=4*pi/alam0*sin(0.5*atan(bklen*(i-xk0)/detdis))
           endif
           if(modeqc.eq.2) then
             xwerte(i,ib)=4*pi/alam0*sin(1.e-6*xwerte(i,ia)/2)
           endif
           ywerte(i,ib)=ywerte(i,ia)
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+3*numpls
 2002    continue
         xname(ib)='q      '
         write(6,*)'q1=',xwerte(1,ib),'  qend=',xwerte(nwert(ib),ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'sp      '.or.comand.eq.'spline  ') then
!                    ---> spline smoothing
          do 2799 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'iequal  ') iequal = rpar(j) + 0.1
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'auto    ') iauto  = 1
            if(vname(i).eq.'noauto  ') iauto  = 0
            if(vname(i).eq.'smpar   ') smpar  = rpar(j)
 2799     continue
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia     = isels(1)
        nn     = nwert(ia)
        numspl = numor(ia)
        nwspl  = nwert(ia)
        if(nneu.le.0.or.nneu.gt.mwert) nneu = nn
        if(iequal.le.0) iequal = 1
        write(6,*)'nneu  = ',nneu
        if(nbuf.lt.mbuf)then
           nbuf=nbuf+1
           ib=nbuf
        endif
!
        if(iauto.eq.1) then
           write(6,*)'iequal= ',iequal
           call csscv(nn,xwerte(1,ia),ywerte(1,ia),iequal,break,cscoef)
        else
           write(6,*)'smpar = ',smpar
           do 3777 i=1,nn
             weight(i) = sqrt(abs(ywerte(i,ia)))
 3777      continue
           call cssmh(nn,xwerte(1,ia),ywerte(1,ia),weight,smpar,break,  &
     &                cscoef)
        endif
         if(iout.gt.1) then
            write(6,*)'break         cscoeff'
            do 2776 i=1,nn
              write(6,2777)i,break(i),(cscoef(j,i),j=1,4)
 2777         format(1x,i3,':',5e14.6)
 2776       continue
         endif
         x01 = break(1)
         x02 = break(nn)
         dx  = (x02-x01)/(nneu-1)
         do 2702 i=1,nneu
            xn = x01 + (i-1)*dx
            xwerte(i,ib) = xn
            ywerte(i,ib) = fsplin(xn)
 2702    continue
         call txfpar(ia,ib)
         numor(ib)= - numor(ia)
         nwert(ib)= nneu
         isels(1) = ia
         ifits(1) = ib
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'fftmx   '.or.comand.eq.'fft-ms  ') then
!                    ---> multiple scattering by fft
           t1 = 0
           if(numspl.eq.0.or.nwspl.eq.0) then
             write(6,*)'spline must be made prior to fft'
             goto 2000
           endif
! -------- identify the original to get parameters ---
           do 3705 ia=1,nbuf
           if(numspl.eq.numor(ia)) goto 3707
 3705     continue
           write(6,*)'original data not found. unable to set parameters'
           goto 2000
 3707     continue
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          if(nx.gt.mwert) nx = mwert
          call txfpar(ia,ib)
          isels(1) = ia
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do 3799 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'mx      ') imx    = 0
            if(vname(i).eq.'rmx     ') imx    = 1
            if(vname(i).eq.'rmx1    ') imx    = 2
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') then
               xmax   = rpar(j)
               ymax   = 10.0/xmax
            endif
            if(vname(i).eq.'ymax    ') ymax   = rpar(j)
            if(vname(i).eq.'nx      ') nx     = rpar(j) +.01
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
 3799     continue
!
         if(imx.eq.0) then
           write(6,*)'evaluate multiple scattering ...'
           call fftmx(xwerte(1,ib),ywerte(1,ib),trans,xmax,dx,nx,nfft)
           numor(ib)= mod(numor(ia),numpls)+10*numpls
           yname(ib)= 'mx('//yname(ia)(1:5)
         else
           write(6,*)'evaluate inverse multiple scattering ...'
           if(imx.eq.1) then
             call fftrmx(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,     &
     &                                                  xmax,dx,nx,nfft)
           else
             call rmx1d(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,      &
     &                                                xmax,ymax,nx,nfft)
           endif
           numor(ib)= mod(numor(ia),numpls)+11*numpls
           yname(ib)= 'imx('//yname(ia)(1:4)
           call parset('i0      ',ai0,ib)
         endif
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nx
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'dmx     ') then
!                    ---> multiple scattering by fft
           t1 = 0
          ia = isels(1)
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          call txfpar(ia,ib)
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do i=1,inames
            j = inapa(i)
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') xmax   = rpar(j)
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
          enddo
!
         call demux(xwerte(1,ia),ywerte(1,ia),nwert(ia),                &
     &              trans,xmax,nfft,                                    &
     &              iout,xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'dmx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'mux     ') then
!                    ---> multiple scattering by fft
           t1 = 0
          ia = isels(1)
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          call txfpar(ia,ib)
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do i=1,inames
            j = inapa(i)
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') xmax   = rpar(j)
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
          enddo
!
         call   mux(xwerte(1,ia),ywerte(1,ia),nwert(ia),                &
     &              trans,xmax,nfft,                                    &
     &              iout,xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'mx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'tracorr ') then
          ia = isels(1)
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)
          do i=1,nwert(ia)
            ywerte(i,ia) = ywerte(i,ia)*thick*trans
            yerror(i,ia) = yerror(i,ia)*thick*trans
          enddo
          write(6,*)'selected item:',ia,' multiplied by:',              &
     &              thick,' * ',trans,' !'
         goto 2000
        endif
!
       if(comand.eq.'des     '.or.comand.eq.'desmear ') then
!                    ---> infinite slit desmearing
          do 2899 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'qmax    ') qmax   = rpar(j)
            if(vname(i).eq.'errabs  ') errabs = rpar(j)
            if(vname(i).eq.'errrel  ') errel  = rpar(j)
 2899     continue
         if(numspl.eq.0.or.nwspl.eq.0) then
           write(6,*)'spline must be made prior to desmearing'
           goto 2000
         endif
! --- identify the original to get parameters ---
         do 3005 ia=1,nbuf
          if(numspl.eq.numor(ia)) goto 3007
 3005    continue
          write(6,*)'original data not found. unable to set parameters'
          goto 2000
 3007    continue
         call parget('delqv   ',delqv,ia,ier)
         if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
         endif
         ic=nbuf
!
! ---- determine limits and steps of the new-data ----
         qstart = break(1)
         qend   = qmax
         if(qend.ge.break(nwspl)) qend = break(nwspl-1)
         if(nneu.eq.0)     nneu = 100
         if(nneu.gt.mwert) nneu = mwert
         dq  = (qend - qstart) / ( nneu-1 )
         write(6,*)'qstart = ',qstart
         write(6,*)'nwspl=',nwspl,'  break(nwspl) =',break(nwspl)
         write(6,*)'qend   = ',qend,'  dq=',dq
! ---- desmearing of all points ----
         do 3003 i=1,nneu
           qziel  = qstart + (i-1)*dq
           gunten = 0
           goben  = sqrt( break(nwspl-1)**2 - qziel**2 )
!       ---> integrate now:
           call erset (0,0,0)
           call qdng(fdes,gunten,goben,errabs,errel,result,errret)
!ccccc     call qdags(fdes,gunten,goben,errabs,errel,result,errret)
           call erset (0,1,0)
           if(iout.gt.0)write(6,*)i,' errret=',errret
           xwerte(i,ic) = qziel
           ywerte(i,ic) = -delqv/pi * result
 3003    continue
         call txfpar(ia,ic)
         nwert(ic)= nneu
         numor(ic)= mod(numspl,numpls) + 8*numpls
         isels(1) = ic
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'ia      '.or.comand.eq.'i-abso  ') then
!                    ---> correction absolut i <-------
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia = isels(1)
           if(nbuf.lt.mbuf)then
           nbuf=nbuf+1
           ib=nbuf
        endif
         if(yname(ia).ne.'i      ') then
           write(6,*)'date are no intensity',yname(ia)
           goto 2000
         endif
         call parget('i0      ',ai0  ,ia,ier)
         call parget('d       ',d    ,ia,ier)
         call parget('delqv   ',delqv,ia,ier)
         call parget('delqh   ',delqh ,ia,ier)
         call parget('lambda  ',alam0 ,ia,ier)
         write(6,*)'i-absolut  i0 =',ai0
         write(6,*)'           d  =',d, ' angstroem'
         write(6,*)'       lambda =',alam0,'  "  '
         write(6,*)'       delqv  =',delqv,' angstroem**-1'
         write(6,*)'       delqh  =',delqh,'    "         '
         fia = (2*pi/alam0)**2 / (ai0*d*delqv*delqh)
         write(6,*)'       normalizing factor =',fai
         do 2021 i=1,nwert(ia)
             xwerte(i,ib)=xwerte(i,ia)
             ywerte(i,ib)= fia*ywerte(i,ia)
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+5*numpls
 2021    continue
         yname(ib)='i-abso  '
         write(6,*)'q1=',ywerte(1,ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'clip      ') then
!                    ---> remove points
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          ifrom = intval('from    ',1,inew)
          ito   = intval('to      ',0,inew)
          if(ifound('last    ').ne.0)then
           ilas = intval('last    ',1,inew)
           ito  = nwert(ia)
           ifrom= ito-ilas+1
          endif
          if(ito   .lt.1) ito   = 1
          if(ifrom .lt.1) ifrom = 1
          if(ito   .gt.nwert(ia) ) ito  = nwert(ia)
          if(ifrom .gt.nwert(ia) ) ifrom= nwert(ia)
          write(6,*)'removing points ',ifrom,ito,' from ',i
          errmax = getval('errmax  ',1d33,inew)
! flag the points to be removed !
          if(found('from    ') .or.                                     &
     &       found('to      ') .or.                                     &
     &       found('last    ')) then
           do j=ifrom,ito
             yerror(j,ia) = 3e33
           enddo
          endif
! and remove the points
          nnw = 0
          do j=1,nwert(ia)
           if(found('rel     ')) then
             errmx = errmax*abs(ywerte(j,ia))
           else
             errmx = errmax
           endif
           if(yerror(j,ia).lt.errmx) then
            nnw = nnw+1
            xwerte(nnw,ia) = xwerte(j,ia)
            ywerte(nnw,ia) = ywerte(j,ia)
            yerror(nnw,ia) = yerror(j,ia)
           endif
          enddo
          nwert(ia) = max(nnw,1)
          write(6,*)'number of points left (min=1) : ',nwert(ia)
         enddo

         goto 2000
        endif
!
!
       if(comand.eq.'fun     '.or.comand.eq.'function') then
!                    ---> y --> fun(y)
         ifunx = 0
         ifuny = 0
         do 2701 i=1,inames
           j = inapa(i)
           if(vname(i).eq.'op       ') fonpla = .true.
           if(vname(i).eq.'np       ') fonpla = .false.
           if(vname(i).eq.'x        ') ifunx=0
           if(vname(i).eq.'re(x)    ') ifunx=1
           if(vname(i).eq.'log(x)   ') ifunx=2
           if(vname(i).eq.'exp(x)   ') ifunx=3
           if(vname(i).eq.'x**2     ') ifunx=4
           if(vname(i).eq.'sqrt(x)  ') ifunx=5
           if(vname(i).eq.'x*       ') then
                                       ifunx=6
                    if(inpar(i).ne.0)  fx   = rpar(j)
           endif
           if(vname(i).eq.'x+       ') then
                                       ifunx=7
                    if(inpar(i).ne.0)  sx   = rpar(j)
           endif
           if(vname(i).eq.'rouse    ') then 
                                       ifunx=8
                                       ifuny=-1
           endif
           if(vname(i).eq.'zimm     ') ifunx=9

           if(vname(i).eq.'y        ') ifuny=0
           if(vname(i).eq.'re(y)    ') ifuny=1
           if(vname(i).eq.'log(y)   ') ifuny=2
           if(vname(i).eq.'exp(y)   ') ifuny=3
           if(vname(i).eq.'y**2     ') ifuny=4
           if(vname(i).eq.'sqrt(y)  ') ifuny=5
           if(vname(i).eq.'y*       ') then
                                       ifuny=6
                    if(inpar(i).ne.0)  fy   = rpar(j)
           endif
           if(vname(i).eq.'y+       ') then
                                       ifuny=7
                    if(inpar(i).ne.0)  sy   = rpar(j)
           endif
           if(vname(i).eq.'deff     ') ifuny=8
 2701    continue
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
! --- treat all selected items ---
         do 25021 ii=1,nsel
           ia = isels(ii)
           if(fonpla) then
!                     ----> treat data on place
             ib = ia
           else
!          ----> look for new place for the treated data ---
             if(nbuf.lt.mbuf)then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
!          -----(of if(fonpla...

           call txfpar(ia,ib)

           do 2502 i=1,nwert(ia)
             if(ifunx.eq.0)          xwerte(i,ib)=xwerte(i,ia)
             if(ifunx.eq.1) then
               if(xwerte(i,ia).ne.0) xwerte(i,ib)=1 / xwerte(i,ia)
             endif
             if(ifunx.eq.2) then
               if(xwerte(i,ia).gt.0) xwerte(i,ib)= alog(xwerte(i,ia))
             endif
             if(ifunx.eq.3) then
                                     xwerte(i,ib)= exp(xwerte(i,ia))
             endif
             if(ifunx.eq.4) then
                                     xwerte(i,ib)=  xwerte(i,ia)**2
             endif
             if(ifunx.eq.5) then
               if(xwerte(i,ib).ge.0.) then
                                     xwerte(i,ib)=  sqrt(xwerte(i,ia))
               else
                                   xwerte(i,ib)=  0
               endif
             endif
             if(ifunx.eq.6) then
                                     xwerte(i,ib)=  fx * xwerte(i,ia)
             endif
             if(ifunx.eq.7) then
                                     xwerte(i,ib)=  sx + xwerte(i,ia)
             endif

                                 ! Rouse Scaling !
             if(ifunx.eq.8) then
              if(i.eq.1) then
               diffrouse = 0
               wl4 = getval('wl4     ',dble(wl4),inew)
               call  parget('wl4     ',wl4,ia,ier)
               call  parget('q       ',qrs,ia,ier)
               if(ier.ne.0) then
                 call errsig(1111,'q-value not found ...$')
               endif
               call parset('wl4     ',wl4,ib)
               cbuffer = chrval('usediff  ','none    ',inew)
               if(cbuffer(1:4) .ne. 'none' ) then
                  diffscal  = getval('unit    ', 1d7, ier)   ! default diffusion in cm**2/s ==> scale to A**2/ns
                 call parget(cbuffer(1:8),diffrouse,ia,ier)
                 diffrouse = diffrouse * diffscal
                 write(6,'("scaling:",i6,"  at q= ",f12.6," with diff ",f12.6," A**2/ns  from:",a)')ia,qrs,diffrouse,cbuffer(1:8)
               endif
              endif
!!              write(6,*)i,ia,diffrouse,qrs,xwerte(i,ia), exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              ywerte(i,ib) = ywerte(i,ia) * exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              yerror(i,ib) = yerror(i,ia) * exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              xwerte(i,ib) = (qrs**2)*sqrt(wl4*xwerte(i,ia))
             endif

                                 ! Zimm Scaling !
             if(ifunx.eq.9) then
              if(i.eq.1) then
               eta_s = getval('eta     ',dble(eta_s),inew)
               temp_z= getval('temp    ',dble(temp_z),inew)
               call  parget('eta_solv  ',eta_s,ia,ier)
               call  parget('temp      ',temp_z,ia,ier)
               call  parget('q       ',qrs,ia,ier)
               if(ier.ne.0) then
                 call errsig(1111,'q-value not found ...$')
               endif
               call parset('eta_solv  ',eta_s,ib)
               call parset('temp      ',temp_z,ib)
              endif
              ywerte(i,ib) = ywerte(i,ia)
              yerror(i,ib) = yerror(i,ia)
              xwerte(i,ib) =                                            &
     &                ((qrs**3)*akb_z*temp_z/(6*pi*eta_s)*              &
     &                xwerte(i,ia))**(2.0/3.0)
             endif

             if(ifuny.eq.0) then
                                     ywerte(i,ib)=ywerte(i,ia)
                                     yerror(i,ib)=yerror(i,ia)
             endif
             if(ifuny.eq.1) then
               if(ywerte(i,ia).ne.0) then
                  ywerte(i,ib)=    1         / ywerte(i,ia)
                  yerror(i,ib)= yerror(i,ia) / ywerte(i,ia)**2
               endif
             endif
             if(ifuny.eq.2) then
               if(ywerte(i,ia).gt.0) then
                  ywerte(i,ib)= alog(ywerte(i,ia))
                  yerror(i,ib)= yerror(i,ia)/abs(ywerte(i,ia))
               endif
             endif
             if(ifuny.eq.3) then
                  ywerte(i,ib)= exp(ywerte(i,ia))
                  yerror(i,ib)= exp(ywerte(i,ia))*yerror(i,ia)
             endif
             if(ifuny.eq.4) then
                  ywerte(i,ib)= ywerte(i,ia) **2
                  yerror(i,ib)= 2 * yerror(i,ia) * abs(ywerte(i,ia))
             endif
             if(ifuny.eq.5) then
               if(ywerte(i,ib).ge.0.) then
                  ywerte(i,ib)=  sqrt(ywerte(i,ia))
                  if(ywerte(i,ia).gt.0) then
                   yerror(i,ib) = yerror(i,ia) / (2*ywerte(i,ib))
                  else
                   yerror(i,ib) = yerror(i,ia)
                  endif
               else
                  ywerte(i,ib)=  0
                  yerror(i,ib)=  0
               endif
             endif
             if(ifuny.eq.6) then
                ywerte(i,ib)=  fy * ywerte(i,ia)
                yerror(i,ib)=  fy * yerror(i,ia)
             endif
             if(ifuny.eq.7) then
                ywerte(i,ib)=  sy + ywerte(i,ia)
                yerror(i,ib)=       yerror(i,ia)
             endif
             if(ifuny.eq.8) then
               ywerte(i,ib) = 0
               if(ywerte(i,ia).gt.0) then
                if(xwerte(i,ia).ne.0) then
                  ywerte(i,ib) = 1.0/(ywerte(i,ia)*xwerte(i,ia)**2)
                  yerror(i,ib) = yerror(i,ia)                           &
     &                          /(ywerte(i,ia)*xwerte(i,ia))**2
                endif
               endif
             endif
 2502      continue
!          call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+(3+ifuny+10*ifunx)*numpls
           lxx = laenge(xname(ia),80,' ')
           lyy = laenge(yname(ia),80,' ')
           if(ifunx.eq.0) xname(ib) = xname(ia)
           if(ifunx.eq.1) xname(ib) = 're('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.2) xname(ib) = 'ln('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.3) xname(ib) = 'exp'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.4) xname(ib) = 'sqr('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.5) xname(ib) = 'sqrt'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.6) xname(ib) = 'fx*'//xname(ia)(1:lxx)
           if(ifunx.eq.7) xname(ib) = 'sx+'//xname(ia)(1:lxx)
           if(ifunx.eq.8) xname(ib) =                                   &
     &                      'q**2*sqrt(Wl4*'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.9) xname(ib) =                                   &
     &               '(q**3*kT/(6*eta*pi)*'//xname(ia)(1:lxx)//')**2/3'
           if(ifuny.eq.0) yname(ib) = yname(ia)
           if(ifuny.eq.1) yname(ib) = 're('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.2) yname(ib) = 'ln('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.3) yname(ib) = 'exp'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.4) yname(ib) = 'sqr('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.5) yname(ib) = 'sqrt'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.6) yname(ib) = 'fy*'//yname(ia)(1:lyy)
           if(ifuny.eq.7) yname(ib) = 'sy+'//yname(ia)(1:lyy)
           if(ifuny.eq.8) yname(ib) =                                   &
     &          '1/('//xname(ia)(1:lxx)//'**2 *'//yname(ia)(1:lyy)//')'
           call parget('ifunx   ',val,ib,ier)
           ival = val+0.1
           val  = ival*10 + ifunx
           call parset('ifunx   ',val,ib)
           call parget('ifuny   ',val,ib,ier)
           ival = val+0.1
           val  = ival*10 + ifuny
           call parset('ifuny   ',val,ib)
           if(ifunx.eq.6) then
             call parget('fx      ',val,ib,ier)
             if(ier.eq.1) val = 1
             val = val * fx
             call parset('fx      ',val,ib)
           endif
           if(ifuny.eq.6) then
             call parget('fy      ',val,ib,ier)
             if(ier.eq.1) val = 1
             val = val * fy
             call parset('fy      ',val,ib)
           endif
           if(ifunx.eq.7) then
             call parget('sx      ',val,ib,ier)
             if(ier.eq.1) val = 0
             val = val + sx
             call parset('sx      ',val,ib)
           endif
           if(ifuny.eq.7) then
             call parget('sy      ',val,ib,ier)
             if(ier.eq.1) val = 0
             val = val + sy
             call parset('sy      ',val,ib)
           endif
           isels(ii) = ib
           ifits(ii) = 0
           write(6,*)'item:',ia,' treated by fun: ',yname(ib),' --> ',ib
25021    continue
         goto 2000
        endif
!
!
       if(comand.eq.'funfun  ') then
           if(.not.found('immediate')) then
             open(15,file='formdat',status='UNKNOWN')
             read(15,'(a)',end=15999) xformel
             read(15,'(a)',end=15999) yformel
             close(15)
15999        continue
           endif
           write(6,*)'treat x by :',xformel
           write(6,*)'treat y by :',yformel

           fonpla = found('op      ')

         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
! --- treat all selected items ---
         do   ii=1,nsel
           ia = isels(ii)
           iibuf = ia
           if(fonpla) then
!                     ----> treat data on place
             ib = ia
           else
!          ----> look for new place for the treated data ---
             if(nbuf.lt.mbuf)then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
           iibuf = ia
           iadda = ia
!          -----(of if(fonpla...
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               yyee = yerror(i,ia)
               call evaluate(xformel,val8x,ierx)
               call evaluate(yformel,val8y,iery)
               xwerte(i,ib) = val8x
               ywerte(i,ib) = val8y
! -- evaluate error --
               yyyy = ywerte(i,ia)+0.5*yerror(i,ia)
               call evaluate(yformel,val8y,iery)
               yyyy = ywerte(i,ia)-0.5*yerror(i,ia)
               call evaluate(yformel,val8x,iery)
               yerror(i,ib) = dabs( val8y - val8x)
           enddo
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+9*numpls
           xname(ib)  = '!'//xformel(1:7)
           yname(ib)  = '!'//yformel(1:7)
           coment(ib) = 'x='//xformel(1:30)//'y='//yformel(1:30)
           isels(ii) = ib
           ifits(ii) = 0
          enddo
         goto 2000
        endif
!
!
       if(comand.eq.'seterr  ') then
!           write(6,*)'treat x by :',xformel
           write(6,*)'create y-errors by :',yformel

          if(nsel.eq.0) then
            write(6,*)'no curve selected'
            goto 2000
          endif
! --- treat all selected items ---
          do   ii=1,nsel
           ia = isels(ii)
           iadda = ia
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               yyee = yerror(i,ia)
               call evaluate(yformel,val8y,iery)
               yerror(i,ia) = val8y
           enddo
          enddo
         goto 2000
        endif
!
       if(comand.eq.'invers  ') then
!                    ---> data ---> data**-1 vs q**2
         if(nsel.eq.0) then
           write(6,*)'no curve seleted'
           goto 2000
         endif
        ia = isels(1)
        if(xname(ia).ne.'q      '.or.numor(ia).lt.3*numpls) then
           write(6,*)'treat data with q !!'
           goto 2000
         endif
         do 2051 i=1,inames
           j = inapa(i)
           if(vname(i).eq.'bkgr     ') bkgr  = rpar(j)
 2051    continue
! ---- look for an additional data-slot ---
           if(nbuf.lt.mbuf)then
           nbuf=nbuf+1
           ib=nbuf
        endif
!
         do 2052 i=1,nwert(ia)
           xwerte(i,ib) = xwerte(i,ia)**2
           if(ywerte(i,ia).eq.0.or.ywerte(i,ia)-bkgr.le.0) then
              yerror(i,ib) = 1e33
              ywerte(i,ib) = 0
            else
              ywerte(i,ib) = 1 / (ywerte(i,ia)-bkgr)
              yerror(i,ib) = yerror(i,ia)/(ywerte(i,ia)-bkgr)**2
            endif
 2052    continue
!
         call txfpar(ia,ib)
         numor(ib)=mod(numor(ia),numpls)+4*numpls
!
         xname(ib)='q**2   '
         yname(ib)='i**-1  '
         write(6,*)'ok (bkgr=',bkgr,')'
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
!
       if(comand.eq.'dir     ') then
!                    ---
         len_comm = intval('clength ',len_comm,inew)
         if(len_comm.lt.1 ) len_comm=1
         if(len_comm.gt.80) len_comm=80
         if(inew.ne.0) then
           jl=2
         else
           jl=1
         endif

!
         do 170 i=1,nbuf
           csel = ' '
           do j=1,nsel
             if(i.eq.isels(j)) csel = '!'
             if(i.eq.ifits(j)) csel = '-'
           enddo

           ip = 0
           iy = 1
           pastring = ' '
           fostring = ' '
           if(inames.ge.jl) then
           do j=jl,min(inames,20)
             call parget(vname(j),pxxd,i,ier)
             if(ier.eq.0) then
               ip = ip+1
               dirpan(ip) = vname(j)
               dirpav(ip) = pxxd
             else
              iy=0
             endif
           enddo
            if(ip.gt.0) then
            write(fostring,'(a,i4,a)')'(',ip,'(1x,a8,f12.6))'
            write(pastring,fostring)(dirpan(j),dirpav(j),j=1,ip)
            endif
           endif


           if(iy.ne.0) then
           write(6,171) csel,i,numor(i),name(i),xname(i),yname(i)       &
     &    ,'"'//coment(i)(1:len_comm)//'"',pastring(1:ip*21)
  171    format(1x,a1,i4,':#',i14,' : ',a8,':',a8,' vs ',a8,'>',a,' ',a)
           endif
  170    continue

         if(numspl.ne.0)write(6,1711)numspl
 1711    format(/'     active spline generated from #',i14)
         goto 2000
       endif
!
!
       if(comand.eq.'purge   ') then
!                    -----
         if(ipars.eq.0 .and. inames.eq.0) then
           write(6,*)'nothing removed use arguments no1 no2 .. or all'
           goto 2000
         endif
         if(found('all     ')) then
           nbuf = 0
           write(6,*)'purging all ...'
           do ia=1,mbuf
            do j=1,mwert
             xwerte(j,ia) = 0
             ywerte(j,ia) = 0
             yerror(j,ia) = 0
            enddo
            nwert(ia) = 0
            numor(ia) = 0
            xname(ia) = ' '
            yname(ia) = ' '
            coment(ia) = ' '
            nopar(ia) = 0
           enddo
           goto 18051
         endif
         do 1801 i=1,ipars
           ia = Nint(rpar(i))
           write(6,*)'purging no.: ',ia,' ...'
            do j=1,mwert
             xwerte(j,ia) = 0
             ywerte(j,ia) = 0
             yerror(j,ia) = 0
            enddo
            nwert(ia) = 0
            numor(ia) = 0
            xname(ia) = ' '
            yname(ia) = ' '
            coment(ia) = ' '
            nopar(ia) = 0
 1801    continue
! ---- compress data ----
!   --- find the first gap in the data ...
         do 1803 i=1,nbuf
            if(numor(i).eq.0) then
              m = i-1
              goto 1804
            endif
 1803    continue
 1804    continue
!   --- ok here it is ...
!     --- and shift the rest down now ---
              do 1805 j=i,nbuf
                if(numor(j).ne.0) then
                  m = m + 1
                  call txfera(j,m)
                endif
 1805         continue
          nbuf = m
18051     continue
          write(6,*)'nbuf is now ',nbuf
          do 1807 i=1,msel
            isels(i) = 0
            ifits(i) = 0
 1807     continue
          nsel = 0
          write(6,*)'selections are removed'
         goto 2000
       endif
!
!
       if(comand.eq.'sel     ') then
!                    ---

      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= sel                                                                        ='
       write(6,*)'= select items from the list of loaded data records                          ='
       write(6,*)'=    sel  <i1> <i2> ... <in>                                                 ='
       write(6,*)'=            selects those with dir adresses i1, i2,.. in                    ='
       write(6,*)'=    sel  <i1>  -<i2>                                                        ='
       write(6,*)'=            selects all addresses between i1 and i2                         ='
       write(6,*)'=    sel  all  <parname> <value>  band <value>  | <i1> ..                    ='
       write(6,*)'=            selects all records with parameter parnam value +- band         ='
       write(6,*)'=    sel  fit+                                                               ='
       write(6,*)'=            adds fit curves to the selection                                ='
       write(6,*)'=    sel next <parname> <value>  band <value>                                ='
       write(6,*)'=            selects the next record matching the parameter value            ='
       write(6,*)'=    sel add ......                                     (OR)                 ='
       write(6,*)'=         adds the selctions done via ... to the present list                ='
       write(6,*)'=    sel narrow <parname> <value>  band <value>                              ='
       write(6,*)'=         applies parameter selction to the prsent list (AND)                ='
       write(6,*)'=    sel exclude <parname> <value>  band <value>                             ='
       write(6,*)'=         removes matching records from the list                             ='
       write(6,*)'=    sel exclude numor mod <n>                                               ='
       write(6,*)'=         removes records with mod(#,n)=0 from the list                      ='
       write(6,*)'=    sel                                                                     ='
       write(6,*)'=         clear the list                                                     ='
       write(6,*)'=                                                                            ='
       write(6,*)'=    use dsl  or dir to check selections!                                    ='
       write(6,*)'=============================================================================='
       goto 2000
      endif


          
narrow:   if(found('narrow  ')) then
             if(nsel <= 0) goto 2000
                selparna  = chrval('narrow  ',selparna  ,inew)
                selparval = getval( selparna ,dble(selparval) ,inew)
                selpartol = getval('band    ',dble(selpartol) ,inew)
                m = 0
                do i=1,nsel
                   call parget(selparna, parval_x , isels(i),ier)
                   if(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                   endif
                   nsel = m
                enddo   
             goto 2000
          endif narrow
          
exclude:   if(found('exclude  ')) then
                if(nsel <= 0) goto 2000

                if(found('numor   ')) then
                   l = intval('mod     ',2,inew)
                   if(l==0) l=100
                   m = 0
                   do i=1,nsel
                    if(.not.(mod(numor(isels(i)),l) == 0)) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                    endif
                    nsel = m
                   enddo   
                 goto 2000
                endif   
                
                selparna  = chrval('exclude  ',selparna  ,inew)
                selparval = getval( selparna ,dble(selparval) ,inew)
                selpartol = getval('band    ',dble(selpartol) ,inew)
                m = 0
                do i=1,nsel
                   call parget(selparna, parval_x , isels(i),ier)
                   if(.not.(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol)) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                   endif
                   nsel = m
                enddo   
             goto 2000
          endif exclude
          
          
         if(inames.eq.0 .or. found('add     ')                          &
     &                  .or. found('fit+    ')                          &
     &                  .or. found('next    ')                          &
     &                  .or. found('all     ')                          &
     &                                       ) then
! keep current selection if add
           if(.not.found('add     ')) then
             m = 0
            else
             m = nsel
           endif

! go and look for data records with a given parameter
           if( found('next    ') .or. found('all     ') )then
             selparna  = chrval('next    ','        ',inew)
             selparna  = chrval('all     ',selparna  ,inew)
             selparval = getval( selparna ,dble(selparval) ,inew)
             selpartol = getval('band    ',dble(selpartol) ,inew)
                                   id = 1
             if(found('next    ')) id = isels(nsel)+1
             if(nsel.le.0)         id = 1
             do i=id,nbuf
                  iadda = i
                  call parget(selparna, parval_x , iadda,ier)
                  if(iout.gt.0) write(6,*) i , selparna, parval_x , ier
                  if(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol) then
                    if(found('next    ')) then
                      if(nsel.le.0) nsel = 1
                      if(found('add     ')) then
                        if(nsel.lt.mbuf) nsel = nsel+1
                      endif
                      isels(nsel) = i
                      if(iout.gt.0) write(6,*)'next:',nsel,isels(nsel)
                      goto 17801
                    endif
                    m = m+1
                    isels(m) = i
                    write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                          isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                    if(iout.gt.0) write(6,*)'all ',m,i
                  endif
             enddo
             nsel = m
17801       continue

           else

            do i=1,ipars
             iss      = Nint(rpar(i))
             if(iss > 0 .and. iss <= nbuf ) then
               write(6,*)'select adress   ',iss
             else
               write(6,*)"selected=",iss,"  nbuf=", nbuf
               call errsig(999,"ERROR: sel selcted address is out of range$")
               goto 2000
             endif
             m = m + 1
             isels(m) = iss
             ifits(m) = 0
             if(iss.lt.0.and.i.gt.1) then
               m   = m - 1
               ias = isels(m)
               ies = -iss
               write(6,*)'select adresses from ',ias,' to ',ies
               do l = ias,ies
                isels(m) = l
                ifits(m) = 0
                m = m + 1
                                        enddo
               m = m - 1
             endif
                                enddo
            nsel = m
         endif

          if(found('fit+    ')) then
             nfsel = 0
             do i=1,nsel
                numr = -numor(isels(i))
                do j=1,nbuf
                 if(numr.eq.numor(j)) then
                   write(6,*)i,' matches with: ',j
                   ifits(i) = j
                   isfits(i) = j
                   nfsel = nfsel+1
                 endif
                enddo
             enddo
           endif
           write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
           goto 2000
         endif !   ende if(inames.eq.0 ...........&
         do 7077 i=1,inames
           j = inapa(i)
           k = inpar(i)
           if(vname(i).eq.'fits    ') then
             do 70771 l=1,nsel
               numr = numor(isels(l))
               do 70772 ll=1,nbuf
                 if(-numr.eq.numor(ll)) then
                    ifits(l) = ll
                    write(6,*)'select ',ll,' as fit of scan: ',numr
                    goto 70773
                 endif
70772          continue
70773          continue
70771        continue
           endif
           if(vname(i).eq.'last    ') then
             l = 1
             nnumi(1) = numor(nbuf)
             call search(nnumi,l)
           endif
           if(vname(i).eq.'sc      ') then
             do 181 l=1,k
               nnumi(l) = Nint(rpar(l-1+j))
  181        continue
             call search(nnumi,k)
           endif
           if(vname(i).eq.'sc+     ') then
              do 182 l=1,nsel
               nnumi(l) = numor(isels(l))
  182         continue
             do 183 l=1,k
               kk = l+nsel
               if(kk.gt.minc) goto 183
                 nnumi(kk) = Nint(rpar(l-1+j))
                 iprs      = kk
  183        continue
             call search(nnumi,iprs)
           endif
 7077    continue

          write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
         goto 2000
       endif
!
       if(comand.eq.'putpar  ') then
!                    ------
         do i=1,nsel
            iaddp = isels(i)
            call parset (vname(1),sngl(rpar(1)),iaddp)
         enddo
         goto 2000
       endif
!
       if(comand.eq.'rename  ') then
!                    ------
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rename  xaxis <xaxisname>                                                   '
           write(6,*)'= rename  yaxis <yaxisname>                                                   '
           write(6,*)'= rename  name  <shortname>                                                   '
           write(6,*)'=                                                                             '
           write(6,*)'= max. length of names: ', len(xname(1)), len(yname(1)), len(name(1))
           write(6,*)'= for compatibility reasons try to stay with length of 8                      '
           write(6,*)'=============================================================================='
        endif

             do j=1,inames
               do i=1,nsel
                 iaddp = isels(i)
                 if(vname(j).eq.'xaxis   ') xname(iaddp) = argvals(j+1)(1:len(xname(1)))
                 if(vname(j).eq.'yaxis   ') yname(iaddp) = argvals(j+1)(1:len(yname(1)))
                 if(vname(j).eq.'name    ')  name(iaddp) = argvals(j+1)(1:len( name(1)))
               enddo
               do i=1,nfsel
                 iaddp = isfits(i)
                 if(vname(j).eq.'xaxis   ') xname(iaddp) = argvals(j+1)(1:len(xname(1)))
                 if(vname(j).eq.'yaxis   ') yname(iaddp) = argvals(j+1)(1:len(yname(1)))
                 if(vname(j).eq.'name    ')  name(iaddp) = argvals(j+1)(1:len( name(1)))
               enddo
             enddo
             goto 2000
           endif
!
!
       if(comand.eq.'dispsel '.or.comand.eq.'dsl     ') then
!                    --------                ---
!          write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
   390     format(' selected items: ',i5/                                &
      &           ' scan-address:  numor:      fit-address:'/            &
      &           (1x,i8,5x,i8,5x,i8))


!
!!============ this section is nearly identical to dir an should be cast into a subroutine

         len_comm = intval('clength ',len_comm,inew)
         if(len_comm.lt.1 ) len_comm=1
         if(len_comm.gt.80) len_comm=80
         if(inew.ne.0) then
           jl=2
         else
           jl=1
         endif

         do  l=1,nsel
             i = isels(l)
            csel = ' '
!           do j=1,nsel
!             if(i.eq.isels(j)) csel = '!'
!             if(i.eq.ifits(j)) csel = '-'
!           enddo

           ip = 0
           iy = 1
           pastring = ' '
           fostring = ' '
           if(inames.ge.jl) then
           do j=jl,min(inames,20)
             call parget(vname(j),pxxd,i,ier)
             if(ier.eq.0) then
               ip = ip+1
               dirpan(ip) = vname(j)
               dirpav(ip) = pxxd
             else
              iy=0
             endif
           enddo
            if(ip.gt.0) then
            write(fostring,'(a,i4,a)')'(',ip,'(1x,a8,f12.6))'
            write(pastring,fostring)(dirpan(j),dirpav(j),j=1,ip)
            endif
           endif


           if(iy.ne.0) then
           write(6,'(1x,a1,i4,2H [,i4,1H],2H:#,i14,3H : ,a8,1H:,a8,4H vs ,a8,1H>,a,1x,a)') &
     &     csel,i,ifits(l),numor(i),name(i),xname(i),yname(i)       &
     &    ,'"'//coment(i)(1:len_comm)//'"',pastring(1:ip*21)
           endif
         enddo

         if(numspl.ne.0)write(6,1711)numspl
 

         goto 2000
       endif
!
!
!
       if(comand.eq.'edit    ') then
!                    ----

         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         if(vname(1).eq.'sc      '.or.vname(1).eq.'numor   ') then
          do 499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 498
           endif
  499     continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
  498     continue
         endif
!
         if(vname(1).eq.'n       '.or.vname(1).eq.'name    ') then
          do 599 i=1,nbuf
           if(vname(2).eq.name(i)) then
             ispc = i
             goto 598
           endif
  599     continue
          write(6,*)'data with name ',vname(2),' not found!'
          ierrs = 201
          goto 2000
  598     continue
         endif
         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
! ----------> write data onto buffer <---------------------------------
         call savdat('datbuf  ',ispc)
! ----------> enter the system editor with that file <-----------------
!         call system('emacs datbuf')
         call system('open -eW datbuf')
! --- reread it to the same place ---
         nbuff = nbuf
         nbuf = ispc - 1
         inames = 3
         ipars = 0
         argvals(1) = 'datbuf  '
         call input
         nbuf = nbuff
         goto 2000
       endif
!
!
!
       if(comand.eq.'save    ') then
!                    ----
         fsname = 'lastsave'
         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         do i=1,inames
           if(trim(vname(i)).eq.'to'.and.i.lt.inames) then
              fsname = vname(i+1)
              exit
           endif
         enddo
         do j=1,inames
         if(vname(1).eq.'sc      '.or.vname(j).eq.'numor   ') then
          do 4499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 4498
           endif
 4499     continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
 4498     continue
         endif
         enddo
!
         if(vname(1).eq.'n       '.or.vname(i).eq.'name    ') then
           do i=1,nbuf
            if(vname(2).eq.name(i)) then
             ispc = i
             goto 4598
            endif
           enddo
          write(6,*)'data with name ',vname(2),' not found!'
          ierrs = 201
          goto 2000
 4598     continue
         endif
         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
! ----------> write data onto buffer <---------------------------------
         call savdat(fsname,ispc)
         goto 2000
       endif
!
       if(comand.eq.'msave   ') then
!                    -----
         fsname = 'lastsave'
         if(inames.gt.0)   fsname = argvals(1)
         call msavdat(fsname)
         goto 2000
       endif
!
       if(comand.eq.'thc     ') then
!                    ---
         npoint = rpar(1) + 0.001
         call couple(1)
         call thc(npoint)
         goto 2000
       endif
!
       if(comand.eq.'fit     ') then
         call fit
         goto 2000
       endif

       if(comand.eq.'ga_fit  ') then
          call ga_fit
          goto 2000
       endif

       if(comand.eq.'paraout     ') then
         call outputparams
         goto 2000
       endif
!
        if(comand.eq.'th_init  ') then
                                write(*,*) 'not yet implemented, sorry will coming soon'
        !call init_theories(thenam,thparn,nthpar,thrapar,thramin,thramax,mth,mtpar)
         goto 2000
       endif
!
       if(comand.eq.'yfitform ') then
         if(ioldc.ne.0) then
           yfitform = reslin
           ioldc   = 0
         endif
         write(6,*)'y-fit-formel=',trim(yfitform)
         goto 2000
       endif
!
       if(comand.eq.'yformel  ') then
!                    -------
!cc      yformel = title
         if(ioldc.ne.0) then
           yformel  = reslin
           ioldc   = 0
         endif
         write(6,*)'y-formel=',trim(yformel )
         goto 2000
       endif
!
       if(comand.eq.'xformel  ') then
!                    -------
!cc      xformel = title
         if(ioldc.ne.0) then
           xformel  = reslin
           ioldc   = 0
         endif
         write(6,*)'x-formel=',trim(xformel )
         goto 2000
       endif
!
!
       if(comand.eq.'open    ') then
!                    ----
         npax = 1
         npay = 1
         fpnam = vname(1)
         xpnam = vname(2)
         if(inpar(2).ne.0) npax  = rpar(inapa(2)) + 0.0001
         ypnam = vname(3)
         if(inpar(3).ne.0) npay  = rpar(inapa(3)) + 0.0001
         open(33,file=fpnam,status='UNKNOWN')
         write(6,*)'file ',fpnam,' a has been opened for write'
         write(6,*)'data stored on write: ',xpnam,'(',npax,')  vs ',    &
     &              ypnam,'(',npay,')'
         write(33,*)title
         write(33,*)fpnam,'    ',ypnam,'  vs   ',xpnam,'     111111'
         write(33,*)
         write(33,*)'values'
         goto 2000
       endif
!
!
       if(comand.eq.'write   ') then
!                    -----
!        --- take it from the first selected file ---
         if(nsel.eq.0) then
           write(6,*)'no files selected !'
           goto 2000
         endif
         iadd = isels(1)
         iaddf= ifits(1)
!          --- look first for parameters ---
            call parget(xpnam,xval,iadd,ier)
            if(ier.ne.0) call fpaget(xpnam,xval,npax,ier)
            if(ier.ne.0) then
              write(6,*)'parameter xpnam(',npax,'):',xpnam,' not found'
              goto 2000
            endif
            call parget(ypnam,yval,iadd,ier)
            if(ier.ne.0) call fpaget(ypnam,yval,npay,ier)
            if(ier.ne.0) then
              write(6,*)'parameter ypnam(',npay,'):',ypnam,' not found'
              goto 2000
            endif
            write(6,*)'xpnam=',xval,'   ypnam=',yval
            write(33,*)xval,yval
         goto 2000
       endif
!
       if(comand.eq.'close   ') then
!                    -----
         write(33,*)
         write(33,'(4h#eod)')
         close(33)
         write(6,*)'file ',vname(1),' a has been closed'
         goto 2000
       endif
!
!
       if(comand.eq.'title   '.or.comand.eq.'tit     ') then
!                    -----
         title = trim(inline(6:))
         goto 2000
       endif
!
       if(comand.eq.'theos   ') then
!                    -----> list available theories
                        write(6,*)' ***** theories available *****'
                        ileng=len_trim(vname(1))
                        if (ileng.gt.0) then
                                do i=1,mth
                                        if (thenam(i)(:ileng).eq.vname(1)(:ileng)) then
                                                write(6,*)'______________________________________'
                                                write(6,'(i3,": ",a8,i3)') i,thenam(i),nthpar(i)
                                                write(6,*)(trim(thparn(j,i))//' ',j=1,mtpar)
                                        endif
                                enddo
                        else
                                do i=1,mth
                                        if (thenam(i).ne.' ') then  !   kein name!!!
                                                write(6,*)'-------------------------------------'
                                                write(6,'(i3,": ",a8,i3)') i,thenam(i),nthpar(i)
                                                write(6,*)(trim(thparn(j,i))//' ',j=1,mtpar)
                                        endif
                                enddo
                        endif

                        goto 2000
      endif
!
       if(comand.eq.'activate'.or.comand.eq.'ac      ') then
!                    -----> activate a theory
         call activa(0)
         call activa(2)
         goto 2000
       endif
!
      if(comand.eq.'label   ') then
!                   -----> assign a label to a parameter
         call lsearch(jpar,itcal,ierr)
         thpala(jpar,itcal) = vname(3)(1:4)
         goto 2000
      endif

      if(comand.eq.'chgthpar') then
!                   --------> change single theory parameters
         call lsearch(jpar,itcal,ierr)
         if(ierr .eq. 0) then
           ithc  = nthtab(itcal)
           thparx(jpar,itcal) =                                         &
     &        getval('par     ',dble(thparx(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)')                          &
     &                'change parameter ',thparn(jpar,ithc),            &
     &                ' of th(',itcal,'):',thenam(ithc),' to ',         &
     &                thparx(jpar,itcal)
           endif
           thpsca(jpar,itcal) =                                         &
     &        getval('scale   ',dble(thpsca(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)')                          &
     &                'change scale of  ',thparn(jpar,ithc),            &
     &                ' of th(',itcal,'):',thenam(ithc),' to ',         &
     &                thpsca(jpar,itcal)
           endif
         endif
         goto 2000
      endif
!
      if(comand.eq.'couple  ') then
!                   ------> couple label to a parameter
         call lsearch(jpar,itcal,ierr)
         if(ierr.ne.0) goto 2000
         if(inames.lt.3) then
           write(6,*)'label must be specified'
           ierrs =804
           goto 2000
         endif
         ncpl = ncoup(jpar,itcal)
         do 8011 i=1,ncpl
            if(thpalc(i,jpar,itcal).eq.vname(3)(1:4)) then
               thpafc(i,jpar,itcal) = rpar(inapa(3))
               goto 2000
            endif
 8011     continue
          if(ncpl.eq.mcoup) then
            write(6,*)'too many coupled items'
            ierrs = 805
           else
            ncpl = ncpl+1
            ncoup(jpar,itcal) = ncpl
            thpalc(ncpl,jpar,itcal) = vname(3)(1:4)
            thpafc(ncpl,jpar,itcal) = rpar(inapa(3))
           endif
           call couple(0)
         goto 2000
      endif
!
!
!
       if(comand.eq.'acl     '.or.comand.eq.'aclast  ') then
!                    -----> reactivate theories as stored in lastth
         call activa(3)
         if(ierrr.eq.0. .and. iout.ge.0) call activa(2)
         goto 2000
       endif
!
       if(comand.eq.'desactiv'.or.comand.eq.'dac     ') then
!                    -----> desactivate all theories
         call activa(1)
         goto 2000
       endif
!
       if(comand.eq.'activlst'.or.comand.eq.'al      ') then
!                    -----> list activated theories
         call activa(2)
         goto 2000
       endif
!
       if(comand.eq.'get_th'.or.comand.eq.'gth      ') then
!                    -----> get theory appended to an file created by msave and copy it to lassth
!                           and load it
         if(inames == 1) then
!           i = LEN_TRIM(argvals(1))
!           write(6,*)"argvals: ",argvals(1)(1:i)
!           call extract_th(argvals(1)(1:i)//" ")
           call extract_th(argvals(1))
           if(ierrr==0) then 
             call activa(3)
             if(ierrr==0) call activa(2)
           endif
         else
           call errsig(999,"need exactly one filename as argument!$")
         endif
         goto 2000
       endif
!

       if(comand.eq.'gplot    ') then
!                    -----> plot selected curves
         call gplot()
         goto 2000
       endif
!
       if(comand.eq.'plot    '.or.comand.eq.'p       '.or.comand.eq.'gp      ') then
!                    -----> plot selected curves
         call splot(.true.)
!         ibild1 = ibild
         goto 2000
       endif
!
!
!
       if(comand.eq.'plot0   '.or.comand.eq.'p0      ') then
!                    -----> set parameters for plot
         call splot(.false.)
!         ibild1 = ibild
         goto 2000
       endif
!
!
!
       if(comand.eq.'zero   ') then
!                    ----      ------> logical zero the common storage
         nbuf= 0
         goto 2000
       endif
!
!
       if (comand.eq.'numorpls') then
!                    -----------> change offset between files
         numpls = Nint(rpar(1)) 
         goto 2000
       endif

       if (comand.eq.'numorchg') then
!                    -----------> change numor
         if(nsel > 0) then
           write(6,*)"chnage numor of record:",isels(1)," from ",numor(isels(1))," to ", Nint(rpar(1)) 
           numor(isels(1)) = Nint(rpar(1))
         endif
 
         goto 2000
       endif
!
!
       if (comand.eq.'cailtab ') then
!                    -----------> create a table for caille scattering
!                                 to enable complete resolution fitting

          if(nsel.le.0) then
              write(6,*)'need a dummy selection!'
              ierrs = 999
              goto 2000
           endif
          if(ntheos.ne.1) then
              write(6,*)'only one theory spec. is allowed'
              ierrs = 999
              goto 2000
           endif

           ith = nthtab(1)
           if(thenam(ith)(1:4).ne.'cail') then
              write(6,*)'theory(',ith,'): ',thenam(ith),                &
     &        '  does not belong to the cail-family'
              ierrs = 999
              goto 2000
           endif

          write(6,*)'Using theory(',ith,'): ',thenam(ith)

          qortho_0 = getval('qortho_0 ',qortho_0,inew)
          write(6,*)qortho_0, inew
          qortho_1 = getval('qortho_1 ',qortho_1,inew)
          write(6,*)qortho_1, inew
          northo   = intval('northo   ',northo,  inew)
          write(6,*)northo, inew
          qz_0     = getval('qz_0     ',qz_0,    inew)
          write(6,*)qz_0, inew
          qz_1     = getval('qz_1     ',qz_1,    inew)
          write(6,*)qz_1, inew
          nz       = intval('nz       ',nz,      inew)
          write(6,*)nz, inew
          filext   = chrval('to       ',filext,  inew)
          write(6,*)filext, inew

          filnam   = 'cailtab.'//filext

          write(6,*)'Writing Results to: ',filnam
          open(85,file=filnam,form='formatted',status='unknown')

          write(6,'(A,F18.7)')'qortho_0   =  ',qortho_0
          write(6,'(A,F18.7)')'qortho_1   =  ',qortho_1
          write(6,'(A,F18.7)')'qz_0       =  ',qz_0
          write(6,'(A,F18.7)')'qz_1       =  ',qz_1
          write(6,'(A,I11)')  'northo     =  ',northo
          write(6,'(A,I11)')  'nz         =  ',nz

          write(85,'(A)') filnam
          write(85,'(A,F18.7)')' ',qortho_0
          write(85,'(A,F18.7)')' ',qortho_1
          write(85,'(A,F18.7)')' ',qz_0
          write(85,'(A,F18.7)')' ',qz_1
          write(85,'(A,I11)')  ' ',northo
          write(85,'(A,I11)')  ' ',nz

          call couple(1)



         do iortho=1,northo

           qortho=qortho_0+(iortho-1)*((qortho_1-qortho_0)/(northo-1))

           do isel=1,nsel
             iadda = isels(isel)
             write(6 ,*)'qortho = ',qortho,'  PARSET to sel: ',iadda
             call parset('qortho  ' ,qortho,iadda)
                 !! isel
           enddo
           write(85,*)'qortho = ',qortho


           do iz=1,nz
             if(isigint().gt.0) then
               write(6,*)'Table computation stopped by Crtl-C'
               goto 20009
             endif
             qz = qz_0+(iz-1)*((qz_1-qz_0)/(nz-1))
             x = qz
             yctab = thval(x)
             write(85,*) yctab

                 !! iz
           enddo

               !! iortho
         enddo

20009    continue

         call theo_out(85)
         close(85)
         write(6,*)'Ready, Data written to: ',filnam

         goto 2000
       endif
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TOF-Stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(comand.eq.'tofdos  ') then
!                    ---> remove points
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif

         Ommax  = getval('omegamax',Ommax,inew)
         Nomega = intval('nomega  ',Nomega,inew)       

         write(6,*)'Preparing Pseudo Density of states on grid'
         write(6,*)'with ', Nomega,' points'
         write(6,*)'spanning to: ',Ommax,' GHz '

         do i=1,nsel
          ia = isels(i)
          if(nbuf.lt.mbuf) then 
            nbuf = nbuf+1
          else
            Write(6,*)'max. number of buffers reached .. '
            goto 2000   
          endif

          call  tof_dos(ia,nbuf,Ommax,Nomega,ierrr)
          if(ierrr.ne.0) goto 2000
          
          isels(i) = nbuf
           
         enddo

         goto 2000
        endif
!

       if(comand.eq.'detsens  ') then
!                    ---> sensitivity correction
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif


         He3_pressure_bar  = getval('p       ',dble(He3_pressure_bar) ,inew)
         Tube_diameter_m   = getval('dtube   ',dble(Tube_diameter_m),inew)       
         lenfrac           = getval('lenfrac ',dble(lenfrac),inew) 

         write(6,'(a)')'Performing wavelength dependent sensitivy correction on Tof data on lambda scale'
         write(6,*)' p       = He3-pressure/bar        = ',  He3_pressure_bar 
         write(6,*)' dtube   = tube diameter/m         = ',  Tube_diameter_m 
         write(6,*)' lenfrac = projected length factor = ',  lenfrac 

         do i=1,nsel
          ia      = isels(i)
          call parget('det_phe3',xx,ia,ier)
          if(ier .ne. 0 ) then
            do j=1,nwert(ia)
              lambdaA =  xwerte(j,ia)
              if(lambdaA > 0.01 ) then
                 detsens = detector_sensitivity(lambdaA,lenfrac,He3_pressure_bar,Tube_diameter_m)
                 ywerte(j,ia) = ywerte(j,ia)/detsens
                 yerror(j,ia) = yerror(j,ia)/detsens
 !                write(6,*)xwerte(j,ia), detsens
              endif
            enddo
            call parset ('det_phe3 ', He3_pressure_bar ,ia)
            call parset ('dtube    ', Tube_diameter_m  ,ia)
          else
            write(6,*)'Record: ',ia,' will not be treated a second time, skipped!'
          endif
         enddo

         goto 2000
        endif
!
!
       if(comand.eq.'gdos_q  ') then
!                    ------> see e.g. J.C. Smith et al, J.Chem. Phys. 85, 3636 (1986)
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
         ia      = isels(1)

         k1                = intval('k1      ', k1,      inew)
         k2                = intval('k2      ', k2,      inew)
         recstep           = intval('recstep ', recstep, inew)

         write(6,'(a)')'Creating q-exptrapolation curve for g(omega) deterimation '
         write(6,'(a)')'see e.g. J.C. Smith et al, J.Chem. Phys. 85, 3636 (1986)  '
         write(6,*)' from channel ',k1,' to ', k2,'   recstep ',recstep
         write(6,*)' lambda       ',xwerte(k1,ia),' to ',xwerte(k2,ia)

! use first record as template for energy scaling, parameter extraction etc.
         call parget('det_phe3',xx,ia,ier)
         if(ier .ne. 0 ) then
           write(6,*)'do detector efficiency correction first!'
           goto 2000
         endif
  
         call parget ('lambda  ',lambdaA ,ia, ier)
         if(ier.ne.0) then
            write(6,*)'ERROR: Parameter lambda is missing in parameterlist of record!'
            goto 2000
         endif

         call parget ('temp    ', temp ,ia, ier)
         if(ier.ne.0) then
            write(6,*)'ERROR: Parameter temp is missing in parameterlist of record!'
            goto 2000
         endif

         lambda0      = (xwerte(k1,ia)+xwerte(k2,ia))/2     * 1d-10                      ! das ist lambda_final
         lambda1      = lambdaA                             * 1d-10                      ! das ist lambda_in
         delta_lambda = (xwerte(k2+1,ia)-xwerte(k1,ia))     * 1d-10
         e0           = NeutronEnergy_fromLambda(dble(lambda0))
         e1           = NeutronEnergy_fromLambda(dble(lambda1))
         dE           = e1-e0                                    ! hquer*omega = E_in - E_final
         de0     = NeutronEnergy_fromLambda(dble(lambda1-delta_lambda/2))-NeutronEnergy_fromLambda(dble(lambda1+delta_lambda/2))
         dee     = NeutronEnergy_fromLambda(dble(lambda0-delta_lambda/2))-NeutronEnergy_fromLambda(dble(lambda0+delta_lambda/2))



         omega   = dE*2*Pi/Planckkonstante / 1d9
         om_cm   = dE/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1

         kinetic_factor = lambda0/lambda1      ! to be multiplied with the cross setion data to come towards s(q,om)
         channel_factor = dee/de0
         channel_width  = dee/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1
         channel_width0 = de0/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1

          
         if(nbuf.lt.mbuf) then 
            nbuf = nbuf+1
         else
            Write(6,*)'max. number of buffers reached .. '
            goto 2000   
         endif
         call txfera(ia,nbuf)
         numor(nbuf)    = numor(ia)+100000
         xname(nbuf)    = 'q      '
         yname(nbuf)    = 'gdos'
         coment(nbuf)   = coment(nbuf)(1:40)//'gdos extrapolation data '
         nwert(nbuf) = 0

         call parset ('omega    ', omega            ,nbuf)
         call parset ('om_cm    ', om_cm            ,nbuf)
         call parset ('rwidth   ', channel_factor   ,nbuf)
         call parset ('cwidth   ', channel_width    ,nbuf)
         call parset ('cwidth0  ', channel_width0   ,nbuf)
         call parset ('lambda1  ', xwerte(k1,ia)    ,nbuf)
         call parset ('lambda2  ', xwerte(k2,ia)    ,nbuf)
         
         
         do i=0,nsel-recstep,recstep
          sumq = 0
          sum  = 0
          sume = 0
 irs:     do irstp=1,recstep
           ia      = isels(i+irstp)
           call parget ('angle   ',angle_2tht, ia, ier)  
           if(ier.ne.0) then
             write(6,*)'ERROR: Parameter angle is missing in parameterlist of record!',isels(ia)
           endif

          do j=k1,k2
            lambda0      = xwerte(k1,ia)                       * 1d-10                      ! das ist lambda_final
            lambda1      = lambdaA                             * 1d-10                      ! das ist lambda_in
            delta_lambda = (xwerte(k2+1,ia)-xwerte(k1,ia))     * 1d-10
            e0           = NeutronEnergy_fromLambda(dble(lambda0))
            e1           = NeutronEnergy_fromLambda(dble(lambda1))
            dE           = e1-e0                                    ! hquer*omega = E_in - E_final
            q  = sqrt(((2*Pi/lambda0)**2+(2*Pi/lambda1)**2 - 2*COS(angle_2tht*Pi/180d0)*(2*Pi/lambda0)*(2*Pi/lambda1)))  
            qa = q * 1d-10
            sumq = sumq+qa

            omega   = dE*2*Pi/Planckkonstante / 1d9
            om_cm   = dE/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1
            kinetic_factor = lambda0/lambda1      ! to be multiplied with the cross setion data to come towards s(q,om)
!                                                   kinetic_factor=1 if already applied
            dbb = exp(dble(-dE/(temp*Boltzmannkonstante))) - 1.0
!                         ! definition of omega in paper
            sum  = sum  +  6*(-omega)/(qa**2)*dbb* kinetic_factor * ywerte(j,ia)
            sume = sume + (6*(-omega)/(qa**2)*dbb* kinetic_factor * yerror(j,ia))**2
!                             ! definition of omega in paper          
           enddo
          enddo irs

           nwert(nbuf) =  nwert(nbuf) + 1
           xwerte(nwert(nbuf), nbuf) = sumq       / (recstep*(k2-k1+1))
           ywerte(nwert(nbuf), nbuf) = sum        / (channel_factor*recstep)
           yerror(nwert(nbuf), nbuf) = sqrt(sume) / (channel_factor*recstep)
           write(6,'(2i4,3e14.7)')i, nwert(nbuf), xwerte(nwert(nbuf), nbuf), ywerte(nwert(nbuf), nbuf),  yerror(nwert(nbuf), nbuf)
         enddo

         nsel     = 1
         isels(1) = nbuf
         ifits    = 0
        
         goto 2000
        endif
!

!============================================== Checkings =================
       if(comand.eq.'argvals  ') then
         write(6,*)' inmames   = ', inames
         write(6,*)' argvals(1)= ', argvals(1)(1:60), ' ... '
         write(6,*)' argvals(2)= ', argvals(2)(1:60), ' ... '
         write(6,*)' argvals(3)= ', argvals(3)(1:60), ' ... '
         write(6,*)' argvals(4)= ', argvals(4)(1:60), ' ... '
         write(6,*)' argvals(5)= ', argvals(5)(1:60), ' ... '
         goto 2000
       endif








!
       call makro(comand)
!
       goto 2000
!      ---------> read now commandlines from makro file
!
!
!


      END ! Ende der Hauptschleife







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_pair_of_points(xin1,yin1,yer1,n1, xin2,yin2,yer2,n2, xout, yout1, yerout1, yout2, yerout2, nout)
!      ===============================================================================================================
!
       implicit none

       real, intent(in)     :: xin1(n1), yin1(n1), yer1(n1)
       real, intent(in)     :: xin2(n2), yin2(n1), yer2(n2)
       integer, intent(in)  :: n1, n2

       real, intent(out)    :: xout(nout), yout1(nout), yerout1(nout), yout2(nout), yerout2(nout) 
       integer              :: nout
         

       integer              :: i, j, jc, ii
!?       real                 :: x, xx, y1, y2, e1, e2, dist,p,y,ye
       real                 :: x, xx, dist,p,y,ye
 
       nout = 0
! assuming nontrivial length check
       if(n1.lt.2 .or. n2.lt.2) then
         write(6,*)'at least one of the vectors has too few components:',n1,n2
         call errsig(990,"ERROR: vectors dont length-match$")
         return
       endif
! assuming equal length check
       if(n1.ne.n2) then
         write(6,*)'vectors have different number of components:',n1,n2
         call errsig(991,"ERROR: vectors dont length-match$")
         return
       endif
! assuming ordered vectors CHECK
       do i=1,n1-1
        if(xin1(i).ge.xin1(i+1)) then
          write(6,*)'vector 1 is not ordered ',i,xin1(i),xin1(i+1)
          call errsig(992,"ERROR: vectors must be x-sequentially ordered $")
          return
        endif
       enddo
       do i=1,n2-1
        if(xin2(i).ge.xin2(i+1)) then
          write(6,*)'vector 2 is not ordered ',i,xin2(i),xin2(i+1)
         call errsig(992,"ERROR: vectors must be x-sequentially ordered $")
          return
        endif
       enddo

! do the assignments !
       do i=1,n1
!      search for matching points       
!      1. closest mate
         dist = abs(xin1(i)-xin2(1))
         jc   = 1
         do j=2,n2
           xx = abs(xin1(i)-xin2(j))
           if(xx.lt.dist) then
              dist = xx
              jc   = j
           endif
         enddo


!        now our closest neighbour pair is i<-->jc
!        if they are equal we are done
         if(xin1(i).eq.xin2(jc)) then
           nout = nout + 1
           xout(nout)    = xin1(i)
           yout1(nout)   = yin1(i)
           yerout1(nout) = yer1(i)
           yout2(nout)   = yin2(jc) 
           yerout2(nout) = yer2(jc)
           cycle
        endif

        if(xin1(i).lt.xin2(jc)) then
          do ii=1,n1
            if(xin1(ii).gt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)  !! this would be good if full statistical independeence
                                                           !! is assumed, however interpolation creates correlation
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif

        if(xin1(i).gt.xin2(jc)) then
          do ii=n1,1,-1
            if(xin1(ii).lt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif
         
       enddo
       

       return
       end subroutine get_pair_of_points



 
   subroutine average_data()
!     -------------------------
      use dimensions
      use cincom
      use cincoc
      use xoutxx
      use cdata
      use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
      use constants
      use PhysicalConstantsPlus
 
      implicit none

      double precision, save        :: xcatch = 0.05d0
      integer                       :: inew, ier
!?      integer                       :: i,j,k,n
      integer                       :: i,j,n
      integer                       :: number_of_data_points
      integer                       :: n_result_point
      character(len=8)              :: cbuf
      logical                       :: found
      double precision              :: getval
      double precision              :: xaver, yaver, yerra
      double precision, allocatable :: xva(:), yva(:), yvaer(:)
      integer         , allocatable :: iperm(:)
      logical         , allocatable :: va_used(:)


      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= average_data                                                               ='
       write(6,*)'= combine error weighted close data points from records in the selection list='
       write(6,*)'= parameter is       xcatch  the relative width of the x-window              ='
       write(6,*)'= current default is xcatch  ',xcatch
       write(6,*)'=============================================================================='
       return
      endif





      xcatch = getval("xcatch  ",xcatch,inew)

      write(6,*)'xcatch=', xcatch

      if(nsel <= 1) then
        write(6,*)"Select more than one curve to average!"
        ierrs = 1000
        return
      endif

      if(nbuf >= mbuf) then
        write(6,*)"..too many buffers", nbuf, mbuf
        ierrs = 2000
        return
      endif

! prepare destination
      nbuf = nbuf + 1
      call txfera(isels(1),nbuf)
      numor(nbuf)    = numor(isels(1))+200000
      call parset("xcatch   ",sngl(xcatch),nbuf)
      do i=1, nsel
        write(cbuf,'("c",i0)')i
        call parset(cbuf,real(numor(isels(i))),nbuf)
      enddo
! check how many data_points are there
      number_of_data_points = 0
      do i=1,nsel
       number_of_data_points = number_of_data_points + nwert(isels(i)) 
      enddo

      if(allocated(xva))     deallocate(xva)
      if(allocated(yva))     deallocate(yva)
      if(allocated(yvaer))   deallocate(yvaer)
      if(allocated(va_used)) deallocate(va_used)
      if(allocated(iperm))   deallocate(iperm)
      allocate(xva    (number_of_data_points))
      allocate(yva    (number_of_data_points))
      allocate(yvaer  (number_of_data_points))
      allocate(va_used(number_of_data_points))
      allocate(  iperm(number_of_data_points))

! distribute data
      n = 0      
      do i=1,nsel
        do j=1,nwert(isels(i))
         n=n+1
         xva(n)    =  xwerte(j,isels(i))
         yva(n)    =  ywerte(j,isels(i))
         yvaer(n)  =  yerror(j,isels(i))
         va_used(n)= .false.
        enddo 
      enddo
      number_of_data_points = n
! error-sorting (using SLATEC dpsort) 
      iperm = [(i,i=1,number_of_data_points)]
      call  DPSORT (yvaer, number_of_data_points, IPERM,  1 , ier)
!C      KFLAG - input -- control parameter:---------------^
!C            =  2  means return the permutation vector resulting from
!C                  sorting DX in increasing order and sort DX also.
!C            =  1  means return the permutation vector resulting from
!C                  sorting DX in increasing order and do not sort DX.
!C            = -1  means return the permutation vector resulting from
!C                  sorting DX in decreasing order and do not sort DX.
!C            = -2  means return the permutation vector resulting from
!C                  sorting DX in decreasing order and sort DX also.

       write(6,*)            " All related points after error sorting: "
       write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,number_of_data_points)

! now proceed with the lowest error point and get all neibours in xcatch distance
       n_result_point = 0
d1:    do i=1,number_of_data_points
          if(va_used(iperm(i))) cycle d1
          va_used(iperm(i)) = .true.
          xaver = xva(iperm(i)) * 1d0/yvaer(iperm(i))**2
          yaver = yva(iperm(i)) * 1d0/yvaer(iperm(i))**2
          yerra = 1d0/yvaer(iperm(i))**2
d2:       do j=1,number_of_data_points
            if(va_used(iperm(j))) cycle d2
            if(abs(2*(xva(iperm(j))-xva(iperm(i)))/(xva(iperm(j))+xva(iperm(j)))) > xcatch ) cycle d2 ! relative match
              xaver = xaver + xva(iperm(j)) * 1d0/yvaer(iperm(j))**2
              yaver = yaver + yva(iperm(j)) * 1d0/yvaer(iperm(j))**2
              yerra = yerra + 1d0/yvaer(iperm(j))**2   
              va_used(iperm(j)) = .true.           
          enddo d2 
          n_result_point = n_result_point + 1
          if(n_result_point > mwert) then
            write(6,*)'too many points created.'
            return
          endif
          xwerte(n_result_point,nbuf) = xaver / yerra
          ywerte(n_result_point,nbuf) = yaver / yerra
          yerror(n_result_point,nbuf) = sqrt(1d0 / yerra)

!          write(6,*)"test: ", n_result_point, xwerte(n_result_point,nbuf), &
!                     ywerte(n_result_point,nbuf), yerror(n_result_point,nbuf)

       enddo d1
       nwert(nbuf) =  n_result_point
!sort the result according to x
       iperm = [(i,i=1, n_result_point)]

!       write(6,*)            " intermediate result: "
!       write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,n_result_point)

       call  SPSORT (xwerte(:,nbuf),n_result_point, IPERM, -1 , ier)
       xva  (1:n_result_point) = xwerte(iperm(1:n_result_point),nbuf)
       yva  (1:n_result_point) = ywerte(iperm(1:n_result_point),nbuf)
       yvaer(1:n_result_point) = yerror(iperm(1:n_result_point),nbuf)
       xwerte(1:n_result_point,nbuf) =  xva  (1:n_result_point)
       ywerte(1:n_result_point,nbuf) =  yva  (1:n_result_point) 
       yerror(1:n_result_point,nbuf) = yvaer (1:n_result_point)

       write(6,'(a)'    , advance='no') "final result from: "
       write(6,'(20i4)' , advance='no') (isels(i),i=1,nsel)
       write(6,'(a,i4,a,f12.6)')" --> ",nbuf,"    xcatch=",xcatch  
       write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,n_result_point)

       nsel     = 1
       isels(1) = nbuf



   end subroutine average_data



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! help functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine command_list()
  write(6,*) " command name               " , " parameters: "
  write(6,*) "?         (synonym) help    " , " " 
  write(6,*) "acl       (synonym) aclast  " , " "
  write(6,*) "activate  (synonym) ac      " , " <theory_name> "
  write(6,*) "activlst  (synonym) al      " , " "
  write(6,*) "addsels                     " , " "
  write(6,*) "addval                      " , " "
  write(6,*) "argvals                     " , " "
  write(6,*) "arit                        " , " f1 <factor1>  f2 <factor2> to <numor> [div | mult | (no)norm ]"
  write(6,*) "arit2                       " , " f1 <factor1>  f2 <factor2> to <numor> [div | mult | (no)norm ]"
  write(6,*) "average                     " , " "
  write(6,*) "chgthpar                    " , " "
  write(6,*) "clip                        " , " "
  write(6,*) "close                       " , " "
  write(6,*) "combine                     " , " "
  write(6,*) "copy                        " , " "
  write(6,*) "couple                      " , " "
  write(6,*) "cs        (synonym) clearsel" , " "
  write(6,*) "des       (synonym) desmear " , " "
  write(6,*) "desactiv  (synonym) dac     " , " "
  write(6,*) "detsens                     " , " "
  write(6,*) "dir                         " , " "
  write(6,*) "dispsel   (synonym) dsl     " , " "
  write(6,*) "dmx                         " , " "
  write(6,*) "echocurv                    " , " "
  write(6,*) "edit                        " , " "
  write(6,*) "ercorrc   (synonym) ecc     " , " "
  write(6,*) "fft                         " , " "
  write(6,*) "fftmx     (synonym) fft-ms  " , " "
  write(6,*) "fftrand                     " , " "
  write(6,*) "fit                         " , " "
  write(6,*) "fun       (synonym) function" , " "
  write(6,*) "funfun                      " , " "
  write(6,*) "ga_fit                      " , " "
  write(6,*) "gaiscale                    " , " "
  write(6,*) "gdos_q                      " , " "
  write(6,*) "gen.res                     " , " "
  write(6,*) "get_th  (synonym) gth       " , " "
  write(6,*) "gplot      (synonym) gp     " , " "
  write(6,*) "hiqextra                    " , " "
  write(6,*) "ia        (synonym) i-abso  " , " "
  write(6,*) "in        (synonym) input   " , " "
  write(6,*) "inscn                       " , " "
  write(6,*) "integrat                    " , " "
  write(6,*) "interpol                    " , " "
  write(6,*) "invers                      " , " "
  write(6,*) "kz                          " , " "
  write(6,*) "label                       " , " "
  write(6,*) "m         (synonym) mirror  " , " "
  write(6,*) "msave                       " , " "
  write(6,*) "mux                         " , " "
  write(6,*) "noise                       " , " "
  write(6,*) "numorchg                    " , " "
  write(6,*) "numorpls                    " , " "
  write(6,*) "open                        " , " "
  write(6,*) "out_gli   (synonym) gli     " , " "
  write(6,*) "paraout                     " , " "
  write(6,*) "parextra                    " , " "
  write(6,*) "plot      (synonym) p       " , " "
  write(6,*) "plot0     (synonym) p0      " , " "
  write(6,*) "purge                       " , " "
  write(6,*) "putpar                      " , " "
  write(6,*) "qc        (synonym) q-conv  " , " "
  write(6,*) "rebin                       " , " "
  write(6,*) "recsum                      " , " "
  write(6,*) "rename                      " , " "
  write(6,*) "rerange                     " , " "
  write(6,*) "save                        " , " "
  write(6,*) "scale                       " , " "
  write(6,*) "sel                         " , " "
  write(6,*) "sequence                    " , " "
  write(6,*) "seterr                      " , " "
  write(6,*) "sp        (synonym) spline  " , " "
  write(6,*) "swapxy                      " , " "
  write(6,*) "sym                         " , " "
  write(6,*) "th_init                     " , " "
  write(6,*) "thc                         " , " "
  write(6,*) "theos                       " , " "
  write(6,*) "title     (synonym) tit     " , " "
  write(6,*) "tofdos                      " , " "
  write(6,*) "tracorr                     " , " "
  write(6,*) "uni_ft                      " , " "
  write(6,*) "write                       " , " "
  write(6,*) "xformel                     " , " "
  write(6,*) "yfitform                    " , " "
  write(6,*) "yformel                     " , " "
  write(6,*) "z         (synonym) zero    " , " "
  write(6,*) "zero                        " , " "
end subroutine command_list
