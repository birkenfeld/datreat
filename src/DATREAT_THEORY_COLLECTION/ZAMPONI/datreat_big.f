c*ds
c*ds
c **********************************************************************
c *    treating any data                                               *
c *    imsl version 10 fitprogram test                                 *
c *                                                                    *
c *    datreat version from 18.5.88                                    *
c *    m. monkenbusch, iff, kfa-juelich                                *
c *    neu 3.11.93: multiply flag fuer theorien                        *
c *    neu 28.1.99: thc -n x1 <x1> x2 <x2> on log scale                *
c *    neu 09.3.99: fun zimm, rouse etc + error                        *
c *                 plot noparplo , achsen lange strings               *
c *                 plot error bars bei log ok                         *
c *                 purge loescht tatsaechlich (auch errors)           *
c *                 sel add                                            *
c *    neu 04.02.00:seterr                                             *
c *                 generiere Fehler aus Daten (Fitweighting)          *
c *                 ueber yformel                            (RCS:1.3) *
c *                                                                    *
c *    neu 23.02.00:fit map ..                                         *
c *    neu 19.05.00:dparam('name   ') fuer th-programmierung verfuegbar*
c *    neu         :msave             fuer Mehrfachsave                *
c **********************************************************************
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c --- mwert  = max. no. of x-y-values in one buffer
c     mbuf   = max. no. of different buffers
c     mpar   = max. no. of parameters associated with one buffer
c ---  maximum scan length ....
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)
c ---  fit dimensions ---
       parameter (mfit=40,msmpl=4000)
c  -- mfit = max no. of fitted parameters
c     msmpl= max no. of datapoints in fit
c
c --- incom common-section ---
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar, rparf, getval, valnxt, get1, get2, get3, dble
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       real*8 xyorig, rotvec
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c --- xwerte(i,j)   x-values on buffer j
c     ywerte(i,j)   y-values on buffer j
c     yerror(i,j)   error of y-values   (only supported by some fktn)
c     xname(j)      name of x-values (x-axis) for buffer j
c     yname(j)       "   "  y-   "    y-  "    "    "    "
c     name(j)       short text-identifier for data on buffer j
c     nwert(j)      no. of valid points on buffer j
c     nbuf          no. of filled buffers
c     numor(j)      numerical idenfication of data on buffer j
c     coment(j)     one line of comment describing data on buffer j
c     params(l,j)   set of parameters associated with data on buffer j
c     napar(l,j)    names of these parameters
c     nopar(j)      no. of valid parameters
c
c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c  isels(i)  = address (/spectr/) of selected scan
c  ifits(i)  = adress of fitted spectrum (isels(i))
c  nsel      = length of this table
c
       common/fslist/isfits(mbuf),nfsel
c  isfits    = address (/spectr/) of selected fits
c  nfsel     = length of this table
c
c ----- theories common block and definitions ----
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
c  thenam(i)  = name of i-th theory
c  thparn(j,i)= name of j-th parameter of i-th theory
c  nthpar(i)  = no of parameters required for i-th theory
c  thparx(j,l)= parameters for the l-th activated theory
c  thpsca(j,l)= corresponding fit scales   "
c  nthtab(l)  = no. of l-th activated theory
c  ntheos     = total no. of activated theories
c  multflg    = flag indicating a multiplicative theory
c
c ------ coupling of theory-parameters -----------------------
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
c  --- thpala(i,j) ... label of i-th parameter , j-th theory
c      thpalc(l,i,j) . labels (l) that point to parameters that are
c                      coupleded to the i-th parameter of the j-th theo.
c      thpafc(l,i,j) . associated proportionality factor
c      thpaco(i,j) ... offset value for a parameter that is derived by
c                      couplede other parameters
c      ncoup(i,j) .... number of parameters that are coupled to i,j-par.
c
c ------ errors of fit --------------------------------------------
       common/therrc/therro(mtpar,mtcal)
c      ------ estimates of 1 sigma errros -----------------------
 
       common/cfc/qziel,cscoef(4,mwert),break(mwert),weight(mwert),
     *            numspl,nwspl
c --- parameters of last spline smoothing + value qziel
c     qziel is the value at witch the spline should be evaluated
c     numspl numor of spline fitted data
c     nwspl  length of splined data vectors
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c xxxx,yyyy = aktuelle addersse fuer wertextraktion xx,yy
c
c ---- communication with subr. func ---
      logical sqwght,sqwbuf
      logical autox1,autox2
      common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)
c  --- iprt     : printing during func calculation (set by fit)
c      x1       : lower limit of fit range (if autox1 = .false.)
c      x2       : upper limit of fit range (if autox2 = .false.)
c      autox1/2 : if true the genuine limits of the datafiled is taken
c                 else x1 /x2 values are taken (set by thc)
c      ferror   : fehler des zu vergleichenden datenfeldes
c ----------------------------------------------------------------------
       real    fcssq
       logical lerrel, lwrtfitdat
       common/cfunce/ lerrel, lwrtfitdat, fcssq
c ---  if lerrel=.true. the fiterrors are take as relative errors
c
c --- echoform related parameters ( nse-programs ) ----
      real*4 j1echo,j2echo,j0delta
 
      common/partran/j1echo,j2echo,j0delta,cdelta
c     j1, j2  = feldintegrale in gauss*meter
c     j0delta = feldintegralvariation unabh&a.ngig von j1,j2
c     cdelta  = koeffizient der feldintegralabh&a.ngigkeit von j1,j2
 
      common/wlntran/alam0, dalam
c     alam0   = mittlere wellenlaenge in a
c     dalam   = wellenla.ngenbreite in a
 
      common/sqtran/tau
c     tau in sekunden, relaxationszeit des streuers
 
      common/ergsum/bcut1,rcoef1,slope1,bcut2,rcoef2,slope2,yint
c --- ergebnisse der sum-auswertung
c     bcut1,2 --- y-achsschittpunkte der high u. low-freq. geraden
c     slope1,2 -- steigungen              "          "
c     rcoef1,2 -- r-parameter
c     yint     -- abgeleitete intensitaet
 
      common/ergfil/ymx,sigma,ymaxf,yminf,
     *              imx,imx2,imaxf,iminf,ic1,ic2,ic3
c --- filter-sektion (parameter & ergebnisse) ---
c  -- ymx : wert des maximums im spektrum
c     imx : kanal des spektralen max.
c     imx2: spiegel kanal (imx)
c     ymaxf: maximum in der gefilterten huellkurve
c     imaxf: kanallage des maximums ymaxf
c     yminf: minimum in der gefilterten huellkurve
c     iminf: kanallage des maximums ymaxf
c     ic1,2,3 : die letzten 3 angesprochenen kurven
 
c ------------ imsl-version 10 workspace -------------------------------
c
c      parameter (iwklen=11*mfit+2*msmpl-1)
       parameter (iwklen=200000)
c      wrklen for unlsf
c      the system must be told about this length by calling iwkin(wrklen
       common/worksp/ rwksp(iwklen)
c ----------------------------------------------------------------------
c
       character fname*25,fnam*8,ftyp*8,fmode*1,fpnam*8,xpnam*8,ypnam*8
       character fsname*8
       character*1 csel
       character*132 xformel,yformel,yfitform
       common/formul/xformel,yformel,yfitform
       real*8 val8x,val8y
 
       complex ca(mwert), cb(mwert), cmplx, conjg
       dimension alim(2), blim(2)
c
       dimension qq(3),par(10),qt(3),sqt(3),iqt(3)
       dimension nnumi(minc)
       logical withmo,fonpla,found,folgt,compare
 
       parameter( mcut=200 )
       dimension rrv(0:mcut)

       common /thiadd/iadda
c ---> transfer der addresse der gerade bearbeiteten kurve nach thXX
c      damit koennen dann ggf. weitere parameter ueber parget gewonnen
c      werden
c

       double precision qortho_0, qortho_1,qz_0, qz_1
       integer          nz, northo, iortho, iz, ith
       character*40     filnam
       character*8      filext, chrval

       integer          len_comm      ! dir-command --> comment-length
c
       external fdes
       external f
c
c
       data fmode/'a'/,ftyp/'scn     '/,withmo/.true./,fonpla/.false./
       data errabs/1.0/,errel/1.e-2/
c
c ---- initialisations ----
!!!!   call iwkin(iwklen)
!!!!   ------------------> inform imsl programs about the workspacesize
c ---- error-set ----------
       call erset (0,1,0)

       call sigset(-1)  !! activates Crtl-C signal handler
       call sig_Reset() !! Clears Signals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(6,*)
       write(6,*)'====================================================='
       write(6,*)'=   datreat_big     Version: 19.05.00               ='
       write(6,*)'=   -----------     -----------------               ='
       write(6,*)'=   new feature: fit map                            ='
       write(6,*)'=   new                                             ='
       write(6,*)'=   new USE dparam( NAM ) in th-programs to replace ='
       write(6,*)'=           parget                                  ='
       write(6,*)'====================================================='
       write(6,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
       pi   = 4 * atan(1.0)
       iout = 0
       inka = 5
       ibild= 0
       numpls = 10000
       facto1 = 1
       facto2 = 2
c ---- zero the commons ----
       nbuf = 0
       ntyp  = 0
       numspl= 0
       nwspl = 0
c ---- fit and calculation limit determination -------------------------
       autox1 = .true.
       autox2 = .true.
       lerrel = .false.
       lwrtfitdat = .false.
c ---- dummy formulae ----
       xformel = '(xx)'
       yformel = '(yy)'
c ---- noise amplitude ----
       ampli   = 1.0
       iseed   = 12347897
c ---- nse-echoform ----
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
 
c ---- ecc -------------
      n = 200
      ri= 5.0
      r0= 18.9
      d0= 0.2
      a = 1.2e-3
      ncut = 7
!!!!!!!!!!
       msel = 0
!!!!!!!!!!
!-------- cailtab -------------
       northo   = 11
       nz       = 11
       qortho_0 = 0.0001d0
       qortho_1 = 0.1001d0
       qz_0     = 0.01d0
       qz_1     = 0.15d0
       filext   = 'last'
c
!-------- fun Rouse, Zimm Scaling -----
       wl4      = 1.0
       eta_s    = 1e-3
       akb_z    = 1e21*1.380662e-23  ! fuer A**-3ns
c
   
!--------------------------------------
       len_comm = 60                 ! comment-length in dir (display)

       title='Datreat Plot'


c
c ------ mask some of the runtime errors -------------------------------
!      call errdev
c      -----------> this routine reads the error-numbers etc. from
c                   file errmsk a
c ----------------------------------------------------------------------
c
c ---- initialize theory tables ----
       ini = 0
       nthpar(1) = mtpar
       dummy = thdatr1(x,par,thenam(1),thparn(1,1),nthpar(1),ini)
       nthpar(2) = mtpar
       dummy = th2(x,par,thenam(2),thparn(1,2),nthpar(2),ini)
 
       nthpar(3) = mtpar
       dummy = th3(x,par,thenam(3),thparn(1,3),nthpar(3),ini)
       nthpar(4) = mtpar
       dummy = th4(x,par,thenam(4),thparn(1,4),nthpar(4),ini)
       nthpar(5) = mtpar
       dummy = th5(x,par,thenam(5),thparn(1,5),nthpar(5),ini)
       nthpar(6) = mtpar
       dummy = th6(x,par,thenam(6),thparn(1,6),nthpar(6),ini)
       nthpar(7) = mtpar
       dummy = th7(x,par,thenam(7),thparn(1,7),nthpar(7),ini)
       nthpar(8) = mtpar
       dummy = th8(x,par,thenam(8),thparn(1,8),nthpar(8),ini)
       nthpar(9) = mtpar
       dummy = th9(x,par,thenam(9),thparn(1,9),nthpar(9),ini)
       nthpar(10) = mtpar
       dummy = th10(x,par,thenam(10),thparn(1,10),nthpar(10),ini)
       nthpar(11) = mtpar
       dummy = th11(x,par,thenam(11),thparn(1,11),nthpar(11),ini)
       nthpar(12) = mtpar
       dummy = th12(x,par,thenam(12),thparn(1,12),nthpar(12),ini)
       nthpar(13) = mtpar
       dummy = th13(x,par,thenam(13),thparn(1,13),nthpar(13),ini)
       nthpar(14) = mtpar
       dummy = th14(x,par,thenam(14),thparn(1,14),nthpar(14),ini)
       nthpar(15) = mtpar
       dummy = th15(x,par,thenam(15),thparn(1,15),nthpar(15),ini)
       nthpar(16) = mtpar
       dummy = th16(x,par,thenam(16),thparn(1,16),nthpar(16),ini)
       nthpar(17) = mtpar
       dummy = th17(x,par,thenam(17),thparn(1,17),nthpar(17),ini)
       nthpar(18) = mtpar
       dummy = th18(x,par,thenam(18),thparn(1,18),nthpar(18),ini)
       nthpar(19) = mtpar
       dummy = th19(x,par,thenam(19),thparn(1,19),nthpar(19),ini)
       nthpar(20) = mtpar
       dummy = th20(x,par,thenam(20),thparn(1,20),nthpar(20),ini)
       nthpar(21) = mtpar
       dummy = th21(x,par,thenam(21),thparn(1,21),nthpar(21),ini)
       nthpar(22) = mtpar
       dummy = th22(x,par,thenam(22),thparn(1,22),nthpar(22),ini)
       nthpar(23) = mtpar
       dummy = th23(x,par,thenam(23),thparn(1,23),nthpar(23),ini)
       nthpar(24) = mtpar
       dummy = th24(x,par,thenam(24),thparn(1,24),nthpar(24),ini)
       nthpar(25) = mtpar
       dummy = th25(x,par,thenam(25),thparn(1,25),nthpar(25),ini)
       nthpar(26) = mtpar
       dummy = th26(x,par,thenam(26),thparn(1,26),nthpar(26),ini)
       nthpar(27) = mtpar
       dummy = th27(x,par,thenam(27),thparn(1,27),nthpar(27),ini)
       nthpar(28) = mtpar
       dummy = th28(x,par,thenam(28),thparn(1,28),nthpar(28),ini)
       nthpar(29) = mtpar
       dummy = th29(x,par,thenam(29),thparn(1,29),nthpar(29),ini)
       nthpar(30) = mtpar
       dummy = th30(x,par,thenam(30),thparn(1,30),nthpar(30),ini)
       nthpar(31) = mtpar
       dummy = th31(x,par,thenam(31),thparn(1,31),nthpar(31),ini)
       nthpar(32) = mtpar
       dummy = th32(x,par,thenam(32),thparn(1,32),nthpar(32),ini)
       nthpar(33) = mtpar
       dummy = th33(x,par,thenam(33),thparn(1,33),nthpar(33),ini)
       nthpar(34) = mtpar
       dummy = th34(x,par,thenam(34),thparn(1,34),nthpar(34),ini)
       nthpar(35) = mtpar
       dummy = th35(x,par,thenam(35),thparn(1,35),nthpar(35),ini)
       nthpar(36) = mtpar
       dummy = th36(x,par,thenam(36),thparn(1,36),nthpar(36),ini)
       nthpar(37) = mtpar
       dummy = th37(x,par,thenam(37),thparn(1,37),nthpar(37),ini)
       nthpar(38) = mtpar
       dummy = th38(x,par,thenam(38),thparn(1,38),nthpar(38),ini)
       nthpar(39) = mtpar
       dummy = th39(x,par,thenam(39),thparn(1,39),nthpar(39),ini)
       nthpar(40) = mtpar
       dummy = th40(x,par,thenam(40),thparn(1,40),nthpar(40),ini)
       ntheos= 0
c
!!     write(6,3001)(i,thenam(i),nthpar(i),
!!   *                   (thparn(j,i),j=1,mtpar),i=1,mth)
c
c
2000   continue
       call sig_Reset()
       if(ierrs.ne.0) then
         write(6,*)' error return code:',ierrs
         ioldc = 0
c        ----------> if error forget rest of commandline (reslin)
         if(inka.ne.5) then
            close(inka)
            inka = 5
c           --------> go back to input from keyboard
         endif
       endif
       ierrs = 0
c

       iibuf= isels(1)
       call incom(comand)
c      ------------------
       iout = iot
c
c
       if(comand.eq.'?       '.or.comand.eq.'help    ') then
c                    -                       ----
         call info
         goto 2000
       endif
c
c
       if(comand.eq.'in      '.or.comand.eq.'input   ') then
c                    --                      -----
         call input
         nsel = 1
         isels(1) = nbuf
         ifits(1) = 0
         goto 2000
       endif
c
       if(comand.eq.'out_gli '.or.comand.eq.'gli     ') then
c                    --                      -----
         call out_gli
         goto 2000
       endif
c
       if(comand.eq.'inscn   ') then
c                    -----
         call inscn
         nsel = 1
         isels(1) = nbuf
         ifits(1) = 0
         goto 2000
       endif
c
c
       if(comand.eq.'cs      '.or.comand.eq.'clearsel') then
c                    --                      --------
2020     continue
c   <------ dieser sprung ist unelegant und nicht systemkonform
c           in der naechsten version subroutine cs verwenden !!
         nsel = 0
         do 333 i=1,msel
          isels(i) = 0
          ifits(i) = 0
333      continue
         write(6,*)'selections are removed ....'
         goto 2000
       endif
c
c
       if(comand.eq.'z       '.or.comand.eq.'zero    ') then
c                    -                       ----
         nbuf = 0
         write(6,*)'nbuf has been resetted ....'
         nsel = 0
         do  i=1,msel
          isels(i) = 0
          ifits(i) = 0
         enddo
         write(6,*)'selections are removed ....'
       endif
c
c
       if(comand.eq.'noise   ') then
c                    -----
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
ccccc    call rnopt(6)
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
c
c
       if(comand.eq.'gen.res ') then
c                    -------
c ---> create a file in in11 spin-echo result format
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
         write(33,'(a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')
     *   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',
     *   bkgr
         write(6,*)'----------------------------------'
         write(6 ,'(a,a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')
     *   ' file=',
     *   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',
     *   bkgr
         write(6,'(1x,a)')coment(iadd)
         do i=1,nn
          write(33,*)xwerte(i,iadd),ywerte(i,iadd)/ynno,
     *                              yerror(i,iadd)/ynno
         enddo
         close(33)
         goto 2000
       endif
c
c
 
       if(comand.eq.'echocurv') then
c                    ---------
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
 
         call  qand(f,2,alim,blim,erraec,errrec,maxint,
     *              result, errest )
         write(6,*)'-------------------------------------------'
         write(6,*)'result =',result
         write(6,*)'errest =',errest
         write(6,*)'-------------------------------------------'
         n0 =     nech/2
         n0     = intval('n0      ',n0    ,inew)
         call     parset('n0      ',(real(n0)) ,ia1)   !! aix.sp
 
         nwert(ia1) = nech
         do i=1,nech
           j2echo = j1echo + (i-n0) * ddj + dj
           xwerte(i,ia1) = j2echo-j1echo
c --- imsl mehrdimensionale integration ---
           call  qand(f,2,alim,blim,erraec,errrec,maxint,
     *                result, errest )
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
c                    -------
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
c --- und filtern ---
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
c --- speziell fuer echo-auswertung ----
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
c --- echoauswertung konventionelle methode ---
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
c                    -------
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
c --- speziell fuer echo-auswertung ----
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
c ---> fuer den filter
         do i=1,n1/2
           ywerte(i,ia4) = 2*atan(famp*exp(-((i-1-n1/4)/fwid)**2))/pi
           if(i.ge.2) then
             ywerte(n1+2-i,ia4) = ywerte(i,ia4)
           endif
         enddo
         ywerte(n1/2+1,ia4) = 0
c <---
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
c --- und filtern ---
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
         yfamp =
     *         ywerte(imaxf,ia5)**2 + ywerte(imaxf+1,ia5)**2
     *    +    ywerte(imaxf,ia5)**2 + ywerte(imaxf-1,ia5)**2
         yfamp = sqrt( yfamp/2.0)
c        ------> filtern: bis hier
c ---- echoauswertung ----
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
c --- 3-pt-konventionell ---
           y1 = igrand(ywerte(n0-2,ia1)*amp*n1/3.0)
           y2 = igrand(ywerte(n0-1,ia1)*amp*n1/3.0)
           y3 = igrand(ywerte(n0,ia1)*amp*n1/3.0 )
           aver = (y1+y3)*0.5
           echo = sqrt( ((y3-y1)*0.5)**2 + (y2-aver)**2 )
           if(iout.gt.0)write(6,*)ira,'  peak:',sump,'  noise:',sumb,
     *                            ' peakf:',yfamp
c
           bkave = bkave + sumb
           bkqave= bkqave+ sumb**2
           pkave = pkave + sump
           pkqave= pkqave+ sump**2
           p3ave = p3ave + echo
           p3qave= p3qave+ echo**2
           pfave = pfave + yfamp
           pfqave= pfqave+ yfamp**2
c
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
 
 
c
       if(comand.eq.'arit    ') then
c                    ----
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
          if(vname(i).eq.'to      ') newnum = rpar(j) * 1.0000001
          if(vname(i).eq.'sc      ') then
            num1 = rpar(j) * 1.0000001
            num2 = rpar(j+1)*1.0000001
          endif
          if(vname(i).eq.'factor1 '.or.vname(i).eq.'f1      ')
     *                                                 facto1 = rpar(j)
          if(vname(i).eq.'factor2 '.or.vname(i).eq.'f2      ')
     *                                                 facto2 = rpar(j)
4711     continue
c --- figuer out the adresses ---
         if(nsel.eq.2) then
           iad1 = isels(1)
           iad2 = isels(2)
         else
           iad1 = 0
           iad2 = 0
         endif
         do 4713 i=1,nbuf
          if(numor(i).eq.num1) iad1 = i
          if(numor(i).eq.num2) iad2 = i
4713     continue
         if(iad1.eq.0) then
           write(6,*)'file :',num1,' not found'
           goto 2000
         endif
         if(iad2.eq.0) then
           write(6,*)'file :',num2,' not found'
           goto 2000
         endif
c --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4715 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
4715     continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.mbuf) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
c ---- normalize to monitor-values ----
         if(withmo) then
           write(6,*)'normalizing to monitor-values ....'
           write(6,*)'use option nonorm to switch off, norm to swit. on'
c ---- extract monitorvalues ----
            xmon1 = 0
           do 4717 i=1,nopar(iad1)
             if(napar(i,iad1).eq.'monitor  ') xmon1 = params(i,iad1)
4717       continue
           xmon2 = 0
           do 4719 i=1,nopar(iad2)
             if(napar(i,iad2).eq.'monitor  ') xmon2 = params(i,iad2)
4719       continue
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
c        ------> of if(withmo ..
c --- start the computation ---
         n1 = nwert(iad1)
         n2 = nwert(iad2)
         n3 = 0
         do 4721 i=1,n1
          x1 = xwerte(i,iad1)
c         -----> look for appropriate values in the seconf file
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
4723      continue
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
           if(idimux.eq.0 ) y  = facto1 * y0/xmon1 +
     *          facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2
           if(idimux.eq.1 ) y  = (facto1 * y0/xmon1) /
     *        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           if(idimux.eq.2 ) y  = (facto1 * y0/xmon1) *
     *        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           ywerte(n3,newadd) = y
           yyy=y
c --- and the errors ! ----
           y0 = yerror(i ,iad1)
           y1 = yerror(j1,iad2)
           y2 = yerror(j2,iad2)
           x1 = xwerte(j1,iad2)
           x2 = xwerte(j2,iad2)
           if(idimux.eq.0) y = (facto1 * y0/xmon1)**2 +
     *          (facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)**2
           if(idimux.eq.1) y = yyy * ((facto1 * y0/xmon1)/y00) *
     *          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           if(idimux.eq.2) y = yyy * ((facto1 * y0/xmon1)/y00) *
     *          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           yerror(n3,newadd) = sqrt(y)
4721     continue
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
4728     continue
         write(6,*)n3,' points resulting from arit are stored as sc',
     *                  newnum
         if(idimux.eq.0)
     *   write(6,*)'<===',facto1,' * ',name(iad1),' + ',facto2,' * ',
     *             name(iad2)
         if(idimux.eq.1)
     *   write(6,*)'<===',facto1,' * ',name(iad1),' / ',facto2,' * ',
     *             name(iad2)
         if(idimux.eq.2)
     *   write(6,*)'<===',facto1,' * ',name(iad1),' * ',facto2,' * ',
     *             name(iad2)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
c
c
c
c
       if(comand.eq.'combine ') then
c                    -------
         newnum = 777777
         do 4811 i=1,inames
          j = inapa(i)
          k = inpar(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'raster  ') then
           x1raster = rpar(j)
           x2raster = rpar(j+1)
           nraster  = rpar(j+2) * 1.00000001d0
          endif
          if(vname(i).eq.'to      ') newnum = rpar(j) * 1.0000001d0
          if(vname(i).eq.'sc      ') then
             do 48181 l=1,k
               nnumi(l) = rpar(l-1+j)*1.0000001d0
48181        continue
             call search(nnumi,k)
          endif
4811     continue
c --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4815 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
4815     continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.mbuf) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
c --- create the raster ----------------------------------------------
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
c --- loop over the selectionlist ----------------------------
         nwert(newadd) = nraster
         do isl=1,nsel
           iadd = isels(isl)
           xmonitor = 1.d0
           call parget ('monitor ',xmonitor,iadd,ier)
           n1 = nwert(newadd)
           n2 = nwert(iadd)
           do 4821 i=1,n1
            x1 = xwerte(i,newadd)
c         -----> look for appropriate values in the second file
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
4823        continue
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
c           ----- temorary misuse of the error-filed to collect
c                 monitor-values
4821       continue
         enddo
c ---- normalize the data and create errors --------
         do i=1,n1
           if(yerror(i,newadd).gt.0.0) then
             yer     = sqrt(ywerte(i,newadd))/ yerror(i,newadd)
             ywerte(i,newadd) = ywerte(i,newadd) / yerror(i,newadd)
             yerror(i,newadd) = yer
           endif
         enddo
c --- transfer parameters ---
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
c
c
c
c
       if(comand.eq.'m       '.or.comand.eq.'mirror  ') then
c                    -                       ------
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
           write(6,*) numor(ia),'=numor > ',numpls,'; warning no origina
     *                data'
         endif
         if(xname(ia).ne.'kanal   '.and.xname(ia).ne.'kanal   ') then
           write(6,*)'warning x-values are not kanal !!!'
         endif
         nw = nwert(ia)
c --- determine the center ---
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
1702     continue
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
c
         call parset ('xk0     ',xk0,ia)
         nopar(nbuf)= npp+1
         do 2007 i=1,npp+1
          params(i,nbuf) = params(i,ia)
          napar (i,nbuf) = napar (i,ia)
2007     continue
         goto 2000
       endif
c
c
c
       if(comand.eq.'ercorrc '.or.comand.eq.'ecc     ') then
c                    -------                 ---
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
c
c
       if(comand.eq.'sym     ') then
c                    ---->  symmetrize
c ---- build now the symmetric average on the right side ---
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
c --- look for xk0 ---
         do 2077 i=1,nopar(ib)
            if(napar(i,ib).eq.'xk0     ') then
              xk0 = params(i,ib)
              npp = i
            endif
2077     continue
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
2010     continue
         nwert(ib) = nw-kan1+1
         xk0            = xk0 - kan1 + 1
         params(npp,ib) = xk0
         numor(ib)=mod(numor(ib),numpls)+2*numpls
         write(6,*)'ok symm-ds has ',nwert(ib),' points'
         goto 2000
       endif
c
c
       if(comand.eq.'kz      ') then
c                    -->    kanalzusammenfassung
c ---- build now the symmetric average on the right side ---
         kz = rpar(1)+0.1
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
         write(6,*)'channels contracted for ',ia,' new no. of ch=',j,
     *             ' (',ir,')'
         enddo
         goto 2000
       endif
c
c
       if(comand.eq.'interpol') then
c                    ---->  interpolate
c ---- build now the symmetric average on the right side ---
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
c
c
       if(comand.eq.'parextra') then
c                    ---> parabola extrapolation to q--> 0
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
c
c
       if(comand.eq.'hiqextra') then
c                    ---> extrapolation for high q .
c                         data == a*q^z
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
c        ---------------
         xcut  = rpar(1)
c        ---------------
         ia    = isels(1)
         do i=1,nwert(ia)
            if(xwerte(i,ia).gt.xcut) goto 77771
         enddo
77771    continue
         icut = i-1
         slope = (ywerte(icut,ia)-ywerte(icut-2,ia))/
     *           (xwerte(icut,ia)-xwerte(icut-2,ia))
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
c
       if(comand.eq.'qc      '.or.comand.eq.'q-conv  ') then
c                    ---> channels to q-values
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
2001     continue
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
2002     continue
         xname(ib)='q      '
         write(6,*)'q1=',xwerte(1,ib),'  qend=',xwerte(nwert(ib),ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
c
       if(comand.eq.'sp      '.or.comand.eq.'spline  ') then
c                    ---> spline smoothing
          do 2799 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'iequal  ') iequal = rpar(j) + 0.1
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'auto    ') iauto  = 1
            if(vname(i).eq.'noauto  ') iauto  = 0
            if(vname(i).eq.'smpar   ') smpar  = rpar(j)
2799      continue
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
c
        if(iauto.eq.1) then
           write(6,*)'iequal= ',iequal
           call csscv(nn,xwerte(1,ia),ywerte(1,ia),iequal,break,cscoef)
        else
           write(6,*)'smpar = ',smpar
           do 3777 i=1,nn
             weight(i) = sqrt(abs(ywerte(i,ia)))
3777       continue
           call cssmh(nn,xwerte(1,ia),ywerte(1,ia),weight,smpar,break,
     *                cscoef)
        endif
         if(iout.gt.1) then
            write(6,*)'break         cscoeff'
            do 2776 i=1,nn
              write(6,2777)i,break(i),(cscoef(j,i),j=1,4)
2777          format(1x,i3,':',5e14.6)
2776        continue
         endif
         x01 = break(1)
         x02 = break(nn)
         dx  = (x02-x01)/(nneu-1)
         do 2702 i=1,nneu
            xn = x01 + (i-1)*dx
            xwerte(i,ib) = xn
            ywerte(i,ib) = fsplin(xn)
2702     continue
         call txfpar(ia,ib)
         numor(ib)= - numor(ia)
         nwert(ib)= nneu
         isels(1) = ia
         ifits(1) = ib
         nsel     = 1
         goto 2000
        endif
c
       if(comand.eq.'fftmx   '.or.comand.eq.'fft-ms  ') then
c                    ---> multiple scattering by fft
           t1 = second()
           if(numspl.eq.0.or.nwspl.eq.0) then
             write(6,*)'spline must be made prior to fft'
             goto 2000
           endif
c -------- identify the original to get parameters ---
           do 3705 ia=1,nbuf
           if(numspl.eq.numor(ia)) goto 3707
3705      continue
           write(6,*)'original data not found. unable to set parameters'
           goto 2000
3707      continue
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
c --- transfer all data to be kept to the new buffer item ---
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
3799      continue
c
         if(imx.eq.0) then
           write(6,*)'evaluate multiple scattering ...'
           call fftmx(xwerte(1,ib),ywerte(1,ib),trans,xmax,dx,nx,nfft)
           numor(ib)= mod(numor(ia),numpls)+10*numpls
           yname(ib)= 'mx('//yname(ia)(1:5)
         else
           write(6,*)'evaluate inverse multiple scattering ...'
           if(imx.eq.1) then
             call fftrmx(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,
     *                                                  xmax,dx,nx,nfft)
           else
             call rmx1d(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,
     *                                                xmax,ymax,nx,nfft)
           endif
           numor(ib)= mod(numor(ia),numpls)+11*numpls
           yname(ib)= 'imx('//yname(ia)(1:4)
           call parset('i0      ',ai0,ib)
         endif
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nx
         isels(1) = ib
         t2 = second()
         write(6,*)'done needed ',t2-t1,' sec cpu'
         goto 2000
        endif
c
       if(comand.eq.'dmx     ') then
c                    ---> multiple scattering by fft
           t1 = second()
          ia = isels(1)
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
c --- transfer all data to be kept to the new buffer item ---
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
c
         call demux(xwerte(1,ia),ywerte(1,ia),nwert(ia),
     *              trans,xmax,nfft,
     *              iout,xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'dmx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = second()
         write(6,*)'done needed ',t2-t1,' sec cpu'
         goto 2000
        endif
c
       if(comand.eq.'mux     ') then
c                    ---> multiple scattering by fft
           t1 = second()
          ia = isels(1)
          if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
            ib=nbuf
          endif
c --- transfer all data to be kept to the new buffer item ---
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
c
         call   mux(xwerte(1,ia),ywerte(1,ia),nwert(ia),
     *              trans,xmax,nfft,
     *              iout,xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'mx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = second()
         write(6,*)'done needed ',t2-t1,' sec cpu'
         goto 2000
        endif
c
       if(comand.eq.'tracorr ') then
          ia = isels(1)
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)
          do i=1,nwert(ia)
            ywerte(i,ia) = ywerte(i,ia)*thick*trans
            yerror(i,ia) = yerror(i,ia)*thick*trans
          enddo
          write(6,*)'selected item:',ia,' multiplied by:',
     *              thick,' * ',trans,' !'
         goto 2000
        endif
c
       if(comand.eq.'des     '.or.comand.eq.'desmear ') then
c                    ---> infinite slit desmearing
          do 2899 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'qmax    ') qmax   = rpar(j)
            if(vname(i).eq.'errabs  ') errabs = rpar(j)
            if(vname(i).eq.'errrel  ') errel  = rpar(j)
2899      continue
         if(numspl.eq.0.or.nwspl.eq.0) then
           write(6,*)'spline must be made prior to desmearing'
           goto 2000
         endif
c --- identify the original to get parameters ---
         do 3005 ia=1,nbuf
          if(numspl.eq.numor(ia)) goto 3007
3005     continue
          write(6,*)'original data not found. unable to set parameters'
          goto 2000
3007     continue
         call parget('delqv   ',delqv,ia,ier)
         if(nbuf.lt.mbuf)then
            nbuf=nbuf+1
         endif
         ic=nbuf
c
c ---- determine limits and steps of the new-data ----
         qstart = break(1)
         qend   = qmax
         if(qend.ge.break(nwspl)) qend = break(nwspl-1)
         if(nneu.eq.0)     nneu = 100
         if(nneu.gt.mwert) nneu = mwert
         dq  = (qend - qstart) / ( nneu-1 )
         write(6,*)'qstart = ',qstart
         write(6,*)'nwspl=',nwspl,'  break(nwspl) =',break(nwspl)
         write(6,*)'qend   = ',qend,'  dq=',dq
c ---- desmearing of all points ----
         do 3003 i=1,nneu
           qziel  = qstart + (i-1)*dq
           gunten = 0
           goben  = sqrt( break(nwspl-1)**2 - qziel**2 )
c       ---> integrate now:
           call erset (0,0,0)
           call qdng(fdes,gunten,goben,errabs,errel,result,errret)
cccccc     call qdags(fdes,gunten,goben,errabs,errel,result,errret)
           call erset (0,1,0)
           if(iout.gt.0)write(6,*)i,' errret=',errret
           xwerte(i,ic) = qziel
           ywerte(i,ic) = -delqv/pi * result
3003     continue
         call txfpar(ia,ic)
         nwert(ic)= nneu
         numor(ic)= mod(numspl,numpls) + 8*numpls
         isels(1) = ic
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
c
       if(comand.eq.'ia      '.or.comand.eq.'i-abso  ') then
c                    ---> correction absolut i <-------
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
2021     continue
         yname(ib)='i-abso  '
         write(6,*)'q1=',ywerte(1,ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
c
       if(comand.eq.'clip      ') then
c                    ---> remove points
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          ifrom = intval('from    ',1,inew)
          ito   = intval('to      ',0,inew)
          if(ifound('last    '))then
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
          if(found('from    ') .or.
     *       found('to      ') .or.
     *       found('last    ')) then
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
c
c
       if(comand.eq.'fun     '.or.comand.eq.'function') then
c                    ---> y --> fun(y)
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
           if(vname(i).eq.'rouse    ') ifunx=8
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
2701     continue
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
c --- treat all selected items ---
         do 25021 ii=1,nsel
           ia = isels(ii)
           if(fonpla) then
c                     ----> treat data on place
             ib = ia
           else
c          ----> look for new place for the treated data ---
             if(nbuf.lt.mbuf)then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
c          -----(of if(fonpla...

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

             if(ifunx.eq.8) then ! Rouse Scaling !
              if(i.eq.1) then
               wl4 = getval('wl4     ',dble(wl4),inew)
               call  parget('wl4     ',wl4,ia,ier)
               call  parget('q       ',qrs,ia,ier)
               if(ier.ne.0) then
                 call errsig(1111,'q-value not found ...$')
               endif
               call parset('wl4     ',wl4,ib)
              endif       
              ywerte(i,ib) = ywerte(i,ia)
              yerror(i,ib) = yerror(i,ia)
              xwerte(i,ib) = (qrs**2)*sqrt(wl4*xwerte(i,ia))
             endif

             if(ifunx.eq.9) then ! Zimm Scaling !
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
              xwerte(i,ib) = 
     *                ((qrs**3)*akb_z*temp_z/(6*pi*eta_s)*
     *                xwerte(i,ia))**(2.0/3.0)
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
                  yerror(i,ib) = yerror(i,ia)
     *                          /(ywerte(i,ia)*xwerte(i,ia))**2
                endif
               endif
             endif
2502       continue
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
           if(ifunx.eq.8) xname(ib) = 
     *                      'q**2*sqrt(Wl4*'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.9) xname(ib) = 
     *               '(q**3*kT/(6*eta*pi)*'//xname(ia)(1:lxx)//')**2/3'
           if(ifuny.eq.0) yname(ib) = yname(ia)
           if(ifuny.eq.1) yname(ib) = 're('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.2) yname(ib) = 'ln('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.3) yname(ib) = 'exp'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.4) yname(ib) = 'sqr('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.5) yname(ib) = 'sqrt'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.6) yname(ib) = 'fy*'//yname(ia)(1:lyy)
           if(ifuny.eq.7) yname(ib) = 'sy+'//yname(ia)(1:lyy)
           if(ifuny.eq.8) yname(ib) = 
     *          '1/('//xname(ia)(1:lxx)//'**2 *'//yname(ia)(1:lyy)//')'
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
c
c
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
           call lowcase(xformel,132)
           call lowcase(yformel,132)
 
           fonpla = found('op      ')
 
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
c --- treat all selected items ---
         do   ii=1,nsel
           ia = isels(ii)
           iibuf = ia
           if(fonpla) then
c                     ----> treat data on place
             ib = ia
           else
c          ----> look for new place for the treated data ---
             if(nbuf.lt.mbuf)then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
           iibuf = ia
c          -----(of if(fonpla...
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               call evaluate(xformel,val8x,ierx)
               call evaluate(yformel,val8y,iery)
               xwerte(i,ib) = val8x
               ywerte(i,ib) = val8y
c -- evaluate error --
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
c
c
       if(comand.eq.'seterr  ') then
           write(6,*)'treat x by :',xformel
           write(6,*)'treat y by :',yformel
           call lowcase(xformel,132)
           call lowcase(yformel,132)
 
          if(nsel.eq.0) then
            write(6,*)'no curve selected'
            goto 2000
          endif
c --- treat all selected items ---
          do   ii=1,nsel
           ia = isels(ii)

           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               call evaluate(yformel,val8y,iery)
               yerror(i,ia) = val8y
           enddo
          enddo
         goto 2000
        endif
c
       if(comand.eq.'invers  ') then
c                    ---> data ---> data**-1 vs q**2
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
2051     continue
c ---- look for an additional data-slot ---
           if(nbuf.lt.mbuf)then
           nbuf=nbuf+1
           ib=nbuf
        endif
c
         do 2052 i=1,nwert(ia)
           xwerte(i,ib) = xwerte(i,ia)**2
           if(ywerte(i,ia).eq.0.or.ywerte(i,ia)-bkgr.le.0) then
              ywerte(i,ib) = 0
            else
              ywerte(i,ib) = 1 / (ywerte(i,ia)-bkgr)
            endif
2052     continue
c
         call txfpar(ia,ib)
         numor(ib)=mod(numor(ia),numpls)+4*numpls
c
         xname(ib)='q**2   '
         yname(ib)='i**-1  '
         write(6,*)'ok (bkgr=',bkgr,')'
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
c
c
       if(comand.eq.'dir     ') then
c                    ---
!         if(ipars.gt.0) then
!c                       ----> this is a nasty trick
!           nnnn = rpar(1) + 0.1
!           write(6,*)'warning: change nbuf from ',nbuf,' to ',nnnn
!           nbuf = nnnn
!         endif

         len_comm = intval('clength ',len_comm,inew)
         if(len_comm.lt.1 ) len_comm=1
         if(len_comm.gt.80) len_comm=80
c
         do 170 i=1,nbuf
           csel = ' '
           do j=1,nsel
             if(i.eq.isels(j)) csel = '!'
             if(i.eq.ifits(j)) csel = '-'
           enddo
           write(6,171)csel,i,numor(i),name(i),xname(i),yname(i)
     *    ,coment(i)(1:len_comm)
171        format(1x,a1,i3,':#',i14,' : ',a8,':',a8,' vs ',a8,'>',a)
170      continue
         if(numspl.ne.0)write(6,1711)numspl
1711     format(/'     active spline generated from #',i14)
         goto 2000
       endif
c
c
       if(comand.eq.'purge   ') then
c                    -----
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
           ia = rpar(i) * 1.000001
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
1801     continue
c ---- compress data ----
c   --- find the first gap in the data ...
         do 1803 i=1,nbuf
            if(numor(i).eq.0) then
              m = i-1
              goto 1804
            endif
1803     continue
1804     continue
c   --- ok here it is ...
c     --- and shift the rest down now ---
              do 1805 j=i,nbuf
                if(numor(j).ne.0) then
                  m = m + 1
                  call txfera(j,m)
                endif
1805          continue
          nbuf = m
18051     continue
          write(6,*)'nbuf is now ',nbuf
          do 1807 i=1,msel
            isels(i) = 0
            ifits(i) = 0
1807      continue
          nsel = 0
          write(6,*)'selections are removed'
         goto 2000
       endif
c
c
       if(comand.eq.'sel     ') then
c                    ---
         if(inames.eq.0 .or. found('add     ') 
     *                  .or. found('fit+    ')) then
           if(.not.found('add     ')) m = 0
           do 179 i=1,ipars
             iss      = rpar(i) * 1.0000001
             if(iss.gt.0) write(6,*)'select adress   ',iss
             m = m + 1
             isels(m) = iss
             ifits(m) = 0
             if(iss.lt.0.and.i.gt.1) then
               m   = m - 1
               ias = isels(m)
               ies = -iss
               write(6,*)'select adresses from ',ias,' to ',ies
               do 17911 l = ias,ies
                isels(m) = l
                ifits(m) = 0
                m = m + 1
17911          continue
               m = m - 1
             endif
179        continue
           nsel = m

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

           goto 2000
         endif
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
               nnumi(l) = rpar(l-1+j)*1.0000001
181          continue
             call search(nnumi,k)
           endif
           if(vname(i).eq.'sc+     ') then
              do 182 l=1,nsel
               nnumi(l) = numor(isels(l))
182           continue
             do 183 l=1,k
               kk = l+nsel
               if(kk.gt.minc) goto 183
                 nnumi(kk) = rpar(l-1+j)*1.0000001
                 iprs      = kk
183          continue
             call search(nnumi,iprs)
           endif
7077     continue



         goto 2000
       endif
c
       if(comand.eq.'putpar  ') then
c                    ------
         do i=1,nsel
            iaddp = isels(i)
            call parset (vname(1),sngl(rpar(1)),iaddp)
         enddo
         goto 2000
       endif
c
       if(comand.eq.'rename  ') then
c                    ------
         do j=1,inames
           do i=1,nsel
             iaddp = isels(i)
             if(vname(j).eq.'xaxis   ') xname(iaddp) = vname(j+1)
             if(vname(j).eq.'yaxis   ') yname(iaddp) = vname(j+1)
             if(vname(j).eq.'name    ')  name(iaddp) = vname(j+1)
           enddo
           do i=1,nfsel
             iaddp = isfits(i)
             if(vname(j).eq.'xaxis   ') xname(iaddp) = vname(j+1)
             if(vname(j).eq.'yaxis   ') yname(iaddp) = vname(j+1)
             if(vname(j).eq.'name    ')  name(iaddp) = vname(j+1)
           enddo
         enddo
         goto 2000
       endif
c
c
       if(comand.eq.'dispsel '.or.comand.eq.'dsl     ') then
c                    --------                ---
          write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
390       format(' selected items: ',i5/
     *           ' scan-address:  numor:      fit-address:'/
     *           (1x,i8,5x,i8,5x,i8))
         goto 2000
       endif
c
c
c
       if(comand.eq.'edit    ') then
c                    ----
         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         if(vname(1).eq.'sc      '.or.vname(1).eq.'numor   ') then
          do 499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 498
           endif
499       continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
498       continue
         endif
c
         if(vname(1).eq.'n       '.or.vname(1).eq.'name    ') then
          do 599 i=1,nbuf
           if(vname(2).eq.name(i)) then
             ispc = i
             goto 598
           endif
599       continue
          write(6,*)'data with name ',vname(2),' not found!'
          ierrs = 201
          goto 2000
598       continue
         endif
         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
c ----------> write data onto buffer <---------------------------------
         call savdat('datbuf  ',ispc)
c ----------> enter the system editor with that file <-----------------
         call system('dxnotepad datbuf')
c --- reread it to the same place ---
         nbuff = nbuf
         nbuf = ispc - 1
         inames = 3
         ipars = 0
         vname(1) = 'datbuf  '
         call input
         nbuf = nbuff
         goto 2000
       endif
c
c
c
       if(comand.eq.'save    ') then
c                    ----
         fsname = 'lastsave'
         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         do i=1,inames
           if(vname(i).eq.'to      ') fsname = vname(i+1)
         enddo
         do j=1,inames
         if(vname(1).eq.'sc      '.or.vname(j).eq.'numor   ') then
          do 4499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 4498
           endif
4499      continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
4498      continue
         endif
         enddo
c
         if(vname(1).eq.'n       '.or.vname(i).eq.'name    ') then
          do 4599 i=1,nbuf
           if(vname(2).eq.name(i)) then
             ispc = i
             goto 4598
           endif
4599      continue
          write(6,*)'data with name ',vname(2),' not found!'
          ierrs = 201
          goto 2000
4598      continue
         endif
         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
c ----------> write data onto buffer <---------------------------------
         if(fsname(1:1).eq.'&') then
           in = 0
           do j=2,8
            if(fsname(j:j).ne.' ') then
             in = in + 1
             name(ispc)(in:in) = fsname(j:j)
            endif
           enddo
           do j=1,8-in
             in = in+1
             name(ispc)(in:in) = coment(ispc)(j:j)
           enddo
           fsname = name(ispc)
         endif
         call savdat(fsname,ispc)
         goto 2000
       endif
c


       if(comand.eq.'msave   ') then
c                    -----
         fsname = 'lastsave'
         if(inames.gt.0)   fsname = vname(1)
         call msavdat(fsname)
         goto 2000
       endif
c
c
       if(comand.eq.'thc     ') then
c                    ---
         npoint = rpar(1) + 0.001
         call couple(1)
         call thc(npoint)
         goto 2000
       endif
c
c
       if(comand.eq.'fit     ') then
c                    ---
         call fit
         goto 2000
       endif
c
       if(comand.eq.'yfitform ') then
c                    -------
         if(ioldc.ne.0) then
           yfitform = reslin
           ioldc   = 0
         endif
         write(6,*)'y-fit-formel=',yfitform
         goto 2000
       endif
c
       if(comand.eq.'yformel  ') then
c                    -------
ccc      yformel = title
         if(ioldc.ne.0) then
           yformel  = reslin
           ioldc   = 0
         endif
         write(6,*)'y-formel=',yformel
         goto 2000
       endif
c
       if(comand.eq.'xformel  ') then
c                    -------
ccc      xformel = title
         if(ioldc.ne.0) then
           xformel  = reslin
           ioldc   = 0
         endif
         write(6,*)'x-formel=',xformel
         goto 2000
       endif
c
c
       if(comand.eq.'open    ') then
c                    ----
         npax = 1
         npay = 1
         fpnam = vname(1)
         xpnam = vname(2)
         if(inpar(2).ne.0) npax  = rpar(inapa(2)) + 0.0001
         ypnam = vname(3)
         if(inpar(3).ne.0) npay  = rpar(inapa(3)) + 0.0001
         open(33,file=fpnam,status='UNKNOWN')
         write(6,*)'file ',fpnam,' a has been opened for write'
         write(6,*)'data stored on write: ',xpnam,'(',npax,')  vs ',
     *              ypnam,'(',npay,')'
         write(33,*)title
         write(33,*)fpnam,'    ',ypnam,'  vs   ',xpnam,'     111111'
         write(33,*)
         write(33,*)'values'
         goto 2000
       endif
c
c
       if(comand.eq.'write   ') then
c                    -----
c        --- take it from the first selected file ---
         if(nsel.eq.0) then
           write(6,*)'no files selected !'
           goto 2000
         endif
         iadd = isels(1)
         iaddf= ifits(1)
c          --- look first for parameters ---
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
c
c
       if(comand.eq.'close   ') then
c                    -----
         write(33,*)
         write(33,'(4h#eod)')
         close(33)
         write(6,*)'file ',vname(1),' a has been closed'
         goto 2000
       endif
c
c
       if(comand.eq.'title   '.or.comand.eq.'tit     ') then
c                    -----
         title = inline(6:80)
         if(comand.eq.'tit     ') title = inline(4:80)
         goto 2000
       endif
c
       if(comand.eq.'theos   ') then
c                    -----> list available theories
         write(6,3001)(i,thenam(i),nthpar(i),
     *                   (thparn(j,i),j=1,mtpar),i=1,mth)
3001     format(' ***** theories available *****'/
     *   (1x,i3,': ',a8,i3,4(/10(1x,a8))))
         goto 2000
       endif
c
       if(comand.eq.'activate'.or.comand.eq.'ac      ') then
c                    -----> activate a theory
         call activa(0)
         call activa(2)
         goto 2000
       endif
c
      if(comand.eq.'label   ') then
c                   -----> assign a label to a parameter
         call lsearch(jpar,itcal,ierr)
         thpala(jpar,itcal) = vname(3)(1:4)
         goto 2000
      endif

      if(comand.eq.'chgthpar') then
c                   --------> change singel theory parameters 
         call lsearch(jpar,itcal,ierr)

         ithc  = nthtab(itcal)

         if(ierr .eq. 0) then
           thparx(jpar,itcal) = 
     *        getval('par     ',dble(thparx(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)') 
     *                'change parameter ',thparn(jpar,ithc),
     *                ' of th(',itcal,'):',thenam(ithc),' to ',
     *                thparx(jpar,itcal) 
           endif
           thpsca(jpar,itcal) = 
     *        getval('scale   ',dble(thpsca(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)') 
     *                'change scale of  ',thparn(jpar,ithc),
     *                ' of th(',itcal,'):',thenam(ithc),' to ',
     *                thpsca(jpar,itcal) 
           endif
         endif
         goto 2000
      endif
c
      if(comand.eq.'couple  ') then
c                   ------> couple label to a parameter
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
8011      continue
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
c
c
c
       if(comand.eq.'acl     '.or.comand.eq.'aclast  ') then
c                    -----> reactivate theories as stored in lastth
         call activa(3)
         if(iout.ge.0) call activa(2)
         goto 2000
       endif
c
       if(comand.eq.'desactiv'.or.comand.eq.'dac     ') then
c                    -----> desactivate all theories
         call activa(1)
         goto 2000
       endif
c
       if(comand.eq.'activlst'.or.comand.eq.'al      ') then
c                    -----> list activated theories
         call activa(2)
         goto 2000
       endif
c
c
       if(comand.eq.'plot    '.or.comand.eq.'p       ') then
c                    -----> plot selected curves
         call splot(.true.)
         ibild1 = ibild
         goto 2000
       endif
c
c
c
       if(comand.eq.'plot0   '.or.comand.eq.'p0      ') then
c                    -----> set parameters for plot
         call splot(.false.)
         ibild1 = ibild
         goto 2000
       endif
c
c
c
       if(comand.eq.'zero   ') then
c                    ----      ------> logical zero the common storage
         nbuf= 0
         goto 2000
       endif
c
c
       if (comand.eq.'numorpls') then
c                    -----------> change offset between files
         numpls = rpar(1) + 0.001
         goto 2000
       endif
c
c
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
              write(6,*)'theory(',ith,'): ',thenam(ith),
     *        '  does not belong to the cail-family'
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
           enddo !! isel
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
 
           enddo !! iz

         enddo !! iortho

20009    continue

         call theo_out(85)
         close(85)
         write(6,*)'Ready, Data written to: ',filnam

         goto 2000
       endif
c
c
       call makro
c
       goto 2000
c      ---------> read now commandlines from makro file
c
c
       end


c*ds
       subroutine splot (doplo)
c      ================  scan-plotting
c
c
       parameter(mth=40,mtpar=40,mtcal=40)
c      parameter(mkurv=minc)
c                ---------> max. no of curves to be plotted
c ---  maximum scan length
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c ---- common containing all the scans ----
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
       common/fslist/isfits(mbuf),nfsel
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
       common/therrc/therro(mtpar,mtcal)
c
c
       character*6 cnum
c
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*8 comand,vname
       character*20 arglst,pmlist
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       dimension x(mwert),y(mwert),yl(mwert),yh(mwert),irecv(minc),
     *          isymb(minc),irecn(minc),icolo(minc),ifrec(minc)
       dimension e(mwert)
 
       logical found 
       logical ptex
       logical paplo
       logical doplo
       logical fitplo
       logical errplo
       logical paxis
       logical taxis

       logical log_x
       logical log_y
       real    log10 

c      --- doplo = false  means: set parameters only ---
       character*80 option,xtext,ytext,tbuf
       character*8  dirnam(4),opart(8)

       character*12 tag, stunde
       character*12 tx,sx
c
      data xmin/0./,xmax/4./,ymin/0./,ymax/300./,nkurv/0/
      data isymb/4,5,23,6,16,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21
     *          ,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41
     *          ,42,43/
      data ixmode/1/,txsize/.45/,sysize/.1/
      data framx/23./,framy/23./,yensca/0./,frlux/4./,frluy/4./
      data dirnam/'qx      ','qy      ','qz      ','en / thz'/
      data epscl/0.001/
      data ptex/.true./
      data icolo/minc * 0/,ifont/1/
      data icol0/1/
      data paplo/.true./
      data paxis/.true./
      data taxis/.true./
      data fitplo/.true./
      data errplo/.false./
      data log_x/.false./
      data log_y/.false./
      data opart/'x=1     ','y=1     ','i       ','a       ',
     *           'f=(3,1) ','m=2.0   ','u=1     ','        '/
      data txsizt/.23/,xtshft/0./,ytshft/0./
c
c
c
c ----- parameter retrieving from stack -----
      nkurv  = 0
      nfkurv = 0
      if(inames.eq.0) then
c                     ----> assume that a list of number numors
       if(ipars.gt.0) then
          nkurv = 0
          nsel  = 0
          if(ipars.gt.minc) ipars = minc
          do 2 i=1,ipars
           irecv(i) = rpar(i) + 0.0001
2         continue
          nkurv = ipars
        endif
       else
c      -----> decode by names
        do 3 i=1,inames
          j = inapa(i)
          if(vname(i).eq.'xmin    ') xmin = rpar(j)
          if(vname(i).eq.'xmax    ') xmax = rpar(j)
          if(vname(i).eq.'ymin    ') ymin = rpar(j)
          if(vname(i).eq.'ymax    ') ymax = rpar(j)
          if(vname(i).eq.'framx   ') framx= rpar(j)
          if(vname(i).eq.'framy   ') framy= rpar(j)
          if(vname(i).eq.'frlux   ') frlux= rpar(j)
          if(vname(i).eq.'frluy   ') frluy= rpar(j)
          if(vname(i).eq.'text    ') ptex = .true.
          if(vname(i).eq.'txon    ') ptex = .true.
          if(vname(i).eq.'notext  ') ptex = .false.
c   ---> txon, txoff are from old version and not documented any more <-
          if(vname(i).eq.'txoff   ') ptex  = .false.
          if(vname(i).eq.'epscl   ') epscl = rpar(j)
          if(vname(i).eq.'txsize  ') txsize= rpar(j)
          if(vname(i).eq.'legsize ') txsizt= rpar(j)
          if(vname(i).eq.'legx    ') xtshft= rpar(j)
          if(vname(i).eq.'legy    ') ytshft= rpar(j)
          if(vname(i).eq.'font    ') ifont = rpar(j) +0.0001
          if(vname(i).eq.'color   ') icol0 = rpar(j) +0.0001
          if(vname(i).eq.'sysize  ') sysize= rpar(j)
          if(vname(i).eq.'parplo  ') paplo=.true.
          if(vname(i).eq.'noparplo') paplo=.false.
          if(vname(i).eq.'errplo  ') errplo=.true.
          if(vname(i).eq.'noerrplo') errplo=.false.
          if(vname(i).eq.'axis    ') paxis =.true.
          if(vname(i).eq.'noaxis  ') paxis =.false.
          if(vname(i).eq.'txaxis  ') taxis =.true.
          if(vname(i).eq.'notxaxis') taxis =.false.
          if(vname(i).eq.'ox      ') opart(1) = vname(i+1)
          if(vname(i).eq.'oy      ') opart(2) = vname(i+1)
          if(vname(i).eq.'oi      ') opart(3) = 'a'
          if(vname(i).eq.'oa      ') opart(3) = 'i'
          if(vname(i).eq.'of      ') opart(5) = vname(i+1)
          if(vname(i).eq.'om      ') opart(6) = vname(i+1)
          if(vname(i).eq.'ou      ') opart(7) = vname(i+1)
          if(vname(i).eq.'fits    ') fitplo = .true.
          if(vname(i).eq.'log_x   ') log_x  = .true.
          if(vname(i).eq.'log_y   ') log_y  = .true.
          if(vname(i).eq.'lin_x   ') log_x  = .false.
          if(vname(i).eq.'lin_y   ') log_y  = .false.
c
          if(vname(i).eq.'symb    ') then
            nsy = 0
            do 49 l=1,inpar(i)
             nsy   = nsy   + 1
             if(nsy.gt.minc) goto 29
             isymb(nsy) = rpar(j) + 0.0001
             j = j + 1
49          continue
          endif
29        continue
c
          if(vname(i).eq.'icolo   ') then
            nco = 0
            do 59 l=1,inpar(i)
             nco   = nco   + 1
             if(nco.gt.minc) goto 39
             icolo(nco) = rpar(j) + 0.0001
             j = j + 1
59          continue
          endif
39        continue
c
          if(vname(i).eq.'sc      ') then
            nsel = 0
            nkurv= 0
            do 5 l=1,inpar(i)
             nkurv = nkurv + 1
             if(nkurv.gt.minc) goto 31
             irecv(nkurv) = rpar(j) * 1.000001
             j = j + 1
5           continue
          endif
31        continue
c
          if(vname(i).eq.'fsc     ') then
            nfsel = 0
            nfkurv= 0
            do 641 l=1,inpar(i)
             nfkurv = nfkurv + 1
             if(nfkurv.gt.minc) goto 64
             ifrec(nfkurv) = rpar(j) * 1.000001
             j = j + 1
641         continue
          endif
64        continue
          if(nfsel.eq.0) then
            if(nfkurv.gt.0) then
               call fsrch(ifrec,nfkurv)
            endif
          endif
c
3       continue
      endif
c
      if (.not.doplo) then
        return
      endif
c
      if(nsel.eq.0) then
        if(nkurv.gt.0) then
          call search(irecv,nkurv)
c       -----------------------------> select the items to be plotted
        else
          write(6,*)'no scan selected, no plot possible !'
          return
        endif
      endif
c
      if (.not.doplo) then
        return
      endif
c ---->                                             <-----------
      nkurv = nsel
c ----> ich denke, das muss so sein !?!?!?!  (m.s.) <-----------
      if (nkurv.eq.0) then
        write(6,*)'there is no curve selected => no plot !'
        return
      endif
c
c
c      write(*,*)' entering plot section .....'
c
c ----- initialize plot (if first) -----
       if(ibild.eq.0) then
         call grstrt(35,8)
       endif
c
       ibild = ibild + 1
c
c ---- set frame & scales ----
c      write(*,*)' set frames & scales .....'
       frxx = frlux + framx
       fryy = frluy + framy
       call grsclc(frlux,frluy,frxx,fryy)
       call grsclv(xmin,ymin,xmax,ymax)
       call grfont(ifont)
       call grnwpn(icol0)
       call grchrc(txsize,0.,16 )     !! aix
c      --------------------------> set size of signs
c ---- prepare axes ----
c
c ---- identify scantyp ----
c
       if(log_x) then
         lxx   = laenge(xname(isels(1)),80,' ')
         xtext = 'log10('//xname(isels(1))(1:lxx)//')'
         xmi_s = log10(abs(xmin)+1e-33)
         xma_s = log10(abs(xmax)+1e-33)
         if(xma_s .lt. xmi_s) then
           xh    = xma_s
           xma_s = xmi_s
           xmi_s = xh
         endif
       else
         xtext = xname(isels(1))
         xmi_s = xmin
         xma_s = xmax
       endif
       if(log_y) then
         lyy   = laenge(yname(isels(1)),80,' ')
         ytext = 'log10('//yname(isels(1))(1:lyy)//')'
         ymi_s = log10(abs(ymin)+1e-33)
         yma_s = log10(abs(ymax)+1e-33)
         if(yma_s .lt. ymi_s) then
           xh    = yma_s
           yma_s = ymi_s
           ymi_s = xh
         endif
       else
         ytext = yname(isels(1))
         ymi_s = ymin
         yma_s = ymax
       endif
       call grsclv(xmi_s,ymi_s,xma_s,yma_s)

       ltext = 25
       option = opart(1)//','//opart(2)//','//opart(3)//','//
     *                         opart(5)//','//opart(6)//','//opart(7)
       call upcase(option,80)
       write(6,*)option
       lopt = 80
c      write(*,*)' make axes .....'
       if(.not.taxis) then 
         lxx = 0
         lyy = 0
       else
         lxx = laenge(xtext,80,' ')
         lyy = laenge(ytext,80,' ')
       endif 
       if(paxis) call graxs(lopt,option,lxx,xtext,lyy,ytext)
c
c
c ---- plot the selected fit-kurves ----
c
       nfkurv = nfsel
       if(fitplo) then
       do 70 i=1,nfkurv
        irfcu = isfits(i)
        npicf = nwert(irfcu)
        nnpi = 0
        do 70010 j=1,npicf
          if(xwerte(j,irfcu).ge.xmin.and.xwerte(j,irfcu).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,irfcu)
!!          if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).lt.ymin) y(nnpi) = ymin
            if(y(nnpi).gt.ymax) y(nnpi) = ymax
            x(nnpi) = xwerte(j,irfcu)
            
            if(log_x) then
              if(x(nnpi).ne.0.0) then
                x(nnpi) = log10(abs(x(nnpi)))
              else
                x(nnpi) = -33.
              endif
            endif
            if(log_y) then
              if(y(nnpi).ne.0.0) then
                y(nnpi) = log10(abs(y(nnpi)))
              else
                y(nnpi) = -33.
              endif
            endif

          endif
70010     continue
          call grln(x,y,nnpi)
70      continue
       endif
c
c
c ---- plot datarecords ----
c
       nkurv = nsel
       do 20 i=1,nkurv
        ircu = isels(i)
        ircf = ifits(i)
        npic = nwert(ircu)
        if(fitplo) then
        if(ircf.ne.0) then
c                     ----> plot the fitted data automatically
c                           this is after a fit-command until you
c                           select new curves
        npicf = nwert(ircf)
        nnpi  = 0
        do 20010 j=1,npicf
          if(xwerte(j,ircf).ge.xmin.and.xwerte(j,ircf).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,ircf)
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).gt.ymax) y(nnpi) = ymax
            x(nnpi) = xwerte(j,ircf)

            if(log_x) then
              if(x(nnpi).ne.0.0) then
                x(nnpi) = log10(abs(x(nnpi)))
              else
                x(nnpi) = -33.
              endif
            endif
            if(log_y) then
              if(y(nnpi).ne.0.0) then
                y(nnpi) = log10(abs(y(nnpi)))
              else
                y(nnpi) = -33.
              endif
            endif

          endif
20010     continue
          call grln(x,y,nnpi)
        endif
       endif
c
c                  ---> prepare data to plot
        nnpi = 0
        do 2001 j=1,npic
          if(xwerte(j,ircu).ge.xmin.and.xwerte(j,ircu).le.xmax) then
            nnpi = nnpi + 1
            y(nnpi) = ywerte(j,ircu)
            x(nnpi) = xwerte(j,ircu)
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02
            if(y(nnpi).gt.ymax) y(nnpi) = ymax
            e(nnpi) = yerror(j,ircu)
            
            if(log_x) then
              if(x(nnpi).ne.0.0) then
                x(nnpi) = log10(abs(x(nnpi)))
              else
                x(nnpi) = -33.
              endif
            endif
            if(log_y) then
              if(y(nnpi).ne.0.0) then
                y(nnpi) = log10(abs(y(nnpi)))
              else
                y(nnpi) = -33.
              endif
            endif

          endif
2001    continue
c
c
c
c --- plot ---
c       if (numor(ircu).gt.0) then
           icco=mod(icolo(i),7) + 1
 
           call grnwpn(icco)
c          ----------- plot a dataline -------
           if(isymb(i).eq.0) then
             call grln(x,y,nnpi)
           else
             sysiz = sysize
             call grchrc(sysiz,0.,16)
             do 2003 ik=1,nnpi
              call grjmps(x(ik),y(ik),isymb(i))
2003         continue
             call grchrc(txsize,0.,16)
           endif
 
           if(errplo) then
             do ik=1,nnpi
              if(log_y) then
               yepl = log10(10.0**(y(ik))+e(ik))
               yeml = 10.0**(y(ik))-e(ik)
               if(yeml.gt.0) then
                 yeml = log10(yeml) 
               else
                 yeml = -33.0
               endif
              else
               yepl = y(ik)+e(ik)
               yeml = y(ik)-e(ik)
              endif
!              call grjmp(x(ik),y(ik)-e(ik))
!              call grdrw(x(ik),y(ik)+e(ik))
              call grjmp(x(ik),yeml)
              call grdrw(x(ik),yepl)
             enddo
           endif
 
             call grnwpn(1)
c       endif
c
20     continue
c
c

c ---- textpart ----
       if(ptex) then
c      write(*,*)' enter textplotting....'
       call grsclv(xmin,ymin,xmax,ymax)
c
c --- title ---
         xtx = xmin + 0.1 * (xmax-xmin)
         ytx = ymax
         ltext = 20./txsize
         if(ltext.gt.74) ltext=74
         call grtxt(xtx,ytx,ltext,title)
c ---- set textwindow ----
c
         call grsclc(0.,0.,39.5,28.7)
         call grsclv(0.,0.,39.5,28.7)
c        txsizt = 0.65 * txsize
         call grchrc(txsizt,0.,16 )   !! aix
         xtx = 28.3 + xtshft
         ytx = 28.  + ytshft
!
! - plot date and time info:
!

       call date_and_time(tag,stunde)
  
       tx = tag(7:8)//'-'//tag(5:6)//'-'//tag(1:4)
       sx = stunde(1:2)//':'//stunde(3:4)//':'//stunde(5:6)
       xtext = tx//'  '//sx
       call grtxt(0.5,0.5,26,xtext)

c ---- plot theory parameters ----
         if(ntheos.ne.0) then
           do 115 it = 1,ntheos
             ith = nthtab(it)
             if(multflg(it).eq.1) then
                write(xtext,'(8htheory* ,a8)')thenam(ith)
             else
                write(xtext,'(8htheory+ ,a8)')thenam(ith)
             endif
             call grtxt(xtx,ytx,16,xtext)
             ytx = ytx - 1.7 * txsizt
             npar = nthpar(ith)
             if(npar.ne.0) then
               do 117 ip = 1,npar
               write(xtext,'(a8,1h=,1e12.4,2h+-,e9.2,e8.1)')
     *         thparn(ip,ith),thparx(ip,it),therro(ip,it),thpsca(ip,it)
                 call grtxt(xtx,ytx,33,xtext)
                 ytx = ytx - 1.7*txsizt
 117           continue
            endif
 115      continue
         endif
c ---- plotted items ----
         do 101 i=1,nkurv
           ircu = isels(i)
           write(xtext,'(a8,i14)') name(ircu),numor(ircu)
           xtxs = xtx - 2*txsizt
           ytxs = ytx + txsizt / 2
           icco=mod(icolo(i),7) + 1
           call grnwpn(icco)
           if(isymb(i).ne.0.and.numor(ircu).gt.0) then
             call grjmps(xtxs,ytxs,isymb(i))
           else
             call grtxt(xtxs,ytxs,1,'-')
           endif
           call grtxt(xtx,ytx,22,xtext)
           ytx = ytx - 1.7 * txsizt
           call grtxt(xtx,ytx,80,coment(ircu))
           ytx = ytx - 1.7 * txsizt
           if(paplo) then
           do 1012 l=1,nopar(ircu)
           write(xtext,'(a8,e14.5)')napar(l,ircu),params(l,ircu)
           call grtxt(xtx,ytx,22,xtext)
           ytx = ytx - 1.7 * txsizt
1012       continue
           else
           do l=1,nopar(ircu)
            if(found(napar(l,ircu)//' ')) then
              write(xtext,'(a8,e14.5)')napar(l,ircu),
     *                                 params(l,ircu)
              call grtxt(xtx,ytx,22,xtext)
              ytx = ytx - 1.7 * txsizt
            endif
           enddo
           endif

           ircu = ifits(i)
           if(ircu.gt.0) then
             write(xtext,'(a8,i14)') name(ircu),numor(ircu)
             xtxs = xtx - 2*txsizt
             ytxs = ytx + txsizt / 2
             icco=mod(icolo(i),7) + 1
             call grnwpn(icco)
             if(isymb(i).ne.0.and.numor(ircu).gt.0) then
               call grjmps(xtxs,ytxs,isymb(i))
             else
               call grtxt(xtxs,ytxs,1,'-')
             endif
             call grtxt(xtx,ytx,22,xtext)
             ytx = ytx - 1.7 * txsizt
             call grtxt(xtx,ytx,80,coment(ircu))
             ytx = ytx - 1.7 * txsizt
             if(paplo) then
             do  l=1,nopar(ircu)
             write(xtext,'(a8,e14.5)')napar(l,ircu),params(l,ircu)
             call grtxt(xtx,ytx,22,xtext)
             ytx = ytx - 1.7 * txsizt
             enddo
             else
             do l=1,nopar(ircu)
              if(found(napar(l,ircu)//' ')) then
                write(xtext,'(a8,e14.5)')napar(l,ircu),
     *                                   params(l,ircu)
                call grtxt(xtx,ytx,22,xtext)
                ytx = ytx - 1.7 * txsizt
              endif
             enddo
             endif
           endif
101      continue
         call grnwpn(1)
      endif
c
c     write(*,*)' ok ready call nextframe ....'
      call grnxtf
c
      return
      end
c
c
c
c
c*ds
c*ed
       subroutine lsearch( jpar, itcal, ierr)
c      ======================================
c
c  ---- search for a theory parameter specification by the commandline
c       cmd  theoryname <n-th occ> parametername
c  output: jpar = parameter-adress
c          itcal= theory adress
c          ierr = errorindicator ( 0=ok  1=not found)
c
c
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)
c
c
 
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
 
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
 
!!     common/outlev/iout,ibild,ierrs,inka   !! aix
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c

       character*8      vnameF
       double precision rparF
       integer          inpaF

c ----------------------------------------------------------------------
c --- which occurence of the theory
       ict = 0
       if(inpaF(1).eq.0) then
          icz = 1
       else
          icz = rparF(1) + 0.1
       endif
c
       do 1 i=1,ntheos
          ith = nthtab(i)
          if(vnameF(1).eq.thenam(ith)) then
            ict = ict + 1
            if(ict.eq.icz) goto 2
c                          ------> ok we have it
          endif
1      continue
          write(6,*)'search: theory=',vnameF(1),' could not be found'
          ierr = 1
          ierrs = 801
          return
2      continue
       itcal= i
c      ========
c ---- ok theory is identified. look now for parameter.....
       np = nthpar(ith)
       do 4 j=1,np
          if(vnameF(2).eq.thparn(j,ith)) goto 6
4      continue
         write(6,*)'search: parameter=',vnameF(2),
     *             ' of theory=',vnameF(1),
     *             ' could not be found'
         ierr = 1
         ierrs = 802
         return
6      continue
       ierr = 0
       jpar = j
c      ========
c
       return
       end
c*ds
c*ds
      subroutine search(irecv,nkurv)
c     ------------------------------------
c
c --- the adresses(/spectr/) of the nkurv  numors given in irecv are
c     located and put to irecn1
c     nkurv is set to the number of actually found items
c  !  selection table is updated
c     imx returns the no. of biggest gradient contribution
c
c
c ----- search for the scans to be selected -----
       parameter(minc=40)
       parameter(mth=40,mtpar=40,mtcal=40)
c      parameter(mkurv=minc)
c                ---------> max. no of curves to be selected
c ---  maximum scan length
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c ---- common containing all the scans ----
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
       common/fslist/isfits(mbuf),nfsel
c
      common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
      dimension irecv(minc),irecn1(minc)
c
      l = 0
      do 1001 ic = 1,nkurv
        num = irecv(ic)
        if(iout.gt.3) write(6,*)'search for scan: ',num
          do 1002 i=1,nbuf
            if(num.eq.numor(i)) then
               l = l+1
               irecn1(l) = i
               goto 1001
            endif
1002      continue
        if(iout.ge.0) write(6,*)'scan: ',num,' has not been found!'
1001   continue
c
       nkurv = l
c      ---------> only found scans may be selected ...
       if(nkurv.le.0) then
        if(iout.ge.0) write(6,*)'no scan has been found'
        return
       endif
c
c
c --- check for the other scans ---
       if(nkurv.ge.2) then
         inum1 = irecn1(1)
         do 1020 i=2,nkurv
           inum2 = irecn1(i)
           if(xname(inum2).ne.xname(inum1)) then
             write(6,*)'warning abscissas may be incompatible:'
             write(6,*)name(inum1),' x: ',xname(inum1),
     *                 name(inum2),' x: ',xname(inum2)
           endif
1020     continue
       endif
c      ------(of if(nkurv.ge.2)....)
c
c --- update selection table ---
       do 2000 i=1,nkurv
        isels(i) = irecn1(i)
        ifits(i) = 0
2000   continue
       nsel = nkurv
       if(iout.ge.0) write(6,*)'selection table updated. ',l,' items'
c
       return
       end
c
c
c*ds
c*ds
      subroutine fsrch(ifrec,nfkurv)
c     ------------------------------------
c this subroutine is analogous to search, it will select the scans of
c the fit curves
c attention: the fsc-command requires positive inputs, therefore this
c            routine searches for -num !!
c
       parameter(minc=40)
       parameter(mth=40,mtpar=40,mtcal=40)
c      parameter(mkurv=minc)
c                ---------> max. no of curves to be selected
c ---  maximum scan length
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c ---- common containing all the scans ----
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
      common/fslist/isfits(mbuf),nfsel
c
      common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
      dimension ifrec(minc),irecn1(minc)
c
      l = 0
      do 1001 ic = 1,nfkurv
        num = -ifrec(ic)
        if(iout.gt.3) write(6,*)'search for fscan: ',num
          do 1002 i=1,nbuf
            if(num.eq.numor(i)) then
               l = l+1
               irecn1(l) = i
               goto 1001
            endif
1002      continue
        if(iout.ge.0) write(6,*)'scan: ',num,' has not been found!'
1001   continue
c
       nkurv = l
c      ---------> only found scans may be selected ...
       if(nfkurv.le.0) then
        if(iout.ge.0) write(6,*)'no scan has been found'
        return
       endif
c
c
c --- check for the other scans ---
       if(nfkurv.ge.2) then
         inum1 = irecn1(1)
         do 1020 i=2,nfkurv
           inum2 = irecn1(i)
           if(xname(inum2).ne.xname(inum1)) then
             write(6,*)'warning abscissas may be incompatible:'
             write(6,*)name(inum1),' x: ',xname(inum1),
     *                 name(inum2),' x: ',xname(inum2)
           endif
1020     continue
       endif
c      ------(of if(nfkurv.ge.2)....)
c
c --- update f-selection table ---
       do 2000 i=1,nfkurv
        isfits(i) = irecn1(i)
2000   continue
       nfsel = nfkurv
       if(iout.ge.0) write(6,*)'selection table updated. ',l,' items'
c
       return
       end
c
c
c*ds
c*ds
       subroutine fit
c      ==============
c
c ---- fit of activated theories to data on /fil2/ ----
c
       parameter (mfit=40,msmpl=4000)
       parameter(mwert=1024,mbuf=200,mpar=200)
c      parameter(mkurv=minc)
c  mfit = max-no. of fitvariables
c  msmpl= max-no. of data points where a comparison is made
       parameter(mth=40,mtpar=40,mtcal=40)
       character*8 ci
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar, getval
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
 
       common/therrc/therro(mtpar,mtcal)
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
      dimension iparam(6),rparam(7),x(mfit),f(msmpl),xjac(msmpl,mfit),
     *          xguess(mfit),xscale(mfit),fscale(msmpl)
      dimension numv(minc),numn(minc)
      dimension ermat(mfit,msmpl), gmat(mfit,mfit), ginv(mfit,mfit)
c
      dimension   xcenter(mfit), xstepsc(mfit)
      character*8 xpname(mfit)  
c
      external func
c
      logical sqwght,sqwbuf
      logical autox1,autox2
      logical found, folgt
      logical final_thc, lbuffer
      common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)

       real    fcssq
       logical lerrel, lwrtfitdat
       common/cfunce/ lerrel, lwrtfitdat, fcssq
c ---- communication with subr. func ---

       real*8 ssq 
c
c
c ---- defaults for parameters -----
       data nsig/6/,eps/1.d-5/,delta/1.d-5/,maxfn/50/,iopt/0/
       data ngood/0/,maxit/0/,stpsz/0.0/,trure/0.0/
       data ifalt/0/,iprint/1/,irecse/0/
       data sqwbuf/.false./,x1i/-1.e30/,x2i/1.e30/
       data final_thc/.false./
c
c ---- take parameters from stack ----
       sqwght= sqwbuf
      igo = 0
      if(inames.eq.0) goto 1000
      do 1011 i=1,inames
       j = inapa(i)
       ci= vname(i)
       if(ci.eq.'go      ') then
         igo = 1
         goto 1011
       endif
       if(ci.eq.'wlin    '.or.ci.eq.'linweigh') then
         sqwght= .false.
         goto 1011
       endif
       if(ci.eq.'thc     ') then
         final_thc= .true.
         goto 1011
       endif
       if(ci.eq.'nothc   ') then
         final_thc= .false.
         goto 1011
       endif
ccc    if(j.eq.0) then
ccc      write(6,*)' numerical parameter expected & no found ...'
ccc      ierrs = 700
ccc      return
ccc    endif
       if(vname(i).eq.'scans   '.or.vname(i).eq.'sc      ') then
         nspf = inpar(i)
         do 7001 l=1,nspf
          numv(l) = rpar(j+l-1) * 1.0000001
7001     continue
       endif
       if(ci.eq.'maxfn   '                    ) maxfn = rpar(j) + 0.001
       if(ci.eq.'maxit   '                    ) maxit = rpar(j) + 0.001
       if(ci.eq.'ngood   '                    ) ngood = rpar(j) + 0.001
       if(ci.eq.'maxstep '                    ) stpsz = rpar(j)
       if(ci.eq.'trustreg'                    ) trure = rpar(j)
       if(ci.eq.'relerr  '                    ) lerrel= .true.
       if(ci.eq.'abserr  '                    ) lerrel= .false.
       if(ci.eq.'wrtfit  '                    ) lwrtfitdat = .true.
       if(ci.eq.'nowrtfit'                    ) lwrtfitdat = .false.
c ----> the parameters  x1 x2 resp. the option autox1 autox2
c       are to be decoded during the first call of thc by thc !!
1011  continue
       sqwbuf= sqwght
c --- search scans & prepare selection table
       if(nsel.eq.0) call search(numv,nspf)
c      -----------------------------------------
c
c ---- write all parameters ----
      if(iout.ge.0) then
        write(6,1)iprint,maxfn,maxit,ngood,
     *           (numor(isels(i)),i=1,nsel)
1       format(' ***** fit-profile in effect ***** '/
     *  ' printing    ',i8/
     *  ' max funcalls',i8,/
     *  ' max iterat .',i8,/
     *  ' n-good-dig .',i8,/
     *  ' scans :'/(1x,10i7))
        if(sqwght)write(6,*)' ***** weight differences by 1/sqrt *****'
        if(lerrel)write(6,*)' ***** errors weigthed relative:relerr **'
        if(lwrtfitdat)write(6,*)' ***** writing fitdat.tmp  **'
        if(final_thc) then
          write(6,*)' ***** opt: thc   ->    final fine mesh call**'
        else
          write(6,*)' ***** opt: nothc -> no final fine mesh call**'
        endif
       endif
c
      if(igo.eq.0) return
1000  continue
c
          iprt = iprint
c
          write(6,*)' startparameters : '
          call activa(2)
c
c ---- prepare startvalues -----
          if(ntheos.eq.0) then
            write(6,*)' fit ==> no theories activated ...'
            return
          endif
c
          nfit = 0
          do 10 it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do 20 ip=1,npar
             if(thpsca(ip,it).ne.0) then
               nfit = nfit + 1
               if(nfit.gt.mfit) then
                 write(6,*)' too many fitparameters for current mfit=',
     *                     mfit
                 return
               endif
               x(nfit) = thparx(ip,it) / thpsca(ip,it)
               xstepsc(nfit) =  thpsca(ip,it)
               xpname(nfit)  =  thparn(ip,ith)
             endif
20        continue
10      continue
c
        n     = nfit
        ixjac = msmpl
        call func(m,n,x,f)
c       -------------------> determine m & di a first test calculation
        write(6,*)' no. of sample points m = ',m
        if(.not.autox1) write(6,*)'set lower limit of x: x1=',x1
        if(.not.autox2) write(6,*)'set upper limit of x: x2=',x2
c
c ------ imsl version 10 :  setup of some new vectors ------------------
c        to meet the function of the old zxssq approximately
        do 10010 i=1,n
           xguess(i) = x(i)
           xscale(i) = 1.0
10010   continue
c
        do 10020 i=1,m
           fscale(i) = 1.0
10020   continue
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(found('map     ')) goto 20000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c --- set the desired execution-parameters ---
        call u4lsf(iparam,rparam)
        if(iout.gt.0) write(6,*)'initial iparam=',iparam
        if(iout.gt.0) write(6,*)'initial rparam=',rparam
        if(ngood.ne.0) iparam(2) = ngood
        if(maxit.ne.0) iparam(3) = maxit
        if(maxfn.ne.0) iparam(4) = maxfn
        if(stpsz.ne.0) rparam(6) = stpsz
        if(trure.ne.0) rparam(7) = trure
        if(iout.gt.0) write(6,*)'used    iparam=',iparam
        if(iout.gt.0) write(6,*)'used    rparam=',rparam
c ---------------------------------------------------------------------
c ---- do the fitting with imsl routine unlsf (old:zxssq) ----
        call unlsf(func,m,n,xguess,xscale,fscale,iparam,rparam,x,f,xjac,
     *             ixjac)
        if(iout.gt.0) write(6,*)'final   iparam=',iparam
        if(iout.gt.0) write(6,*)'final   rparam=',rparam
c
        icode = iercd()
c       ---------------> get the errorcode
c ----- output -----
        write(6,*)'icode = ',icode
        if(icode.eq.1) write(6,*)' **-1- : no longer converging   *****'
        if(icode.eq.2) write(6,*)' **-2- : conv. to a noncritical pnt *'
        if(icode.eq.3) write(6,*)' **-3- : maxfn exceeded *************'
        if(icode.eq.4) write(6,*)' **-4- : maxit exceeded *************'
        if(icode.eq.5) write(6,*)' **-5- : 5 steps with max steplength*'
c
        write(6,501)iparam(4),iparam(3),iparam(2)
501     format(
     *   ' no. of func calls ..... ',i6/
     *   ' est. no. of sig. dig. . ',i6/
     *   ' no of iterations ...... ',i6/)
c
c
 
c ----- error-determination ?? test ?? -----------------------------
c
      if(.not.found('noerrors ')) then
         write(6,*)'error determination...'
         do i=1,n
          do j=1,n
           sum = 0
           do l=1,m
            sum = sum + xjac(l,i)*xjac(l,j)
           enddo
           gmat(i,j) = sum
          enddo
         enddo
c --- invertiere gmat ---
         call linrg(n,gmat,mfit,ginv,mfit)
c --- bilde g**-1 f  (f=xjac) ---
         do i=1,m
          do l=1,n
           sum = 0
           do j=1,n
            sum = sum + ginv(l,j)*xjac(i,j)
           enddo
           ermat(l,i) = sum
          enddo
         enddo
c --- bilde fehlerquadratmatrix ---
         do l=1,n
          do ll=1,n
           sum = 0
           do i=1,m
!!linwght!! sum = sum + ermat(l,i)*ermat(ll,i)*(ferror(i)**2)
            sum = sum + ermat(l,i)*ermat(ll,i)
           enddo
           gmat(l,ll) = sum
          enddo
         enddo
c --- extrahiere fehler und ausgabe ---
          nff = 0
          do    it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do    ip=1,npar
             therro(ip,it) = 0.0
             if(thpsca(ip,it).ne.0) then
               nff  = nff  + 1
               xerr = sqrt(gmat(nff,nff))
               thperr        = xerr    * thpsca(ip,it)
               write(6,*)'error for ', thenam(ith),'(',thparn(ip,ith),
     *                   ') = +- ',thperr
               therro(ip,it) = thperr
             endif
          enddo
        enddo
      else
          do    it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do    ip=1,npar
             therro(ip,it) = 0.0
            enddo
          enddo
      endif
c ------error-determination ?? end  ?? -----------------------------
 
        call func(m,n,x,f)
c       ------------------> final call including elastic lines
        write(6,*)' ******* final parameters ********'
        if(irecse.ne.0) write(6,*)'       record ',iru
        call activa(2)

        call extract('ssq0    ' ,ssq,ier)
        if(ier.ne.0) ssq = 0
        do i=1,nsel
          call parset ('ssq     ',sngl(ssq),ifits(i))
        enddo           
c
        call couple(1)
        if(final_thc) call thc(mwert)
c       ----------------> fitted curve with finest mesh
c
!!??    if(ier.ne.0) then
!!??      ierrs = ier
!!??      return
!!??    endif
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  create a map                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
20000 continue  
      if(found('map     ')) then
        map   = NINT(getval('map     ',0.0d0     ,inew))
        divis =      getval('div     ',20.0d0    ,inew)
     
        if(map.le.0) then
          write(6,*)'no map created, give no of points!'
          return
        endif

! check parameters:
        if(n.gt.3 .or. n.lt.0) then
          write(6,*)'fit map: the number of parameters=',n
          write(6,*)'is out of range 1..3'
        endif

        do i=1,n
         xcenter(i) = x(i)
        enddo
        

! ------> open-file <-------
        open(14,file='map.xy' ,status='unknown')

        write(6,*)'start mapping in ',n,' dimensions ...'

        write(14,'(a,i4,a)')'MAPPING in ',n,' dimensions'
        write(14,'(a)')     'Data selected: '
        do i=1,nsel
          j = isels(i)
          write(14,'(a,1x,i9,a,a)')name(j)(1:12),
     *                   numor(j),' :',coment(j)(1:60)
        enddo

        write(14,'(3a20)')(xpname(i),i=1,n)
        do j=-map,map
          do i=1,n
           x(i) = (xcenter(i) + j/divis)*xstepsc(i)
          enddo
          write(14,'(3E15.6)')(x(i),i=1,n)
        enddo

        close(14)

        open(15,file='map.ssq',status='unknown')

        write(15,'(a,2i8)')'ssq | map,n : ',map,n
        write(15,'(a)') name(isels(1))(1:20)
        do i=1,n
         write(15,'(a)') xpname(i)
        enddo
        if       (n.eq.1) then
          do j=-map,map
            x(1) = xcenter(1) + j/divis
            call func(m,n,x,f)
            write(15,'(E15.6)') fcssq
          enddo
        elseif   (n.eq.2) then
          do j1=-map,map
           x(1) = xcenter(1) + j1/divis
           do j2=-map,map
            x(2) = xcenter(2) + j2/divis
            call func(m,n,x,f)
            write(15,'(E15.6)') fcssq
           enddo
          enddo 
        elseif   (n.eq.3) then
          do j1=-map,map
           x(1) = xcenter(1) + j1/divis
           do j2=-map,map
            x(2) = xcenter(2) + j2/divis
             do j3=-map,map
               x(3) = xcenter(3) + j3/divis
               call func(m,n,x,f)
               write(15,'(E15.6)') fcssq
             enddo
           enddo
          enddo 
         
        endif    

        close(15)

! restore old setting !
        do i=1,n
         x(i) = xcenter(i)
        enddo             
        call func(m,n,x,f)
        call activa(2)
     
        write(6,*)'MAP data written to: map.xy and map.ssq'
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      return
      end
c
c
c*ds
c*ds
       subroutine func(m,nff,x,f)
c      ==========================
c
c ----- calculates deviation data theory for zxssq ----
c
       parameter (mfit=40,msmpl=4000)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*8 dum
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
       logical sqwght,sqwbuf
       logical autox1,autox2
       common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)

       real fcssq
       logical lerrel, lwrtfitdat
       common/cfunce/ lerrel, lwrtfitdat, fcssq
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       real      x(mfit),f(msmpl)
       dimension xp(mtpar),q2(3)
c
       data icall/0/
c
c
c ---- restore startvalues & parameters ----
c
          nfit = 0
          do 10 it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do 20 ip=1,npar
             if(thpsca(ip,it).ne.0) then
               nfit = nfit + 1
               thparx(ip,it) = x(nfit) * thpsca(ip,it)
             endif
20        continue
10      continue
c
c ----- calculate theory-values ------
c
       call couple(1)
       call thc(0)
c      -----------
c
c ----- compare with data ----
       if(lwrtfitdat) then
        open(90,file='fitdat.tmp',status='unknown')
       endif

       m   = 0
       ssq = 0
       do 100 isel=1,nsel
        iad1 = isels(isel)
        iad2 = ifits(isel)
        n    = nwert(iad1)
        mn   = 0
c       --------> counts valid values in the thc-computed data
        do 103 ipt=1,n
          xx=xwerte(ipt,iad1)
          if(((xx.ge.x1).or.autox1).and.((xx.le.x2).or.autox2)) then
c ---------> this if stmnt must correspond to the equival. in thc ! --
             m = m+1
             mn= mn+1
             f(m) = (ywerte(ipt,iad1)-ywerte(mn,iad2))

            if(lwrtfitdat) then
             write(90,*)xx,ywerte(ipt,iad1),ywerte(mn,iad2)
            endif

             if(lerrel) f(m) = f(m) / ywerte(ipt,iad1)
             ferr      = yerror(ipt,iad1)
             ferror(m) = ferr
             if(ferr.gt.0.0d0) then
                f(m) = f(m)/ferr
                ssq  = ssq + f(m)**2                
             else
                ssq  = ssq + f(m)**2
             endif
             if(m.ge.msmpl) then
               write(6,*)' func ==> no. of compared datapts ',m,
     *                   ' exhausts current dimension msmpl=',msmpl
               goto 111
             endif
           endif
c          -----> of if(xx.ge...
103      continue
         if(mn.ne.nwert(iad2)) then
           write(6,*)'func - error(programming!) exp and thc dont match'
           if(lwrtfitdat) close(90)
           return
         endif
100    continue
111    continue
          if(m.eq.0) then
            write(6,*)' no point in fit window !!'
            m=1
            f(1)=9999
          endif
c
       icall = icall + 1
c ---- output if option is set ---
       ssq = ssq/m
       if(iprt.gt.0) write(6,200)icall,ssq,(x(i),i=1,nfit)
200    format(' ',i4,': ssq=',5e12.4/(23x,4e12.4))
       call setudf('ssq0 ',dble(ssq),ier)
       fcssq = ssq
c
c
       if(lwrtfitdat) then
        call theo_out(90)
        close(90)
       endif

       return
       end
c
c
c*ds
c*ds
       subroutine thc(npunkt)
c      ----------------------
c
c ---- compute theoretical scans with a resolution of npoint points
c      of the scans given in /selist/isels(1..nsel)
c      the addresses of the computed scans are put to /selist/ifits
c      if npoint= 0 the sample points match exactly those of the
c      corresponding isels-scan
c
       parameter (mfit=40,msmpl=4000)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*8 dum
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       real*8 getval, dble
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
       logical autox1,autox2
       common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       common /thiadd/iadda
c ---> transfer der addresse der gerade bearbeiteten kurve nach thdatr1
c      damit koennen dann ggf. weitere parameter ueber parget gewonnen
c      werden
c
       dimension dq(4),q2(4)
       logical convolute/.false./, found, folgt
       logical imu
c
       data extup/0./,extlow/0./
 
c
c
       npoint = IABS(npunkt)
       npk    = npunkt
c
       if(found('convolut'))           convolute = .true.
       if(folgt('off     ','convolut'))convolute = .false.
c
       do 2001 i=1,inames
        j = inapa(i)
        if(vname(i).eq.'x1      ') then
            x1     = rpar(j)
            autox1 = .false.
        endif
        if(vname(i).eq.'x2      ') then
            x2     = rpar(j)
            autox2 = .false.
        endif
        if(vname(i).eq.'auto    ') then
            autox1 = .true.
            autox2 = .true.
        endif
2001   continue
c
c
       if(nbuf+nsel.gt.mbuf) then
         write(6,*)'not enough space to store computed spectra!'
         ierrs = 301
         return
       endif
c
       if(iout.gt.2) write(6,*)'thc called with ',nsel,' selected items'
       if(nsel.le.0) then
          write(6,*)'thc: error  no items selected'
          ierrs = 303
          return
       endif
       if(npoint.gt.mwert) then
          write(6,*)'thc warning: no. of points maximum is ',mwert
          npoint = mwert
       endif
c

       if(isigint().gt.0) then
         write(6,*)'Thc Return withut calc due to Sigint:',isigint()
         return
       endif

c ---- compute ....
       do 1001 isel=1,nsel
        iadd = isels(isel)
        iadda= iadd
        if(ifits(isel).eq.0) then
          nbuf        = nbuf+1
          ifits(isel) =  nbuf
        endif
        inbuf = ifits(isel)
c       --------------------> where the data are stored
        if(iout.gt.3) write(6,*)'isel=',isel,'  iadd=',iadd
        if((npoint.le.1).and..not.convolute) then
c                                 ----> take old scan structure
!!!!!!!!
!! erst mal kommentiert lassen, da in diesem Fall auch parget geaendert werden muss
!! eventuell x1 und x2 sowie ggf xc1 und xc2 als standard Parameter mit in die Daten-
!! definition aufnehmen
!!!!!! neu !!!! 
!!          call parget('x1     ',x1,iadd,ier)
!!          x1  = getval('x1      ',dble(x1) ,inew)
!!          if(inew.ne.0 .and. ier.eq.0) then
!!            write(6,*)'commandline spec. of x1 overruled file param.'
!!          endif             
!!          call parget('x2     ',x2,iadd,ier)
!!          x2  = getval('x2      ',dble(x2) ,inew)
!!          if(inew.ne.0 .and. ier.eq.0) then
!!            write(6,*)'commandline spec. of x2 overruled file param.'
!!          endif             
!!!!!! neu !!!! 


          n = nwert(iadd)
          m = 0
          do 9901 i=1,n
       
          if(isigint().gt.2) then
            write(6,*)'THC interupt calc due to Sigint:',isigint()
            return
          endif


           xxx = xwerte(i,iadd)
           if(((xxx.ge.x1).or.autox1).and.((xxx.le.x2).or.autox2)) then
c ---- ! ---- this if condition must meet the corresponding if in the
c             do loop 103 of func ------ ! -----
             m = m+1
             xwerte(m,inbuf) = xxx
           endif
9901      continue
          n = m
        else
c       ----> rescale the gradient
         npoint = intval('n       ',npoint, inew)
         npk    = npoint
         npoint = IABS(npoint)
         if(convolute.and.inew.ne.0) nconvo = npoint
         xmin = xwerte(1,iadd)
         xmax = xwerte(1,iadd)
         do 20001 im=2,nwert(iadd)
          if(xwerte(im,iadd).lt.xmin)   xmin = xwerte(im,iadd)
          if(xwerte(im,iadd).gt.xmax)   xmax = xwerte(im,iadd)
20001    continue
          if(.not.autox1) xmin = x1
          if(.not.autox2) xmax = x2
c --- convolution may require overlapping data --
          if(convolute) then
!!!! neu !!!! 
           call parget  ('xc1     ',xc1, iadd, ier)
           xc1  = getval('xc1     ',dble(xc1) ,inew)
           if(inew.ne.0 .and. ier.eq.0) then
             write(6,*)'commandline spec. of xc1 overruled file param.'             
           endif
           call parget  ('xc2     ',xc2, iadd, ier)
           xc2  = getval('xc2     ',dble(xc2) ,inew)
           if(inew.ne.0 .and. ier.eq.0) then
             write(6,*)'commandline spec. of xc2 overruled file param.'             
           endif
           call parget  ('nconvo  ',xnconvo, iadd,ier)
           nconvo = xnconvo
           if(ier.eq.0) then
             write(6,*)'Param:nconvo(',iadd,') overuled spec. of n'
             write(6,*)'nconvo = ',nconvo
           endif
!!!! neu !!!! 

           xmin = xc1
           xmax = xc2
           npoint = nconvo
           if(npoint.le.0) npoint = nwert(iadd)
c        --- no partial fitting allowed in this version
c            if convolute is used ---> set automatic range on ---
           autox1 = .true.
           autox2 = .true.
          endif
c ------------------------------------------------
          if(npoint.gt.mwert) npoint = mwert
          n = npoint
 
          if(npk.ge.0 
     *      .or. xmin.le.0.0d0 
     *      .or. xmax .le. 0.0d0 
     *      .or. xmax.le.xmin     ) then
           dx = (xmax-xmin)/(n-1)
           do im=1,n
            xwerte(im,inbuf) = xmin + dx * (im-1)
            enddo
          else    ! a neg. n flags that thc computes a even intervals on a log scale
           dx = (xmax-xmin)/(n-1)
           da = alog(xmax/xmin)/(n-1)          
           do im=1,n
            xwerte(im,inbuf) = xmin*exp(da*(im-1))
            enddo
          endif
        endif
    
        nwert(inbuf) = n

c --- compute now...
          do 1003 i=1,n
             x = xwerte(i,inbuf)
c

            sum = thval(x)

            xwerte(i,inbuf) = x
            ywerte(i,inbuf) = sum
            if(iout.gt.5) write(6,*)'i: ',i,'  ywerte=',sum
 
1003     continue

c --- transfer the scanparameters ---
         call txfpar(iadd,inbuf)
         numor(inbuf) =-numor(iadd)
         name(inbuf)  = 'fit'//name(iadd)(1:5)
         np           = nopar(iadd) + 2
         call parset('x1      ',x1,inbuf)
         call parset('x2      ',x2,inbuf)
c ---- vorlaeufig !!! -----
         call parset('x1      ',x1,iadd)
         call parset('x2      ',x2,iadd)
c -------------------------
         if(iout.gt.3)write(6,*)'nbuf=',inbuf,'  numor=',numor(inbuf)
         nwert(inbuf) = n
c ----- convolution ------ ( subroutine has to be supplied ! ) ------
         if(convolute) then
           call datconv(inbuf,xwerte(1,inbuf),ywerte(1,inbuf),n,
     *                  xwerte(1,iadd) ,nwert(iadd))
           nwert(inbuf) = n
         endif
c -------------------------------------------------------------------
c
1001   continue
c
       return
       end
c


       real function thval(x)
c      ----------------------
c
c ---- compute the value of the active set of theories at value x
c
       parameter (mfit=40,msmpl=4000)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*8 dum
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       real*8 getval, dble
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
       logical autox1,autox2
       common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       LOGICAL imu
c
c
             sum = 0
c ----> sum over selected theories ----
             do 10032 it=1,ntheos
               ith = nthtab(it)
               imu = (multflg(it).eq.1)
 
               goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     *              21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     *              38,39,40),
     *              ith
1              continue
               if(imu) then
                 sum=sum*thdatr1(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+thdatr1(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
2              continue
               if(imu) then
                 sum=sum*th2 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th2 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
3              continue
               if(imu) then
                 sum=sum*th3 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th3 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
4              continue
               if(imu) then
                 sum=sum*th4 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th4 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
5              continue
               if(imu) then
                 sum=sum*th5 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th5 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
6              continue
               if(imu) then
                 sum=sum*th6 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th6 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
7              continue
               if(imu) then
                 sum=sum*th7 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th7 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
8              continue
               if(imu) then
                 sum=sum*th8 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th8 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
9              continue
               if(imu) then
                 sum=sum*th9 (x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th9 (x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
10             continue
               if(imu) then
                 sum=sum*th10(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th10(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
11             continue
               if(imu) then
                 sum=sum*th11(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th11(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
12             continue
               if(imu) then
                 sum=sum*th12(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th12(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
13             continue
               if(imu) then
                 sum=sum*th13(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th13(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
14             continue
               if(imu) then
                 sum=sum*th14(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th14(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
15             continue
               if(imu) then
                 sum=sum*th15(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th15(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
16             continue
               if(imu) then
                 sum=sum*th16(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th16(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
17             continue
               if(imu) then
                 sum=sum*th17(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th17(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
18             continue
               if(imu) then
                 sum=sum*th18(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th18(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
19             continue
               if(imu) then
                 sum=sum*th19(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th19(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
20             continue
               if(imu) then
                 sum=sum*th20(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th20(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
21             continue
               if(imu) then
                 sum=sum*th21(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th21(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
22             continue
               if(imu) then
                 sum=sum*th22(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th22(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
23             continue
               if(imu) then
                 sum=sum*th23(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th23(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
24             continue
               if(imu) then
                 sum=sum*th24(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th24(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
25             continue
               if(imu) then
                 sum=sum*th25(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th25(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
26             continue
               if(imu) then
                 sum=sum*th26(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th26(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
27             continue
               if(imu) then
                 sum=sum*th27(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th27(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
28             continue
               if(imu) then
                 sum=sum*th28(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th28(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
29             continue
               if(imu) then
                 sum=sum*th29(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th29(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
30             continue
               if(imu) then
                 sum=sum*th30(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th30(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
31             continue
               if(imu) then
                 sum=sum*th31(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th31(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
32             continue
               if(imu) then
                 sum=sum*th32(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th32(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
33             continue
               if(imu) then
                 sum=sum*th33(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th33(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
34             continue
               if(imu) then
                 sum=sum*th34(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th34(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
35             continue
               if(imu) then
                 sum=sum*th35(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th35(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
36             continue
               if(imu) then
                 sum=sum*th36(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th36(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
37             continue
               if(imu) then
                 sum=sum*th37(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th37(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
38             continue
               if(imu) then
                 sum=sum*th38(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th38(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
39             continue
               if(imu) then
                 sum=sum*th39(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th39(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
40             continue
               if(imu) then
                 sum=sum*th40(x ,thparx(1,it),dum,dum,idum,1)
               else
                 sum=sum+th40(x ,thparx(1,it),dum,dum,idum,1)
               endif
               goto 250
250          continue
10032       continue
c
       thval = sum 
       return
       end
c
c*ds
c*ds
       subroutine parset (pname,pvalue,iadd)
c
c ---- this routine changes or adds a parameter ----
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       character*8 pname
c
c
       np = nopar(iadd)
       do 100 i=1,np
         if (napar(i,iadd).eq.pname) then
           params(i,iadd) = pvalue
           return
         endif
100    continue
       if (np.ge.mpar) then
         write(6,*)'impossible to add ',pname,' to parameterlist of '
     *      ,numor(iadd)
         return
       endif
c
       np = np + 1
       nopar(iadd) = np
       params(np,iadd) = pvalue
       napar(np,iadd) = pname
       return
c
       end
c
c
c*ds
c*ds
       subroutine parget (pname,pvalue,iadd,ier)
c      =========================================
c
c ---- this routine gets the value of a parameter ----
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       character*8 pname
c
c ----------- look for the parameter -----------------------------------
       np = nopar(iadd)
       do 100 i=1,np
         if (napar(i,iadd).eq.pname) then
           pvalue = params(i,iadd)
           ier = 0
           return
         endif
100    continue
c ------------- if not found .. ---------------------------------------
         if(ier.ge.0) then
         if(iout().gt.0)
     *    write(6,*)'impossible to find ',pname,' in parameterlist of '
     *      ,numor(iadd)
          ier = 1
         endif
         return
c
c
       end
c*ds
c*ds
       real*4 function dparam (pname)
c      ==============================
c
c ---- this routine gets the value of a parameter ----
c      besser gekapselte Routine zur Kommunikation mit th-Routinen
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       common /thiadd/iadda       
       character*8 pname
c
c ----------- look for the parameter -----------------------------------
       np = nopar(iadda)
       do i=1,np
         if (napar(i,iadda).eq.pname) then
           dparam = params(i,iadda)
           ier = 0
           return
         endif
       enddo
c ------------- if not found .. ---------------------------------------

       
       dparam = 0.0

       return
c
c
       end
c
c
c*ds
c*ds
       subroutine fpaget(pname,value,nothe,ier)
c      ========================================
c
c ----- get value of the parameter pname of the nothe-th theory -----
c
       parameter(mth=40,mtpar=40,mtcal=40)
c ----- theories common block and definitions ----
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
c  thenam(i)  = name of i-th theory
c  thparn(j,i)= name of j-th parameter of i-th theory
c  nthpar(i)  = no of parameters required for i-th theory
c  thparx(j,l)= parameters for the l-th activated theory
c  thpsca(j,l)= corresponding fit scales   "
c  nthtab(l)  = no. of l-th activated theory
c  ntheos     = total no. of activated theoties
c
       character*8 pname
c
       if(nothe.gt.ntheos) then
         write(6,*)'only ',ntheos,' theories activated !!'
         ier = 1
         return
       endif
c
       nth  = nthtab(nothe)
       npar = nthpar(nth)
       do 100 i=1,npar
        if(thparn(i,nth).eq.pname) then
           value = thparx(i,nothe)
           ier   = 0
           return
        endif
100    continue
c ------------------- not found ----------------------------------------
       write(6,*)'theory parameter: ',pname,' (',nothe,' ) not found'
       ier = 1
       return
       end
c
c
c*ds
c*ds
       subroutine activa(iopt)
c      =======================
c --- activate a theory
c     iopt = 0 ===> activate
c     iopt = 1 ===> desactivate all or selected entries
c     iopt = 2 ===> list actvated theories
c     iopt = 3 ===> reactivate theories from lastth
c
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)
c
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
 
       character*132 xformel,yformel,yfitform
       common/formul/xformel,yformel,yfitform
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
 
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
 
       common/therrc/therro(mtpar,mtcal)
 
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       character*132  rlibuf
       integer       iocbuf
       character*4   label
       character*8   pname
       real*8        dble
       logical       found
c
c
       if(iopt.lt.0.or.iopt.gt.3) then
         write(6,*)' activate invalid option : iopt = ',iopt
         ierrs = 800
         return
       endif
c
       if(iopt.eq.0) then
c      -------------------> activate
         if(inames.eq.0) then
           write(6,*)' theoryname & parameters required'
           ierrs = 800
           return
         endif
         do 10 i=1,mth
          if(vname(1).eq.thenam(i)) then
            ith = i
            goto 11
          endif
10       continue
         write(6,*)' wrong theoryname use theos to list available names'
         ierrs = 800
         return
11       continue
         npar = nthpar(ith)
c        if(ipars.ne.npar*2) then
c          write(6,*)' wrong no. of parameters & scales'
c          write(6,*)' for each parameter a scale must be given!'
c          ierrs = 800
c          return
c        endif
c
         ntheos         = ntheos + 1
         nthtab(ntheos) = ith
         if(found('multiply')) then
           multflg(ntheos) = 1
           inaml           = 2
         else
           multflg(ntheos) = 0
           inaml           = 1
         endif
c
         if(inames.eq.inaml) then
c                        -----> assign noname parameterlist by sequence
           do 20 i=1,npar
              j1 = 2*i - 1
              j2 = j1 + 1
              thparx(i,ntheos) = rpar(j1)
              thpsca(i,ntheos) = rpar(j2)
              thpala(i,ntheos) = '    '
              ncoup(i,ntheos)  = 0
20         continue
           return
         endif
c
        if(inames.ne.npar+inaml) then
          write(6,*)' none or all parameter scale pairs must be given'
          write(6,*)' by names : parnam x s parnam x s ... npar=',npar
          ntheos = ntheos - 1
          ierrs = 800
        endif
        do 30 j=2,inames
          j1 = 2*(j-1) - 1
          j2 = j1 + 1
          do 31 i=1,npar
           if(vname(j).eq.thparn(i,ith)) then
             thparx(i,ntheos) = rpar(j1)
             thpsca(i,ntheos) = rpar(j2)
             ncoup(i,ntheos) = 0
             thpala(i,ntheos)= '    '
             goto 30
           endif
31        continue
          write(6,*)vname(j),' is no parametername ..'
          ntheos = ntheos - 1
          ierrs = 800
          return
30      continue
        write(6,*)vname(1),' is the ',ntheos,'-th  activated entry'
c
        return
       endif
c
       if(iopt.eq.1) then
c                    ----> clear activations
         if(ipars.eq.0.or.ntheos.le.1) then
           ntheos = 0
           write(6,*)' theory activation stack counter set to zero'
           return
         endif
c ---- desactivate only some theories ---
         do 1501 i=1,ipars
           j = rpar(i)+0.001
           if(j.lt.1.or.j.gt.ntheos) then
             write(6,*)' theory to be activated does not exist ..'
             ierrs = 700
             return
           endif
           nthtab(j) = 0
1501     continue
c ---- update the stack ----
         nthn = ntheos
         it   = 1
1511     continue
           if(nthtab(it).eq.0) then
            if(it.eq.nthn)then
              nthn = nthn - 1
              goto 1530
            endif
             do 1512 itt = it+1,nthn
               ith = nthtab(itt)
               if(ith.eq.0) then
                 npars = 1
               else
                 npars = nthpar(ith)
               endif
               nthtab(itt-1)  = nthtab(itt)
               multflg(itt-1) = multflg(itt)
               do 1522 ip=1,npars
                 thparx(ip,itt-1) = thparx(ip,itt)
                 thpsca(ip,itt-1) = thpsca(ip,itt)
                 thpala(ip,itt-1) = thpala(ip,itt)
                 thpaco(ip,itt-1) = thpaco(ip,itt)
                 ncc              = ncoup(ip,itt)
                 ncoup(ip,itt-1)  = ncc
                 do 1522 lp = 1,ncc
                    thpalc(lp,ip,itt-1) = thpalc(lp,ip,itt-1)
                    thpafc(lp,ip,itt-1) = thpafc(lp,ip,itt-1)
1522           continue
1512         continue
             nthn = nthn - 1
             if(it.le.nthn) goto 1511
            endif
            if(it.lt.nthn) then
               it = it + 1
               goto 1511
            endif
1530       continue
           ntheos = nthn
         return
       endif
c
       if(iopt.eq.2) then
c                    ----> output of activated theories
c                    ----> and writing the makro lastth
         write(6,*)' ***** activated theories *****'
         if(ntheos.eq.0) then
          write(6,*)'      none, use  acl to restore lastth setting'
          return
c         -------> preserve the contents of file lastth in this case
         else
c
c ---- > write the file file lastth to restore the theory-params
          kk = 6
c         ------> write first to std-output...
          write(kk,*)'current parametersetting:'
6666      continue
          

          if(kk.ne.6) open(kk,file='lastth',
     *                        form='formatted',status='UNKNOWN')
c                                   ------
          call theo_out(kk)
       
!          do 150 i=1,ntheos
!           ith = nthtab(i)
!           if(multflg(i).eq.1) then
!              write(kk,*)'theory   ',thenam(ith),'   multiply'
!           else
!              write(kk,*)'theory   ',thenam(ith)
!           endif
!           do 150 j=1,nthpar(ith)
!             if(kk.ne.6) then
!               call setudf(thparn(j,ith)//' ',dble(thparx(j,i)),ier)
!             endif
!             write(kk,170)thpala(j,i),thparn(j,ith),thparx(j,i),
!     *                    thpsca(j,i),therro(j,i),
!     *                    (thpalc(l,j,i),thpafc(l,j,i),l=1,ncoup(j,i))
!170        format(1x,a4,1x,a8,1x,e12.4,1x,e8.1,1x,e9.2,
!     *                                              10(1x,a4,1x,f5.2))
!150       continue
!          write(kk,*)'end'
!          if(thenam(nthtab(1)).eq.'eval    ') write(kk,*)yfitform

          if(kk.eq.6) then
            kk = 30
            goto 6666
c           ---------> no save this to file
          else
            close(kk)
          endif
        endif
        return
       endif
c  ------> reactivate theories as stored in lastth
       if(iopt.eq.3) then
c ---- > read the file file lastth to restore the theory-params
          open(30,file='lastth',form='formatted',
     *                          status='UNKNOWN',err=999)
c                       ------
c ---- save the old incom status ----
          rlibuf = reslin
          iocbuf = ioldc
c ---- initialize ----
          ntheos = 0
c ---- read an decode the parameter-file --------
7777      continue
          read(30,'(a)') reslin
c         write(6,*)reslin
          ioldc = 1
          call lowcase(reslin,80)
          call incom(comand)
c      write(6,*)'inames=',inames,' ipars=',ipars,' inpar(1)=',inpar(1),
c    *           ' inapa(1)=',inapa(1),' vname(1..3)=',vname(1),vname(2)
c    *           ,vname(3)
          if(comand.eq.'end     ') goto 8888
          if(comand.eq.'theory  ') then
            ntheos = ntheos + 1
            if(ntheos.gt.mtcal) then
               write(6,*)'too many theories ntheos=',ntheos,' max=',mtca
               ierrs = 703
               goto 999
            endif
c           ------- identify the theory ------
            do 7801 ith=1,mth
               if(vname(1).eq.thenam(ith)) goto 7802
7801        continue
              write(6,*)'theoryname: ',vname(1),' is not known!'
              ierrs = 700
              goto 999
7802        continue
c           --------> ok name is known
            nthtab(ntheos) = ith
            if(found('multiply')) then
              multflg(ntheos) = 1
            else
              multflg(ntheos) = 0
            endif
          else
c           -------- look now for a label -------
            if(ntheos.eq.0) then
               write(6,*)'theory must be specified prior to params'
               goto 999
            endif
            if(iparn(1).ne.0) then
              label = comand(1:4)
              pname = vname(1)
              ilook = 2
            else
              label = '    '
              pname = comand
              ilook = 1
            endif
c ---------- look for the right entry ------------
            do 7803 i=1,nthpar(ith)
              if(pname.eq.thparn(i,ith)) then
                 thpala(i,ntheos) = label
                 thparx(i,ntheos) = rpar(1)
                 thpsca(i,ntheos) = rpar(2)
                 goto 7805
              endif
7803        continue
            write(6,*) pname ,' is an invalid parametername'
            ierrs = 701
            goto 999
7805        continue
c --------- look for coupleded items --------
            ncoup(i,ntheos) = 0
            if(inames.ge.ilook) then
              do 7807 j=ilook,inames
                 ncoup(i,ntheos) = ncoup(i,ntheos) + 1
                 if(ncoup(i,ntheos).gt.mcoup) then
                    write(6,*)'too many couplings max is ',mcoup
                    goto 999
                 endif
                 thpalc(ncoup(i,ntheos),i,ntheos) = vname(j)
                 thpafc(ncoup(i,ntheos),i,ntheos) = rpar(inapa(j))
7807          continue
            endif
          endif
          goto 7777
c ---- restore the incom resline buffers --------
8888      continue
          reslin = rlibuf
          ioldc  = iocbuf
          close(30)
          call couple(0)
          write(6,*)' theory-setting restored'
          return
c
999       continue
          close(30)
          ierrs = 99
          write(6,*)' theories could not be restored ....'
          ntheos = 0
          return
       endif
c
       end
c
c
       subroutine theo_out(kk)
!      -----------------------           
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)

       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       character*8 combinam,cha*1
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
 
       character*132 xformel,yformel,yfitform
       common/formul/xformel,yformel,yfitform
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
 
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
 
       common/therrc/therro(mtpar,mtcal)
 
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)

       integer i, ith, kk

          if(ntheos.le.0) return

          do 150 i=1,ntheos
           ith = nthtab(i)
           if(multflg(i).eq.1) then
              write(kk,*)'theory   ',thenam(ith),'   multiply'
           else
              write(kk,*)'theory   ',thenam(ith)
           endif
           do 150 j=1,nthpar(ith)
            if(kk.ne.6) then
             call setudf(thparn(j,ith)//' ',dble(thparx(j,i)),ier)
             call setudf('e.'//thparn(j,ith)//' ',dble(therro(j,i)),ier)
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! bis zur Erweiterung auf lange Strings hier eine Notloesung zum Setzen von th-params !!!
             write(cha,'(i1)')i
             combinam = thenam(ith)(1:2)//cha//thparn(j,ith)(1:5)
             do lf=1,nsel
              if(ifits(lf).gt.0) then
               call parset(combinam,thparx(j,i),ifits(lf))
               call parset('e'//combinam(1:7),therro(j,i),ifits(lf))
              endif
             enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            endif
             write(kk,170)thpala(j,i),thparn(j,ith),thparx(j,i),
     *                    thpsca(j,i),therro(j,i),
     *                    (thpalc(l,j,i),thpafc(l,j,i),l=1,ncoup(j,i))
170        format(1x,a4,1x,a8,1x,e12.4,1x,e8.1,1x,e9.2,
     *                                              10(1x,a4,1x,f5.2))
150       continue
          write(kk,*)'end'
          if(thenam(nthtab(1)).eq.'eval    ') write(kk,*)yfitform
       
       return
       end

c
       block data bcoup
c      ==========
c ----------------------------------------------------------------------
c ---- this blockdata is used to replace the standard filling of
c      character variables with nulls ! by blanks !
c      this is importatnt for the incom parser, which is only sensitiv
c      to blanks, whereas nulls are written and read if default
c      character varaibles are printed!
c ----------------------------------------------------------------------
c
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)
c
c ---- dervived auxiliary parameters only for block data ----
       parameter(l1=mtpar*mtcal, l2=mcoup*mtpar*mtcal)
       parameter(l3=mth*mtpar)
 
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
c
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
c
c
c --- set the relevant character-varaibles ---
       data thpala/l1*'    '/,thpalc/l2*'    '/
       data thparn/l3*'        '/
       end
c
c
c
c*ds
c*ed
       subroutine couple(iopt)
c      =======================
c ----- couple theory-parameters -----
c
c ----------------------------------------------------------------------
c      the theoryparameters are coupled linearily to other parameters
c      according to the label settings
c      the formula is:
c               p[l]    <== sum(i<>l) { f[l,i] p[i] }  + p0[l]
c               p0[l]   <== p[l|0] - sum(i<>l) { f[l,i] p[i|0] }
c      p[j|0] designates the startvalue of the j-th parameter and
c      f[l,i] is the linearcoupling factor for parameters l with i.
c
c      for iopt=0   the p[j|0] are evaluated and stored on thpaco
c          iopt=1   the actual p[l] are computed --> thparx
c      for a couplede parameter the fitscale thpsca is forced to 0
c
c ----------------------------------------------------------------------
c
c
       parameter(mth=40,mtpar=40,mtcal=40,mcoup=10)
c
       character*8 thenam,thparn,thenax,thpanx(mtpar)
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
c
        character*4 thpala, thpalc
        common /theorc/ thpala(mtpar,mtcal), thpalc( mcoup,mtpar,mtcal),
     *        thpafc(mcoup,mtpar,mtcal), thpaco(mtpar,mtcal),
     *        ncoup(mtpar,mtcal)
c
!!     common/outlev/iout,ibild,ierrs,inka   !! aix
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       character*4   label
c
c ----------------------------------------------------------------------
c
c ----- go through all couplings ----
       do 10 i=1,ntheos
          ith = nthtab(i)
          np  = nthpar(ith)
          do 20 j=1,np
             n = ncoup(j,i)
             y = 0
             if(n.gt.0) then
                do 30 l=1,n
                   label = thpalc(l,j,i)
c                  ---- search now for a parameter with this label ---
                   do 40 ii=1,ntheos
                      ithh = nthtab(ii)
                      npp  = nthpar(ith)
                      do 40 jj=1,npp
                         if(label.eq.thpala(jj,ii)) goto 50
40                 continue
                      write(6,*)'couple: label=',label,' not found'
                      ierrs = 800
50                 continue
                   if((j.eq.jj).and.(i.eq.ii)) then
                      write(6,*)'couple: label=',label,' points to ',
     *                          'itself'
                      ierrs = 801
                   endif
c           ---- do the coupling ----
                   if(ierrs.eq.0) y = y + thpafc(l,j,i) * thparx(jj,ii)
30              continue
                if(iopt.eq.0) then
                   thpaco(j,i) = thparx(j,i) - y
                   if(iout.gt.5) write(6,*)'couple i,j,y =',i,j,y
                   if(thpsca(j,i).ne.0) then
                      write(6,*)'couple: force scale (',j,',',i,') to 0'
                      thpsca(j,i) = 0
                   endif
                else
                   thparx(j,i) = thpaco(j,i) + y
                endif
             endif
c            ----- ( if(n.gt.0).... )
20        continue
10     continue
c
       return
       end
c
c
c*ds
c*ds
       subroutine info
c      ===============
c ---- reads info file & presets informations ----
c
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       open(22,file='datinf',status='old',form='formatted',err=999)
c
       if(inames.eq.0) then
10       continue
         read(22,1,end=999) inline
1        format(a)
         if(inline(1:1).ne.' ') then
           close(22)
           return
         endif
         write(6,2)inline(1:79)
2        format(1x,a)
         goto 10
       endif
c
c ---- else look for special item ---
100    continue
       read(22,1,end=999) inline
         if(inline(1:8).eq.vname(1)) then
200        continue
           write(6,2) inline
           read(22,1,end=999)inline
           if(inline(1:1).ne.' ') then
             close(22)
             return
           endif
           goto 200
        endif
        goto 100
c
999     continue
        close(22)
        write(6,*)' info not available in file datinf'
c
       return
       end
c
c
c*ds
c*ds
       subroutine out_gli
c      ==================
c
       parameter(mth=40,mtpar=40,mtcal=40)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
       common/fslist/isfits(mbuf),nfsel
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       character*12 infile
       character*80 rline
       logical*4    fileda
       dimension sq(3),ihkl(3),xh(3)
c
c -- open the file -----------------------------------------------------
       if(inames.eq.0) then
         ierrs = 1
         write(6,*)'input(e1): filename is lacking!'
         return
       endif
c
       infile = vname(1)        
       do i=1,8
         if(infile(i:i).eq.' ') infile(i:i)='_'
       enddo
       infile(9:12)='.dat'
       inquire(file=infile,exist=fileda)
!!     if(.not.fileda) then
!!       write(6,*)'input(e2): file ',infile,' does not exist!'
!!       ierrs = 3
!!       return
!!     endif
 
       if(nsel.gt.18 ) then
         write(6,*)'gli_out: too many selected items!'           
         ierrs = 3
         return
       endif
 
       open(20,file=infile,status='UNKNOWN')
       write(20,'(A)')(name(isels(i)),i=1,nsel)
       do j=1,nwert(isels(1))
         write(20,'(18(1x,E13.6))')xwerte(j,isels(1)),
     *        (ywerte(j,isels(i)),i=1,nsel)
         xxx = xwerte(j,isels(1))
         do i=1,nsel
           if( xxx-xwerte(j,isels(i)) .ne. 0.0 ) then
             write(6,*) 'Warning: x-value of ',isels(i),' differs!'
           endif
         enddo
       enddo
 
       close(20)
       write(6,*)'File: ',infile,' written!'
       return
       end
c*ds
c*ds
       subroutine input
c      ================
c
       parameter(mth=40,mtpar=40,mtcal=40)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c

       character*140  infile
c       character*8  infile
       character*80 rline
       logical*4    fileda
       dimension sq(3),ihkl(3),xh(3)

       real*8       xshift, yshift, getval


       LOGICAL      plain_numbers, jump_flag

c  neu fuer pfad 
       integer iz
       character*19 uspfad
       character*132 tempfad
       character*8 shortinf

       data tempfad /'./'/
       save tempfad
       save uspfad

c       do iz=1,132
c          if (reslin(iz:iz).eq.' ') then
c             goto 123
c          endif
c       enddo

c 123   continue
c       iz = iz-1
c       ioldc = 0
       uspfad='./'
       

c 
c
c -- open the file -----------------------------------------------------
       if(inames.eq.0) then
         ierrs = 1
         write(6,*)'input(e1): filename is lacking!'
         return
       endif

       xshift = 0.0
       xshift = getval('xshift  ',xshift,inew)
       yshift = 0.0
       yshift = getval('yshift  ',yshift,inew)

       if(xshift .ne. 0.0d0) then
         write(6,'(a,F12.6)') 'modify x-values by adding ',xshift
       endif

       if(yshift .ne. 0.0d0) then
         write(6,'(a,F12.6)') 'modify y-values by adding ',yshift
       endif

c pfad!!
c       tempfad=uspfad//reslin(1:iz)

       if (ioldc.gt.0) then
c          if(iz.gt.0)then
             infile = uspfad//trim(reslin)//vname(1)
             tempfad=trim(uspfad)//trim(reslin)
c        infile = reslin(1:iz)//vname(1)
             write(*,*) infile
             write(*,*) trim(tempfad)
          else 
             infile = trim(tempfad)//vname(1)
c             infile = vname(1)
             write(*,*) infile
             write(*,*) tempfad
c          endif
       endif
       ioldc = 0
c

       inquire(file=infile,exist=fileda)
       if(.not.fileda) then
         write(6,*)'input(e2): file ',infile,' does not exist!'
         ierrs = 3
         return
       endif
       open(20,file=infile,status='UNKNOWN')
c -- try to read the first line --
20000  continue
       plain_numbers = .false.
       jump_flag     = .false.
       read(20,'(a80)',end=10) rline
       goto 20
10     continue
         write(6,*)'input(e2): file ',infile,'is empty or nonexistent'
         ierrs = 2
         close(20)
         return
20     continue
c
       nbuf = nbuf + 1
       if(nbuf.gt.mbuf) then
         write(6,*)'cdata buffer is full!'
         close(20)
         return
       endif
c -- output the identification & header-line
       write(6,*) 'Comment: ',rline
       coment(nbuf)= rline
c
c -- read the items --
2000   continue
       read(20,'(a80)',end=40) rline
c -- detect end of data -----------------------------------------------
       if(rline(1:4).eq.'#nxt') goto 20000
       if(rline(1:4).eq.'#eod') then
         close(20)
         write(6,*)'input: ready'
         return
       endif
c ----------------------------------------------------------------------
c
       call decode(rline)
c      ------------------
       if(ipars.eq.0.and.inames.eq.0) goto 2000
c                                     ---------> skip blank lines
       if(ipars.ne.0.and.inames.eq.0) then  !! assume plain numbers
           name(nbuf) = infile
          yname(nbuf) = 'y-data'
          xname(nbuf) = 'x-data'
          numor(nbuf) = 100*nbuf
          nopar(nbuf) = 1
          params(1,nbuf) = 0.0
          napar(1,nbuf)  = 'padummy'
          lx = 0
          ly = 0
          le = 0
          plain_numbers = .true.
          jump_flag     = .true.
          goto 202       !! -> and try to read the numbers ....
       endif
c                                
c -- distribute datainputs --
c
c -- read the lattice data --
       if(vname(3).eq.'vs      '.or.vname(3).eq.'versus  ') then
c                      --                        ------
           name(nbuf) = vname(1)
          yname(nbuf) = vname(2)
          xname(nbuf) = vname(4)
          numor(nbuf) = rpar(1) * 1.000001
         goto 2000
       endif
c -- atomtypes --
       if(vname(1).eq.'paramete') then
c                      --------
         l = 0
         nopar(nbuf) = 0

!!         l = l+1     
!!         params(l,nbuf) = xshift
!!         napar (l,nbuf) = 'xshift'
!!         l = l+1     
!!         params(l,nbuf) = yshift
!!         napar (l,nbuf) = 'yshift'


101      continue
         read(20,'(a80)',end=40) rline
         call decode(rline)
         if(inames.eq.0) goto 2000
         if(inames.ne.ipars) then
           write(6,*)'inconsistent data (parameters) !'
           close(20)
           return
         endif
         do 103 i=1,inames
          l = l + 1
          if(l.gt.mpar) then
            write(6,*)'too many parameters'
            close(20)
            ierrs = 100
            return
          endif
          params(l,nbuf) = rpar(inapa(i))
          napar (l,nbuf) = vname(i)
103      continue
         nopar(nbuf) = l
         goto 101
       endif


202    continue

       if(vname(1).eq.'values  '.or.jump_flag) then
c                      ------
        jump_flag = .false.

        if(.not.plain_numbers)then
         lx = 0
         ly = 0
         le = 0
           read(20,'(a80)',end=40) rline
           call decode(rline)
        endif 
201     continue
           if(ipars.eq.0.and.inames.eq.0) then
             if(lx.ne.ly) then
               write(6,*)'number of x-values=',lx,' does not match ly=',
     *                    ly
             endif
             nwert(nbuf) = lx
             goto 2000
           endif
c ---- simple list without any text ! ----
           if(inames.eq.0) then
             lx = lx + 1
             ly = ly + 1
             le = le + 1
             if(lx.gt.mwert.or.ly.gt.mwert) then
               write(6,*)'too many x-y values'
               close(20)
               return
             endif
             xwerte(lx,nbuf) = rpar(1) + xshift
             ywerte(ly,nbuf) = rpar(2) + yshift
             yerror(le,nbuf) = rpar(3)
           else
c ----- the more explicit lists -----
             do 203 j=1,inames
              la = inapa(j)
              ln = inpar(j)
              if(vname(j).eq.'x       ') then
                do 205 i=1,ln
                 lx = lx + 1
                 if(lx.gt.mwert) then
                   write(6,*)'too many x-y values'
                   close(20)
                   return
                 endif
                 xwerte(lx,nbuf) = rpar(la+i-1)
205             continue
              endif
              if(vname(j).eq.'y       ') then
                do 207 i=1,ln
                 ly = ly + 1
                 ywerte(ly,nbuf) = rpar(la+i-1)
207             continue
              endif
              if(vname(j).eq.'e       ') then
                do 209 i=1,ln
                 le = le + 1
                 yerror(le,nbuf) = rpar(la+i-1)
209             continue
              endif
203         continue
c
         endif
c

         read(20,'(a80)',end=40) rline
         call decode(rline)
         goto 201
       endif
c
c -- report nonidentified lines --
       write(6,*)'input(w): nonidentified, ignored line:'
       write(6,*)rline
       goto 2000
c
c -- error returns --
40     continue
        write(6,*)'input(e3): unexpected end of file ',
     *             '(not closed by: #eod)'
        if(lx.ne.ly) then
               write(6,*)'number of x-values=',lx,' does not match ly=',
     *                    ly
        endif
        nwert(nbuf) = lx
        
!!
        if(xshift.ne.0.0) then
          call parset ('xshift  ',xshift,nbuf)
        endif
        if(yshift.ne.0.0) then
          call parset ('yshift  ',yshift,nbuf)
        endif
!!
        ierrs = 0
        close(20)
        return
999    continue
        close(20)
       return
       end
c*ds
c*ds
       subroutine inscn
c      ================ input of sv4-scantype data
c
       parameter(mth=40,mtpar=40,mtcal=40)
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
       character*8  infile,intyp,inmode
       character*131 rline
       dimension sq(3),ihkl(3),xh(3)
       dimension ah(mwert),ak(mwert),al(mwert),anue(mwert),atemp(mwert),
     *           afeld(mwert),acounts(mwert),acorrcnt(mwert),
     *           amonitor(mwert),extramon(mwert),atime(mwert)
c
c -- open the file -----------------------------------------------------
       if(inames.lt.2) then
         ierrs = 1
         write(6,*)'inscn(e1): filename and/or ...type is lacking!'
         return
       endif
c
       infile = vname(1)
c      intyp  = vname(2)
c      inmode = 'f'
c      rline  = 'filedef 20 disk '//infile//' '//intyp//' '//inmode//
c    *          ' ( recfm  u '
c      write(6,*) rline(1:80)
c      call obey(rline,80,ier)
       write(6,*)'opening: ',infile,' ...'
       open(20,file=infile,status='UNKNOWN')
c
c -- read the items --
 
2000   continue
 
       nbuf = nbuf + 1
       if(nbuf.gt.mbuf) then
         write(6,*)'cdata buffer is full!'
         close(20)
         return
       endif
 
c -- look for a header --
2001   continue
         read(20,'(a80)',end=998)rline
2002     continue
         write(6,*) rline
       if(rline(1:2).ne.'sc') goto 2001
 
       coment(nbuf)= rline(1:80)
       rline = ' '
       j = 0
       do i=3,20
         if(coment(nbuf)(i:i).eq.'.') goto 1011
         j = j + 1
         rline(j:j) = coment(nbuf)(i:i)
       enddo
1011   continue
       read(rline,*)nscn
 
c --- search for the beginning of the data-block ---
2003   continue
         read(20,'(a131)',end=999)rline
         if(rline(1:12).ne.'  h       k ') goto 2003
         write(6,*)'1:',rline(1:70)
c --- read the data ---
       ndd = 0
2004   continue
         read(20,'(a131)',end=999)rline
ccc      write(6,*)'*:',rline
         if(rline(1:1).ne.' ') goto 2010
         ndd = ndd + 1
         if(ndd.gt.mwert) then
           write(6,*)'warning too many data !'
           ndd = ndd-1
           goto 2010
         endif
         read(rline,*,end=2004)
     *                ah(ndd),ak(ndd),al(ndd),anue(ndd),atemp(ndd),
     *                afeld(ndd),acounts(ndd),acorrcnt(ndd),
     *                amonitor(ndd),extramon(ndd),atime(ndd)
       goto 2004
c ----------------------------------------------------------------------
c
2010   continue
 
c ---- uebertragen der werte -----
       tavera = 0.
       do i=1,ndd
          xwerte(i,nbuf) = ah(i)
          ywerte(i,nbuf) = acounts(i)
          yerror(i,nbuf) = sqrt(acounts(i))
          tavera = tavera + atemp(i)
       enddo
       tavera = tavera / ndd
 
        name(nbuf) = coment(nbuf)(1:8)
        yname(nbuf) = 'counts'
        xname(nbuf) = 'q'
        numor(nbuf) = nscn
        nwert(nbuf) = ndd
 
c --- setzen von parametern ----
       napar (1,nbuf)  = 'monitor'
       params(1,nbuf)  = amonitor(1)
       napar (2,nbuf)  = 'temp'
       params(2,nbuf)  = tavera
       nopar(nbuf)     = 2
 
       nbuf = nbuf + 1
       if(nbuf.gt.mbuf) then
         write(6,*)'cdata buffer is full!'
         close(20)
         return
       endif
       goto 2002
c      ---------> read the next scan
 
998    continue
       nbuf = nbuf-1
999    continue
       close(20)
       return
       end
c
c
c*ds
c*ds
       subroutine decode(tline)
c      ========================
c
c ---- decode provides a means to treat a 80 characters line by the
c      incom facilities (comand is always set to '&')
c      the results are stored in the usual fashion in /cincom/
c
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c
       character*80 tline
c
c
       if(iout.gt.2) then
         write(6,*)'Decode :',tline
       endif

       if(tline(79:80).ne.'  ') then
         write(6,*)'decode(w): truncating ',tline(79:80),'  of:'
         write(6,*)tline
       endif
c
       ioldc  = 1
       reslin = '& '//tline(1:78)

       if(iout.gt.2) then
         write(6,*)'Decode2:',reslin
       endif

       call incom(comand)

       if(iout.gt.2) then
         write(6,*)'Decode3:',comand
       endif
c
c -- write control info if iout is greater 3 ---------------------------
       if(iout.gt.2) then
         write(6,*)tline
         write(6,*)'inames=',inames,'  ipars=',ipars
         if(inames.ne.0) then
          do 100 i=1,inames
           write(6,10)vname(i),(rpar(j),j=inapa(i),inapa(i)+inpar(i)-1)
10         format(1x,a8,2x,10f12.6)
100       continue
         else
          write(6,200)(i,rpar(i),i=1,ipars)
200       format(' decode numerical parameters:'/
     *          (1x,i2,':  ',f12.6))
         endif
       endif
c-----------------------------------------------------------------------
c
       return
       end
c
c
       function smirro(y,y2,xk0,k1,k2,n)
c      =================================
c --- computes the mirrorimage on the k1 to k2 region mirrored at xk0
c     onto the k1 to k2 region of y2 ( y2 is set to zero for the rest)
c     the squared error difference of y1 and y2 is retuned as function
c     value
c     n is the total number of valid points
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c
       dimension y(mwert),y2(mwert)
c
        do 10 i=1,n
         y2(i) = 0
10      continue
c
       sum = 0
       do 100 i=k1,k2
        xmirro = 2*xk0 - i
        km     = xmirro
        kp     = km + 1
        t      = xmirro - km
        if(kp.gt.n.or.km.lt.1) then
          sum = sum + 1.e20
          goto 100
        endif
          y2(i) = y(km) + (y(kp)-y(km))*t
          sum   = (y(i)-y2(i))**2 + sum
100    continue
       smirro = sum
       return
       end
c
c
       subroutine detxk1(y,y2,xk0,k1,k2,n)
c      ===================================
c
c ---- determine the xk0 - value with the minimum error to miirorsymm.
c      for the description of parameters see function smirro
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       dimension y(mwert),y2(mwert)
c
c
       m  = 2
       f1 = 1.e30
       step = 5
       x1   = xk0
       do 100 j=1,20
        step = step/m
        xk00 = x1
         do 100 i=-m,m
           x = xk00+ step*i
           f = smirro(y,y2,x,k1,k2,n)
           if(f.lt.f1) then
             f1 = f
             x1 = x
           endif
          write(6,*)j,i,' x x1 ',x,x1
          write(6,*)' f f1 ',f,f1
100    continue
        xk0 = x1
c
        return
        end
c
c
       subroutine detxk0(y,y2,xk0,k1,k2,n)
c      ===================================
c
c ---- determine the xk0 - value with the minimum error to miirorsymm.
c      for the description of parameters see function smirro
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       dimension y(mwert),y2(mwert)
c
        f0 = smirro(y,y2,xk0,k1,k2,n)
        do 100 i=1,10
          xk00 = xk0
          x1 = xk0 + 0.5
          x2 = xk0
          x3 = xk0 - 0.5
          f1 = smirro(y,y2,x1,k1,k2,n)
          f2 = smirro(y,y2,x2,k1,k2,n)
          f3 = smirro(y,y2,x3,k1,k2,n)
          f12 = (f1-f2)/(x1-x2)
          f23 = (f2-f3)/(x2-x3)
          a   = (f12-f23)/(x1-x3)
          xk0 = -0.5* (f12/a - x1 - x2)
          write(6,*)'xk00=',xk00,'  xk0=',xk0,'  f2=',f2
          if(abs(xk0-xk00).lt.0.001) goto 200
100     continue
200     continue
        if(f2.gt.f0) xk0 = 0
        return
        end
c
c
       subroutine txfpar(ia,ib)
c      ========================
c --- copies all accompanying infos from data slot ia to slot ib
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       np=nopar(ia)
       nopar(ib)=np
        do 10 i=1,np
          params(i,ib)=params(i,ia)
          napar(i,ib)=napar(i,ia)
10      continue
       name(ib)  = name(ia)
       coment(ib)= coment(ia)
       xname(ib) = xname(ia)
       yname(ib) = yname(ia)
       nwert(ib) = nwert(ia)
       numor(ib) = numor(ia)
c
       return
       end
c
c
       subroutine txfera(ia,ib)
c      ========================
c ---- copy the contens of data slot ia into slot ib -----
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       np=nopar(ia)
       nopar(ib)=np
       do 5 i=1,nwert(ia)
          xwerte(i,ib) = xwerte(i,ia)
          ywerte(i,ib) = ywerte(i,ia)
5      continue
       do 10 i=1,np
         params(i,ib)=params(i,ia)
         napar(i,ib)=napar(i,ia)
10     continue
       name(ib)  = name(ia)
       coment(ib)= coment(ia)
       xname(ib) = xname(ia)
       yname(ib) = yname(ia)
       nwert(ib) = nwert(ia)
       numor(ib) = numor(ia)
c
       return
       end
c
c
c*ds
c*ds
       subroutine savdat(fname,ispc)
c      =============================
c -- save data on adress ispc onto a file named file fname a --
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       character*8 fname,finame
c
c ----------> transfer new name if not scratchfile for edit <-----------
         if(fname.eq.'datbuf  ') then
            finame = name(ispc)
         else
            finame = fname
         endif
c ----------> write data onto buffer <---------------------------------
         open(18,file=fname,status='UNKNOWN')
         write(18,'(a)')coment(ispc)
         write(18,'(a,a,a,a,a,i14)')finame,' ',yname(ispc)(1:20),
     *                                   ' vs ',xname(ispc)(1:20)
     *             ,numor(ispc)
         write(18,*)'parameter'
         write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i=
     *                                 1,nopar(ispc))
         write(18,*)' '
         write(18,*)'values'
         write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),
     *                 i=1,nwert(ispc))
! 501      format(2x,'x  ',e14.7,5x,'y  ',e14.7,5x,'e  ',e14.7)
501      format(2x,'   ',e14.7,5x,'   ',e14.7,5x,'   ',e14.7)
         write(18,*)' '
         write(18,'(a4)')'#eod'

         call theo_out(18)

         close(18)
c ----------------------------------------------------------------------
         write(6,*)'file(',ispc,') : ',name(ispc),
     *             ' has been saved to  ',
     *              fname
       return
       end
c*ds
c*ds
       subroutine msavdat(fname)
c      =========================
c -- MULTI save data on adress ispc onto a file named file fname a --
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c  isels(i)  = address (/spectr/) of selected scan
c  ifits(i)  = adress of fitted spectrum (isels(i))
c  nsel      = length of this table
c
       common/fslist/isfits(mbuf),nfsel
c  isfits    = address (/spectr/) of selected fits
c  nfsel     = length of this table
c
       character*8 fname,finame
c
c ----------> write data onto buffer <---------------------------------

! first write data selected:

         open(18,file=fname,status='UNKNOWN')

         do l=1,nsel
         
         ispc = isels(l)

         write(18,'(a)')coment(ispc)
         write(18,'(a,a,a,a,a,i14)')name(ispc)(1:8)
     *             ,' ',yname(ispc)(1:20),
     *                                   ' vs ',xname(ispc)(1:20)
     *             ,numor(ispc)
         write(18,*)'parameter'
         write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i=
     *                                 1,nopar(ispc))
         write(18,*)' '
         write(18,*)'values'
         write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),
     *                 i=1,nwert(ispc))
! 501      format(2x,'x  ',e14.7,5x,'y  ',e14.7,5x,'e  ',e14.7)
501      format(2x,'   ',e14.7,5x,'   ',e14.7,5x,'   ',e14.7)
         write(18,*)' '
         if(l.lt.nsel .or. ifits(l).gt.0) write(18,'(a4)')'#nxt'

         if(ifits(l).gt.0) then
         ispc = ifits(l)
         write(18,'(a)')coment(ispc)
         write(18,'(a,a,a,a,a,i14)')name(ispc)(1:8)
     *             ,' ',yname(ispc)(1:20),
     *                                   ' vs ',xname(ispc)(1:20)
     *             ,numor(ispc)
         write(18,*)'parameter'
         write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i=
     *                                 1,nopar(ispc))
         write(18,*)' '
         write(18,*)'values'
         write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),
     *                 i=1,nwert(ispc))
         write(18,*)' '
         if(l.lt.nsel) write(18,'(a4)')'#nxt'

         endif

         enddo

         write(18,'(a4)')'#eod'

         call theo_out(18)

         close(18)
c ----------------------------------------------------------------------
         write(6,*)'data and fits have been saved to: ',fname
       return
       end
c
c
c*ds
c*ds
       function fsplin(x)
c      ==================
c --- spline interpolated data evaluation ---
c     spline coefficients are taken from the last call of spline
       parameter(mwert=1024,mbuf=200,mpar=200)
       common/cfc/qziel,cscoef(4,mwert),break(mwert),weight(mwert),
     *            numspl,nwspl
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
         nn  = nwspl
         do 2705 j=1,nn
           if(break(j).ge.x) goto 2707
2705     continue
           if(iout.ge.0) then
             write(6,*)'fsplin: x-value is out of range x =',x,xx
             write(6,*)'range:',break(1),break(nn)
           endif
           fsplin = 0
           return
2707     continue
         zn = x - break(j)
c  ---- compute the interpolated data ----
         yd = ((cscoef(4,j)*zn+cscoef(3,j))*zn+cscoef(2,j))*zn+
     *          cscoef(1,j)
c
         fsplin = yd
c        -----------
       return
       end
c
c
c*ds
c*ds
       function fdes(x)
c      ================
c --- integrand for desmearing ---
       parameter(mwert=1024,mbuf=200,mpar=200)
       common/cfc/qziel,cscoef(4,mwert),break(mwert),weight(mwert),
     *            numspl,nwspl
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
         nn  = nwspl
         xx = sqrt(qziel**2 + x**2)
         do 2705 j=1,nn
           if(break(j).ge.xx) goto 2707
2705     continue
             write(6,*)'fdes x-value is out of range x =',x,xx
             write(6,*)'range:',break(1),break(nn)
             fdes = 0
             return
2707     continue
         zn = xx - break(j)
c  ---- compute the first derivative ----
         yd =  (3*cscoef(4,j)*zn+2*cscoef(3,j))*zn+cscoef(2,j)
c
         fdes = yd / xx
c        --------------
       return
       end
c
c
c*ds
c*ed
       subroutine fftmx(x,y,t,xmax,dx,nx,nfft)
c      =======================================
c  --- small angle multiple scattering evaluator (test-version) ---
c
c   -- it is assumed that the input data (derived form the spline
c      coefficients by fsplin) represent some scattering law in the
c      form:
c      yi <---  i0 * sigma(x) * [unit sample thickness:d]   =
c               lim[d-->0] ( i-measured(x,d) / d )
c      the output data on y will the be scaled such that they represent
c      yo <---  i0 * i-measured(x,d) / d      (note: no limit !!)
c      ---------------------------------
c      i.e. the scattering intensity including all orders of multiple
c      elastic small angle scattering will be calculated for a sample
c      of unit thickness assuming the same primary intensity factor,i0,
c      as for the input data.
c ----------------------------------------------------------------------
c      input variables:
c      t ....... : transmission reduction factor due to small angle sc.
c      xmax .... : largest scattering angle (or q) to be considered
c      dx ...... : increment of x
c      nx ...... : no. of points to be generated
c      nfft .... : fft no. of points (optimal choice nfft=2**m)
c      output variables:
c      x(1..nfft/2) : output x-values
c      y(1..nfft/2) : output y-values
c ----------------------------------------------------------------------
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c --- fft-dimensioning -------------------------------------------------
       parameter(mdim=1024)
       parameter(lda = mdim+1)
       complex*8  ca,cexp,clog,cabs
       common/fftwrk/ca(lda,lda)
c             ---> this common saves storage by use of ca in fftrmx also
       dimension x(*),y(*)
c ----------------------------------------------------------------------
       data pi/3.1415926535897/
c ----------------------------------------------------------------------
c
       if(nfft.eq.0)    nfft = mdim
       if(nfft.gt.mdim) nfft = mdim
       if(nx.gt.nfft/2) nx   = nfft/2
       dd = xmax / (nx-1)
       if(dx.lt.dd) dx = dd
       write(6,*)'fftmx uses the following parameters:'
       write(6,*)'nfft ........... :',nfft
       write(6,*)'nx ............. :',nx
       write(6,*)'xmax ........... :',xmax
       write(6,*)'dx ............. :',dx
       write(6,*)'transmission ... :',t
       if(t.le.0. .or. t.ge.1. ) then
         write(6,*)'fftmx: transmission value is out of range !'
         ierrs = 100
         return
       endif
c ----- period of fft is:
       xl = nfft * dx
c ----- reciprocal fourier increment:
       dy = 1 / xl
c ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
100    continue
c ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
c ---- init data fields ----
       do 110 i=1,nfft
        do 110 j=1,nfft
         ca(i,j) = (0.,0.)
110    continue
       ioold= iout
       iout = -1
       do 200 i=1,nfft/2
        do 200 j=1,nfft/2
         x1 = (i-1) * dx
         x2 = (j-1) * dx
         xx = sqrt(x1**2 + x2**2)
         if(xx.gt.xmax) goto 200
         yy = fsplin(xx) * dx * dx
         ca(       i,       j) = cmplx(yy,0.)
         ca(nfft+2-i,       j) = ca(i,j)
         ca(nfft+2-i,nfft+2-j) = ca(i,j)
         ca(       i,nfft+2-j) = ca(i,j)
200    continue
       iout = ioold
c
       write(6,*)'starting fft....'
       call  fft2d(nfft,nfft,ca,lda,ca,lda)
       write(6,*)'ca(1,1)=',ca(1,1)
       ai1= cabs(ca(1,1))
       ai0= ai1          /(1.-t)
c ---  === : i0 intensity factor provided the original data do not
c            contain any primary beam components. otherwise the
c            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
       xn = -at / cabs(ca(1,1))
c ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
c      xn = normalisation factor of the fourier transform to meet the
c           absolute scattering cross-section scale
       do 150 i=1,nfft
        do 150 j=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i,j)
c
         ca(i,j) = cexp( xn*ca(i,j) + at ) - t
c ----                                       = <- no primary beam int.
c                                              will be in the results !
c                                              omitt this t if primary
c                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i,j)
150    continue
c
c ---- reverse transformation back to scattering angle field ----
       call  fft2b(nfft,nfft,ca,lda,ca,lda)
c
c ---- return values ---------------------------------------------------
c ---- normalize to the same total scattering --------------------------
       xn = ai1 * dy * dy / (-at)
       do 1000 i=1,nfft/2
        if(iout.gt.0) write(6,*)'result(',i,')= ',ca(i,1) * xn
        y(i) = real(ca(i,1)) * xn
1000   continue
c ----------------------------------------------------------------------
c
       return
       end
c*ds
c*ed
       subroutine fftrmx(x,y,ai0,t,thick,xmax,dx,nx,nfft)
c      ==============================================
c  --- small angle multiple scattering inverse evaluator ----
c   -- it is assumed that the spline coefficients that are to be
c      generated by the datreatx command spline from the data to be
c      cleaned from multiple scattering will be used by fsplin to
c      provide equally spaced interpolated (and smoothed) original
c      data for the first fourier transform.
c      the output will then be:
c      y <--- i0 * sigma(x)
c      provided the appropriate value for the sample thickness has been
c      supplied.
c      the primary intensity factor has been included to faciliate the
c      comparison of original data and cleaned data on the same scale
c      (in that case thickness = 1 may be appropriate).
c      i0 will be printed to allow for a subsequent computation of sigma
c      sigma corresponds the cross-section per unit area and thickness
c      times the solid angle represented by the spacing dx of x-values.
c      solid angle ---> dx * dx /(ki**2)
c ----------------------------------------------------------------------
c      input variables:
c      t ....... : transmission reduction factor due to small angle sc.
c      thick ... : sample thickness (units will be the units of the res.
c      xmax .... : largest scattering angle (or q) to be considered
c      dx ...... : increment of x
c      nx ...... : no. of points to be generated
c      nfft .... : fft no. of points (optimal choice nfft=2**m)
c      output variables:
c      x(1..nfft/2) : output x-values
c      y(1..nfft/2) : output y-values
c      ai0 ........ : primary intensity varaible as assumed (i0)
c ----------------------------------------------------------------------
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c --- fft-dimensioning -------------------------------------------------
       parameter(mdim=1024)
       parameter(lda = mdim+1)
ccc    complex*8  ca,cexp,clog
       common/fftwr2/ca(lda,lda)  !! aix
c             ---> this common saves storage by use of ca in fftrmx also
       dimension x(*),y(*)
c ----------------------------------------------------------------------
       data pi/3.1415926535897/
c ----------------------------------------------------------------------
c
       if(nfft.eq.0)       nfft = mdim
       if(nfft.gt.mdim)    nfft = mdim
       if(nx.gt.nfft) nx        = nfft
       dd = xmax / (nx-1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   if(dx.lt.dd) dx = dd
       dx = dd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(6,*)'fftrmx(real) uses the following parameters:'
       write(6,*)'nfft ........... :',nfft
       write(6,*)'nx ............. :',nx
       write(6,*)'xmax ........... :',xmax
       write(6,*)'dx ............. :',dx
       write(6,*)'transmission ... :',t
       write(6,*)'thickness(sample):',thick
       if(t.le.0. .or. t.ge.1. ) then
         write(6,*)'fftmx: transmission value is out of range !'
         ierrs = 100
         return
       endif
       if(thick.le.0) then
         write(6,*)'sample thickness of :',thick,' is replaced'
         thick= 1
         write(6,*)'                 by :',thick
       endif
c ----- period of fft is:
ccccc  xl =     nfft * dx
       xl = 2 * nfft * dx
c ----- reciprocal fourier increment:
       dy = 1 / xl
c ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
100    continue
c ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
c ---- init data fields ----
       do 110 i=1,nfft
        do 110 j=1,nfft
         ca(i,j) = 0.
110    continue
       ioold= iout
       iout = -1
       do 200 i=1,nfft
        do 200 j=1,nfft
         x1 = (i-1) * dx
         x2 = (j-1) * dx
         xx = sqrt(x1**2 + x2**2)
         if(xx.gt.xmax) goto 200
ccccc    yy = fsplin(xx) * dx * dx
         yy = yinterp(xx) * dx * dx
         ca(       i,       j) = yy
200    continue
       iout = ioold
c
       write(6,*)'starting fft....'
       call  rfft2d(nfft,nfft,ca,lda,ca,lda)
c!!!!  ai0= cabs(ca(1,1))/(1.-t)
       ai0=     (ca(1,1))/(1.-t)
c!!!!
c ---  === : i0 inttensity factor provided the original data do not
c            contain any primary beam components. otherwise the
c            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
ccccc  xn = -at / cabs(ca(1,1))
c ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
c      xn = normalisation factor of the fourier transform to meet the
c           absolute scattering cross-section scale
       do 150 i=1,nfft
        do 150 j=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i,j)/ai0 + t
c
c!!!!c   ca(i,j) = clog( ca(i,j)/ai0       + t) - at
         rca     = ca(i,j)
         ca(i,j) = alog( rca          /ai0 + t) - at
c!!!!c
c ----                                       = <- no primary beam int.
c                                              will be in the results !
c                                              omitt this t if primary
c                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i,j)
150    continue
c
c ---- reverse transformation back to scattering angle field ----
       call  rfft2d(nfft,nfft,ca,lda,ca,lda)
c
c ---- return values ---------------------------------------------------
       do 1000 i=1,nfft/2
        if(iout.gt.0) write(6,*)'result(',i,')=',ca(i,1)*dy*dy*ai0/thick
cc      y(i) = real(ca(i,1)) * dy * dy * ai0 / thick
        y(i) =     (ca(1,i)) * dy * dy * ai0 / thick
1000   continue
       write(6,*)'intensity factor i0 =',ai0
       write(6,*)'                     ','=================='
c ----------------------------------------------------------------------
c
c
       return
       end
c*ds
c*ed
       subroutine rmx1d(x,y,ai0,t,thick,xmax,ymax,nx,nfft)
c      =================================================
c  ---> like fftrmx but using 1d-fourier-transforms <---
c  --- small angle multiple scattering inverse evaluator ----
c   -- it is assumed that the spline coefficients that are to be
c      generated by the datreatx command spline from the data to be
c      cleaned from multiple scattering will be used by fsplin to
c      provide equally spaced interpolated (and smoothed) original
c      data for the first fourier transform.
c      the output will then be:
c      y <--- i0 * sigma(x)
c      provided the appropriate value for the sample thickness has been
c      supplied.
c      the primary intensity factor has been included to faciliate the
c      comparison of original data and cleaned data on the same scale
c      (in that case thickness = 1 may be appropriate).
c      i0 will be printed to allow for a subsequent computation of sigma
c      sigma corresponds the cross-section per unit area and thickness
c      times the solid angle represented by the spacing dx of x-values.
c      solid angle ---> dx * dx /(ki**2)
c ----------------------------------------------------------------------
c      input variables:
c      t ....... : transmission reduction factor due to small angle sc.
c      thick ... : sample thickness (units will be the units of the res.
c      xmax .... : largest scattering angle (or q) to be considered
c      ymax .... : fouriertransform-range
c      nx ...... : no. of points to be generated
c      nfft .... : fft no. of points (optimal choice nfft=2**m)
c      output variables:
c      x(1..nfft/2) : output x-values
c      y(1..nfft/2) : output y-values
c      ai0 ........ : primary intensity varaible as assumed (i0)
c ----------------------------------------------------------------------
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c --- fft-dimensioning -------------------------------------------------
       parameter(mdim=1024)
       dimension  ca(mdim)   ,cb(mdim)
c             ---> this common saves storage by use of ca in fftrmx also
       dimension x(*),y(*)
c ----------------------------------------------------------------------
       data pi/3.1415926535897/
c ----------------------------------------------------------------------
c
       if(nfft.eq.0)         nfft = mdim
       if(nfft.gt.mdim     ) nfft = mdim
       if(nx.gt.mdim  )      nx   = mdim
       dx = xmax / (nx-1)
       write(6,*)'fftrmx uses the following parameters:'
       write(6,*)'nfft ........... :',nfft
       write(6,*)'nx ............. :',nx
       write(6,*)'xmax ........... :',xmax
       write(6,*)'ymax ........... :',ymax
       write(6,*)'dx ............. :',dx
       write(6,*)'transmission ... :',t
       write(6,*)'thickness(sample):',thick
       if(t.le.0. .or. t.ge.1. ) then
         write(6,*)'fftmx: transmission value is out of range !'
         ierrs = 100
         return
       endif
       if(thick.le.0) then
         write(6,*)'sample thickness of :',thick,' is replaced'
         thick= 1
         write(6,*)'                 by :',thick
       endif
c ----- reciprocal fourier increment:
       dy = ymax / (nfft-1)
c ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
100    continue
c ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
c ---- init data fields ----
!!??   iout = ioold
c
       write(6,*)'starting transform ....'
       do i=1,nfft
        qq  = dy*(i-1)
        sum = 0
        do j=1,nx
          xx = dx*(j-1)
cccc      sum = sum + fsplin(xx)*xx*bsj0(xx*qq)
          sum = sum + yinterp(xx)*xx*bsj0(xx*qq)
        enddo
        cb(i) = sum*dx
        if(iout.gt.0)write(6,*)i,qq,cb(i)
       enddo
 
       write(6,*)'cb(1) = ',cb(1)
       ai0=     (cb(1))/(1.-t)
c ---  === : i0 inttensity factor provided the original data do not
c            contain any primary beam components. otherwise the
c            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
cccc   xn = -at /      cb(1)
c ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
c      xn = normalisation factor of the fourier transform to meet the
c           absolute scattering cross-section scale
       do 150 i=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i)/ai0 + t
c
         ca(i  ) = alog( cb(i  )/ai0       + t) - at
c ----                                       = <- no primary beam int.
c                                              will be in the results !
c                                              omitt this t if primary
c                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i  )
150    continue
c
c ---- reverse transformation back to scattering angle field ----
       write(6,*)'reverse transform...'
       do i=1,nx
        xx  = dx*(i-1)
        sum = 0
        do j=1,nfft
          qq = (j-1)*dy
          sum= sum + ca(j)*qq*bsj0(xx*qq)
        enddo
        cb(i) = sum*dy
        if(iout.gt.0) write(6,*)i,xx,cb(i)
       enddo
c
c ---- return values ---------------------------------------------------
       do 1000 i=1,nx
        if(iout.gt.0) write(6,*)'result(',i,')=',cb(i  )*ai0/thick
        y(i) = (cb(i  ))* ai0 / thick
1000   continue
       write(6,*)'intensity factor i0 =',ai0
       write(6,*)'                     ','=================='
c ----------------------------------------------------------------------
c
c
       return
       end
c
c
       subroutine rfft2d(n1,n2,a,lda,b,ldb)
c      ==================>==>==>==>==<==>==
c                              a === b ok.
c
c 2d-reelle cos-fouriertransformation
c -----------------------------------
       parameter (mdim=1024)
       dimension a(lda,lda)
       dimension b(ldb,ldb)
       dimension y(mdim)
 
 
       if(n1.gt.mdim.or.n2.gt.mdim) then
          write(6,*)'rfft2d: n1=',n1,' or n2=',n1,' > mdim!'
          if(n1.gt.mdim) n1=mdim
          if(n2.gt.mdim) n2=mdim
       endif
 
       do i=1,n2
          call fcost(n1,a(1,i),b(1,i))
       enddo
 
       do j=1,n1
          do i=1,n2
             y(i) = b(j,i)
          enddo
          call fcost(n2,y,y)
          do i=1,n2
             b(j,i) = y(i)
          enddo
       enddo
 
       return
       end
c*ds
c*ed
c
       function yinterp(x)
c      ---------------------
c --> interpolates y-value at x for the file on top of the
c     selection list
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c --- mwert  = max. no. of x-y-values in one buffer
c     mbuf   = max. no. of different buffers
c     mpar   = max. no. of parameters associated with one buffer
c ---  maximum scan length ....
       parameter(mth=40,mtpar=40,mtcal=40)
c ---  fit dimensions ---
       parameter (mfit=40,msmpl=4000)
c  -- mfit = max no. of fitted parameters
c     msmpl= max no. of datapoints in fit
c
c --- incom common-section ---
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c --- xwerte(i,j)   x-values on buffer j
c     ywerte(i,j)   y-values on buffer j
c     yerror(i,j)   error of y-values   (only supported by some fktn)
c     xname(j)      name of x-values (x-axis) for buffer j
c     yname(j)       "   "  y-   "    y-  "    "    "    "
c     name(j)       short text-identifier for data on buffer j
c     nwert(j)      no. of valid points on buffer j
c     nbuf          no. of filled buffers
c     numor(j)      numerical idenfication of data on buffer j
c     coment(j)     one line of comment describing data on buffer j
c     params(l,j)   set of parameters associated with data on buffer j
c     napar(l,j)    names of these parameters
c     nopar(j)      no. of valid parameters
c
c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c  isels(i)  = address (/spectr/) of selected scan
c  ifits(i)  = adress of fitted spectrum (isels(i))
c  nsel      = length of this table
c
       common/fslist/isfits(mbuf),nfsel
c  isfits    = address (/spectr/) of selected fits
c  nfsel     = length of this table
c
c ----------------------------------------------------------------
       ia = isels(1)
       n  = nwert(ia)
       if(x.lt.xwerte(1,ia).or.x.gt.xwerte(n,ia)) then
          yinterp = 0.0
          return
       endif
       do i=2,n
        if(xwerte(i,ia).gt.x) then
           yinterp = ywerte(i,ia)*(x           -xwerte(i-1,ia))/
     *                            (xwerte(i,ia)-xwerte(i-1,ia))
     *             + ywerte(i-1,ia)*(xwerte(i,ia)-x)/
     *                              (xwerte(i,ia)-xwerte(i-1,ia))
           return
        endif
       enddo
       yinterp = 0.0
       return
       end
c*ds
c*ed
       subroutine demux(xin,yin,nin,t,xmax,nfft,iout,xout,yout,ier)
c      ============================================================
c  --- small angle multiple scattering inverse evaluator ----
c      input:
c      xin(1..nin) ... : inputvector containig scatteringangles
c                        or equivalent in ascending order
c      yin(1..nin) ... : inputvector containing scattering intensities
c      nin ........... : active length of xin and yin
c      t ............. : effective transmission factor associated
c                        with the small angle scattering
c      xmax .......... : x (scattering vector/angle) range
c                        for the result
c      nfft .......... : number of points of the internal
c                        fouriertransforms = the number of
c                        points in the output vectors.
c                        preferred value: 2**n.
c      iout .......... : variable to control additional printoutput
c                        from this subroutine
c      output:
c      xout(1..nfft).. : output x-values
c      yout(1..nfft).. : output y-values
c      ier ........... : errorindicator (0=ok)
c
c      the output will be the hypothetical scattering data of the
c      same sample using the same apparatus
c      but without any effects of multiple scattering,
c      e.g. no multiple sacttered intensity and no transmission
c      reduction.
c ----------------------------------------------------------------------
c
c ---- dimensioning ----------------------------------------------------
       parameter(pi=3.1415926535897)
       parameter(mdim=1024)
       real xin (nin) ,yin (nin)
       real xout(nfft),yout(nfft)
       common/fftwr1/  ca(mdim,mdim), yinter(mdim)     !! k--> 1 aix !!
c ----------------------------------------------------------------------
       ier = 0
       if(nfft.eq.0)       nfft = mdim
       if(nfft.gt.mdim)    nfft = mdim
 
       if(iout.gt.0) then
         write(6,*)'demux uses the following parameters:'
         write(6,*)'nfft ........... :',nfft
         write(6,*)'xmax ........... :',xmax
         write(6,*)'transmission ... :',t
       endif
 
       if(t.le.0. .or. t.ge.1. ) then
         write(6,*)'demux: transmission value is out of range !'
         ier = 1
         return
       endif
 
c------- data preparation ----------------------------------------------
c check for the first and last nonzero element
         ifirst = 0
         ilast  = nin
         do i=1,nin
           if(yin(i).ne.0.0) then
              ifirst = i
              goto 1
           endif
         enddo
1        continue
         do i=nin,1,-1
           if(yin(i).ne.0.0) then
              ilast  = i
              goto 2
           endif
         enddo
2        continue
         if(ifirst.eq.0) then
           write(6,*)'demux: data are identical zero !'
           ier = 3
           return
         endif
c extrapolate data to x --> 0 assuming a parabola
         slope = (yin(ifirst+2)-yin(ifirst))
     *          /(xin(ifirst+2)-xin(ifirst))
         apar  = slope / (2.0*xin(ifirst+1))
         bpar  = yin(ifirst+1) - slope*0.5*xin(ifirst+1)
         dxin  = xin(ifirst+1) - xin(ifirst)
         npl   = xin(ifirst) / dxin + 1
         if(iout.gt.0) then
            write(6,*)'parabola extrapolation to x-->0 : a*x^2 + b'
            write(6,*)'a       =  ',apar
            write(6,*)'b       =  ',bpar
            write(6,*)'+points =  ',npl
         endif
         do i=1,npl
           xxin      = (i-1)*dxin
           yinter(i) = apar*xxin**2 + bpar
         enddo
c
c -- high-q-extrapolation --
c
         slope = (yin(ilast)-yin(ilast-2))/
     *           (xin(ilast)-xin(ilast-2))
         zparhq= xin(ilast-1)*slope / yin(ilast-1)
         aparhq= yin(ilast-1)*xin(ilast-1)**(-zparhq)
         if(iout.gt.0) then
           write(6,*)'high-q extrapolation :  a*q^z'
           write(6,*)'a =  ',aparhq
           write(6,*)'z =  ',zparhq
         endif
c
c --- fill the rest of the interpolation vector with equidistant data --
c
         ninter = xmax / dxin + 1
         if(ninter.gt.mdim) then
           write(6,*)'not enough space to store intermediate',
     *               ' interpolation vector'
           ier = 2
           return
         endif
         jl = ifirst+1
         do i=npl+1,ninter
           xxin = (i-1)*dxin
           do j=jl,ilast
             if(xin(j).gt.xxin) goto 10
           enddo
10         continue
           jl = j
           if(jl.gt.ilast) then
             yinter(i) = aparhq * xxin**zparhq
           else
             p = (xxin-xin(jl-1))/(xin(jl)-xin(jl-1))
             yinter(i) = yin(jl)*p + yin(jl-1)*(1.0-p)
           endif
         enddo
cdebug
         call pushda(yinter,0.0,dxin,ninter)
cdebug
c
c ----- prepare the 1d-data for 2d-fft -------------------------------
       dx = xmax / (nfft-1)
c period of fft is:
       xl = 2 * nfft * dx
c reciprocal fourier increment:
       dy = 1.0 / xl
c ---- generate the 2d field of y-values ----
       if(iout.gt.1)write(6,*)'generating data field to be transformed'
c ---- init 2d-intermediate data field for fft :
       do i=1,nfft
        x1 = (i-1) * dx
        do j=1,nfft
         x2 = (j-1) * dx
         xx = sqrt(x1**2 + x2**2)
         nn = xx/dxin + 1
         if(nn.lt.ninter) then
           p  = xx/dxin - (nn-1)
           yy = yinter(nn+1)*p + yinter(nn)*(1.0-p)
           ca(i,j) = yy * dx * dx
         else
           ca(i,j) = 0.0
         endif
        enddo
       enddo
c
cdebug
         do i=1,nfft
            yinter(i) = ca(1,i)/dx/dx
         enddo
         call pushda(yinter,0.0,dx,nfft)
cdebug
       if(iout.gt.1)write(6,*)'starting 1. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
cdebug
       call pushda(ca(1,1),0.0,dy,nfft)
cdebug
c
       ai0= ca(1,1)/(1.-t)
c ---  === : i0 intensity factor provided the original data do not
c            contain any primary beam components. otherwise the
c            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
c ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
c
c removing of multiple scattered events and transmission correction
       do i=1,nfft
        do j=1,nfft
         ca(i,j) = alog( ca(i,j)/ai0 + t) - at
c ----                                 = <- no primary beam int.
c                                           will be in the results !
c                                           omitt this t if primary
c                                           beam is in the data.
        enddo
       enddo
c
cdebug
       call pushda(ca(1,1),0.0,dy,nfft)
cdebug
c ---- reverse transformation back to scattering angle field ----
       if(iout.gt.1)write(6,*)'starting 2. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
c
c ---- return values:
       do i=1,nfft
        xout(i) = (i-1) * dx
        yout(i) = ca(1,i) * dy * dy * ai0
       enddo
       if(iout.gt.0) then
          write(6,*)'intensity factor i0 =',ai0
          write(6,*)'                     ','=================='
       endif
c ----------------------------------------------------------------------
c
       return
       end
c*ds
c*ed
       subroutine mux(xin,yin,nin,t,xmax,nfft,iout,xout,yout,ier)
c      ==========================================================
c  --- small angle multiple scattering evaluator ----
c      input:
c      xin(1..nin) ... : inputvector containig scatteringangles
c                        or equivalent in ascending order
c      yin(1..nin) ... : inputvector containing scattering intensities
c      nin ........... : active length of xin and yin
c      t ............. : effective transmission factor associated
c                        with the small angle scattering
c      xmax .......... : x (scattering vector/angle) range
c                        for the result
c      nfft .......... : number of points of the internal
c                        fouriertransforms = the number of
c                        points in the output vectors.
c                        preferred value: 2**n.
c      iout .......... : variable to control additional printoutput
c                        from this subroutine
c      output:
c      xout(1..nfft).. : output x-values
c      yout(1..nfft).. : output y-values
c      ier ........... : errorindicator (0=ok)
c
c      the output will be the hypothetical scattering data of the
c      same sample using the same apparatus
c      but with any effects of multiple scattering,
c      e.g. no multiple sacttered intensity and transmission
c      reduction.
c ----------------------------------------------------------------------
c
c ---- dimensioning ----------------------------------------------------
       parameter(pi=3.1415926535897)
       parameter(mdim=1024)
       real xin (nin) ,yin (nin)
       real xout(nfft),yout(nfft)
       common/fftwr1/  ca(mdim,mdim), yinter(mdim)     !! k--> 1 aix !!
c ----------------------------------------------------------------------
       ier = 0
       if(nfft.eq.0)       nfft = mdim
       if(nfft.gt.mdim)    nfft = mdim
 
       if(iout.gt.0) then
         write(6,*)'mux uses the following parameters:'
         write(6,*)'nfft ........... :',nfft
         write(6,*)'xmax ........... :',xmax
         write(6,*)'transmission ... :',t
       endif
 
       if(t.le.0. .or. t.ge.1. ) then
         write(6,*)'mux: transmission value is out of range !'
         ier = 1
         return
       endif
 
c------- data preparation ----------------------------------------------
c
c --- fill the rest of the interpolation vector with equidistant data --
c
         dx = xmax / (nfft-1)
         ninter = xmax / dx + 1
         if(ninter.gt.mdim) then
           write(6,*)'not enough space to store intermediate',
     *               ' interpolation vector'
           ier = 2
           return
         endif
         jl = 2
         do i=1,ninter
           xxin = (i-1)*dx
           do j=jl,nin
             if(xin(j).gt.xxin) goto 10
           enddo
10         continue
           jl = j
           if(jl.gt.nin) then
             yinter(i) = 0.0
           else
             p = (xxin-xin(jl-1))/(xin(jl)-xin(jl-1))
             yinter(i) = yin(jl)*p + yin(jl-1)*(1.0-p)
           endif
         enddo
cdebug
         call pushda(yinter,0.0,dxin,ninter)
cdebug
c
c ----- prepare the 1d-data for 2d-fft -------------------------------
       dx = xmax / (nfft-1)
c period of fft is:
       xl = 2 * nfft * dx
c reciprocal fourier increment:
       dy = 1.0 / xl
c ---- generate the 2d field of y-values ----
       if(iout.gt.1)write(6,*)'generating data field to be transformed'
c ---- init data 2d-intermediate data field for fft :
       do i=1,nfft
        x1 = (i-1) * dx
        do j=1,nfft
         x2 = (j-1) * dx
         xx = sqrt(x1**2 + x2**2)
         nn = xx/dx  + 1
         if(nn.lt.ninter) then
           p  = xx/dx   - (nn-1)
           yy = yinter(nn+1)*p + yinter(nn)*(1.0-p)
           ca(i,j) = yy * dx * dx
         else
           ca(i,j) = 0.0
         endif
        enddo
       enddo
c
cdebug
         do i=1,nfft
            yinter(i) = ca(1,i)/dx/dx
         enddo
         call pushda(yinter,0.0,dx,nfft)
cdebug
       if(iout.gt.1)write(6,*)'starting 1. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
cdebug
       call pushda(ca(1,1),0.0,dy,nfft)
cdebug
c
       ai1= ca(1,1)
       at = alog(t)
       xn = -at / ai1
c
c multiple convolution series:
       do i=1,nfft
        do j=1,nfft
         ca(i,j) = exp(xn*ca(i,j)+at)- t
c ----                                 = <- no primary beam int.
c                                           will be in the results !
c                                           omitt this t if primary
c                                           beam is to be generated.
        enddo
       enddo
cdebug
       call pushda(ca(1,1),0.0,dy,nfft)
cdebug
c
c ---- reverse transformation back to scattering angle field ----
       if(iout.gt.1)write(6,*)'starting 2. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
c
c ---- return values:
       do i=1,nfft
        xout(i) = (i-1) * dx
        yout(i) = ca(1,i) * dy * dy / xn
       enddo
c ----------------------------------------------------------------------
c
       return
       end
c
c
       subroutine pushda(y,x0,dx,n)
c      ----------------------------
c --> pushes data y-values to the datafiles
c
       parameter(mwert=1024,mbuf=200,mpar=200)
c --- mwert  = max. no. of x-y-values in one buffer
c     mbuf   = max. no. of different buffers
c     mpar   = max. no. of parameters associated with one buffer
c ---  maximum scan length ....
       parameter(mth=40,mtpar=40,mtcal=40)
c ---  fit dimensions ---
       parameter (mfit=40,msmpl=4000)
c  -- mfit = max no. of fitted parameters
c     msmpl= max no. of datapoints in fit
c
c --- incom common-section ---
       parameter(minc=40)
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
      logical cray
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c --- xwerte(i,j)   x-values on buffer j
c     ywerte(i,j)   y-values on buffer j
c     yerror(i,j)   error of y-values   (only supported by some fktn)
c     xname(j)      name of x-values (x-axis) for buffer j
c     yname(j)       "   "  y-   "    y-  "    "    "    "
c     name(j)       short text-identifier for data on buffer j
c     nwert(j)      no. of valid points on buffer j
c     nbuf          no. of filled buffers
c     numor(j)      numerical idenfication of data on buffer j
c     coment(j)     one line of comment describing data on buffer j
c     params(l,j)   set of parameters associated with data on buffer j
c     napar(l,j)    names of these parameters
c     nopar(j)      no. of valid parameters
c
c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c  isels(i)  = address (/spectr/) of selected scan
c  ifits(i)  = adress of fitted spectrum (isels(i))
c  nsel      = length of this table
c
       common/fslist/isfits(mbuf),nfsel
c  isfits    = address (/spectr/) of selected fits
c  nfsel     = length of this table
c
       dimension y(n)
       data  num/9000/
c ----------------------------------------------------------------
       if(nbuf.ge.mbuf) then
         write(6,*)'pushda unsuccessful: no room on stack !'
         ierrs = 100
         return
       endif
 
       num  = num  + 1
       nbuf = nbuf + 1
       nwert(nbuf) = n
       numor(nbuf) = num
       name(nbuf)  = 'check'
       xname(nbuf) = 'x'
       yname(nbuf) = 'y'
       coment(nbuf)= 'pushda data, checking purpose only!'
       nopar(nbuf) = 0
       do i=1,n
          xwerte(i,nbuf) = x0 + (i-1)*dx
          ywerte(i,nbuf) = y(i)
          yerror(i,nbuf) = 0.0
       enddo
       write(6,*)'pushda: created item:',nbuf,' with numor:',num
       return
       end
c*ds
c*ed
      subroutine usrfun(nam,x,nx,ier)
c     -------------------------------
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       parameter (mfit=40,msmpl=4000)
       common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
       character*8 thenam,thparn
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),
     *  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos
     * ,multflg(mtcal)
c ------ errors of fit --------------------------------------------
       common/therrc/therro(mtpar,mtcal)

c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
       common/fslist/isfits(mbuf),nfsel
c
       common/cfc/qziel,cscoef(4,mwert),break(mwert),weight(mwert),
     *            numspl,nwspl
c --- parameters of last spline smoothing + value qziel
c     qziel is the value at witch the spline should be evaluated
c     numspl numor of spline fitted data
c     nwspl  length of splined data vectors
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c ---- communication with subr. func ---
      logical sqwght,sqwbuf
      logical autox1,autox2
      common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)
c  --- iprt     : printing during func calculation (set by fit)
c      x1       : lower limit of fit range (if autox1 = .false.)
c      x2       : upper limit of fit range (if autox2 = .false.)
c      autox1/2 : if true the genuine limits of the datafiled is taken
c                 else x1 /x2 values are taken (set by thc)
c ----------------------------------------------------------------------
      character*20 nam
      character*8 pname
      real*8 x(*)
      logical compare
 
      ier = 0
      if(compare(nam,'myfun3 ')) then
        if(nx.lt.3) then
          ier =1
          return
        endif
        x1   = x(nx-2)
        x2   = x(nx-1)
        x3   = x(nx)
        x(nx-2) = myfun(x1,x2,x3)
        nx = nx-2
        return
      endif
 
 
      if(compare(nam,'p  ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        iaddr= x(nx)  +0.01
        x(nx) = ptxf(iaddr)
        return
      endif
 
      if(compare(nam,'yv ')) then
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf = x(nx-1)+0.01
        iwert= x(nx)  +0.01
        x(nx-1) = ywerte(iwert,ibuf)
        nx = nx-1
        return
      endif
 
      if(compare(nam,'xv ')) then
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf = x(nx-1)+0.01
        iwert= x(nx)  +0.01
        x(nx-1) = xwerte(iwert,ibuf)
        nx = nx-1
        return
      endif
 
      if(compare(nam,'ye ')) then
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf = x(nx-1)+0.01
        iwert= x(nx)  +0.01
        x(nx-1) = yerror(iwert,ibuf)
        nx = nx-1
        return
      endif
 
      if(compare(nam,'sc ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        inumr= x(nx)*1.000001
        iadr = 0
        do l=1,nbuf
          if(numor(l).eq.inumr) then
            iadr = l
            goto 22
          endif
        enddo
22      continue
        x(nx) = iadr+0.0000001
        return
      endif
 
      if(compare(nam,'num ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        inumr= x(nx)+0.001
        iadr = numor(inumr)
        x(nx) = iadr*1.00000001
        return
      endif
 
      if(compare(nam,'sumy ')) then
        if(nx.lt.3) then
          ier =1
          return
        endif
        ibuf  = x(nx-2)+0.01
        iwert = x(nx-1)+0.01
        iwert2= x(nx  )+0.01
        sum = 0
        do i=iwert,iwert2
          sum = sum + ywerte(i,ibuf)
        enddo
        x(nx-2) = sum
        nx = nx-2
        return
      endif
 
      if(compare(nam,'sumx ')) then
        if(nx.lt.3) then
          ier =1
          return
        endif
        ibuf  = x(nx-2)+0.01
        iwert = x(nx-1)+0.01
        iwert2= x(nx  )+0.01
        sum = 0
        do i=iwert,iwert2
          sum = sum + xwerte(i,ibuf)
        enddo
        x(nx-2) = sum
        nx = nx-2
        return
      endif

 
      if(compare(nam,'th_par ')) then
        if(nx.lt.2) then
          ier =1
          return
        endif
        iparn = x(nx-1)+0.01
        itheo = x(nx  )+0.01
        x(nx-1) = thparx(iparn,itheo)
        nx = nx-1
        return
      endif
 
 
      if(compare(nam,'th_err ')) then
        if(nx.lt.2) then
          ier =1
          return
        endif
        iparn = x(nx-1)+0.01
        itheo = x(nx  )+0.01
        x(nx-1) = therro(iparn,itheo)
        nx = nx-1
        return
      endif
 
 
      pname = ' '
      do i=1,8
       if(nam(i:i).eq.' ') goto 111
       pname(i:i) = nam(i:i)
      enddo
111   continue
      iiadr = x(nx) + 0.001
      call parget (pname,pvalue,iiadr,ier)
      x(nx) = pvalue
 
      return
      end
c*ds
c*ed
      subroutine usrextr(nam,val,ier)
c
       parameter(mwert=1024,mbuf=200,mpar=200)
       parameter(mth=40,mtpar=40,mtcal=40)
       parameter (mfit=40,msmpl=4000)
       common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
       character*80 name,xname,yname,napar,coment*80
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),
     *        yerror(mwert,mbuf),
     *        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),
     *        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),
     *        napar(mpar,mbuf),nopar(mbuf)
c
c ---- common containing a selected list of spectra ----
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
c
       common/fslist/isfits(mbuf),nfsel
c
       common/cfc/qziel,cscoef(4,mwert),break(mwert),weight(mwert),
     *            numspl,nwspl
c --- parameters of last spline smoothing + value qziel
c     qziel is the value at witch the spline should be evaluated
c     numspl numor of spline fitted data
c     nwspl  length of splined data vectors
c
c ---- outputlevel
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
c
c ---- communication with subr. func ---
      logical sqwght,sqwbuf
      logical autox1,autox2
      common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl)
c  --- iprt     : printing during func calculation (set by fit)
c      x1       : lower limit of fit range (if autox1 = .false.)
c      x2       : upper limit of fit range (if autox2 = .false.)
c      autox1/2 : if true the genuine limits of the datafiled is taken
c                 else x1 /x2 values are taken (set by thc)
c ----------------------------------------------------------------------
      character*1 nam(*)
      character*8 pname
      logical compare
      real*8 val
      ier = 0
      if(compare(nam,'sel ')) then
        val = isels(1)
        return
      endif
      if(compare(nam,'maxx ')) then
        ibuf = isels(1)
        val = xwerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(xwerte(i,ibuf).gt.val) val = xwerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'minx ')) then
        ibuf = isels(1)
        val = xwerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(xwerte(i,ibuf).lt.val) val = xwerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'maxy ')) then
        ibuf = isels(1)
        val = ywerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(ywerte(i,ibuf).gt.val) val = ywerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'miny ')) then
        ibuf = isels(1)
        val = ywerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(ywerte(i,ibuf).lt.val) val = ywerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'yy   ')) then
        val = yyyy
        return
      endif
      if(compare(nam,'xx   ')) then
        val = xxxx
        return
      endif
      if(compare(nam,'iout ')) then
        val = iout
        return
      endif
      if(compare(nam,'nbuf ')) then
        val = nbuf
        return
      endif
      pname = ' '
      do i=1,8
       if(nam(i).eq.' ') goto 111
       pname(i:i) = nam(i)
      enddo
111   continue
      
      call parget (pname,pvalue,isels(1),ier)
      
        if(ier.ne.0) then
         if(ifits(1).gt.0.and.nsel.gt.0) then
          call parget (pname,pvalue,ifits(1),ier)
         endif
        endif

      val = pvalue
      return
      end
 
 
       function wlambda( alam )
c      ------------------------
c --- repraesentiert die wellenlaengenverteilung
       implicit real*4 (a-h,o-z)
       common/wlntran/alam0, dalam
          arg     = ( (alam-alam0)/dalam )**2
          if(arg.lt.50.e0) then
            wlambda =  exp( -arg )
          else
            wlambda = 0.e0
          endif
       return
       end
c*ds
c*ed
       function sofqom( omega )
c      ------------------------
c --- repraesentiert die streufunktion, omega in s**-1
       implicit real*4 (a-h,o-z)
       common/sqtran/tau
c      tau in sekunden
 
       x      = omega * tau
       sofqom = tau /(1.e0+x*x)
 
       return
       end
c*ds
c*ds
       function f(n,x)
c      ---------------
c --- integrand ---
       implicit real*4 (a-h,o-z)
       real*4 j1echo,j2echo,j0delta
       common/partran/j1echo,j2echo,j0delta,cdelta
       dimension x(n)
c
c      x(1) = lambda
c      x(2) = omega
c
       data gamma/18303.33e0/
c ---- larmorkonstante in rad/s/gauss
       data xmh  /2.50607e-4/
c ---- neutronenmasse durch h in 2/m/angstroem
       data zpi  /6.283185307e0/
 
       a  = gamma*xmh*x(1)/(1.d0+xmh*x(1)*x(1)*x(2)*1d-10/zpi)
       b  = gamma*xmh*x(1)
       del= j0delta**2 + ( (j1echo+j2echo)*0.5e0*cdelta )**2
       f  = ( 1.d0 +  exp(-(a**2+b**2)*del*0.25e0 )*
     *               cos( a*j1echo - b*j2echo )
     *      ) * sofqom(x(2)) * wlambda(x(1))
       return
       end
 
c*ds
c*ed
       function yintp2 (x,y,n,t)
c      ------------------------ einfache tabelleninterpolation ----
       implicit real*4 (a-h,o-z)
       dimension x(n),y(n)
c
       if(t.lt.x(1) .or. t.gt.x(n-1) ) then
         yintp2  = 0.e0
         return
       endif
c
       do i=2,n
         if(x(i).gt.t) then
           p = (t-x(i-1))/(x(i)-x(i-1))
           yintp2  = (1.e0-p)*y(i-1)+p*y(i)
           return
         endif
       enddo
       yintp2  = 0.e0
       return
       end
 
       function igrand(an)
c      ------------------
c ---- gaussian approximation to poisson statistics with
c      expectation value = an
c      second moment     = an
c      integer result
c      to mock up counting statistics
c ---------------------------------------------------------------------
       implicit real*4 (a-h,o-z)
1      continue
        y =  rnunf()*2.0-1.0
        z =  erfi(y)
        igrand = an + sqrt(2.0*an)*z + 0.5
        if(igrand.ge.0) return
       goto 1
       end
      
