       include 'GAmodule.f90'


       subroutine lsearch( jpar, itcal, ierr)
!      ======================================
       use outlev
       use theory
       implicit none
!  ---- search for a theory parameter specification by the commandline
!       cmd  theoryname <n-th occ> parametername
!  output: jpar = parameter-adress
!          itcal= theory adress
!          ierr = errorindicator ( 0=ok  1=not found)
!

       integer jpar, itcal, ierr
       integer i, j, ict, icz, ith, np
       character*8      vnameF
       double precision rparF
       integer          inpaF

! ----------------------------------------------------------------------
! --- which occurence of the theory
       ict = 0
       if(inpaF(1).eq.0) then
          icz = 1
       else
          icz = rparF(1) + 0.1
       endif
!
       do 1 i=1,ntheos
          ith = nthtab(i)
          if(vnameF(1).eq.thenam(ith)) then
            ict = ict + 1
            if(ict.eq.icz) goto 2
!                          ------> ok we have it
          endif
    1  continue
          write(6,*)'search: theory=',vnameF(1),' could not be found'
          write(6,*)'check spelling or if it is in theos directory'
          ierr = 1
          ierrs = 801
          return
    2  continue
       itcal= i
!      ========
! ---- ok theory is identified. look now for parameter.....
       np = nthpar(ith)
       do 4 j=1,np
          if(vnameF(2).eq.thparn(j,ith)) goto 6
    4  continue
         write(6,*)'search: parameter=',vnameF(2),                      &
     &             ' of theory=',vnameF(1),                             &
     &             ' could not be found'
         ierr = 1
         ierrs = 802
         return
    6  continue
       ierr = 0
       jpar = j
!      ========
!
       return
      END


      subroutine search(irecv,nkurv)
!     ------------------------------------
!
! --- the adresses(/spectr/) of the nkurv  numors given in irecv are
!     located and put to irecn1
!     nkurv is set to the number of actually found items
!  !  selection table is updated
!     imx returns the no. of biggest gradient contribution
!
!
! ----- search for the scans to be selected -----

       use constants
       use cdata
       use outlev
       use selist
       implicit none

      integer irecv,nkurv
      integer i, ic, inum1, inum2, irecn1, l, num
      dimension irecv(minc),irecn1(minc)
!
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
 1002     continue
        if(iout.ge.0) write(6,*)'scan: ',num,' has not been found!'
 1001  continue
!
       nkurv = l
!      ---------> only found scans may be selected ...
       if(nkurv.le.0) then
        if(iout.ge.0) write(6,*)'no scan has been found'
        return
       endif
!
!
! --- check for the other scans ---
       if(nkurv.ge.2) then
         inum1 = irecn1(1)
         do 1020 i=2,nkurv
           inum2 = irecn1(i)
           if(xname(inum2).ne.xname(inum1)) then
             write(6,*)'warning abscissas may be incompatible:'
             write(6,*)name(inum1),' x: ',xname(inum1),                 &
     &                 name(inum2),' x: ',xname(inum2)
           endif
 1020    continue
       endif
!      ------(of if(nkurv.ge.2)....)
!
! --- update selection table ---
       do 2000 i=1,nkurv
        isels(i) = irecn1(i)
        ifits(i) = 0
 2000  continue
       nsel = nkurv
       if(iout.ge.0) write(6,*)'selection table updated. ',l,' items'
!
       return
      END


      subroutine fsrch(ifrec,nfkurv)
!     ------------------------------------
! this subroutine is analogous to search, it will select the scans of
! the fit curves
! attention: the fsc-command requires positive inputs, therefore this
!            routine searches for -num !!
!
       use cdata
       use outlev
       use fslist
       use constants
       implicit none

      integer ifrec,nfkurv
      integer i, ic, inum1, inum2, irecn1, l, nkurv, num
      dimension ifrec(minc),irecn1(minc)
!
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
 1002     continue
        if(iout.ge.0) write(6,*)'scan: ',num,' has not been found!'
 1001  continue
!
       nkurv = l
!      ---------> only found scans may be selected ...
       if(nfkurv.le.0) then
        if(iout.ge.0) write(6,*)'no scan has been found'
        return
       endif
!
!
! --- check for the other scans ---
       if(nfkurv.ge.2) then
         inum1 = irecn1(1)
         do 1020 i=2,nfkurv
           inum2 = irecn1(i)
           if(xname(inum2).ne.xname(inum1)) then
             write(6,*)'warning abscissas may be incompatible:'
             write(6,*)name(inum1),' x: ',xname(inum1),                 &
     &                 name(inum2),' x: ',xname(inum2)
           endif
 1020    continue
       endif
!      ------(of if(nfkurv.ge.2)....)
!
! --- update f-selection table ---
       do 2000 i=1,nfkurv
        isfits(i) = irecn1(i)
 2000  continue
       nfsel = nfkurv
       if(iout.ge.0) write(6,*)'selection table updated. ',l,' items'
!
       return
      END
!
!


       real function thval(x)
!      ----------------------
!
! ---- compute the value of the active set of theories at value x 
!
       use cdata
       use theory
       use thparc
       use constants
       implicit none

       real x
       integer iadda
       common /thiadd/iadda

       integer sum
       real theory_x

       sum = 0
!      here we call the function value with all neccessary to calc it   it only needs params from current dataset not from all
!       this will make it easier for future releases to incorporate C without all these commonblocks RB
       thval =  theory_x(x,ntheos,nthtab,multflg,thparx,thrapar,thramin,thramax,nopar(iadda),params(:,iadda),napar(:,iadda),mbuf,mth,mtpar,mtcal)
       return
      END


       subroutine fit
!      ==============
!
! ---- fit of activated theories to data on /fil2/ ----
!
       use cincom
       use cincoc
       use cdata
       use outlev
       use theory
       use selist
       use therrc
       use cfunc
       use cfunce
       use constants

       implicit none
       integer i, it, icode, ier, igo, inew, iparam, ip, iru,l, map
       integer iercd, ll, j, j1, j2, j3, ith, ixjac
       integer npar, numv, nfit, nff, n, m, nspf
       real thperr, sum, rparam
       real x, xscale, xjac, xcenter, xguess, xerr, xstepsc
       real divis, ermat, f, gmat, fscale, ginv
       real*8 ssq

       character*8 ci
       real*8 getval
      dimension iparam(6),rparam(7),x(mfit),f(msmpl),xjac(msmpl,mfit),  &
     &          xguess(mfit),xscale(mfit),fscale(msmpl)
      dimension numv(minc)
      dimension ermat(mfit,msmpl), gmat(mfit,mfit), ginv(mfit,mfit)

      dimension   xcenter(mfit), xstepsc(mfit)
      character*8 xpname(mfit)

      external func

      logical found
      logical :: final_thc=.false., sqwbuf=.false.

! ---- defaults for parameters -----
       integer :: maxfn = 500, ngood = 0, maxit = 0, iprint = 1, irecse = 0
       real :: stpsz = 0.0, trure = 0.0

! ---- take parameters from stack ----
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
       if(vname(i).eq.'scans   '.or.vname(i).eq.'sc      ') then
         nspf = inpar(i)
         do 7001 l=1,nspf
          numv(l) = rpar(j+l-1) * 1.0000001
 7001    continue
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
! ----> the parameters  x1 x2 resp. the option autox1 autox2
!       are to be decoded during the first call of thc by thc !!
 1011 continue
       sqwbuf= sqwght
! --- search scans & prepare selection table
       if(nsel.eq.0) call search(numv,nspf)
!      -----------------------------------------
!
! ---- write all parameters ----
      if(iout.ge.0) then
        write(6,1)iprint,maxfn,maxit,ngood,                             &
     &           (numor(isels(i)),i=1,nsel)
    1   format(' ***** fit-profile in effect ***** '/                   &
     &  ' printing    ',i8/                                             &
     &  ' max funcalls',i8,/                                            &
     &  ' max iterat .',i8,/                                            &
     &  ' n-good-dig .',i8,/                                            &
     &  ' scans :'/(1x,10i7))
        if(sqwght)write(6,*)' ***** weight differences by 1/sqrt *****'
        if(lerrel)write(6,*)' ***** errors weigthed relative:relerr **'
        if(lwrtfitdat)write(6,*)' ***** writing fitdat.tmp  **'
        if(final_thc) then
          write(6,*)' ***** opt: thc   ->    final fine mesh call**'
        else
          write(6,*)' ***** opt: nothc -> no final fine mesh call**'
        endif
       endif
!
      if(igo.eq.0) return
 1000 continue
!
          iprt = iprint
!
          write(6,*)' startparameters : '
          call activa(2)
!
! ---- prepare startvalues -----
          if(ntheos.eq.0) then
            write(6,*)' fit ==> no theories activated ...'
            return
          endif
!
          nfit = 0
          do 10 it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do 20 ip=1,npar
             if(thpsca(ip,it).ne.0) then
               nfit = nfit + 1
               if(nfit.gt.mfit) then
                 write(6,*)' too many fitparameters for current mfit=', &
     &                     mfit
                 return
               endif
               x(nfit) = thparx(ip,it) / thpsca(ip,it)
               xstepsc(nfit) =  thpsca(ip,it)
               xpname(nfit)  =  thparn(ip,ith)
             endif
   20     continue
   10   continue
!
        n     = nfit
        ixjac = msmpl
        call func(m,n,x,f)
!       -------------------> determine m & di a first test calculation
        write(6,*)' no. of sample points m = ',m
        if(.not.autox1) write(6,*)'set lower limit of x: x1=',x1
        if(.not.autox2) write(6,*)'set upper limit of x: x2=',x2
!
! ------ imsl version 10 :  setup of some new vectors ------------------
!        to meet the function of the old zxssq approximately
        do 10010 i=1,n
           xguess(i) = x(i)
           xscale(i) = 1.0
10010   continue
!
        do 10020 i=1,m
           fscale(i) = 1.0
10020   continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(found('map     ')) goto 20000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- set the desired execution-parameters ---
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
! ---------------------------------------------------------------------
! ---- do the fitting with imsl routine unlsf (old:zxssq) ----
        call unlsf(func,m,n,xguess,xscale,fscale,iparam,rparam,x,f,xjac,ixjac)
        if(iout.gt.0) write(6,*)'final   iparam=',iparam
        if(iout.gt.0) write(6,*)'final   rparam=',rparam
!
        icode = iercd()
!       ---------------> get the errorcode
! ----- output -----
        write(6,*)'icode = ',icode
        if(icode.eq.1) write(6,*)' **-1- : no longer converging   *****'
        if(icode.eq.2) write(6,*)' **-2- : conv. to a noncritical pnt *'
        if(icode.eq.3) write(6,*)' **-3- : maxfn exceeded *************'
        if(icode.eq.4) write(6,*)' **-4- : maxit exceeded *************'
        if(icode.eq.5) write(6,*)' **-5- : 5 steps with max steplength*'
!
        write(6,501)iparam(4),iparam(3),iparam(2)
  501   format(                                                         &
     &   ' no. of func calls ..... ',i6/                                &
     &   ' est. no. of sig. dig. . ',i6/                                &
     &   ' no of iterations ...... ',i6/)
!
!

! ----- error-determination ?? test ?? -----------------------------
!
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
! --- invertiere gmat ---
         call linrg(n,gmat,mfit,ginv,mfit)
! --- bilde g**-1 f  (f=xjac) ---
         do i=1,m
          do l=1,n
           sum = 0
           do j=1,n
            sum = sum + ginv(l,j)*xjac(i,j)
           enddo
           ermat(l,i) = sum
          enddo
         enddo
! --- bilde fehlerquadratmatrix ---
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
! --- extrahiere fehler und ausgabe ---
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
               write(6,*)'error for ', thenam(ith),'(',thparn(ip,ith),  &
     &                   ') = +- ',thperr
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
! ------error-determination ?? end  ?? -----------------------------

        call func(m,n,x,f)
!       ------------------> final call including elastic lines
        write(6,*)' ******* final parameters ********'
        if(irecse.ne.0) write(6,*)'       record ',iru
        call activa(2)

        call extract('ssq0    ' ,ssq,ier)
        if(ier.ne.0) ssq = 0
        do i=1,nsel
          call parset ('ssq     ',sngl(ssq),ifits(i))
        enddo
!
        call couple(1)
        if(final_thc) call thc(mwert)
!       ----------------> fitted curve with finest mesh
!
!!??    if(ier.ne.0) then
!!??      ierrs = ier
!!??      return
!!??    endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  create a map
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
          write(14,'(a,1x,i9,a,a)')name(j)(1:12),                       &
     &                   numor(j),' :',coment(j)(1:60)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      return
      END

       subroutine ga_fit
!      ==================
!
! ---- fit of activated theories to data on /fil2/ ----
!      using genetic algorithm
!
       use cincom
       use cincoc
       use cdata
       use outlev
       use theory
       use selist
       use therrc
       use cfunc
       use cfunce
       use constants
       use GAmodule

       implicit none
       integer i, it, icode, ier, igo, inew, iparam, ip, iru,l, map
       integer iercd, ll, j, j1, j2, j3, ith, ixjac
       integer npar, numv, nfit, nff, n, m, nspf
       real thperr, sum, rparam
   
       real*8 ssq

       character*8 ci
       real*8 getval

       real       ::   x(mfit)
       real       ::   xguess(mfit),xscale(mfit),f(msmpl)
       real       ::   xcenter(mfit), xstepsc(mfit)
       dimension numv(minc)

       integer, save    :: npop = 30          ! default population size
       integer, save    :: ngen = 50          ! default number of generations
       real,    save    :: mutrate = 0.002    ! default mutation rate
       integer, save    :: ncodeBits = 12     ! default bitcoding

       real             :: best_fitness
      

      character*8 xpname(mfit)

      external func

      logical found
      logical :: final_thc=.false., sqwbuf=.false.

! ---- defaults for parameters -----
       integer :: maxfn = 500, ngood = 0, maxit = 0, iprint = 1, irecse = 0
       real :: stpsz = 0.0, trure = 0.0



      GAtrace = 0                            ! default is restored, must be given explicit 'trace' 
                                             ! for each call

! ---- take parameters from stack ----
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
       if(vname(i).eq.'scans   '.or.vname(i).eq.'sc      ') then
         nspf = inpar(i)
         do 7001 l=1,nspf
          numv(l) = rpar(j+l-1) * 1.0000001
 7001    continue
       endif
       if(ci.eq.'npop    '                    ) npop      = nint(rpar(j))
       if(ci.eq.'ngen    '                    ) ngen      = nint(rpar(j))
       if(ci.eq.'mutation'                    ) mutrate   = rpar(j) 
       if(ci.eq.'trace   '                    ) GAtrace   = nint(rpar(j))
       if(ci.eq.'bits    '                    ) ncodeBits = nint(rpar(j))
       if(ci.eq.'wrtfit  '                    ) lwrtfitdat = .true.
       if(ci.eq.'nowrtfit'                    ) lwrtfitdat = .false.
! ----> the parameters  x1 x2 resp. the option autox1 autox2
!       are to be decoded during the first call of thc by thc !!
 1011 continue
       sqwbuf= sqwght
! --- search scans & prepare selection table
       if(nsel.eq.0) call search(numv,nspf)
!      -----------------------------------------
!
! ---- write all parameters ----
      if(iout.ge.0) then
        write(6,1)iprint,npop,ngen,ncodeBits,mutrate,                   &
     &           (numor(isels(i)),i=1,nsel)
    1   format(' ***** fit-profile in effect ***** '/                   &
     &  ' printing              ',i8/                                   &
     &  ' population    npop    ',i8,/                                  &
     &  ' generations   ngen    ',i8,/                                  &
     &  ' n-code-bits   bits    ',i8,/                                  &
     &  ' mutationrate  mutation',F12.9,/                               &
     &  ' scans :'/(1x,10i7))
        if(sqwght)write(6,*)' ***** weight differences by 1/sqrt *****'
        if(lerrel)write(6,*)' ***** errors weigthed relative:relerr **'
        if(lwrtfitdat)write(6,*)' ***** writing fitdat.tmp  **'
        if(final_thc) then
          write(6,*)' ***** opt: thc   ->    final fine mesh call**'
        else
          write(6,*)' ***** opt: nothc -> no final fine mesh call**'
        endif
       endif
!
      if(igo.eq.0) return
 1000 continue
!
          iprt = iout
!
          write(6,*)' startparameters : '
          call activa(2)
!
! ---- prepare startvalues -----
          if(ntheos.eq.0) then
            write(6,*)' fit ==> no theories activated ...'
            return
          endif
!
          nfit = 0
          do 10 it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do 20 ip=1,npar
             if(thpsca(ip,it).ne.0) then
               nfit = nfit + 1
               if(nfit.gt.mfit) then
                 write(6,*)' too many fitparameters for current mfit=', &
     &                     mfit
                 return
               endif
               x(nfit) = thparx(ip,it) / thpsca(ip,it)
               xstepsc(nfit) =  thpsca(ip,it)
               xpname(nfit)  =  thparn(ip,ith)
             endif
   20     continue
   10   continue
!
        n     = nfit
        ixjac = msmpl
        call func(m,n,x,f)
!       -------------------> determine m & di a first test calculation
        write(6,*)' no. of sample points m = ',m
        if(.not.autox1) write(6,*)'set lower limit of x: x1=',x1
        if(.not.autox2) write(6,*)'set upper limit of x: x2=',x2
!
        do 10010 i=1,n
           xguess(i) = x(i)
           xscale(i) = 1.0
10010   continue
!
! ---------------------------------------------------------------------
        call ga_evolve( xguess ,xscale, n, ncodeBits, mutrate, npop, ngen, x, best_fitness)

!
        write(6,'(a,f12.6)') 'best_fitness = ',best_fitness


        call func(m,n,x,f)
!       ------------------> final call 
        write(6,*)' ******* final parameters ********'
        if(irecse.ne.0) write(6,*)'       record ',iru
        call activa(2)

        call extract('ssq0    ' ,ssq,ier)
        if(ier.ne.0) ssq = 0
        do i=1,nsel
          call parset ('ssq     ',sngl(ssq),ifits(i))
        enddo
!
        call couple(1)
        if(final_thc) call thc(mwert)
!       ----------------> fitted curve with finest mesh
!
! 
       return
       END subroutine ga_fit


        function ga_func(x,n)
!       =====================
        use cdata
        use outlev
        use cfunc
        use cfunce
        use constants

        implicit none
        real                                :: ga_func
        real             , intent(inout)    :: x(n)
        integer          , intent(inout)    :: n


        real             ::   f(msmpl)


    
        integer          :: i,m 
        real             :: y
    

        call func(m,n,x,f)

!        write(6,*)'f: m,n,x ',m,n,x(1:n)

        y = 0
        do i=1,m
          y = y + f(i)**2
        enddo

        y = y/m 
        
        ga_func = 1/y 
    
        end function ga_func
    


       subroutine func(m,nff,x,f)
!      ==========================
!
! ----- calculates deviation data theory for zxssq ----
!
! compiler complains about nff unsused!!!!

       use cdata
       use theory
       use selist
       use cfunc
       use cfunce
       use constants

       implicit none
       integer m, nff
       real      x(mfit),f(msmpl)
       integer npar, nfit, n, mn, ith, isel, ipt, it, ip, ier, icall, iad1, iad2, i
       real xx, ssq, ferr

!
! ---- restore startvalues & parameters ----
!
          icall = 0
          nfit = 0
          do 10 it =1,ntheos
            ith = nthtab(it)
            npar= nthpar(ith)
            do 20 ip=1,npar
             if(thpsca(ip,it).ne.0) then
               nfit = nfit + 1
               thparx(ip,it) = x(nfit) * thpsca(ip,it)
             endif
   20     continue
   10   continue

! ==============================================================
!
!  initialize random number generator to the same value each time
!  in order to enable fitting with gradients
!
       call srand(75251)


!
! ----- calculate theory-values ------
!
       call couple(1)
       call thc(0)
!      -----------
!
! ----- compare with data ----
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
!       --------> counts valid values in the thc-computed data
        do 103 ipt=1,n
          xx=xwerte(ipt,iad1)
          if(((xx.ge.x1).or.autox1).and.((xx.le.x2).or.autox2)) then
! ---------> this if stmnt must correspond to the equival. in thc ! --
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
               write(6,*)' func ==> no. of compared datapts ',m,        &
     &                   ' exhausts current dimension msmpl=',msmpl
               goto 111
             endif
           endif
!          -----> of if(xx.ge...
  103    continue
         if(mn.ne.nwert(iad2)) then
           write(6,*)'func - error(programming!) exp and thc dont match'
           if(lwrtfitdat) close(90)
           return
         endif
  100  continue
  111  continue
          if(m.eq.0) then
            write(6,*)' no point in fit window !!'
            m=1
            f(1)=9999
          endif
!
       icall = icall + 1
! ---- output if option is set ---
       ssq = ssq/m
       if(iprt.gt.0) write(6,200)icall,ssq,(x(i),i=1,nfit)
  200  format(' ',i4,': ssq=',5e12.4/(23x,4e12.4))
       call setudf('ssq0 ',dble(ssq),ier)
       fcssq = ssq
!
!
       if(lwrtfitdat) then
        call theo_out(90)
        close(90)
       endif

       return
      END



       subroutine thc(npunkt)
!      ----------------------
!
! ---- compute theoretical scans with a resolution of npoint points
!      of the scans given in /selist/isels(1..nsel)
!      the addresses of the computed scans are put to /selist/ifits
!      if npoint= 0 the sample points match exactly those of the
!      corresponding isels-scan
!
       use cincom
       use cincoc
       use outlev
       use cdata
       use selist
       use cfunc
       use constants

       implicit none

       real*8 getval, dble
!
       integer :: iadda
       common /thiadd/iadda
! ---> transfer der addresse der gerade bearbeiteten kurve nach thdatr1
!      damit koennen dann ggf. weitere parameter ueber parget gewonnen
!      werden
!
       logical :: convolute = .false., found, folgt
!       real :: extup = 0. ,extlow = 0.
       integer :: iadd, i, j, ier, m, n, np, npk, npoint, intval, isigint, nconvo, isel, inew
       integer :: inbuf, im, npunkt
       real :: da, dx, sum, x, xmax, xmin, xc1, xc2, thval, xxx, xnconvo

       npoint = IABS(npunkt)
       npk    = npunkt
!
       if(found('convolut'))           convolute = .true.
       if(folgt('off     ','convolut'))convolute = .false.
!
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
 2001  continue
!
!
       if(nbuf+nsel.gt.mbuf) then
         write(6,*)'not enough space to store computed spectra!'
         ierrs = 301
         return
       endif
!
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
!

       if(isigint().gt.0) then
         write(6,*)'Thc Return withut calc due to Sigint:',isigint()
         return
       endif

! ---- compute ....
       do 1001 isel=1,nsel
        iadd = isels(isel)
        iadda= iadd
        if(ifits(isel).eq.0) then
          nbuf        = nbuf+1
          ifits(isel) =  nbuf
        endif
        inbuf = ifits(isel)
!       --------------------> where the data are stored
        if(iout.gt.3) write(6,*)'isel=',isel,'  iadd=',iadd
        if((npoint.le.1).and..not.convolute) then
!                                 ----> take old scan structure
!!!!!!!!
!! erst mal kommentiert lassen, da in diesem Fall auch parget geaendert
!! eventuell x1 und x2 sowie ggf xc1 und xc2 als standard Parameter mit
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
! ---- ! ---- this if condition must meet the corresponding if in the
!             do loop 103 of func ------ ! -----
             m = m+1
             xwerte(m,inbuf) = xxx
           endif
 9901     continue
          n = m
        else
!       ----> rescale the gradient
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
! --- convolution may require overlapping data --
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
!        --- no partial fitting allowed in this version
!            if convolute is used ---> set automatic range on ---
           autox1 = .true.
           autox2 = .true.
          endif
! ------------------------------------------------
          if(npoint.gt.mwert) npoint = mwert
          n = npoint

          if(npk.ge.0                                                   &
     &      .or. xmin.le.0.0d0                                          &
     &      .or. xmax .le. 0.0d0                                        &
     &      .or. xmax.le.xmin     ) then
           dx = (xmax-xmin)/(n-1)
           do im=1,n
            xwerte(im,inbuf) = xmin + dx * (im-1)
            enddo
                  ! a neg. n flags that thc computes a even intervals on
          else
           dx = (xmax-xmin)/(n-1)
           da = alog(xmax/xmin)/(n-1)
           do im=1,n
            xwerte(im,inbuf) = xmin*exp(da*(im-1))
            enddo
          endif
        endif

        nwert(inbuf) = n
! --- compute now...

          do 1003 i=1,n
             x = xwerte(i,inbuf)
!

            sum = thval(x)

            xwerte(i,inbuf) = x
            ywerte(i,inbuf) = sum
            if(iout.gt.5) write(6,*)'i: ',i,'  ywerte=',sum

 1003    continue

! --- transfer the scanparameters ---
         call txfpar(iadd,inbuf)
         numor(inbuf) =-numor(iadd)
         name(inbuf)  = 'fit'//name(iadd)(1:5)
         np           = nopar(iadd) + 2
         call parset('x1      ',x1,inbuf)
         call parset('x2      ',x2,inbuf)
! ---- vorlaeufig !!! -----
         call parset('x1      ',x1,iadd)
         call parset('x2      ',x2,iadd)
! -------------------------
         if(iout.gt.3)write(6,*)'nbuf=',inbuf,'  numor=',numor(inbuf)
         nwert(inbuf) = n
! ----- convolution ------ ( subroutine has to be supplied ! ) ------
         if(convolute) then
           call datconv(inbuf,xwerte(1,inbuf),ywerte(1,inbuf),n,        &
     &                  xwerte(1,iadd) ,nwert(iadd))
           nwert(inbuf) = n
         endif
! -------------------------------------------------------------------
!
 1001  continue

       return
      END



       subroutine parset (pname,pvalue,iadd)
!      =====================================
!
! ---- this routine changes or adds a parameter ----
!
       use cdata
       use constants
       implicit none

       character*8 pname
       real pvalue

       integer i,np, iadd

       np = nopar(iadd)
       do 100 i=1,np
         if (napar(i,iadd).eq.pname) then
           params(i,iadd) = pvalue
           return
         endif
  100  continue
       if (np.ge.mpar) then
         write(6,*)'impossible to add ',pname,' to parameterlist of '   &
     &      ,numor(iadd)
         return
       endif
!
       np = np + 1
       nopar(iadd) = np
       params(np,iadd) = pvalue
       napar(np,iadd) = pname
       return
!
      END  subroutine parset



       subroutine parscale (pname,svalue,iadd)
!      =======================================
!
! ---- this routine scales a parameter ----
!
       use cdata
       use constants
       implicit none

       character*8 pname
       real svalue

       integer i,np, iadd

       np = nopar(iadd)
       do 100 i=1,np
         if (napar(i,iadd).eq.pname) then
           params(i,iadd) = params(i,iadd) * svalue
           return
         endif
  100  continue

       return
!
      END subroutine parscale





       subroutine parget (pname,pvalue,iadd,ier)
!      =========================================
!
! ---- this routine gets the value of a parameter ----
!
       use cdata
       implicit none

       character*8 pname
       real pvalue
       integer iadd,ier

       integer i, np, iout
!
! ----------- look for the parameter -----------------------------------
       if (iadd .lt. 1) then
          write(*,*)'parget: invalid argument'
	  ier = 1
          return
       end if
       np = nopar(iadd)
       do 100 i=1,np
         if (napar(i,iadd).eq.pname) then
           pvalue = params(i,iadd)
           ier = 0
           return
         endif
  100  continue
! ------------- if not found .. ---------------------------------------
         if(ier.ge.0) then
         if(iout().gt.0)  write(6,*)'impossible to find ',pname,' in parameterlist of ',numor(iadd)
          ier = 1
         endif
         return
!
!
      END subroutine parget


       subroutine get_full_xvec (iadd,n,xv)
!      ===================================
!
! ---- this routine gets the value of a parameter ----
!
       use cdata
       use outlev  
       implicit none

       integer, intent(in)   :: iadd
       integer, intent(inout):: n
       real,    intent(out)  :: xv(*)

       if( iadd < 1 .or. iadd > nbuf) then
          write(6,*)'get_full_xvec: iadd = ',iadd,' out of range !'
          n = 0
          ierrs = 800
          return
       endif     
       if(nwert(iadd) > n) then
          write(6,*)'get_full_xvec: too many tau points',n,nwert(iadd)
          n = 0
          ierrs = 800
          return
       endif

       n = nwert(iadd)
       
       if(n <= 0) return
       xv(1:n) = xwerte(1:n,iadd)
!
      END subroutine get_full_xvec



       subroutine par_gaint_scale(svalue, iadd)
!      ----------------------------------------
!
! scale gauss-convolution parameters intensity values
! to cope with data scaing 
!
       use cdata
       implicit none

       real    :: svalue, xsc
       integer :: iadd, ier

       call parscale('ga1inten',svalue,iadd)
       call parscale('ga2inten',svalue,iadd)
       call parscale('ga3inten',svalue,iadd)
       call parscale('ga4inten',svalue,iadd)
       call parscale('ga5inten',svalue,iadd)
       call parscale('ga6inten',svalue,iadd)
       call parscale('ga7inten',svalue,iadd)
       call parscale('ga8inten',svalue,iadd)
       call parscale('ga9inten',svalue,iadd)
! we assume 9 is the maximum of gaussians !
       xsc = 1.0
       call parget('_gai_scal',xsc   ,iadd,ier)
       xsc = xsc * svalue
       call parset('_gai_scal',xsc   ,iadd)
              
       return

       end subroutine par_gaint_scale






       subroutine outputparams
!      =========================================
!
! ---- this outputs all parameters ----
!
       use cdata
       implicit none

       integer iadda
       common/thiadd/iadda

       integer i, np
!
! ----------- look for the parameter -----------------------------------
       np = nopar(iadda)
       do 100 i=1,np
           write(*,*) napar(i,iadda), params(i,iadda)
  100  continue

      END

       real*4 function dparam (pname)
!      ==============================
!
! ---- this routine gets the value of a parameter ----
!      besser gekapselte Routine zur Kommunikation mit th-Routinen
!
       use cdata
       implicit none

       integer iadda
       common /thiadd/iadda

       character*8 pname

       integer i, np, ier

! ----------- look for the parameter -----------------------------------
       np = nopar(iadda)
       do i=1,np
         if (napar(i,iadda).eq.pname) then
           dparam = params(i,iadda)
           ier = 0
           return
         endif
       enddo
! ------------- if not found .. ---------------------------------------
       dparam = 0.0

       return

      END


       subroutine fpaget(pname,value,nothe,ier)
!      ========================================

! ----- get value of the parameter pname of the nothe-th theory -----

       use theory
       implicit none

       character*8 pname
       real value
       integer nothe,ier

       integer i, npar, nth

       if(nothe.gt.ntheos) then
         write(6,*)'only ',ntheos,' theories activated !!'
         ier = 1
         return
       endif
!
       nth  = nthtab(nothe)
       npar = nthpar(nth)
       do 100 i=1,npar
        if(thparn(i,nth).eq.pname) then
           value = thparx(i,nothe)
           ier   = 0
           return
        endif
  100  continue
! ------------------- not found ----------------------------------------
       write(6,*)'theory parameter: ',pname,' (',nothe,' ) not found'
       ier = 1
       return
      END



       subroutine activa(iopt)
!      =======================
! --- activate a theory
!     iopt = 0 ===> activate
!     iopt = 1 ===> desactivate all or selected entries
!     iopt = 2 ===> list actvated theories
!     iopt = 3 ===> reactivate theories from lastth
!
!
       use cincom
       use cincoc
       use outlev
       use theory
       use theorc
       use therrc
       use thparc
       use constants
       implicit none

       integer iopt 

       character*1024  rlibuf
       integer       i,j,iocbuf,inew,inaml,ilook,ip,itt,it,ith
       integer       kk,npar,npars,mtca,ncc,nthn,lp
       character*4   label
       character*8   pname, chrval
       real*8        getval
       logical       found


       if(iopt.lt.0.or.iopt.gt.3) then
         write(6,*)' activate invalid option : iopt = ',iopt
         ierrs = 800
         return
       endif
!
       if(iopt.eq.0) then
!      -------------------> activate
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
   10    continue
         write(6,*)' wrong theoryname use theos to list available names'
         ierrs = 800
         return
   11    continue
         npar = nthpar(ith)
!        if(ipars.ne.npar*2) then
!          write(6,*)' wrong no. of parameters & scales'
!          write(6,*)' for each parameter a scale must be given!'
!          ierrs = 800
!          return
!        endif
!
         ntheos         = ntheos + 1
         nthtab(ntheos) = ith
         if(found('multiply')) then
           multflg(ntheos) = 1
           inaml           = 2
         else
           multflg(ntheos) = 0
           inaml           = 1
         endif

!! check whether Fit-Range is specified !!
         if(found('range   ')) then
            thrapar(ntheos) = chrval('range   ','        ',inew)
            thramin(ntheos) = getval('min     ',-1d30, inew)
            thramax(ntheos) = getval('max     ', 1d30, inew)
            write(6,*)'theo :',ntheos,' valid for ',thrapar(ntheos),    &
     &                 thramin(ntheos), thramax(ntheos)
         else
            thrapar(ntheos) = '        '
         endif

           do  i=1,npar
              thparx(i,ntheos) = 0
              thpsca(i,ntheos) = 0
              therro(i,ntheos) = 0
              thpala(i,ntheos) = '    '
              ncoup(i,ntheos)  = 0
           enddo

!
!         if(inames.eq.inaml) then
!c                        -----> assign noname parameterlist by sequence
!           do 20 i=1,npar
!              j1 = 2*i - 1
!              j2 = j1 + 1
!              thparx(i,ntheos) = rpar(j1)
!              thpsca(i,ntheos) = rpar(j2)
!              thpala(i,ntheos) = '    '
!              ncoup(i,ntheos)  = 0
!20         continue
!           return
!         endif
!c
!        if(inames.ne.npar+inaml) then
!          write(6,*)' none or all parameter scale pairs must be given'
!          write(6,*)' by names : parnam x s parnam x s ... npar=',npar
!          ntheos = ntheos - 1
!          ierrs = 800
!        endif
!        do 30 j=2,inames
!          j1 = 2*(j-1) - 1
!          j2 = j1 + 1
!          do 31 i=1,npar
!           if(vname(j).eq.thparn(i,ith)) then
!             thparx(i,ntheos) = rpar(j1)
!             thpsca(i,ntheos) = rpar(j2)
!             ncoup(i,ntheos) = 0
!             thpala(i,ntheos)= '    '
!             goto 30
!           endif
!31        continue
!          write(6,*)vname(j),' is no parametername ..'
!          ntheos = ntheos - 1
!          ierrs = 800
!          return
!30      continue

        write(6,*)vname(1),' is the ',ntheos,'-th  activated entry'
!
        return
       endif
!
       if(iopt.eq.1) then
!                    ----> clear activations
         if(ipars.eq.0.or.ntheos.le.1) then
           ntheos = 0
           write(6,*)' theory activation stack counter set to zero'
           return
         endif
! ---- desactivate only some theories ---
         do 1501 i=1,ipars
           j = rpar(i)+0.001
           if(j.lt.1.or.j.gt.ntheos) then
             write(6,*)' theory to be activated does not exist ..'
             ierrs = 700
             return
           endif
           nthtab(j) = 0
 1501    continue
! ---- update the stack ----
         nthn = ntheos
         it   = 1
 1511    continue
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
 1522          continue
 1512        continue
             nthn = nthn - 1
             if(it.le.nthn) goto 1511
            endif
            if(it.lt.nthn) then
               it = it + 1
               goto 1511
            endif
 1530      continue
           ntheos = nthn
         return
       endif
!
       if(iopt.eq.2) then
!                    ----> output of activated theories
!                    ----> and writing the makro lastth
         write(6,*)' ***** activated theories *****'
         if(ntheos.eq.0) then
          write(6,*)'      none, use  acl to restore lastth setting'
          return
!         -------> preserve the contents of file lastth in this case
         else
!
! ---- > write the file file lastth to restore the theory-params
          kk = 6
!         ------> write first to std-output...
          write(kk,*)'current parametersetting:'
 6666     continue


          if(kk.ne.6) open(kk,file='lastth',                            &
     &                        form='formatted',status='UNKNOWN')
!                                   ------
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
!           ---------> no save this to file
          else
            close(kk)
          endif
        endif
        return
       endif
!  ------> reactivate theories as stored in lastth
       if(iopt.eq.3) then
! ---- > read the file file lastth to restore the theory-params
          open(30,file='lastth',form='formatted',                       &
     &                          status='UNKNOWN',err=999)
!                       ------
! ---- save the old incom status ----
          rlibuf = reslin
          iocbuf = ioldc
! ---- initialize ----
          ntheos = 0
! ---- read an decode the parameter-file --------
 7777     continue
          read(30,'(a)') reslin
!         write(6,*)reslin
          ioldc = 1
!          call lowcase(reslin,80)
          call incom(comand)
!      write(6,*)'inames=',inames,' ipars=',ipars,' inpar(1)=',inpar(1),
!    *           ' inapa(1)=',inapa(1),' vname(1..3)=',vname(1),vname(2)
!    *           ,vname(3)
          if(comand.eq.'end     ') goto 8888
          if(comand.eq.'theory  ') then
            ntheos = ntheos + 1
            if(ntheos.gt.mtcal) then
               write(6,*)'too many theories ntheos=',ntheos,' max=',mtca
               ierrs = 703
               goto 999
            endif
!           ------- identify the theory ------
            do 7801 ith=1,mth
               if(vname(1).eq.thenam(ith)) goto 7802
 7801       continue
              write(6,*)'theoryname: ',vname(1),' is not known!'
              ierrs = 700
              goto 999
 7802       continue
!           --------> ok name is known
            nthtab(ntheos) = ith
            if(found('multiply')) then
              multflg(ntheos) = 1
            else
              multflg(ntheos) = 0
            endif
!! check whether Fit-Range is specified !!
            if(found('range   ')) then
              thrapar(ntheos) = chrval('range   ','        ',inew)
              thramin(ntheos) = getval('min     ',-1d30, inew)
              thramax(ntheos) = getval('max     ', 1d30, inew)
              write(6,*)'theo :',ntheos,' valid for ',thrapar(ntheos),  &
     &                   thramin(ntheos), thramax(ntheos)
            else
             thrapar(ntheos) = '       '
            endif

          else
!           -------- look now for a label -------
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
! ---------- look for the right entry ------------
            do 7803 i=1,nthpar(ith)
              if(pname.eq.thparn(i,ith)) then
                 thpala(i,ntheos) = label
                 thparx(i,ntheos) = rpar(1)
                 thpsca(i,ntheos) = rpar(2)
                 goto 7805
              endif
 7803       continue
            write(6,*) pname ,' is an invalid parametername'
            ierrs = 701
            goto 999
 7805       continue
! --------- look for coupleded items --------
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
 7807         continue
            endif
          endif
          goto 7777
! ---- restore the incom resline buffers --------
 8888     continue
          reslin = rlibuf
          ioldc  = iocbuf
          close(30)
          call couple(0)
          write(6,*)' theory-setting restored'
          return
!
  999     continue
          close(30)
          ierrs = 99
          write(6,*)' theories could not be restored ....'
          ntheos = 0
          return
       endif
!
      END
!
!
       subroutine theo_out(kk)
!      -----------------------

       use theory
       use selist
       use theorc
       use therrc
       use thparc
       use formul

       implicit none
       character*8 combinam,cha*1
       character*30 buf1
       character*60 buf2
       integer i, ith, kk, ier,l,j,lf

          if(ntheos.le.0) return

          do 150 i=1,ntheos
           ith = nthtab(i)
           if(multflg(i).eq.1) then
              write(buf1,*)'theory   ',thenam(ith),'   multiply'
           else
              write(buf1,*)'theory   ',thenam(ith)
           endif
           if(thrapar(i).ne.'        ') then
              write(buf2,'(2a8,a5,e13.6,a5,e13.6)')' range ',           &
     &             thrapar(i),' min ',thramin(i),' max ', thramax(i)
           else
              buf2 = ' '
           endif
           write(kk,'(a30,a60)')buf1,buf2


           do 150 j=1,nthpar(ith)
            if(kk.ne.6) then
             call setudf(thparn(j,ith)//' ',dble(thparx(j,i)),ier)
             call setudf('e.'//thparn(j,ith)//' ',dble(therro(j,i)),ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! bis zur Erweiterung auf lange Strings hier eine Notloesung zum Setze
             write(cha,'(i1)')i
             combinam = thenam(ith)(1:2)//cha//thparn(j,ith)(1:5)
             do lf=1,nsel
              if(ifits(lf).gt.0) then
               call parset(combinam,thparx(j,i),ifits(lf))
               call parset('e'//combinam(1:7),therro(j,i),ifits(lf))
              endif
              if(isels(lf).gt.0) then
               call parset(combinam,thparx(j,i),isels(lf))
!!??               call parset('e'//combinam(1:7),therro(j,i),isels(lf))
              endif
             enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            endif
             write(kk,170)thpala(j,i),thparn(j,ith),thparx(j,i),        &
     &                    thpsca(j,i),therro(j,i),                      &
     &                    (thpalc(l,j,i),thpafc(l,j,i),l=1,ncoup(j,i))
  170      format(1x,a4,1x,a8,1x,e12.4,1x,e8.1,1x,e9.2,                 &
     &                                              10(1x,a4,1x,f5.2))
  150     continue
          write(kk,*)'end'
          if(thenam(nthtab(1)).eq.'eval    ') write(kk,*)yfitform

       return
      END

       subroutine couple(iopt)
!      =======================
! ----- couple theory-parameters -----
!
! ----------------------------------------------------------------------
!      the theoryparameters are coupled linearily to other parameters
!      according to the label settings
!      the formula is:
!               p[l]    <== sum(i<>l) { f[l,i] p[i] }  + p0[l]
!               p0[l]   <== p[l|0] - sum(i<>l) { f[l,i] p[i|0] }
!      p[j|0] designates the startvalue of the j-th parameter and
!      f[l,i] is the linearcoupling factor for parameters l with i.
!
!      for iopt=0   the p[j|0] are evaluated and stored on thpaco
!          iopt=1   the actual p[l] are computed --> thparx
!      for a couplede parameter the fitscale thpsca is forced to 0
!
! ----------------------------------------------------------------------
!
!
       use outlev
       use theory
       use theorc

       implicit none

       character*4   label
       integer iopt, i, ii, n, np, npp, j, jj, ith, ithh, l
       real y

! ----------------------------------------------------------------------
! ----- go through all couplings ----
       do 10 i=1,ntheos
          ith = nthtab(i)
          np  = nthpar(ith)
          do 20 j=1,np
             n = ncoup(j,i)
             y = 0
             if(n.gt.0) then
                do 30 l=1,n
                   label = thpalc(l,j,i)
!                  ---- search now for a parameter with this label ---
                   do 40 ii=1,ntheos
                      ithh = nthtab(ii)
                      npp  = nthpar(ith)
                      do 40 jj=1,npp
                         if(label.eq.thpala(jj,ii)) goto 50
   40              continue
                      write(6,*)'couple: label=',label,' not found'
                      ierrs = 800
   50              continue
                   if((j.eq.jj).and.(i.eq.ii)) then
                      write(6,*)'couple: label=',label,' points to ',   &
     &                          'itself'
                      ierrs = 801
                   endif
!           ---- do the coupling ----
                   if(ierrs.eq.0) y = y + thpafc(l,j,i) * thparx(jj,ii)
   30           continue
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
!            ----- ( if(n.gt.0).... )
   20     continue
   10  continue
!
       return
      END
!
!

       subroutine input
!      ================
!
       use cincom
       use cincoc
       use cmargs
       use cdata
       use outlev
       use constants

       implicit none
       integer ln, le, la, lx, ly, j, inew, i,ier
       character*1024  infile
       character*80 rline
       logical*4    fileda
       real*8       xshift, yshift, getval

       logical      :: is_inx, is_ins 


!  neu fuer pfad
       character*1024 uspfad

       character*1024 X_datpath
       character*8 ignoreline,nextset
       integer     il

       save uspfad

        uspfad = X_datpath()     ! this is simply the datapath   Why so complicate?
		  ignoreline='#'
		  nextset='#nxt'
!
! -- some arguments
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

! pfad!!
		infile =argvals(1)
		!if local path by '.' or global path by / is indicated omit the datapath in front of content
      if (argvals(1)(1:1).ne.'.' .and. argvals(1)(1:1).ne.'/') infile = trim(uspfad)//infile
!
         ioldc = 0
         il = LEN_TRIM(infile)
         write(6,*)'opening: >',infile(1:il),'<'
         is_inx = (infile(il-2:il).eq.'inx').or.(infile(il-2:il).eq.'INX')
         is_ins = (infile(il-2:il).eq.'ins').or.(infile(il-2:il).eq.'INS')
         if(is_inx) write(6,*)' assuming inx format.... '
         if(is_ins) write(6,*)' assuming ins = SPHERES inx format.... '
 
!
       inquire(file=infile(1:il),exist=fileda)
       if(.not.fileda) then
         write(6,*)'input(e2): file >',trim(infile),'< does not exist!'
         ierrs = 3
         return
       endif
! open the file----------------------------------------------------------------------
       open(20,file=infile(1:il),status='UNKNOWN')
! -- main input loop !!!!!!!!!!!!!!!!!


!! --> Save the start address to a user-variable
     
      call setudf('read1    ',dble(nbuf+1),ier)                        
      call setudf('readlast ',dble(0),ier)



! -- and use a new subroutine for inx reading (clearer and more efficient) ---
      if(is_inx) then
        call inx_reading(20,infile)
        close(20)
        call setudf('readlast ',dble(nbuf),ier)
        return
      endif

      if(is_ins) then  ! spheres inx format
        call ins_reading(20,infile)
        close(20)
        call setudf('readlast ',dble(nbuf),ier)
        return
      endif


20000  continue                         !  new dataset
      nbuf = nbuf + 1
      nopar(nbuf) = 0                          ! Anzahl der parameter in dataset nbuf
      lx = 0   ! anzahl x y ey werte
      ly = 0
      le = 0

	if(nbuf.gt.mbuf) then
        write(6,*)'cdata buffer is full!'
        close(20)        
        call setudf('readlast ',dble(mbuf),ier)
        return
      endif
      name(nbuf)  = infile
      yname(nbuf) = 'y-data'
      xname(nbuf) = 'x-data'
      numor(nbuf) = 100*nbuf
      coment(nbuf)= ''
 2000 continue  !new line
		read(20,'(a)',end=40) rline            ! einziges lesecommando im loop    end=999  ist ausgang
		rline=ADJUSTL(rline)
		!write(6,*)trim(rline)
		call decode(rline)
		 !write(6,*)trim(rline)
		if(iout.gt.2) then                 ! only to show output
			write(6,*)'Decoded line: ',trim(rline),'#end rline'
			write(6,*)'ipars',ipars,'inames',inames,'iparn(1)',iparn(1),'vname(1)',vname(1),'nwert'
			write(6,*)'lx',lx
			do i=1,inames
        		write(6,"(a,1x,a8,2x,10f12.6)")"=",vname(i),(rpar(j),j=inapa(i),inapa(i)+inpar(i)-1)
        	enddo
			write(6,*)'############Decoding ende'
		endif
		!write(6,*)trim(rline)
 		! start a new set if an empty line follows the data so already datavalues exist
		if(     (ipars.eq.0.and.inames.eq.0)                         & !  an empty line
			.or.	trim(rline).eq.trim(nextset))       then
			if (lx .gt.0 ) then    !if we have already some data
				if(lx.ne.ly) then
!                		write(6,*)'number of x-values=',lx,' does not match ly=',ly
				endif
		  		nwert(nbuf)=lx       ! we tell the nwert buf how much
				write(6,"(I4,':lx',I4,' para',I4,' Comment: ',a)") nbuf,lx,nopar(nbuf),trim(coment(nbuf))
				goto 20000             !and start a nes set
			else
				goto 2000      !it was only an empty line btween params and data next line
			endif
		endif
		if(rline(1:len_trim(ignoreline)).eq.trim(ignoreline))    goto 2000  ! ignore lines starting with #  (will be a command argument later)
		if(trim(rline).eq.'values')    goto 2000
		if(trim(rline).eq.'parameter')    goto 2000
 		! start a new set if non data start the line after ly .gt.0
		if( (lx.gt.0 .and. .not.(ipars.gt.0.and.iparn(1).eq.0)) )   then   !non data line after data
				if(lx.ne.ly) then
               write(6,*)'number of x-values=',lx,' does not match ly=',ly
				endif
		  		nwert(nbuf)=lx       ! we tell the nwert buf how much data
				!######################
				write(6,"(I4,':lx',I4,' para',I4,' Comment: ',a)") nbuf,lx,nopar(nbuf),trim(coment(nbuf))
				nbuf = nbuf + 1
				nopar(nbuf) = 0                          ! Anzahl der parameter in dataset nbuf
				lx = 0   ! anzahl x y ey werte
				ly = 0
				le = 0
				if(nbuf.gt.mbuf) then
				write(6,*)'cdata buffer is full!'
				close(20)        
                                call setudf('readlast ',real(nbuf),ier)
				return
				endif
				name(nbuf) = infile
				yname(nbuf) = 'y-data'
				xname(nbuf) = 'x-data'
				numor(nbuf) = 100*nbuf
				coment(nbuf)= ''
				!#######################
				!   !and go on with actual rline
		endif


!        data have no names only values (ipars .ne. 0)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  		if (ipars.gt.0 .and. iparn(1).eq.0) then        !  data  line started with values
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
			if (ipars.gt.2) then
				yerror(le,nbuf) = rpar(3)
			!else
			!	yerror(le,nbuf) = 0
			endif
			goto 2000
      !endif  ! only data lines
		elseif(ipars.gt.0 .and. iparn(1).eq.1 .and.(trim(vname(1)).eq.trim('x').or.trim(vname(1)).eq.trim('y').or.trim(vname(1)).eq.trim('e')) ) then    ! we found parameters  starting with a name (starting with number is iparn(1)=0
			! !       data like x 1 2 3 4 5 6 y 9 8 7 6 5 e 1 1 1 1 1 1 ! ----- the more explicit lists -----with numbers
        	do j=1,inames
				la = inapa(j)
				ln = inpar(j)
				if(trim(vname(j)).eq.trim('x')) then ! found parameterlike x line
					do i=1,ln
					lx = lx + 1
					if(lx.gt.mwert) then
						write(6,*)'too many x-y values'
						close(20)
						return
					endif
					xwerte(lx,nbuf) = rpar(la+i-1)
					enddo
				elseif(trim(vname(j)).eq.trim('y')) then ! found parameterlike y line
					do i=1,ln
						ly = ly + 1
						ywerte(ly,nbuf) = rpar(la+i-1)
					enddo
				elseif(trim(vname(j)).eq.trim('e')) then ! found parameterlike er line
					do i=1,ln
						le = le + 1
						yerror(le,nbuf) = rpar(la+i-1)
					enddo
				endif
         	enddo
			goto 2000
      !endif ! parameter like data
		elseif(ipars.gt.0 .and. iparn(1).eq.1 .and. inames.eq.1 ) then    ! we found parameters  starting with a name (starting with number is iparn(1)=0
			nopar(nbuf) =   nopar(nbuf) + 1
			if(nopar(nbuf).gt.mpar) then
				write(6,*)'too many parameters'
				close(20)
				ierrs = 100
				return
			endif
			params(nopar(nbuf),nbuf) = rpar(inapa(1))
			napar (nopar(nbuf),nbuf) = vname(1)
			goto 2000
       !endif ! parameters
		elseif(vname(3).eq.'vs      '.or.vname(3).eq.'versus  ') then !!1   a line to set x ,y axis names
          name(nbuf) = vname(1)
          yname(nbuf) = vname(2)
          xname(nbuf) = vname(4)
          numor(nbuf) = rpar(1) * 1.000001
          goto 2000
		else ! --  identification & comment-line        evrything else
			if ( coment(nbuf).eq. ''  ) then
				coment(nbuf)= rline
			else
				write(6,*) 'Comment ignor: ',trim(rline)
			endif
       		goto 2000  ! end of main loop !    now next line
		endif
!
! -- error returns --
   40  continue
        if(lx.ne.ly) then
               write(6,*)'number of x-values=',lx,' does not match ly=',ly
        endif
        nwert(nbuf) = lx

       if (nwert(nbuf).eq.0) then ! keine Daten in letztem datensatz
	  nbuf=nbuf -1             ! nbuf  =>reset
   	 write (6,*) 'appended comments found and ignored '
       endif

        if(xshift.ne.0.0) then
          call parset ('xshift  ',sngl(xshift),nbuf)
        endif
        if(yshift.ne.0.0) then
          call parset ('yshift  ',sngl(yshift),nbuf)
        endif
      

        ierrs = 0
        close(20)
        call setudf('readlast ',dble(nbuf),ier)
        return
  999  continue
        close(20)
        call setudf('readlast ',dble(nbuf),ier)
       return
      END   ! input



       subroutine inx_reading(ikan,filnam)
!      ===================================
!
       use cincom
       use cincoc
       use cmargs
       use cdata
       use outlev
       use constants
       use PhysicalConstantsPlus

       implicit none
       integer ln, le, la, lx, ly, ji, inew, i, ier
       character*1024  infile
       character*80 rline, cline, filnam
       logical*4    fileda
       real*8       xshift, yshift, getval

       integer      :: ikan            ! Fortran Lesekanal
       integer      :: nk_tot, idum1, idum2, idum3, idum4,idum5,idum6, nk_act
       double precision :: angle, E_incident_meV, q_Ain, temperature, xdum1, xdum2
       double precision :: channel_width_musec
       double precision :: xvector(mwert), yvector(mwert), yervector(mwert)

       double precision :: Ein_Joule, Lambda_in_m, velocity_in_m, Q_elastic, lambda_in_A
       integer          :: numor0
       double precision :: xnum
       logical          :: found

       double precision :: t0, tx, kell, ejoule, lam_kina, ymx, velocx
       double precision :: flpath=4d0
!       double precision :: NeutronWavelength_fromE


! --und nun hier das Einlesen 

     xnum = 0
     call extract('numor0 ',xnum,ier)
     numor0 = NINT(xnum)
     Write(6,*)'adding ',numor0,' to autogenerated numors..'

     do ji=1,mbuf

      write(6,*)'inx reading record # ',ji, ' (',mbuf,')'

!     Nb tot channels, -, -, -, -, nb channels actually used
      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline
      read(rline,*) nk_tot,idum1,idum2,idum3,idum4,idum5,idum6,nk_act
      write(6,*)'  actual no of channels = ',nk_act

!     COMMENTLINE
      read(ikan,'(a)',end=999,err=999) rline
      cline = rline
      write(6,*)'  comment               = ',cline(1:60)
!     Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -,
      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline
      read(rline,*) angle, E_incident_meV, q_Ain, temperature, xdum1, xdum2
      write(6,*)'  angle                 = ',angle
      write(6,*)'  incident energy       = ',E_incident_meV,' meV'
      write(6,*)'  Q                     = ',q_Ain
      write(6,*)'  Temperature           = ',temperature,' K' 

      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline 
      read(rline,*) xdum1, channel_width_musec, xdum2 
      write(6,*)'  Channel width         = ', channel_width_musec,' musec' 
      
!
!     und nun die eigentlichen Daten
!
      if(nk_act.gt.mwert) then
        write(6,*)'too many channels...'
        return
      endif
       

      kell = 0
      ymx  = 0
      do lx=1,nk_act
        read(ikan,'(a)',end=999,err=999) rline
!        write(6,*)':',rline
        read(rline,*) xvector(lx),yvector(lx),yervector(lx)
        if(yvector(lx).gt.ymx) then
           ymx  = yvector(lx)
           kell = lx
        endif
      enddo
      write(6,*)'  Maximum channel kell    = ', kell
!
      write(6,*)'------------------------------------------------------'
! und nun auf den internen Speicher ...
      nbuf = nbuf + 1
      if(nbuf.gt.mbuf) then
         nbuf = nbuf - 1
         write(6,*)'cdata buffer is full!'
        return
      endif
      nopar(nbuf)   = 0                          ! Anzahl der parameter in dataset nbuf
      name(nbuf)    = filnam(1:8)
      infile        = trim(cline)//':'//trim(filnam)
      coment(nbuf)  = infile(1:80)
      numor(nbuf)   = ji + numor0
      call parset('angle   ',sngl(angle),nbuf)
      call parset('q_par   ',sngl(Q_Ain),nbuf)
      call parset('temp    ',sngl(temperature),nbuf)
      call parset('Ein     ',sngl(E_incident_meV),nbuf)
      call parset('binmusec',sngl(channel_width_musec),nbuf)
 
      Ein_Joule     =  E_incident_meV*1d-3*Elektronenladung
      lambda_in_m   =  NeutronWavelength_fromE(Ein_Joule)
      lambda_in_A   =  lambda_in_m*1d10
      velocity_in_m =  NeutronVelocity_fromLambda(lambda_in_m)
      Q_elastic     =  4*Pi/(lambda_in_A)*sin(Pi*angle/180.0d0/2)
      call parset('q       ',sngl(Q_elastic),nbuf)
      call parset('lambda  ',sngl(lambda_in_A),nbuf)
      call parset('veloc   ',sngl(velocity_in_m),nbuf)

      kell = getval('kell    ',kell,inew)
      call parset('kell    ',sngl(kell),nbuf)


      if(abs(xvector(2)-xvector(1)-1.0).lt.1e-5) then
        if(found('lambda  ')) then
          xname(nbuf)   = 'lambda'
          flpath = getval('flpath  ',flpath,inew)
          call parset('flpath  ',sngl(flpath),nbuf)
          do lx=1,nk_act
             t0       = flpath/velocity_in_m
             tx       = (lx-kell)*channel_width_musec*1d-6+t0
             velocx   = flpath/tx
             ejoule   = 0.5d0*Neutronenmasse*velocx**2
             lam_kina = NeutronWavelength_fromE(ejoule)*1d10
             xvector(lx) = lam_kina
          enddo
        else
          xname(nbuf)   = 'tof_bin#'
        endif
        yname(nbuf)   = 'tof_int '
      else
        xname(nbuf)   = 'meV     '
        if(found('GHz     '))then    ! omega = Kreisfrequenz
          xvector = xvector *1d-3*Elektronenladung*2*Pi/Planckkonstante/1d9
          xname(nbuf) = 'GHz     '
        endif
        yname(nbuf)   = 'sqw     '
      endif



! und noch die eigentlichen Daten
      if(nk_act.gt.mwert) then
        nk_act = mwert
        write(6,*)'number of channels too high: truncated'
      endif
      do lx=1,nk_act
        xwerte(lx,nbuf) = xvector(lx)
        ywerte(lx,nbuf) = yvector(lx)
        yerror(lx,nbuf) = yervector(lx)
      enddo
      nwert(nbuf) = nk_act


      call parset('_xwidth ',xwerte(kell+1,nbuf)-xwerte(kell,nbuf),nbuf)


      write(6,*)'sum of record ',nbuf,' = ',sum(yvector(1:nk_act))
      write(6,*)'next record ....'
      
     enddo
     write(6,*)'end of loop'
     return

999  continue
     write(6,*)'end of reading'
       return
      END   ! inx_reading


       subroutine ins_reading(ikan,filnam)
!      ===================================
!
! read that what SPHERES produces as inx 
!
       use cincom
       use cincoc
       use cmargs
       use cdata
       use outlev
       use constants
       use PhysicalConstantsPlus

       implicit none
       integer ln, le, la, lx, ly, ji, inew, i, ier
       integer :: idum
       character*1024  infile
       character*80 rline, cline, filnam
       logical*4    fileda
       real*8       xshift, yshift, getval

       integer      :: ikan            ! Fortran Lesekanal
       integer      :: nk_tot, idum1, idum2, idum3, idum4,idum5,idum6, nk_act
       double precision :: angle, E_incident_meV, q_Ain, temperature, xdum1, xdum2
       double precision :: channel_width_musec
       double precision :: xvector(mwert), yvector(mwert), yervector(mwert)

       double precision :: Ein_Joule, Lambda_in_m, velocity_in_m, Q_elastic, lambda_in_A
       integer          :: numor0
       double precision :: xnum
       logical          :: found

       double precision :: t0, tx, kell, ejoule, lam_kina, ymx, velocx
       double precision :: flpath=4d0
!       double precision :: NeutronWavelength_fromE


! --und nun hier das Einlesen 

     xnum = 0
     call extract('numor0 ',xnum,ier)
     numor0 = NINT(xnum)
     Write(6,*)'adding ',numor0,' to autogenerated numors..'

     do ji=1,mbuf

      write(6,*)'inx reading record # ',ji, ' (',mbuf,')'

!     Nb tot channels, -, -, -, -, nb channels actually used
      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline
      read(rline,*) nk_tot,idum1,idum2,idum3,idum4,idum5,idum6,nk_act
      write(6,*)'  actual no of channels = ',nk_act

!     COMMENTLINE
      read(ikan,'(a)',end=999,err=999) rline
      cline = rline
      write(6,*)'  comment               = ',cline(1:60)
!     Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -,
      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline
!     read(rline,*) angle, E_incident_meV, q_Ain, temperature, xdum1, xdum2    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      read(rline,*) q_Ain, E_incident_meV, angle, temperature, xdum1, xdum2    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if(angle.eq.0) then   ! try to reconstruct it from the other data
        Ein_Joule     =  E_incident_meV*1d-3*Elektronenladung
        lambda_in_m   =  NeutronWavelength_fromE(Ein_Joule)
        lambda_in_A   =  lambda_in_m*1d10
        velocity_in_m =  NeutronVelocity_fromLambda(lambda_in_m)
        angle         = 2*180/Pi*asin(lambda_in_A*q_Ain/4/Pi)
      endif

      write(6,*)'  angle                 = ',angle
      write(6,*)'  incident energy       = ',E_incident_meV,' meV'
      write(6,*)'  Q                     = ',q_Ain
      write(6,*)'  Temperature           = ',temperature,' K' 

      read(ikan,'(a)',end=999,err=999) rline
!      write(6,*)':',rline 
      read(rline,*) xdum1, channel_width_musec, xdum2 
      write(6,*)'  Channel width         = ', channel_width_musec,' musec' 
      
!
!     und nun die eigentlichen Daten
!
      if(nk_act.gt.mwert) then
        write(6,*)'too many channels...'
        return
      endif
       

      kell = 0
      ymx  = 0
      do lx=1,nk_act
        read(ikan,'(a)',end=999,err=999) rline
!        write(6,*)':',rline
!        read(rline,*) xvector(lx),yvector(lx),yervector(lx)
        read(rline,*) idum, xvector(lx),yvector(lx),yervector(lx)
        if(yvector(lx).gt.ymx) then
           ymx  = yvector(lx)
           kell = lx
        endif
      enddo
      write(6,*)'  Maximum channel kell    = ', kell
!
      write(6,*)'------------------------------------------------------'
! und nun auf den internen Speicher ...
      nbuf = nbuf + 1
      if(nbuf.gt.mbuf) then
         nbuf = nbuf - 1
         write(6,*)'cdata buffer is full!'
        return
      endif
      nopar(nbuf)   = 0                          ! Anzahl der parameter in dataset nbuf
      name(nbuf)    = filnam(1:8)
      infile        = trim(cline)//':'//trim(filnam)
      coment(nbuf)  = infile(1:80)
      numor(nbuf)   = ji + numor0
      call parset('angle   ',sngl(angle),nbuf)
      call parset('q_par   ',sngl(Q_Ain),nbuf)
      call parset('temp    ',sngl(temperature),nbuf)
      call parset('Ein     ',sngl(E_incident_meV),nbuf)
      call parset('binmusec',sngl(channel_width_musec),nbuf)
 
      Ein_Joule     =  E_incident_meV*1d-3*Elektronenladung
      lambda_in_m   =  NeutronWavelength_fromE(Ein_Joule)
      lambda_in_A   =  lambda_in_m*1d10
      velocity_in_m =  NeutronVelocity_fromLambda(lambda_in_m)
      Q_elastic     =  4*Pi/(lambda_in_A)*sin(Pi*angle/180.0d0/2)
      call parset('q       ',sngl(Q_elastic),nbuf)
      call parset('lambda  ',sngl(lambda_in_A),nbuf)
      call parset('veloc   ',sngl(velocity_in_m),nbuf)

      kell = getval('kell    ',kell,inew)
      call parset('kell    ',sngl(kell),nbuf)


      if(abs(xvector(2)-xvector(1)-1.0).lt.1e-5) then
        if(found('lambda  ')) then
          xname(nbuf)   = 'lambda'
          flpath = getval('flpath  ',flpath,inew)
          call parset('flpath  ',sngl(flpath),nbuf)
          do lx=1,nk_act
             t0       = flpath/velocity_in_m
             tx       = (lx-kell)*channel_width_musec*1d-6+t0
             velocx   = flpath/tx
             ejoule   = 0.5d0*Neutronenmasse*velocx**2
             lam_kina = NeutronWavelength_fromE(ejoule)*1d10
             xvector(lx) = lam_kina
          enddo
        else
          xname(nbuf)   = 'tof_bin#'
        endif
        yname(nbuf)   = 'tof_int '
      else
        xname(nbuf)   = 'micro_eV'
        if(found('GHz     '))then    ! omega = Kreisfrequenz
          xvector = xvector *1d-6*Elektronenladung*2*Pi/Planckkonstante/1d9
          xname(nbuf) = 'GHz     '
        endif
        yname(nbuf)   = 'sqw     '
      endif



! und noch die eigentlichen Daten
      if(nk_act.gt.mwert) then
        nk_act = mwert
        write(6,*)'number of channels too high: truncated'
      endif
      do lx=1,nk_act
        xwerte(lx,nbuf) = xvector(lx)
        ywerte(lx,nbuf) = yvector(lx)
        yerror(lx,nbuf) = yervector(lx)
      enddo
      nwert(nbuf) = nk_act


      call parset('_xwidth ',xwerte(kell+1,nbuf)-xwerte(kell,nbuf),nbuf)


      write(6,*)'sum of record ',nbuf,' = ',sum(yvector(1:nk_act))
      write(6,*)'next record ....'
      
     enddo
     write(6,*)'end of loop'
     return

999  continue
     write(6,*)'end of reading'
       return
      END   ! ins_reading



       subroutine inscn
!      ================ input of sv4-scantype data
!
       use cincom
       use cincoc
       use cdata
       use outlev
       use constants

       implicit none
!
       integer nscn, ndd, j, i
       real tavera, extramon, atemp, atime, ah, amonitor, al, anue, acounts, afeld, ak, acorrcnt
       character*8  infile
       character*131 rline
       dimension ah(mwert),ak(mwert),al(mwert),anue(mwert),atemp(mwert),&
     &           afeld(mwert),acounts(mwert),acorrcnt(mwert),           &
     &           amonitor(mwert),extramon(mwert),atime(mwert)
!
! -- open the file -----------------------------------------------------
       if(inames.lt.2) then
         ierrs = 1
         write(6,*)'inscn(e1): filename and/or ...type is lacking!'
         return
       endif
!
       infile = vname(1)
!      intyp  = vname(2)
!      inmode = 'f'
!      rline  = 'filedef 20 disk '//infile//' '//intyp//' '//inmode//
!    *          ' ( recfm  u '
!      write(6,*) rline(1:80)
!      call obey(rline,80,ier)
       write(6,*)'opening: ',infile,' ...'
       open(20,file=infile,status='UNKNOWN')
!
! -- read the items --

 2000  continue

       nbuf = nbuf + 1
       if(nbuf.gt.mbuf) then
         write(6,*)'cdata buffer is full!'
         close(20)
         return
       endif

! -- look for a header --
 2001  continue
         read(20,'(a80)',end=998)rline
 2002    continue
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
 1011  continue
       read(rline,*)nscn

! --- search for the beginning of the data-block ---
 2003  continue
         read(20,'(a131)',end=999)rline
         if(rline(1:12).ne.'  h       k ') goto 2003
         write(6,*)'1:',rline(1:70)
! --- read the data ---
       ndd = 0
 2004  continue
         read(20,'(a131)',end=999)rline
!cc      write(6,*)'*:',rline
         if(rline(1:1).ne.' ') goto 2010
         ndd = ndd + 1
         if(ndd.gt.mwert) then
           write(6,*)'warning too many data !'
           ndd = ndd-1
           goto 2010
         endif
         read(rline,*,end=2004)                                         &
     &                ah(ndd),ak(ndd),al(ndd),anue(ndd),atemp(ndd),     &
     &                afeld(ndd),acounts(ndd),acorrcnt(ndd),            &
     &                amonitor(ndd),extramon(ndd),atime(ndd)
       goto 2004
! ----------------------------------------------------------------------
!
 2010  continue

! ---- uebertragen der werte -----
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

! --- setzen von parametern ----
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
!      ---------> read the next scan

  998  continue
       nbuf = nbuf-1
  999  continue
       close(20)
       return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TOF routines                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine tof_dos(isource,idest,OMmax,Nom,ierrr)
!      =================================================
!
       use cincom
       use cincoc
       use cmargs
       use cdata
       use outlev
       use constants
       use PhysicalConstantsPlus


       implicit none

       integer                :: isource       ! record number of source data
       integer                :: idest         ! record number of destination data (must be different)
       double precision       :: OMmax         ! maximum Frequency of linear scale in result
       integer                :: Nom           ! number of bins in result
       integer                :: ierrr         ! error return value


       integer       :: ierr, inew, i, ier

       real*4        ::  xh

       integer      :: ikan            ! Fortran Lesekanal

       double precision :: angle, E_incident_meV, q_Ain, temperature, xdum1, xdum2

       double precision :: Ein_Joule, Lambda_in_m, velocity_in_m, Q_elastic, lambda_in_A
       double precision :: lambda0_angstroem, lambda_scale, omega_scale, E_elastic
       double precision :: lambda_1, lambda_2, lam_inc
       integer          :: numor0
       double precision :: xnum
       logical          :: found


      double precision :: k_i, k_f, qq
      double precision :: p_index_1, p_index_2, p1, p2, yval, yerr
      double precision :: Eom, bet, nbose, EomGHz, Eomega_start, Eomega_end
      double precision :: dosconv
      integer          :: index_1, index_2, ji



      write(6,*)'tofdos ', isource, idest, OMmax, Nom



      ierrr = 0
!!    first check   !!
      if(idest.eq.isource) then
        write(6,*)'source = destination is not possible', idest
        ierrr = 10
        return
      endif    

      if(idest.lt.1 .or. idest.gt.mbuf) then
        write(6,*)'destination is out of range', idest
        ierrr = 1
        return
      endif    

      if(isource.lt.1 .or. isource.gt.mbuf) then
        write(6,*)'source data address is out of range', isource
        ierrr = 2
        return
      endif    

      if(xname(isource).ne.'lambda  ') then
        write(6,*)'x-values are not lambda but',xname(isource)
        ierrr = 3
        return
      endif    

!! try to extract necessary parameters
      call parget('lambda  ', xh ,isource ,ier)
      if(ier.ne.0) then
       write(6,*)'Parameter lambda (incident lambda) not found'
       ierrr =4
       return 
      endif
      lambda0_angstroem = xh   ! assuming value is in Angstroem

      call parget('temp    ', xh ,isource ,ier)
      if(ier.ne.0) then
       write(6,*)'Parameter temp not found'
       ierrr =5
       return 
      endif
      temperature = xh   

      call parget('q       ', xh ,isource ,ier)
      if(ier.ne.0) then
       write(6,*)'Parameter q    not found'
       ierrr =6
       return
      endif
      Q_Elastic   = xh   

!! create destination
      call txfera(isource, idest)
      call txfpar(isource, idest)

      numor(idest) = numor(idest) + 2000000
      xname(idest) = 'GHz     '
      yname(idest) = 'DOS(om) '
      if(Nom.lt.mwert .and. Nom.gt.0) then
        nwert(idest) = Nom
      else
        nwert(idest) = mwert
      endif
! ..omega scale
      omega_scale = 1d9
      do i=1,nwert(idest)
        xwerte(i,idest) = i*OMmax/nwert(idest)     ! ok here we create a GHz scale
      enddo 
! .. and we assume lambda in A
      lambda_scale = 1d-10

! .. fill the histogram from tof data with energy gain
!    determine the lambda limits of the omega bin
!
     E_elastic =  NeutronEnergy_fromLambda(lambda0_angstroem*lambda_scale)
!     write(6,*)':4  E=',E_elastic
     Eomega_start = E_elastic
     do i=1,nwert(idest)
       Eomega_end = xwerte(i,idest)*omega_scale*Planckkonstante/2d0/Pi + E_elastic
       lambda_1 = NeutronWavelength_fromE(Eomega_end)   / lambda_scale               ! lambda and Q
       lambda_2 = NeutronWavelength_fromE(Eomega_start) / lambda_scale               ! in Angstroem units ! 
       k_i      = 2*Pi/(lambda0_angstroem)     
       k_f      = 4*Pi/(lambda_1+lambda_2)     
       qq       = k_i**2+k_f**2-2*k_i*k_f*(1d0-Q_Elastic**2/(2*k_i**2))
!
!      .. and now get the data range in the TOF spectrum ..
!
       lam_inc = (xwerte(nwert(isource),isource)-xwerte(1,isource))/(nwert(isource)-1)      ! also here we assume A units
       p_index_1 = (lambda_1-xwerte(1,isource))/lam_inc
       p_index_2 = (lambda_2-xwerte(1,isource))/lam_inc
       index_1   = INT(p_index_1)
       index_2   = INT(p_index_2)
       p1        = p_index_1-index_1
       p2        = p_index_2-index_2

!       write(6,*)': Eomega_start = ',Eomega_start  
!       write(6,*)': Eomega_en    = ',Eomega_end  
!       write(6,*)': lambda_1     = ',lambda_1  
!       write(6,*)': lambda_2     = ',lambda_2  
!       write(6,*)':6 ',i,p_index_1,p_index_2,index_1,index_2
     
       if(index_1.eq.index_2) then
         yval = (p2-p1)*ywerte(index_1,isource)
         yerr = (p2-p1)*yerror(index_1,isource)**2   ! this assumes counting error for a fraction of a bin!
       else
         yval = (1d0-p1)*ywerte(index_1,isource)+p2*ywerte(index_2,isource)
         yerr = (1d0-p1)*yerror(index_1,isource)**2+p2*yerror(index_2,isource)**2
         if(index_2-index_1 .gt. 1) then
          do ji=index_1+1, index_2-1
           yval = yval + ywerte(ji,isource)
           yerr = yerr + yerror(ji,isource)**2
          enddo
         endif
       endif
       yerr = sqrt(yerr)
    
       Eom      = E_elastic-(Eomega_start+Eomega_end)/2
       bet      = Eom/(temperature*Boltzmannkonstante)
       nbose    = 1d0/(exp(bet)-1d0)

       EomGHz   = 2*Pi*Eom/Planckkonstante/omega_scale

       dosconv = (k_i/k_f)/qq*EomGHz/(nbose+1d0)  ! * exp(2*W)

       ywerte(i,idest) = yval * dosconv
       yerror(i,idest) = yerr * dosconv

       Eomega_start = Eomega_end
     enddo

      return

      END












       subroutine decode(tline)
!      ========================
!
! ---- decode provides a means to treat a 80 characters line by the
!      incom facilities (comand is always set to '&')
!      the results are stored in the usual fashion in /cincom/
!
       use cincom
       use cincoc
       use outlev

       implicit none
       character*80 tline
       integer i,j
!
!
       if(iout.gt.2)  write(6,*)'Decode :',trim(tline),'#-------------'

       if(tline(79:80).ne.'  ') then
         write(6,*)'decode(w): truncating ',tline(79:80),'  of:'
         write(6,*)tline
       endif
!
       ioldc  = 1
       reslin = '& '//tline(1:78)

       if(iout.gt.2) then
         write(6,*)'Decode2:',trim(reslin)//'#-------------'
       endif

	call incom(comand)

       if(iout.gt.2) then
         write(6,*)'Decode3:',trim(comand)
       endif
!
! -- write control info if iout is greater 3 ---------------------------
       if(iout.gt.2) then
         write(6,*)tline
         write(6,*)'inames=',inames,'  ipars=',ipars
         if(inames.ne.0) then
          do i=1,inames
           write(6,"(1x,a8,2x,10f12.6)")vname(i),(rpar(j),j=inapa(i),inapa(i)+inpar(i)-1)
          enddo
         else
          write(6,200)(i,rpar(i),i=1,ipars)
  200     format(' decode numerical parameters:'(1x,i2,':  ',f12.6))
         endif
       endif
!-----------------------------------------------------------------------
!
       return
      END
!
!
       real function smirro(y,y2,xk0,k1,k2,n)
!      =================================
! --- computes the mirrorimage on the k1 to k2 region mirrored at xk0
!     onto the k1 to k2 region of y2 ( y2 is set to zero for the rest)
!     the squared error difference of y1 and y2 is retuned as function
!     value
!     n is the total number of valid points
!
       use constants
       implicit none
!
       real y,y2,xk0
       integer k1,k2,n

       integer i,kp, km
       real xmirro, sum, t
       dimension y(mwert),y2(mwert)
!
        do 10 i=1,n
         y2(i) = 0
   10   continue
!
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
  100  continue
       smirro = sum
       return
      END
!
!
       subroutine detxk1(y,y2,xk0,k1,k2,n)
!      ===================================
!
! ---- determine the xk0 - value with the minimum error to miirorsymm.
!      for the description of parameters see function smirro
!
       use constants
       implicit none

       real y,y2,xk0
       integer k1,k2,n

       integer i, j, m
       real f, f1, step, x, x1, xk00, smirro
       dimension y(mwert),y2(mwert)

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
  100  continue
        xk0 = x1
!
        return
      END
!
!
       subroutine detxk0(y,y2,xk0,k1,k2,n)
!      ===================================
!
! ---- determine the xk0 - value with the minimum error to miirorsymm.
!      for the description of parameters see function smirro
!
       use constants
       implicit none

       real y,y2,xk0,k1,k2
       integer n

       dimension y(mwert),y2(mwert)
       real a,x1,x2,x3,smirro,xk00,f0,f1,f12,f2,f3,f23
       integer i

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
  100   continue
  200   continue
        if(f2.gt.f0) xk0 = 0
        return
      END
!
!
       subroutine txfpar(ia,ib)
!      ========================
! --- copies all accompanying infos from data slot ia to slot ib
! !
       use cdata
       implicit none

       integer ia, ib
       integer i, np
!
       np=nopar(ia)
       nopar(ib)=np
        do 10 i=1,np
          params(i,ib)=params(i,ia)
          napar(i,ib)=napar(i,ia)
   10   continue
       name(ib)  = name(ia)
       coment(ib)= coment(ia)
       xname(ib) = xname(ia)
       yname(ib) = yname(ia)
       nwert(ib) = nwert(ia)
       numor(ib) = numor(ia)
!
       return
      END
!
!
       subroutine txfera(ia,ib)
!      ========================
! ---- copy the contens of data slot ia into slot ib -----
!
       use cdata
       implicit none

       integer ia, ib
       integer i, np

       np=nopar(ia)
       nopar(ib)=np
       do 5 i=1,nwert(ia)
          xwerte(i,ib) = xwerte(i,ia)
          ywerte(i,ib) = ywerte(i,ia)
          yerror(i,ib) = yerror(i,ia) 
    5  continue
       do 10 i=1,np
         params(i,ib)=params(i,ia)
         napar(i,ib)=napar(i,ia)
   10  continue
       name(ib)  = name(ia)
       coment(ib)= coment(ia)
       xname(ib) = xname(ia)
       yname(ib) = yname(ia)
       nwert(ib) = nwert(ia)
       numor(ib) = numor(ia)
!
       return
      END



       subroutine savdat(fname,ispc)
!      =============================
! -- save data on adress ispc onto a file named file fname a --
!
       use cdata
       implicit none

       character*1024 fname,finame,savepath, pathbuf, outfile
       integer i, ispc

		 pathbuf = savepath()
!
! ----------> transfer new name if not scratchfile for edit <-----------
         if(trim(fname).eq.trim('datbuf   ')) then
            finame = name(ispc)
         else
            finame = fname
         endif
! ----------> write data onto buffer <---------------------------------
			if (fname(1:1).ne.'.' .and. fname(1:1).ne.'/') then
				outfile = trim(pathbuf)//trim(fname)
			else
				outfile = trim(fname)
			endif
         write(6,*)'opening:',trim(outfile),'....'
         open(18,file=trim(outfile),status='UNKNOWN',err=999)
         write(18,'(a)')trim(coment(ispc))
         write(18,'(a,a,a,a,a,i14)')trim(finame(index(finame,'/',back=.true.)+1:)),' ',yname(ispc)(1:20),' vs ',xname(ispc)(1:20),numor(ispc)
         write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i=1,nopar(ispc))
         write(18,*)' '
         write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),i=1,nwert(ispc))
! 501      format(2x,'x  ',e14.7,5x,'y  ',e14.7,5x,'e  ',e14.7)
  501    format(2x,'   ',e14.7,5x,'   ',e14.7,5x,'   ',e14.7)
         write(18,*)' '

         call theo_out(18)

         close(18)
! ----------------------------------------------------------------------
         write(6,*)'file(',ispc,') : ',name(ispc),                      &
     &             ' has been saved to  ',                              &
     &              trim(outfile)
       return
999    continue
       call errsig(999,'file open failed...$')
       return
      END



       subroutine msavdat(fname)
!      =========================
! -- MULTI save data on adress ispc onto a file named file fname a --
!
       use cdata
       use selist
       implicit none

       character*1024 fname,savepath, pathbuf, outfile

       integer i, ispc, l

! first write data selected:
         pathbuf = savepath()
!                                                                                                                                                outfile = trim(pathbuf)//fname
         if (fname(1:1).ne.'.' .and. fname(1:1).ne.'/') then
				outfile = trim(pathbuf)//trim(fname)
			else
				outfile = trim(fname)
			endif
         write(6,*)'opening:',trim(outfile),'....'
         open(18,file=trim(outfile),status='UNKNOWN',err=999)
         do l=1,nsel
				ispc = isels(l)
				write(18,'(a)')trim(coment(ispc) )
				write(18,'(a,a,a,a,a,i14)')trim(name(ispc)(index(name(ispc),'/',back=.true.)+1:)),' ',yname(ispc)(1:20), ' vs ',xname(ispc)(1:20),numor(ispc)
				write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i= 1,nopar(ispc))
				write(18,*)' '
				write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),i=1,nwert(ispc))
501	      format(2x,'   ',e14.7,5x,'   ',e14.7,5x,'   ',e14.7)
				write(18,*)' '
				if(ifits(l).gt.0) then
					ispc = ifits(l)
					write(18,'(a)')coment(ispc)
					write(18,'(a,a,a,a,a,i14)')name(ispc)(1:8) ,' ',yname(ispc)(1:20),' vs ',xname(ispc)(1:20) ,numor(ispc)
					write(18,'(2x,a8,10x,e14.7)')(napar(i,ispc),params(i,ispc),i= 1,nopar(ispc))
					write(18,*)' '
					write(18,501)(xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),i=1,nwert(ispc))
					write(18,*)' '
				endif
         enddo
         call theo_out(18)
         close(18)
! ----------------------------------------------------------------------
         write(6,*)'data and fits have been saved to: ',trim(outfile)
       return
999    continue
       call errsig(999,'file open failed...$')
       return
      END


       subroutine sumseldat
!      ====================
! -- Summation of data in sellist
!
       use cdata
       use selist
       use constants

       implicit none

       character*1024 fname,savepath, pathbuf, outfile

       integer i, ispc, l, j, nn

       double precision :: xh, deltax, sumy, sumyerq

       if(nsel.le.0) then
         call errsig(999,'no data selected..$')
         return
       endif

       if(nbuf.ge.mbuf) then
         call errsig(999,'maximum number of data records ...$')
       endif
 
       write(6,*)'Summation of ',nsel,' selected records...'
! first write data selected:
       nbuf = nbuf + 1

       ispc = isels(1)
       nn   = nwert(ispc)
!      (xwerte(i,ispc),ywerte(i,ispc),yerror(i,ispc),i=1,nwert(ispc))
       do j=1,nn
         sumy    = 0
         sumyerq = 0
         xh      = xwerte(j,ispc)
         deltax  = abs(xwerte(nn,ispc)-xwerte(1,ispc))/(5.0*nn)
         do i=1,nsel
           sumy     =     sumy    + ywerte(j,isels(i))
           sumyerq  =     sumyerq + yerror(j,isels(i))**2
           if(abs(xwerte(j,isels(i))-xwerte(j,ispc)).gt.deltax) then
              write(6,*)i,':',isels(i),'  xpoint:',j
              write(6,*)xwerte(j,isels(i)), xwerte(j,ispc)
              call errsig(999,'x-values are different...$')
              nbuf = nbuf-1
              return
           endif
          enddo
          xwerte(j,nbuf) = xh
          ywerte(j,nbuf) = sumy
          yerror(j,nbuf) = sqrt(sumyerq)
       enddo
     
! transfer of parameters
       call txfpar(ispc,nbuf)

! set average parameters where appropriate
! ... yet to be done properly
! ... for now take the parameters from center of range
!
      if(nsel.gt.1) then
        j = nsel/2
      else
        j = 1
      endif
      call txfpar(isels(j),nbuf)
  
     
      Write(6,*)' sum stored at position: ', nbuf
      Write(6,*)' attention: parameters are copied from element',j,' at position ',isels(j)
     
      nsel     = 1
      isels(1) = nbuf

     

       return
      END



       real function fsplin(x)
!      ==================
! --- spline interpolated data evaluation ---
!     spline coefficients are taken from the last call of spline
       use outlev
       use cfc

       implicit none 

       real x, xx, yd, zn
       integer j, nn

         nn  = nwspl
         do 2705 j=1,nn
           if(break(j).ge.x) goto 2707
 2705    continue
           if(iout.ge.0) then
             write(6,*)'fsplin: x-value is out of range x =',x,xx
             write(6,*)'range:',break(1),break(nn)
           endif
           fsplin = 0
           return
 2707    continue
         zn = x - break(j)
!  ---- compute the interpolated data ----
         yd = ((cscoef(4,j)*zn+cscoef(3,j))*zn+cscoef(2,j))*zn+         &
     &          cscoef(1,j)
!
         fsplin = yd
!        -----------
       return
      END



       real function fdes(x)
!      ================
! --- integrand for desmearing ---
       use cfc
       implicit none

       real x, xx, yd, zn
       integer nn, j

         nn  = nwspl
         xx = sqrt(qziel**2 + x**2)
         do 2705 j=1,nn
           if(break(j).ge.xx) goto 2707
 2705    continue
             write(6,*)'fdes x-value is out of range x =',x,xx
             write(6,*)'range:',break(1),break(nn)
             fdes = 0
             return
 2707    continue
         zn = xx - break(j)
!  ---- compute the first derivative ----
         yd =  (3*cscoef(4,j)*zn+2*cscoef(3,j))*zn+cscoef(2,j)
!
         fdes = yd / xx
!        --------------
       return
      END


       subroutine fftmx(x,y,t,xmax,dx,nx,nfft)
!      =======================================
!  --- small angle multiple scattering evaluator (test-version) ---
!
!   -- it is assumed that the input data (derived form the spline
!      coefficients by fsplin) represent some scattering law in the
!      form:
!      yi <---  i0 * sigma(x) * [unit sample thickness:d]   =
!               lim[d-->0] ( i-measured(x,d) / d )
!      the output data on y will the be scaled such that they represent
!      yo <---  i0 * i-measured(x,d) / d      (note: no limit !!)
!      ---------------------------------
!      i.e. the scattering intensity including all orders of multiple
!      elastic small angle scattering will be calculated for a sample
!      of unit thickness assuming the same primary intensity factor,i0,
!      as for the input data.
! ----------------------------------------------------------------------
!      input variables:
!      t ....... : transmission reduction factor due to small angle sc.
!      xmax .... : largest scattering angle (or q) to be considered
!      dx ...... : increment of x
!      nx ...... : no. of points to be generated
!      nfft .... : fft no. of points (optimal choice nfft=2**m)
!      output variables:
!      x(1..nfft/2) : output x-values
!      y(1..nfft/2) : output y-values
! ----------------------------------------------------------------------
! --- fft-dimensioning -------------------------------------------------
       use outlev
       use constants
       use fftwrk
       implicit none

       real x,y,t,xmax,dx
       integer nx,nfft

       complex*8  cexp,cabs
       dimension x(*),y(*)

       real xx, yy , x1, x2, xl, xn, fsplin, dy, dd, ai1, ai0, at
       integer ioold, i, j
! ----------------------------------------------------------------------

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
! ----- period of fft is:
       xl = nfft * dx
! ----- reciprocal fourier increment:
       dy = 1 / xl
! ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
  100  continue
! ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
! ---- init data fields ----
       do 110 i=1,nfft
        do 110 j=1,nfft
         ca(i,j) = (0.,0.)
  110  continue
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
  200  continue
       iout = ioold
!
       write(6,*)'starting fft....'
       call  fft2d(nfft,nfft,ca,lda,ca,lda)
       write(6,*)'ca(1,1)=',ca(1,1)
       ai1= cabs(ca(1,1))
       ai0= ai1          /(1.-t)
! ---  === : i0 intensity factor provided the original data do not
!            contain any primary beam components. otherwise the
!            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
       xn = -at / cabs(ca(1,1))
! ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
!      xn = normalisation factor of the fourier transform to meet the
!           absolute scattering cross-section scale
       do 150 i=1,nfft
        do 150 j=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i,j)
!
         ca(i,j) = cexp( xn*ca(i,j) + at ) - t
! ----                                       = <- no primary beam int.
!                                              will be in the results !
!                                              omitt this t if primary
!                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i,j)
  150  continue
!
! ---- reverse transformation back to scattering angle field ----
       call  fft2b(nfft,nfft,ca,lda,ca,lda)
!
! ---- return values ---------------------------------------------------
! ---- normalize to the same total scattering --------------------------
       xn = ai1 * dy * dy / (-at)
       do 1000 i=1,nfft/2
        if(iout.gt.0) write(6,*)'result(',i,')= ',ca(i,1) * xn
        y(i) = real(ca(i,1)) * xn
 1000  continue
! ----------------------------------------------------------------------
!
       return
      END



       subroutine fftrmx(x,y,ai0,t,thick,xmax,dx,nx,nfft)
!      ==============================================
!  --- small angle multiple scattering inverse evaluator ----
!   -- it is assumed that the spline coefficients that are to be
!      generated by the datreatx command spline from the data to be
!      cleaned from multiple scattering will be used by fsplin to
!      provide equally spaced interpolated (and smoothed) original
!      data for the first fourier transform.
!      the output will then be:
!      y <--- i0 * sigma(x)
!      provided the appropriate value for the sample thickness has been
!      supplied.
!      the primary intensity factor has been included to faciliate the
!      comparison of original data and cleaned data on the same scale
!      (in that case thickness = 1 may be appropriate).
!      i0 will be printed to allow for a subsequent computation of sigma
!      sigma corresponds the cross-section per unit area and thickness
!      times the solid angle represented by the spacing dx of x-values.
!      solid angle ---> dx * dx /(ki**2)
! ----------------------------------------------------------------------
!      input variables:
!      t ....... : transmission reduction factor due to small angle sc.
!      thick ... : sample thickness (units will be the units of the res.
!      xmax .... : largest scattering angle (or q) to be considered
!      dx ...... : increment of x
!      nx ...... : no. of points to be generated
!      nfft .... : fft no. of points (optimal choice nfft=2**m)
!      output variables:
!      x(1..nfft/2) : output x-values
!      y(1..nfft/2) : output y-values
!      ai0 ........ : primary intensity varaible as assumed (i0)
! ----------------------------------------------------------------------
       use outlev
       use fftwr2
       use constants
       implicit none

       real x,y,ai0,t,thick,xmax,dx
       integer nx,nfft
! --- fft-dimensioning -------------------------------------------------
       dimension x(*),y(*)

       real xx,yy, x1,x2,xl,yinterp,rca, dd, at, dy
       integer ioold, i, j

! ----------------------------------------------------------------------
!
       if(nfft.eq.0)       nfft = mdim
       if(nfft.gt.mdim)    nfft = mdim
       if(nx.gt.nfft) nx        = nfft
       dd = xmax / (nx-1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   if(dx.lt.dd) dx = dd
       dx = dd
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
! ----- period of fft is:
!cccc  xl =     nfft * dx
       xl = 2 * nfft * dx
! ----- reciprocal fourier increment:
       dy = 1 / xl
! ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
  100  continue
! ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
! ---- init data fields ----
       do 110 i=1,nfft
        do 110 j=1,nfft
         ca(i,j) = 0.
  110  continue
       ioold= iout
       iout = -1
       do 200 i=1,nfft
        do 200 j=1,nfft
         x1 = (i-1) * dx
         x2 = (j-1) * dx
         xx = sqrt(x1**2 + x2**2)
         if(xx.gt.xmax) goto 200
!cccc    yy = fsplin(xx) * dx * dx
         yy = yinterp(xx) * dx * dx
         ca(       i,       j) = yy
  200  continue
       iout = ioold
!
       write(6,*)'starting fft....'
       call  rfft2d(nfft,nfft,ca,lda,ca,lda)
!!!!!  ai0= cabs(ca(1,1))/(1.-t)
       ai0=     (ca(1,1))/(1.-t)
!!!!!
! ---  === : i0 inttensity factor provided the original data do not
!            contain any primary beam components. otherwise the
!            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
!cccc  xn = -at / cabs(ca(1,1))
! ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
!      xn = normalisation factor of the fourier transform to meet the
!           absolute scattering cross-section scale
       do 150 i=1,nfft
        do 150 j=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i,j)/ai0 + t
!
!!!!!c   ca(i,j) = clog( ca(i,j)/ai0       + t) - at
         rca     = ca(i,j)
         ca(i,j) = alog( rca          /ai0 + t) - at
!!!!!c
! ----                                       = <- no primary beam int.
!                                              will be in the results !
!                                              omitt this t if primary
!                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i,j)
  150  continue
!
! ---- reverse transformation back to scattering angle field ----
       call  rfft2d(nfft,nfft,ca,lda,ca,lda)
!
! ---- return values ---------------------------------------------------
       do 1000 i=1,nfft/2
        if(iout.gt.0) write(6,*)'result(',i,')=',ca(i,1)*dy*dy*ai0/thick
!c      y(i) = real(ca(i,1)) * dy * dy * ai0 / thick
        y(i) =     (ca(1,i)) * dy * dy * ai0 / thick
 1000  continue
       write(6,*)'intensity factor i0 =',ai0
       write(6,*)'                     ','=================='
! ----------------------------------------------------------------------
!
!
       return
      END



       subroutine rmx1d(x,y,ai0,t,thick,xmax,ymax,nx,nfft)
!      =================================================
!  ---> like fftrmx but using 1d-fourier-transforms <---
!  --- small angle multiple scattering inverse evaluator ----
!   -- it is assumed that the spline coefficients that are to be
!      generated by the datreatx command spline from the data to be
!      cleaned from multiple scattering will be used by fsplin to
!      provide equally spaced interpolated (and smoothed) original
!      data for the first fourier transform.
!      the output will then be:
!      y <--- i0 * sigma(x)
!      provided the appropriate value for the sample thickness has been
!      supplied.
!      the primary intensity factor has beened to faciliate the
!      comparison of original data and cleaned data on the same scale
!      (in that case thickness = 1 may be appropriate).
!      i0 will be printed to allow for a subsequent computation of sigma
!      sigma corresponds the cross-section per unit area and thickness
!      times the solid angle represented by the spacing dx of x-values.
!      solid angle ---> dx * dx /(ki**2)
! ----------------------------------------------------------------------
!      input variables:
!      t ....... : transmission reduction factor due to small angle sc.
!      thick ... : sample thickness (units will be the units of the res.
!      xmax .... : largest scattering angle (or q) to be considered
!      ymax .... : fouriertransform-range
!      nx ...... : no. of points to be generated
!      nfft .... : fft no. of points (optimal choice nfft=2**m)
!      output variables:
!      x(1..nfft/2) : output x-values
!      y(1..nfft/2) : output y-values
!      ai0 ........ : primary intensity varaible as assumed (i0)
! ----------------------------------------------------------------------
       use outlev
       use constants
       implicit none

       real x,y,ai0,t,thick,xmax,ymax
       integer nx,nfft

       real xx, qq, sum, dx, dy, ca, cb, bsj0, at, yinterp
       integer i,j

! --- fft-dimensioning -------------------------------------------------
       dimension  ca(mdim)   ,cb(mdim)
!             ---> this common saves storage by use of ca in fftrmx also
       dimension x(*),y(*)
! ----------------------------------------------------------------------

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
! ----- reciprocal fourier increment:
       dy = ymax / (nfft-1)
! ---- generate the x-vector --------
       do 100 i=1,nx
        x(i) = (i-1) * dx
  100  continue
! ---- generate the 2d field of y-values ----
       write(6,*)'generating data field to be transformed ...'
! ---- init data fields ----
!!??   iout = ioold
!
       write(6,*)'starting transform ....'
       do i=1,nfft
        qq  = dy*(i-1)
        sum = 0
        do j=1,nx
          xx = dx*(j-1)
!ccc      sum = sum + fsplin(xx)*xx*bsj0(xx*qq)
          sum = sum + yinterp(xx)*xx*bsj0(xx*qq)
        enddo
        cb(i) = sum*dx
        if(iout.gt.0)write(6,*)i,qq,cb(i)
       enddo

       write(6,*)'cb(1) = ',cb(1)
       ai0=     (cb(1))/(1.-t)
! ---  === : i0 inttensity factor provided the original data do not
!            contain any primary beam components. otherwise the
!            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
!ccc   xn = -at /      cb(1)
! ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
!      xn = normalisation factor of the fourier transform to meet the
!           absolute scattering cross-section scale
       do 150 i=1,nfft
         if(iout.gt.2) write(6,*)'1:',ca(i)/ai0 + t
!
         ca(i  ) = alog( cb(i  )/ai0       + t) - at
! ----                                       = <- no primary beam int.
!                                              will be in the results !
!                                              omitt this t if primary
!                                              beam is to be generated.
         if(iout.gt.2) write(6,*)'2:',ca(i  )
  150  continue
!
! ---- reverse transformation back to scattering angle field ----
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
!
! ---- return values ---------------------------------------------------
       do 1000 i=1,nx
        if(iout.gt.0) write(6,*)'result(',i,')=',cb(i  )*ai0/thick
        y(i) = (cb(i  ))* ai0 / thick
 1000  continue
       write(6,*)'intensity factor i0 =',ai0
       write(6,*)'                     ','=================='
! ----------------------------------------------------------------------
!
!
       return
      END
!
!
       subroutine rfft2d(n1,n2,a,ldaa,b,ldb)
!      ==================>==>==>==>==<==>==
!                              a === b ok.
!
! 2d-reelle cos-fouriertransformation
! -----------------------------------
       use constants
       implicit none

       integer n1,n2,ldaa,ldb
       real a,b

       integer i, j
       real y

       dimension a(ldaa,ldaa)
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
      END


       real function yinterp(x)
!      ---------------------
! --> interpolates y-value at x for the file on top of the
!     selection list

       use cdata
       use selist
       implicit none

       real x
       integer i, ia, n

       ia = isels(1)
       n  = nwert(ia)
       if(x.lt.xwerte(1,ia).or.x.gt.xwerte(n,ia)) then
          yinterp = 0.0
          return
       endif
       do i=2,n
        if(xwerte(i,ia).gt.x) then
           yinterp = ywerte(i,ia)*(x           -xwerte(i-1,ia))/        &
     &                            (xwerte(i,ia)-xwerte(i-1,ia))         &
     &             + ywerte(i-1,ia)*(xwerte(i,ia)-x)/                   &
     &                              (xwerte(i,ia)-xwerte(i-1,ia))
           return
        endif
       enddo
       yinterp = 0.0
       return
      END


       subroutine demux(xin,yin,nin,t,xmax,nfft,iout,xout,yout,ier)
!      ============================================================
!  --- small angle multiple scattering inverse evaluator ----
!      input:
!      xin(1..nin) ... : inputvector containig scatteringangles
!                        or equivalent in ascending order
!      yin(1..nin) ... : inputvector containing scattering intensities
!      nin ........... : active length of xin and yin
!      t ............. : effective transmission factor associated
!                        with the small angle scattering
!      xmax .......... : x (scattering vector/angle) range
!                        for the result
!      nfft .......... : number of points of the internal
!                        fouriertransforms = the number of
!                        points in the output vectors.
!                        preferred value: 2**n.
!      iout .......... : variable to control additional printoutput
!                        from this subroutine
!      output:
!      xout(1..nfft).. : output x-values
!      yout(1..nfft).. : output y-values
!      ier ........... : errorindicator (0=ok)
!
!      the output will be the hypothetical scattering data of the
!      same sample using the same apparatus
!      but without any effects of multiple scattering,
!      e.g. no multiple sacttered intensity and no transmission
!      reduction.
! ----------------------------------------------------------------------
!
! ---- dimensioning ----------------------------------------------------
       use fftwr1
       use constants
       parameter(pi=3.1415926535897)
       real xin (nin) ,yin (nin)
       real xout(nfft),yout(nfft)
                                                       !! k--> 1 aix !!
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

!------- data preparation ----------------------------------------------
! check for the first and last nonzero element
         ifirst = 0
         ilast  = nin
         do i=1,nin
           if(yin(i).ne.0.0) then
              ifirst = i
              goto 1
           endif
         enddo
    1    continue
         do i=nin,1,-1
           if(yin(i).ne.0.0) then
              ilast  = i
              goto 2
           endif
         enddo
    2    continue
         if(ifirst.eq.0) then
           write(6,*)'demux: data are identical zero !'
           ier = 3
           return
         endif
! extrapolate data to x --> 0 assuming a parabola
         slope = (yin(ifirst+2)-yin(ifirst))                            &
     &          /(xin(ifirst+2)-xin(ifirst))
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
!
! -- high-q-extrapolation --
!
         slope = (yin(ilast)-yin(ilast-2))/                             &
     &           (xin(ilast)-xin(ilast-2))
         zparhq= xin(ilast-1)*slope / yin(ilast-1)
         aparhq= yin(ilast-1)*xin(ilast-1)**(-zparhq)
         if(iout.gt.0) then
           write(6,*)'high-q extrapolation :  a*q^z'
           write(6,*)'a =  ',aparhq
           write(6,*)'z =  ',zparhq
         endif
!
! --- fill the rest of the interpolation vector with equidistant data --
!
         ninter = xmax / dxin + 1
         if(ninter.gt.mdim) then
           write(6,*)'not enough space to store intermediate',          &
     &               ' interpolation vector'
           ier = 2
           return
         endif
         jl = ifirst+1
         do i=npl+1,ninter
           xxin = (i-1)*dxin
           do j=jl,ilast
             if(xin(j).gt.xxin) goto 10
           enddo
   10      continue
           jl = j
           if(jl.gt.ilast) then
             yinter(i) = aparhq * xxin**zparhq
           else
             p = (xxin-xin(jl-1))/(xin(jl)-xin(jl-1))
             yinter(i) = yin(jl)*p + yin(jl-1)*(1.0-p)
           endif
         enddo
!debug
         call pushda(yinter,0.0,dxin,ninter)
!debug
!
! ----- prepare the 1d-data for 2d-fft -------------------------------
       dx = xmax / (nfft-1)
! period of fft is:
       xl = 2 * nfft * dx
! reciprocal fourier increment:
       dy = 1.0 / xl
! ---- generate the 2d field of y-values ----
       if(iout.gt.1)write(6,*)'generating data field to be transformed'
! ---- init 2d-intermediate data field for fft :
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
!
!debug
         do i=1,nfft
            yinter(i) = ca(1,i)/dx/dx
         enddo
         call pushda(yinter,0.0,dx,nfft)
!debug
       if(iout.gt.1)write(6,*)'starting 1. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
!debug
       call pushda(ca(1,1),0.0,dy,nfft)
!debug
!
       ai0= ca(1,1)/(1.-t)
! ---  === : i0 intensity factor provided the original data do not
!            contain any primary beam components. otherwise the
!            1/(1-t) factor has to be omitted or modified.
       at = alog(t)
! ---  at = d * sigma(total-kws) as derived from the transmission(kws) -
!
! removing of multiple scattered events and transmission correction
       do i=1,nfft
        do j=1,nfft
         ca(i,j) = alog( ca(i,j)/ai0 + t) - at
! ----                                 = <- no primary beam int.
!                                           will be in the results !
!                                           omitt this t if primary
!                                           beam is in the data.
        enddo
       enddo
!
!debug
       call pushda(ca(1,1),0.0,dy,nfft)
!debug
! ---- reverse transformation back to scattering angle field ----
       if(iout.gt.1)write(6,*)'starting 2. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
!
! ---- return values:
       do i=1,nfft
        xout(i) = (i-1) * dx
        yout(i) = ca(1,i) * dy * dy * ai0
       enddo
       if(iout.gt.0) then
          write(6,*)'intensity factor i0 =',ai0
          write(6,*)'                     ','=================='
       endif
! ----------------------------------------------------------------------
!
       return
      END


       subroutine mux(xin,yin,nin,t,xmax,nfft,iout,xout,yout,ier)
!      ==========================================================
!  --- small angle multiple scattering evaluator ----
!      input:
!      xin(1..nin) ... : inputvector containig scatteringangles
!                        or equivalent in ascending order
!      yin(1..nin) ... : inputvector containing scattering intensities
!      nin ........... : active length of xin and yin
!      t ............. : effective transmission factor associated
!                        with the small angle scattering
!      xmax .......... : x (scattering vector/angle) range
!                        for the result
!      nfft .......... : number of points of the internal
!                        fouriertransforms = the number of
!                        points in the output vectors.
!                        preferred value: 2**n.
!      iout .......... : variable to control additional printoutput
!                        from this subroutine
!      output:
!      xout(1..nfft).. : output x-values
!      yout(1..nfft).. : output y-values
!      ier ........... : errorindicator (0=ok)
!
!      the output will be the hypothetical scattering data of the
!      same sample using the same apparatus
!      but with any effects of multiple scattering,
!      e.g. no multiple sacttered intensity and transmission
!      reduction.
! ----------------------------------------------------------------------
!
! ---- dimensioning ----------------------------------------------------
       use constants
       use fftwr1
       parameter(pi=3.1415926535897)

       real xin (nin) ,yin (nin)
       real xout(nfft),yout(nfft)
                                                       !! k--> 1 aix !!
! ----------------------------------------------------------------------
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

!------- data preparation ----------------------------------------------
!
! --- fill the rest of the interpolation vector with equidistant data --
!
         dx = xmax / (nfft-1)
         ninter = xmax / dx + 1
         if(ninter.gt.mdim) then
           write(6,*)'not enough space to store intermediate',          &
     &               ' interpolation vector'
           ier = 2
           return
         endif
         jl = 2
         do i=1,ninter
           xxin = (i-1)*dx
           do j=jl,nin
             if(xin(j).gt.xxin) goto 10
           enddo
   10      continue
           jl = j
           if(jl.gt.nin) then
             yinter(i) = 0.0
           else
             p = (xxin-xin(jl-1))/(xin(jl)-xin(jl-1))
             yinter(i) = yin(jl)*p + yin(jl-1)*(1.0-p)
           endif
         enddo
!debug
         call pushda(yinter,0.0,dxin,ninter)
!debug
!
! ----- prepare the 1d-data for 2d-fft -------------------------------
       dx = xmax / (nfft-1)
! period of fft is:
       xl = 2 * nfft * dx
! reciprocal fourier increment:
       dy = 1.0 / xl
! ---- generate the 2d field of y-values ----
       if(iout.gt.1)write(6,*)'generating data field to be transformed'
! ---- init data 2d-intermediate data field for fft :
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
!
!debug
         do i=1,nfft
            yinter(i) = ca(1,i)/dx/dx
         enddo
         call pushda(yinter,0.0,dx,nfft)
!debug
       if(iout.gt.1)write(6,*)'starting 1. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
!debug
       call pushda(ca(1,1),0.0,dy,nfft)
!debug
!
       ai1= ca(1,1)
       at = alog(t)
       xn = -at / ai1
!
! multiple convolution series:
       do i=1,nfft
        do j=1,nfft
         ca(i,j) = exp(xn*ca(i,j)+at)- t
! ----                                 = <- no primary beam int.
!                                           will be in the results !
!                                           omitt this t if primary
!                                           beam is to be generated.
        enddo
       enddo
!debug
       call pushda(ca(1,1),0.0,dy,nfft)
!debug
!
! ---- reverse transformation back to scattering angle field ----
       if(iout.gt.1)write(6,*)'starting 2. fft ....'
       call  rfft2d(nfft,nfft,ca,mdim,ca,mdim)
!
! ---- return values:
       do i=1,nfft
        xout(i) = (i-1) * dx
        yout(i) = ca(1,i) * dy * dy / xn
       enddo
! ----------------------------------------------------------------------
!
       return
      END


       subroutine pushda(y,x0,dx,n)
!      ----------------------------
! --> pushes data y-values to the datafiles

       use cdata
       use constants
       use outlev
       implicit none

       real y,x0,dx
       integer n

       dimension y(n)
       integer :: num=9000
       integer i
! ----------------------------------------------------------------
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
      END



      double precision function tof_omega(irec,ikan)
!     ----------------------------------------------
!  Extract Omega (in Rad/sec) for a inx tof record
!
       use cdata
       use outlev
       use theory
       use selist
       use therrc
       use cfunc
       use PhysicalConstantsPlus
       implicit none

       integer            :: irec 
       integer            :: ikan

       integer            :: i, j, ier
       real*4             :: xh
       double precision   :: ymax, kell, br, flpath, lambda
       double precision   :: veloc0, velocx, ejoul, omega, t0, tx


       call parget('binmusec', xh ,irec ,ier)
       if(ier.eq.0) then
         br = xh*1d-6
       else
         write(6,*)'tof_omega: channelwidth binmusec not found! ',irec
         br = 1
       endif

       call parget('kell    ', xh ,irec ,ier)
       if(ier.eq.0) then
         kell = xh
       else      
         ymax = ywerte(1,irec)
         kell = 1
         do i= 1, nwert(irec)
          if(ywerte(i,irec).gt.ymax) then
             ymax = ywerte(i,irec)
             kell = i
          endif
         enddo
         call  parset('kell    ',sngl(kell),irec)
       endif

       call parget('flpath  ', xh ,irec ,ier)
       if(ier.eq.0) then
         flpath = xh
       else      
         flpath = 4d0
         call  parset('flpath  ',sngl(flpath),irec)
       endif

       call parget('lambda  ', xh ,irec ,ier)
       if(ier.eq.0) then
         lambda = xh*1d-10
       else      
         write(6,*)' >>>>>>>> enter lambda/A >>>>>> '
         read(5,*) lambda
         call  parset('lambda   ',sngl(lambda),irec)
         lambda = lambda*1d-10
       endif

       veloc0 = NeutronVelocity_fromLambda(Lambda)
       t0     = flpath/veloc0
       tx     = (ikan-kell)*br+t0
       velocx = flpath/tx
       ejoul  = 0.5d0*Neutronenmasse*(velocx**2-veloc0**2)
       omega  = 2*Pi*ejoul/Planckkonstante
       tof_omega = omega

       return
       end function tof_omega





      subroutine usrfun(nam,x,nx,ier)
!     -------------------------------

       use cdata
       use outlev
       use theory
       use selist
       use therrc
       use cfunc
       use PhysicalConstantsPlus
       implicit none

      character*20 nam
      real*8 x(*)
      integer ier, nx

      real sum, pvalue, x3
      integer i, iwert, iwert2, itheo, iparn, ibuf, iadr, iaddr, iiadr, inumr, l, myfun
      character*8 pname

      real*8 xh, dx
      logical compare

      double precision :: tof_omega

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

      if(compare(nam,'om ')) then  ! TOF channel frequency
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf = x(nx-1)+0.01
        iwert= x(nx)  +0.01
        x(nx-1) = tof_omega(iwert,ibuf)
        nx = nx-1
        return
      endif

      if(compare(nam,'en ')) then   ! TOF channel energy in meV
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf = x(nx-1)+0.01
        iwert= x(nx)  +0.01
        x(nx-1) = tof_omega(iwert,ibuf)*Planckkonstante/2/Pi/Elektronenladung*1d3
        nx = nx-1
        return
      endif

      if(compare(nam,'sc ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        inumr= NINT(x(nx))
        iadr = 0
        do l=1,nbuf
          if(numor(l).eq.inumr) then
            iadr = l
            goto 22
          endif
        enddo
   22   continue
        x(nx) = iadr+0.0000001
        return
      endif

      if(compare(nam,'nv ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        inumr= NINT(x(nx))
        iadr = nwert(inumr)
        x(nx) = iadr*1.00000001
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

      if(compare(nam,'isel ')) then
        if(nx.lt.1) then
          ier =1
          return
        endif
        inumr= x(nx)+0.001
        iadr = isels(inumr)
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
!        write(6,*) iwert,min(iwert,nwert(ibuf))
!        write(6,*) iwert2,min(iwert2,nwert(ibuf))
        do i=min(iwert,nwert(ibuf)),min(iwert2,nwert(ibuf))
          sum = sum + ywerte(i,ibuf)
        enddo
        x(nx-2) = sum
        nx = nx-2
        return
      endif

      if(compare(nam,'sumyerr ')) then
        if(nx.lt.3) then
          ier =1
          return
        endif
        ibuf  = x(nx-2)+0.01
        iwert = x(nx-1)+0.01
        iwert2= x(nx  )+0.01
        sum = 0
!        write(6,*) iwert,min(iwert,nwert(ibuf))
!        write(6,*) iwert2,min(iwert2,nwert(ibuf))
        do i=min(iwert,nwert(ibuf)),min(iwert2,nwert(ibuf))
          sum = sum + yerror(i,ibuf)**2
        enddo
        x(nx-2) = sqrt(sum)
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
        do i=min(iwert,nwert(ibuf)),min(iwert2,nwert(ibuf))
          sum = sum + xwerte(i,ibuf)
        enddo
        x(nx-2) = sum
        nx = nx-2
        return
      endif

      if(compare(nam,'indxval ')) then   ! index of x-value
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf  = x(nx-1)+0.01
        xh    = x(nx  )

        if(xh.lt.xwerte(1,ibuf) .or. xh.gt.xwerte(nwert(ibuf),ibuf))then
          x(nx-1) = 0
          nx = nx-1
          return
        endif

        do i=1,nwert(ibuf)
          if(xwerte(i,ibuf).ge.xh) then
           x(nx-1) = i
           nx = nx-1
           return
          endif
        enddo
        x(nx-1) = 0
        nx = nx-1
        return
      endif

      if(compare(nam,'intxval ')) then   ! index of x-value
        if(nx.lt.2) then
          ier =1
          return
        endif
        ibuf  = x(nx-1)+0.01
        xh    = x(nx  )

        if(xh.lt.xwerte(1,ibuf) .or. xh.gt.xwerte(nwert(ibuf),ibuf))then
          x(nx-1) = 0
          nx = nx-1
          return
        endif

        do i=1,nwert(ibuf)
          if(xwerte(i,ibuf).ge.xh) then
           if(i.le.1) then
            dx=xwerte(i+1,ibuf)-xwerte(i,ibuf)
            if(dx.ne.0) dx=(xh-xwerte(i,ibuf))/dx
            x(nx-1) = dx
           else
            dx=xwerte(i,ibuf)-xwerte(i-1,ibuf)
            if(dx.ne.0) dx=(xh-xwerte(i,ibuf))/dx
            x(nx-1) = i+dx
           endif
           nx = nx-1
           return
          endif
        enddo
        x(nx-1) = 0
        nx = nx-1
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
  111 continue
      iiadr = x(nx) + 0.001
      call parget (pname,pvalue,iiadr,ier)
      x(nx) = pvalue

      return
      END


      subroutine usrextr(nam,val,ier)
!     -------------------------------
       use cdata
       use outlev
       use selist
       use cfunc
       implicit none

       integer iadda
       common /thiadd/iadda

      character*1 nam(*)
      real*8 val
      integer ier

      character*8 pname
      logical compare
      integer i, ibuf
      real pvalue, sum, sumy, sumyy, xh

      ier = 0
      if(compare(nam,'sel ')) then
        val = isels(1)
        return
      endif
      if(compare(nam,'nsel ')) then
        val = nsel
        return
      endif
      if(compare(nam,'maxx ')) then
        ibuf = iadda
        val = xwerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(xwerte(i,ibuf).gt.val) val = xwerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'minx ')) then
        ibuf = iadda
        val = xwerte(1,ibuf)
        do i=2,nwert(ibuf)
          if(xwerte(i,ibuf).lt.val) val = xwerte(i,ibuf)
        enddo
        return
      endif
      if(compare(nam,'maxy ')) then
        ibuf = iadda
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
!! newnewnewnewnew>
!! center of mass type from subrange limited by x1 x2 fit range
      if(compare(nam,'centerx ')) then
        ibuf = isels(1)
        sumy  = 0
        sum   = 0
        do i=1,nwert(ibuf)
          xh   = xwerte(i,ibuf)
          if(((xh.ge.x1) .or. autox1).and.((xh.le.x2).or.autox2))then
            sum = sum + ywerte(i,ibuf)*xwerte(i,ibuf)
            sumy= sumy+ ywerte(i,ibuf)
          endif
        enddo
        if(sumy.eq.0.0) sumy = 1.0
        val = sum/sumy
        return
      endif

      if(compare(nam,'widthx ')) then
        ibuf = isels(1)
        sumy  = 0
        sumyy = 0
        sum   = 0
        do i=1,nwert(ibuf)
          xh   = xwerte(i,ibuf)
          if(((xh.ge.x1).or.autox1).and.((xh.le.x2).or.autox2))then
            sum = sum + ywerte(i,ibuf)*xwerte(i,ibuf)
            sumy= sumy+ ywerte(i,ibuf)
          endif
        enddo
        if(sumy.eq.0.0) sumy = 1.0
        val = sum/sumy

        sumy  = 0
        sumyy = 0
        sum   = 0
        do i=1,nwert(ibuf)
          xh   = xwerte(i,ibuf)
          if(((xh.ge.x1).or.autox1).and.((xh.le.x2).or.autox2))then
            sum = sum + ywerte(i,ibuf)*(xwerte(i,ibuf)-val)**2
            sumy= sumy+ ywerte(i,ibuf)
          endif
        enddo
        if(sumy.eq.0.0) sumy = 1.0
        if(sum.lt.0) write(6,*)'widthx determination: sum is negative'
        val = sqrt(abs(sum/sumy))
        return
      endif

!! <newnewnewnewnew
      if(compare(nam,'yy   ')) then
        val = yyyy
        return
      endif
      if(compare(nam,'ye   ')) then
        val = yyee
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
  111 continue

      call parget (pname,pvalue,iadda,ier)

        if(ier.ne.0) then
         if(ifits(1).gt.0.and.nsel.gt.0) then
          call parget (pname,pvalue,ifits(1),ier)
         endif
        endif

      val = pvalue
      return
      END


       function wlambda( alam )
!      ------------------------
! --- repraesentiert die wellenlaengenverteilung
	use wlntran
        implicit none
	real*4 arg, alam, wlambda
          arg     = ( (alam-alam0)/dalam )**2
          if(arg.lt.50.e0) then
            wlambda =  exp( -arg )
          else
            wlambda = 0.e0
          endif
       return
      END


       function sofqom( omega )
! --- repraesentiert die streufunktion, omega in s**-1
       use sqtran
       implicit none

       real*4 x, omega, sofqom
       x      = omega * tau
       sofqom = tau /(1.e0+x*x)
       return
      END


       function f(n,x)
! --- integrand ---
       use partran
       implicit none

       integer n
       real*4 x
       dimension x(n)
!
!      x(1) = lambda
!      x(2) = omega
!
! ---- larmorkonstante in rad/s/gauss
! ---- neutronenmasse durch h in 2/m/angstroem

       real*4 :: gamma=18303.33e0, xmh = 2.50607e-4, zpi = 6.283185307e0
       real*4 a,b,f,del, wlambda, sofqom
       a  = gamma*xmh*x(1)/(1.d0+xmh*x(1)*x(1)*x(2)*1d-10/zpi)
       b  = gamma*xmh*x(1)
       del= j0delta**2 + ( (j1echo+j2echo)*0.5e0*cdelta )**2
       f  = ( 1.d0 +  exp(-(a**2+b**2)*del*0.25e0 )*                    &
     &               cos( a*j1echo - b*j2echo )                         &
     &      ) * sofqom(x(2)) * wlambda(x(1))
       return
      END


       function yintp2 (x,y,n,t)
!      ------------------------ einfache tabelleninterpolation ----
       implicit none 
       real*4 x,y,t
       integer n

       dimension x(n),y(n)
       integer i
       real*4 p,yintp2

       if(t.lt.x(1) .or. t.gt.x(n-1) ) then
         yintp2  = 0.e0
         return
       endif
!
       do i=2,n
         if(x(i).gt.t) then
           p = (t-x(i-1))/(x(i)-x(i-1))
           yintp2  = (1.e0-p)*y(i-1)+p*y(i)
           return
         endif
       enddo
       yintp2  = 0.e0
       return
      END

       integer function igrand(an)
!      ------------------
! ---- gaussian approximation to poisson statistics with
!      expectation value = an
!      second moment     = an
!      integer result
!      to mock up counting statistics
! update to use Sim_Add_Noise to fix inconsistencies
! ---------------------------------------------------------------------
       implicit none 
       real*4 an

       real*4 y, z, erfi, rnunf
       double precision :: Sim_Add_Noise
       double precision :: averd

       averd  = an
       igrand = Sim_Add_Noise( averd )

      END


 double precision function Sim_Add_Noise( average ) 
!--------------------------------------------------

  implicit none
  double precision, intent(in)   :: average 
 

  double precision ::  r1, r2, r, f, grand1, grand2, now
  real   rnunf
 
    do 
      r1 = rnunf()*2.d0-1.d0
      r2 = rnunf()*2.d0-1.d0
      r  = r1*r1 +r2*r2
      if(( (r<1.d0) .and. (r .ne. 0.d0) )) exit
    enddo
 
    f = sqrt(-2.d0*log(r)/r)
    grand1 = r1*f
!   grand2 = r2*f;  */ /* verwerfe 2te um keine Korrelationen zu haben */
 
    now    = average + grand1*sqrt(average)
 
    if(now < 0.d0) now = 0.d0
 
    Sim_Add_Noise = now
 
  end function Sim_Add_Noise 






       FUNCTION BSJN(M,X)                                               
!      ==================                                               
!                                                                       
! ----- SPHERICAL BESSELFUNCTION OF ORDER M>=0 -------------------------
!                                                                       
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
        DOUBLE PRECISION :: bsjn, x, eps, f, fmn, fmnm1, fmnp1, fn, fnm1, fnp1, z, zzt
!                                                                       
        INTEGER          :: i, m, maxm, n
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
       DATA EPS/1.0D-6/                                        
!                          
!
       z = x
       if(z.eq.0d0) z = 1d-8
!                                                   
! --- LOOK FOR BAD M VALUE ---                                          
       IF(M.LT.-2) THEN                                                  
        WRITE(6,*)' ERROR: BSJN CALLED M =',M,' < -2 --> BSJN = 0'      
        BSJN = 0.0D0
        RETURN   
       else if(m.eq.-2) then
         BSJN = -cos(z)/(z*z)-sin(z)/z
         return
       else if(m.eq.-1) then          
         BSJN =  cos(z)/z
         return                                                      
       ENDIF                                                            
!                                                                       
       IF(Z.EQ.0.0D0) THEN 
         BSJN = 1.00D0
         IF(M.GT.0) BSJN = 0.0D0
         RETURN                                                         
       ENDIF                                                            
! --- PREPARE FACTORIAL FOR LIMITING FORMULA ---                        
       F = 1.00D0
       DO 101 I=1,2*M+1,2                                               
         F = F * I                                                      
101    CONTINUE                                                         
!                                                                       
! ---- TREAT CASE Z SMALL ----                                          
       ZZT = (Z**M) / F                                                 
       IF(ABS(ZZT).LT.EPS) THEN                                         
!        -----> USE LIMITING FORMULA                                    
           BSJN = (Z**M / F ) * (1.0D0-Z**2/(2*(2*M+3)) +     &
                                 Z**4/(8*(2*M+3)*(2*M+5)) )             
!        WRITE(6,*)' ASYMPTOTIC FORMULA'                                
         RETURN                                                         
       ENDIF                                                            
!                                                                       
! --- TREAT CASES M=0 AND M=1 ---                                       
       IF(M.LT.2) THEN                                                  
        IF(M.EQ.0) THEN                                                 
          BSJN = SIN(Z) / Z                                             
          RETURN                                                        
        ENDIF                                                           
        IF(M.EQ.1) THEN                                                 
          BSJN = ( SIN(Z) / Z  - COS(Z) ) / Z                           
          RETURN                                                        
        ENDIF                                                           
      ENDIF                                                             
!                                                                       
! --- TREAT CASES M >= 2 ---                                            
!                                                                       
      FNM1 = 1.00D0/Z   
      FN   = FNM1**2                                                    
      FMN  = 0.0D0
      FMNP1= FNM1                                                       
      MAXM = M - 1                                                      
!                                                                       
      DO 1 N=1,MAXM                                                     
        FNP1   = (2*N+1) * FN/Z - FNM1                                  
        FMNM1  = (1-2*N) * FMN/Z- FMNP1                                 
        FNM1   = FN                                                     
        FN     = FNP1                                                   
        FMNP1  = FMN                                                    
        FMN    = FMNM1                                                  
1     CONTINUE                                                          
      FMNM1 = (1-2*M) * FMN/Z - FMNP1                                   
      BSJN  = FN*SIN(Z) - (-1)**M * FMNM1*COS(Z)                        
!                                                                       
      RETURN                                                            
      END  function bsjn                                                



      function radial_bsjn2_integral(n,q,r1,r2)
!     -----------------------------------------
!
!     Integral over a shell from r1..r2 of bsjn(n,q*r)**2  as is needed in the rotational diffusion
!     modelling for a rotating protein in incoherent approximation
!     The integral is normalized to the shell volume
!
!     i.e. integral_[r1..r2](bsjn(x,q*r)**2 *r**2 * dr) / integral_[r1..r2](r**2 *dr) 
!
      implicit none
      double precision              :: radial_bsjn2_integral
      double precision, intent(in)  :: q, r1, r2
      integer         , intent(in)  :: n
    
      double precision              :: bsjn

      radial_bsjn2_integral = -(0.5d0)*R1**3*bsjn(n  , R1*q)**2                 &
                              +(0.5d0)*R1**3*bsjn(n-1, R1*q)  *bsjn(n+1, R1*q)  &
                              +(0.5d0)*R2**3*bsjn(n  , R2*q)**2                 &
                              -(0.5d0)*R2**3*bsjn(n-1, R2*q)  *bsjn(n+1, R2*q) 


!     normalizing to volume
      
      radial_bsjn2_integral =  radial_bsjn2_integral*3/(r2**3-r1**3)

      return
      end function  radial_bsjn2_integral


!                                                                       
!                                                                       
!                                                                       
! ---- SCATTERING LAW OF A DIFFUSION INSIDE A SPHERE ----    
!      time domain: intermediate scattering function           
!      ACCORDING TO F.VOLINO & A.J. DIANOUX MOL.PHYS (1980),VOL41,271   
!                                                                       
       FUNCTION DINSPH_t(Q,t,D,A,LX,ZX)                                  
!      --------------------------------                                  
!                                                                       
!   Q   = MOMENTUM TRANSFER                                             
!   OM  = FREQUENCY (ENERGY TRANSFER)                                   
!   D   = DIFFUSION CONSTANT                                            
!   A   = RADIUS OF SPHERE                                              
!   LX  = MAX L TO BE USED                                              
!   ZX  = MAX VALUE OF ZERO TO BE USED                                  
!                                                                       
!   ZEROS(I) = VALUE OF I-TH ZERO = XNL                                 
!   LOFZ (I) = L-VALUE FOR THIS ZERO                                    
!   NOFZ (I) = NO OF THIS ZERO                                          
!   NZERO    = TOTAL NO OF ZEROS                                        
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!        INCLUDE  'th1_vdzero_comdef.f'

        DOUBLE PRECISION        zeros
        INTEGER*4               lofz, nofz, nzero, maxz
        PARAMETER               (maxz = 5000)
        DIMENSION               zeros(maxz), lofz(maxz), nofz(maxz)
!                                                                           !
        COMMON /vdzero/         zeros, lofz, nofz, nzero

        DOUBLE PRECISION        zmx
        INTEGER*4               lmx
!                                                                           !
        COMMON /dinbuf/         zmx, lmx


!        INCLUDE  'th1_dinbuf_comdef.f'
!        INCLUDE 'natur_konstanten_defs.f'
!                                                                       
        DOUBLE PRECISION a, aq, d, da, delta, dinsph_t, dom, om, oml, &
                         q, qa, sum, x, xx, zx, t
        DOUBLE PRECISION bsjn

        double precision, parameter :: pi=3.141592654d0
!                                                                       
        INTEGER*4        i, l, lx, n
!                                                                       
!        DATA             OML/1.0D50/
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
       DELTA = 1.0D-6
!      -------------> LIMIT FOR AVOIDING SINGULARITY                    
! ---- OML IS NEEDED TO IDENTIFY THE ZERO POSITION & GET THE FACTOR FOR 
!      THE DELTA FUNCTION                                               
!       IF(OML.GT.OM) OML = OM                                           
!                    --------> RESETS OML IF A NEW DETECTOR IS TAKEN    
!                                                                       
! ---- GENERATE ZEROS IF NECESSARY ----                                 
!                                                                       
       IF(LX.NE.LMX.OR.DABS(ZX-ZMX).GT.0.2D0.OR.NZERO.EQ.0) THEN        
         LMX = LX                                                       
         ZMX = ZX                                                       
         CALL VDCOEF(ZMX,LMX)                                           
         WRITE(6,*)' .....',NZERO,' NEW VDCOEFFS GENERATED'             
       ENDIF                                                            
!                                                                       
! ---- DO THE SUMMATION ----                                            
!                                                                       
       QA = Q*A                                                         
       DA = D / A**2                                                    
       SUM= 0.0D0
!                                                                       
       DO 10 I=1,NZERO                                                  
        N = NOFZ(I)                                                     
        L = LOFZ(I)  
        if(l.eq.0 .and.n.eq.0) cycle                                                   
        X = ZEROS(I)                                                    
        XX= X**2                                                        
        IF(DABS(QA-ZEROS(I)).GT.DELTA) THEN                             
          AQ = ( QA * BSJN(L+1,QA) - L * BSJN(L,QA) ) / (QA**2-XX)
          AQ = AQ**2 * 6.0D0 * XX / ( XX - L*(L+1) )
        ELSE                                                            
          AQ = 1.5D0 * BSJN(L,X)**2 * (XX - L*(L+1) ) / XX 
        ENDIF                                                           
!                                                                       
!        SUM = SUM + (2*L+1)*AQ *    XX*DA / ( (XX*DA)**2 + OM**2 )       ! hier Zeitfunktion einsetzen  
        SUM = SUM + (2*L+1)*AQ *    Pi*exp(-XX*DA*abs(t))
!
!   write(6,*)'1: ',i,n,l,qa,aq,xx,sum                                                                       
10     CONTINUE                                                         
       SUM = SUM / PI  
!   write(6,*)'1a: ',sum                                                 
!                                                                       
! ---- ADD THE DELTA FUNCTION ----                                      
!                                                                       
!       IF(OM*OML.LT.0.0D0.OR.OM.EQ.0.0D0) THEN 
!        DOM = DABS(OM-OML)                                              
!        SUM = SUM + (3.0D0*BSJN(1,QA)/QA)**2 / DOM 
!       ENDIF                                                            
!
        SUM = SUM + (3.0D0*BSJN(1,QA)/QA)**2 
!    write(6,*)'2: ',qa,aq,xx,sum   
       DINSPH_t = SUM                           
!      --------------  
!   write(6,*)'2a: ',sum                                  
       RETURN                                
       END                                  



       DOUBLE PRECISION FUNCTION BLJ(Z,L)  
!      ==================================
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
        DOUBLE PRECISION z 
        DOUBLE PRECISION bsjn
!                                                                       
        INTEGER*4       l
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
       IF(L.EQ.0) THEN                                                  
         BLJ = BSJN(1,Z)  
         RETURN                                                         
       ENDIF                                                            
!                                                                       
       BLJ = L*BSJN(L,Z) - Z * BSJN(L+1,Z)    
!      -------------------------------------
       RETURN                                                           
       END                                                              
!                                                                       
!                                                                       
!                                                                       

!                                                                       
! ---- DETERMINE COEFFICIENTS XNL FOR VOLINO DIANOUX THEORY FOR         
!      THE SCATTERING OF DIFFUSION IN A SPHERICAL POTENTIAL WELL ----   
!                                                                       
       SUBROUTINE VDCOEF(ZMX,LMX)                                       
!      --------------------------                                       
!                                                                       
! ---- ZMX = UPPER LIMIT FOR ZEROES ----                                
!      LMX = UPPER LIMIT FOR L                                          
!                                                                       
!                                                                       
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
	IMPLICIT NONE
!                                                                       
!       INCLUDE  'outlev_comdef.f'
!       INCLUDE  'th1_vdzero_comdef.f'
        DOUBLE PRECISION        zeros
        INTEGER*4               lofz, nofz, nzero, maxz
        PARAMETER               (maxz = 5000)
        DIMENSION               zeros(maxz), lofz(maxz), nofz(maxz)
!                                                                           !
        COMMON /vdzero/         zeros, lofz, nofz, nzero
!                                                                       
        DOUBLE PRECISION dz, dz0, eps, f, f0, f1, z0, zmx
        DOUBLE PRECISION blj
!                                                                       
        INTEGER*4        i, j, l, lmx, maxit, ncnt, nl
        integer  :: ierrs
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
!                                                                       
!                                                                       
! ---- SET PARAMETER-VALUES ---                                         
       MAXIT = 50                                                       
!      -----------> MAX NO. OF NEWTON ITERATIONS                        
       DZ    = 0.001D0
!      -----------> SETP FOR NEWTON DIFFERENTIATION                     
       EPS   = 1.0D-9  
!      -----------> REQUIRED ACCURACY FOR NEWTON ITERATION              
       DZ0   = 0.2D0
!      -----------> INCREMENT FOR EXHAUSTIVE ZERO SEARCH                
!                                                                       
       WRITE(6,*)' VDCOEF GENERATING COEFFICIENTS FOR ',LMX,ZMX,' ...'  
!                                                                       
       NL = ZMX / DZ0                                                   
       NZERO = 0                                                        
!                                                                       
       DO 5 L = 0,LMX                                                   
         NCNT = 0                                                       
         DO 10 J = 1,NL                                                 
          Z0 = (J-1) * DZ0                                              
          F  = BLJ(Z0,L) * BLJ(Z0-DZ0,L)                                
          IF(F.LT.0) THEN                                               
!                                                                       
! ---- NEWTON-ITERATION ----                                            
            DO 100 I=1,MAXIT                                            
             F0 = BLJ(Z0,L)                                             
             F1 = BLJ(Z0+DZ,L)                                          
             Z0 = Z0 - F0 * DZ / (F1-F0)                                
!            WRITE(6,*)I,F0,Z0                                          
             IF(DABS(F0).LT.EPS) GOTO 101                               
100         CONTINUE                                                    
            WRITE(6,*)' WARNING BAD CONVERGENCE , RESIDUAL=',F0         
            IERRS = 55                                                  
!                                                                       
101         CONTINUE                                                    
            NZERO = NZERO + 1                                           
            IF(NZERO.GT.MAXZ) THEN                                      
              WRITE(6,*)' NO OF ZEROS EXCEEDS MAX DIMENSION = ',MAXZ    
              IERRS = 777                                               
              NZERO = MAXZ                                              
            ENDIF                                                       
!                                                                       
            ZEROS(NZERO) = Z0                                           
            LOFZ(NZERO)  = L                                            
            NOFZ(NZERO)  = NCNT                                         
!                                                                       
!       WRITE(6,*)L,NCNT, Z0                          
            NCNT = NCNT + 1                                             
          ENDIF                                                         
10       CONTINUE                                                       
5      CONTINUE                                                         
!                                                                       
       RETURN                                                           
!                                                                       
       END                                                              








!..............................................................................................
!. functions pertaining the detector sensitivity: parameters for LET ..........................
!..............................................................................................

  real function detector_sensitivity(lambdaA,lenfrac,p,dtube)
! -----------------------------------------------------------
  implicit none
  real, intent(in)    :: lambdaA             ! wavelength in Angstroems
  real, intent(in)    :: lenfrac             ! here the projection strech factor given by distance(i)/distance0
  real, intent(in)    :: p                   ! detector pressure in bar 
  real, intent(in)    :: dtube               ! detector diameter in m
  

  real, parameter     :: sigmaA_He3 = 5333e-28
  real, parameter     :: sigma_Iron = 20.0                ! absorption per meter (estimated)
  real, parameter     :: effective_tube_thickness = 0.001 ! default assumption 1mm
  real, parameter     :: NA         = 6.022045e23
  real, parameter     :: Vmol       = 0.02241383
  real, parameter     :: lambda0    = 1.798
  double precision    :: alpha
  double precision    :: wcyldet

  alpha = NA/Vmol * p * dtube * lenfrac * sigmaA_He3 * lambdaA/lambda0

  detector_sensitivity = max(0.01, wcyldet(alpha) * exp(-lenfrac*effective_tube_thickness*sigma_Iron) )
!                                                                                         -> here we may needs some sort of lambda dependence
  
!  write(6,*)' p = ',p
!  write(6,*)' d = ',dtube
!  write(6,*)' a=  ',alpha
!  write(6,*)' w=  ',wcyldet(alpha) 

  end function detector_sensitivity



  double precision function wcyldet(alpha)
! ----------------------------------------
!
! series expansion of integral to estimate the detection efficiency of a cylindrical
! (elliptical) detector tube
! alpha is = Diameter x Sigma-abs
! alpha is   tubediameter * Sigma(He3) * (distance(j)/distance0)
!
! Sigma(He3) = numberdensity * absorption-crossection * lambda/lambda0
!                  = Na / (Vmol) * (p/p0) *  5333barn * lambda/1.798A
!
! the expression returned is an approximation to 
! 
! W := (1/2)*(int(1-exp(-alpha*sqrt(1-x^2)), x = -1 .. 1))
! 
! obtained by the taylor expansion of the integral
!
! the exact result would be:
!
! W = 1-(-(1/2)*alpha*Pi*BesselI(1, alpha)+(1/2)*alpha*Pi*StruveL(1, alpha)+alpha)/alpha
!

    implicit none
    double precision, intent(in)   :: alpha
    double precision               :: t1, w


     t1 = exp(-alpha)
     w  = 1-t1+alpha*(-0.148562771061958D18/0.6977792621150208D18*t1+ &
          alpha*(-0.3226622737612121D17/0.6977792621150208D18*t1+     &
          alpha*(-72378131847533.D0/0.872224077643776D16*t1+          &
          alpha*(-0.2549005820162231D16/0.2093337786345062D19*t1+     &
          alpha*(-102518622323707.D0/0.6977792621150208D18*t1+        &
          alpha*(-32516961685057.D0/0.2220206743093248D19*t1+         &
          alpha*(-77335696019.D0/0.64044425281536D17*t1+              &
          alpha*(-44134000963.D0/0.5379731723649024D18*t1+            &
          alpha*(-1439093561.D0/0.316454807273472D18*t1+              &
          alpha*(-30269161.D0/0.149899645550592D18*t1+                &
          alpha*(-4097321.D0/0.58889146466304D18*t1+                  &
          alpha*(-18013.D0/0.10241590689792D18*t1+                    &
          alpha*(-467.D0/0.1597688147607552D18*t1-                    &
          alpha*t1/0.414215445676032D17                   )))))))))))))

          wcyldet = w

   end function wcyldet   


