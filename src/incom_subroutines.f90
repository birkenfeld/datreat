!***********************************************************************
!       incom-system  subroutines                                      *
!       ------------                                                   *
!       by  michael monkenbusch                                        *
!       institut fuer festkoerperforschung des                         *
!       forschungszentrum juelich, d-5170 juelich, po-box 1913, frg    *
!       iff115 @ djukfa11                                              *
!                                                                      *
!***********************************************************************
!cc+++ noch optimierungsbeduerftig !!! ++++++++++++++++++<<<!!!!!!
       subroutine creplace(a,b,c,s)
!---------------------------------------------------
! ersetzten von string b durch c in string a
! strings werden durch den character s begrenzt !
!----------------------------------------------------
       implicit none
       integer, parameter :: mclen=1024
       character*(mclen) a,b,c,d
       character*1 s

       integer i, j, l, jj, ld, istart, iend, len
!
! ---  look for occurence of substring b in a ---
       i = 1
       d = a
       len = 0
       do l=1,mclen
         if(a(l:l).eq.s) then
           len = l-1
           goto 1
         endif
       enddo
    1  continue
         a = d
         if(a(i:i).eq.b(1:1)) then
           istart = i
           do j=2,mclen
             if(b(j:j).eq.s  ) goto 10
             l=i-1+j
             if(b(j:j).ne.a(l:l)) then
               istart = 0
               goto 10
             endif
           enddo
   10     continue
          if(istart.ne.0) then
            iend = l
            do jj=1,mclen
              ld= istart-1+jj
              if(c(jj:jj).eq.s  ) goto 20
              d(ld:ld)= c(jj:jj)
            enddo
            i = ld
   20       continue
            do jj=iend+1,mclen
             if(l.gt.mclen) goto 30
             d(ld:ld)= a(jj:jj)
             if(a(jj:jj).eq.s  ) goto 30
             ld = ld+1
            enddo
   30       continue
          endif
         endif
       i = i+1
       if(i.le.mclen) then
        if( a(i:i).ne.s  ) goto 1
       endif

       do l=1,len
         if(d(l:l).eq.s)then
            do i=l,len
              d(i:i)=' '
            enddo
         endif
       enddo


       a = d

       return
      END

      subroutine scan(ctext,r,ierr)
! ----------------------------------------------------------------------
!  scan - scannen einer real-zahl (version fuer ibm pc prof. fortran)
!  autor: g. egerer, datum: 26. 9.85, letzte aenderung: 31. 7.87
! ----------------------------------------------------------------------
!                           * eingabeparameter:
!                           * ctext  - textzeile (64 zeichen)
!                           * ausgabeparameter:
!                           * r      - real-zahl
!                           * ierr   - fehlercode
      implicit none 
      real*8 r
      integer ierr
      integer, parameter :: ctxtle = 64
      character*(ctxtle) ctext
      character*10       form


!                           * benamte programmkonstanten:
      integer, parameter :: cnil = -32766

!                           * programmvariablen:
      integer i,j
      integer   icur, iscntb(9,0:7), istate, isymbl


!                           * initialisierungen:
! tabelle zur syntaxdefinition:
!
! symbol   s0   " "   "-"   "+"   "."   "e"   "d"   ziff
! ------+--------------------------------------------------
! zu-   |
! stand |
! --> 1 |         1     2     2     3                  8
!     2 |                           3                  8
!     3 |                                              9
!     4 |               5     5                        6
!     5 |                                              6
! end-  |
! zust. |
!     6 |         7                                    7
!     7 |         7
!     8 |         7                 9     4     4      8
!     9 |         7                       4     4      9
!
! s0 (symbolklasse 0): enthaelt alle symbole, die nicht gesondert aufge-
!                      fuehrt sind
! ziff               : enthaelt alle ziffern von "0" - "9"
!
      data ((iscntb(i,j), j = 0,7), i = 1,9)                            &
     &   / cnil,    1,    2,    2,    3, cnil, cnil,    8,              &
     &     cnil, cnil, cnil, cnil,    3, cnil, cnil,    8,              &
     &     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    9,              &
     &     cnil, cnil,    5,    5, cnil, cnil, cnil,    6,              &
     &     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    6,              &
     &     cnil,    7, cnil, cnil, cnil, cnil, cnil,    7,              &
     &     cnil,    7, cnil, cnil, cnil, cnil, cnil, cnil,              &
     &     cnil,    7, cnil, cnil,    9,    4,    4,    8,              &
     &     cnil,    7, cnil, cnil, cnil,    4,    4,    9 /
!                           * ------------------------------------------

      icur = 1
      istate = 1
 1100 continue
      if ( icur .gt. ctxtle ) then
         if ( istate .ge. 6 ) then
!                           * istate ist ein endzustand
            ierr = 0
            form = '(BN,Eww.0)'
            write (form(6:7),'(i2)') ctxtle
            read(ctext,form) r
         else
!                           * istate ist kein endzustand
            ierr = icur
         endif
      else
!                           * ermittle index des aktuellen symbols
         isymbl = index(' -+.EDed0123456789',ctext(icur:icur))
         if (isymbl .ge. 7 ) then
            if (isymbl .ge. 9 ) then
               isymbl = 7
            else
              isymbl = isymbl - 2
            endif
         endif
         istate = iscntb(istate,isymbl)
         if ( istate .ge. 0 ) then
            icur = icur+1
            go to 1100
         endif
         ierr = icur
      endif
      END
!
!
!*ds
!*ed
      subroutine parse( zeile )
!     =========================
! ----------------------------------------------------------------------
! zerlegung der zeile nach incom manier
! die daten der letzten eingabedekodierung werden dabei ueberschrieben
! eine ggf. vorhandene restzeile bleibt erhalten
! restzeilen (; ...) beim parsen werden vollstaendig ignoriert !
! ----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       implicit none 

       character*1024 zeile, reslbf
       character*8   cmd
       integer ioldbf
!-----------------------------------------------------------------------
!
       reslbf = reslin
       ioldbf = ioldc
       reslin = zeile
!       call lowcase(reslin,)
       ioldc  = 1
       call incom (cmd)

       reslin = reslbf
       ioldc  = ioldbf

      return
      END
!
!
      character*1024 function X_datpath()
!     =================================
      use icpathes
      implicit none

      X_datpath(1:1024) = data_path

      return
      end
!
      character*1024 function savepath()
!     =================================
      use icpathes
      implicit none

      savepath = save_path

      return
      end

      character*1024 function makropath()
!     =================================
      use icpathes
      implicit none

      makropath = makro_path

      return
      end

      character*1024 function homepath()
!     =================================
      use icpathes
      implicit none

      homepath = home

      return
      end

	character*1024 function get_argvals(ii)
!     =================================
        use cmargs
        implicit none
        integer ii

! --- variables for makro-parameter-passing ---  see incom for details
      get_argvals =argvals(ii)(:)

      return
      end
!

      integer function iout()
!     ===============
      use xoutxx
      implicit none
         iout = iot
      return
      END

      integer function ierrs()
!     ===============
      use xoutxx
      implicit none
         ierrs = ierrr
      return
      END
!
!
!*ds
!*ed
      subroutine errsig(ierr,say)
!     ================================
      use xoutxx
      implicit none
      integer ierr, laenge, lsay
      character*128 say, sayit
! ----------------------------------------------------------------------
!  error signalisierung
! ----------------------------------------------------------------------
      ierrr = ierr
      lsay = laenge(say,80,'$')
      sayit = say(1:lsay)
      write(6,*)'error:',ierr,' ',trim(sayit)
      return
      END


      subroutine clean
!     ================
      use xoutxx
      implicit none
          if(ierrr.ne.0) then
            write(6,*)'error condition:',ierrr,' cleaned'
            ierrr = 0
          endif
      return
      END


      subroutine pushn(pname)
!-----------------------------------------------------------------------
!  wert eintrag in den namensstack
!-----------------------------------------------------------------------
       use cincom
       use cincoc
       use constants
       implicit none 

       character*8 pname
       real*8 rx

       if(inames.ge.minc) then
         call errsig(-1,'pushn failed, too many items! '//pname//' $')
         return
       endif
!
       inames = inames + 1
       vname(inames) = pname

        return

      entry pushr(rx)     ! not used as far qas i can see ->RB
!     ===============
!----------------------------------------------------------------------
!     zufuegen eines real-wertes auf den stack
!----------------------------------------------------------------------
      if(ipars.ge.minc) then
        call errsig(-1,'pushr failed, too many items ! $')
        return
      endif

      ipars = ipars+1
      rpar(ipars)   =  rx
      if(inpar(inames).eq.0) inapa(inames) =  ipars
      inpar(inames) =  inpar(inames) + 1
      iparn(ipars)  =  inames
!
      END
!
!
!*ds
!*ed
      function getval(pname,defval,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants

       implicit real*8 (a-h,o-z)
       character*8 pname
!-----------------------------------------------------------------------
       inew      = 0
       wert      = defval
       getval    = defval
       call setudf(pname//' ',wert,ier)
       if(inames.le.0) return
!
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            if(inpar(i).gt.0) wert = rpar(inapa(i))
            inew = i
            lstpar = inapa(i)
            lstnam = i
            goto 20
         endif
   10  continue
   20  continue
!
      getval = wert

!c-setudf--
!c-setudf--

      return
!
      END
!
!
!*ds
!*ed
      function valnxt(defval,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       inew      = 0
       valnxt    = defval
       if(ipars.le.lstpar) return
!
      valnxt = rpar(lstpar+1)
      inew   = lstpar + 1
      lstpar = inew

      return
!
      END
!*ds
!*ed
      subroutine get1(p1)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)


       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       lstpar = 1
       return
       entry get2(p1,p2)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       lstpar = 2
       return
       entry get3(p1,p2,p3)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       lstpar = 3
       return
       entry get4(p1,p2,p3,p4)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       lstpar = 4
       return
       entry get5(p1,p2,p3,p4,p5)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       if(ipars.ge.5) p5        = rpar(5)
       lstpar = 5
       return
       entry get6(p1,p2,p3,p4,p5,p6)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       if(ipars.ge.5) p5        = rpar(5)
       if(ipars.ge.6) p6        = rpar(6)
       lstpar = 6
       return
       entry get7(p1,p2,p3,p4,p5,p6,p7)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       if(ipars.ge.5) p5        = rpar(5)
       if(ipars.ge.6) p6        = rpar(6)
       if(ipars.ge.7) p7        = rpar(7)
       lstpar = 7
       return
       entry get8(p1,p2,p3,p4,p5,p6,p7,p8)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       if(ipars.ge.5) p5        = rpar(5)
       if(ipars.ge.6) p6        = rpar(6)
       if(ipars.ge.7) p7        = rpar(7)
       if(ipars.ge.8) p8        = rpar(8)
       lstpar = 8
       return
       entry get9(p1,p2,p3,p4,p5,p6,p7,p8,p9)
       if(inames.ne.0) return
       if(ipars.ge.1) p1        = rpar(1)
       if(ipars.ge.2) p2        = rpar(2)
       if(ipars.ge.3) p3        = rpar(3)
       if(ipars.ge.4) p4        = rpar(4)
       if(ipars.ge.5) p5        = rpar(5)
       if(ipars.ge.6) p6        = rpar(6)
       if(ipars.ge.7) p7        = rpar(7)
       if(ipars.ge.8) p8        = rpar(8)
       if(ipars.ge.9) p9        = rpar(9)
       lstpar = 9
       return
!
      END
!*ds
!*ed
      function intval(pname,idef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  integer
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*8 pname
!-----------------------------------------------------------------------
       inew      = 0
       intval    = idef
       call setudf(pname//' ',wert,ier)
       if(inames.le.0) return
!
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            intval = rpar(inapa(i))*1.0000000001d0
            inew = i
            lstpar = inapa(i)
            lstnam = i
            goto 20
         endif
   10  continue
   20  continue
!
!c-setudf--
      wert = intval
!c-setudf--
      return
!
      END
!*ds
!*ed
      function intnxt(idef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  integer
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       inew      = 0
       intnxt    = idef
       if(ipars.le.lstpar) return
!
       intnxt = rpar(lstpar+1)*1.0000000001d0
       inew   = lstpar+1
       lstpar = lstpar

      return
!
      END
!*ds
!*ed
      character*8 function chrval(pname,cdef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  name (character)
!-----------------------------------------------------------------------
!
!ray -------------------------------                                    \
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*8 pname, cdef
!-----------------------------------------------------------------------
!
       inew      = 0
       chrval    = cdef
       if(inames.le.0) return
!
       do 10 i=1,inames-1
         if(vname(i).eq.pname) then
            chrval = vname(i+1)
            inew = i
            lstnam = i+1
            goto 20
         endif
   10  continue
   20  continue
!
      return
!
      END
!*ds
!*ed
      character*8 function chrnxt(cdef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  name (character)
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*8 cdef
!-----------------------------------------------------------------------
!
       inew      = 0
       chrnxt    = cdef
       if(inames.le.lstnam) return

       chrnxt = vname(lstnam+1)
       inew   = lstnam + 1
       lstnam = inew

      return
!
      END
      character*8 function vnamef(iaddr)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  name (character)
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       vnamef = vname(iaddr)
       lstnam = iaddr
!
      return
!
      END
!*ds
!*ed
      subroutine settit(newtit)
!-----------------------------------------------------------------------
!  setzen des titels
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*80 newtit
!-----------------------------------------------------------------------
!
       ltit  = laenge(newtit,80,'$')
       title = newtit(1:ltit)
!
      return
!
      END
!*ds
!*ed
      character*80 function titlef()
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  titel
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       titlef = title(1:80)
!
      return
!
      END
!*ds
!*ed
      function rparf(iaddr)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  real
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       rparf = rpar(iaddr)
       lstpar= iaddr
!
      return
!
      END

      function ifound(pname)
!-----------------------------------------------------------------------
!  logische optionserkennung  ifound wird=fundstelle sonst 0
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*8 pname
!-----------------------------------------------------------------------
!
       ifound    = 0
       if(inames.le.0) return
!
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            ifound=i
            lstnam=i
            goto 20
         endif
   10  continue
   20  continue
!
      return
!
      END
!*ds
!*ed
      logical function found(pname)
!-----------------------------------------------------------------------
!  logische optionserkennung
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       character*8 pname

      found = ifound(pname) .ne. 0
      return

      END
!
!
!*ds
!*ed
      logical function folgt(popt,pname)
!-----------------------------------------------------------------------
!  logische optionserkennung, wahr wenn popt pname folgt
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       character*8 popt, pname
!-----------------------------------------------------------------------
!
       folgt = .false.
       j     = ifound(pname)
       if(j.gt.0 .and. j.lt.inames) then
          folgt = popt .eq. vname(j+1)
       endif
       if(folgt) lstnam = j+1
!
      return
!
      END
!*ds
!*ed
      function intvno(ipnum,idef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  integer nach addr
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       inew      = 0
       intvno    = idef
       if(ipars.lt.ipnum) return
!
       intvno = rpar(ipnum) * 1.0000000001d0
       inew      = ipnum
       lstpar    = ipnum

      return
!
      END
!*ds
!*ed
      function getvno(ipnum,adef,inew)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  real    nach addr
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use constants
       implicit real*8 (a-h,o-z)

       inew      = 0
       getvno    = adef
       if(ipars.lt.ipnum) return
!
       getvno = rpar(ipnum)
       inew   = ipnum
       lstpar = ipnum

      return
!
      END
!*ds
!*ed
      subroutine getvec(pname,def1,def2,def3,vec ,inew)
!-----------------------------------------------------------------------
!  vektor-wert extraktion aus dem incom parameterstack mit rotation
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use xroxxx
       use constants
       implicit real*8 (a-h,o-z)
!-----------------------------------------------------------------------
       character*16 pnamc
       character*8 pname
       dimension vec(3)
!-----------------------------------------------------------------------
!
       inew      = 0
       vec(1)    = def1
       vec(2)    = def2
       vec(3)    = def3
       call unrot(vec,xyorig,rotvec)
       call cappend(pname//' ','.1 ',pnamc)
       call setudf(pnamc//' ',vec(1),ier)
       call cappend(pname//' ','.2 ',pnamc)
       call setudf(pnamc//' ',vec(2),ier)
       call cappend(pname//' ','.1 ',pnamc)
       call setudf(pnamc//' ',vec(1),ier)
       call cappend(pname//' ','.3 ',pnamc)
       call setudf(pnamc//' ',vec(3),ier)
       if(inames.le.0) goto 20
!
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            if(inpar(i).gt.0) then
               vec(1) = rpar(inapa(i))
               lstpar = inapa(i)
            endif
            if(inpar(i).gt.1) then
               vec(2) = rpar(inapa(i)+1)
               lstpar = inapa(i)+1
            endif
            if(inpar(i).gt.2) then
               vec(3) = rpar(inapa(i)+2)
               lstpar = inapa(i)+2
            endif
            lstnam = i
            inew = i

            goto 20
         endif
   10  continue
   20  continue
!
!c-setudf--
      call cappend(pname//' ','.1 ',pnamc)
      call setudf(pnamc//' ',vec(1),ier)
      call cappend(pname//' ','.2 ',pnamc)
      call setudf(pnamc//' ',vec(2),ier)
      call cappend(pname//' ','.1 ',pnamc)
      call setudf(pnamc//' ',vec(1),ier)
      call cappend(pname//' ','.3 ',pnamc)
      call setudf(pnamc//' ',vec(3),ier)
!c-setudf--

       call rotavc(vec,xyorig,rotvec)
      return
      END
!
!
!*ds
!*ed
      subroutine getnva(pname ,pardef ,parout, np, inew)
!-----------------------------------------------------------------------
!  extraktion meherere einem namen zugeordneter werte
!-----------------------------------------------------------------------
!
!ray -------------------------------
       use cincom
       use cincoc
       use xroxxx
       use constants
       implicit real*8 (a-h,o-z)
!ray -------------------------------
!-----------------------------------------------------------------------
       character*8 pname
       character*8 napp, pnamc*16
       dimension pardef(np), parout(np)
!-----------------------------------------------------------------------
!
       inew      = 0
       do 1 i=1,np
         parout(i) = pardef(i)
    1  continue
       if(inames.le.0) goto 20
!
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            napr = inpar(i)
            if(napr.gt.np) napr = np
            if(napr.gt.0) then
            do 11 j=1,inpar(i)
               parout(j) = rpar(inapa(i)-1+j)
   11       continue
            inew = napr
            lstpar = inapa(i)+inpar(i)-1
            lstnam = i
            endif
            goto 20
         endif
   10  continue
   20  continue

       do 30 i=1,np
        write(napp,'(i3,5x)')100+i
        napp(1:1)='.'
        call cappend(pname//' ',napp,pnamc)
        call setudf(pnamc//' ',parout(i),ier)
   30  continue

!
      return
      END
!*ds
!*ed
      subroutine rotavc( x, xyorig, rotvec )
!     ======================================
!     unterwirft den koordinatenvector x einer drehung um xyorigin
!     um den winkel rotang
!-----------------------------------------------------------------------
      implicit real*8   (a-h,o-z)
!
      dimension  x(3), y(3)
      dimension  xyorig(3), rotvec(3), r(3,3)
      data    degtor /0.017453292d0/
!
      y(1) = x(1) - xyorig(1)
      y(2) = x(2) - xyorig(2)
      y(3) = x(3) - xyorig(3)
      call rotax(rotvec, degtor, r, angdeg )
      call matvec( r, y, y)
      x(1) = y(1) + xyorig(1)
      x(2) = y(2) + xyorig(2)
      x(3) = y(3) + xyorig(3)
!
      return
      END
!*ds
!*ed
      subroutine unrot( x, xyorig, rotvec )
!     =====================================
!     unterwirft den koordinatenvector x einer drehung um xyorigin
!     um die drehung         rotang**-1
!-----------------------------------------------------------------------
      implicit real*8   (a-h,o-z)
!
      dimension x(3), y(3)
      dimension  xyorig(3), rotvec(3), r(3,3)
      data    degtor /0.017453292d0/
!
      y(1) = x(1) - xyorig(1)
      y(2) = x(2) - xyorig(2)
      y(3) = x(3) - xyorig(3)

      call rotax(rotvec, degtor, r, angdeg )
      call mattxp( r, r)
      call matvec( r, y, y)

      x(1) = y(1) + xyorig(1)
      x(2) = y(2) + xyorig(2)
      x(3) = y(3) + xyorig(3)
!
      return
      END
!
!
!*ds
!*ed
       subroutine rotax( ax, scale, rmat, angdeg )
!      ===========================================
!                       -->  ---->  <---   <---
!-----------------------------------------------------------------------
! drehmatrix rmat wird als drehung um die achse ax berechnet
! die laenge wird zur bestimmung des drehwinkels mit scale multipliziert
!-----------------------------------------------------------------------

       implicit real*8 (a-h,o-z)
       dimension ax(3), rmat(3,3)
       dimension en(3)

!-----------------------------------------------------------------------

       do 1 i=1,3
        do 1 j=1,3
           rmat(j,i) = 0.0d0
    1  continue

       s    = dsqrt( ax(1)**2 + ax(2)**2 + ax(3)**2 )
       if( s.gt.0.d0) then
         ang  = s * scale
         angdeg = ang * 57.29577951d0
         co   = dcos(ang)
         si   = dsin(ang)
         s    = 1.d0 / s
         en(1)= ax(1) * s
         en(2)= ax(2) * s
         en(3)= ax(3) * s
         do 2 i=1,3
          do 2 j=1,3
            rmat(i,j) = en(i)*en(j) * ( 1.d0 - co )
    2    continue
         rmat(1,1) = rmat(1,1) + co
         rmat(2,2) = rmat(2,2) + co
         rmat(3,3) = rmat(3,3) + co

         rmat(1,2) = rmat(1,2) - si * en(3)
         rmat(2,1) = rmat(2,1) + si * en(3)

         rmat(3,1) = rmat(3,1) - si * en(2)
         rmat(1,3) = rmat(1,3) + si * en(2)

         rmat(2,3) = rmat(2,3) - si * en(1)
         rmat(3,2) = rmat(3,2) + si * en(1)
       else
         rmat(1,1) = 1.d0
         rmat(2,2) = 1.d0
         rmat(3,3) = 1.d0
         ang       = 0.d0
         angdeg    = 0.d0
       endif

       return
      END
!
!
!*ds
!*ed
       subroutine matmux( a, b, c )
!      ============================
!-----------------------------------------------------------------------
!      3x3 matrixmultiplikation
!      a * b ---> c
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), b(3,3), c(3,3)
       dimension h(3,3)

       do 10 i=1,3
          do 10 j=1,3
             s = 0.d0
             do 20 l=1,3
                s = s + a(i,l)*b(l,j)
   20        continue
             h(i,j) = s
   10  continue

       do 30 i=1,3
         do 30 j=1,3
            c(j,i) = h(j,i)
   30  continue

       return
      END
!
!
!*ds
!*ed
       subroutine matvec( a, x, y )
!      ============================
!-----------------------------------------------------------------------
!      3x3 matrix - vector - multiplikation
!      a * x ---> y
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), x(3), y(3)
       dimension h(3)

       do 10 i=1,3
             s = 0.d0
             do 20 l=1,3
                s = s + a(i,l)*x(l)
   20        continue
             h(i) = s
   10  continue

       do 30 i=1,3
            y(i) = h(i)
   30  continue

       return
      END
!
!
!*ds
!*ed
       subroutine mattxp( a, at)
!      =========================
!-----------------------------------------------------------------------
!      3x3 matrix - transposition
!      a(t)  ---> at
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), at(3,3)

       at(1,1) = a(1,1)
       at(2,2) = a(2,2)
       at(3,3) = a(3,3)

       xu      = a(1,2)
       xd      = a(2,1)
       at(1,2) = xd
       at(2,1) = xu

       xu      = a(1,3)
       xd      = a(3,1)
       at(1,3) = xd
       at(3,1) = xu

       xu      = a(2,3)
       xd      = a(3,2)
       at(2,3) = xd
       at(3,2) = xu

       return
      END
!
!
!*ds
!*ed
       subroutine matcpy( a, b )
!      =========================
!-----------------------------------------------------------------------
!      3x3 matrix - kopie
!      a     --->   b
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), b(3,3)

       do 10 i=1,3
         do 10 j=1,3
           b(j,i) = a(j,i)
   10  continue

       return
      END
!
!
!*ds
!*ed
       subroutine matone( a )
!      ======================
!-----------------------------------------------------------------------
!      3x3 matrix ---> einheitsmatrix
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3)

       do 10 i=1,3
         do 20 j=1,3
           a(j,i) = 0.d0
   20    continue
         a(i,i) = 1.d0
   10  continue

       return
      END
!
!
!*ds
!*ed
       subroutine axrot( rmat, ax)
!      ===========================
!-----------------------------------------------------------------------
!      rotationsachsendarstellung aus rmat
!      rmat ==~==> ax
!-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension rmat(3,3), ax(3), rtest(3,3)
       dimension d1(3), d2(3), d3(3) , x1(3), x2(3), x3(3)
       dimension tx(3), ty(3)

       d1(1) = rmat(1,1) - 1.d0
       d1(2) = rmat(2,1)
       d1(3) = rmat(3,1)

       d2(1) = rmat(1,2)
       d2(2) = rmat(2,2) - 1.d0
       d2(3) = rmat(3,2)

       d3(1) = rmat(1,3)
       d3(2) = rmat(2,3)
       d3(3) = rmat(3,3) - 1.d0

       x1(1) = d1(2)*d2(3) - d1(3)*d2(2)
       x1(2) = d1(3)*d2(1) - d1(1)*d2(3)
       x1(3) = d1(1)*d2(2) - d1(2)*d2(1)

       s1    =  x1(1)**2 + x1(2)**2 + x1(3)**2
       imax =   1
       smax =   s1

       x2(1) = d3(2)*d1(3) - d3(3)*d1(2)
       x2(2) = d3(3)*d1(1) - d3(1)*d1(3)
       x2(3) = d3(1)*d1(2) - d3(2)*d1(1)

       s2    =  x2(1)**2 + x2(2)**2 + x2(3)**2
       if(s2.gt.smax) then
         imax = 2
         smax = s2
       endif

       x3(1) = d2(2)*d3(3) - d2(3)*d3(2)
       x3(2) = d2(3)*d3(1) - d2(1)*d3(3)
       x3(3) = d2(1)*d3(2) - d2(2)*d3(1)

       s3    =  x3(1)**2 + x3(2)**2 + x3(3)**2
       if(s3.gt.smax) imax = 3

       goto (100,200,300),imax
  100  continue
         s1 = 1.0d0 / dsqrt(s1)
         ax(1) = x1(1) * s1
         ax(2) = x1(2) * s1
         ax(3) = x1(3) * s1
         sn    = 1.d0 / dsqrt( ax(2)**2+ ax(3)**2)
         tx(1) = 0
         tx(2) = ax(3) * sn
         tx(3) =-ax(2) * sn
         call matvec( rmat, tx, ty)
         co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
       goto 999
  200  continue
         s2 = 1.0d0 / dsqrt(s2)
         ax(1) = x2(1) * s2
         ax(2) = x2(2) * s2
         ax(3) = x2(3) * s2
         sn    = 1.d0 / dsqrt( ax(2)**2+ ax(3)**2)
         tx(1) = 0
         tx(2) = ax(3) * sn
         tx(3) =-ax(2) * sn
         call matvec( rmat, tx, ty)
         co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
       goto 999
  300  continue
         s3 = 1.0d0 / dsqrt(s3)
         ax(1) = x3(1) * s3
         ax(2) = x3(2) * s3
         ax(3) = x3(3) * s3
         sn    = 1.d0 / dsqrt( ax(1)**2+ ax(3)**2)
         tx(1) = ax(3) * sn
         tx(2) = 0.0
         tx(3) =-ax(1) * sn
         call matvec( rmat, tx, ty)
         co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
  999  continue

       ang = dacos(co)
       ax(1) = ax(1) * ang
       ax(2) = ax(2) * ang
       ax(3) = ax(3) * ang

! ------------- check sign ! ----------
       call rotax( ax, 1.0d0, rtest, angdeg)
       err =                                                            &
     &       (rtest(1,2)-rmat(1,2))**2+                                 &
     &       (rtest(1,3)-rmat(1,3))**2+                                 &
     &       (rtest(2,3)-rmat(2,3))**2

       if(err.gt.1d-10) then
!        write(6,*)'++++++ rotate sign +++++'
         ang = -ang
         ax(1) =-ax(1)
         ax(2) =-ax(2)
         ax(3) =-ax(3)
       endif

       return
      END
!
!
!*ds
!*ed
       subroutine cross ( v1, v2, vc )
!      ===============================
!
! ** cross product of   v1 x v2  ==> vc
!
       real*8 v1,v2,vc, vch
       dimension v1(3),v2(3),vc(3),vch(3)
!
       vch(1) = v1(2)*v2(3) - v1(3)*v2(2)
       vch(2) = v1(3)*v2(1) - v1(1)*v2(3)
       vch(3) = v1(1)*v2(2) - v1(2)*v2(1)
!
       vc(1)  = vch(1)
       vc(2)  = vch(2)
       vc(3)  = vch(3)

       return
      END
!*ds
!*ed
       subroutine norm( xin, xout, n)
!      ==============================
! --- vektor norm nach skalarprodukt ---
       implicit real*8 (a-h,o-z)
       dimension xin(n), xout(n)

       s = 0.d0
       do 1 i=1,n
         s = s + xin(i)**2
    1  continue
       if(s.gt.0.d0) s = 1.d0/dsqrt(s)

       do 2 i=1,n
         xout(i) = s * xin(i)
    2  continue

       return
      END
!
!
!*ds
!*ed
      function laenge( string, maxlen, delim)
!     =======================================
!     laenge von string, wobei das ende durch delim(1) gekennzeichnet is
!     ------------------------------------------------------------------
      character*1024 string
      character*1  delim

      if(maxlen.gt.1024) then
         n = 1024
      else
         n = maxlen
      endif

      do 1 i=1,n
       if(string(i:i).eq.delim) then
         laenge = i-1
         return
       endif
    1 continue

      laenge = n

      return
      END
!
!
!=======================================================================
!
!=================== io-add-ons for kfa-ibm-use =======================*
!
!*ds
!*ed
      subroutine openkfa(file,ikan,mode,irc)
!     ======================================
!
!    subroutine for cms filedef command
!    file  =    filename filetyp   oder filename.filetyp
!    ikan  =    fortran io-kanal
!    mode  =    'new' oder 'old' oder 'mod'
!    irc   =    return code
!
      character*24  file, fili, cnum*4, line*80, mode*3
      line = ' '
!
      fili = file
      do 10 i=1,24
        if (fili(i:i).eq.'.') fili(i:i)=' '
   10 continue
      ll    = laenge(fili,24,'$')
      fili  = fili(1:ll)
!
      write(cnum,'(i4)') ikan
!
      if(mode.ne.'new') then
         line = 'estate '//fili
         write(6,*)line
         call system(line)
         if (irc .ne. 0) then
            write(6,*)'openkfa file:',file,' not found!'
            return
         endif
      endif
!
      line = 'filedef '//cnum//' '//'disk '//fili
      if (mode.eq.'new') then
          line = 'filedef '//cnum//' '//'disk '//fili                   &
     &           //' ( recfm f lrecl 132 '
      endif
      if (mode.eq.'mod') then
          line = 'filedef '//cnum//' '//'disk '//fili//' ( disp mod'
      endif
!
      write(6,*)'mode=',mode
      write(6,*)line
      call system(line)
      open (ikan)
!
      return
      END
!
!***********************************************************************
!---------- formula-interpreter-section --------------------------------
!***********************************************************************
       subroutine foinit
!      ----------------
!====================================================================
! formelinterpreter kern
! -----------------
!====================================================================
       use formnu
       use formch
       use constants

! ---- internal use ---
       character*4 op
       integer     prio
       real*8      val
       logical     yesno

! ---- inits -----------------------------------------------------------
       ok          = .true.
       error       = .false.
       typ         = '    '
       topnumstack = 0
       topopstack  = 0
       klammerprio = 0
       actchar     = 0
       len         = lclen(formula,maxformlength)
!ccw   write(6,*) 'len=',len
       delims(0)='('
       delims(1)=')'
       delims(2)='+'
       delims(3)='-'
       delims(4)='*'
       delims(5)='/'
       delims(6)='^'
       delims(7)=','
!      do i=0,maxitemlength
!         item(i) = ' '
!      enddo
       return
! ---------------------------------------------------------------------
       entry sayit(yesno)
!      ------------------
       say = yesno
       return
! ---------------------------------------------------------------------
       entry setdeg
!      ------------
       degree = datan(1.d0)/45.d0
       write(6,*)'evaluator set to degree'
       return
! ---------------------------------------------------------------------
       entry setrad
!      ------------
       degree = 1.d0
       write(6,*)'evaluator set to rad'
       return
! ---------------------------------------------------------------------
       entry show
!      ----------
       if(.not.say) return
       write(6,*) (formula(i),i=0,len)
       write(6,*) 'actchar=',actchar, '  klammerprio=',klammerprio
       write(6,*) 'item   =',item
       do i=1,topopstack
         write(6,*)i,'  ',opstack(i),'   ',priostack(i)
       enddo
       write(6,*)
       do i=1,topnumstack
         write(6,*)i,'  ',numstack(i)
       enddo
       return
! ---------------------------------------------------------------------
       entry putopstack( op, prio)
!      -----------------------------

!ccw   write(6,*)'putopstack:',op,' ',prio
       if(topopstack.le.maxopstack) then
          topopstack = topopstack + 1
          opstack(topopstack)   = op
          priostack(topopstack) = prio+klammerprio
!ccw   write(6,*)'topopstack=',topopstack
       else
          if(say) write(6,*)'error: opstack overflow!'
          error = .true.
       endif

       return
! ---------------------------------------------------------------------
       entry enternumstack( val )
!      --------------------------

       if(topnumstack.le.maxnumstack) then
          topnumstack = topnumstack + 1
          numstack(topnumstack) = val
       else
          if(say) write(6,*)'error: numstack overflow!'
          error = .true.
       endif

       return
! ---------------------------------------------------------------------
       entry checknum( n )
         if(topnumstack.lt.n) then
            if(say) write(6,*)'error: too few num operands!'
            error = .true.
            do i=topnumstack,n
              numstack(i) = 0.d0
            enddo
         endif
       return
! ---------------------------------------------------------------------
       entry getword
       litem =-1
       do i = actchar,len
         do j=0,nodelims
          if(formula(i).eq.delims(j)) goto 100
         enddo
         litem = litem+1
         if(litem.gt.maxitemlength) goto 999
         item(litem) = formula(i)
       enddo
       actchar = len+1
  100  continue
       actchar = i
       if(litem.gt.maxitemlength) goto 999
       item(litem+1) = ' '
!ccw   write(6,*)'getitem: item=',(item(j),j=0,maxitemlength)
       return

  999  continue
       if(say) write(6,*)'getitem: item too long'
       error = .true.
       return

      END


!*ds
!*ed
       subroutine getitem
!      ------------------


       use formnu
       use formch
       use constants

! --- internal use ---
       implicit none
       integer j, intn, ier, ierr, l, ll, lit
       logical compare
       logical anklam

       character*(maxitemlength+1) citem

       !needed since common-block conversion - maybe this could be simpliefied
       !equivalence is the problem here
       character*1 tempitem(0:maxitemlength)

       equivalence(citem,tempitem(0))

       !the assignment hast to follow equivalence or else the compiler will explode!
       tempitem=item
       if(actchar.gt.len) then
          typ = 'end '
          call putopstack('end ',0)
          return
       endif

       if(formula(actchar).eq.'(') then
          typ = 'klam'
          klammerprio = klammerprio+klammerinc
          actchar     = actchar + 1
          goto 100
       endif
       if(formula(actchar).eq.')') then
          typ = 'klam'
          klammerprio = klammerprio-klammerinc
          actchar     = actchar + 1
          goto 100
       endif
       if(formula(actchar).eq.',') then
          typ = 'list'
          call putopstack('nop ',komprio  )
          actchar     = actchar + 1
          goto 100
       endif
       if(formula(actchar).eq.'+') then
          anklam = (actchar.eq.0)
          if(.not.anklam) anklam = (formula(actchar-1).eq.'(')
          if(anklam) then
            actchar     = actchar + 1
          else
            typ = 'bin '
            call putopstack('+   ',iplusprio)
            actchar     = actchar + 1
          endif
          goto 100
       endif
       if(formula(actchar).eq.'-') then
          anklam = (actchar.eq.0)
          if(.not.anklam) anklam = (formula(actchar-1).eq.'(')
          if(anklam) then
            call putopstack('neg ',iuprio)
            actchar     = actchar + 1
          else
            typ = 'bin '
            call putopstack('-   ',minusprio)
            actchar     = actchar + 1
          endif
          goto 100
       endif
       if(formula(actchar).eq.'*') then
            typ = 'bin '
            call putopstack('*   ',multprio)
            actchar     = actchar + 1
          goto 100
       endif
       if(formula(actchar).eq.'/') then
            typ = 'bin '
            call putopstack('/   ',idivprio)
            actchar     = actchar + 1
          goto 100
       endif
       if(formula(actchar).eq.'^') then
            typ = 'bin '
            call putopstack('^   ',iexpprio)
            actchar     = actchar + 1
          goto 100
       endif
! ---  suche bis zum naechsten delimiter ---
       citem = ' '
       item=tempitem
	call getword
! --- ist item ein unaerer operator ? ---
       if( compare(item,'sin ') ) then
         typ = 'unae'
         call putopstack('sin ',iuprio)
         goto 100
       endif
       if( compare(item,'cos ') ) then
         typ = 'unae'
         call putopstack('cos ',iuprio)
         goto 100
       endif
       if( compare(item,'tan ') ) then
         typ = 'unae'
         call putopstack('tan ',iuprio)
         goto 100
       endif
       if( compare(item,'atan ') ) then
         typ = 'unae'
         call putopstack('atan',iuprio)
         goto 100
       endif
       if( compare(item,'asin ') ) then
         typ = 'unae'
         call putopstack('asin',iuprio)
         goto 100
       endif
       if( compare(item,'acos ') ) then
         typ = 'unae'
         call putopstack('acos',iuprio)
         goto 100
       endif
       if( compare(item,'ln ') ) then
         typ = 'unae'
         call putopstack('ln  ',iuprio)
         goto 100
       endif
       if( compare(item,'exp ') ) then
         typ = 'unae'
         call putopstack('exp ',iuprio)
         goto 100
       endif
       if( compare(item,'sqrt ') ) then
         typ = 'unae'
         call putopstack('sqrt',iuprio)
         goto 100
       endif
       if( compare(item,'abs ') ) then
         typ = 'unae'
         call putopstack('abs ',iuprio)
         goto 100
       endif
       if( compare(item,'int ') ) then
         typ = 'unae'
         call putopstack('int ',iuprio)
         goto 100
       endif
! --- ist item eine userfunction ? -----
       if( formula(actchar).eq.'(') then
         typ = 'usrf'
         call putopstack('usrf',iuprio)
         if(tusrfstack.lt.musrfstack) then
            tusrfstack=tusrfstack+1
         else
            error = .true.
            write(6,*)'userfuncstack exceeded'
         endif
         usrfstack(tusrfstack) = citem(1:19)//' '
         goto 100
       endif
! --- ist item eine zahl ? ----
       call scan(item ,valnum,ierr)
! --- try exp-num-representation --
       if(ierr.ne.0) then
!ccw   write(6,*)'check exp-num item(litem)=',item(litem),'  ',litem
!ccw   write(6,*)'formula(actchar)=',formula(actchar)
          lit = litem
          if((item(litem).eq.'e').or.(item(litem).eq.'d') .or.(item(litem).eq.'E').or.(item(litem).eq.'D')) then
            if((formula(actchar).eq.'+').or.(formula(actchar).eq.'-')) then
               lit = lit+1
               item(lit)=formula(actchar)
               do l=actchar+1,len
                  do ll=0,nodelims
                    if(formula(l).eq.delims(ll)) goto 400
                  enddo
                  lit  =lit  +1
                  item(lit  ) = formula(l)
               enddo
               l = len+1
  400          continue
               call scan(item ,valnum,ierr)
               if(ierr.eq.0) then
                 actchar = l
               else
!cc              do l=litem+1,maxitemlength
!cc                item(l)=' '
!cc              enddo
                 item(litem+1)=' '
               endif
            endif
          endif
       endif
       if (ierr.eq.0) then
         typ = 'num '
         call enternumstack(valnum)
         goto 100
       endif
! --- ist item eine referenz ? ----
       call extract(item,valnum,ierr)
       if (ierr.eq.0) then
         typ = 'num '
         call enternumstack(valnum)
         goto 100
       endif
! --- fehler:
       error = .true.
       if(say) write(6,*)'item:',item,' not decodable'

  100  continue

       return

! ----------------------------------------------------------------------
       entry stackevaluate
!      -------------------
       j = topopstack-1
       if(j.lt.0) return

    1  continue
       if((j.gt.0).and.(priostack(j).ge.priostack(topopstack))) then

         if(opstack(j).eq.'nop ') then
            goto 200
         endif
         if(opstack(j).eq.'usrf') then
            call checknum(1)
            call usrfun(usrfstack(tusrfstack),numstack,topnumstack,ier)
            if(ier.ne.0) error = .true.
            call checknum(1)
            tusrfstack=tusrfstack-1
            goto 200
         endif
         if(opstack(j).eq.'+   ') then
            call checknum(2)
            numstack(topnumstack-1)=                                    &
     &      numstack(topnumstack-1)+numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'-   ') then
            call checknum(2)
            numstack(topnumstack-1)=                                    &
     &      numstack(topnumstack-1)-numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'*   ') then
            call checknum(2)
            numstack(topnumstack-1)=                                    &
     &      numstack(topnumstack-1)*numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'/   ') then
            call checknum(2)
            numstack(topnumstack-1)=                                    &
     &      numstack(topnumstack-1)/numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'^   ') then
            call checknum(2)
            numstack(topnumstack-1)=                                    &
     &      numstack(topnumstack-1)**numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif

         if(opstack(j).eq.'neg ') then
            call checknum(1)
            numstack(topnumstack)=    -numstack(topnumstack)
            goto 200
         endif
         if(opstack(j).eq.'sin ') then
            call checknum(1)
            numstack(topnumstack)=dsin(degree*numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'cos ') then
            call checknum(1)
            numstack(topnumstack)=dcos(degree*numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'tan ') then
            call checknum(1)
            numstack(topnumstack)=dtan(degree*numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'asin') then
            call checknum(1)
            numstack(topnumstack)=dasin(numstack(topnumstack))          &
     &                            /degree
            goto 200
         endif
         if(opstack(j).eq.'acos') then
            call checknum(1)
            numstack(topnumstack)=dacos(numstack(topnumstack))          &
     &                            /degree
            goto 200
         endif
         if(opstack(j).eq.'atan') then
            call checknum(1)
            numstack(topnumstack)=datan(numstack(topnumstack))          &
     &                            /degree
            goto 200
         endif
         if(opstack(j).eq.'ln  ') then
            call checknum(1)
            numstack(topnumstack)=dlog(numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'exp ') then
            call checknum(1)
            numstack(topnumstack)=dexp(numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'sqrt') then
            call checknum(1)
            numstack(topnumstack)=dsqrt(numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'abs ') then
            call checknum(1)
            numstack(topnumstack)=dabs(numstack(topnumstack))
            goto 200
         endif
         if(opstack(j).eq.'int ') then
            call checknum(1)
            intn = numstack(topnumstack)
            numstack(topnumstack)=intn
            goto 200
         endif
  200  continue
       j = j-1
       else
        j = j+1
        opstack(j) = opstack(topopstack)
        priostack(j)=priostack(topopstack)
        topopstack = j
        return
       endif
       goto 1

      END

!*ds
!*ed
       subroutine evaluate( f, val, ierr)
!      ----------------------------------

       use xoutxx
       use formnu
       use formch
       use outlev !for iout
       use constants
! --- internal use ---

       implicit none
       real*8 val
       integer i, ierr
       character*(maxformlength+1) f
!
       do i=0,maxformlength
         formula(i) = ' '
       enddo
       do i=1,maxformlength
        if(f(i:i).eq.' ') goto 1
        formula(i-1)= f(i:i)
       enddo
    1  continue

       call foinit
       do i=1,maxformlength
         call getitem
         call show
         if(typ.ne.'num ') then
            call stackevaluate
            call show
         endif
         if((typ.eq.'end ').or.error) goto 10
       enddo
   10  continue
       if(error) then
         ierr = 1000
       else
         ierr = 0
       endif
       if(klammerprio.ne.0) then
         if(say)write(6,*)'unmatched parentheses'
         ierr = ierr+1
       endif
       if(opstack(1).ne.'end ') then
         if(say)write(6,*)'opstack items left!'
         ierr = ierr+10
       endif
       if(topnumstack.ne.1) then
         if(say)write(6,*)'numstack not ok'
         ierr = ierr+100
       endif

       val = numstack(1)

       call scan(formula,valnum,ierrs)
       if(ierrs.ne.0 .and. iout.ge.1) write(6,*)'evaluate: ',(formula(i),i=0,len),' to ',val
			 return
      END



       subroutine cappend( a, b, c)
!      ----------------------------
       implicit none
       character*16 a,b,c

       integer i,j,l
       c = ' '
       do i=1,16
         if(a(i:i).ne.' ') then
          c(i:i) = a(i:i)
         else
          goto 1
         endif
       enddo
    1  continue
       j = i
       do i=j,16
         l=i-j+1
         if(b(l:l).ne.' ') then
          c(i:i) = b(l:l)
         else
          goto 2
         endif
       enddo
    2  continue
       j = i
       do i=j,16
        c(j:j) = ' '
       enddo
      END

       integer function lclen( string, m)
!      ----------------------------------
      implicit none
      integer m 
      character*1 string(0:m)
       do lclen=m,0,-1
         if(string(lclen).ne.' ') return
       enddo
      END

       logical function compare( s1, s2 )
!      ----------------------------------
       implicit none
       character*1 s1(*),s2(*)
       
       integer i
       compare = .true.
       do i=1,255
         if(s1(i).ne.s2(i)) then
           compare = .false.
           return
         endif
         if((s1(i).eq.' ').or.(s2(i).eq.' ')) return
       enddo
       compare = .false.
       return
      END


!*ds
!*ed
      subroutine extract(nam,val,ier)
!     -------------------------------
      use usevar
      use imargs

      implicit none
      integer i, ier
      character*1 nam(*)
      real*8 val
      logical compare

      ier = 0
! --- look in uservars ---
      do i=1,nousev
        if(compare(nam,usenam(i))) then
          val = useval(i)
          return
        endif
      enddo
! --- number of doargs ---
      if(compare(nam,'noargs ')) then
        val = ipmlst(ktop)
        return
      endif
! --- user defined vals ---
      call usrextr(nam,val,ier)
      return
      END

!*ds
!*ed
      subroutine setudf(nam,val,ier)
!     ------------------------------
      use xoutxx
      use usevar
      use constants
      implicit none
      character*16 nam
      real*8 val
      integer ier

      logical compare
      integer i,j

      ier = 0
! --- look in uservars ---
!ccw  write(6,*)'setudf: nam=',nam,' val=',val
      do i=1,nousev
        if(compare(nam,usenam(i))) then
          useval(i) = val
          if(iot.ge.1)                                                  &
     &    write(6,*)'uservar:',usenam(i),' is set to:',useval(i)
          return
        endif
      enddo
!ccw  write(6,*)'create new udf'
! --- create new var ---
      if(nousev.lt.musevar) then
        nousev = nousev+1
        do i=1,16
          usenam(nousev)(i:i)= nam(i:i)
          if(nam(i:i).eq.' ') then
           do j=i,16
            usenam(nousev)(j:j) = ' '
           enddo
           goto 1
          endif
        enddo
    1   continue
        useval(nousev) = val
        i = nousev
        if(iot.ge.1)                                                    &
     &  write(6,*)'uservar:',usenam(i),' is created:',useval(i)
        return
      else
        write(6,*)'uservar set: buffer is full ! error'
        ier = 200
      endif
      return
      END

!*ds
!*ed
      subroutine clrudf(nam)
!     ----------------------
      use usevar
      implicit none
      integer i,j
      character*16 nam

      logical compare

      if(nam(1:3).eq.'all') then
        nousev = 0
        write(6,*)'uservar stack cleared'
      else
      do i=1,nousev
        if(compare(nam,usenam(i))) then
          write(6,*)'uservar:',usenam(i),'=',useval(i),' removed!'
          do j=i+1,nousev
            usenam(j-1)=usenam(j)
            useval(j-1)=useval(j)
          enddo
          nousev = nousev-1
          return
        endif
      enddo
      endif
      return
      END

!*ds
!*ed
      subroutine shwudf
!     -----------------
      use usevar
      implicit none
      integer i

      if(nousev.eq.0) then
        write(6,*)'no uservars defined!'
      else
        write(6,*)'defined uservars:'
        do i=1,nousev
           write(6,*)usenam(i),useval(i)
        enddo
      endif
      return
      END

      function inpaf(iaddr)
!-----------------------------------------------------------------------
!  wert extraktion aus dem incom parameterstack  zuordnung inpar
!  still needed for one routine in datreat_subroutines
!-----------------------------------------------------------------------
       use cincom
       implicit none
       integer iaddr, inpaf
       inpaf = inpar(iaddr)
      return
      END
