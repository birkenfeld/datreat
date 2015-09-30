c***********************************************************************
c       incom-system                                                   *
c       ------------                                                   *
c       by  michael monkenbusch                                        *
c       institut fuer festkoerperforschung des                         *
c       forschungszentrum juelich, d-5170 juelich, po-box 1913, frg    *
c       iff115 @ djukfa11                                              *
c                                                                      *
c***********************************************************************
c*ds
c*ed
       subroutine incom( cmd )
c      -----------------------
c
c ---- command line decoder ----
c
c
c ---- separates commands parameternames & parameters ----
c
c ---- commandline format may be :
c      command <name parameter name parameter < ;command ....>>
c
c ---- parameter separator is blank !!!
c ---- command separator is ;
c
c ---- the current command is stored on comand
c      the names on stack vname(*)
c      the numbers on stack rpar(*)
c      in common/cincom/
c
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
 
c ---- communication common block conainig the analysed inputline
c      comand   = actual commoand keyword
c      vname(*) = names stack
c      rpar(*)  = number stack
c      inames   = length of names stack
c      ipars    = length of number stack
c      ioldc    = flag to indicate second command in line
c      inline   = actual inputline
c      reslin   = residual line conaining all still nonprocessed command
c      inpar(*) = no. of parameters associated with the i-th name
c      iparn(*) = pointer parameters to names
c      inapa(*) = adress of first parameter following the *-th name
c      arglst(*)= character*20 list of arguments
c      iargs    = length of arglst = no. of parameters
c      pmlist(*,i) => dummy, left only for compatibility reasons
c      ipmls    =       "                    "
c      iolbuf   = buffer for saving ioldc during a makro call
c      rlbuf    = buffer to save reslin during a makro call
c      lstpar   = number of the last decoded parameter
c      lstnam   = number of the last decoded name
 
      parameter(mdepth=20)
      character*80 argvals
      character*80 pmlst
      common/cmargs/argvals(minc),pmlst(mdepth,minc,2)
      common/imargs/iargvs,ipmlst(mdepth),kanal(0:mdepth),ktop
c --- variables for makro-parameter-passing ---
c     argvals(i)  = parameterlist at a call of a makro, only temp.
c     pmlst(k,i,1..2) = replace list, in the k-th makro-layer
c                       all substrings=pmlst(k,*,1) will be
c                       replaced by pmlst(k,*,2)
c     iargvs      = number of arguments
c     ipmlst(k)   = number of given replacing items for makro k
c     kanal(k)    = fortran io-number associated to makro layer k
c     ktop        = current makro layer (0=keyboard-input)
c
c
c ---- outputlevel
      logical cray
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
       common/xroxxx/  xyorig(3), rotvec(3)
c
c
       character*8   glabel,xlab,itypc,vnamef
       character*132 buf,bla132
       character*1   blank,csep,numtst(14),numeqv*14,inread*80
       character*8   cmd
       equivalence  (numeqv,numtst(1))
       logical       name, lmakro, found
 
       data          lmakro /.false./
c
c ---- initialisations & conventions ----
c
       csep = ';'
c              = --> command separator
       blank= ' '
c              = --> parameter separator
       numeqv = '(.+-0123456789'
c                --------------> to reckognize the beginning of a number
       bla132 = ' '
c
 
8888   continue
c      --------> reentry-point
       do 5 i=1,minc
        vname(i) = ' '
        rpar(i)  = 0
        iparn(i) = 0
        inpar(i) = 0
        inapa(i) = 0
5      continue
       comand = ' '
 
       lstpar    = 0
       lstnam    = 0
c
c ---- error-response -----------------------------------------
       if(ierrr.ne.0) then
         write(6,*)' error return code:',ierrr
c
cray -------------------------------------------
         if(cray) stop
cray -------------------------------------------
c
         ioldc = 0
c        ----------> if error forget rest of commandline (reslin)
         if(ktop.ne.0) then
           do k=ktop,1,-1
            close(kanal(k))
           enddo
           ktop = 0
c          --------> go back to input from keyboard
         endif
       endif
       ierrr = 0
c
c
c ---- input a line or take nodecode rest of previous line ----
1000   continue
       if(ioldc.eq.0) then
         goto 1002
1001     continue
         if(iot.gt.-3) write(6,*)'end of data '
         close(kanal(ktop))
         if(ktop.eq.0) then
            open(kanal(ktop))
         else
           ktop = ktop-1
c                      --------> restore old reslin state
           ioldc = iolbuf
           reslin= rlbuf
           ipmls = 0
           if(ktop.eq.0)write(6,*)' command input back to keybord'
cccweg     inka = 5
           if(ioldc.ne.0) then
             inline = reslin
             goto 1011
           endif
         endif
1002     continue
 
         if(ktop.eq.0) then
            write(6,*)' #-->'
            ipmls = 0
         endif
         read(kanal(ktop),1,end=1001) inread
         call lowcase(inread,80)
         if(ktop.ne.0.and.iot.gt.-2) write(6,*)'! ',inread(1:78)
         if(ktop.eq.0)               write(6,*)'> ',inread(1:78)
         inline = inread
         inlbuf = inline
1        format(a)
       else
         inline = reslin
       endif
1011   continue
       ioldc = 0
       ioldna =0
c
c ---- look for a commnd separator ----
       do 10 i=1,132
        if(inline(i:i).eq.csep) then
          ioldc = 1
          j = i+1
          reslin = inline(j:132)
c         ---------------------> store second command on reslin
          j = i-1
          inline = inline(1:j)
          goto 11
        endif
10     continue
11     continue
c
c
c ---- remove any leading blanks ----
       do 20 i=1,132
        if(inline(i:i).ne.blank) then
          inline = inline(i:132)
          goto 21
        endif
20     continue
21     continue
c
c ---- separate the command ----
       if(inline(1:1).eq.blank) then
         comand = ' '
       else
         do 30 i=1,132
          if(inline(i:i).eq.blank) then
            j = i-1
            comand = inline(1:j)
            goto 31
          endif
30       continue
31       continue
       endif 
 
       cmd = comand

c mz
c replace '.' with 'dir' so that the programm doesn't crash
       if(comand.eq.'.       ') then
          ioldc = 1
          reslin = 'dir    '
          goto 8888
       endif
c
 
c ---- dont analyse further if it is a title ----
         if(comand.eq.'tit     '.or.comand.eq.'title   ') then
           if(comand.eq.'tit     ') title = inline(5:132)
           if(comand.eq.'title   ') title = inline(7:132)
           goto 8888
         endif
c ---- dont analyse further if it is a message ---
         if(comand.eq.'msg     ') then
           write(6,*)inline
           goto 8888
         endif
c ---- dont analyse further if it is a blank line ---
         if(comand.eq.'        ') then
           write(6,*)inline
           goto 8888
         endif
c ---- dont analyse further if it is a comment ---
         if(comand.eq.'c       ') then
           goto 8888
         endif
c ---- dont analyse further if it is a label ----
         if(comand(1:1).eq.':') then
           if(iot.gt.2)write(6,*)'label encountered ',comand
           goto 8888
         endif
c
c ---- look for keywords & parameters ----
c
       j = j+1
       inames = 0
       ipars  = 0
       iargs  = 0
       iargvs = 0
c
400    continue
       i = j
       l = i+1
       if(inline(i:i).eq.blank.and.inline(l:l).ne.blank) then
c ---- look for the end of the currend item -----
         do 60 k=l,132
          if(inline(k:k).eq.blank) then
            j = k-1
            buf = inline(l:j)
            goto 61
          endif
60       continue
61       continue
 
c --- prepare the textlist of formal makro arguments --
         iargs = iargs + 1
         arglst(iargs) = buf(1:20)
 
c --- decode for makro argument values replace strings ! ----
         if(ktop.ne.0) then
         if(ipmlst(ktop).ne.0) then
          if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',
     *                 ipmlst(ktop)
           do k=1,ipmlst(ktop)
        if(iot.gt.3)write(6,*)buf
            call creplace(buf,pmlst(ktop,k,1)//' ',
     *                           pmlst(ktop,k,2)//' ',' ')
        if(iot.gt.3)write(6,*)buf
           enddo
          if(buf.eq.bla132) goto 410
         endif
         endif
c
c ----   discriminate between name & parameter ----
         name = .true.
         do 50 k=1,14
           if(buf(1:1).eq.numtst(k)) then
             name = .false.
             goto 51
           endif
50       continue
51       continue
         if(.not.name) then
           call evaluate( buf//' ', val, ierr)
           if(ierr.eq.0) then
              name = .false.
           else
              comand = 'f-error!'
              cmd    = comand
              write(6,*)'syntax error in:',buf
           endif
         endif
c
         iargvs = iargvs+1
         if(name) then
          inames = inames + 1
          if(inames.gt.minc) goto 999
          vname(inames) = buf(1:8)
          argvals(iargvs) = buf(1:80)
          inapa(inames) = 0
         else
          ipars = ipars + 1
          if(ipars.gt.minc) goto 999
cccc      read(buf,*)rpar(ipars)
cccc      call scan(buf,rpar(ipars),ierr)
          rpar(ipars) = val
          iival = val+1d-11
          if(dabs(val-iival).lt.3d-11) then
            write(argvals(iargvs),'(i80)') iival
          else
            write(argvals(iargvs),'(e80.12)') rpar(ipars)
          endif
          iparn(ipars) = inames
          if(inames.ne.ioldna) then
            inapa(inames) = ipars
            ioldna = inames
          endif
         endif
       endif
410    continue
       j = j+1
       if(j.lt.132-1) goto 400
c
c ---- prepare inpar (no. of parameters / name) ----
       do 700 i=1,inames
         isum = 0
         if(ipars.ne.0) then
           do 702 j=1,ipars
             if(iparn(j).eq.i) isum = isum + 1
702        continue
         endif
         inpar(i) = isum
700    continue
c
       inline = inlbuf
c
       if(lmakro) goto 9888
 
c --------------- decoding of opaque commands -------------------------
c
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c cms - kommando absetzen                             | cms <cmd >
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'cms     ') then
c                    ---
         inline = inline(5:80)
ccccc
         do i=1,75
          if(inline(i:i).eq.' ') inline(i:i) = '~'
         enddo
         inline(75:75) = ' '
         if(ktop.ne.0) then
         if(ipmlst(ktop).ne.0) then
          if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',
     *                 ipmlst(ktop)
           do k=1,ipmlst(ktop)
            call creplace(inline,pmlst(ktop,k,1)//' ',
     *                           pmlst(ktop,k,2)//' ',' ')
           enddo
         endif
         endif
         do i=1,75
          if(inline(i:i).eq.'~') inline(i:i) = ' '
         enddo
         inline(75:75) = ' '
ccccc
         call system(inline)
         if(irc.ne.0) then
           call errsig(1,'wrong cms command ...$')
         endif
         goto 8888
       endif
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c output"level" setzen:                               | iout <iout>
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'iout    ') then
c                    ----
         ioold= iot
         iot  = rpar(1)*1.0000000001d0
         if(vname(1).eq.'old     ') iot = ioold
         write(6,*)'outputlevel is now ',iot
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c cray-flag setzen: (toggle)                          | cray
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'cray    ') then
c                    ----
         cray = .not. cray
         if(cray) then
           write(6,*)'cray execution mode'
         else
           write(6,*)'cms  execution mode'
         endif
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c if                                                  |
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'if      ') then
c                    ----------------> if-bedingung (numerisch)   !
         itypc = vnamef(1)
         if(itypc.eq.'=       ') then
           if(rpar(1).eq.rpar(2)) then
             j = index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'>       ') then
           if(rpar(1).gt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<       ') then
           if(rpar(1).lt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<=      ') then
           if(rpar(1).le.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'>=      ') then
           if(rpar(1).ge.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<>      ') then
           if(rpar(1).ne.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:132)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         write(6,*)'if operator ',itypc, ' unknown'
         ierrs = 1
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c realisation goto                                    | goto
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'goto    ') then
c                    ----------------> comment feature for makros !
         if(ktop.gt.0) then
           glabel = vnamef(1)
           if(glabel(1:1).ne.':') then
              write(6,*)glabel,' is not a valid label'
              write(6,*)'use >>> :',glabel,' <<< !'
              ierrs = 1
              goto 8888
           endif
           rewind(kanal(ktop))
88881      continue
            read(kanal(ktop),'(a8)',end=88883)xlab
            if(glabel.eq.xlab) goto 8888
            goto 88881
         else
           write(6,*)'goto in keyboardmode ignored'
         endif
88883    continue
         close(kanal(ktop))
         ktop = ktop-1
         write(6,*)'label:',glabel,' not found!'
         ierrs = 1
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c return from makro                                   | return
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'return  ') then
c                    ----------------> return from makro
         if(ktop.gt.0) then
           close(kanal(ktop))
           ktop = ktop - 1
         endif
         goto 8888
       endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c exit all makros                                     | exit
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'exit    ') then
c                    ----------------> exit   from makro
         if(ktop.gt.0) then
           do k=ktop,1,-1
             close(kanal(ktop))
           enddo
           ktop = 0
         endif
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c evaluation info                                     | ??
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'??      ') then
c                    ----------------> evaluation info !
         do i=1,ipars
          write(6,*) rpar(i)
         enddo
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c evaluation and  print                               | print
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'print   ') then
c                    ----------------> print evaluated numbers !
    
         if(inames.lt.1) then
           write(6,*)'print: needs filename'
           goto 8888
         endif 
         open(99,file=vname(1),status='unknown',position='append')  
         write(6,'(10(1x,e13.6))') (rpar(i),i=1,ipars)
         write(99,'(10(1x,e13.6))') (rpar(i),i=1,ipars)
         close(99)
         goto 8888
       endif
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'f-error!') then
c                    ----------------> fehler-reaktion !
         write(6,*)'!! commandline in error --> ignored !!'
         ierrr = 111
         goto 8888
       endif
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c (user)vars zeigen                                   | vars?
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'vars?   ') then
c                    ----------------> uservars display
         call shwudf
         goto 8888
       endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c formelauswertung protokoll                          | say
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'ftrace  ') then
c                    ----------> set deg-mode
         call sayit(.not.found('off     '))
         goto 8888
       endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c degmode setzen                                      | setdeg
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'setdeg  ') then
c                    ----------------> set deg-mode
         call setdeg
         goto 8888
       endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c radmode  setzen                                     | setrad
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'setrad  ') then
c                    ----------------> set deg-mode
         call setrad
         goto 8888
       endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c uservariable setzen                                 | set
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'set     ') then
c                    ----------------> set uservar
         do i=1,inames
           call setudf(vname(i)//' ',rpar(inapa(i)),ier)
         enddo
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c uservariable entefrnen                              | clr
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'clr     ') then
c                    ----------------> set uservar
         do i=1,inames
           call clrudf(vname(i)//' ')
         enddo
         goto 8888
       endif
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c rotationsparameter setzen                           | setrot <xyz> <r>
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'setrot  ') then
c                    ----------------> set rotation-parameters !
           call get6( rotvec(1), rotvec(2), rotvec(3),
     *                xyorig(1), xyorig(2), xyorig(3))
           xyorig(1) = getval('x       ',xyorig(1),inew)
           xyorig(2) = getval('y       ',xyorig(2),inew)
           xyorig(3) = getval('z       ',xyorig(3),inew)
           rotvec(1) = getval('rx      ',rotvec(1),inew)
           rotvec(2) = getval('ry      ',rotvec(2),inew)
           rotvec(3) = getval('rz      ',rotvec(3),inew)
           write(6,*)'------------------------------------------------'
           write(6,'(a,3f12.6/a,3f12.6)')' rotation origin =',xyorig,
     *                                   ' rotation vector =',rotvec
           write(6,*)'------------------------------------------------'
         goto 8888
       endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c ausgang                                             | q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
       if(comand.eq.'q       ') then
c                    ----------------> exit the program
         if(ibild.gt.0) then
           write(6,*)'calling grende ...'
           call grende
         endif
         write(6,*)'exit by q, bye...'
         stop
       endif
c
       return
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c ---- stack overfow error trap ----
999    continue
       write(6,*)' incom(',comand,') :  stack overflow   reenter'
       ioldc = 0
       ierrr = 2
       goto 8888
c
c
c
c --- makro - aktivator - entry ----------------------------------------
c
       entry makro
c      ===========
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-
c --- try to interpret non-known commands as makro input files ---
c ------> any makro input file must have the keyword  <makro>
c         at the first position
       if(ktop.ge.mdepth) then
         write(6,*)' too many stacked makros : ',comand
         ierrr = 3
         return
       endif
       ktop = ktop + 1

!!! Modification for two macro sources: global, individual !!!
!!! try first the individual makro as we had it already    !!!

       open(kanal(ktop),
     *          file=comand,status='old',form='formatted',err=99991)
       read(kanal(ktop),9998,end=99991) inline
9998   format(a80)
       call lowcase( inline , 80)
       if(inline(1:6).eq.'makro ') goto 19999
       
!!! ok this was no makro !!!
       close(kanal(ktop))
99991  continue
c -------> try whether in the global makro section file is a makro file
       open(kanal(ktop),
     *          file='/usr/users/iff_ns/monk/datreat/makros/'//comand,
     *          status='old',form='formatted',err=9999)
       read(kanal(ktop),9998,end=9999) inline
!9998   format(a80)
       call lowcase( inline , 80)
       if(inline(1:5).ne.'makro') goto 9999
c -------> try whether file is a makro file
c
19999  continue
c --- ok this is a makro file ! save first the parametervalues ---
         do i=1,minc
           pmlst(ktop,i,2)=' '
         enddo
         do 9005 i=1,iargvs
           k = 0
           do l=1,80
             if(argvals(i)(l:l).ne.' ') then
               k = k + 1
               pmlst(ktop,i,2)(k:k) = argvals(i)(l:l)
             endif
           enddo
           if(iot.gt.3) write(6,*)'argvals : ',pmlst(ktop,i,2)
9005     continue
         ipmlst(ktop) = iargvs
         do 9007 i=iargvs+1,minc
           pmlst(ktop,i,2) = ' '
9007     continue
c --- save reslin form the current keyboard command ---
         iolbuf = ioldc
         rlbuf  = reslin
c --- decode the first makro line to get the parameternames
         ioldc = 1
         reslin = inline
         lmakro = .true.
         if(iot.gt.-3)write(6,*)'this is a makro ...'
         goto 8888
9888     continue
         lmakro = .false.
         if(iot.gt.-3)write(6,*)'go on with it ...'
         if(iargs.ne.0) then
           do 9010 i=1,iargs
            if(iot.gt.3) write(6,*)'arglst:',arglst(i)
            pmlst(ktop,i,1) = arglst(i)
9010       continue
         endif
         ipmlst(ktop) = iargs
c        ipmls = iargs
c
       return
c      ---------> read now commandlines from makro file
9999   continue
       do k=ktop,1,-1
         close(kanal(k))
       enddo
       ktop = 0
       write(6,777)comand
  777  format(' program name ''',a8,''' not known')
       ierrr = 4
       return
 
       end
c
       block data kanass
c      -----------------
       parameter(minc=40)
       parameter(mdepth=20)
       character*80 argvals
       character*80 pmlst
       common/cmargs/argvals(minc),pmlst(mdepth,minc,2)
       common/imargs/iargvs,ipmlst(mdepth),kanal(0:mdepth),ktop
       data kanal/ 5,40,41,42,43,44,45,46,47,48,
     *            49,50,51,52,53,54,55,56,57,58, 59/
       data ktop /0/
       end
 
c*ds
c*ed
ccc+++ noch optimierungsbeduerftig !!! ++++++++++++++++++<<<!!!!!!
       subroutine creplace(a,b,c,s)
c---------------------------------------------------
c ersetzten von string b durch c in string a
c strings werden durch den character s begrenzt !
c----------------------------------------------------
       parameter(mclen=132)
       character*(mclen) a,b,c,d
       character*1 s
c
c ---  look for occurence of substring b in a ---
       i = 1
       d = a
       len = 0
       do l=1,mclen
         if(a(l:l).eq.s) then
           len = l-1
           goto 1
         endif
       enddo
1      continue
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
10        continue
          if(istart.ne.0) then
            iend = l
            do jj=1,mclen
              ld= istart-1+jj
              if(c(jj:jj).eq.s  ) goto 20
              d(ld:ld)= c(jj:jj)
            enddo
            i = ld
20          continue
            do jj=iend+1,mclen
             if(l.gt.mclen) goto 30
             d(ld:ld)= a(jj:jj)
             if(a(jj:jj).eq.s  ) goto 30
             ld = ld+1
            enddo
30          continue
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
       end
c
c
c*ds
c*ed
       subroutine oldlowcase( string , l)
c      ==============================
c ---> reverse action for unix system <----
c
       character*132 string
       character*26  upper, lower , c*1
       DATA LOWER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
       data upper/'abcdefghijklmnopqrstuvwxyz'/
c
       do 1 i=1,l
        c = string(i:i)
        do 2 j=1,26
          if(c.eq.lower(j:j)) then
             c = upper(j:j)
             goto 3
          endif
2       continue
3       continue
        string(i:i) = c
1      continue
       return
       end
c*ds
c*ed
       subroutine oldupcase( string , l)
c      ==============================
c ---> uppercase translation <---
c
       character*132 string
       character*26  upper, lower , c*1
       DATA upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
       data lower/'abcdefghijklmnopqrstuvwxyz'/
c
       do 1 i=1,l
        c = string(i:i)
        do 2 j=1,26
          if(c.eq.lower(j:j)) then
             c = upper(j:j)
             goto 3
          endif
2       continue
3       continue
        string(i:i) = c
1      continue
       return
       end
c
c
c*ds
c*ed
      subroutine scan(ctext,r,ierr)
c ----------------------------------------------------------------------
c  scan - scannen einer real-zahl (version fuer ibm pc prof. fortran)
c  autor: g. egerer, datum: 26. 9.85, letzte aenderung: 31. 7.87
c ----------------------------------------------------------------------
 
c                           * eingabeparameter:
c                           * ctext  - textzeile (64 zeichen)
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
      integer   ctxtle
      parameter (ctxtle = 64)
      character*(ctxtle) ctext
      character*10       form
 
c                           * ausgabeparameter:
c                           * r      - real-zahl
c                           * ierr   - fehlercode
c     real    r
      integer ierr
 
c                           * benamte programmkonstanten:
      integer   cnil
      parameter (cnil = -32766)
 
c                           * programmvariablen:
      integer   icur, iscntb(9,0:7), istate, isymbl
 
 
c                           * initialisierungen:
c tabelle zur syntaxdefinition:
c
c symbol   s0   " "   "-"   "+"   "."   "e"   "d"   ziff
c ------+--------------------------------------------------
c zu-   |
c stand |
c --> 1 |         1     2     2     3                  8
c     2 |                           3                  8
c     3 |                                              9
c     4 |               5     5                        6
c     5 |                                              6
c end-  |
c zust. |
c     6 |         7                                    7
c     7 |         7
c     8 |         7                 9     4     4      8
c     9 |         7                       4     4      9
c
c s0 (symbolklasse 0): enthaelt alle symbole, die nicht gesondert aufge-
c                      fuehrt sind
c ziff               : enthaelt alle ziffern von "0" - "9"
c
      data ((iscntb(i,j), j = 0,7), i = 1,9)
     .   / cnil,    1,    2,    2,    3, cnil, cnil,    8,
     .     cnil, cnil, cnil, cnil,    3, cnil, cnil,    8,
     .     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    9,
     .     cnil, cnil,    5,    5, cnil, cnil, cnil,    6,
     .     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    6,
     .     cnil,    7, cnil, cnil, cnil, cnil, cnil,    7,
     .     cnil,    7, cnil, cnil, cnil, cnil, cnil, cnil,
     .     cnil,    7, cnil, cnil,    9,    4,    4,    8,
     .     cnil,    7, cnil, cnil, cnil,    4,    4,    9 /
c                           * ------------------------------------------
 
      icur = 1
      istate = 1
 1100 continue
      if ( icur .gt. ctxtle ) then
         if ( istate .ge. 6 ) then
c                           * istate ist ein endzustand
            ierr = 0
            form = '(BN,Eww.0)'
            write (form(6:7),'(i2)') ctxtle
            read(ctext,form) r
         else
c                           * istate ist kein endzustand
            ierr = icur
         endif
      else
c                           * ermittle index des aktuellen symbols
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
      end
c
c
c*ds
c*ed
      subroutine parse( zeile )
c     =========================
c ----------------------------------------------------------------------
c zerlegung der zeile nach incom manier
c die daten der letzten eingabedekodierung werden dabei ueberschrieben
c eine ggf. vorhandene restzeile bleibt erhalten
c restzeilen (; ...) beim parsen werden vollstaendig ignoriert !
c ----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c
       character*132 zeile, reslbf
       character*8   cmd
c-----------------------------------------------------------------------
c
       reslbf = reslin
       ioldbf = ioldc
       reslin = zeile
       call lowcase(reslin,132)
       ioldc  = 1
       call incom (cmd)
 
       reslin = reslbf
       ioldc  = ioldbf
 
      return
      end
c
c
c*ds
c*ed
      function iout()
c     ===============
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
         iout = iot
      return
      end
c
c*ds
c*ed
      function ierrs()
c     ===============
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
         ierrs = ierrr
      return
      end
c
c
c*ds
c*ed
      subroutine errsig(ierr,say)
c     ================================
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
      character*132 say, sayit
c ----------------------------------------------------------------------
c  error signalisierung
c ----------------------------------------------------------------------
      ierrr = ierr
      lsay = laenge(say,80,'$')
      sayit = say(1:lsay)
      write(6,*)'error:',ierr,' ',sayit
      return
      end
c
c
c*ds
c*ed
      subroutine clean
c     ================
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
          if(ierrr.ne.0) then
            write(6,*)'error condition:',ierrr,' cleaned'
            ierrr = 0
          endif
      return
      end
c
c
c*ds
c*ed
      subroutine pushn(pname)
c-----------------------------------------------------------------------
c  wert eintrag in den namensstack
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       if(inames.ge.minc) then
         call errsig(-1,'pushn failed, too many items! '//pname//' $')
         return
       endif
c
       inames = inames + 1
       vname(inames) = pname
 
        return
 
      entry pushr(rx)
c     ===============
c----------------------------------------------------------------------
c     zufuegen eines real-wertes auf den stack
c----------------------------------------------------------------------
      if(ipars.ge.minc) then
        call errsig(-1,'pushr failed, too many items ! $')
        return
      endif
 
      ipars = ipars+1
      rpar(ipars)   =  rx
      if(inpar(inames).eq.0) inapa(inames) =  ipars
      inpar(inames) =  inpar(inames) + 1
      iparn(ipars)  =  inames
c
      end
c
c
c*ds
c*ed
      function getval(pname,defval,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
       character*16 pnamc
c-----------------------------------------------------------------------
c
       inew      = 0
       wert      = defval
       getval    = defval
       call setudf(pname//' ',wert,ier)
       if(inames.le.0) return
c
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            if(inpar(i).gt.0) wert = rpar(inapa(i))
            inew = i
            lstpar = inapa(i)
            lstnam = i
            goto 20
         endif
10     continue
20     continue
c
      getval = wert
 
cc-setudf--
cc-setudf--
 
      return
c
      end
c
c
c*ds
c*ed
      function valnxt(defval,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       inew      = 0
       valnxt    = defval
       if(ipars.le.lstpar) return
c
      valnxt = rpar(lstpar+1)
      inew   = lstpar + 1
      lstpar = inew
 
      return
c
      end
c*ds
c*ed
      subroutine get1(p1)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
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
c
      end
c*ds
c*ed
      function intval(pname,idef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  integer
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
       character*16 pnamc
c-----------------------------------------------------------------------
c
       inew      = 0
       intval    = idef
       call setudf(pname//' ',wert,ier)
       if(inames.le.0) return
c
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            intval = rpar(inapa(i))*1.0000000001d0
            inew = i
            lstpar = inapa(i)
            lstnam = i
            goto 20
         endif
10     continue
20     continue
c
cc-setudf--
      wert = intval
cc-setudf--
      return
c
      end
c*ds
c*ed
      function intnxt(idef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  integer
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c

       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       inew      = 0
       intnxt    = idef
       if(ipars.le.lstpar) return
c
       intnxt = rpar(lstpar+1)*1.0000000001d0
       inew   = lstpar+1
       lstpar = lstpar
 
      return
c
      end
c*ds
c*ed
      character*8 function chrval(pname,cdef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  name (character)
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname, cdef
c-----------------------------------------------------------------------
c
       inew      = 0
       chrval    = cdef
       if(inames.le.0) return
c
       do 10 i=1,inames-1
         if(vname(i).eq.pname) then
            chrval = vname(i+1)
            inew = i
            lstnam = i+1
            goto 20
         endif
10     continue
20     continue
c
      return
c
      end
c*ds
c*ed
      character*8 function chrnxt(cdef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  name (character)
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname, cdef
c-----------------------------------------------------------------------
c
       inew      = 0
       chrnxt    = cdef
       if(inames.le.lstnam) return
 
       chrnxt = vname(lstnam+1)
       inew   = lstnam + 1
       lstnam = inew
 
      return
c
      end
c*ds
c*ed
      character*8 function cmdf()
c-----------------------------------------------------------------------
c  kommando extraktion aus dem incom parametercommon (character)
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       cmdf = comand
c
      return
c
      end
c*ds
c*ed
      character*8 function vnamef(iaddr)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  name (character)
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       vnamef = vname(iaddr)
       lstnam = iaddr
c
      return
c
      end
c*ds
c*ed
      subroutine settit(newtit)
c-----------------------------------------------------------------------
c  setzen des titels
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*80 newtit
c-----------------------------------------------------------------------
c
       ltit  = laenge(newtit,80,'$')
       title = newtit(1:ltit)
c
      return
c
      end
c*ds
c*ed
      character*80 function titlef()
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  titel
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
c
       titlef = title(1:80)
c
      return
c
      end
c*ds
c*ed
      function rparf(iaddr)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  real
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       rparf = rpar(iaddr)
       lstpar= iaddr
c
      return
c
      end
c*ds
c*ed
      function inapf(iaddr)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  zuordnung inapa
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       inapf = inapa(iaddr)
c
      return
c
c
      end
c
c
c*ds
c*ed
      function inpaf(iaddr)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  zuordnung inpar
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       inpaf = inpar(iaddr)
c
      return
      end
c
c
c*ds
c*ed
      function iparf()
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  zuordnung ipars
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       iparf = ipars
c
      return
      end
c
c
c*ds
c*ed
      function inamf()
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  zuordnung inames
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       inamf = inames
c
      return
c
      end
c*ds
c*ed
      function ifound(pname)
c-----------------------------------------------------------------------
c  logische optionserkennung  ifound wird=fundstelle sonst 0
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
       ifound    = 0
       if(inames.le.0) return
c
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            ifound=i
            lstnam=i
            goto 20
         endif
10     continue
20     continue
c
      return
c
      end
c*ds
c*ed
      logical function found(pname)
c-----------------------------------------------------------------------
c  logische optionserkennung
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c-----------------------------------------------------------------------
       character*8 pname
c-----------------------------------------------------------------------
c
      found = ifound(pname) .ne. 0
      return
c
      end
c
c
c*ds
c*ed
      logical function folgt(popt,pname)
c-----------------------------------------------------------------------
c  logische optionserkennung, wahr wenn popt pname folgt
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       character*8 popt, pname
c-----------------------------------------------------------------------
c
       folgt = .false.
       j     = ifound(pname)
       if(j.gt.0 .and. j.lt.inames) then
          folgt = popt .eq. vname(j+1)
       endif
       if(folgt) lstnam = j+1
c
      return
c
      end
c*ds
c*ed
      function intvno(ipnum,idef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  integer nach addr
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
c
       inew      = 0
       intvno    = idef
       if(ipars.lt.ipnum) return
c
       intvno = rpar(ipnum) * 1.0000000001d0
       inew      = ipnum
       lstpar    = ipnum
 
      return
c
      end
c*ds
c*ed
      function getvno(ipnum,adef,inew)
c-----------------------------------------------------------------------
c  wert extraktion aus dem incom parameterstack  real    nach addr
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
c
       inew      = 0
       getvno    = adef
       if(ipars.lt.ipnum) return
c
       getvno = rpar(ipnum)
       inew   = ipnum
       lstpar = ipnum
 
      return
c
      end
c*ds
c*ed
      subroutine getvec(pname,def1,def2,def3,vec ,inew)
c-----------------------------------------------------------------------
c  vektor-wert extraktion aus dem incom parameterstack mit rotation
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*16 pnamc
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       common/xroxxx/  xyorig(3), rotvec(3)
c-----------------------------------------------------------------------
       character*8 pname
       dimension defvec(3), vec(3)
c-----------------------------------------------------------------------
c
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
c
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
10     continue
20     continue
c
cc-setudf--
      call cappend(pname//' ','.1 ',pnamc)
      call setudf(pnamc//' ',vec(1),ier)
      call cappend(pname//' ','.2 ',pnamc)
      call setudf(pnamc//' ',vec(2),ier)
      call cappend(pname//' ','.1 ',pnamc)
      call setudf(pnamc//' ',vec(1),ier)
      call cappend(pname//' ','.3 ',pnamc)
      call setudf(pnamc//' ',vec(3),ier)
cc-setudf--
 
       call rotavc(vec,xyorig,rotvec)
      return
      end
c
c
c*ds
c*ed
      subroutine getnva(pname ,pardef ,parout, np, inew)
c-----------------------------------------------------------------------
c  extraktion meherere einem namen zugeordneter werte
c-----------------------------------------------------------------------
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
c-----------------------------------------------------------------------
       common/xroxxx/  xyorig(3), rotvec(3)
c-----------------------------------------------------------------------
       character*8 pname
       character*8 napp, pnamc*16
       dimension pardef(np), parout(np)
c-----------------------------------------------------------------------
c
       inew      = 0
       do 1 i=1,np
         parout(i) = pardef(i)
1      continue
       if(inames.le.0) goto 20
c
       do 10 i=1,inames
         if(vname(i).eq.pname) then
            napr = inpar(i)
            if(napr.gt.np) napr = np
            if(napr.gt.0) then
            do 11 j=1,inpar(i)
               parout(j) = rpar(inapa(i)-1+j)
11          continue
            inew = napr
            lstpar = inapa(i)+inpar(i)-1
            lstnam = i
            endif
            goto 20
         endif
10     continue
20     continue
 
       do 30 i=1,np
        write(napp,'(i3,5x)')100+i
        napp(1:1)='.'
        call cappend(pname//' ',napp,pnamc)
        call setudf(pnamc//' ',parout(i),ier)
30     continue
 
c
      return
      end
c*ds
c*ed
      subroutine rotavc( x, xyorig, rotvec )
c     ======================================
c     unterwirft den koordinatenvector x einer drehung um xyorigin
c     um den winkel rotang
c-----------------------------------------------------------------------
      implicit real*8   (a-h,o-z)
c
      dimension  x(3), y(3)
      dimension  xyorig(3), rotvec(3), r(3,3)
      logical cray
      data    degtor /0.017453292d0/
c
      y(1) = x(1) - xyorig(1)
      y(2) = x(2) - xyorig(2)
      y(3) = x(3) - xyorig(3)
      call rotax(rotvec, degtor, r, angdeg )
      call matvec( r, y, y)
      x(1) = y(1) + xyorig(1)
      x(2) = y(2) + xyorig(2)
      x(3) = y(3) + xyorig(3)
c
      return
      end
c*ds
c*ed
      subroutine unrot( x, xyorig, rotvec )
c     =====================================
c     unterwirft den koordinatenvector x einer drehung um xyorigin
c     um die drehung         rotang**-1
c-----------------------------------------------------------------------
      implicit real*8   (a-h,o-z)
c
      dimension x(3), y(3)
      dimension  xyorig(3), rotvec(3), r(3,3)
      data    degtor /0.017453292d0/
c
      y(1) = x(1) - xyorig(1)
      y(2) = x(2) - xyorig(2)
      y(3) = x(3) - xyorig(3)
 
      call rotax(rotvec, degtor, r, angdeg )
      call mattxp( r, r)
      call matvec( r, y, y)
 
      x(1) = y(1) + xyorig(1)
      x(2) = y(2) + xyorig(2)
      x(3) = y(3) + xyorig(3)
c
      return
      end
c
c
c*ds
c*ed
       subroutine rotax( ax, scale, rmat, angdeg )
c      ===========================================
c                       -->  ---->  <---   <---
c-----------------------------------------------------------------------
c drehmatrix rmat wird als drehung um die achse ax berechnet
c die laenge wird zur bestimmung des drehwinkels mit scale multipliziert
c-----------------------------------------------------------------------
 
       implicit real*8 (a-h,o-z)
       dimension ax(3), rmat(3,3)
       dimension en(3)
 
c-----------------------------------------------------------------------
 
       do 1 i=1,3
        do 1 j=1,3
           rmat(j,i) = 0.0d0
1      continue
 
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
2        continue
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
       end
c
c
c*ds
c*ed
       subroutine matmux( a, b, c )
c      ============================
c-----------------------------------------------------------------------
c      3x3 matrixmultiplikation
c      a * b ---> c
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), b(3,3), c(3,3)
       dimension h(3,3)
 
       do 10 i=1,3
          do 10 j=1,3
             s = 0.d0
             do 20 l=1,3
                s = s + a(i,l)*b(l,j)
20           continue
             h(i,j) = s
10     continue
 
       do 30 i=1,3
         do 30 j=1,3
            c(j,i) = h(j,i)
30     continue
 
       return
       end
c
c
c*ds
c*ed
       subroutine matvec( a, x, y )
c      ============================
c-----------------------------------------------------------------------
c      3x3 matrix - vector - multiplikation
c      a * x ---> y
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), x(3), y(3)
       dimension h(3)
 
       do 10 i=1,3
             s = 0.d0
             do 20 l=1,3
                s = s + a(i,l)*x(l)
20           continue
             h(i) = s
10     continue
 
       do 30 i=1,3
            y(i) = h(i)
30     continue
 
       return
       end
c
c
c*ds
c*ed
       subroutine mattxp( a, at)
c      =========================
c-----------------------------------------------------------------------
c      3x3 matrix - transposition
c      a(t)  ---> at
c-----------------------------------------------------------------------
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
       end
c
c
c*ds
c*ed
       subroutine matcpy( a, b )
c      =========================
c-----------------------------------------------------------------------
c      3x3 matrix - kopie
c      a     --->   b
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3), b(3,3)
 
       do 10 i=1,3
         do 10 j=1,3
           b(j,i) = a(j,i)
10     continue
 
       return
       end
c
c
c*ds
c*ed
       subroutine matone( a )
c      ======================
c-----------------------------------------------------------------------
c      3x3 matrix ---> einheitsmatrix
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       dimension a(3,3)
 
       do 10 i=1,3
         do 20 j=1,3
           a(j,i) = 0.d0
20       continue
         a(i,i) = 1.d0
10     continue
 
       return
       end
c
c
c*ds
c*ed
       subroutine axrot( rmat, ax)
c      ===========================
c-----------------------------------------------------------------------
c      rotationsachsendarstellung aus rmat
c      rmat ==~==> ax
c-----------------------------------------------------------------------
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
100    continue
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
200    continue
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
300    continue
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
999    continue
 
       ang = dacos(co)
       ax(1) = ax(1) * ang
       ax(2) = ax(2) * ang
       ax(3) = ax(3) * ang
 
c ------------- check sign ! ----------
       call rotax( ax, 1.0d0, rtest, angdeg)
       err =
     *       (rtest(1,2)-rmat(1,2))**2+
     *       (rtest(1,3)-rmat(1,3))**2+
     *       (rtest(2,3)-rmat(2,3))**2
 
       if(err.gt.1d-10) then
c        write(6,*)'++++++ rotate sign +++++'
         ang = -ang
         ax(1) =-ax(1)
         ax(2) =-ax(2)
         ax(3) =-ax(3)
       endif
 
       return
       end
c
c
c*ds
c*ed
       subroutine cross ( v1, v2, vc )
c      ===============================
c
c ** cross product of   v1 x v2  ==> vc
c
       real*8 v1,v2,vc, vch
       dimension v1(3),v2(3),vc(3),vch(3)
c
       vch(1) = v1(2)*v2(3) - v1(3)*v2(2)
       vch(2) = v1(3)*v2(1) - v1(1)*v2(3)
       vch(3) = v1(1)*v2(2) - v1(2)*v2(1)
c
       vc(1)  = vch(1)
       vc(2)  = vch(2)
       vc(3)  = vch(3)
 
       return
       end
c*ds
c*ed
       subroutine norm( xin, xout, n)
c      ==============================
c --- vektor norm nach skalarprodukt ---
       implicit real*8 (a-h,o-z)
       dimension xin(n), xout(n)
 
       s = 0.d0
       do 1 i=1,n
         s = s + xin(i)**2
1      continue
       if(s.gt.0.d0) s = 1.d0/dsqrt(s)
 
       do 2 i=1,n
         xout(i) = s * xin(i)
2      continue
 
       return
       end
c
c
c*ds
c*ed
       block data bdinco
c      =================
c
cray -------------------------------
       implicit real*8 (a-h,o-z)
cray -------------------------------
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
 
       logical cray
       common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
 
       common/xroxxx/  xyorig(3), rotvec(3)
 
       data iot/0/, inka/5/, ierrr/0/, ibild/0/, ioold/0/
       data cray/.false./
       data xyorig/3*0.d0/, rotvec/3*0.d0/
 
       end
c*ds
c*ed
       subroutine preplo
c      ================= prepare plotting (plotmachine)
c
       parameter(minc=40)
c --- minc = incom stack depth
c
       character*132 inline,reslin,title,rlbuf,inlbuf
       character*20 arglst,pmlist
       character*8 comand,vname
c
       real*8 rpar
       common/cincom/rpar(minc),inames,ipars
     *  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf
     *  ,lstpar, lstnam
       common/cincoc/comand,vname(minc),title,reslin,inline
     *  ,arglst(minc),pmlist(minc,2),rlbuf
 
      logical cray
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
 
       logical ptex
       logical paplo
       logical doplo
       logical fitplo
       logical paxis
       logical taxis
c      --- doplo = false  means: set parameters only ---
       character*80 option,xtext,ytext,tbuf,ttext,txtv,fmt,ttout
       character*8  dirnam(4),opart(8)
       character*8  codena, codefn, chrval
       logical      found, folgt
       dimension isymb(minc), icolo(minc)
       real*8 xdbl(*), ydbl(*), xcode(*)
       real*8 dble, getval, val
c
      data xmin/0./,xmax/4./,ymin/0./,ymax/2./,nkurv/0/
      data isymb/4,5,23,6,16,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21
     *          ,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41
     *          ,42,43/
      data ixmode/1/,txsize/.45/,sysize/.1/
      data framx/23./,framy/23./,yensca/0./,frlux/4./,frluy/4./
      data epscl/0.001/
      data ptex/.true./
      data paxis/.true./
      data taxis/.true./
      data icolo/minc * 1/,ifont/1/
      data icol0/1/
      data paplo/.true./
      data fitplo/.true./
      DATA OPART/'X=1     ','Y=1     ','I       ','A       ',
     *           'F=(3,1) ','M=2.0   ','U=1     ','        '/
      data txsizt/.3/,xtshft/0./,ytshft/0./
      data ttscrx/39.5/, ttscry/28.7/
c
c
c
c ----- parameter retrieving from stack -----
c      -----> decode by names
        xmin = getval('xmin    ',dble(xmin),inew)
        xmax = getval('xmax    ',dble(xmax),inew)
        ymin = getval('ymin    ',dble(ymin),inew)
        ymax = getval('ymax    ',dble(ymax),inew)
        framx= getval('framx   ',dble(framx),inew)
        framy= getval('framy   ',dble(framy),inew)
        frlux= getval('frlux   ',dble(frlux),inew)
        frluy= getval('frluy   ',dble(frluy),inew)
        if(found('text    ')) ptex = .true.
        if(found('txon    ')) ptex = .true.
        if(found('notext  ')) ptex = .false.
        if(found('txoff   ')) ptex = .false.
        if(found('noaxis  ')) paxis= .false.
        if(found('axis    ')) paxis= .true.
        if(found('notxaxis')) taxis= .false.
        if(found('txaxis  ')) taxis= .true.
        epscl =  getval('epscl   ',dble(epscl),inew)
        txsize=  getval('txsize  ',dble(txsize),inew)
        txsizt = getval('legsize ',dble(txsizt) ,inew)
        xtshtf = getval('legx    ',dble(xtshft) ,inew)
        ytshtf = getval('legy    ',dble(ytshft) ,inew)
        ifont  = intval('font    ',ifont  ,inew)
        icol0  = intval('color   ',icol0  ,inew)
        sysize = getval('sysize  ',dble(sysize) ,inew)
        do 3 i=1,inames
          j = inapa(i)
          if(vname(i).eq.'parplo  ') paplo=.true.
          if(vname(i).eq.'noparplo') paplo=.false.
          if(vname(i).eq.'ox      ') opart(1) = vname(i+1)
          if(vname(i).eq.'oy      ') opart(2) = vname(i+1)
          if(vname(i).eq.'oi      ') opart(3) = 'a'
          if(vname(i).eq.'oa      ') opart(3) = 'i'
          if(vname(i).eq.'of      ') opart(5) = vname(i+1)
          if(vname(i).eq.'om      ') opart(6) = vname(i+1)
          if(vname(i).eq.'ou      ') opart(7) = vname(i+1)
          if(vname(i).eq.'fits    ') fitplo = .true.
          if(vname(i).eq.'nofits  ') fitplo = .false.
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
3       continue
c
c
c ----- initialize plot (if first) -----
       if(ibild.eq.0) then
         call grstrt(35,8)
         write(6,*)'incom: preplo: grstrt(35,8) <======='
       endif
c
       ibild = ibild + 1
       isycnt= 0
c
c ---- set frame & scales ----
c      write(*,*)' set frames & scales .....'
       frxx = frlux + framx
       fryy = frluy + framy
       write(6,*)'xmin=',xmin
       write(6,*)'xmax=',xmax
       write(6,*)'ymin=',ymin
       write(6,*)'ymax=',ymax
       call grsclc(frlux,frluy,frxx,fryy)
       call grsclv(xmin,ymin,xmax,ymax)
       call grfont(ifont)
       call grnwpn(icol0)
       call grchrc(txsize,0.,16.)
c      --------------------------> set size of signs
 
c --- title ---
         xtx = xmin + 0.1 * (xmax-xmin)
         ytx = ymax
         ltext = 20./txsize
         if(ltext.gt.74) ltext=74
         call grtxt(xtx,ytx,ltext,title)
c
         xtx = 28.3 + xtshft
         ytx = 28.  + ytshft
 
       return
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
       entry setplt( codena, codefn, xcode, ncode )
c      --------------------------------------------
c      eingabe von plotparametern durch unterprogrammaufruf statt
c      von der kommandoebene aus
c      --------------------------------------------
          if(iout.gt.0) then
            write(6,*)'plot default setplt: ',codena,'  ',codefn
            write(6,'(1x,5f12.6)')(xcode(i),i=1,ncode)
          endif
 
          if(codena.eq.'xmin    ') xmin = xcode(1)
          if(codena.eq.'xmax    ') xmax = xcode(1)
          if(codena.eq.'ymin    ') ymin = xcode(1)
          if(codena.eq.'ymax    ') ymax = xcode(1)
          if(codena.eq.'framx   ') framx= xcode(1)
          if(codena.eq.'framy   ') framy= xcode(1)
          if(codena.eq.'frlux   ') frlux= xcode(1)
          if(codena.eq.'frluy   ') frluy= xcode(1)
          if(codena.eq.'text    ') ptex = .true.
          if(codena.eq.'txon    ') ptex = .true.
          if(codena.eq.'notext  ') ptex = .false.
c   ---> txon, txoff are from old version and not documented any more <-
          if(codena.eq.'txoff   ') ptex  = .false.
          if(codena.eq.'epscl   ') epscl = xcode(1)
          if(codena.eq.'txsize  ') txsize= xcode(1)
          if(codena.eq.'legsize ') txsizt= xcode(1)
          if(codena.eq.'legx    ') xtshft= xcode(1)
          if(codena.eq.'legy    ') ytshft= xcode(1)
          if(codena.eq.'font    ') ifont = xcode(1) +0.0001
          if(codena.eq.'sysize  ') sysize= xcode(1)
          if(codena.eq.'parplo  ') paplo=.true.
          if(codena.eq.'noparplo') paplo=.false.
          if(codena.eq.'ox      ') opart(1) = codefn
          if(codena.eq.'oy      ') opart(2) = codefn
          if(codena.eq.'oi      ') opart(3) = 'a'
          if(codena.eq.'oa      ') opart(3) = 'i'
          if(codena.eq.'of      ') opart(5) = codefn
          if(codena.eq.'om      ') opart(6) = codefn
          if(codena.eq.'ou      ') opart(7) = codefn
          if(codena.eq.'fits    ') fitplo = .true.
          if(codena.eq.'nofits  ') fitplo = .false.
c
          if(codena.eq.'symb    ') then
            nsy = 0
            do 3049 l=1,ncode
             nsy   = nsy   + 1
             if(nsy.gt.minc) goto 3029
             isymb(nsy) = xcode(l) + 0.0001
3049        continue
          endif
3029      continue
c
          if(codena.eq.'icolo   ') then
            nco = 0
            do 3059 l=1,ncode
             nco   = nco   + 1
             if(nco.gt.minc) goto 3039
             icolo(nco) = xcode(l) + 0.0001
3059        continue
          endif
3039      continue
c
       return
c      ======
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
       entry preaxe( xtext, ytext)
c      ---------------------------
 
c ---- prepare axes ----
c
c
       option = opart(1)//','//opart(2)//','//opart(3)//','//
     *                         opart(5)//','//opart(6)//','//opart(7)
       write(6,*)option
       lopt = 80
c      write(*,*)' make axes .....'
       if(taxis) then
          lxtxt = laenge(xtext,80,'$')
          lytxt = laenge(ytext,80,'$')
       else
          lxtxt = 0
          lytxt = 0
       endif
       if(paxis) call graxs(lopt,option,lxtxt,xtext,lytxt,ytext)
c
       return
c      ======
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
       entry pldata( xdbl, ydbl, ndbl)
c      ===============================
c --- resetscaling
       call grsclc(frlux,frluy,frxx,fryy)
       call grsclv(xmin,ymin,xmax,ymax)
       call grfont(ifont)
       call grchrc(txsize,0.,16.)
c      --------------------------> set size of signs
c ---  plot daten
       isycnt = mod(isycnt,minc) + 1
       isy    = isymb(isycnt)
       ico    = icolo(isycnt)
       call grnwpn(ico)
       if(isy.eq.0) then
         call grjmp( sngl(xdbl(1)), sngl(ydbl(1)) )
         do 200 i=2,ndbl
            call grdrw( sngl(xdbl(i)), sngl(ydbl(i)) )
200      continue
       else
         call grjmp( sngl(xdbl(1)), sngl(ydbl(1)) )
         do 201 i=2,ndbl
            call grjmps( sngl(xdbl(i)), sngl(ydbl(i)), isy )
201      continue
       endif
 
       return
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
       entry grcolo( icoln )
c      =====================
         isycc  = mod(isycnt,minc) + 1
         icolo(isycc ) = icoln
       return
c
c
       entry grwrit( ttext )
c      =====================
c
c ---- textpart ----
       if(.not.ptex) then
         write(*,*)'no textplotting due to notext'
         return
       endif
c
c ---- set textwindow ----
         call grsclc(0.,0.,ttscrx, ttscry)
         call grsclv(0.,0.,ttscrx, ttscry)
c
         call grchrc(txsizt,0.,16.)
         ltt = laenge(ttext,80,'$')
         call grtxt(xtx,ytx,ltt,ttext)
         ytx = ytx - 1.7 * txsizt
 
         return
c
       entry grwrva( txtv , val, fmt )
c      ===============================
c
c ---- textpart ----
       if(.not.ptex) then
         write(*,*)'no textplotting due to notext'
         return
       endif
c
c ---- set textwindow ----
         call grsclc(0.,0.,ttscrx, ttscry)
         call grsclv(0.,0.,ttscrx, ttscry)
c
         call grchrc(txsizt,0.,16.)
         tbuf = ' '
         lfmt = laenge(fmt,80,'$')
         tbuf='(a,'//fmt(1:lfmt)//',a)'
         ltxtv= laenge(txtv,80,'$')
         write(ttout,tbuf)txtv(1:ltxtv),val,'$'
         ltt = laenge(ttout,80,'$')
         call grtxt(xtx,ytx,ltt,ttout)
         ytx = ytx - 1.7 * txsizt
 
         return
c
      end
c
c
c*ds
c*ed
      function laenge( string, maxlen, delim)
c     =======================================
c     laenge von string, wobei das ende durch delim(1) gekennzeichnet is
c     ------------------------------------------------------------------
      character*132 string
      character*1  delim
 
      if(maxlen.gt.132) then
         n = 132
      else
         n = maxlen
      endif
 
      do 1 i=1,n
       if(string(i:i).eq.delim) then
         laenge = i-1
         return
       endif
1     continue
 
      laenge = n
 
      return
      end
c
c
c=======================================================================
c
c=================== io-add-ons for kfa-ibm-use =======================*
c
c*ds
c*ed
      subroutine openkfa(file,ikan,mode,irc)
c     ======================================
c
c    subroutine for cms filedef command
c    file  =    filename filetyp   oder filename.filetyp
c    ikan  =    fortran io-kanal
c    mode  =    'new' oder 'old' oder 'mod'
c    irc   =    return code
c
      character*24  file, fili, fm*1,cnum*4, line*80, mode*3
      line = ' '
c
      fili = file
      do 10 i=1,24
        if (fili(i:i).eq.'.') fili(i:i)=' '
10    continue
      ll    = laenge(fili,24,'$')
      fili  = fili(1:ll)
c
      write(cnum,'(i4)') ikan
c
      if(mode.ne.'new') then
         line = 'estate '//fili
         write(6,*)line
         call system(line)
         if (irc .ne. 0) then
            write(6,*)'openkfa file:',file,' not found!'
            return
         endif
      endif
c
      line = 'filedef '//cnum//' '//'disk '//fili
      if (mode.eq.'new') then
          line = 'filedef '//cnum//' '//'disk '//fili
     *           //' ( recfm f lrecl 132 '
      endif
      if (mode.eq.'mod') then
          line = 'filedef '//cnum//' '//'disk '//fili
     *         //' ( disp mod'
      endif
c
      write(6,*)'mode=',mode
      write(6,*)line
      call system(line)
      open (ikan)
c
      return
      end
c
c***********************************************************************
c---------- formula-interpreter-section --------------------------------
c***********************************************************************
c
 
 
       block data formb1
c      -----------------
       parameter(maxformlength=132)
       parameter(maxitemlength=80)
       parameter(maxnumstack=50)
       parameter(maxopstack=50)
       parameter(musrfstack=50)
       parameter(klammerinc=10)
       parameter(iplusprio=1, minusprio=1, multprio=2, idivprio=2)
       parameter(iexpprio=3, iuprio=7,komprio=0)
 
       parameter(nodelims=7)
       character*1 delims(0:nodelims)
       character*1 formula(0:maxformlength)
       character*1 item(0:maxitemlength)
       real*8      numstack(maxnumstack),valnum,degree
       character*4 opstack(maxopstack)
       integer     priostack(0:maxopstack)
       character*4 typ
       character*20 usrfstack(musrfstack)
       integer     tusrfstack
       integer     topnumstack,topopstack,klammerprio,actchar,len
       logical     ok, error, say
       common/formnu/numstack,valnum,priostack,topnumstack,topopstack,
     *      degree,
     *      tusrfstack,  klammerprio,actchar,len,litem,ok,error, say
       common/formch/formula,item,delims,typ,opstack,usrfstack
 
       data degree  /1.d0/
       data say     /.false./
       end
 
 
 
 
c*ds
c*ed
       subroutine foinit
c      ----------------
c====================================================================
c formelinterpreter kern
c -----------------
c====================================================================
 
       parameter(maxformlength=132)
       parameter(maxitemlength=80)
       parameter(maxnumstack=50)
       parameter(maxopstack=50)
       parameter(musrfstack=50)
       parameter(klammerinc=10)
       parameter(iplusprio=1, minusprio=1, multprio=2, idivprio=2)
       parameter(iexpprio=3, iuprio=7,komprio=0)
 
       parameter(nodelims=7)
       character*1 delims(0:nodelims)
       character*1 formula(0:maxformlength)
       character*1 item(0:maxitemlength)
       real*8      numstack(maxnumstack),valnum,degree
       character*4 opstack(maxopstack)
       integer     priostack(0:maxopstack)
       character*4 typ
       character*20 usrfstack(musrfstack)
       integer     tusrfstack
       integer     topnumstack,topopstack,klammerprio,actchar,len
       logical     ok, error, say
       common/formnu/numstack,valnum,priostack,topnumstack,topopstack,
     *      degree,
     *      tusrfstack,  klammerprio,actchar,len,litem,ok,error, say
       common/formch/formula,item,delims,typ,opstack,usrfstack
 
c ---- internal use ---
       character*4 op
       integer     prio
       real*8      val
       logical     yesno,compare
 
c ---- inits -----------------------------------------------------------
       ok          = .true.
       error       = .false.
       typ         = '    '
       topnumstack = 0
       topopstack  = 0
       klammerprio = 0
       actchar     = 0
       len         = lclen(formula,maxformlength)
cccw   write(6,*) 'len=',len
       delims(0)='('
       delims(1)=')'
       delims(2)='+'
       delims(3)='-'
       delims(4)='*'
       delims(5)='/'
       delims(6)='^'
       delims(7)=','
c      do i=0,maxitemlength
c         item(i) = ' '
c      enddo
       return
c ---------------------------------------------------------------------
       entry sayit(yesno)
c      ------------------
       say = yesno
       return
c ---------------------------------------------------------------------
       entry setdeg
c      ------------
       degree = datan(1.d0)/45.d0
       write(6,*)'evaluator set to degree'
       return
c ---------------------------------------------------------------------
       entry setrad
c      ------------
       degree = 1.d0
       write(6,*)'evaluator set to rad'
       return
c ---------------------------------------------------------------------
       entry show
c      ----------
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
c ---------------------------------------------------------------------
       entry putopstack( op, prio)
c      -----------------------------
 
cccw   write(6,*)'putopstack:',op,' ',prio
       if(topopstack.le.maxopstack) then
          topopstack = topopstack + 1
          opstack(topopstack)   = op
          priostack(topopstack) = prio+klammerprio
cccw   write(6,*)'topopstack=',topopstack
       else
          if(say) write(6,*)'error: opstack overflow!'
          error = .true.
       endif
 
       return
c ---------------------------------------------------------------------
       entry enternumstack( val )
c      --------------------------
 
       if(topnumstack.le.maxnumstack) then
          topnumstack = topnumstack + 1
          numstack(topnumstack) = val
       else
          if(say) write(6,*)'error: numstack overflow!'
          error = .true.
       endif
 
       return
c ---------------------------------------------------------------------
       entry checknum( n )
         if(topnumstack.lt.n) then
            if(say) write(6,*)'error: too few num operands!'
            error = .true.
            do i=topnumstack,n
              numstack(i) = 0.d0
            enddo
         endif
       return
c ---------------------------------------------------------------------
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
100    continue
       actchar = i
       if(litem.gt.maxitemlength) goto 999
       item(litem+1) = ' '
cccw   write(6,*)'getitem: item=',(item(j),j=0,maxitemlength)
       return
 
999    continue
       if(say) write(6,*)'getitem: item too long'
       error = .true.
       return
 
       end
 
 
c*ds
c*ed
       subroutine getitem
c      ------------------
 
       parameter(maxformlength=132)
       parameter(maxitemlength=80)
       parameter(maxnumstack=50)
       parameter(maxopstack=50)
       parameter(musrfstack=50)
       parameter(klammerinc=10)
       parameter(iplusprio=1, minusprio=1, multprio=2, idivprio=2)
       parameter(iexpprio=3, iuprio=7,komprio=0)
 
       parameter(nodelims=7)
       character*1 delims(0:nodelims)
       character*1 formula(0:maxformlength)
       character*1 item(0:maxitemlength)
       real*8      numstack(maxnumstack),valnum,degree
       character*4 opstack(maxopstack)
       integer     priostack(0:maxopstack)
       character*4 typ
       character*20 usrfstack(musrfstack)
       integer     tusrfstack
       integer     topnumstack,topopstack,klammerprio,actchar,len
       logical     ok, error, say
       common/formnu/numstack,valnum,priostack,topnumstack,topopstack,
     *      degree,
     *      tusrfstack,  klammerprio,actchar,len,litem,ok,error, say
       common/formch/formula,item,delims,typ,opstack,usrfstack
 
c --- internal use ---
       logical compare
       logical anklam
 
       character*(maxitemlength+1) citem
       character*(maxformlength+1) cform
       equivalence(citem,item(0))
       equivalence(cform,formula(0))
 
ccc    citem = ' '
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
c ---  suche bis zum naechsten delimiter ---
       citem = ' '
       call getword
c --- ist item ein unaerer operator ? ---
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
c --- ist item eine userfunction ? -----
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
c --- ist item eine zahl ? ----
cccw   write(6,*)'scan:',citem
       call scan(item ,valnum,ierr)
c --- try exp-num-representation --
       if(ierr.ne.0) then
cccw   write(6,*)'check exp-num item(litem)=',item(litem),'  ',litem
cccw   write(6,*)'formula(actchar)=',formula(actchar)
          lit = litem
          if((item(litem).eq.'e').or.(item(litem).eq.'d')
     *    .or.(item(litem).eq.'E').or.(item(litem).eq.'D')) then
            if((formula(actchar).eq.'+')
     *         .or.(formula(actchar).eq.'-')) then
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
400            continue
               call scan(item ,valnum,ierr)
               if(ierr.eq.0) then
                 actchar = l
               else
ccc              do l=litem+1,maxitemlength
ccc                item(l)=' '
ccc              enddo
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
c --- ist item eine referenz ? ----
       call extract(item,valnum,ierr)
       if (ierr.eq.0) then
         typ = 'num '
         call enternumstack(valnum)
         goto 100
       endif
c --- fehler:
       error = .true.
       if(say) write(6,*)'item:',item,' not decodable'
 
100    continue
 
       return
 
c ----------------------------------------------------------------------
       entry stackevaluate
c      -------------------
       j = topopstack-1
       if(j.lt.0) return
 
1      continue
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
            numstack(topnumstack-1)=
     *      numstack(topnumstack-1)+numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'-   ') then
            call checknum(2)
            numstack(topnumstack-1)=
     *      numstack(topnumstack-1)-numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'*   ') then
            call checknum(2)
            numstack(topnumstack-1)=
     *      numstack(topnumstack-1)*numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'/   ') then
            call checknum(2)
            numstack(topnumstack-1)=
     *      numstack(topnumstack-1)/numstack(topnumstack)
            topnumstack = topnumstack-1
            goto 200
         endif
         if(opstack(j).eq.'^   ') then
            call checknum(2)
            numstack(topnumstack-1)=
     *      numstack(topnumstack-1)**numstack(topnumstack)
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
            numstack(topnumstack)=dasin(numstack(topnumstack))
     *                            /degree
            goto 200
         endif
         if(opstack(j).eq.'acos') then
            call checknum(1)
            numstack(topnumstack)=dacos(numstack(topnumstack))
     *                            /degree
            goto 200
         endif
         if(opstack(j).eq.'atan') then
            call checknum(1)
            numstack(topnumstack)=datan(numstack(topnumstack))
     *                            /degree
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
200    continue
       j = j-1
       else
        j = j+1
        opstack(j) = opstack(topopstack)
        priostack(j)=priostack(topopstack)
        topopstack = j
        return
       endif
       goto 1
 
       end
 
c*ds
c*ed
       subroutine evaluate( f, val, ierr)
c      ----------------------------------
 
       parameter(maxformlength=132)
       parameter(maxitemlength=80)
       parameter(maxnumstack=50)
       parameter(maxopstack=50)
       parameter(musrfstack=50)
       parameter(klammerinc=10)
       parameter(iplusprio=1, minusprio=1, multprio=2, idivprio=2)
       parameter(iexpprio=3, iuprio=7,komprio=0)
 
       parameter(nodelims=7)
       character*1 delims(0:nodelims)
       character*1 formula(0:maxformlength)
       character*1 item(0:maxitemlength)
       real*8      val
       real*8      numstack(maxnumstack),valnum,degree
       character*4 opstack(maxopstack)
       integer     priostack(0:maxopstack)
       character*4 typ
       character*20 usrfstack(musrfstack)
       integer     tusrfstack
       integer     topnumstack,topopstack,klammerprio,actchar,len
       logical     ok, error, say
       common/formnu/numstack,valnum,priostack,topnumstack,topopstack,
     *      degree,
     *      tusrfstack,  klammerprio,actchar,len,litem,ok,error, say
       common/formch/formula,item,delims,typ,opstack,usrfstack
 
      logical cray
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
 
c --- internal use ---
       character*(maxitemlength+1) citem
       character*(maxformlength+1) cform ,f
       equivalence(citem,item(0))
       equivalence(cform,formula(0))
 
       do i=0,maxformlength
         formula(i) = ' '
       enddo
       do i=1,maxformlength
        if(f(i:i).eq.' ') goto 1
        formula(i-1)= f(i:i)
       enddo
1      continue
 
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
10     continue
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
       if(ierrs.ne.0 .and. iot.ge.1)
     *   write(6,*)'evaluate: ',(formula(i),i=0,len),' to ',val
 
       return
       end
 
 
 
       subroutine cappend( a, b, c)
c      ----------------------------
       character*16 a,b,c
       c = ' '
       do i=1,16
         if(a(i:i).ne.' ') then
          c(i:i) = a(i:i)
         else
          goto 1
         endif
       enddo
1      continue
       j = i
       do i=j,16
         l=i-j+1
         if(b(l:l).ne.' ') then
          c(i:i) = b(l:l)
         else
          goto 2
         endif
       enddo
2      continue
       j = i
       do i=j,16
        c(j:j) = ' '
       enddo
       end
 
       integer function lclen( string, m)
c      ----------------------------------
       character*1 string(0:m)
       do lclen=m,0,-1
         if(string(lclen).ne.' ') return
       enddo
       end
 
       logical function compare( s1, s2 )
c      ----------------------------------
       character*1 s1(*),s2(*)
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
       end
 
 
c*ds
c*ed
      subroutine extract(nam,val,ier)
c     -------------------------------
      character*1 nam(*)
      parameter (musevar=100)
      character*16 usenam
      real*8 val, useval
      common/usevar/useval(musevar),nousev
      common/usevac/usenam(musevar)
 
      parameter(mdepth=20)
      common/imargs/iargvs,ipmlst(mdepth),kanal(0:mdepth),ktop
 
      logical compare
 
      ier = 0
c --- look in uservars ---
      do i=1,nousev
        if(compare(nam,usenam(i))) then
          val = useval(i)
          return
        endif
      enddo
c --- number of doargs ---
      if(compare(nam,'noargs ')) then
        val = ipmlst(ktop)
        return
      endif
c --- user defined vals ---
      call usrextr(nam,val,ier)
      return
      end
 
c*ds
c*ed
      subroutine setudf(nam,val,ier)
c     ------------------------------
      character*16 nam
      parameter (musevar=100)
      character*16 usenam
      real*8 val, useval
      common/usevar/useval(musevar),nousev
      common/usevac/usenam(musevar)
      logical cray
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
      logical compare
 
      ier = 0
c --- look in uservars ---
cccw  write(6,*)'setudf: nam=',nam,' val=',val
      do i=1,nousev
        if(compare(nam,usenam(i))) then
          useval(i) = val
          if(iot.ge.1)
     *    write(6,*)'uservar:',usenam(i),' is set to:',useval(i)
          return
        endif
      enddo
cccw  write(6,*)'create new udf'
c --- create new var ---
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
1       continue
        useval(nousev) = val
        i = nousev
        if(iot.ge.1)
     *  write(6,*)'uservar:',usenam(i),' is created:',useval(i)
        return
      else
        write(6,*)'uservar set: buffer is full ! error'
        ier = 200
      endif
      return
      end
 
c*ds
c*ed
      subroutine clrudf(nam)
c     ----------------------
      character*16 nam
      parameter (musevar=100)
      character*16 usenam
      real*8 val, useval
      common/usevar/useval(musevar),nousev
      common/usevac/usenam(musevar)
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
      end
 
c*ds
c*ed
      subroutine shwudf
c     -----------------
      character*16 nam
      parameter (musevar=100)
      character*16 usenam
      real*8 val, useval
      common/usevar/useval(musevar),nousev
      common/usevac/usenam(musevar)
      logical compare
 
      if(nousev.eq.0) then
        write(6,*)'no uservars defined!'
      else
        write(6,*)'defined uservars:'
        do i=1,nousev
           write(6,*)usenam(i),useval(i)
        enddo
      endif
      return
      end
c*ds
c*ed
      subroutine adpicn
c     -----------------
      logical cray
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray
 
      ibild = ibild + 1
c     -----------------
      return
      end
