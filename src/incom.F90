!***********************************************************************
!       incom-system                                                   *
!       ============                                                   *
!       by  michael monkenbusch                                        *
!       institut fuer festkoerperforschung des                         *
!       forschungszentrum juelich, d-5170 juelich, po-box 1913, frg    *
!       iff115 @ djukfa11                                              *
!                                                                      *
!***********************************************************************
       subroutine incom( cmd )
       use cincom
       use cincoc
       use icpathes
       use imargs
       use cmargs
       use xoutxx
       use xroxxx
       use constants
       use outlev
       implicit none
!      -----------------------
!
! ---- command line decoder ----
!
!
! ---- separates commands parameternames & parameters ----
!
! ---- commandline format may be :
!      command <name parameter name parameter < ;command ....>>
!
! ---- parameter separator is blank !!!
! ---- command separator is ;
!
! ---- outputlevel
!
!
       character*8   glabel,xlab,itypc,vnamef
       character*1024 buf,bla132
       character*1   blank,csep,numtst(14),numeqv*14,inread*1024
       character*8   cmd
       equivalence  (numeqv,numtst(1))
       logical :: name, lmakro = .false., found, fileda

       character*1024 ma_fil
       integer       ilma, i, j, k, l, ii, irc, ipmlen, isum, ioold
       integer       ier, ierr, ioldna, inew, iival

       integer :: init_run = 1
       real*8 val
       real getval
!
! ---- initialisations & conventions ----
!
       if(init_run.eq.1) then
			data_path =  './'
			save_path =  './'
			makro_path = './'
!			call getenv('HOME', home)
            call getenv('PWD', PWD)
			datreat_path = PATH_TO_LOCAL_DATREAT//'/'
!           datreat_path = PWD(:index(PWD,'/',BACK = .TRUE.))
			do i=20-1,0,-1
			history(i)=''
			enddo

			! to load a initialisation macro if it exists
			inquire(file=trim(datreat_path)//'makros/'//'initdatr',exist=fileda)
			if(fileda) then      ! for a init file
				ioldc=1
				reslin='initdatr  '
			endif
			init_run = 0
       endif

!
       csep = ';'  !              = --> command separator
       blank= ' '  !              = --> parameter separator
       numeqv = '(.+-0123456789' !=> to reckognize the beginning of a number
       bla132 = ' '
!

 8888  continue ! startingpoint of incom loop
!      --------> reentry-point
       do i=1,minc
        vname(i) = ' '
        rpar(i)  = 0
        iparn(i) = 0
        inpar(i) = 0
        inapa(i) = 0
       enddo
       comand = ' '

       lstpar    = 0
       lstnam    = 0
!
! ---- error-response -----------------------------------------
       if(ierrr.ne.0) then
         write(6,*)' error return code:',ierrr
         ioldc = 0 !-----> if error forget rest of commandline (reslin)
           if(ktop.ne.0) then
             do k=ktop,1,-1
               close(kanal(k))
             enddo
             ktop = 0
!            --------> go back to input from keyboard
           endif
       endif
       ierrr = 0
!
!
! ---- input a line or take nodecode rest of previous line ----
       if(ioldc.eq.0) then
         goto 1002
 1001  continue

         if(iot.gt.-3) write(6,*)'end of data '
         close(kanal(ktop))
         if(ktop.eq.0) then
            open(kanal(ktop))
         else
           ktop = ktop-1
!                      --------> restore old reslin state
           ioldc = iolbuf
           reslin= rlbuf
           ipmls = 0
           if(ktop.eq.0)write(6,*)' command input back to keybord'
!			write(6,*)' #-###->'
           if(ioldc.ne.0) then
             inline = reslin

             goto 1011
           endif
         endif
 1002    continue
        if(ktop.eq.0) then
!            write(6,*)' #-->'
            write(6,'(a)',advance='no') prompt
            ipmls = 0
		endif
        read(kanal(ktop),"(a)",end=1001) inread
        ! ---- last command from history ----
       	ii=len_trim(inread(2:))
		if(inread(1:1).eq.'_') then
			do i=0,20
				if (index(history(i),inread(2:ii+1)).eq.1) then
					inread=history(i)
					exit
				endif
			enddo
		elseif(ktop.eq.0 .and. index('history',inread(:4)).ne.1) then
			do i=20-1,0,-1
			  history(i+1)=history(i)
			enddo
			history(0)=inread
		endif  ! go on with new command if found

		if(ktop.ne.0.and.iot.gt.-1)write(6,*)'! '//trim(inread)
        if(ktop.eq.0.and.iot.gt.-2)write(6,*)'> '//trim(inread)
			inline = inread

		else
			inline = trim(reslin)
			ioldc=0
		endif
 1011   continue


		ioldc = 0
		ioldna =0
!
! ---- look for a commnd separator ----
       if(SCAN(inline,csep).gt.0) then  ! falls csep vorkommt ist scan >0 , falls nicht =0
          ioldc = 1
!         ---------------------> store second command on reslin
          reslin = inline(SCAN(inline,csep)+1:len(inline))
		  inline = inline(1:SCAN(inline,csep)-1)//' '

       endif
!
! ---- remove any leading blanks ----
       inline = ADJUSTL(inline)
		 !
!
! ---- separate the command ----   SCAN ('FORTRAN', 'R')=3
       j=SCAN(inline,' ')-1        ! findet erstes blank j ein davor
       if(j.gt.len(comand))then
       	 comand = inline(1:8)
		else
			comand = inline(1:j)
		endif
!
		 cmd = comand

! ---- dont analyse further if it is a pathdefinition ----
        if(comand.eq.'path') then
           if((inline(6:6).eq.'?') .or. (len_trim(inline).eq.len_trim('path') )) then
                write(6,*)'read data from : '//TRIM(data_path)
				write(6,*)'save data to   : '//TRIM(save_path)
				write(6,*)'load makro from: '//TRIM(makro_path)
				write(6,*)'additional to path: '//TRIM(datreat_path)//'makros/'
				write(6,*)'change with datapath/savepath/macrpath  "PATH"'
             goto 8888
           endif

           goto 8888
         endif
		if(comand.eq.'datapath') then
            if((inline(10:10).eq.'?') .or. (len_trim(inline).eq.len_trim('datapath') )) then
              write(6,*)TRIM(data_path)
              goto 8888
            endif
		    if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
				data_path=trim(inline(10:))//'/'
			else
				data_path=trim(inline(10:))
			endif
        if(iot.ge.0) write(6,*)'datapath: ',TRIM(data_path)
           goto 8888
         endif
         if(comand.eq.'savepath') then
           if((inline(10:10).eq.'?').or. (len_trim(inline).eq.len_trim('savepath') )) then
             write(6,*)TRIM(save_path)
             goto 8888
           endif
           if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
					save_path=trim(inline(10:))//'/'
				else
					save_path=trim(inline(10:))
				endif
           if(iot.ge.0) write(6,*)'savepath: ',TRIM(save_path)
           goto 8888
         endif
         if(comand.eq.'macrpath') then
           if((inline(10:10).eq.'?').or. (len_trim(inline).eq.len_trim('macrpath') )) then
             write(6,*)TRIM(makro_path)
             goto 8888
           endif
            if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
					makro_path=trim(inline(10:))//'/'
				else
					makro_path=trim(inline(10:))
				endif
           if(iot.ge.0) write(6,*)'macrpath: ',TRIM(makro_path)
           goto 8888
         endif

! ---- dont analyse further if it is a title ----
         if(comand.eq.'tit     '.or.comand.eq.'title   ') then
           if(comand.eq.'tit     ') title = trim(inline(5:))
           if(comand.eq.'title   ') title = trim(inline(7:))
           goto 8888
         endif
! ---- dont analyse further if it is a message ---
         if(comand.eq.'msg     ') then
           write(6,*)trim(inline)//'                 '
           goto 8888
         endif
! ---- dont analyse further if it is a blank line ---
         if(comand.eq.'        ') then
           !write(6,*)inline
           goto 8888
         endif
! ---- dont analyse further if it is a comment ---
         if(comand.eq.'c       ') then
           goto 8888
         endif
! ---- dont analyse further if it is a label ----
         if(comand(1:1).eq.':') then
           if(iot.gt.2)write(6,*)'label encountered ',comand
           goto 8888
         endif
!
! ---- look for keywords & parameters ----
!
       inames = 0
       ipars  = 0
       iargs  = 0
       iargvs = 0
       j = j+1    !erstes blank nach comand da j letzter char von comand
       do while (len_trim(inline(j+1:len(inline))).ne.0 )   ! schleife parameterersetzung
			i = j ! i = aktuelles erstes blank
			do while (inline(i+1:i+1).eq.blank)
				 i=i+1       ! ! if there are more than single blank
			enddo
			l = i+1
!   ---- look for the end of the currend item -----
			j=i+SCAN(inline(i+1:len(inline)),blank)      ! j = next blank
			if (j.eq.i) j=len(inline)+1       ! no further blanks till end j=ende+1 virtuell blank
			buf = inline(i+1:j-1)
	! --- prepare the textlist of formal makro arguments --
			iargs = iargs + 1
			arglst(iargs) = buf(1:20)
	! --- decode for makro argument values replace strings ! ----
			if(ktop.ne.0) then
				if(ipmlst(ktop).ne.0) then
					if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
					do k=1,ipmlst(ktop)
						if(iot.gt.3)write(6,*)buf
						call creplace(buf,pmlst(ktop,k,1)//' ',pmlst(ktop,k,2)//' ',' ')
						if(iot.gt.3)write(6,*)buf
					enddo
					if(len_trim(buf).eq.0) cycle
				endif
			endif
! ----   discriminate between name & parameter ----
			name = .true.
			! starts with

			if((scan('(+-0123456789',buf(1:1)).gt.0) .or. &         !number starts with number or +-
		      ('.'.eq.buf(1:1) .and. scan('(0123456789',buf(2:2)).gt.0)  .or. &  ! number omitted 0 like .534
				((scan(buf(1:len(buf)),'()+-*/^').gt.0 .and.buf(1:2).ne.'./'.and.&
						 scan('/^*.',buf(1:1)).eq.0.and. scan('/^*.',buf(len_trim(buf):len_trim(buf))).eq.0 ))  & ! this is a formula   !
				)   name = .false.       ! then it is a number or should be evaluated as formula
			if(.not.name) then ! .not.name = zahl oder formel
				call evaluate( buf//' ', val, ierr)
				if(ierr.eq.0) then
					name = .false.
				else
					name = .true.
!					comand = 'f-error!'
!					cmd    = comand  (inline(1:1).ne.'&') .and.
!					write(6,*)'syntax error in: '//trim(buf)
!					write(6,*)'       in line : '//trim(inline)
				endif
			endif
!
			iargvs = iargvs+1
			if(name) then
				inames = inames + 1
				if(inames.gt.minc) goto 999
				vname(inames) = buf(1:8)
				argvals(iargvs) = buf(1:1024)
				inapa(inames) = 0
			else
				ipars = ipars + 1
				if(ipars.gt.minc) goto 999
				rpar(ipars) = val
				iival = val+1d-11
				if(dabs(val-iival).lt.3d-11) then
					write(argvals(iargvs),'(i132)') iival
				else
					write(argvals(iargvs),'(e132.12)') rpar(ipars)
				endif
				iparn(ipars) = inames
				if(inames.ne.ioldna) then
					inapa(inames) = ipars
					ioldna = inames
				endif
			endif
      enddo ! ende schleife parameterersetzung

!
! ---- prepare inpar (no. of parameters / name) ----
       do i=1,inames
         isum = 0
         if(ipars.ne.0) then
            do j=1,ipars
             if(iparn(j).eq.i) isum = isum + 1
         	enddo
			endif
         inpar(i) = isum
		  enddo
!
!
       if(lmakro) goto 9888

! --------------- decoding of opaque commands -------------------------
!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! cms - kommando absetzen                             | cms <cmd >
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'cms     '.or.comand.eq.'sys     ') then
!        same code as  bottom shell command            ---
			buf=ADJUSTL(inline(5:))
         if(ktop.ne.0) then
         if(ipmlst(ktop).ne.0) then
          if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
           	do k=1,ipmlst(ktop)
					ipmlen=0
					do while (index(buf,trim(pmlst(ktop,k,1))).gt.ipmlen )
						ipmlen=index(buf,trim(pmlst(ktop,k,1)))
						buf=buf(1:ipmlen-1)//trim(pmlst(ktop,k,2))//trim(buf(ipmlen+len_trim(pmlst(ktop,k,1)):))
						ipmlen=ipmlen+len_trim(pmlst(ktop,k,2))
					enddo
				enddo
          endif
         endif
			call system(trim(buf))
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! history anzeigen                             | history    hist
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(index('history',comand(:4)).eq.1) then
!                    ----
         if (index(vname(1),'clear').gt.0) then
				do i=20-1,0,-1
				history(i)=''
				enddo
				goto 8888
			endif
			open(10,file=trim(makro_path)//'history')
			do i=20-1,0,-1
				if (len_trim(history(i)).gt.0) then
					write(*,*) i,' ',trim(history(i))
					write(10,*) ,trim(history(i))
				endif
			enddo
			close(10)
         goto 8888
       endif
!

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! cd - kommando absetzen                             | cd <dir >
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'cd     ') then
!                    ---
         call chdir(trim(inline(4:)))
         call system('pwd')
			goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'..      ') then
!                    ---
         call chdir('..')
         call system('pwd')
			goto 8888
       endif
       if(comand.eq.'datreat ') then
!                    ---
         call chdir(trim(datreat_path))
         call system('pwd')
			goto 8888
       endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! output"level" setzen:                               | iout <iout>
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'iout    ') then
!                    ----
         ioold= iot
         iot  = rpar(1)*1.0000000001d0
         if(vname(1).eq.'old     ') iot = ioold
         write(6,*)'outputlevel is now ',iot
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! if                                                  |
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'if      ') then
!                    ----------------> if-bedingung (numerisch)   !
         itypc = vnamef(1)
         if(itypc.eq.'=       ') then
           if(rpar(1).eq.rpar(2)) then
             j = index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'>       ') then
           if(rpar(1).gt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<       ') then
           if(rpar(1).lt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<=      ') then
           if(rpar(1).le.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'>=      ') then
           if(rpar(1).ge.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         if(itypc.eq.'<>      ') then
           if(rpar(1).ne.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
           endif
           goto 8888
         endif
         write(6,*)'if operator ',itypc, ' unknown'
         ierrs = 1
         goto 8888
       endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! realisation goto                                    | goto
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'goto    ') then
!                    ----------------> comment feature for makros !
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
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! return from makro                                   | return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'return  ') then
!                    ----------------> return from makro
         if(ktop.gt.0) then
           close(kanal(ktop))
           ktop = ktop - 1
         endif
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! exit all makros                                     | exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'exit    ') then
!                    ----------------> exit   from makro
         if(ktop.gt.0) then
           do k=ktop,1,-1
             close(kanal(ktop))
           enddo
           ktop = 0
         endif
         goto 8888
       endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! evaluation info                                     | ??
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'??      ') then
!                    ----------------> evaluation info !
         do i=1,ipars
          write(6,*) rpar(i)
         enddo
         goto 8888
       endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! evaluation and  print                               | print
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'print   ') then
!                    ----------------> print evaluated numbers !

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
!
!

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! (user)vars zeigen                                   | vars?
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'vars?   ') then
!                    ----------------> uservars display
         call shwudf
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! formelauswertung protokoll                          | say
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'ftrace  ') then
!                    ----------> set deg-mode
         call sayit(.not.found('off     '))
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! degmode setzen                                      | setdeg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'setdeg  ') then
!                    ----------------> set deg-mode
         call setdeg
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! radmode  setzen                                     | setrad
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'setrad  ') then
!                    ----------------> set deg-mode
         call setrad
         goto 8888
       endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! uservariable setzen                                 | set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'set     ') then
!                    ----------------> set uservar

			do i=1,inames
				if (inpar(i).eq.1) then
           			call setudf(vname(i)//' ',rpar(inapa(i)),ier)
			  	else
					write(6,*)'error in value cannot be evaluated     '//trim(vname(i))
			  endif
         enddo
         goto 8888
       endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! uservariable entefrnen                              | clr
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'clr     ') then
!                    ----------------> set uservar
         do i=1,inames
           call clrudf(vname(i)//' ')
         enddo
         goto 8888
       endif
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! rotationsparameter setzen                           | setrot <xyz> <r>
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'setrot  ') then
!                    ----------------> set rotation-parameters !
           call get6( rotvec(1), rotvec(2), rotvec(3),                  &
     &                xyorig(1), xyorig(2), xyorig(3))
           xyorig(1) = getval('x       ',xyorig(1),inew)
           xyorig(2) = getval('y       ',xyorig(2),inew)
           xyorig(3) = getval('z       ',xyorig(3),inew)
           rotvec(1) = getval('rx      ',rotvec(1),inew)
           rotvec(2) = getval('ry      ',rotvec(2),inew)
           rotvec(3) = getval('rz      ',rotvec(3),inew)
           write(6,*)'------------------------------------------------'
           write(6,'(a,3f12.6/a,3f12.6)')' rotation origin =',xyorig,   &
     &                                   ' rotation vector =',rotvec
           write(6,*)'------------------------------------------------'
         goto 8888
       endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ausgang                                             | q
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'q       ') then
!         old gr stuff
!                    ----------------> exit the program
!         if(ibild.gt.0) then
!            write(6,*)'calling grende ...'
!           !call grende
!         endif
         vname(1)='cleanio '
         call gplot
         write(6,*)'exit by q, bye...'
         stop
       endif
!
       return
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! ---- stack overfow error trap ----
  999  continue
       write(6,*)' incom(',comand,') :  stack overflow   reenter'
       ioldc = 0
       ierrr = 2
       goto 8888

!
!
!
! --- makro - aktivator - entry ----------------------------------------
!
       entry makro (cmd)
!      ===========
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- try to interpret non-known commands as makro input files ---
! ------> any makro input file must have the keyword  <makro>
!         at the first position
       if(ktop.ge.mdepth) then
         write(6,*)' too many stacked makros : ',comand
         ierrr = 3
         return
       endif
       ktop = ktop + 1

!!! Modification for two makro sources: global, individual !!!
!!! try first the individual makro as we had it already    !!!

       open(kanal(ktop),file=comand,status='old',form='formatted',err=99991)
       read(kanal(ktop),9998,end=99991) inline
 9998  format(a80)
       if(inline(1:6).eq.'makro ') goto 19999

!!! ok this was no makro !!!
       close(kanal(ktop))
99991  continue

! maybe on user-std-makropath
       ma_fil = trim(makro_path)//comand
       ilma = LEN_TRIM(ma_fil)
       open(kanal(ktop),file=ma_fil(1:ilma),status='old',form='formatted',err=99992)
       read(kanal(ktop),9998,end=99992) inline
       if(inline(1:6).eq.'makro ') goto 19999

!!! ok this was no makro !!!
       close(kanal(ktop))
99992  continue
!
! -------> try whether in the global makro section file is a makro file
       open(kanal(ktop),file=trim(datreat_path)//'makros/'//comand,status='old',form='formatted',err=9999)
       read(kanal(ktop),9998,end=9999) inline
       if(inline(1:5).ne.'makro') goto 9999
! -------> try whether file is a makro file
!
19999  continue
! --- ok this is a makro file ! save first the parametervalues ---
         do i=1,minc
           pmlst(ktop,i,2)=' '
         enddo
         do i=1,iargvs
				pmlst(ktop,i,2) = ADJUSTL(argvals(i))
				if(iot.gt.4) write(6,*)'argvals : '//'##'//trim(ADJUSTL(argvals(i)))//'##'//pmlst(ktop,i,2)//'##'
         enddo
			ipmlst(ktop) = iargvs
         do i=iargvs+1,minc
           pmlst(ktop,i,2) = ' '
			enddo
! --- save reslin form the current keyboard command ---
         iolbuf = ioldc
         rlbuf  = reslin
! --- decode the first makro line to get the parameternames
         ioldc = 1
         reslin = inline
         lmakro = .true.
         if(iot.gt.-3)write(6,*)'this is a makro ...'
         goto 8888
 9888    continue
         lmakro = .false.
         if(iot.gt.-3)write(6,*)'go on with it ...'
         if(iargs.ne.0) then
           	do i=1,iargs
            	if(iot.gt.3) write(6,*)'arglst:',arglst(i)
            	pmlst(ktop,i,1) = arglst(i)
         	enddo
			endif
         ipmlst(ktop) = iargs
!        ipmls = iargs
!
       return
!      ---------> read now commandlines from makro file

 9999  continue
		ktop = ktop - 1 ! back to old ktop macrolevel if it was no macro
 !###################################################
!        if command not known till now we try to sent it to the shell
			buf=inline
         if(ktop.ne.0) then
         if(ipmlst(ktop).ne.0) then
          if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
           	do k=1,ipmlst(ktop)
					ipmlen=0
					do while (index(buf,trim(pmlst(ktop,k,1))).gt.ipmlen )
						ipmlen=index(buf,trim(pmlst(ktop,k,1)))
						buf=buf(1:ipmlen-1)//trim(pmlst(ktop,k,2))//trim(buf(ipmlen+len_trim(pmlst(ktop,k,1)):))
						ipmlen=ipmlen+len_trim(pmlst(ktop,k,2))
					enddo
				enddo
          endif
         endif
		call system(trim(buf))
		If (irc .ne. -1) return

!        back from (s)hell
!###################################################
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(comand.eq.'f-error!') then
!                    ----------------> fehler-reaktion !
         write(6,*)'!! commandline in error --> ignored !!'
         ierrr = 111
         goto 8888
       endif
!
       do k=ktop,1,-1
         close(kanal(k))
       enddo
       ktop = 0
       write(6,"(' program name ''',a8,''' not known')")comand
       ierrr = 4

		 return

      END
!     END of incom!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


