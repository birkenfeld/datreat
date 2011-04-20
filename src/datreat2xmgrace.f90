subroutine gplot ()
  use cincom
  use cincoc
  use constants
  !      ================  scan-plotting
  !
  !
  !       implicit none
  parameter(mth=40,mtpar=40,mtcal=40)
  !      parameter(mkurv=minc)
  !                ---------> max. no of curves to be plotted
  ! ---  maximum scan length
  !
  parameter(mwert=1024,mbuf=200,mpar=200)
  ! ---- common containing all the scans ----
  character*80 name,xname,yname,napar,coment*80
  common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),yerror(mwert,mbuf), xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),&
       numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),napar(mpar,mbuf),nopar(mbuf)
  !
  common/selist/isels(mbuf),ifits(mbuf),nsel,numpls
  common/fslist/isfits(mbuf),nfsel
  !
  character*8 thenam,thparn
  common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),thparx(mtpar,mtcal),thpsca(mtpar,mtcal),&
       nthtab(mtcal),ntheos,multflg(mtcal)
  common/therrc/therro(mtpar,mtcal)
  !
  character*6 cnum
  !

  character*1024 inlbuf
  logical cray
  common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray
  common/xroxxx/  xyorig(3), rotvec(3)
  !
  !
  common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
  !
  dimension x(mwert),y(mwert),yl(mwert),yh(mwert),irecv(minc),isymb(minc),irecn(minc),icolo(minc),ifrec(minc)
  dimension e(mwert)
  ! ###################################################Anfang RB
  logical fitplo,errplo,legend,closeGrace,noautoscale,gsave,clearplo
  logical comon,writefit
  logical found
  real*8 params_all(mpar),paramsf_all(mpar)
  character*80 napar_all(mpar),naparf_all(mpar)
  character*80 commentar*80, commentarfit*80
  !  prepare grace #############################
  CHARACTER*1900 buf,buf2
  character*128 savename          ! Buffer for writing
  character*8 ii,jj               !buffers for colum and row
  character*8 g0s,g0si           ! short for graph and set in grace g0s1
  character*1 quot
  integer sets,fb                ! Anzahl der geplotteten sets in grace
  character*1024 startupsetting  ! startbegfehle für grace layout
  !
  INTEGER GraceOpenf, GraceIsOpenf
  EXTERNAL GraceOpenf, GraceIsOpenf
  EXTERNAL GraceRegistErerrorFunctionf
  EXTERNAL GraceCommandf,GraceClosef,GraceClosePipef
  EXTERNAL MyError, GraceFlushf
  !
  !CALL GraceRegistErerrorFunctionf (MyError)

  ! ende prepare grace #############################
  quot=achar(34)            ! this is "
  !       set some startoptions
  sets=0            !  number of sets in grace set to first (write in 0 and move to next) 0 is empty
  data fitplo /.true./
  data errplo /.true./
  data legend /.true./
  data fb /0/
  data closeGrace /.false./
  data gsave /.false./
  data comon /.false./
  data writefit / .false./
  data noautoscale / .false./
  data g0s /'g0.s'/
  data startupsetting /' '/
  data clearplo /.false./
  data savename /'lastgrace'/
  if(vname(1).eq.'start   ' ) then
     write(6,*) 'setting grace startupoptions'
     startupsetting=inline(SCAN(inline,' ')+1:len_trim(inline))//';'//trim(reslin)                       ! first blank in inline
     reslin=''  ! rest of line should be startupsetting
     startupsetting=startupsetting(SCAN(startupsetting,' ')+1:len(startupsetting))  ! second blank in inline
     return
  endif
  if(vname(1).eq.'cleanio ' ) then
     write(6,*) 'delete temporary .iofile.tmp'
     call system('rm -f .iofile.tmp*')
     return
  endif
  if(vname(1).eq.'comand  ' .or. vname(1).eq.'cm      ' ) then
     !direct command to grace
     !full command is in inline, so we need to strip gplot and comand all after is sent to grace as comand
     ! erreicht durch 2 mal suchen nach space da die commandos dadurch getrennt sind
     buf=inline(SCAN(inline,' ')+1:len(inline))       ! gplot wegschneiden
     if ( GraceIsOpenf () .NE. 0 ) then
        CALL MyGraceCommandf (trim(buf(SCAN(buf,' ')+1:len(buf))))
        return
     else
        write(*,*)'no grace pipe found ; is there already a plot opened?'
     endif
     return
  endif
  if(vname(1).eq.'save    ' .or. vname(1).eq.'sv      ' ) then
     if (inames.gt.1) then
        savename = vname(1+1)
        gsave=.true.
     else
        savename='lastgrace'
        gsave=.not.gsave
     endif
  endif
  !       -----> decode by names
  do i=1,inames
     j = inapa(i)
     if(vname(i).eq.'help    ' .or. vname(i).eq.'h       ' .or. vname(i).eq.'?       ' ) then
        ! Helptext   #################################
        write(*,*)' help ,h,? for helptext'
        write(*,*)' gplot,gp to plot selected files in connected xmgrace'
        write(*,*)'          gp plots into same grace until connection is closed'
        write(*,*)'          without open connection a new grace process is created'
        write(*,*)'          common parameter were extracted and a legend is shown'
        write(*,*)'          in grace: C-t write text; C-d delete text; '
        write(*,*)'                    C-M moves objekt; C-L moves legend '
        write(*,*)'Options:  switches between on/off '
        write(*,*)'         fits,fi         hide fits (but transfers data) '
        write(*,*)'         fitpara,fp      write last fit parameters'
        write(*,*)'         error,er        no error bars (but transfered) '
        write(*,*)'         legend,le       hide legend; '
        write(*,*)'                         parameters were appended, they build up shorter legend'
        write(*,*)'         noplot,np       no plot; to set switches'
        write(*,*)'         jointpar,jp     write joint parameters'
        write(*,*)'         close,cl        cloce connection to grace after plotting'
        write(*,*)'         clearplt,clear  deletes all previous data from plot '
        write(*,*)'                          and plots into a clean graph'
        write(*,*)'         autoscal,as      autoscale switched off/on'
        write(*,*)'Options:  '
        write(*,*)'         setting         shows current settings '
        write(*,*)'         closenow,cn     closes connection without plotting of new data'
        write(*,*)'         kill,ki         kills given sets in grace '
        write(*,*)'                   k,l deltes sets k l;  n -m deletes from n to m'
        write(*,*)'         comand,cm       rest of input line as direct command to grace,'
        write(*,*)'                           see xmgrace help or grace files as example'
        write(*,*)'         graph,gr        plot data in graph i, new graph i inside grace is build'
        write(*,*)'                                arrange them in EDIT/ARRANGE Graphs'
        write(*,*)'         save,sv         save graceplot with following name (.agr is added)'
        write(*,*)'                         name "on" saves each new plot to lastgrace.agr'
        write(*,*)'         setoffset,so    offset for numbering new data in grace, '
        write(*,*)'                            data with smaller numbers are preserved'
        write(*,*)'         cleanio         deletes .iofile if exists'
        write(*,*)'         start          the rest of the line will be sent directly to grace'
        write(*,*)'                        used to customise the startupbehaviour of xmgrace'
        write(*,*)'                    e.g. title "NSE data";subtitle "\xt\N is fouriertime";xaxis label "\xt\N / ns"'
        write(*,*)'                          its a semicolon separated list of xmgrace commands '
        write(*,*)'         fitcol,fc     fits coloured in i   i=1 for black, 0 same color as data'
        !############################################################    
        return
     endif

     if(vname(i).eq.'status  ' .or. vname(i).eq.'sta     ' ) then
        write(*,*)'Behaviour        Setting         change with'
        write(*,*)'hide fits        ',fitplo        ,'         fits     '
        write(*,*)'plot errors      ',errplo        ,'         error    '
        write(*,*)'plot legend      ',legend        ,'         legend   '
        write(*,*)'fit color        ',fb            ,'    fitcol   '
        write(*,*)'closeGrace       ',closeGrace    ,'         close    '
        write(*,*)'grace save       ',gsave         ,'         save     '
        write(*,*)'plot joint para  ',comon         ,'         jointpar '
        write(*,*)'write fit para   ',writefit      ,'         fitpara  '
        write(*,*)'autoscale        ',noautoscale   ,'         autoscal '
        write(*,*)'current graph    ',g0s           ,'      graph    '
        write(*,*)'clearplo         ',clearplo      ,'         clearplo '
        write(*,*)'startupsetting   ','"',trim(startupsetting),'"','         start '
        write(*,*)'connected Grace Pipe  ',GraceIsOpenF()
        return
     endif
     if(vname(i).eq.'fits    ' .or. vname(i).eq.'fi      ' ) then
        fitplo = .not. fitplo
        write(*,*)'fitplo set to->',fitplo
     endif
     if(vname(i).eq.'clearplo' .or. vname(i).eq.'clear   ' ) then
        clearplo = .not.clearplo
        write(*,*)'clearplot set to->',clearplo
     endif
     if(vname(i).eq.'fitpara ' .or. vname(i).eq.'fp      ' ) then
        writefit = .not.writefit
        write(*,*)'writefit',writefit
     endif
     if(vname(i).eq.'error   ' .or. vname(i).eq.'er      ' ) then
        errplo = .not.errplo
        write(*,*)'errplot set to->',errplo
     endif
     if(vname(i).eq.'legend  ' .or. vname(i).eq.'le      ' ) then
        legend = .not.legend
        write(*,*)'legend set to->',legend
     endif
     if(vname(i).eq.'close   ' .or. vname(i).eq.'cl      ' ) then
        closeGrace = .not.closeGrace
     endif
     if(vname(i).eq.'autoscal' .or. vname(i).eq.'as      ' ) then
        noautoscale = .not.noautoscale
        write(*,*)'noautoscale set to->',noautoscale
     endif
     if(vname(i).eq.'jointpar' .or. vname(i).eq.'jp      ' ) then
        comon = .not.comon
        write(*,*)'jointpar set to->',comon
     endif
     if(vname(i).eq.'noplot  ' .or. vname(i).eq.'np      ' ) then
        write(*,*)'no plot'
        return
     endif

     if(vname(i).eq.'kill    ' .or. vname(i).eq.'ki      ' ) then

        if ( inpar(i) .gt. 0 ) then
           do l=1,inpar(i)
              if(l.gt.minc) exit
              istart=0
              if (NINT(rpar(j)).lt.0 .and. NINT(-rpar(j)).gt.istart) then
                 iend=NINT(-rpar(j))
                 do ll=istart,iend,1
                    write(buf,"(I2)")ll
                    CALL MyGraceCommandf('kill g0.s'//TRIM(adjustl(buf)))
                 enddo
              else
                 istart=NINT(rpar(j))
                 write(buf,"(I2)")istart
                 CALL MyGraceCommandf('kill g0.s'//TRIM(adjustl(buf)))
              endif
              j = j + 1   ! j is adress parameter
           enddo
           CALL MyGraceCommandf ('redraw')
           CALL MyGraceCommandf ('updateall')
        endif
        return
     endif

     if(vname(i).eq.'graph   ' .or. vname(i).eq.'gr      ' ) then
        if ( inpar(i) .gt. 0 ) then
           ! gj.si    fuer grace hier wird j geschrieben i kommt spaeter
           write(buf,"(I2)") NINT(100*rpar(j)/100)
           g0s='g'//TRIM(ADJUSTL(buf))//'.s'
           CALL MyGraceCommandf ('with g'//TRIM(ADJUSTL(buf)))
           CALL MyGraceCommandf ('g'//TRIM(ADJUSTL(buf))//' on')
        else 
           g0s='g0.s'
           CALL MyGraceCommandf ('with g0')
           CALL MyGraceCommandf ('g0 on')
        endif
        CALL MyGraceCommandf ('updateall')
        write(*,*)'to arrange use Menue EDIT/ARRANGE Graphs'
     endif
     if(vname(i).eq.'fitcol  ' .or. vname(i).eq.'fc      ' ) then
        if ( inpar(i) .gt. 0 ) then
           fb=int(rpar(j))
        else
           fb=0
        endif
     endif

     if(vname(i).eq.'setoffs ' .or. vname(i).eq.'so      ' ) then
        sets = rpar(j)
     endif
     if(vname(i).eq.'closenow '.or.vname(i).eq.'cn      ') then
        call system('rm -f .iofile.tmp*')
        if ( GraceIsOpenf ().NE. 0 ) then
           CALL GraceClosePipef ()
           return
        endif
     endif
  enddo

  write(*,*)'  gplot help,h,? for helptext' 
  ! ----- ende parameter retrieving from stack -----
  !       is there something to plot?
  nkurv = nsel       ! number of selected sets , length of table of selected data
  if (nkurv.eq.0) then
     write(6,*)'there is no curve selected => no plot, even no close !'
     return
  endif


  !       equal params list #####################################
  !       build a list with for all sets equall parmas (parameters) (extra comment for fits
  !       fits)   fits have same params as data with added fit parameters
  commentar=''
  inallend=0
  ircf1=0
  commentarfit=''
  if (nkurv.gt.1 .or. ifits(1).ne.0) then
     ircu1 = isels(1)
     do l=1,nopar(ircu1)  ! initialisierung des Kommentars einfach der erste
        napar_all(l)=napar(l,ircu1)
        params_all(l)=params(l,ircu1)
        commentar=coment(1)
     enddo
     do i=2,nkurv
        ircu = isels(i)  ! adress of selected set i
        if (commentar.ne.coment(ircu)) then 
           commentar='unequal#comment_found'  
        endif
        do l=1,nopar(ircu)   ! loop for all params of set 
           do ij=1,nopar(ircu1) ! loop for common params
              if (napar_all(ij).eq.napar(l,ircu)) then     ! equal names
                 if (params_all(ij).ne.params(l,ircu)) then ! unequal values
                    napar_all(ij)='' ! all which are different in a parameterset where marked by a space
                 endif
              endif
           enddo
        enddo
     enddo
  endif
  !shrink list
  if (ircf1 .gt. 0) then
    inallend=nopar(ircu1) !needed later as length of parameterlist
    ij=1
    do while (ij<=inallend)  ! sort
        if (napar_all(ij).eq.'') then
            do while (napar_all(inallend).eq.'' .and. inallend .gt. 1)
            inallend=inallend-1
            enddo
            napar_all(ij)=napar_all(inallend)
            params_all(ij)=params_all(inallend)
            inallend=inallend-1
        endif
        ij=ij+1
    enddo  !erstellen der napar_all fertig
  endif
  ! now find all sets with fit  = name starts with fit
  ! first look in data then additional in corresponding fits
  inallendf=0
  nofitfoundflag=1 !flags no fit found ;should be a space (different from empty as indicator of different parameter values)
  ircf1 = 0
  if (nkurv.gt.1 .or. ifits(1).ne.0) then ! now only for fits
     do i=1,nkurv
        ircf = isels(i) ! adress of selected set i
        if (name(ircf)(1:3).eq.'fit') then
           if (nofitfoundflag.eq.1) then  ! first fit  as init 
              ircf1 = isels(i)
              do l=1,nopar(ircf1)  ! initialisierung des Kommentars einfach der erste
                 naparf_all(l)=napar(l,ircf1)
                 paramsf_all(l)=params(l,ircf1)
                 commentarfit=coment(1)
              enddo
              nofitfoundflag=0
           else   
              do l=1,nopar(ircf)   ! loop for all params of set 
                 do ij=1,nopar(ircf1) ! loop for common params
                    if (naparf_all(ij).eq.napar(l,ircf)) then     ! equal names
                       if (paramsf_all(ij).ne.params(l,ircf)) then ! unequal values
                          naparf_all(ij)='' ! all which are different in a parameterset where marked by a space
                       endif
                    endif
                 enddo
              enddo
           endif
        endif
        !corresponding fits  
        if (ifits(i) .gt. 0) then
           ircf = ifits(i)       ! adress of fit corresponding to selected set i
           if (name(ircf)(1:3).eq.'fit') then ! sollte immer erfüllt sein aber gut...
              if (nofitfoundflag.eq.1) then  ! first fit  as init
                 ircf1 = ifits(i)
                 do l=1,nopar(ircf1)  ! initialisierung des Kommentars einfach der erste
                    naparf_all(l)=napar(l,ircf1)
                    paramsf_all(l)=params(l,ircf1)
                    commentarfit=coment(1)
                 enddo
                 nofitfoundflag=0
              else   
                 do l=1,nopar(ircf)   ! loop for all params of set 
                    do ij=1,nopar(ircf1) ! loop for common params
                       if (naparf_all(ij).eq.napar(l,ircf)) then     ! equal names
                          if (paramsf_all(ij).ne.params(l,ircf)) then ! unequal values
                             naparf_all(ij)='' ! all which are different in a parameterset where marked by a space
                          endif
                       endif
                    enddo
                 enddo
              endif
           endif
        endif
     enddo !i=1,nkurv
  endif     ! nkurv.gt.1

  !			ende    equal params list ##############                ?????
  !shrink list by filling equal elsments with unequal up from bach to front ij
  if (ircf1 .gt. 0) then
     inallendf=nopar(ircf1)    ! needed later as length of parameterlist fits
     ij=1
     do while (ij<inallendf)  ! sort
        if (naparf_all(ij).eq.'') then
           do while (naparf_all(inallendf).eq.'')
              inallendf=inallendf-1
           enddo
           naparf_all(ij)=naparf_all(inallendf)
           paramsf_all(ij)=paramsf_all(inallendf)
           inallendf=inallendf-1
        endif
        ij=ij+1
     enddo
  endif
  !			
  !          now we have alist of equal parameters from data also present in fits in 
  !          and a list of equal fit parameters only present in fits
  !#############################################################################
  !      Start Grace with a buffer size of 2048 and open the pipe
  CALL MyGraceCommandf('updateall') ! only to see if it is realy open
  if ( GraceIsOpenf() .eq. 0) then
     write(*,*) 'open grace'
     write(*,*) '     Ctrl - d   to delete text/objects in xmgrace '
     write(*,*) '     Ctrl - m   to move   text/objects'
     write(*,*) '     Ctrl - l   to move   legend    '
     IF (GraceOpenf(2048) .EQ. -1) THEN
        WRITE (*,*) 'Can not run grace.'
        CALL EXIT (1)
     ENDIF
     g0s='g0.s'  ! to start in first graph
     CALL MyGraceCommandf (trim(startupsetting(:1024)))
     write(*,*)    trim(startupsetting(:1024))
     CALL GraceRegistErerrorFunctionf (MyError)
  else
     if (clearplo .eqv. .true.) then ! to start always with a clean plot
        do ll=0,99,1
           write(buf,"(I2)")ll
           CALL MyGraceCommandf('kill g0.s'//TRIM(adjustl(buf)))
        enddo
     endif
  endif
  CALL MyGraceCommandf ("default symbol size 0.500000")

  if ( noautoscale .eqv. .false. ) then
     CALL MyGraceCommandf('AUTOSCALE ONREAD XYAXES')
  else
     CALL MyGraceCommandf('AUTOSCALE ONREAD NONE')
  endif
  !     	Entering plot section
  ! 	---- plot datarecords ----#####################
  !
  !CALL MyGraceCommandf (TRIM(g0s)//'0 HIDDEN TRUE')
  do i=1,nkurv
     sets=sets+1
     write(ii,"(I8)") sets-1
     g0si=TRIM(g0s)//TRIM(ADJUSTL(ii))
     ircu = isels(i)                             ! adress of selected set i
     npic = nwert(ircu)                    ! number of points in ircu
     CALL MyGraceCommandf (g0si//' on')
     CALL MyGraceCommandf ('kill '//g0si)        ! Falls noch Daten drinstanden

     if  (npic < 50) then ! !write data to grace 
        CALL MyGraceCommandf (g0si//' type xydy')
        do  j=1,npic       
           write(jj,"(I8)") j-1
           WRITE (buf, "(E12.4,' , ',E12.4)") xwerte(j,ircu),ywerte(j,ircu)
           CALL MyGraceCommandf (g0si//' point '//buf)
           WRITE (buf, "(E12.4)") yerror(j,ircu)
           CALL MyGraceCommandf (TRIM(g0si)//'.y1['//TRIM(ADJUSTL(jj))//']='//buf)
        enddo
     else    ! writing complete set and read from file is faster for big sets
        write(buf,"(I4.4)") i-1
        open(10,file='.iofile.tmp'//trim(g0si))
        do  j=1,npic
           write(10,"(E12.4,' ',E12.4,' ',E12.4)") xwerte(j,ircu),ywerte(j,ircu),yerror(j,ircu)
        enddo
        close(10)
        CALL MyGraceCommandf (g0si//' point 0,0')
        CALL MyGraceCommandf ('READ xydy "'//'.iofile.tmp'//trim(g0si)//'"')
        CALL MyGraceCommandf ('MOVE S_ TO '//TRIM(g0si)//' ')
        CALL MyGraceCommandf ('s_ HIDDEN FALSE')
        ! datreat is faster than xmgrace so dont delete the files
        ! clean files with q (exit grace))

     endif! !write data to grace 				

     write(buf,"(I8)") sets-10*(sets/10)+1
     if (name(ircu)(1:3).eq.'fit') then
        CALL MyGraceCommandf (g0si//' symbol 0')
        CALL MyGraceCommandf (g0si//' line type 1')
     else 
        CALL MyGraceCommandf (g0si//' symbol '//buf)          
        CALL MyGraceCommandf (g0si//' line type 0')
     endif
     write(buf,"(I8)") sets-12*(sets/12)+1
     CALL MyGraceCommandf (g0si//' line color '//buf)
     CALL MyGraceCommandf (g0si//' symbol color '//buf)
     CALL MyGraceCommandf (g0si//' symbol fill color '//buf)
     CALL MyGraceCommandf (g0si//' errorbar color '//buf)
     write(buf,"(I1)") sets-2*(sets/2)
     CALL MyGraceCommandf (g0si//' symbol fill pattern '//buf)


     CALL MyGraceCommandf (g0si//' errorbar size 0.3')
     if(errplo) then
        CALL MyGraceCommandf (g0si//' errorbar on')
     else
        CALL MyGraceCommandf (g0si//' errorbar off')
     endif
     !          writing parameters as legend and always as comment
     buf=""
     if (commentar.eq.'unequal#comment_found') then
        write(buf,"('comment ',A,'\n')") TRIM(coment(ircu))
     endif
     buf2=""
     do l=1,nopar(ircu)  ! write buffer of unequal params
        innapar_all=0
        do ij=1,inallend
           if(napar_all(ij).eq.napar(l,ircu)) then
              innapar_all=1
              exit
           endif
        enddo
        do ij=1,inallendf
           if(naparf_all(ij).eq.napar(l,ircu)) then
              innapar_all=1
              exit
           endif
        enddo
        if (innapar_all.eq.0) then
           write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf(:1800)),TRIM(napar(l,ircu)),params(l,ircu)
        endif
        write(buf2,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf2(:1800)),TRIM(napar(l,ircu)),params(l,ircu)
     enddo
     CALL MyGraceCommandf(g0si//'comment '//quot//TRIM(buf2(:1900))//quot)
     if(legend) then
        CALL MyGraceCommandf(g0si//'legend '//quot//TRIM(buf(:1900))//quot)
     else
        buf=""
        do l=1,nopar(ircu)
           if(found(napar(l,ircu)//' ')) then
              write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf),TRIM(napar(l,ircu)),params(l,ircu)
           endif
        enddo
        CALL MyGraceCommandf(g0si//'legend '//quot//TRIM(buf(:1900))//quot)
     endif
     !           ende writing parameters as legend and always as comment


     !###################add fits to plot####################################
     !        if a new corresponding fit exists plot it
     if(fitplo) then  !c  ----> plot the fitted data automatically immediatly after fit
        ircf = ifits(i)		               ! adress of corresponding fit (to ircu)
        if(ircf.ne.0) then         ! =0 is not fitted !=0 it is fitted
           sets=sets+1
           write(ii,"(I8)") sets-1
           g0si=TRIM(g0s)//TRIM(ADJUSTL(ii))
           CALL MyGraceCommandf (g0si//' on')
           CALL MyGraceCommandf ('kill '//g0si)        ! Falls noch Daten drinstanden

           npicf = nwert(ircf)   ! number of fittpointsc
           if  (npicf < 30) then ! !write data to grace 
              CALL MyGraceCommandf (g0si//' type xy')
              do  j=1,npicf       
                 write(jj,"(I8)") j-1
                 WRITE (buf, "(E12.4,' , ',E12.4)") xwerte(j,ircf),ywerte(j,ircf)
                 CALL MyGraceCommandf (TRIM(g0si)//' point '//buf)
              enddo
           else    ! writing complete set and read from file is faster for big sets
              write(buf,"(I4.4)") i-1
              open(10,file='.iofile.tmpf'//trim(buf))
              do  j=1,npicf
                 write(10,"(E12.4,' ',E12.4)") xwerte(j,ircf),ywerte(j,ircf)
              enddo
              close(10)
              CALL MyGraceCommandf (g0si//' point 0,0')
              CALL MyGraceCommandf ('READ xy "./.iofile.tmpf'//trim(buf)//'"')
              CALL MyGraceCommandf ('MOVE S_ TO '//TRIM(g0si)//' ')
              CALL MyGraceCommandf ('S_ HIDDEN FALSE')
              ! see above comment datreat to fast for grace
           endif! !write data to grace 

           CALL MyGraceCommandf (g0si//' symbol 0')
           CALL MyGraceCommandf (g0si//' line type 1')
           if (fb.eq.0) then
              write(buf,"(I8)") (sets-1)-12*((sets-1)/12)+1
              CALL MyGraceCommandf (g0si//' line color '//buf)
              CALL MyGraceCommandf (g0si//' symbol color '//buf)            
              CALL MyGraceCommandf (g0si//' symbol fill color '//buf)
           else
              write(buf,"(I8)") int(fb)
              CALL MyGraceCommandf (g0si//' line color '//buf)
              CALL MyGraceCommandf (g0si//' symbol color '//buf)            
              CALL MyGraceCommandf (g0si//' symbol fill color '//buf)
           endif
           CALL MyGraceCommandf (g0si//' errorbar off')

           !				      writing parameters as legend and always as comment
           buf=""
           buf2=""
           do l=1,nopar(ircf)  ! write buffer of unequal params
              innapar_all=0
              do ij=1,inallend
                 if(napar_all(ij).eq.napar(l,ircf)) then
                    innapar_all=1
                    exit
                 endif
              enddo
              do ij=1,inallendf
                 if(naparf_all(ij).eq.napar(l,ircf)) then
                    innapar_all=1
                    exit
                 endif
              enddo
              if (innapar_all.eq.0) then
                 write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf),TRIM(napar(l,ircf)),params(l,ircf)
              endif
              write(buf2,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf2),TRIM(napar(l,ircf)),params(l,ircf)
           enddo
           CALL MyGraceCommandf(g0si//'comment '//quot//TRIM(buf2(:1900))//quot)
           if(legend) then
              CALL MyGraceCommandf(g0si//'legend '//quot//TRIM(buf(:1900))//quot)
           else
              buf=""
              do l=1,nopar(ircf)
                 if(found(napar(l,ircf)//' ')) then
                    write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf),TRIM(napar(l,ircf)),params(l,ircf)
                 endif
              enddo
              CALL MyGraceCommandf(g0si//'legend '//quot//TRIM(buf(:1900))//quot)
           endif
           !                 ende writing parameters as legend and always as comment
        endif
     endif
     !###################add fit ende####################################
  enddo ! ende plotten der Daten

  ! write all equal params

  buf=''
  do ij=1,inallend
     write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf),TRIM(napar_all(ij)),params_all(ij)
  enddo
  if (LEN_TRIM(buf).gt.0 .and. comon ) then
     CALL MyGraceCommandf ('with string')
     CALL MyGraceCommandf ('string on')
     CALL MyGraceCommandf ('string loctype view ')
     CALL MyGraceCommandf ('string 0.30, 0.4 ')
     CALL MyGraceCommandf ('string char size 0.600000')
     CALL MyGraceCommandf ('string def '//quot//'common params\n'//TRIM(buf(:1900))//quot )
  endif
  ! write all equal fit params
  buf=''
  do ij=1,inallendf
     iflag=0
     do ik=1,inallend
        if (napar_all(ik).eq.naparf_all(ij)) then
           iflag=1 ! kommt schon in params vor also nicht schreiben
           exit
        endif
     enddo
     if (iflag.eq.0) then              
        write(buf,"(A,' ',A,'=',E11.5,'\n')") TRIM(buf(:1800)),TRIM(naparf_all(ij)),paramsf_all(ij)
     endif
  enddo
  if (LEN_TRIM(buf).gt.0 .and. comon ) then
     CALL MyGraceCommandf ('with string')
     CALL MyGraceCommandf ('string on')
     CALL MyGraceCommandf ('string loctype view ')
     CALL MyGraceCommandf ('string 0.50, 0.4 ')
     CALL MyGraceCommandf ('string char size 0.600000')
     CALL MyGraceCommandf ('string def '//quot//'common fitparam:\n'//TRIM(buf(:1900))//quot )
  endif

  ! #############
  ! ---- plot theory parameters ----
  buf=""
  if(ntheos.ne.0 .and. writefit.eqv..true. ) then
     do it = 1,ntheos
        ith = nthtab(it)
        if(multflg(it).eq.1) then
           write(buf,"(8htheory* ,a8,'\n')") thenam(ith)
        else
           write(buf,"(8htheory+ ,a8,'\n')") thenam(ith)
        endif
        npar = nthpar(ith)
        if(npar.ne.0) then
           do ip = 1,npar
              write(buf,"(A,a8,1h=,1e12.4,2h+-,e9.2,e8.1,'\n')") TRIM(buf),thparn(ip,ith),thparx(ip,it),&
                   therro(ip,it),thpsca(ip,it)
           enddo
        endif
     enddo
     CALL MyGraceCommandf ('with string')
     CALL MyGraceCommandf ('string on')
     CALL MyGraceCommandf ('string loctype view ')
     CALL MyGraceCommandf ('string 0.80, 0.13 ')
     CALL MyGraceCommandf ('string def '//quot//buf(:1900)//quot)
  endif
  ! ---- plot the selected fit-kurves ----#################################
  !      This is not used and not functional
  nfkurv = nfsel                       ! length of selected fits list isfits
  if(fitplo) then                       ! should fits be plotted if yes do.......
     do i=1,nfkurv                      ! for all sel fits
        sets=sets+1
        g0si=TRIM(g0s)//TRIM(ADJUSTL(ii))
        CALL MyGraceCommandf (g0si//' on')
        CALL MyGraceCommandf ('kill '//g0si)        ! Falls noch Daten drinstanden
        CALL MyGraceCommandf (g0si//' symbol 0')
        CALL MyGraceCommandf (g0si//' line type 1')
        write(buf,"(I8)") sets-12*(sets/12)+1
        CALL MyGraceCommandf (g0si//' line color '//buf)
        !CALL MyGraceCommandf (g0si//' symbol color '//buf)
        CALL MyGraceCommandf (g0si//' type xy')
        CALL MyGraceCommandf (g0si//' errorbar off')
        irfcu = isfits(i)                   ! adress of fitted spectrum
        npicf = nwert(irfcu)         !  number of valid points on buffer
        do j=1,npicf
           write(jj,"(I8)") j
           WRITE (buf, "(A,' point ',E12.4,' , ',E12.4)") g0si,xwerte(j,irfcu),ywerte(j,irfcu)
           CALL MyGraceCommandf (buf)
        enddo
     enddo
  endif
  !   ende plot selected fit curves #########################################
  !

  !c ---- ende   plot datarecords ----#####################
  !       close pipe  and leave Grace alone			

  CALL MyGraceCommandf ('redraw')
  CALL MyGraceCommandf ('updateall')
  if (gsave.eqv..true.) then
     if (trim(savename).eq.trim('on')) then
        CALL MyGraceCommandf('saveall '//quot//'lastgrace.agr'//quot)
     else
        CALL MyGraceCommandf('saveall '//quot//TRIM(savename)//'.agr'//quot)
        gsave=.false.
     endif
  endif
  if ( closeGrace.eqv..true. )  CALL GraceClosePipef ()

  return
end subroutine gplot



!     for Grace error
SUBROUTINE MyError (str)
  IMPLICIT NONE
  CHARACTER*(*) str
  WRITE (0, '(''library message : "'', a, ''"'')') str
  RETURN
END SUBROUTINE MyError


SUBROUTINE MyGraceCommandf(str)
  common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)
  EXTERNAL GraceCommandf
  CHARACTER*(*) str
  if(iout.gt.0) WRITE(*,*) '>>>'//trim(str)
! CALL GraceCommandf(trim(str)//CHAR(0))
  CALL GraceCommandf(trim(str))
END SUBROUTINE MyGraceCommandf

