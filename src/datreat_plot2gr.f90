!*ds                                                                    
       subroutine splot (doplo) 
!      ================  scan-plotting                                  
!                                                                       
!                                                                       
       use cincom
       use cincoc
       use constants
       parameter(mth=40,mtpar=40,mtcal=40) 
!      parameter(mkurv=minc)                                            
!                ---------> max. no of curves to be plotted             
! ---  maximum scan length                                              
!                                                                       
       parameter(mwert=1024,mbuf=200,mpar=200) 
! ---- common containing all the scans ----                             
       character*80 name,xname,yname,napar,coment*80 
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),             &
     &        yerror(mwert,mbuf),                                       &
     &        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),           &
     &        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),          &
     &        napar(mpar,mbuf),nopar(mbuf)                              
!                                                                       
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls 
       common/fslist/isfits(mbuf),nfsel 
!                                                                       
       character*8 thenam,thparn 
       common/theory/thenam(mth),thparn(mtpar,mth),nthpar(mth),         &
     &  thparx(mtpar,mtcal),thpsca(mtpar,mtcal),nthtab(mtcal),ntheos    &
     & ,multflg(mtcal)                                                  
       common/therrc/therro(mtpar,mtcal) 
       character*8 thrapar 
       real*4      thramin, thramax 
       common /thparc/ thrapar(mth),thramin(mth),thramax(mth) 
!                                                                       
!                                                                       
       character*6 cnum 
!                                                                       
       character*1024 inlbuf 
      logical cray 
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray 
       common/xroxxx/  xyorig(3), rotvec(3) 
!                                                                       
!                                                                       
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20),yyee 
!                                                                       
       dimension x(mwert),y(mwert),yl(mwert),yh(mwert),irecv(minc),     &
     &          isymb(minc),irecn(minc),icolo(minc),ifrec(minc)         
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
                                                                        
!      --- doplo = false  means: set parameters only ---                
       character*80 option,xtext,ytext,tbuf 
       character*8  dirnam(4),opart(8) 
                                                                        
       character*12 tag, stunde 
       character*12 tx,sx 
!                                                                       
      data xmin/0./,xmax/1./,ymin/0./,ymax/1./,nkurv/0/ 
      data isymb/4,5,23,6,16,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21  &
     &          ,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41  &
     &          ,42,43/                                                 
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
      data opart/'x=1     ','y=1     ','i       ','a       ',           &
     &           'f=(3,1) ','m=2.0   ','u=1     ','        '/           
      data txsizt/.23/,xtshft/0./,ytshft/0./ 
!                                                                       
!                                                                       
!                                                                       
! ----- parameter retrieving from stack -----                           
      nkurv  = 0 
      nfkurv = 0 
      if(inames.eq.0) then 
!                     ----> assume that a list of number numors         
       if(ipars.gt.0) then 
          nkurv = 0 
          nsel  = 0 
          if(ipars.gt.minc) ipars = minc 
          do 2 i=1,ipars 
           irecv(i) = rpar(i) + 0.0001 
    2     continue 
          nkurv = ipars 
        endif 
       else 
!      -----> decode by names                                           
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
!   ---> txon, txoff are from old version and not documented any more <-
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
!                                                                       
          if(vname(i).eq.'symb    ') then 
            nsy = 0 
            do 49 l=1,inpar(i) 
             nsy   = nsy   + 1 
             if(nsy.gt.minc) goto 29 
             isymb(nsy) = rpar(j) + 0.0001 
             j = j + 1 
   49       continue 
          endif 
   29     continue 
!                                                                       
          if(vname(i).eq.'icolo   ') then 
            nco = 0 
            do 59 l=1,inpar(i) 
             nco   = nco   + 1 
             if(nco.gt.minc) goto 39 
             icolo(nco) = rpar(j) + 0.0001 
             j = j + 1 
   59       continue 
          endif 
   39     continue 
!                                                                       
          if(vname(i).eq.'sc      ') then 
            nsel = 0 
            nkurv= 0 
            do 5 l=1,inpar(i) 
             nkurv = nkurv + 1 
             if(nkurv.gt.minc) goto 31 
             irecv(nkurv) = rpar(j) * 1.000001 
             j = j + 1 
    5       continue 
          endif 
   31     continue 
!                                                                       
          if(vname(i).eq.'fsc     ') then 
            nfsel = 0 
            nfkurv= 0 
            do 641 l=1,inpar(i) 
             nfkurv = nfkurv + 1 
             if(nfkurv.gt.minc) goto 64 
             ifrec(nfkurv) = rpar(j) * 1.000001 
             j = j + 1 
  641       continue 
          endif 
   64     continue 
          if(nfsel.eq.0) then 
            if(nfkurv.gt.0) then 
               call fsrch(ifrec,nfkurv) 
            endif 
          endif 
!                                                                       
    3   continue 
      endif 
!                                                                       
      if (.not.doplo) then 
        return 
      endif 
!                                                                       
      if(nsel.eq.0) then 
        if(nkurv.gt.0) then 
          call search(irecv,nkurv) 
!       -----------------------------> select the items to be plotted   
        else 
          write(6,*)'no scan selected, no plot possible !' 
          return 
        endif 
      endif 
!                                                                       
      if (.not.doplo) then 
        return 
      endif 
! ---->                                             <-----------        
      nkurv = nsel 
! ----> ich denke, das muss so sein !?!?!?!  (m.s.) <-----------        
      if (nkurv.eq.0) then 
        write(6,*)'there is no curve selected => no plot !' 
        return 
      endif 
!                                                                       
!                                                                       
!      write(*,*)' entering plot section .....'                         
!                                                                       
! ----- initialize plot (if first) -----                                
       if(ibild.eq.0) then 
         call grstrt(35,8) 
       endif 
!                                                                       
       ibild = ibild + 1 
!                                                                       
! ---- set frame & scales ----                                          
!      write(*,*)' set frames & scales .....'                           
       frxx = frlux + framx 
       fryy = frluy + framy 
       call grsclc(frlux,frluy,frxx,fryy) 
       call grsclv(xmin,ymin,xmax,ymax) 
       call grfont(ifont) 
       call grnwpn(icol0) 
                                      !! aix                            
       call grchrc(txsize,0.,16 ) 
!      --------------------------> set size of signs                    
! ---- prepare axes ----                                                
!                                                                       
! ---- identify scantyp ----                                            
!                                                                       
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
       option = opart(1)//','//opart(2)//','//opart(3)//','//           &
     &                         opart(5)//','//opart(6)//','//opart(7)   
       call upcase(option,80) 
       write(6,*)option 
       lopt = 80 
!      write(*,*)' make axes .....'                                     
       if(.not.taxis) then 
         lxx = 0 
         lyy = 0 
       else 
         lxx = laenge(xtext,80,' ') 
         lyy = laenge(ytext,80,' ') 
       endif 
       if(paxis) call graxs(lopt,option,lxx,xtext,lyy,ytext) 
!                                                                       
!                                                                       
! ---- plot the selected fit-kurves ----                                
!                                                                       
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
   70   continue 
       endif 
!                                                                       
!                                                                       
! ---- plot datarecords ----                                            
!                                                                       
       nkurv = nsel 
       do 20 i=1,nkurv 
        ircu = isels(i) 
        ircf = ifits(i) 
        npic = nwert(ircu) 
        if(fitplo) then 
        if(ircf.ne.0) then 
!                     ----> plot the fitted data automatically          
!                           this is after a fit-command until you       
!                           select new curves                           
        npicf = nwert(ircf) 
        nnpi  = 0 
        do 20010 j=1,npicf 
          if(xwerte(j,ircf).ge.xmin.and.xwerte(j,ircf).le.xmax) then 
            nnpi = nnpi + 1 
            y(nnpi) = ywerte(j,ircf) 
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02 
            if(y(nnpi).gt.ymax) y(nnpi) = ymax+(ymax-ymin)*0.02 
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
!                                                                       
!                  ---> prepare data to plot                            
        nnpi = 0 
        do 2001 j=1,npic 
          if(xwerte(j,ircu).ge.xmin.and.xwerte(j,ircu).le.xmax) then 
            nnpi = nnpi + 1 
            y(nnpi) = ywerte(j,ircu) 
            x(nnpi) = xwerte(j,ircu) 
            if(y(nnpi).lt.ymin) y(nnpi) = ymin-(ymax-ymin)*0.02 
            if(y(nnpi).gt.ymax) y(nnpi) = ymax+(ymax-ymin)*0.02 
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
 2001   continue 
!                                                                       
!                                                                       
!                                                                       
! --- plot ---                                                          
!       if (numor(ircu).gt.0) then                                      
           icco=mod(icolo(i),7) + 1 
                                                                        
           call grnwpn(icco) 
!          ----------- plot a dataline -------                          
           if(isymb(i).eq.0) then 
             call grln(x,y,nnpi) 
           else 
             sysiz = sysize 
             call grchrc(sysiz,0.,16) 
             do 2003 ik=1,nnpi 
              call grjmps(x(ik),y(ik),isymb(i)) 
 2003        continue 
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
!       endif                                                           
!                                                                       
   20  continue 
!                                                                       
!                                                                       
                                                                        
! ---- textpart ----                                                    
       if(ptex) then 
!      write(*,*)' enter textplotting....'                              
       call grsclv(xmin,ymin,xmax,ymax) 
!                                                                       
! --- title ---                                                         
         xtx = xmin + 0.1 * (xmax-xmin) 
         ytx = ymax 
         ltext = 20./txsize 
         if(ltext.gt.74) ltext=74 
         call grtxt(xtx,ytx,ltext,title) 
! ---- set textwindow ----                                              
!                                                                       
         call grsclc(0.,0.,39.5,28.7) 
         call grsclv(0.,0.,39.5,28.7) 
!        txsizt = 0.65 * txsize                                         
                                      !! aix                            
         call grchrc(txsizt,0.,16 ) 
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
                                                                        
! ---- plot theory parameters ----                                      
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
             if(thrapar(it).ne.'        ') then 
               write(xtext,'(4hfor ,a8,2f12.6)')thrapar(it),thramin(it),thramax(it)                            
               call grtxt(xtx,ytx,36,xtext) 
               ytx = ytx - 1.7 * txsizt 
             endif 
             npar = nthpar(ith) 
             if(npar.ne.0) then 
               do 117 ip = 1,npar 
               write(xtext,'(a8,1h=,1e12.4,2h+-,e9.2,e8.1)')thparn(ip,ith),thparx(ip,it),therro(ip,it),thpsca(ip,it) 
                 call grtxt(xtx,ytx,33,xtext) 
                 ytx = ytx - 1.7*txsizt 
  117          continue 
            endif 
  115     continue 
         endif 
! ---- plotted items ----                                               
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
 1012      continue 
           else 
           do l=1,nopar(ircu) 
            if(found(napar(l,ircu)//' ')) then 
              write(xtext,'(a8,e14.5)')napar(l,ircu),                   &
     &                                 params(l,ircu)                   
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
                write(xtext,'(a8,e14.5)')napar(l,ircu),                 &
     &                                   params(l,ircu)                 
                call grtxt(xtx,ytx,22,xtext) 
                ytx = ytx - 1.7 * txsizt 
              endif 
             enddo 
             endif 
           endif 
  101    continue 
         call grnwpn(1) 
      endif 
!                                                                       
!     write(*,*)' ok ready call nextframe ....'                         
      call grnxtf 
!                                                                       
      return 
      END                                           
!                                                                       
!                                                                       
       subroutine out_gli 
!      ==================                                               
!                                                                       
       use cincom
       use cincoc
       use constants
       parameter(mth=40,mtpar=40,mtcal=40) 
       parameter(mwert=1024,mbuf=200,mpar=200) 
       character*1024 inlbuf 
      logical cray 
      common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray 
       common/xroxxx/  xyorig(3), rotvec(3) 
!                                                                       
!                                                                       
       character*80 name,xname,yname,napar,coment*80 
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),             &
     &        yerror(mwert,mbuf),                                       &
     &        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),           &
     &        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),          &
     &        napar(mpar,mbuf),nopar(mbuf)                              
!                                                                       
       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls 
       common/fslist/isfits(mbuf),nfsel 
!                                                                       
!                                                                       
       common/outlev/iout,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20),yyee 
!                                                                       
       character*12 infile 
       character*80 rline 
       logical*4    fileda 
       dimension sq(3),ihkl(3),xh(3) 
!                                                                       
! -- open the file -----------------------------------------------------
       if(inames.eq.0) then 
         ierrs = 1 
         write(6,*)'input(e1): filename is lacking!' 
         return 
       endif 
!                                                                       
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
         write(20,'(18(1x,E13.6))')xwerte(j,isels(1)),                  &
     &        (ywerte(j,isels(i)),i=1,nsel)                             
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
      END                                           
!*ds                                                                    
!*ds
       subroutine preplo 
!      ================= prepare plotting (plotmachine)                 
!                                                                       
! --- minc = incom stack depth                                          
!                                                                       
       use cincom
       use cincoc
       use constants
       character*1024 inlbuf 
      logical cray 
      common/xoutxx/iot,ioold,ibild,ierrr,inka, cray 
                                                                        
       logical ptex 
       logical paplo 
       logical doplo 
       logical fitplo 
       logical paxis 
       logical taxis 
!      --- doplo = false  means: set parameters only ---                
       character*80 option,xtext,ytext,tbuf,ttext,txtv,fmt,ttout 
       character*8  dirnam(4),opart(8) 
       character*8  codena, codefn, chrval 
       logical      found, folgt 
       dimension isymb(minc), icolo(minc) 
       real*8 xdbl(*), ydbl(*), xcode(*) 
       real*8 dble, getval, val 
!                                                                       
      data xmin/0./,xmax/4./,ymin/0./,ymax/2./,nkurv/0/ 
      data isymb/4,5,23,6,16,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21  &
     &          ,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41  &
     &          ,42,43/                                                 
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
      DATA OPART/'X=1     ','Y=1     ','I       ','A       ',           &
     &           'F=(3,1) ','M=2.0   ','U=1     ','        '/           
      data txsizt/.3/,xtshft/0./,ytshft/0./ 
      data ttscrx/39.5/, ttscry/28.7/ 
!                                                                       
!                                                                       
!                                                                       
! ----- parameter retrieving from stack -----                           
!      -----> decode by names                                           
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
!                                                                       
          if(vname(i).eq.'symb    ') then 
            nsy = 0 
            do 49 l=1,inpar(i) 
             nsy   = nsy   + 1 
             if(nsy.gt.minc) goto 29 
             isymb(nsy) = rpar(j) + 0.0001 
             j = j + 1 
   49       continue 
          endif 
   29     continue 
!                                                                       
          if(vname(i).eq.'icolo   ') then 
            nco = 0 
            do 59 l=1,inpar(i) 
             nco   = nco   + 1 
             if(nco.gt.minc) goto 39 
             icolo(nco) = rpar(j) + 0.0001 
             j = j + 1 
   59       continue 
          endif 
   39     continue 
!                                                                       
    3   continue 
!                                                                       
!                                                                       
! ----- initialize plot (if first) -----                                
       if(ibild.eq.0) then 
         call grstrt(35,8) 
         write(6,*)'incom: preplo: grstrt(35,8) <=======' 
       endif 
!                                                                       
       ibild = ibild + 1 
       isycnt= 0 
!                                                                       
! ---- set frame & scales ----                                          
!      write(*,*)' set frames & scales .....'                           
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
!      --------------------------> set size of signs                    
                                                                        
! --- title ---                                                         
         xtx = xmin + 0.1 * (xmax-xmin) 
         ytx = ymax 
         ltext = 20./txsize 
         if(ltext.gt.74) ltext=74 
         call grtxt(xtx,ytx,ltext,title) 
!                                                                       
         xtx = 28.3 + xtshft 
         ytx = 28.  + ytshft 
                                                                        
       return 
                                                                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
                                                                        
       entry setplt( codena, codefn, xcode, ncode ) 
!      --------------------------------------------                     
!      eingabe von plotparametern durch unterprogrammaufruf statt       
!      von der kommandoebene aus                                        
!      --------------------------------------------                     
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
!   ---> txon, txoff are from old version and not documented any more <-
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
!                                                                       
          if(codena.eq.'symb    ') then 
            nsy = 0 
            do 3049 l=1,ncode 
             nsy   = nsy   + 1 
             if(nsy.gt.minc) goto 3029 
             isymb(nsy) = xcode(l) + 0.0001 
 3049       continue 
          endif 
 3029     continue 
!                                                                       
          if(codena.eq.'icolo   ') then 
            nco = 0 
            do 3059 l=1,ncode 
             nco   = nco   + 1 
             if(nco.gt.minc) goto 3039 
             icolo(nco) = xcode(l) + 0.0001 
 3059       continue 
          endif 
 3039     continue 
!                                                                       
       return 
!      ======                                                           
                                                                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
                                                                        
       entry preaxe( xtext, ytext) 
!      ---------------------------                                      
                                                                        
! ---- prepare axes ----                                                
!                                                                       
!                                                                       
       option = opart(1)//','//opart(2)//','//opart(3)//','//           &
     &                         opart(5)//','//opart(6)//','//opart(7)   
       write(6,*)option 
       lopt = 80 
!      write(*,*)' make axes .....'                                     
       if(taxis) then 
          lxtxt = laenge(xtext,80,'$') 
          lytxt = laenge(ytext,80,'$') 
       else 
          lxtxt = 0 
          lytxt = 0 
       endif 
       if(paxis) call graxs(lopt,option,lxtxt,xtext,lytxt,ytext) 
!                                                                       
       return 
!      ======                                                           
                                                                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                        
       entry pldata( xdbl, ydbl, ndbl) 
!      ===============================                                  
! --- resetscaling                                                      
       call grsclc(frlux,frluy,frxx,fryy) 
       call grsclv(xmin,ymin,xmax,ymax) 
       call grfont(ifont) 
       call grchrc(txsize,0.,16.) 
!      --------------------------> set size of signs                    
! ---  plot daten                                                       
       isycnt = mod(isycnt,minc) + 1 
       isy    = isymb(isycnt) 
       ico    = icolo(isycnt) 
       call grnwpn(ico) 
       if(isy.eq.0) then 
         call grjmp( sngl(xdbl(1)), sngl(ydbl(1)) ) 
         do 200 i=2,ndbl 
            call grdrw( sngl(xdbl(i)), sngl(ydbl(i)) ) 
  200    continue 
       else 
         call grjmp( sngl(xdbl(1)), sngl(ydbl(1)) ) 
         do 201 i=2,ndbl 
            call grjmps( sngl(xdbl(i)), sngl(ydbl(i)), isy ) 
  201    continue 
       endif 
                                                                        
       return 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                       
       entry grcolo( icoln ) 
!      =====================                                            
         isycc  = mod(isycnt,minc) + 1 
         icolo(isycc ) = icoln 
       return 
!                                                                       
!                                                                       
       entry grwrit( ttext ) 
!      =====================                                            
!                                                                       
! ---- textpart ----                                                    
       if(.not.ptex) then 
         write(*,*)'no textplotting due to notext' 
         return 
       endif 
!                                                                       
! ---- set textwindow ----                                              
         call grsclc(0.,0.,ttscrx, ttscry) 
         call grsclv(0.,0.,ttscrx, ttscry) 
!                                                                       
         call grchrc(txsizt,0.,16.) 
         ltt = laenge(ttext,80,'$') 
         call grtxt(xtx,ytx,ltt,ttext) 
         ytx = ytx - 1.7 * txsizt 
                                                                        
         return 
!                                                                       
       entry grwrva( txtv , val, fmt ) 
!      ===============================                                  
!                                                                       
! ---- textpart ----                                                    
       if(.not.ptex) then 
         write(*,*)'no textplotting due to notext' 
         return 
       endif 
!                                                                       
! ---- set textwindow ----                                              
         call grsclc(0.,0.,ttscrx, ttscry) 
         call grsclv(0.,0.,ttscrx, ttscry) 
!                                                                       
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
!                                                                       
      END                                           
!                                                                       
!                                                                       
!*ds                                                                    
!*ed                                                                    