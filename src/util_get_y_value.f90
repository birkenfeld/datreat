        function getyval(xvalue, specnum) 
!       ===================================================             
!                                                                       
! -------> returns the y value of spectrum number specnum  <--------                  
!          which corresponds to xvalue                                           
!          used in backscattering background determination  
!           	                                                            
!       implicit none                 
	   real*4            :: getyval
	   real*4            :: xvalue
	   integer           :: specnum, i
	   integer, parameter:: mwert=1024, mbuf=200, mpar=200
! --- mwert  = max. no. of x-y-values in one buffer                     
!     mbuf   = max. no. of different buffers                            
!     mpar   = max. no. of parameters associated with one buffer        
! ---  maximum scan length ....                                         
	   integer, parameter:: mth=40, mtpar=40,mtcal=40,mcoup=10
! ---  fit dimensions ---                                               
	   integer, parameter:: mfit=40,msmpl=4000
!  -- mfit = max no. of fitted parameters                               
!     msmpl= max no. of datapoints in fit                               
!                                                                       
! --- incom common-section ---                                          
       integer, parameter:: minc=40
	   
	   common/cincom/rpar(minc),inames,ipars                            &
     &  ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf   &
     &  ,lstpar, lstnam                                                 
       common/cincoc/comand,vname(minc),title,reslin,inline             &
     &  ,arglst(minc),pmlist(minc,2),rlbuf                 
       logical cray              
       common/xoutxx/iot,ioold,ibild1,ierrr,inka1, cray 
       common/xroxxx/  xyorig(3), rotvec(3) 
!                                                                       
!                                                                       
!
       character(len=80) name,xname,yname,napar,coment*80        
       common/cdata/ xwerte(mwert,mbuf),ywerte(mwert,mbuf),             &
     &        yerror(mwert,mbuf),                                       &
     &        xname(mbuf),yname(mbuf),name(mbuf),nwert(mbuf),           &
     &        numor(mbuf),nbuf,coment(mbuf),params(mpar,mbuf),          &
     &        napar(mpar,mbuf),nopar(mbuf)                              
! --- xwerte(i,j)   x-values on buffer j                                
!     ywerte(i,j)   y-values on buffer j                                
!     yerror(i,j)   error of y-values   (only supported by some fktn)   
!     xname(j)      name of x-values (x-axis) for buffer j              
!     yname(j)       "   "  y-   "    y-  "    "    "    "              
!     name(j)       short text-identifier for data on buffer j          
!     nwert(j)      no. of valid points on buffer j                     
!     nbuf          no. of filled buffers                               
!     numor(j)      numerical idenfication of data on buffer j          
!     coment(j)     one line of comment describing data on buffer j     
!     params(l,j)   set of parameters associated with data on buffer j  
!     napar(l,j)    names of these parameters                           
!     nopar(j)      no. of valid parameters                             


       
	   if (specnum.gt.nbuf) then
		   write(6,*)'Spectrum number too large...', nbuf
		   write(6,*)'nbuf = ', nbuf, ' < specnum =', specnum
		   getyval = 0.0
		   return
	   endif

       getyval = 1.0                                                                    
	   xsearch: do i=1,mwert
	       if (xvalue.le.xwerte(i,specnum)) then
               getyval = ywerte(i,specnum)
			   exit 
		   endif
		   getyval = 0.0	   
	   enddo xsearch
       return 
	   end                                                

	   
	   
