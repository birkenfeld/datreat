        function getyval(xvalue, specnum) 
!       ===================================================             
!                                                                       
! -------> returns the y value of spectrum number specnum  <--------                  
!          which corresponds to xvalue                                           
!          used in backscattering background determination  
!           	                                                            
!       implicit none                 
      use cincom
      use cincoc
      use xoutxx
      use xroxxx
      use constants
	   real*4            :: getyval
	   real*4            :: xvalue
	   integer           :: specnum, i
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

	   
	   
