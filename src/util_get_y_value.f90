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
      use cdata
      use constants
	   real*4            :: getyval
	   real*4            :: xvalue
	   integer           :: specnum, i


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
