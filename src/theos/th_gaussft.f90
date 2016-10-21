      FUNCTION thgaussft (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!      gaussian fit + ft (point fits or area fit)
!      change doFT parameter for analytic FT
!      Y=intensit*exp((x-center)/width)->Y=intensit*exp(-1/4*x*width)*width*exp(-center*x)*sqrt(Pi)                     
!      implicit none                                           
  
      double precision, parameter :: Pi=3.141592653589793D0   
      CHARACTER(8) thnam, parnam (29) 
      DIMENSION pa (29)
      integer     , intent(inout) :: mbuf
      integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
      real        , intent(inout) :: params(mbuf)             ! value des parameters n
      real*8    arg, erraccu, width, ftreal, ftimag, d

!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'gaussft' 
         nparx = 29
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           thgaussft = 0.0d0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1)  = 'inten1' 
         parnam (2)  = 'width1' 
         parnam (3)  = 'cente1'
         parnam (4)  = 'inten2' 
         parnam (5)  = 'width2' 
         parnam (6)  = 'cente2'
         parnam (7)  = 'inten3' 
         parnam (8)  = 'width3' 
         parnam (9)  = 'cente3'
         parnam (10) = 'inten4' 
         parnam (11) = 'width4' 
         parnam (12) = 'cente4'
         parnam (13) = 'inten5' 
         parnam (14) = 'width5' 
         parnam (15) = 'cente5'
         parnam (16) = 'inten6' 
         parnam (17) = 'width6' 
         parnam (18) = 'cente6'
         parnam (19) = 'inten7' 
         parnam (20) = 'width7' 
         parnam (21) = 'cente7'
         parnam (22) = 'inten8' 
         parnam (23) = 'width8' 
         parnam (24) = 'cente8'
         parnam (25) = 'wmin'            ! minimal width for gaussians
         parnam (26) = 'level'
         parnam (27) = 'slope'
         parnam (28) = 'doFT'            ! do fourier transform? (1 or 0)
         parnam (29) = 'area'            ! do histogram area fitting ala M.M. (energy bin width or 0)

!                                                                       
         thgaussft = 0.0d0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----
       thgaussft = 0.0d0
       if (pa(28).eq.0.0d0) then                        ! do fit in w space

        if (pa(29).eq.0.0d0) then                        ! do point fit
         do i=0,7
          if (pa(i*3+1).ne.0.0d0) then
           width = abs(pa(i*3+2))
           if (width.le.pa(25)) then
            width = abs(pa(25))
           endif
           arg = ( (x-pa(i*3+3)) / width ) **2.0d0
           thgaussft = thgaussft + pa(i*3+1) * exp(-arg)
          endif
         enddo
         thgaussft = thgaussft + abs(pa(26)) + pa(27)*x
         return

        elseif (pa(29).ne.0.0d0) then                   ! do area fit of histogram data
         d = pa(29)/2.0d0                              
         do i=0,7
          if (pa(i*3+1).ne.0.0d0) then
           width = abs(pa(i*3+2))
           if (width.le.pa(25)) then
             width = abs(pa(25))
           endif
           arg = erf((x-pa(i*3+3)+d)/width) - erf((x-pa(i*3+3)-d)/width) 
           thgaussft=thgaussft + pa(i*3+1) * 0.5d0 * sqrt(Pi) * width / abs(pa(29)) * arg
          endif
         enddo
         thgaussft = thgaussft + abs(pa(26)) + pa(27)*x
         return
        endif


       elseif (pa(28).eq.1.0d0) then                   ! do FT of gaussians
        ftreal = 0.0d0
        ftimag = 0.0d0                             
        do i=0,7
         if (pa(i*3+1).ne.0.0d0) then
          width = abs(pa(i*3+2))
          if (width.le.pa(25)) then
           width = abs(pa(25))
          endif
          arg = (width * x ) **2.0d0
          thgaussftt = pa(i*3+1) * exp(-1.0d0/4.0d0*arg) * width * sqrt(Pi)
          ftreal = ftreal + thgaussftt * cos(-pa(i*3+3)*x)
          ftimag = ftimag + thgaussftt * sin(-pa(i*3+3)*x)         
         endif
        enddo
        thgaussft = sqrt(ftreal**2.0d0 + ftimag**2.0d0)
        return
       endif
      write(6,*) 'No method given!'
      return
      END FUNCTION thgaussft


!     example snippet for cutoff determination, fails because of diverging erf (better implementation needed)     
!!       do cutoff determination using erf
!           erfarg = cmplx((pa(30)-pa(i*3+3))/width, width*x)
!           call CERROR(erfarg, erfresult)
!           erfarg = cmplx((-pa(30)-pa(i*3+3))/width, width*x)
!           call CERROR(erfarg, erfresult2)
!           erfreal = real(erfresult-erfresult2)
!           erfimag = aimag(erfresult-erfresult2)
!           ftreal=ftreal+0.5d0*thgaussftt*(cos(-pa(i*3+3)*x)*erfreal-sin(-pa(i*3+3)*x)*erfimag)
!           ftimag=ftimag+0.5d0*thgaussftt*(sin(-pa(i*3+3)*x)*erfreal+cos(-pa(i*3+3)*x)*erfimag)
!          endif                               ! do normal ft




!       complex erf implementation
	SUBROUTINE CERROR(Z,CER)
!
!       ====================================================
!       Purpose: Compute error function erf(z) for a complex
!                argument (z=x+iy)
!       Input :  z   --- Complex argument
!       Output:  CER --- erf(z)
!       ====================================================
!
	IMPLICIT COMPLEX *16 (C,Z)
	DOUBLE PRECISION A0,PI
	A0=CDABS(Z)
	C0=CDEXP(-Z*Z)
	PI=3.141592653589793D0
	Z1=Z
	IF (REAL(Z).LT.0.0) THEN
	   Z1=-Z
	ENDIF
	IF (A0.LE.5.8D0) THEN    
	   CS=Z1
	   CR=Z1
	   DO 10 K=1,120
	      CR=CR*Z1*Z1/(K+0.5D0)
	      CS=CS+CR
	      IF (CDABS(CR/CS).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CER=2.0D0*C0*CS/DSQRT(PI)
	ELSE                              
	   CL=1.0D0/Z1              
	   CR=CL
	   DO 20 K=1,13
	      CR=-CR*(K-0.5D0)/(Z1*Z1)
	      CL=CL+CR
	      IF (CDABS(CR/CL).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CER=1.0D0-C0*CL/DSQRT(PI)
	ENDIF
	IF (REAL(Z).LT.0.0) THEN
	   CER=-CER
	ENDIF
	RETURN
	END
