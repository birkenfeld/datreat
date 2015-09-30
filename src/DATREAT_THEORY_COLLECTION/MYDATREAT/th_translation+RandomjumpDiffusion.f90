       FUNCTION tradiff(x, pa, thnam, parnam, npar,ini)
!     ========================================
!     name of this function in quotes  to call in datreat
!
!     "transdiff"
!     # description
!	#start description
!	Translational Motion described by random jump diffusion model
!      S(q,w)= 1/pi)*Gamma(q)/(Gamma^2+w^2) with Gamma =D^2q^2/(1+Dq^2tau0)
!	Gamma  hwhm
!      D Translational Diffusion  A^2/ns
!     tau0   residence time ns
! 	#end description

!	enforce of  explicit variable definition
	implicit none
	real*4 tradiff                   ! function return value real*4 for datreat
!	some constants
	double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)
       character*8 thnam,parnam(20)       !      # theoriename and parameter names passed over (returned to datreat) in   parnam(i) and thname for ini=0
!      variables from datreat were all real*4
       real*4 pa(20)		!      # theorie parameter values in pa(i)    # definition real*4
       real*4 x                	!      x values data
       real*4 q_dat		!      # needed parameter from data
       integer npar,ini,ier    ! number of datreat parameters, initialization flag  0 or !=0  , errorflag
!       locale variables  can be real*8   ;   by local_var=datreat_var  datreat real*4 is cast to local real*8
       real*8    tau0    ! local parameters
	real*8 w,amp,Dt,l_jump,ga
	real*8 q
       real*8 Sqw                   ! result
       integer nparx                ! number of local parameters (excluding result)
!      Declaration of used functions
       real*8 gamma
       real*8 trandiff

!      thiadda actuall datreat dataset address in iadda
	integer iadda
       common/thiadd/ iadda

!
! ----- initialisation ----- #Init from datreat called with ini=0 to tell datreat all about this function
       if(ini.eq.0) then
!	Theoriename
         thnam = 'transdiff'
!	number of theorie parameter
	 nparx = 1
         if(npar.lt.nparx) then            ! check it and give a hint
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           tradiff = 0
           return
         endif
         npar = nparx
!        --------------> set the name of parameters
	 parnam(1)  = 'Amplitu'           !Amplitude of function
	 parnam(2)  = 'tranDiff'        ! Dt
	 parnam(3)   ='l_jump'		! mean square jump length
	 parnam(4)  = 'resTime'         ! residence time between jumps


!
!        default value for return
         trandiff = 0
         return
       endif
!  ------------end initialization ---------------ier=0


!      put  values  to some significant variablenames to make it easier additional this casts datreat real*4 to local real*8
       w      = x      !  x values
       amp      = pa(1)  !
       Dt          = pa(2)
       l_jump  = pa(2)

!      to get the value of parameter eg 'q' from dataset
! 	returnvalue ier=0 if succeded, value from 'q' in real*4 from datreat  (min a space at end as end of string, max 8 char)
       call        parget('q ',q_dat,iadda,ier)
       q   = q_dat   !  if necessary  here cast to definition of of q eg real*8
!     default value or evaluation dependent on ier
       if(ier.ne.0) then
       		q   = 1.0
        endif
! ---- calculate theory here -----
	Dt=trandiff(l_jump,tau0)
	ga=gamma(Dt,q,tau0)
	Sqw=1/pi *  (ga/(ga**2+w**2))


	! pass the result here at the end of calculation (cast to real*4)
	tradiff     =sngl( amp * Sqw )
!-------end of calculation value returned in tradiff                -------------

!!      to pass additional values from here to the dataset parameters  (sngl()  casts real*8 to real*4 for datreat)
!       call        parset('l       ',sngl(l),iadda,ier)
!
       return
       end
!##########end######################################################################end
!subroutine and functions
!      ==================================================
       real*8 function gamma(D,q,tau)

       implicit none
       real*8 D,q,tau
	gamma=D*q**2/(1+D*q**2*tau)
       return
       end
!      ==================================================
       real*8 function tradiff(l,tau)

       implicit none
       real*8 l,tau
		tradiff=l**2/6/tau
       return
       end

