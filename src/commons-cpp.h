#define __MINC 40
#define __MDEPTH 20
#define __MUSEVAR 100
#define __MBUF 1200
#define __MWERT 2048
#define __MPAR 400
#define __MTH 40
#define __MTPAR 40
#define __MTCAL 40
#define __MCOUP 10
#define __MAXFORMLENGTH 1024
#define __MAXITEMLENGTH 80
#define __MAXNUMSTACK 50
#define __MAXOPSTACK 50
#define __MUSRFSTACK 500
#define __NODELIMS 7
#define __MSMPL 10000
#define __MDIM 2048
! ---- communication common block containig the analysed inputline       
!      comand   = actual command keyword                               
!      vname(*) = names stack                                           
!      rpar(*)  = number stack                                          
!      inames   = length of names stack                                 
!      ipars    = length of number stack                                
!      ioldc    = flag to indicate second command in line               
!      inline   = actual inputline                                      
!      reslin   = residual line conaining all still nonprocessed command
!      inpar(*) = no. of parameters associated with the i-th name       
!      iparn(*) = pointer parameters to names                           
!      inapa(*) = adress of first parameter following the *-th name     
!      arglst(*)= character*20 list of arguments                        
!      iargs    = length of arglst = no. of parameters                  
!      pmlist(*,i) => dummy, left only for compatibility reasons        
!      ipmls    =       "                    "                          
!      iolbuf   = buffer for saving ioldc during a makro call           
!      rlbuf    = buffer to save reslin during a makro call             
!      lstpar   = number of the last decoded parameter                  
!      lstnam   = number of the last decoded name                       
                                                                        

	module cincom
	    save
	    real*8 rpar(__MINC)
	    integer inames
	    integer ipars
	    integer ioldc
	    integer inpar(__MINC)
	    integer iparn(__MINC)
	    integer inapa(__MINC)
	    integer iargs
	    integer ipmls
	    integer iolbuf !Can we delete this and put in a local var?
	    integer lstpar
	    integer lstnam
	end module cincom

! ---- the current command is stored on comand                          
!      the names on stack vname(*)                                      
!      the numbers on stack rpar(*)                                     
!      in common/cincom/                                                
!
	module cincoc
		save
		character*8 comand
		character*8 vname(__MINC)
		character*1024 title
		character*1024 reslin
		character*1024 inline
		character*20 arglst(__MINC)
		character*20 pmlist(__MINC,2)
		character*1024 rlbuf
	end module cincoc

	module icpathes
		save
		character*1024 data_path
		character*1024 save_path
		character*1024 makro_path
		character*1024 datreat_path,PWD
		character*1024 history(0:20)
	end module icpathes


! --- variables for makro-parameter-passing ---                         
!     argvals(i)  = parameterlist at a call of a makro, only temp.      
!     pmlst(k,i,1..2) = replace list, in the k-th makro-layer           
!                       all substrings=pmlst(k,*,1) will be             
!                       replaced by pmlst(k,*,2)                        
!     iargvs      = number of arguments                                 
!     ipmlst(k)   = number of given replacing items for makro k         
!     kanal(k)    = fortran io-number associated to makro layer k       
!     ktop        = current makro layer (0=keyboard-input)              
!                                                                       
 

	module cmargs
		save
		character*1024 argvals(__MINC)
		character*80 pmlst(__MDEPTH,__MINC,2)
	end module cmargs
	
	module imargs
		save
		integer iargvs
		integer ipmlst(__MDEPTH)
		integer kanal(0:__MDEPTH)
		integer ktop
		data kanal/ 5,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58, 59/
       		data ktop /0/ 
	end module imargs
 
  

       	module xoutxx
		save
		integer :: iot=0
		integer :: ierrr=0
	end module xoutxx

	module xroxxx
		save
		real*8 xyorig(3)
		real*8 rotvec(3)
		data xyorig/3*0.d0/, rotvec/3*0.d0/
        end module xroxxx 

! usenum and usevar are now combined
	module usevar
		save
		real*8 useval(__MUSEVAR)
		integer nousev 
		character*16 usenam(__MUSEVAR)
	end module usevar

! ---- common containing all the scans ----                             
!                                                                       
! --- xwerte(i,j)   x-values on buffer j                                
!     ywerte(i,j)   y-values on buffer j                                
!     yerror(i,j)   error of y-values   (only supported by some fktn)   
!     xname(j)      name of x-values (x-axis) for buffer j              
!     yname(j)              y-        y-                            
!     name(j)       short text-identifier for data on buffer j          
!     nwert(j)      no. of valid points on buffer j                     
!     nbuf          no. of filled buffers                               
!     numor(j)      numerical identification of data on buffer j          
!     coment(j)     one line of comment describing data on buffer j     
!     params(l,j)   set of parameters associated with data on buffer j  
!     napar(l,j)    names of these parameters                           
!     nopar(j)      no. of valid parameters                             
!                                                                       


	MODULE cdata
		save
		real xwerte(__MWERT,__MBUF)
		real ywerte(__MWERT,__MBUF)
		real yerror(__MWERT,__MBUF)
     		character*80 xname(__MBUF)
		character*80 yname(__MBUF)
		character*80 name(__MBUF)
		integer nwert(__MBUF)
		integer numor(__MBUF)
		integer :: nbuf = 0
		character*80 coment(__MBUF)*80
		real params(__MPAR,__MBUF)
		character*80 napar(__MPAR,__MBUF)
		integer nopar(__MBUF)


        CONTAINS

        subroutine DataCopy( isource, idestination)
!       -------------------------------------------
        implicit none
        integer, intent(inout) :: isource, idestination
        integer                :: i
        
        if(isource.le.0 .or. isource.gt.__MBUF) then
          write(6,*)'cdata:DataCopy: isource out of range: ',isource,' NO ACTION '
          return
        endif
        if(idestination.eq.0) idestination = nbuf+1
        if(idestination.le.0 .or. idestination.gt.__MBUF) then
          write(6,*)'cdata:DataCopy: idestination out of range: ',idestination,' NO ACTION '
          return
        endif

        xwerte(1:__MWERT,idestination) =  xwerte(1:__MWERT,isource)
        ywerte(1:__MWERT,idestination) =  ywerte(1:__MWERT,isource)
        yerror(1:__MWERT,idestination) =  yerror(1:__MWERT,isource)
        xname(idestination)            =  xname(isource)
        yname(idestination)            =  yname(isource)
        name(idestination)             =  name(isource)
        nwert(idestination)            =  nwert(isource)
        numor(idestination)            =  numor(isource)
        coment(idestination)           =  coment(isource)
        params(1:__MPAR,idestination)  =  params(1:__MPAR,isource)
        napar(1:__MPAR,idestination)   =  napar(1:__MPAR,isource)
        nopar(idestination)            =  nopar(isource)

        if(idestination.gt.nbuf) nbuf=idestination

        end subroutine DataCopy
!       -------------------------------------------


        subroutine DataGet( isource, x, y, dy, n)
!       -----------------------------------------
        implicit none
        integer, intent(in)           :: isource
        integer, intent(out)          :: n
        double precision, intent(out) :: x(__MWERT), y(__MWERT), dy(__MWERT)
        
        if(isource.le.0 .or. isource.gt.nbuf) then
          write(6,*)'cdata:DataGet: isource out of range: ',isource,' NO ACTION '
          return
        endif

        n             =  nwert(isource)
        x (1:n)       =  xwerte(1:n,isource)
        y (1:n)       =  ywerte(1:n,isource)
        dy(1:n)       =  yerror(1:n,isource)
 
        return

        end subroutine DataGet
!       -------------------------------------------


        subroutine DataPut( idestination, x, y, dy, n)
!       ----------------------------------------------
        implicit none
        integer, intent(inout)        :: idestination
        integer, intent(in)           :: n
        double precision, intent(out) :: x(__MWERT), y(__MWERT), dy(__MWERT)
        
        if(idestination.eq.0 .and. nbuf .lt. __MBUF) then
          nbuf         = nbuf+1 
          idestination = nbuf
        endif

        if(idestination.le.0 .or. idestination.gt.__MBUF) then
          write(6,*)'cdata:DataPut: idestination out of range: ',idestination,' NO ACTION '
          return
        endif

        xwerte(1:__MWERT,idestination)  =  x (1:__MWERT) 
        ywerte(1:__MWERT,idestination)  =  y (1:__MWERT) 
        yerror(1:__MWERT,idestination)  =  dy(1:__MWERT)

        nwert(idestination)             =  n
          
        return

        end subroutine DataPut
!       -------------------------------------------


        !! call parset (vname(1),sngl(rpar(1)),iaddp)
        !! call parget('thick   ',thick,ia,ier)

	END MODULE cdata


! ---- outputlevel                                                      
! xxxx,yyyy = aktuelle addersse fuer wertextraktion xx,yy               
!
	module outlev 
		save
		integer :: iout = 0 
		integer :: ibild = 0
		integer ierrs
		integer :: inka = 5
		integer iibuf
		real xxxx
		real yyyy
		real ptxf(20)
		real yyee 
	end module outlev
! ----- theories common block and definitions ----                      
!                                                                       
!  thenam(i)  = name of i-th theory                                     
!  thparn(j,i)= name of j-th parameter of i-th theory                   
!  nthpar(i)  = no of parameters required for i-th theory               
!  thparx(j,l)= parameters for the l-th activated theory                
!  thpsca(j,l)= corresponding fit scales                              
!  nthtab(l)  = no. of l-th activated theory                            
!  ntheos     = total no. of activated theories                         
!  multflg    = flag indicating a multiplicative theory                 
!                                                                       

	module theory
		save
		character*8 thenam(__MTH)
		character*8 thparn(__MTPAR,__MTH)
		integer nthpar(__MTH)
		real thparx(__MTPAR,__MTCAL)
		real thpsca(__MTPAR,__MTCAL)
		integer nthtab(__MTCAL)
		integer ntheos
		integer multflg(__MTCAL)
	        integer, parameter :: l3=__MTH*__MTPAR
                data thparn/l3*'        '/
	end module theory

	! ---- common containing a selected list of spectra ----                
	!  isels(i)  = address (/spectr/) of selected scan                      
	!  ifits(i)  = adress of fitted spectrum (isels(i))                     
	!  nsel      = length of this table                                     
	module selist
		save
		integer isels(__MBUF)
		integer ifits(__MBUF)
		integer nsel
		integer :: numpls = 10000
	end module selist
 
	!  isfits    = address (/spectr/) of selected fits                      
	!  nfsel     = length of this table                                     
	module fslist
		save
		integer isfits(__MBUF)
		integer nfsel
	end module fslist

! ------ coupling of theory-parameters -----------------------          
!  --- thpala(i,j) ... label of i-th parameter , j-th theory            
!      thpalc(l,i,j) . labels (l) that point to parameters that are     
!                      coupleded to the i-th parameter of the j-th theo.
!      thpafc(l,i,j) . associated proportionality factor                
!      thpaco(i,j) ... offset value for a parameter that is derived by  
!                      couplede other parameters                        
!      ncoup(i,j) .... number of parameters that are coupled to i,j-par.
!

        module theorc 
		save
		character*4 thpala(__MTPAR,__MTCAL)
		character*4 thpalc(__MCOUP,__MTPAR,__MTCAL)
		real thpafc(__MCOUP,__MTPAR,__MTCAL)
		real thpaco(__MTPAR,__MTCAL)
		integer ncoup(__MTPAR,__MTCAL)

		! ---- this blockdata is used to replace the standard filling of        
		!      character variables with nulls ! by blanks !                     
		!      this is importatnt for the incom parser, which is only sensitiv  
		!      to blanks, whereas nulls are written and read if default         
		!      character varaibles are printed!                                 
		! ---- dervived auxiliary parameters only for block data ----
       		integer, parameter :: l1=__MTPAR*__MTCAL, l2=__MCOUP*__MTPAR*__MTCAL
		data thpala/l1*'    '/,thpalc/l2*'    '/
        end module theorc

! ------ errors of fit --------------------------------------------     
!      ------ estimates of 1 sigma errros -----------------------       
	module therrc
		save
		real therro(__MTPAR,__MTCAL)
	end module therrc	

! range definition of theories (only to be evaluated if parameter thrap
!                               given range)
	module thparc
		save
		character*8 thrapar(__MTH)
		real*4 thramin(__MTH)
		real*4 thramax(__MTH)
	end module thparc

	module formnu
		save
		real*8 numstack(__MAXNUMSTACK)
		real*8 :: degree=1.d0
		real*8 valnum
		integer priostack(0:__MAXOPSTACK)
		integer topnumstack
		integer topopstack
		integer tusrfstack
		integer klammerprio
		integer actchar
		integer len
		integer litem
		logical ok
		logical error
		logical :: say=.false.
	end module formnu

	module formch
		save
		character*1 formula(0:__MAXFORMLENGTH)
		character*1 item(0:__MAXITEMLENGTH)
		character*1 delims(0:__NODELIMS)
		character*4 typ
		character*4 opstack(__MAXOPSTACK) 
		character*20 usrfstack(__MUSRFSTACK)
                !needed in subroutine getitem - dirty!
                character*(__MAXITEMLENGTH+1) citemx
                equivalence(citemx,item(0))
	end module formch


	module formul
		save
		character*1024 :: xformel = '(xx)'
		character*1024 :: yformel = '(yy)'
		character*1024 yfitform
	end module formul

! --- parameters of last spline smoothing + value qziel                 
!     qziel is the value at witch the spline should be evaluated        
!     numspl numor of spline fitted data                                
!     nwspl  length of splined data vectors                             

	module cfc
		save
		real qziel
		real cscoef(4,__MWERT)
		real break(__MWERT)
		real weight(__MWERT)
		integer :: numspl = 0
		integer :: nwspl = 0
	end module cfc

! ---- communication with subr. func ---                                
!  --- iprt     : printing during func calculation (set by fit)         
!      x1       : lower limit of fit range (if autox1 = .false.)        
!      x2       : upper limit of fit range (if autox2 = .false.)        
!      autox1/2 : if true the genuine limits of the datafiled is taken  
!                 else x1 /x2 values are taken (set by thc)             
!      ferror   : fehler des zu vergleichenden datenfeldes              
! ----------------------------------------------------------------------

	module cfunc
		save
		integer iprt
		logical sqwght
		real x1
		real x2
		logical :: autox1 = .true.
		logical :: autox2 = .true.
		real ferror(__MSMPL)
	end module cfunc


! ---  if lerrel=.true. the fiterrors are take as relative errors 
	module cfunce
		save
		logical :: lerrel = .false.
		logical :: lwrtfitdat = .false.
		real fcssq
	end module cfunce 

! --- echoform related parameters ( nse-programs ) ----                 
!     j1, j2  = feldintegrale in gauss*meter                            
!     j0delta = feldintegralvariation unabh&a.ngig von j1,j2            
!     cdelta  = koeffizient der feldintegralabh&a.ngigkeit von j1,j2    
 
	module partran
		save
		real*4 j1echo
		real*4 j2echo
		real*4 j0delta
		real*4 cdelta
	end module partran 

!     alam0   = mittlere wellenlaenge in a                              
!     dalam   = wellenla.ngenbreite in a                                

	module wlntran
		save
		real alam0
		real dalam 
	end module wlntran

!     tau in sekunden, relaxationszeit des streuers
	module sqtran
		save
		real tau
	end module sqtran 


!             ---> this common saves storage by use of ca in fftrmx also
	module fftwrk
		save
		complex*8 ca(__MDIM+1,__MDIM+1) 
	end module fftwrk

	module fftwr1
		save
		real ca(__MDIM,__MDIM) 
                real yinter(__MDIM)
	end module fftwr1

	module fftwr2
		save
		real ca(__MDIM+1,__MDIM+1) 
	end module fftwr2


	module constants
		save
		! --- minc = incom stack depth                                          
		integer, parameter :: minc=__MINC
      		integer, parameter :: mdepth=__MDEPTH
	        integer, parameter :: mwert=__MWERT, mbuf=__MBUF, mpar=__MPAR
		! --- mwert  = max. no. of x-y-values in one buffer                     
		!     mbuf   = max. no. of different buffers                            
		!     mpar   = max. no. of parameters associated with one buffer        
		! ---  maximum scan length ....                                         
       		integer, parameter:: mth=__MTH, mtpar=__MTPAR,mtcal=__MTCAL,mcoup=__MCOUP
		! ---  fit dimensions ---                                               
		!  -- mfit = max no. of fitted parameters                               
		!     msmpl= max no. of datapoints in fit                               
       		integer, parameter :: mfit=100,msmpl=__MSMPL
		integer, parameter :: musevar=__MUSEVAR
		integer, parameter :: maxformlength=__MAXFORMLENGTH
		integer, parameter :: maxitemlength=__MAXITEMLENGTH
		integer, parameter :: maxnumstack=__MAXNUMSTACK 
		integer, parameter :: maxopstack=__MAXOPSTACK
		integer, parameter :: musrfstack=__MUSRFSTACK
		integer, parameter :: nodelims=__NODELIMS
		integer, parameter :: mdim=__MDIM 
	        integer, parameter :: lda=__MDIM+1
	end module constants


