#define __MINC 40
#define __MDEPTH 20
#define __MUSEVAR 100
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
	    real*8 rpar(__MINC)
	    integer inames
	    integer ipars
	    integer ioldc
	    integer inpar(__MINC)
	    integer iparn(__MINC)
	    integer inapa(__MINC)
	    integer iargs
	    integer ipmls
	    integer iolbuf
	    integer lstpar
	    integer lstnam
	end module cincom

! ---- the current command is stored on comand                          
!      the names on stack vname(*)                                      
!      the numbers on stack rpar(*)                                     
!      in common/cincom/                                                
!
	module cincoc
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
		character*1024 data_path
		character*1024 save_path
		character*1024 makro_path
		character*1024 home
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
		character*1024 argvals(__MINC)
		character*80 pmlst(__MDEPTH,__MINC,2)
	end module cmargs
	
	module imargs
		integer iargvs
		integer ipmlst(__MDEPTH)
		integer kanal(0:__MDEPTH)
		integer ktop
		data kanal/ 5,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58, 59/                    
       		data ktop /0/ 
	end module imargs
 
  

       	module xoutxx
		integer iot
		integer ioold
		integer ibild1
		integer ierrr
		integer inka1
		logical cray 
		data cray/.false./
		data iot/0/, inka1/5/, ierrr/0/, ibild1/0/, ioold/0/
	end module xoutxx
	
	module xroxxx
		real*8 xyorig(3)
		real*8 rotvec(3)
		data xyorig/3*0.d0/, rotvec/3*0.d0/
        end module xroxxx 

! usenum and usevar are now combined
	module usevar
		real*8 useval(__MUSEVAR)
		integer nousev 
		character*16 usenam(__MUSEVAR)
	end module usevar



! --- minc = incom stack depth                                          
	module constants
		save
		integer, parameter :: minc=__MINC
      		integer, parameter :: mdepth=__MDEPTH
	        integer, parameter :: mwert=1024, mbuf=200, mpar=200
		! --- mwert  = max. no. of x-y-values in one buffer                     
		!     mbuf   = max. no. of different buffers                            
		!     mpar   = max. no. of parameters associated with one buffer        
		! ---  maximum scan length ....                                         
       		integer, parameter:: mth=40, mtpar=40,mtcal=40,mcoup=10
		! ---  fit dimensions ---                                               
		!  -- mfit = max no. of fitted parameters                               
		!     msmpl= max no. of datapoints in fit                               
       		integer, parameter :: mfit=40,msmpl=4000
		integer, parameter :: musevar=__MUSEVAR
	end module constants
