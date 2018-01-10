
  module dimensions
   save
   integer, parameter :: c_MINC        = 40
   integer, parameter :: c_MDEPTH      = 20
   integer, parameter :: c_MUSEVAR     = 300
   integer, parameter :: c_MBUF        = 1200
   integer, parameter :: c_MWERT       = 20000  ! was 10000
   integer, parameter :: c_MPAR        = 400
   integer, parameter :: c_MTH         = 80
   integer, parameter :: c_MTPAR       = 40
   integer, parameter :: c_MTCAL       = 40
   integer, parameter :: c_MCOUP       = 10
   integer, parameter :: c_MAXFORMLENGTH = 1024
   integer, parameter :: c_MAXITEMLENGTH = 1024 !80
   integer, parameter :: c_MAXNUMSTACK   = 50
   integer, parameter :: c_MAXOPSTACK    = 50
   integer, parameter :: c_MUSRFSTACK    = 500
   integer, parameter :: c_NODELIMS      = 7
   integer, parameter :: c_MSMPL         = 20000 ! was 10000
   integer, parameter :: c_MDIM          = 2048
   integer, parameter :: c_MFIT          = 100
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
                                                                        
  end module dimensions


	module cincom
            use dimensions
	    save
	    real*8 rpar(c_MINC)
	    integer inames
	    integer ipars
	    integer ioldc
	    integer inpar(c_MINC)
	    integer iparn(c_MINC)
	    integer inapa(c_MINC)
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
                use dimensions
		save
		character*8 comand
		character*8 vname(c_MINC)
		character*1024 title
		character*1024 reslin
		character*1024 inline
		character*20 arglst(c_MINC)
		character*20 pmlist(c_MINC,2)
		character*1024 rlbuf
                character(len=20) :: prompt = "--> "
	end module cincoc

	module icpathes
                use dimensions
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
                use dimensions
		save
		character*1024 argvals(c_MINC)
		character*80 pmlst(c_MDEPTH,c_MINC,2)
                logical :: argquoted(c_MINC)
		integer iargvs
	end module cmargs
	
	module imargs
                use dimensions
		save
!		integer iargvs
		integer ipmlst(c_MDEPTH)
		integer kanal(0:c_MDEPTH)
		integer ktop
		data kanal/ 5,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58, 59/
       		data ktop /0/ 
	end module imargs
 
  

       	module xoutxx
		save
		integer :: iot=0
		integer :: ierrr=0
                logical :: mask_err = .false.
	end module xoutxx

	module xroxxx
		save
		real*8 xyorig(3)
		real*8 rotvec(3)
		data xyorig/3*0.d0/, rotvec/3*0.d0/
        end module xroxxx 

! usenum and usevar are now combined
	module usevar
                use dimensions
		save
		real*8 useval(c_MUSEVAR)
		integer nousev 
		character*16 usenam(c_MUSEVAR)
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
                use dimensions
		save
		real xwerte(c_MWERT,c_MBUF)
		real ywerte(c_MWERT,c_MBUF)
		real yerror(c_MWERT,c_MBUF)
     		character*80 xname(c_MBUF)
		character*80 yname(c_MBUF)
		character*80 name(c_MBUF)
		integer nwert(c_MBUF)
		integer numor(c_MBUF)
		integer :: nbuf = 0
		character*80 coment(c_MBUF)*80
		real params(c_MPAR,c_MBUF)
		character*80 napar(c_MPAR,c_MBUF)
		integer nopar(c_MBUF)
                integer :: params_display_level(c_MPAR,c_MBUF) = 0


        CONTAINS

        subroutine DataCopy( isource, idestination)   
!       -------------------------------------------   
        implicit none
        integer, intent(inout) :: isource, idestination
        integer                :: i, j, np
        
        if(isource.le.0 .or. isource.gt.c_MBUF) then
          write(6,*)'cdata:DataCopy: isource out of range: ',isource,' NO ACTION '
          return
        endif
        if(idestination.eq.0) idestination = nbuf+1
        if(idestination.le.0 .or. idestination.gt.c_MBUF) then
          write(6,*)'cdata:DataCopy: idestination out of range: ',idestination,' NO ACTION '
          return
        endif

        xwerte(:,idestination)         =  xwerte(:,isource)
        ywerte(:,idestination)         =  ywerte(:,isource)
        yerror(:,idestination)         =  yerror(:,isource)
        xname(idestination)            =  xname(isource)
        yname(idestination)            =  yname(isource)
        name(idestination)             =  name(isource)
        nwert(idestination)            =  nwert(isource)
        numor(idestination)            =  numor(isource)
        coment(idestination)           =  coment(isource)

        np = nopar(isource)
        j  = 0
        do i=1,np
         if(index(napar(i,isource),"!") == 0 ) then
          j = j+1
          params(j,idestination)                =  params(i,isource)
          params_display_level(j,idestination)  =  params_display_level(i,isource)
          napar(j,idestination)                 =  napar(i,isource)
         endif
        enddo
        nopar(idestination) = j

        if(idestination.gt.nbuf) nbuf=idestination

        end subroutine DataCopy
!       -------------------------------------------


        subroutine DataGet( isource, x, y, dy, n)
!       -----------------------------------------
        implicit none
        integer, intent(in)           :: isource
        integer, intent(out)          :: n
        double precision, intent(out) :: x(c_MWERT), y(c_MWERT), dy(c_MWERT)
        
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
        double precision, intent(out) :: x(c_MWERT), y(c_MWERT), dy(c_MWERT)
        
        if(idestination.eq.0 .and. nbuf .lt. c_MBUF) then
          nbuf         = nbuf+1 
          idestination = nbuf
        endif

        if(idestination.le.0 .or. idestination.gt.c_MBUF) then
          write(6,*)'cdata:DataPut: idestination out of range: ',idestination,' NO ACTION '
          return
        endif

        xwerte(1:c_MWERT,idestination)  =  x (1:c_MWERT) 
        ywerte(1:c_MWERT,idestination)  =  y (1:c_MWERT) 
        yerror(1:c_MWERT,idestination)  =  dy(1:c_MWERT)

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
                use dimensions
		save
		character*8 thenam(c_MTH)
		character*8 thparn(c_MTPAR,c_MTH)
		integer nthpar(c_MTH)
		real thparx(c_MTPAR,c_MTCAL)
		real thpsca(c_MTPAR,c_MTCAL)
		integer nthtab(c_MTCAL)
		integer ntheos
		integer multflg(c_MTCAL)
	        integer, parameter :: l3=c_MTH*c_MTPAR
                data thparn/l3*'        '/
	end module theory

	! ---- common containing a selected list of spectra ----                
	!  isels(i)  = address (/spectr/) of selected scan                      
	!  ifits(i)  = adress of fitted spectrum (isels(i))                     
	!  nsel      = length of this table                                     
	module selist
                use dimensions
		save
		integer isels(c_MBUF)
		integer ifits(c_MBUF)
		integer nsel
		integer :: numpls = 10000
	end module selist
 
	!  isfits    = address (/spectr/) of selected fits                      
	!  nfsel     = length of this table                                     
	module fslist
                use dimensions
		save
		integer isfits(c_MBUF)
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
                use dimensions
		save
		character*4 thpala(c_MTPAR,c_MTCAL)
		character*4 thpalc(c_MCOUP,c_MTPAR,c_MTCAL)
		real thpafc(c_MCOUP,c_MTPAR,c_MTCAL)
		real thpaco(c_MTPAR,c_MTCAL)
		integer ncoup(c_MTPAR,c_MTCAL)

		! ---- this blockdata is used to replace the standard filling of        
		!      character variables with nulls ! by blanks !                     
		!      this is importatnt for the incom parser, which is only sensitiv  
		!      to blanks, whereas nulls are written and read if default         
		!      character varaibles are printed!                                 
		! ---- dervived auxiliary parameters only for block data ----
       		integer, parameter :: l1=c_MTPAR*c_MTCAL, l2=c_MCOUP*c_MTPAR*c_MTCAL
		data thpala/l1*'    '/,thpalc/l2*'    '/
        end module theorc

! ------ errors of fit --------------------------------------------     
!      ------ estimates of 1 sigma errros -----------------------       
	module therrc
                use dimensions
		save
		real therro(c_MTPAR,c_MTCAL)
	end module therrc	

! range definition of theories (only to be evaluated if parameter thrap
!                               given range)
	module thparc
                use dimensions
		save
		character*8 thrapar(c_MTH)
		real*4 thramin(c_MTH)
		real*4 thramax(c_MTH)
	end module thparc

	module formnu
                use dimensions
		save
		real*8 numstack(c_MAXNUMSTACK)
		real*8 :: degree=1.d0
		real*8 valnum
		integer priostack(0:c_MAXOPSTACK)
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
                use dimensions
		save
		character*1 formula(0:c_MAXFORMLENGTH)
		character*1 item(0:c_MAXITEMLENGTH)
		character*1 delims(0:c_NODELIMS)
		character*4 typ
		character*4 opstack(c_MAXOPSTACK) 
		character*20 usrfstack(c_MUSRFSTACK)
                !needed in subroutine getitem - dirty!
                character*(c_MAXITEMLENGTH+1) citemx
                character*(c_MAXFORMLENGTH+1) cformulax
                equivalence(citemx,item(0))
                equivalence(cformulax,formula(0))
	end module formch


	module formul
                use dimensions
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
                use dimensions
		save
		real qziel
		real cscoef(4,c_MWERT)
		real break(c_MWERT)
		real weight(c_MWERT)
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
                use dimensions
		save
		integer iprt
                integer :: icall
		logical sqwght
		real x1
		real x2
		logical :: autox1 = .true.
		logical :: autox2 = .true.
		real ferror(c_MSMPL)
                real    :: xinitial(c_MFIT)
                double precision :: pardev_scale = 0d0
	end module cfunc


! ---  if lerrel=.true. the fiterrors are take as relative errors 
	module cfunce
                use dimensions
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
                use dimensions
		save
		complex*8 ca(c_MDIM+1,c_MDIM+1) 
	end module fftwrk

	module fftwr1
                use dimensions
		save
		real ca(c_MDIM,c_MDIM) 
                real yinter(c_MDIM)
	end module fftwr1

	module fftwr2
                use dimensions
		save
		real ca(c_MDIM+1,c_MDIM+1) 
	end module fftwr2


	module constants
                use dimensions
		save
		! --- minc = incom stack depth                                          
		integer, parameter :: minc=c_MINC
      		integer, parameter :: mdepth=c_MDEPTH
	        integer, parameter :: mwert=c_MWERT, mbuf=c_MBUF, mpar=c_MPAR
		! --- mwert  = max. no. of x-y-values in one buffer                     
		!     mbuf   = max. no. of different buffers                            
		!     mpar   = max. no. of parameters associated with one buffer        
		! ---  maximum scan length ....                                         
       		integer, parameter:: mth=c_MTH, mtpar=c_MTPAR,mtcal=c_MTCAL,mcoup=c_MCOUP
		! ---  fit dimensions ---                                               
		!  -- mfit = max no. of fitted parameters                               
		!     msmpl= max no. of datapoints in fit                               
       		integer, parameter :: mfit=c_MFIT,msmpl=c_MSMPL
		integer, parameter :: musevar=c_MUSEVAR
		integer, parameter :: maxformlength=c_MAXFORMLENGTH
		integer, parameter :: maxitemlength=c_MAXITEMLENGTH
		integer, parameter :: maxnumstack=c_MAXNUMSTACK 
		integer, parameter :: maxopstack=c_MAXOPSTACK
		integer, parameter :: musrfstack=c_MUSRFSTACK
		integer, parameter :: nodelims=c_NODELIMS
		integer, parameter :: mdim=c_MDIM 
	        integer, parameter :: lda=c_MDIM+1
	end module constants




 module theory_description
   use dimensions
   use theory
   save
        integer, parameter         :: M_recin_par  = 100 
        integer, parameter         :: M_recout_par = 100 
        character(len=1),parameter :: cr = char(10)
        character(len=16),parameter:: parspace = "                "
	character(len=8)           :: th_identifier(c_MTH) = " "
	character(len=1024)        :: th_explanation(c_MTH)
	character(len=1024)        :: th_citation(c_MTH) = " "
	character(len=1024)         :: th_param_desc(c_MTPAR,c_MTH)
	character(len=1024)         :: th_file_param(M_recin_par,c_MTH)
	character(len=1024)         :: th_out_param(M_recout_par,c_MTH)
	integer, private           :: nthdesc = 0
        integer                    :: idesc
   
   contains

        integer function next_th_desc() 
          implicit none
          if(nthdesc < c_MTH) then
            nthdesc = nthdesc+1
          else
            call errsig(9999,"Theory description fault!")
          endif
          next_th_desc = nthdesc
        end function next_th_desc

       logical function output_th_explanation( thn ) result(ok)
          implicit none
          character(len=8), intent(in) :: thn
          
          integer :: i, ith, ipa
          logical :: parout
       
!          write(6,*)'...searching for : ',thn, 'among: ',nthdesc


nt:       do i=1,nthdesc
            if(thn == th_identifier(i)) then
              write(6,'(a)')"-----------------------------------------------------------------------------------"
              write(6,'("theory id-name: ",a)')thn
              write(6,'(a)')"-----------------------------------------------------------------------------------"
              write(6,'(a)')trim(th_explanation(i))
              write(6,'(a)')"-----------------------------------------------------------------------------------"
dt:          do ith=1, c_MTH
                if(thenam(ith) == thn) then
                  write(6,'(a,i2,a)') "Parameters(",nthpar(ith),"): "
                  do ipa=1,nthpar(ith)
                    write(6,'(i3,": ",a8," > ",a)') ipa,thparn(ipa,ith), trim(th_param_desc(ipa,i))
                  enddo
                  exit dt
                endif
              enddo dt  

             parout = .false.
             do ipa = 1, 1,M_recin_par
               parout = (len_trim(th_file_param(ipa,i)) > 0) .or. parout
             enddo
             if(parout) then
               write(6,'(a)')"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               write(6,'(a)')"INPUT: Parameters that are extracted from the actual considered data records:"
                write(6,'(a)')"..there may be default assumptions, but better make sure that these parameters are set properly!"
  dt2:          do ipa = 1,M_recin_par
                  if(len_trim(th_file_param(ipa,i)) > 0)  write(6,'(i3,": ",a)') ipa, trim(th_file_param(ipa,i))
                enddo dt2
             endif

             parout = .false.
             do ipa = 1, 1,M_recout_par
               parout = (len_trim(th_out_param(ipa,i)) > 0) .or. parout
             enddo
             if(parout) then
              write(6,'(a)')"==================================================================================="
              write(6,'(a)')"OUTPUT: Parameters that are computed and added to the records parameters as information:"
  dt3:         do ipa = 1,M_recout_par
                 if(len_trim(th_out_param(ipa,i)) > 0)  write(6,'(i3,": ",a)') ipa, trim(th_out_param(ipa,i))
               enddo dt3
            endif
          
              if(len_trim(th_citation(i)) > 1) then
                write(6,'(a)')"..................................................................................."
                write(6,'("cite: ",a," !")') trim(th_citation(i))
              endif
              write(6,'(a)')"-----------------------------------------------------------------------------------"
              ok = .true.
              return
            endif
          enddo nt

          ok = .false.
!          write(6,'(a,a)')"... no further info available for: ",thn

       end function output_th_explanation
 

 end module theory_description
