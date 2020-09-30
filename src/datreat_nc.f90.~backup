! INCLUDE "commons.h"

! **********************************************************************
! *    treating any data                                               *
! *    imsl version 10 fitprogram test                                 *
! *                                                                    *
! *    datreat version from 18.5.88                                    *
! *    m. monkenbusch, iff, kfa-juelich                                *
! *    neu 3.11.93: multiply flag fuer theorien                        *
! *    neu 28.1.99: thc -n x1 <x1> x2 <x2> on log scale                *
! *    neu 09.3.99: fun zimm, rouse etc + error                        *
! *                 plot noparplo , achsen lange strings               *
! *                 plot error bars bei log ok                         *
! *                 purge loescht tatsaechlich (auch errors)           *
! *                 sel add                                            *
! *    neu 04.02.00:seterr                                             *
! *                 generiere Fehler aus Daten (Fitweighting)          *
! *                 ueber yformel                            (RCS:1.3) *
! *                                                                    *
! *    neu 23.02.00:fit map ..                                         *
! *    neu 19.05.00:dparam('name   ') fuer th-programmierung verfuegbar*
! *    neu         :msave             fuer Mehrfachsave                *
! *                                                                    *
! *                                                                    *
! *    BSS-Version datreat05_1                                         *
! *    1. sel with new options: next, all, range                       *
! *                                                                    *
! **********************************************************************
!
  program datreat
      use dimensions
      use new_com
      ! use cincoc
      ! use cmargs
      use xoutxx
      use cdata
      !  use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use therrc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
!      use constants
      use PhysicalConstantsPlus
      use unift
      use theory_description
      use gr_interface
      use lmfit
      implicit none

      integer iadda
      common /thiadd/iadda


       integer isigint
       real yyy, zpar, ywr, yw2, yy, ysum, ynno, yval, ymm2, ymaxf,ymax
       real yw1, yminf, ym1, yfamp, yerr, yer, yampf, y3, yctab, y11, y1
       real y00, y2, xz, y, y0, xx, xval, xsum, yinterp, xr, xnc, xmon2
       real xmonitor, xmon1, xn, xmax, xk00, fsplin, xm, ximaxf, xi, ximaxq
       real xcut, x1raster, xk0, x2raster, x, val, tt, wl4, thval, x01, x02
       real theta2, thick, temp_z, trans, t1, sump, sumb, sx, slope, smpar
       real t2, sy, sigma, selpartol, r0, selparval, ri, result, qend, qstart
       real pkave, pfkave, pfave, qz, qrs, qmax, pxxd, pkqave, qortho
       real pfqave, p3qave, p3ave, parval_x, gunten, fy, fwid, goben
       real fia, famp, facto2, eta_s, f1, errmx, facto1, errrec, errmax
       real erraec, echo, dxx, fx, smirro, fai, esum, errest, errret, df
       real delqv, ddj, delqh, dj, draster, d0, d1, d, d2, blim, bklen, a2
       real a, alam, aa, akb_z, ai0, ampli, amp, apar, alim, an0, bfave, bkgr
       real bkqave, bkave, dq, dt, dx, bpar, bfkave, aver, bqkave, detdis, rrv
       integer nw, nx, numr,  num2, num1, nra, npp, npoint, nraster
       integer npk, npay, npl, npax, np, nold, nnpar, nnn, nnumi, nn4, nn3
       integer nnw, nn1, nn, nn2, nfits, newadd, nneu, ncp, ncpl, ncut, nech
       integer nbuff, newnum, nb, n1, n2, n3, n0, nfft, n, msel, mpk, modeqc
       integer maxint, m, lyy, lxx, ll, kk, mcut, k1, k2, kz, k, jpar, jl, l
       integer j, j2, ival, kan1, ithc, iy, j1, iss, ispc, isl, ira, iseed, ir
       integer ito, isel, ip, imx, imaxf, iminf, ikz, ilas, inew, ii
       integer ifunx, ifuny, ik, iery, ies, idimux, ier, iprs, ifrom, ierr
       integer ierx, iequal, ides, icut, id, ib, ileng, ic, iauto, iaddp 
       integer iaddf, iad2, iadd, iad1, ia1, ia2, ia3, ia4, ia5, ia6, i, ia
       integer itcal, ias, igrand

       real :: xh  

       character :: fpnam*8,xpnam*8,ypnam*8
       character fsname*1024
       character(len=1024) :: charbuf
       character*1 csel
       real*8 val8x,val8y

       complex ca(mwert), cb(mwert), cmplx, conjg
       dimension alim(2), blim(2)
!
       dimension nnumi(size(inpar))
       logical :: withmo = .true.,fonpla = .false.

       parameter( mcut=200 )
       dimension rrv(0:mcut)

! ---> transfer der addresse der gerade bearbeiteten kurve nach thXX
!      damit koennen dann ggf. weitere parameter ueber parget gewonnen
!      werden

! Hilfsvariablen fuer echocurv
       double precision            :: fouriertime
       double precision, parameter :: unit_Gauss = 1d-4



       double precision qortho_0, qortho_1,qz_0, qz_1
       integer          nz, northo, iortho, iz, ith
       character*40     filnam
       character*8      filext
       character*8      selparna

                                      ! dir-command --> comment-length
       integer          len_comm

       character*420    pastring
       character*80     fostring
       character*8      dirpan(20)
       real*4           dirpav(20)

       real          :: y_scaling
       real          :: diffscal, diffrouse



       character(len=80)         :: cbuffer
       character(len=len(title)) :: tbuffer
 
       character(len=len(xname(1)))         :: cbufferx, nbufferx
       character(len=len(yname(1)))         :: cbuffery, nbuffery




!
       real :: errabs=1.0,errel=1.e-2
       real :: xcopy1, xcopy2

!
! lokale Hilfsvariablen fuer arit2 (geht noch besser) !
!
       integer :: nout
       integer, parameter :: mwork=c_MWERT
       real :: yout1(mwork), yerout1(mwork), yout2(mwork), yerout2(mwork)



       real    :: x1int, x2int, sum, sumer
       integer :: i1int, i2int

!
       double precision :: Ommax = 1000d0
       integer          :: nomega = 100

       real             :: He3_pressure_bar = 4.75  ! (IN5 value in 2013)
       real             :: Tube_diameter_m  = 0.025 ! (IN5 1 inch tubes)
       real             :: lenfrac          = 1.0
       real             :: detector_sensitivity, detsens, lambdaA


       double precision :: unift_range_expand     = 1d0
       double precision :: unift_resolution_limit = 0.1d0
       

         real  :: lambda0, lambda1, delta_lambda, temp,angle_2tht 
!?         real  :: e0, e1, dE, de0, dee, dbb, omega, om_cm, qrc, qa, q
         real  :: e0, e1, dE, de0, dee, dbb, omega, om_cm, qa, q
         real  :: kinetic_factor, channel_factor,channel_width,channel_width0
         real  :: sume, sumq
  
       integer :: irstp, recstep = 1   


       real    ::   lower_range 
       real    ::   upper_range
       logical ::   range_is_y  

       logical ::   keep
       logical ::   content_da, usv_da, sel_da
 
       character(len=len(thenam(1))) :: hwork(mth)
       integer                       :: iperm(mth)


       CHARACTER(len=255) :: farg

       integer           :: length, status
       character(len=80) :: argval


!       character(len=cmd_len)  :: mycommand
!       integer                 :: ier
       integer                 :: istat, ixf, iyf
       integer                 :: iselstep
!
       external fdes
       external f
!
       external                :: usrextr
       external                :: usrfun 



! ---- initialisations ----
! ---- error-set ----------
       call erset (0,1,0)

                        !! activates Crtl-C signal handler
       call sigset(-1)
                        !! Clears Signals
       call sig_Reset()



  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! this is for installation of the usr_extract program needed by new_com
  call set_usrextr( usrextr )
  call set_usrfunc( usrfun  )

  ! ------------------------------------------------------------------------

  !! add tab expansions as one adds commands
  istat = add_tab_expansion('help')
  istat = add_tab_expansion('clear')
  istat = add_tab_expansion('quit')
  !




     CALL get_command_argument(1,farg)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        write(6,*)
                        write(6,*)'======================================================='
                        write(6,*)'=   datreat     Version: mm-develop 3.0               ='
                        write(6,*)'=   -------     --------                              ='
                        write(6,*)'=   Author: M.Monkenbusch  R. Biehl, O.Holderer, JCNS ='
                        write(6,*)'======================================================='


                        prompt       = "#datreat => " 
                        history_file = ".last_dtr_commands"


                        write(6,*)
                        write(6,*)
                        write(6,*)
                        write(6,*)'=================================================================='
                        write(6,*)'=  NEW NEW NEW NEW NEW                                           ='
                        write(6,*)'=  sel : advanced functions, type "sel help" to learn more       ='
                        write(6,*)'=  average : combine data,   type "average help"                 ='
                        write(6,*)'=  copy : copy selected, if x1 <x1> x2 <x2> are given            ='
                        write(6,*)'=         copy only that range                                   ='
                        write(6,*)'=  sequence: replace x-values by sequence nubers 1..n            ='
                        write(6,*)'=  swapxy:   exchange x and y values, discard errors             ='
                        write(6,*)'=  numorchg: change numor of selected record                     ='
                        write(6,*)'=  arit:     back to old one                                     ='
                        write(6,*)'=  arit2:    note needs sorted data...                           ='
                        write(6,*)'=  get_th:   activates theory of msave files:: gth <filenam>     ='
                        write(6,*)'=  msave:    unlimited filename length                           ='
                        write(6,*)'=  .... up to 80 theories, 40 th-params, 300 uservars            ='
                        write(6,*)'=  fit: maxfn , ngood  parameters revitalized                    ='
                        write(6,*)'=  .... more commands with help option, try <cmd> help           ='
                        write(6,*)'=  .... improved output formats in plots and lastth              ='
                        write(6,*)'=  .... observes ERROR (and stops makros in expression eval)     ='
                        write(6,*)'=  HINT: use sys to issue system commands including -,+ ...      ='
                        write(6,*)'=        do this always in makros                                ='
                        write(6,*)'=  In parameter list quotes are respected to protect strings     ='
                        write(6,*)'=  (internal incom: we only assume evaluation if 1st char is     ='
                        write(6,*)'=  .(+-1..9)                                                     =' 
                        write(6,*)'=  with fit: parameter parwght  , go, help, couple               =' 
                        write(6,*)'=  new parameter display level ; parlev <dl> ; plot parlev <dl>  =' 
                        write(6,*)'=  range reset if not specified                                  =' 
                        write(6,*)'=  range reset if not specified                                  =' 
                        write(6,*)'=  parameters with potected names (i.e. containing !)            =' 
                        write(6,*)'=   will not be copied to processed data records                 =' 
                        write(6,*)'=   thereby allowing selection of originals                      =' 
                        write(6,*)'=                                                                =' 
                        write(6,*)'=  NEW (9/2018)                                                  =' 
                        write(6,*)'=  serfit      : fits the selecte items one after the other      =' 
                        write(6,*)'=                and creates data records with the parameters    =' 
                        write(6,*)'=  creatser    : creates data record series with parameter var.  =' 
                        write(6,*)'=                as templates                                    =' 
                        write(6,*)'=  fxy         : immediate function fxy & x=(..) y=(..)          =' 
                        write(6,*)'=  csep        : set resline separator, default now is: &        =' 
                        write(6,*)'=  average     : combine points for one or several records       =' 
                        write(6,*)'=  NEW (5/2019)                                                  =' 
                        write(6,*)'=  aligny      : determine scaling to match (e.g SANS col1,2..   =' 
                        write(6,*)'=  NEW (9/2019)                                                  =' 
                        write(6,*)'=  mexp        : automatic matching of sum(n=1..N, an*exp(-t/tn))=' 
                        write(6,*)'=  NEW (2/2020)                                                  =' 
                        write(6,*)'=  Feature: autosave at quit and restore at start !              =' 
                        write(*,*)'= See manual DTRman.pdf !                                        ='
                        write(*,*)'=                                                                ='
                        write(*,*)'=  datreat help   opens Manual                                              ='
                        write(*,*)'=  datreat ni     OR                                             ='
                        write(*,*)'=  datreat 0      starts without reloading previous content      ='
                        write(*,*)'=                 (old behavior)                                 ='
                        write(6,*)'=================================================================='
                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       pi   = 4 * atan(1.0)
       call DataClean(1)

       facto1 = 1
       facto2 = 2

! ---- noise amplitude ----
       ampli   = 1.0
       iseed   = 12347897
! ---- nse-echoform ----
       nech    = 64
       nra     = 10
       amp     = 10
       ncp     = 1
       mpk     = 3
       ic      = 1
       j1echo  = 100.e0
       j2      = 100.e0
       j0delta = 0.e0
       cdelta  = 1.e-6
       alam0   = 8.e0
       dalam   = 0.5e0
       tau     = 1d0
       dj      = 0
       ddj     = 1.358237/(4*alam0)
       famp    = 20.0
       fwid    = 2.0

       erraec  = 1.e-4
       errrec  = 1.e-4
       maxint  = 2000

       sigma   = 1d0

! ---- ecc -------------
      n = 200
      ri= 5.0
      r0= 18.9
      d0= 0.2
      a = 1.2e-3
      ncut = 7
!!!!!!!!!!
       msel = 0
       selpartol = 0.001
       selparval = 0
!!!!!!!!!!
!-------- cailtab -------------
       northo   = 11
       nz       = 11
       qortho_0 = 0.0001d0
       qortho_1 = 0.1001d0
       qz_0     = 0.01d0
       qz_1     = 0.15d0
       filext   = 'last'
!
!-------- fun Rouse, Zimm Scaling -----
       wl4      = 1.0
       eta_s    = 1e-3
                                     ! fuer A**-3ns
       akb_z    = 1e21*1.380662e-23
!

!--------------------------------------
       ! comment-length in dir (display)
       len_comm = 30

       title='Datreat Plot'

!       makro_path='/Users/monk/Desktop/datreat/makros/'     ! here better solution needed ...!
       makro_path=PATH_TO_MAKROS                            ! set in Makefile
       editor    =DEFAULT_EDIT                              !  "        "

       write(6,*)"Actual makro path: ",trim(makro_path)

       call init_theories(thenam,thparn,nthpar,thrapar,thramin,thramax,mth,mtpar)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,2
      call get_command_argument(i,argval,length, status)
      if(status == 0 ) then
        if(argval(1:2) == "he") then   !! TBD nach vorne!
          write(*,*)"Data treatment program datreat"
          write(*,*)"See manual DTRman.pdf !"
          write(*,*)
          write(*,*)"datreat ni    starts without reloading previous content (old behavior)"
          write(*,*)
          call execute_command_line(OPEN_MANUAL)
        endif
        if(argval(1:2) == "ni" .or. argval(1:2) == "0" ) goto 2000
      endif
    enddo

       inquire(file = "last_datreat_content", exist= content_da)
       inquire(file = "lastusv", exist= usv_da)
       inquire(file = "lastselections", exist= sel_da)
       if(content_da .and. usv_da .and. sel_da  .and.   trim(farg) .ne. 'purge') then
         call push_cmd_line("inplus last_datreat_content & lastusv & sel restore & sel add fit+ & dsl ")
! ===>                                                   ! this is the default command separator
!                                                        ! as used in new_com
       endif
      

!!   Anfang der Hauptschleife
 2000  continue
       call sig_Reset()
       if(ierrs.ne.0) then
         write(6,*)' Final: error return code:',ierrs
         ioldc = 0
!        ----------> if error forget rest of commandline (reslin)
!         write(*,*)" check how to deal with those errors when usingh new_com"
!         write(*,*)" check how to deal with those errors when usingh new_com"
!         write(*,*)" check how to deal with those errors when usingh new_com"
!!         if(inka.ne.5) then
!!            close(inka)
!!            inka = 5
!!!           --------> go back to input from keyboard
!!         endif
       endif
       ierrs = 0
!       mask_err = .false.
!

       iibuf= isels(1)
       iadda= isels(1)
       mask_err = .false.   ! expression evaluation error detection is activated
                            ! since evaluation is done in an early stage of line interpretation
                            ! this error should be masked if comment lines etc. are read
                            ! e.g. in input in order to avoid erroneous error messages
                            ! when constructions that resemble expression are contained.

!       call incom(comand)
!      ------------------

       call get_newcom(comand, istat)  

       if(istat < 0) then   ! exit
 !         write(*,*)"istat = ",istat
          select case(istat)
          case (-1)
            write(*,*)'..to exit datreat enter "quit" or "exit"!' 
          case (-2)
            write(*,*)'..to exit datreat enter "quit" or "exit"!' 
          case (-3)
            call SaveStatus
!            open(19,file="lastselections")
!              write(19,'(i8)')nsel,(isels(i),i=1,nsel)            
!              write(19,'(i8)')nsel,(ifits(i),i=1,nsel)            
!              write(19,'(i8)')nsel,(isfits(i),i=1,nsel)            
!            close(19)
!            open(19,file="lastusv")
!             write(19,'(a)')"makro"
!             if(nousev > 0) write(19,'("set ",a,"  ",e16.8)')(usenam(i),useval(i),i=1,nousev)
!            close(19)
!     
!            isels  = 0
!            ifits  = 0
!            isfits = 0 
!            nsel = nbuf
!            isels(1:nsel) = [(i,i=1,nsel)]
!            fsname = "last_datreat_content"
!            call msavdat(fsname)
!            write(*,*) "QUIT: content saved to last_datreat_content ... "
            write(*,*) " exit datreat BYE.... "
            stop
          case (-4)
            write(*,*) "EXIT: content not saved, state of previous section is still valid "
            write(*,*) " exit datreat BYE.... "
            stop
          case default
             write(*,*)"PROGRAMMING ERROR (new_com): unknown istatus ", istat
          end select 
          goto 2000
       endif
!       iout() = iot
!
!
       if(comand.eq.'?       '.or.comand.eq.'help    ') then
!                    -                       ----
!         call system('firefox '//datreat_path//'/doc/DatreatManual6.html &')
         call command_list()
         write(*,*)"Please check datreat../doc/DatreatManual6.html for more help"
         goto 2000
       endif
!
!
       if(comand.eq.'prompt  ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= prompt <name>                                                               '
           write(6,*)'=    replaeces chars 2:9 of prompt by the given <name>                        '
           write(6,*)'=    thus a datreat session can be given an easy identifier if several are run'
           write(6,*)'=============================================================================='
           goto 2000
        endif
                   
         if(inames .ne. 1) then
           call errsig(999,"ERROR: prompt modification requires exactly 1 name$")
         else
           prompt = "# "//trim(vname(1))//":"//trim(prompt(len_trim(prompt)-7:))//" "
         endif
         goto 2000
       endif

!
       if(comand.eq.'freeze  ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= freeze                                                                      '
           write(6,*)'=    freezes current state (as quit would do, but stays active)               '
           write(6,*)'=    if THIS state shall be kept leave datreat by exit                        '
           write(6,*)'=============================================================================='
           goto 2000
        endif
                   
         call SaveStatus

         goto 2000
       endif

!
       if(comand.eq.'prompt  ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= prompt <name>                                                               '
           write(6,*)'=    replaeces chars 2:9 of prompt by the given <name>                        '
           write(6,*)'=    thus a datreat session can be given an easy identifier if several are run'
           write(6,*)'=============================================================================='
           goto 2000
        endif
                   
         if(inames .ne. 1) then
           call errsig(999,"ERROR: prompt modification requires exactly 1 name$")
         else
           prompt = "# "//trim(vname(1))//":"//trim(prompt(len_trim(prompt)-7:))//" "
         endif
         goto 2000
       endif

       if(comand.eq.'editor  ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= editor <editor>                                                             '
           write(6,*)'=    sets the system command to invoke the default deitor to be use           '
           write(6,*)'=    on mac: open -W -a <mypreferrededitor>                                   '
           write(6,*)'=============================================================================='
           write(6,*)'= present setting for editor cmd: ',trim(editor)
           goto 2000
        endif
    
        i = index(inline,"editor ")
        if(i>0) then
          editor = inline(8:87)
          write(*,'(a)')"Editor now will be invoked by: ",trim(editor)
        endif                
   
         goto 2000
       endif

       if(comand.eq.'ed      ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= ed <file>                                                                   '
           write(6,*)'=    edit file                                                                '
           write(6,*)'=============================================================================='
           goto 2000
        endif
        i = index(inline,"ed ")
        if(i>0) then
          if(iargvs == 0) then
            call execute_command_line(trim(editor)//" "//trim(inline(4:)))
          else
            call execute_command_line(trim(editor)//" "//trim(argvals(1)))       
          endif
        endif                
   
         goto 2000
       endif

       if(comand.eq.'title   ') then
!                    ------ 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= title "some string"                                                         '
           write(6,*)'=    other than tit with evaluation                                           '
           write(6,*)'=============================================================================='
           goto 2000
        endif
                   
 !        if(inames < 1) then
 !          call errsig(999,"ERROR: title modification requires  1 (long)name$")
 !        else
          title = " "
          do i=1,iargvs
           title = trim(title)//" "// argvals(i)(1:len(title))
          enddo
 !        endif
         write(6,'(a,a)')"current plot title: ",trim(title)
         goto 2000
       endif

       if(comand.eq.'in      '.or.comand.eq.'input   ') then
!                    --                      -----
         mask_err = .true.
         call input
 !?        nsel = 1
 !?        isels(1) = nbuf
 !?        ifits(1) = 0
         ifits = 0
         mask_err = .false.
         goto 2000
       endif
!
       if(comand.eq.'out_gli '.or.comand.eq.'gli     ') then
!                    --                      -----
         call out_gli
         goto 2000
       endif
!
       if(comand.eq.'inscn   ') then
!                    -----
         call inscn
         nsel = 1
         isels(1) = nbuf
         ifits(1) = 0
         goto 2000
       endif
!
!
       if(comand.eq.'cs      '.or.comand.eq.'clearsel') then
!                    --                      --------
! 2020    continue
!   <------ dieser sprung ist unelegant und nicht systemkonform
!           in der naechsten version subroutine cs verwenden !!
         nsel = 0
         nfits = 0
         isels = 0
         ifits = 0
         isfits = 0

         write(6,*)'selections are removed ....'
         goto 2000
       endif
!
!
       if(comand.eq.'z       '.or.comand.eq.'zero    ') then
!                    -                       ----
         nbuf = 0
         write(6,*)'nbuf has been resetted ....'
         nsel  = 0
         isels = 0
         ifits = 0
         isfits = 0
         write(6,*)'selections are removed ....'
       endif
!
!
       if(comand.eq.'noise   ') then
!                    -----
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= noise                                                                       '
           write(6,*)'=    add gaussian noise to the selected data record                           '
           write(6,*)'=    parameters:                                                              '
           write(6,*)'=      a   <val>        : factor scaling ydata to original counts for errdet  '
           write(6,*)'=      iseed  <val>     : seed for random number generator (optional)         '
           write(6,*)'=      errors           : option, if there: only estimate errors from y-values'
           write(6,*)'=                         if not: add noise to y-valuse and set errors        '
           write(6,*)'=   RESULT: a copy on a new record is generated     xs                        '
           write(6,*)'=============================================================================='
           goto 2000
        endif

         if(nsel <= 0) then
           write(6,*)'no records selected !'
           call errsig(999,"ERROR: noise no records selected $")
           goto 2000
         endif

         if(nbuf+nsel > ubound(ywerte,dim=2)) then
           write(6,*)'too many records selected !'
           call errsig(999,"ERROR: noise too many records selected $")
           goto 2000
         endif

         ampli = getval('a       ',dble(ampli),inew)
         iseed = intval('iseed   ',iseed,inew)
         if(inew.ne.0) then
            call rnset(iseed)
            write(6,*)'new seed=',iseed
         endif

       
 
         do j=1,nsel
            iadd = isels(j)
            nbuf = nbuf + 1
            if(nbuf.gt.size(nwert)) nbuf = size(nwert)
            ides = nbuf
            call txfera(iadd,ides)
            nn   = nwert(iadd)
            if(found('errors  ')) then
              do i=1,nn
                yerror(i,ides) = sqrt(ampli*ywerte(i,ides))/ampli
              enddo
            else
              do i=1,nn
               ywerte(i,ides) = igrand(ywerte(i,iadd)*ampli)
               yerror(i,ides) = sqrt(ywerte(i,ides))
              enddo
            endif
            xh = 1
            call     parget('noise_a ',xh    ,iadd, ier)
            call     parset('noise_a ',ampli*xh ,ides)
            if(xh .ne. 1) write(*,*)"WARNING: noise applied more than once!"
         enddo
 
         if(found('errors  ')) then
           write(*,*) "y-values have been scaled with factor a (noise_a) = ",ampli
           write(*,*) "         in order to have counts equivalent. "
           write(*,*) "y-errors then are set to the expecation of count statistics"
         else
           write(*,*) "y-values have been scaled with factor a (noise_a) = ",ampli
           write(*,*) "  and according GAUSSIAN NOISE HAS BEEN ADDED."
           write(*,*) "         in order to have counts equivalent. "
           write(*,*) "y-errors then are set to the expecation of count statistics"
         endif

         write(*,'(a)',advance='no')"Originally selected:"
         write(*,'(40(1x,i0))') isels(1:nsel)
   
         isels(1:nsel) = [(i,i=nbuf-nsel+1,nbuf)]  ! select the newly created records
 
         write(*,'(a)',advance='no')"Now the results are selected:"
         write(*,'(40(1x,i0))') isels(1:nsel)
         
         goto 2000
       endif
!
!
       if(comand.eq.'gen.res ') then
!                    -------
! ---> create a file in in11 spin-echo result format
         if(nsel.eq.0) then
           write(6,*)'no item selected'
           goto 2000
         endif
         iadd = isels(1)
         nn   = nwert(iadd)
         do i=1,8
           if(name(iadd)(i:i).eq.' ') name(iadd)(i:i)='x'
         enddo
         open(33,file=name(iadd)//'.res',status='UNKNOWN')
         ynno = ywerte(1,iadd)
         call parget('theta   ',theta2,iadd,ier)
         call parget('alam    ',alam  ,iadd,ier)
         call parget('bkgr    ',bkgr  ,iadd,ier)
         alam = alam*1.0e10
         write(33,'(a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')                &
     &   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',               &
     &   bkgr
         write(6,*)'----------------------------------'
         write(6 ,'(a,a8,1x,a,f6.2,1x,a,f6.2,1x,a,f12.5)')              &
     &   ' file=',                                                      &
     &   name(iadd),'2tht=',2*theta2,'lam=',alam,'bkgr=',               &
     &   bkgr
         write(6,'(1x,a)')coment(iadd)
         do i=1,nn
          write(33,*)xwerte(i,iadd),ywerte(i,iadd)/ynno,                &
     &                              yerror(i,iadd)/ynno
         enddo
         close(33)
         goto 2000
       endif
!
!

       if(comand.eq.'echocurv') then
!                    ---------
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= echocurv j <j> ......                                                       '
           write(6,*)'=    computation of simulated phase scan of NSE for a simple Lorenzian spect. '
           write(6,*)'=    ft by integration assuming gaussian distribution functions               '
           write(6,*)'=    parameters:                                                              '
           write(6,*)'=      j   <val>        : field integral in Gauss m                           '
           write(6,*)'=      dj  <val>        : asymmetry (phase0) in Gauss m                       '
           write(6,*)'=      j0delta <val>    : inhomogeneity offset (external fields..) in Gm      '
           write(6,*)'=      ddj <val>        : phase scan step (field integral) in Gm              '
           write(6,*)'=      n <val>          : number of points in scan                            '
           write(6,*)'=      cdelta <val>     : relative inhomogeneity                              '
           write(6,*)'=      lambda0 <val>    : wavelength     in Angstroem                         '
           write(6,*)'=      dlambda <val>    : wavelength width                                    '
           write(6,*)'=      tau <val>        : relaxation time --> spectrum  in ns                 '
           write(6,*)'=      errabs <> errel <> maxint  : integration parameters                    '
           write(6,*)'=============================================================================='
           goto 2000
        endif

         if(nbuf.lt.size(nwert)) then
            ia1  =  nbuf+1
            nbuf =  ia1
         else
            call errsig(100,'buffer full!$')
            goto 2000
         endif

         j1echo = getval('j       ',dble(j1echo)    ,inew)
         dj     = getval('dj      ',dble(dj )       ,inew)
         j0delta= getval('j0delta ',dble(j0delta)   ,inew)
         ddj    = getval('ddj     ',dble(ddj)       ,inew)
         nech   = intval('n       ',nech  ,inew)
         cdelta = getval('cdelta  ',dble(cdelta)    ,inew)
         alam0  = getval('lambda0 ',dble(alam0)     ,inew)
         dalam  = getval('dlambda ',dble(dalam)     ,inew)
         tau    = getval('tau     ',dble(tau  )     ,inew) * unit_ns
         erraec = getval('errabs  ',dble(erraec)    ,inew)
         errrec = getval('errrel  ',dble(errrec)    ,inew)
         maxint = intval('maxint  ',maxint,inew)

         call     parset('j       ',j1echo     ,ia1)
         call     parset('dj      ',dj         ,ia1)
         call     parset('j0delta ',j0delta    ,ia1)
         call     parset('ddj     ',ddj        ,ia1)
         call     parset('cdelta  ',cdelta     ,ia1)
         call     parset('lambda0 ',alam0      ,ia1)
         call     parset('dlambda ',dalam      ,ia1)
         call     parset('tau     ',tau        ,ia1)
         call     parset('errabs  ',erraec     ,ia1)
         call     parset('errrel  ',errrec     ,ia1)


! Fouriertime
! J \lambda^3 \, \frac{\gamma_n m_n^2}{2\pi h^2}}
    
 
         fouriertime = j1echo*unit_Gauss * (alam0 * unit_Angstroem)**3 *  Larmorkonstante * Neutronenmasse**2 / &
                       ( Planckkonstante**2 ) /  unit_ns 


         alim(1) = alam0 - 3*dalam
         blim(1) = alam0 + 3*dalam

         alim(2) = -5.d0/tau
         blim(2) = -alim(2)

         alim(1)= getval('a1      ',dble(alim(1)),inew)
         alim(2)= getval('a2      ',dble(alim(2)),inew)
         blim(1)= getval('b1      ',dble(blim(1)),inew)
         blim(2)= getval('b2      ',dble(blim(2)),inew)


         write(6,'(a,f18.9," ns        ")')      '-> fouriertime =  ',fouriertime
         write(6,'(a,f18.9," Gauss*m   ")')      'j                 ',j1echo
         write(6,'(a,f18.9," Gauss*m   ")')      'dj                ',dj
         write(6,'(a,f18.9," Gauss*m   ")')      'ddj               ',ddj
         write(6,'(a,f18.9," Gauss*m   ")')      'j0delta           ',j0delta
         write(6,'(a,f18.9,"           ")')      'cdelta            ',cdelta
         write(6,'(a,f18.9," Angstroem ")')      'lambda0           ',alam0
         write(6,'(a,f18.9," Angstroem ")')      'dlambda           ',dalam
         write(6,'(a,f18.9," ns        ")')      'tau               ',tau / unit_ns
         write(6,'(a,f18.9,"           ")')      'errabs            ',erraec
         write(6,'(a,f18.9,"           ")')      'errrel            ',errrec
         write(6,'(a,i18  ,"           ")')      'maxint            ',maxint
         write(6,'(a,f18.9," Angstroem ")')      'a1                ',alim(1)
         write(6,'(a,f18.9," 1/ns      ")')      'a2                ',alim(2) * unit_ns
         write(6,'(a,f18.9," Angstroem ")')      'b1                ',blim(1)
         write(6,'(a,f18.9," 1/ns      ")')      'b2                ',blim(2) * unit_ns

         j2echo = j1echo + dj

         call  qand(f,2,alim,blim,erraec,errrec,maxint,                 &
     &              result, errest )


         write(6,*)'-------------------------------------------'
         write(6,*)'result =',result
         write(6,*)'errest =',errest
         write(6,*)'-------------------------------------------'
         n0 =     nech/2
         n0     = intval('n0      ',n0    ,inew)
                                                       !! aix.sp
         call     parset('n0      ',(real(n0)) ,ia1)
         call     parset('fouriert',sngl(fouriertime),ia1)

         nwert(ia1) = nech
         do i=1,nech
           j2echo = j1echo + (i-n0) * ddj + dj
           xwerte(i,ia1) = j2echo-j1echo
! --- imsl mehrdimensionale integration ---
           call  qand(f,2,alim,blim,erraec,errrec,maxint,               &
     &                result, errest )
!           write(6,*)'result =',result,'     errest=',errest
           ywerte(i,ia1) = result
           yerror(i,ia1) = 0
         enddo
         xname(ia1) = 'dj/gm'
         yname(ia1) = 'counts'
         name(ia1)  = 'echo'
         write(coment(ia1),'(a,f18.9, "ns")')"echo ft=", fouriertime
         numor(ia1) = intval('numor   ',1,inew)
         write(6,*)'ok'

         nsel = 1
         isels(1) = ia1

         goto 2000
       endif


       if(comand.eq.'fft     ') then
!                    -------
         ia1    = isels(1)
         call parget ('n0      ',an0     ,ia1 ,ier)
         if(ier.eq.0) then
           n0 = Nint(an0)
         else
           n0 = nwert(ia1)/2
         endif
         mpk = intval('mpeak   ',mpk,inew)
         famp= getval('famp    ',dble(famp),inew)
         fwid= getval('fwid    ',dble(fwid),inew)
         write(6,*)'n0 ....... =',n0
         write(6,*)'mpeak..... =',mpk
         write(6,*)'famp ..... =',famp
         write(6,*)'fwid ..... =',fwid
         if(nbuf.le.size(nwert)-5) then
           ia2 = nbuf+1
           ia3 = nbuf+2
           ia4 = nbuf+3
           ia5 = nbuf+4
           ia6 = nbuf+5
           nbuf= ia6
         else
           call errsig(100,'buffer full!$')
           goto 2000
         endif

         write(6,*)'fft    ',ia1,' to  ',ia2,' phase -->',ia3
         write(6,*)'filter ',ia4,' filtered: ',ia5
         call txfera(ia1,ia2)
         call txfera(ia1,ia3)
         call txfera(ia1,ia4)
         call txfera(ia1,ia5)
         call txfera(ia1,ia6)
         n1 = nwert(ia1)
         dt = xwerte(2,ia1)-xwerte(1,ia1)
         do i=1,n1
          ca(i)    =  cmplx( ywerte(i,ia1)*dt, 0.0 )
         enddo
         call  fftcf( n1, ca, cb )
         tt = xwerte(n1,ia1) - xwerte(1,ia1)
         df = 6.283185307/tt
         do i=1,n1
          xwerte(i,ia2) = (i-1)*df
          xwerte(i,ia3) = xwerte(i,ia2)
          ywerte(i,ia2) = sqrt(real(cb(i)*conjg(cb(i))))
          ywerte(i,ia3) = atan2(real(cb(i)),aimag(cb(i)))
         enddo
! --- und filtern ---
         do i=1,n1/2
           ywerte(i,ia4) = 2*atan(famp*exp(-((i-1-n1/4)/fwid)**2))/pi
           xwerte(i,ia4) = xwerte(i,ia2)
           cb(i) = cb(i) * ywerte(i,ia4)
           xwerte(n1+1-i,ia4) = xwerte(n1+1-i,ia2)
           if(i.ge.2) then
             cb(n1+2-i)= cb(n1+2-i) * ywerte(i,ia4)
             ywerte(n1+2-i,ia4) = ywerte(i,ia4)
           endif
         enddo
         ywerte(n1/2+1,ia4) = 0
         cb(n1/2+1) = (0.0,0.0)
         call  fftcb( n1, cb, ca )
         ymaxf = 0
         yminf = 1.e30
         do i=1,n1
          ywerte(i,ia5) = real(ca(i))/(dt*n1)
          xwerte(i,ia5) = xwerte(i,ia1)
          ywerte(i,ia6) = aimag(ca(i))/(dt*n1)
          xwerte(i,ia6) = xwerte(i,ia1)
          if(ymaxf.lt.ywerte(i,ia5)) then
             ymaxf = ywerte(i,ia5)
             imaxf = i
             ymm2  = ywerte(i,ia5)-ywerte(i+2,ia5)
             ym1   = ywerte(i+1,ia5)
          endif
          if(yminf.gt.ywerte(i,ia5)) then
             yminf = ywerte(i,ia5)
             iminf = i
          endif
         enddo
         yampf = sqrt( (ymm2/2.0)**2 + ym1**2)

         name(ia2) = 'fftabs '
         name(ia3) = 'fftphas'
         name(ia4) = 'filter '
         name(ia5) = 'filtrd.r'
         name(ia6) = 'filtrd.i'
         xname(ia2)= '_1/dj*gm'
         xname(ia3)= '_1/dj*gm'
         xname(ia4)= '_1/dj*gm'
         xname(ia5)= '   dj/gm'
         xname(ia6)= '   dj/gm'
         numor(ia2)= numor(ia1)+1000
         numor(ia3)= numor(ia1)+2000
         numor(ia4)= numor(ia1)+3000
         numor(ia5)= numor(ia1)+4000
         numor(ia6)= numor(ia1)+5000
         nwert(ia4)= n1
         nwert(ia5)= n1
         nwert(ia6)= n1
! --- speziell fuer echo-auswertung ----
         nnn = nwert(ia2)
         nn1 = 3
         nn2 = nnn/4-mpk
         nn3 = nnn/4+mpk
         nn4 = nnn/2
         sumb= 0
         nb  = 0
         sump= 0
         npk = 0
         do i=nn1,nn2
          sumb = sumb + ywerte(i,ia2)**2
          nb   = nb+1
         enddo
         do i=nn3,nn4
          sumb = sumb + ywerte(i,ia2)**2
          nb   = nb+1
         enddo
         sumb = sumb/nb
         do i=nn2,nn3
           yyy  = ywerte(i,ia2)**2 - sumb
           if(yyy.ge.0.0) then
              sump = sump + sqrt(yyy)
           else
              sump = sump - sqrt(-yyy)
           endif
           npk = npk+1
         enddo
         sumb  = sqrt(sumb)
         sump  = sump*df/pi
         write(6,*)'peak ---> ',sump
         write(6,*)'noise --> ',sumb
         call parset('peak    ',sump,ia2)
         call parset('mpk     ',float(mpk),ia2)
         call parset('noisef  ',sumb,ia2)
! --- echoauswertung konventionelle methode ---
         y1 = ywerte(n0-2,ia1)
         y2 = ywerte(n0-1,ia1)
         y3 = ywerte(n0  ,ia1)
         aver = (y1+y3)*0.5
         echo = sqrt( ((y3-y1)*0.5)**2 + (y2-aver)**2 )
         write(6,*)'three point echo = ',echo
         call parset('echo3   ',echo,ia2)
         write(6,*)'filter maximum   = ',ymaxf,'   (min=',yminf,')'
         write(6,*)'filter amplitude = ',yampf,'   (min=',yminf,')'

         goto 2000
       endif


       if(comand.eq.'fftrand ') then
!                    -------
         ia1    = isels(1)
         call parget ('n0      ',an0     ,ia1 ,ier)
         if(ier.eq.0) then
           n0 = an0+0.1
         else
           n0 = nwert(ia1)/2
         endif
         if(nbuf.le.size(nwert)-5) then
           ia2 = nbuf+1
           ia3 = nbuf+2
           ia4 = nbuf+3
           ia5 = nbuf+4
           ia6 = nbuf+5
           nbuf= ia6
         else
           call errsig(100,'buffer full!$')
           goto 2000
         endif
         mpk = intval('mpeak   ',mpk,inew)
         amp = getval('a       ',dble(amp),inew)
         nra = intval('nrand   ',nra,inew)
         famp= getval('famp    ',dble(famp),inew)
         fwid= getval('fwid    ',dble(fwid),inew)

         write(6,*)'fft  ',ia1,' to  ',ia2,' phase -->',ia3
         write(6,*)'filter ',ia4,' filtered: ',ia5
         write(6,*)'(a)mplitude=',amp
         write(6,*)'nrand      =',nra
         write(6,*)'mpeak      =',mpk
         write(6,*)'n0         =',n0
         write(6,*)'famp       =',famp
         write(6,*)'fwid       =',fwid
         call txfera(ia1,ia2)
         call txfera(ia1,ia3)
         call txfera(ia1,ia4)
         call txfera(ia1,ia5)
         call txfera(ia1,ia6)
         n1 = nwert(ia1)
         dt = xwerte(2,ia1)-xwerte(1,ia1)
! --- speziell fuer echo-auswertung ----
         nnn = nwert(ia2)
         nn1 = 3
         nn2 = nnn/4-mpk
         nn3 = nnn/4+mpk
         nn4 = nnn/2
         do i=1,n1
          ywerte(i,ia3) = 0
          yerror(i,ia3) = 0
          ywerte(i,ia6) = 0
          yerror(i,ia6) = 0
         enddo
! ---> fuer den filter
         do i=1,n1/2
           ywerte(i,ia4) = 2*atan(famp*exp(-((i-1-n1/4)/fwid)**2))/pi
           if(i.ge.2) then
             ywerte(n1+2-i,ia4) = ywerte(i,ia4)
           endif
         enddo
         ywerte(n1/2+1,ia4) = 0
! <---
         bkave = 0
         bkqave= 0
         pkave = 0
         pkqave= 0
         p3ave = 0
         p3qave= 0
         pfave = 0
         pfqave= 0
         ximaxf = 0
         ximaxq = 0

         do ira=1,nra
           do i=1,n1
             ywr = igrand(ywerte(i,ia1)*amp)*dt
             ca(i)    =  cmplx( ywr, 0.0 )
           enddo
           call  fftcf( n1, ca, cb )
           tt = xwerte(n1,ia1) - xwerte(1,ia1)
           df = 6.283185307/tt
           do i=1,n1
            ywerte(i,ia2) = sqrt(real(cb(i)*conjg(cb(i))))
            ywerte(i,ia3) = ywerte(i,ia3) + ywerte(i,ia2)
            yerror(i,ia3) = yerror(i,ia3) + ywerte(i,ia2)**2
           enddo
! --- und filtern ---
         do i=1,n1
           cb(i) = cb(i) * ywerte(i,ia4)
         enddo
         call  fftcb( n1, cb, ca )
         ymaxf = 0
         yminf = 1.e30
         do i=1,n1
          ywerte(i,ia5) = ca(i)/(dt*n1)
          ywerte(i,ia6) = ywerte(i,ia6) + ywerte(i,ia5)
          yerror(i,ia6) = yerror(i,ia6) + ywerte(i,ia5)**2
          if(ymaxf.lt.ywerte(i,ia5)) then
             ymaxf = ywerte(i,ia5)
             imaxf = i
          endif
          if(yminf.gt.ywerte(i,ia5)) then
             yminf = ywerte(i,ia5)
             iminf = i
          endif
         enddo
         imaxf = intval('imaxf   ',imaxf,inew)
         ximaxf = ximaxf + imaxf
         ximaxq = ximaxq + imaxf**2
         yfamp =                                                        &
     &         ywerte(imaxf,ia5)**2 + ywerte(imaxf+1,ia5)**2            &
     &    +    ywerte(imaxf,ia5)**2 + ywerte(imaxf-1,ia5)**2
         yfamp = sqrt( yfamp/2.0)
!        ------> filtern: bis hier
! ---- echoauswertung ----
           sumb= 0
           nb  = 0
           sump= 0
           npk = 0
           do i=nn1,nn2
            sumb = sumb + ywerte(i,ia2)
            nb   = nb+1
           enddo
           do i=nn3,nn4
            sumb = sumb + ywerte(i,ia2)
            nb   = nb+1
           enddo
           sumb = sumb/nb
           do i=nn2,nn3
             yyy  = ywerte(i,ia2)**2 - sumb
             if(yyy.ge.0.0) then
                sump = sump + sqrt(yyy)
             else
                sump = sump - sqrt(-yyy)
             endif
             npk = npk+1
           enddo
           sump  = sump*df/pi
! --- 3-pt-konventionell ---
           y1 = igrand(ywerte(n0-2,ia1)*amp*n1/3.0)
           y2 = igrand(ywerte(n0-1,ia1)*amp*n1/3.0)
           y3 = igrand(ywerte(n0,ia1)*amp*n1/3.0 )
           aver = (y1+y3)*0.5
           echo = sqrt( ((y3-y1)*0.5)**2 + (y2-aver)**2 )
           if(iout().gt.0)write(6,*)ira,'  peak:',sump,'  noise:',sumb,   &
     &                            ' peakf:',yfamp
!
           bkave = bkave + sumb
           bkqave= bkqave+ sumb**2
           pkave = pkave + sump
           pkqave= pkqave+ sump**2
           p3ave = p3ave + echo
           p3qave= p3qave+ echo**2
           pfave = pfave + yfamp
           pfqave= pfqave+ yfamp**2
!
         enddo
         bkave = bkave/nra
         pkave = pkave/nra
         p3ave = p3ave/nra
         pfave = pfave/nra
         bkqave= bkqave/nra
         pkqave= pkqave/nra
         p3qave= p3qave/nra
         pfqave= pfqave/nra
         ximaxf = ximaxf/nra
         ximaxq = ximaxq/nra
         bkqave= sqrt(bkqave - bkave**2)
         pkqave= sqrt(pkqave - pkave**2)
         p3qave= sqrt(p3qave - p3ave**2)*3.0/n1
         p3ave = p3ave *3.0/n1
         pfqave= sqrt(pfqave - pfave**2)
         ximaxq= sqrt(ximaxq-ximaxf**2)
         write(6,*)'noise average=',bkave,'  deviation=',bkqave
         write(6,*)'peak  average=',pkave,'  deviation=',pkqave
         write(6,*)'peak3 average=',p3ave,'  deviation=',p3qave
         write(6,*)'peakf average=',pfave,'  deviation=',pfqave
         write(6,*)'peakf posit. =',ximaxf,'  deviation=',ximaxq
         call parset('peak    ',pkave,ia3)
         call parset('peak3   ',p3ave,ia3)
         call parset('peakdev ',pkqave,ia3)
         call parset('peak3dev',p3qave,ia3)
         call parset('mpk     ',float(mpk),ia3)
         call parset('noisef  ',bkave,ia3)
         call parset('noisedev',bqkave,ia3)
         call parset('peakf   ',pfave,ia5)
         call parset('peakfdev',pfkave,ia5)
         call parset('famp    ',pfave,ia5)
         call parset('fwid    ',bfkave,ia5)
         call parset('famp    ',bfave,ia4)
         call parset('fwid    ',pfkave,ia4)
         call parset('peakf   ',pfave,ia6)
         call parset('peakfdev',pfkave,ia6)
         call parset('famp    ',pfave,ia6)
         call parset('fwid    ',pfkave,ia6)


         do i=1,n1
            xwerte(i,ia2) = (i-1)*df
            xwerte(i,ia3) = xwerte(i,ia2)
            xwerte(i,ia4) = xwerte(i,ia2)
            xwerte(i,ia5) = xwerte(i,ia1)
            xwerte(i,ia6) = xwerte(i,ia1)
            ywerte(i,ia3) = ywerte(i,ia3)/nra
            ywerte(i,ia6) = ywerte(i,ia6)/nra
            yerror(i,ia3) = yerror(i,ia3)/nra
            yerror(i,ia6) = yerror(i,ia6)/nra
            yerror(i,ia3) = sqrt(yerror(i,ia3)-ywerte(i,ia3)**2)
            yerror(i,ia6) = sqrt(yerror(i,ia6)-ywerte(i,ia6)**2)
         enddo

         name(ia2) = 'fftabs '
         name(ia3) = 'fftave '
         name(ia4) = 'filter '
         name(ia5) = 'filtered'
         name(ia6) = 'filtaver'
         xname(ia2)= '_1/dj*gm'
         xname(ia3)= '_1/dj*gm'
         xname(ia4)= '_1/dj*gm'
         xname(ia5)= '   dj/gm'
         xname(ia6)= '   dj/gm'
         numor(ia2)= numor(ia1)+1000
         numor(ia3)= numor(ia1)+2000
         numor(ia4)= numor(ia1)+3000
         numor(ia5)= numor(ia1)+4000
         numor(ia6)= numor(ia1)+5000
         nwert(ia4)= n1
         nwert(ia5)= n1
         nwert(ia6)= n1
         write(6,*)'peak ---> ',sump
         write(6,*)'noise --> ',sumb
         call parset('peak    ',sump,ia2)
         call parset('mpk     ',float(mpk),ia2)
         call parset('noisef  ',sumb,ia2)

         goto 2000
       endif

       if(comand.eq.'uni_ft   ') then
!                    -------
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= uni_ft [store_at <n>] [resnorm]  [rexpand <val>]  [reslim <val>]            '
           write(6,*)'=      performs fourier transform of the first selected records               '
           write(6,*)'=      the second selected record (may) hold the resolution                   '
           write(6,*)'=   options:                                                                  '
           write(6,*)'=            resnorm : normalisation to resolution intensity                  '
           write(6,*)'=            rexpand : multiplies (instead of adding)                         '
           write(6,*)'=            reslim:   limit for resolution at deconvolution                  '
           write(6,*)'=   spectrum x-axis unit given by valid names: micro-eV, meV, GHz, omega      '
           write(6,*)'=   intermediate results (Ft(data), Ft(res), ..) are stored on etra records   '
           write(6,*)'=============================================================================='
           goto 2000
        endif

         if(nsel.le.0) then
           write(6,*)'error: no data records selected: no action !'
           write(6,*)
           write(6,*)'select the data record and the resolution record '
           write(6,*)'using the sel command; sequence of selcted items matters! '
           write(6,*)'if only one record is selected, only simple FT without '
           write(6,*)'deconvolution is performed'
           ierrs = 1
           goto 2000
         endif 
  
         if(nsel.gt.2) then
           write(6,*)'error: too many data records selected: no action !'
           write(6,*)
           write(6,*)'select only the data record and the resolution record '
           write(6,*)' ... for the rest its not clear what to do... '
           write(6,*)' ... in order to process several q-values reiterate '
           write(6,*)' ... directly or by using a makro '
           ierrs = 2
           goto 2000
         endif   

         ia1 = nbuf+1
         ia1 = intval('store_at',ia1,inew)

         if(found('resnorm  ') ) then
            write(6,*)"Option:   resnorm    ==> use resolution for intensity normalisation also"
         else
            write(6,*)"NO option: resnorm  specified ==> use resolution for shape deconvolution ONLY"
         endif

         unift_range_expand     = getval('rexpand ',unift_range_expand,ier)
         unift_resolution_limit = getval('reslim  ',unift_resolution_limit,ier)

         
         write(6,*)'uni_ft: selected records and store data beginning from recordnr.:',ia1
         if (nsel.eq.1) then 
           isels(2) = 0
         endif
         call uni_ft(isels(1), isels(2), ia1, found('resnorm  '), unift_range_expand, unift_resolution_limit )

         goto 2000
       endif

!
       if(comand.eq.'aritold ') then ! alte version von arit 
!                    ----
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= arit f1 <factor1> f2 <factor2> [to <numor>]  [options]                      '
           write(6,*)'=      adds the two selected records with factors 1,2 and stores as numor     '
           write(6,*)'=      for nonmatchig x-values interpolation is applied                       '
           write(6,*)'=   options:                                                                  '
           write(6,*)'=            div     : divides (instead of adding)                            '
           write(6,*)'=            mult    : multiplies (instead of adding)                         '
           write(6,*)'=            sc <numor1> <numor2>  : selection by numor match instead of sel  '
           write(6,*)'=   old version                                                               '
           write(6,*)'=============================================================================='
           goto 2000
        endif
         newnum = 777777
         idimux = 0

         num1 = 0
         num2 = 0
         do 4711 i=1,inames
          j = inapa(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'div     ') idimux = 1
          if(vname(i).eq.'mult    ') idimux = 2
          if(vname(i).eq.'to      ') newnum = Nint(rpar(j))
          if(vname(i).eq.'sc      ') then
            num1 = Nint(rpar(j))
            num2 = Nint(rpar(j+1))
          endif
          if(vname(i).eq.'factor1 '.or.vname(i).eq.'f1      ')          &
     &                                                 facto1 = rpar(j)
          if(vname(i).eq.'factor2 '.or.vname(i).eq.'f2      ')          &
     &                                                 facto2 = rpar(j)
 4711    continue
! --- figure out the adresses ---
         if(nsel.eq.2) then
           iad1 = isels(1)
           iad2 = isels(2)
         else
           iad1 = 0
           iad2 = 0
         endif
da1:     do i=1,nbuf
          if(numor(i) == 0 ) cycle da1
          if(numor(i).eq.num1) iad1 = i
          if(numor(i).eq.num2) iad2 = i
         enddo  da1
         if(iad1.eq.0) then
           write(6,*)'file :',num1,' not found'
           goto 2000
         endif
         if(iad2.eq.0) then
           write(6,*)'file :',num2,' not found'
           goto 2000
         endif
! --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4715 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
 4715    continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.size(nwert)) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
! ---- normalize to monitor-values ----
         if(withmo) then
           write(6,*)'normalizing to monitor-values ....'
           write(6,*)'use option nonorm to switch off, norm to swit. on'
! ---- extract monitorvalues ----
            xmon1 = 0
           do 4717 i=1,nopar(iad1)
             if(napar(i,iad1).eq.'monitor  ') xmon1 = params(i,iad1)
 4717      continue
           xmon2 = 0
           do 4719 i=1,nopar(iad2)
             if(napar(i,iad2).eq.'monitor  ') xmon2 = params(i,iad2)
 4719      continue
           if(xmon1.eq.0) then
             write(6,*)'monitor of first file is lacking take 1'
             xmon1 = 1
           endif
           if(xmon2.eq.0) then
             write(6,*)'monitor of 2nd   file is lacking take 1'
             xmon2 = 1
           endif
         else
           xmon1 =1
           xmon2 =1
         endif
!        ------> of if(withmo ..
! --- start the computation ---
         n1 = nwert(iad1)
         n2 = nwert(iad2)
         n3 = 0
         do 4721 i=1,n1
          x1 = xwerte(i,iad1)
!         -----> look for appropriate values in the seconf file
          j1 = 0
          j2 = 0
          d1 = 1.e30
          d2 = 1.e30
          do 4723 j=1,n2
           x2 = xwerte(j,iad2)
           if(abs(x1-x2).lt.abs(d1).and.x1.ge.x2) then
            j1 = j
            d1 = x1-x2
           endif
           if(abs(x1-x2).lt.abs(d2).and.x1.lt.x2) then
            j2 = j
            d2 = x1-x2
           endif
 4723     continue
          if(j1.eq.0.or.j2.eq.0) goto 4721
           n3 = n3+1
           xwerte(n3,newadd) = x1
           xx = x1
           y0 = ywerte(i ,iad1)
           y00= y0
           y1 = ywerte(j1,iad2)
           y2 = ywerte(j2,iad2)
           x1 = xwerte(j1,iad2)
           x2 = xwerte(j2,iad2)
           y11= y2 + (y1-y2)/(x1-x2)*(xx-x2)
           if(idimux.eq.0 ) y  = facto1 * y0/xmon1 +                    &
     &          facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2
           if(idimux.eq.1 ) y  = (facto1 * y0/xmon1) /                  &
     &        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           if(idimux.eq.2 ) y  = (facto1 * y0/xmon1) *                  &
     &        (  facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)
           ywerte(n3,newadd) = y
           yyy=y
! --- and the errors ! ----
           y0 = yerror(i ,iad1)
           y1 = yerror(j1,iad2)
           y2 = yerror(j2,iad2)
           x1 = xwerte(j1,iad2)
           x2 = xwerte(j2,iad2)
           if(idimux.eq.0) y = (facto1 * y0/xmon1)**2 +                 &
     &          (facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)**2
           if(idimux.eq.1) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           if(idimux.eq.2) y = yyy * ((facto1 * y0/xmon1)/y00) *        &
     &          ((facto2 * (y2 + (y1-y2)/(x1-x2)*(xx-x2)) / xmon2)/y11)
           yerror(n3,newadd) = sqrt(y)
 4721    continue
         call txfpar(iad1,newadd)
         nwert(newadd) = n3
         numor(newadd) = newnum
         xname(newadd) = xname(iad1)
         yname(newadd) = yname(iad1)
          name(newadd) =  name(iad1)
         coment(newadd)= coment(iad1)(1:40)//coment(iad2)(1:40)
         nnpar=nopar(iad1)
       if(nnpar.le.mpar-2) then
         nnpar=nnpar+2
         nopar(newadd) = nnpar
         params(1,newadd) = facto1
         params(2,newadd) = facto2
         napar(1,newadd)  = name(iad1)
         napar(2,newadd)  = name(iad2)
         np=3
       endif
         do 4728 i=np,nnpar
           params(i,newadd)= params(i+1-np,iad1)
           napar(i,newadd) = napar (i+1-np,iad1)
 4728    continue
         write(6,*)n3,' points resulting from arit are stored as sc',   &
     &                  newnum
         if(idimux.eq.0)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' + ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.1)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' / ',facto2,' * ',   &
     &             name(iad2)
         if(idimux.eq.2)                                                &
     &   write(6,*)'<===',facto1,' * ',name(iad1),' * ',facto2,' * ',   &
     &             name(iad2)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
!
!
       if(comand.eq.'arit     ') then  ! deals better with interpolation at boundaries of range (vorher arit2)!
!                    ----
        call arit()

        goto 2000
       endif

!
!
       if(comand.eq.'average ') then
!                    -------
         call average_data()

         goto 2000
       endif
!
!
       if(comand.eq.'aligny ') then
!                    -------
         call align_y()

         goto 2000
       endif

       if(comand.eq.'aligna ') then
!                    -------
         call align_a()

         goto 2000
       endif
!
       if(comand.eq.'addval  ') then
!                    ------
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= addval <x> <y> <err>                                                        '
           write(6,*)'=    appends one data point to selected record                                '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
           goto 2000
        endif

         if(nsel.eq.1) then
           iad1 = isels(1)
         else
           write(6,*)'select one and only one record'
           ierrs= 1
           goto 2000
         endif
! --- start the computation ---
         
         nwert(iad1) = nwert(iad1)+1
         xwerte( nwert(iad1), iad1 ) = rpar(1)
         ywerte( nwert(iad1), iad1 ) = rpar(2)
         yerror( nwert(iad1), iad1 ) = rpar(3)

         goto 2000
       endif
!
       if(comand.eq.'integrat') then
!                    --------           approximates integral from binned data
!                                       by simple summation of data times bin-width
!                                       there may be better interpolation based schemes
!                                       but this is clear an simple                        
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= integrat x1 <x1> x2 <x2>                                                    '
           write(6,*)'=        approximates integral from binned data over the given range          '
           write(6,*)'=        by simple summation of data times bin-width                          '
           write(6,*)'=        there may be better interpolation based schemes                      '
           write(6,*)'=        but this is clear an simple                                          '
           write(6,*)'=============================================================================='
           goto 2000
        endif
         if(nsel<1) then
           call errsig(200,"ERROR: integrat slect at least one record$")
           goto 2000
         endif
! --- start the computation --- these are the integration limits
 
         x1int        = getval('x1      ',dble(x1int),inew)
         x2int        = getval('x2      ',dble(x2int),inew)
 
         do j=1,nsel
            iad1 = isels(j)         
            i1int = 1
            i2int = nwert(iad1)
            
            do i=1,nwert(iad1)-1
             if(xwerte(i,iad1).le.x1int .and.xwerte(i+1,iad1).gt.x1int) i1int=i
             if(xwerte(i,iad1).le.x2int .and.xwerte(i+1,iad1).gt.x2int) i2int=i
            enddo
       
            
            sum   = (xwerte(i1int,iad1)-x1int)*ywerte(i1int,iad1)
            sumer = ((xwerte(i1int,iad1)-x1int)*yerror(i1int,iad1))**2       
       
            sum   = sum + (x2int-xwerte(i2int,iad1))*ywerte(i2int,iad1)
            sumer = sumer + ((x2int-xwerte(i2int,iad1))*yerror(i2int,iad1))**2
            do i=i1int+1,i2int-1
              sum   = sum   + (xwerte(i+1,iad1)-xwerte(i,iad1))*ywerte(i+1,iad1)
              sumer = sumer + ((xwerte(i+1,iad1)-xwerte(i,iad1))*yerror(i+1,iad1))**2
            enddo
           
            sumer = sqrt(sumer)
       
            write(6,*)"================= Record: ",iad1,"  ===================="
            write(6,'(a,i9,1x,e13.6,a,e13.6,a)')'integral from: ',i1int,x1int,' (',xwerte(i1int,iad1),' )'
            write(6,'(a,i9,1x,e13.6,a,e13.6,a)')'           to: ',i2int,x2int,' (',xwerte(i2int,iad1),' )'
          
            write(6,*                  )'         ===>: ',sum, '+-', sumer
!            xwerte( nwert(iad1), iad1 ) = rpar(1)
!            ywerte( nwert(iad1), iad1 ) = rpar(2)
!            yerror( nwert(iad1), iad1 ) = rpar(2)
       
            call parset ('integral',sum        ,iad1 ) 
            call parset ('x1integ ',x1int      ,iad1 ) 
            call parset ('x2integ ',x2int      ,iad1 ) 
         enddo 
         goto 2000
       endif
!

!
!
       if(comand.eq.'combine ') then
!                    -------
         newnum = 777777
         do 4811 i=1,inames
          j = inapa(i)
          k = inpar(i)
          if(vname(i).eq.'norm    ') withmo = .true.
          if(vname(i).eq.'nonorm  ') withmo = .false.
          if(vname(i).eq.'raster  ') then
           x1raster = rpar(j)
           x2raster = rpar(j+1)
           nraster  = Nint(rpar(j+2))
          endif
          if(vname(i).eq.'to      ') newnum = Nint(rpar(j))
          if(vname(i).eq.'sc      ') then
             do 48181 l=1,k
               nnumi(l) = Nint(rpar(l-1+j))
48181        continue
             call search(nnumi,k)
          endif
 4811    continue
! --- if result is already there , replace it ---
         newadd = nbuf+1
         do 4815 i=1,nbuf
          if(numor(i).eq.newnum) newadd = i
 4815    continue
         if(newadd.gt.nbuf) nbuf = newadd
         if(nbuf.ge.size(nwert)) then
           write(6,*)'no storage for the result'
           nbuf = nbuf - 1
           goto 2000
         endif
! --- create the raster ----------------------------------------------
         if(nraster.gt.mwert) nraster = mwert
         if(nraster.lt.2    ) nraster = 2
         draster = (x2raster-x1raster)/(nraster-1)
         do  i=1,nraster
           xwerte(i,newadd) = x1raster + (i-1)*draster
           ywerte(i,newadd) = 0.d0
           yerror(i,newadd) = 0.d0
         enddo
         write(6,*)'raster generated from x1=',x1raster
         write(6,*)'                 to   x2=',x2raster
         write(6,*)'with a number of points of n = ',nraster
         write(6,*)'the newaddress is at ......... ',newadd
! --- loop over the selectionlist ----------------------------
         nwert(newadd) = nraster
         do isl=1,nsel
           iadd = isels(isl)
           xmonitor = 1.d0
           call parget ('monitor ',xmonitor,iadd,ier)
           n1 = nwert(newadd)
           n2 = nwert(iadd)
           do 4821 i=1,n1
            x1 = xwerte(i,newadd)
!         -----> look for appropriate values in the second file
            j1 = 0
            j2 = 0
            d1 = 1.e30
            d2 = 1.e30
            do 4823 j=1,n2
             x2 = xwerte(j,iadd)
             if(abs(x1-x2).lt.abs(d1).and.x1.ge.x2) then
              j1 = j
              d1 = x1-x2
             endif
             if(abs(x1-x2).lt.abs(d2).and.x1.lt.x2) then
              j2 = j
              d2 = x1-x2
             endif
 4823       continue
            if(j1.eq.0.or.j2.eq.0) goto 4821
            xx = x1
            y0 = ywerte(i ,newadd)
            y1 = ywerte(j1,iadd)
            y2 = ywerte(j2,iadd)
            x1 = xwerte(j1,iadd)
            x2 = xwerte(j2,iadd)
            y  = y0 + (y2 + (y1-y2)/(x1-x2)*(xx-x2))
            ywerte(i,newadd) = y
            yerror(i,newadd) = yerror(i,newadd) + xmonitor
!           ----- temorary misuse of the error-filed to collect
!                 monitor-values
 4821      continue
         enddo
! ---- normalize the data and create errors --------
         do i=1,n1
           if(yerror(i,newadd).gt.0.0) then
             yer     = sqrt(ywerte(i,newadd))/ yerror(i,newadd)
             ywerte(i,newadd) = ywerte(i,newadd) / yerror(i,newadd)
             yerror(i,newadd) = yer
           endif
         enddo
! --- transfer parameters ---
         call txfpar(isels(1),newadd)
         nwert(newadd) = n1
         numor(newadd) = newnum
         xname(newadd) = xname(isels(1))
         yname(newadd) = yname(isels(1))
         name(newadd)  =  name(isels(1))
         coment(newadd)= coment(isels(1))(1:40)//coment(isels(2))(1:40)
         nnpar=nopar(isels(1))
         call  parset ('monitor ',1.0,newadd)
         isels(1) = newadd
         ifits(1) = 0
         nsel     = 1
         goto 2000
       endif
!
!
!
!
       if(comand.eq.'m       '.or.comand.eq.'mirror  ') then
!                    -                       ------
         xk0= rpar(1)
         k1 = rpar(2)+0.0001
         k2 = rpar(3)+0.0001
         if(ipars.le.0) xk0 = 32
         if(ipars.le.1) k1  = xk0 - 20
         if(ipars.le.1) k2  = xk0 - 10
         if(nsel.eq.0) then
           write(6,*)'no scan selected ; plot one to select it'
           goto 2000
         endif
         ia = isels(1)
         if(numor(ia).ge.numpls) then
           write(6,*) numor(ia),'=numor > ',numpls,'; warning no origina&
     &                data'
         endif
         if(xname(ia).ne.'kanal   '.and.xname(ia).ne.'kanal   ') then
           write(6,*)'warning x-values are not kanal !!!'
         endif
         nw = nwert(ia)
! --- determine the center ---
         xk00 = xk0
         call detxk0(ywerte(1,ia),ywerte(1,nbuf+1),xk0,k1,k2,nw)
         if(xk0.le.0)then
            xk0 = xk00
            call detxk1(ywerte(1,ia),ywerte(1,nbuf+1),xk0,k1,k2,nw)
         endif
         write(6,*)'center channel is = ',xk0
         write(6,*)'generate the mirror-image ...'
         f1 = smirro(ywerte(1,ia),ywerte(1,nbuf+1),xk0,1,nw,nw)
         do 1702 i=1,nw
           xwerte(i,nbuf+1) = xwerte(i,ia)
 1702    continue
         nbuf       = nbuf + 1
         name(nbuf) = name(ia)
         xname(nbuf)= xname(ia)
         yname(nbuf)= 'sym-'//yname(ia)(1:4)
         numor(nbuf)=mod(numor(ia),numpls)+numpls
         coment(nbuf)=coment(ia)
         nwert(nbuf)= nw
         write(6,*)'ok'
         isels(2) = ia
         isels(1) = nbuf
         ifits(1)=0
         ifits(2)=0
         nsel     = 2
         npp=nopar(ia)
         if(npp.ge.mpar) then
           write(6,*)'xk0 cannot be stored on params because too many p'
           goto 2000
         endif
!
         call parset ('xk0     ',xk0,ia)
         nopar(nbuf)= npp+1
         do 2007 i=1,npp+1
          params(i,nbuf) = params(i,ia)
          napar (i,nbuf) = napar (i,ia)
 2007    continue
         goto 2000
       endif
!
!
!
       if(comand.eq.'ercorrc '.or.comand.eq.'ecc     ') then
!                    -------                 ---
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= ercorrc n <n> a <a> ri <ri> r0 <r0> d0 <d0> ncut <nc>                       '
           write(6,*)'=    NSE correction coil design special function                              '
           write(6,*)'=    check what this reall does (remove?)xs                                     '
           write(6,*)'=============================================================================='
           goto 2000
        endif

          if(nbuf.eq.size(nwert)) then
            write(6,*)' no space for an additional item'
            call errsig(999,"ERROR: ercorrc $")
            goto 2000
          endif
         nbuf = nbuf+1
         nsel = 1
         isels(1)= nbuf
         nfits        = 0
         numor(nbuf)  = 100*nbuf
         xname(nbuf)  = 'x'
         yname(nbuf)  = 'err'
         name(nbuf)   = 'cc_err'
         coment(nbuf)= 'cc_err'
         n            = intval('n       ',n,inew)
         aa           = getval('a       ',dble(aa),inew)
         ri           = getval('ri      ',dble(ri),inew)
         r0           = getval('r0      ',dble(r0),inew)
         d0           = getval('d0      ',dble(d0),inew)
         ncut         = intval('ncut    ',ncut,inew)

         if(n.gt.mwert) then
           write(6,*) ' n is too large'
           call errsig(999,"ERROR: ercorrc $")
           goto 2000
         endif

         nwert(nbuf)  = n

         if(ncut.gt.mcut) then
           write(6,*)' ncut is too large max=',mcut
           call errsig(999,"ERROR: ercorrc $")
           goto 2000
         endif

         a = aa !??
         write(6,*)' a ........ = ',a
         write(6,*)' ri ....... = ',ri
         write(6,*)' r0 ....... = ',r0
         write(6,*)' d0 ....... = ',d0
         write(6,*)' ncut ..... = ',ncut
         write(6,*)' n ........ = ',n


         call parset ('a       ',a  ,nbuf)
         call parset ('ri      ',ri ,nbuf)
         call parset ('r0      ',r0 ,nbuf)
         call parset ('d0      ',d0 ,nbuf)
         xnc = ncut
         call parset ('ncut    ',xnc ,nbuf)

         a2 = a/2.0
         x1 = 0.0
         xm = r0*sqrt(xnc)
         dx = xm  / n
         do i=1,ncut+1
           xi = i
           rrv(i) = r0 *sqrt( xi )
         enddo
         rrv(0) = ri
         do i=0,n
           xr = x1 + i*dx
           j  = 0
           y  = 0.0
           do j=0,ncut
            if( xr .lt. rrv(j)) goto 20022
           enddo
20022      continue
           if(j.gt.0) then
             xn = a2*(rrv(j)**2-rrv(j-1)**2)+d0*alog(rrv(j)/rrv(j-1))
             xz = a2*(xr**2    -rrv(j-1)**2)+d0*alog(xr/rrv(j-1))
             y  =  xz/xn + j-1
           endif
20032      continue

           ywerte(i,nbuf) = y
           yerror(i,nbuf) = 0
           xwerte(i,nbuf) = xr

         enddo

         goto 2000
       endif
!
!
       if(comand.eq.'sym     ') then
!                    ---->  symmetrize
! ---- build now the symmetric average on the right side ---
         if(nsel.lt.2) then
           write(6,*)'invalid isel-list, start with mirro'
           goto 2000
         endif
         ib = isels(1)
         ia = isels(2)
         nn1 = numor(ia)
         nn2 = numor(ib)
         if(nn1+numpls.ne.nn2) then
            write(6,*)'no pair present, use mirro'
            goto 2000
         endif
! --- look for xk0 ---
         do 2077 i=1,nopar(ib)
            if(napar(i,ib).eq.'xk0     ') then
              xk0 = params(i,ib)
              npp = i
            endif
 2077    continue
            if(vname(1).eq.'xk0     ') xk0=rpar(1)
         write(6,*)'xk0 = ',xk0
         nw = nwert(ia)
         kan1 = xk0 + 1
         do 2010 i=kan1,nw
           yw1 = ywerte(i,ia)
           yw2 = ywerte(i,ib)
           yy  = ( yw1 + yw2 ) / 2
           if(yw1.eq.0) yy = yw2
           if(yw2.eq.0) yy = yw1
           ywerte(i+1-kan1,ib) = yy
 2010    continue
         nwert(ib) = nw-kan1+1
         xk0            = xk0 - kan1 + 1
         params(npp,ib) = xk0
         numor(ib)=mod(numor(ib),numpls)+2*numpls
         write(6,*)'ok symm-ds has ',nwert(ib),' points'
         goto 2000
       endif
!
!
       if(comand.eq.'kz      ') then
!                    -->    kanalzusammenfassung
! ---- build now the symmetric average on the right side ---
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= kz <nk>                                                                     '
           write(6,*)'=    contraction of spectra by adding <nk> adjacent channele                  '
           write(6,*)'=    channel width paramezters etc. are update in the parameterlist           '
           write(6,*)'=    such that spectra fitting stays with the proper _xwidth                  '
           write(6,*)'=    The contratction is on place (i.e.)                                      '
           write(6,*)'=    All selected items are treated                                           '
           write(6,*)'=============================================================================='
        endif

         kz = NINT(rpar(1))
         if(kz.lt.2) goto 2000
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          ikz=kz
          if(ikz.gt.n)ikz=n
          j  = 0
          m  = n/ikz
          ir = n-ikz*m
          do l=1,m
            xsum=0.0
            ysum=0.0
            esum=0.0
            do k=1,ikz
              ik=(l-1)*ikz+k
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ikz
            ywerte(j,ia) = ysum/ikz
            yerror(j,ia) = sqrt(esum)/ikz
          enddo
         if(ir.gt.0) then
            xsum=0.0
            ysum=0.0
            esum=0.0
            do ik=m*ikz+1,n
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ir
            ywerte(j,ia) = ysum/ir
            yerror(j,ia) = sqrt(esum)/ir
         endif
         nwert(ia)=j
         call parget ('_rebin  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * ikz
         call parset ('_rebin  ',xx,ia)

         call parget ('_xwidth ',xx,ia,ier)
         if(ier.ne.0) then
           xx = xwerte(m/2+1,ia)-xwerte(m/2+1,ia)
         else
           xx = xx * ikz
         endif
         call parset ('_xwidth ',xx,ia)

         write(6,*)'channels contracted for ',ia,' new no. of ch=',j,   &
     &             ' (',ir,')'
         enddo
         goto 2000
       endif



       if(comand.eq.'rebin   ') then
!                    ----->    kanalzusammenfassung mit Kopie auf
!                              neuen Platz
! 
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rebin <nk>                                                                  '
           write(6,*)'=    contraction of spectra by adding <nk> adjacent channele                  '
           write(6,*)'=    channel width paramezters etc. are update in the parameterlist           '
           write(6,*)'=    such that spectra fitting stays with the proper _xwidth                  '
           write(6,*)'=    The contracted spectra are stored on new records.                        '
           write(6,*)'=    All selected items are treated   (similar kz)                            '
           write(6,*)'=============================================================================='
        endif
         kz = NINT(rpar(1))
         if(kz.lt.2) goto 2000
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

!     >  copy first  
         if(nsel+nbuf.gt.size(nwert)) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          ikz=kz
          if(ikz.gt.n)ikz=n
          j  = 0
          m  = n/ikz
          ir = n-ikz*m
          do l=1,m
            xsum=0.0
            ysum=0.0
            esum=0.0
            do k=1,ikz
              ik=(l-1)*ikz+k
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ikz
            ywerte(j,ia) = ysum/ikz
            yerror(j,ia) = sqrt(esum)/ikz
          enddo
         if(ir.gt.0) then
            xsum=0.0
            ysum=0.0
            esum=0.0
            do ik=m*ikz+1,n
              xsum = xsum + xwerte(ik,ia)
              ysum = ysum + ywerte(ik,ia)
              esum = esum + yerror(ik,ia)**2
            enddo
            j = j+1
            xwerte(j,ia) = xsum/ir
            ywerte(j,ia) = ysum/ir
            yerror(j,ia) = sqrt(esum)/ir
         endif
         nwert(ia)=j
         numor(ia)=numor(ia)+100000
         call parget ('_rebin  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * ikz
         call parset ('_rebin  ',xx,ia)

         call parget ('_xwidth ',xx,ia,ier)
         if(ier.ne.0) then
           xx = xwerte(m/2+1,ia)-xwerte(m/2+1,ia)
         else
           xx = xx * ikz
         endif
         call parset ('_xwidth ',xx,ia)

         

         write(6,'(a,i4,a,i4)')'channels of[',ia,'] contracted by ',ikz 
         enddo
         goto 2000
       endif
!


       if(comand.eq.'rerange ') then
!                    ----->    neuer x-range durch selektion
!                              neuen Platz
!           
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rerange <min> <max>     [y]                                                 '
           write(6,*)'=    copy of selected records such that                                       '
           write(6,*)'=    only values betwen min < x|y < max  are retained , y: yrange             '
           write(6,*)'=    without option x-range is pertained                                      '
           write(6,*)'=============================================================================='
        endif


         lower_range = rpar(1)
         upper_range = rpar(2)
         range_is_y  = found('y       ')

!     >  copy first  
         if(nsel+nbuf.gt.size(nwert)) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          j  = 0
drer1:    do ik=1,n
             if(range_is_y) then
               if(ywerte(ik,ia) < lower_range .or. ywerte(ik,ia) > upper_range) cycle drer1
             else
               if(xwerte(ik,ia) < lower_range .or. xwerte(ik,ia) > upper_range) cycle drer1
             endif
             j = j+1
             xwerte(j,ia) = xwerte(ik,ia)
             ywerte(j,ia) = ywerte(ik,ia)
             yerror(j,ia) = yerror(ik,ia)
          enddo drer1

          nwert(ia)=j
          numor(ia)=numor(ia)+10000
          if(range_is_y) write(6,'(a)',advance='no') "Y-"
          write(6,'(a,i4,a,2e14.6)')'range of[',ia,'] limited to ',lower_range, upper_range 
         enddo
         goto 2000
       endif
!




       if(comand.eq.'swapxy   ') then
!                    ----->    vertauschen von x und y
!                              
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= swapxy                                                                      '
           write(6,*)'=    role of x amd y-values are exchanged                                     '
           write(6,*)'=    NOTE: y-errors are lost                                                  '
           write(6,*)'=    operation is on place no new copies!                                     '
           write(6,*)'=============================================================================='
        endif
! 
         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do l=1,n
              xsum = xwerte(l,ia)
              ysum = ywerte(l,ia)
              xwerte(l,ia) = ysum
              ywerte(l,ia) = xsum             
              yerror(l,ia) = 0
           enddo
           cbuffer   = xname(ia)
           xname(ia) = yname(ia)
           yname(ia) = cbuffer
           write(6,'(a,i5,a)')" x and y values of ",ia," have been swapped, errors set to zero!" 
         enddo
         goto 2000
       endif
!

       if(comand.eq.'sequence   ') then
!                    ----->  sequence  
!                                       
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= sequence                                                                    '
           write(6,*)'=    replace x-values by sequence numbers 1..n                                '
           write(6,*)'=    operation is on place, previous x-values are lost                        '
           write(6,*)'=============================================================================='
        endif

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do l=1,n
             xwerte(l,ia) = l
          enddo
           xname(ia) = "#"
           write(6,'(a,i5,a)')" x values of ",ia," have been replace by sequence numbers" 
         enddo
         goto 2000
       endif
!
!
       if(comand.eq.'copy   ') then
!                    ----->  copy 
!                                       
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= copy [x1 <min> x2 <max>                                                     '
           write(6,*)'=    copy of selected records to new record buffers                           '
           write(6,*)'=    if x1, x2 are given only x-values within the interval are copied         '
           write(6,*)'=============================================================================='
        endif

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

         xcopy1 = getval('x1      ',dble(xcopy1),ier)
         xcopy2 = getval('x2      ',dble(xcopy2),ier)

         do i=1,nsel
          ia = isels(i)
          if(nbuf < size(nwert)) then
             nbuf = nbuf + 1
          else
             Write(6,*)"ERROR: too many records !"
             goto 2000
          endif

           call DataCopy(ia,nbuf)
           write(6,'(a,i5,a,i5)')" copy ", ia," to ",nbuf 
           isels(i) = nbuf 

           if(found('x1      ') .and. found('x2      ')) then
             Write(6,'(a,2e14.6)')" ...restricted to the interval: ",xcopy1, xcopy2
             nwert(nbuf) = 0
             do j=1,nwert(ia)
               if(xwerte(j,ia) >= xcopy1 .and. xwerte(j,ia) <= xcopy2) then
                  nwert(nbuf) = nwert(nbuf) + 1
                  xwerte(nwert(nbuf) ,nbuf)  = xwerte(j,ia)
                  ywerte(nwert(nbuf) ,nbuf)  = ywerte(j,ia)
                  yerror(nwert(nbuf) ,nbuf)  = yerror(j,ia)
               endif
             enddo
           endif

         enddo
         goto 2000
       endif
!
       

       if(comand.eq.'scale   ') then
!                    ----->    Skalierung der Intesnitaet          

        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= scale <factor>                                                              '
           write(6,*)'=    scales selected spectra with a factor, bookkeeping in parameter section  '
           write(6,*)'=    also for gaussian resolution parameters and errors                       '
           write(6,*)'=============================================================================='
        endif
                              
         if(ipars.ne.1) then
           write(6,*)'ERROR: need just one value as parameter..'
           ierrs = 1
           goto 2000
         endif

         y_scaling =  rpar(1)

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif

!     >  copy first  
         if(nsel+nbuf.gt.size(nwert)) then
           write(6,*)'ERROR: copying selection would exceed max records'
           goto 2000
         endif 
         do i=1,nsel
           ia = isels(i)
           ib = nbuf+1
           write(6,*)'copy record: ',ia,' to record: ',ib
           call DataCopy(ia,ib)
           isels(i) = ib
         enddo
!     <  end copy
 
         do i=1,nsel
          ia = isels(i)
          n  = nwert(ia)
          do ik=1,n
              ywerte(ik,ia) = ywerte(ik,ia) *y_scaling
              yerror(ik,ia) = yerror(ik,ia) *y_scaling
          enddo
         numor(ia)=numor(ia)+200000

         call par_gaint_scale(y_scaling,ia)


         call parget ('_scale  ',xx,ia,ier)
         if(ier.ne.0) xx = 1
         xx = xx * y_scaling
         call parset ('_scale  ',xx,ia)

         write(6,'(a,i4,a,e13.6)')'values of[',ia,'] scaled by ',y_scaling 
         enddo
         goto 2000
       endif
!
!
!
!
       if(comand.eq.'gaiscale') then
!                    -------->    Skalierung der Aufloesungsparameter
!                              
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= gaiscale <factor>                                                           '
           write(6,*)'=    scaling of the gaussian resolution parameters associated                 '
           write(6,*)'=    with the selected spectra                                                '
           write(6,*)'=============================================================================='
        endif

         if(ipars.ne.1) then
           write(6,*)'ERROR: need just one value as parameter..'
           ierrs = 1
           goto 2000
         endif

         y_scaling =  rpar(1)

         if(nsel.lt.1) then
           write(6,*)'no items selected, use sel !'
           goto 2000
         endif
 
         do i=1,nsel
           call par_gaint_scale(y_scaling,isels(i))
           write(6,'(a,i4,a,e13.6)')'ga#inten values of[',isels(i),'] scaled by ',y_scaling 
         enddo
         goto 2000
       endif
!




!
       if(comand.eq.'interpol') then
!                    ---->  interpolate
! ---- build now the symmetric average on the right side ---
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= interpol <n>                                                                '
           write(6,*)'=    creates a clone of the firts (only) item of the selection list           '
           write(6,*)'=    by interpolating y-values such that the x-range is covered by <n>        '
           write(6,*)'=    equidistant points.                                                      '
           write(6,*)'=    Errors are NOT treated (yet)                                             '
           write(6,*)'=============================================================================='
        endif

         ia = isels(1)
         ib = nbuf+1
         if(ib.gt.size(nwert)) ib=size(nwert)
         nbuf = ib
         call txfpar(ia,ib)
         nn = rpar(1)
         if(nn.le.1) nn=mwert
         dxx= (xwerte(nwert(ia),ia)-xwerte(1,ia))/(nn-1)
         do i=1,nn
           xx = xwerte(1,ia)+(i-1)*dxx
           xwerte(i,ib) = xx
           ywerte(i,ib) = yinterp(xx)
         enddo
         nwert(ib) = nn
         isels(1)  = ib
         nsel      = 1
         goto 2000
       endif
!
!
       if(comand.eq.'addsels ') then
!                    ---->  summation over selected records : will become obsolete due to name
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= addsels                                                                     '
           write(6,*)'=    creates record containing the sum of all selected (compatible)           '
           write(6,*)'=    records.      (similar recsum)                                           '
           write(6,*)'=============================================================================='
        endif
         call sumseldat
         goto 2000
       endif
!
!
       if(comand.eq.'recsum  ') then
!                    ---->  summation over selected records: the new one
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= recsum from <n1> to <n2>                                                    '
           write(6,*)'=    creates record containing the sum of all records from a continuous       '
           write(6,*)'=    address range <n1> .. <n2>      (similar addsels)                        '
           write(6,*)'=============================================================================='
        endif
         if(found('from    ')) then
           ia = intval('from    ',isels(1),inew)
           ib = intval('to      ',isels(nsel),inew)
           nsel = ib-ia+1
           do i=1,nsel
             isels(i) = ia-1+i
           enddo
           write(6,'(a,i4,a,i4)')'sum records from ',ia,' ... ',ib
         endif
         call sumseldat
         goto 2000
       endif
!
!
       if(comand.eq.'parextra') then
!                    ---> parabola extrapolation to q--> 0

        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= parextra                                                                    '
           write(6,*)'=    adds a point with x = 0 by parabola extrapolation of data starting with  '
           write(6,*)'=    finite x-values (e.g. extrapolation of SANS data to q=0)                 '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
        endif

         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
         ia    = isels(1)
         slope = (ywerte(3,ia)-ywerte(1,ia))/(xwerte(3,ia)-xwerte(1,ia))
         apar  = slope / (2.0*xwerte(2,ia))
         bpar  = ywerte(2,ia) - slope*0.5*xwerte(2,ia)
         dx    = xwerte(2,ia)-xwerte(1,ia)
         npl   = xwerte(1,ia) / dx + 1
         write(6,*)'parabola: a*x^2 + b'
         write(6,*)'a = ',apar
         write(6,*)'b = ',bpar
         write(6,*)'dx= ',dx
         write(6,*)'n+= ',npl
         nold = nwert(ia)
         if(nwert(ia)+npl.gt.mwert) then
            write(6,*)'too many points'
            ierrs = 1
            goto 2000
         endif
         nwert(ia) = nwert(ia) + npl
         yerr      = yerror(1,ia)
         do i=nwert(ia),npl+1,-1
            ywerte(i,ia) = ywerte(i-npl,ia)
            xwerte(i,ia) = xwerte(i-npl,ia)
            yerror(i,ia) = yerror(i-npl,ia)
         enddo
         do i=1,npl
           xwerte(i,ia) = (i-1)*dx
           ywerte(i,ia) = apar*xwerte(i,ia)**2 + bpar
           yerror(i,ia) = yerr
         enddo
         write(6,*)'done'
         goto 2000
        endif
!
!
       if(comand.eq.'hiqextra') then
!                    ---> extrapolation for high q .
!                         data == a*q^z
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= hiqextra <xcut>                                                             '
           write(6,*)'=    adds points from hiq x (q) extrapolation by a*q^z  up to x=xcut          '
           write(6,*)'=    only the first selected record is treated                                '
           write(6,*)'=============================================================================='
        endif
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
!        ---------------
         xcut  = rpar(1)
!        ---------------
         ia    = isels(1)
         do i=1,nwert(ia)
            if(xwerte(i,ia).gt.xcut) goto 77771
         enddo
77771    continue
         icut = i-1
         slope = (ywerte(icut,ia)-ywerte(icut-2,ia))/                   &
     &           (xwerte(icut,ia)-xwerte(icut-2,ia))
         zpar  = xwerte(icut-1,ia)*slope / ywerte(icut-1,ia)
         apar  = ywerte(icut-1,ia)*xwerte(icut-1,ia)**(-zpar)
         dx    = xwerte(icut,ia)-xwerte(icut-1,ia)
         do i=icut,mwert
           xx = xwerte(icut,ia)+(i-icut)*dx
           ywerte(i,ia) = apar*xx**zpar
           yerror(i,ia) = yerror(icut,ia)
           xwerte(i,ia) = xx
         enddo
         write(6,*)'high-q extrapolation :  a*q^z'
         write(6,*)'a =  ',apar
         write(6,*)'z =  ',zpar
         write(6,*)'dx=  ',dx
         write(6,*)'ncut=',icut
         nwert(ia) = mwert
         write(6,*)'done'
         goto 2000
        endif
!
       if(comand.eq.'qc      '.or.comand.eq.'q-conv  ') then
!                    ---> channels to q-values
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia = isels(1)
        modeqc=0
        if(xname(ia)(1:1).eq.'k'.or.xname(ia)(1:1).eq.'c')modeqc=1
        if(xname(ia).eq.'rad*1.e-6')modeqc=2
        if(modeqc.eq.0) then
           write(6,*)'data do not have the proper x-values:',xname(ia)
           ierrs=100
           goto 2000
         endif
         do 2001 i=1,nopar(ia)
           if(napar(i,ia).eq.'xk0      ') xk0   = params(i, ia)
           if(napar(i,ia).eq.'lambda   ') alam0 = params(i, ia)
           if(napar(i,ia).eq.'detdis   ') detdis= params(i, ia)
           if(napar(i,ia).eq.'bklen    ') bklen = params(i, ia)
 2001    continue
           if(nbuf.lt.size(nwert))then
           nbuf=nbuf+1
           ib=nbuf
        endif
         do 2002 i=1,nwert(ia)
           if(modeqc.eq.1) then
             xwerte(i,ib)=4*pi/alam0*sin(0.5*atan(bklen*(i-xk0)/detdis))
           endif
           if(modeqc.eq.2) then
             xwerte(i,ib)=4*pi/alam0*sin(1.e-6*xwerte(i,ia)/2)
           endif
           ywerte(i,ib)=ywerte(i,ia)
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+3*numpls
 2002    continue
         xname(ib)='q      '
         write(6,*)'q1=',xwerte(1,ib),'  qend=',xwerte(nwert(ib),ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'sp      '.or.comand.eq.'spline  ') then
!                    ---> spline smoothing
          do 2799 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'iequal  ') iequal = rpar(j) + 0.1
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'auto    ') iauto  = 1
            if(vname(i).eq.'noauto  ') iauto  = 0
            if(vname(i).eq.'smpar   ') smpar  = rpar(j)
 2799     continue
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia     = isels(1)
        nn     = nwert(ia)
        numspl = numor(ia)
        nwspl  = nwert(ia)
        if(nneu.le.0.or.nneu.gt.mwert) nneu = nn
        if(iequal.le.0) iequal = 1
        write(6,*)'nneu  = ',nneu
        if(nbuf.lt.size(nwert))then
           nbuf=nbuf+1
           ib=nbuf
        endif
!
        if(iauto.eq.1) then
           write(6,*)'iequal= ',iequal
           call csscv(nn,xwerte(1,ia),ywerte(1,ia),iequal,break,cscoef)
        else
           write(6,*)'smpar = ',smpar
           do 3777 i=1,nn
             weight(i) = sqrt(abs(ywerte(i,ia)))
 3777      continue
           call cssmh(nn,xwerte(1,ia),ywerte(1,ia),weight,smpar,break,  &
     &                cscoef)
        endif
         if(iout().gt.1) then
            write(6,*)'break         cscoeff'
            do 2776 i=1,nn
              write(6,2777)i,break(i),(cscoef(j,i),j=1,4)
 2777         format(1x,i3,':',5e14.6)
 2776       continue
         endif
         x01 = break(1)
         x02 = break(nn)
         dx  = (x02-x01)/(nneu-1)
         do 2702 i=1,nneu
            xn = x01 + (i-1)*dx
            xwerte(i,ib) = xn
            ywerte(i,ib) = fsplin(xn)
 2702    continue
         call txfpar(ia,ib)
         numor(ib)= - numor(ia)
         nwert(ib)= nneu
         isels(1) = ia
         ifits(1) = ib
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'fftmx   '.or.comand.eq.'fft-ms  ') then
!                    ---> multiple scattering by fft
           t1 = 0
           if(numspl.eq.0.or.nwspl.eq.0) then
             write(6,*)'spline must be made prior to fft'
             goto 2000
           endif
! -------- identify the original to get parameters ---
           do 3705 ia=1,nbuf
           if(numspl.eq.numor(ia)) goto 3707
 3705     continue
           write(6,*)'original data not found. unable to set parameters'
           goto 2000
 3707     continue
          if(nbuf.lt.size(nwert))then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          if(nx.gt.mwert) nx = mwert
          call txfpar(ia,ib)
          isels(1) = ia
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do 3799 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'mx      ') imx    = 0
            if(vname(i).eq.'rmx     ') imx    = 1
            if(vname(i).eq.'rmx1    ') imx    = 2
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') then
               xmax   = rpar(j)
               ymax   = 10.0/xmax
            endif
            if(vname(i).eq.'ymax    ') ymax   = rpar(j)
            if(vname(i).eq.'nx      ') nx     = rpar(j) +.01
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
 3799     continue
!
         if(imx.eq.0) then
           write(6,*)'evaluate multiple scattering ...'
           call fftmx(xwerte(1,ib),ywerte(1,ib),trans,xmax,dx,nx,nfft)
           numor(ib)= mod(numor(ia),numpls)+10*numpls
           yname(ib)= 'mx('//yname(ia)(1:5)
         else
           write(6,*)'evaluate inverse multiple scattering ...'
           if(imx.eq.1) then
             call fftrmx(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,     &
     &                                                  xmax,dx,nx,nfft)
           else
             call rmx1d(xwerte(1,ib),ywerte(1,ib),ai0,trans,thick,      &
     &                                                xmax,ymax,nx,nfft)
           endif
           numor(ib)= mod(numor(ia),numpls)+11*numpls
           yname(ib)= 'imx('//yname(ia)(1:4)
           call parset('i0      ',ai0,ib)
         endif
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nx
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'dmx     ') then
!                    ---> multiple scattering by fft
           t1 = 0
          ia = isels(1)
          if(nbuf.lt.size(nwert))then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          call txfpar(ia,ib)
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do i=1,inames
            j = inapa(i)
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') xmax   = rpar(j)
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
          enddo
!
         call demux(xwerte(1,ia),ywerte(1,ia),nwert(ia),                &
     &              trans,xmax,nfft,                                    &
     &              iout(),xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'dmx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'mux     ') then
!                    ---> multiple scattering by fft
           t1 = 0
          ia = isels(1)
          if(nbuf.lt.size(nwert))then
            nbuf=nbuf+1
            ib=nbuf
          endif
! --- transfer all data to be kept to the new buffer item ---
          call txfpar(ia,ib)
          nsel     = 1
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)

          do i=1,inames
            j = inapa(i)
            if(vname(i).eq.'trans   ') trans  = rpar(j)
            if(vname(i).eq.'thick   ') thick  = rpar(j)
            if(vname(i).eq.'xmax    ') xmax   = rpar(j)
            if(vname(i).eq.'nfft    ') nfft   = rpar(j) +.01
          enddo
!
         call   mux(xwerte(1,ia),ywerte(1,ia),nwert(ia),                &
     &              trans,xmax,nfft,                                    &
     &              iout(),xwerte(1,ib),ywerte(1,ib),ier)
         numor(ib)= mod(numor(ia),numpls)+11*numpls
         yname(ib)= 'mx('//yname(ia)(1:4)
         call parset('i0      ',ai0,ib)
         call parset('thick   ',thick,ib)
         call parset('trans   ',trans,ib)
         nwert(ib) = nfft
         isels(1) = ib
         t2 = 0
         write(6,*)'done '
         goto 2000
        endif
!
       if(comand.eq.'tracorr ') then
          ia = isels(1)
          call parget('thick   ',thick,ia,ier)
          call parget('trans   ',trans,ia,ier)
          do i=1,nwert(ia)
            ywerte(i,ia) = ywerte(i,ia)*thick*trans
            yerror(i,ia) = yerror(i,ia)*thick*trans
          enddo
          write(6,*)'selected item:',ia,' multiplied by:',              &
     &              thick,' * ',trans,' !'
         goto 2000
        endif
!
       if(comand.eq.'des     '.or.comand.eq.'desmear ') then
!                    ---> infinite slit desmearing
          do 2899 i=1,inames
            j = inapa(i)
            if(vname(i).eq.'nneu    ') nneu   = rpar(j) + 0.1
            if(vname(i).eq.'qmax    ') qmax   = rpar(j)
            if(vname(i).eq.'errabs  ') errabs = rpar(j)
            if(vname(i).eq.'errrel  ') errel  = rpar(j)
 2899     continue
         if(numspl.eq.0.or.nwspl.eq.0) then
           write(6,*)'spline must be made prior to desmearing'
           goto 2000
         endif
! --- identify the original to get parameters ---
         do 3005 ia=1,nbuf
          if(numspl.eq.numor(ia)) goto 3007
 3005    continue
          write(6,*)'original data not found. unable to set parameters'
          goto 2000
 3007    continue
         call parget('delqv   ',delqv,ia,ier)
         if(nbuf.lt.size(nwert))then
            nbuf=nbuf+1
         endif
         ic=nbuf
!
! ---- determine limits and steps of the new-data ----
         qstart = break(1)
         qend   = qmax
         if(qend.ge.break(nwspl)) qend = break(nwspl-1)
         if(nneu.eq.0)     nneu = 100
         if(nneu.gt.mwert) nneu = mwert
         dq  = (qend - qstart) / ( nneu-1 )
         write(6,*)'qstart = ',qstart
         write(6,*)'nwspl=',nwspl,'  break(nwspl) =',break(nwspl)
         write(6,*)'qend   = ',qend,'  dq=',dq
! ---- desmearing of all points ----
         do 3003 i=1,nneu
           qziel  = qstart + (i-1)*dq
           gunten = 0
           goben  = sqrt( break(nwspl-1)**2 - qziel**2 )
!       ---> integrate now:
           call erset (0,0,0)
           call qdng(fdes,gunten,goben,errabs,errel,result,errret)
!ccccc     call qdags(fdes,gunten,goben,errabs,errel,result,errret)
           call erset (0,1,0)
           if(iout().gt.0)write(6,*)i,' errret=',errret
           xwerte(i,ic) = qziel
           ywerte(i,ic) = -delqv/pi * result
 3003    continue
         call txfpar(ia,ic)
         nwert(ic)= nneu
         numor(ic)= mod(numspl,numpls) + 8*numpls
         isels(1) = ic
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'ia      '.or.comand.eq.'i-abso  ') then
!                    ---> correction absolut i <-------
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
        ia = isels(1)
           if(nbuf.lt.size(nwert))then
           nbuf=nbuf+1
           ib=nbuf
        endif
         if(yname(ia).ne.'i      ') then
           write(6,*)'date are no intensity',yname(ia)
           goto 2000
         endif
         call parget('i0      ',ai0  ,ia,ier)
         call parget('d       ',d    ,ia,ier)
         call parget('delqv   ',delqv,ia,ier)
         call parget('delqh   ',delqh ,ia,ier)
         call parget('lambda  ',alam0 ,ia,ier)
         write(6,*)'i-absolut  i0 =',ai0
         write(6,*)'           d  =',d, ' angstroem'
         write(6,*)'       lambda =',alam0,'  "  '
         write(6,*)'       delqv  =',delqv,' angstroem**-1'
         write(6,*)'       delqh  =',delqh,'    "         '
         fia = (2*pi/alam0)**2 / (ai0*d*delqv*delqh)
         write(6,*)'       normalizing factor =',fai
         do 2021 i=1,nwert(ia)
             xwerte(i,ib)=xwerte(i,ia)
             ywerte(i,ib)= fia*ywerte(i,ia)
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+5*numpls
 2021    continue
         yname(ib)='i-abso  '
         write(6,*)'q1=',ywerte(1,ib)
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
       if(comand.eq.'clip      ') then
!                    ---> remove points
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif

         do i=1,nsel
          ia = isels(i)
          ifrom = intval('from    ',1,inew)
          ito   = intval('to      ',0,inew)
          if(ifound('last    ').ne.0)then
           ilas = intval('last    ',1,inew)
           ito  = nwert(ia)
           ifrom= ito-ilas+1
          endif
          if(ito   .lt.1) ito   = 1
          if(ifrom .lt.1) ifrom = 1
          if(ito   .gt.nwert(ia) ) ito  = nwert(ia)
          if(ifrom .gt.nwert(ia) ) ifrom= nwert(ia)
          write(6,*)'removing points ',ifrom,ito,' from ',i
          errmax = getval('errmax  ',1d33,inew)
! flag the points to be removed !
          if(found('from    ') .or.                                     &
     &       found('to      ') .or.                                     &
     &       found('last    ')) then
           do j=ifrom,ito
             yerror(j,ia) = 3e33
           enddo
          endif
! and remove the points
          nnw = 0
          do j=1,nwert(ia)
           if(found('rel     ')) then
             errmx = errmax*abs(ywerte(j,ia))
           else
             errmx = errmax
           endif
           if(yerror(j,ia).lt.errmx) then
            nnw = nnw+1
            xwerte(nnw,ia) = xwerte(j,ia)
            ywerte(nnw,ia) = ywerte(j,ia)
            yerror(nnw,ia) = yerror(j,ia)
           endif
          enddo
          nwert(ia) = max(nnw,1)
          write(6,*)'number of points left (min=1) : ',nwert(ia)
         enddo

         goto 2000
        endif
!
!
       if(comand.eq.'fun     '.or.comand.eq.'function') then
!                    ---> y --> fun(y)
         ifunx = 0
         ifuny = 0
         do 2701 i=1,inames
           j = inapa(i)
           if(vname(i).eq.'op       ') fonpla = .true.
           if(vname(i).eq.'np       ') fonpla = .false.
           if(vname(i).eq.'x        ') ifunx=0
           if(vname(i).eq.'re(x)    ') ifunx=1
           if(vname(i).eq.'log(x)   ') ifunx=2
           if(vname(i).eq.'exp(x)   ') ifunx=3
           if(vname(i).eq.'x**2     ') ifunx=4
           if(vname(i).eq.'sqrt(x)  ') ifunx=5
           if(vname(i).eq.'x*       ') then
                                       ifunx=6
                    if(inpar(i).ne.0)  fx   = rpar(j)
           endif
           if(vname(i).eq.'x+       ') then
                                       ifunx=7
                    if(inpar(i).ne.0)  sx   = rpar(j)
           endif
           if(vname(i).eq.'rouse    ') then 
                                       ifunx=8
                                       ifuny=-1
           endif
           if(vname(i).eq.'zimm     ') ifunx=9

           if(vname(i).eq.'y        ') ifuny=0
           if(vname(i).eq.'re(y)    ') ifuny=1
           if(vname(i).eq.'log(y)   ') ifuny=2
           if(vname(i).eq.'exp(y)   ') ifuny=3
           if(vname(i).eq.'y**2     ') ifuny=4
           if(vname(i).eq.'sqrt(y)  ') ifuny=5
           if(vname(i).eq.'y*       ') then
                                       ifuny=6
                    if(inpar(i).ne.0)  fy   = rpar(j)
           endif
           if(vname(i).eq.'y+       ') then
                                       ifuny=7
                    if(inpar(i).ne.0)  sy   = rpar(j)
           endif
           if(vname(i).eq.'deff     ') ifuny=8
 2701    continue
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
! --- treat all selected items ---
         do 25021 ii=1,nsel
           ia = isels(ii)
           if(fonpla) then
!                     ----> treat data on place
             ib = ia
           else
!          ----> look for new place for the treated data ---
             if(nbuf.lt.size(nwert))then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
!          -----(of if(fonpla...

           call txfpar(ia,ib)

           do 2502 i=1,nwert(ia)
             if(ifunx.eq.0)          xwerte(i,ib)=xwerte(i,ia)
             if(ifunx.eq.1) then
               if(xwerte(i,ia).ne.0) xwerte(i,ib)=1 / xwerte(i,ia)
             endif
             if(ifunx.eq.2) then
               if(xwerte(i,ia).gt.0) xwerte(i,ib)= alog(xwerte(i,ia))
             endif
             if(ifunx.eq.3) then
                                     xwerte(i,ib)= exp(xwerte(i,ia))
             endif
             if(ifunx.eq.4) then
                                     xwerte(i,ib)=  xwerte(i,ia)**2
             endif
             if(ifunx.eq.5) then
               if(xwerte(i,ib).ge.0.) then
                                     xwerte(i,ib)=  sqrt(xwerte(i,ia))
               else
                                   xwerte(i,ib)=  0
               endif
             endif
             if(ifunx.eq.6) then
                                     xwerte(i,ib)=  fx * xwerte(i,ia)
             endif
             if(ifunx.eq.7) then
                                     xwerte(i,ib)=  sx + xwerte(i,ia)
             endif

                                 ! Rouse Scaling !
             if(ifunx.eq.8) then
              if(i.eq.1) then
               diffrouse = 0
               wl4 = getval('wl4     ',dble(wl4),inew)
               call  parget('wl4     ',wl4,ia,ier)
               call  parget('q       ',qrs,ia,ier)
               if(ier.ne.0) then
                 call errsig(1111,'q-value not found ...$')
               endif
               call parset('wl4     ',wl4,ib)
               cbuffer = chrval('usediff  ','none    ',inew)
               if(cbuffer(1:4) .ne. 'none' ) then
                  diffscal  = getval('unit    ', 1d7, ier)   ! default diffusion in cm**2/s ==> scale to A**2/ns
                 call parget(cbuffer(1:8),diffrouse,ia,ier)
                 diffrouse = diffrouse * diffscal
                 write(6,'("scaling:",i6,"  at q= ",f12.6," with diff ",f12.6," A**2/ns  from:",a)')ia,qrs,diffrouse,cbuffer(1:8)
               endif
              endif
!!              write(6,*)i,ia,diffrouse,qrs,xwerte(i,ia), exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              ywerte(i,ib) = ywerte(i,ia) * exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              yerror(i,ib) = yerror(i,ia) * exp( +diffrouse *qrs**2 * xwerte(i,ia) )
              xwerte(i,ib) = (qrs**2)*sqrt(wl4*xwerte(i,ia))
             endif

                                 ! Zimm Scaling !
             if(ifunx.eq.9) then
              if(i.eq.1) then
               eta_s = getval('eta     ',dble(eta_s),inew)
               temp_z= getval('temp    ',dble(temp_z),inew)
               call  parget('eta_solv  ',eta_s,ia,ier)
               call  parget('temp      ',temp_z,ia,ier)
               call  parget('q       ',qrs,ia,ier)
               if(ier.ne.0) then
                 call errsig(1111,'q-value not found ...$')
               endif
               call parset('eta_solv  ',eta_s,ib)
               call parset('temp      ',temp_z,ib)
              endif
              ywerte(i,ib) = ywerte(i,ia)
              yerror(i,ib) = yerror(i,ia)
              xwerte(i,ib) =                                            &
     &                ((qrs**3)*akb_z*temp_z/(6*pi*eta_s)*              &
     &                xwerte(i,ia))**(2.0/3.0)
             endif

             if(ifuny.eq.0) then
                                     ywerte(i,ib)=ywerte(i,ia)
                                     yerror(i,ib)=yerror(i,ia)
             endif
             if(ifuny.eq.1) then
               if(ywerte(i,ia).ne.0) then
                  ywerte(i,ib)=    1         / ywerte(i,ia)
                  yerror(i,ib)= yerror(i,ia) / ywerte(i,ia)**2
               endif
             endif
             if(ifuny.eq.2) then
               if(ywerte(i,ia).gt.0) then
                  ywerte(i,ib)= alog(ywerte(i,ia))
                  yerror(i,ib)= yerror(i,ia)/abs(ywerte(i,ia))
               endif
             endif
             if(ifuny.eq.3) then
                  ywerte(i,ib)= exp(ywerte(i,ia))
                  yerror(i,ib)= exp(ywerte(i,ia))*yerror(i,ia)
             endif
             if(ifuny.eq.4) then
                  ywerte(i,ib)= ywerte(i,ia) **2
                  yerror(i,ib)= 2 * yerror(i,ia) * abs(ywerte(i,ia))
             endif
             if(ifuny.eq.5) then
               if(ywerte(i,ib).ge.0.) then
                  ywerte(i,ib)=  sqrt(ywerte(i,ia))
                  if(ywerte(i,ia).gt.0) then
                   yerror(i,ib) = yerror(i,ia) / (2*ywerte(i,ib))
                  else
                   yerror(i,ib) = yerror(i,ia)
                  endif
               else
                  ywerte(i,ib)=  0
                  yerror(i,ib)=  0
               endif
             endif
             if(ifuny.eq.6) then
                ywerte(i,ib)=  fy * ywerte(i,ia)
                yerror(i,ib)=  fy * yerror(i,ia)
             endif
             if(ifuny.eq.7) then
                ywerte(i,ib)=  sy + ywerte(i,ia)
                yerror(i,ib)=       yerror(i,ia)
             endif
             if(ifuny.eq.8) then
               ywerte(i,ib) = 0
               if(ywerte(i,ia).gt.0) then
                if(xwerte(i,ia).ne.0) then
                  ywerte(i,ib) = 1.0/(ywerte(i,ia)*xwerte(i,ia)**2)
                  yerror(i,ib) = yerror(i,ia)                           &
     &                          /(ywerte(i,ia)*xwerte(i,ia))**2
                endif
               endif
             endif
 2502      continue
!          call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+(3+ifuny+10*ifunx)*numpls
           lxx = laenge(xname(ia),80,' ')
           lyy = laenge(yname(ia),80,' ')
           if(ifunx.eq.0) xname(ib) = xname(ia)
           if(ifunx.eq.1) xname(ib) = 're('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.2) xname(ib) = 'ln('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.3) xname(ib) = 'exp'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.4) xname(ib) = 'sqr('//xname(ia)(1:lxx)//')'
           if(ifunx.eq.5) xname(ib) = 'sqrt'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.6) xname(ib) = 'fx*'//xname(ia)(1:lxx)
           if(ifunx.eq.7) xname(ib) = 'sx+'//xname(ia)(1:lxx)
           if(ifunx.eq.8) xname(ib) =                                   &
     &                      'q**2*sqrt(Wl4*'//xname(ia)(1:lxx)//')'
           if(ifunx.eq.9) xname(ib) =                                   &
     &               '(q**3*kT/(6*eta*pi)*'//xname(ia)(1:lxx)//')**2/3'
           if(ifuny.eq.0) yname(ib) = yname(ia)
           if(ifuny.eq.1) yname(ib) = 're('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.2) yname(ib) = 'ln('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.3) yname(ib) = 'exp'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.4) yname(ib) = 'sqr('//yname(ia)(1:lyy)//')'
           if(ifuny.eq.5) yname(ib) = 'sqrt'//yname(ia)(1:lyy)//')'
           if(ifuny.eq.6) yname(ib) = 'fy*'//yname(ia)(1:lyy)
           if(ifuny.eq.7) yname(ib) = 'sy+'//yname(ia)(1:lyy)
           if(ifuny.eq.8) yname(ib) =                                   &
     &          '1/('//xname(ia)(1:lxx)//'**2 *'//yname(ia)(1:lyy)//')'
           call parget('ifunx   ',val,ib,ier)
           ival = val+0.1
           val  = ival*10 + ifunx
           call parset('ifunx   ',val,ib)
           call parget('ifuny   ',val,ib,ier)
           ival = val+0.1
           val  = ival*10 + ifuny
           call parset('ifuny   ',val,ib)
           if(ifunx.eq.6) then
             call parget('fx      ',val,ib,ier)
             if(ier.eq.1) val = 1
             val = val * fx
             call parset('fx      ',val,ib)
           endif
           if(ifuny.eq.6) then
             call parget('fy      ',val,ib,ier)
             if(ier.eq.1) val = 1
             val = val * fy
             call parset('fy      ',val,ib)
           endif
           if(ifunx.eq.7) then
             call parget('sx      ',val,ib,ier)
             if(ier.eq.1) val = 0
             val = val + sx
             call parset('sx      ',val,ib)
           endif
           if(ifuny.eq.7) then
             call parget('sy      ',val,ib,ier)
             if(ier.eq.1) val = 0
             val = val + sy
             call parset('sy      ',val,ib)
           endif
           isels(ii) = ib
           ifits(ii) = 0
           write(6,*)'item:',ia,' treated by fun: ',yname(ib),' --> ',ib
25021    continue
         goto 2000
        endif
!
!
       if(comand.eq.'fxy     ') then
        if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= fxy x=(expression) y=(expression)                                           '
           write(6,*)'=    expression=formula containing no blanks!                                 '
           write(6,*)'=    in expression the values of the selected data are referred as:  X, Y, ERR'
           write(6,*)'=    NOTE: variable names X, Y, ERR must be in capital letters   !!!          '
           write(6,*)'=    the processed records are newly generated and are selected               '
           write(6,*)'=    except if option op  (on place) is given (in that case the selected      '
           write(6,*)'=    records are modified                                                     '
           write(6,*)'=============================================================================='
           goto 2000
        endif
  
         if(ioldc==0) then
           reslin = inline(4:)
           ioldc  = 1
         endif

         if(ioldc.ne.0) then
           ixf=index(reslin,"x=")
           iyf=index(reslin,"y=")
           if(iyf > 0 .and. iyf < ixf) then
             call errsig(2999,"needs x-formula prior to y-formula$")
             goto 2000
           endif
           if(ixf==0) then 
!?             xformel="(X)"
           else
             if(iyf > 0) then
                xformel=reslin(ixf+2:iyf-1)
             else
                xformel=reslin(ixf+2:)
             endif 
           endif   
           if(iyf==0) then
!?             yformel="(Y)"
           else
             yformel=reslin(iyf+2:)
           endif               
 
           ioldc   = 0
         else
           write(*,*)"no contiunuation with formula found!"
           write(*,*)"use command csep to check and/or change spearator!"
           write(*,*)"default changed form ; to &"
           call errsig(999,"needs a formula after separator!$")
         endif
!
           write(6,*)'treat x by :',trim(xformel)
           write(6,*)'treat y by :',trim(yformel)

           fonpla = found('op      ')

         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
! --- treat all selected items ---
         do   ii=1,nsel
           ia = isels(ii)
           iibuf = ia
           if(fonpla) then
!                     ----> treat data on place
             ib = ia
           else
!          ----> look for new place for the treated data ---
             if(nbuf.lt.size(nwert))then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
           iibuf = ia
           iadda = ia
!          -----(of if(fonpla...
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               yyee = yerror(i,ia)
               call evaluate(xformel,val8x,ierx)
               call evaluate(yformel,val8y,iery)
               xwerte(i,ib) = val8x
               ywerte(i,ib) = val8y
! -- evaluate error --
               yyyy = ywerte(i,ia)+0.5*yerror(i,ia)
               call evaluate(yformel,val8y,iery)
               yyyy = ywerte(i,ia)-0.5*yerror(i,ia)
               call evaluate(yformel,val8x,iery)
               yerror(i,ib) = dabs( val8y - val8x)
           enddo
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+9*numpls

            if(xname(ib)(1:2) == "$$") then
              cbufferx     = xname(ib)(1:2) 
              nbufferx     = xname(ib)(3:)
           else
              nbufferx = xname(ib)
              cbufferx = " "
           endif
           if(yname(ib)(1:2) == "$$") then
              cbuffery     = yname(ib)(1:2) 
              nbuffery     = yname(ib)(3:)
           else
              nbuffery = yname(ib)
              cbuffery = " "
           endif
 
          charbuf = trim(xformel)//" "
! write(*,*)"TC1:", len(charbuf), trim(charbuf)," -- ",trim(xname(ib))
           if(cbufferx(1:2)=="$$")  call replace_str(charbuf,"/",":")
           call replace_str(charbuf,"X",trim(nbufferx))
           call replace_str(charbuf,"Y",trim(nbuffery))           

! write(*,*)"TC2:", len(charbuf), trim(charbuf)
           xname(ib)  = trim(cbufferx)//charbuf(1:80-len_trim(cbufferx))

           charbuf = trim(yformel)//" "
           if(cbuffery(1:2)=="$$")  call replace_str(charbuf,"/",":")
 !write(*,*)"TC1:", len(charbuf), trim(charbuf)," -- ",trim(yname(ib))
           call replace_str(charbuf,"X",trim(nbufferx))
           call replace_str(charbuf,"Y",trim(nbuffery))
 !write(*,*)"TC2:", len(charbuf), trim(charbuf)
           yname(ib)  = trim(cbuffery)//charbuf(1:80-len_trim(cbuffery))

           coment(ib) = 'x='//trim(xformel(1:36))//' y='//trim(yformel(1:36))
           isels(ii) = ib
           ifits(ii) = 0
          enddo
         goto 2000
        endif
!

       if(comand.eq.'funfun  ') then
           if(.not.found('immediate')) then
             open(15,file='formdat',status='UNKNOWN')
             read(15,'(a)',end=15999) xformel
             read(15,'(a)',end=15999) yformel
             close(15)
15999        continue
           endif
           write(6,*)'treat x by :',xformel
           write(6,*)'treat y by :',yformel

           fonpla = found('op      ')

         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
! --- treat all selected items ---
         do   ii=1,nsel
           ia = isels(ii)
           iibuf = ia
           if(fonpla) then
!                     ----> treat data on place
             ib = ia
           else
!          ----> look for new place for the treated data ---
             if(nbuf.lt.size(nwert))then
               nbuf=nbuf+1
               ib=nbuf
             endif
           endif
           iibuf = ia
           iadda = ia
!          -----(of if(fonpla...
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               yyee = yerror(i,ia)
               call evaluate(xformel,val8x,ierx)
               call evaluate(yformel,val8y,iery)
               xwerte(i,ib) = val8x
               ywerte(i,ib) = val8y
! -- evaluate error --
               yyyy = ywerte(i,ia)+0.5*yerror(i,ia)
               call evaluate(yformel,val8y,iery)
               yyyy = ywerte(i,ia)-0.5*yerror(i,ia)
               call evaluate(yformel,val8x,iery)
               yerror(i,ib) = dabs( val8y - val8x)
           enddo
           call txfpar(ia,ib)
           numor(ib)=mod(numor(ia),numpls)+9*numpls
           xname(ib)  = '!'//xformel(1:7)
           yname(ib)  = '!'//yformel(1:7)
           coment(ib) = 'x='//xformel(1:30)//'y='//yformel(1:30)
           isels(ii) = ib
           ifits(ii) = 0
          enddo
         goto 2000
        endif
!

!
       if(comand.eq.'seterr  ') then
!           write(6,*)'treat x by :',xformel
           write(6,*)'create y-errors by :',yformel

          if(nsel.eq.0) then
            write(6,*)'no curve selected'
            goto 2000
          endif
! --- treat all selected items ---
          do   ii=1,nsel
           ia = isels(ii)
           iadda = ia
           do i=1,nwert(ia)
               xxxx = xwerte(i,ia)
               yyyy = ywerte(i,ia)
               yyee = yerror(i,ia)
write(*,'(a,a,4f12.6)')"TEST: form1=",trim(yformel),xxxx,yyyy,yyee
               call evaluate(yformel,val8y,iery)
write(*,*)"evaluate: ",trim(yformel), val8y, iery
               yerror(i,ia) = val8y
write(*,'(a,a,4f12.6)')"TEST: form2=",trim(yformel),xxxx,yyyy,yyee,val8y
           enddo
          enddo
         goto 2000
        endif
!
       if(comand.eq.'invers  ') then
!                    ---> data ---> data**-1 vs q**2
         if(nsel.eq.0) then
           write(6,*)'no curve seleted'
           goto 2000
         endif
        ia = isels(1)
        if(xname(ia).ne.'q      '.or.numor(ia).lt.3*numpls) then
           write(6,*)'treat data with q !!'
           goto 2000
         endif
         do 2051 i=1,inames
           j = inapa(i)
           if(vname(i).eq.'bkgr     ') bkgr  = rpar(j)
 2051    continue
! ---- look for an additional data-slot ---
           if(nbuf.lt.size(nwert))then
           nbuf=nbuf+1
           ib=nbuf
        endif
!
         do 2052 i=1,nwert(ia)
           xwerte(i,ib) = xwerte(i,ia)**2
           if(ywerte(i,ia).eq.0.or.ywerte(i,ia)-bkgr.le.0) then
              yerror(i,ib) = 1e33
              ywerte(i,ib) = 0
            else
              ywerte(i,ib) = 1 / (ywerte(i,ia)-bkgr)
              yerror(i,ib) = yerror(i,ia)/(ywerte(i,ia)-bkgr)**2
            endif
 2052    continue
!
         call txfpar(ia,ib)
         numor(ib)=mod(numor(ia),numpls)+4*numpls
!
         xname(ib)='q**2   '
         yname(ib)='i**-1  '
         write(6,*)'ok (bkgr=',bkgr,')'
         isels(1) = ib
         ifits(1) = 0
         nsel     = 1
         goto 2000
        endif
!
!
       if(comand.eq.'dir     ') then
!                    ---

      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= dir                                                                        ='
       write(6,*)'= list loaded records in memory                                              ='
       write(6,*)'=    dir clength <n>                                                         ='
       write(6,*)'=            sets display length of comment                                  ='
       write(6,*)'=    dir  <parname1> <parname2> ....                                         ='
       write(6,*)'=            displays also values of specified parameters                    ='                  
       write(6,*)'=                                                                            ='
       write(6,*)'=    dir  unique#                                                            ='
       write(6,*)'=            creates unique numors                                           ='                  
       write(6,*)'=                                                                            ='
       write(6,*)'=    use dsl to check selections!                                            ='
       write(6,*)'=============================================================================='
       goto 2000
      endif


         if(found('unique# ') .or. found('u#      ')) call UniqueNumors


         len_comm = intval('clength ',len_comm,inew)
         if(len_comm.lt.1 ) len_comm=1
         if(len_comm.gt.80) len_comm=80
         if(inew.ne.0) then
           jl=2
         else
           jl=1
         endif

!
         do 170 i=1,nbuf
           csel = ' '
           do j=1,nsel
             if(i.eq.isels(j)) csel = '!'
             if(i.eq.ifits(j)) csel = '-'
           enddo

           ip = 0
           iy = 1
           pastring = ' '
           fostring = ' '
           if(inames.ge.jl) then
           do j=jl,min(inames,20)
             call parget(vname(j),pxxd,i,ier)
             if(ier.eq.0) then
               ip = ip+1
               dirpan(ip) = vname(j)
               dirpav(ip) = pxxd
             else
              iy=0
             endif
           enddo
            if(ip.gt.0) then
            write(fostring,'(a,i4,a)')'(',ip,'(1x,a,f12.6))'
            write(pastring,fostring)(trim(dirpan(j)),dirpav(j),j=1,ip)
            endif
           endif


!!           if(iy.ne.0) then
           write(6,171) csel,i,numor(i),name(i),yname(i),xname(i)       &
     &    ,'"'//coment(i)(1:len_comm)//'"',pastring(1:ip*21)
  171    format(1x,a1,i4,':#',i10,' : ',a8,' :',a8,' vs ',a8,'>',a,' ',a)
!!           endif
  170    continue

         if(numspl.ne.0)write(6,1711)numspl
 1711    format(/'     active spline generated from #',i14)
         goto 2000
       endif
!
!
       if(comand.eq.'purge   ') then
!                    -----
       if(found('help    ')) then 
        write(6,*)'=============================================================================='
        write(6,*)'= purge                                                                      ='
        write(6,*)'= remove items from the list of loaded data records                          ='
        write(6,*)'=    purge   rest                                                            ='
        write(6,*)'=            clears all records that are NOT SELECTED                        ='
        write(6,*)'=    purge sel                                                               ='
        write(6,*)'=            removes all SELECTED records                                    ='
        write(6,*)'=    purge all                                                               ='
        write(6,*)'=            clears record list completely                                   ='
        write(6,*)'=    Hint:                                                                   ='
        write(6,*)'=    use sel ....    to select records                                       ='
        write(6,*)'=    use dsl  or dir to check selections!                                    ='
        write(6,*)'=============================================================================='
        goto 2000
       endif

! 
          if(.not.(found('all     ').or.found('sel     ').or. found('rest    ' )) &
              .or.   (iparf() > 0 ) ) then
              call errsig(999,"NEED all, sel OR rest AS OPTION => nothing purged! $")
              goto 2000
          endif

          if(found('all     ')) then 
!             nwert = 0
             nsel  = 0
             isels = 0
             nfits = 0
             ifits = 0
!             xwerte = 0
!             ywerte = 0
!             yerror = 0
              call DataClean(1)        
              goto 2000    
          endif

          if(found('sel     ')) then
            if(nsel > 0) then
              nwert(isels(1:nsel)) = 0
            else
              call errsig(999,"NOTHING SELECTED => nothing purged! $")
            endif
          else
            do i=1,nbuf
             keep = .false.
             do j=1,nsel
               if(i.eq.isels(j)) then
                 keep = .true.                 
               endif
             enddo
             if(.not.keep) nwert(i) = 0   
           enddo
          endif

          m = 0
          do  j=1,nbuf
            if(nwert(j).ne.0) then
              m = m + 1
              if(j>m) call txfera(j,m)
            endif
          enddo
          nbuf = m
          nwert(nbuf+1:) = 0
18051     continue
          write(6,*)'nbuf is now ',nbuf
          ifits = 0
          if(found('all     ') .or. found('sel     ') ) then
            nsel  = 0
            isels = 0
            ifits = 0
            write(6,*)'selections are removed'
          else
            nsel  = m
            isels = [(i,i=1,m)]
          endif

          call DataClean(nbuf+1)   !! clean all residual data from freeded records 

         goto 2000
       endif
!
!
       if(comand.eq.'sel     ') then
!                    ---

      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= sel                                                                        ='
       write(6,*)'= select items from the list of loaded data records                          ='
       write(6,*)'=    sel  <i1> <i2> ... <in>                                                 ='
       write(6,*)'=            selects those with dir adresses i1, i2,.. in                    ='
       write(6,*)'=    sel  <i1>  -<i2>                                                        ='
       write(6,*)'=            selects all addresses between i1 and i2                         ='
       write(6,*)'=    sel  <i1>  -<i2> -<is>                                                  ='
       write(6,*)'=            selects all <is>-th addresses between i1 and i2                 ='
       write(6,*)'=    sel  all  <parname> <value>  band <value>  | <i1> ..                    ='
       write(6,*)'=            selects all records with parameter parnam value +- band         ='
       write(6,*)'=    sel  fit+                                                               ='
       write(6,*)'=            adds fit curves to the selection                                ='
       write(6,*)'=    sel next <parname> <value>  band <value>                                ='
       write(6,*)'=            selects the next record matching the parameter value            ='
       write(6,*)'=    sel add ......                                     (OR)                 ='
       write(6,*)'=         adds the selctions done via ... to the present list                ='
       write(6,*)'=    sel narrow <parname> <value>  band <value>                              ='
       write(6,*)'=         applies parameter selction to the prsent list (AND)                ='
       write(6,*)'=    sel exclude <parname> <value>  band <value>                             ='
       write(6,*)'=         removes matching records from the list                             ='
       write(6,*)'=    sel exclude numor mod <n>                                               ='
       write(6,*)'=         removes records with mod(#,n)=0 from the list                      ='
       write(6,*)'=    sel                                                                     ='
       write(6,*)'=         clear the list                                                     ='
       write(6,*)'=                                                                            ='
       write(6,*)'=    use dsl  or dir to check selections!                                    ='
       write(6,*)'=============================================================================='
       goto 2000
      endif


          
restore:   if(found('restore  ')) then
             open(20,file="lastselections",iostat=i)
               if(i .ne. 0) then
                 write(*,*)'WARNING: could not open lastselections'
               endif
               isels = 0
               ifits = 0
               read(20,*)nsel            
               if(nsel > 0) then
                  read(20,*)isels(1:nsel)    
               endif        
               read(20,*)nfits
               if(nfits > 0 ) read(20,*)ifits(1:nfits)                  
               read(20,*)nfsel
               if(nfsel > 0 ) read(20,*)isfits(1:nfsel)            
             close(20)
 
             goto 2000
          endif restore


narrow:   if(found('narrow  ')) then
             if(nsel <= 0) goto 2000
                selparna  = chrval('narrow  ',selparna  ,inew)
                selparval = getval( selparna ,dble(selparval) ,inew)
                selpartol = getval('band    ',dble(selpartol) ,inew)
                m = 0
                do i=1,nsel
                   call parget(selparna, parval_x , isels(i),ier)
                   if(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                   endif
                   nsel = m
                enddo   
             goto 2000
          endif narrow
          
exclude:   if(found('exclude  ')) then
                if(nsel <= 0) goto 2000

                if(found('numor   ')) then
                   l = intval('mod     ',2,inew)
                   if(l==0) l=100
                   m = 0
                   do i=1,nsel
                    if(.not.(mod(numor(isels(i)),l) == 0)) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                    endif
                    nsel = m
                   enddo   
                 goto 2000
                endif   
                
                selparna  = chrval('exclude  ',selparna  ,inew)
                selparval = getval( selparna ,dble(selparval) ,inew)
                selpartol = getval('band    ',dble(selpartol) ,inew)
                m = 0
                do i=1,nsel
                   call parget(selparna, parval_x , isels(i),ier)
                   if(.not.(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol)) then
                      m = m+1
                      isels(m) = isels(i)
                      ifits(m) = ifits(i)
                      write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                           isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                   endif
                   nsel = m
                enddo   
             goto 2000
          endif exclude
          
          
         if(inames.eq.0 .or. found('add     ')                          &
     &                  .or. found('fit+    ')                          &
     &                  .or. found('next    ')                          &
     &                  .or. found('all     ')                          &
     &                                       ) then
! keep current selection if add
           if(.not.found('add     ')) then
             m = 0
            else
             m = nsel
           endif

! go and look for data records with a given parameter
           if( found('next    ') .or. found('all     ') )then
             selparna  = chrval('next    ','        ',inew)
             selparna  = chrval('all     ',selparna  ,inew)
             selparval = getval( selparna ,dble(selparval) ,inew)
             selpartol = getval('band    ',dble(selpartol) ,inew)
                                   id = 1
             if(found('next    ')) id = isels(nsel)+1
             if(nsel.le.0)         id = 1
             do i=id,nbuf
                  iadda = i
                  call parget(selparna, parval_x , iadda,ier)
                  if(iout().gt.0) write(6,*) i , selparna, parval_x , ier
                  if(ier.eq.0 .and. abs(parval_x-selparval).lt.selpartol &
                              .and. numor(iadda)>=0                      ) then  ! this line experimental
                    if(found('next    ')) then
                      if(nsel.le.0) nsel = 1
                      if(found('add     ')) then
                        if(nsel.lt.size(nwert)) nsel = nsel+1
                      endif
                      isels(nsel) = i
                      if(iout().gt.0) write(6,*)'next:',nsel,isels(nsel)
                      goto 17801
                    endif
                    m = m+1
                    isels(m) = i
                    write(6,'(i4,"[",i4,"]:#",i9,3x,a,2x,e14.7)') &
                          isels(m),ifits(m),numor(isels(m)),trim(selparna),parval_x
                    if(iout().gt.0) write(6,*)'all ',m,i
                  endif
             enddo
             nsel = m
17801       continue

           else
            ifits  = 0
            isfits = 0
ipl:        do i=1,ipars
             iss      = Nint(rpar(i))
             if(abs(iss) <= nbuf) then
               write(6,*)'select adress   ',iss
             else
               write(6,*)"selected=",iss,"  nbuf=", nbuf
               call errsig(999,"ERROR: sel selcted address is out of range$")
               goto 2000
             endif
             m = m + 1
             isels(m) = iss
!             ifits(m) = 0
             if(iss.lt.0.and.i.gt.1) then
               m   = m - 1
               ias = isels(m)
               ies = -iss
               iselstep = 1
               if(ipars == 3 .and. Nint(rpar(ipars)) < 0) then
                 iselstep = abs(nint(rpar(ipars)))
               endif
               write(6,*)'select adresses from ',ias,' to ',ies,' step ',iselstep
               do l = ias,ies,iselstep
                isels(m) = l
!                ifits(m) = 0
                m = m + 1
               enddo
               m = m - 1
               if(iselstep .ne. 1) exit ipl
             endif
            enddo ipl
            nsel = m
         endif

          if(found('fit+    ')) then
             nfsel = 0
             do i=1,nsel
                numr = -numor(isels(i))
                do j=1,nbuf
                 if(numr.eq.numor(j)) then
                   write(6,*)i,' matches with: ',j
                   ifits(i) = j
                   isfits(i) = j
                   nfsel = nfsel+1
                 endif
                enddo
             enddo
           endif
           write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
           goto 2000
         endif !   ende if(inames.eq.0 ...........&


         do 7077 i=1,inames
           j = inapa(i)
           k = inpar(i)
           if(vname(i).eq.'fits    ') then
             do 70771 l=1,nsel
               numr = numor(isels(l))
               do 70772 ll=1,nbuf
                 if(-numr.eq.numor(ll)) then
                    ifits(l) = ll
                    write(6,*)'select ',ll,' as fit of scan: ',numr
                    goto 70773
                 endif
70772          continue
70773          continue
70771        continue
           endif
           if(vname(i).eq.'last    ') then
             l = 1
             nnumi(1) = numor(nbuf)
             call search(nnumi,l)
           endif
           if(vname(i).eq.'sc      ') then
             do 181 l=1,k
               nnumi(l) = Nint(rpar(l-1+j))
  181        continue
             call search(nnumi,k)
           endif
           if(vname(i).eq.'sc+     ') then
              do 182 l=1,nsel
               nnumi(l) = numor(isels(l))
  182         continue
             do 183 l=1,k
               kk = l+nsel
               if(kk.gt.size(inpar)) goto 183
                 nnumi(kk) = Nint(rpar(l-1+j))
                 iprs      = kk
  183        continue
             call search(nnumi,iprs)
           endif
 7077    continue

          write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
         goto 2000
       endif
!
       if(comand.eq.'putpar  ') then
!                    ------
         do i=1,nsel
            iaddp = isels(i)
            call parset (vname(1),sngl(rpar(1)),iaddp)
         enddo
         goto 2000
       endif

       if(comand.eq.'parlev  ') then
!                    ------   set display level
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= parlev <parname> <level>                                                    '
           write(6,*)'=        change display level of parameter                                    '
           write(6,*)'=        only parammeters with level < 1 are displayed on plots with default  '
           write(6,*)'=        -> plot parlev <level>    may shift the plot display level           '
           write(6,*)'=============================================================================='
           goto 2000
        endif

           do i=1,nsel
              iaddp = isels(i)
              call parset_display (vname(1),Nint(rpar(1)),iaddp)
           enddo

         goto 2000
       endif
!
       if(comand.eq.'rename  ') then
!                    ------
          if(found('help    ')) then 
           write(6,*)'=============================================================================='
           write(6,*)'= rename  xaxis <xaxisname>                                                   '
           write(6,*)'= rename  yaxis <yaxisname>                                                   '
           write(6,*)'= rename  name  <shortname>                                                   '
           write(6,*)'=                                                                             '
           write(6,*)'= max. length of names: ', len(xname(1)), len(yname(1)), len(name(1))
           write(6,*)'= for compatibility reasons try to stay with length of 8                      '
           write(6,*)'=============================================================================='
           goto 2000
        endif

             do j=1,inames
               do i=1,nsel
                 iaddp = isels(i)
                 if(vname(j).eq.'xaxis   ') xname(iaddp) = argvals(j+1)(1:len(xname(1)))
                 if(vname(j).eq.'yaxis   ') yname(iaddp) = argvals(j+1)(1:len(yname(1)))
                 if(vname(j).eq.'name    ')  name(iaddp) = argvals(j+1)(1:len( name(1)))
               enddo
               do i=1,nfsel
                 iaddp = isfits(i)
                 if(iaddp <= 0) cycle
                 if(vname(j).eq.'xaxis   ') xname(iaddp) = argvals(j+1)(1:len(xname(1)))
                 if(vname(j).eq.'yaxis   ') yname(iaddp) = argvals(j+1)(1:len(yname(1)))
                 if(vname(j).eq.'name    ')  name(iaddp) = argvals(j+1)(1:len( name(1)))
               enddo
             enddo
             goto 2000
           endif
!
!
       if(comand.eq.'dispsel '.or.comand.eq.'dsl     ') then
!                    --------                ---
!          write(6,390)nsel,(isels(i),numor(isels(i)),ifits(i),i=1,nsel)
   390     format(' selected items: ',i5/                                &
      &           ' scan-address:  numor:      fit-address:'/            &
      &           (1x,i8,5x,i12,5x,i12))


!
!!============ this section is nearly identical to dir an should be cast into a subroutine

         len_comm = intval('clength ',len_comm,inew)
         if(len_comm.lt.1 ) len_comm=1
         if(len_comm.gt.80) len_comm=80
         if(inew.ne.0) then
           jl=2
         else
           jl=1
         endif

         do  l=1,nsel
             i = isels(l)
            csel = ' '
!           do j=1,nsel
!             if(i.eq.isels(j)) csel = '!'
!             if(i.eq.ifits(j)) csel = '-'
!           enddo

           ip = 0
           iy = 1
           pastring = ' '
           fostring = ' '
           if(inames.ge.jl) then
           do j=jl,min(inames,20)
             call parget(vname(j),pxxd,i,ier)
             if(ier.eq.0) then
               ip = ip+1
               dirpan(ip) = vname(j)
               dirpav(ip) = pxxd
             else
              iy=0
             endif
           enddo
            if(ip.gt.0) then
            write(fostring,'(a,i4,a)')'(',ip,'(1x,a,f12.6))'
            write(pastring,fostring)(trim(dirpan(j)),dirpav(j),j=1,ip)
            endif
           endif


!!           if(iy.ne.0) then
           write(6,'(1x,a1,i4,2H [,i4,1H],2H:#,i10,3H : ,a8,1H:,a8,4H vs ,a8,1H>,a,1x,a)') &
     &     csel,i,ifits(l),numor(i),name(i),xname(i),yname(i)       &
     &    ,'"'//coment(i)(1:len_comm)//'"',pastring(1:ip*21)
!!           endif
         enddo

         if(numspl.ne.0)write(6,1711)numspl
 

         goto 2000
       endif
!
!
!
       if(comand.eq.'edit    ') then
!                    ----

         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         if(vname(1).eq.'sc      '.or.vname(1).eq.'numor   ') then
          do 499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 498
           endif
  499     continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
  498     continue
         endif
!
!??obsolete??         if(vname(1).eq.'n       '.or.vname(1).eq.'name    ') then
!??obsolete??            do 599 i=1,nbuf
!??obsolete??            if(vname(2).eq.name(i)) then
!??obsolete??             ispc = i
!??obsolete??             goto 598
!??obsolete??           endif
!??obsolete??  599     continue
!??obsolete??          write(6,*)'data with name ',vname(2),' not found!'
!??obsolete??          ierrs = 201
!??obsolete??          goto 2000
!??obsolete??  598     continue
!??obsolete??         endif
!??obsolete??

         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
! ----------> write data onto buffer <---------------------------------
         call savdat('datbuf  ',ispc)
! ----------> enter the system editor with that file <-----------------
!!!         call system('emacs datbuf')
!!         call get_environment_variable("OSTYPE",cbuffer)
!!           if(index(cbuffer,"linux") > 0) then
!!           call execute_command_line("emacs datbuf")
!!         else
!!!           call execute_command_line('open -eW datbuf')
!!            call execute_command_line('vim datbuf')
!!!           call execute_command_line('open -W -n -a MacVim datbuf')
!!         endif
          call execute_command_line(trim(editor)//" datbuf")
! --- reread it to the same place ---
!         nbuff = nbuf
!         nbuf = ispc - 1
!         inames = 1
!         ipars = 0
!         iargvs = 1
!         argvals(1) = 'datbuf  '
!
!        call input
!         nbuf = nbuff
          ioldc  = 1
          reslin = 'in datbuf '
     
         goto 2000
       endif
!
!
!
       if(comand.eq.'save    ') then
!                    ----
         fsname = 'lastsave'
         ispc = rpar(1) + 0.001
         if(ispc.eq.0) ispc = isels(1)
         do i=1,inames
           if(trim(vname(i)).eq.'to'.and.i.lt.inames) then
              fsname = vname(i+1)
              exit
           endif
         enddo
         do j=1,inames
         if(vname(1).eq.'sc      '.or.vname(j).eq.'numor   ') then
          do 4499 i=1,nbuf
           if(numor(i).eq.ispc) then
             ispc = i
             goto 4498
           endif
 4499     continue
          write(6,*)'numor: ',ispc,' not found !'
          ierrs = 201
          goto 2000
 4498     continue
         endif
         enddo
!
         if(vname(1).eq.'n       '.or.vname(i).eq.'name    ') then
           do i=1,nbuf
            if(vname(2).eq.name(i)) then
             ispc = i
             goto 4598
            endif
           enddo
          write(6,*)'data with name ',vname(2),' not found!'
          ierrs = 201
          goto 2000
 4598     continue
         endif
         if(ispc.gt.nbuf.or.ispc.lt.1) then
           write(6,*)'scan-address ispc=',ispc,' is out of range!'
           ierrs = 202
           goto 2000
         endif
! ----------> write data onto buffer <---------------------------------
         call savdat(fsname,ispc)
         goto 2000
       endif
!
       if(comand.eq.'msave   ') then
!                    -----
         fsname = 'lastsave'
         if(inames.gt.0)   fsname = argvals(1)
 !        title = trim(fsname)  !! experimental... 
         call msavdat(fsname)
         if(found("pl      ")) then   ! also plotting by option "pl"
           tbuffer = title
           title = trim(title)//"("//trim(fsname)//")"
           call splot(.true.)
           write(charbuf,'(a)')'cp "'//trim(LAST_DTR_PLOT)//'" '//'plot_'//trim(fsname)//'.pdf' 
           write(*,*)"last plot: ",trim(LAST_DTR_PLOT)
           write(*,*)"executing: ",trim(charbuf)
           call execute_command_line(trim(charbuf)) 
           title = tbuffer
         endif
         goto 2000
       endif
!
       if(comand.eq.'thc     ') then
!                    ---
         npoint = rpar(1) + 0.001
         call couple(1)
         call thc(npoint)
         goto 2000
       endif
!
       if(comand.eq.'fit     ') then
         call fit
         goto 2000
       endif

       if(comand.eq.'serfit   ') then
         call serial_fit
         goto 2000
       endif


       if(comand.eq.'creatser   ') then
         call create_ser
         goto 2000
       endif

       if(comand.eq.'fit_plus ') then
         call fit_plus
         goto 2000
       endif

       if(comand.eq.'ga_fit  ') then
          call ga_fit
          goto 2000
       endif

       if(comand.eq.'paraout     ') then
         call outputparams
         goto 2000
       endif
!
        if(comand.eq.'th_init  ') then
                                write(*,*) 'not yet implemented, sorry will coming soon'
        !call init_theories(thenam,thparn,nthpar,thrapar,thramin,thramax,mth,mtpar)
         goto 2000
       endif
!
       if(comand.eq.'yfitform ') then
        if(found('help    ')) then 
         write(6,*)'=============================================================================='
         write(6,*)'= yfitform y= <myformula>                                                    ='
         write(6,*)'=     specify an expression = <myformula> that is evaluate with the special  ='
         write(6,*)'=     theory eval                                                            ='
         write(6,*)'=                                                                            ='
         write(6,*)'=     the expression may not contain blanks                                  ='
         write(6,*)'=     the fit parameters can be specified as p(1) ... p(10)  (max no. is 10) ='
         write(6,*)'=     the indepependent variable is specified as X                           ='
         write(6,*)'=     prameters assoiated to the records are referred by their name          ='
         write(6,*)'=                                                                            ='
         write(6,*)'= example a parabola: yfitform y=p(1)+p(2)*X+p(3)*X*X                        ='
         write(6,*)'=============================================================================='
         goto 2000
      endif

         if(ioldc.ne.0) then
           yfitform = trim(reslin)
           ioldc   = 0
         else
           i = index(inline,"y=")
           yfitform = trim(inline(i+2:))
           if(i>0) then
           else
             write(*,*)"no formula found!"
             write(*,*)"use yfitform y=<myformula> !"
             write(*,*)"eval parameter names are p(1) ... p(9)"
             write(*,*)"name for x-variable: X"
             call errsig(999,"needs a formula !$")
           endif
         endif
         write(6,*)'y-fit-formel=',trim(yfitform)
         goto 2000
       endif
!
       if(comand.eq.'yformel  ') then
!                    -------
!cc      yformel = title
         if(ioldc.ne.0) then
           yformel  = adjustl(reslin)
           ioldc   = 0
         else
           write(*,*)"no contiunuation with formula found!"
           write(*,*)"use command csep to check and/or change spearator!"
           write(*,*)"default changed form ; to &"
           call errsig(999,"needs a formula after separator!$")
         endif
         write(6,*)'y-formel=',trim(yformel )
         goto 2000
       endif
!
       if(comand.eq.'xformel  ') then
!                    -------
!cc      xformel = title
         if(ioldc.ne.0) then
           xformel  = adjustl(reslin)
           ioldc   = 0
         else
           write(*,*)"no contiunuation with formula found!"
           write(*,*)"use command csep to check and/or change spearator!"
           write(*,*)"default changed form ; to &"
           call errsig(999,"needs a formula after separator!$")
         endif
         write(6,*)'x-formel=',trim(xformel )
         goto 2000
       endif
!
!
       if(comand.eq.'fopen   ') then
!                    -----
         npax = 1
         npay = 1
         fpnam = vname(1)
         xpnam = vname(2)
         if(inpar(2).ne.0) npax  = rpar(inapa(2)) + 0.0001
         ypnam = vname(3)
         if(inpar(3).ne.0) npay  = rpar(inapa(3)) + 0.0001
         open(33,file=fpnam,status='UNKNOWN')
         write(6,*)'file ',fpnam,' a has been opened for write'
         write(6,*)'data stored on write: ',xpnam,'(',npax,')  vs ',    &
     &              ypnam,'(',npay,')'
         write(33,*)title
         write(33,*)fpnam,'    ',ypnam,'  vs   ',xpnam,'     111111'
         write(33,*)
         write(33,*)'values'
         goto 2000
       endif
!
!
       if(comand.eq.'write   ') then
!                    -----
!        --- take it from the first selected file ---
         if(nsel.eq.0) then
           write(6,*)'no files selected !'
           goto 2000
         endif
         iadd = isels(1)
         iaddf= ifits(1)
!          --- look first for parameters ---
            call parget(xpnam,xval,iadd,ier)
            if(ier.ne.0) call fpaget(xpnam,xval,npax,ier)
            if(ier.ne.0) then
              write(6,*)'parameter xpnam(',npax,'):',xpnam,' not found'
              goto 2000
            endif
            call parget(ypnam,yval,iadd,ier)
            if(ier.ne.0) call fpaget(ypnam,yval,npay,ier)
            if(ier.ne.0) then
              write(6,*)'parameter ypnam(',npay,'):',ypnam,' not found'
              goto 2000
            endif
            write(6,*)'xpnam=',xval,'   ypnam=',yval
            write(33,*)xval,yval
         goto 2000
       endif
!
       if(comand.eq.'close   ') then
!                    -----
         write(33,*)
         write(33,'(4h#eod)')
         close(33)
         write(6,*)'file ',vname(1),' a has been closed'
         goto 2000
       endif
!
!
       if(comand.eq.'title   '.or.comand.eq.'tit     ') then
!                    -----
         title = trim(inline(6:))
         goto 2000
       endif
!
       if(comand.eq.'theos   ' .or. comand.eq.'th      ') then
!                    -----> list available theories
        if(found('help    ')) then 
         write(6,*)'=============================================================================='
         write(6,*)'= theos [thnam]            (short th)                                        ='
         write(6,*)'=     optional parameter theory name (only ckeck for this theory and list it ='
         write(6,*)'=     otherwise give a list of all available (linked) theories               ='
         write(6,*)'=                                                                            ='
         write(6,*)'= HINT: the theory impementations are stored as *.f or *.f90 sourecs in      ='
         write(6,*)'=       subdirectory ~/datreat/src/theos                                     ='
         write(6,*)'=       other theroies may be found in ~/datreat/src/unused_theos            ='
         write(6,*)'=       moving of the sources form unused_theos to theos and subseqeunt      ='
         write(6,*)'=       make clean  and make (from subdir src) will install the new config.  ='
         write(6,*)'=       available theories in theos may serve as templates for own creations ='
         write(6,*)'=                                                                            ='
         write(6,*)'=  th[eos] <name> will check for a theory and set the uservar theoryok = 1   ='
         write(6,*)'=               if the theory is available                                   ='
         write(6,*)'=                                                                            ='
         write(6,*)'=============================================================================='
         goto 2000
      endif

         ileng=len_trim(vname(1))

         if (ileng.gt.0) then
            call setudf('theoryok ',0d0,ier)
            if(ier.ne.0) call errsig(999,"ERROR: cannot create/set uservar theoryok!$")

            
            if( output_th_explanation(vname(1)) ) then
               call setudf('theoryok ',1d0,ier)
              goto 2000
            else
            do i=1,mth
               if (thenam(i).eq.vname(1)) then
                       write(6,*)'______________________________________'
                       write(6,'(i3,": ",a8,i3)') i,thenam(i),nthpar(i)
                       write(6,'(i3,": ",a)')(j,trim(thparn(j,i))//' ',j=1,nthpar(i))
                       write(6,*)'______________________________________'
                       call setudf('theoryok ',1d0,ier)
               endif
            enddo            
           endif

         else
           iperm = [(i,i=1,mth)]
           call HPSORT (thenam, mth,1,8, IPERM, 1, HWORK, IER) 
           write(6,'(a)')" ***** available theories available *****"
           write(6,*)'---------------------------------------------------------------------------------'
 !          do i=1,mth
 !              if (thenam(i).ne.' ') then  !   kein name!!!
 !                      write(6,'(i3,": ",a8,i3," :")',advance='no') i,thenam(i),nthpar(i)
 !                      write(6,'(30a)')(trim(thparn(j,i))//' ',j=1,mtpar)
 !              endif
 !          enddo
!!             write(6,'(8(2x,a8))') thenam(1:mth)
!!          write(6,'(8(2x,a8))') thenam(iperm(1:ntheos))
          j = 0
          do i=1,mth
             if(thenam(iperm(i)) == "        ") cycle
             j = j+1
             if( mod(j,8) .ne. 0 .and. i < mth) then
               write(6,'(2x,a8)',advance='no') thenam(iperm(i))
             else
               write(6,'(2x,a8)',advance='yes')thenam(iperm(i))
             endif
          enddo
          write(6,*)'----------------------------------------------------------------------------------'
          write(6,*)
          write(6,*)' .... for more details on one of these: --> th <theoryname> ! '
          write(6,*)
         endif

   
         goto 2000
      endif
!
       if(comand.eq.'activate'.or.comand.eq.'ac      ') then
!                    -----> activate a theory
         call activa(0)
         call activa(2)
         goto 2000
       endif
!
      if(comand.eq.'label   ') then
!                   -----> assign a label to a parameter
         call lsearch(jpar,itcal,ierr)
         thpala(jpar,itcal) = vname(3)(1:4)
         goto 2000
      endif

      if(comand.eq.'extthpar') then
!                   --------> extract a thory parameter
         call lsearch(jpar,itcal,ierr)
         if(ierr .eq. 0) then
           ithc  = nthtab(itcal)
           call setudf('lastpar  ',dble( thparx(jpar,itcal) ), ierr)
           call setudf('lastperr ',dble( therro(jpar,itcal) ), ierr)
           write(*,'(a,i0,a,i0,a,e15.7,a)')"Parameter value(",jpar,",",itcal,") = ",thparx(jpar,itcal),"  put to variable lastpar"
           write(*,'(a,i0,a,i0,a,e15.7,a)')"Parameter error(",jpar,",",itcal,") = ",therro(jpar,itcal),"  put to variable lastperr"
         else
           call errsig(1999,"theory parameter not found! $") 
         endif
         goto 2000
      endif
!

      if(comand.eq.'chgthpar') then
!                   --------> change single theory parameters
         call lsearch(jpar,itcal,ierr)
         if(ierr .eq. 0) then
           ithc  = nthtab(itcal)
           thparx(jpar,itcal) =                                         &
     &        getval('par     ',dble(thparx(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)')                          &
     &                'change parameter ',thparn(jpar,ithc),            &
     &                ' of th(',itcal,'):',thenam(ithc),' to ',         &
     &                thparx(jpar,itcal)
           endif
           thpsca(jpar,itcal) =                                         &
     &        getval('scale   ',dble(thpsca(jpar,itcal)),inew)
           if(inew.ne.0) then
             write(6,'(a,a,a,i3,a,a,a,e13.6)')                          &
     &                'change scale of  ',thparn(jpar,ithc),            &
     &                ' of th(',itcal,'):',thenam(ithc),' to ',         &
     &                thpsca(jpar,itcal)
           endif
         endif
         goto 2000
      endif
!
      if(comand.eq.'couple  ') then
!                   ------> couple label to a parameter
         call lsearch(jpar,itcal,ierr)
         if(ierr.ne.0) goto 2000
         if(inames.lt.3) then
           write(6,*)'label must be specified'
           ierrs =804
           goto 2000
         endif
         ncpl = ncoup(jpar,itcal)
         do 8011 i=1,ncpl
            if(thpalc(i,jpar,itcal).eq.vname(3)(1:4)) then
               thpafc(i,jpar,itcal) = rpar(inapa(3))
               goto 2000
            endif
 8011     continue
          if(ncpl.eq.mcoup) then
            write(6,*)'too many coupled items'
            ierrs = 805
           else
            ncpl = ncpl+1
            ncoup(jpar,itcal) = ncpl
            thpalc(ncpl,jpar,itcal) = vname(3)(1:4)
            thpafc(ncpl,jpar,itcal) = rpar(inapa(3))
           endif
           call couple(0)
         goto 2000
      endif
!
!
!
       if(comand.eq.'acl     '.or.comand.eq.'aclast  ') then
!                    -----> reactivate theories as stored in lastth
         call activa(3)
         if(ierrr.eq.0. .and. iout().ge.0) call activa(2)
         goto 2000
       endif
!
       if(comand.eq.'desactiv'.or.comand.eq.'dac     ') then
!                    -----> desactivate all theories
         call activa(1)
         goto 2000
       endif
!
       if(comand.eq.'activlst'.or.comand.eq.'al      ') then
!                    -----> list activated theories
         call activa(2)
         goto 2000
       endif
!
       if(comand.eq.'get_th'.or.comand.eq.'gth      ') then
!                    -----> get theory appended to an file created by msave and copy it to lassth
!                           and load it
         if(inames == 1) then
!           i = LEN_TRIM(argvals(1))
!           write(6,*)"argvals: ",argvals(1)(1:i)
!           call extract_th(argvals(1)(1:i)//" ")
           call extract_th(argvals(1))
           if(ierrr==0) then 
             call activa(3)
             if(ierrr==0) call activa(2)
           endif
         else
           call errsig(999,"need exactly one filename as argument!$")
         endif
         goto 2000
       endif
!

       if(comand.eq.'gplot    ') then
!                    -----> plot selected curves
         call gplot()
         goto 2000
       endif
!
 
      if(comand.eq.'plot    '.or.comand.eq.'p       '.or.comand.eq.'gp      ') then
!                    -----> plot selected curves
         call splot(.true.)
!         ibild1 = ibild
         goto 2000
       endif
!
!
!
       if(comand.eq.'plot0   '.or.comand.eq.'p0      ') then
!                    -----> set parameters for plot
         call splot(.false.)
!         ibild1 = ibild
         goto 2000
       endif
!
!
!
       if(comand.eq.'zero   ') then
!                    ----      ------> logical zero the common storage
         nbuf= 0
         goto 2000
       endif
!
!
       if (comand.eq.'numorpls') then
!                    -----------> change offset between files
         numpls = Nint(rpar(1)) 
         goto 2000
       endif

       if (comand.eq.'numorchg') then
!                    -----------> change numor
         if(nsel > 0) then
           write(6,*)"chnage numor of record:",isels(1)," from ",numor(isels(1))," to ", Nint(rpar(1)) 
           numor(isels(1)) = Nint(rpar(1))
         endif
 
         goto 2000
       endif
!
!
       if (comand.eq.'cailtab ') then
!                    -----------> create a table for caille scattering
!                                 to enable complete resolution fitting

          if(nsel.le.0) then
              write(6,*)'need a dummy selection!'
              ierrs = 999
              goto 2000
           endif
          if(ntheos.ne.1) then
              write(6,*)'only one theory spec. is allowed'
              ierrs = 999
              goto 2000
           endif

           ith = nthtab(1)
           if(thenam(ith)(1:4).ne.'cail') then
              write(6,*)'theory(',ith,'): ',thenam(ith),                &
     &        '  does not belong to the cail-family'
              ierrs = 999
              goto 2000
           endif

          write(6,*)'Using theory(',ith,'): ',thenam(ith)

          qortho_0 = getval('qortho_0 ',qortho_0,inew)
          write(6,*)qortho_0, inew
          qortho_1 = getval('qortho_1 ',qortho_1,inew)
          write(6,*)qortho_1, inew
          northo   = intval('northo   ',northo,  inew)
          write(6,*)northo, inew
          qz_0     = getval('qz_0     ',qz_0,    inew)
          write(6,*)qz_0, inew
          qz_1     = getval('qz_1     ',qz_1,    inew)
          write(6,*)qz_1, inew
          nz       = intval('nz       ',nz,      inew)
          write(6,*)nz, inew
          filext   = chrval('to       ',filext,  inew)
          write(6,*)filext, inew

          filnam   = 'cailtab.'//filext

          write(6,*)'Writing Results to: ',filnam
          open(85,file=filnam,form='formatted',status='unknown')

          write(6,'(A,F18.7)')'qortho_0   =  ',qortho_0
          write(6,'(A,F18.7)')'qortho_1   =  ',qortho_1
          write(6,'(A,F18.7)')'qz_0       =  ',qz_0
          write(6,'(A,F18.7)')'qz_1       =  ',qz_1
          write(6,'(A,I11)')  'northo     =  ',northo
          write(6,'(A,I11)')  'nz         =  ',nz

          write(85,'(A)') filnam
          write(85,'(A,F18.7)')' ',qortho_0
          write(85,'(A,F18.7)')' ',qortho_1
          write(85,'(A,F18.7)')' ',qz_0
          write(85,'(A,F18.7)')' ',qz_1
          write(85,'(A,I11)')  ' ',northo
          write(85,'(A,I11)')  ' ',nz

          call couple(1)



         do iortho=1,northo

           qortho=qortho_0+(iortho-1)*((qortho_1-qortho_0)/(northo-1))

           do isel=1,nsel
             iadda = isels(isel)
             write(6 ,*)'qortho = ',qortho,'  PARSET to sel: ',iadda
             call parset('qortho  ' ,qortho,iadda)
                 !! isel
           enddo
           write(85,*)'qortho = ',qortho


           do iz=1,nz
             if(isigint().gt.0) then
               write(6,*)'Table computation stopped by Crtl-C'
               goto 20009
             endif
             qz = qz_0+(iz-1)*((qz_1-qz_0)/(nz-1))
             x = qz
             yctab = thval(x)
             write(85,*) yctab

                 !! iz
           enddo

               !! iortho
         enddo

20009    continue

         call theo_out(85)
         close(85)
         write(6,*)'Ready, Data written to: ',filnam

         goto 2000
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! multi exp fits !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(comand.eq.'mexp    ') then
!                    ---> estimate a multi exp representation 
         call dtr_mexp()

         goto 2000
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TOF-Stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(comand.eq.'tofdos  ') then
!                    ---> remove points
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif

         Ommax  = getval('omegamax',Ommax,inew)
         Nomega = intval('nomega  ',Nomega,inew)       

         write(6,*)'Preparing Pseudo Density of states on grid'
         write(6,*)'with ', Nomega,' points'
         write(6,*)'spanning to: ',Ommax,' GHz '

         do i=1,nsel
          ia = isels(i)
          if(nbuf.lt.size(nwert)) then 
            nbuf = nbuf+1
          else
            Write(6,*)'max. number of buffers reached .. '
            goto 2000   
          endif

          call  tof_dos(ia,nbuf,Ommax,Nomega,ierrr)
          if(ierrr.ne.0) goto 2000
          
          isels(i) = nbuf
           
         enddo

         goto 2000
        endif
!

       if(comand.eq.'detsens  ') then
!                    ---> sensitivity correction
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif


         He3_pressure_bar  = getval('p       ',dble(He3_pressure_bar) ,inew)
         Tube_diameter_m   = getval('dtube   ',dble(Tube_diameter_m),inew)       
         lenfrac           = getval('lenfrac ',dble(lenfrac),inew) 

         write(6,'(a)')'Performing wavelength dependent sensitivy correction on Tof data on lambda scale'
         write(6,*)' p       = He3-pressure/bar        = ',  He3_pressure_bar 
         write(6,*)' dtube   = tube diameter/m         = ',  Tube_diameter_m 
         write(6,*)' lenfrac = projected length factor = ',  lenfrac 

         do i=1,nsel
          ia      = isels(i)
          call parget('det_phe3',xx,ia,ier)
          if(ier .ne. 0 ) then
            do j=1,nwert(ia)
              lambdaA =  xwerte(j,ia)
              if(lambdaA > 0.01 ) then
                 detsens = detector_sensitivity(lambdaA,lenfrac,He3_pressure_bar,Tube_diameter_m)
                 ywerte(j,ia) = ywerte(j,ia)/detsens
                 yerror(j,ia) = yerror(j,ia)/detsens
 !                write(6,*)xwerte(j,ia), detsens
              endif
            enddo
            call parset ('det_phe3 ', He3_pressure_bar ,ia)
            call parset ('dtube    ', Tube_diameter_m  ,ia)
          else
            write(6,*)'Record: ',ia,' will not be treated a second time, skipped!'
          endif
         enddo

         goto 2000
        endif
!
!
       if(comand.eq.'gdos_q  ') then
!                    ------> see e.g. J.C. Smith et al, J.Chem. Phys. 85, 3636 (1986)
         if(nsel.eq.0) then
           write(6,*)'no curve selected'
           goto 2000
         endif
         ia      = isels(1)

         k1                = intval('k1      ', k1,      inew)
         k2                = intval('k2      ', k2,      inew)
         recstep           = intval('recstep ', recstep, inew)

         write(6,'(a)')'Creating q-exptrapolation curve for g(omega) deterimation '
         write(6,'(a)')'see e.g. J.C. Smith et al, J.Chem. Phys. 85, 3636 (1986)  '
         write(6,*)' from channel ',k1,' to ', k2,'   recstep ',recstep
         write(6,*)' lambda       ',xwerte(k1,ia),' to ',xwerte(k2,ia)

! use first record as template for energy scaling, parameter extraction etc.
         call parget('det_phe3',xx,ia,ier)
         if(ier .ne. 0 ) then
           write(6,*)'do detector efficiency correction first!'
           goto 2000
         endif
  
         call parget ('lambda  ',lambdaA ,ia, ier)
         if(ier.ne.0) then
            write(6,*)'ERROR: Parameter lambda is missing in parameterlist of record!'
            goto 2000
         endif

         call parget ('temp    ', temp ,ia, ier)
         if(ier.ne.0) then
            write(6,*)'ERROR: Parameter temp is missing in parameterlist of record!'
            goto 2000
         endif

         lambda0      = (xwerte(k1,ia)+xwerte(k2,ia))/2     * 1d-10                      ! das ist lambda_final
         lambda1      = lambdaA                             * 1d-10                      ! das ist lambda_in
         delta_lambda = (xwerte(k2+1,ia)-xwerte(k1,ia))     * 1d-10
         e0           = NeutronEnergy_fromLambda(dble(lambda0))
         e1           = NeutronEnergy_fromLambda(dble(lambda1))
         dE           = e1-e0                                    ! hquer*omega = E_in - E_final
         de0     = NeutronEnergy_fromLambda(dble(lambda1-delta_lambda/2))-NeutronEnergy_fromLambda(dble(lambda1+delta_lambda/2))
         dee     = NeutronEnergy_fromLambda(dble(lambda0-delta_lambda/2))-NeutronEnergy_fromLambda(dble(lambda0+delta_lambda/2))



         omega   = dE*2*Pi/Planckkonstante / 1d9
         om_cm   = dE/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1

         kinetic_factor = lambda0/lambda1      ! to be multiplied with the cross setion data to come towards s(q,om)
         channel_factor = dee/de0
         channel_width  = dee/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1
         channel_width0 = de0/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1

          
         if(nbuf.lt.size(nwert)) then 
            nbuf = nbuf+1
         else
            Write(6,*)'max. number of buffers reached .. '
            goto 2000   
         endif
         call txfera(ia,nbuf)
         numor(nbuf)    = numor(ia)+100000
         xname(nbuf)    = 'q      '
         yname(nbuf)    = 'gdos'
         coment(nbuf)   = coment(nbuf)(1:40)//'gdos extrapolation data '
         nwert(nbuf) = 0

         call parset ('omega    ', omega            ,nbuf)
         call parset ('om_cm    ', om_cm            ,nbuf)
         call parset ('rwidth   ', channel_factor   ,nbuf)
         call parset ('cwidth   ', channel_width    ,nbuf)
         call parset ('cwidth0  ', channel_width0   ,nbuf)
         call parset ('lambda1  ', xwerte(k1,ia)    ,nbuf)
         call parset ('lambda2  ', xwerte(k2,ia)    ,nbuf)
         
         
         do i=0,nsel-recstep,recstep
          sumq = 0
          sum  = 0
          sume = 0
 irs:     do irstp=1,recstep
           ia      = isels(i+irstp)
           call parget ('angle   ',angle_2tht, ia, ier)  
           if(ier.ne.0) then
             write(6,*)'ERROR: Parameter angle is missing in parameterlist of record!',isels(ia)
           endif

          do j=k1,k2
            lambda0      = xwerte(k1,ia)                       * 1d-10                      ! das ist lambda_final
            lambda1      = lambdaA                             * 1d-10                      ! das ist lambda_in
            delta_lambda = (xwerte(k2+1,ia)-xwerte(k1,ia))     * 1d-10
            e0           = NeutronEnergy_fromLambda(dble(lambda0))
            e1           = NeutronEnergy_fromLambda(dble(lambda1))
            dE           = e1-e0                                    ! hquer*omega = E_in - E_final
            q  = sqrt(((2*Pi/lambda0)**2+(2*Pi/lambda1)**2 - 2*COS(angle_2tht*Pi/180d0)*(2*Pi/lambda0)*(2*Pi/lambda1)))  
            qa = q * 1d-10
            sumq = sumq+qa

            omega   = dE*2*Pi/Planckkonstante / 1d9
            om_cm   = dE/Planckkonstante/Lichtgeschwindigkeit/100   ! werte in cm**-1
            kinetic_factor = lambda0/lambda1      ! to be multiplied with the cross setion data to come towards s(q,om)
!                                                   kinetic_factor=1 if already applied
            dbb = exp(dble(-dE/(temp*Boltzmannkonstante))) - 1.0
!                         ! definition of omega in paper
            sum  = sum  +  6*(-omega)/(qa**2)*dbb* kinetic_factor * ywerte(j,ia)
            sume = sume + (6*(-omega)/(qa**2)*dbb* kinetic_factor * yerror(j,ia))**2
!                             ! definition of omega in paper          
           enddo
          enddo irs

           nwert(nbuf) =  nwert(nbuf) + 1
           xwerte(nwert(nbuf), nbuf) = sumq       / (recstep*(k2-k1+1))
           ywerte(nwert(nbuf), nbuf) = sum        / (channel_factor*recstep)
           yerror(nwert(nbuf), nbuf) = sqrt(sume) / (channel_factor*recstep)
           write(6,'(2i4,3e14.7)')i, nwert(nbuf), xwerte(nwert(nbuf), nbuf), ywerte(nwert(nbuf), nbuf),  yerror(nwert(nbuf), nbuf)
         enddo

         nsel     = 1
         isels(1) = nbuf
         ifits    = 0
        
         goto 2000
        endif
!

!============================================== Checkings =================
       if(comand.eq.'argvals  ') then
         write(6,*)' inmames   = ', inames
         write(6,*)' argvals(1)= ', argvals(1)(1:60), ' ... '
         write(6,*)' argvals(2)= ', argvals(2)(1:60), ' ... '
         write(6,*)' argvals(3)= ', argvals(3)(1:60), ' ... '
         write(6,*)' argvals(4)= ', argvals(4)(1:60), ' ... '
         write(6,*)' argvals(5)= ', argvals(5)(1:60), ' ... '
         goto 2000
       endif








!
       call makro(comand)
!>> check       call unused( 1, 1, 1, ier)

!
       goto 2000
!      ---------> read now commandlines from makro file
!
!
!
      CONTAINS
      
       subroutine SaveStatus
            integer :: i, j, numx
 d1:        do i=1,nbuf-1
              numx = numor(i)
              do j=i+1,nbuf
               if(numx == numor(j)) then
                 Write(*,*)"Assigning UNIQUE NUMOR prior saving!"
                 call UniqueNumors
                 exit d1
               endif
              enddo
            enddo d1
  
            open(19,file="lastselections")
              write(19,'(i8)')nsel,(isels(i),i=1,nsel)            
              write(19,'(i8)')nsel,(ifits(i),i=1,nsel)            
              write(19,'(i8)')nsel,(isfits(i),i=1,nsel)            
            close(19)
            open(19,file="lastusv")
             write(19,'(a)')"makro"
             if(nousev > 0) write(19,'("set ",a,"  ",e16.8)')(usenam(i),useval(i),i=1,nousev)
            close(19)
     
            isels  = 0
            ifits  = 0
            isfits = 0 
            nsel = nbuf
            isels(1:nsel) = [(i,i=1,nsel)]
            fsname = "last_datreat_content"
            call msavdat(fsname)
            write(*,*) "actual content saved to last_datreat_content etc. ... "


       end subroutine SaveStatus

       subroutine UniqueNumors
           integer :: i, j, numx, numc, numbase
    
           numbase =  log10(dble(maxval(abs(numor(1:nbuf))) + 1)) + 1
           numbase = 10**numbase
                     
           numc    = 1
           call parset("prevnum ",real(numor(nbuf)),nbuf)
           do i = 1, nbuf-1
              if(numor(i) == 0) numor(i) = i
              call parset("prevnum ",real(numor(i)),i)
              numx = numor(i)
              do j = i+1,nbuf
                if(numx == numor(j)) then
                  call parset("prevnum ",real(numor(j)),j)
                  numor(j) = sign((abs(numor(j))+numc*numbase),numor(j))
                  numc = numc+1
                endif
              enddo                                  
           enddo
        end subroutine UniqueNumors



      END ! Ende der Hauptschleife







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_pair_of_points(xin1,yin1,yer1,n1, xin2,yin2,yer2,n2, xout, yout1, yerout1, yout2, yerout2, nout)
!      ===============================================================================================================
!
       use new_com
       implicit none

       real, intent(in)     :: xin1(n1), yin1(n1), yer1(n1)
       real, intent(in)     :: xin2(n2), yin2(n1), yer2(n2)
       integer, intent(in)  :: n1, n2

       real, intent(out)    :: xout(nout), yout1(nout), yerout1(nout), yout2(nout), yerout2(nout) 
       integer              :: nout
         

       integer              :: i, j, jc, ii
!?       real                 :: x, xx, y1, y2, e1, e2, dist,p,y,ye
       real                 :: x, xx, dist,p,y,ye
 
       nout = 0
! assuming nontrivial length check
       if(n1.lt.2 .or. n2.lt.2) then
         write(6,*)'at least one of the vectors has too few components:',n1,n2
         call errsig(990,"ERROR: vectors dont length-match$")
         return
       endif
! assuming equal length check
       if(n1.ne.n2) then
         write(6,*)'vectors have different number of components:',n1,n2
         call errsig(991,"ERROR: vectors dont length-match$")
         return
       endif
! assuming ordered vectors CHECK
       do i=1,n1-1
        if(xin1(i).ge.xin1(i+1)) then
          write(6,*)'vector 1 is not ordered ',i,xin1(i),xin1(i+1)
          call errsig(992,"ERROR: vectors must be x-sequentially ordered $")
          return
        endif
       enddo
       do i=1,n2-1
        if(xin2(i).ge.xin2(i+1)) then
          write(6,*)'vector 2 is not ordered ',i,xin2(i),xin2(i+1)
         call errsig(992,"ERROR: vectors must be x-sequentially ordered $")
          return
        endif
       enddo

! do the assignments !
       do i=1,n1
!      search for matching points       
!      1. closest mate
         dist = abs(xin1(i)-xin2(1))
         jc   = 1
         do j=2,n2
           xx = abs(xin1(i)-xin2(j))
           if(xx.lt.dist) then
              dist = xx
              jc   = j
           endif
         enddo


!        now our closest neighbour pair is i<-->jc
!        if they are equal we are done
         if(xin1(i).eq.xin2(jc)) then
           nout = nout + 1
           xout(nout)    = xin1(i)
           yout1(nout)   = yin1(i)
           yerout1(nout) = yer1(i)
           yout2(nout)   = yin2(jc) 
           yerout2(nout) = yer2(jc)
           cycle
        endif

        if(xin1(i).lt.xin2(jc)) then
          do ii=1,n1
            if(xin1(ii).gt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)  !! this would be good if full statistical independeence
                                                           !! is assumed, however interpolation creates correlation
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif

        if(xin1(i).gt.xin2(jc)) then
          do ii=n1,1,-1
            if(xin1(ii).lt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif
         
       enddo
       

       return
       end subroutine get_pair_of_points



 
   subroutine average_data()
!     -------------------------
      use dimensions
      use new_com
      ! use cincoc
      use xoutxx
      use cdata
      ! use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
!      use constants
      use PhysicalConstantsPlus
 
      implicit none

      double precision, save        :: xcatch = 0.05d0
      double precision              :: yem
      integer                       :: inew, ier
!?      integer                       :: i,j,k,n
      integer                       :: i,j,n
      integer                       :: number_of_data_points
      integer                       :: n_result_point
      character(len=8)              :: cbuf
      
      double precision              :: xaver, yaver, yerra
      double precision, allocatable :: xva(:), yva(:), yvaer(:)
      integer         , allocatable :: iperm(:)
      logical         , allocatable :: va_used(:)
      double precision, allocatable :: weights(:)
      real                          :: xwidth
      logical                       :: relative_catch 
      logical                       :: put_width   
      logical                       :: logy

      real                          :: xpa, xpaav
      integer                       :: ipa   


      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= average                                                                    ='
       write(6,*)'= combine error weighted close data points from records in the selection list='
       write(6,*)'= parameter is       xcatch  the RELATIVE or ABSOLUTE width of the x-window  ='
       write(6,*)'= current default is xcatch  ',xcatch,'  relative '
       write(6,*)'= option (volatile): absolute  catches within absolute distance               ='
       write(6,*)'= option (volatile): putwidth  set _xwidth parameter to min channel distance'
       write(6,*)'= option (volatile): log_y     averfage y on log scale                        '
       write(6,*)'=============================================================================='
       return
      endif





      xcatch = getval("xcatch  ",xcatch,inew)
      if(found("absolute ")) then
         relative_catch = .false.
      else
         relative_catch = .true.
      endif
      if(found("putwidth ")) then 
         put_width      = .true.
      else
         put_width      = .false.
      endif

      logy = found('log_y   ')


      if(nsel <= 0) then
        write(6,*)"No curves selected!"
        ierrs = 1000
        return
      endif

      if(nbuf >= size(nwert)) then
        write(6,*)"..too many buffers", nbuf, size(nwert)
        ierrs = 2000
        return
      endif

      write(*,'(a,f12.6)', advance='no')"Averaging with xcatch = ",xcatch
      if(.not. relative_catch) write(*,'(a)', advance='no')" absolute catch " 
      if(logy)                 write(*,'(a)')"   y-values on log-scale (log_y)"


! prepare destination
      nbuf = nbuf + 1
      call txfera(isels(1),nbuf)
      numor(nbuf)    = numor(isels(1))+200000
      if(relative_catch) then
        call parset("xcatch_r ",sngl(xcatch),nbuf)
      else
        call parset("xcatch_a ",sngl(xcatch),nbuf)
      endif
      if(logy) then
        call parset("log_y    ",1.0,nbuf)
      else
        call parset("log_y    ",0.0,nbuf)
      endif
      do i=1, nsel
        write(cbuf,'("c",i0)')i
        call parset(cbuf,real(numor(isels(i))),nbuf)
      enddo
! check how many data_points are there
      number_of_data_points = 0
      do i=1,nsel
       number_of_data_points = number_of_data_points + nwert(isels(i)) 
      enddo

      if(allocated(xva))     deallocate(xva)
      if(allocated(yva))     deallocate(yva)
      if(allocated(yvaer))   deallocate(yvaer)
      if(allocated(va_used)) deallocate(va_used)
      if(allocated(iperm))   deallocate(iperm)
      if(allocated(weights)) deallocate(weights)
      allocate(xva    (number_of_data_points))
      allocate(yva    (number_of_data_points))
      allocate(yvaer  (number_of_data_points))
      allocate(va_used(number_of_data_points))
      allocate(  iperm(number_of_data_points))
      allocate(weights(nsel))

! distribute data
      n = 0      
      do i=1,nsel
        do j=1,nwert(isels(i))
         n=n+1
         xva(n)    =  xwerte(j,isels(i))
         yva(n)    =  ywerte(j,isels(i))
         yvaer(n)  =  yerror(j,isels(i))
         if(logy) then
           if(yva(n) <= 0d0) then
             n = n - 1
             cycle
           endif
           yvaer(n) = abs(yvaer(n) / yva(n))
           yva(n)    = log(yva(n))  
         endif


         if(yvaer(n) <= 0.0d0) then
           yem = yvaer(max(1,n-1))
           if(yem > 0.0) then
             yvaer(n) = yem
             write(6,*)"WARNING: borrowed error from previous point",i,j
           else
             write(6,*)"Error occured with record:",isels(i)," point:",j
             call errsig(999,"ERROR: Average: Missing y-errors $")
             return
           endif
         endif
         va_used(n)= .false.
        enddo 
      enddo
      number_of_data_points = n
! error-sorting (using SLATEC dpsort) 
      iperm = [(i,i=1,number_of_data_points)]
      call  DPSORT (yvaer, number_of_data_points, IPERM,  1 , ier)
!C      KFLAG - input -- control parameter:---------------^
!C            =  2  means return the permutation vector resulting from
!C                  sorting DX in increasing order and sort DX also.
!C            =  1  means return the permutation vector resulting from
!C                  sorting DX in increasing order and do not sort DX.
!C            = -1  means return the permutation vector resulting from
!C                  sorting DX in decreasing order and do not sort DX.
!C            = -2  means return the permutation vector resulting from
!C                  sorting DX in decreasing order and sort DX also.

       if(iout() > 1 ) then
         write(6,*)            " All related points after error sorting: "
         write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,number_of_data_points)
       endif

! now proceed with the lowest error point and get all neibours in xcatch distance
       n_result_point = 0
d1:    do i=1,number_of_data_points
          if(va_used(iperm(i))) cycle d1
          va_used(iperm(i)) = .true.
          xaver = xva(iperm(i)) * 1d0/yvaer(iperm(i))**2
          yaver = yva(iperm(i)) * 1d0/yvaer(iperm(i))**2
          yerra = 1d0/yvaer(iperm(i))**2
d2:       do j=1,number_of_data_points
            if(va_used(iperm(j))) cycle d2
            if(relative_catch) then
              if(abs(2*(xva(iperm(j))-xva(iperm(i)))/(xva(iperm(j))+xva(iperm(j)))) > xcatch ) cycle d2 ! relative match
            else
              if(abs(2*(xva(iperm(j))-xva(iperm(i))))                               > xcatch ) cycle d2 ! absolute match
           endif
              xaver = xaver + xva(iperm(j)) * 1d0/yvaer(iperm(j))**2
              yaver = yaver + yva(iperm(j)) * 1d0/yvaer(iperm(j))**2
              yerra = yerra + 1d0/yvaer(iperm(j))**2   
              va_used(iperm(j)) = .true.           
          enddo d2 
          n_result_point = n_result_point + 1
          if(n_result_point > size(xwerte(:,1))) then
            write(6,*)'too many points created.'
            return
          endif
          xwerte(n_result_point,nbuf) = xaver / yerra
          ywerte(n_result_point,nbuf) = yaver / yerra
          yerror(n_result_point,nbuf) = sqrt(1d0 / yerra)

!          write(6,*)"test: ", n_result_point, xwerte(n_result_point,nbuf), &
!                     ywerte(n_result_point,nbuf), yerror(n_result_point,nbuf)

       enddo d1
       nwert(nbuf) =  n_result_point
!sort the result according to x
       iperm = [(i,i=1, n_result_point)]

!       write(6,*)            " intermediate result: "
!       write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,n_result_point)

       call  SPSORT (xwerte(:,nbuf),n_result_point, IPERM, 1 , ier)
       xva  (1:n_result_point) = xwerte(iperm(1:n_result_point),nbuf)
       yva  (1:n_result_point) = ywerte(iperm(1:n_result_point),nbuf)
       yvaer(1:n_result_point) = yerror(iperm(1:n_result_point),nbuf)
       if(logy) then
         yva(1:n_result_point)   = exp(yva(1:n_result_point))
         yvaer(1:n_result_point) = yva(1:n_result_point) * yvaer(1:n_result_point) 
       endif
       xwerte(1:n_result_point,nbuf) =  xva  (1:n_result_point)
       ywerte(1:n_result_point,nbuf) =  yva  (1:n_result_point) 
       yerror(1:n_result_point,nbuf) = yvaer (1:n_result_point)

       if(put_width) then
         write(6,'(a)')" ------ OPTION: putwidth ------ "
         xwidth = minval(abs(xva(2:n_result_point)-xva(1:n_result_point-1)))
         call parset ("_xwidth ",xwidth,nbuf)
         write(6,'(a,es14.5)')"the channel width _xwidth used by convoluted theories (kohl..) is ",xwidth
         write(6,'(a)')" the minimum distance between channels! "
         write(6,'(a)')" CARE: consider whether this is adeqaute here (for quaiselastic spectra like kohl...)"
         write(6,'(a)')"       with xcatch relative it is probably ok and good "
         write(6,'(a)')"       with xcatch absolute the new points it must be carefully considered what is good "
         write(6,'(a)')"       average is a kind of interpolation and not an integration, so the original width "
         write(6,'(a)')"       may be closer to what is needed...... "
         write(6,'(a)')" ------ "
      endif

!! determine contributed weights and average parameters accordingly
!! get weights:
! write(*,*)"TEST0: determine weights nsel=", nsel
      weights = 0
      do i=1,nsel
        do j=1,nwert(isels(i))
          if( yerror(j,isels(i)) > 0d0) weights(i) = weights(i) + 1d0/yerror(j,isels(i))**2
        enddo 
      enddo
      weights = weights / sum(weights(1:nsel))

! write(*,*)"TEST0a: weights", weights(1:nsel)

! write(*,*)"TEST0b: nopar", nopar(nbuf), nbuf, nsel

dipa: do ipa=1,nopar(nbuf)
        if(nopar(nbuf) <=0) exit
        xpaav = 0
!TP: write(*,*)"TEST ipa-loop: ",ipa, nopar(nbuf)
dse:    do i=1,nsel
!TP: write(*,*)"TEST1:a ",i, ipa, isels(i), trim(napar(ipa,isels(i))),nopar(isels(i))       
          cbuf = "_nomatch" 
!?          if(nopar(isels(i)) > ipa) cycle dse
          cbuf = trim(napar(ipa,isels(i)))   
!TP: write(*,*)"TEST:     cbuf >",trim(cbuf),"<",i,ipa, len_trim(cbuf) 
          call parget(cbuf, xpa, isels(i), ier )
!TP: write(*,*)"TEST: get cbuf >",trim(cbuf),"<", xpa, isels(i), ier

          if(ier.ne.0) cycle dipa
          xpaav = xpaav + weights(i) * xpa
!TP: write(*,*)"TEST1: ",i, isels(i), ipa, cbuf, xpa, xpaav,  weights(i), ier  
        enddo dse
        call parset(napar(ipa,isels(1)),xpaav,nbuf)
!TP: write(*,*)"TEST SETTING: ",trim(napar(ipa,isels(1))),":",trim(cbuf),":", xpaav, nbuf
      enddo dipa





       write(6,'(a)'    , advance='no') "final result from: "
       write(6,'(20i4)' , advance='no') (isels(i),i=1,nsel)
       if(relative_catch) then
         write(6,'(a,i4,a,f12.6,a)' )" --> ",nbuf,"    xcatch=",xcatch," relative"  
       else
         write(6,'(a,i4,a,es14.6,a)')" --> ",nbuf,"    xcatch=",xcatch," absolute" 
       endif
       if(iout() > 0) write(6,'(2i5,3e14.7)') (i,iperm(i),xva(iperm(i)),yva(iperm(i)),yvaer(iperm(i)),i=1,n_result_point)

       write(*,'(a)') "parameters have been combined with global error weightnfrom the selected records"

  
       nsel     = 1
       isels(1) = nbuf


      deallocate(xva)
      deallocate(yva)
      deallocate(yvaer)
      deallocate(va_used)
      deallocate(iperm)
      deallocate(weights)


   end subroutine average_data

 
   subroutine arit()
!     -------------------------
      use dimensions
      use new_com
      ! use cincoc
      use xoutxx
      use cdata
      ! use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
!      use constants
      use PhysicalConstantsPlus
 
      implicit none

      double precision       :: y2   (size(xwerte(:,1)))  ! selection values interpolated to x1 vals
      double precision       :: y2err(size(xwerte(:,1)))  ! selection values interpolated to x1 vals
      double precision, save :: f1 = 1d0                  ! factor for sel 1
      double precision, save :: f2 = 1d0                  ! factor for sel 2 
      integer                :: mode  
      integer                :: is1, is2  
      integer                :: number_of_datapoints, n

      integer                :: imatch1, imatch2
      integer                :: inew
      double precision       :: p
      integer                :: i, to_numor = 1

      integer, parameter     :: add  = 0
      integer, parameter     :: mult = 1
      integer, parameter     :: div  = 2
      integer, parameter     :: vid  = 3

      character(len=1)       :: op



      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= arit2 f1 <factor1> f2 <factor2> [to <numor>]  [options]                     '
       write(6,*)'=      adds the two selected records with factors 1,2 and stores as numor     '
       write(6,*)'=      for nonmatchig x-values interpolation is applied                       '
       write(6,*)'=      x-values of first selection serve as master                            '
       write(6,*)'=   options:                                                                  '
       write(6,*)'=            div     : divides (1)/(2) (instead of adding)                    '
       write(6,*)'=            vid     : divides (2)/(1) (instead of adding)                    '
       write(6,*)'=            mult    : multiplies (instead of adding)                         '
       write(6,*)'=============================================================================='
       return
      endif





      f1       = getval("f1      ",f1,inew)
      f2       = getval("f2      ",f2,inew)
      to_numor = intval("to      ",to_numor,inew)
                                      mode = add       ! +
      if(         found("mult    "))  mode = mult      ! * 
      if(         found("div     "))  mode = div       ! / 
      if(         found("vid     "))  mode = vid       ! / 

         
      if(nsel .ne. 2) then
        write(6,*)"Needs selection of exactly TWO records!"
        ierrs = 1000
        return
      endif

      is1 = isels(1)      
      is2 = isels(2)      
     


      if(nbuf >= size(nwert)) then
        write(6,*)"..too many buffers", nbuf, size(nwert)
        ierrs = 2000
        return
      endif

! prepare destination by pushing paramters etc.
      nbuf = nbuf + 1
      call txfera(is1,nbuf)
      numor(nbuf)        = numor(is1)+1000000

      if(   found("to      ")) numor(nbuf) = to_numor
      call parset("numor1  ",real(numor(is1)),nbuf)
      call parset("numor2  ",real(numor(is2)),nbuf)
! check how many data_points are there
      number_of_datapoints = nwert(is1)

! prepare selection 2 linear interpolation vector
      do i=1,nwert(is1)
       imatch1 = minloc(abs(xwerte(i,is1)-xwerte(1:nwert(is2),is2)),dim=1)
       imatch2 = minloc(abs(xwerte(i,is1)-xwerte(1:nwert(is2),is2)),dim=1, &
                        mask=(xwerte(1:nwert(is2),is2).ne.xwerte(imatch1,is2)))
       p        = (xwerte(i,is1)  -  xwerte(imatch1,is2)) /   &
                 (xwerte(imatch2,is2) - xwerte(imatch1,is2))
       y2(i)    = ywerte(imatch1,is2) * (1d0-p) + ywerte(imatch2,is2) * p
       y2err(i) = yerror(imatch1,is2) * (1d0-p) + yerror(imatch2,is2) * p 
       if(abs(p) > 1d0) y2err(i) = Huge(p)   
! write(*,*)"T1:",i, imatch1, imatch2, p, xwerte(i,is1), xwerte(imatch1,is2), xwerte(imatch2,is2)
! write(*,*)"  :",i, imatch1, imatch2, ywerte(i,is1), y2(i), ywerte(imatch1,is2), ywerte(imatch2,is2)
      enddo

! do the arithmetic
      select case(mode)
      case(add)
        xwerte(1:number_of_datapoints,nbuf)         = xwerte(1:number_of_datapoints,is1)
        ywerte(1:number_of_datapoints,nbuf)         =   f1 * ywerte(1:number_of_datapoints,is1)  &
                                                      + f2 * y2(1:number_of_datapoints)

        yerror(1:number_of_datapoints,nbuf)   = sqrt( (f1 * yerror(1:number_of_datapoints,is1))**2  &
                                                     +(f2 * y2err(1:number_of_datapoints))**2 )
        op = "+"
      case(mult)
        xwerte(1:number_of_datapoints,nbuf)         = xwerte(1:number_of_datapoints,is1)
        ywerte(1:number_of_datapoints,nbuf)         =   f1 * ywerte(1:number_of_datapoints,is1)  &
                                                      * f2 * y2(1:number_of_datapoints)

        yerror(1:number_of_datapoints,nbuf)   = f1*f2*sqrt( ( y2(1:number_of_datapoints)*yerror(1:number_of_datapoints,is1))**2  &
                                                     +( ywerte(1:number_of_datapoints,is1) * y2err(1:number_of_datapoints))**2 )
        op = "*"
      case(div)
        xwerte(1:number_of_datapoints,nbuf)         = xwerte(1:number_of_datapoints,is1)
        ywerte(1:number_of_datapoints,nbuf)         =     f1 * ywerte(1:number_of_datapoints,is1)  &
                                                      / ( f2 * y2(1:number_of_datapoints))

        yerror(1:number_of_datapoints,nbuf)   =(f1/f2)*sqrt( ( y2(1:number_of_datapoints)*yerror(1:number_of_datapoints,is1))**2  &
                                                     +( ywerte(1:number_of_datapoints,is1) * y2err(1:number_of_datapoints))**2 )  &
                                                     / y2(1:number_of_datapoints)**2
        op = "/"
      case(vid)
        xwerte(1:number_of_datapoints,nbuf)         = xwerte(1:number_of_datapoints,is1)
        ywerte(1:number_of_datapoints,nbuf)         = 1d0/  (  f1 * ywerte(1:number_of_datapoints,is1) )  &
                                                      * ( f2 * y2(1:number_of_datapoints))

        yerror(1:number_of_datapoints,nbuf)   =(f2/f1)*sqrt( ( y2(1:number_of_datapoints)*yerror(1:number_of_datapoints,is1))**2  &
                                                     +( ywerte(1:number_of_datapoints,is1) * y2err(1:number_of_datapoints))**2 )  &
                                                     / ywerte(1:number_of_datapoints,is1)**2
        op = "/"
        i  = is1;  is1 = is2; is2 = i
      case default
        write(6,*)"..unknown arithmetic operation requested", mode
        ierrs = 3000
        return
  
      end select

!  remove invalid points
      n = 0
      do i=1,number_of_datapoints
         if(y2err(i) < Huge(p)) then
           n = n + 1
           xwerte(n,nbuf) = xwerte(i,nbuf) 
           ywerte(n,nbuf) = ywerte(i,nbuf) 
           yerror(n,nbuf) = yerror(i,nbuf) 
         endif
      enddo
      nwert(nbuf) = n
      nsel = 1
      isels(1)    = nbuf

     write(*,'(a,i0,a,i0,a,i0)')"arit: f1 * (",is1,") "//op//" f2 * (",is2,")  ==> ",nbuf
     write(*,'(a,e15.8)')       "      f1 = ",f1
     write(*,'(a,e15.8)')       "      f2 = ",f2
     write(*,'(a,i0   )')       "      new number of points = ",nwert(nbuf)
     write(*,'(a,i0,a,i0)')     "      ==> ",nbuf," now selected, new numor = ",numor(nbuf)

   end subroutine arit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! intensity aligning two data sets by scale an offset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  align_y
 
      use dimensions
      use new_com
      ! use cincoc
      use xoutxx
      use cdata
      ! use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
!      use constants
      use PhysicalConstantsPlus
 
      implicit none
      double precision     :: mpar(2) = [0d0,1d0]
      double precision     :: errmin
      double precision     :: step(2) = 0.05d0
      integer              :: n1, n2, ier


      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= aligny                                                                     ='
       write(6,*)'= determine offset and scale y<=(y+offse)*scale                              ='
       write(6,*)'= of second selected record in order to match the first selected             ='
       write(6,*)'= the scaled result is copied to a new record                                ='
       write(6,*)'= this and the "master" record stay selected such tha one may                ='
       write(6,*)'= use average to merge them                                                  ='
       write(6,*)'= HINT:                                                                      ='
       write(6,*)'=     check whether clip (ing) marginal points in the overlap could          ='
       write(6,*)'=     improve the accuracy of the result !                                   ='
       write(6,*)'=  aligny offset <val>                                                       =' 
       write(6,*)'=  aligny yscale <val>                                                       =' 
       write(6,*)'=  aligny ostep  <val>                                                       ='  
       write(6,*)'=  aligny ystep  <val>                                                       ='  
       write(6,*)'=============================================================================='
       return
      endif


      if(nsel .ne. 2) then
        write(*,*)"ERROR: aligny needs selection of two curves"
        write(*,*)"Select just two curves: master and slave to be scaled!"
        ierrs = 1000
        return
      endif

      if(nbuf >= size(nwert)-1) then
        write(*,*)"..too many buffers", nbuf, size(nwert)
        ierrs = 2000
        return
      endif

      mpar(1) =  getval("offset ",mpar(1))
      mpar(2) =  getval("yscale ",mpar(2))
      step(1) =  getval("ostep  ",step(1)) 
      step(2) =  getval("ystep  ",step(2)) 

! prepare destination
      nbuf = nbuf + 1
      call txfera(isels(2),nbuf)
      numor(nbuf)    = numor(isels(2))+300000

 !     mpar = [0d0, 1d0] 
      n1   = nwert(isels(1))
      n2   = nwert(isels(2))
      call smatch(xwerte(1:n1,isels(1)),ywerte(1:n1,isels(1)),yerror(1:n1,isels(1)),n1, &
                  xwerte(1:n2,isels(2)),ywerte(1:n2,isels(2)),yerror(1:n2,isels(2)),n2, &
                  mpar, step, errmin )
      write(*,'("offset=",f12.6," scale=",f12.6,"   residual error=",f12.6)')mpar(1:2), errmin

      xwerte(:,nbuf) =  xwerte(:,isels(2))
      ywerte(:,nbuf) = (ywerte(:,isels(2))+mpar(1))*mpar(2)
      yerror(:,nbuf) =  yerror(:,isels(2))         *mpar(2)
      nwert(nbuf)    =  n2
      call parset("refnum  ",real(numor(isels(1))),nbuf)
      call parset("aligny_o",sngl(mpar(1)),nbuf)
      call parset("aligny_s",sngl(mpar(2)),nbuf)
      isels(2) = nbuf

      write(*,'(a,i6)') "Scaled data written to record ",nbuf
      write(*,'(a,2i6)')"Selected are now              ",isels(1:2)
      write(*,'(a)')    "you may use  -> average  to merge them!"
      call setudf("aligny_o ",mpar(1),ier)
      call setudf("aligny_s ",mpar(2),ier)

end subroutine  align_y





subroutine smatch(x1,y1,y1err,n1, x2,y2,y2err,n2, mpar, step, errmin )
  use new_com
  implicit none
  real            , intent(in ) :: x1(n1)      ! x-values of master record
  real            , intent(in ) :: y1(n1)      ! y-values of master record
  real            , intent(in ) :: y1err(n1)   ! y-error  of master record
  integer         , intent(in ) :: n1          ! nuber of values in master
  real            , intent(in ) :: x2(n2)      ! x-values of slave record
  real            , intent(in ) :: y2(n2)      ! y-values of slave record
  real            , intent(in ) :: y2err(n2)   ! y-error  of slave record
  integer         , intent(in ) :: n2          ! nuber of values in slave
  double precision, intent(inout) :: mpar(2)     ! offset and scale ....
  double precision, intent(in ) :: step(2)     ! offset and scale step ....
  double precision, intent(out) :: errmin      ! residual matching error
  
  double precision :: reqmin = 1d-7
  integer          :: konvge = 1
  integer          :: kcount = 10000
  integer          :: icount
  integer          :: numres
  integer          :: ifault
  double precision :: mpar_min(size(mpar))

! --> is mainly for testing
  reqmin = getval("reqmin ",reqmin)
  konvge = intval("konvge ",konvge)
  kcount = intval("kcount ",kcount)
  write(*,*)"TESTPARA: ", reqmin, konvge, kcount


  errmin = match_err(mpar, size(mpar))

  call nelmin ( match_err, size(mpar), mpar, mpar_min, errmin, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

  mpar = mpar_min


contains 

function match_err(sp, n ) result(val) 
  implicit none
  double precision, intent(in ) :: sp(n)
  integer         , intent(in ) :: n 

  double precision :: val
  double precision :: yinterp, yinterp_err, xt, yt, yte, p
  integer          :: i, ma, mb, ncompare

  !! go through the points of slave and compare with master
  val      =  0
  ncompare = 0
ds1: do i=1,n2
       xt  = x2(i)
       yt  = (y2(i) + sp(1))* sp(2)
       yte =  y2err(i)      * sp(2)
       ma = minloc(abs(x1(1:n1)-xt),dim=1)   
       mb = min(ma + 1,n1)
       if(mb>n1) mb = ma-1
       p             = (xt - x1(ma))/(x1(mb) - x1(ma))

       if(abs(p) > 1d0) cycle ds1


       yinterp      = (1-p) * y1(ma) + p * y1(mb)
       yinterp_err  = sqrt(  ((1d0-p)*y1err(ma))**2  +((p)*y1err(mb))**2 )
       ncompare = ncompare + 1
       val = val + sqrt( (yinterp-yt)**2 / (yinterp_err**2 + yte**2) )
     enddo ds1
 
     if(ncompare > 0) then
       val = val / ncompare
     else
       val = Huge(val)
     endif
 

end function match_err

end subroutine smatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!>>> TBD: not yet functional align_a  must be reworked !!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>>> TBD: not yet functional align_a  must be reworked !!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>>> TBD: not yet functional align_a  must be reworked !!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>>> TBD: not yet functional align_a  must be reworked !!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  align_a
 
      use dimensions
      use new_com
      ! use cincoc
      use xoutxx
      use cdata
      ! use outlev
      use theory
      use selist
      use fslist
      use theorc
      use thparc
      use formul
      use cfc
      use cfunc
      use cfunce
      use partran
      use wlntran
      use sqtran
!      use constants
      use PhysicalConstantsPlus
 
      implicit none
      double precision     :: errmin
      double precision,save:: mpar(2)=1d0, step(2)=0.1d0
      integer              :: n1, n2, ier


      if(found('help    ')) then 
       write(6,*)'=============================================================================='
       write(6,*)'= aligna        !!! NOT YET FUNCTIONAL, CONTAINS BUGS or ALGORITH MUST BE MODIFIED ='
       write(6,*)'= determine x and y scales                                                   ='
       write(6,*)'= of second selected record in order to match the first selected             ='
       write(6,*)'= the scaled result is copied to a new record                                ='
       write(6,*)'= this and the "master" record stay selected such tha one may                ='
       write(6,*)'= use average to merge them                                                  ='
       write(6,*)'= USAGE: aligna xscale <startval> yscale <startval> xstep <val> ystep <val>  ='
       write(6,*)'=============================================================================='
       return
      endif

 stop "THIS ROUTINE aligna is not yet ready, must be reworked, befor new testsb remove this msg"    

      if(nsel .ne. 2) then
        write(*,*)"ERROR: aligna needs selection of two curves"
        write(*,*)"Select just two curves: master and slave to be scaled!"
        ierrs = 1000
        return
      endif

      if(nbuf >= size(nwert)-1) then
        write(*,*)"..too many buffers", nbuf, size(nwert)
        ierrs = 2000
        return
      endif

      mpar(1) =  getval("xscale ",mpar(1))
      mpar(2) =  getval("yscale ",mpar(2))
      step(1) =  getval("xstep  ",step(1)) 
      step(2) =  getval("ystep  ",step(2)) 

! prepare destination
      nbuf = nbuf + 1
      call txfera(isels(2),nbuf)
      numor(nbuf)    = numor(isels(2))+300000

      n1   = nwert(isels(1))
      n2   = nwert(isels(2))
      call amatch(xwerte(1:n1,isels(1)),ywerte(1:n1,isels(1)),yerror(1:n1,isels(1)),n1, &
                  xwerte(1:n2,isels(2)),ywerte(1:n2,isels(2)),yerror(1:n2,isels(2)),n2, &
                  mpar, step, errmin )
      write(*,'("xscale=",f12.6," yscale=",f12.6,"   residual error=",f12.6)')mpar(1:2), errmin

      xwerte(:,nbuf) =  xwerte(:,isels(2))
      ywerte(:,nbuf) = (ywerte(:,isels(2))+mpar(1))*mpar(2)
      yerror(:,nbuf) =  yerror(:,isels(2))         *mpar(2)
      nwert(nbuf)    =  n2
      call parset("refnum  ",real(numor(isels(1))),nbuf)
      call parset("aligna_o",sngl(mpar(1)),nbuf)
      call parset("aligna_s",sngl(mpar(2)),nbuf)
      isels(2) = nbuf

      write(*,'(a,i6)') "Scaled data written to record ",nbuf
      write(*,'(a,2i6)')"Selected are now              ",isels(1:2)
      write(*,'(a)')    "you may use  -> average  to merge them!"
      call setudf("aligna_o ",mpar(1),ier)
      call setudf("aligna_s ",mpar(2),ier)

end subroutine  align_a


subroutine amatch(x1,y1,y1err,n1, x2,y2,y2err,n2, mpar,step, errmin )
  implicit none
  real            , intent(in ) :: x1(n1)      ! x-values of master record
  real            , intent(in ) :: y1(n1)      ! y-values of master record
  real            , intent(in ) :: y1err(n1)   ! y-error  of master record
  integer         , intent(in ) :: n1          ! nuber of values in master
  real            , intent(in ) :: x2(n2)      ! x-values of slave record
  real            , intent(in ) :: y2(n2)      ! y-values of slave record
  real            , intent(in ) :: y2err(n2)   ! y-error  of slave record
  integer         , intent(in ) :: n2          ! nuber of values in slave
  double precision, intent(inout) :: mpar(2)     ! offset and scale ....
  double precision, intent(in ) :: step(2)     ! offset and scale ....
  double precision, intent(out) :: errmin      ! residual matching error
 
  double precision :: mpar_min(size(mpar))

  double precision :: reqmin = 1d-7
  integer          :: konvge = 1
  integer          :: kcount = 10000
  integer          :: icount
  integer          :: numres
  integer          :: ifault


  errmin = amatch_err(mpar, size(mpar))

  call nelmin ( amatch_err, size(mpar), mpar, mpar_min, errmin, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

  mpar = mpar_min


contains 

function amatch_err(sp, n ) result(val) 
  implicit none
  double precision, intent(in ) :: sp(n)
  integer         , intent(in ) :: n 

  double precision :: val
  double precision :: yinterp, yinterp_err, xt, yt, yte, p
  integer          :: i, ma, mb, ncompare

  !! go through the points of slave and compare with master
  val      =  0
  ncompare = 0
ds1: do i=1,n2
       xt  = x2(i)          *  sp(1)
       yt  = (y2(i))        *  sp(2)
       yte =  y2err(i)      *  sp(2)
       ma = minloc(abs(x1(1:n1)-xt),dim=1)            ; write(*,'(a,i8,2e14.7,2i8)')"ma: ",i,x1(ma),xt,n1,n2
       mb = min(ma + 1,n1)
       if(mb>n1) mb = ma-1
       p             = (xt - x1(ma))/(x1(mb) - x1(ma))   ; write(*,'(a,e14.7)')"p: ",p
       if(abs(p) > 2d0) cycle ds1
       yinterp      = (p-1d0) * y1(ma) + p * y1(mb)      ; write(*,'(a,3e14.7)')"yinterp: ",yinterp, y2(i), yt
 
       yinterp_err  = sqrt(  ((p-1d0)*y1err(ma))**2  +(p*y1err(mb))**2 )
       ncompare = ncompare + 1
       val = val - 1d0 / ((yinterp-yt)**2 / (yinterp_err**2 + yte**2))
write(*,'(a,3i5,7e15.7)')"TP1: ",i,ma,mb,xt,x1(ma),x1(mb),yt,yinterp, yinterp_err,val
     enddo ds1
 
     if(ncompare > 0) then
       val = val / ncompare
     else
       val = Huge(val)
     endif
 
!!TP: write(*,'(3f12.6)') sp, val

end function amatch_err

end subroutine amatch

!<<< TBD: not yet functional align_a  must be reworked !!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Suroutines from external sources                                          !!
!!!! ASA047                                                                    !!
!!!! Nelder-Mead Minimization                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

!*****************************************************************************
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, n )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
  real ( kind = 8 ) del
  real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
  real ( kind = 8 ), parameter :: eps = 0.001D+00
  real ( kind = 8 ), external :: fn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcount
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) l
  integer ( kind = 4 ) numres
  real ( kind = 8 ) p(n,n+1)
  real ( kind = 8 ) p2star(n)
  real ( kind = 8 ) pbar(n)
  real ( kind = 8 ) pstar(n)
  real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) rq
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) y(n+1)
  real ( kind = 8 ) y2star
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ynewlo
  real ( kind = 8 ) ystar
  real ( kind = 8 ) z
!
!  Check the input parameters.
!
  if ( reqmin <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( n < 1 ) then
    ifault = 1
    return
  end if

  if ( konvge < 1 ) then
    ifault = 1
    return
  end if
!
!  Initialization.
!
  icount = 0
  numres = 0
  jcount = konvge
  del = 1.0D+00
  rq = reqmin * real ( n, kind = 8 )
!
!  Initial or restarted loop.
!
  do

    p(1:n,n+1) = start(1:n)
    y(n+1) = fn ( start, n )
    icount = icount + 1
!
!  Define the initial simplex.
!
    do j = 1, n
      x = start(j)
      start(j) = start(j) + step(j) * del
      p(1:n,j) = start(1:n)
      y(j) = fn ( start, n )
      icount = icount + 1
      start(j) = x
    end do
!
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
!
    ilo = minloc ( y(1:n+1), 1 )
    ylo = y(ilo)
!
!  Inner loop.
!
    do while ( icount < kcount )
!
!  YNEWLO is, of course, the HIGHEST value???
!
      ihi = maxloc ( y(1:n+1), 1 )
      ynewlo = y(ihi)
!
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
!
      do i = 1, n
        pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
      end do
!
!  Reflection through the centroid.
!
      pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
      ystar = fn ( pstar, n )
      icount = icount + 1
!
!  Successful reflection, so extension.
!
      if ( ystar < ylo ) then

        p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
        y2star = fn ( p2star, n )
        icount = icount + 1
!
!  Retain extension or contraction.
!
        if ( ystar < y2star ) then
          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
        else
          p(1:n,ihi) = p2star(1:n)
          y(ihi) = y2star
        end if
!
!  No extension.
!
      else

        l = 0
        do i = 1, n + 1
          if ( ystar < y(i) ) then
            l = l + 1
          end if
        end do

        if ( 1 < l ) then

          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
!
!  Contraction on the Y(IHI) side of the centroid.
!
        else if ( l == 0 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
          y2star = fn ( p2star, n )
          icount = icount + 1
!
!  Contract the whole simplex.
!
          if ( y(ihi) < y2star ) then

            do j = 1, n + 1
              p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
              xmin(1:n) = p(1:n,j)
              y(j) = fn ( xmin, n )
              icount = icount + 1
            end do

            ilo = minloc ( y(1:n+1), 1 )
            ylo = y(ilo)

            cycle
!
!  Retain contraction.
!
          else
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          end if
!
!  Contraction on the reflection side of the centroid.
!
        else if ( l == 1 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
          y2star = fn ( p2star, n )
          icount = icount + 1
!
!  Retain reflection?
!
          if ( y2star <= ystar ) then
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          else
            p(1:n,ihi) = pstar(1:n)
            y(ihi) = ystar
          end if

        end if

      end if
!
!  Check if YLO improved.
!
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      end if

      jcount = jcount - 1

      if ( 0 < jcount ) then
        cycle
      end if
!
!  Check to see if minimum reached.
!
      if ( icount <= kcount ) then

        jcount = konvge

        x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
        z = sum ( ( y(1:n+1) - x )**2 )

        if ( z <= rq ) then
          exit
        end if

      end if

    end do
!
!  Factorial tests to check that YNEWLO is a local minimum.
!
    xmin(1:n) = p(1:n,ilo)
    ynewlo = y(ilo)

    if ( kcount < icount ) then
      ifault = 2
      exit
    end if

    ifault = 0

    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn ( xmin, n )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) - del - del
      z = fn ( xmin, n )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) + del
    end do

    if ( ifault == 0 ) then
      exit
    end if
!
!  Restart the procedure.
!
    start(1:n) = xmin(1:n)
    del = eps
    numres = numres + 1

  end do

  return
end


  ! -------------------------------------------------------------------
  !> replace substring by a new string
  subroutine replace_str(string,substr,newstring)
    character(len=*), intent(inout) :: string
    character(len=*), intent(in)    :: substr
    character(len=*), intent(in)    :: newstring
    !
    integer   :: pos1, slen, sublen, cadd, mlen

    sublen  = len_trim(substr)
    cadd    = len_trim(newstring)-sublen
    mlen    = len(string)

    lp: do
        pos1    = index(string, trim(substr))
        if(pos1 <= 0 )       return
        slen    = len_trim(string)
        if(slen+cadd > mlen) return
        string  = string(1:pos1-1)//trim(newstring)//string(pos1+sublen:)
    enddo lp
  end subroutine replace_str


  ! =========================================================================
  !> another one, this time it is a function and trims empty space in between
  function string_replace(string, substring1, substring2) result(new_string)
    character(len=*), intent(in)         :: string
    character(len=*), intent(in)         :: substring1
    character(len=*), intent(in)         :: substring2
    character(len=len(string)+(len(string)/len(substring1)+1)*len(substring2)) :: new_string
    integer  :: i, k, l, i1, i2, len0, len1, len2

    l    = len(string)
    len0 = len_trim(string)
    len1 = len_trim(substring1)
    len2 = len_trim(substring2)

    new_string = ' '
    i1 = 1
    i2 = l

    do k=1,l
        i = index(string(i1:l),substring1(1:len1))
        if(i == 0) then
            new_string = trim(new_string)//string(i1:l)
            exit
        endif

        i = i + i1 - 1
        new_string = trim(new_string)//string(i1:i-1)//substring2(1:len2)
        i1 = i+len1
    enddo

    new_string = trim(new_string)
    return
  end function string_replace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! help functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine command_list()
  write(6,*) " command name               " , " parameters: "
  write(6,*) "?         (synonym) help    " , " " 
  write(6,*) "acl       (synonym) aclast  " , " "
  write(6,*) "activate  (synonym) ac      " , " <theory_name> "
  write(6,*) "activlst  (synonym) al      " , " "
  write(6,*) "addsels                     " , " "
  write(6,*) "addval                      " , " "
  write(6,*) "aligny                      " , " determine scale to match y-values of two records "
  write(6,*) "argvals                     " , " "
  write(6,*) "arit                        " , " f1 <factor1>  f2 <factor2> to <numor> [div | mult | (no)norm ]"
  write(6,*) "arit2                       " , " f1 <factor1>  f2 <factor2> to <numor> [div | mult | (no)norm ]"
  write(6,*) "average                     " , " "
  write(6,*) "chgthpar                    " , " "
  write(6,*) "clip                        " , " "
  write(6,*) "close                       " , " "
  write(6,*) "combine                     " , " "
  write(6,*) "copy                        " , " "
  write(6,*) "couple                      " , " "
  write(6,*) "creatser                    " , " "
  write(6,*) "cs        (synonym) clearsel" , " "
  write(6,*) "des       (synonym) desmear " , " "
  write(6,*) "desactiv  (synonym) dac     " , " "
  write(6,*) "detsens                     " , " "
  write(6,*) "dir                         " , " "
  write(6,*) "dispsel   (synonym) dsl     " , " "
  write(6,*) "dmx                         " , " "
  write(6,*) "echocurv                    " , " "
  write(6,*) "edit                        " , " "
  write(6,*) "ercorrc   (synonym) ecc     " , " "
  write(6,*) "fft                         " , " "
  write(6,*) "fftmx     (synonym) fft-ms  " , " "
  write(6,*) "fftrand                     " , " "
  write(6,*) "fit                         " , " "
  write(6,*) "fit_plus                    " , " "
  write(6,*) "fun       (synonym) function" , " "
  write(6,*) "funfun                      " , " "
  write(6,*) "fxy                         " , " & x=(formel) y=(formel)"
  write(6,*) "ga_fit                      " , " "
  write(6,*) "gaiscale                    " , " "
  write(6,*) "gdos_q                      " , " "
  write(6,*) "gen.res                     " , " "
  write(6,*) "get_th  (synonym) gth       " , " "
  write(6,*) "gplot      (synonym) gp     " , " "
  write(6,*) "hiqextra                    " , " "
  write(6,*) "ia        (synonym) i-abso  " , " "
  write(6,*) "in        (synonym) input   " , " "
  write(6,*) "inscn                       " , " "
  write(6,*) "integrat                    " , " "
  write(6,*) "interpol                    " , " "
  write(6,*) "invers                      " , " "
  write(6,*) "kz                          " , " "
  write(6,*) "label                       " , " "
  write(6,*) "m         (synonym) mirror  " , " "
  write(6,*) "msave                       " , " "
  write(6,*) "mux                         " , " "
  write(6,*) "noise                       " , " "
  write(6,*) "numorchg                    " , " "
  write(6,*) "numorpls                    " , " "
  write(6,*) "fopen                       " , " "
  write(6,*) "out_gli   (synonym) gli     " , " "
  write(6,*) "paraout                     " , " "
  write(6,*) "parextra                    " , " "
  write(6,*) "oplot    plot to ocular     " , " "
  write(6,*) "plot      (synonym) p       " , " "
  write(6,*) "plot0     (synonym) p0      " , " "
  write(6,*) "purge                       " , " "
  write(6,*) "prompt                      " , " "
  write(6,*) "putpar                      " , " "
  write(6,*) "parlev                      " , " "
  write(6,*) "qc        (synonym) q-conv  " , " "
  write(6,*) "rebin                       " , " "
  write(6,*) "recsum                      " , " "
  write(6,*) "rename                      " , " "
  write(6,*) "rerange                     " , " "
  write(6,*) "save                        " , " "
  write(6,*) "scale                       " , " "
  write(6,*) "sel                         " , " "
  write(6,*) "sequence                    " , " "
  write(6,*) "serfit                      " , " "
  write(6,*) "seterr                      " , " "
  write(6,*) "sp        (synonym) spline  " , " "
  write(6,*) "swapxy                      " , " "
  write(6,*) "sym                         " , " "
  write(6,*) "th_init                     " , " "
  write(6,*) "thc                         " , " "
  write(6,*) "theos                       " , " "
  write(6,*) "title     (synonym) tit     " , " "
  write(6,*) "tofdos                      " , " "
  write(6,*) "tracorr                     " , " "
  write(6,*) "uni_ft                      " , " "
  write(6,*) "write                       " , " "
  write(6,*) "xformel                     " , " "
  write(6,*) "yfitform                    " , " "
  write(6,*) "yformel                     " , " "
  write(6,*) "z         (synonym) zero    " , " "
  write(6,*) "zero                        " , " "
end subroutine command_list
