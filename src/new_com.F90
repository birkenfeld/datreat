!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Adaption of incom using newer Fortran ( M.M.)
!! This is to
!! - make a MODULE of incom
!! - hide internals of the parser form the rest of the prog.
!! - allow for long names
!! - stay compatible (possibly then some private declaratuion must be public...)
!! - it is not yet fully tested
!! - it was just a reasonable method for me to learn a bit of the new fortran features...
!!
!! - Piotr Zolnierczuk (10/21/2015) - added linenoise interface
!!   for toplevel input, see https://github.com/antirez/linenoise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! comment the line below if you do not want/have linenoise
#define USE_LINENOISE


MODULE new_com
  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit

!#ifdef __INTEL_COMPILER_BUILD_DATE
!  use ifport            ! ifortran specific to improve system communication (form of system call, open...)
!#endif

#ifdef USE_LINENOISE
  use iso_c_binding, only : c_int, c_char, c_null_char
#endif

  use os_utils

  PUBLIC

  ABSTRACT INTERFACE
     subroutine sub_usr_func (name, x, nx, ier)
       character(len=*), intent(in)    :: name
       double precision, intent(out)   :: x(*)
       integer,          intent(inout) :: nx
       integer,          intent(out)   :: ier
     END subroutine sub_usr_func
  END INTERFACE

  ABSTRACT INTERFACE
     subroutine sub_usr_extr (name, x, ier)
       character(len=*), intent(in)  :: name
       double precision, intent(out) :: x
       integer,          intent(out) :: ier
     END subroutine sub_usr_extr
  END INTERFACE


  interface command
     module procedure  iscmd, iscmdc
  end interface

  interface get_named_value
     module procedure  getval, intval, chrval
  end interface


#ifdef USE_LINENOISE
  !----------------------------------------------------------------------------
  ! linenoise interface
  !----------------------------------------------------------------------------
  interface
     subroutine fclinenoise_init(filename) bind(c, name='fclinenoise_init')
       use iso_c_binding, only : c_char
       character(kind=c_char),intent(in)   ::  filename(*)
     end subroutine fclinenoise_init
  end interface

  interface
     function fclinenoise_add(command) bind(c, name='fclinenoise_add')
       use iso_c_binding, only : c_char, c_int
       integer(kind=c_int):: fclinenoise_add
       character(kind=c_char),intent(in)   ::  command(*)
     end function fclinenoise_add
  end interface

  interface
     ! internal wrapper around read
     function fclinenoise_read(buflen,buf,prompt) bind(c, name='fclinenoise_read')
       use iso_c_binding, only : c_int, c_char
       integer(kind=c_int):: fclinenoise_read
       integer(kind=c_int),intent(in),value ::  buflen
       character(kind=c_char),intent(out)   ::  buf(*)
       character(kind=c_char),intent(in)    ::  prompt(*)
     end function fclinenoise_read
  end interface
  !----------------------------------------------------------------------------
#endif



  SAVE


  procedure(sub_usr_func), pointer, private  :: newcom_usrfun  => null()
  procedure(sub_usr_extr), pointer, private  :: newcom_usrextr => null()


  integer, parameter,public :: cmd_len       = 256
  integer, parameter,public :: cmd_long      = 1024

  integer, parameter :: minc          = 100
  integer, parameter :: mdepth        = 20
  integer, parameter :: musevar       = 500
!  integer, parameter :: mbuf          = cmd_long+200
  integer, parameter :: maxformlength = cmd_long
  integer, parameter :: maxitemlength = cmd_len
  integer, parameter :: maxnumstack   = 50
  integer, parameter :: maxopstack    = 50
  integer, parameter :: musrfstack    = 500
  integer, parameter :: nodelims      = 7


  integer, parameter :: history_depth = 1000
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


  double precision   , public ::  rpar(minc)
  integer         ::  iitems
  integer, public ::  inames
  integer, public ::  ipars
  integer, public ::  ioldc
  integer, public ::  inpar(minc)
  integer, public ::  iparn(minc)
  integer, public ::  inapa(minc)
  integer, public ::  iargs
  integer, public ::  ipmls
  integer, public ::  iolbuf !Can we delete this and put in a local var?
  integer, public ::  lstpar
  integer, public ::  lstnam
  integer, public ::  last_item

  ! ---- the current command is stored on comand
  !      the names on stack vname(*)
  !      the numbers on stack rpar(*)
  !      in common/cincom/
  !
  !character(cmd_long)   , public ::  mybuf  !! fuer re_scan testing ! (PAZ) not-used
  character(cmd_len)    , public ::  comand            ! was 8
  character(cmd_len)    , public ::  vname(minc)       ! was 8
  character(cmd_len)             ::  cmd_item(minc)    ! was 8
  logical                        ::  item_used(minc)   
  integer                        ::  item_nr_name(minc) !  i-th name associated with this item
  integer                        ::  name_nr_item(minc) !  i-th item associated with this (j-th) name
  integer                        ::  item_nr_rpar(minc) !  i-th rpar associated with this item
  integer                        ::  rpar_nr_item(minc) !  i-th item associated with this (r-par)
  character(cmd_long)   , public ::  title
  character(cmd_long)   , public ::  reslin
  character(cmd_long)   , public ::  inline
  character(cmd_len)    , public ::  arglst(minc)
  !character(cmd_len)    , public ::  pmlist(minc,2) ! (PAZ) not-used
  character(cmd_long)   , public ::  rlbuf

   character(cmd_len)   , public  :: top_macro_file = " "
   character(cmd_long)  , public  :: top_macro_call = " "
   character(cmd_long)  , public  :: raw_input      = " "
 
  ! path variables

  character(len=cmd_long) :: data_path  = './'
  character(len=cmd_long) :: save_path  = './'
  character(len=cmd_long) :: makro_path = './'
  character(len=cmd_long) :: mxx_path   = './'
  character(len=cmd_long) :: init_file_path  = ' ' ! default no auto init
#ifdef USE_LINENOISE
  character(cmd_long), public  :: history_file = ".history"
#else
  character(cmd_long), private :: history(0:history_depth)
#endif
  character(cmd_len),  public :: prompt = "drspine --> "

  logical, private :: cused, vused(minc), rused(minc), seten
  integer          :: info
  !
  !       common/cincom/rpar(minc),inames,ipars &
  !        ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf &
  !        ,lstpar, lstnam, cused, vused(minc), rused(minc), seten, info


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



  character(cmd_long) , public ::  argvals(minc)
  character(cmd_len)  , public ::  pmlst(mdepth,minc,2)

  integer       , public ::  iargvs
  integer       , public ::  ipmlst(mdepth)
  integer       , public ::  kanal(0:mdepth)
  integer       , public ::  ktop
  !data kanal/ 5,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58, 59/
  !data ktop /0/


  integer, private :: iot=0
  integer, public :: ierrr=0


  double precision  , private :: xyorig(3)
  double precision  , private :: rotvec(3)
  data xyorig/3*0.d0/, rotvec/3*0.d0/



  double precision  ,  public      ::  useval(musevar)
  integer, public      ::  nousev
  character*16, public :: usenam(musevar)
  ! ---- outputlevel
  ! xxxx,yyyy = aktuelle addersse fuer wertextraktion xx,yy
  !
  !integer, private  :: iout_nc = 0 ! (PAZ) not-used
  integer, public  :: ibild = -2   ! (PAZ) not-used
  integer, public  :: ierrs
  !integer, private  :: inka = 5    ! (PAZ) not-used
  integer, public  :: iibuf       ! (PAZ) not-used
  real, public ::  xxxx     ! (PAZ) not-used
  real, public ::  yyyy     ! (PAZ) not-used
  real, public ::  ptxf(20) ! (PAZ) not-used
  real, public ::  yyee     ! (PAZ) not-used




  !! evaluator data

  double precision   , private ::  numstack(maxnumstack)
  double precision   , private ::  degree=1.d0
  double precision   , private ::  valnum
  integer, private ::  priostack(0:maxopstack)
  integer, private ::  topnumstack
  integer, private ::  topopstack
  integer, private ::  tusrfstack
  integer, private ::  klammerprio
  integer, private ::  actchar
  integer, private ::  len_form
  integer, private ::  litem
  logical, private ::  ok
  logical, private ::  error
  logical, private ::  say=.false.

  character*1 , private ::  formula(0:maxformlength)
  character*1 , private ::  item(0:maxitemlength)
  character*1 , private ::  delims(0:nodelims)
  character*4 , private ::  typ
  character*4 , private ::  opstack(maxopstack)

  character(cmd_len), private ::  usrfstack(musrfstack)

  !needed in subroutine getitem - dirty!
  character*(maxitemlength+1), private :: citemx
  character*(maxformlength+1), private :: cformulax
  equivalence(citemx,item(0))               ! do we really need this maybe fix
  equivalence(cformulax,formula(0))         ! "   --> re_scan


  integer :: isignal  ! utilux communication --> rather make utilux a module
  common/sig/isignal


  ! constants

  !private :: dummy_usrfun  !(PAZ) not-used
  !private :: dummy_usrextr !(PAZ) not-used
  public  :: set_usrfunc
  public  :: set_usrextr
  public  :: set_init_file
  public  :: get_newcom
  public  :: creplace
  private :: re_scan
  public  :: parse
  public  :: errsig
  ! private :: pushn ! (PAZ) not-used
  public  :: get1
  public  :: getvec
  public  :: getnva
  public  :: rotavc
  public  :: unrot
  public  :: rotax
  public  :: matmux
  public  :: matvec
  public  :: mattxp
  public  :: matcpy
  public  :: matone
  public  :: axrot
  public  :: cross
  public  :: norm
  public  :: foinit
  public :: getitem
  public :: evaluate
  public :: cappend
  public  :: extract
  public  :: setudf
  public  :: clrudf
  public  :: shwudf
  public  :: settit

  public ::  copy_chars
  public  ::  compar
  public  ::  X_datpath
  public ::  savepath
  public  ::  iout
  public  ::  getval
  public  ::  valnxt
  public  ::  intval
  public  ::  intnxt
  public  ::  chrval
  public  ::  chrnxt
  public  ::  vnamef
  public  ::  titlef
  public  ::  rparf
  public  ::  ifound
  public  ::  found
  public  ::  folgt
  public ::  intvno
  public ::  getvno
  public  ::  laenge
  public  ::  lclen
  public ::  compare
  public  ::  inpaf
  public  ::  chrval8
  public  ::  cmdf
  public ::  iscmd
  public ::  iscmdc
  public  ::  iparf
  public  ::  inamf
  public  ::  inapf
  !
  public  ::  lvused ! whether the parameter has already been "used"


  public :: close_all_macros

  public  :: intnva

  public  :: u2c_datum
  public  :: set_prompt, get_prompt

CONTAINS


!!$  subroutine dummy_usrfun( name, x, nx, ier )
!!$    !   ===========================================
!!$    implicit none
!!$    character(len=*),          intent(in)    :: name
!!$    double precision,      intent(out)   :: x(*)
!!$    integer,               intent(inout) :: nx
!!$    integer,               intent(out)   :: ier
!!$
!!$    write(6,*)'dummy_usrfun ', name, nx
!!$    x(1) = 0d0
!!$    ier  = 0
!!$
!!$  end subroutine dummy_usrfun
!!$
!!$
!!$
!!$  subroutine dummy_usrextr( name, x, ier )
!!$    !   ========================================
!!$    implicit none
!!$    character(len=*),          intent(in)    :: name
!!$    double precision,      intent(out)   :: x
!!$    integer,               intent(out)   :: ier
!!$
!!$    write(6,*)'dummy_usrextr ', name
!!$    x   = 0d0
!!$    ier = 0
!!$
!!$  end subroutine dummy_usrextr


  subroutine set_usrfunc( fp )
    !   ============================
    implicit none
    procedure(sub_usr_func) :: fp
    newcom_usrfun => fp
  end subroutine set_usrfunc

  subroutine set_usrextr( fp )
    !   ============================
    implicit none
    procedure(sub_usr_extr) :: fp
    newcom_usrextr => fp
  end subroutine set_usrextr

  character(1024) function copy_chars(carray,n)
    !   ===========================
    implicit none
    integer, parameter :: m=1024
    !character(len=m)    ::  buffer
    character(len=1), dimension(30), intent(in) :: carray
    integer :: i,n
    do i=1,n
       !write(6,*)i,carray(i)
       if(carray(i)==' ') exit
       copy_chars(i:i) = carray(i)
    enddo


  end function copy_chars


  logical FUNCTION compar( s1, s2 )

    implicit none
    character(*), intent(in) :: s1, s2

    compar = (trim(s1) .eq. trim(s2))

  END FUNCTION compar


#ifdef USE_LINENOISE
  function add_tab_expansion(cmd) result(istat)
    character(len=*), intent(in) :: cmd
    integer :: istat
    istat = fclinenoise_add(trim(cmd)//c_null_char)
  end function add_tab_expansion
#endif

  !> set initialization file
  subroutine set_init_file(cfilename)
    character(len=*), intent(in) :: cfilename
    init_file_path = cfilename
  end subroutine set_init_file


  !***********************************************************************
  !       newcom-system (clone of incom)                                 *
  !       ------------                                                   *
  !       by  michael monkenbusch                                        *
  !       institut fuer festkoerperforschung des                         *
  !       forschungszentrum juelich, d-5170 juelich, po-box 1913, frg    *
  !       iff115 @ djukfa11                                              *
  !                                                                      *
  !***********************************************************************
  subroutine get_newcom( cmd , istatus)

    implicit none
    !
    character(*) :: cmd
    integer      :: istatus
    !      -----------------------
    !
    ! ---- command line decoder ----
    !
    !
    ! ---- separates commands parameternames & parameters ----
    !
    ! ---- commandline format may be :
    !      command <name parameter name parameter < ;command ....>>
    !
    ! ---- parameter separator is blank !!!
    ! ---- command separator is ;
    !
    ! ---- outputlevel
    !
    !
    !?       character*8   glabel,xlab,itypc
    character(cmd_len)  ::  glabel,xlab,itypc
    !character(cmd_long) ::  buf,bla132, numbuf
    character(cmd_long) ::  buf, numbuf
    !character(1)        ::  blank,csep,numtst(14),numeqv*14,inread*(cmd_long)
    character(1)        ::  blank,csep,inread*(cmd_long)
    character(1)        ::  crem
    !equivalence  (numeqv,numtst(1))
    logical :: name, lmakro = .false.
    logical :: file_exists

    character*1024 ma_fil
    integer       ilma, i, j, k, l, ii, ipmlen, isum, ioold
    integer       ier, ierr, ioldna, inew, iival

    integer :: init_run = 1, irc=0
    double precision   val

    integer, intent(in) ::  icmdus, ivnuse, irpuse
    integer             ::  iretus
    logical             ::  sflag
    logical,save        ::  aufcal=.true.

    save :: blank, csep, crem

    istatus = 0
    !
    ! ---- initialisations & conventions ----
    !
    if(init_run.eq.1) then

       ktop = 0 
       kanal(ktop) = input_unit
       !
       !data_path =  './'
       !save_path =  './'
       !makro_path = './'
       !datreat_path = './'
       csep = '&'  !              = --> command separator
       blank= ' '  !              = --> parameter separator
       crem = '!'  !              = --> comment !! >> mm (1/16) new comment feature
       !numeqv = '(.+-0123456789' !=> to reckognize the beginning of a number
       !bla132 = ' '
       !

!       prompt       = "-->"
!			call getenv('HOME', home)
        !write(prompt,'("drspine(",i0,".",i0,")-> ")')    _VERSION_MAJOR_ , _VERSION_MINOR_
#ifndef USE_LINENOISE
       do i=history_depth-1,0,-1
          history(i)=''
       enddo
#endif

       ! to load a initialisation macro if it exists
       !inquire(file=trim(datreat_path)//'makros/'//'initdatr',exist=fileda)
       if (len_trim(init_file_path)>0) then
          inquire(file=trim(init_file_path),exist=file_exists)
          if(file_exists) then      ! for a init file
             ioldc=1
             reslin=trim(init_file_path)//csep//trim(reslin)
             write(6,*) 'initialization file "'//trim(init_file_path)//'" found'
          else
             write(6,*) 'initialization file "'//trim(init_file_path)//'" not found'
          endif
       end if

#ifdef USE_LINENOISE
       ! initialize linenoise
       call fclinenoise_init(trim(history_file)//c_null_char)
#endif
       init_run = 0
!    endif
    !
    endif



8888 continue ! startingpoint of incom loop
    !      --------> reentry-point
    cused = .false.
    
    vname     = ' '
    rpar      = 0
    iparn     = 0
    inpar     = 0
    inapa     = 0
    vused     = .false.
    rused     = .false.
    item_used = .false.
   
    comand = ' '

    lstpar    = 0
    lstnam    = 0
    last_item = 0

    cmd_item     = " "
    item_nr_name = 0
    name_nr_item = 0
    item_nr_rpar = 0
    rpar_nr_item = 0
 

    !
    ! ---- error-response -----------------------------------------
    if(ierrr.ne.0 .or. isignal >=3 ) then
       if(ierrr   > 0)   write(6,*)' error return code:',ierrr
       if(isignal >=3) write(6,*)' Halt due to 3+ times ctrl-c',isignal
       ioldc = 0 !-----> if error forget rest of commandline (reslin)
       if(ktop.ne.0) then
          do k=ktop,1,-1
             close(kanal(k))
          enddo
          ktop = 0
          !            --------> go back to input from keyboard
       endif
    endif
    ierrr = 0
    !
    !
    ! ---- input a line or take nodecode rest of previous line ----
    if(ioldc.eq.0) then
       goto 1002
1001   continue

       if(iot.gt.-3) write(6,*)'end of data '
       if(ktop /= 0) close(kanal(ktop))
       !if(ktop.eq.0) then
       !   open(kanal(ktop))
       !else
       if(ktop /= 0) then
          ktop = ktop-1
          !                      --------> restore old reslin state
          ioldc = iolbuf
          reslin= rlbuf
          ipmls = 0
          if(ktop.eq.0) then 
            write(6,*)'## command input back to keybord ## '
          else
            write(6,'(a)')'# '
          endif 
          if(ioldc.ne.0) then
             inline = reslin

             goto 1011
          endif
       endif
1002   continue
       if(ktop.eq.0) then
#ifdef USE_LINENOISE
          istatus = fclinenoise_read(len(inread),inread, trim(prompt)//' '//c_null_char)
          raw_input = inread
          if (istatus<0) return
#else
          write(6,'(a)',advance='no')trim(prompt)//' ' 
          ipmls = 0
          read(kanal(ktop),"(a)",end=1001) inread
#endif
       else
          read(kanal(ktop),"(a)",end=1001) inread
          !replace tabs with spaces, bug #2295
          inread = replacetabs(inread)
       endif
       ! ---- last command from history ----
       ii=len_trim(inread(2:))
#ifndef USE_LINENOISE
       if(inread(1:1).eq.'_') then
          do i=0,history_depth
             if (index(history(i),inread(2:ii+1)).eq.1) then
                inread=history(i)
                exit
             endif
          enddo
       elseif(ktop.eq.0 .and. index('history',inread(:4)).ne.1) then
          do i=history_depth-1,0,-1
             history(i+1)=history(i)
          enddo
          history(0)=inread
       endif  ! go on with new command if found
#endif
       if(ktop.ne.0.and.iot.gt.-1)write(6,*)'! '//trim(inread)
       if(ktop.eq.0.and.iot.gt.-2)write(6,*)'> '//trim(inread)
       inline = inread

    else
       inline = trim(reslin)
       ioldc=0
    endif
1011 continue

!!mm>>
!write(*,*)"i1: ",trim(inline)
    inline = markstring(inline)
!write(*,*)"i2: ",trim(inline)
!!mm<<


    ioldc = 0
    ioldna =0
    !
    ! ---- look for a commnd separator ----
    if(SCAN(inline,csep).gt.0) then  ! falls csep vorkommt ist scan >0 , falls nicht =0
       ioldc = 1
       !         ---------------------> store second command on reslin
       reslin = inline(SCAN(inline,csep)+1:len(inline))
       inline = inline(1:SCAN(inline,csep)-1)//' '

    endif
    !
    ! ---- remove any leading blanks ----
    inline = ADJUSTL(inline)
    !

    !!>> mm(1/16): allow more modern comment style in macros
    j = SCAN(inline,crem)     ! finds comment separator
    if(j > 0) then
    !  write(6,*)crem," comment: ",trim(inline(j+1:))
       if( j == 1 ) goto 8888
       inline = inline(1:j-1)
    endif    
    !!<< mm(1/16): allow more modern comment style in macros
    
    !
    ! ---- separate the command ----   SCAN ('FORTRAN', 'R')=3
    j=SCAN(inline,' ')-1        ! findet erstes blank j ein davor
    if(j.gt.len(comand))then
       comand = inline(1:len(comand))
    else
       comand = inline(1:j)
    endif
    !
    cmd = comand

    ! --- ignore macro/makro command on the top-level
    if(ktop==0) then
       if(compar(comand,'makro ').or.compar(comand,'macro ')) goto 8888
    end if

    ! ---- dont analyse further if it is a pathdefinition ----
    if(compar(comand,'path ')) then
       if((inline(6:6).eq.'?') .or. (len_trim(inline).eq.len_trim('path') )) then
          write(6,*)'read data from : '//TRIM(data_path)
          write(6,*)'save data to   : '//TRIM(save_path)
          write(6,*)'load makro from: '//TRIM(makro_path)
          !write(6,*)'additional to path: '//TRIM(datreat_path)//'makros/'
          write(6,*)'change with datapath/savepath/macrpath  "PATH"'
          goto 8888
       endif

       goto 8888
    endif
    if(comand.eq.'datapath') then
       if((inline(10:10).eq.'?') .or. (len_trim(inline).eq.len_trim('datapath') )) then
          write(6,*)TRIM(data_path)
          goto 8888
       endif
       if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
          data_path=trim(inline(10:))//'/'
       else
          data_path=trim(inline(10:))
       endif
       if(iot.ge.0) write(6,*)'datapath: ',TRIM(data_path)
       goto 8888
    endif
    if(comand.eq.'savepath') then
       if((inline(10:10).eq.'?').or. (len_trim(inline).eq.len_trim('savepath') )) then
          write(6,*)TRIM(save_path)
          goto 8888
       endif
       if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
          save_path=trim(inline(10:))//'/'
       else
          save_path=trim(inline(10:))
       endif
       if(iot.ge.0) write(6,*)'savepath: ',TRIM(save_path)
       goto 8888
    endif
    if(comand.eq.'macrpath') then
       if((inline(10:10).eq.'?').or. (len_trim(inline).eq.len_trim('macrpath') )) then
          write(6,*)TRIM(makro_path)
          goto 8888
       endif
       if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
          makro_path=trim(inline(10:))//'/'
       else
          makro_path=trim(inline(10:))
       endif
       if(iot.ge.0) write(6,*)'macrpath: ',TRIM(makro_path)
       goto 8888
    endif

    if(comand.eq.'mxxpath') then
       if((inline(10:10).eq.'?').or. (len_trim(inline).eq.len_trim('mxxpath') )) then
          write(6,*)TRIM(mxx_path)
          goto 8888
       endif
       if (inline(len_trim(inline):len_trim(inline)).ne.'/') then
          mxx_path=trim(inline(10:))//'/'
       else
          mxx_path=trim(inline(10:))
       endif
       if(iot.ge.0) write(6,*)'mxxpath: ',TRIM(mxx_path)
       goto 8888
    endif

    ! ---- dont analyse further if it is a title ----
    if(comand.eq.'tit     '.or.comand.eq.'title   ') then
       if(comand.eq.'tit     ') title = trim(inline(5:))
       if(comand.eq.'title   ') title = trim(inline(7:))
       goto 8888
    endif
    ! ---- dont analyse further if it is a message ---
    if(comand.eq.'msg     ') then
       write(6,*)trim(inline)//'                 '
       goto 8888
    endif
    ! ---- dont analyse further if it is a blank line ---
    if(comand.eq.'        ') then
       !write(6,*)inline
       goto 8888
    endif
    ! ---- dont analyse further if it is a comment ---
    if(comand.eq.'c       ') then
       goto 8888
    endif
!    !!>> mm(1/16): allow more modern comment style in macros
!    if( inline(1:1) == '!' .or. inline(1:1) == '#') then
!       goto 8888
!    endif
!    !!<< mm(1/16): allow more modern comment style in macros

    
    ! ---- dont analyse further if it is a label ----
    if(comand(1:1).eq.':') then
       if(iot.gt.2)write(6,*)'label encountered ',comand
       goto 8888
    endif
    !
    ! ---- look for keywords & parameters ----
    !
    iitems = 0
    inames = 0
    ipars  = 0
    iargs  = 0
    iargvs = 0
    j = j+1    !erstes blank nach comand da j letzter char von comand
    do while (len_trim(inline(j+1:len(inline))).ne.0 )   ! schleife parameterersetzung
       i = j ! i = aktuelles erstes blank
       do while (inline(i+1:i+1).eq.blank)
          i=i+1       ! ! if there are more than single blank
       enddo
       l = i+1
       !   ---- look for the end of the currend item -----
       j=i+SCAN(inline(i+1:len(inline)),blank)      ! j = next blank
       if (j.eq.i) j=len(inline)+1       ! no further blanks till end j=ende+1 virtuell blank
       buf = inline(i+1:j-1)
       ! --- prepare the textlist of formal makro arguments --
       iargs = iargs + 1
       arglst(iargs) = buf(1:20)   !! >>>>>>>>> ??? To be checked MM !!!
       ! --- decode for makro argument values replace strings ! ----
       if(ktop.ne.0) then
          if(ipmlst(ktop).ne.0) then
             if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
             do k=1,ipmlst(ktop)
                if(iot.gt.3)write(6,*)buf
                call creplace(buf,pmlst(ktop,k,1)//' ',pmlst(ktop,k,2)//' ',' ', len(buf))
                if(iot.gt.3)write(6,*)buf
             enddo
             if(len_trim(buf).eq.0) cycle
          endif
       endif
       ! ----   discriminate between name & parameter ----

       iitems = iitems + 1
       cmd_item(iitems) = trim(buf)



       name = .true.
       ! starts with

!?       if((scan('(+-0123456789',buf(1:1)).gt.0) .or. &         !number starts with number or +-
!?            ('.'.eq.buf(1:1) .and. scan('(0123456789',buf(2:2)).gt.0)  .or. &  ! number omitted 0 like .534
!?            ((scan(buf(1:len(buf)),'()+-*/^').gt.0 .and.buf(1:2).ne.'./'.and.&
!?            scan('/^*.',buf(1:1)).eq.0.and. scan('/^*.',buf(len_trim(buf):len_trim(buf))).eq.0 ))  & ! this is a formula   !
!?            )   name = .false.       ! then it is a number or should be evaluated as formula
       if((scan('(+-.0123456789',buf(1:1)).gt.0)  &         !number starts with number or +-
                 .and.(buf(1:2).ne.'./')  & ! 
            )   name = .false.                        ! then it is a number or should be evaluated as formula
       if(.not.name) then ! .not.name = zahl oder formel
          call evaluate( buf//' ', val, ierr)
          !  write(6,*)'evaluate: ',buf(1:len_trim(buf))
          !  write(6,*)' gives ', val, ierr
          if(ierr.eq.0) then
             name = .false.
          else
             name = .true.                                    
!<             comand = 'f-error!'                            !! ??
!<            cmd    = comand  (inline(1:1).ne.'&') .and.
!<             write(6,*)'syntax error in: '//trim(buf)
!<             write(6,*)'       in line : '//trim(inline)
          endif
       endif
       !
       iargvs = iargvs+1
       if(name) then
          inames = inames + 1
          if(inames.gt.minc) goto 999
          vname(inames) = trim(buf(1:len(vname(1))))
!!mm>>
!write(*,*)"v1: ",inames,trim(vname(inames))
          vname(inames) = stripstring(vname(inames))
!write(*,*)"v1: ",inames,trim(vname(inames))
!!mm<<
          argvals(iargvs) = buf(1:1024)
          inapa(inames) = 0

          item_nr_name(iitems)  = inames         
          name_nr_item(inames) = iitems

       else
          ipars = ipars + 1
          if(ipars.gt.minc) goto 999
          rpar(ipars) = val
          iival = int(val+1d-11)
          if(dabs(val-iival).lt.3d-11) then
             write(argvals(iargvs),'(i132)') iival
          else
             write(argvals(iargvs),'(e132.12)') rpar(ipars)
          endif
          iparn(ipars) = inames
          if(inames.ne.ioldna) then
             inapa(inames) = ipars
             ioldna = inames
          endif

          item_nr_rpar(iitems)  = ipars
          rpar_nr_item(ipars)  = iitems
 
       endif
    enddo ! ende schleife parameterersetzung

    !
    ! ---- prepare inpar (no. of parameters / name) ----
    do i=1,inames
       isum = 0
       if(ipars.ne.0) then
          do j=1,ipars
             if(iparn(j).eq.i) isum = isum + 1
          enddo
       endif
       inpar(i) = isum
    enddo
    !
    !
    if(lmakro) goto 9888

    ! --------------- decoding of opaque commands -------------------------
    !

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! cms - kommando absetzen                             | cms <cmd >
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'cms     '.or.comand.eq.'sys     ') then
       !        same code as  bottom shell command            ---
       buf=ADJUSTL(inline(5:))
       if(ktop.ne.0) then
          if(ipmlst(ktop).ne.0) then
             if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
             do k=1,ipmlst(ktop)
                ipmlen=0
                do while (index(buf,trim(pmlst(ktop,k,1))).gt.ipmlen )
                   ipmlen=index(buf,trim(pmlst(ktop,k,1)))
                   buf=buf(1:ipmlen-1)//trim(pmlst(ktop,k,2))//trim(buf(ipmlen+len_trim(pmlst(ktop,k,1)):))
                   ipmlen=ipmlen+len_trim(pmlst(ktop,k,2))
                enddo
             enddo
          endif
       endif

!!       call system(trim(buf),ier)
       call execute_command(trim(buf), exitstat=ier)
       if(ier==0)  then
         cused = .true.
         rused = .true.
         vused = .true.
         item_used = .true.
       endif


       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! history anzeigen                             | history    hist
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef USE_LINENOISE
    if(index('history',comand(:4)).eq.1) then
       !                    ----
       if (index(vname(1),'clear').gt.0) then
          do i=history_depth-1,0,-1
             history(i)=''
          enddo
          goto 8888
       endif
       open(10,file=trim(makro_path)//'history')
       do i=history_depth-1,0,-1
          if (len_trim(history(i)).gt.0) then
             write(*,*) i,' ',trim(history(i))
             write(10,*) ,trim(history(i))
          endif
       enddo
       close(10)
       goto 8888
    endif
#endif
    !

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! cd - kommando absetzen                             | cd <dir >
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'cd     ') then
       !                    ---
       call change_directory(trim(inline(4:)), stat=ier, pwd=.true.)
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'..      ') then
       !                    ---
       call change_directory('..', stat=ier, pwd=.true.)
       goto 8888
    endif

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! output"level" setzen:                               | iout <iout>
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'iout    ') then
       !                    ----
       ioold= iot
       iot  = NINT(rpar(1))
       if(vname(1).eq.'old     ') iot = ioold
       write(6,*)'outputlevel is now ',iot
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! if                                                  |
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'if      ') then
       !                    ----------------> if-bedingung (numerisch)   !
       itypc = vnamef(1)
       

 !!>>>      ! reevaluate rpar to detect syntax errors
       if(.not. found("then    "))then
          call errsig(9009,'ERROR: if needs then as keyword $')
          goto 8888
       endif

      
       call evaluate(trim(cmd_item(1))//' ', val, ierr) 
       if(ierr .ne. 0) then
          call errsig(9010,'ERROR: if needs a number/valid expression as 1st value $')
          goto 8888
       endif      
       rpar(1) = val
       item_used(1) = .true.

       call evaluate(trim(cmd_item(3))//' ', val, ierr) 
       if(ierr .ne. 0) then
          call errsig(9011,'ERROR: if needs a number/valid expression as 2nd value $')
          goto 8888
       endif      
       rpar(2) = val
       item_used(3) = .true.

!!<<< 
       
       if(itypc.eq.'=       ') then
          if(rpar(1).eq.rpar(2)) then
             j = index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       if(itypc.eq.'>       ') then
          if(rpar(1).gt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       if(itypc.eq.'<       ') then
          if(rpar(1).lt.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       if(itypc.eq.'<=      ') then
          if(rpar(1).le.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       if(itypc.eq.'>=      ') then
          if(rpar(1).ge.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       if(itypc.eq.'<>      ') then
          if(rpar(1).ne.rpar(2)) then
             j=index(inline,'then')
             if(j.ne.0)then
                reslin = inline(j+4:)
                ioldc  = 1
             endif
          endif
          goto 8888
       endif
       write(6,*)'if operator ',trim(itypc), ' unknown'
       ierrs = 1
       call errsig(9012,'if construction contains invalid items! Check variables or operator! $ ')
       goto 8888
    endif
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! realisation goto                                    | goto
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'goto    ') then
       !                    ----------------> comment feature for makros !
       if(ktop.gt.0) then
          glabel = vnamef(1)
          if(glabel(1:1).ne.':') then
             write(6,*)glabel,' is not a valid label'
             write(6,*)'use >>> :',glabel,' <<< !'
             ierrs = 1
             goto 8888
          endif
          rewind(kanal(ktop))
88881     continue
          read(kanal(ktop),'(a8)',end=88883)xlab
          if(glabel.eq.xlab) goto 8888
          goto 88881
       else
          write(6,*)'goto in keyboardmode ignored'
       endif
88883  continue
       close(kanal(ktop))
       ktop = ktop-1
       write(6,*)'label:',glabel,' not found!'
       ierrs = 1
       goto 8888
    endif
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! return from makro                                   | return
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'return  ') then
       !                    ----------------> return from makro
       if(ktop.gt.0) then
          close(kanal(ktop))
          ktop = ktop - 1
       endif
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! exit all makros                                     | exit
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'exit    ') then
       !                    ----------------> exit   from makro
       if(ktop.gt.0) then
          do k=ktop,1,-1
             close(kanal(ktop))
          enddo
          ktop = 0
       else ! if we are at the command line level quit program
          istatus = -4
          return
       endif
       goto 8888
    endif
    !
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! evaluation info                                     | ??
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'??      ') then
       !                    ----------------> evaluation info !
       do i=1,ipars
          write(6,*) rpar(i)
       enddo
       goto 8888
    endif
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! evaluation and  print                               | print
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'print   ') then
       !                    ----------------> print evaluated numbers !

       if(inames.lt.1) then
          write(6,*)'print: needs filename'
          goto 8888
       endif
       open(99,file=vname(1),status='unknown',position='append')
       write(6,'(10(1x,e13.6))') (rpar(i),i=1,ipars)
       write(99,'(10(1x,e13.6))') (rpar(i),i=1,ipars)
       close(99)
       goto 8888
    endif

!>>new22
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! writinga tech table                                  | print
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'table   ') then
       !                    ----------------> print evaluated numbers !

       if(inames.lt.1) then
          write(6,*)'table: needs filename'
          goto 8888
       endif

       if(trim(vname(2))=="open") then
         open(99,file=trim(vname(1))//".tex",status='unknown',position='append')
         write(99,'(a)')"\begin{table}[h!]"
         write(99,'(a)',advance="no")"\begin{tabular}{|"
         j = iitems - 2
         do i=1,j 
           write(99,'(a)',advance="no")"l|"
         enddo
         write(99,'(a)')"}"
         write(99,'(a)')"\hline"
         do i=1,j-1
           write(99,'(a)')trim(cmd_item(i+2))//" &" 
         enddo
         write(99,'(a)')trim(cmd_item(i+2))
         write(99,'(a)')" \\ \hline" 
         close(99)
         goto 8888
       endif

       if(trim(vname(2))=="close") then
         open(99,file=trim(vname(1))//".tex",status='unknown',position='append')
             write(99,'(a)')"\hline"
             write(99,'(a)')"\end{tabular}"
             write(99,'(a)')"\caption{ \label{tab:"//trim(vname(1))//"}"
             write(99,'(a)')"Caption of table:"//trim(vname(1))//"}"
             write(99,'(a)')"\end{table}"
         close(99)
         goto 8888
       endif

         open(99,file=trim(vname(1))//".tex",status='unknown',position='append')
         j = iitems - 1
         if(j .ne. ipars) then
           write(*,*)"WARNING TABLE GENERATION number of items and entries mismatch!",i, ipars, trim(vname(1))
         endif
         do i=1,ipars-1 
           write(99,'(f10.3,a)', advance="no")rpar(i)," & "
         enddo
         write(99,'(f10.3,a)')rpar(ipars)," \\ "
         close(99)
       goto 8888
    endif
    !
!<<new22    !

    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! (user)vars zeigen                                   | vars?
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'vars?   ') then
       !                    ----------------> uservars display
       call shwudf
       goto 8888
    endif
     !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! set command line separator                           | csep
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'csep    ') then
       !                    ----------------> uservars display
       write(*,*)"csep actual/previous command line separator is: ",csep
       if(inames > 0) then
         csep = vname(1)(1:1)
         write(*,*)"csep command line separator (rest-line) is set to: ",csep
       endif
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! formelauswertung protokoll                          | say
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'ftrace  ') then
       !                    ----------> set deg-mode
       call sayit(.not.found('off     '))
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! degmode setzen                                      | setdeg
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'setdeg  ') then
       !                    ----------------> set deg-mode
       call setdeg
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! radmode  setzen                                     | setrad
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'setrad  ') then
       !                    ----------------> set deg-mode
       call setrad
       goto 8888
    endif
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! uservariable setzen                                 | set
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'set     ') then
       !                    ----------------> set uservar

!??       do i=1,inames
!??          if (inpar(i).eq.1) then
!??             call setudf(vname(i)//' ',rpar(inapa(i)),ier)
!??          else
!??             call errsig(1002,&
!??                  '"'//trim(vname(i))//'" cannot be evaluated $')
!??             return
!??          endif
!??       enddo


       if(iitems .ne. 2) then
          call errsig(9002,'ERROR: set needs ONE name and ONE number/expression as arguments $')
          goto 8888
       endif

       call evaluate(trim(cmd_item(2))//' ', val, ierr) 
       if(ierr .ne. 0) then
          call errsig(9003,'ERROR: set needs a number/valid expression as 2nd argument $')
          goto 8888
       endif      

       if( index("+-.()0123456789",cmd_item(1)(1:1)) > 0 ) then
          call errsig(9003,'ERROR: invalid variable name $')
          goto 8888
       endif 

       call setudf(trim(cmd_item(1))//' ',val,ier) 


       goto 8888
    endif
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! uservariable entefrnen                              | clr
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'clr     ') then
       !                    ----------------> set uservar
       do i=1,inames
          call clrudf(vname(i)//' ')
       enddo
       goto 8888
    endif
    !
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ausgang                                             | q
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'quit    '.or.comand.eq.'q       ') then
       !         old gr stuff
       !                    ----------------> exit the program
       !         if(ibild.gt.0) then
       !            write(6,*)'calling grende ...'
       !           !call grende
       !         endif
       vname(1)='cleanio '

       write(6,*)'exit by q, bye...'
       !write(6,*)"presently unlike in drspine the stop is performed within new_com"
       !write(6,*)"full stop now!"
       !stop
       istatus = -2
       if(comand.eq.'quit    ') istatus = -3
       return
    endif
    
!!!!! >>>>>> MM: TEST PRINT FOR DEVELOMENT ONLY !!!!!!
!  character(cmd_len)             ::  cmd_item(minc)    ! was 8
!  integer                        ::  item_nr_name(minc) !  i-th name associated with this item
!  integer                        ::  name_nr_item(minc) !  i-th item associated with this (j-th) name
!  integer                        ::  item_nr_rpar(minc) !  i-th rpar associated with this item
!  integer                        ::  rpar_nr_item(minc) !  i-th item associated with this (r-par)
    if(iot >  0) then 
       write(*,*)"NEW_COM -> PARSER ANALYSIS"
       write(*,*)"NEW_COM inline : >",trim(inline),"<"
       write(*,*)"NEW_COM buf    : >",trim(buf),"<"
       write(*,*)
       write(*,*)" #    item               ptnam ptrpar    name                 rpar "
       write(*,*)" ------------------------------------------------------------------"
       do i=1,iitems
         write(*,'(i3,": ",a20,2i6,4x,a20,4x,e14.7)') i, cmd_item(i)(1:20), item_nr_name(i), item_nr_rpar(i), &
                                              vname(max(1,item_nr_name(i)))(1:20),rpar(max(1,item_nr_rpar(i)))
       enddo
       write(*,*)" ------------------------------------------------------------------"
   
       write(*,*)
       write(*,*)" #    name               ptitem inpar   item                 "
       write(*,*)" ------------------------------------------------------------------"
       do i=1,inames
         write(*,'(i3,": ",a20,2i6,4x,a20)') i, vname(i)(1:20), name_nr_item(i),inpar(i), &
                                             cmd_item(max(1,name_nr_item(i)))(1:20)
       enddo
       write(*,*)" ------------------------------------------------------------------"
   
       write(*,*)
       write(*,*)" #    rpar               ptitem iparn inapa   item                 "
       write(*,*)" ------------------------------------------------------------------"
       do i=1,ipars
         write(*,'(i3,": ",e20.7,3i6,4x,a20)') i, rpar(i), rpar_nr_item(i), iparn(i), inapa(i), &
                                             cmd_item(max(1,rpar_nr_item(i)))(1:20)
       enddo
       write(*,*)" ------------------------------------------------------------------"
     endif   
   



    do i=1,inames
       isum = 0
       if(ipars.ne.0) then
          do j=1,ipars
             if(iparn(j).eq.i) isum = isum + 1
          enddo
       endif
       inpar(i) = isum
    enddo
 
!!!!! <<<<<<


!
    return
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! ---- stack overfow error trap ----
999 continue
    write(6,*)' incom(',comand,') :  stack overflow   reenter'
    ioldc = 0
    ierrr = 2
    goto 8888

    !
    !
    !
    ! --- makro - aktivator - entry ----------------------------------------
    !
    entry makro (cmd)
    !      ===========
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! --- try to interpret non-known commands as makro input files ---
    ! ------> any makro input file must have the keyword  <makro>
    !         at the first position
    if(ktop.ge.mdepth) then
       write(6,*)' too many stacked makros : ',comand
       ierrr = 3
       return
    endif
    ktop = ktop + 1

!!! Modification for two makro sources: global, individual !!!
!!! try first the individual makro as we had it already    !!!

9998 format(a80)

    open(newunit=kanal(ktop),file=comand,status='old',form='formatted',err=99991)
    read(kanal(ktop),9998,end=99991) inline
!!  if(inline(1:6).eq.'makro '.or. inline(1:6).eq.'macro ') goto 19999
!! mm>> 10.12.18:
    if(inline(1:6).eq.'makro '.or. inline(1:6).eq.'macro ') then
      if( ktop == 1 ) then
         top_macro_file = comand
         top_macro_call = raw_input
      endif
      goto 19999
    endif
!! mm<<
    !! ok this was no makro !!!
    close(kanal(ktop))

    !! ok this was no makro !!!
    close(kanal(ktop))
99991 continue

    ! maybe on user-std-makropath
    ma_fil = trim(makro_path)//comand
    ilma = LEN_TRIM(ma_fil)
    open(newunit=kanal(ktop),file=ma_fil(1:ilma),status='old',form='formatted',err=99992)
    read(kanal(ktop),9998,end=99992) inline
    if(inline(1:6).eq.'makro '.or. inline(1:6).eq.'macro ') goto 19999
    !! ok this was no makro !!!
    close(kanal(ktop))
99992 continue

    goto  9999
    !
    ! -------> try whether in the global makro section file is a makro file
    !open(newunit=kanal(ktop),file=trim(datreat_path)//'makros/'//comand,status='old',form='formatted',err=9999)
    !read(kanal(ktop),9998,end=9999) inline
    !if(inline(1:5).ne.'makro') goto 9999
    ! -------> try whether file is a makro file
    !
19999 continue
    ! --- ok this is a makro file ! save first the parametervalues ---
    cused = .true.
    do i=1,minc
       pmlst(ktop,i,2)=' '
       rused(i) = .true.
       vused(i) = .true.
    enddo
    do i=1,iargvs
       buf = ADJUSTL(argvals(i)) ! (PAZ) silence truncation warning
       pmlst(ktop,i,2) = buf(:cmd_len)
       if(iot.gt.4) write(6,*)'argvals : '//'##'//trim(ADJUSTL(argvals(i)))//'##'//pmlst(ktop,i,2)//'##'
    enddo
    ipmlst(ktop) = iargvs
    do i=iargvs+1,minc
       pmlst(ktop,i,2) = ' '
    enddo
    ! --- save reslin form the current keyboard command ---
    iolbuf = ioldc
    rlbuf  = reslin
    ! --- decode the first makro line to get the parameternames
    ioldc = 1
    reslin = inline
    lmakro = .true.
    if(iot.gt.3)write(6,*)'this is a makro ...'
    goto 8888
9888 continue
    lmakro = .false.
    cused  = .true.
    if(iot.gt.3)write(6,*)'go on with it ...'
    if(iargs.ne.0) then
       do i=1,iargs
          if(iot.gt.3) write(6,*)'arglst:',arglst(i)
          pmlst(ktop,i,1) = arglst(i)
          rused(i) = .true.
          vused(i) = .true.
       enddo
    endif
    ipmlst(ktop) = iargs
    !        ipmls = iargs
    !
    return
    !      ---------> read now commandlines from makro file

9999 continue
    ktop = ktop - 1 ! back to old ktop macrolevel if it was no macro
    !###################################################
    !        if command not known till now we try to sent it to the shell
    buf=inline
    if(ktop.ne.0) then
       if(ipmlst(ktop).ne.0) then
          if(iot.gt.3)write(6,*)'replacing.. ipmlst(',ktop,')=',ipmlst(ktop)
          do k=1,ipmlst(ktop)
             ipmlen=0
             do while (index(buf,trim(pmlst(ktop,k,1))).gt.ipmlen )
                ipmlen=index(buf,trim(pmlst(ktop,k,1)))
                buf=buf(1:ipmlen-1)//trim(pmlst(ktop,k,2))//trim(buf(ipmlen+len_trim(pmlst(ktop,k,1)):))
                ipmlen=ipmlen+len_trim(pmlst(ktop,k,2))
             enddo
          enddo
       endif
    endif
!!    call  system(trim(buf),ier)
    call execute_command(trim(buf), exitstat=ier) ! F2008
    if(ier == 0 ) then
      cused = .true.
      rused = .true.
      vused = .true.
    endif
    If (irc .ne. -1) return

    !        back from (s)hell
    !###################################################
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(comand.eq.'f-error!') then
       !                    ----------------> fehler-reaktion !
       write(6,*)'!! commandline in error --> ignored !!'
       ierrr = 111
       goto 8888
    endif
    !
    do k=ktop,1,-1
       close(kanal(k))
    enddo
    ktop = 0
    write(6,"(' program name ''',a8,''' not known')")comand
    ierrr = 4

    return



       entry close_all_macros()
!!     ========================
        do k=ktop,1,-1
         close(kanal(k))
       enddo
       ktop = 0
       write(6,*)'all macros closed '
       return



    entry unused( icmdus, ivnuse, irpuse, iretus)
    !      ---------------------------------------------
    !      check whether all command line items have been used
    !      icmdus, ivnuse, irpuse  are message levels
    !      0 : ignore
    !      1 : warning
    !      2 : error ---> sets iretus > 0
    !---------------------------------------------------------------------
    iretus = 0

    !  -- ignore the first call ! maybe this should be handled by the caller ?!
    if(aufcal) then
       aufcal = .false.
       return
    endif

!! for backward compatibility we distribute transfer item usage 
    do i=1,iitems
        if(item_nr_name(i)   > 0 ) vused(item_nr_name(i))   = vused(item_nr_name(i)) .or. item_used(i) 
        if(item_nr_rpar(i)   > 0 ) rused(item_nr_rpar(i))   = rused(item_nr_rpar(i)) .or. item_used(i) 

        if(item_nr_name(i)   > 0 ) item_used(i)             = vused(item_nr_name(i)) .or. item_used(i) 
        if(item_nr_rpar(i)   > 0 ) item_used(i)             = rused(item_nr_rpar(i)) .or. item_used(i) 
    enddo



    if(iot.gt.0) write(6,*)'command: ',trim(comand)
    if(.not.cused) then
!mmnc       cused = .true.
       if(icmdus.gt.0) then
          call errsig(icmdus-1,'command not used:'// &
               trim(comand)//'$')
          if(icmdus.gt.1) iretus = iretus + 10000
       endif
    endif

    do i=1,iitems
       if(iot.gt.0) &
            write(6,*)'item: ',i,cmd_item(i)(1:lclen(vname(i),80)+1)
       if(.not.item_used(i)) then
!mmnc          item_used(i) = .true.
          if(ivnuse.gt.0 .or. irpuse.gt.0) then
             call errsig(max(ivnuse-1,irpuse-1),'item    not used:'// &
                  cmd_item(i)(1:lclen(vname(i),80)+1)//'$')
             if(max(ivnuse,irpuse).gt.1) iretus = iretus + 200*i
          endif
       endif
    enddo


    do i=1,inames
       if(iot.gt.0) &
            write(6,*)'vname: ',i,vname(i)(1:lclen(vname(i),80)+1)
       if(.not.vused(i)) then
!mmnc          vused(i) = .true.
          if(ivnuse.gt.0) then
             call errsig(ivnuse-1,'name    not used:'// &
                  vname(i)(1:lclen(vname(i),80)+1)//'$')
             if(ivnuse.gt.1) iretus = iretus + 100*i
          endif
       endif
    enddo

    do i=1,ipars
       if(iot.gt.0) Write(6,*)'rpar :',i,rpar(i)
       if(.not.rused(i)) then
!mmnc          rused(i) = .true.
          if(irpuse.gt.0) then
             write(numbuf,'(e13.6)') rpar(i)
             call errsig(irpuse-1,'number  not used:'// &
                  numbuf(1:13)//'$')
             if(irpuse.gt.1) iretus = iretus + i
          endif
       endif
    enddo

    return


    entry setena( sflag )
    !      ---------------------
    !      enable/disable creation of uservariables
    !-----------------------------------------------------------------------
    seten = sflag
    return


  END subroutine get_newcom


  !cc+++ noch optimierungsbeduerftig !!! ++++++++++++++++++<<<!!!!!!
  subroutine creplace(a,b,c,s, mclen)
    !---------------------------------------------------
    ! ersetzten von string b durch c in string a
    ! strings werden durch den character s begrenzt !
    !----------------------------------------------------
    implicit none
    integer, intent(in) :: mclen
    !integer, parameter :: mclen=1024
    character(len=*) a,b,c
    character(len=mclen) d
    character*1 s

    integer i, j, l, jj, ld, istart, iend, len
    !
    ! ---  look for occurence of substring b in a ---
    i = 1
    d = a
    len = 0
    do l=1,mclen
       if(a(l:l).eq.s) then
          len = l-1
          goto 1
       endif
    enddo
1   continue
    a = d
    if(a(i:i).eq.b(1:1)) then
       istart = i
       do j=2,mclen
          if(b(j:j).eq.s  ) goto 10
          l=i-1+j
          if(b(j:j).ne.a(l:l)) then
             istart = 0
             goto 10
          endif
       enddo
10     continue
       if(istart.ne.0) then
          iend = l
          do jj=1,mclen
             ld= istart-1+jj
             if(c(jj:jj).eq.s  ) goto 20
             d(ld:ld)= c(jj:jj)
          enddo
          i = ld
20        continue
          do jj=iend+1,mclen
             if(l.gt.mclen) goto 30
             d(ld:ld)= a(jj:jj)
             if(a(jj:jj).eq.s  ) goto 30
             ld = ld+1
          enddo
30        continue
       endif
    endif
    i = i+1
    if(i.le.mclen) then
       if( a(i:i).ne.s  ) goto 1
    endif

    do l=1,len
       if(d(l:l).eq.s)then
          do i=l,len
             d(i:i)=' '
          enddo
       endif
    enddo


    a = d

    return
  END subroutine creplace

  subroutine re_scan(ctext,r,ierr)
    ! ----------------------------------------------------------------------
    !  scan - scannen einer real-zahl (version fuer ibm pc prof. fortran)
    !  autor: g. egerer, datum: 26. 9.85, letzte aenderung: 31. 7.87
    ! ----------------------------------------------------------------------
    !                           * eingabeparameter:
    !                           * ctext  - textzeile (64 zeichen)
    !                           * ausgabeparameter:
    !                           * r      - real-zahl
    !                           * ierr   - fehlercode
    implicit none
    double precision   r
    integer ierr
    integer, parameter       :: ctxtle = 64
    !      character(1), intent(in) :: ctext(*)
    character(*), intent(in) :: ctext
    !character(len=ctxtle)    :: ctxbuf
    character*10       form


    !                           * benamte programmkonstanten:
    integer, parameter :: cnil = -32766

    !                           * programmvariablen:
    integer   i,j
    integer   icur, iscntb(9,0:7), istate, isymbl

    !integer  :: itest


    !                           * initialisierungen:
    ! tabelle zur syntaxdefinition:
    !
    ! symbol   s0   " "   "-"   "+"   "."   "e"   "d"   ziff
    ! ------+--------------------------------------------------
    ! zu-   |
    ! stand |
    ! --> 1 |         1     2     2     3                  8
    !     2 |                           3                  8
    !     3 |                                              9
    !     4 |               5     5                        6
    !     5 |                                              6
    ! end-  |
    ! zust. |
    !     6 |         7                                    7
    !     7 |         7
    !     8 |         7                 9     4     4      8
    !     9 |         7                       4     4      9
    !
    ! s0 (symbolklasse 0): enthaelt alle symbole, die nicht gesondert aufge-
    !                      fuehrt sind
    ! ziff               : enthaelt alle ziffern von "0" - "9"
    !
    data ((iscntb(i,j), j = 0,7), i = 1,9)                            &
         &   / cnil,    1,    2,    2,    3, cnil, cnil,    8,              &
         &     cnil, cnil, cnil, cnil,    3, cnil, cnil,    8,              &
         &     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    9,              &
         &     cnil, cnil,    5,    5, cnil, cnil, cnil,    6,              &
         &     cnil, cnil, cnil, cnil, cnil, cnil, cnil,    6,              &
         &     cnil,    7, cnil, cnil, cnil, cnil, cnil,    7,              &
         &     cnil,    7, cnil, cnil, cnil, cnil, cnil, cnil,              &
         &     cnil,    7, cnil, cnil,    9,    4,    4,    8,              &
         &     cnil,    7, cnil, cnil, cnil,    4,    4,    9 /
    !                           * ------------------------------------------

    icur = 1
    istate = 1
1100 continue
    if ( icur .gt. ctxtle ) then
       if ( istate .ge. 6 ) then
          !                           * istate ist ein endzustand
          ierr = 0
          form = '(BN,Eww.0)'
          write (form(6:7),'(i2)') ctxtle
          !            write (ctxbuf,'(64a)') ctext(1:64)   ! convert character vector to character string !!
          !            read(ctxbuf,form) r
          read(ctext,form) r
       else
          !                           * istate ist kein endzustand
          ierr = icur
       endif
    else
       !                           * ermittle index des aktuellen symbols
       !?        isymbl = index(' -+.EDed0123456789',ctext(icur))
       isymbl = index(' -+.EDed0123456789',ctext(icur:icur))
       if (isymbl .ge. 7 ) then
          if (isymbl .ge. 9 ) then
             isymbl = 7
          else
             isymbl = isymbl - 2
          endif
       endif
       istate = iscntb(istate,isymbl)
       if ( istate .ge. 0 ) then
          icur = icur+1
          go to 1100
       endif
       ierr = icur
    endif
  END subroutine re_scan


  subroutine parse( zeile )
    !     =========================
    ! ----------------------------------------------------------------------
    ! zerlegung der zeile nach incom manier
    ! die daten der letzten eingabedekodierung werden dabei ueberschrieben
    ! eine ggf. vorhandene restzeile bleibt erhalten
    ! restzeilen (; ...) beim parsen werden vollstaendig ignoriert !
    ! ----------------------------------------------------------------------

    implicit none

    character(cmd_long) :: zeile, reslbf
    character(cmd_len)  ::  cmd
    integer ioldbf, istatus
    !-----------------------------------------------------------------------
    !
    reslbf = reslin
    ioldbf = ioldc
    reslin = zeile
    !       call lowcase(reslin,)
    ioldc  = 1
    call get_newcom (cmd, istatus)

    reslin = reslbf
    ioldc  = ioldbf

    return
  END subroutine parse
  !
  !
  character(cmd_long) function X_datpath()
    !     ========================================
    implicit none

    X_datpath(1:cmd_long) = data_path

    return
  end function X_datpath

  character(cmd_long) function savepath()
    !     =======================================
    implicit none

    savepath = save_path

    return
  end function savepath


  integer function iout()
    !     =======================
    implicit none
    iout = iot
    return
  END function iout


  subroutine errsig(ierr,say)
    !     ===========================
    implicit none
    integer, intent(in)      ::ierr
    integer                  ::lsay
    character(*), intent(in) :: say
    ! ----------------------------------------------------------------------
    !  error signalisierung
    ! ----------------------------------------------------------------------
    ierrr = ierr
    ierrs = ierrs+1                 ! check which one we need 
    lsay = laenge(say,cmd_len,'$')
    write(6,*)'error:',ierr,' ',say(1:lsay)
    if(iot.gt.0) then
        write(*,*)
        write(*,*)"raw_input        : ",trim(raw_input)
        write(*,*)"inline           : ",trim(inline)
        write(*,*)"top_macro_file   : ",trim(top_macro_file)
        write(*,*)"top_macro_call   : ",trim(top_macro_call)
    end if
    return
  END subroutine errsig




  double precision   function getval(pname,defval,inew0)
    !-----------------------------------------------------------------------
    !  wert extraktion aus dem incom parameterstack
    !-----------------------------------------------------------------------
    implicit none

    character(len=*), intent(in) :: pname
    double precision             :: defval
    integer, intent(out),optional:: inew0
    !
    integer :: i, inew, ev_err
    double precision :: ev_val
    !-----------------------------------------------------------------------
    inew      = 0
    last_item = 0
    if(present(inew0)) inew0 = inew
    getval    = defval
    if(inames.le.0) return


   search: do i=1,iitems
      if(compar(cmd_item(i),pname//' ')) then
          item_used(i)                                        = .true.
          call evaluate(trim(cmd_item(i+1))//' ',ev_val,ev_err)
          if(ev_err .ne. 0) then
            call errsig(9991,"ERROR: Parameter: "//trim(pname)//" requires a numeric arg $ ")
            ierrs   = ierrs+9991 !??
            inew   = -1
            exit search
          endif
          item_used(i+1) = .true.
          last_item      = i+1
          getval         = ev_val
          inew = i
      endif
    end do search
  

    if(present(inew0)) inew0 = inew

    return

  END function getval




  integer function intval(pname,idef,inew0)
    !-----------------------------------------------------------------------
    !  wert extraktion aus dem incom parameterstack  integer
    !-----------------------------------------------------------------------
    implicit none
    character(len=*) :: pname
    integer          :: idef
    integer,optional :: inew0
    !
    integer :: i, inew, ev_err
    double precision :: ev_val
    !-----------------------------------------------------------------------
    inew      = 0
    last_item = 0
    if(present(inew0)) inew0 = inew
    intval    = idef
    if(inames.le.0) return
    !

   search: do i=1,iitems
      if(compar(cmd_item(i),pname//' ')) then
          item_used(i)                                        = .true.
          call evaluate(trim(cmd_item(i+1))//' ',ev_val,ev_err)
          if(ev_err .ne. 0) then
            call errsig(9991,"ERROR: Parameter: "//trim(pname)//" requires a numeric arg $ ")
            ierrs   = ierrs+9991 !??
            inew   = -1
            exit search
          endif
          item_used(i+1) = .true.
          last_item      = i+1
          intval         = NINT(ev_val)
          inew = i
      endif
    end do search
  

    if(present(inew0)) inew0 = inew
    return


  END function intval



  double precision   function valnxt(defval,inew0)     
    !-----------------------------------------------------------------------
    !  wert extraktion aus dem incom parameterstack
    !-----------------------------------------------------------------------
    implicit none

    double precision, intent(in)              ::  defval
    integer,          intent(out), optional   ::  inew0

    double precision ::  inew, ev_val
    integer          ::  ev_err

    inew      = 0
    valnxt    = defval
    if(present(inew0)) inew0 = inew


!??    if(ipars.le.lstpar) return
!??
!??    valnxt = rpar(lstpar+1)
!??    rused(lstpar+1) = .true.
!??    inew   = lstpar + 1
!??    lstpar = inew

    if(iitems <= last_item) return
    last_item = last_item+1
    call evaluate(trim(cmd_item(last_item))//' ',ev_val,ev_err)
    if(ev_err .ne. 0) then
       call errsig(9995,"ERROR: Parameter: exected a numeric arg but could not evaluate $ ")
       ierrs   = ierrs+9995 !??
       inew   = -1
    else
       valnxt               = ev_val
       item_used(last_item) = .true.
       inew                 = last_item   !! check usage of valnxt whether this is ok
    endif

    if(present(inew0)) inew0 = inew
     

  END function valnxt

  integer function ifound(pname)
    !-----------------------------------------------------------------------
    !  logische optionserkennung  ifound wird=fundstelle sonst 0
    !-----------------------------------------------------------------------
    implicit none

    character(*) pname

    integer i
    !-----------------------------------------------------------------------

    ifound    = 0
    if(inames.le.0) return

    do i=1,inames
       if(compar(vname(i),pname//' ')) then
          vused(i) = .true.
          ifound=i
          lstnam=i
          last_item = name_nr_item(i)  
          goto 20
       endif
    end do

20  continue

    return

  END function ifound


  logical function found(pname)
    !-----------------------------------------------------------------------
    !  logische optionserkennung
    !-----------------------------------------------------------------------
    implicit none
    character(*) pname



    found = ifound(pname) .ne. 0
    return

  END function found



  logical function folgt(popt,pname)
    !-----------------------------------------------------------------------
    !  logische optionserkennung, wahr wenn popt pname folgt
    !-----------------------------------------------------------------------
    implicit none

    character(*) popt, pname

    integer j
    !-----------------------------------------------------------------------

    folgt = .false.
    j     = ifound(pname)
    if(j.gt.0 .and. j.lt.inames) then
       folgt = compar(popt//' ', vname(j+1))
    endif
    if(folgt) then
       lstnam = j+1
       last_item = name_nr_item(j+1)  
       vused(j+1) = .true.
    endif

    return

  END function folgt


  integer function intvno(ipnum,idef,inew)
    !-----------------------------------------------------------------------
    !  wert extraktion aus dem incom parameterstack  integer nach addr
    !-----------------------------------------------------------------------
    implicit none

    integer ipnum,idef,inew

    inew      = 0
    intvno    = idef
    if(ipars.lt.ipnum) return

    intvno = NINT(rpar(ipnum))
    rused(ipnum) = .true.
    inew      = ipnum
    lstpar    = ipnum
    last_item = rpar_nr_item(ipnum) 
    return

  END function intvno

  double precision   function getvno(ipnum,adef,inew)
    !-----------------------------------------------------------------------
    !  wert extraktion aus dem incom parameterstack  real    nach addr
    !-----------------------------------------------------------------------
    implicit none

    double precision   adef
    integer ipnum, inew

    inew      = 0
    getvno    = adef
    if(ipars.lt.ipnum) return

    getvno = rpar(ipnum)
    inew   = ipnum
    lstpar = ipnum
    last_item = rpar_nr_item(ipnum) 

    return

  END function getvno




  subroutine getnva(pname ,pardef ,parout, np, inew)
    !-----------------------------------------------------------------------
    !  extraktion meherere einem namen zugeordneter werte
    !-----------------------------------------------------------------------
    implicit none

    !-----------------------------------------------------------------------
    character(*) pname
    double precision   pardef, parout
    integer np, inew
    dimension pardef(np), parout(np)

    character(cmd_len)   :: napp
    character(2*cmd_len) :: pnamc
    integer :: i=0, j=0, ier, napr
    !-----------------------------------------------------------------------

    inew      = 0
    do i=1,np
       parout(i) = pardef(i)
    end do
    if(inames.le.0) goto 20

    do i=1,inames
       if(compar(vname(i),pname//' ')) then
          vused(i)            = .true.
          napr = inpar(i)
          if(napr.gt.np) napr = np
! if(napr<=0)  call errsig(9991,"Parameter: "//trim(pname)//" requires a numeric (numbre or expression) arg.$ ")
          if(napr.gt.0) then
             do j=1,min(inpar(i),np)
                parout(j) = rpar(inapa(i)-1+j)
                rused(inapa(i)-1+j) = .true.
             end do
             inew = napr
             lstpar = inapa(i)+inpar(i)-1
             lstnam = i
          endif
          goto 20
       endif
    end do
20  continue

    do i=1,np
       write(napp,'(i3,5x)')100+i
       napp(1:1)='.'
       call cappend(pname//' ',napp,pnamc)
       call setudf(pnamc//' ',parout(i),ier)
    end do


    !
    return
  END subroutine getnva


!!mm>>
   subroutine push_cmd_line(cmdl, ierr)
     implicit none
     character(len=*), intent(in)   :: cmdl
     integer, intent(out), optional :: ierr

     if(present(ierr)) ierr = 0
     if(ioldc > 0) then
       write(*,*)"ERROR: previous residual command line part is lost! ", trim(reslin) 
       if(present(ierr)) ierr = 1
     endif
     
     reslin  = trim(cmdl)
     ioldc   = 1
     
   end subroutine push_cmd_line
!!mm<<


  subroutine rotavc( x, xyorig, rotvec )
    !     ======================================
    !     unterwirft den koordinatenvector x einer drehung um xyorigin
    !     um den winkel rotang
    !-----------------------------------------------------------------------
    implicit none
    double precision   :: x(3), y(3), xyorig(3), rotvec(3), r(3,3)
    double precision   :: degtor=0.017453292d0
    double precision   :: angdeg

    y(1) = x(1) - xyorig(1)
    y(2) = x(2) - xyorig(2)
    y(3) = x(3) - xyorig(3)
    call rotax(rotvec, degtor, r, angdeg )
    call matvec( r, y, y)
    x(1) = y(1) + xyorig(1)
    x(2) = y(2) + xyorig(2)
    x(3) = y(3) + xyorig(3)

    return
  END subroutine rotavc


  subroutine unrot( x, xyorig, rotvec )
    !     =====================================
    !     unterwirft den koordinatenvector x einer drehung um xyorigin
    !     um die drehung         rotang**-1
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: x(3), y(3), xyorig(3), rotvec(3), r(3,3)
    double precision   :: degtor=0.017453292d0
    double precision   :: angdeg

    y(1) = x(1) - xyorig(1)
    y(2) = x(2) - xyorig(2)
    y(3) = x(3) - xyorig(3)

    call rotax(rotvec, degtor, r, angdeg )
    call mattxp( r, r)
    call matvec( r, y, y)

    x(1) = y(1) + xyorig(1)
    x(2) = y(2) + xyorig(2)
    x(3) = y(3) + xyorig(3)
    !
    return
  END subroutine unrot


  subroutine rotax( ax, scale, rmat, angdeg )
    !      ===========================================
    !                       -->  ---->  <---   <---
    !-----------------------------------------------------------------------
    ! drehmatrix rmat wird als drehung um die achse ax berechnet
    ! die laenge wird zur bestimmung des drehwinkels mit scale multipliziert
    !-----------------------------------------------------------------------

    implicit none
    double precision   :: ax(3), rmat(3,3), en(3)
    integer :: i, j
    double precision   :: angdeg, scale, ang, co, si, s
    !-----------------------------------------------------------------------

    do i=1,3
       do j=1,3
          rmat(j,i) = 0.0d0
       end do
    end do

    s    = dsqrt( ax(1)**2 + ax(2)**2 + ax(3)**2 )
    if( s.gt.0.d0) then
       ang  = s * scale
       angdeg = ang * 57.29577951d0
       co   = dcos(ang)
       si   = dsin(ang)
       s    = 1.d0 / s
       en(1)= ax(1) * s
       en(2)= ax(2) * s
       en(3)= ax(3) * s
       do i=1,3
          do  j=1,3
             rmat(i,j) = en(i)*en(j) * ( 1.d0 - co )
          end do
       end do
       rmat(1,1) = rmat(1,1) + co
       rmat(2,2) = rmat(2,2) + co
       rmat(3,3) = rmat(3,3) + co

       rmat(1,2) = rmat(1,2) - si * en(3)
       rmat(2,1) = rmat(2,1) + si * en(3)

       rmat(3,1) = rmat(3,1) - si * en(2)
       rmat(1,3) = rmat(1,3) + si * en(2)

       rmat(2,3) = rmat(2,3) - si * en(1)
       rmat(3,2) = rmat(3,2) + si * en(1)
    else
       rmat(1,1) = 1.d0
       rmat(2,2) = 1.d0
       rmat(3,3) = 1.d0
       ang       = 0.d0
       angdeg    = 0.d0
    endif

    return
  END subroutine rotax


  subroutine matmux( a, b, c )
    !      ============================
    !-----------------------------------------------------------------------
    !      3x3 matrixmultiplikation
    !      a * b ---> c
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: a(3,3), b(3,3), c(3,3)
    double precision   :: h(3,3)
    double precision   :: s
    integer :: i,j,l

    do i=1,3
       do j=1,3
          s = 0.d0
          do l=1,3
             s = s + a(i,l)*b(l,j)
          end do
          h(i,j) = s
       end do
    end do

    do  i=1,3
       do  j=1,3
          c(j,i) = h(j,i)
       end do
    end do

    return
  END subroutine matmux


  subroutine matvec( a, x, y )
    !      ============================
    !-----------------------------------------------------------------------
    !      3x3 matrix - vector - multiplikation
    !      a * x ---> y
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: a(3,3), x(3), y(3)
    double precision   :: h(3)
    double precision   :: s
    integer :: i,l

    do i=1,3
       s = 0.d0
       do l=1,3
          s = s + a(i,l)*x(l)
       end do

       h(i) = s
    end do

    do i=1,3
       y(i) = h(i)
    end do

    return
  END subroutine matvec


  subroutine mattxp( a, at)
    !      =========================
    !-----------------------------------------------------------------------
    !      3x3 matrix - transposition
    !      a(t)  ---> at
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: a(3,3), at(3,3)
    double precision   :: xd, xu

    at(1,1) = a(1,1)
    at(2,2) = a(2,2)
    at(3,3) = a(3,3)

    xu      = a(1,2)
    xd      = a(2,1)
    at(1,2) = xd
    at(2,1) = xu

    xu      = a(1,3)
    xd      = a(3,1)
    at(1,3) = xd
    at(3,1) = xu

    xu      = a(2,3)
    xd      = a(3,2)
    at(2,3) = xd
    at(3,2) = xu

    return
  END subroutine mattxp




  subroutine matone( a )
    !      ======================
    !-----------------------------------------------------------------------
    !      3x3 matrix ---> einheitsmatrix
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: a(3,3)
    integer :: i, j

    do i=1,3
       do j=1,3
          a(j,i) = 0.d0
       end do
       a(i,i) = 1.d0
    end do

    return
  END subroutine matone


  subroutine axrot( rmat, ax)
    !      ===========================
    !-----------------------------------------------------------------------
    !      rotationsachsendarstellung aus rmat
    !      rmat ==~==> ax
    !-----------------------------------------------------------------------
    implicit none

    double precision   :: rmat(3,3), ax(3), rtest(3,3)
    double precision   :: d1(3), d2(3), d3(3) , x1(3), x2(3), x3(3)
    double precision   :: tx(3), ty(3)
    double precision   :: angdeg, ang, err, co, s1, s2, s3, sn, smax
    integer :: imax = 1

    d1(1) = rmat(1,1) - 1.d0
    d1(2) = rmat(2,1)
    d1(3) = rmat(3,1)

    d2(1) = rmat(1,2)
    d2(2) = rmat(2,2) - 1.d0
    d2(3) = rmat(3,2)

    d3(1) = rmat(1,3)
    d3(2) = rmat(2,3)
    d3(3) = rmat(3,3) - 1.d0

    x1(1) = d1(2)*d2(3) - d1(3)*d2(2)
    x1(2) = d1(3)*d2(1) - d1(1)*d2(3)
    x1(3) = d1(1)*d2(2) - d1(2)*d2(1)

    s1    =  x1(1)**2 + x1(2)**2 + x1(3)**2
    imax =   1
    smax =   s1

    x2(1) = d3(2)*d1(3) - d3(3)*d1(2)
    x2(2) = d3(3)*d1(1) - d3(1)*d1(3)
    x2(3) = d3(1)*d1(2) - d3(2)*d1(1)

    s2    =  x2(1)**2 + x2(2)**2 + x2(3)**2
    if(s2.gt.smax) then
       imax = 2
       smax = s2
    endif

    x3(1) = d2(2)*d3(3) - d2(3)*d3(2)
    x3(2) = d2(3)*d3(1) - d2(1)*d3(3)
    x3(3) = d2(1)*d3(2) - d2(2)*d3(1)

    s3    =  x3(1)**2 + x3(2)**2 + x3(3)**2
    if(s3.gt.smax) imax = 3

    goto (100,200,300),imax
100 continue
    s1 = 1.0d0 / dsqrt(s1)
    ax(1) = x1(1) * s1
    ax(2) = x1(2) * s1
    ax(3) = x1(3) * s1
    sn    = 1.d0 / dsqrt( ax(2)**2+ ax(3)**2)
    tx(1) = 0
    tx(2) = ax(3) * sn
    tx(3) =-ax(2) * sn
    call matvec( rmat, tx, ty)
    co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
    goto 999
200 continue
    s2 = 1.0d0 / dsqrt(s2)
    ax(1) = x2(1) * s2
    ax(2) = x2(2) * s2
    ax(3) = x2(3) * s2
    sn    = 1.d0 / dsqrt( ax(2)**2+ ax(3)**2)
    tx(1) = 0
    tx(2) = ax(3) * sn
    tx(3) =-ax(2) * sn
    call matvec( rmat, tx, ty)
    co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
    goto 999
300 continue
    s3 = 1.0d0 / dsqrt(s3)
    ax(1) = x3(1) * s3
    ax(2) = x3(2) * s3
    ax(3) = x3(3) * s3
    sn    = 1.d0 / dsqrt( ax(1)**2+ ax(3)**2)
    tx(1) = ax(3) * sn
    tx(2) = 0.0
    tx(3) =-ax(1) * sn
    call matvec( rmat, tx, ty)
    co    = tx(1)*ty(1) + tx(2)*ty(2) + tx(3)*ty(3)
999 continue

    ang = dacos(co)
    ax(1) = ax(1) * ang
    ax(2) = ax(2) * ang
    ax(3) = ax(3) * ang

    ! ------------- check sign ! ----------
    call rotax( ax, 1.0d0, rtest, angdeg)
    err =                                                            &
         &       (rtest(1,2)-rmat(1,2))**2+                                 &
         &       (rtest(1,3)-rmat(1,3))**2+                                 &
         &       (rtest(2,3)-rmat(2,3))**2

    if(err.gt.1d-10) then
       !        write(6,*)'++++++ rotate sign +++++'
       ang = -ang
       ax(1) =-ax(1)
       ax(2) =-ax(2)
       ax(3) =-ax(3)
    endif

    return
  END subroutine axrot


  subroutine cross ( v1, v2, vc )
    !      ===============================
    !
    ! ** cross product of   v1 x v2  ==> vc
    !
    implicit none

    double precision   :: v1(3),v2(3),vc(3),vch(3)
    !
    vch(1) = v1(2)*v2(3) - v1(3)*v2(2)
    vch(2) = v1(3)*v2(1) - v1(1)*v2(3)
    vch(3) = v1(1)*v2(2) - v1(2)*v2(1)
    !
    vc(1)  = vch(1)
    vc(2)  = vch(2)
    vc(3)  = vch(3)

    return
  END subroutine cross


  subroutine norm( xin, xout, n)
    !      ==============================
    ! --- vektor norm nach skalarprodukt ---
    implicit none

    integer :: n
    double precision   :: xin(n), xout(n)
    double precision   :: s
    integer :: i

    s = 0.d0
    do i=1,n
       s = s + xin(i)**2
    end do
    if(s.gt.0.d0) s = 1.d0/dsqrt(s)

    do i=1,n
       xout(i) = s * xin(i)
    end do

    return
  END subroutine norm


  integer function laenge( string, maxlen, delim)
    !     ===============================================
    !     laenge von string, wobei das ende durch delim(1) gekennzeichnet is
    !     ------------------------------------------------------------------
    implicit none
    character(*) :: string
    integer      :: maxlen
    character(1) :: delim
    integer      :: i, n

    if(maxlen.gt.cmd_long) then
       n = cmd_long
    else
       n = maxlen
    endif

    do i=1,n
       if(string(i:i).eq.delim) then
          laenge = i-1
          return
       endif
    end do

    laenge = n

    return
  END function laenge

  !***********************************************************************
  !---------- formula-interpreter-section --------------------------------
  !***********************************************************************
  subroutine foinit
    !      ----------------
    !====================================================================
    ! formelinterpreter kern
    ! -----------------
    !====================================================================
    implicit none
    ! ---- internal use ---
    character*4 op
    integer     prio
    double precision        val
    logical     yesno
    integer     i,j,n

    ! ---- inits -----------------------------------------------------------
    ok          = .true.
    error       = .false.
    typ         = '    '
    topnumstack = 0
    topopstack  = 0
    klammerprio = 0
    actchar     = 0
    len_form    = lclen(formula,maxformlength)
    !ccw   write(6,*) 'len=',len_form
    delims(0)='('
    delims(1)=')'
    delims(2)='+'
    delims(3)='-'
    delims(4)='*'
    delims(5)='/'
    delims(6)='^'
    delims(7)=','
    !      do i=0,maxitemlength
    !         item(i) = ' '
    !      enddo
    return
    ! ---------------------------------------------------------------------
    entry sayit(yesno)
    !      ------------------
    say = yesno
    return
    ! ---------------------------------------------------------------------
    entry setdeg
    !      ------------
    degree = datan(1.d0)/45.d0
    write(6,*)'evaluator set to degree'
    return
    ! ---------------------------------------------------------------------
    entry setrad
    !      ------------
    degree = 1.d0
    write(6,*)'evaluator set to rad'
    return
    ! ---------------------------------------------------------------------
    entry show
    !      ----------
    if(.not.say) return
    write(6,*) (formula(i),i=0,len_form)
    write(6,*) 'actchar=',actchar, '  klammerprio=',klammerprio
    write(6,*) 'item   =',item
    do i=1,topopstack
       write(6,*)i,'  ',opstack(i),'   ',priostack(i)
    enddo
    write(6,*)
    do i=1,topnumstack
       write(6,*)i,'  ',numstack(i)
    enddo
    return
    ! ---------------------------------------------------------------------
    entry putopstack( op, prio)
    !      -----------------------------

    !ccw   write(6,*)'putopstack:',op,' ',prio
    if(topopstack.le.maxopstack) then
       topopstack = topopstack + 1
       opstack(topopstack)   = op
       priostack(topopstack) = prio+klammerprio
       !ccw   write(6,*)'topopstack=',topopstack
    else
       if(say) write(6,*)'error: opstack overflow!'
       error = .true.
    endif

    return
    ! ---------------------------------------------------------------------
    entry enternumstack( val )
    !      --------------------------

    if(topnumstack.le.maxnumstack) then
       topnumstack = topnumstack + 1
       numstack(topnumstack) = val
    else
       if(say) write(6,*)'error: numstack overflow!'
       error = .true.
    endif

    return
    ! ---------------------------------------------------------------------
    entry checknum( n )
    if(topnumstack.lt.n) then
       if(say) write(6,*)'error: too few num operands!'
       error = .true.
       do i=max(1,topnumstack),n
          numstack(i) = 0.d0
       enddo
    endif
    return
    ! ---------------------------------------------------------------------
    entry getword
    litem =-1
    do i = actchar,len_form
       do j=0,nodelims
          if(formula(i).eq.delims(j)) goto 100
       enddo
       litem = litem+1
       if(litem.gt.maxitemlength) goto 999
       item(litem) = formula(i)
    enddo
    actchar = len_form+1
100 continue
    actchar = i
    if(litem.gt.maxitemlength) goto 999
    item(litem+1) = ' '
    !ccw   write(6,*)'getitem: item=',(item(j),j=0,maxitemlength)
    return

999 continue
    if(say) write(6,*)'getitem: item too long'
    error = .true.
    return

  END subroutine foinit


  subroutine getitem
    !      ------------------

    ! --- internal use ---
    implicit none

    integer, parameter :: klammerinc=10
    integer, parameter :: iplusprio=1, minusprio=1, multprio=2, idivprio=2
    integer, parameter :: iexpprio=3, iuprio=7,komprio=0

    integer j, intn, ier, ierr, l, ll, lit

    logical anklam

    if(actchar.gt.len_form) then
       typ = 'end '
       call putopstack('end ',0)
       return
    endif

    if(formula(actchar).eq.'(') then
       typ = 'klam'
       klammerprio = klammerprio+klammerinc
       actchar     = actchar + 1
       goto 100
    endif
    if(formula(actchar).eq.')') then
       typ = 'klam'
       klammerprio = klammerprio-klammerinc
       actchar     = actchar + 1
       goto 100
    endif
    if(formula(actchar).eq.',') then
       typ = 'list'
       call putopstack('nop ',komprio  )
       actchar     = actchar + 1
       goto 100
    endif
    if(formula(actchar).eq.'+') then
       anklam = (actchar.eq.0)
       if(.not.anklam) anklam = (formula(actchar-1).eq.'(')
       if(anklam) then
          actchar     = actchar + 1
       else
          typ = 'bin '
          call putopstack('+   ',iplusprio)
          actchar     = actchar + 1
       endif
       goto 100
    endif
    if(formula(actchar).eq.'-') then
       anklam = (actchar.eq.0)
       if(.not.anklam) anklam = (formula(actchar-1).eq.'(')
       if(anklam) then
!?          call putopstack('neg ',iuprio)
          call putopstack('neg ',iuprio-1)    !! this is to fix -ln(x), and other fknt. TBD check with Klammerprio..
          actchar     = actchar + 1
       else
          typ = 'bin '
          call putopstack('-   ',minusprio)
          actchar     = actchar + 1
       endif
       goto 100
    endif
    if(formula(actchar).eq.'*') then
       typ = 'bin '
       call putopstack('*   ',multprio)
       actchar     = actchar + 1
       goto 100
    endif
    if(formula(actchar).eq.'/') then
       typ = 'bin '
       call putopstack('/   ',idivprio)
       actchar     = actchar + 1
       goto 100
    endif
    if(formula(actchar).eq.'^') then
       typ = 'bin '
       call putopstack('^   ',iexpprio)
       actchar     = actchar + 1
       goto 100
    endif
    ! ---  suche bis zum naechsten delimiter ---
    citemx = ' '
    call getword
    ! --- ist item ein unaerer operator ? ---
    if( compare(item,'sin ') ) then
       typ = 'unae'
       call putopstack('sin ',iuprio)
       goto 100
    endif
    if( compare(item,'cos ') ) then
       typ = 'unae'
       call putopstack('cos ',iuprio)
       goto 100
    endif
    if( compare(item,'tan ') ) then
       typ = 'unae'
       call putopstack('tan ',iuprio)
       goto 100
    endif
    if( compare(item,'atan ') ) then
       typ = 'unae'
       call putopstack('atan',iuprio)
       goto 100
    endif
    if( compare(item,'asin ') ) then
       typ = 'unae'
       call putopstack('asin',iuprio)
       goto 100
    endif
    if( compare(item,'acos ') ) then
       typ = 'unae'
       call putopstack('acos',iuprio)
       goto 100
    endif
    if( compare(item,'ln ') ) then
       typ = 'unae'
       call putopstack('ln  ',iuprio)
       goto 100
    endif
    if( compare(item,'exp ') ) then
       typ = 'unae'
       call putopstack('exp ',iuprio)
       goto 100
    endif
    if( compare(item,'sqrt ') ) then
       typ = 'unae'
       call putopstack('sqrt',iuprio)
       goto 100
    endif
    if( compare(item,'abs ') ) then
       typ = 'unae'
       call putopstack('abs ',iuprio)
       goto 100
    endif
    if( compare(item,'int ') ) then
       typ = 'unae'
       call putopstack('int ',iuprio)
       goto 100
    endif
    ! --- ist item eine userfunction ? -----
    if( formula(actchar).eq.'(') then
       typ = 'usrf'
       call putopstack('usrf',iuprio)
       if(tusrfstack.lt.musrfstack) then
          tusrfstack=tusrfstack+1
       else
          error = .true.
          write(6,*)'userfuncstack exceeded'
       endif
       usrfstack(tusrfstack) = citemx(1:19)//' '
       goto 100
    endif
    ! --- ist item eine zahl ? ----


    call re_scan(citemx ,valnum,ierr)

    ! --- try exp-num-representation --
    if(ierr.ne.0) then
       !ccw   write(6,*)'check exp-num item(litem)=',item(litem),'  ',litem
       !ccw   write(6,*)'formula(actchar)=',formula(actchar)
       lit = litem
       if((item(litem).eq.'e').or.(item(litem).eq.'d') .or.(item(litem).eq.'E').or.(item(litem).eq.'D')) then
          if((formula(actchar).eq.'+').or.(formula(actchar).eq.'-')) then
             lit = lit+1
             item(lit)=formula(actchar)
             do l=actchar+1,len_form
                do ll=0,nodelims
                   if(formula(l).eq.delims(ll)) goto 400
                enddo
                lit  =lit  +1
                item(lit  ) = formula(l)
             enddo
             l = len_form+1
400          continue

             call re_scan(citemx ,valnum,ierr)
             if(ierr.eq.0) then
                actchar = l
             else
                !cc              do l=litem+1,maxitemlength
                !cc                item(l)=' '
                !cc              enddo
                item(litem+1)=' '
             endif
          endif
       endif
    endif
    if (ierr.eq.0) then
       typ = 'num '
       call enternumstack(valnum)
       goto 100
    endif
    ! --- ist item eine referenz ? ----
    call extract(citemx,valnum,ierr)
    if (ierr.eq.0) then
       typ = 'num '
       call enternumstack(valnum)
       goto 100
    endif
    ! --- fehler:
    error = .true.
    if(say) write(6,*)'item:',item,' not decodable'

100 continue

    return

    ! ----------------------------------------------------------------------
    entry stackevaluate
    !      -------------------
    j = topopstack-1
    if(j.lt.0) return

1   continue
    if((j.gt.0).and.(priostack(j).ge.priostack(topopstack))) then

       if(opstack(j).eq.'nop ') then
          goto 200
       endif
       if(opstack(j).eq.'usrf') then
          call checknum(1)
          call newcom_usrfun(usrfstack(tusrfstack),numstack,topnumstack,ier)
          if(ier.ne.0) error = .true.
          call checknum(1)
          tusrfstack=tusrfstack-1
          goto 200
       endif
       if(opstack(j).eq.'+   ') then
          call checknum(2)

          if(topnumstack-1 < 1) then
            error = .true.
            goto 200
          endif

          numstack(topnumstack-1)=                                    &
               &      numstack(topnumstack-1)+numstack(topnumstack)
          topnumstack = topnumstack-1
          goto 200
       endif
       if(opstack(j).eq.'-   ') then
          call checknum(2)

          if(topnumstack-1 < 1) then
            error = .true.
            goto 200
          endif

          numstack(topnumstack-1)=                                    &
               &      numstack(topnumstack-1)-numstack(topnumstack)
          topnumstack = topnumstack-1
          goto 200
       endif
       if(opstack(j).eq.'*   ') then
          call checknum(2)

          if(topnumstack-1 < 1) then
            error = .true.
            goto 200
          endif

          numstack(topnumstack-1)=                                    &
               &      numstack(topnumstack-1)*numstack(topnumstack)
          topnumstack = topnumstack-1
          goto 200
       endif
       if(opstack(j).eq.'/   ') then
          call checknum(2)

          if(topnumstack-1 < 1) then
            error = .true.
            goto 200
          endif

          numstack(topnumstack-1)=                                    &
               &      numstack(topnumstack-1)/numstack(topnumstack)
          topnumstack = topnumstack-1
          goto 200
       endif
       if(opstack(j).eq.'^   ') then
          call checknum(2)

          if(topnumstack-1 < 1) then
            error = .true.
            goto 200
          endif

          numstack(topnumstack-1)=                                    &
               &      numstack(topnumstack-1)**numstack(topnumstack)
          topnumstack = topnumstack-1
          goto 200
       endif

       if(opstack(j).eq.'neg ') then
          call checknum(1)
          numstack(topnumstack)=    -numstack(topnumstack)
          goto 200
       endif
       if(opstack(j).eq.'sin ') then
          call checknum(1)
          numstack(topnumstack)=dsin(degree*numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'cos ') then
          call checknum(1)
          numstack(topnumstack)=dcos(degree*numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'tan ') then
          call checknum(1)
          numstack(topnumstack)=dtan(degree*numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'asin') then
          call checknum(1)
          numstack(topnumstack)=dasin(numstack(topnumstack))          &
               &                            /degree
          goto 200
       endif
       if(opstack(j).eq.'acos') then
          call checknum(1)
          numstack(topnumstack)=dacos(numstack(topnumstack))          &
               &                            /degree
          goto 200
       endif
       if(opstack(j).eq.'atan') then
          call checknum(1)
          numstack(topnumstack)=datan(numstack(topnumstack))          &
               &                            /degree
          goto 200
       endif
       if(opstack(j).eq.'ln  ') then
          call checknum(1)
          numstack(topnumstack)=dlog(numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'exp ') then
          call checknum(1)
          numstack(topnumstack)=dexp(numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'sqrt') then
          call checknum(1)
          numstack(topnumstack)=dsqrt(numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'abs ') then
          call checknum(1)
          numstack(topnumstack)=dabs(numstack(topnumstack))
          goto 200
       endif
       if(opstack(j).eq.'int ') then
          call checknum(1)
          intn = int(numstack(topnumstack))
          numstack(topnumstack)=intn
          goto 200
       endif
200    continue
       j = j-1
    else
       j = j+1
       opstack(j) = opstack(topopstack)
       priostack(j)=priostack(topopstack)
       topopstack = j
       return
    endif
    goto 1

  END subroutine getitem



  subroutine evaluate( f, val, ierr)
    !      ----------------------------------

    implicit none

    integer :: ierrx

    ! --- internal use ---
    double precision   val
    integer i, ierr
    character(len=*) f
    !
    do i=0,maxformlength
       formula(i) = ' '
    enddo
    do i=1,maxformlength
       if(f(i:i).eq.' ') goto 1
       formula(i-1)= f(i:i)
    enddo
1   continue


    call foinit
    do i=1,maxformlength
       call getitem
       call show
       if(typ.ne.'num ') then
          call stackevaluate
          call show
       endif
       if((typ.eq.'end ').or.error) goto 10
    enddo
10  continue
    if(error) then
       ierr = 1000
    else
       ierr = 0
    endif
    if(klammerprio.ne.0) then
       if(say)write(6,*)'unmatched parentheses'
       ierr = ierr+1
    endif
    if(opstack(1).ne.'end ') then
       if(say)write(6,*)'opstack items left!'
       ierr = ierr+10
    endif
    if(topnumstack.ne.1) then
       if(say)write(6,*)'numstack not ok'
       ierr = ierr+100
    endif

    val = numstack(1)


    call re_scan(cformulax,valnum,ierrx)
    if(ierrx.ne.0 .and. iout().ge.1) write(6,*)'evaluate: ',(formula(i),i=0,len_form),' to ',val
    return
  END subroutine evaluate



  subroutine cappend( a, b, c)
    !      ----------------------------
    implicit none
    character(*), intent(in)  ::  a,b
    character(len(a)+len(b)), intent(out) :: c

    integer i,j,l,n
    n = len(a)+len(b)

    c = ' '
    do i=1,n
       if(a(i:i).ne.' ') then
          c(i:i) = a(i:i)
       else
          goto 1
       endif
    enddo
1   continue
    j = i
    do i=j,n
       l=i-j+1
       if(b(l:l).ne.' ') then
          c(i:i) = b(l:l)
       else
          goto 2
       endif
    enddo
2   continue
    j = i
    do i=j,n
       c(j:j) = ' '
    enddo
  END subroutine cappend

  integer function lclen( string, m)
    !      ----------------------------------
    implicit none
    integer m
    character*1 string(0:m)
    do lclen=m,0,-1
       if(string(lclen).ne.' ') return
    enddo
  END function lclen

  logical function compare( s1, s2 )
    !      ----------------------------------
    implicit none
    character(1) ::  s1(*),s2(*)

    integer i
    compare = .true.
    do i=1,cmd_len-1
       if(s1(i).ne.s2(i)) then
          compare = .false.
          return
       endif
       if((s1(i).eq.' ').or.(s2(i).eq.' ')) return
    enddo
    compare = .false.
    return
  END function compare


  subroutine extract(nam,val,ier)
    !     -------------------------------

    implicit none
    integer i, ier
    character(*) :: nam
    double precision   val


    ier = 0
    ! --- look in uservars ---
    do i=1,nousev
       if(compare(nam,usenam(i))) then
          val = useval(i)
          return
       endif
    enddo
    ! --- number of doargs ---
    if(compare(nam,'noargs ')) then
       val = ipmlst(ktop)
       return
    endif
    ! --- user defined vals ---
    call newcom_usrextr(nam,val,ier)
    return
  END subroutine extract


  subroutine setudf(nam,val,ier)
    !     ------------------------------
    implicit none
    character(len=*), intent(in) :: nam
    double precision, intent(in) :: val
    integer, intent(out) :: ier

    integer i,j

    ier = 0
    ! --- look in uservars ---
    !ccw  write(6,*)'setudf: nam=',nam,' val=',val
    do i=1,nousev
       if(compare(nam,usenam(i))) then
          useval(i) = val
          if(iot.ge.1)                                                  &
               &    write(6,*)'uservar:',usenam(i),' is set to:',useval(i)
          return
       endif
    enddo
    !ccw  write(6,*)'create new udf'
    ! --- create new var ---
    if(nousev.lt.musevar) then
       nousev = nousev+1
       do i=1,len(usenam(nousev))
          usenam(nousev)(i:i)= nam(i:i)
          if(nam(i:i).eq.' ') then
             do j=i,len(usenam(nousev)) !16
                usenam(nousev)(j:j) = ' '
             enddo
             goto 1
          endif
       enddo
1      continue
       useval(nousev) = val
       i = nousev
       if(iot.ge.1)                                                    &
            &  write(6,*)'uservar:',usenam(i),' is created:',useval(i)
       return
    else
       write(6,*)'uservar set: buffer is full ! error'
       ier = 200
    endif
    return
  END subroutine setudf


  subroutine clrudf(nam)
    !     ----------------------
    implicit none
    character(len=*), intent(in) :: nam

    integer i,j
    !character(2*cmd_len) nam

    if(nam(1:3).eq.'all') then
       nousev = 0
       write(6,*)'uservar stack cleared'
    else
       do i=1,nousev
          if(compare(nam,usenam(i))) then
             write(6,*)'uservar:',usenam(i),'=',useval(i),' removed!'
             do j=i+1,nousev
                usenam(j-1)=usenam(j)
                useval(j-1)=useval(j)
             enddo
             nousev = nousev-1
             return
          endif
       enddo
    endif
    return
  END subroutine clrudf


  subroutine shwudf
    !     -----------------
    implicit none
    integer i

    if(nousev.eq.0) then
       write(6,*)'no uservars defined!'
    else
       write(6,*)'defined uservars:'
       do i=1,nousev
          write(6,*)usenam(i),useval(i)
       enddo
    endif
    return
  END subroutine shwudf

  !:::      integer function inpaf(iaddr)
  !:::!-----------------------------------------------------------------------
  !:::!  wert extraktion aus dem incom parameterstack  zuordnung inpar
  !:::!  still needed for one routine in datreat_subroutines
  !:::!-----------------------------------------------------------------------
  !:::       implicit none
  !:::       integer iaddr
  !:::       inpaf = inpar(iaddr)
  !:::      return
  !:::      END function inpaf


  FUNCTION intnxt(idef,inew)
    !     --------------------------
    !     extract value immediately following the last numerical value
    !     as integer, idef and inew like for intval
    !-----------------------------------------------------------------------
    implicit none
    integer     :: intnxt
    integer     :: idef
    integer     :: inew
    !-----------------------------------------------------------------------
    !
    inew      = 0
    intnxt    = idef
    if(ipars.le.lstpar) return
    !
    intnxt = NINT(rpar(lstpar+1))
    rused(lstpar+1) = .true.
    inew   = lstpar+1
    lstpar = lstpar

  END FUNCTION intnxt


  character(cmd_len) FUNCTION chrval(pname,cdef,inew0)
    !     --------------------------------------------
    !     extract a long name following (associated to) pname
    !     cdef=default, inew analogous to getval and intval
    !-----------------------------------------------------------------------
    !
    implicit none
    character(*), intent(in)          :: pname, cdef
    integer, intent(out) ,optional    :: inew0
    !
    !       logical cused, vused, rused, seten, compar
    !
    !       common/cincom/rpar(minc),inames,ipars &
    !        ,ioldc,inpar(minc),iparn(minc),inapa(minc),iargs,ipmls,iolbuf &
    !        ,lstpar, lstnam, cused, vused(minc), rused(minc), seten, info

    !       integer :: lclen
    integer :: i, inew
    !-----------------------------------------------------------------------
    !
    inew      = 0
    if(present(inew0)) inew0 = inew
    chrval    = cdef
    if(inames.le.0) return
    !
    search: do  i=1,inames
       if(compar(vname(i),pname//' ')) then
          vused(i) = .true.
          if ( len_trim(vname(i+1)) > 0 ) then
             chrval = vname(i+1)
             vused(i+1) = .true.
             inew = i
          else
             inew   = -i ! argument not found
          end if
          lstnam = i+1
          exit search
       endif
    end do search
    !
    if(present(inew0)) inew0 = inew
    return
    !
    !      end
  END FUNCTION chrval



  character(8) FUNCTION chrval8(pname,cdef,inew)
    !     --------------------------------------------
    !     extract a character*8 name following (associated to) pname
    !     cdef=default, inew analogous to getval and intval
    !-----------------------------------------------------------------------
    !
    implicit none
    character(*), intent(in) :: pname
    character(*), intent(in) :: cdef
    integer, intent(out)     :: inew
    !
    !       integer lclen
    integer :: i
    !-----------------------------------------------------------------------
    !
    inew      = 0
    chrval8    = cdef(1:8)
    if(inames.le.0) return
    !
    do i=1,inames-1
       if(compar(vname(i),pname//' ')) then
          vused(i) = .true.
          chrval8 = vname(i+1)(1:8)
          vused(i+1) = .true.
          inew = i
          lstnam = i+1
          goto 20
       endif
    end do
20  continue
    !
  END FUNCTION chrval8




  character(cmd_len) FUNCTION chrnxt(cdef,inew)
    !     --------------------------------------
    !     extract the name that follows the last extracted name immediately
    !     cdef=default, inew analogous to valnxt, intnxt.
    !-----------------------------------------------------------------------
    implicit none

    character(*), intent(in) :: cdef
    integer, intent(out)     :: inew
    !

    !       integer lclen
    !-----------------------------------------------------------------------
    !
    inew      = 0
    chrnxt    = cdef
    if(inames.le.lstnam) return

    chrnxt = vname(lstnam+1)
    vused(lstnam+1) = .true.
    inew   = lstnam + 1
    lstnam = inew

    !      end
  END FUNCTION chrnxt



  character(cmd_len) FUNCTION cmdf()
    !     ---------------------------
    !     returns the last command as function value.
    !
    implicit none

    cmdf = comand
    cused= .true.

  END FUNCTION cmdf



  logical FUNCTION iscmd(cname)
    !     -----------------------------
    !     returns true if the last command matches cname.
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(in) :: cname
    !
    !       integer clen, lclen, laenge

    if(info.ne.0) then
       iscmd = .false.
       write(6,'(a)') trim(cname)
       return
    endif

    if ( compar(comand,cname) ) then
       iscmd = .true.
       cused = .true.
    else
       iscmd = .false.
    endif

  END FUNCTION iscmd



  logical FUNCTION iscmdc(cname,explan)
    !     -------------------------------------
    !     returns true if the last command matches cname and gives
    !     explanation, explan is assumed to be terminated by $
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(in) :: cname
    character(*), intent(in) :: explan

    integer :: lexpl

    if(info.ne.0 ) then
       if (inames==0 .or. found(cname)) then
       iscmdc = .false.
       cused  = .true.
       !lexpl =laenge(explan,cmd_long,'$')
       lexpl = index(explan, '$')-1
       write(6,*)'-------------------command-------------------------'
       write(6,'(a)') trim(cname)
       write(6,*)'USAGE:'
       write(6,'(a)') trim(explan(1:lexpl))
       write(6,*)
       return
       endif
    endif

    if ( compar(comand,cname) ) then
       iscmdc= .true.
       cused = .true.
       if(iout().gt.0) then
          write(*,*)'performing: ',trim(cname)
       endif
    else
       iscmdc= .false.
    endif

  END FUNCTION iscmdc



  character(cmd_len) FUNCTION vnamef(iaddr)
    !     ----------------------------------
    !     returns the iaddr-th name that occured in the last command line
    !-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: iaddr
    !

    vnamef = vname(iaddr)
    lstnam = iaddr
    vused(iaddr) = .true.
    !       write(6,*)'vnamef,lstnam is: ', vnamef, lstnam

  END FUNCTION vnamef



  SUBROUTINE settit(newtit)
    !     -------------------------
    !     set the title string ( for plot etc.)
    !     the string should be terminatred by a '$' character !
    !-----------------------------------------------------------------------
    implicit none
    character(*) newtit
    integer :: ltit
    !
    !
    ltit  = laenge(newtit,cmd_len,'$')
    title = newtit(1:ltit)
    !
  END SUBROUTINE settit



  character FUNCTION titlef()
    !     ------------------------------
    !     returns the current title-string
    !
    implicit none


    titlef = title(1:cmd_len)

  END FUNCTION titlef



  FUNCTION rparf(iaddr)
    !     ---------------------
    !     returns the iaddr-th numerical value that occured in the last
    !     command line.
    !-----------------------------------------------------------------------
    implicit none
    double precision    :: rparf
    integer, intent(in) :: iaddr
    !

    if(iaddr > ipars) then
      rparf = 0d0                 ! may be replace by an "unknown" constant
      return
    endif

    rparf = rpar(iaddr)
    rused(iaddr) = .true.
    lstpar= iaddr

    return
    !      end
  END FUNCTION rparf



  FUNCTION inapf(iaddr)
    !     ---------------------
    !     returns the table-address of the numerical parameter that
    !     follows immediately to the iaddr-th name-parameter
    !-----------------------------------------------------------------------
    implicit none
    integer             :: inapf
    integer, intent(in) :: iaddr

    inapf = inapa(iaddr)

  END FUNCTION inapf



  FUNCTION inpaf(iaddr)
    !     ---------------------
    !     returns the number of numerical parameters following the iaddr-th
    !     name-parameter
    !-----------------------------------------------------------------------
    implicit none
    integer             :: inpaf
    integer, intent(in) :: iaddr

    inpaf = inpar(iaddr)

  END FUNCTION inpaf




  FUNCTION iparf()
    !     ----------------
    !     returns the total number of numerical parameters in the last
    !     command line
    !-----------------------------------------------------------------------
    implicit none
    integer :: iparf

    iparf = ipars

  END FUNCTION iparf



  FUNCTION inamf()
    !     ----------------
    !     returns the total number of name-parameters in the last command line
    !-----------------------------------------------------------------------
    implicit none
    integer :: inamf

    inamf = inames

  END FUNCTION inamf





 subroutine intnva(pname ,intdef ,intout, np, inew)
!     --------------------------------------------------
!     extract a secified number of numerical values that follow
!     the name pname as integers.
!     np values are extracted and returned on intout(1..np).
!     default values are to be given in intdef(1..np).
!     inew is analogous to inew in getval.
!-----------------------------------------------------------------------
       implicit none
       integer,          intent(in)        :: np
       character(len=*), intent(in)        :: pname
       integer,          intent(in)        :: intdef(np)
       integer,          intent(out)       :: intout(np)
       integer,          intent(out)       :: inew

       integer  :: i, j, ier, napr
       character(len=cmd_len) :: pnamc, napp
!
       inew      = 0
       do 1 i=1,np
         intout(i) = intdef(i)
1      continue
       if(inames.le.0) goto 20
!
       do 10 i=1,inames
         if(compar(vname(i),pname//' ')) then
            vused(i) = .true.
            napr = inpar(i)
            if(napr.gt.np) napr = np
            if(napr.gt.0) then
            do 11 j=1,inpar(i)
               intout(j) = nint(rpar(inapa(i)-1+j))
               rused(inapa(i)-1+j) = .true.
11          continue
            inew = napr
            lstpar = inapa(i)+inpar(i)-1
            lstnam = i
            endif
            goto 20
         endif
10     continue
20     continue

       do 30 i=1,np
        write(napp,'(i3,5x)')100+i
        napp(1:1)='.'
!        call capnd(pname//' ',napp,pnamc)
        pnamc = trim(pname)//trim(napp)
        call setudf(pnamc//' ',real(intout(i), kind=8),ier)
30     continue

!
      return

END SUBROUTINE intnva



!##################################################################################
!*+unix2c Converts Unix system time to date/time integer array.
      subroutine unix2c(utime, idate)
      implicit none
      integer, intent(in)   :: utime
      integer, intent(out)  :: idate(6)
!*utime  input  Unix system time, seconds since 1970.0
!*idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
!*-Author  Clive Page, Leicester University, UK.   1995-MAY-2
! -- 
! -- 
! Clive Page,                         e-mail:  cgp
! Dept of Physics & Astronomy,                 (at) le
! University of Leicester.                     (dot) ac (dot) uk 
      integer  :: mjday, nsecs
      real     :: day
!*Note the MJD algorithm only works from years 1901 to 2099.
      mjday    = int(utime/86400 + 40587)
      idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
      day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
      idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
      idate(3) = 1 + int(mod(day,30.6))
      nsecs    = mod(utime, 86400)
      idate(6) = mod(nsecs, 60)
      nsecs    = nsecs / 60
      idate(5) = mod(nsecs, 60)
      idate(4) = nsecs / 60

      end subroutine unix2c


  character(len=30) function u2c_datum(utime)
! ===========================================
   implicit none
   integer, intent(in)  :: utime   ! unix time: seconds since 1.1 1970
   integer              :: idate(6)
   character(len=30)    :: buffer

   call unix2c(utime, idate)
   write(buffer,'(i0.2,":",i0.2,":",i0.2," am ",i0.2,".",i0.2,".",i4," UTC ")') & 
               idate(4:6), idate(3), idate(2), idate(1)
   u2c_datum = buffer
  end function u2c_datum


  !> set command prompt
  subroutine set_prompt( s )
   implicit none
   character(len=*), intent(in) :: s
   prompt = trim(s)
  end subroutine set_prompt

  !> get command prompt
  function get_prompt() result( s )
   implicit none
   character(len=cmd_len) :: s
   s = prompt
  end function get_prompt




  ! (paz) whether "name" has been used
  function lvused(iname) result(res)
    integer, intent(in) :: iname
    logical :: res
    res = vused(iname)
  end function lvused







 ! =========================================================================
 !> convert blanks between quoites to ~ (aux function to enable newcom to deal 
 !> with quote strings as nametype items

    function markstring(str) result(marked)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str))      :: marked
         
      integer :: i
      logical :: instring 
      
      instring = .false.
      marked = str
      do i=1,len_trim(str)
         if(str(i:i) == '"') instring = .not. instring
         if(str(i:i) == ' '  .and.       instring) marked(i:i) = '~'
      enddo
    end function markstring      


 ! =========================================================================
 !> remove quotes and change  ~'s  to blanks
 !> (aux function to enable newcom to deal 
 !> with quote strings as nametype items

    function stripstring(str) result(stripped)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str))      :: stripped
         
      integer :: i, j

      stripped = " "
      j        = 0
      do i=1,len_trim(str)
         if( .not. str(i:i) == '"') then
           j = j+1
           stripped(j:j) = str(i:i)
           if(stripped(j:j) == '~') stripped(j:j) = ' '
         endif 
      enddo
    end function stripstring      

    function replacetabs(inpstr) result(outstr)
       implicit none
       character(len=1), parameter :: ch_tab = char(9)
       character(len=*), intent(in) :: inpstr
       character(len=len(inpstr))   :: outstr
       !
       integer :: i
       !
       outstr = repeat(' ', len(inpstr))
       do i=1,len_trim(inpstr)
          if(inpstr(i:i)/=ch_tab) outstr(i:i)=inpstr(i:i)
       end do
    end function replacetabs



END MODULE new_com
