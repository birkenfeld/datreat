!> DrSPINE os-dependend system utilities
module os_utils

#ifdef __INTEL_COMPILER_BUILD_DATE
  use ifport            ! ifortran specific to improve system communication (form of system call, open...)
#endif

  private
  public :: execute_command
  public :: create_directory, change_directory
  public :: get_osname, os_openfile

contains
  !! =====================================================================================================
  !> execute system command
  !! command resurrected because ifort and gfortran seem to differ how to handle non-existing commands
  subroutine execute_command(cmd, exitstat, cmdstat, cmdmsg)
    character(len=*),  intent(in)    :: cmd
    integer, optional, intent(inout) :: exitstat
    integer, optional, intent(out)   :: cmdstat
    character(len=*), optional, intent(inout) :: cmdmsg
    !
    integer, parameter :: msg_len = 1024
    integer :: iexit
    integer :: istat
    character(len=msg_len) :: cmsg
    istat = 0
    cmsg  = ' '
#ifdef __INTEL_COMPILER_BUILD_DATE
    iexit = system(trim(cmd))
    if(iexit/=0) then
        istat=ierrno()
        write(cmsg,'(a,i0,a,i0)') 'exit code:', iexit, ', error code:', istat
    end if
#else
    call execute_command_line(trim(cmd), .true., iexit, istat, cmsg)
#endif
    if (present(exitstat)) exitstat= iexit
    if (present(cmdstat )) cmdstat = istat
    if (present(cmdmsg  ) .and. istat/=0) cmdmsg  = cmsg
  end subroutine execute_command


  !! =====================================================================================================
  !> change directory
  subroutine change_directory(cdirname, stat, pwd)
    character(len=*), intent(in)   :: cdirname
    integer, optional, intent(out) :: stat
    logical, optional, intent(in)  :: pwd
    !
    integer :: istat
    istat = chdir(cdirname)
    if (present(stat)) stat=istat
    if (istat.ne.0) return

    if (present(pwd)) then
       if (pwd) call execute_command_line('pwd', cmdstat=istat) ! F2008
       if (present(stat)) stat=istat
    end if
  end subroutine change_directory

  !! =====================================================================================================
  !> create directory
  subroutine create_directory(cdirname, writable, clean, stat)
    character(len=*), intent(in)   :: cdirname
    logical, optional, intent(in)  :: writable
    logical, optional, intent(in)  :: clean
    integer, optional, intent(out) :: stat
    !
    logical :: is_writable
    logical :: do_clean
    integer :: iex, ier
    ! get options
    is_writable = .true.
    do_clean    = .false.
    if (present(writable)) is_writable = writable
    if (present(clean))  do_clean = clean
    ! create directory
    call execute_command_line('mkdir -p '//trim(cdirname), exitstat=iex, cmdstat=ier)
    if (iex/=0) then
       if(present(stat)) stat = iex
       return
    endif
    if (ier/=0) then
       if(present(stat)) stat = ier
       return
    endif

    ! check if it is writable
    if(is_writable) then
       is_writable = is_writable_directory(trim(cdirname), stat=ier)
       if(ier/=0 .or. .not.is_writable) then
          if(present(stat)) stat = -1
          return
       end if
    end if
    if(do_clean) then
       call clean_directory(trim(cdirname), ier)
       if (ier/=0) then
          if(present(stat)) stat = ier
          return
       end if
    end if

    if(present(stat)) stat = 0
  end subroutine create_directory

  !! =====================================================================================================
  !> test is the directory is writable
  function is_writable_directory(cdirname, stat) result(res)
    character(len=*), intent(in)   :: cdirname
    integer, optional, intent(out) :: stat
    logical :: res
    !
    integer :: ier, iex
    res = .true.
    if(present(stat)) stat = 0

    call execute_command_line('touch '//trim(cdirname)//'/_is_writable', exitstat=iex, cmdstat=ier)
    if(ier/=0 .or. iex/=0) then
       res = .false.
       if(present(stat)) stat = -1
       return
    end if

    call execute_command_line('rm -f '//trim(cdirname)//'/_is_writable', exitstat=iex, cmdstat=ier)
    if(ier/=0 .or. iex/=0) then
       res = .false.
       if(present(stat)) stat = -1
       return
    end if
  end function is_writable_directory

  !! =====================================================================================================
  !> clean directory
  subroutine clean_directory(cdirname, stat)
    character(len=*), intent(in)   :: cdirname
    integer, optional, intent(out) :: stat
    !
    integer :: ier
    if(present(stat)) stat = 0
    call execute_command_line('rm -f '//trim(cdirname)//'/*', cmdstat=ier)
    if(ier/=0) then
       if(present(stat)) stat = -1
       return
    end if
  end subroutine clean_directory


  !! =====================================================================================================
  !> clean directory
  function get_osname() result(cname)
    character(len=32) :: cname
    !
    character(len=32), save :: os_name = ' '
    integer :: ier, iex
    !! FIXME - need a better way of capturing operating system (paz)
    if (len_trim(os_name)<=0) then ! not cached value
        call execute_command_line('uname -a|grep -iq darwin', exitstat=iex, cmdstat=ier)
        if (iex/=0 .or. ier/=0) then ! command failed, assume linux
            os_name = "linux"
        else
            os_name = "darwin"
        end if
    endif
    cname = os_name
  end function get_osname


  !! =====================================================================================================
  !> OS open file
  subroutine os_openfile(cfilename)
    character(len=*), intent(in) :: cfilename

#if   defined(__linux__)
    call execute_command_line('xdg-open '//trim(cfilename))
#elif defined(__APPLE__) || defined(__apple__)
    call execute_command_line('open     '//trim(cfilename))
#else
    character(len=32) :: osname

    osname = get_osname()
    if ( trim(osname) == "linux" ) then
        call execute_command_line('xdg-open '//trim(cfilename))
    else
        call execute_command_line('open     '//trim(cfilename))
    end if
#endif
  end subroutine  os_openfile


end module os_utils
