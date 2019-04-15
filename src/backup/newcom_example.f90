!! ===========================================================
!! DrSPINE MAIN PROGRAM
!!
!! ===========================================================
program newcom_example

  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

  use new_com
  !
  character(len=cmd_len)  :: mycommand
  integer                 :: ier
  !
  external ::  myextract 
  external ::  myusrfun 


  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! this is for installation of the usr_extract program needed by new_com
  call set_usrextr( myextract )
  call set_usrfunc( myusrfun  )

  ! ------------------------------------------------------------------------

  !! add tab expansions as one adds commands
  istat = add_tab_expansion('help')
  istat = add_tab_expansion('clear')
  istat = add_tab_expansion('datapath')
  !

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ------ THE LOOP -----------------------------------------!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  info = 0
  mycommand = "  "
  commandloop: do

     !!>>>>
     !!> TO BE DONE: check here for errors issued in the last command and perform appropriate action if available
     !!> we may introduce an extraction function 'get_last_error_codes' reporting the error code accumulated in new_com either
     !!> during comand parsing/interpretation or obtained during calls of 'unused' and 'errsig'
     !!>>>>

     !------------------------------------------------------------
     !> get the next command line from keyboard or makro file
     !------------------------------------------------------------
     if(info == 0) then
        call get_newcom(mycommand, istat)  !!! read and parse command line
        if (istat<0) exit commandloop
        if(command("help ")) then
           info = 1
           cycle commandloop
        endif
     endif
     !-------------------------------------------------------------
     !> COMMAND: clear
     !! clear
     !-------------------------------------------------------------
     if(command('clear ','clear $')) then
       write(*,*)"Kommando CLEAR recieved "
        cycle commandloop
     endif


     !-------------------------------------------------------------
     !> COMMAND: .....
     !-------------------------------------------------------------
     if(command('quit   ','exit program $')) then
        exit commandloop
     endif


     ! ... and more entries in this sequence


     !-------------------------------------------------------------
     !> COMMAND: help
     !-------------------------------------------------------------
     !> switch off help mode
     if(info /= 0) then
        info = 0
        cycle commandloop
     endif

     !-------------------------------------------------------------
     !> COMMAND: check for the possibility that mycommad is a makro
     !-------------------------------------------------------------
     !> finally check for makro files
     call makro(mycommand)
     call unused( 1, 1, 1, ier)

  enddo commandloop

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ------ THE LOOP ENDS HERE--------------------------------!!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! The rest CONTAINS the implementation of commands         !!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program newcom_example


! CONTAINS


  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! this is a simple example for the program
  !! that allows to extract data
  !! that can be installed by the program to
  !! get the role of the former usrextract
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine myextract(name,x,ier)
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! if in an expression name matches, a corresponding value
    !! is returned on x
    !! trivial example: name='pi',  return x=3.14159
    !! if no match: ier = -ier
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    character(*),     intent(in)  :: name
    double precision, intent(out) :: x
    integer,          intent(out) :: ier

    !integer                       :: i

    if(trim(name).eq.'myextract') then
       write(6,*)'MYEXTRACT :',trim(name)
       x = 777d0
       ier = 0
    else
       ier = -1
    endif

  end subroutine myextract

  subroutine myusrfun( name, x, nx, ier )
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! if in an expression name(x1,x2,...x(nx)) name matches,
    !! a corresponding value which may depend on x(1)..x(nx)
    !! is returned on x(1)
    !! trivial example: name='sin',  return x(1)=sin(x(1))
    !! trivial example: name='add',  return x(1)=x(1)+x(2)
    !! if no match: ier = -ier
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    character(*),          intent(in)    :: name
    double precision,      intent(out)   :: x(*)
    integer,               intent(inout) :: nx
    integer,               intent(out)   :: ier

    integer                              :: i, j

    ! :: a dummy ::
    if(trim(name).eq.'myusrfun') then
       write(6,*)'MYUSEFUN :',trim(name), nx
       x(1) = 777d0
       ier = 0
    else
       ier = -1
    endif

 end subroutine myusrfun




