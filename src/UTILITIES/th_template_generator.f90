program th_template_generator
! author: michael monkenbusch, JCNS, FZ-Juelich
! create a template for a datreat theory
implicit none

 integer, parameter  :: namelen = 8     ! length of parameter an theory names
 integer, parameter  :: commlen = 256   ! length of help texts


 integer, parameter  :: mthpar  = 40    ! max number of parameters
 integer, parameter  :: mrecpar = 400   ! max number of record parameters input
 integer, parameter  :: moutpar = 400   ! max number of record parameters output

 integer             :: npar
 integer             :: nrecpar
 integer             :: noutpar

 
 character(len=namelen)    :: theoryname                      =" "
 character(len=4*commlen)  :: theory_description              =" "
 character(len=commlen)    :: theory_citation                 =" "
 character(len=namelen)    :: parametername(mthpar)           =" "
 character(len=commlen)    :: parameter_description(mthpar)   =" "
 character(len=namelen)    :: recparin(mrecpar)               =" "
 character(len=commlen)    :: recparin_default(mrecpar)       =" "
 character(len=commlen)    :: recparin_description(mrecpar)   =" "
 character(len=namelen)    :: recparout(mrecpar)              =" "
 character(len=commlen)    :: recparout_description(moutpar)  =" "

 character(len=commlen+namelen+1) :: inline, buf

 integer, parameter   :: in=5, out=10
 integer :: i,j

 !
 ! open output file
 !
   open(out,file="out.f90")
 !
 ! collect information on the theory !
 !
 ! we have to keep a given sequence starting with the name
 !
 
 ! read until keyword THEORY
 d1: do 
       read(in,'(a)',end=999,err=999) buf; inline = adjustl(buf)
       i = index(inline,"#THEORY ")
       if( index(inline,"#THEORY ") > 0) then
         theoryname = adjustl(trim(inline(i+7:len_trim(inline))))
         write(6,*)inline(1:80)
         write(6,*)"theory name = ",theoryname
         exit d1
       endif
     enddo d1
 !
 ! read now the theory description until the keyword #PARAMETERS occurs
 !
    theory_description = " "
d2: do
       read(in,'(a)',end=999,err=999)  buf; inline = adjustl(buf)
       i = index(inline,"#CITE")
       if(i>0) exit d2
       theory_description = trim(theory_description)//" "//trim(inline)
    enddo d2

    write(6,*)"thdesc: ",trim(theory_description)

d2a: do
       read(in,'(a)',end=999,err=999)   buf; inline = adjustl(buf)
       i = index(inline,"#PARAMETERS")
       if(i>0) exit d2a
       theory_citation = trim(theory_citation)//" "//trim(inline)
    enddo d2a

    write(6,*)"thcite: ",trim(theory_citation)
!
! read now the parameters until recin parameters
!
    npar = 0
d3: do
       read(in,'(a)',end=999,err=999)  buf; inline = adjustl(buf)
       i = index(inline,"#RECIN-PARAMETERS")
       if(i>0) exit d3
       npar = npar+1
       buf = adjustl(trim(inline))
       j=index(buf,"!")
       if(j==0) j=len_trim(buf)+1
       parametername(npar)         = buf(1:j-1)
       parameter_description(npar) = adjustl(trim(buf(j+1:len(buf))))  
       write(6,'(i3,a,3x,a)') npar,  parametername(npar),  trim(parameter_description(npar))    
    enddo d3

!
! read now the parameters until recoutparameters
!
    nrecpar = 0
d4: do
       read(in,'(a)',end=999,err=999)  buf; inline = adjustl(buf)
       i = index(inline,"#RECOUT-PARAMETERS")
       if(i>0) exit d4
       nrecpar = nrecpar+1
       buf = adjustl(trim(inline))
       recparin(nrecpar)         = buf(1:8)
       j=index(buf,"!")
       if(j==0) j=len_trim(buf)+1
       recparin_default(nrecpar) = buf(9:j-1)
       recparin_description(nrecpar) = adjustl(trim(buf(j+1:)))  
       write(6,'(i5," << ",a,3x,a)') nrecpar,  recparin(nrecpar),  trim(recparin_description(nrecpar))    
    enddo d4
!
! read now the parameters until recoutparameters
!
    noutpar = 0
d5: do
       read(in,'(a)',end=999,err=999)  buf; inline = adjustl(buf)
       i = index(inline,"#VARIABLES")
       if(i>0) exit d5
       noutpar = noutpar+1
       buf = adjustl(trim(inline))
       recparout(noutpar)         = buf(1:8)
       j=index(buf,"!")
       if(j==0) j=len_trim(buf)+1
       recparout_description(noutpar) = adjustl(trim(buf(j+1:)))  
       write(6,'(i5," << ",a,3x,a)') noutpar,  recparout(noutpar),  trim(recparout_description(noutpar))    
    enddo d5

!
! so far so good now we can write the header section
!
   call make_the_header

!
! read now 1 to 1 VARIABLE definitions
!
d6: do
      read(in,'(a)',end=999,err=999) inline
      i = index(inline,"#IMPLEMENTATION")
      if(i>0) exit d6
      write(out,'(a)') trim(inline)
    enddo d6
!
! start implementation with the init section
!
   call make_init_section
!
! read now 1 to 1 implementation
! 
d7: do
      read(in,'(a)',end=999,err=999) inline
      i = index(inline,"#END") 
      if(i>0) then
        write(out,'(3a)')"     th_",trim(adjustl(theoryname))," = th"
        write(out,'(a)')" "
        write(out,'(a)')"! ---- writing computed parameters to the record >>>  "
        do i=1,noutpar
            write(out,'(5a)')"      call parset('",recparout(i),"',sngl(",trim(recparout(i)),"),iadda,ier)"
        enddo 
        goto 1000
      endif
      i = index(inline,"#SUBROUTINES")
      if(i>0) exit d7
      write(out,'(a)') trim(inline)
    enddo d7

   write(out,'(3a)')"     th_",trim(adjustl(theoryname))," = th"
   write(out,'(a)')" "
   write(out,'(a)')"! ---- writing computed parameters to the record >>>  "
o1: do i=1,noutpar
      write(out,'(5a)')"      call parset('",recparout(i),"',sngl(",trim(recparout(i)),"),iadda,ier)"
    enddo o1

!
! read now 1 to 1 subroutines impelemns
!
  write(out,'(a)')" "
  write(out,'(a)')" CONTAINS "
  write(out,'(a)')" "
  write(out,'(a)')"! subroutines and functions entered here are private to this theory and share its variables "
  write(out,'(a)')" "
  
d8: do
      read(in,'(a)',end=999,err=999) inline
      i = index(inline,"#END")
      if(i>0) exit d8
      write(out,'(a)') trim(inline)
    enddo d8

1000 continue

 write(out,'(a,a)')" end function th_",trim(adjustl(theoryname))

 
 close(out)

 write(6,*)
 write(6,*)"The theory template is written to out.f90, complete it by editing and copy to th_",trim(theoryname),".f90"
 stop
999 continue
 write(6,*)"ERROR: premature end"
 close(out)


CONTAINS


  subroutine make_the_header
   implicit none
   integer :: i

   write(out,'(a,a,a)')" FUNCTION th_",trim(adjustl(theoryname)),"(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)"
   write(out,'(a)')    "!================================================================================"
   write(out,'(8("! ",128a))') trim(theory_description)
   write(out,'(2("! ",128a))') trim(theory_citation)
   write(out,'(4x,a)') "  use theory_description "
   write(out,'(4x,a)') "  implicit none "
   write(out,'(4x,2a)')"  real    :: th_",trim(adjustl(theoryname))
   write(out,'(4x,a)') "  character(len=8) :: thnam, parnam (*) "
   write(out,'(4x,a)') "  real    :: pa (*) "
   write(out,'(4x,a)') "  real    :: x , xh"
   write(out,'(4x,a)') "  integer :: mbuf, nparx, ier, ini, npar, iadda"
   write(out,'(4x,a)') "  integer, intent(inout) :: nopar       "                ! Anzahl der Parameter data
   write(out,'(4x,a)') "  character(len=80), intent(inout) :: napar(mbuf) "      ! name des parameters n
   write(out,'(4x,a)') "  real, intent(inout) :: params(mbuf) "                  ! value des parameters n
   write(out,'(4x,a)') " "
   write(out,'(4x,a)') "  double precision, parameter :: Pi = 4*atan(1d0)"
   write(out,'(4x,a)') "  integer                     :: actual_record_address"
   write(out,'(4x,a)') " "
   write(out,'(a)'   ) "! the internal parameter representation "
   do i=1,npar
      write(out,'(5x,a,a,a,a)') "double precision :: ",parametername(i),"   ! ",parameter_description(i)(1:80)
   enddo
   write(out,'(a)') "! the recin parameter representation "
   do i=1,nrecpar
      write(out,'(5x,a,a,a,a)') "double precision :: ",recparin(i),"   ! ",recparin_description(i)(1:80)
   enddo
   write(out,'(a)') "! the reout parameter representation "
   do i=1,noutpar
      write(out,'(5x,a,a,a,a)') "double precision :: ",recparout(i),"   ! ",recparout_description(i)(1:80)
   enddo
   write(out,'(a)') " "
   write(out,'(5x,a)') "double precision :: th"
   write(out,'(a)') " "


  end subroutine make_the_header




  subroutine make_init_section
    implicit none
    integer :: i
      write(out,'(a)')     "!"
      write(out,'(a)')     "! ----- initialisation ----- "
      write(out,'(a)')     "    IF (ini.eq.0) then     "
      write(out,'(a,a,a)') "       thnam = '",trim(adjustl(theoryname)),"'"
      write(out,'(a,i8)')  "       nparx = ",npar
      write(out,'(a)')     "       IF (npar.lt.nparx) then"
      write(out,*)         "          WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar"
      write(out,'(a,a,a)') "          th_",trim(adjustl(theoryname))," = 0" 
      write(out,'(a)')     "          RETURN"
      write(out,'(a)')     "       ENDIF"
      write(out,'(a)')     "       npar = nparx"
      write(out,'(a)')     "! >>>>> describe theory with >>>>>>> "
      write(out,'(a)')     "       idesc = next_th_desc()"
      write(out,'(a)')     "       th_identifier(idesc)   = thnam"
      write(out,'(4a)')    "       th_explanation(idesc)  = ",'"',trim(theory_description),'"'
      write(out,'(4a)')    "       th_citation(idesc)     = ",'"',trim(theory_citation),'"'
      write(out,'(a)')     "!       --------------> set the parameter names --->"
  
l1:   do i=1,npar
       write(out,'(a,i2,a,a,a,a)') "        parnam (",i,") = '",parametername(i),"'  ! ",parameter_description(i)(1:80)
      enddo l1     
      write(out,'(a)')     "! >>>>> describe parameters >>>>>>> "
l2:   do i=1,npar
        write(out,'(a,i2,3a)')     '        th_param_desc(',i,',idesc) = "',trim(parameter_description(i)),'" !//cr//parspace//&'
      enddo l2
      write(out,'(a)')     "! >>>>> describe record parameters used >>>>>>>"
      write(out,'(a)')     '        th_file_param(:,idesc) = " " '
l3:   do i=1,nrecpar
       write(out,'(a,i3,5a)')    '        th_file_param(',i,',idesc) = "',recparin(i),' > ',trim(recparin_description(i)),'"'
      enddo l3
      write(out,'(a)')     "! >>>>> describe record parameters creaqted by this theory >>>>>>> "
      write(out,'(a)')     '        th_out_param(:,idesc)  = " "'  
l4:   do i=1,noutpar
       write(out,'(a,i3,5a)')    '        th_out_param(',i,',idesc) = "',recparout(i),' > ',trim(recparout_description(i)),'"'
      enddo l4
      write(out,'(a)')     "! "
      write(out,'(3a)')    "        th_",trim(adjustl(theoryname))," = 0.0"
      write(out,'(a)')     " "
      write(out,'(a)')     "        RETURN"
      write(out,'(a)')     "     ENDIF"
      write(out,'(a)')     "!"
      write(out,'(a)')     "! ---- transfer parameters -----"

l5:   do i=1,npar
        write(out,'(6x,a,a,i2,a)') parametername(i)," =      pa(",i,")" 
      enddo l5

      write(out,'(a)') "! ---- extract parameters that are contained in the present record under consideration by fit or thc ---"
      write(out,'(a)')     "      iadda = actual_record_address()"
l6:   do i=1,nrecpar
        write(out,'(2a)')  "! >>> extract: ",trim(recparin_description(i))
        write(out,'(2a)')  "      xh = ",trim(recparin_default(i))
        write(out,'(5a)')  "      call parget('",recparin(i),"',xh,iadda,ier)"
        write(out,'(3a)')  "      ",recparin(i)," = xh"
      enddo l6



      write(out,'(a)')     "! "
      write(out,'(a)')     "! ------------------------------------------------------------------"
      write(out,'(a)')     "! ----------------------- implementation ---------------------------"
      write(out,'(a)')     "! ------------------------------------------------------------------"
      write(out,'(a)')     "! "



  end subroutine make_init_section


end program th_template_generator
