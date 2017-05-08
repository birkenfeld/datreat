subroutine extract_th(filename)
! get theory section appended on msave generated files
! and write it to lastth
! useage: exctract_th <filename>
!
implicit none
character(len=*), intent(in) :: filename

character(len=256) :: line
integer            :: inunit, outunit
integer            :: i, length, status
character(len=256) :: val
  
open(newunit=inunit,file=trim(filename),iostat=status)
if(status /= 0) then
  call errsig(999,"ERROR #### extract_th: cannot open:"//trim(val)//"$")
  stop
else
  call errsig(999,"EXTRACTING lastth from appendix of file: "//trim(val)//"$")
endif

l1: do 
       read(inunit,'(a)',end=999) line
       if(line(1:7) == " theory") exit l1
    enddo l1

open(newunit=outunit,file="lastth")

write(outunit,'(a)') trim(line)
l2: do
       read(inunit,'(a)',end=998) line
       write(outunit,'(a)') trim(line)      
       if(adjustl(trim(line)) == "end") exit l2
    enddo l2

close(outunit)
close(inunit)
stop

998 continue
  call errsig(999,"ERROR #### extract_th: theory end not found! $")
  close(outunit)
  close(inunit)
  stop
999 continue
  call errsig(6,*)"ERROR #### extract_th: theory section not found! $")
  close(inunit)
end subroutine extract_th
