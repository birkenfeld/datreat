program add_help_to_th
! small utility to enter necessary structure to add help info into datreat theos
! reads program source (f90) and outputs augmented one
implicit none
integer :: npar, i


character (len=256) :: inline, buf
integer, parameter  :: in=5, out=6

! assume first line is the function th_definition
!
     read(in,'(a)') inline
     write(out,'(a)') trim(inline)


d1: do 
     read(in,'(a)') inline
     if(inline(1:1) .ne. '!' .and. len_trim(inline) > 0) then
        write(out,'(a)')"      use theory_description"
        write(out,'(a)') trim(inline)
       exit d1
    else
       write(out,'(a)') trim(inline)
    endif
  
    enddo d1


d2: do 
     read(in,'(a)') inline
     write(out,'(a)') trim(inline)
     if(index(inline,"nparx")>0  .and. index(inline,"=")>0 ) exit d2
   enddo d2

   read(inline(index(inline,"=")+1:len_trim(inline)),*) npar


d2a: do 
     read(in,'(a)') inline
     write(out,'(a)') trim(inline)
!     if(index(inline,"npar")>0 .and. index(inline,"nparx")>0 .and. index(inline,"=")>0 ) exit d2a
!     if(index(inline,"npar = nparx")>0  ) exit d2a
      if(index(inline,"npar")>0 .and. index(inline,"nparx")>0 .and. index(inline,"=")>0 &
          .and. index(inline,")")==0 &
          .and. (index(inline,"npar") < index(inline,"=")) &
          .and. (index(inline,"=") < index(inline,"nparx"))) exit d2a
   enddo d2a

write(out,'(a,i4,a)')"! >>>>> describe theory with ",npar," parameters >>>>>>>"
   
write(out,'(a)')'        idesc = next_th_desc()'
write(out,'(a)')'        th_identifier(idesc)   = thnam'
write(out,'(a)')'        th_explanation(idesc)  = "DESCRIBR THEORY HERE                               "//cr//&'
write(out,'(a)')'                                 "CONTINUE HERE (max 1024 chars)                     "!'
write(out,'(a)')' !'
write(out,'(a)')'        th_citation(idesc) = "CITATIONS OF LIT HERE"'
 
! now look for the last parnam entry
     write(buf,'(i8)') npar
d3:  do 
       read(in,'(a)') inline
       write(out,'(a)') trim(inline)
       if(index(inline,"parnam") .ne.0 .and. index(inline,trim(adjustl(buf(1:8)))).ne.0) exit d3
     enddo d3   
! assuming that that was the last parnme
write(out,'(a)')"! >>>>> describe parameters >>>>>>>"
do i=1,npar
   write(out,'(a,i2,a)')'        th_param_desc(',i,',idesc) = " RARAMETER DESCRIPTION " !//cr//parspace//&'
enddo
write(out,'(a)')"! >>>>> describe record parameters used >>>>>>>"
write(out,'(a)')'       th_file_param(:,idesc) = " " '                                                                     
write(out,'(a)')'       th_file_param(1,idesc) = "FIRST PARAMETER DESC  " '                                                                     
write(out,'(a)')'       th_file_param(2,idesc) = "2nd   PARAMETER DESC  " '                                                                     
write(out,'(a)')'       th_file_param(2,idesc) = " ...  PARAMETER DESC  " '                                                                     
write(out,'(a)')"! >>>>> describe record parameters creaqted by this theory >>>>>>>"
write(out,'(a)')'       th_out_param(:,idesc)  = " "  '  
write(out,'(a)')'       th_out_param(1,idesc)  = "PARMETER CREATED "  '  
write(out,'(a)')'       th_out_param(2,idesc)  = "PARMETER CREATED "  '  
write(out,'(a)')'!'


! and add the rest 
 d4:  do 
       read(in,'(a)',end=999, err=999) inline
       write(out,'(a)') trim(inline)
     enddo d4   

999 continue

end program add_help_to_th
