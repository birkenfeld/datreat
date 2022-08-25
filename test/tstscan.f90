program testscan
implicit none
 
character(len=1) :: cn(64)
character(len=size(cn)) :: cnc
equivalence(cn,cnc)
double precision  :: x = 0
integer           :: ierr =0

  read(5,*) cnc
  write(6,*) cn
  call scan(cn,x,ierr)
  
  write(6,'(i8,e14.7,"  ",a)')ierr, x, cnc

  read(cnc,*,iostat=ierr) x
  
  write(6,'(i8,e14.7,"  ",a)')ierr, x, cnc


end program testscan



      subroutine scan(ctext,r,ierr)
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
      real*8 r
      integer ierr
      integer, parameter :: ctxtle = 64
      character*(ctxtle) ctext
      character*10       form


!                           * benamte programmkonstanten:
      integer, parameter :: cnil = -32766

!                           * programmvariablen:
      integer i,j
      integer   icur, iscntb(9,0:7), istate, isymbl


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
            read(ctext,form) r
         else
!                           * istate ist kein endzustand
            ierr = icur
         endif
      else
!                           * ermittle index des aktuellen symbols
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
      END

