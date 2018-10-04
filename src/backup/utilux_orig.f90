       subroutine upcase(cx,lx)
!      ------------------------
!      Translate Uppercase Characters in String into lowercase
!      -------------------------------------------------------
       character cx(*)
       integer lx
       character*1   c(0:255), test
       character*26  upper,lower
       integer*1     byte
       equivalence(test,byte)
       save c
       data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
       data lower/'abcdefghijklmnopqrstuvwxyz'/
       data c(0)/'N'/

! -- create the translation-table - only once at first call ----
       if(c(0).eq.'N') then
         write(6,*)'Translate-Table'
         do i=0,255
          c(i) = char(i)
          do j=1,26
            if(c(i).eq.lower(j:j)) c(i) = upper(j:j)
          enddo
         enddo
       endif

! -- perform the translation (always) --
       do j=1,lx
        test   = cx(j)
        cx(j)  = c(byte)
       enddo

       return
      END



       subroutine lowcase(cx,lx)
!      -------------------------
!      Translate lowercase Characters in String into Uppercase
!      -------------------------------------------------------
       character cx(*)
       integer lx
       character*1   c(0:255), test
       character*26  upper,lower
       integer*1     byte
       equivalence(test,byte)
       save c
       data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
       data lower/'abcdefghijklmnopqrstuvwxyz'/
       data c(0)/'N'/

! -- create the translation-table - only once at first call ----
       if(c(0).eq.'N') then
         write(6,*)'Translate-Table'
         do i=0,255
          c(i) = char(i)
          do j=1,26
            if(c(i).eq.upper(j:j)) c(i) = lower(j:j)
          enddo
         enddo
       endif

! -- perform the translation (always) --
       do j=1,lx
        test   = cx(j)
        cx(j)  = c(byte)
       enddo

       return
      END



!----------------------------------------------------------------------
!  Group of Signal-routines that reckognize Crtl-C (SIGINT)
!
!  call sigset(-1)   to establish customized Crtl-C action
!
!  call sigset(0)    to restore system value
!
!  each time Crtl-C is pressed a counter is incremented
!  the value of this counter may be obtained via
!  function:
!            isigint
!
!  each time isigint is called, the counter is set to zero again
!
!----------------------------------------------------------------------


       subroutine sigset(iflag)
!      ------------------------
!      establish SIGINT (Crtl-C) for use in own fortran programs
!      if iflag = -1 --> SIGINT is flagged in common/sig/
!      if iflga =  0 --> Systemaction is restore

       integer SIGINT,ret, iflag, iflag0, signal
       parameter(SIGINT=2)
       common/sig/isignal
       external signal, IntHandler

       isignal = 0

       if(iflag.eq.-1) then
          iflag0 = -1
          write(6,*)'Signal with iflag0=',iflag0
!--          ret = signal( SIGINT, IntHandler, iflag )
          write(6,*)'Signal:',SIGINT,'  redefined, ret=',ret
       else
          iflag0 = 0
!--          ret = signal( SIGINT, IntHandler, iflag0)
          write(6,*)'Signal:',SIGINT,'  is set to system default',      &
     &              ', ret=',ret
       endif

       return
      END

       subroutine IntHandler(isig)
!      --------------------------
       common/sig/isignal

       isignal = isignal + 1
       write(6,*)'Sigint counter is now: ',isignal
       if(isignal.gt.5) then
         write(6,*)'Stop due to multiple Crlt-C'
         stop
       endif
       return
      END

       integer function isigint()
!      --------------------------
       common/sig/isignal
       isigint = isignal
       return
      END

       subroutine sig_reset()
!      --------------------------
       common/sig/isignal
       isignal = 0
       return
      END


       block data mysigi
!      -----------------
       common/sig/isignal
       data isignal/0/
      END
