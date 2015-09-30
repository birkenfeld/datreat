       subroutine upcase(cx,lx)
c      ------------------------
c      Translate Uppercase Characters in String into lowercase
c      -------------------------------------------------------
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
 
c -- create the translation-table - only once at first call ----
       if(c(0).eq.'N') then
         write(6,*)'Translate-Table'
         do i=0,255
          c(i) = char(i)
          do j=1,26
            if(c(i).eq.lower(j:j)) c(i) = upper(j:j)
          enddo
         enddo
       endif
 
c -- perform the translation (always) --
       do j=1,lx
        test   = cx(j)
        cx(j)  = c(byte)
       enddo
 
       return
       end
 
 
 
       subroutine lowcase(cx,lx)
c      -------------------------
c      Translate lowercase Characters in String into Uppercase
c      -------------------------------------------------------
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
 
c -- create the translation-table - only once at first call ----
       if(c(0).eq.'N') then
         write(6,*)'Translate-Table'
         do i=0,255
          c(i) = char(i)
          do j=1,26
            if(c(i).eq.upper(j:j)) c(i) = lower(j:j)
          enddo
         enddo
       endif
 
c -- perform the translation (always) --
       do j=1,lx
        test   = cx(j)
        cx(j)  = c(byte)
       enddo
 
       return
       end
 
 
 
c----------------------------------------------------------------------
c  Group of Signal-routines that reckognize Crtl-C (SIGINT)
c
c  call sigset(-1)   to establish customized Crtl-C action
c
c  call sigset(0)    to restore system value
c
c  each time Crtl-C is pressed a counter is incremented
c  the value of this counter may be obtained via
c  function:
c            isigint
c
c  each time isigint is called, the counter is set to zero again
c
c----------------------------------------------------------------------
 
       subroutine sigset(iflag)
c      ------------------------
c      establish SIGINT (Crtl-C) for use in own fortran programs
c      if iflag = -1 --> SIGINT is flagged in common/sig/
c      if iflga =  0 --> Systemaction is restore
 
       integer SIGINT,ret, iflag, iflag0, signal
       parameter(SIGINT=2)
       common/sig/isignal
       external signal, IntHandler
 
       isignal = 0

       if(iflag.eq.-1) then
          iflag0 = -1
          write(6,*)'Signal with iflag0=',iflag0
          ret = signal( SIGINT, IntHandler, iflag )
          write(6,*)'Signal:',SIGINT,'  redefined, ret=',ret
       else
          iflag0 = 0
          ret = signal( SIGINT, IntHandler, iflag0)
          write(6,*)'Signal:',SIGINT,'  is set to system default',
     *              ', ret=',ret
       endif
 
       return
       end
 
       subroutine IntHandler(isig)
c      --------------------------
       common/sig/isignal

       isignal = isignal + 1
       write(6,*)'Sigint counter is now: ',isignal
       if(isignal.gt.5) then
         write(6,*)'Stop due to multiple Crlt-C'
         stop
       endif
       return
       end
 
       integer function isigint()
c      --------------------------
       common/sig/isignal
       isigint = isignal
       return
       end

       subroutine sig_reset()
c      --------------------------
       common/sig/isignal
       isignal = 0
       return
       end
 
 
       block data mysigi
c      -----------------
       common/sig/isignal
       data isignal/0/
       end
