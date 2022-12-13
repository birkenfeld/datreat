      function th_radbrush(x,pa,thnam,parnam,npar,idum,ini)
c     =====================================================
c
c -------> see PRL 71, 1015 (1993) <--------
c
c
         implicit none

         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)

       character*8 thnam,parnam(20)
       REAL        th_radbrush, x, pa, qq, zpi, xh, vol_frac
       dimension pa(20),qq(3)
       data zpi/6.283185/
       real*8 q,tau,r0,sigmar0,rho_innen,ratio,aiqt,aiqt0,sq
       real(kind=4) :: aint, tskala, rskala, diffus, qp

       integer :: ini , nparx, npar, kmax, lmax, idum
       integer :: iflcorr, ninterp, isqflg
 
 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'brush'
!       nparx = 14
         nparx = 13
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'tskala'
         parnam(3) = 'rskala'
         parnam(4) = 'r0'
         parnam(5) = 'sigmar0'
         parnam(6) = 'rho_inne'
         parnam(7) = 'ratio*ts'
         parnam(8) = 'kmax'
         parnam(9) = 'lmax'
         parnam(10)= 'korrflag'
         parnam(11)= 'ninterp'
         parnam(12)= 'd/cm^2/s'
         parnam(13)= 'sq_flag'
!         parnam(14)= 'resolut'
c
         th_radbrush = 0
         return
       endif
c
c ---- calculate theory here -----
       aint      = pa(1)
       tskala    = pa(2)
       rskala    = pa(3)
       r0        = pa(4)
       sigmar0   = pa(5)
       rho_innen = pa(6)
       ratio     = pa(7)
       kmax      = pa(8) + 0.1
       lmax      = pa(9) + 0.1
       iflcorr   = pa(10)+ 0.1
       ninterp   = pa(11)+ 0.1
       diffus    = pa(12)
       isqflg    = pa(13)+ 0.1
 !      resol     = pa(14)
 
       ratio     = ratio / tskala
 
       if(isqflg.eq.0) then
          call        parget('q       ',qp,iadda,ier)
          if(ier.ne.0) then
            write(6,*)'th_radbrush-brush: q-not found! set q=0.1'
            qp = 0.1
          endif
          q   = qp * rskala
          tau = x / tskala
       else
          q   = x * rskala
          tau = 0.0d0
       endif
 
 
       call thbrush( q, tau, r0, sigmar0, rho_innen, ratio,
     *               kmax, lmax, iflcorr, ninterp,
     *               aiqt, aiqt0, sq, ier)
 
c --- einbeziehen einer globalen diffusion ?? ---
       if(diffus.ne.0.0) then
!          aiqt = aiqt * exp(-x  * (qp *qp * diffus * 1e7 + 1.0/resol) )
           aiqt = aiqt * exp(-x  * (qp *qp * diffus * 1e7 ) )
       endif
 
       if( isqflg.eq.0) then
 
          if(aint.eq.0.0) then
            th_radbrush = aiqt / aiqt0
          else
            th_radbrush = aint * aiqt
          endif
 
       else
          if(isqflg.eq.1) then
             th_radbrush = aint * aiqt0
          else
             th_radbrush = aint * sq
          endif
       endif
 
c
       return
       end
