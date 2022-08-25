      function thrmo(x,pa,thnam,parnam,npar,ini,nopar,params,napar,mbuf)
c     ===================================================================
c
c -------> rousemod <--------
c
c
       implicit none
       character*8 thnam,parnam(20)
       real     thrmo, x, pa(20),qq(3)
       integer  mbuf,npar,ini,nparx,ier,i
       integer nopar ! Anzahl der Parameter data
       character*80 napar(mbuf)
       real*8 params(mbuf)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, l, Sq0, Sq, Sqt
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff, dr
       real*4    qget, tget, an
       integer   n, nqw
       logical   sqof0

       double precision Rp, Rw, f0, wl4

       integer iadda
       common/thiadd/iadda

       double precision zpi
       data zpi/6.283185d0/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'rousemod'
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           thrmo = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '
         parnam(3) = 'n_segmnt'   ! : if <0 return non-normalized s(q,t)
         parnam(4) = 're      '   ! : re**2 = N * b**2
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width '
         parnam(8) = 'rmode   '   ! Uebergangslaenge
         parnam(9) = 'rmode_wd'   ! Schaerfe des Uebergangs
         parnam(10)= 'flowlim '   ! Factor (< 1) fuer langwellige Moden (f=1--> keine Modenanalyse) 


c
         thrmo = 0
         return
       endif
c
       sqof0 = .false.
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       wl4      = abs(pa(2))
       an       = pa(3)
       N        = nint(abs(an))
       R        = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       Rp       = pa(8)
       Rw       = pa(9)
       f0       = pa(10)

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) then 
         write(6,*)'Assume Q as independend variable!'
         qz   = x
         tau = 0.0d0
         call parget('tau     ',qget,iadda,ier)
         if(ier.eq.0) then
          tau = qget
         endif
         write(6,*)'tau = ',tau 
         sqof0 = .true.
       endif
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif


      thrmo     = 0
      sum     = 0
      sumnorm = 0
      nqw     = 15
      dqw     = 4*q_width/nqw
      if(q_width.eq.0) then
        nqw = 0
        q_width = 1.0d0
      endif 

      do i=-nqw,nqw
       qzz = qz + i*dqw
       if(qzz.gt.0) then 
         fn   = dexp( -(i*dqw/q_width)**2)

c --- include center of mass diffusion ---
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
          dr = diff
       
          call Mrouse(qzz,tau,temp,Dr,Wl4,N,R,Rp,Rw,f0,W,l,xi,Sq,Sqt)

          sum = sum + fn*Sqt       
          sumnorm = sumnorm + Sq*fn

        endif
       enddo

       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            thrmo = sum
          else
            if(an.gt.0) then
              thrmo = sum/sumnorm
            else
              thrmo = sum
            endif
          endif
       else
          thrmo = 0
       endif 
       thrmo = thrmo*a0

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('xi      ',sngl(xi),iadda,ier)     !    SI
       dr        = dr /( 1d-9 / 1d-16 )                      ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
c
       return
       end




       Subroutine Mrouse(q,t,temp,Dr,Wl4,N,R,Rp,Rw,f0,W,l,xi,Sq,Sqt)
!      =============================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Wl4   ----> Rouse parameter in nsA units
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule   in A
!    Rp    ----> transition length for mode restriction        in A
!    Rpw   ----> width of transition                           in A
!    f0    ----> long mode rate-factor
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    xi    ----> friction coefficient (Si: Ns/m=kg/s, hier: kg/ns)
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt
       double precision Wl4, Rp,Rw,f0
       double precision pt, pw, f_mode

       integer N, nn,mm,p

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2

       integer iout
       
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2

!      W   = 3*kbt/(xi*(l**2))  ! in 1/ns
       W   = Wl4/l**4           ! in 1/ns    
       xi  = 3*kbt/(W*(l**2))   ! in SI


!      transition mode numbers
       if(Rp.gt.0d0 .and. Rw.gt.0d0) then
         pt = (R/Rp)**2
         pw = (R**2)/(2*Rw*Rp)
       else
         pt = 0d0
         pw = 1d0
       endif
!       write(6,'(a,2f12.6)')' pt,pw= ',pt,pw


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)
    
          arg2 = 0
          arg20= 0
          do p = 1,N
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))        ! this is the rate !!!
            tau_p = tau_p * f_mode(dfloat(p),f0,pt,pw)   ! mode rate modificator
            a0    = -t*tau_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)

            arg2  = arg2  + ffc*e0
            arg20 = arg20 + ffc

          enddo   
          arg2  = arg2  * ff2
          arg20 = arg20 * ff2

          aa1 = arg10
          if(aa1.lt.-300.0d0) aa1 = -300.0d0
          if(aa1.gt. 300.0d0) aa1 =  300.0d0

          aa2 = arg1+arg2
          if(aa2.lt.-300.0d0) aa2 = -300.0d0
          if(aa2.gt. 300.0d0) aa2 =  300.0d0


          Sq  = Sq  + exp(aa1)
          Sqt = Sqt + exp(aa2)

        enddo
       enddo

       Sq  = Sq /N
       Sqt = Sqt/N

       if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end



       double precision function f_mode(p,f0,pt,pw)
!      --------------------------------------------
!      Uebergangsfunktion 
!
       implicit none

       double precision f0, pt, pw, p
 
       if(pw.ne.0) then
        f_mode=(1d0-f0)/(1d0+exp((pt-p)/pw))+f0
       else
        if(p.lt.pt) f_mode = f0
        if(p.eq.pt) f_mode = f0+0.5*(1-f0)
        if(p.gt.pt) f_mode = 1d0
       endif
 
       return
       end
       
