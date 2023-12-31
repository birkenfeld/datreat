      function th_nrouse3(x,pa,thnam,parnam,npar,idum,ini)
c     ========================================
c
c -------> nrouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, l, Sq0, Sq, Sqt, wl4
       real*8    Rn
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff, dr
       real*4    qget, tget
       integer   n, ifix, nmin, Nn
       logical   sqof0

 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'nrousefn'
         nparx = 9
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_nrouse3 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4'
         parnam(3) = 'n_segmnt'
         parnam(4) = 're      '   ! : re**2 = N * b**2
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width '
         parnam(8) = 'fixend  '   ! : 0 normal Rouse, 1 fixed segment (N)
         parnam(9) = 'nmin    '   ! : minimum N for averaging lengths

c
         th_nrouse3 = 0
         return
       endif
c
       sqof0 = .false.
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       wl4      = abs(pa(2))  ! in A**4/ns
       N        = nint(pa(3))
       R        = abs(pa(4))  ! in A
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       ifix     = NINT(pa(8))
       nmin     = NINT(pa(9))

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


      th_nrouse3     = 0
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
    

   
          do Nn=nmin, N      ! loop over partial chains

            Rn = R*sqrt(dble(Nn)/dble(N))
            call NrouseY(qzz,tau,temp,dr,wl4,Nn,Rn, W, l,ifix, Sq,Sqt)
            sum = sum + fn*Sqt       
            sumnorm = sumnorm + Sq*fn

          enddo



        endif
       enddo

       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            th_nrouse3 = sum
          else
            th_nrouse3 = sum/sumnorm
          endif
       else
          th_nrouse3 = 0
       endif 
       th_nrouse3 = th_nrouse3*a0

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
c
       return
       end




       subroutine NrouseY(q,t,temp,Dr,wl4,N,R, W, l,ifix, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    ifix  ----> if 0: normal Rouse, else segment(N) fix
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

!       integer iout
       
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
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

    
       if(ifix.eq.0) then
         p0fix = 0
         pfac  = 1
       else
         p0fix = -0.5d0
         pfac  = 0.5d0   ! this is wrong ? !
         pfac=1         !! I think this has to be like this !!
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
          do ip=1,N
            p = ip+p0fix
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))
            a0    = -t*tau_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)
            ffc   = ffc*pfac

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

 !      if(iout().gt.1)write(6,'(1x,a,6E14.6)')
 !    *        'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end

