      function th19(x,pa,thnam,parnam,npar,ini)
c     ========================================
c
c -------> nrouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, l, Sq0, Sq, Sqt
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff, dr
       real*4    qget, tget
       integer   n
       logical   sqof0

       double precision taue,d,tau0,wloc,wfix

 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'cr'
!         nparx = 7
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th19 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'xi_frict'
         parnam(3) = 'b_segmnt'
         parnam(4) = 'epsilon '
!         parnam(4) = 'w_loc '
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width '
         parnam(8) = 'tau_e   '
         parnam(9) = 'd       '
         parnam(10) = 'wfix    '
c
         th19 = 0
         return
       endif
c
       sqof0 = .false.
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       xi       = abs(pa(2))
       l        = abs(pa(3))
       epsilon  = abs(pa(4))
!       wloc     = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       taue     = pa(8)
       d        = pa(9)
       wfix     = pa(10)

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

      th19     = 0
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
            a=fn*a0*dexp(-qzz*qzz*diff*tau)
            
!           call Nrousecr(qzz,tau,temp,dr,xi,W,l,Sq,Sqt,taue,d,tau0,wloc)
           call Nrousecr(qzz,tau,temp,dr,xi,W,l,Sq,Sqt,taue,d,tau0,wfix)
!            sum = sum + fn*Sqt       
            sumnorm = sumnorm + Sq*fn            
            sum = sum + a*Sqt       

         endif
      enddo
      
      if(sumnorm.gt.0.0d0) then
         if(sqof0) then
            th19 = sum
         else
            th19 = sum/sumnorm
         endif
      else
         th19 = 0
      endif 
      th19 = th19*a0

      call        parset('l       ',sngl(l),iadda,ier) ! in ns A units
      call        parset('w       ',sngl(W*1d31),iadda,ier) !     "
      call        parset('wl4     ',sngl(W*l**4*1d31),iadda,ier) !     "
      dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
      call        parset('diff    ',sngl(dr),iadda,ier)
      call parset('tau_0    ',sngl(tau0),iadda,ier)

c
      return
      end




!       subroutine Nrousecr(q,t,temp,Dr,xi,W,b,Sq,Sqt,taue,d,tau0,wloc)
       subroutine Nrousecr(q,t,temp,Dr,xi,W,b,Sq,Sqt,taue,d,tau0,wfix)
!      ==================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    xi    ----> friction coefficient (Si: Ns/m=kg/s, hier: kg/ns)
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
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

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt
       integer N, nn,mm,p

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2

       integer iout

       double precision epsilon,b,xifix,sq_rouset,taue,nu
       double precision tau0,t0,eqd,d,x,tx,derfc,loc,wfix
       double precision tnorm,txnorm,norm1,norm2,norm3,sqrou,sqtrou,wloc

       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
!       kbt = kbt*100          ! in         kg*A**2/ns**2
!       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

! ---- init sums ----
       Sq  = 0
       Sqt = 0
       
!       b=3.612d0
!       xifix=0.3925d8 ! wrong l^2, W
       xifix=1.3*0.3925d8
       epsilon=1.0d-3
!       wfix=411.d0 ! wrong l^2, W
!       wfix=242.d0 !parameter

!       if(t.lt.taue) then
!          sqt=sq_rouset(t,q,temp,xifix,b,epsilon)
!          sq=1.
!       else
          sqt=sq_rouset(t,q,temp,xi,b,epsilon)

          w   = 3*kbt/(xi*(b**2))
!          tau0 =  36.d0/(Wfix*((q*b)**4))
!          tau0 =  36.d0/(Wloc*((q*b)**4))
          tau0 =  36.d0/((wfix-w*1d31 )*((q*b)**4))
          t0 = t/tau0
          eqd = exp(-((q*d/6.d0)**2))
          
          if(t0.gt.50.d0) then
             x = t0 
      tx = 1/sqrt(0.3141592653589793D1)*sqrt(1/x)-1/sqrt(0.3141592653589
     #793D1)*sqrt(1/x)**3/2+3.D0/4.D0/sqrt(0.3141592653589793D1)*sqrt(1/
     #x)**5-15.D0/8.D0/sqrt(0.3141592653589793D1)*sqrt(1/x)**7
          else
             tx = exp(t0)*derfc(sqrt(t0))
          endif
        
          loc=(1.d0-eqd)*tx +eqd

! normierung
          norm1=sq_rouset(taue,q,temp,xifix,b,epsilon)

          tnorm=taue/tau0
          if(tnorm.gt.50.d0) then
             x = tnorm 
             txnorm = 1/sqrt(0.3141592653589793D1)*sqrt(1/x)-
     #          1/sqrt(0.3141592653589793D1)*sqrt(1/x)**3/2
     #          +3.D0/4.D0/sqrt(0.3141592653589793D1)*sqrt(1/x)**5
     #          -15.D0/8.D0/sqrt(0.3141592653589793D1)*sqrt(1/x)**7
          else
             txnorm = exp(tnorm)*derfc(sqrt(tnorm))
          endif
          
          norm2=(1.d0-eqd)*txnorm +eqd
!        
          norm3=sq_rouset(taue,q,temp,xi,b,epsilon)
!        
!          sqt=sqt*loc*norm1/(norm2*norm3)
          sqt=sqt*loc
          sq=1.
          
!       endif

       if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,Sq,Sqt,Sqt/Sq,w 

       return
       end

