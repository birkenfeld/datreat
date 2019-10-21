       function th_vrouse(x,pa,thnam,parnam,npar,idum,ini)
C      function th_nrouse(x,pa,thnam,parnam,npar,ini)
c     =============================================
c
c -------> nrouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, Sq0, Sq, Sqt,xi2,pc,rmesh2
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff, dr
       real*4    qget, tget
 !     integer   n,switch
       integer   n,switch,switch2



 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'vilgis'
         nparx = 12
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_vrouse = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'xi_frict'
         parnam(3) = 'n_segmnt'
         parnam(4) = 're      '   ! : re**2 = N * b**2
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width '
         parnam(8) = 'labi    '
         parnam(9) = 'labf    '
         parnam(10)= 'xi2     '
         parnam(11)= 'pc      '
         parnam(12)= 'rmesh2  '

c

         th_vrouse = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       xi       = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       switch   = nint(pa(8))
       switch2  = nint(pa(9))
       xi2      = pa(10)
       pc       = pa(11)
       rmesh2   = pa(12)

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) write(6,*)'Warning q not found' 
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif

      th_vrouse     = 0

      rou     = 0
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
       
          call Vrouse(qzz,tau,temp,dr,xi,N,R, W,Sq,
     _               Sqt,switch,switch2,xi2,pc,rmesh2)
          sum = sum + fn*Sqt       
          sumnorm = sumnorm + Sq*fn

        endif
       enddo

       if(sumnorm.gt.0.0d0) then
          th_vrouse = a0*sum/sumnorm
       else
          th_vrouse = 0
c          rou = 0
       endif 

       call        parset('w       ',sngl(W),iadda,ier)
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
c
       return
       end




       subroutine Vrouse(q,t,temp,Dr,xi,N,R, W,
     _                   Sq,Sqt,switch,switch2,xi2,pc,rmesh2)
!      =============================================================
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
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt,xi2,pc,rmesh2
       integer N, nn,mm,p

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20,derf
       double precision aa1 , aa2,W2,ffct0

       integer iout,switch,switch2
       
       

c *****************************************
       double precision star(199),sqinc,pinc
c *****************************************


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
       W    = 3*kbt/(xi*(l**2))  ! in 1/ns
       W2   = 3*kbt/(xi2*(l**2))


c ******************************************       
c       do p=1,199
c       star(p) = 0.d0
c       enddo


c       open(10,file='straube',status='old')
c       do p=1,199
c       read(10,*) star(p)
c       enddo
c       close(10)
c ******************************************



! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----

!       do nn = 1,N
!          do mm = 1,N
       do nn = switch,switch2
          do mm = switch,switch2
c          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
c          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
C          ff2  = -2*N*(l*q)**2/(3*pi**2)   
             ff2 = -(q**2)/6.d0
             arg2 = 0
             arg20= 0


C          do p = 1,N   

C        if(switch.eq.0) then
C           tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))
C            tau_p = p**2*pi**2*W/dfloat(N)**2  
C        else

 
C read from file (output from Straube program)             
C            tau_p=xi2*star(p)

c or put in test-function

           

c             tau_p = (W2*(atan((p*pi/dfloat(N)-pc*pi/dfloat(N))/step)
c     _               +pi/2.d0)+W*pi)/pi*(1.d0-cos(p*pi/dfloat(N)))

C        endif

C            a0    = -t*tau_p
C            if(a0.lt.-200.0d0) a0 = -200.0d0
             e0    = 1.0d-4
          
  
           
C            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
C            ffc   = ffc / (p**2)

             ffc = 3.d0*rmesh2*(1.d0-0.5*(2.*cosh(l**2/(3.*rmesh2)
     _            *abs(nn-mm))-exp(-abs(nn-mm)*l**2/(3.*rmesh2))*
     _            derf(1./(3.*rmesh2)*sqrt(W*l**4*t)-abs(nn-mm)*l**2
     _            /(2.*sqrt(W*l**4*t)))
     _            -exp(abs(nn-mm)*l**2/(3.*rmesh2))*
     _            derf(1./(3.*rmesh2)*sqrt(W*l**4*t)+abs(nn-mm)*l**2
     _            /(2.*sqrt(W*l**4*t)))))
             
             ffct0 = 3.d0*rmesh2*(1.d0-0.5*(2.*cosh(l**2/(3.*rmesh2)
     _            *abs(nn-mm))-exp(-abs(nn-mm)*l**2/(3.*rmesh2))*
     _            derf(1./(3.*rmesh2)*sqrt(W*l**4*e0)-abs(nn-mm)*l**2
     _            /(2.*sqrt(W*l**4*e0)))
     _            -exp(abs(nn-mm)*l**2/(3.*rmesh2))*
     _            derf(1./(3.*rmesh2)*sqrt(W*l**4*e0)+abs(nn-mm)*l**2
     _            /(2.*sqrt(W*l**4*e0)))))
             
             
             arg2  =  ffc
             arg20 = ffct0
    
C     enddo   
             arg2  = arg2  * ff2
             arg20 = arg20 * ff2
             
             aa1 = arg20
             if(aa1.lt.-300.0d0) aa1 = -300.0d0
             if(aa1.gt. 300.0d0) aa1 =  300.0d0
             
C     aa2 = arg1+arg2
             
             aa2=arg2
             
             if(aa2.lt.-300.0d0) aa2 = -300.0d0
             if(aa2.gt. 300.0d0) aa2 =  300.0d0
             
             
             Sq  = Sq  + exp(aa1)
             Sqt = Sqt + exp(aa2)
             
          enddo
       enddo

c************************************************************************
c calculate the incoherent scattering function

c init the p-sum
       
       pinc =  0.d0
       
       do p = 1,N
          
          
!          if(switch.eq.0) then
             tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))
C     tau_p = p**2*pi**2*W/dfloat(N)**2  
!          else
             
             
C     read from file (output from Straube program)    
!             tau_p=xi2*star(p)
             
c     or put in test-function
             
c     tau_p = (W2*(atan((p*pi/dfloat(N)-pc*pi/dfloat(N))/step)
c     _               +pi/2.d0)+W*pi)/pi*(1.d0-cos(p*pi/dfloat(N)))
             
!          endif
          
          pinc = pinc + 1.d0/p**2*(cos(p*pi/dfloat(N)))**2
     _         *(1.d0-exp(-t*tau_p))
          
          
       enddo   
       
       
       sqinc = exp(-2.d0/3.d0*q**2*dfloat(N)*l**2/pi**2*pinc)
       
       
c *************************************************************************

C put the coherent (switch2 = 0) or incoherent (switch2 = 1) function

!       if(switch2.eq.0) then

          Sq  = Sq /N
c       Sq = 1.d0
          Sqt = Sqt/N

!       else
          
!          Sqt = sqinc
!          Sq = 1.d0
          
!       endif
       
       if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,Sq,Sqt, Sqt/Sq, w 
       
       return
       end

