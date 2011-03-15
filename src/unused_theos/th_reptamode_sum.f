

      function th_reptamode_sum(x,pa,thnam,parnam,npar,idum,ini)
c     ==========================================================
c
c -------> reptamode_sum <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,w
       real*8    R, Sq0, Sq, Sqt
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff, dr
       real*4    qget, tget
       integer   n,idum
       real*8    ptrans, pwidth, beta, alpha 
 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'reptamod'
         nparx = 11
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_reptamode_sum = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'    ! Intensität tau=0
         parnam(2) = 'w'           ! Rate 
         parnam(3) = 'n_segmnt'    ! number Chain elements
         parnam(4) = 're      '    ! : re**2 = N * b**2 in Angstroem
         parnam(5) = 'temp    '    ! temp 
         parnam(6) = 'com_diff'    ! Diffusion coefficient
         parnam(7) = 'q_width '    ! q width of detector (width for an averaging in q )
         parnam(8) = 'p_trans '    ! : transition between slow and fast regime
         parnam(9) = 'p_width '    ! : width of transition regime f=1/(1+exp((p-pmax)/p_width))
         parnam(10)= 'beta    '    ! : exponent steepness of transistion
         parnam(11)= 'alpha   '    ! : fraction of slow rate/p**2 to fast rate/p**2

c
         th_reptamode_sum = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       W        = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       ptrans   = abs(pa(8))
       pwidth   = abs(pa(9))
       beta     = abs(pa(10))
       alpha    = abs(pa(11))


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


      th_reptamode_sum     = 0
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

! --- get center of mass diffusion only effectve if diff .ne. 0  ---
         dr = diff
       
         call ReptaSum(qzz,tau,temp,Dr,W,N,R,
     *                     ptrans,pwidth,beta,alpha,Sq,Sqt)


         sum = sum + fn*Sqt       
         sumnorm = sumnorm + Sq*fn

        endif
       enddo

       if(sumnorm.gt.0.0d0) then
          th_reptamode_sum = a0*sum/sumnorm
       else
          th_reptamode_sum = 0
       endif 

       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
c
       return
       end



       subroutine ReptaSum(q,t,temp,Dr,W,N,R,
     *                     ptrans,pwidth,beta,alpha,Sq,Sqt)
!      ===============================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    W     ----> rate constant
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    ptrans----> mode nr. (transition)
!    pwidth----> width of modenr. transition
!    beta  ----> transition exponent
!    alpha ----> ratio slow/fast limit  (a+(a-1)) 
! Output parameters:
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,W,R,Sq,Sqt
       double precision ptrans, pwidth, beta, alpha
       integer N, nn,mm,p

       double precision l, tau_p, kbt, Sq0, arg1, arg2, nue
       double precision a0,e0, ff2, ffc,    arg10          ! ,arg20
       double precision aa1 , aa2, rsi
       double precision rate_p, p_spectrum
       double precision transf_p
       double precision fdf
  
       double precision ptrans_c, pwidth_c, fract_c, beta_c
       common /cpspec/ ptrans_c, pwidth_c, fract_c, beta_c



          
       ptrans_c = ptrans
       pwidth_c = pwidth
       fract_c  = alpha
       beta_c   = beta


       if(N.le.0) then
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

       
       nue = 0.5d0
! ---- determine the segment length l ----
       l = R/(dfloat(N)**nue)       
       
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       Rsi = R*1d-10            ! Re in meter                   

! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         ! do a linear interpolation between values for Gaussian and good solvent chain
         ! limting factors give in Doi-Edwards
         !        fdf= 0.196d0 + (0.203d0-0.196d0)*(nue-0.5d0)/(0.6d0-0.5d0)
         !        Dr = fdf*kbt/(Rsi*etasolv)      ! in m**2/s
         !        Dr = Dr*1d20/1d9                ! in A**2/ns
       endif

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0


!       write(6,*)'l=',l
!       write(6,*)'W=',w

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + (abs(nn-mm)**(2*nue))*(l**2)/6.0d0)
          arg10= -(q**2)*(       (abs(nn-mm)**(2*nue))*(l**2)/6.0d0)
!         ff2  = -2*N*(l*q)**2/(3*pi**2)
!         ff2  = -2*(dfloat(N)**(2*nue))*(l*q)**2/(3*pi**2)
          ff2  = -2*(R*q)**2/(3*pi**2)
    
          arg2 = 0
!         arg20= 0
          do p = 1,N
! rouse:    
!           rate_p = 2*W*(1-cos((pi*p)/dfloat(N)))                      ! hier rate in s**-1 nicht tau
            rate_p = W*pi/dfloat(N) * p_spectrum(dfloat(p))             ! hier rate in s**-1 nicht tau

! zimm:     tau_p = 0.325*etasolv*Rsi**3/kbt / 
!     *              (dfloat(p)**(3*nue)+alpha*dfloat(p)**(4-nue))      ! in s (?? 4-nue ok? fuer nue=1/2 ja)
!                                 --- since we use R-end-to-end durectly
!                                     here the chain statistics is implicit
!                   ----- this is the Gaussian chain factor (hope that it still is ok!
!           tau_p = tau_p * 1d9                                        ! in ns
 
            rate_p = rate_p*1d-9                                       ! in ns**-1 

!            write(6,*)'p r',p, rate_p                   

            a0    = -t*rate_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
! zimm      ffc   = ffc/(dfloat(p)**(2*nue+1.0d0)+alpha*dfloat(p)**(4))
            ffc   = ffc/(dfloat(p)**2)

            arg2  = arg2  + ffc*e0
!           arg20 = arg20 + ffc

          enddo   
          arg2  = arg2  * ff2
!         arg20 = arg20 * ff2

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

       return
       end


!!        double precision function transf_p(x,p,w)
!! !      -----------------------------------------
!!        implicit none
!!        double precision x, p, w, f, arg
!! 
!!        if(w.le.0.0d0) then
!!          if(x.gt.p) then
!!           f = 0.0d0
!!          else
!!           f = 1.0d0
!!          endif
!!        else
!!           arg = (x-p)/w
!!           if(arg.lt.-200d0) arg = -200d0
!!           if(arg.gt. 200d0) arg =  200d0
!!           f   = 1.0d0/(1.0d0+exp(arg))
!!        endif      
!! 
!!        transf_p = f
!! 
!!        return
!!        end
!! 
!! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Modifiziertes Modenspektrum                                       !!
!! motiviert durch Computersimulationen zur Reptation                !!
!! Limit p--> inf  --> p**2                                          !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       double precision function p_spectrum(p)
!      ---------------------------------------
!
       implicit none
       double precision p, y

       double precision ptrans_c, pwidth_c, fract_c, beta_c
       common /cpspec/ ptrans_c, pwidth_c, fract_c, beta_c

       y=fract_c+
     * (1-fract_c)*(0.5d0*(1.0d0+tanh((p-pwidth_c)/ptrans_c)))**beta_c
       
       p_spectrum = p*p*y

       return
       end

