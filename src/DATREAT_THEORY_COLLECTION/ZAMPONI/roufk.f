      function th18(x,pa,thnam,parnam,npar,ini)
c     ======================================== c
c
c -------> nrouse(t<taue)+fatkullin-kimmich(t>taue) <-------- 
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

      integer   labi,labf,pmod
      double precision  d,r2loc,sinc,taue,dx     
      
      common/thiadd/iadda

      data zpi/6.283185/
c     
c ----- initialisation -----       
      if(ini.eq.0) then          
         thnam = 'roufk' 
!         nparx = 12         
         nparx = 13         
         if(npar.lt.nparx) then            
            write(6,1)thnam,nparx,npar
 1          format(' theory: ',a8,' no of parametrs=',i8,      
     *           ' exceeds current max. = ',i8)            
            th18 = 0           
            return          
         endif          
         npar = nparx
c        --------------> set the number of parameters          
         parnam(1) = 'intensit'          
         parnam(2) = 'xi_frict'          
         parnam(3) = 'n_segmnt'          
         parnam(4) = 're      ' ! : re**2 = N * b**2          
         parnam(5) = 'temp    '          
         parnam(6) = 'diff    ' ! umfunktioniert         
         parnam(7) = 'q_width '          
         parnam(8) = 'lab_i   '          
         parnam(9) = 'lab_f   '          
         parnam(10)= 'd       '
         parnam(11)= 'tau_e   '
         parnam(12)= 'pmod    '
         parnam(13)= 'dx      ' !roehrendurchmesser fuer stern
c         
         th18 = 0          
         return        
      endif 
c       
      sqof0 = .false. 
c ---- calculate theory here -----       
      tau      = x        
      a0       = pa(1)       
      xi       = abs(pa(2))        
      N        = nint(pa(3))        
      R        = abs(pa(4))
      temp     = pa(5)       
      diff     = abs(pa(6))     ! in cm**2/sec
      q_width  = pa(7)       
      labi     = nint(pa(8))        
      labf     = nint(pa(9))        
      d        = pa(10)
      taue     = pa(11)
      pmod     = nint(pa(12))
      dx       = pa(13)

!      diff     = diff * 1d-9 / 1d-16 ! in A**2/ns        
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

      th18     = 0      
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
!            dr = diff                  

            call Nroufk(qzz,tau,temp,diff,xi,N,R,W,l,Sq,Sqt,labi,labf,d,
     &           taue,pmod,dx)
            sum = sum + fn*Sqt                  
            sumnorm = sumnorm + Sq*fn         
         endif        
      enddo

      if(sumnorm.gt.0.0d0) then           
         if(sqof0) then             
            th18 = sum
         else             
            th18 = sum/sumnorm           
         endif        
      else           
         th18 = 0       
      endif         
      th18 = th18*a0            

      call        parset('l       ',sngl(l),iadda,ier) ! in ns A units        
      call        parset('w       ',sngl(W),iadda,ier) !     "        
      call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "        
!      dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s        
!      call        parset('diff    ',sngl(dr),iadda,ier) 
c        
      return        
      end



      subroutine Nroufk(q,t,temp,diff,xi,N,R,W,l,Sq,Sqt,labi,labf,d,taue
     &     ,pmod,dx)
!      =============================================================== 
!
! Rouse expression for a chain of finite length: 
! Input parameters:
!    q     ----> momentum transfer in A**-1 
!    t     ----> time in nano-sec
!    temp  ----> temperature in K 
!!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation!    xi    ----> friction coefficient (Si: Ns/m=kg/s, hier: kg/ns)
!    N     ----> number of chain segments 
!    R     ----> end-to-end distance of the polymer molecule
! Output parameters: 
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l" 
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t) 
! ------------------------------------------------------------
!
!---mz---
!     Rouse modified for label structure: 
!     label starting at labi ending at labf
!     no cm diffusion
!     for t>taue Sinc after Fatkullin & Kimmich,
!     renormalized to Rouse


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

      integer   labi,labf,pmod
      double precision d,r2loc,r2,sinc,derfc,r2norm,sincnorm,dx     
      double precision lfk,Wfk,rounorm,taue,sqrou,sqtrou
! 
      double precision num, st, sumlik, sumlik0, st0, mue
      double precision x,tx
      double precision Hs,lab
      double precision nreal,ne,diff

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
      kbt = temp*kb             ! in Joule = kg*m**2/s**2        
      kbt = kbt * 100           ! in         kg*A**2/ns**2        
      W   = 3*kbt/(xi*(l**2))   ! in 1/ns 

! ---- set the diffusion constant if input is zero --- !
!       if(Dr.eq.0.0d0) then 
!         Dr = kbt/(N*xi)
!       endif
! mz-- no cm-diffusion        
      Dr=0.0d0                  

! ---- init sums ----       
      Sq0 = 0        
      Sq  = 0       
      Sqt = 0
      
!      open (10,file='test',status='unknown')
                     
      if(t.lt.taue) then
! ---- Do the sums -----
! mz-- sum only over labeled part 
!       do nn = 1,N        
!        do mm = 1,N        
         if(diff.eq.0) then
            do nn = labi,labf        
               do mm = labi,labf           
                  arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)           
                  arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)           
                  ff2  = -2*N*(l*q)**2/(3*pi**2)               
                  
                  arg2 = 0          
                  arg20= 0           
                  do p = 1,N             
!                  do p =n/(labf-labi),n
                     tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))             
                     a0    = -t*tau_p            
                     if(a0.lt.-200.0d0) a0 = -200.0d0             
                     e0    = 1.0d0-exp(a0)                          
                   ffc=cos((pi*p*nn)/dfloat(N))*cos((pi*p*mm)/dfloat(N))
! suppressing the larger modes    
!                     ffc   = ffc / (p**2)
!                     if(p.lt.(n/(labf-labi))) then
                     if(p.lt.pmod) then
                        ffc   = ffc / (3*p**2)
                     else
                        ffc   = ffc / (p**2)
                     endif
!     
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
         
         else  ! end label
            do nn = 1,N        
               do mm = 1,N
                  if((nn.lt.labi.and.mm.lt.labi).or.
     &                 (nn.gt.labf.and.mm.gt.labf)) then
!            do nn = labi,labf        
!               do mm = labi,labf           
                     arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)           
                     arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)           
                     ff2  = -2*N*(l*q)**2/(3*pi**2)               
                     
                     arg2 = 0          
                     arg20= 0           
                     do p = 1,N             
                        tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))             
                        a0    = -t*tau_p            
                        if(a0.lt.-200.0d0) a0 = -200.0d0             
                        e0    = 1.0d0-exp(a0)                          
                   ffc=cos((pi*p*nn)/dfloat(N))*cos((pi*p*mm)/dfloat(N))
                        if(p.lt.pmod) then
                           ffc   = ffc / (3*p**2)
                        else
                           ffc   = ffc / (p**2)
                        endif
!     
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
!                     write(10,*) nn,mm,sq,sqt,sqt/sq
                  endif
               enddo        
            enddo        
            
            Sq  = Sq /N       
            Sqt = Sqt/N               
            
         endif
      else
! mz-- fatkullin&kimmich term       
!     PE, 509K
!         lfk=3.612d0
!         Wfk=411.d0
!     l^2=17A^2! nicht 13A^2
         lfk=4.12d0
         Wfk=242.d0


         r2=2*sqrt(Wfk*lfk**4*t/pi)
         if(dx.ne.0) then
            if(r2.gt.dx**2) r2=dx**2 !beschraenkung fuer stern
         endif
         sinc=dexp(q**4*d**2*r2/216.d0)*derfc(q**2*d*sqrt(r2)/      
     &        (6.d0*sqrt(real(6))))        

! locale reptation, wrong!
!      r2loc=d*sqrt(real(2)/real(3))*(Wfk*lfk**4*t/pi)**0.25               
!      sinc=dexp(q**4*d**2*r2loc/216.d0)*derfc(q**2*d*sqrt(r2loc)/      
!     &     (6.d0*sqrt(real(6))))        

! normierung
         r2norm=2*sqrt(Wfk*lfk**4*taue/pi)
         sincnorm=dexp(q**4*d**2*r2norm/216.d0)*derfc(q**2*d*
     &        sqrt(r2norm)/(6.d0*sqrt(real(6))))        

         sqrou=0.
         sqtrou=0.
        
         if(diff.eq.0) then
            do nn = labi,labf        
               do mm = labi,labf           
                  arg1 = -(q**2)*(Dr*taue + abs(nn-mm)*(l**2)/6.0d0)           
                  arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)           
                  ff2  = -2*N*(l*q)**2/(3*pi**2)               
               
                  arg2 = 0          
                  arg20= 0           
                  do p = 1,N             
!                  do p =n/(labf-labi),n
                     tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))             
                     a0    = -taue*tau_p            
                     if(a0.lt.-200.0d0) a0 = -200.0d0             
                     e0    = 1.0d0-exp(a0)                          
                   ffc=cos((pi*p*nn)/dfloat(N))*cos((pi*p*mm)/dfloat(N))  
! suppressing the larger modes    
!                     ffc   = ffc / (p**2)
!                     if(p.lt.(n/(labf-labi))) then
                     if(p.lt.pmod) then
                        ffc   = ffc / (3*p**2)
                     else
                        ffc   = ffc / (p**2)
                     endif
!     
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
               
                  Sqrou  = Sqrou  + exp(aa1)           
                  Sqtrou = Sqtrou + exp(aa2)         
!                  Sq  = Sq  + exp(aa1)           
!                  Sqt = Sqt + exp(aa2)         
               enddo        
            enddo        
         else ! end label
            do nn = 1,N        
               do mm = 1,N
                  if((nn.lt.labi.and.mm.lt.labi).or.
     &                 (nn.gt.labf.and.mm.gt.labf)) then
                     arg1 = -(q**2)*(Dr*taue + abs(nn-mm)*(l**2)/6.0d0)           
                     arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)           
                     ff2  = -2*N*(l*q)**2/(3*pi**2)               
                     
                     arg2 = 0          
                     arg20= 0           
                     do p = 1,N             
                        tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))             
                        a0    = -taue*tau_p            
                        if(a0.lt.-200.0d0) a0 = -200.0d0             
                        e0    = 1.0d0-exp(a0)                          
                   ffc=cos((pi*p*nn)/dfloat(N))*cos((pi*p*mm)/dfloat(N))
                        if(p.lt.pmod) then
                           ffc   = ffc / (3*p**2)
                        else
                           ffc   = ffc / (p**2)
                        endif
!     
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
               
                     Sqrou  = Sqrou  + exp(aa1)           
                     Sqtrou = Sqtrou + exp(aa2)         
                  endif
               enddo        
            enddo        
         endif

         rounorm=sqtrou/sqrou     
            
         Sq  = 1.
         
         if(diff.eq.0) then
            Sqt = sinc*rounorm/sincnorm
         else
! end label, hdh12
!            nreal=925.d0
            nreal=932.d0
            lab=0.33d0
            ne=(d/lfk)**2
            mue=q**2*nreal*lfk**2/12.d0
            num = nreal/ne
            st = 0.5d0*((1.5d0/num)*(t/taue)**0.25)
            st0 = 0.5d0*(1.5d0/num)

            sumlik=(1.d0+dexp(-2*mue)+dexp(-4*mue*st)-4.*dexp(-2*mue*st)
     #           -2.*dexp(2*mue*(lab-1))+dexp(2*mue*(2*lab-1))
     #           +2.*dexp(-2*mue*lab)+4.*mue*(lab-st)+Hs(st-lab)*(-1.d0
     #           +4.*dexp(-2*mue*st)-2.*dexp(2*mue*(-2*st+lab))
     #           +dexp(4*mue*(lab-st))-2.*dexp(-2*mue*lab)
     #           +4.*mue*(st-lab)))/(2.*mue**2)
            
            sumlik0=(1.d0+dexp(-2*mue)+dexp(-4*mue*st0)-4.*dexp(-2*mue*
     #           st0)-2.*dexp(2*mue*(lab-1))+dexp(2*mue*(2*lab-1))
     #           +2.*dexp(-2*mue*lab)+4.*mue*(lab-st0)+Hs(st0-lab)*
     #           (-1.d0+4.*dexp(-2*mue*st0)-2.*dexp(2*mue*(-2*st0+lab))
     #           +dexp(4*mue*(lab-st0))-2.*dexp(-2*mue*lab)
     #           +4.*mue*(st0-lab)))/(2.*mue**2)
         
            sqt = sumlik/sumlik0*rounorm
         endif
!     
      endif
!      close(10)
         
      if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,Sq,Sqt, Sqt/Sq, w         
      return        
      end
