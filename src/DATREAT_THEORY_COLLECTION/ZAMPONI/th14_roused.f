      function th14(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> rouse <--------
c
c
       character*8 thnam,parnam(20)
       dimension pa(20),qq(3)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8                            SQ_rouseT
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    epsilon, diff
       real*4    qget, tget
       real*4    kbolz
       Parameter (kbolz=1.380662e-23)
 
       common/thiadd/iadda

       data zpi/6.283185/
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'roused'
         nparx = 7
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th14 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'xi_frict'
         parnam(3) = 'b_segmnt'
         parnam(4) = 'epsilon '
         parnam(5) = 'temp    '
         parnam(6) = 'n       '
         parnam(7) = 'q_width '

c
         th14 = 0
         return
       endif
c
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       xi       = abs(pa(2))
       b        = abs(pa(3))
       epsilon  = abs(pa(4))
       temp     = pa(5)
       anmol    = abs(pa(6))  
       q_width  = pa(7)

     
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


       diff     = kbolz*temp*1e4/(anmol*xi*1e-20)  
       diff     = diff * 1d-9 / 1d-16          ! in A**2/ns
       call parset('diff    ',sngl(diff),iadda)

      th14    = 0
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
         sumnorm = sumnorm + fn

c --- include center of mass diffusion ---
         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
       
         if(pa(4).lt.0) then
           write(6,*)'Full rouse computation from scratch!'
           sum = sum + a * SQ_rouse(tau,qzz,temp,xi,b,epsilon)
         else
           sum = sum + a * SQ_rouseT(tau,qzz,temp,xi,b,epsilon)
         endif
        endif
       enddo

       th14 = sum/sumnorm
c
       return
       end
