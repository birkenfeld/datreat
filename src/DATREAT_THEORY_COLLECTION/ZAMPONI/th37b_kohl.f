      function th37(x,pa,thnam,parnam,npar,ini)
c     ===================================================
c
c -------> kohl <--------
c
c
       implicit none
       
       character*8 thnam,parnam(20)
   
       real*4    x, pa, qq, th37
       dimension pa(20),qq(3)
       integer   npar, ini, nparx

       integer mgaussians
       parameter(mgaussians=10)
       double precision Pi
       parameter(Pi=3.141592654d0)
 
       

       real*8    tau, tau0, beta, a0, epsilon
       real*8    qz
       real*4    tget, temp, qget

       real*4    gampli (mgaussians)
       real*4    gwidth (mgaussians)
       real*4    gcenter(mgaussians)
       integer   ng


       integer   maxit
       integer   ier
       integer   inew

       integer   i

       real*4    kbolz
       Parameter (kbolz=1.380662e-23)
 
       integer iadda
       common/thiadd/iadda


       real*8 omega0, u_sqr, dwf
       real*8 qcenter, qband, qexpt0, qexpbeta, bkgr        

       real*8 a,b, domega, o0
       real*8 adapint, sum, result, result2, err, erraccu
       real*8 strex_kernel
       external strex_kernel

! communication with Kernel
 

       real*8 Omega, str_beta, str_tau0, str_delta
       common /cstrex1/Omega, str_beta, str_tau0, str_delta



c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'kohl_q'
         nparx = 11
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th37 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          ! prefactor
         parnam(2) = 'tau0'              ! KWW time-constant, prefactor in front of q-dependence
         parnam(3) = 'beta'              ! streched exp, prefactor in front of-q-dependence
         parnam(4) = 'epsilon '          ! accuracy parameter for FT-integrations (DO NOT FIT)
         parnam(5) = 'omega0'            ! omega scale zero shift
         parnam(6) = 'u_sqr'             ! <u**2> value for Debye-Waller-Factor
         parnam(7) = 'q_center'          ! center of q_range for switching
         parnam(8) = 'q_band'            ! +-width of q_range for switching
         parnam(9) = 'qexp_t0'           ! q-exponent for tau0
         parnam(10)= 'qexp_bet'          ! beta-exponent 
         parnam(11)= 'bkgr'              ! constant background 



c
         th37 = 0
         return
       endif
c
c ---- calculate theory here -----
       o0           = x   -   pa(5)
       a0           = pa(1)
       str_tau0     = abs(pa(2))
       str_beta     = abs(pa(3))
       epsilon      = abs(pa(4))
       u_sqr        = abs(pa(6))
       qcenter      = pa(7)
       qband        = pa(8)
       qexpt0       = pa(9)
       qexpbeta     = pa(10)         
       bkgr         = pa(11)

     
       if(epsilon.eq.0.0d0) epsilon = 1.0d-8
       maxit = 1000

        qget = 0.0
        call        parget('q       ',qget,iadda,ier)
        qz   = qget
        if(ier.ne.0) write(6,*)'Warning q not found' 

       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif


!
! Check whether this instance of the theorie matches the q-range
!
       if(qz .lt. (qcenter-qband) .or. 
     *    qz .ge. (qcenter+qband)     ) then
!         Outside range: no contribution !       
!          write(6,*)'z ',qz,qcenter,qband
          th37 = 0.0
          return
       endif

! Apply the q-scalings
       if(qz.gt.0.0d0) then
         str_tau0 = str_tau0 * qz**qexpt0
         str_beta = str_beta * qz**qexpbeta
       endif
!
! push these into the parameter section
       call parset('tau0_q  ',sngl(str_tau0),iadda)
       call parset('beta_q  ',sngl(str_beta),iadda)



! extract gaussian parameters 
!
! assuming Gaussian parameters according to the form of th:gauss
!
!      g = ampli*exp(-((x-center)/width)**2)
!      -------------------------------------


       gampli(1)   = 1d0
       gwidth(1)   = 1d0
       gcenter(1)  = 0.0
       ng          = 0
       call        parget('ga1inten',gampli(1) ,iadda,ier)
       call        parget('ga1width',gwidth(1) ,iadda,ier)
       call        parget('ga1cente',gcenter(1),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif

       call        parget('ga2inten',gampli(2) ,iadda,ier)
       call        parget('ga2width',gwidth(2) ,iadda,ier)
       call        parget('ga2cente',gcenter(2),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga3inten',gampli(3) ,iadda,ier)
       call        parget('ga3width',gwidth(3) ,iadda,ier)
       call        parget('ga3cente',gcenter(3),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga4inten',gampli(4) ,iadda,ier)
       call        parget('ga4width',gwidth(4) ,iadda,ier)
       call        parget('ga4cente',gcenter(4),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga5inten',gampli(5) ,iadda,ier)
       call        parget('ga5width',gwidth(5) ,iadda,ier)
       call        parget('ga5cente',gcenter(5),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga6inten',gampli(6) ,iadda,ier)
       call        parget('ga6width',gwidth(6) ,iadda,ier)
       call        parget('ga6cente',gcenter(6),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga7inten',gampli(7) ,iadda,ier)
       call        parget('ga7width',gwidth(7) ,iadda,ier)
       call        parget('ga7cente',gcenter(7),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

       call        parget('ga8inten',gampli(8) ,iadda,ier)
       call        parget('ga8width',gwidth(8) ,iadda,ier)
       call        parget('ga8cente',gcenter(8),iadda,ier)
       if(ier.ne.0) goto 100
       if(ng.lt.mgaussians) then
         ng = ng + 1
       else
         write(6,*)'too many gaussians '
         goto 100
       endif
              

 100   continue

        sum = 0 
        do i=1,ng
    
         Omega     = o0 - gcenter(i)
         str_delta = abs(gwidth(i))

         a = 0
         b = 5.0d0/str_delta
         result  = adapint(strex_kernel,a,b,epsilon,maxit,erraccu)*2

         sum = sum + gampli(i)*result/(2*Pi)*sqrt(Pi)         

!         write(6,'(f12.6,i2,4f12.6)')
!     *   o0,i,gcenter(i),gwidth(i),gampli(i),sum

        enddo

       dwf  = exp(-u_sqr*qz*qz/3.0d0)
 
       th37 = a0*dwf*sum + bkgr
c
       return
       end




       function strex_kernel(t)
!      ------------------------
       implicit none

       real*8 strex_kernel, t

       real*8 Omega, str_beta, str_tau0, str_delta
       common /cstrex1/Omega, str_beta, str_tau0, str_delta

       strex_kernel= 
     *  exp(-(t/str_tau0)**str_beta) * 
     *  exp(-1d0/4d0*(str_delta*t)**2) * cos(t*Omega)*str_delta    

       return
       end



