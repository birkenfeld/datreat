      function th_grotmod2(x,pa,thnam,parnam,npar,idum,ini)
c     ===================================================
!
!     model to describe the combined effect of rotational diffusion
!     1D local motion  (domain)
!     and a 3D Gaussian confined motion  
!     
!     
!     
!
!    
!   
!    
       implicit none
       
       character*8 thnam,parnam(20)
   
       real*4    x, pa, qq, th_grotmod2
       dimension pa(20),qq(3)
       integer   npar, ini, nparx,idum

       integer mgaussians
       parameter(mgaussians=10)
       double precision Pi
       parameter(Pi=3.141592654d0)
 
       

       real*8    tau, tau0, beta, a0, epsilon
       real*8    qz, dummy
       real*4    tget, temp, qget

       real*4    bgr_level, bgr_slope

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


       real*8 omega0, u_sqr, dwf, eisf
       real*8 qcenter, qband, qexpt0, qexpbeta, bkgr, nmg 
       real*8 jlen, beta0

       real*8 a,b, domega, o0
       real*8 adapint, sum, result, result2, err, erraccu
       real*8 strex_kernel_gm2
       external strex_kernel_gm2

! communication with Kernel
 

       real*8 Omega, xwidth, str_delta, qc, 
     *       diffcm, rotdiff,r1i,r1o,amp1d,tau1d,rconfine,diffconf 
     *       , frac_ch3, tau_ch3, beta_ch3
       integer lmax
       common /cstgrot2/Omega, xwidth, str_delta, qc,
     *       diffcm, rotdiff,r1i,r1o,amp1d,tau1d,rconfine,diffconf  
     *       , frac_ch3, tau_ch3, beta_ch3




c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'grotmod2'
         nparx = 16
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_grotmod2 = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          ! prefactor
         parnam(2) = 'diffcm'            ! diffusion constant  
         parnam(3) = 'rotdiff'           ! rotational diffusion constant
         parnam(4) = 'r1i'               ! first shell inner radius
         parnam(5) = 'r1o'               ! frist shell outer radius
         parnam(6) = 'amp1d'             ! amplitude of 1D motion
         parnam(7) = 'tau1d'             ! time constant fo 1D motion
         parnam(8) = 'rconfine'          ! Gaussian confinement radius
         parnam(9) = 'diffconf'          ! effective Diffusion in confined Volume
         parnam(10)= 'frac_ch3'          ! fraction of CH3 groups
         parnam(11)= 'tau_ch3'           ! tau of CH3 groups relaxation
         parnam(12)= 'beta_ch3'          ! stretching exponentof CH3 groups relaxation
        


         parnam(13) = 'epsilon '          ! accuracy parameter for FT-integrations (DO NOT FIT)
         parnam(14) = 'omega0'            ! omega scale zero shift
         parnam(15) = 'u_sqr'             ! <u**2> value for Debye-Waller-Factor
         parnam(16) = 'xwidth'            ! channel integration witdth (see gauss2)

c
         th_grotmod2 = 0
         return
       endif
c
c ---- calculate theory here -----
       a0           = pa(1)
       diffcm       = abs(pa(2))
       rotdiff      = abs(pa(3))
       r1i          = abs(pa(4))   
       r1o          = abs(pa(5))   
       amp1d        = abs(pa(6))   
       tau1d        = abs(pa(7))    
       rconfine     = abs(pa(8))   
       diffconf     = abs(pa(9))   
       frac_ch3     = abs(pa(10))   
       tau_ch3      = abs(pa(11))   
       beta_ch3     = abs(pa(12))   

       epsilon      = abs(pa(13))
       o0   = x   -   pa(14)
       u_sqr        = abs(pa(15))
       xwidth       = pa(16)

       if(epsilon.eq.0.0d0) epsilon = 1.0d-10
       maxit = 10000

        qget = 0.0
        call        parget('q       ',qget,iadda,ier)
        qz   = qget
        if(ier.ne.0) write(6,*)'Warning q not found' 

       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif


       if(xwidth.eq.0.0d0) then
         call        parget('_xwidth ',tget,iadda,ier)
         xwidth = tget
       endif

       if(xwidth.eq.0.0d0) then
         write(6,*)'ATTENTION XWIDTH = 0 set 0.1 to avoud crash'
         xwidth = 0.1
       endif

!



! extract gaussian parameters 
!
! assuming Gaussian parameters according to the form of th:gauss
!
!      g = ampli*exp(-((x-center)/width)**2)
!      -------------------------------------

       bgr_level = 0
       bgr_slope = 0
 
       call        parget('bk1level',bgr_level ,iadda,ier)
       call        parget('bk1slope',bgr_slope ,iadda,ier)

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

        qc   = qz
        sum  = 0 
        do i=1,ng
    
         Omega     = o0 - gcenter(i)
         str_delta = abs(gwidth(i))

         a = 0
         b = 9.0d0/str_delta
         result  = adapint(strex_kernel_gm2,a,b,epsilon,maxit,erraccu)*2

         sum = sum + gampli(i)*result/(2*Pi)*sqrt(Pi)         

!         write(6,'(f12.6,i2,4f12.6)')
!     *   o0,i,gcenter(i),gwidth(i),gampli(i),sum

        enddo

       dwf  = exp(-u_sqr*qz*qz/3.0d0)

       th_grotmod2 = a0*dwf*sum

c
       return
       end




       function strex_kernel_gm2(t)
!      ---------------------------
       implicit none

       real*8 strex_kernel_gm2, t

       double precision :: u1d, u2d, u3d, uuz, uux, sq, ddaws, sumtf
       double precision :: eisf, sch3
       double precision, parameter :: rhh=1.2d0
       double precision, parameter :: Pi=3.141592654d0
  

       real*8 Omega, xwidth, str_delta, qc, 
     *       diffcm, rotdiff,r1i,r1o,amp1d,tau1d,rconfine,diffconf 
     *       , frac_ch3, tau_ch3, beta_ch3

       integer lmax
       common /cstgrot2/Omega, xwidth, str_delta, qc,
     *       diffcm, rotdiff,r1i,r1o,amp1d,tau1d,rconfine,diffconf  
     *       , frac_ch3, tau_ch3, beta_ch3




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! computing the time function                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       u1d = amp1d**2*(1-exp(-t/tau1d))
       u2d = rotdiff*(r1i/2+r1o/2)**2*t
       u3d = (rconfine**2)/5*(1-exp(-t*diffconf/((rconfine**2)/5)))

!!! hier die Dimesionsfaktoren ueberpruefen!!!!! !!!! 

       uuz = u1d+u3d
       uux = u2d+u3d
     

       if(uuz.gt.uux) then
         sq = qc*sqrt(uuz-uux)
         if(sq.lt.1d-4) sq=1d-4
!         write(6,*)'sq=',sq
!         th_unrou_inc = a0*S15AEF(sq,ifail)*sqrt(Pi)*exp(-q*q*uux/2d0)/(2*sq)
          sumtf = derf(sq)*sqrt(Pi)*exp(-qc*qc*uux)/(2*sq)
!                         Erf
       else
         sq = qc*sqrt(-uuz+uux)
         if(sq.lt.1d-4) sq=1d-4
!         th_unrou_inc = a0* S15AFF(sq,ifail)*2*exp(sq**2-q*q*uux/2d0)/(2*sq)
         sumtf = ddaws(sq)*2*exp(sq**2-qc*qc*uux)/(2*sq)
!                          Dawson, wegen erf(i*x)=2*i*exp(x**2)*Dawson(x)/sqrt(Pi)
       endif

       
        sumtf = sumtf * exp(-diffcm*qc*qc*t)

!      include the effect of ch3-groups
        if(frac_ch3.gt.0d0) then
          eisf  = (1d0+2d0*sin(qc*rhh)/(qc*rhh))/3d0
        
          sch3  =(eisf+(1-eisf)*exp(-(t/tau_ch3)**beta_ch3))
          sumtf = sumtf * ((1-frac_ch3)+frac_ch3*sch3)
        endif


        strex_kernel_gm2= 
     *  sumtf * 
     *  exp(-1d0/4d0*(str_delta*t)**2) * 
     *  (sin(-t*Omega+0.5d0*t*xwidth)+sin(t*Omega+0.5d0*t*xwidth))/    !! this replaces cos(t*Omega)
     *  (t*xwidth) * str_delta                                         !! in order to yield the integral
                                                                       !! over one channel width in omega


       return
       end



