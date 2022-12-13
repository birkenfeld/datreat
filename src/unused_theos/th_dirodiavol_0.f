      function th_ddiavo(x,pa,thnam,parnam,npar,idum,ini)
c     ===================================================
!
!      rotationl diffusion and translational diffusion
!      see quasielastic models and paper: PHYSICAL REVIEW A VOLUME 38, NUMBER 5 SEPTEMBER 1, 1988
!      Effect of rotational diffusion on quasielastic light scattering from fractal colloid aggregates
!      H. M. Lindsay R. Klein et al
!
!      There is one zone of defined undrgoing rotational diffusion
!      This zone may contain groups that diffuse according to dianoux volino inside a sphere
!      if fraction is set to zero only dionoux volino is computed


       implicit none
       
       character*8 thnam,parnam(20)
   
       real*4    x, pa, qq, th_ddiavo
       dimension pa(20),qq(3)
       integer   npar, ini, nparx,idum

       integer mgaussians
       parameter(mgaussians=10)
       double precision Pi
       parameter(Pi=3.141592654d0)
 
       

       real*8    tau, tau0, beta, a0, epsilon
       real*8    qz
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
       real*8 strex_kernel_dd
       external strex_kernel_dd

! communication with Kernel
 

       real*8 Omega, xwidth, str_delta, qc, 
     *       diffcm, rotdiff, 
     *       r1i,r1o,frac1, r2i,r2o,frac2, r3i,r3o,frac3
       integer lmax
       common /cstkd0/Omega, xwidth, str_delta, qc,
     *       diffcm, rotdiff, 
     *       r1i,r1o,frac1, r2i,r2o,frac2, r3i,r3o,frac3, lmax

       double precision rdivo, ddivo, zmxdivo
       integer          lmxdivo
       common /cstdivo/rdivo, ddivo, zmxdivo, lmxdivo



c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'ddiavo'
         nparx = 15
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_ddiavo = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          ! prefactor
         parnam(2) = 'diffcm'            ! diffusion constant  
         parnam(3) = 'rotdiff'           ! rotational diffusion constant
         parnam(4) = 'r1i'               ! first shell inner radius
         parnam(5) = 'r1o'               ! frist shell outer radius
         parnam(6) = 'frac1'             ! proton fraction in first shell
         parnam(7) = 'rdivo'             ! radius of diffuion inside sphere
         parnam(8) = 'ddivo'             ! diffusion constant inside sphere
         parnam(9) = 'lmxdivo'           ! lmax for dianous volino
         parnam(10)= 'zmxdivo'           ! zmax (for zeroes in Dianoux Volino

         parnam(11) = 'epsilon '          ! accuracy parameter for FT-integrations (DO NOT FIT)
         parnam(12) = 'omega0'            ! omega scale zero shift
         parnam(13) = 'u_sqr'             ! <u**2> value for Debye-Waller-Factor
         parnam(14) = 'xwidth'            ! channel integration witdth (see gauss2)
         parnam(15) = 'lmax'              ! max l-summation limit for rotdiff

c
         th_ddiavo = 0
         return
       endif
c
c ---- calculate theory here -----
       o0           = x   -   pa(12)
       a0           = pa(1)
       diffcm       = abs(pa(2))
       rotdiff      = abs(pa(3))
       r1i          = abs(pa(4))   
       r1o          = abs(pa(5))   
       frac1        = abs(pa(6))   
       rdivo        = abs(pa(7))    
       ddivo        = abs(pa(8))   
       lmxdivo      = nint(pa(9))   
       zmxdivo      = abs(pa(10))   

       epsilon      = abs(pa(11))
       u_sqr        = abs(pa(13))
       xwidth       = pa(14)
       lmax         = NINT(pa(15))

       if(epsilon.eq.0.0d0) epsilon = 1.0d-8
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
         b = 5.0d0/str_delta
         result  = adapint(strex_kernel_dd,a,b,epsilon,maxit,erraccu)*2

         sum = sum + gampli(i)*result/(2*Pi)*sqrt(Pi)         

!         write(6,'(f12.6,i2,4f12.6)')
!     *   o0,i,gcenter(i),gwidth(i),gampli(i),sum

        enddo

       dwf  = exp(-u_sqr*qz*qz/3.0d0)

       th_ddiavo = a0*dwf*sum

c
       return
       end




       function strex_kernel_dd(t)
!      ---------------------------
       implicit none

       real*8 strex_kernel_dd, t
       double precision dinsph_t


       real*8 Omega, xwidth, str_delta, qc,
     *       diffcm, rotdiff, 
     *       r1i,r1o,frac1, r2i,r2o,frac2, r3i,r3o,frac3
       integer lmax
       common /cstkd0/Omega, xwidth, str_delta, qc,
     *       diffcm, rotdiff, 
     *       r1i,r1o,frac1, r2i,r2o,frac2, r3i,r3o,frac3, lmax


       double precision rdivo, ddivo, zmxdivo
       integer          lmxdivo
       common /cstdivo/rdivo, ddivo, zmxdivo, lmxdivo


       double precision  sumtf, rri(3), rro(3), fr(3)
       double precision  radial_bsjn2_integral
       integer           i, l

!       strex_kernel_k0= 
!     *  exp(-(t/str_tau0)**str_beta) * 
!     *  exp(-1d0/4d0*(str_delta*t)**2) * cos(t*Omega)*str_delta    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! computing the time function                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        sumtf = 0

        if(frac1.gt.0d0) then
!         ---> rotational diffusion of protons in a zone from r1i to r1o  <----    
        do l=0,lmax
           sumtf=sumtf + 
     *             (2*l+1)*radial_bsjn2_integral(l,qc,r1i,r1o)*
     *             exp(-l*(l+1)*rotdiff*t)
        enddo
        endif

        sumtf = (1d0-frac1)+frac1*sumtf
!       -------------------------------
!
        sumtf = sumtf * dinsph_t(qc,t,ddivo,rdivo,lmxdivo,zmxdivo)
!       ---------------------------------------------------------- ! Multiply with Dianoux Volino
!            assuming that local bound protons can diffuse inside a limited
!            'sphere' of radisu rdivo
        

       
        sumtf = sumtf * exp(-diffcm*qc*qc*t)


        strex_kernel_dd= 
     *  sumtf * 
     *  exp(-1d0/4d0*(str_delta*t)**2) * 
     *  (sin(-t*Omega+0.5d0*t*xwidth)+sin(t*Omega+0.5d0*t*xwidth))/    !! this replaces cos(t*Omega)
     *  (t*xwidth) * str_delta                                         !! in order to yield the integral
                                                                       !! over one channel width in omega


       return
       end



