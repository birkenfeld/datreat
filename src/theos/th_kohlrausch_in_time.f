      function th_kohl_q(x,pa,thnam,parnam,npar,idum,ini)
c     ===================================================
!
!      kohlrausch in time
!      kww-function with DWF and EISF
!      Perez Aparcio, Arbe, Colmenro, Macromeleules 2006, 1060
!


       implicit none
       
       character*8 thnam,parnam(20)
   
       real*4    x, pa, qq, th_kohl_q
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
       real*8 strex_kernel
       external strex_kernel

! communication with Kernel
 

       real*8 Omega, str_beta, str_tau0, str_delta
       common /cstrex1/Omega, str_beta, str_tau0, str_delta



c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'kohl_q'
         nparx = 12
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_kohl_q = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          ! prefactor
         parnam(2) = 'tau0'              ! KWW time-constant, prefactor in front of q-dependence
         parnam(3) = 'beta'              ! streched exp, prefactor in front of-q-dependence
         parnam(4) = 'epsilon '          ! accuracy parameter for FT-integrations (DO NOT FIT)
         parnam(5) = 'omega0'            ! omega scale zero shift
         parnam(6) = 'u_sqr'             ! < u^2> value for Debye-Waller-Factor
         parnam(7) = 'j0'                ! jump length (if applicable)
         parnam(8) = 'beta0'             ! beta offset
         parnam(9) = 'qexp_t0'           ! q-exponent for tau0
         parnam(10)= 'qexp_bet'          ! beta-exponent 
         parnam(11)= 'bkgr'              ! constant background 
	 parnam(12)= 'n_mg'              ! fraction of hydrogen in side chains


c
         th_kohl_q = 0
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
       jlen         = pa(7)
       beta0        = pa(8)
       qexpt0       = pa(9)
       qexpbeta     = pa(10)
       bkgr         = pa(11)
       nmg          = pa(12)

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


! Apply the q-scalings
       if(qz.gt.0.0d0) then
         str_tau0 = str_tau0 * (1+jlen*qz**qexpt0)/(qz**qexpt0)
         str_beta = str_beta * (qz**qexpbeta + beta0)
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

       eisf = (1.0/3.0)*(1.0+2.0*(sin(qz*0.178)/qz*1.78))
 
       th_kohl_q = a0*dwf*(1-nmg+(nmg*eisf))*sum + bkgr

       th_kohl_q =  th_kohl_q + bgr_level + bgr_slope*o0
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



