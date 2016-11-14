
      function th_kohlfas0(x,pa,thnam,parnam,npar,idum,ini)
c     ===================================================
!
!     polymer chain with ch3-groups: description accordinmg to 
!     paper:R. 
!     Perez Aparicio,A. Arbe,J. Colmenero,B. Frick,L. Willner,D. Richter,L. J. Fetters
!     Macromolecules, 2006, 39,1060
!  
!     Further modifed by adding a fast local motion to the model with an on eisf2 
!
c
       use theory_description

       implicit none
       
       character*8 thnam,parnam(20)
   
       real*4    x, pa, qq, th_kohlfas0
       dimension pa(20),qq(3)
       integer   npar, ini, nparx, idum

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
 !      integer   :: idesc

       real*4    kbolz
       Parameter (kbolz=1.380662e-23)
 
       integer iadda
       common/thiadd/iadda


       real*8 omega0, u_sqr, dwf
       real*8 qcenter, qband, qexpt0, qexpbeta, bkgr 
       real*8 jlen, beta0       

       real*8 a,b, domega, o0
       real*8 adapint, sum, result, result2, err, erraccu
!       real*8 kch3_kernel
!       external kch3_kernel

! communication with Kernel
       double precision  ::  q, nch, nmg, rhh, betakww, taukww, uu
       double precision  ::  betafast, taufast, eisf2
       double precision  ::  qexpkww, qexpfast, taukww0, taufast0
       double precision  ::  rg_eisf
!       common /ckohlch3/ q, nch, nmg, rhh, betakww, taukww, uu,
!     *                   betafast, taufast, eisf2


       real*8 Omega, str_beta, str_tau0, str_delta, xwidth
!       common /cstrex1/Omega, str_beta, str_tau0, str_delta, xwidth

       real*8 ln_tau0, ln_beta, ln_width, cl_t
!       common /crlognor/ ln_tau0, ln_beta, ln_width, cl_t


c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'kohlfas0'
         nparx = 18
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th_kohlfas0 = 0
           return
         endif

        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = 
     *  "kohlrausch function augmented with ch3 motion:         "//cr//
     *  "  and an aditional fast component                      "
        th_citation(idesc)     = 
     *  "Perez Aparicio,A. Arbe,J. Colmenero,B. Frick,L. Willner"//
     *  ",D. Richter,L. J. Fetters"//cr//"Macromolecules, 2006, 39,1060"
 


         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          ! prefactor (should be 1)
         parnam(2) = 'n_chain'           ! chain protons fraction
         parnam(3) = 'n_ch3'             ! ch3-protons fraction
         parnam(4) = 'tau_kww'           ! main-chain tau
         parnam(5) = 'beta_kww'          ! main-chain beta
         parnam(6) = 'lntau_mg'          ! log of center of methyl tau's (here we enter ln(tau0) !! )
         parnam(7) = 'sigma_mg'          ! width of lognormal distribution (sigma of the paper)
         parnam(8) = 'beta_mg'           ! beta of mg rotation (paper ==> 1)
         parnam(9) = 'rhh'               ! radius of ch3-rotation (paper ==> 1.78 Angstroem)
         parnam(10)= 'u_square'          ! <u**2> value for Debye-Waller-Factor
         parnam(11)= 'omega0'            ! omega scale zero shift
         parnam(12)= 'epsilon'           ! accuracy parameter for FT-integrations (DO NOT FIT)
         parnam(13)= 'xwidth'            ! channelwidth
         parnam(14)= 'rgeisf2'           ! "eisf" for the fast local motion as exp(-(1/3)(Rg*q)**2)
         parnam(15)= 'tau_fast'          ! tau of fast component
         parnam(16)= 'beta_fast'         ! beta of fast 
         parnam(17)= 'qexpkww'           ! q exponent of tau kww
         parnam(18)= 'qexpfast'          ! q exponent of fast component
 
!
        th_param_desc(1,idesc) = "prefactor (default 1) " 
        th_param_desc(2,idesc) = "fraction of (non-ch3) chain protons " 
        th_param_desc(3,idesc) = "fraction of -ch3 protons " 
        th_param_desc(4,idesc) = "tau_kww (alpha relaxation) "//
     *                            "time in units 1/units(xaxis) " 
        th_param_desc(5,idesc) = "streching exponent of the KWW fkt." 
        th_param_desc(6,idesc) = "log of center of methyl taus" 
        th_param_desc(7,idesc) = "log normal width methyl taus" 
        th_param_desc(8,idesc) = "stretching exp of methyl fkt." 
        th_param_desc(9,idesc) = "radius of ch3-rotation (==> 1.78A)" 
        th_param_desc(10,idesc) = "u**2 value for Debye-Waller-Factor" 
        th_param_desc(11,idesc) = "omega0 zero adjustment xaxis" 
        th_param_desc(12,idesc) = "accuracy FT-integ. (DO NOT FIT)" 
        th_param_desc(13,idesc) = "channel width, set 0 if _xwidth "//
     *   "parameter is present in the data records"//cr//parspace// 
     *   "   takes care of proper integration of sharp features over"//
     *   " the experimental box width"
        th_param_desc(14,idesc) = "eisf2= exp(-(1/3)( Rgeisf2 *q)**2) "       
        th_param_desc(15,idesc) = "characteristic time of fast comp. "       
        th_param_desc(16,idesc) = "stretching exponent of fast comp. "       
        th_param_desc(17,idesc) = "exponent of tau_KWW dispersion  "       
        th_param_desc(18,idesc) = "exponent of tau_fast dispersion  "       

! 
        th_file_param(:,idesc) = " "  
        th_file_param(1,idesc) = 
     *  " q     =  q value of the S(q,t), units THE LENGTH UNITS"
        th_file_param(2,idesc) = 
     *  " temp  =  T value in K (currently not used)"
       th_file_param(3,idesc) = 
     *  " _xwidth  =  channel width (x-channel)"
       th_file_param(4,idesc) = 
     *  " ga#inten, ga#witdh, ga#cente   #=1..8"//cr//parspace//
     *  " Gaussian resolution parameter NEEDED at least #=1"

         th_out_param(:,idesc)  = " "    

c
         th_kohlfas0 = 0
         return
       endif
c
c ---- calculate theory here -----
       o0           = x   -   pa(11)

       a0           = pa(1)
       nch          = pa(2)
       nmg          = pa(3)
       taukww0      = abs(pa(4))
       betakww      = abs(pa(5))
       ln_tau0      = pa(6)
       ln_width     = pa(7)
       ln_beta      = abs(pa(8))
       rhh          = pa(9)
       uu           = pa(10)

       epsilon      = abs(pa(12))
       xwidth       = abs(pa(13))

       rg_eisf      = pa(14)
       taufast0     = abs(pa(15))
       betafast     = pa(16)
       qexpkww      = pa(17)
       qexpfast     = pa(18)
       
      

       if(xwidth.eq.0d0) xwidth =1d-6


     
       if(epsilon.eq.0.0d0) epsilon = 1.0d-7
       maxit = 10000

        qget = 0.0
        call        parget('q       ',qget,iadda,ier)
        qz   = qget
        q    = qget
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


       taukww   = taukww0  * q**(-qexpkww)
       taufast  = taufast0 * q**(-qexpfast)
       eisf2    = exp(-0.33333d0 *(rg_eisf * q)**2)    
 
      ! write(6,*)"taukww   ",taukww,  qexpkww
      ! write(6,*)"taufast  ",taufast, qexpfast
      ! write(6,*)"eisf     ",eisf2, rg_eisf, q




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
         b = 9.0d0/str_delta
         result  = adapint(kch3_kernel,a,b,epsilon,maxit,erraccu)*2

         sum = sum + gampli(i)*result/(2*Pi)*sqrt(Pi)         

!         write(6,'(f12.6,i2,4f12.6)')
!     *   o0,i,gcenter(i),gwidth(i),gampli(i),sum

        enddo

 
       th_kohlfas0 = a0*sum 
c
       return
!       end





       contains


       function kch3_kernel(t)
!      -----------------------
       implicit none

       double precision :: kch3_kernel
       real*8 t
      

!       real*8 Omega, str_beta, str_tau0, str_delta, xwidth
!       common /cstrex1/Omega, str_beta, str_tau0, str_delta, xwidth

       kch3_kernel= 
     *  ftsq(t) * 
     *  exp(-1d0/4d0*(str_delta*t)**2) *    
     *  (sin(-t*Omega+0.5d0*t*xwidth)+sin(t*Omega+0.5d0*t*xwidth))/    !! this replaces cos(t*Omega)
     *  (t*xwidth) * str_delta                                         !! in order to yield the integral

       return
       end


       function ftsq(t)
!      ---------------- time domain scattering function
       implicit none
       double precision :: ftsq       


!       double precision  ::  q, nch, nmg, rhh, betakww, taukww, uu
!       double precision  ::  betafast, taufast, eisf2
!       common /ckohlch3/ q, nch, nmg, rhh, betakww, taukww, uu,
!     *                   betafast, taufast, eisf2

       double precision ::  t

       double precision :: eisf, skww_inc, s_rot, dwf
       double precision :: skww_fast

       eisf       = (1d0+2d0*sin(q*rhh)/(q*rhh))/3d0
       dwf        = exp(-uu*q*q/3d0)

       skww_inc   = exp(-(t/taukww)**betakww)
       skww_fast  = exp(-(t/taufast)**betafast)
       s_rot      = eisf+(1-eisf)*f_rlognor(t)

       ftsq = dwf*(nch+nmg*s_rot)*skww_inc*(eisf2+(1-eisf2)*skww_fast)


       return
       end
       

       function f_rlognor(t)
!      --------------------
       implicit none

       real*8 f_rlognor, t

!       real*8 Omega, str_beta, str_tau0, str_delta, xwidth
!       common /cstrex1/Omega, str_beta, str_tau0, str_delta, xwidth

!       real*8 ln_tau0, ln_beta, ln_width, cl_t
!       common /crlognor/ ln_tau0, ln_beta, ln_width, cl_t

       real*8 a,b
       real*8 a2dapint, result, err, erraccu, epsilon

       integer maxit

!      real*8   rlnor_kernel
!      external rlnor_kernel


       maxit   =  5000
       epsilon =  1d-5


       cl_t = t
       a    = ln_tau0-4*ln_width
       b    = ln_tau0+4*ln_width
        
       f_rlognor  = a2dapint(rlnor_kernel,a,b,epsilon,maxit,erraccu)

       return
       end



       function rlnor_kernel(lt)
!      -----------------------
       implicit none

       double precision Pi
       parameter(Pi=3.141592654d0)


       real*8 rlnor_kernel, lt

!      real*8 Omega, str_beta, str_tau0, str_delta, xwidth
!       common /cstrex1/Omega, str_beta, str_tau0, str_delta, xwidth

!       real*8 ln_tau0, ln_beta, ln_width, cl_t
!       common /crlognor/ ln_tau0, ln_beta, ln_width, cl_t


       real*8 arg1, arg2, arg


       arg1 = -(((lt-ln_tau0)/ln_width)**2)/2d0
       arg2 = -(cl_t/exp(lt))**ln_beta
 
       arg  = arg1+arg2

       if(arg.lt.-50d0) arg = -50d0
    

       rlnor_kernel = exp(arg)/(sqrt(2d0*Pi)*ln_width)


      
       return
       end

       end function th_kohlfas0
