 FUNCTION th_rpalir3(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!============================================================================
!  applies rpa mixture tor the scattering function for a ring-linear(log) polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))
! 
      use theory_description 
      use lmfit
      use rpa_laplace
      use polynom
      use reptation3       ! the module for reptation we use in reptest6 as use in the JCP üaüer
      implicit none
!! 
      SAVE
!!    ====


      real    :: th_rpalir3
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
 
!      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address

 
! the internal parameter representation 
     double precision :: ampli          ! prefactor

!! ring diffusion related:
     double precision :: ring_diff       ! limiting diffusion (NMR value) in [A**2/ns]                                     
     double precision :: ring_r021       ! reference mean squared displacement at first transition point to medium timed di
     double precision :: ring_alpha      ! sub-diffusion  exponent of SHORT time diffusion                                 
     double precision :: ring_r022       ! reference mean squared displacement at second transition point to D0 in [A**2]  
     double precision :: ring_beta       ! sub-diffusion  exponent of medium time diffusion                                
     double precision :: ring_a_cross    ! transition exponent between short and long time diffusion (sharper kink for larg
!!                                     
!!     double precision :: l          ! effective segment length     == aring                                                     
     double precision :: ring_nue        ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
     double precision :: ring_nue0       ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
     double precision :: ring_ne0        !                      
     double precision :: ring_wl4        ! Rouse rate in [A**4/ns]                                                         
     double precision :: ring_pmin       ! transition mode number between simple ring-Rouse and large p modification       
     double precision :: ring_pwidth     ! sharpness of transition     

     double precision :: ring_ntran      ! structure transitions from nu0 to nue
     double precision :: ring_nwid       ! structure transitions width

                                                    
     double precision :: ring_f0         ! prefactor f(p) limit for small p values (default 1)                             
     double precision :: ring_finf       ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
     double precision :: ring_pexinf     ! large p tau(p) = tauinf/p**pexinf  

!     double precision :: diffrous   ! centre-of-mass diffusino star ignore diffusion of the lonh linear chains !?!
!?!     double precision :: diff       
!?!     double precision :: r02   , rr     
!?!     double precision :: tc         
!?!     double precision :: nu_subdiff 
!?!     double precision :: a_cross    
!?!     double precision :: beta_q     
     double precision :: diff,alpha,beta,r021,r022,a_cross


     double precision :: rep_l
     double precision :: rep_ne
     double precision :: rep_wl4
     double precision :: rep_re
     double precision :: rep_re_fac
     double precision :: rep_alpha0
     double precision :: rep_talphamax
     double precision :: rep_talphawd
     integer          :: rep_clf
  

! the recin parameter representation 
     double precision :: q          ! scattering wavevector
     double precision :: temp       ! temperature
     integer          :: nring      ! number of segments in ring polymer
     double precision :: lring      ! segment length ring 
     double precision :: llin       ! segment length linear
     double precision :: phiring    ! volume fraction of star
     integer          :: nlin       ! number of segments of linear polymer
     integer          :: nring_cc   ! number of segments of ring polymer matrix 
     double precision :: philin     ! volume fraction of linear polymer

     double precision :: blin       ! scattering contrast of linear polymer omponent 
     double precision :: bring      ! scattering contrast of start polymer

     integer          :: nxpoints   ! number of points in sum(exp) modelling of S00 functions 

     integer          :: mode       ! mode select 
     integer          :: modeex     ! mode select exp representation
     integer          :: modeex1      ! mode select exp representation
     integer          :: modeex2      ! mode select exp representation
     integer          :: modeexcc     ! mode select exp representation


     double precision :: Sqt_result(2)
     double precision :: tau_e


     double precision, parameter :: tmin = 0.001d0
     double precision            :: t0  ! virtual 0 > 0 to allow the same numerics as for finite t to compute s(q,0)



     double precision :: tau        ! for modes > 0
 
     double precision   :: t

     double precision   :: sqt, sqt0, fqt
     double precision   :: pring0, pring
     double precision   :: plin0, plin
!     double precision   :: plin0cc, plincc
     double precision   :: pring0cc, pringcc
!?!     double precision   :: dr, wx, lx
     double precision   :: Re_ring
     double precision   :: Rgr
!     double precision   :: Sinv0, Sinv
     double precision   :: repta_fqt

     double precision, parameter  :: t_table_spacing1 = 300d0 
     double precision, parameter  :: t_table_spacing2 = 100d0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     integer, parameter      :: mparams=60
     double precision        :: last_params(mparams) = 0d0
     double precision        :: allparams(mparams)   = 0d0
     integer, parameter      :: mexp=10
     double precision        :: aexp11(mexp), rexp11(mexp)
     integer                 :: nex11
     double precision        :: aexp22(mexp), rexp22(mexp)
     integer                 :: nex12
     double precision        :: aexpcc(mexp), rexpcc(mexp)
     integer                 :: nexcc

     integer, parameter      :: mxpoints = 10000
     double precision        :: t_samples(mxpoints)
     double precision        :: s_samples(mxpoints)
     double precision        :: exp_samples(mxpoints)
     double precision        :: tmax = 1000d0
     logical                 :: newcomp_required
     integer                 :: i
     double precision        :: ts, rmsdev, rmsdev_limit = 9.9d-3

     double precision        :: ss11, ss110, ss12, ss120, ss22, ss220
     double precision        :: sqtres(2)

     double precision        :: a1, a2, r1, r2, r3, b1, b2, g1, g2, g3
     integer                 :: iout = 0


     complex(kind=XPREC)     :: il_coeffs11(3*mexp)
     complex(kind=XPREC)     :: il_coeffs12(3*mexp)
     complex(kind=XPREC)     :: il_coeffs22(3*mexp)
     complex(kind=XPREC)     :: il_alphas(3*mexp)
     integer                 :: nnsum
     integer                 :: j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ----- initialisation ----- 
    IF (ini.eq.0) then 
       thnam = 'rpalir3 '
       nparx =       32
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpalir3 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " applies danamic rpa to the scattering function for a "//cr//parspace//&
                                " ring-lin(long)  polymer mixture   "//cr//parspace//&
                                " with ringmodel ring2am7 (without com diffusion) "//cr//parspace//&
                                " and long linear according to reptest6 " //cr//parspace//&
                                " TEST WITH CC=RING instead CC=lin "

       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli'     ! prefactor

        parnam ( 2) = 'rep_ne'
        parnam ( 3) = 'rep_wl4'
        parnam ( 4) = 'rep_f'      ! refac
        parnam ( 5) = 'rep_alp0'
        parnam ( 6) = 'rep_tamx'
        parnam ( 7) = 'rep_tawd'
        parnam ( 8) = 'rep_clf'

        parnam ( 9) = 'diff    '  ! limiting diffusion (NMR value) in [A**2/ns]                                     
        parnam (10) = 'r021    '  ! reference mean squared displacement at first transition point to medium timed di
        parnam (11) = 'alpha   '  ! sub-diffusion  exponent of SHORT time diffusion                                 
        parnam (12) = 'r022    '  ! reference mean squared displacement at second transition point to D0 in [A**2]  
        parnam (13) = 'beta    '  ! sub-diffusion  exponent of medium time diffusion                                
        parnam (14) = 'a_cross '  ! transition exponent between short and long time diffusion (sharper kink for larg
                    
        parnam (15) = 'ri_nue  '  ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
        parnam (16) = 'ri_wl4  '  ! Rouse rate in [A**4/ns]                                                         
        parnam (17) = 'ri_pmin '  ! transition mode number between simple ring-Rouse and large p modification       
        parnam (18) = 'ri_pwid '  ! sharpness of transition                                                         
        parnam (19) = 'ri_f0   '  ! prefactor f(p) limit for small p values (default 1)                             
        parnam (20) = 'ri_finf '  ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
        parnam (21) = 'ri_ne0  '  ! experimental ne0                                              
        parnam (22) = 'ri_nue0 '  ! nue0 (Gaussian=0.5)                                              
        parnam (23) = 'ri_ntran'  ! structure n transition                                              
        parnam (24) = 'ri_nwid '  ! width of that    
        parnam (25) = 'ri_pexin'  ! large p tau(p) = tauinf/p**pexinf     


        parnam (26) = 'mode    '  ! what to compute
        parnam (27) = 'modeex  '  ! selects way of representation in terms of nexp
        parnam (28) = 'nxpoints'  ! number of points to determine representation in terms of nexp
        parnam (29) = 'tzero   '  ! preliminary: effective zero time for S(Q,t=0=tzero)
        parnam (30) = 'tmaxrng '  ! maximum range that shall be spanned by creating the n-exp model
        parnam (31) = 'minrate '  ! minimum rate allowed in the exp models
        parnam (32) = 'rmsdev  '  ! n-exp-fit error limit (some 1e-3)
  
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of entanglement segements of linear long chain" !//cr//parspace//&
        th_param_desc( 3,idesc) = "rouse rate of linear long chain" !//cr//parspace//&
        th_param_desc( 4,idesc) = "pseudo anisotropy factor of tube dimensions of linear long chain" !//cr//parspace//&
        th_param_desc( 5,idesc) = "Non-Gaussianty amplitude (for blobs of long chain) " !//cr//parspace//&
        th_param_desc( 6,idesc) = "t-maximum of NG alpha(t) function" !//cr//parspace//&
        th_param_desc( 7,idesc) = "t-width (log) of NG alpha(t) function" !//cr//parspace//&
        th_param_desc( 8,idesc) = "switch (off/on) clf contribution" !//cr//parspace//&
        th_param_desc( 9,idesc) = "limiting diffusion (NMR value) in [A**2/ns]                                     "
        th_param_desc(10,idesc) = "reference mean squared displacement at first transition point to medium timed di"
        th_param_desc(11,idesc) = "sub-diffusion  exponent of SHORT time diffusion                                 "
        th_param_desc(12,idesc) = "reference mean squared displacement at second transition point to D0 in [A**2]  "
        th_param_desc(13,idesc) = "sub-diffusion  exponent of medium time diffusion                                "
        th_param_desc(14,idesc) = "transition exponent between short and long time diffusion (sharper kink for larg"

        th_param_desc(15,idesc) = "chain statistics exponent for rings"  !  
        th_param_desc(16,idesc) = "rouse rate of rings "  !   
        th_param_desc(17,idesc) = "transition mode number between simple ring-Rouse and large p modification "  !
        th_param_desc(18,idesc) = "width of the transition"  !    

        th_param_desc(19,idesc) = "prefactor f(p) limit for small p values (default 1)  " !//cr//parspace//&
        th_param_desc(20,idesc) = "prefactor f(p) limit for large p values (default F=0.9??) transition width" !//cr//parspace//&
        th_param_desc(21,idesc) = "experimental ne0   " !//cr//parspace//&
        th_param_desc(22,idesc) = "nue0 (Gaussian=0.5) " !//cr//parspace//&
        th_param_desc(23,idesc) = "structure n transition " !//cr//parspace//&
        th_param_desc(24,idesc) = "width of structure n transition " !//cr//parspace//&
        th_param_desc(25,idesc) = "large p tau(p) = tauinf/p**pexinf   " !//cr//parspace//&
       
        th_param_desc(26,idesc) = "mode: 0 normal, >0 S(x=Q,tau) 1 lin, 2 star 3 rpa" !//cr//parspace//&
        th_param_desc(27,idesc) = "modex: number of exps to represent the undistrurbed S(q,t) relaxation curves" !//cr//parspace//&
        th_param_desc(28,idesc) = "number of points to determine representation in terms of nexp" !//cr//parspace//&
        th_param_desc(29,idesc) = "preliminary: effective zero time for S(Q,t=0=tzero)" //cr//parspace//&
                                "for the numerical integration we have a smoothing of the 0-->1 jump " //cr//parspace//&
                                "at t=0 due to the pratically less than infinite width of the s-integration"  //cr//parspace//&
                                "for the analytic=1 method this may be zero, for the analytic=0 method something" //cr//parspace//&
                                "like 0.001 * tmax could be a good start to try"
        th_param_desc(30,idesc) = " maximum range that shall be spanned by creating the n-exp model" !//cr//parspace//&
        th_param_desc(31,idesc) = " minimum rate allowed in the exp models" !//cr//parspace//&
        th_param_desc(32,idesc) = " n-exp fit error-limit rmsdev_limit" !//cr//parspace//&

! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  3,idesc) = "nring    > number of segments per star arm"
        th_file_param(  4,idesc) = "lring    > segment length ring"
        th_file_param(  5,idesc) = "llin     > segment length lin"
        th_file_param(  6,idesc) = "phiring  > volume fraction of star"
        th_file_param(  7,idesc) = "nlin     > number of segments of linear polymer"
        th_file_param(  8,idesc) = "nring_cc > number of segments of linear polymer matrix"
        th_file_param(  9,idesc) = "philin   > volume fraction of linera polymer"
        th_file_param( 10,idesc) = "tau      > tau val for mode > 0 calc"
        th_file_param( 11,idesc) = "blin     > scattering contrast linear"
        th_file_param( 12,idesc) = "bring    > scattering contrast star"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpalir3 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
     ampli         =      pa( 1)
     rep_ne        = abs( pa( 2) ) 
     rep_wl4       = abs( pa( 3) ) 
     rep_re_fac    = abs( pa( 4) )    ! refac
     rep_alpha0    = abs( pa( 5) ) 
     rep_talphamax = abs( pa( 6) ) 
     rep_talphawd  = abs( pa( 7) ) 
     rep_clf       = nint(pa( 8) ) 

     diff          = abs(  pa( 9) )
     r021          = abs(  pa(10) )
     alpha         = abs(  pa(11) )
     r022          = abs(  pa(12) )
     beta          = abs(  pa(13) )
     a_cross       = abs(  pa(14) )
                                                    
     ring_nue      =  abs( pa(15) )   ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
     ring_wl4      =  abs( pa(16) )   ! Rouse rate in [A**4/ns]                                                         
     ring_pmin     =     ( pa(17) )   ! transition mode number between simple ring-Rouse and large p modification       
     ring_pwidth   =  abs( pa(18) )   ! sharpness of transition                                                         
     ring_f0       =  abs( pa(19) )   ! prefactor f(p) limit for small p values (default 1)                             
     ring_finf     =  abs( pa(20) )   ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
     ring_ne0      =  abs( pa(21) )   ! experimental ne0                                              
     ring_nue0     =  abs( pa(22) )   ! nue0 (Gaussian=0.5)                                              
     ring_ntran    =  abs( pa(23) )   ! structure n transition                                              
     ring_nwid     =  abs( pa(24) )   ! width of that    
     ring_pexinf   =  abs( pa(25) )   ! large p tau(p) = tauinf/p**pexinf   

      mode         =  nint(pa(26))
      modeex       =  nint(pa(27))
      nxpoints     =  nint(pa(28))
      t0           =  abs( pa(29))
      tmax         =  abs( pa(30))
      rlow         =  abs( pa(31))  ! parameter in rpa_laplace
      rmsdev_limit = max(1e-4,abs( pa(32)))  ! parameter in rpa_laplace

      maximum_rouse_n = 300         ! sets summation limit in module reptation3


      t0       = max(t0      ,tmin)
      nxpoints = max(nxpoints  ,16)
      modeex   = max(modeex    , 8)

      if( modeex > mexp ) then
         write(6,*)"INFORMATION: Modeex limited to: ",mexp
         modeex = mexp
      endif

      if( nxpoints > mxpoints ) then
         write(6,*)"INFORMATION: nxpoints limited to: ",mxpoints
         nxpoints = mxpoints
      endif

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: scattering wavevector
      xh = 0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh = 400d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: number of segments per star arm
      xh = 1
      call parget('nring    ',xh,iadda,ier)
      nring     = xh
! >>> extract: segment length
      xh = 5.68d0   !! the PEO value as defaul   !! the PEO value as defaultt
      call parget('lring    ',xh,iadda,ier)
      lring      = xh
! >>> extract: segment length
      xh = 5.68d0
      call parget('llin     ',xh,iadda,ier)
      llin        = xh
! >>> extract: volume fraction of rings
      xh = 0.01d0
      call parget('phiring ',xh,iadda,ier)
      phiring  = xh
! >>> extract: number of segments of linear polymer
      xh = 100
      call parget('nlin    ',xh,iadda,ier)
      nlin     = nint(xh)
!! >>> extract: number of segments of linear polymer
      xh = 100
      call parget('nring_cc ',xh,iadda,ier)
      nring_cc  = nint(xh)
! >>> extract: volume fraction of linear polymer
      xh = 0.9d0
      call parget('philin  ',xh,iadda,ier)
      philin   = xh
!! >>> extract: tau for mode > 0
      xh = 0.0d0
      call parget('tau     ',xh,iadda,ier)
      tau   = xh
!! >>> extract: contratst lin
      xh = 0.0d0
      call parget('blin     ',xh,iadda,ier)
      blin   = xh
!!! >>> extract: contrats ring
      xh = 1.0d0
      call parget('bring    ',xh,iadda,ier)
      bring  = xh
!!!!!!!! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
i1:  if( mode == 0 ) then     ! normal spin-echo
       t   = x
     else
       t   = tau
       q   = x
     endif i1

      allparams(1:42) = [ q, &
                          lring, llin ,  dble(nring), dble(nlin), dble(nring_cc), phiring, philin, blin, bring, &
                          rep_ne, rep_ne, rep_wl4, rep_re_fac, rep_alpha0,  rep_talphamax ,rep_talphawd, dble(rep_clf), &
                          ring_nue, ring_wl4,  ring_pmin, ring_pwidth, ring_f0, ring_finf, &
                          ring_ne0, ring_nue0, ring_ntran, ring_nwid, ring_pexinf, &
                          dble(mode), dble(modeex), dble(nxpoints), t0, tmax, rlow, rmsdev_limit, &
                          diff,alpha,beta,r021,r022,a_cross]  


      newcomp_required = ( sum(abs(allparams-last_params)) .ne. 0 ) .and. (modeex > 0)
    
!!?? testing if(analytic < 0) newcomp_required = .true.

!! testing only
! goto 3333
!!!!!!!!!!!!!!!



      if(newcomp_required) then
        last_params = allparams

        ss11   = 0
        ss110  = 0
        ss12   = 0
        ss120  = 0
        ss22   = 0
        ss220  = 0

      endif 


! Preparation of the n-exp representation of "pure" scattering functions
ilr: if( newcomp_required ) then

      call prepare_ttable1(t0, tmax, nxpoints, t_samples)

      ! linear polymer componenet: here labelled (blin) part with concentration philin
      ! modelling of the long chain linear componente by the empirical locrep scheme  
      ! linear chain function according to the reptation interpolation model locrep2
      ! HERE: parametes nlin, blin, philin --> s11
      ! 
      do i=1,nxpoints
             ts     =  t_samples(i)
             rep_Re = rep_re_fac * sqrt(rep_Ne) * llin
             tau_e  = rep_Re**4/(Pi**2*rep_Wl4)

             if(rep_clf == 0) then 
               sqtres =  reptation2_sqt  (q,ts, dble(Nlin), llin, rep_Ne, rep_Re, rep_wl4, rep_alpha0, tau_e * rep_talphamax,rep_talphawd) 
             else
               sqtres =  reptationlik_sqt(q,ts, dble(Nlin), llin, rep_Ne, rep_Re, rep_wl4, rep_alpha0, tau_e * rep_talphamax,rep_talphawd) 
             endif
             fqt = sqtres(2) / sqtres(1)
             plin0   = nlin  * Debye_qnl(q, nlin, llin)  !?!  sqtres(1) ??
             plin    = plin0 * fqt
             s_samples(i) = fqt   

!! write(*,'(a,5(a,f12.6))')"T1 lin sq: ","sqtres(1)=",sqtres(1)," nlin=",dble(nlin),"  plin0=",plin0, &
!!                         "   t=", ts, "  fqt=", fqt


          enddo
               
          rmsdev = rmsdev_limit
          call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexp1,aexp11,rexp11,rmsdev,iout)

!! checking mexpfit >>>>
         rexp11(1:nexp1) = max( rexp11(1:nexp1) , rlow )
         do j=1,nxpoints
            exp_samples(j) = sum(aexp11(1:nexp1)*exp(-rexp11(1:nexp1)*t_samples(j)))
          enddo
          call write_mexp_results(q, t_samples,s_samples, exp_samples, nxpoints, &
                                  aexp11, rexp11, nexp1, 1, iadda)
!! <<< end checking

          



          if(rmsdev > rmsdev_limit) then
             write(*,'(a)')" LOW ACC WARN"
!            write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
!            write(*,*)"         rpa_test exp model bad match 11", rmsdev , rmsdev_limit
!            write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
          else
            write(*,*)
          endif

! ring polymer component, here: matrix:  Scc (or linear??)
! modelling of the long chain linear componente by the empirical locrep scheme same as component 11
! implicit assumption: chain length ncc = nlin, contrast acc = 0, volumen fraction phicc = 1-philin-phiring 

      do i=1,nxpoints
             ts     =  t_samples(i)
             call  ring2am7(q,ts, nring, ring_ne0, lring, ring_nue, ring_nue0, ring_wl4,  &
                      ring_pmin, ring_pwidth, ring_f0, ring_finf, ring_pexinf, &
                      ring_ntran, ring_nwid, &   
                      diff,alpha,beta,r021,r022,a_cross, &
                      sqtres)                       !??
            fqt = sqtres(2)/sqtres(1)

!             plin0cc   = nring_cc  * Debye_qnl(q, nring_cc, l) 
             pring0cc   = sqtres(1)   !! TBD faktor ???? 
             pringcc    = pring0cc * fqt
             s_samples(i) = fqt 

!!write(*,'(a,5(a,f12.6))')"T3 rin sq: ","sqtres(1)=",sqtres(1)," nring=",dble(nring),"  plin0=",&
!!             nring_debye(q, lring, ring_nue, ring_nwid, dble(nring), ring_ntran, Rgr ) * nring, &
!!             "   t=", ts, "  fqt=", fqt

       enddo


      rmsdev = rmsdev_limit
      call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexpcc,aexpcc,rexpcc,rmsdev,iout)
!! checking mexpfit >>>>
         rexpcc(1:nexpcc) = max( rexpcc(1:nexpcc) , rlow )
         do j=1,nxpoints
            exp_samples(j) = sum(aexpcc(1:nexpcc)*exp(-rexpcc(1:nexpcc)*t_samples(j)))
          enddo
          call write_mexp_results(q, t_samples,s_samples, exp_samples, nxpoints, &
                                  aexpcc, rexpcc, nexpcc, 2, iadda)
!! <<< end checking




      if(rmsdev > rmsdev_limit) then
        write(*,'(a)')" LOW ACC WARN"
!        write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
!        write(*,*)"         rpa_test exp model bad match CC", rmsdev , rmsdev_limit
!        write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
      else
        write(*,*)
      endif



      ! ring polymer componnent: here labelled (bring) part with concentration phiring
      ! modelling in terms of N-Rouse model (explicit summation of Roise modes) multiplied by
      ! a center of mass diffsuion with sublinear initial phase and simple diffusion beyond r**2 = r02
      ! HERE: parameters nring, bring, phiring --> s22

      do i=1,nxpoints
             ts     =  t_samples(i)
             call  ring2am7(q,ts, nring, ring_ne0, lring, ring_nue, ring_nue0, ring_wl4,  &
                      ring_pmin, ring_pwidth, ring_f0, ring_finf, ring_pexinf, &
                      ring_ntran, ring_nwid, &   
                      diff,alpha,beta,r021,r022,a_cross, &
                      sqtres)                       !??
            fqt = sqtres(2)/sqtres(1)

!             plin0cc   = nring_cc  * Debye_qnl(q, nring_cc, l) 
             pring0   = sqtres(1)   !! TBD faktor ???? 
             pring    = pring0 * fqt
             s_samples(i) = fqt 

!!write(*,'(a,5(a,f12.6))')"T3 rin sq: ","sqtres(1)=",sqtres(1)," nring=",dble(nring),"  plin0=",&
!!             nring_debye(q, lring, ring_nue, ring_nwid, dble(nring), ring_ntran, Rgr ) * nring, &
!!             "   t=", ts, "  fqt=", fqt
 



       enddo

       call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexp2,aexp22,rexp22,rmsdev,iout)
!! checking mexpfit >>>>
         rexp22(1:nexp1) = max( rexp22(1:nexp1) , rlow )
         do j=1,nxpoints
            exp_samples(j) = sum(aexp22(1:nexp2)*exp(-rexp22(1:nexp2)*t_samples(j)))
          enddo
          call write_mexp_results(q, t_samples,s_samples, exp_samples, nxpoints, &
                                  aexp22, rexp22, nexp2, 3, iadda)
!! <<< end checking




       if(rmsdev > rmsdev_limit) then
         write(*,'(a)')" LOW ACC WARN"
!         write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
!         write(*,*)"         rpa_test exp model bad match 22", rmsdev , rmsdev_limit
!         write(*,*)"WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING"
       else
         write(*,*)
       endif


  !! communication with the rpa_laplace module 
  !! transfer of the pure function n-exp representation parameters to global variables
  !! of module rpa_laplace
   ! componente c    = lineraes polymer (matrix) gleiche Parameter wie comp 1
   ! componente 1    = lineares polymer (markierter Teil)
   ! componente 2    = ring-polymer
        phi1              =   philin   ! volume fraction of polymer component 1
        phi2              =   phiring  ! volume fraction of polymer component 2
        
        Scc00             =   pring0cc ! unperturbed structure factor S(Q) of "matrix" polymers
        S0011             =   plin0    ! unperturbed structure factor S(Q) of polymer 1
        S0022             =   pring0   ! unperturbed structure factor S(Q) of polymer 2

write(*,'(a,f10.4,a,3f12.6,a,2f12.6,a,2f12.6)')"Q=",q," Scc00, S0011, S0022 =",Scc00,S0011,S0022,&
                                               "  phi1,2=",phi1,phi2,"  blin, bring=",blin,bring

!        nexpcc            =   modeexcc   ! number of exp-functions to describe background
!        nexp1             =   modeex1   ! number of exp-functions to describe component1
!        nexp2             =   modeex2   ! number of exp-functions to describe component2
        aexp_cc(1:nexpcc) =   aexpcc(1:nexpcc)   ! amplitude coeffs for laplace-exp representation of "matrix"
        rexp_cc(1:nexpcc) =   rexpcc(1:nexpcc)   ! rate      coeffs for laplace-exp representation of "matrix"
        aexp_s1(1:nexp1)  =   aexp11(1:nexp1)    ! amplitude coeffs for laplace-exp representation of polymer 1
        rexp_s1(1:nexp1)  =   rexp11(1:nexp1)    ! rate      coeffs for laplace-exp representation of polymer 1
        aexp_s2(1:nexp2)  =   aexp22(1:nexp2)    ! amplitude coeffs for laplace-exp representation of polymer 2
        rexp_s2(1:nexp2)  =   rexp22(1:nexp2)    ! rate      coeffs for laplace-exp representation of polymer 2

!        write(6,*)"numeric integration of log spaced oscillating integral.."
        if(phi1 < 1d-5) phi1 = 1d-5
        if(phi2 < 1d-5) phi2 = 1d-5

        call prepare_intervals(900,1d-8,1d7)  !! <<<<<< Parameter ggf. noch optimiern bzw. nach aussen reichen

        call sel_scomp(1) 
        call create_coefficients(1, Ssfu)
        call sel_scomp(2)
        call create_coefficients(2, Ssfu)
        call sel_scomp(4)  
        call create_coefficients(3, Ssfu)

     endif ilr

     ss11 =  get_integral_value(1, max(t0,t))
     ss12 =  get_integral_value(2, max(t0,t))
     ss22 =  get_integral_value(3, max(t0,t))

     if(newcomp_required) then
       ss110 =  get_integral_value(1, t0)
       ss120 =  get_integral_value(2, t0)
       ss220 =  get_integral_value(3, t0)
     endif


     sqt  = ss11  * blin*blin + 2*blin*bring * ss12  +  ss22  * bring*bring 
     sqt0 = ss110 * blin*blin + 2*blin*bring * ss120 +  ss220 * bring*bring 


 select case(mode)  ! be sure that in the beginning the assignment of t, q are properly made (if i1:)
   case (0)
       th_rpalir3 = Sqt/Sqt0

  case (1) ! linear contribution at tau as function of q
     th_rpalir3 = plin

  case(2)  ! star contribution at tau  as function of q
     th_rpalir3 = pring
 
  case(3)  ! rpa contribution at tau as function of q
      th_rpalir3 = Sqt

  case default
     write(6,*)'invalid mode !'
      th_rpalir3 = 0

 end select

!!! just for testing !!!!
goto 2222
1111     continue
             call  ring2am7(q,t, nring, ring_ne0, lring, ring_nue, ring_nue0, ring_wl4,  &
                      ring_pmin, ring_pwidth, ring_f0, ring_finf, ring_pexinf, &
                      ring_ntran, ring_nwid, &   
                      diff,alpha,beta,r021,r022,a_cross, &
                      sqtres)                       !??
            th_rpalir3 = sqtres(2)/sqtres(1)

            write(*,*)"T1:",q,t,sqtres, th_rpalir3 
            write(*,*)nring, ring_ne0
            goto 2222

3333 continue
             rep_Re = rep_re_fac * sqrt(rep_Ne) * llin
             tau_e  = rep_Re**4/(Pi**2*rep_Wl4)
             sqtres = reptationlik_sqt(q,t, dble(Nlin), llin, rep_Ne, rep_Re, rep_wl4, &
                      rep_alpha0, tau_e * rep_talphamax,rep_talphawd) 
            th_rpalir3 = sqtres(2)/sqtres(1)
2222 continue
!!!!!!!!!!!!!!!!!!!!!!!!!


  th_rpalir3 =  th_rpalir3 * ampli

  call parset('nlin    ' ,sngl(1.d0*nlin),iadda) 
  call parset('nring_cc ',sngl(1.d0*nring_cc),iadda) 
  call parset('nring   ' ,sngl(1.d0*nring),iadda) 
  call parset('rep_f   ' ,sngl(1.d0*rep_re_fac),iadda) 

  call parset('llin    ',sngl(llin),iadda) 
  call parset('philin  ',sngl(philin),iadda) 
  call parset('phiring ',sngl(phiring),iadda) 
  call parset('blin    ',sngl(blin),iadda) 
  call parset('bring   ',sngl(bring),iadda) 


 CONTAINS 

! Debugging aid to assess the quality of the n-exp fit by independent datreat evaluation 
! by writing it to a file
 subroutine write_mexp_results(q, times, sqtvalues, sqtexp, nvals, aexp, rexp, mexp, ityp, irec)
   implicit none
   double precision, intent(in) :: q
   double precision, intent(in) :: times(nvals)
   double precision, intent(in) :: sqtvalues(nvals)
   double precision, intent(in) :: sqtexp(nvals)
   integer         , intent(in) :: nvals
   double precision, intent(in) :: aexp(mexp)
   double precision, intent(in) :: rexp(mexp)
   integer         , intent(in) :: mexp
   integer         , intent(in) :: ityp
   integer         , intent(in) :: irec

   character(len=80) :: fnam, buf
   integer           :: iunit, ios, i, j

   write(fnam,'(a,i0,a,i0,a)')"sqtexptst",ityp,"rec",irec,".dtr"
   open(newunit=iunit,file=trim(fnam))
   write(iunit,'(a)')trim(fnam)
   write(iunit,'(a,i0)')fnam(8:16)//" fqt  vs  t ",1000*ityp+irec
   write(iunit,'(a,e14.7)')"q    ",q
   write(iunit,'(a,i5)')   "ityp ",ityp
   write(iunit,'(a)')"values"
   do i=1,nvals
     write(iunit,'(3e16.7)') times(i), sqtvalues(i), 0d0
   enddo
   write(iunit,'(a)')"#nxt"

   write(iunit,'(a)')trim(fnam)//" expfit"
   write(iunit,'(a,i0)')fnam(8:16)//" fqt  vs  t ",10000+1000*ityp+irec
   write(iunit,'(a,e14.7)')"q    ",q
   write(iunit,'(a,i5)')   "ityp ",ityp
   write(iunit,'(a,i3)')"mexp    ",mexp
   do i=1,mexp
     write(buf,'("aexp",i0)') i
     write(iunit,'(a,e14.7)') trim(buf)//"  ",aexp(i)
     write(buf,'("rexp",i0)') i
     write(iunit,'(a,e14.7)') trim(buf)//"  ",rexp(i)
   enddo
   write(iunit,'(a)')"values"
   do i=1,nvals
     write(iunit,'(3e16.7)') times(i), sqtexp(i), 0d0
   enddo
   write(iunit,'(a)')"#nxt"

   close(iunit)
 
 end subroutine write_mexp_results 
 
! subroutines and functions entered here are private to this theory and share its variables 
      function Debye_qnl( q, nlin, l) result(p0)
        implicit none
        double precision, intent(in)   :: q
        integer         , intent(in)   :: nlin
        double precision, intent(in)   :: l
        double precision               :: p0

        double precision               :: x

        x   = (q*l)**2 * (nlin / 6d0)
        p0  = 2d0/x**2 *((exp(-x)-1d0)+x)

      end function  Debye_qnl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





function nring_debye(q, l, nu, nuwidth, n, n1, Rg ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: l           ! effective segment length
   double precision ,intent(in) :: nu          ! exponent
   double precision ,intent(in) :: nuwidth
   double precision ,intent(in) :: n           ! number of segments
   double precision ,intent(in) :: n1
   double precision ,intent(out):: Rg          ! radius of gyration
   double precision             :: val

   integer          :: i, j
   double precision :: eterms(0:nint(n)-1)
   double precision :: rterms(0:nint(n)-1)
   double precision :: Rg0, rr, mu=0.5d0


!!$OMP PARALLEL DO PRIVATE(rr)   
 do i=0,nint(n)-1
  rr = l**2 * (i*(1d0-i/N))**(nucross(dble(i),n1,nuwidth,nu,mu)*2)
  rterms(i) = 0.5d0 * rr 
  eterms(i) = exp(-q*q*rr/6d0)
 enddo
 val = nint(n) * eterms(0)
 rg0 = nint(n) * rterms(0)
!!$OMP PARALLEL DO REDUCTION(+:val,rg0)
   do i=1,nint(n)-1
      val = val + eterms(i)*(nint(n)-i) * 2
      rg0 = rg0 + rterms(i)*(nint(n)-i) * 2
   enddo
!!$OMP END PARALLEL DO

   val = val / (nint(n)**2)
   rg  = sqrt(rg0 / (nint(n)**2))

end function nring_debye

function nucross(x,x0,width,nu,mu) result(nuc)
  double precision, intent(in) :: x, x0, width, nu, mu
  double precision :: nuc
  nuc = fermi(x,x0,width)*mu + (1d0-  fermi(x,x0,width))*nu
end function nucross


function fermi(x,x0,width) result(ff)
  double precision, intent(in) :: x, x0, width
  double precision :: ff
  ff = 1d0/(1d0+exp((x-x0)/width))
  
end function fermi






end function th_rpalir3



!ring2am7(q,ts, nring, ring_ne0, lring, ring_nue, ring_nue0, ring_wl4,  &
!                      ring_pmin, ring_pwidth, ring_f0, ring_finf, ring_pexinf, &
!                      ring_ntran, ring_nwid,sqtres)        


!?! subroutine ring2am7(q,t, N, ne00, l, nue, nu0, wl4,  pmin, pwidth, f0, finf, pexinf, &
!?!                   ntrans, nwidth, &
!?!                   diff,alpha,beta,r021,r022,a_cross, &
!?!                   sqtres )    
!?! 
!?!   implicit none
!?!   double precision, intent(in)  :: q
!?!   double precision, intent(in)  :: t
!?!   integer,          intent(in)  :: N
!?!   double precision, intent(in)  :: l
!?!   double precision, intent(in)  :: ne00
!?!   double precision, intent(in)  :: nue
!?!   double precision, intent(in)  :: nu0
!?!   double precision, intent(in)  :: wl4
!?!   double precision, intent(in)  :: pmin
!?!   double precision, intent(in)  :: pwidth
!?!   double precision, intent(in)  :: f0
!?!   double precision, intent(in)  :: finf
!?!   double precision, intent(in)  :: pexinf
!?!   double precision, intent(in)  :: ntrans
!?!   double precision, intent(in)  :: nwidth
!?! 
!?!  double precision, intent(in) :: diff
!?!  double precision, intent(in) :: alpha
!?!  double precision, intent(in) :: beta
!?!  double precision, intent(in) :: r021
!?!  double precision, intent(in) :: r022
!?!  double precision, intent(in) :: a_cross
!?! 
!?!   double precision, intent(out)           :: sqtres(2)
!?!  
!?! 
!?!   double precision :: rr, sqr, sqtr
!?!   double precision :: tau_R, W, taue
!?! 
!?!   double precision :: ne0
!?! 
!?!   ne0 = ne00
!?!   rr  = 0
!?! 
!?! 
!?! 
!?!      rr=((exp(((-alpha+beta)*log(r021/r022)+alpha*beta*(log(2d0)+log(3d0)-       &
!?!          log(r022/diff)))/beta)*r022*t**alpha)**a_cross+(exp((log(2d0)+log(3d0)- &
!?!          log(r022/diff))*beta)*r022*t**beta)**a_cross+(6d0*diff*t)**a_cross)**(1d0/a_cross)
!?! 
!?!   call Nrouse_ring2am7(q,t,n, sqr, sqtr)
!?!   sqtres =  [sqr , exp(-rr * q*q / 6d0) * sqtr]
!?! 
!?! !?! write(*,*)"TESTx:",q,t,sqr,sqtr,sqtres
!?!  
!?! ! ---- writing computed parameters to the record >>>  
!?! !?!      call parset('Rg      ',sngl(Rg)   ,iadda,ier)
!?! !?!      call parset('tau_R   ',sngl(tau_R), iadda, ier)
!?! !?!      call parset('taue    ',sngl(taue) , iadda, ier)
!?! !?!      call parset('W       ',sngl(W)    , iadda, ier)
!?! !?!      call parset('Ne0     ',sngl(Ne0)  , iadda, ier)
!?! !?!      call parset('tau_pmin',sngl(tau_R/pmin**2)  , iadda, ier)
!?! !?!      call parset('tau_px  ',sngl(1/( (W*pi**2) /Ne0**2)  )  , iadda, ier)
!?! !?!      call parset('tauring ',sngl(tauring)  , iadda, ier)
!?!  
!?!  CONTAINS 
!?! 
!?!        subroutine Nrouse_ring2am7(q,t,N,sqx,sqtx)
!?! !      ================================================================
!?! !
!?! ! Rouse expression for a ring polymer
!?! ! with Nb segments each, the ideal Ree of the linear version of the polymer is given by R.
!?! ! amod contains mode modifiers
!?! ! ------------------------------------------------------------
!?! !
!?!        implicit none
!?! 
!?!        double precision kb, pi
!?!        parameter(kb=1.380662d-23)
!?!        parameter(pi=3.141592654d0)
!?! 
!?! 
!?! 
!?!        double precision, intent(in)     ::  q
!?!        double precision, intent(in)     ::  t
!?!        integer, intent(in)              ::  N
!?!        double precision, intent(out)    :: sqx
!?!        double precision, intent(out)    :: sqtx
!?! 
!?!        integer                          ::  nn,mm,p
!?! 
!?!        double precision :: rate, traf, nun ! , tau_R, W, Ne0, taue
!?! 
!?! 
!?!        double precision :: cosarray(0:N,N/2), ewfac(N/2), ff2(N/2), taus(N/2)
!?! 
!?!        double precision :: nu(0:N)
!?!        double precision :: tram(0:N)
!?! 
!?!        double precision :: msd_im, Rg, tauring
!?! 
!?!        double precision :: Sq , Sqt
!?! 
!?! 
!?!        integer          :: ip, N2
!?! 
!?! 
!?!        if(N.le.0) then
!?!          W  = 999
!?!          Sq = 999
!?!          Sqt= 999
!?!          write(6,*)'Error Number of chain segments is <= 0!',N
!?!          return
!?!        endif
!?! 
!?! ! ---- and the Rousetime ----
!?!        W     = Wl4 / l**4
!?!        tau_R = N**2/(W * Pi**2)
!?!        if(Ne0 <= 0d0)  Ne0   = N/pmin
!?!        taue  = Ne0**2 * l**4/(Wl4*pi**2)
!?! !!       tauinf= taue * pmin**pexinf !! taue * (N/Ne0)**pexinf 
!?! !!       rp=(Wl4*pi^2)*(p/pmin)**pexp/Ne0^2.
!?! 
!?! ! ---- init sums ----
!?!        N2 = N/2
!?! ! p(even) = 2*ip in the following...a
!?! 
!?! !$OMP PARALLEL DO
!?!   do nn=0,N
!?!     tram(nn) = 1d0/(1.0d0+exp((nn-ntrans)/nwidth))
!?! !                                 -----   ------
!?!   enddo
!?! !$OMP END PARALLEL DO
!?! 
!?! 
!?! !$OMP PARALLEL DO
!?!   do nn=0,N
!?!     nu(nn) = 2*(nu0*tram(nn) + (1-tram(nn))*nue)
!?!   enddo
!?! !$OMP END PARALLEL DO
!?! 
!?! 
!?! 
!?! !$OMP PARALLEL DO PRIVATE(traf,nun)
!?!        do nn=0,N
!?!         do ip=1,N2
!?!          traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
!?!          nun  = traf * (2*nue) + (1d0-traf) * (2*nu0)
!?!          ff2(ip) =  -2d0* dble(N)**nun *(l*q)**2/(3d0*pi**2)
!?!          cosarray(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))               &
!?!                           /  dble(2*ip)**(2d0)                         &
!?!                           *  (f0*(1d0-traf)+finf*(traf))
!?!         enddo
!?!        enddo
!?! !$OMP END PARALLEL DO
!?! 
!?! 
!?! !       ff2  = -2d0* dble(N)**nu *(l*q)**2/(3d0*pi**2)
!?! 
!?!        msd_im = 0
!?! !$OMP PARALLEL DO PRIVATE(rate,traf) REDUCTION(+:msd_im)
!?! !!$OMP PARALLEL DO PRIVATE(rate,traf)
!?!        do ip=1,N/2
!?!          traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
!?! 
!?!          rate      = dble(2d0*ip)**(2d0)/tau_R*(1d0-traf) + &
!?!                      (W*pi**2)* (2d0*ip/pmin)**pexinf/Ne0**2 * traf
!?! 
!?!          taus(ip)  = 1d0/rate
!?! 
!?!          ewfac(ip) = 1.0d0-exp(-t * rate )
!?! 
!?! 
!?!          msd_im = msd_im - ff2(ip)  /  dble(2*ip)**(1d0+1.d0)                      &
!?!                           *  (f0*(1d0-traf)+finf*(traf))                       &
!?!                           *  ewfac(ip)
!?!        enddo
!?! !$OMP END PARALLEL DO
!?!        msd_im = 6 * msd_im / q**2
!?! 
!?! 
!?!            Sq  = 0
!?!            Sqt = 0
!?!            Rg  = 0
!?! 
!?! 
!?! ! ---- Do the sums -----
!?! !$OMP PARALLEL DO REDUCTION(+:Sq,Sqt,Rg)
!?!            do nn = 1,N
!?!             do mm = 1,N
!?! 
!?! 
!?!               Sq  = Sq  + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) )
!?!               Sqt = Sqt + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) &
!?!                               + sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)*ff2(1:N2)) )
!?! 
!?!               Rg  = Rg  + l**2  * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) 
!?! 
!?!             enddo
!?!            enddo
!?! !$OMP END PARALLEL DO
!?! 
!?!            Sqx  = Sq /N
!?!            Sqtx = Sqt/N
!?! 
!?!            Rg  = sqrt(Rg/2d0/N**2)
!?! 
!?! !?! write(*,*)"Test2:",Sqx,Sqtx,Rg,W,l,taue
!?! 
!?!            tauring = taus(1)
!?! 
!?!        end subroutine Nrouse_ring2am7
!?! 
!?! 
!?!  end subroutine ring2am7
