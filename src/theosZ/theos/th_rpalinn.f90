 FUNCTION th_rpalinn(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!============================================================================
!  applies rpa mixture tor the scattering function for a star-linear polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))
! 
      use theory_description 
      use lmfit
      use rpa_laplace
      use polynom
      implicit none
!! 
      SAVE
!!    ====


      real    :: th_rpalinn
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
     double precision :: ampli      ! prefactor
     double precision :: wl4rous    ! rouse rate star 
     integer          :: nrouseff    ! effective segment star arm 
!     double precision :: diffrous   ! centre-of-mass diffusino star
     double precision :: diff       
     double precision :: r02   , rr     
     double precision :: tc         
     double precision :: nu_subdiff 
     double precision :: a_cross    
     double precision :: beta_q     

     double precision :: rep_b
     double precision :: rep_ne
     double precision :: rep_wl4
     double precision :: rep_re
     double precision :: rep_f_re
     double precision :: rep_wd_mod
     double precision :: rep_alpha
!?!     double precision :: locr2_b    ! locrep2 parameter b 
!?!     double precision :: locr2_a    ! locrep2 parameter a
!?!     double precision :: locr2_ta   ! locrep2 parameter tau 
!?!     double precision :: locr2_lz   ! locrep2 parameter lz 
!?!     double precision :: locr2_te   ! locrep2 parameter teps 
!?!     integer          :: nro_me     ! effective segments for Me-Modelling 
!?!     double precision :: re_me      ! Re of Me (entangle strand Re)
!?!     double precision :: wl4        ! rouse rate
!?!
!?!     double precision :: lr2_b_c    ! locrep2 parameter b 
!?!     double precision :: lr2_a_c    ! locrep2 parameter a
!?!     double precision :: lr2_ta_c   ! locrep2 parameter tau 
!?!     double precision :: lr2_lz_c   ! locrep2 parameter lz 
!?!     double precision :: lr2_te_c   ! locrep2 parameter teps 
!?!     integer          :: nro_me_c   ! effective segments for Me-Modelling 
!?!     double precision :: re_me_c    ! Re of Me (entangle strand Re)
!?!     double precision :: wl4_c      ! rouse rate
!?!
!?!     double precision :: difflin    ! diffusion linear
!?!     double precision :: betadifl   ! beta for that
!?!
!?!     double precision :: diffmatc   ! diffusion linear matrix
!?!     double precision :: betadifc   ! beta for that
  
     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset


! the recin parameter representation 
     double precision :: q          ! scattering wavevector
     double precision :: temp       ! temperature
     integer          :: nrous       ! number of segments per star arm
     double precision :: l          ! segment length
     integer          :: f          ! number of arms per star
     double precision :: phirous    ! volume fraction of star
     integer          :: nlin       ! number of segments of linear polymer
     integer          :: nlin_cc    ! number of segments of linear polymer matrix 
     double precision :: philin     ! volume fraction of linera polymer
     double precision :: betadif

     double precision :: alin       ! scattering contrast of linear polymer omponent 
     double precision :: arous      ! scattering contrast of start polymer

     integer          :: nxpoints   ! number of points in sum(exp) modelling of S00 functions 

     integer          :: mode       ! mode select 
     integer          :: modeex     ! mode select exp representation
     integer          :: modeex1      ! mode select exp representation
     integer          :: modeex2      ! mode select exp representation
     integer          :: modeexcc     ! mode select exp representation


     double precision, parameter :: tmin = 0.001d0
     double precision            :: t0  ! virtual 0 > 0 to allow the same numerics as for finite t to compute s(q,0)



     double precision :: tau        ! for modes > 0
 
     double precision   :: t

     double precision   :: sqt, sqt0, fqt
!?!     double precision   :: locrep2
     double precision   :: prous0, prous
     double precision   :: plin0, plin
     double precision   :: plin0cc, plincc
     double precision   :: dr, wx, lx
     double precision   :: Re_rous
     double precision   :: Sinv0, Sinv
!?!     double precision   :: local_reptation2
     double precision   :: repta_fqt

     double precision, parameter  :: t_table_spacing1 = 300d0 
     double precision, parameter  :: t_table_spacing2 = 100d0 
! exp(i*log(tmax*300d0)/nxpoints)/100d0
! log t table, first point at exp(log(tmax*t_table_spacing1)/nxpoints)/t_table_spacing2 
!              last  point at exp(log(tmax*t_table_spacing1))/t_table_spacing2 = tmax*(t_table_spacing1/t_table_spacing2) 



     integer            :: ifix

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

     double precision        :: a1, a2, r1, r2, r3, b1, b2, g1, g2, g3
     integer                 :: iout = 0


     complex(kind=XPREC)     :: il_coeffs11(3*mexp)
     complex(kind=XPREC)     :: il_coeffs12(3*mexp)
     complex(kind=XPREC)     :: il_coeffs22(3*mexp)
     complex(kind=XPREC)     :: il_alphas(3*mexp)
     integer                 :: nnsum
     integer                 :: j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(6,*)"mbuf=",mbuf
!
! ----- initialisation ----- 
    IF (ini.eq.0) then 
       thnam = 'rpalinn '
       nparx =       29
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpalinn = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " applies danamic rpa to the scattering function for a "//cr//parspace//&
                                " lin-lin(long)  polymer mixture and extends this into  "//cr//parspace//&
                                " with non-Gauss alpha (Guenza2014) "//cr//parspace//&
                                " alpha(1) =alpha0*exp(-(((log(t+eps)-log(talpmax))/talpwd)**2)/2)+alpoff" !//cr//parspace//&

       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor
        parnam ( 2) = 'wl4rous '  ! rouse rate star
        parnam ( 3) = 'nrouseff '  ! effective segment star arm
        parnam ( 4) = 'd0'
        parnam ( 5) = 'r02'
        parnam ( 6) = 'alpha'
        parnam ( 7) = 'a_cross'
        parnam ( 8) = 'beta_q'

        parnam ( 9) = 'alpha0     ' 
        parnam (10) = 'talphamax  ' 
        parnam (11) = 'talphawd   ' 
        parnam (12) = 'aoffset    ' 
        parnam (13) = 'rep_b      ' 
        parnam (14) = 'rep_ne     ' 
        parnam (15) = 'rep_re     ' 
        parnam (16) = 'rep_wl4    ' 
        parnam (17) = 'rep_f_re   '
        parnam (18) = 'rep_wd_mod '
        parnam (19) = 'rep_alpha  '

        parnam (20) = 'mode    '  ! what to compute
        parnam (21) = 'modeex  '  ! selects way of representation in terms of nexp
        parnam (22) = 'nxpoints'  ! number of points to determine representation in terms of nexp
        parnam (23) = 'tzero   '  ! preliminary: effective zero time for S(Q,t=0=tzero)
        parnam (24) = 'tmaxrng '  ! maximum range that shall be spanned by creating the n-exp model
        parnam (25) = 'minrate '  ! minimum rate allowed in the exp models
        parnam (26) = 'rmsdev  '  ! n-exp-fit error limit (some 1e-3)
        parnam (27) = 'xil     '  ! distance of path for inv-laplace integration 0.01
        parnam (28) = 'epap    '  ! = 1d-8      ! Apodisationfactor for numerical integrand
        parnam (29) = 'epsrpa  '  ! = 1d-5      ! accuracy parameter


  
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse short liner chains rate (phirous)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment number in summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "limiting diffusion (NMR value)" !//cr//parspace//&
        th_param_desc( 5,idesc) = " reference mean squared displacement at transition point to D0 " !//cr//parspace//&
        th_param_desc( 6,idesc) = " = nu_subdiff  exponent " !//cr//parspace//&
        th_param_desc( 7,idesc) = " transition exponent between short and long time diffusuion " !//cr//parspace//&
        th_param_desc( 8,idesc) = " heterogeneity exponent  " !//cr//parspace//&
        th_param_desc( 9,idesc) = " amlitude of non-Gaussian term in Rouse  "  !  
        th_param_desc(10,idesc) = " tau of maximum of non-Gaussiabity alpha(t) "  !   
        th_param_desc(11,idesc) = " width of maximum of non Gaussianity "  !
        th_param_desc(12,idesc) = " offset for non-Gaussianity alpha "  !    

        th_param_desc(13,idesc) = "longchain locrep parameter b" !//cr//parspace//&
        th_param_desc(14,idesc) = "longchain estim. entangelment Ne" !//cr//parspace//&
        th_param_desc(15,idesc) = "longchain Wl4" !//cr//parspace//&
        th_param_desc(16,idesc) = "longchain = entanglement Re (transversal diameter)" !//cr//parspace//&
        th_param_desc(17,idesc) = "longchain dblob = f_re * rep_re modificator" !//cr//parspace//&
        th_param_desc(18,idesc) = "transition time (t0) WIDTH (tw) modificator tw=wd_mod*t0" !//cr//parspace//&
        th_param_desc(19,idesc) = "longchain non-Gaussianity strength alpha" !//cr//parspace//&
        th_param_desc(20,idesc) = "mode: 0 normal, >0 S(x=Q,tau) 1 lin, 2 star 3 rpa" !//cr//parspace//&
        th_param_desc(21,idesc) = "modex: number of exps to represent the undistrurbed S(q,t) relaxation curves" !//cr//parspace//&
        th_param_desc(22,idesc) = "number of points to determine representation in terms of nexp" !//cr//parspace//&
        th_param_desc(23,idesc) = "preliminary: effective zero time for S(Q,t=0=tzero)" //cr//parspace//&
                                "for the numerical integration we have a smoothing of the 0-->1 jump " //cr//parspace//&
                                "at t=0 due to the pratically less than infinite width of the s-integration"  //cr//parspace//&
                                "for the analytic=1 method this may be zero, for the analytic=0 method something" //cr//parspace//&
                                "like 0.001 * tmax could be a good start to try"
        th_param_desc(24,idesc) = " maximum range that shall be spanned by creating the n-exp model" !//cr//parspace//&
        th_param_desc(25,idesc) = " minimum rate allowed in the exp models" !//cr//parspace//&
        th_param_desc(26,idesc) = " n-exp fit error-limit rmsdev_limit" !//cr//parspace//&
        th_param_desc(27,idesc) = " distance of path for inv-laplace integration 0.01"
        th_param_desc(28,idesc) = " Apodisationfactor for numerical integrand = 1d-8 "
        th_param_desc(29,idesc) = " accuracy parameter 1e-5"


! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  3,idesc) = "nrous     > number of segments per star arm"
        th_file_param(  4,idesc) = "l        > segment length"
        th_file_param(  5,idesc) = "f        > number of arms per star"
        th_file_param(  6,idesc) = "phirous  > volume fraction of star"
        th_file_param(  7,idesc) = "nlin     > number of segments of linear polymer"
        th_file_param(  8,idesc) = "nlin_cc  > number of segments of linear polymer matrix"
        th_file_param(  9,idesc) = "philin   > volume fraction of linera polymer"
        th_file_param( 10,idesc) = "tau      > tau val for mode > 0 calc"
        th_file_param( 11,idesc) = "alin     > scattering contrast linear"
        th_file_param( 12,idesc) = "arous    > scattering contrast star"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpalinn = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli       =      pa( 1)
      wl4rous     = abs( pa( 2))
      nrouseff    = nint(pa( 3))
      diff        = abs( pa (4))    ! in A**2/ns
      r02         = abs( pa (5))
      nu_subdiff  = abs( pa (6))
      a_cross     = abs( pa (7))
      beta_q      = abs( pa (8))
      alpha0      =      pa( 9)
      talphamax   =  abs(pa(10))
      talphawd    =  abs(pa(11))
      aoffset     =      pa(12)

      rep_b       = abs( pa(13))
      rep_ne      = abs( pa(14))
      rep_re      = abs( pa(15))
      rep_wl4     = abs( pa(16))
      rep_f_re    = abs( pa(17))
      rep_wd_mod  = abs( pa(18))
      rep_alpha   = abs( pa(19))

      mode        = nint(pa(20))
      modeex      = nint(pa(21))
      nxpoints    = nint(pa(22))
      t0          = abs( pa(23))
      tmax        = abs( pa(24))
      rlow        = abs( pa(25))  ! parameter in rpa_laplace
      rmsdev_limit= max(1e-4,abs( pa(26)))  ! parameter in rpa_laplace
      xil         = pa(26)
      epap        = pa(27)
      epsrpa      = pa(28)
      
      if(xil    == 0d0) xil    = 0.01d0
      if(epap   == 0d0) epap   = 1d-8
      if(epsrpa == 0d0) epsrpa = 1d-5
      



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
      call parget('nrous    ',xh,iadda,ier)
      nrous     = xh
! >>> extract: segment length
      xh = 3d0
      call parget('l       ',xh,iadda,ier)
      l        = xh
! >>> extract: number of arms per star
      xh = 1
      call parget('f       ',xh,iadda,ier)
      f        = nint(xh)
! >>> extract: volume fraction of star
      xh = 0.01d0
      call parget('phirous ',xh,iadda,ier)
      phirous  = xh
! >>> extract: number of segments of linear polymer
      xh = 100
      call parget('nlin    ',xh,iadda,ier)
      nlin     = nint(xh)
!! >>> extract: number of segments of linear polymer
      xh = 100
      call parget('nlin_cc ',xh,iadda,ier)
      nlin_cc  = nint(xh)
! >>> extract: volume fraction of linera polymer
      xh = 0.9d0
      call parget('philin  ',xh,iadda,ier)
      philin   = xh
!! >>> extract: tau for mode > 0
      xh = 0.0d0
      call parget('tau     ',xh,iadda,ier)
      tau   = xh
!! >>> extract: tau for mode > 0
      xh = 0.0d0
      call parget('alin     ',xh,iadda,ier)
      alin   = xh
!!! >>> extract: tau for mode > 0
      xh = 1.0d0
      call parget('arous    ',xh,iadda,ier)
      arous  = xh
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


      allparams(1:38) = [ wl4rous, dble(nrouseff), rep_b, rep_ne, rep_re, rep_wl4, rep_f_re, &
                         rep_wd_mod, rep_alpha , q , dble(nrous), l, dble(f),   &
                         phirous, dble(nlin), philin, alin, arous, &
                         dble(nlin_cc),  &
                         diff, r02, nu_subdiff, a_cross, beta_q , &
                         dble(mode), dble(modeex), dble(nxpoints), t0, tmax, rlow, &
                         alpha0, talphamax, talphawd, aoffset, rmsdev_limit, xil, epap, epsrpa ]  


      newcomp_required = ( sum(abs(allparams-last_params)) .ne. 0 ) .and. (modeex > 0)
    
!!?? testing if(analytic < 0) newcomp_required = .true.

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

      ! linear polymer componenet: here labelled (alin) part with concentration philin
      ! modelling of the long chain linear componente by the empirical locrep scheme  
      ! linear chain function according to the reptation interpolation model locrep2
      ! HERE: parametes nlin, alin, philin --> s11
      ! 
      do i=1,nxpoints
             ts  =  t_samples(i)
             fqt =  repta_fqt(q,ts,rep_wl4,l,rep_b,rep_ne,dble(nlin),rep_re,rep_f_re,rep_wd_mod,rep_alpha,  iadda) 

             plin0   = nlin  * Debye_qnl(q, nlin, l) 
             plin    = plin0 * fqt
             s_samples(i) = fqt   
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

      ! linear polymer componenet, here: matrix:  Scc
      ! modelling of the long chain linear componente by the empirical locrep scheme same as component 11
      ! implicit assumption: chain length ncc = nlin, contrast acc = 0, volumen fraction phicc = 1-philin-phirous 

      do i=1,nxpoints
      
             ts  =  t_samples(i)
             fqt =  repta_fqt(q,ts,rep_wl4,l,rep_b,rep_ne,dble(nlin),rep_re,rep_f_re,rep_wd_mod,rep_alpha,  iadda) 
             plin0cc   = nlin_cc  * Debye_qnl(q, nlin_cc, l) 
             plincc    = plin0cc * fqt
             s_samples(i) = fqt   
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



      ! (short) linear polymer componnent: here labelled (alin) part with concentration phirous
      ! modelling in terms of N-Rouse model (explicit summation of Roise modes) multiplied by
      ! a center of mass diffsuion with sublinear initial phase and simple diffusion beyond r**2 = r02
      ! HERE: parametes nrous, arous, phirous --> s22

      do i=1,nxpoints
             ts     =  t_samples(i)
             Re_rous = sqrt(nrous * l**2)

             dr      = 1d-20
             if(alpha0 == 0d0) then
               ifix    = 0
               call  NrousePX(q,ts,temp,Dr,wl4rous,nrouseff,Re_rous, Wx, lx,ifix, sqt0,sqt)
             else
!               call nrousalpha(q,t ,temp,Dr,wl4    ,N,       R,       W,  l, pmin,pmax, Sq,Sqt)
               call nrousalpha(q,ts,temp,Dr,wl4rous,nrouseff,Re_rous, Wx, lx, 1d0,dble(nrouseff), sqt0,sqt)
             endif

             prous0 =  nrous  * Debye_qnl(q, nrous, l) 

             rr = ((exp(-log(r02/diff/6)*nu_subdiff)*r02*ts** nu_subdiff)** a_cross +&
                  (6*diff*ts)**a_cross)**(1d0/a_cross)

             s_samples(i) = sqt/sqt0 * exp( -(q*q*rr/6d0)**beta_q )
             prous        = prous0 * s_samples(i)

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
   ! componente 2    = start-polymer
        phi1              =   philin   ! volume fraction of polymer component 1
        phi2              =   phirous  ! volume fraction of polymer component 2
        
        Scc00             =   plin0cc  ! unperturbed structure factor S(Q) of "matrix" polymers
        S0011             =   plin0    ! unperturbed structure factor S(Q) of polymer 1
        S0022             =   prous0   ! unperturbed structure factor S(Q) of polymer 2
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


     sqt  = ss11  * alin*alin + 2*alin*arous * ss12  +  ss22  * arous*arous 
     sqt0 = ss110 * alin*alin + 2*alin*arous * ss120 +  ss220 * arous*arous 


 select case(mode)  ! be sure that in the beginning the assignment of t, q are properly made (if i1:)
   case (0)
       th_rpalinn = Sqt/Sqt0

  case (1) ! linear contribution at tau as function of q
     th_rpalinn = plin

  case(2)  ! star contribution at tau  as function of q
     th_rpalinn = prous
 
  case(3)  ! rpa contribution at tau as function of q
      th_rpalinn = Sqt

  case default
     write(6,*)'invalid mode !'
      th_rpalinn = 0

 end select

  th_rpalinn =  th_rpalinn * ampli

  call parset('nlin    ',(1.0*nlin),iadda) 
  call parset('nlin_cc ',(1.0*nlin_cc),iadda) 
  call parset('nrous   ',(1.0*nrous),iadda) 
  call parset('f       ',(1.0*f),iadda) 

  call parset('l       ',sngl(l),iadda) 
  call parset('philin  ',sngl(philin),iadda) 
  call parset('phirous ',sngl(phirous),iadda) 
  call parset('alin    ',sngl(alin),iadda) 
  call parset('arous   ',sngl(arous),iadda) 


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

!        subroutine nrousalpha(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)


       subroutine NrousePX(q,t,temp,Dr,wl4,N,R, W, l,ifx , Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    ifx   ----> fixend or not
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4
       integer N, nn,mm,ifx,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix

       double precision :: cosarray(N,N), ewfac(N)

       integer :: ipmin, ipmax, i

!       integer iout
       
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
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

       if(ifx.eq.0) then
         p0fix = 0
       else
         p0fix = -0.5d0
       endif

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*(ip+p0fix)*nn)/dfloat(N)) / (ip+p0fix)
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO PRIVATE(p)    
       do i=1,N 
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*(i+p0fix))/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    

       ipmin = 1
       ipmax = N

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0


! ---- Do the sums -----
! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt) 
       do nn = 1,N
        do mm = 1,N
          Sq  = Sq  + exp( -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp( -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)  &
                          + (-2*N*(l*q)**2/(3*pi**2)* &
                  sum(ewfac(ipmin:ipmax)*cosarray(nn,ipmin:ipmax)*cosarray(mm,ipmin:ipmax))))
        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end subroutine NrousePX
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       subroutine nrousalpha(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    pmin  ----> minimum p
!    pmax  ----> maximum p
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4, pmin, pmax
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)
       double precision :: rmm, rmm0=0, fqq, fqq0

       integer :: ipmin, ipmax, i

!       integer iout
       
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! write(*,*)"Talpha:",t,alpha(t),alpha0,talphamax,talphawd

! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip 
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO    
       do i=1,N
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*i)/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    

       ipmin = max(1,nint(pmin))
       ipmax = min(N,nint(pmax))

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

       rmm  = 0
!$OMP PARALLEL DO REDUCTION(+:rmm)
       do mm = 1,N
             rmm = rmm + 4d0*N*l**2/(pi**2) * &
                   sum(cosarray(mm,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )
       enddo
!$OMP END PARALLEL DO
       rmm = rmm/N

       fqq  = 1d0 -q**2 * rmm/12d0  * alpha(t)       !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
!       fqq0 = 1d0 -q**2 * rmm0/12d0 * alpha(0d0)     !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)


!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

!          Sq  = Sq  + exp(- fqq0*(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
!          Sqt = Sqt + exp(- fqq* (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
!                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))


          Sq  = Sq  + exp(- (q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(- (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
   double precision, intent(in) :: t
   double precision             :: a

! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
!                                        are shared (common) with those of the th-function
! the model
   a = alpha0 * exp(-(((log(t+1d-8)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   
end function alpha


 end function th_rpalinn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!?!  function local_reptation2(q, t, B, L) result(val)
!?!    implicit none
!?!    double precision, intent(in)   :: q, t, B, L
!?!    double precision, parameter    :: A = 1d0
!?!    double precision, parameter    :: Pi = 3.141592653589793d0 
!?!    double precision               :: val
!?!    double precision, parameter    :: sqp = sqrt(4*atan(1d0))
!?!
!?!    double precision               :: ec1, ec2, dec, edec, z
!?!
!?!
!?!    z   = sqrt(t) * q**2
!?!
!?!    ec1 = erfc((q ** 2 * t + 3 * L) * t ** (-0.5d0) / 6d0)
!?!    ec2 = erfc(sqrt(t) * q ** 2 / 6d0)
!?!
!?!    if( isnan(ec1) ) ec1 = 1d0 - erf((q ** 2 * t + 3 * L) * t ** (-0.5d0) / 6d0)
!?!    if( isnan(ec2) ) ec2 = 1d0 - erf(sqrt(t) * q ** 2 / 6d0)
!?!
!?!    dec  =   (-ec1 + ec2 )
!?!    edec = exp(t * q ** 4 / 36d0) *  ( dec ) 
!?!
!?!    if( isnan(edec) .or. z > 50d0 ) then
!?!          edec = - ( -6/(sqrt(Pi)*z)+108/(sqrt(Pi)*z**3)-5832/(sqrt(Pi)*z**5) )
!?!    endif
!?!
!?!      val= 72d0 * (sqrt(t) * q ** 4 * exp((-2d0 * L * q ** 2 * t - 3 * L ** 2) / t / 12d0) / 36d0 + &
!?!           sqrt(Pi) * (q ** 2 * t / 3 + L) * q ** 4 * edec / 72d0 - sqrt(t) * q ** 4 / 36d0) &
!?!           * B / q ** 4 * Pi ** (-0.5d0) / L +  72d0 * (A * exp(-q ** 2 * L / 6d0) * sqrt(Pi) &
!?!           + (A * L * q ** 2 / 6d0 - A) * sqrt(Pi)) * Pi ** (-0.5d0) / q ** 4 / L
!?!
!?!  end function local_reptation2
!?!


  function repta_fqt(q, t, wl4, a, b,  ne, n, re, f_re, wd_mod, alpha0,  iadda) result(fqt)
      implicit none
      double precision, parameter :: Pi = 4*atan(1d0)

      double precision, intent(in)  :: q       !
      double precision, intent(in)  :: t       !
      double precision, intent(in)  :: wl4     ! Rouse rate
      double precision, intent(in)  :: a       ! (true) segment length of reptating polymer
      double precision, intent(in)  :: b       ! fluctaution prefactor dnse = a * sqrt(b*ne)
      double precision, intent(in)  :: ne      ! number of segments in an entanglement strand
      double precision, intent(in)  :: n       ! total number of segments
      double precision, intent(in)  :: re      ! Re of entanglement blob
      double precision, intent(in)  :: f_re    ! factor to relate Re and dblob (default=1)
      double precision, intent(in)  :: wd_mod  ! modifcation factor relating width and location of t_blob
      double precision, intent(in)  :: alpha0  ! non-Gaussian correction strength
      integer         , intent(in)  :: iadda   ! address of records used for parset

      double precision              :: fqt
     
      double precision :: sqt0, sqt, sqdeb0, sqdebt
      double precision :: dnse        ! step length (local reptation) related tube dimension
      double precision :: deff        ! effective tude diameter from lateral and logitudial dimensions
      double precision :: dblob       ! effective (scaled with f_re) blob dimension 
      double precision :: f_dnse = 0  ! admixture of dnse to dblob
      double precision :: f_deff = 0  ! admixture of deff to dblob
      double precision :: t0_locrep   ! transition time from blob to local reptation dynamics 
      double precision :: tw_locrep   ! width of transition (here: tw = wd_mod * t0)
      double precision :: talphamax   ! time of maximum non_Gaussianity  
      double precision :: talphawd    ! width of Gaussiantity maximum 

      double precision :: pmin, pmax, t0, tww, aoffset
      double precision :: w0, ww0, wr, Dr, l
      double precision :: temp = 400d0   ! virtual temperature (has no effect) 
      integer          :: n_segmen    ! effective N for summation 
      integer          :: ier                                                     
     
      W0   = wl4/a**4

      dnse = sqrt(b*ne)*a
      deff = sqrt( (dnse**2 + 2*re**2)/3d0 )

      dblob = f_re * re + f_dnse * dnse + f_deff * deff

      t0_locrep = dblob**4 /(Pi**2 * wl4)  
      tw_locrep = t0_locrep * wd_mod

      talphamax = t0_locrep
      talphawd  = log( talphamax * wd_mod )

      n_segmen = ne
      pmin     = 1
      pmax     = n_segmen
      aoffset  = 0
      temp     = 400d0     ! ist irrelevant

   


     Dr  = 1d-20
     call nrousalpha_corr(q,t,temp,Dr,wl4,N_segmen,Re, Wr, l,pmin,pmax, Sqdeb0,Sqdebt) 

     t0  = t0_locrep
     tww = tw_locrep 
     sqt0  = 1d0
     sqt   =  local_reptationdr2( q, t,         a, W0, n, ne,  b, 0)

     fqt   = sqt/sqt0  * Sqdebt/Sqdeb0


      if(iadda > 0) then
         call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
         call        parset('w0      ',sngl(W0),iadda,ier)      !     "
         call        parset('wr      ',sngl(Wr),iadda,ier)      !     "
         call        parset('wl4     ',sngl(W0*l**4),iadda,ier) !     "
         call        parset('deff    ',sngl(deff),iadda,ier) !     "
         call        parset('dnse    ',sngl(dnse),iadda,ier) !     "
         call        parset('dblob   ',sngl(dblob),iadda,ier) !     "
         call        parset('t0_locrep     ',sngl(t0),iadda,ier) !     "
         call        parset('tw_locrep     ',sngl(tww),iadda,ier) !     "
      endif

CONTAINS 


 
  function local_reptationdr2(q, t, a, W0, n, ne, b, modelr) result(val)   !! with series expansion
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: a        ! step (entanglement?) length
    double precision, intent(in)   :: W0        ! (Rouse?) rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments/entanglement
    double precision, intent(in)   :: b        ! b-fluctuation intensity
    integer         , intent(in)   :: modelr   ! full deGennes 0, or only T2 (locrep) 1
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)

    double precision :: T1, T2, T3, ec1, ec2, dec
    double precision :: x, xa
    double precision, parameter :: xlimser = 15d0
    double precision :: W, ww0

    ww0 = 1d0/ne
    W = w0*(1d0 - (1d0-ww0)/(1+exp((t-t0)/tww)))

    T1 = 72d0 * (exp(-q**2 * n * a**2 / 6d0) + q**2 * n * a**2 / 6d0 - 1d0) / q**4 / a**4
!    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
!        ne*(n/3+a**2*q**2*w*t/9)*exp(a**4*q**4*w*t/36)* &
!         (erfc(a**2*q**2*sqrt(w*t)/6)-erfc(n/(2*sqrt(w*t))+a**2*q**2*sqrt(w*t)/6)))
    T3 = 72d0/((q**4*a**4))*(exp(-(q**2*n*a**2)/6d0)+(q**2*n*a**2)/6d0-1d0)+b*ne*n/3d0
!
!    val = (T1+T2)/T3

    x = sqrt(a**4*q**4*w*t/36)
    xa= n/(2*sqrt(w*t))
    
    if( abs(x) < xlimser ) then
      dec = exp(x**2) *(erfc(x)-erfc(x+xa))
    else
      dec =  1d0/(sqrt(Pi)*x)-1d0/(2*sqrt(Pi)*x**3)+3d0/(4*sqrt(Pi)*x**5) &
          - (x**4-xa*x**3+(xa**2-1d0/2d0)*x**2+(-xa**3+3d0/2d0*xa)*x+xa**4-3*xa**2+3d0/4d0) &
           * exp(-xa**2)*exp(-2*xa*x)/(sqrt(Pi)*x**5)
    endif


    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
        ne*(n/3+a**2*q**2*w*t/9) &
        * dec )

    if(modelr == 1) then
      val = T2
    else
      val = (T1+T2)/T3
    endif


  end function local_reptationdr2




function fermi(t, t0, w) result(f)
  double precision, intent(in) :: t
  double precision, intent(in) :: t0
  double precision, intent(in) :: w
  double precision             :: f

  f = 1d0/(1d0+exp((t-t0)/w))

 end function fermi



 

       subroutine nrousalpha_corr(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    pmin  ----> minimum p
!    pmax  ----> maximum p
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4, pmin, pmax
       integer :: N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)
       double precision :: rmm, fqq, fqq0

       integer :: ipmin, ipmax, i

!       integer iout
       
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
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip 
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO    
       do i=1,N
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*i)/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    

       ipmin = max(1,nint(pmin))
       ipmax = min(N,nint(pmax))

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

       rmm = 0
!$OMP PARALLEL DO REDUCTION(+:rmm)
       do mm = 1,N
             rmm = rmm + 4d0*N*l**2/(pi**2) * &
                   sum(cosarray(mm,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )
       enddo
!$OMP END PARALLEL DO
       rmm = rmm/N

       fqq  = 1d0 -q**2 * rmm/12d0 * alpha(t)       !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
       fqq0 = 1d0                                   !! since rmm(0)=0


!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(- fqq0*(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(- fqq* (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

!          Sq  = Sq  + exp(- (q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
!          Sqt = Sqt + exp(- (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
!                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
   double precision, intent(in) :: t
   double precision             :: a

! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
!                                        are shared (common) with those of the th-function
! the model
   a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   
end function alpha


end function repta_fqt 





 
