 FUNCTION th_rpalind(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
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


      real    :: th_rpalind
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


     double precision :: locr2_b    ! locrep2 parameter b 
     double precision :: locr2_a    ! locrep2 parameter a
     double precision :: locr2_ta   ! locrep2 parameter tau 
     double precision :: locr2_lz   ! locrep2 parameter lz 
     double precision :: locr2_te   ! locrep2 parameter teps 
     integer          :: nro_me     ! effective segments for Me-Modelling 
     double precision :: re_me      ! Re of Me (entangle strand Re)
     double precision :: wl4        ! rouse rate

     double precision :: lr2_b_c    ! locrep2 parameter b 
     double precision :: lr2_a_c    ! locrep2 parameter a
     double precision :: lr2_ta_c   ! locrep2 parameter tau 
     double precision :: lr2_lz_c   ! locrep2 parameter lz 
     double precision :: lr2_te_c   ! locrep2 parameter teps 
     integer          :: nro_me_c   ! effective segments for Me-Modelling 
     double precision :: re_me_c    ! Re of Me (entangle strand Re)
     double precision :: wl4_c      ! rouse rate

     double precision :: difflin    ! diffusion linear
     double precision :: betadifl   ! beta for that

     double precision :: diffmatc   ! diffusion linear matrix
     double precision :: betadifc   ! beta for that
  
     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset
     double precision :: alphaq4    ! q**6 non Gaussian prefactor
     double precision :: fdiffal    ! non_Gauss in diff
     double precision :: xnexp      ! alpha0 N exp


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

     double precision   :: sqt, sqt0
     double precision   :: locrep2
     double precision   :: prous0, prous
     double precision   :: plin0, plin
     double precision   :: plin0cc, plincc
     double precision   :: dr, wx, lx
     double precision   :: Re_rous
     double precision   :: Sinv0, Sinv
     double precision   :: local_reptation2

     double precision   :: fqalpha

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

     integer, parameter      :: mxpoints = 1000
     double precision        :: t_samples(mxpoints)
     double precision        :: s_samples(mxpoints)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(6,*)"mbuf=",mbuf
!
! ----- initialisation ----- 
    IF (ini.eq.0) then 
       thnam = 'rpalind '
       nparx =       41
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpalind = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " applies danamic rpa to the scattering function for a "//cr//parspace//&
                                " lin-lin(long)  polymer mixture and extends this into  "//cr//parspace//&
                                " with non-Gauss alpha (Guenza2014) on internal modes AND *fdiffal in com diffusion "//cr//parspace//&
                                " alpha(1) =alpha0*(nrous/100)**xnexp*exp(-(((log(t+eps)-log(talpmax))/talpwd)**2)/2)+alpoff" !//cr//parspace//&

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

 
        parnam ( 7+2) = 'locr2_b '  ! locrep2 parameter b
        parnam ( 8+2) = 'locr2_a '  ! locrep2 parameter a
        parnam ( 9+2) = 'locr2_ta'  ! locrep2 parameter tau
        parnam (10+2) = 'locr2_lz'  ! locrep2 parameter lz
        parnam (11+2) = 'locr2_te'  ! locrep2 parameter teps
        parnam (12+2) = 'nro_me  '  ! effective segments for Me-Modelling
        parnam (13+2) = 're_me   '  ! Re of Me (entangled strand Re)
        parnam (14+2) = 'wl4     '  ! Wl4 (local entanglement strand Re)
        parnam (15+2) = 'difflin '  ! center of mass diffusion
        parnam (16+2) = 'betadifl'  ! beta of linear polymer diff

        parnam (17+2) = 'lr2_b_c '  ! locrep2 parameter b for matrix c
        parnam (18+2) = 'lr2_a_c '  ! locrep2 parameter a for matrix c
        parnam (19+2) = 'lr2_ta_c'  ! locrep2 parameter tau for matrix c 
        parnam (20+2) = 'lr2_lz_c'  ! locrep2 parameter lz for matrix c 
        parnam (21+2) = 'lr2_te_c'  ! locrep2 parameter teps for matrix c 
        parnam (22+2) = 'nro_me_c'  ! effective segments for Me-Modelling for matrix c 
        parnam (23+2) = 're_me_c '  ! Re of Me (entangled strand Re) for matrix c 
        parnam (24+2) = 'wl4_c   '  ! Wl4 (local entanglement strand Re) for matrix c 
        parnam (25+2) = 'diffmatc'  ! center of mass diffusion for matrix c 
        parnam (26+2) = 'betadifc'  ! beta of linear polymer diff for matrix c 

        parnam (27+2) = 'mode    '  ! what to compute
        parnam (28+2) = 'modeex  '  ! selects way of representation in terms of nexp
        parnam (29+2) = 'nxpoints'  ! number of points to determine representation in terms of nexp
        parnam (30+2) = 'tzero   '  ! preliminary: effective zero time for S(Q,t=0=tzero)
        parnam (31+2) = 'tmaxrng '  ! maximum range that shall be spanned by creating the n-exp model
        parnam (32+2) = 'minrate '  ! minimum rate allowed in the exp models
        parnam (33+2) = 'alpha0  '  !                                                     
        parnam (34+2) = 'talpmax '  !                                                     
        parnam (35+2) = 'talpwd  '  !                                                    
        parnam (36+2) = 'alpoffs '  !    
        parnam (36+3) = 'alphaq4 '  !    
        parnam (36+4) = 'fdiffal '  ! 
        parnam (36+5) = 'xnexp   '  !    
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse short liner chains rate (phirous)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment number in summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "limiting diffusion (NMR value)" !//cr//parspace//&
        th_param_desc( 5,idesc) = " reference mean squared displacement at transition point to D0 " !//cr//parspace//&
        th_param_desc( 6,idesc) = " = nu_subdiff  exponent " !//cr//parspace//&
        th_param_desc( 7,idesc) = " transition exponent between short and long time diffusuion " !//cr//parspace//&
        th_param_desc( 8,idesc) = " heterogeneity exponent  " !//cr//parspace//&


        th_param_desc( 7+2,idesc) = "locrep2 parameter b" !//cr//parspace//&
        th_param_desc( 8+2,idesc) = "locrep2 parameter a" !//cr//parspace//&
        th_param_desc( 9+2,idesc) = "locrep2 parameter tau" !//cr//parspace//&
        th_param_desc(10+2,idesc) = "locrep2 parameter lz" !//cr//parspace//&
        th_param_desc(11+2,idesc) = "locrep2 parameter teps" !//cr//parspace//&
        th_param_desc(12+2,idesc) = "effective segments for Me-Modelling" !//cr//parspace//&
        th_param_desc(13+2,idesc) = "Re of Me (entangle strand Re)" !//cr//parspace//&
        th_param_desc(14+2,idesc) = "Wl4 entaglement strand" !//cr//parspace//&
        th_param_desc(15+2,idesc) = "centre-of-mass diffusion linear" !//cr//parspace//&
        th_param_desc(16+2,idesc) = "betadiff linear" !//cr//parspace//&

        th_param_desc(17+2,idesc) = "locrep2 parameter b  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(18+2,idesc) = "locrep2 parameter a  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(19+2,idesc) = "locrep2 parameter tau for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(20+2,idesc) = "locrep2 parameter lz  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(21+2,idesc) = "locrep2 parameter teps  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(22+2,idesc) = "effective segments for Me-Modelling  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(23+2,idesc) = "Re of Me (entangle strand Re)  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(24+2,idesc) = "Wl4 entaglement strand  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(25+2,idesc) = "centre-of-mass diffusion linear  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(26+2,idesc) = "betadiff linear  for matrix polymer (c)" !//cr//parspace//&

        th_param_desc(27+2,idesc) = "mode: 0 normal, >0 S(x=Q,tau) 1 lin, 2 star 3 rpa" !//cr//parspace//&
        th_param_desc(28+2,idesc) = "modex: number of exps to represent the undistrurbed S(q,t) relaxation curves" !//cr//parspace//&
        th_param_desc(29+2,idesc) = "number of points to determine representation in terms of nexp" !//cr//parspace//&
        th_param_desc(30+2,idesc) = "preliminary: effective zero time for S(Q,t=0=tzero)" //cr//parspace//&
                                  "for the numerical integration we have a smoothing of the 0-->1 jump " //cr//parspace//&
                                  "at t=0 due to the pratically less than infinite width of the s-integration"  //cr//parspace//&
                                  "for the analytic=1 method this may be zero, for the analytic=0 method something" //cr//parspace//&
                                  "like 0.001 * tmax could be a good start to try"
        th_param_desc(31+2,idesc) = " maximum range that shall be spanned by creating the n-exp model" !//cr//parspace//&
        th_param_desc(32+2,idesc) = " minimum rate allowed in the exp models" !//cr//parspace//&
        th_param_desc(33+2,idesc) = " amlitude of non-Gaussian term in Rouse  "  !  
        th_param_desc(34+2,idesc) = " tau of maximum of non-Gaussiabity alpha(t) "  !   
        th_param_desc(35+2,idesc) = " width of maximum of non Gaussianity "  !
        th_param_desc(36+2,idesc) = " offset for non-Gaussianity alpha "  !    
        th_param_desc(36+3,idesc) = " q**6 non-Gaussian prefcator "  !    
        th_param_desc(36+4,idesc) = " admixture of f(q**2) non-Gauss to diffusion "  !    
        th_param_desc(36+5,idesc) = " exponent x  N**x of alpha0 dependence  "  !    
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  3,idesc) = "nrous    > number of segments in short (Rouse) chain"
        th_file_param(  4,idesc) = "l        > segment length"
        th_file_param(  5,idesc) = "f        > number of arms per star OBSOLETE"
        th_file_param(  6,idesc) = "phirous  > volume fraction of short Rouse chains"
        th_file_param(  7,idesc) = "nlin     > number of segments of linear polymer"
        th_file_param(  8,idesc) = "nlin_cc  > number of segments of linear polymer matrix"
        th_file_param(  9,idesc) = "philin   > volume fraction of linera polymer"
        th_file_param( 10,idesc) = "tau      > tau val for mode > 0 calc"
        th_file_param( 11,idesc) = "alin     > scattering contrast linear"
        th_file_param( 12,idesc) = "arous    > scattering contrast short Rouse chains"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpalind = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4rous  = abs( pa( 2))
      nrouseff = nint(pa( 3))
      diff        = abs( pa (4) )    ! in A**2/ns
      r02         = abs( pa (5) )
      nu_subdiff  = abs( pa (6) )
      a_cross     = abs( pa (7) )
      beta_q      = abs( pa (8) )

      locr2_b  = abs( pa( 7+2))
      locr2_a  = abs( pa( 8+2))
      locr2_ta = abs( pa( 9+2))
      locr2_lz = abs( pa(10+2))
      locr2_te = abs( pa(11+2))
      nro_me   = nint(pa(12+2))
      re_me    = abs( pa(13+2))
      wl4      = abs( pa(14+2))
      difflin  = abs( pa(15+2))
      betadifl = abs( pa(16+2))

      lr2_b_c  = abs( pa(17+2))
      lr2_a_c  = abs( pa(18+2))
      lr2_ta_c = abs( pa(19+2))
      lr2_lz_c = abs( pa(20+2))
      lr2_te_c = abs( pa(21+2))
      nro_me_c = nint(pa(22+2))
      re_me_c  = abs( pa(23+2))
      wl4_c    = abs( pa(24+2))
      diffmatc = abs( pa(25+2))
      betadifc = abs( pa(26+2))

      mode     = nint(pa(27+2))
      modeex   = nint(pa(28+2))
      nxpoints = nint(pa(29+2))
      t0       = abs( pa(30+2))
      tmax     = abs( pa(31+2))
      rlow     = abs( pa(32+2))  ! parameter in rpa_laplace
      alpha0   =      pa(33+2)
      talphamax=  abs(pa(34+2))
      talphawd =  abs(pa(35+2))
      aoffset  =      pa(36+2)
      alphaq4  =      pa(36+3)

      fdiffal  =      pa(36+4)
      xnexp    =      pa(36+5)

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


      allparams(1:51) = [ wl4rous, dble(nrouseff), locr2_b, locr2_a, locr2_ta, locr2_lz, locr2_te, &
                         dble(nro_me) , re_me,   wl4     , betadif , q , dble(nrous), l, dble(f),   &
                         phirous, dble(nlin), philin, alin, arous, &
                         dble(nlin_cc),  &
                         lr2_b_c, lr2_a_c, lr2_ta_c, lr2_lz_c, lr2_te_c, &
                         dble(nro_me_c) , re_me_c,   wl4_c     , diffmatc, betadifc, difflin, betadifl, &
                         diff, r02, nu_subdiff, a_cross, beta_q , &
                         dble(mode), dble(modeex), dble(nxpoints), t0, tmax, rlow, &
                         alpha0, talphamax, talphawd, aoffset, alphaq4, fdiffal, xnexp ]  

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
             ts      =  t_samples(i)

             sqt0    =  local_reptation2( q*locr2_a, locr2_te   , locr2_b, locr2_lz)
             sqt     =  local_reptation2( q*locr2_a, ts/locr2_ta, locr2_b, locr2_lz)
             locrep2 =  sqt/sqt0
             dr      = 1d-20
             ifix    = 0
             call  NrousePX(q,ts,temp,Dr,wl4,nro_me,re_me, Wx, lx,ifix, sqt0,sqt)
             plin0   = nlin  * Debye_qnl(q, nlin, l) 
             plin    = plin0 * locrep2 * sqt / sqt0
             s_samples(i) = locrep2 * sqt / sqt0 *  exp( -difflin * q*q * (ts)**betadifl )  
          enddo
               
          rmsdev = rmsdev_limit
          call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexp1,aexp11,rexp11,rmsdev,iout)

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
             ts      =  t_samples(i)

             sqt0    =  local_reptation2( q*lr2_a_c, lr2_te_c   , lr2_b_c, lr2_lz_c)
             sqt     =  local_reptation2( q*lr2_a_c, ts/lr2_ta_c, lr2_b_c, lr2_lz_c)
             locrep2 =  sqt/sqt0
             dr      = 1d-20
             ifix    = 0
             call  NrousePX(q,ts,temp,Dr,wl4_c,nro_me_c,re_me_c, Wx, lx,ifix, sqt0,sqt)
             plin0cc   = nlin_cc  * Debye_qnl(q, nlin_cc, l) 
             plincc    = plin0cc * locrep2 * sqt / sqt0
             s_samples(i) = locrep2 * sqt / sqt0 *  exp( -diffmatc * q*q * (ts)**betadifc )    
      enddo

      rmsdev = rmsdev_limit
      call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexpcc,aexpcc,rexpcc,rmsdev,iout)
 

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
               call nrousalpha(q,ts,temp,Dr,wl4rous,nrouseff,Re_rous, Wx, lx, 1d0,dble(nrouseff), sqt0,sqt, fqalpha)
             endif

             prous0 =  nrous  * Debye_qnl(q, nrous, l) 

             rr = ((exp(-log(r02/diff/6)*nu_subdiff)*r02*ts** nu_subdiff)** a_cross +&
                  (6*diff*ts)**a_cross)**(1d0/a_cross)

             rr = (1d0 + fdiffal * (fqalpha-1d0)) * rr

             s_samples(i) = sqt/sqt0 * exp( -(q*q*rr/6d0)**beta_q )
             prous        = prous0 * s_samples(i)

       enddo
       rmsdev = rmsdev_limit
       call match_exp_auto(t_samples,s_samples,nxpoints,modeex,nexp2,aexp22,rexp22,rmsdev,iout)

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
       th_rpalind = Sqt/Sqt0

  case (1) ! linear contribution at tau as function of q
     th_rpalind = plin

  case(2)  ! star contribution at tau  as function of q
     th_rpalind = prous
 
  case(3)  ! rpa contribution at tau as function of q
      th_rpalind = Sqt

  case default
     write(6,*)'invalid mode !'
      th_rpalind = 0

 end select

  th_rpalind =  th_rpalind * ampli

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



       subroutine nrousalpha(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt, fqq)
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
       double precision :: rmm, fqq, fqq0, q2rmm

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

       rmm = 0
!$OMP PARALLEL DO REDUCTION(+:rmm)
       do mm = 1,N
             rmm = rmm + 4d0*N*l**2/(pi**2) * &
                   sum(cosarray(mm,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )
       enddo
!$OMP END PARALLEL DO
       rmm = rmm/N

!       fqq  = 1d0 -q**2 * rmm/12d0 * alpha(t)       !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
!       fqq0 = 1d0 -q**2 * rmm/12d0 * alpha(0d0)     !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
       q2rmm = rmm * q**2
       fqq  = 1d0 -(q2rmm-alphaq4*q2rmm**2)/12d0 * alpha(t)   !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
!       fqq0 = 1d0 -(q2rmm-alphaq4*q2rmm**2)/12d0 * alpha(0d0)  !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
! with a heuristic q**6 contribution with factor alphaq4 included


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
!   a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   a = alpha0 * (nrous/100d0)**xnexp * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   
end function alpha








 end function th_rpalind


