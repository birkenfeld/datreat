 FUNCTION th_rpalin(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
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


      real    :: th_rpalin
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
     double precision :: diffrous   ! centre-of-mass diffusino star
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
     double precision :: plimit     ! limit of star modes (internal EV value) which are frozen

     double precision :: alin       ! scattering contrast of linear polymer omponent 
     double precision :: arous      ! scattering contrast of start polymer

     integer          :: nxpoints   ! number of points in sum(exp) modelling of S00 functions 

     integer          :: mode       ! mode select 
     integer          :: modeex     ! mode select exp representation


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
     double precision        :: ts, ssq

     double precision        :: ss11, ss110, ss12, ss120, ss22, ss220

     double precision        :: a1, a2, r1, r2, r3, b1, b2, g1, g2, g3
     integer                 :: analytic = 0


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
       thnam = 'rpalin  '
       nparx =       32
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpalin = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " applies rpa mixture tor the scattering function for a star-linear polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor
        parnam ( 2) = 'wl4rous '  ! rouse rate star
        parnam ( 3) = 'nrouseff '  ! effective segment star arm
        parnam ( 4) = 'plimit  '  ! mode ev limit for star contribution
        parnam ( 5) = 'diffrous'  ! centre-of-mass diffusino star
        parnam ( 6) = 'betadifs'  ! beta of stardiff
 
        parnam ( 7) = 'locr2_b '  ! locrep2 parameter b
        parnam ( 8) = 'locr2_a '  ! locrep2 parameter a
        parnam ( 9) = 'locr2_ta'  ! locrep2 parameter tau
        parnam (10) = 'locr2_lz'  ! locrep2 parameter lz
        parnam (11) = 'locr2_te'  ! locrep2 parameter teps
        parnam (12) = 'nro_me  '  ! effective segments for Me-Modelling
        parnam (13) = 're_me   '  ! Re of Me (entangled strand Re)
        parnam (14) = 'wl4     '  ! Wl4 (local entanglement strand Re)
        parnam (15) = 'difflin '  ! center of mass diffusion
        parnam (16) = 'betadifl'  ! beta of linear polymer diff

        parnam (17) = 'lr2_b_c '  ! locrep2 parameter b for matrix c
        parnam (18) = 'lr2_a_c '  ! locrep2 parameter a for matrix c
        parnam (19) = 'lr2_ta_c'  ! locrep2 parameter tau for matrix c 
        parnam (20) = 'lr2_lz_c'  ! locrep2 parameter lz for matrix c 
        parnam (21) = 'lr2_te_c'  ! locrep2 parameter teps for matrix c 
        parnam (22) = 'nro_me_c'  ! effective segments for Me-Modelling for matrix c 
        parnam (23) = 're_me_c '  ! Re of Me (entangled strand Re) for matrix c 
        parnam (24) = 'wl4_c   '  ! Wl4 (local entanglement strand Re) for matrix c 
        parnam (25) = 'diffmatc'  ! center of mass diffusion for matrix c 
        parnam (26) = 'betadifc'  ! beta of linear polymer diff for matrix c 

        parnam (27) = 'mode    '  ! what to compute
        parnam (28) = 'modeex  '  ! selects way of representation in terms of nexp
        parnam (29) = 'nxpoints'  ! number of points to determine representation in terms of nexp
        parnam (30) = 'tzero   '  ! preliminary: effective zero time for S(Q,t=0=tzero)
        parnam (31) = 'tmaxrng '  ! maximum range that shall be spanned by creating the n-exp model
        parnam (32) = 'minrate '  ! minimum rate allowed in the exp models

! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse short liner chains rate (phirous)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment number in summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "NOT USED !          " !//cr//parspace//&
        th_param_desc( 5,idesc) = "centre-of-mass diffusino star" !//cr//parspace//&
        th_param_desc( 6,idesc) = "betadiff" !//cr//parspace//&

        th_param_desc( 7,idesc) = "locrep2 parameter b" !//cr//parspace//&
        th_param_desc( 8,idesc) = "locrep2 parameter a" !//cr//parspace//&
        th_param_desc( 9,idesc) = "locrep2 parameter tau" !//cr//parspace//&
        th_param_desc(10,idesc) = "locrep2 parameter lz" !//cr//parspace//&
        th_param_desc(11,idesc) = "locrep2 parameter teps" !//cr//parspace//&
        th_param_desc(12,idesc) = "effective segments for Me-Modelling" !//cr//parspace//&
        th_param_desc(13,idesc) = "Re of Me (entangle strand Re)" !//cr//parspace//&
        th_param_desc(14,idesc) = "Wl4 entaglement strand" !//cr//parspace//&
        th_param_desc(15,idesc) = "centre-of-mass diffusion linear" !//cr//parspace//&
        th_param_desc(16,idesc) = "betadiff linear" !//cr//parspace//&

        th_param_desc(17,idesc) = "locrep2 parameter b  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(18,idesc) = "locrep2 parameter a  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(19,idesc) = "locrep2 parameter tau for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(20,idesc) = "locrep2 parameter lz  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(21,idesc) = "locrep2 parameter teps  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(22,idesc) = "effective segments for Me-Modelling  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(23,idesc) = "Re of Me (entangle strand Re)  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(24,idesc) = "Wl4 entaglement strand  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(25,idesc) = "centre-of-mass diffusion linear  for matrix polymer (c)" !//cr//parspace//&
        th_param_desc(26,idesc) = "betadiff linear  for matrix polymer (c)" !//cr//parspace//&

        th_param_desc(27,idesc) = "mode: 0 normal, >0 S(x=Q,tau) 1 lin, 2 star 3 rpa" !//cr//parspace//&
        th_param_desc(28,idesc) = "modex: number of exps to represent the undistrurbed S(q,t) relaxation curves" !//cr//parspace//&
        th_param_desc(29,idesc) = "number of points to determine representation in terms of nexp" !//cr//parspace//&
        th_param_desc(30,idesc) = "preliminary: effective zero time for S(Q,t=0=tzero)" //cr//parspace//&
                                  "for the numerical integration we have a smoothing of the 0-->1 jump " //cr//parspace//&
                                  "at t=0 due to the pratically less than infinite width of the s-integration"  //cr//parspace//&
                                  "for the analytic=1 method this may be zero, for the analytic=0 method something" //cr//parspace//&
                                  "like 0.001 * tmax could be a good start to try"
        th_param_desc(31,idesc) = " maximum range that shall be spanned by creating the n-exp model" !//cr//parspace//&
        th_param_desc(32,idesc) = " minimum rate allowed in the exp models" !//cr//parspace//&
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
        th_file_param( 13,idesc) = "analytic > if=1 use polynomial method"
        th_file_param( 14,idesc) = "dss      > polynom sampling distance used for polynomial method only" //cr//parspace//&
                                   "           this is critical for numerical accuracy/stability try other" //cr//parspace//&
                                   "           values if the computation fails "
        th_file_param( 15,idesc) = "sfak0    > polynom scaling factor used for polynomial method only" //cr//parspace//&
                                   "           this is critical for numerical accuracy/stability try other" //cr//parspace//&
                                   "           values if the computation fails "
        th_file_param( 16,idesc) = "npp_plus > polynom sampling extra points used for polynomial method only" //cr//parspace//&
                                   "           increasing the number npp in excess of the minimum needed" //cr//parspace//&
                                   "           may improve accuracy, eventually decrease dss when increasing npp  "
        th_file_param( 17,idesc) = "xil      > for analytic=0 (or not present): integration for invlaplace" //cr//parspace//&
                                   "           xil is the displacement of the -I*inf...I*Inf integration path" //cr//parspace//&
                                   "           along the real axis, a very small 1e-3 value shoul usually be good  "
        th_file_param( 18,idesc) = "epap     > for analytic=0 (or not present): integration for invlaplace" //cr//parspace//&
                                   "           epap is an optional apodisation factor" //cr//parspace//&
                                   "           a very small 1e-9 value shoul usually be good  "
        th_file_param( 19,idesc) = "epsrpa   > for analytic=0 (or not present): integration for invlaplace" //cr//parspace//&
                                   "           epsrpa is an accuracy target for the numerical integratio" //cr//parspace//&
                                   "           a 1e-5 value should usually be good  "
    
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpalin = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4rous  = abs( pa( 2))
      nrouseff = nint(pa( 3))
      plimit   = abs( pa( 4))
      diffrous = abs( pa( 5))
      betadif  = abs( pa( 6))

      locr2_b  = abs( pa( 7))
      locr2_a  = abs( pa( 8))
      locr2_ta = abs( pa( 9))
      locr2_lz = abs( pa(10))
      locr2_te = abs( pa(11))
      nro_me   = nint(pa(12))
      re_me    = abs( pa(13))
      wl4      = abs( pa(14))
      difflin  = abs( pa(15))
      betadifl = abs( pa(16))

      lr2_b_c  = abs( pa(17))
      lr2_a_c  = abs( pa(18))
      lr2_ta_c = abs( pa(19))
      lr2_lz_c = abs( pa(20))
      lr2_te_c = abs( pa(21))
      nro_me_c = nint(pa(22))
      re_me_c  = abs( pa(23))
      wl4_c    = abs( pa(24))
      diffmatc = abs( pa(25))
      betadifc = abs( pa(26))

      mode     = nint(pa(27))
      modeex   = nint(pa(28))
      nxpoints = nint(pa(29))
      t0       = abs( pa(30))
      tmax     = abs( pa(31))
      rlow     = abs( pa(32))  ! parameter in rpa_laplace


      t0       = max(t0      ,tmin)
      nxpoints = max(nxpoints  , 7)
      modeex   = max(modeex    , 2)

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
!!! >>> tentative only testing
      xh = 0d0
      call parget('analytic ',xh,iadda,ier)
      analytic  = nint(xh) 
!!! >>> tentative only testing
      xh =  0.0001d0
      call parget('xil      ',xh,iadda,ier)
      xil = xh
!!! >>> tentative only testing
      xh =  1d-9
      call parget('epap      ',xh,iadda,ier)
      epap = xh
!!! >>> tentative only testing
      xh =  1d-5
      call parget('epsrpa    ',xh,iadda,ier)
      epsrpa = xh
!!!! >>> in the polynomial analytic this determines the spcing of the polynomial testpoints
!!      this parameter is critical for the numerical accuracy 
!!      TBD determine it automatically
!!      just now: try it out....
      xh =  0.001
      call parget('dss       ',xh,iadda,ier)
      dss = xh
!!! >>> in the polynomial analytic this determines the scaling of the polynomial testpoints
!!      this parameter is critical for the numerical accuracy 
!!      TBD determine it automatically
!!      just now: try it out....
      xh =  0.001
      call parget('sfak0     ',xh,iadda,ier)
      sfak0 = xh
!!! >>> in the polynomial analytic this determines the excess point number of tables for polyfit
!!      this parameter is critical for the numerical accuracy 
!!      this and the two preceeding parameters are in the global section of module polynom
!!      just now: try it out....
      xh =  30
      call parget('npp_plus  ',xh,iadda,ier)
      npp_plus = nint(xh)
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


      allparams(1:42) = [ wl4rous, dble(nrouseff), diffrous, locr2_b, locr2_a, locr2_ta, locr2_lz, locr2_te, &
                         dble(nro_me) , re_me,   wl4     , betadif, plimit , q , dble(nrous), l, dble(f),   &
                         phirous, dble(nlin), philin, alin, arous, &
                         dble(nlin_cc), dble(analytic), dss, sfak0, dble(npp_plus),xil,epap,epsrpa, &
                         lr2_b_c, lr2_a_c, lr2_ta_c, lr2_lz_c, lr2_te_c, &
                         dble(nro_me_c) , re_me_c,   wl4_c     , diffmatc, betadifc, difflin, betadifl ]

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



! linear chain function according to the reptation interpolation model locrep2
ilr: if( newcomp_required ) then

      ! linear polymer componenet: here labelled part with concentration philin
      ! modelling of the long chain linear componente by the empirical locrep scheme  

      do i=1,nxpoints
             ts      =  exp(i*log(tmax*t_table_spacing1)/nxpoints)/t_table_spacing2

             sqt0    =  local_reptation2( q*locr2_a, locr2_te   , locr2_b, locr2_lz)
             sqt     =  local_reptation2( q*locr2_a, ts/locr2_ta, locr2_b, locr2_lz)
             locrep2 =  sqt/sqt0
             dr      = 1d-20
             ifix    = 0
             call  NrouseY(q,ts,temp,Dr,wl4,nro_me,re_me, Wx, lx,ifix, sqt0,sqt)
             plin0   = nlin  * Debye_qnl(q, nlin, l) 
             plin    = plin0 * locrep2 * sqt / sqt0
             t_samples(i) = ts
             s_samples(i) = locrep2 * sqt / sqt0 *  exp( -difflin * q*q * (ts)**betadifl )  
          enddo
           
          call nexp_match(t_samples,s_samples,nxpoints,modeex,aexp11,rexp11,ssq)

          if(ssq > 1d-4) then
            write(6,*)"rpa_test exp model bad match 11", ssq
          endif


      ! linear polymer componenet, here: matrix
      ! modelling of the long chain linear componente by the empirical locrep scheme  

      do i=1,nxpoints
             ts      =  exp(i*log(tmax*t_table_spacing1)/nxpoints)/t_table_spacing2

             sqt0    =  local_reptation2( q*lr2_a_c, lr2_te_c   , lr2_b_c, lr2_lz_c)
             sqt     =  local_reptation2( q*lr2_a_c, ts/lr2_ta_c, lr2_b_c, lr2_lz_c)
             locrep2 =  sqt/sqt0
             dr      = 1d-20
             ifix    = 0
             call  NrouseY(q,ts,temp,Dr,wl4_c,nro_me_c,re_me_c, Wx, lx,ifix, sqt0,sqt)
             plin0cc   = nlin_cc  * Debye_qnl(q, nlin_cc, l) 
             plincc    = plin0cc * locrep2 * sqt / sqt0
             t_samples(i) = ts
             s_samples(i) = locrep2 * sqt / sqt0 *  exp( -diffmatc * q*q * (ts)**betadifc )    
      enddo
           
          call nexp_match(t_samples,s_samples,nxpoints,modeex,aexpcc,rexpcc,ssq)
          if(ssq > 1d-4) then
            write(6,*)"rpa_test exp model bad match cc", ssq
          endif


     ! star ( primary sample component) here the component with index 2

      do i=1,nxpoints
             ts     =  exp(i*log(tmax*t_table_spacing1)/nxpoints)/t_table_spacing2

             Re_rous = sqrt(nrous * l**2)


             dr      = 1d-20
             ifix    = 0
             call  NrouseY(q,ts,temp,Dr,wl4rous,nrouseff,Re_rous, Wx, lx,ifix, sqt0,sqt)
             prous0 =  nrous  * Debye_qnl(q, nrous, l) 

             prous  =  sqt/sqt0 * prous0                                         !! osolete ?
             prous  =  prous * exp( -diffrous * q*q * (ts)**betadif )            !! osolete ? add diffusion expilicitly

             t_samples(i) = ts
             s_samples(i) = sqt/sqt0 *  exp( -diffrous * q*q * (ts)**betadif )    
       enddo
          call nexp_match(t_samples,s_samples,nxpoints,modeex,aexp22,rexp22,ssq)

           if(ssq > 1d-4) then
            write(6,*)"rpa_test exp model bad match 22", ssq
           endif


   
   !
   ! componente c    = lineraes polymer (matrix) gleiche Parameter wie comp 1
   ! componente 1    = lineares polymer (markierter Teil)
   ! componente 2    = start-polymer
        phi1              =   philin   ! volume fraction of polymer component 1
        phi2              =   phirous  ! volume fraction of polymer component 2
        
        Scc00             =   plin0cc  ! unperturbed structure factor S(Q) of "matrix" polymers
        S0011             =   plin0    ! unperturbed structure factor S(Q) of polymer 1
        S0022             =   prous0   ! unperturbed structure factor S(Q) of polymer 2
        nexpcc            =   modeex   ! number of exp-functions to describe background
        nexp1             =   modeex   ! number of exp-functions to describe component1
        nexp2             =   modeex   ! number of exp-functions to describe component2
        aexp_cc(1:nexpcc) =   aexpcc(1:nexpcc)   ! amplitude coeffs for laplace-exp representation of "matrix"
        rexp_cc(1:nexpcc) =   rexpcc(1:nexpcc)   ! rate      coeffs for laplace-exp representation of "matrix"
        aexp_s1(1:nexp1)  =   aexp11(1:nexp1)    ! amplitude coeffs for laplace-exp representation of polymer 1
        rexp_s1(1:nexp1)  =   rexp11(1:nexp1)    ! rate      coeffs for laplace-exp representation of polymer 1
        aexp_s2(1:nexp2)  =   aexp22(1:nexp2)    ! amplitude coeffs for laplace-exp representation of polymer 2
        rexp_s2(1:nexp2)  =   rexp22(1:nexp2)    ! rate      coeffs for laplace-exp representation of polymer 2
    
        if(analytic >= 1) then
 !! AUS NUMERISCHEN GRUENDEN DARF Phi1 bzw Phi2 nicht NULL sein
 !! erstmal default Aktion    
          if(phi1 < 1d-5) phi1 = 1d-5
          if(phi2 < 1d-5) phi2 = 1d-5
          call get_invlaplace_coeffs( il_coeffs11, il_coeffs12, il_coeffs22 , il_alphas, nnsum)  
        endif

 !!! testouptput
 !!  write(6,'(a,d14.7)') "Scc00          = ", Scc00  
 !!  write(6,'(a,d14.7)') "S0011          = ", S0011  
 !!  write(6,'(a,d14.7)') "S0022          = ", S0022  
 !!  write(6,'(a,d14.7)') "phi1           = ", phi1   
 !!  write(6,'(a,d14.7)') "phi2           = ", phi2   
 !!  write(6,'(a,i8   )') "nexp1          = ", nexp1  
 !!  write(6,'(a,i8   )') "nexp2          = ", nexp2  
 !!  write(6,'(a,i8   )') "nexpcc         = ", nexpcc 
 !!  write(6,'(a,10(d14.7,","))') "aexp_s1(1:nexp1)  = [", aexp_s1(1:nexp1)  
 !!  write(6,'(a,10(d14.7,","))') "aexp_s2(1:nexp2)  = [", aexp_s2(1:nexp2)  
 !!  write(6,'(a,10(d14.7,","))') "aexp_cc(1:nexpcc) = [", aexp_cc(1:nexpcc) 
 !!  write(6,'(a,10(d14.7,","))') "rexp_s1(1:nexp1)  = [", rexp_s1(1:nexp1)  
 !!  write(6,'(a,10(d14.7,","))') "rexp_s2(1:nexp2)  = [", rexp_s2(1:nexp2)  
 !!  write(6,'(a,10(d14.7,","))') "rexp_cc(1:nexpcc) = [", rexp_cc(1:nexpcc) 


       
   

     endif ilr

   

! RPA VARIABLE TRANSFER


!!  xil         =  0.0001d0    ! distance of path for inv-laplace integration
!!  epap        =  1d-9
!!  epsrpa      =  1d-5

 

    if( analytic >= 1 ) then  
!     try the polynomial approach
      call compute_invlaplace(  t, il_coeffs11, il_coeffs12, il_coeffs22 , il_alphas, nnsum, ss11, ss12, ss22) 
      call compute_invlaplace(0d0, il_coeffs11, il_coeffs12, il_coeffs22 , il_alphas, nnsum, ss110, ss120, ss220) 

     else

  write(6,*) "St_rpa: "

       if(alin .ne. 0d0) then
  write(6,*)"s11"
                                ss11  = St_rpa(t ,1,1)
  write(6,*)"1:", ss11, t, t0

           if(newcomp_required) ss110 = St_rpa(t0,1,1)
  write(6,*) ss11, ss110
        endif
        if(arous .ne. 0d0) then
  write(6,*)"s22"
           ss22                       = St_rpa(t ,2,2)
           if(newcomp_required) ss220 = St_rpa(t0,2,2)
  write(6,*) ss22, ss220
        endif
        if( alin*arous .ne. 0d0 ) then
  write(6,*)"s12"
           ss12                       = St_rpa(t ,1,2) 
           if(newcomp_required) ss120 = St_rpa(t0,1,2) 
  write(6,*) ss12, ss120

        endif 

!
     endif


     sqt  = ss11  * alin*alin + 2*alin*arous * ss12  +  ss22  * arous*arous 
     sqt0 = ss110 * alin*alin + 2*alin*arous * ss120 +  ss220 * arous*arous 


 select case(mode)  ! be sure that in the beginning the assignment of t, q are properly made (if i1:)
   case (0)
       th_rpalin = Sqt/Sqt0

  case (1) ! linear contribution at tau as function of q
     th_rpalin = plin

  case(2)  ! star contribution at tau  as function of q
     th_rpalin = prous
 
  case(3)  ! rpa contribution at tau as function of q
      th_rpalin = Sqt

  case default
     write(6,*)'invalid mode !'
      th_rpalin = 0

 end select

  call parset('nlin    ',(1.0*nlin),iadda) 
  call parset('nlin_cc ',(1.0*nlin_cc),iadda) 
  call parset('nrous   ',(1.0*nrous),iadda) 
  call parset('f       ',(1.0*f),iadda) 

  call parset('l       ',sngl(l),iadda) 
  call parset('philin  ',sngl(philin),iadda) 
  call parset('phirous ',sngl(phirous),iadda) 
  call parset('alin    ',sngl(alin),iadda) 
  call parset('arous   ',sngl(arous),iadda) 


  call parset('dss     ',sngl(dss),iadda) 
  call parset('npp_plus',(1.0*npp_plus),iadda) 
  call parset('sfak0   ',sngl(sfak0),iadda) 

  call parset('xil     ',sngl(xil),iadda) 
  call parset('epap    ',sngl(epap),iadda) 
  call parset('epsrpa  ',sngl(epsrpa),iadda) 


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 end function th_rpalin


