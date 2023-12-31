 FUNCTION th_rpatest(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  applies rpa mixture tor the scattering function for a star-linear polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))
! 
      use theory_description 
      use lmfit
      use rpa_laplace
      implicit none
!! 
      SAVE
!!    ====


      real    :: th_rpatest
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
     double precision :: wl4star    ! rouse rate star 
     integer          :: narmeff    ! effective segment star arm 
     double precision :: diffstar   ! centre-of-mass diffusino star
     double precision :: locr2_b    ! locrep2 parameter b 
     double precision :: locr2_a    ! locrep2 parameter a
     double precision :: locr2_ta   ! locrep2 parameter tau 
     double precision :: locr2_lz   ! locrep2 parameter lz 
     double precision :: locr2_te   ! locrep2 parameter teps 
     integer          :: nro_me     ! effective segments for Me-Modelling 
     double precision :: re_me      ! Re of Me (entangle strand Re)
     double precision :: wl4        ! rouse rate
! the recin parameter representation 
     double precision :: q          ! scattering wavevector
     double precision :: temp       ! temperature
     integer          :: narm       ! number of segments per star arm
     double precision :: l          ! segment length
     integer          :: f          ! number of arms per star
     double precision :: phistar    ! volume fraction of star
     integer          :: nlin       ! number of segments of linear polymer
     double precision :: philin     ! volume fraction of linera polymer
     double precision :: betadif
     double precision :: plimit     ! limit of star modes (internal EV value) which are frozen

     double precision :: alin       ! scattering contrast of linear polymer omponent 
     double precision :: astar      ! scattering contrast of start polymer

     integer          :: nxpoints   ! number of points in sum(exp) modelling of S00 functions 

     integer          :: mode       ! mode select 
     integer          :: modeex     ! mode select exp representation


     double precision, parameter :: tmin = 0.01d0
     double precision            :: t0  ! virtual 0 > 0 to allow the same numerics as for finite t to compute s(q,0)



     double precision :: tau        ! for modes > 0
 
     double precision   :: t

     double precision   :: sqt, sqt0
     double precision   :: locrep2
     double precision   :: pstar0, pstar
     double precision   :: plin0, plin
     double precision   :: dr, wx, lx
     double precision   :: Re_arm
     double precision   :: Sinv0, Sinv
     double precision   :: local_reptation2
     double precision   :: pericostar_sqt3
 

     integer            :: ifix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     integer, parameter      :: mparams=30
     double precision        :: last_params(mparams) = 0d0
     double precision        :: allparams(mparams)   = 0d0
     integer, parameter      :: mexp=10
     double precision        :: aexp11(mexp), rexp11(mexp)
     integer                 :: nex11
     double precision        :: aexp22(mexp), rexp22(mexp)
     integer                 :: nex12
     double precision        :: aexpcc(mexp), rexpcc(mexp)
     integer                 :: nexcc

     integer, parameter      :: mxpoints = 30
     double precision        :: t_samples(mxpoints)
     double precision        :: s_samples(mxpoints)
     double precision        :: tmax = 1000d0
     logical                 :: newcomp_required
     integer                 :: i
     double precision        :: ts, ssq

     double precision        :: ss11, ss110, ss12, ss120, ss22, ss220

     double precision        :: a1, a2, r1, r2, r3, b1, b2, g1, g2, g3
     integer                 :: analytic = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! ----- initialisation ----- 
    IF (ini.eq.0) then 
       thnam = 'rpatest'
       nparx =       18
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpatest = 0
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
        parnam ( 2) = 'wl4star '  ! rouse rate star
        parnam ( 3) = 'narmeff '  ! effective segment star arm
        parnam ( 4) = 'diffstar'  ! centre-of-mass diffusino star
        parnam ( 5) = 'locr2_b '  ! locrep2 parameter b
        parnam ( 6) = 'locr2_a '  ! locrep2 parameter a
        parnam ( 7) = 'locr2_ta'  ! locrep2 parameter tau
        parnam ( 8) = 'locr2_lz'  ! locrep2 parameter lz
        parnam ( 9) = 'locr2_te'  ! locrep2 parameter teps
        parnam (10) = 'nro_me  '  ! effective segments for Me-Modelling
        parnam (11) = 're_me   '  ! Re of Me (entangled strand Re)
        parnam (12) = 'wl4     '  ! Wl4 (local entanglement strand Re)
        parnam (13) = 'betadif '  !
        parnam (14) = 'plimit  '  ! mode ev limit for star contribution
        parnam (15) = 'mode    '  ! what to compute
        parnam (16) = 'modeex  '  ! selects way of representation in terms of nexp
        parnam (17) = 'nxpoints'  ! number of points to determine representation in terms of nexp
        parnam (18) = 'tzero   '  ! preliminary: effective zero time for S(Q,t=0=tzero)

! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate star" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment star arm" !//cr//parspace//&
        th_param_desc( 4,idesc) = "centre-of-mass diffusino star" !//cr//parspace//&
        th_param_desc( 5,idesc) = "locrep2 parameter b" !//cr//parspace//&
        th_param_desc( 6,idesc) = "locrep2 parameter a" !//cr//parspace//&
        th_param_desc( 7,idesc) = "locrep2 parameter tau" !//cr//parspace//&
        th_param_desc( 8,idesc) = "locrep2 parameter lz" !//cr//parspace//&
        th_param_desc( 9,idesc) = "locrep2 parameter teps" !//cr//parspace//&
        th_param_desc(10,idesc) = "effective segments for Me-Modelling" !//cr//parspace//&
        th_param_desc(11,idesc) = "Re of Me (entangle strand Re)" !//cr//parspace//&
        th_param_desc(12,idesc) = "Wl4 entaglement strand" !//cr//parspace//&
        th_param_desc(13,idesc) = "betadiff" !//cr//parspace//&
        th_param_desc(14,idesc) = "plimit EV limit below modes are frozen           " !//cr//parspace//&
        th_param_desc(15,idesc) = "mode: 0 normal, >0 S(x=Q,tau) 1 lin, 2 star 3 rpa" !//cr//parspace//&
        th_param_desc(16,idesc) = "mode: number of exps to represent" !//cr//parspace//&
        th_param_desc(17,idesc) = "number of points to determine representation in terms of nexp" !//cr//parspace//&
        th_param_desc(18,idesc) = "preliminary: effective zero time for S(Q,t=0=tzero)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  3,idesc) = "narm     > number of segments per star arm"
        th_file_param(  4,idesc) = "l        > segment length"
        th_file_param(  5,idesc) = "f        > number of arms per star"
        th_file_param(  6,idesc) = "phistar  > volume fraction of star"
        th_file_param(  7,idesc) = "nlin     > number of segments of linear polymer"
        th_file_param(  8,idesc) = "philin   > volume fraction of linera polymer"
        th_file_param(  9,idesc) = "tau      > tau val for mode > 0 calc"
        th_file_param( 10,idesc) = "alin     > scattering contrast linear"
        th_file_param( 11,idesc) = "astar    > scattering contrast star"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpatest = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4star  = abs( pa( 2))
      narmeff  = nint(pa( 3))
      diffstar = abs( pa( 4))
      locr2_b  = abs( pa( 5))
      locr2_a  = abs( pa( 6))
      locr2_ta = abs( pa( 7))
      locr2_lz = abs( pa( 8))
      locr2_te = abs( pa( 9))
      nro_me   = nint(pa(10))
      re_me    = abs( pa(11))
      wl4      = abs( pa(12))
      betadif  = abs( pa(13))
      plimit   = abs( pa(14))
      mode     = nint(pa(15))
      modeex   = nint(pa(16))
      nxpoints = nint(pa(17))
      t0       = abs( pa(18))


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
      call parget('narm    ',xh,iadda,ier)
      narm     = xh
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
      call parget('phistar ',xh,iadda,ier)
      phistar  = xh
! >>> extract: number of segments of linear polymer
      xh = 10
      call parget('nlin    ',xh,iadda,ier)
      nlin     = nint(xh)
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
      call parget('astar    ',xh,iadda,ier)
      astar  = xh
!!! >>> tentative onl testing
      xh = 0d0
      call parget('analytic ',xh,iadda,ier)
      analytic  = nint(xh) 
!!! 
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


      allparams(1:22) = [ wl4star, dble(narmeff), diffstar, locr2_b, locr2_a, locr2_ta, locr2_lz, locr2_te, &
                         dble(nro_me) , re_me,   wl4     , betadif, plimit , q , dble(narm), l, dble(f),   &
                         phistar, dble(nlin), philin, alin, astar ]

      newcomp_required = ( sum(abs(allparams-last_params)) .ne. 0 ) .and. (modeex > 0)
    
! newcomp_required = .true.
!!!!!!
!!!!! Aus bisher unbekannter Ursache reagiert die Rechnung (analytic vs numeric) auf diesen Flag
!!!!! eigentlich sollte alles einmal berechnet so bleiben
!!!!! was ist da los, der Effekt is zwar klein aber sichtbar !!! ????
!!!!! parameter vergessen ??
!!!!!!

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
      do i=1,nxpoints
             ts      =  exp(i*log(tmax*300d0)/nxpoints)/100d0
             sqt0    =  local_reptation2( q*locr2_a, locr2_te   , locr2_b, locr2_lz)
             sqt     =  local_reptation2( q*locr2_a, ts/locr2_ta, locr2_b, locr2_lz)
             locrep2 =  sqt/sqt0
             dr      = 1d-20
             ifix    = 0
             call  NrouseY(q,ts,temp,Dr,wl4,nro_me,re_me, Wx, lx,ifix, sqt0,sqt)
             plin0   = nlin  * Debye_qnl(q, nlin, l) 
             plin    = plin0 * locrep2 * sqt / sqt0
             t_samples(i) = ts
             s_samples(i) = locrep2 * sqt / sqt0
          enddo
          call nexp_match(t_samples,s_samples,nxpoints,modeex,aexp11,rexp11,ssq)
          aexpcc = aexp11
          rexpcc = rexp11
          if(ssq > 1d-4) then
            write(6,*)"rpa_test exp model bad match 11", ssq
          endif
     elseif(modeex==0) then ! diese Berechnung und die analogen jeweils in einen separate Prozedur packen ! --> F(q,t), S(q)
       sqt0    =  local_reptation2( q*locr2_a, locr2_te  , locr2_b, locr2_lz)
       sqt     =  local_reptation2( q*locr2_a, t/locr2_ta, locr2_b, locr2_lz)
       locrep2 =  sqt/sqt0
       dr      = 1d-20
       ifix    = 0
       call  NrouseY(q,t,temp,Dr,wl4,nro_me,re_me, Wx, lx,ifix, sqt0,sqt)
       plin0   = nlin  * Debye_qnl(q, nlin, l) 
       plin    = plin0 * locrep2 * sqt / sqt0
     endif ilr


 ! star
ist: if( newcomp_required ) then
      do i=1,nxpoints
             ts     =  exp(i*log(tmax*300d0)/nxpoints)/100d0
             Re_arm = sqrt(narm * l**2)
             sqt0   =  pericostar_sqt3(q,0d0,f,narmeff,Re_arm,Wl4star,plimit)    ! automatic diffusion if diffstar==0 
             sqt    =  pericostar_sqt3(q,ts ,f,narmeff,Re_arm,Wl4star,plimit)
 ! fix normalization
             pstar0 =  SQstar(q,narm,f,l)
             pstar  =  sqt/sqt0 * pstar0
             pstar  =  pstar * exp( -diffstar * q*q * (ts)**betadif )            !! add diffusion expilicitly
             t_samples(i) = ts
             s_samples(i) = sqt/sqt0 *  exp( -diffstar * q*q * (ts)**betadif )    
           enddo
          call nexp_match(t_samples,s_samples,nxpoints,modeex,aexp22,rexp22,ssq)
          if(ssq > 1d-4) then
            write(6,*)"rpa_test exp model bad match 22", ssq
          endif
     elseif(modeex==0) then     

        Re_arm = sqrt(narm * l**2)
        sqt0   =  pericostar_sqt3(q,0d0,f,narmeff,Re_arm,Wl4star,plimit)    ! automatic diffusion if diffstar==0 
        sqt    =  pericostar_sqt3(q,t  ,f,narmeff,Re_arm,Wl4star,plimit)
 ! fix normalization
        pstar0 =  SQstar(q,narm,f,l)
        pstar  =  sqt/sqt0 * pstar0
        pstar  =  pstar * exp( -diffstar * q*q * (t)**betadif )                       !! add diffusion expilicitly

     endif ist

   
!!     if( modeex > 0 ) then
!!       plin   = 0
!!       pstar  = 0
!!       do i=1,modeex
!!         plin  = plin  + aexp11(i) * exp( -t * rexp11(i) ) 
!!         pstar = pstar + aexp22(i) * exp( -t * rexp22(i) ) 
!!       enddo
!!       plin  = plin  * plin0
!!       pstar = pstar * pstar0
!!     endif

! RPA VARIABLE TRANSFER


 xil         =  0.0001d0    ! distance of path for inv-laplace integration
 epap        =  1d-9
 epsrpa      =  1d-5

 
!
! componente c    = lineraes polymer (matrix) gleiche Parameter wie comp 1
! componente 1    = lineares polymer (markierter Teil)
! componente 2    = start-polymer
     phi1              =   philin   ! volume fraction of polymer component 1
     phi2              =   phistar  ! volume fraction of polymer component 2
     Scc00             =   plin0    ! unperturbed structure factor S(Q) of "matrix" polymers
     S0011             =   plin0    ! unperturbed structure factor S(Q) of polymer 1
     S0022             =   pstar0   ! unperturbed structure factor S(Q) of polymer 2
     nexpcc            =   modeex   ! number of exp-functions to describe background
     nexp1             =   modeex   ! number of exp-functions to describe component1
     nexp2             =   modeex   ! number of exp-functions to describe component2
     aexp_cc(1:nexpcc) =   aexpcc(1:nexpcc)    ! amplitude coeffs for laplace-exp representation of "matrix"
     rexp_cc(1:nexpcc) =   rexpcc(1:nexpcc)   ! rate      coeffs for laplace-exp representation of "matrix"
     aexp_s1(1:nexp1)  =   aexp11(1:nexp1)    ! amplitude coeffs for laplace-exp representation of polymer 1
     rexp_s1(1:nexp1)  =   rexp11(1:nexp1)    ! rate      coeffs for laplace-exp representation of polymer 1
     aexp_s2(1:nexp2)  =   aexp22(1:nexp2)    ! amplitude coeffs for laplace-exp representation of polymer 2
     rexp_s2(1:nexp2)  =   rexp22(1:nexp2)    ! rate      coeffs for laplace-exp representation of polymer 2
 

    if( analytic == 1 .and. modeex == 3 ) then  
       a1 = aexp_s2(1)
       a2 = aexp_s2(2)
       r1 = rexp_s2(1)
       r2 = rexp_s2(2)
       r3 = rexp_s2(3)

       b1 = aexp_cc(1)
       b2 = aexp_cc(2)
       g1 = rexp_cc(1)
       g2 = rexp_cc(2)
       g3 = rexp_cc(3)

!! we have to exchange the roles of 11 and 22 since the analytic function imply S22=Scc
       if(alin .ne. 0d0) then
         ss11   = InvLaplace2d22(t,   S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
         if(newcomp_required) ss110  = InvLaplace2d22(0d0, S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
write(6,*)"i11:",t,t0,ss11,ss110
       endif
       if(astar .ne. 0d0) then
         ss22   = InvLaplace2d11(t,   S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
         if(newcomp_required) ss220  = InvLaplace2d11(0d0, S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
write(6,*)"i22:",t,t0,ss22,ss220

       endif
        if( alin*astar .ne. 0d0 ) then
         ss12   = InvLaplace2d12(t,   S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
         if(newcomp_required) ss120  = InvLaplace2d12(0d0, S0022, a1, a2, r1, r2, r3, S0011, Scc00, b1, b2, g1, g2, g3, phi2, phi1 )
write(6,*)"i12:",t,t0,ss12,ss120
        endif 

     else

       if(alin .ne. 0d0) then
                                ss11  = St_rpa(t ,1,1)
           if(newcomp_required) ss110 = St_rpa(t0,1,1)
        endif
        if(astar .ne. 0d0) then
           ss22                       = St_rpa(t ,2,2)
           if(newcomp_required) ss220 = St_rpa(t0,2,2)
        endif
        if( alin*astar .ne. 0d0 ) then
           ss12                       = St_rpa(t ,1,2) 
           if(newcomp_required) ss120 = St_rpa(t0,1,2) 
        endif 

!
     endif
!!     sqt  = St_rpa(t ,1,1)*alin*alin + 2*alin*astar*St_rpa(t ,1,2) + St_rpa(t ,2,2)*astar*astar 
!!     sqt0 = St_rpa(t0,1,1)*alin*alin + 2*alin*astar*St_rpa(t0,1,2) + St_rpa(t0,2,2)*astar*astar 



     sqt  = ss11  * alin*alin + 2*alin*astar * ss12  +  ss22  * astar*astar 
     sqt0 = ss110 * alin*alin + 2*alin*astar * ss120 +  ss220 * astar*astar 


! Rpa (assuming contrast only between h-lin polymer and the rest DAS MUSS NOCH GEAENDER WERDEN !

!      Sqt0   = philin*((phistar+philin-1)*Plin0-Pstar0*phistar)*Plin0/(Plin0*(phistar-1)-Pstar0*phistar)
!      Sqt    = philin*((phistar+philin-1)*Plin -Pstar *phistar)*Plin /(Plin *(phistar-1)-Pstar *phistar)




 select case(mode)  ! be sure that in the beginning the assignment of t, q are properly made (if i1:)
   case (0)
       th_rpatest = Sqt/Sqt0

  case (1) ! linear contribution at tau as function of q
     th_rpatest = plin

  case(2)  ! star contribution at tau  as function of q
     th_rpatest = pstar
 
  case(3)  ! rpa contribution at tau as function of q
      th_rpatest = Sqt

  case default
     write(6,*)'invalid mode !'
      th_rpatest = 0

 end select

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

      function SQstar( q, narm, f, l) result(p0)
        implicit none
        double precision, intent(in)   :: q
        integer         , intent(in)   :: narm
        integer         , intent(in)   :: f
        double precision, intent(in)   :: l

        double precision               :: p0

        double precision               :: x

        x   = (q*l)**2 * (narm / 6d0)
        p0  = narm * (2d0/x**2)
        p0  = p0 * ( x - (1d0-exp(-x) ) + (f-1)/2d0 * (1d0-exp(-x))**2) 

      end function  SQstar





 end function th_rpatest


 
double precision function pericostar_sqt3(q,t,f_arm,n_arm,Re_arm,Wl4,plim)
!-------------------------------------------------------------------------
!! star following the derivation of:
!!      M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!
!! with experimental mode suppression
!!
!! this module written by: michael Monkenbusch, Jcns
!!

  implicit none

  double precision, intent(in)  :: q
  double precision, intent(in)  :: t
  integer,          intent(in)  :: f_arm
  integer,          intent(in)  :: n_arm 
  double precision, intent(in)  :: Re_arm
  double precision, intent(in)  :: Wl4
  double precision, intent(in)  :: plim


  double precision, parameter   :: Pi=3.141592654d0
  integer                       :: N
  integer                       :: i, j, ip, l, i1, i2, job, info
  double precision              :: l0         ! effective segment length
  double precision              :: tauR       ! basic Rouse-time of a chain of 2 arm lengths

  double precision              :: p, sum, s11, s12, dec
  double precision              :: rmodifier

  integer, parameter   :: real_kind = 4
 
  real(kind=real_kind)         :: U(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: H(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: X(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: M(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: A(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: phi(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: work((1+n_arm*f_arm)*2)
  real(kind=real_kind), allocatable, save   :: Eigenvalues(:)
  real(kind=real_kind), allocatable, save   :: Eigenvectors(:,:)

  logical, save                :: first_run = .true.
  integer, save                :: f_arm0, n_arm0 




  N      = f_arm * n_arm + 1
  l0     = Re_arm/sqrt(dble(n_arm))      ! must be changed if a nontrivial U-Matrix will be used

!! compute Eigenvectors and Eigenvalues for the problem:
!! see: M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!

    N=1+n_arm*f_arm
!
! the first index pertains to the center of the star
!


 first_run = .true.    ! allocated are not save in all cases thus we need this as a workaround 


eig:  if(first_run .or. n_arm .ne. n_arm0 .or. f_arm .ne. f_arm0) then

    first_run = .false.
    f_arm0    = f_arm
    n_arm0    = n_arm

    if(allocated(Eigenvalues)) then
      deallocate(Eigenvalues)
      deallocate(Eigenvectors)
    endif

    allocate(Eigenvalues(N))
    allocate(Eigenvectors(N,N))

    U = 0            ! inverse of static bond correlation matrix <l_i*l_j>/l**2  (here  1)
    M = 0            ! bead-to-bond vector transformation matrix
    A = 0            ! A = Mt(U)M   (were the first col and row of U are 0)
    H = 0            ! H = hydrodynamic matrix (so far here 1) friction is lumped in the Wl4 parameter
 
    forall (i=2:N)   U(i,i) = 1.0
    forall (i=1:N)   M(1,i) = 1.0/N
    forall (i=2:N)   M(i,i) = 1.0
    forall (i=1:N)   H(i,i) = 1.0
    do j=1,f_arm
      i1 = 2+n_arm*(j-1)
      i2 = j*n_arm
      forall (i=i1:i2)   M(i+1,i) = -1_real_kind
    enddo
    do j=1,f_arm
      i = 2+(j-1)*n_arm
      M(i,1) = -1
    enddo
 
   !
   ! compute A = M^t (U) M
   !
   do i=1,N
    do j=1,N
     X(i,j) = dot_product(U(i,1:N),M(1:N,j))
    enddo
   enddo
   do i=1,N
    do j=1,N
     A(i,j) = dot_product(M(1:N,i),X(1:N,j))
    enddo
   enddo
 
   ! and for completeness HA !
   do i=1,N
    do j=1,N
     X(i,j) = dot_product(H(i,1:N),A(1:N,j))
    enddo
   enddo

   ! check symmetry
   do i=1,N
    do j=i+1,N
     if(abs(X(i,j)-X(j,i)) > 1e-5) then
       write(6,*)'ERROR nonsymmetric matrix HA '
       write(6,*)'USE star2 to deal with cases where H is nontrivial'
       pericostar_sqt3 = 0
       return
     endif
    enddo
   enddo

 
   Eigenvectors = X
   job = 1            ! signal to ssiev that the matrix is to be replaced by the eigenvectors
   !
   ! Slatec:
    call SSIEV (Eigenvectors, N , N, Eigenvalues, WORK, JOB, INFO)

 !!   write(6,*)"Eigenvalues:"
 !!   do i=1,N
 !!    write(6,'(i8,f12.6)')i, Eigenvalues(i)
 !!   enddo
 
 endif eig



! -- fill the (symmetric) phi(n,m) matrix -- ( the L**2 * d_ij(t) values of the papaers) 
philp: do i=1,N
        do j=1,N 
          sum = 0d0
dl:       do l = 2 ,N  ! sjkip the first eigenvalue = 0 
           p = Eigenvalues(l) 

           if( p < plim ) then
               rmodifier = 0d0
           else
               rmodifier = 1d0
           endif

           sum = sum + (  Eigenvectors(i,l)**2+Eigenvectors(j,l)**2- &
                        2*Eigenvectors(i,l)*Eigenvectors(j,l)*exp(- rmodifier * wl4/(l0**4)*p*t))/p
          enddo dl
          phi(j,i) = (l0**2) * sum 
        enddo
       enddo philp
!
! --> ok phi matrix filled (check p-summation role, fixed middle point --> p odd)
!

!
! -- and now the scattering contibutions of the arms
!

        s11 = 0         ! the "self term"
        do i=2,n_arm+1
           s11 = s11 + 2 * exp(-(q*q/6d0) * phi(1,i)) + exp(-(q*q/6d0)* phi(i,i)) 
        enddo
        do i=2,n_arm
         do j=i+1,n_arm+1
           s11 = s11 + 2* exp(-(q*q/6d0)* phi(i,j))
         enddo
        enddo

        s12 = 0         ! the "cross term"
        if(f_arm > 1) then
          do i=2,n_arm+1
             s12 = s12 + exp(-(q*q/6d0) * phi(i+n_arm,i))
          enddo
          do i=2,n_arm
             do j=i+1,n_arm+1
               s12 = s12 + 2 * exp(-(q*q/6d0) * phi(i+n_arm,j))
             enddo
          enddo
        endif

        sum = (1+f_arm*(s11+(f_arm-1)*s12))/N**2

     pericostar_sqt3 = sum


!!! tentatively add diffusion !
!!    if(diff == 0d0) then
!!       pericostar_sqt3 = pericostar_sqt3 * exp( -q*q* wl4/(3*f_arm*n_arm*l0*l0) * t)
!!    endif
!!


end function pericostar_sqt3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
