 FUNCTION th_sinhafra(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  fractal scattering assuming a massfractal with dimension df and correlation length xi direct integration for r0 .ne.0
!  J.K. Kjems, T. Freltoft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_sinhafra
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor                                                                       
     double precision :: c          ! sq = 1 + c* ...                                                               
     double precision :: df         ! fractal dimension                                                               
     double precision :: xi         ! correlation length                                                              
     double precision :: r0         ! correlation hole   ??                                                              
     double precision :: va2        ! a2 type ww parameter                                                               
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision   :: adapint, eps=1d-9, erracc
     integer            :: maxit = 1000
     integer            :: nper  = 20
     double precision   :: q
     double precision   :: acc = 1d-9
     double precision   :: g1, g2, gm, dg
     double precision   :: sq_molecule_dendrimer
     integer            :: withsq
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sinhafra'
       nparx =        7
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sinhafra = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " fractal scattering assuming a massfractal with dimension df and correlation length xi direct integration for r0 .ne.0"
       th_citation(idesc)     = " J.K. Kjems, T. Freltoft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'cf      '  ! sq = 1 +c * ...                                                              
        parnam ( 3) = 'df      '  ! fractal dimension                                                               
        parnam ( 4) = 'xi      '  ! correlation length                                                              
        parnam ( 5) = 'r0      '  ! correlation hole                                                                
        parnam ( 6) = 'va2     '  ! a2 type                                                                 
        parnam ( 7) = 'withsq  '  ! multiplied by molecular formfactor as in sq_molecule_dendirmer.f90                                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "factor c in S(Q) = 1 +c * fractal expression" !//cr//parspace//&
        th_param_desc( 3,idesc) = "fractal dimension" !//cr//parspace//&
        th_param_desc( 4,idesc) = "correlation length" !//cr//parspace//&
        th_param_desc( 5,idesc) = "correlation hole: caveat r0 .ne. 0 --> neg. S(Q) needs fix!" !//cr//parspace//&
        th_param_desc( 6,idesc) = "a2 type correction parameter" !//cr//parspace//&
        th_param_desc( 7,idesc) = " selector (0/1) multiply with molecular sq as in sq_molecule_dendirmer.f90 "                                                                
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_sinhafra = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      c        =      abs( pa( 2) )
      df       =      abs( pa( 3) )
      xi       =      abs( pa( 4) )
      r0       =      abs( pa( 5) )
      va2      =      abs( pa( 6) )
      withsq   =     nint( pa( 7) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
!      call parget('xh      ',xh,iadda,ier)
!      xh         = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x

     if(r0 == 0d0) then
         th_sinhafra = GAMMA(-1+df)*sqrt(q**2 *xi**2 +1) * sin( (-1+df)*atan(q*xi)) / ( (q**2 * xi**2+1)**(0.5d0*df) *q*xi )
     else
          th_sinhafra = adapint(gr_kernel_innen,0d0,2*r0, eps, maxit, erracc)
          g1 = r0*2
          gm = -log(acc) * xi
          dg = nper*2*pi/q
          if(dg > gm) dg = gm
 d1:      do 
            g2 = g1+dg
            th_sinhafra =  th_sinhafra + adapint(gr_kernel,g1,g2, eps, maxit, erracc) 
            if(g2 > gm) exit d1
            g1 = g2
          enddo d1
  
          q = 1d-4
          th_sinhafra = th_sinhafra / ((xi**df*GAMMA(-1d0+df)*(-1d0+df))  ) ! &
 !                                     + adapint(gr_kernel_innen,0d0,2*r0, eps, maxit, erracc)  ! normalisiert so dass S(Q) = 1 
 ! NOTE: die vorgeschlagene Realisierung des Korrelation hole wie hier durch r0 gegeben fuehrt bei groesserem q zu oszillationen ins Negative !!
 !       hier muss modifiziert werden !!
  
!
! !! DIESE FORMEL (folgende) (6) AUS DEM PAPER IST FEHLERHAFT midestens muss xi**df --> xi**df bzw xi**(df/2-1)
! !! aber auch dann ist das resultat unphysikalisch
!
!        th_sinhafra = xi**df * sqrt(1+q*q*xi*xi) / ( (q*xi) * (1+xi*xi*q*q)**(0.5d0*df) ) * &
!             Imagpart(cgamm_bf(df-1d0,2*r0/xi*cmplx(1d0,-q*xi))*exp(cmplx(0d0,(df-1d0)*atan(q*xi)))) 
!        
     endif

! write(*,*)"T1:",q,th_sinhafra

     th_sinhafra = 1d0 + c * th_sinhafra

     th_sinhafra = 1d0 / ( 1d0/th_sinhafra + va2)

     th_sinhafra = ampli * th_sinhafra

! write(*,*)"T2:",q,th_sinhafra, ampli, c, va2, r0

     if(withsq == 1) then
      th_sinhafra = th_sinhafra * sq_molecule_dendrimer(q)
     endif
! write(*,*)"T3:",q,th_sinhafra
 
! ---- writing computed parameters to the record >>>  
      rg = xi * sqrt(df*(1+df))
      call parset('rg      ',sngl(rg),iadda,ier)
 
 CONTAINS 
 
function cgamm_bf(a,z) result(val) 
! incomplete Gamma function with complex second argument z
! brute force integral THIS CAN BE MADE BETTER, but takes time...
implicit none
double precision, intent(in) :: a
 complex(kind=8), intent(in) :: z
 complex(kind=8) :: val
 complex(kind=8) :: y
 
 double precision :: u, du = 0.01d0
 integer          :: i, imax
 
 imax = min(nint(-log(1d-8)/(du*Realpart(z))), 100000)

 u   = 1d0 
 y = (z*u)**(a-1d0) * exp(-z*u)
 val = y/2d0
!$OMP PARALLEL DO REDUCTION(+:val) PRIVATE(u,y) 
 do i=1,imax
   u = i * du  + 1d0
   y = (z*u)**(a-1d0) * exp(-z*u)
   val = val + y
 enddo
!$OMP END PARALLEL DO

 val = val * du * z
 
end function cgamm_bf



! subroutines and functions entered here are private to this theory and share its variables 
 
  double precision function gr_kernel(r) 
    implicit none
    double precision, intent(in) :: r

    gr_kernel = sin(q*r) * r**(-2d0+df) * exp(-r/xi) / q

  end function gr_kernel

 
  double precision function gr_kernel_innen(r) 
    implicit none
    double precision, intent(in) :: r
    double precision             :: gri

!!! ??????????????  !!!

    gri = (1d0/3d0)*Pi*((1d0/8d0)*(2*r0-r)**2*(4*r0+r)+(1d0/2d0*(r0-(1/4)*r))*(2*r0+r)**2-(3*(r0**2-(1d0/12d0)*r**2))*r)
    gri = 0 !! gri ist Kugel-Kugel ueberlapp-funktion (FT gibt Kugelformfaktor Sp(q)
    !! Ansatz scheint nicht zu kunktionieren........
    !! so wie es im Paper sthet geht e auch ins negative, also...

    gr_kernel_innen = sin(q*r) * r**2 * (gri+(2*r0)**(-3d0+df)) * exp(-(2*r0)/xi) / q

  end function gr_kernel_innen

 end function th_sinhafra
