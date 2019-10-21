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
     double precision :: nseg       ! nseg prefactor                                                              
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
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sinhafra = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " fractal scattering assuming a massfractal with dimension df and correlation length xi "//cr//&
       " normalized such that P(q=0)=1 and final form "//cr//&
       " S(q) = ampli * 1/(1/(nseg*P(q) + va2)   * MolFromfactor "

       th_citation(idesc)     = " J.K. Kjems, T. Freltoft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'nseg    '  ! number segments (prefactor) ...  
        parnam ( 3) = 'df      '  ! fractal dimension                                                               
        parnam ( 4) = 'xi      '  ! correlation length                                                              
        parnam ( 5) = 'va2     '  ! a2 type                                                                 
        parnam ( 6) = 'withsq  '  ! multiplied by molecular formfactor as in sq_molecule_dendirmer.f90                                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "nseg as prefactor to P(q) for comparability of a, va2 parameters with ndebye " !//cr//parspace//&
        th_param_desc( 3,idesc) = "fractal dimension" !//cr//parspace//&
        th_param_desc( 4,idesc) = "correlation length" !//cr//parspace//&
        th_param_desc( 5,idesc) = "a2 type correction parameter" !//cr//parspace//&
        th_param_desc( 6,idesc) = " selector (0/1) multiply with molecular sq as in sq_molecule_dendirmer.f90 "                                                                
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
      nseg     =      abs( pa( 2) )
      df       =      abs( pa( 3) )
      xi       =      abs( pa( 4) )
      va2      =      pa( 5) 
      withsq   =     nint( pa( 6) )
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

     th_sinhafra = 1d0/(-1+df)*sqrt(q**2 *xi**2 +1) * sin( (-1+df)*atan(q*xi)) / ( (q**2 * xi**2+1)**(0.5d0*df) *q*xi )

     th_sinhafra = nseg * th_sinhafra

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
 end function th_sinhafra
