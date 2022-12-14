 FUNCTION th_debstar(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  fromfactor of an f-arm star with arms with Gaussian chain statistics
!  H. Benoit, J. Polym. Sci. (1953), 11, 507-510
      use theory_description 
      implicit none 
      real    :: th_debstar
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
     double precision :: f          ! number of arms                                                                  
     double precision :: rgstar     ! star radius of gyration                                                         
     double precision :: nulin      ! exponent of linear                                                              
     double precision :: rglin      ! radius of gyration linear polymer                                               
     double precision :: chipar     ! potential chi-parameter                                                         
! the recin parameter representation 
     double precision :: phistar    ! volume fraction of start                                                        
     double precision :: philin     ! volume fraction of linear                                                       
     double precision :: mwstar     ! molecular weight star polymer  (g/mol)                                          
     double precision :: mwlin      ! molecular weight linear polymer(g/mol)                                          
     double precision :: mwseg      ! molecular weight of a segment  (g/mol)                                          
     double precision :: rhostar    ! density star polymer                                                            
     double precision :: rholin     ! linear polymer (g/cm**2)                                                        
     double precision :: astar      ! scattering length density contrast star (cm**-2)                                         
     double precision :: alin       ! scattering length density contrast lin  (cm**-2)                                        
     integer          :: linstar    ! switch to enable mixture of lin/d stars 
     integer          :: fstar      ! number of arms 
 
 
     double precision :: th
 
     double precision, parameter :: Navogadro = 6.022140857d23
     double precision   :: q
     double precision   :: df, rg
     double precision   :: sqstar, sqlin, sq, s11, s22, s12, s33
     double precision   :: v11, v22, v12
     double precision   :: vseg
     double precision   :: nseg_star, nseg_lin

     double precision   :: sqstar0, sqlin0
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'debstar'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_debstar = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " RPA scattering intensity for a ternary system " //cr//&
                                " d-linear polymer (background) " //cr// &
                                " h/d-linear polymer (volume fraction philin) contrast (alin) " //cr//&
                                " h/d-star polymer (volume fraction phistar) contrast (astar) " //cr//&
                                " if linstar parameter =1, phlin/alin also refers to a star (Krücke) " //cr//&
                                " fit parameter f is only active if no fstar parameter is in the data "

       th_citation(idesc)     = " H. Benoit, J. Polym. Sci. (1953), 11, 507-510"//cr//parspace//&
                                " G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146"//cr//parspace//&
                                " H. Benoit, Compt. Rend. (1957) 245, 2244-2247 "//cr//parspace//&

                                " A.Z. Akcasu, M. Tombakoglu, Macromolecules 23, 607-612 (1990) RPA "//cr//parspace//&
                                " B. Hammouda, Journal Of Chemical Physics 98, 3439-3444 (1993) "

!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'f       '  ! number of arms                                                                  
        parnam ( 3) = 'rgstar  '  ! star radius of gyration                                                         
        parnam ( 4) = 'nulin   '  ! exponent of linear                                                              
        parnam ( 5) = 'rglin   '  ! radius of gyration linear polymer                                               
        parnam ( 6) = 'chipar  '  ! potential chi-parameter                                                         
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of arms" !//cr//parspace//&
        th_param_desc( 3,idesc) = "star radius of gyration" !//cr//parspace//&
        th_param_desc( 4,idesc) = "exponent of linear" !//cr//parspace//&
        th_param_desc( 5,idesc) = "radius of gyration linear polymer" !//cr//parspace//&
        th_param_desc( 6,idesc) = "potential chi-parameter" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "phistar  > volume fraction of star"
        th_file_param(  2,idesc) = "philin   > volume fraction of labeled linear"
        th_file_param(  3,idesc) = "mwstar   > molecular weight star polymer  (g/mol)"
        th_file_param(  4,idesc) = "mwlin    > molecular weight linear polymer(g/mol)"
        th_file_param(  5,idesc) = "mwseg    > molecular weight of a segment  (g/mol)"
        th_file_param(  6,idesc) = "rhostar  > density star polymer"
        th_file_param(  7,idesc) = "rholin   > linear polymer (g/cm**2)"
        th_file_param(  8,idesc) = "astar    > scattering length density contrast star (cm**-2)"
        th_file_param(  9,idesc) = "alin     > scattering length density contrast lin  (cm**-2)"
        th_file_param( 10,idesc) = "linstar  > if = 1 the labelled lin comp is star ???? check "
        th_file_param( 11,idesc) = "fstar    > number of arms"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "         > "
! 
        th_debstar = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      f        = abs( pa( 2) ) 
      rgstar   = abs( pa( 3) )
      nulin    = abs( pa( 4) )
      rglin    = abs( pa( 5) )
      chipar   =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: volume fraction of start
      xh = 0
      call parget('phistar ',xh,iadda,ier)
      phistar  = xh
!! >>> extract: volume fraction of start
      xh = 0
      call parget('philin  ',xh,iadda,ier)
      philin  = xh
! >>> extract: molecular weight star polymer  (g/mol)
      xh = 0
      call parget('mwstar  ',xh,iadda,ier)
      mwstar   = xh
! >>> extract: molecular weight linear polymer(g/mol)
      xh = 0
      call parget('mwlin   ',xh,iadda,ier)
      mwlin    = xh
! >>> extract: molecular weight of a segment  (g/mol)
      xh = 0
      call parget('mwseg   ',xh,iadda,ier)
      mwseg    = xh
! >>> extract: density star polymer
      xh = 1
      call parget('rhostar ',xh,iadda,ier)
      rhostar  = xh
! >>> extract: linear polymer (g/cm**2)
      xh = 1
      call parget('rholin  ',xh,iadda,ier)
      rholin   = xh
! >>> extract: scattering length density star (cm**-2)
      xh = 0
      call parget('astar   ',xh,iadda,ier)
      astar    = xh
! >>> extract: scattering length density lin  (cm**-2)
      xh = 0
      call parget('alin    ',xh,iadda,ier)
      alin     = xh
!! >>> extract: switch variable (lin--> star)
      xh = 0
      call parget('linstar ',xh,iadda,ier)
      linstar     = nint(xh)
!!! >>> extract: arm number
      xh = 0
      call parget('fstar   ',xh,iadda,ier)
      fstar     = (xh)
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x

     nseg_star = mwstar/mwseg
     nseg_lin  = mwlin/mwseg
     vseg      = mwseg / (rholin * Navogadro)

     if( fstar > 0d0 ) f = fstar

     sqstar0  =  debyestar(q, rgstar, f)
     sqlin0   =  benoitmf( q, rglin, nulin )

     if(linstar == 1 ) sqlin0 = sqstar0     ! dies ist eine Krücke, ???

     sqstar  = sqstar0 * nseg_star * phistar            * vseg
     sqlin   = sqlin0  * nseg_lin  * philin             * vseg

     s33     = sqlin0  * nseg_lin  * (1-phistar-philin)  * vseg

     v11     = 1/s33 - 2*chipar/vseg
     v22     = 1/s33 - 2*chipar/vseg
     v12     = 1/s33 -   chipar/vseg

     s11     = sqstar * (1+ v22 * sqlin)  / ((1+v11*sqstar)*(1+v22*sqlin)-v12**2 *sqstar*sqlin)
     s22     = sqlin  * (1+ v11 * sqstar) / ((1+v11*sqstar)*(1+v22*sqlin)-v12**2 *sqstar*sqlin)
     s12     = -sqlin  * v11 * sqstar     / ((1+v11*sqstar)*(1+v22*sqlin)-v12**2 *sqstar*sqlin)
    


     sq = astar**2 * s11 + alin**2 * s22 + 2*astar*alin * s12

     th = ampli * sq


     th_debstar = th

!! write(6,*)"===================================================================="
!! write(6,*)" q         = ",   q
!! write(6,*)" nseg_star = ",   nseg_star 
!! write(6,*)" nseg_lin  = ",   nseg_lin 
!! write(6,*)" astar     = ",   astar 
!! write(6,*)" alin      = ",   alin 
!! write(6,*)" vseg      = ",   vseg
!! write(6,*)" f         = ",   f
!! write(6,*)" chipar    = ",   chipar
!! write(6,*)" sqstar    = ",   sqstar
!! write(6,*)" sqlin     = ",   sqlin
!! write(6,*)" phistar   = ",   phistar
!! write(6,*)" philin    = ",   philin
!! write(6,*)" s33       = ",   s33
!! write(6,*)" v11-12    = ",   v11, v22, v12
!! write(6,*)" s11-12    = ",   s11, s22, s12
!! write(6,*)" v11*sqstar= ",   v22 * sqstar
!! write(6,*)" v22*sqlin = ",   v22 * sqlin
!! write(6,*)" v12**2*s*s= ",   v12**2 *sqstar*sqlin
!! write(6,*)" sq        = ",   sq
!!  
! ---- writing computed parameters to the record >>>  
!      call parset('        ',sngl(),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
     function debyestar( qs, rgs, fs ) result(sqd)
       implicit none
       double precision, intent(in) :: qs
       double precision, intent(in) :: rgs
       double precision, intent(in) :: fs

       double precision             :: sqd
       double precision             :: gam

       gam = rgs**2 * qs**2 * fs / (3*fs-2d0)

       sqd = 2/(fs*gam**2) * ( gam - (1-exp(-gam)) + (fs-1)/2d0 * (1-exp(-gam))**2 )

     end function debyestar


     function benoitmf( q, rgb, nu ) result(bmf)
       implicit none
       double precision, intent(in) :: q
       double precision, intent(in) :: rgb
       double precision, intent(in) :: nu
       double precision             :: bmf

       double precision, external   :: adapint
       double precision             :: eps, erracc
       integer                      :: maxit

       df = 1d0/nu
       rg = rgb 

       maxit = 1000
       eps   = 1d-8
       bmf = df/(q*rg)**df *  adapint(benmf_kernel,0d0, (q*rg)**2, eps, maxit, erracc)

     end function benoitmf



    function benmf_kernel(y) result(val)
      implicit none
      double precision, intent(in) :: y
      double precision             :: val

      val = (1d0-y**(df/2) / (q*rg)**df) * exp(-y) * y**(df/2-1d0)

     end function benmf_kernel

 end function th_debstar
