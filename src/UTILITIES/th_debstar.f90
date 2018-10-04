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
     double precision :: mwstar     ! molecular weight star polymer  (g/mol)                                          
     double precision :: mwlin      ! molecular weight linear polymer(g/mol)                                          
     double precision :: mwseg      ! molecular weight of a segment  (g/mol)                                          
     double precision :: rhostar    ! density star polymer                                                            
     double precision :: rholin     ! linear polymer (g/cm**2)                                                        
     double precision :: astar      ! scattering length density star (cm**-2)                                         
     double precision :: alin       ! scattering length density lin  (cm**-2)                                         
! the reout parameter representation 
     double precision ::            !                                                                                 
 
     double precision :: th
 
     double precision, parameter :: Navogadro = 6.022140857d23
     double precision   :: q
     double precision   :: df
     double precision   :: sqstar, sqlin, sq
     double precision   :: vseg
     double precision   :: nse_star, nseg_lin
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
       th_explanation(idesc)  = " fromfactor of an f-arm star with arms with Gaussian chain statistics"
       th_citation(idesc)     = " H. Benoit, J. Polym. Sci. (1953), 11, 507-510"//cr//parspace//&
                                " G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146"//cr//parspace//&
                                " H. Benoit, Compt. Rend. (1957) 245, 2244-2247 "
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
        th_file_param(  1,idesc) = "phistar  > volume fraction of start"
        th_file_param(  2,idesc) = "mwstar   > molecular weight star polymer  (g/mol)"
        th_file_param(  3,idesc) = "mwlin    > molecular weight linear polymer(g/mol)"
        th_file_param(  4,idesc) = "mwseg    > molecular weight of a segment  (g/mol)"
        th_file_param(  5,idesc) = "rhostar  > density star polymer"
        th_file_param(  6,idesc) = "rholin   > linear polymer (g/cm**2)"
        th_file_param(  7,idesc) = "astar    > scattering length density star (cm**-2)"
        th_file_param(  8,idesc) = "alin     > scattering length density lin  (cm**-2)"
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
      f        =      pa( 2)
      rgstar   =      pa( 3)
      nulin    =      pa( 4)
      rglin    =      pa( 5)
      chipar   =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: volume fraction of start
      xh = 
      call parget('phistar ',xh,iadda,ier)
      phistar  = xh
! >>> extract: molecular weight star polymer  (g/mol)
      xh = 
      call parget('mwstar  ',xh,iadda,ier)
      mwstar   = xh
! >>> extract: molecular weight linear polymer(g/mol)
      xh = 
      call parget('mwlin   ',xh,iadda,ier)
      mwlin    = xh
! >>> extract: molecular weight of a segment  (g/mol)
      xh = 
      call parget('mwseg   ',xh,iadda,ier)
      mwseg    = xh
! >>> extract: density star polymer
      xh = 
      call parget('rhostar ',xh,iadda,ier)
      rhostar  = xh
! >>> extract: linear polymer (g/cm**2)
      xh = 
      call parget('rholin  ',xh,iadda,ier)
      rholin   = xh
! >>> extract: scattering length density star (cm**-2)
      xh = 
      call parget('astar   ',xh,iadda,ier)
      astar    = xh
! >>> extract: scattering length density lin  (cm**-2)
      xh = 
      call parget('alin    ',xh,iadda,ier)
      alin     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x

     nseg_star = mwstar/mwseg
     nseg_lin  = mwlin/mwseg
     vseg      = mwseg / (rholin * Navogadro)
     sqstar  =  debyestar(q, rgstar, f)
     sqlin   =  benoitmf( q, rglin, nulin )

     sqstar  = sqstar * nseg_star * phistar
     sqslin  = sqlin  * nseg_lin  * (phistar-1)

     sq = 1d0/ ( 1d0/sqstar + 1d0/sqlin - 2*chipar )

     th = ampli * ((astar-alin))**2 * vseg * sq


     th_debstar = th
 
! ---- writing computed parameters to the record >>>  
!      call parset('        ',sngl(),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
     function debyestar( qs, rgs, fs ) result(sq)
       implicit none
       double precision, intent(in) :: qs
       double precision, intent(in) :: rgs
       double precision, intent(in) :: fs

       double precision             :: sq
       double precision             :: gam

       gam = rgs**2 * q**2 * fs / (3*fs-2d0)

       sq = 2/(f*gam**2) * ( gam - (1-exp(-gam)) + (f-1)/2d0 * (1-exp(-gam))**2 )

     end function debyestar


     function benoitmf( q, rg, nu ) result(bmf)
       implicit none
       double precision, intent(in) :: q
       double precision, intent(in) :: rg
       double precision, intent(in) :: nu
       double precision             :: bmf

       double precision             :: eps, erracc
       integer                      :: maxit

       df = 1d0/nu

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
