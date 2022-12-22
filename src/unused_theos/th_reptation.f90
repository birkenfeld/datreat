 FUNCTION th_reptatio(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      use reptation

      implicit none 
      real    :: th_reptatio
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
     double precision :: b          ! fluctuation intensity (relative)   
     double precision :: a          ! length scale                      
     double precision :: tau        ! timescale                         
     double precision :: lz         ! total length    

     double precision :: ampdebye    ! prefactor                          
     double precision :: wl4        ! rouse rate                          
     integer          :: n_segmen   ! effective N for summation            
     double precision :: re         ! end-to-end radius of Gaussian coil     
     double precision :: temp       ! temperature if not in parameters       
     double precision :: com_diff   ! center of mass diffusion (if 0, use Rouse default) 
     double precision :: pmin       ! minimum mode to be included 
     double precision :: pmax       ! maximum mode to be included  
     double precision :: blocrep
                                                              
! the reout parameter representation 
     double precision :: W0, wr      ! rouse     W
     double precision :: l          ! sqgment length    
     double precision :: ne, n

     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset

     double precision :: dummy
!
 
! the recin parameter representation 
     double precision :: q          ! q-value    

     double precision :: sqt(2) 

     double precision              :: imod
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'reptatio'
       nparx =        20
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_reptatio = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " fermi W modificatio and b,w,ne coupling," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters  locrep t=0 limit exp(-t) "
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor   
        parnam ( 2) = 'a       '  ! length scale      
        parnam ( 3) = 'ne      '  ! entanglemen ne length  
        parnam ( 4) = 'n       '  ! n                     
        parnam ( 5) = 'wl4     '  ! rouse rate           
        parnam ( 6) = 'n_segmen'  ! effective N for summation 
        parnam ( 7) = 're      '  ! end-to-end radius of Gaussian coil  
        parnam ( 8) = 'temp    '  ! temperature if not in parameters     
        parnam ( 9) = 'com_diff'  ! center of mass diffusion (if 0, use Rouse default) 
        parnam (10) = 'pmin    '  ! minimum mode to be included                 
        parnam (11) = 'pmax    '  ! minimum mode to be included           
        parnam (12) = 't0locr  '  ! cutoff of locrep
        parnam (13) = 'twlocr  '  ! width of fermi    
        parnam (14) = 'alpha0  '  !                                                     
        parnam (15) = 'talpmax '  !                                                     
        parnam (16) = 'talpwd  '  !                                                    
        parnam (17) = 'alpoffs '  !                                                                 
        parnam (18) = 'blocrep '  !                                                                 
        parnam (19) = 'dtlocrep'  !  
        parnam (20) = 'imod    '                                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "totla n" !//cr//parspace//&
        th_param_desc( 5,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 6,idesc) = "effective N for summation" !//cr//parspace//&
        th_param_desc( 7,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 8,idesc) = "temperature if not in parameters" !//cr//parspace//&
        th_param_desc( 9,idesc) = "center of mass diffusion (if 0, use Rouse default)" !//cr//parspace//&
        th_param_desc(10,idesc) = "OBSOLETE no EFFECT minimum mode to be included" !//cr//parspace//&
        th_param_desc(11,idesc) = "OBSOLETE no EFFECT maximum mode to be included" !//cr//parspace//&
        th_param_desc(12,idesc) = "OBSOLETE NO EFFECT " !//cr//parspace//&
        th_param_desc(13,idesc) = "OBSOLETE NO EFFECT " !//cr//parspace//&
        th_param_desc(14,idesc) = "prefactor of non Gaussianity alpha(t) expression" !//cr//parspace//&
        th_param_desc(17,idesc) = "t of max in log distribution alpha(t)" !//cr//parspace//&
        th_param_desc(16,idesc) = "width of Gaussion max. in log t of alpha(t)" !//cr//parspace//&
        th_param_desc(17,idesc) = "offset in alpha(t)" !//cr//parspace//&
        th_param_desc(18,idesc) = "prefactor to locrep DEFAULT is 1 !" !//cr//parspace//&
        th_param_desc(19,idesc) = "OBSOLETE NO EFFECT tzero shift of locrep" !//cr//parspace//&
        th_param_desc(20,idesc) = "OBSOLETE NO EFFECT modulation between eq 43/46 for Lrep" !/ !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
                                             
! >>>>> describe parameters >>>>>>> 
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "diff     > diffusion in cm**2/s"
        th_out_param(  2,idesc) = "l        > effective segment length (Rouse summation)"
        th_out_param(  3,idesc) = "w0       > rouse rate: W for WL4/a**4"
        th_out_param(  4,idesc) = "wl4      > inferred value of Wl4"
        th_out_param(  5,idesc) = "rga      > sqrt(ne) * a"
!
        th_reptatio = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =        pa( 1)
      a        =   abs( pa( 2) )
      ne       =   abs( pa( 3) )
      n        =   abs( pa( 4) )
      wl4      =   abs( pa( 5))
      n_segmen =  nint( pa( 6))
      re       =   abs( pa( 7))
      temp     =        pa( 8)
      com_diff =   abs( pa( 9))
      pmin     =        pa(10)
      pmax     =        pa(11)
      dummy    =      ( pa(12))
      dummy    =   abs( pa(13))
      alpha0   =        pa(14)
      talphamax=   abs(pa(15))
      talphawd =   abs(pa(16))
      aoffset  =        pa(17)
      blocrep  =   abs( pa(18) )
      dummy    =   abs( pa(19) )
      imod     =        pa(20)

      b   = blocrep
      W0  = wl4/a**4

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 

      sqt =  reptation_sqt(q,t, N, a, Ne, Re, wl4, alpha0, talphamax,talphawd) 

      th_reptatio = ampli * sqt(2)/sqt(1)


! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(a),iadda)      ! in ns A units
      call        parset('w0      ',sngl(W0),iadda)      !     "
      call        parset('wl4     ',sngl(Wl4),iadda) !     "
      call        parset('rga     ',sngl(sqrt(ne)*a),iadda) !     "

                   
END FUNCTION th_reptatio

