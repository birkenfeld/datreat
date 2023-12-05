 FUNCTION th_reptestt(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptesttn expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      use reptation2

      implicit none 
      real    :: th_reptestt
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
     double precision :: a          ! segment length                      

     double precision :: wl4        ! rouse rate                          
     double precision :: re         ! end-to-end radius of Gaussian coil     
     double precision :: re_fac     ! end-to-end radius of Gaussian coil factor with respect to sqrt(ne)*l    
!     double precision :: blocrep
     double precision :: W0         ! rouse     W
     double precision :: l          ! sqgment length    
     double precision :: ne, n      ! length (number of segments) of entanglement strand, full chain
     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in NG = talphamax * tau_e
     double precision :: tau_e      ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset

! the recin parameter representation 
     double precision :: q          ! q-value    
     double precision :: sqt(2) 
     double precision :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'reptestt'
       nparx =        11
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_reptestt = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptesttn expression along the lines of deGennes, " //cr//parspace//&
                                " modification as decribed in m.monkenbusch et al, Macromolecules 2023 "! //cr//parspace//
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>" //cr//parspace//&
                                " Macrolmolecules (2023)?? (M.Monkenbusch, M.Krutyeva, D.Richter) "
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor   
        parnam ( 2) = 'a       '  ! segment length     
        parnam ( 3) = 'ne      '  ! entanglement ne length  
        parnam ( 4) = 'n       '  ! n                     
        parnam ( 5) = 'wl4     '  ! rouse rate           
        parnam ( 6) = 're_fac  '  ! end-to-end radius of Gaussian coil  
        parnam ( 7) = 'alpha0  '  ! Non-Gaussianity amplitude 
        parnam ( 8) = 'talpmax '  !                                                    
        parnam ( 9) = 'talpwd  '  !                                                    
        parnam (10) = 'bfluct  '  ! Locrep fluctuation amplitude    DEFAULT=1 !
        parnam (11) = 'bfloor  '  ! Locrep floor (1/Z) modificator  DEFAULT=1 !
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "segemnt length" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "total number of segments n" !//cr//parspace//&
        th_param_desc( 5,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 6,idesc) = "Ne to end-to-end R fcator" !//cr//parspace//&
        th_param_desc( 7,idesc) = "amplitude of Non_Gaussianity of (Rouse) entanglement blob" !//cr//parspace//&
        th_param_desc( 8,idesc) = "tau_e*talphamax = t of max in log distribution alpha(t)" !//cr//parspace//&
        th_param_desc( 9,idesc) = "width of Gaussion max. in log t of alpha(t)" !//cr//parspace//&
        th_param_desc(10,idesc) = "prefactor to fluctuation term DEFAULT is 1 !" !//cr//parspace//&
        th_param_desc(10,idesc) = "prefactor to (1/Z) intesity   DEFAULT is 1 !" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
                                            
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "l        > effective segment length (Rouse summation)"
        th_out_param(  2,idesc) = "w0       > rouse rate: W for WL4/a**4"
        th_out_param(  3,idesc) = "wl4      > inferred value of Wl4"
        th_out_param(  4,idesc) = "rga      > sqrt(ne) * a"
!
        th_reptestt = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =        pa( 1)
      a        =   abs( pa( 2) )
      ne       =   abs( pa( 3) )
      n        =   abs( pa( 4) )
      wl4      =   abs( pa( 5) )
      re_fac   =   abs( pa( 6) )
      alpha0   =        pa( 7)
      talphamax=   abs( pa( 8) )
      talphawd =   abs( pa( 9) )
      b_fluctuation = abs(pa(10))     ! shared with reptation2 enabling preliminary tests 
      b_z_floor     = abs(pa(11))     !   "     "      "

      if(b_fluctuation == 0d0) b_fluctuation = 1d0
      if(b_z_floor     == 0d0) b_z_floor     = 1d0


     
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
      Re = re_fac * sqrt(Ne) * a

      tau_e = Re**4/(Pi**2*Wl4)

      sqt =  reptation2_sqt(q,t, N, a, Ne, Re, wl4, alpha0, talphamax * tau_e ,talphawd) 

      th_reptestt = ampli * sqt(2)/sqt(1)


! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(a),iadda)      ! in ns A units
      call        parset('w0      ',sngl(W0),iadda)      !     "
      call        parset('wl4     ',sngl(Wl4),iadda) !     "
      call        parset('rga     ',sngl(sqrt(ne)*a),iadda) !     "
      call        parset('bmcdef  ',sngl(b_mcdefault),iadda) !     "  from module reptation2
      call        parset('tau_e   ',sngl(tau_e),iadda) !     "  from module reptation2
      call        parset('talmax    ',sngl(tau_e*talphamax),iadda) !     "  from module reptation2

                   
END FUNCTION th_reptestt

