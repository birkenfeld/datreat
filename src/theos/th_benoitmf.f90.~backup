 FUNCTION th_benoitmf(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  general expression allowing the description of the small angle scattering of a mass fractal according to benoit (1957)
!  ..a generalisation of the Debye chain formula for chasin statistyics deviationg from nu=0.5 (df=2)
!  G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146; H. Benoit, Compt. Rend. (1957) 245, 2244-2247
      use theory_description 
      implicit none 
      real    :: th_benoitmf
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
     double precision :: rg         ! radius of gyration                                                              
     double precision :: df         ! fractal dimension, see eq (6) in Beaucages Paper  
     double precision :: va2        ! A2 like excluded vol parmeter                              
! the recin parameter representation 
     double precision :: thick      ! Dicke der Probe                                                                 
! the reout parameter representation 
 
     double precision   :: q

     double precision   :: adapint, eps=1d-8, erracc
     integer            :: maxit  = 1000
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'benoitmf'
       nparx =        4
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_benoitmf = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " general expression allowing the description of the small angle scattering of a"//cr//&
                                " chain with satistics deviating from nu=0.5 (df=2) (Debye formula, mass fractal)"//cr//&
                                " according to benoit (1957)"
       th_citation(idesc)     = " G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146; H. Benoit, Compt. Rend. (1957) 245, 2244-2247"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'rg      '  ! radius of gyration                                                              
        parnam ( 3) = 'df      '  ! fractal dimension, see eq (6) in Beaucages Paper                                
        parnam ( 4) = 'va2     '  ! a2 type excluded volume parameter                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "radius of gyration" !//cr//parspace//&
        th_param_desc( 3,idesc) = "fractal dimension, see eq (6) in Beaucages Paper" !//cr//parspace//&
        th_param_desc( 4,idesc) = "a2 type exclude vol param:   sq <- 1(1/sq + va2)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "thick    > Dicke der Probe"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_benoitmf = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      rg       =      pa( 2)
      df       =      pa( 3)
      va2      =      pa( 4)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: Dicke der Probe
      xh =  0.1
      call parget('thick   ',xh,iadda,ier)
      thick    = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     th_benoitmf = df/(q*rg)**df * &
                   adapint(benmf_kernel,0d0, (q*rg)**2, eps, maxit, erracc)
 

     th_benoitmf = 1d0 / (1d0/th_benoitmf + va2)

     th_benoitmf = ampli * th_benoitmf

! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function benmf_kernel(y) result(val)
    implicit none
    double precision, intent(in) :: y
    double precision             :: val

    val = (1d0-y**(df/2) / (q*rg)**df) * exp(-y) * y**(df/2-1d0)

  end function benmf_kernel
 end function th_benoitmf
