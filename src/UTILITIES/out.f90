<<<<<<< HEAD
 FUNCTION th_benoitmf(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  general expression allowing the description of the small angle scattering of a mass fractal according to benoit (1957)
!  G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146; H. Benoit, Compt. Rend. (1957) 245, 2244-2247
      use theory_description 
      implicit none 
      real    :: th_benoitmf
=======
 FUNCTION th_havneg1(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Havriliak Negami function with tau following Vogel-Fulcher-Tamman tau=tau0 * exp(B/(T-Tvft))
!  Havriliak and Negami
      use theory_description 
      implicit none 
      real    :: th_havneg1
>>>>>>> mm-develop
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
<<<<<<< HEAD
     double precision :: ampli      ! prefactor                                                                       
     double precision :: rg         ! radius of gyration                                                              
     double precision :: df         ! fractal dimension, see eq (6) in Beaucages Paper                                
! the recin parameter representation 
     double precision :: thick      ! Dicke der Probe                                                                 
! the reout parameter representation 
 
     double precision   :: q
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'benoitmf'
       nparx =        3
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_benoitmf = 0
=======
     double precision :: amplitud   ! prefactor                                                                       
     double precision :: tau0       ! characteristic time                                                             
     double precision :: B          ! Vogel-Fulcher B-Parameter                                                       
     double precision :: Tvft       ! Vogel-Fulcher Temperature                                                       
     double precision :: alpha      ! Havriliak Negami alpha                                                          
     double precision :: beta       ! Havriliak Negami beta                                                           
! the recin parameter representation 
     double precision :: temp       ! temperature                                                                     
     double precision :: re_im      ! 0=real (epsilon) 1 =-imag (epsilon)                                             
! the reout parameter representation 
     double precision :: dummy      !                                                                                 
 
     double precision :: th
 
     double precision   :: om
     double precision   :: tau
     complex(kind=8)    :: ci=(0d0,1d0)
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'havneg1'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_havneg1 = 0
>>>>>>> mm-develop
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
<<<<<<< HEAD
       th_explanation(idesc)  = " general expression allowing the description of the small angle scattering of a mass fractal according to benoit (1957)"
       th_citation(idesc)     = " G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146; H. Benoit, Compt. Rend. (1957) 245, 2244-2247"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'rg      '  ! radius of gyration                                                              
        parnam ( 3) = 'df      '  ! fractal dimension, see eq (6) in Beaucages Paper                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "radius of gyration" !//cr//parspace//&
        th_param_desc( 3,idesc) = "fractal dimension, see eq (6) in Beaucages Paper" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "thick    > Dicke der Probe"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_benoitmf = 0.0
=======
       th_explanation(idesc)  = " Havriliak Negami function with tau following Vogel-Fulcher-Tamman tau=tau0 * exp(B/(T-Tvft))"
       th_citation(idesc)     = " Havriliak and Negami"
!       --------------> set the parameter names --->
        parnam ( 1) = 'amplitud'  ! prefactor                                                                       
        parnam ( 2) = 'tau0    '  ! characteristic time                                                             
        parnam ( 3) = 'B       '  ! Vogel-Fulcher B-Parameter                                                       
        parnam ( 4) = 'Tvft    '  ! Vogel-Fulcher Temperature                                                       
        parnam ( 5) = 'alpha   '  ! Havriliak Negami alpha                                                          
        parnam ( 6) = 'beta    '  ! Havriliak Negami beta                                                           
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "characteristic time" !//cr//parspace//&
        th_param_desc( 3,idesc) = "Vogel-Fulcher B-Parameter" !//cr//parspace//&
        th_param_desc( 4,idesc) = "Vogel-Fulcher Temperature" !//cr//parspace//&
        th_param_desc( 5,idesc) = "Havriliak Negami alpha" !//cr//parspace//&
        th_param_desc( 6,idesc) = "Havriliak Negami beta" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "temp     > temperature"
        th_file_param(  2,idesc) = "re_im    > 0=real (epsilon) 1 =-imag (epsilon)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "dummy    > "
! 
        th_havneg1 = 0.0
>>>>>>> mm-develop
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
<<<<<<< HEAD
      ampli    =      pa( 1)
      rg       =      pa( 2)
      df       =      pa( 3)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: Dicke der Probe
      xh =  0.1
      call parget('thick   ',xh,iadda,ier)
      thick    = xh
=======
      amplitud =      pa( 1)
      tau0     =      pa( 2)
      B        =      pa( 3)
      Tvft     =      pa( 4)
      alpha    =      pa( 5)
      beta     =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: temperature
      xh =          300
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: 0=real (epsilon) 1 =-imag (epsilon)
      xh =          1
      call parget('re_im   ',xh,iadda,ier)
      re_im    = xh
>>>>>>> mm-develop
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
<<<<<<< HEAD
     q   = x
     th_benoitmf = ! INSERT RESULTING VALUE OF TH EVALUATION HERE
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function benmf_kernel(y) result(val)
    implicit none
    double precision, intent(in) :: y
    double precision             :: val

    val = (1d0-y**(df/2) / (q*rg)**df) * exp(-y) * y**(df/2-1d0)

  end function bc_sum
 end function th_benoitmf
=======
     om   = x
     tau  = tau0 * exp(B/(temp-Tvft))
     th   = Real( amplitud * ci**(nint(re_im)) / ( 1d0 + (ci*om*tau)**alpha )**beta )
     th_havneg1 = th
 
! ---- writing computed parameters to the record >>>  
      call parset('dummy   ',sngl(dummy),iadda,ier)
 end function th_havneg1
>>>>>>> mm-develop
