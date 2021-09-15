 FUNCTION th_spatcorr(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  reptation approach , Spatial correlations of entangled polymer dynamics, Jihong Ma
!  PHYSICAL REVIEW E 104, 024503 (2021)
      use theory_description 
      implicit none 
      real    :: th_spatcorr
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
     double precision :: r2         ! msd-scale, prefactor to Kt2                                                     
     double precision :: tau        ! timescale                                                                       
     double precision :: beta       ! spatial streching exponent                                                      
     double precision :: a          ! transition exponent in K(t)                                                     
! the recin parameter representation 
     double precision :: q          ! q-value    default value                                                        
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'spatcorr'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_spatcorr = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " reptation approach , Spatial correlations of entangled polymer dynamics, Jihong Ma"
       th_citation(idesc)     = " PHYSICAL REVIEW E 104, 024503 (2021)"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'r2      '  ! msd-scale, prefactor to Kt2                                                     
        parnam ( 3) = 'tau     '  ! timescale                                                                       
        parnam ( 4) = 'beta    '  ! spatial streching exponent                                                      
        parnam ( 5) = 'a       '  ! transition exponent in K(t)                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "msd-scale, prefactor to Kt2" !//cr//parspace//&
        th_param_desc( 3,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "spatial streching exponent" !//cr//parspace//&
        th_param_desc( 5,idesc) = "transition exponent in K(t)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_spatcorr = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      r2       =      pa( 2)
      tau      =      pa( 3)
      beta     =      pa( 4)
      a        =      pa( 5)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t  = x              ! since we prefer to call the independent variable t, x must be copied to t
     th = ampli * exp(-(1d0/6d0)*(q*q*Kt2)**beta)

     th_spatcorr = th
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function Kt2(t,a) result(val)
    implicit none
    double precision, intent(in)   :: t, a
    double precision               :: val


    val = r2/((1/(t/tau)**(1/4))**a+(1/sqrt(t/tau))**a)**(1/a)


  end function Kt2
 end function th_spatcorr
