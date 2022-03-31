 FUNCTION th_spatcorr2(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  reptation approach , Spatial correlations of entangled polymer dynamics, Jihong Ma
!  PHYSICAL REVIEW E 104, 024503 (2021)
      use theory_description 
      implicit none 
      real    :: th_spatcorr2
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
     double precision :: exp1, exp2
     double precision :: betagt0, betagwd, betat
     integer          :: mode
! the recin parameter representation 
     double precision :: q          ! q-value    default value                                                        
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'spatcor2'
       nparx =        10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_spatcorr2 = 0
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
        parnam ( 6) = 'exp1    '  ! t-exp 1 (1/4)                                                    
        parnam ( 7) = 'exp2    '  ! t-exp 2 (1/2)                                                    
        parnam ( 8) = 'betagt0 '  !                                                     
        parnam ( 9) = 'betagwd '  !                                                     
        parnam (10) = 'mode    '  !                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "msd-scale, prefactor to Kt2" !//cr//parspace//&
        th_param_desc( 3,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "spatial streching exponent" !//cr//parspace//&
        th_param_desc( 5,idesc) = "transition exponent in K(t)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "small t exp (1/4)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "large t exp (1/2)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "gaussian center beta" !//cr//parspace//&
        th_param_desc( 9,idesc) = "gaussian width beta " !//cr//parspace//&
        th_param_desc(10,idesc) = "select 0=Sqt, 1=xi2(t), 2=beta(t) " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_spatcorr2 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      r2       =  abs(pa( 2))
      tau      =  abs(pa( 3))
      beta     =  abs(pa( 4))
      a        =  abs(pa( 5))
      exp1     =  abs(pa( 6))
      exp2     =  abs(pa( 7))
      betagt0  =  abs(pa( 8))
      betagwd  =  abs(pa( 9))
      mode     = nint(pa(10))
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

     betat = 1 + (beta-1)*exp(-(log(t/betagt0)/betagwd)**2)

     select case(mode)
     case(0)
       th = ampli * exp(-(1d0/6d0)*(q*q*Kt2(t,a))**betat)
     case(1)
       th = Kt2(t,a)
     case(2) 
       th = betat
     case default 
       th = ampli * exp(-(1d0/6d0)*(q*q*Kt2(t,a))**betat)
     end select

!      th = ampli * exp(-(1d0/6d0)*(q*q*Kt2(t,a))**betat)

     th_spatcorr2 = th
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function Kt2(t,a) result(val)
    implicit none
    double precision, intent(in)   :: t, a
    double precision               :: val


    val = r2/((1d0/(t/tau)**(exp1))**a+(1d0/(t/tau)**exp2)**a)**(1d0/a)


  end function Kt2


 end function th_spatcorr2
