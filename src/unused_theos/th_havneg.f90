 FUNCTION th_havneg(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Havriliak Negami function with tau following Vogel-Fulcher-Tamman tau=tau0 * exp(B/(T-Tvft)
!  Havriliak and Negami
      use theory_description 
      implicit none 
      real    :: th_havneg
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
     double precision :: amplitud   ! prefactor                                                                       
     double precision :: tau0       ! characteristic time                                                             
     double precision :: B          ! Vogel-Fulcher B-Parameter                                                       
     double precision :: Tvft       ! Vogel-Fulcher Temperature                                                       
     double precision :: alpha      ! Havriliak Negami alpha                                                          
     double precision :: beta       ! Havriliak Negami beta                                                           
! the recin parameter representation 
     double precision :: temp       ! temperature                                                                     
     integer          :: re_im      ! 1=real (epsilon) else =-imag (epsilon)                                          
                                                                          
 
     double precision   :: th
 
     double precision   :: om
     double precision   :: tau
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'havneg'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_havneg = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Havriliak Negami function with tau following Vogel-Fulcher-Tamman tau=tau0 * exp(B/(T-Tvft)"
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
        th_file_param(  2,idesc) = "re_im    > 1=real (epsilon) else =-imag (epsilon)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "         > "
! 
        th_havneg = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      amplitud =      pa( 1)
      tau0     =  abs(pa( 2)) 
      B        =  abs(pa( 3))
      Tvft     =  abs(pa( 4))
      alpha    =  abs(pa( 5))
      beta     =  abs(pa( 6))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: temperature
      xh = 300
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: 1=real (epsilon) else =-imag (epsilon)
      xh = 2
      call parget('re_im   ',xh,iadda,ier)
      re_im    = nint(xh)
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     tau  = tau0 * exp(B*(temp-Tvft))
     om   = x

     th = amplitud * havneg(om,tau,alpha,beta,re_im)


     th_havneg = th
 
! ---- writing computed parameters to the record >>>  
!      call parset('        ',sngl(),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
     function  havneg(om,tau,alpha,beta,re_im) result(hn)
       implicit none
       double precision, intent(in) :: om
       double precision, intent(in) :: tau
       double precision, intent(in) :: alpha
       double precision, intent(in) :: beta
       integer         , intent(in) :: re_im

       double precision             :: hn

       complex(kind=8) :: y

       y = (1d0,0d0) / ( (1d0,0d0) + ((0d0,1d0)*om*tau)**alpha)**beta

       if(re_im == 1) then
        hn =  RealPart(y)
       else
        hn = -ImagPart(y)
       endif

     end function havneg

 end function th_havneg
