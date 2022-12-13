 FUNCTION th_fjchainp(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  freely joint chain statistics test P(N)= Anorm * N**eta1**(-2/3)**(3/2)*exp(-r**2*3d0/(2*a**2*N**eta2))
! 
      use theory_description 
      implicit none 
      real    :: th_fjchainp
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
     double precision :: ampli
     double precision :: r, a                                                                       
     double precision :: eta, eta2                                                                                     
! the recin parameter representation 
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: n
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'fjchainp'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_fjchainp = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " freely joint chain statistics test P(N)= (3/(2a**2*N))**(3/2)*exp(-r**2*3d0/(2*a**2*N))"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'r       '  ! radius of gyration                                                              
        parnam ( 3) = 'a       '  ! l                                                                               
        parnam ( 4) = 'eta     '  ! l                                                                               
        parnam ( 5) = 'eta2    '  ! l                                                                               
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "end-to-end distance" !//cr//parspace//&
        th_param_desc( 3,idesc) = "segment length" !//cr//parspace//&
        th_param_desc( 4,idesc) = "modiying exponen normalizing term" !//cr//parspace//&
        th_param_desc( 5,idesc) = "modifying exponen exp term" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_fjchainp = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      r        =      pa( 2)
      a        =      pa( 3)
      eta      =      pa( 4)
      eta2     =      pa( 5)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     n   = x
     th_fjchainp = ampli * p_of_n(n)
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function p_of_n(n) result(val)
    implicit none
    double precision, intent(in) :: n
    double precision             :: val

    if(abs(n)  < epsilon(n) ) then
      val = 0d0
    else    !! this is the normalize version with int(P(N)dN)=1 !!
      val = (0.2D1 / 0.3D1 * R ** 2 / a ** 2) ** ((eta - 1) / eta2) & 
           * dble(eta2) * N ** (-eta) * exp(-0.2D1 / 0.3D1 * R ** 2 / &
           a** 2 * N ** (-eta2)) / GAMMA((eta - 1) / eta2)  
    endif



  end function p_of_n
 end function th_fjchainp
