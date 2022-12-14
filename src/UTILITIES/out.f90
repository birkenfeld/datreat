 FUNCTION th_fjchainp(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  freely joint chain statistics test P(N)= (3/(2a**2*N))**(3/2)*exp(-r**2*3d0/(2*a**2*N))
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
     double precision :: ampli      ! prefactor                                                                       
     double precision :: r          ! radius of gyration                                                              
     double precision :: a          ! l                                                                               
! the recin parameter representation 
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: n
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'fjchainp'
       nparx =        3
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
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "radius of gyration" !//cr//parspace//&
        th_param_desc( 3,idesc) = "l" !//cr//parspace//&
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

    val = (3d0/(2d0*a**2*n))**(3d0/2d0) * exp(-r**2 * 3d0/(2d0*a**2*n))

  end function p_of_n
 end function th_fjchainp
