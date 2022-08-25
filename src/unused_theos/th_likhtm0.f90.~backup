 FUNCTION th_likhtm0(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  just the (nrmalized) likhtman factor
! 
      use theory_description 
      implicit none 
      real    :: th_likhtm0
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
     double precision :: N          ! N                                                                               
     double precision :: Ne         ! Ne                                                                              
     double precision :: l          ! l segment                                                                       
     double precision :: alpha      ! prefactor of the t**-1/4 CLF progression                                        
     double precision :: taue       ! time scale (tau_e) of the CLF progression                                       
! the recin parameter representation 
     double precision :: q          ! q                                                                               
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: t, mu, st
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'likhtm0'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_likhtm0 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " just the (nrmalized) likhtman factor"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'N       '  ! N                                                                               
        parnam ( 3) = 'Ne      '  ! Ne                                                                              
        parnam ( 4) = 'l       '  ! l segment                                                                       
        parnam ( 5) = 'alpha   '  ! prefactor of the t**-1/4 CLF progression                                        
        parnam ( 6) = 'taue    '  ! time scale (tau_e) of the CLF progression                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "N" !//cr//parspace//&
        th_param_desc( 3,idesc) = "Ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "l segment" !//cr//parspace//&
        th_param_desc( 5,idesc) = "prefactor of the t**-1/4 CLF progression" !//cr//parspace//&
        th_param_desc( 6,idesc) = "time scale (tau_e) of the CLF progression" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_likhtm0 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      N        =      pa( 2)
      Ne       =      pa( 3)
      l        =      pa( 4)
      alpha    =      pa( 5)
      taue     =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q
      xh =  0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     mu  =  N/12d0 * (q*l)**2
     st  =  0.5d0 * alpha*(Ne/N) *(t/taue)**0.25d0
     th  =  (2*mu+exp(-2*mu)+2-4*mu*st-4*exp(-2*mu*st)+exp(-4*mu*st))/(2*mu+exp(-2*mu)-1)
     th_likhtm0 = th
 
! ---- writing computed parameters to the record >>>  

 end function th_likhtm0
