 FUNCTION th_trigfun(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  a * cos(fx*(x-dx) + a2*cos(fx2*(x-dx2))**2
! 
      use theory_description 
      implicit none 
      real    :: th_trigfun
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
     double precision :: a          ! prefactor                                                                       
     double precision :: dx         ! x-offset                                                                        
     double precision :: fx         ! x factor                                                                        
     double precision :: a2         ! quadratic part amplitude                                                        
     double precision :: fx2        ! quadratic part x factor                                                         
     double precision :: dx2        ! quadratic part x-offset                                                         
! the recin parameter representation 
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'trigfun'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_trigfun = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " a * cos(fx*(x-dx) + a2*cos(fx2*(x-dx2))**2"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'a       '  ! prefactor                                                                       
        parnam ( 2) = 'dx      '  ! x-offset                                                                        
        parnam ( 3) = 'fx      '  ! x factor                                                                        
        parnam ( 4) = 'a2      '  ! quadratic part amplitude                                                        
        parnam ( 5) = 'fx2     '  ! quadratic part x factor                                                         
        parnam ( 6) = 'dx2     '  ! quadratic part x-offset                                                         
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "x-offset" !//cr//parspace//&
        th_param_desc( 3,idesc) = "x factor" !//cr//parspace//&
        th_param_desc( 4,idesc) = "quadratic part amplitude" !//cr//parspace//&
        th_param_desc( 5,idesc) = "quadratic part x factor" !//cr//parspace//&
        th_param_desc( 6,idesc) = "quadratic part x-offset" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_trigfun = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      a        =      pa( 1)
      dx       =      pa( 2)
      fx       =      pa( 3)
      a2       =      pa( 4)
      fx2      =      pa( 5)
      dx2      =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th =  a * cos(fx*(x-dx) + a2*cos(fx2*(x-dx2))**2

     th_trigfun = th
 
! ---- writing computed parameters to the record >>>  
 

 end function th_trigfun
