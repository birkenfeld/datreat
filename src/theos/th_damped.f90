 FUNCTION th_damped(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  a * cos(fx*(x-dx) + a2*cos(fx2*(x-dx2))**2
! 
      use theory_description 
      implicit none 
      real    :: th_damped
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
     double precision :: dt         ! x-offset                                                                        
     double precision :: f          ! x factor                                                                        
     double precision :: damp        ! quadratic part amplitude   

     double precision :: L,C, omega0         ! induct                                                     
                                                     
! the recin parameter representation 
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'damped'
       nparx =        4
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_damped = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " a * cos(f*(x-dt)*exp(-t/damping)"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'a       '  ! prefactor                                                                       
        parnam ( 2) = 'dt      '  ! t-offset                                                                        
        parnam ( 3) = 'f       '  ! frequency                                                                       
        parnam ( 4) = 'damping '  ! quadratic part amplitude                                                        
                                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "f-offset" !//cr//parspace//&
        th_param_desc( 3,idesc) = "frequency" !//cr//parspace//&
        th_param_desc( 4,idesc) = "damping " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_damped = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      a        =      pa( 1)
      dt       =      pa( 2)
      f        =      pa( 3)
      damp     =      pa( 4)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th =  a * cos(2*Pi*f*(x-dt)) * exp(-(x-dt)/damp)

     th_damped = th

     omega0 = sqrt( (2*Pi*f)**2 - 1d0/damp**2)

     C = 4.7d-9
     L = 1d0/(C*omega0**2)

     call parset('L       ',sngl(l),iadda)  
 
! ---- writing computed parameters to the record >>>  
 

 end function th_damped
