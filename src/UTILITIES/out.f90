 FUNCTION th_zcoil1(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Absolute value of the impedance of a coil coupled to a conductive support as function of the frequency
! 
      use theory_description 
      implicit none 
      real    :: th_zcoil1
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
     double precision :: R          ! coil resistance                                                                 
     double precision :: Rb         ! coil body resistance                                                            
     double precision :: L          ! coil inductance                                                                 
     double precision :: Lb         ! coil body inductance                                                            
     double precision :: Lcb        ! mutual inductance (coil body)                                                   
     double precision :: C          ! effective parallel capacitance                                                  
! the recin parameter representation 
! the reout parameter representation 
     double precision :: k          ! coupling factor                                                                 
 
     double precision :: th
 
     double precision :: omega
     double precision :: Zabs
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'zcoil1'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_zcoil1 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Absolute value of the impedance of a coil coupled to a conductive support as function of the frequency"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'R       '  ! coil resistance                                                                 
        parnam ( 2) = 'Rb      '  ! coil body resistance                                                            
        parnam ( 3) = 'L       '  ! coil inductance                                                                 
        parnam ( 4) = 'Lb      '  ! coil body inductance                                                            
        parnam ( 5) = 'Lcb     '  ! mutual inductance (coil body)                                                   
        parnam ( 6) = 'C       '  ! effective parallel capacitance                                                  
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "coil resistance" !//cr//parspace//&
        th_param_desc( 2,idesc) = "coil body resistance" !//cr//parspace//&
        th_param_desc( 3,idesc) = "coil inductance" !//cr//parspace//&
        th_param_desc( 4,idesc) = "coil body inductance" !//cr//parspace//&
        th_param_desc( 5,idesc) = "mutual inductance (coil body)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "effective parallel capacitance" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "k        > coupling factor"
! 
        th_zcoil1 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      R        =      pa( 1)
      Rb       =      pa( 2)
      L        =      pa( 3)
      Lb       =      pa( 4)
      Lcb      =      pa( 5)
      C        =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     omega   = 2*Pi*x
     th =  ampli *  exp(-rr * q*q / 6d0) * Sqt/Sq
     th_zcoil1 = th
 
! ---- writing computed parameters to the record >>>  
      call parset('k       ',sngl(k),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
     function Zabs(om) result(Z)
       double precision, intent(in) :: om

       Z=((L*Lb-Lcb**2)**2*om**4+(L**2*Rb**2+2*Lcb**2*R*Rb+Lb**2*R**2)*om**2+R**2*Rb**2) /  &
          (C**2*(L*Lb-Lcb**2)**2*om**6+((L**2*Rb**2+2*Lcb**2*R*Rb+Lb**2*R**2)*C- &
           2*L*Lb**2+2*Lb*Lcb**2)*C*om**4+(C**2*R**2*Rb**2-2*C*L*Rb**2+Lb**2)*om**2+Rb**2)

       Z = sqrt(Z)

     end function Zabs
 end function th_zcoil1
