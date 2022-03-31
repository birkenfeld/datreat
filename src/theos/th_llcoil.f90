 FUNCTION th_llcoil(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  compute the frequency dependent impedance of the sc-maincoil
! 
      use theory_description 
      use lcalc_Zomega_cc
      implicit none 
      real    :: th_llcoil
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
     double precision :: dlengtha    ! extralength of winding cores                                                    
     double precision :: wwidtha     ! thickness of winding cores                                                      
     double precision :: rsepa       ! separation radius winding-winding core                                          
     double precision :: cper       ! capacitance per squaremetzer winding/core                                       
     double precision :: cpa        ! overll parallel capacitance                                                     
     double precision :: rhoal      ! specific resistity of Alcore (ca. 2e-9)                                         
     double precision :: r0         ! series resistance (of connection)                                               
! the recin parameter representation 
     double precision :: abspha     ! 1=abs, 2=phase of Z                                                             
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)                                                       
     double precision :: Rwc_mod(10)
 
     double precision :: th
 
     double precision   :: frequenz
     complex            :: zx
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'llcoil'
       nparx =        18
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_llcoil = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " compute the frequency dependent impedance of the sc-maincoil" //cr//parspace//&
        " NOTE: module lcalc_Zomega_cc must be manually compiled after make clean " //cr//parspace//&
        "     ...some not resulved clash with the integrated adapint functions, or rename them .... TBD"

       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'dlength '  ! extralength of winding cores                                                    
        parnam ( 3) = 'wwidth  '  ! thickness of winding cores                                                      
        parnam ( 4) = 'rsep    '  ! separation radius winding-winding core                                          
        parnam ( 5) = 'cpersqm '  ! capacitance per squaremetzer winding/core                                       
        parnam ( 6) = 'cparalle'  ! overll parallel capacitance                                                     
        parnam ( 7) = 'rhoal   '  ! specific resistity of Alcore (ca. 2e-9)                                         
        parnam ( 8) = 'r0      '  ! series resistance (of connection)
                                               
        parnam ( 9) = 'r_mod_1 '  ! winding core resistance modifier                                               
        parnam (10) = 'r_mod_2 '  ! winding core resistance modifier
        parnam (11) = 'r_mod_3 '  ! winding core resistance modifier
        parnam (12) = 'r_mod_4 '  ! winding core resistance modifier
        parnam (13) = 'r_mod_5 '  ! winding core resistance modifier
        parnam (14) = 'r_mod_6 '  ! winding core resistance modifier
        parnam (15) = 'r_mod_7 '  ! winding core resistance modifier
        parnam (16) = 'r_mod_8 '  ! winding core resistance modifier
        parnam (17) = 'r_mod_9 '  ! winding core resistance modifier
        parnam (18) = 'r_mod_10'  ! winding core resistance modifier
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "extralength of winding cores" !//cr//parspace//&
        th_param_desc( 3,idesc) = "thickness of winding cores" !//cr//parspace//&
        th_param_desc( 4,idesc) = "separation radius winding-winding core" !//cr//parspace//&
        th_param_desc( 5,idesc) = "capacitance per squaremetzer winding/core" !//cr//parspace//&
        th_param_desc( 6,idesc) = "overll parallel capacitance" !//cr//parspace//&
        th_param_desc( 7,idesc) = "specific resistity of Alcore (ca. 2e-9)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc( 9,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(10,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(11,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(12,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(13,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(14,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(15,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(16,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(17,idesc) = "winding core resistance modifier" !//cr//parspace//&
        th_param_desc(18,idesc) = "winding core resistance modifier" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "abspha   > 1=abs, 2=phase of Z"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "intens0  > intensity factor (cm**-1)"
! 
        th_llcoil = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli     =      pa( 1)
      dlengtha  =     abs( pa( 2) )
      wwidtha   =     abs( pa( 3) )
      rsepa     =     abs( pa( 4) )
      cper      =     abs( pa( 5) )
      cpa       =     abs( pa( 6) )
      rhoal     =     abs( pa( 7) )
      r0        =     abs( pa( 8) )
      Rwc_mod   =     abs( pa(9:18) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 1=abs, 2=phase of Z
      xh = 1
      call parget('abspha  ',xh,iadda,ier)
      abspha   = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     frequenz   = x
     zx         =  Z_llcoil(frequenz, rhoAl, cpa, cper, dlengtha, wwidtha, rsepa, Rwc_mod) + R0
     zx         = zx * ampli
     if(nint(abspha) <= 1)  then
        th = abs(zx)
     else
        th = atan2(aimag(zx),real(zx))*180d0/Pi
     endif



     th_llcoil = th
 
! ---- writing computed parameters to the record >>>  
      call parset('abspha ',sngl(abspha),iadda,ier)
 


 end function th_llcoil
