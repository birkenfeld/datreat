 FUNCTION th_llcoil(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  compute the frequency dependent impedance of the sc-maincoil
! 
      use theory_description 
      use lcalc_20_Zomega_cc
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
     double precision :: dlength    ! extralength of winding cores                                                    
     double precision :: wwidth     ! thickness of winding cores                                                      
     double precision :: rsep       ! separation radius winding-winding core                                          
     double precision :: cpersqm    ! capacitance per squaremetzer winding/core                                       
     double precision :: cparalle   ! overll parallel capacitance                                                     
     double precision :: rhoal      ! specific resistity of Alcore (ca. 2e-9)                                         
     double precision :: r0         ! series resistance (of connection)                                               
! the recin parameter representation 
     double precision :: abspha     ! 1=abs, 2=phase of Z                                                             
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)                                                       
 
     double precision :: th
 
     double precision   :: frequenz
     complex            :: zx
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'llcoil'
       nparx =        8
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_llcoil = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " compute the frequency dependent impedance of the sc-maincoil"
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
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "extralength of winding cores" !//cr//parspace//&
        th_param_desc( 3,idesc) = "thickness of winding cores" !//cr//parspace//&
        th_param_desc( 4,idesc) = "separation radius winding-winding core" !//cr//parspace//&
        th_param_desc( 5,idesc) = "capacitance per squaremetzer winding/core" !//cr//parspace//&
        th_param_desc( 6,idesc) = "overll parallel capacitance" !//cr//parspace//&
        th_param_desc( 7,idesc) = "specific resistity of Alcore (ca. 2e-9)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "series resistance (of connection)" !//cr//parspace//&
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
      ampli    =      pa( 1)
      dlength  =      pa( 2)
      wwidth   =      pa( 3)
      rsep     =      pa( 4)
      cpersqm  =      pa( 5)
      cparalle =      pa( 6)
      rhoal    =      pa( 7)
      r0       =      pa( 8)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 1=abs, 2=phase of Z
      xh = 
      call parget('abspha  ',xh,iadda,ier)
      abspha   = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     frequenz   = x
     zx         =  Z_llcoil(frequenz, rhoAl, cparallel, cpersqm, dlength, wwidth, rsep) + R0
     zx         = zx * ampli
     if(nint(abspha) <= 1)  then
        th = abs(zx)
     else
        th = atan2(aimag(zx),real(zx))*180d0/Pi
     endif



     th_llcoil = th
 
! ---- writing computed parameters to the record >>>  
      call parset('intens0 ',sngl(intens0),iadda,ier)
 


 end function th_llcoil
