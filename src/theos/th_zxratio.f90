 FUNCTION th_zxratio(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  impedance or RC bridge ratio for Quench detection
!  QC
      use theory_description 
      implicit none 
      real    :: th_zxratio
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
     double precision :: R10        ! base resistor branch 1                                                          
     double precision :: R20        ! base resistor branch 2                                                          
     double precision :: R11        ! series resistor for C1 in branch 1                                              
     double precision :: R11p       ! parallel resistor for C1 in branch 1                                            
     double precision :: C11        ! capacitance C1 in branch 2                                                      
     double precision :: R21        ! series resistor for C1 in branch 2                                              
     double precision :: R21p       ! parallel resistor for C1 in branch 2                                            
     double precision :: C21        ! capacitance C2 in branch 2                                                      
     double precision :: R12        ! series resistor for C2 in branch 1                                              
     double precision :: R12p       ! parallel resistor for C2 in branch 1                                            
     double precision :: C12        ! capacitance C2 in branch 2                                                      
     double precision :: R22        ! series resistor for C2 in branch 2                                              
     double precision :: R22p       ! parallel resistor for C2 in branch 2                                            
     double precision :: C22        ! capacitance C2 in branch 2                                                      
     double precision :: R13        ! series resistor for C3 in branch 1                                              
     double precision :: R13p       ! parallel resistor for C3 in branch 1                                            
     double precision :: C13        ! capacitance C3 in branch 2                                                      
     double precision :: R23        ! series resistor for C3 in branch 2                                              
     double precision :: R23p       ! parallel resistor for C1 in branch 2                                            
     double precision :: C23        ! capacitance C3 in branch 2                                                      
! the recin parameter representation 
     double precision :: re_im      ! real imag selectopr                                                             
! the reout parameter representation 
     double precision :: dummy      ! dummy                                                                           
 
     double precision :: th
 
     double precision   :: freq
     integer, parameter :: n = 3
     integer            :: reim
     double precision   :: R1(0:n), C1(1:n), R1p(1:n), R2(0:n), C2(1:n), R2p(1:n)
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'zxratio'
       nparx =       20
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_zxratio = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " impedance or RC bridge ratio for Quench detection"
       th_citation(idesc)     = " QC"
!       --------------> set the parameter names --->
        parnam ( 1) = 'R10     '  ! base resistor branch 1                                                          
        parnam ( 2) = 'R20     '  ! base resistor branch 2                                                          
        parnam ( 3) = 'R11     '  ! series resistor for C1 in branch 1                                              
        parnam ( 4) = 'R11p    '  ! parallel resistor for C1 in branch 1                                            
        parnam ( 5) = 'C11     '  ! capacitance C1 in branch 2                                                      
        parnam ( 6) = 'R21     '  ! series resistor for C1 in branch 2                                              
        parnam ( 7) = 'R21p    '  ! parallel resistor for C1 in branch 2                                            
        parnam ( 8) = 'C21     '  ! capacitance C2 in branch 2                                                      
        parnam ( 9) = 'R12     '  ! series resistor for C2 in branch 1                                              
        parnam (10) = 'R12p    '  ! parallel resistor for C2 in branch 1                                            
        parnam (11) = 'C12     '  ! capacitance C2 in branch 2                                                      
        parnam (12) = 'R22     '  ! series resistor for C2 in branch 2                                              
        parnam (13) = 'R22p    '  ! parallel resistor for C2 in branch 2                                            
        parnam (14) = 'C22     '  ! capacitance C2 in branch 2                                                      
        parnam (15) = 'R13     '  ! series resistor for C3 in branch 1                                              
        parnam (16) = 'R13p    '  ! parallel resistor for C3 in branch 1                                            
        parnam (17) = 'C13     '  ! capacitance C3 in branch 2                                                      
        parnam (18) = 'R23     '  ! series resistor for C3 in branch 2                                              
        parnam (19) = 'R23p    '  ! parallel resistor for C1 in branch 2                                            
        parnam (20) = 'C23     '  ! capacitance C3 in branch 2                                                      
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "base resistor branch 1" !//cr//parspace//&
        th_param_desc( 2,idesc) = "base resistor branch 2" !//cr//parspace//&
        th_param_desc( 3,idesc) = "series resistor for C1 in branch 1" !//cr//parspace//&
        th_param_desc( 4,idesc) = "parallel resistor for C1 in branch 1" !//cr//parspace//&
        th_param_desc( 5,idesc) = "capacitance C1 in branch 2" !//cr//parspace//&
        th_param_desc( 6,idesc) = "series resistor for C1 in branch 2" !//cr//parspace//&
        th_param_desc( 7,idesc) = "parallel resistor for C1 in branch 2" !//cr//parspace//&
        th_param_desc( 8,idesc) = "capacitance C2 in branch 2" !//cr//parspace//&
        th_param_desc( 9,idesc) = "series resistor for C2 in branch 1" !//cr//parspace//&
        th_param_desc(10,idesc) = "parallel resistor for C2 in branch 1" !//cr//parspace//&
        th_param_desc(11,idesc) = "capacitance C2 in branch 2" !//cr//parspace//&
        th_param_desc(12,idesc) = "series resistor for C2 in branch 2" !//cr//parspace//&
        th_param_desc(13,idesc) = "parallel resistor for C2 in branch 2" !//cr//parspace//&
        th_param_desc(14,idesc) = "capacitance C2 in branch 2" !//cr//parspace//&
        th_param_desc(15,idesc) = "series resistor for C3 in branch 1" !//cr//parspace//&
        th_param_desc(16,idesc) = "parallel resistor for C3 in branch 1" !//cr//parspace//&
        th_param_desc(17,idesc) = "capacitance C3 in branch 2" !//cr//parspace//&
        th_param_desc(18,idesc) = "series resistor for C3 in branch 2" !//cr//parspace//&
        th_param_desc(19,idesc) = "parallel resistor for C1 in branch 2" !//cr//parspace//&
        th_param_desc(20,idesc) = "capacitance C3 in branch 2" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "re_im    > real imag selectopr"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "dummy    > dummy"
! 
        th_zxratio = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      R10      =     abs( pa( 1))  
      R20      =     abs( pa( 2))
      R11      =     abs( pa( 3))
      R11p     =     abs( pa( 4))
      C11      =     abs( pa( 5))
      R21      =     abs( pa( 6))
      R21p     =     abs( pa( 7))
      C21      =     abs( pa( 8))
      R12      =     abs( pa( 9))
      R12p     =     abs( pa(10))
      C12      =     abs( pa(11))
      R22      =     abs( pa(12))
      R22p     =     abs( pa(13))
      C22      =     abs( pa(14))
      R13      =     abs( pa(15))
      R13p     =     abs( pa(16))
      C13      =     abs( pa(17))
      R23      =     abs( pa(18))
      R23p     =     abs( pa(19))
      C23      =     abs( pa(20))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: real imag selectopr
      xh = 0 
      call parget('re_im   ',xh,iadda,ier)
      re_im    = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     freq   = x
     reim   = nint(re_im)
     R1     = [R10,R11,R12,R13]
     R2     = [R20,R21,R22,R23]
     C1     = [C11,C12,C13]
     C2     = [C21,C22,C23]
     R1p    = [R11p, R12p, R13p]
     R2p    = [R21p, R22p, R23p]

     th = zxratio(freq,reim , R1, C1, R1p, R2, C2, R2p, n )

     th_zxratio = th
 
! ---- writing computed parameters to the record >>>  
      call parset('dummy   ',sngl(dummy),iadda)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

function zxratio(freq, re_im, R1, C1, R1p, R2, C2, R2p, n ) result(zr)
  implicit none
  double precision, intent(in)  :: freq
  integer         , intent(in)  :: re_im      ! select 1=Real part , 2=Imaginary, 0=abs
  double precision, intent(in)  :: R1(0:)     ! Array with resistivity values of lower branch
  double precision, intent(in)  :: C1(1:)     ! Array with capacitance values of lower branch
  double precision, intent(in)  :: R1p(1:)    ! Array with R parallel capacitance values of lower branch
  double precision, intent(in)  :: R2(0:)     ! Array with resistivity values of upper branch
  double precision, intent(in)  :: C2(1:)     ! Array with capacitance values of upper branch
  double precision, intent(in)  :: R2p(1:)    ! Array with R parallel capacitance values of upper branch
  integer         , intent(in)  :: n          ! active parallel elements 1

  double precision :: zr
  double precision :: omega
  double precision, parameter :: Pi=4d0*atan(1d0)
  complex(kind=8) , parameter :: CI=(0d0,1d0)
  complex(kind=8)  :: X1, X2, X12r
  complex(kind=8)  :: Ya, Yb, YCa, YCb
  integer          :: i



  omega = 2*Pi*freq

  X1 = R1(0)
  X2 = R2(0)

  do i=1,n
    YCa = 1d0/(CI*omega*C1(i) + 1d0/R1p(i))
    Ya = R1(i) + Yca
    YCb = 1d0/(CI*omega*C2(i) + 1d0/R2p(i))
    Yb = R2(i) + YCb
    X1 = 1d0/(1d0/X1+1d0/Ya)
    X2 = 1d0/(1d0/X2+1d0/Yb)
  enddo

  X12r = X1/(X1+X2)

  select case(re_im)
  case(0)
    zr = abs(X12r)
  case(1)
    zr = Real(X12r)
  case(2)
    zr = Imag(X12r)
  case default
    write(*,*)"Error in zxratio: re_im = ",re_im," is not valid, return abs value as default"
    zr = abs(X12r)
  end select


end function zxratio


 end function th_zxratio
