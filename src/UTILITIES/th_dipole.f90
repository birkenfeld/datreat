 FUNCTION th_dipole(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  dipole field if dipole position is x
!  dipole
      use theory_description 
      implicit none 
      real    :: th_dipole
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
     double precision :: ypos       ! y-position                                                                      
     double precision :: zpos       ! z-position                                                                      
     double precision :: mx         ! moment_x                                                                        
     double precision :: my         ! moment_y                                                                        
     double precision :: mz         ! moment_z                                                                        
! the recin parameter representation 
     double precision :: bcompone   ! component of field                                                              
     double precision :: xsense     ! sensor coordinate x                                                             
     double precision :: ysense     ! sensor coordinate y                                                             
     double precision :: zsense     ! sensor coordinate z                                                             
! the reout parameter representation 
     double precision :: mx         ! moment_x                                                                        
     double precision :: my         ! moment_y                                                                        
     double precision :: mz         ! moment_z                                                                        
     double precision :: ypos       ! y-position                                                                      
     double precision :: zpos       ! z-position                                                                      
 
     double precision :: th
 
   double precision :: xpos, rdipole(3), rauf(3), Bfeld(3), momemt(3)
   integer          :: bcomp
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'dipole'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_dipole = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " dipole field if dipole position is x"
       th_citation(idesc)     = " dipole"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'ypos    '  ! y-position                                                                      
        parnam ( 3) = 'zpos    '  ! z-position                                                                      
        parnam ( 4) = 'mx      '  ! moment_x                                                                        
        parnam ( 5) = 'my      '  ! moment_y                                                                        
        parnam ( 6) = 'mz      '  ! moment_z                                                                        
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "y-position" !//cr//parspace//&
        th_param_desc( 3,idesc) = "z-position" !//cr//parspace//&
        th_param_desc( 4,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc( 5,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc( 6,idesc) = "moment_z" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "bcompone > component of field"
        th_file_param(  2,idesc) = "xsense   > sensor coordinate x"
        th_file_param(  3,idesc) = "ysense   > sensor coordinate y"
        th_file_param(  4,idesc) = "zsense   > sensor coordinate z"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "mx       > moment_x"
        th_out_param(  2,idesc) = "my       > moment_y"
        th_out_param(  3,idesc) = "mz       > moment_z"
        th_out_param(  4,idesc) = "ypos     > y-position"
        th_out_param(  5,idesc) = "zpos     > z-position"
! 
        th_dipole = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      ypos     =      pa( 2)
      zpos     =      pa( 3)
      mx       =      pa( 4)
      my       =      pa( 5)
      mz       =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: component of field
      xh =     1
      call parget('bcompone',xh,iadda,ier)
      bcompone = xh
! >>> extract: sensor coordinate x
      xh =  19.0
      call parget('xsense  ',xh,iadda,ier)
      xsense   = xh
! >>> extract: sensor coordinate y
      xh =   7.
      call parget('ysense  ',xh,iadda,ier)
      ysense   = xh
! >>> extract: sensor coordinate z
      xh =   0.
      call parget('zsense  ',xh,iadda,ier)
      zsense   = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     xpos     = x
     bcomp    =  nint(bcomponen)
     rauf     =  [xsense, ysense, zsense]
     rdipole  =  [ x, ypos, zpos]
     moment   =  [ mx, my, mz ]
     Bfeld    =  bdipol( rauf, rdipole, moment )
     th  = ampli * Bfeld(bcomp)


     th_dipole = th
 
! ---- writing computed parameters to the record >>>  
      call parset('mx      ',sngl(mx),iadda,ier)
      call parset('my      ',sngl(my),iadda,ier)
      call parset('mz      ',sngl(mz),iadda,ier)
      call parset('ypos    ',sngl(ypos),iadda,ier)
      call parset('zpos    ',sngl(zpos),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

!***********************************************************************
!* bdipol  Magnetfeld eines Dipols m am Ort  x
!***********************************************************************

function bdipol( r, x, m ) result (B)
  implicit none
  double precision, intent(in)    :: r(3)       ! Aufpunkt
  double precision, intent(in)    :: x(3)       ! Dipolzentrum
  double precision, intent(in)    :: m(3)       ! Richtung des Einheitsdipols
  double precision                :: B(3)

  double precision, dimension(3)  :: rs, mm
  double precision                :: rr, m0

  double precision, parameter     :: mu_4pi = 1d-7

  rs  = r - x
  rr  = sqrt(dot_product(rs,rs))
!  m0  = sqrt(dot_product(m,m))       ! to normalize m , this may or may not be desired ..
!  if(m0 == 0d0) then
!     B = 0
!     return
!  endif

  mm = m

  B = (mu_4pi/rr**5) * (3*rs*dot_product(mm,rs)-mm*rr**2)

end function bdipol




 end function th_dipole
