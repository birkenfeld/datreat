 FUNCTION th_dipolep(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  dipole field if dipole position is x +pol
!  dipole
      use theory_description 
      implicit none 
      real    :: th_dipolep
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
     double precision :: xx1        !                                                                                 
     double precision :: mx1        ! moment_x                                                                        
     double precision :: my1        ! moment_y                                                                        
     double precision :: mz1        ! moment_z                                                                        
     double precision :: xx2        !                                                                                 
     double precision :: mx2        ! moment_x                                                                        
     double precision :: my2        ! moment_y                                                                        
     double precision :: mz2        ! moment_z                                                                        
     double precision :: xx3        !                                                                                 
     double precision :: mx3        ! moment_x                                                                        
     double precision :: my3        ! moment_y                                                                        
     double precision :: mz3        ! moment_z                                                                        
     double precision :: xx4        !                                                                                 
     double precision :: mx4        ! moment_x                                                                        
     double precision :: my4        ! moment_y                                                                        
     double precision :: mz4        ! moment_z                                                                        
     double precision :: xx5        !                                                                                 
     double precision :: mx5        ! moment_x                                                                        
     double precision :: my5        ! moment_y                                                                        
     double precision :: mz5        ! moment_z                                                                        
     double precision :: px1        ! polarizable dipole1 pos x                                                       
     double precision :: px2        ! polarizable dipole2 pos x                                                       
     double precision :: px3        ! polarizable dipole3 pos x                                                       
     double precision :: pp         ! polarizability factor                                                           
     double precision :: p11        ! polmatrix                                                                       
     double precision :: p12        !                                                                                 
     double precision :: p13        !                                                                                 
     double precision :: p21        !                                                                                 
     double precision :: p22        !                                                                                 
     double precision :: p23        !                                                                                 
     double precision :: p31        !                                                                                 
     double precision :: p32        !                                                                                 
     double precision :: p33        !                                                                                 
! the recin parameter representation 
     double precision :: bcompone   ! component of field                                                              
     double precision :: xsense     ! sensor coordinate x                                                             
     double precision :: ysense     ! sensor coordinate y                                                             
     double precision :: zsense     ! sensor coordinate z                                                             
! the reout parameter representation 
 
 
     double precision :: th
 
   double precision :: xpos, rdipole(3), rauf(3), Bfeld(3), moment(3), P(3,3), BfeldP(3)
   integer          :: bcomp
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'dipolep'
       nparx =       36
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_dipolep = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " dipole field if dipole position is x +pol"
       th_citation(idesc)     = " dipole"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'ypos    '  ! y-position                                                                      
        parnam ( 3) = 'zpos    '  ! z-position                                                                      
        parnam ( 4) = 'xx1     '  !                                                                                 
        parnam ( 5) = 'mx1     '  ! moment_x                                                                        
        parnam ( 6) = 'my1     '  ! moment_y                                                                        
        parnam ( 7) = 'mz1     '  ! moment_z                                                                        
        parnam ( 8) = 'xx2     '  !                                                                                 
        parnam ( 9) = 'mx2     '  ! moment_x                                                                        
        parnam (10) = 'my2     '  ! moment_y                                                                        
        parnam (11) = 'mz2     '  ! moment_z                                                                        
        parnam (12) = 'xx3     '  !                                                                                 
        parnam (13) = 'mx3     '  ! moment_x                                                                        
        parnam (14) = 'my3     '  ! moment_y                                                                        
        parnam (15) = 'mz3     '  ! moment_z                                                                        
        parnam (16) = 'xx4     '  !                                                                                 
        parnam (17) = 'mx4     '  ! moment_x                                                                        
        parnam (18) = 'my4     '  ! moment_y                                                                        
        parnam (19) = 'mz4     '  ! moment_z                                                                        
        parnam (20) = 'xx5     '  !                                                                                 
        parnam (21) = 'mx5     '  ! moment_x                                                                        
        parnam (22) = 'my5     '  ! moment_y                                                                        
        parnam (23) = 'mz5     '  ! moment_z                                                                        
        parnam (24) = 'px1     '  ! polarizable dipole1 pos x                                                       
        parnam (25) = 'px2     '  ! polarizable dipole2 pos x                                                       
        parnam (26) = 'px3     '  ! polarizable dipole3 pos x                                                       
        parnam (27) = 'pp      '  ! polarizability factor                                                           
        parnam (28) = 'p11     '  ! polmatrix                                                                       
        parnam (29) = 'p12     '  !                                                                                 
        parnam (30) = 'p13     '  !                                                                                 
        parnam (31) = 'p21     '  !                                                                                 
        parnam (32) = 'p22     '  !                                                                                 
        parnam (33) = 'p23     '  !                                                                                 
        parnam (34) = 'p31     '  !                                                                                 
        parnam (35) = 'p32     '  !                                                                                 
        parnam (36) = 'p33     '  !                                                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "y-position" !//cr//parspace//&
        th_param_desc( 3,idesc) = "z-position" !//cr//parspace//&
        th_param_desc( 4,idesc) = "" !//cr//parspace//&
        th_param_desc( 5,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc( 6,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc( 7,idesc) = "moment_z" !//cr//parspace//&
        th_param_desc( 8,idesc) = "" !//cr//parspace//&
        th_param_desc( 9,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc(10,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc(11,idesc) = "moment_z" !//cr//parspace//&
        th_param_desc(12,idesc) = "" !//cr//parspace//&
        th_param_desc(13,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc(14,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc(15,idesc) = "moment_z" !//cr//parspace//&
        th_param_desc(16,idesc) = "" !//cr//parspace//&
        th_param_desc(17,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc(18,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc(19,idesc) = "moment_z" !//cr//parspace//&
        th_param_desc(20,idesc) = "" !//cr//parspace//&
        th_param_desc(21,idesc) = "moment_x" !//cr//parspace//&
        th_param_desc(22,idesc) = "moment_y" !//cr//parspace//&
        th_param_desc(23,idesc) = "moment_z" !//cr//parspace//&
        th_param_desc(24,idesc) = "polarizable dipole1 pos x" !//cr//parspace//&
        th_param_desc(25,idesc) = "polarizable dipole2 pos x" !//cr//parspace//&
        th_param_desc(26,idesc) = "polarizable dipole3 pos x" !//cr//parspace//&
        th_param_desc(27,idesc) = "polarizability factor" !//cr//parspace//&
        th_param_desc(28,idesc) = "polmatrix" !//cr//parspace//&
        th_param_desc(29,idesc) = "" !//cr//parspace//&
        th_param_desc(30,idesc) = "" !//cr//parspace//&
        th_param_desc(31,idesc) = "" !//cr//parspace//&
        th_param_desc(32,idesc) = "" !//cr//parspace//&
        th_param_desc(33,idesc) = "" !//cr//parspace//&
        th_param_desc(34,idesc) = "" !//cr//parspace//&
        th_param_desc(35,idesc) = "" !//cr//parspace//&
        th_param_desc(36,idesc) = "" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "bcompone > component of field"
        th_file_param(  2,idesc) = "xsense   > sensor coordinate x"
        th_file_param(  3,idesc) = "ysense   > sensor coordinate y"
        th_file_param(  4,idesc) = "zsense   > sensor coordinate z"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "ypos     > y-position"
        th_out_param(  2,idesc) = "zpos     > z-position"
! 
        th_dipolep = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      ypos     =      pa( 2)
      zpos     =      pa( 3)
      xx1      =      pa( 4)
      mx1      =      pa( 5)
      my1      =      pa( 6)
      mz1      =      pa( 7)
      xx2      =      pa( 8)
      mx2      =      pa( 9)
      my2      =      pa(10)
      mz2      =      pa(11)
      xx3      =      pa(12)
      mx3      =      pa(13)
      my3      =      pa(14)
      mz3      =      pa(15)
      xx4      =      pa(16)
      mx4      =      pa(17)
      my4      =      pa(18)
      mz4      =      pa(19)
      xx5      =      pa(20)
      mx5      =      pa(21)
      my5      =      pa(22)
      mz5      =      pa(23)
      px1      =      pa(24)
      px2      =      pa(25)
      px3      =      pa(26)
      pp       =      pa(27)
      p11      =      pa(28)
      p12      =      pa(29)
      p13      =      pa(30)
      p21      =      pa(31)
      p22      =      pa(32)
      p23      =      pa(33)
      p31      =      pa(34)
      p32      =      pa(35)
      p33      =      pa(36)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: component of field
      xh =    1
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
     bcomp    =  nint(bcompone)
     rauf     =  [xsense, ysense, zsense]
     Bfeld    = 0
     Bfeld    = bdipol( rauf,[ xpos, ypos+xx1, zpos] ,[ mx1, my1, mz1 ]) &
              + bdipol( rauf,[ xpos, ypos+xx2, zpos] ,[ mx2, my2, mz2 ]) &
              + bdipol( rauf,[ xpos, ypos+xx3, zpos] ,[ mx3, my3, mz3 ]) &
              + bdipol( rauf,[ xpos, ypos+xx4, zpos] ,[ mx4, my4, mz4 ]) &
              + bdipol( rauf,[ xpos, ypos+xx5, zpos] ,[ mx5, my5, mz5 ])

     Bfeld    = ampli * Bfeld

     BfeldP   = bdipol([ px1, px2, px3],[ xpos, ypos+xx1, zpos] ,[ mx1, my1, mz1 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx2, zpos] ,[ mx2, my2, mz2 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx3, zpos] ,[ mx3, my3, mz3 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx4, zpos] ,[ mx4, my4, mz4 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx5, zpos] ,[ mx5, my5, mz5 ])

     BfeldP    = ampli * BfeldP



     P(1,1:3) = [p11,p12,p13]
     P(2,1:3) = [p21,p22,p23]
     P(3,1:3) = [p31,p32,p33]

     P        = P * pp
     moment   = matmul(P,BfeldP )

     Bfeld  =  Bfeld + bdipol( rauf,[ px1, px2, px3] , moment)

     th = Bfeld(bcomp)



     xpos     = 50d0
     bcomp    =  nint(bcompone)
     rauf     =  [xsense, ysense, zsense]
     Bfeld    = 0
     Bfeld    = bdipol( rauf,[ xpos, ypos+xx1, zpos] ,[ mx1, my1, mz1 ]) &
              + bdipol( rauf,[ xpos, ypos+xx2, zpos] ,[ mx2, my2, mz2 ]) &
              + bdipol( rauf,[ xpos, ypos+xx3, zpos] ,[ mx3, my3, mz3 ]) &
              + bdipol( rauf,[ xpos, ypos+xx4, zpos] ,[ mx4, my4, mz4 ]) &
              + bdipol( rauf,[ xpos, ypos+xx5, zpos] ,[ mx5, my5, mz5 ])

     Bfeld    = ampli * Bfeld

     BfeldP   = bdipol([ px1, px2, px3],[ xpos, ypos+xx1, zpos] ,[ mx1, my1, mz1 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx2, zpos] ,[ mx2, my2, mz2 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx3, zpos] ,[ mx3, my3, mz3 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx4, zpos] ,[ mx4, my4, mz4 ]) &
              + bdipol([ px1, px2, px3],[ xpos, ypos+xx5, zpos] ,[ mx5, my5, mz5 ])

     BfeldP    = ampli * BfeldP



     P(1,1:3) = [p11,p12,p13]
     P(2,1:3) = [p21,p22,p23]
     P(3,1:3) = [p31,p32,p33]

     P        = P * pp
     moment   = matmul(P,BfeldP )

     Bfeld  =  Bfeld + bdipol( rauf,[ px1, px2, px3] , moment)


     th_dipolep = th -Bfeld(bcomp)
 
! ---- writing computed parameters to the record >>>  
      call parset('ypos    ',sngl(ypos),iadda)
      call parset('zpos    ',sngl(zpos),iadda)
 
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




 end function th_dipolep
