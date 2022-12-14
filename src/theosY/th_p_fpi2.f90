 FUNCTION th_p_fpi2(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  a * P(Ipi=x,Bx,By,Bz,lambda) Pi/2 flipper polarisation function
! 
      use theory_description 
      implicit none 
      real    :: th_p_fpi2
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
     double precision :: bx         ! longitudinal field in Gauss                                                     
     double precision :: by         ! transversal  field in Gauss                                                     
     double precision :: bz         ! vertical field in Gauss                                                         
     double precision :: turnpm     ! winding density in turns per m                                                  
     double precision :: d          ! thickness in m                                                                  
     double precision :: lambda     ! wavelength in Angstroem                                                        
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'p_fpi2'
       nparx =        7
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_p_fpi2 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " a * P(Ipi=x,Bx,By,Bz,lambda) Pi/2 flipper polarisation function"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'a       '  ! prefactor                                                                       
        parnam ( 2) = 'bx      '  ! longitudinal field in Gauss                                                     
        parnam ( 3) = 'by      '  ! transversal  field in Gauss                                                     
        parnam ( 4) = 'bz      '  ! vertical field in Gauss                                                         
        parnam ( 5) = 'turnpm  '  ! winding density in turns per m                                                  
        parnam ( 6) = 'd       '  ! thickness in m                                                                  
        parnam ( 7) = 'lambda  '  ! wavelength in Angstroem                                                         
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "longitudinal field in Gauss" !//cr//parspace//&
        th_param_desc( 3,idesc) = "transversal  field in Gauss" !//cr//parspace//&
        th_param_desc( 4,idesc) = "vertical field in Gauss" !//cr//parspace//&
        th_param_desc( 5,idesc) = "winding density in turns per m" !//cr//parspace//&
        th_param_desc( 6,idesc) = "thickness in m" !//cr//parspace//&
        th_param_desc( 7,idesc) = "wavelength in Angstroem" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_p_fpi2 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      a        =      pa( 1)
      bx       =      pa( 2)
      by       =      pa( 3)
      bz       =      pa( 4)
      turnpm   =      pa( 5)
      d        =      pa( 6)
      lambda   =      pa( 7)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th =  a * pol_of_fpi2(dble(x))

     th_p_fpi2 = th
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
   function pol_of_fpi2(Ipi) result(P)
     implicit none
     double precision, intent(in) :: Ipi
     double precision             :: P

     double precision, parameter :: Pi                  =     3.141592653589793238462643d0
     double precision, parameter :: Larmorkonstante     =     2916.469323d+4! (older: 2916.46943d+4) -> nist new
     double precision, parameter :: Mue0                =     1.25663706d-6
     double precision, parameter :: Neutronenmasse      =     1.674927351d-27
     double precision, parameter :: Planckkonstante     =     6.62606957d-34
     double precision, parameter :: Gauss_Tesla         =     1d-4
     double precision, parameter :: mu0                 =     4*Pi*1d-7

     double precision :: sigma, gam, lam, pref, Bpi2, arg, b2

     gam   = 2*Pi* Larmorkonstante
     sigma = d * gam * Neutronenmasse / Planckkonstante
     lam   = lambda * 1d-10
     pref  = lam * sigma * Gauss_Tesla
     Bpi2  = Ipi*turnpm*mu0 / Gauss_Tesla
     b2    = bx**2+by**2+bz**2
     arg   = (b2+2*bz*Bpi2+Bpi2**2)
     P     = (cos(pref*sqrt(arg))*(Bpi2**2)*(b2-bz**2) +(b2+bz*Bpi2)**2)/(b2*arg)

   end function pol_of_fpi2





 end function th_p_fpi2
