 FUNCTION th_beaucage(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  general expression allowing the description of the small angle scattering of a hierarchically structured fractal system
!  G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146
      use theory_description 
      implicit none 
      real    :: th_beaucage
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
     double precision :: g1         ! gaussian(1) amplitude                                                           
     double precision :: rg1        ! radius of gyration of first contribution (Gaussian)                             
     double precision :: b1         ! power law component amplitiude of first contribution                            
     double precision :: pow1       ! exponent of first power law contribution                                        
     double precision :: g2         ! gaussian(2) amplitude                                                           
     double precision :: rg2        ! radius of gyration of 2nd  contribution (Gaussian)                              
     double precision :: b2         ! power law component amplitiude of 2nd  contribution                             
     double precision :: pow2       ! exponent of 2nd  power law contribution                                         
     double precision :: g3         ! gaussian(3) amplitude                                                           
     double precision :: rg3        ! radius of gyration of 3rd contribution (Gaussian)                               
     double precision :: b3         ! power law component amplitiude of 3rd contribution                              
     double precision :: pow3       ! exponent of 3rd power law contribution                                          
     double precision :: g4         ! gaussian(4) amplitude                                                           
     double precision :: rg4        ! radius of gyration of 4th contribution (Gaussian)                               
     double precision :: b4         ! power law component amplitiude of 4th contribution                              
     double precision :: pow4       ! exponent of 4th power law contribution                                          
     double precision :: g5         ! gaussian(5) amplitude                                                           
     double precision :: rg5        ! radius of gyration of 5th contribution (Gaussian)                               
     double precision :: b5         ! power law component amplitiude of 5th contribution                              
     double precision :: pow5       ! exponent of 5th power law contribution                                          
! the recin parameter representation 
     double precision :: thick      ! Dicke der Probe                                                                 
! the reout parameter representation 
 
     integer, parameter :: maxcomp=5
     double precision   :: q
     double precision   :: g(maxcomp), rg(maxcomp+1), b(maxcomp), pow(maxcomp)
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'beaucage'
       nparx =       21
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_beaucage = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " general expression allowing the description of the small angle scattering of a hierarchically structured fractal system"
       th_citation(idesc)     = " G. Beaucage, J. Appl. Cryst. (1996) 29, 134-146"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'g1      '  ! gaussian(1) amplitude                                                           
        parnam ( 3) = 'rg1     '  ! radius of gyration of first contribution (Gaussian)                             
        parnam ( 4) = 'b1      '  ! power law component amplitiude of first contribution                             
        parnam ( 5) = 'pow1    '  ! exponent of first power law contribution                                        
        parnam ( 6) = 'g2      '  ! gaussian(2) amplitude                                                           
        parnam ( 7) = 'rg2     '  ! radius of gyration of 2nd  contribution (Gaussian)                              
        parnam ( 8) = 'b2      '  ! power law component amplitiude of 2nd  contribution                             
        parnam ( 9) = 'pow2    '  ! exponent of 2nd  power law contribution                                         
        parnam (10) = 'g3      '  ! gaussian(3) amplitude                                                           
        parnam (11) = 'rg3     '  ! radius of gyration of 3rd contribution (Gaussian)                               
        parnam (12) = 'b3      '  ! power law component amplitiude of 3rd contribution                              
        parnam (13) = 'pow3    '  ! exponent of 3rd power law contribution                                          
        parnam (14) = 'g4      '  ! gaussian(4) amplitude                                                           
        parnam (15) = 'rg4     '  ! radius of gyration of 4th contribution (Gaussian)                               
        parnam (16) = 'b4      '  ! power law component amplitiude of 4th contribution                              
        parnam (17) = 'pow4    '  ! exponent of 4th power law contribution                                          
        parnam (18) = 'g5      '  ! gaussian(5) amplitude                                                           
        parnam (19) = 'rg5     '  ! radius of gyration of 5th contribution (Gaussian)                               
        parnam (20) = 'b5      '  ! power law component amplitiude of 5th contribution                              
        parnam (21) = 'pow5    '  ! exponent of 5th power law contribution                                          
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "gaussian(1) amplitude" !//cr//parspace//&
        th_param_desc( 3,idesc) = "radius of gyration of first contribution (Gaussian)" !//cr//parspace//&
        th_param_desc( 4,idesc) = "power law component of first contribution " //cr//parspace//&
                                  "amplitude scaled by Eq(7) such that expectation is 1"
        th_param_desc( 5,idesc) = "exponent of first power law contribution" !//cr//parspace//&
        th_param_desc( 6,idesc) = "gaussian(2) amplitude" !//cr//parspace//&
        th_param_desc( 7,idesc) = "radius of gyration of 2nd  contribution (Gaussian)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "power law component amplitiude of 2nd  contribution" //cr//parspace//&
                                  "amplitude scaled by Eq(7) such that expectation is 1"
        th_param_desc( 9,idesc) = "exponent of 2nd  power law contribution" !//cr//parspace//&
        th_param_desc(10,idesc) = "gaussian(3) amplitude" !//cr//parspace//&
        th_param_desc(11,idesc) = "radius of gyration of 3rd contribution (Gaussian)" !//cr//parspace//&
        th_param_desc(12,idesc) = "power law component amplitiude of 3rd contribution" //cr//parspace//&
                                  "amplitude scaled by Eq(7) such that expectation is 1"
        th_param_desc(13,idesc) = "exponent of 3rd power law contribution" !//cr//parspace//&
        th_param_desc(14,idesc) = "gaussian(4) amplitude" !//cr//parspace//&
        th_param_desc(15,idesc) = "radius of gyration of 4th contribution (Gaussian)" !//cr//parspace//&
        th_param_desc(16,idesc) = "power law component amplitiude of 4th contribution" //cr//parspace//&
                                  "amplitude scaled by Eq(7) such that expectation is 1"
        th_param_desc(17,idesc) = "exponent of 4th power law contribution" !//cr//parspace//&
        th_param_desc(18,idesc) = "gaussian(5) amplitude" !//cr//parspace//&
        th_param_desc(19,idesc) = "radius of gyration of 5th contribution (Gaussian)" !//cr//parspace//&
        th_param_desc(20,idesc) = "power law component amplitiude of 5th contribution" //cr//parspace//&
                                  "amplitude scaled by Eq(7) such that expectation is 1"
        th_param_desc(21,idesc) = "exponent of 5th power law contribution" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "thick    > Dicke der Probe"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_beaucage = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      g1       =      pa( 2)
      rg1      =      pa( 3)
      b1       =      pa( 4)
      pow1     =      pa( 5)
      g2       =      pa( 6)
      rg2      =      pa( 7)
      b2       =      pa( 8)
      pow2     =      pa( 9)
      g3       =      pa(10)
      rg3      =      pa(11)
      b3       =      pa(12)
      pow3     =      pa(13)
      g4       =      pa(14)
      rg4      =      pa(15)
      b4       =      pa(16)
      pow4     =      pa(17)
      g5       =      pa(18)
      rg5      =      pa(19)
      b5       =      pa(20)
      pow5     =      pa(21)

      g  = [  g1,   g2,   g3,   g4,   g5]
      rg = [ rg1,  rg2,  rg3,  rg4,  rg5, 0d0]
      b  = [  b1,   b2,   b3,   b4,   b5]
      pow= [pow1, pow2, pow3, pow4, pow5]
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: Dicke der Probe
      xh =  0.1
      call parget('thick   ',xh,iadda,ier)
      thick    = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     th_beaucage = ampli * bc_sum(q,g,rg,b,pow,maxcomp)
 
! ---- writing computed parameters to the record >>>  
!      call parset('xyz  ',sngl(xyz),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function bc_sum(q,g,rg,b,pow,maxcomp) result(val)
    implicit none
    double precision, intent(in) :: q
    double precision, intent(in) :: g  (maxcomp)
    double precision, intent(in) :: rg (maxcomp+1)
    double precision, intent(in) :: b  (maxcomp)
    double precision, intent(in) :: pow(maxcomp)
    integer         , intent(in) :: maxcomp

    double precision :: qstar, bstar, bpow
    double precision :: k = 1.06d0
    double precision :: val
    integer          :: i
   
    val = 0
    do i=1,maxcomp
      qstar = q / (erf( q * k * rg(i)/sqrt(6d0)))**3
      if(rg(i) > 0) then
         bstar = g(i)*pow(i) / rg(i)**pow(i) * Gamma(pow(i)/2)    ! equation (7) of benoits paper
      else
         bstar = 0
      endif
      if(qstar .ne. 0d0) then
        bpow = bstar*b(i) * exp(-((q*rg(i+1))**2)/3d0) *(1d0/qstar)**pow(i)
      else
        bpow = 0
      endif
      val   = val + g(i) * exp(-((q*rg(i))**2)/3d0) + bpow 
    enddo

  end function bc_sum
 end function th_beaucage
