 FUNCTION th_sinhafra(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  fractal scattering assuming a massfractal with dimension df and correlation length xi direct integration for r0 .ne.0
!  J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_sinhafra
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
     double precision :: df         ! fractal dimension                                                               
     double precision :: xi         ! correlation length                                                              
     double precision :: r0         ! correlation hole                                                                
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision   :: adapint, eps=1d-9, erracc
     integer            :: maxit = 1000
     integer            :: nper  = 5
     double precision   :: q
     double precision   :: acc = 1d-10
     double precision   :: g1, g2, gm, dg
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sinhafra'
       nparx =        4
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sinhafra = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " fractal scattering assuming a massfractal with dimension df and correlation length xi direct integration for r0 .ne.0"
       th_citation(idesc)     = " J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'df      '  ! fractal dimension                                                               
        parnam ( 3) = 'xi      '  ! correlation length                                                              
        parnam ( 4) = 'r0      '  ! correlation hole                                                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fractal dimension" !//cr//parspace//&
        th_param_desc( 3,idesc) = "correlation length" !//cr//parspace//&
        th_param_desc( 4,idesc) = "correlation hole" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_sinhafra = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      df       =      abs( pa( 2) )
      xi       =      abs( pa( 3) )
      r0       =      abs( pa( 4) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
!      call parget('xh      ',xh,iadda,ier)
!      xh         = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x

     if(r0 == 0d0) then
         th_sinhafra = GAMMA(-1+df)*sqrt(q**2 *xi**2 +1) * sin( (-1+df)*atan(q*xi)) / ( (q**2 * xi**2+1)**(0.5d0*df) *q*xi ) 
     else
         th_sinhafra = 0d0
         g1 = r0
         gm = -log(acc) * xi
         dg = nper*2*pi/q
         if(dg > gm) dg = gm
d1:      do 
           g2 = g1+dg
           th_sinhafra =  th_sinhafra + adapint(gr_kernel,g1,g2, eps, maxit, erracc) 
           if(g2 > gm) exit d1
           g1 = g2
         enddo d1

     endif

     th_sinhafra = GAMMA(-1+df)*sqrt(q**2 *xi**2 +1) * sin( (-1+df)*atan(q*xi)) / ( (q**2 * xi**2+1)**(0.5d0*df) *q*xi ) 
     write(6,*)"q = ", q, " th_sinhafra0 = ", th_sinhafra
     th_sinhafra =  th_sinhafra + adapint(gr_kernel,r0,gm, eps, maxit, erracc) 
     write(6,*)"q = ", q, " th_sinhafra1 = ", th_sinhafra

     th_sinhafra = ampli * th_sinhafra
 
! ---- writing computed parameters to the record >>>  
      rg = xi * sqrt(df*(1+df))
      call parset('rg      ',sngl(rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function gr_kernel(r) result(y)

    implicit none
    double precision, intent(in) :: r
    double precision             :: y

    y = sin(q*r)* r**(-2d0+df) * exp(-r/xi)/q

  end function gr_kernel

 end function th_sinhafra
