      FUNCTION th_diffonsphere (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
! anisotropic diffusionincoherent with PFG NMR options 
! ref. Hall P-L. and Ross D.K. , Molecular Physics , 1978, p 1549-1554
!    
      implicit none                                                
      CHARACTER(8) thnam, parnam (20) 
      real    :: pa (20), qq (3) 
      real    :: th_diffonsphere
      real    :: x, xh
      integer :: mbuf, nparx, ier, ini, npar
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n

      double precision, parameter    :: Pi = 4*atan(1d0)
      double precision :: q, tau, diffusion, radius    ! will be valid in the scope of sqt_ker and used there, keep the names!
      integer          :: lmax
      integer          :: mode


!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'd2sphere' 
         nparx = 7
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_diffonsphere = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplit'  ! std prefactor 
         parnam (2) = 'diff  '  ! rotational diffusion constant
         parnam (3) = 'r'       ! radius of sphere
         parnam (4) = 'tau   '
         parnam (5) = 'qval  '  ! q (if not give as parameter)
         parnam (6) = 'lmax  '  ! summation limit
         parnam (7) = 'mode  '  ! Sqt in neutron or nmr mode
         
                                                                        
         th_diffonsphere = 0.0 
                                                                        
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q   = pa(5) 
      CALL getpar ('q       ', xh,nopar ,params,napar,mbuf, ier)  
      if(ier.eq.0) q = xh

      tau = abs(pa(4))
      CALL getpar ('tau     ', xh,nopar ,params,napar,mbuf, ier)  
      if(ier.eq.0) tau = xh
 

      diffusion   = abs(pa(2))/(pa(3)**2)/2
      radius      = abs(pa(3))
      lmax        = nint(pa(6))
      mode        = nint(pa(7))
      
      if(lmax <= 0 ) lmax = 30
      select case(mode)
      case (0)
         tau = x
         th_diffonsphere = pa(1) * sqt_dsp(q,tau,radius,diffusion) 
      case (1)       ! x axis is assumed to be q
         q = x
         th_diffonsphere = pa(1) * sqt_dsp(q,tau,radius,diffusion) 
      case (2)       ! pfg scale x=q*q*tau
         q = sqrt(x / tau)
         th_diffonsphere = pa(1) * sqt_dsp(q,tau,radius,diffusion) 
      case default
        write(6,*)"diffonsphere mode: ",mode, " not supported"
        th_diffonsphere = 0
      end select
       
      

contains

  double precision function sqt_dsp(q,t, R, Dr)
    implicit none
    double precision, intent(in)  :: q
    double precision, intent(in)  :: t
    double precision, intent(in)  :: Dr   ! rotational diffusion constant
    double precision, intent(in)  :: R    ! Radius of the sphere

    double precision, parameter  :: Pi  = 4*atan(1d0)
    double precision             :: eps = 1d-7
    double precision             :: bsjn
!    integer                      :: lmax = 100

    integer :: l

    sqt_dsp = 0
    do l=0,lmax
      sqt_dsp = sqt_dsp + (2*l+1) * bsjn(l,q*R)**2 * exp(-(l*(l+1))*Dr*t)
    enddo

  end function sqt_dsp  



END FUNCTION th_diffonsphere

