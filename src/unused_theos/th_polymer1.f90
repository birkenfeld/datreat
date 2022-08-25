      FUNCTION th_polym1 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                        
!     Conherent+incoherent contribution from linear polymer (currently with approximated Rouse)                                                                  
!       
      use Sample_scattering

      implicit none    
      real :: th_polym1                                                            
      CHARACTER(8) thnam, parnam (20) 
      real    :: x, pa, qq
      integer :: ier, nparx, npar, ini
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) SQ_rouseT 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL qget, tget 
      REAL kbolz 

      double precision  :: scqt, scq0, siqt, siq0, diff, Fx, q, Dcm
                                                                        
     double precision, parameter ::  zpi = 6.283185d0 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'polym1'
         nparx = 14
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_polym1 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'Mw'        ! molecular weight  (g/mol)
         parnam (3) = 'h_volfr' 
         parnam (4) = 'h-atoms'
         parnam (5) = 'c-atoms'
         parnam (6) = 'o-atoms'
         parnam (7) = 'n-atoms'
         parnam (8) = 'si-atoms'
         parnam (9) = 'l0'         ! effective bond length (b) in Angstroem
         parnam (10)= 'density'    ! density in g/cm**3                                                                     
         parnam (11)= 'only_coh'                                                                        
         parnam (12)= 'Wl4'        ! Wl4 in A ns units                                                                        
         parnam (13)= 'dtube'      ! tube diameter (not yet implemented)                                                                  
         parnam (14)= 'diff'       ! in cm**2/s                                                                 
!                                                                       
         th_polym1 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----    

      Polymer_name      = 'linear polymer'
      Polymer_Mw        = abs(pa(2))/1000.0   ! molecular weight  (kg/mol) 
      Polymer_volfrac_h = abs(pa(3))          ! volume fraction of h-chains
      Polymer_no_h      = nint(pa(4))         ! number of hydrogens per segmental unit
      Polymer_no_C      = nint(pa(5))         !   "       carbons     "    "       "
      Polymer_no_O      = nint(pa(6))         !   "       oxygens     "    "       "
      Polymer_no_N      = nint(pa(7))         !   "       nitrogens   "    "       "
      Polymer_no_Si     = nint(pa(8))         !   "       silicons    "    "       "
      Polymer_l0        = abs(pa(9))*1d-10    ! effective bond length (b)
      Polymer_density   = abs(pa(10))/1000.0  ! density
     
      only_coherent     = (nint(pa(11)) > 0) 
 
      Polymer_Wl4       = abs(pa(12))*1d-31   ! Rouse rate in m**4/sec 
      Polymer_dtube     = abs(pa(13))*1d-10   ! tube diameter
      diff              = abs (pa (14) ) * 1d16/1d9  ! conversion to A**2/ns

                                  
      tau      = x 
      a0       = pa (1) 
 

      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier) 
      Q = qget 
      IF (ier.ne.0) write (6, * ) 'Warning q not found' 
      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier)  
         temp = tget 
      ENDIF 

      Polymer_reftemp = temp               
                                                       
      scq0 = S_coh_Polymer(Q*1d10,0d0)
      scqt = S_coh_Polymer(Q*1d10,tau*1d-9)       ! needed to know CM_Diffusion
      
      if(diff == 0d0) then
        Dcm = 0d0                                 ! use  internal diffusuion calculation (rouse default) 
      else   
        Dcm = -CM_Diffusion * 1d-9 / 1d-20 + diff ! undo internal value and replace by diff  
      endif
      CALL setpar ('diff    ', sngl (diff*1d9/1d16) ,nopar ,params,napar,mbuf, ier)        ! to ext we use cm**2/s 
      CALL setpar ('diffrou ', sngl (CM_Diffusion / 1d-4) ,nopar ,params,napar,mbuf, ier)  ! to ext we use cm**2/s


      siq0 = S_inc_Polymer(Q*1d10,0d0)
      siqt = S_inc_Polymer(Q*1d10,tau*1d-9)

      Fx =  (scqt-(1d0/3d0)*siqt)/(scq0-(1d0/3d0)*siq0)

      th_polym1 = a0*Fx
                                                                  
      RETURN 
      END FUNCTION th_polym1                             


!! 
!! 
!! 
!!   MODULE Sample_scattering
!!  
!!   use PhysicalConstantsPlus
!! 
!!   integer, parameter     :: peo = 1
!!   integer, parameter     :: pe  = 2
!! 
!!   double precision       :: Sample_thickness  = 0.002d0
!!   double precision       :: Sample_area       = 0.0009d0
!! 
!! 
!! 
!! 
!! 
!! ! ------------------------------------------------------------
!! ! Parameters for Polymer, all in SI units
!! ! ------------------------------------------------------------
!!   character(len=80)      :: Polymer_name      = 'Polyethylene'
!!   double precision       :: Polymer_Mw        = 10d0     ! molecular weight  (kg/mol) 
!!   double precision       :: Polymer_volfrac_h = 0.5d0    ! volume fraction of h-chains
!!   integer                :: Polymer_no_h      = 2        ! number of hydrogens per segmental unit
!!   integer                :: Polymer_no_C      = 1        !   "       carbons     "    "       "
!!   integer                :: Polymer_no_O      = 0        !   "       oxygens     "    "       "
!!   integer                :: Polymer_no_N      = 0        !   "       nitrogens   "    "       "
!!   integer                :: Polymer_no_Si     = 0        !   "       silicons    "    "       "
!!   double precision       :: Polymer_l0        = 3.7d-10  ! effective bond length (b)
!!   double precision       :: Polymer_density   = 0.79d3   ! density
!! 
!!  ! dynamics parameters
!!   double precision       :: Polymer_Wl4       = 80000d-31 ! Rouse rate
!!   double precision       :: Polymer_dtube     = 48d-10    ! tube diameter
!!   double precision       :: Polymer_reftemp   = 509.0d0   ! temperature where the data are valid    
!!  
!! 
!! 
!! 
!! 
!! 
!! ! ------------------------------------------------------------
!! ! and for checking purpose some visible variables
!! ! ------------------------------------------------------------
!!   double precision :: Rgx
!!   double precision :: betax  
!!   double precision :: Mmx
!!   double precision :: Nx
!!   double precision :: CM_Diffusion = 0d0
!! 
!! !
!!   logical          :: only_coherent = .false.
!! 
!! 
!!   CONTAINS
!! 
!! ! ---------------------------------------------
!!   subroutine report_polymer_data(output_handle)
!! ! ---------------------------------------------
!!   implicit none
!!   integer, intent(in) :: output_handle
!!   integer             :: io_error   
!! 
!!   double precision    :: trans, sip, P, total, inc_1, coh_1
!!   double precision    :: Q=0.1d0*1d10
!! 
!! 
!! 
!!   write(output_handle,*,iostat=io_error)' ============ Sample(polymer) data: ============='
!!   if(io_error .ne. 0) then
!!     write(6,*)'report_polymer_data: cannot write to: ',output_handle
!!     return
!!   endif
!!   write(output_handle,'(a40,a        )')' Polymer_name                           ',  trim(polymer_name) 
!!   write(output_handle,'(a40,f12.6    )')' molecular weight  (kg/mol)             ',  Polymer_Mw       
!!   write(output_handle,'(a40,f12.6    )')' volume fraction of h-chains            ',  Polymer_volfrac_h
!!   write(output_handle,'(a40,i6       )')' number of hydrogens per segmental unit ',  Polymer_no_h     
!!   write(output_handle,'(a40,i6       )')'   "       carbons     "    "       "   ',  Polymer_no_C     
!!   write(output_handle,'(a40,i6       )')'   "       oxygens     "    "       "   ',  Polymer_no_O     
!!   write(output_handle,'(a40,i6       )')'   "       nitrogens   "    "       "   ',  Polymer_no_N     
!!   write(output_handle,'(a40,i6       )')'   "       silicons    "    "       "   ',  Polymer_no_Si    
!!   write(output_handle,'(a40,e14.7    )')' effective bond length (b)              ',  Polymer_l0       
!!   write(output_handle,'(a40,f12.6    )')' density                                ',  Polymer_density  
!!                                                                                   
!!                                                                                   
!!   write(output_handle,'(a40,e14.7    )')' Rouse rate                             ',  Polymer_Wl4      
!!   write(output_handle,'(a40,e14.7    )')' tube diameter                          ',  Polymer_dtube    
!!   write(output_handle,'(a40,f12.6    )')' temperature where the data are valid   ',  Polymer_reftemp  
!!   write(output_handle,*)
!! 
!! 
!!   if(only_coherent) then
!!     write(6,*)'@@@@@@@@@@@ coherent scattering only @@@@@@@@@@@'
!!   endif
!! 
!! 
!!   write(output_handle,*)
!!   write(output_handle,'(a)')' derived data:'
!! 
!!     trans = Transmission_Polymer(Sample_thickness)
!!     sip   = Sigma_inc_Polymer()*4*Pi/100
!!     P     = Polarisation_inc_Polymer(Sample_thickness)
!!     total = sigma_total_inc_Polymer(Q,Sample_thickness, Sample_area)
!!     inc_1 = sigma_1_inc_Polymer(Q, Sample_thickness, Sample_area)
!!     coh_1= sigma_1_coh_Polymer(Q, Sample_thickness, Sample_area)
!! 
!!     write(output_handle,'(a40,f12.6    )')' sample_thickness/m                     ',Sample_thickness
!!     write(output_handle,'(a40,f12.6    )')' sample_area/m**2                       ',Sample_area
!!     write(output_handle,'(a40,f12.6    )')' transmission                           ',trans
!!     write(output_handle,'(a40,f12.6    )')' SIGMA incoherent                       ',sip
!!     write(output_handle,'(a40,e14.7    )')' total incoherent(+multiple scat.)/m**2 ',total
!!     write(output_handle,'(a40,e14.7    )')' inc_1 incoherent                 /m**2 ',inc_1
!!     write(output_handle,'(a40,f12.6    )')' polarisation of inc.                   ',P
!!     write(output_handle,'(a40,e14.7    )')' coh_1   coherent                 /m**2 ',coh_1
!!     write(output_handle,'(a40,e14.7    )')' Rg                               /m    ',Rgx
!!     write(output_handle,'(a40,e14.7    )')' Dcm                             m**2/s ',CM_Diffusion
!!     write(output_handle,'(a40,f12.6    )')' N-momoner                              ',Nx  
!!     write(output_handle,'(a40,e14.7    )')' Mw-momoner                      kg/mol ',Mmx
!!     write(output_handle,'(a40,e14.7    )')' beta                             /m**2 ',betax
!!   
!!   return
!!  
!!   end subroutine report_polymer_data
!! 
!! 
!! 
!! ! ------------------------------------------
!!   double precision function DebyeFunction(x)
!! ! ------------------------------------------
!!     implicit none
!!     double precision, intent(in) :: x
!!     
!!     if(x > 0d0) then
!!       DebyeFunction = 2d0*((exp(-x)-1d0)+x)/(x**2)
!!     else
!!       DebyeFunction = 1d0
!!     endif
!!     return
!!   end function DebyeFunction
!! 
!! 
!! ! -----------------------------------------------------------------------
!!   double precision function Sigma_coh_Polymer(Q)
!! ! -----------------------------------------------------------------------
!! ! yields the volume normalized coherent scatterinng cross section
!! ! in SI units (i.e. in 1/m)
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector
!! 
!!     double precision              :: number_of_segments
!!     double precision              :: volume_per_segment
!!     double precision              :: scattering_length_density_h
!!     double precision              :: scattering_length_density_d
!!     double precision              :: molweight_per_segment
!!     double precision              :: radius_of_gyration
!!     double precision              :: beta               ! contrast
!! 
!! 
!!     molweight_per_segment = Polymer_no_h    * molar_mass_hydrogen    &
!!                           + Polymer_no_C    * molar_mass_carbon      &
!!                           + Polymer_no_O    * molar_mass_oxygen      &
!!                           + Polymer_no_N    * molar_mass_nitrogen    &
!!                           + Polymer_no_Si   * molar_mass_silicon
!! 
!!     number_of_segments    = Polymer_Mw / molweight_per_segment
!! 
!!     volume_per_segment    = molweight_per_segment / (Polymer_density * Avogadrozahl)
!! 
!!     scattering_length_density_h = ( Polymer_no_h    * b_hydrogen    &
!!                                  + Polymer_no_C    * b_carbon      &
!!                                  + Polymer_no_O    * b_oxygen      &
!!                                  + Polymer_no_N    * b_nitrogen    &
!!                                  + Polymer_no_Si   * b_silicon     & 
!!                                  ) /  volume_per_segment 
!!     scattering_length_density_d = ( Polymer_no_h    * b_deuterium   &
!!                                  + Polymer_no_C    * b_carbon      &
!!                                  + Polymer_no_O    * b_oxygen      &
!!                                  + Polymer_no_N    * b_nitrogen    &
!!                                  + Polymer_no_Si   * b_silicon     & 
!!                                  ) /  volume_per_segment 
!! 
!!     beta = scattering_length_density_d -scattering_length_density_h
!! 
!!     radius_of_gyration = sqrt(number_of_segments/6d0) * Polymer_l0
!! 
!!     
!! 
!!     Sigma_coh_Polymer = Polymer_volfrac_h * (1d0-Polymer_volfrac_h  )           &
!!                            * beta**2 * Polymer_Mw / (Polymer_density * Avogadrozahl)  &
!!                            * DebyeFunction( (radius_of_gyration * Q)**2 )
!! 
!! 
!! ! for checking purpose !
!!      Rgx = radius_of_gyration
!!      Mmx = molweight_per_segment
!!      betax = beta
!!      Nx  = number_of_segments    
!!  
!! !!     write(6,*)'c: phi   ',Polymer_volfrac_h
!! !!     write(6,*)'c: Polymer_Mw    ',Polymer_Mw*1d3
!! !!     write(6,*)'c: beta_h',scattering_length_density_h, scattering_length_density_h*1d-4
!! !!     write(6,*)'c: beta_d',scattering_length_density_d, scattering_length_density_d*1d-4
!! !!     write(6,*)'c: Deby  ',DebyeFunction( (radius_of_gyration * Q)**2 )
!! !!     write(6,*)
!! !!
!! 
!! 
!!      return
!!   end function Sigma_coh_Polymer
!! 
!! 
!! 
!! ! ------------------------------------------------------------------
!!   double precision function Sigma_inc_Polymer()
!! ! ------------------------------------------------------------------
!! ! yields the volume normalized incoherent (1) scattering cross section
!! ! in SI units (i.e. in 1/m)
!!     implicit none
!! 
!!     double precision              :: number_of_segments
!!     double precision              :: volume_per_segment
!!     double precision              :: molweight_per_segment
!! 
!! 
!! 
!!     if(only_coherent) then
!!       Sigma_inc_Polymer = 1d-10
!!       return
!!     endif  
!! 
!!   
!! 
!!     molweight_per_segment = Polymer_no_h    * molar_mass_hydrogen    &
!!                           + Polymer_no_C    * molar_mass_carbon      &
!!                           + Polymer_no_O    * molar_mass_oxygen      &
!!                           + Polymer_no_N    * molar_mass_nitrogen    &
!!                           + Polymer_no_Si   * molar_mass_silicon
!! 
!!     volume_per_segment    = molweight_per_segment / (Polymer_density * Avogadrozahl)
!! 
!! 
!!     Sigma_inc_Polymer = Polymer_no_h/(4d0*Pi*volume_per_segment) *            &
!!           (Polymer_volfrac_h * sigma_inc_hydrogen + (1d0-Polymer_volfrac_h) * sigma_inc_deuterium) 
!!  
!! !    write(6,*)'Volfrac      ',Polymer_volfrac_h
!! !    write(6,*)'molwtseg     ',molweight_per_segment*1d3
!! !    write(6,*)'volseg       ',volume_per_segment*1d30
!! !    write(6,*)'Sigma_inc_Polymer ',Sigma_inc_Polymer, Sigma_inc_Polymer*1d-2
!! !!
!! 
!!     return
!!   end function Sigma_inc_Polymer
!! 
!! ! -----------------------------------------------------------------------------------
!! !  differential scattering cross sections 
!! ! -----------------------------------------------------------------------------------
!! !
!! 
!! 
!! ! --------------------------------------------------------------------------------
!!   double precision function Transmission_Polymer(thickness)
!! ! --------------------------------------------------------------------------------
!! ! yields the transmission (only incoherent scattering extinction effect)
!! ! 
!!     implicit none
!!     double precision, intent(in)  :: thickness          ! sample (slab) thickness
!! 
!! 
!!     Transmission_Polymer = exp(-4d0*Pi*Sigma_inc_Polymer() * thickness)
!!     
!!     return
!! 
!!   end function Transmission_Polymer
!! 
!! 
!! ! ------------------------------------------------------------------------------------
!!   double precision function Polarisation_inc_Polymer(thickness)
!! ! ------------------------------------------------------------------------------------
!! ! yields the estimated polarisation of the total incoherent scattering
!! ! EXPERIMENTA
!! ! 
!!     implicit none
!!     double precision, intent(in)  :: thickness          ! sample (slab) thickness
!! 
!!     Polarisation_inc_Polymer =    &
!!               (1d0 -exp(-(1d0/3d0)*4d0*Pi*Sigma_inc_Polymer() * thickness)) &
!!             / (1d0 -exp(          4d0*Pi*Sigma_inc_Polymer() * thickness))
!! 
!!      return
!!   end function Polarisation_inc_Polymer
!! 
!! 
!! ! --------------------------------------------------------------------------------------------
!!   double precision function sigma_total_inc_Polymer(Q, thickness, area)
!! ! --------------------------------------------------------------------------------------------
!! ! yields the the total differential incoherent scattering
!! ! assuming isotropic scattering, no absorption
!! ! 
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector (here dummy for systematic)
!!     double precision, intent(in)  :: thickness          ! sample (slab) thickness
!!     double precision, intent(in)  :: area               ! sample area
!! 
!!     sigma_total_inc_Polymer =  (area/(4d0*Pi)) * (1d0-Transmission_Polymer(thickness)) 
!! 
!!      return
!!   end function sigma_total_inc_Polymer
!! 
!! 
!! ! ----------------------------------------------------------------------------------------
!!   double precision function sigma_1_inc_Polymer(Q, thickness, area)
!! ! ----------------------------------------------------------------------------------------
!! ! yields the the direct differential incoherent scattering
!! ! assuming isotropic scattering, no absorption
!! ! 
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector (here dummy for systematic)
!!     double precision, intent(in)  :: thickness          ! sample (slab) thickness
!!     double precision, intent(in)  :: area               ! sample area
!! 
!! !    double precision              :: Transmission_Polymer
!! !    double precision              :: Sigma_inc_Polymer
!! 
!!     sigma_1_inc_Polymer     =  area * thickness* Transmission_Polymer(thickness) &
!!                                  *  Sigma_inc_Polymer()
!! 
!!      return
!!   end function sigma_1_inc_Polymer
!! 
!! 
!! ! ----------------------------------------------------------------------------------------
!!   double precision function sigma_1_coh_Polymer(Q, thickness, area)
!! ! ----------------------------------------------------------------------------------------
!! ! yields the the differential coherent scattering
!! ! 
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector 
!!     double precision, intent(in)  :: thickness          ! sample (slab) thickness
!!     double precision, intent(in)  :: area               ! sample area
!! 
!! 
!!     sigma_1_coh_Polymer =  area * thickness* Transmission_Polymer(thickness) &
!!                                  *  Sigma_coh_Polymer(Q)
!! 
!!      return
!!   end function sigma_1_coh_Polymer
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!  time functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  
!!   double precision function S_coh_Polymer(Q,t)
!! ! --------------------------------------------
!! ! coherent time function
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector 
!!     double precision, intent(in)  :: t                  ! Fouriertime
!! 
!!     double precision              :: molweight_per_segment
!!     double precision              :: number_of_segments
!!     double precision              :: SQ
!!     double precision              :: X                  ! aux variable (time)
!!     double precision              :: Dcm                ! centre-of-mass diffusion               
!!     double precision, parameter   :: a1=0.5871d0, a2=0.3282d0, a3=0.08475d0   !! approximation parameters for ROUSE F(x)
!!     double precision, parameter   :: x1=4.112d0,  x2=1.780d0,  x3=0.5282d0     
!!     double precision, parameter   :: b1=0.7484d0, b2=0.8973d0, b3=1.0d0
!!     double precision              :: Fx
!! 
!! 
!!     molweight_per_segment = Polymer_no_h    * molar_mass_hydrogen    &
!!                           + Polymer_no_C    * molar_mass_carbon      &
!!                           + Polymer_no_O    * molar_mass_oxygen      &
!!                           + Polymer_no_N    * molar_mass_nitrogen    &
!!                           + Polymer_no_Si   * molar_mass_silicon
!! 
!!     number_of_segments    = Polymer_Mw / molweight_per_segment
!! 
!! !   just for the simple start approximation to Rouse 
!!     X   = Polymer_Wl4 * (Q**4) * t / 36d0
!!     Dcm = Polymer_Wl4/(3d0*Polymer_l0**2)/number_of_segments
!!     if(Dcm.ne.0d0) CM_Diffusion = Dcm
!!     SQ  = Sigma_coh_Polymer(Q)
!!  
!!     if( t > 0d0 ) then
!!       Fx = a1*exp(-(X/x1)**b1) +  a2*exp(-(X/x2)**b2) + a3*exp(-(X/x3)**b3)  
!! !      S_coh_Polymer = SQ * exp(-(X/2.7d0)**(0.7d0)) * exp(-Q*Q*Dcm*t)
!!       S_coh_Polymer = SQ * Fx * exp(-Q*Q*Dcm*t)
!!     else
!!      S_coh_Polymer = SQ
!!     endif  
!! 
!! ! ---- for future use here we may add a fast component ---------------------------------
!! !   S_coh_Polymer = S_coh_Polymer + Sigma_fast_coh_Polymer(Q) * exp(-(t/tau_fast_coh_Polymer(Q))    
!! ! --------------------------------------------------------------------------------------
!! 
!!     return
!! 
!!    end function S_coh_Polymer
!! 
!! ! and the resulting normalized time function
!! 
!! !  --------------------------------------------
!!    double precision function F_coh_Polymer(Q,t)
!! !  --------------------------------------------
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector 
!!     double precision, intent(in)  :: t                  ! Fouriertime
!! 
!!     F_coh_Polymer = S_coh_Polymer(Q,t) / S_coh_Polymer(Q, 0d0)
!! 
!!    end function F_coh_Polymer
!! 
!! 
!! 
!! ! ------------------------------------------------------------------------------------------
!! ! --- and the incoherent part --------------------------------------------------------------
!! ! ------------------------------------------------------------------------------------------
!! 
!! 
!! ! --------------------------------------------
!!   double precision function S_inc_Polymer(Q,t)
!! ! --------------------------------------------
!! ! incoherent time function
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector 
!!     double precision, intent(in)  :: t                  ! Fouriertime
!! 
!!     double precision              :: molweight_per_segment
!!     double precision              :: number_of_segments
!!     double precision              :: SQ
!!     double precision              :: X                  ! aux variable (time)
!!     double precision              :: Dcm                ! centre-of-mass diffusion               
!! 
!! 
!!     molweight_per_segment = Polymer_no_h    * molar_mass_hydrogen    &
!!                           + Polymer_no_C    * molar_mass_carbon      &
!!                           + Polymer_no_O    * molar_mass_oxygen      &
!!                           + Polymer_no_N    * molar_mass_nitrogen    &
!!                           + Polymer_no_Si   * molar_mass_silicon
!! 
!!     number_of_segments    = Polymer_Mw / molweight_per_segment
!! 
!! !   just for the simple start approximation to Rouse 
!!     X   = sqrt(Polymer_Wl4 * (Q**4) * t / (9d0 * Pi))
!!     Dcm = Polymer_Wl4/(3d0*Polymer_l0**2)/number_of_segments
!!     SQ  = Sigma_inc_Polymer()
!!  
!!     if( t > 0d0 ) then
!!       S_inc_Polymer = SQ * exp(-X) * exp(-Q*Q*Dcm*t)
!!     else
!!       S_inc_Polymer = SQ
!!     endif  
!! 
!! ! ---- for future use here we may add a fast component ---------------------------------
!! !   S_inc_Polymer = S_inc_Polymer + Sigma_fast_inc_Polymer(Q) * exp(-(t/tau_fast_inc_Polymer(Q))    
!! ! --------------------------------------------------------------------------------------
!! 
!!     if(Dcm.ne.0d0) CM_Diffusion = Dcm
!!     return
!! 
!!    end function S_inc_Polymer
!! 
!! ! and the resulting normalized time function
!! 
!! !  --------------------------------------------
!!    double precision function F_inc_Polymer(Q,t)
!! !  --------------------------------------------
!!     implicit none
!!     double precision, intent(in)  :: Q                  ! wavevector 
!!     double precision, intent(in)  :: t                  ! Fouriertime
!! 
!!     F_inc_Polymer = S_inc_Polymer(Q,t) / S_inc_Polymer(Q, 0d0)
!! 
!!    end function F_inc_Polymer
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!    subroutine set_polymer(polymer_id, Mw, volfrac_h, Wl4)
!! !  ------------------------------------------------------
!!    implicit none
!! 
!!    integer,          intent(in)  :: polymer_id
!!    double precision, intent(in)  :: Mw
!!    double precision, intent(in)  :: volfrac_h
!!    double precision, intent(in)  :: Wl4
!! 
!!    select case(polymer_id)
!!      case(peo)
!!        Polymer_name      = 'PEO'
!!        Polymer_Mw        = Mw          ! molecular weight  (kg/mol) 
!!        Polymer_volfrac_h = volfrac_h   ! volume fraction of h-chains
!!        Polymer_no_h      = 4           ! number of hydrogens per segmental unit 4
!!        Polymer_no_C      = 2           !   "       carbons     "    "       "
!!        Polymer_no_O      = 1           !   "       oxygens     "    "       "
!!        Polymer_no_N      = 0           !   "       nitrogens   "    "       "
!!        Polymer_no_Si     = 0           !   "       silicons    "    "       "
!!        Polymer_l0        = 7.6d-10     ! effective bond length (b)   ??
!!        Polymer_density   = 1.100d3     ! density                     ??
!!     
!!       
!!        if(Wl4 > 0d0) then
!!            Polymer_Wl4       = Wl4        ! Rouse rate                 
!!        else
!!            Polymer_Wl4       =  30000d-31 ! ??
!!        endif
!!        
!!        Polymer_dtube     = 60d-10    ! tube diameter              ??
!!        Polymer_reftemp   = 400.0d0   ! temperature where the data are valid 
!! 
!!     case (pe)
!!        Polymer_name      = 'Polyethylene'                                      
!!        Polymer_Mw        = Mw           ! molecular weight  (kg/mol)               
!!        Polymer_volfrac_h = volfrac_h    ! volume fraction of h-chains              
!!        Polymer_no_h      = 2        ! number of hydrogens per segmental unit   
!!        Polymer_no_C      = 1        !   "       carbons     "    "       "     
!!        Polymer_no_O      = 0        !   "       oxygens     "    "       "     
!!        Polymer_no_N      = 0        !   "       nitrogens   "    "       "     
!!        Polymer_no_Si     = 0        !   "       silicons    "    "       "     
!!        Polymer_l0        = 3.7d-10  ! effective bond length (b)                
!!        Polymer_density   = 0.79d3   ! density                                  
!!                                                                                                          
!!        if(Wl4 > 0d0) then
!!            Polymer_Wl4       = Wl4        ! Rouse rate                 
!!        else
!!            Polymer_Wl4       =  80000d-31 ! ??
!!        endif
!!                                                                                    
!!        Polymer_dtube     = 48d-10    ! tube diameter                           
!!        Polymer_reftemp   = 509.0d0   ! temperature where the data are valid    
!! 
!!     case default
!!        write(6,*)'unknonw polymer  id: ',polymer_id
!!     end select
!! 
!!        write(6,'(a,a,a,f12.6,a,f12.6)')'Polymer is now: ', trim(Polymer_name),' Mw=', Polymer_Mw, &
!!                                      ' kg/mol  volumefraction(H)= ', Polymer_volfrac_h    
!! 
!!    end subroutine set_polymer
!! 
!! 
!! 
!! 
!! 
!! 
!!   END MODULE Sample_scattering
