      FUNCTION rousei (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rousei <--------                                              
!                                                                       
!                
      implicit none                                                       
      CHARACTER(len=8), intent(inout) :: thnam, parnam (20) 
      real, intent(in)                :: x,  pa (20)
      integer, intent(in)             :: ini, mbuf
      integer, intent(inout)          :: npar
      integer, intent(inout)          :: nopar            ! Anzahl der Parameter data
      character*80, intent(inout)     :: napar(mbuf)      ! name des parameters n
      real, intent(inout)             :: params(mbuf)     ! value des parameters n
      real                            :: rousei

      double precision :: temp, tau, q 
      double precision :: SQ_rouseTW 
      double precision :: a0 
      double precision :: diff 
      REAL             :: qget, tget 

      double precision  :: rho, com_diff, wl4, c_inf, mw, mmon, contrast, phi_h, b
      double precision  :: scoh0, scoh, sinc0, sinc
      double precision  :: Rg, qrg2, debye
      integer           :: n_h_mon                   ! protons per mononmer
      integer           :: ier, nparx

    double precision, parameter :: Pi                  =     3.141592653589793238462643d0
    
    
    double precision, parameter :: Larmorkonstante     =     2916.46943d+4 
    double precision, parameter :: Mue0                =     1.25663706d-6
    double precision, parameter :: Neutronenmasse      =     1.674927351d-27 
    double precision, parameter :: Planckkonstante     =     6.62606957d-34
    double precision, parameter :: Erdbeschleunigung   =     9.80665d0       
    double precision, parameter :: Boltzmannkonstante  =     1.3806488d-23
    double precision, parameter :: Elektronenladung    =     1.602176565d-19
    double precision, parameter :: Avogadrozahl        =     6.02214129d+23
    double precision, parameter :: Molvolumen_Gas      =     0.02241383d0
    double precision, parameter :: Elektronenmasse     =     9.109534d-31
    double precision, parameter :: amu                 =     1.660538921d-27 
    double precision, parameter :: Epsilon0            =     8.854187818d-12  ! Vakuumdielktrizitaetskonstante
    double precision, parameter :: Lichtgeschwindigkeit=     299792458.0d0
    double precision, parameter :: Gaskonstante        =     8.3144621d0
    
    
    double precision, parameter :: unit_Angstroem      =     1d-10
    double precision, parameter :: unit_cm             =     1d-2
    double precision, parameter :: unit_eV             =     1.602176565d-19
    double precision, parameter :: unit_barn           =     1d-28
    double precision, parameter :: unit_ns             =     1d-9
    
    double precision, parameter :: sigma_inc_hydrogen  = 79.9040d-28
                                                                        
    integer iadda
    common/thiadd/iadda
                                                                   
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'rousei' 
         nparx = 11 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            rousei = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters   units: micro ns A , macro g, cm, mol                
         parnam (1)  = 'amplitud' 
         parnam (2)  = 'wl4     ' 
         parnam (3)  = 'b_segmnt' 
         parnam (4)  = 'c_inf   ' 
         parnam (5)  = 'mw      '
         parnam (6)  = 'mmonom  '
         parnam (7)  = 'n_h_mon '
         parnam (8)  = 'phi_h   '
         parnam (9)  = 'rho     '
         parnam (10) = 'contrast'
         parnam (11) = 'com_diff' 
!                                                                       
         rousei = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau     = x 
      a0      = pa (1) 
      wl4     = abs (pa (2) ) 
      b       = abs (pa (3) ) 
      c_inf   = abs (pa (4) ) 
      mw      = pa (5) 
      mmon    = pa(6)
      n_h_mon = nint(pa(7))
      phi_h   = pa(8)
      rho     = pa(9)
      contrast= pa(10)
      com_diff= pa(11)


                                                         
      diff = abs (com_diff )        ! in cm**2/sec
      diff = diff * 1d-9 / 1d-16   ! in A**s/ns
    
      qget = 0.01 
      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ier)  
      q = qget 
      if(ier.ne.0) then
        write(6,*)"ERROR q-vqlue (q)  missing in parameter list "
        write(6,*)"use putpar q <q>    to add it to the selected data records !"
        rousei = 0
        return
      endif
      tget = 400 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier) 
      temp = tget 
      if(ier.ne.0) then
        write(6,*)"ERROR temperature (temp) missing in parameter list "
        write(6,*)"use putpar temp <T>    to add it to the selected data records !"
        rousei = 0
        return
      endif

                                                                        
                                                                        
      rousei = 0 
       
      Rg =b/sqrt(6d0)*sqrt((c_inf*Mw)/Mmon)

      sinc0 = phi_h*(sigma_inc_hydrogen/unit_cm**2)/(4*Pi)*n_h_mon/Mmon*rho*Avogadrozahl
      sinc  = sinc0 * exp(-sqrt((1d0/9d0)*Q**4*Wl4*tau/Pi))

 
      qrg2  = (q*Rg)**2
      debye = (2d0/qrg2**2)*(exp(-qrg2)-1d0+qrg2)                                                    
      scoh0 = phi_h*(1d0-phi_h)*(contrast)**2*Mw/rho/Avogadrozahl * debye         
      scoh  = scoh0*SQ_rouseTW (tau, q, Wl4)    
 
                                                            
!         IF (pa (4) .lt.0) then 
!            WRITE (6, * ) 'Full rouse computation from scratch!' 
!            sum = sum + a * SQ_rouse (tau, qzz, temp, xi, b, 1d-8) 
 
      rousei =  a0 * dexp ( - q*q * diff * tau) * (scoh-sinc/3)/(scoh0-sinc0/3) 
!      

      
      call parset('rg       ',sngl(rg),iadda) 
      call parset('sinc0    ',sngl(sinc0),iadda) 
      call parset('scoh0    ',sngl(scoh0),iadda) 
                                                                 
      RETURN 
      END FUNCTION rousei   
