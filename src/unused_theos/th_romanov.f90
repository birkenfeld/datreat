      FUNCTION th_romanov (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> romanov layer <--------
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
      real                            :: th_romanov
      real                            :: sq, sqt, a0, qget
      integer                         :: i, it, ik

      integer, parameter              :: mtau = 200
      real, save                      :: xvec(mtau), last_xvec(mtau) = 0, last_pa(20) = 0

      integer            :: ier
      integer, save      :: nparx, nxw

      double precision, parameter :: Pi = 4*atan(1d0)
      double precision, parameter :: Boltzmannkonstante  =     1.3806488d-23
    
      double precision, save   :: kip , last_kip = 0d0    !! experimental scattering vector components
      double precision, save   :: kz  , last_kz  = 0d0    !! experimental scattering vector components
    
    
      integer            :: n      = 10            !! number of layers
      double precision   :: B      = 1.0d5         !! compression modulus in SI-units
      double precision   :: gam    = 0.030d0       !! free surface tension in SI-units
      double precision   :: eta3   = 0.002d0       !! layer viscosity in SI-units
      double precision   :: d      = 60d-10        !! layer-layer distance in SI-units
      double precision   :: rho    = 1000d0        !! density (of fluid)
      double precision   :: kappa  = 10d0          !! bending modulus in kT units
      double precision   :: temp   = 300d0         !! temperature in Kelvin
      double precision   :: K      = 0d0           !! (bulk) bending modulus in SI units
    
      
      double precision   :: d_evanescent = 500d-10  !! evanescent waveextension
    
      integer            :: nr = 200                ! number of in plane radius points 
      integer            :: nq = 300                ! number of q integration points
    
      double precision   :: qmax       = 0.3d10     ! maximum for q integration
      double precision   :: rmax       = 2000d-10   ! maximum for r-integration (this is the size of the patterson function carrier)
      double precision   :: rcutfactor = 0.3d0      ! rmax rcutfactor = Gaussian width (plain) of Gaussian Patteron fucntion for Gnm

    
      integer            :: nk = 1
      integer            :: nz = 1
      integer, save      :: nt 


      double precision, save          :: tau(0:mtau) 
      double precision, save          :: skt(1,0:mtau) 
      double precision, save          :: kzv(1:1), kipv(1:1)
    
      logical            :: recompute


      integer          :: iadda 
      common/thiadd/iadda

 
 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'romanov'
         nparx = 15
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_romanov = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'n       '   ! number of layers
         parnam(3) = 'B       '
         parnam(4) = 'gam     '
         parnam(5) = 'eta3    '
         parnam(6) = 'd       '
         parnam(7) = 'rho     '
         parnam(8) = 'kappa   '
         parnam(9) = 'temp    '
         parnam(10)= 'd_evanes'
         parnam(11)= 'nr      '  
         parnam(12)= 'nq      ' 
         parnam(13)= 'qmax    '
         parnam(14)= 'rmax    '
         parnam(15)= 'rcf     '                            
 
         th_romanov =  0
         return      
       endif

        a0              =       pa(1)
        n               =  nint(pa(2))
        B               =   abs(pa(3))
        gam             =   abs(pa(4)) 
        eta3            =   abs(pa(5)) 
        d               =   abs(pa(6))
        rho             =   abs(pa(7)) 
        kappa           =   abs(pa(8))
        temp            =   abs(pa(9))
        d_evanescent    =   abs(pa(10))
        nr              =  nint(pa(11))
        nq              =  nint(pa(12))
        qmax            =   abs(pa(13))
        rmax            =   abs(pa(14))
        rcutfactor      =   abs(pa(15))
 

        call        parget('q       ',qget,iadda,ier)
        kz    = qget            
        call        parget('kip     ',qget,iadda,ier)
        kip   = qget
             

        nt = mtau
        call get_full_xvec(iadda,nt,xvec)
        if(nt <= 0) then
           th_romanov = 0
           return
        endif

!! now decide whether romanov needs a new computation !!
        recompute = .false.

        do i=1,nt   ; if( xvec(i) .ne. last_xvec(i) ) recompute = .true. ; enddo
        if(recompute) last_xvec = xvec

        do i=2,nparx ; if( pa(i)   .ne. last_pa(i)   ) recompute = .true. ; enddo
        if(recompute) last_pa = pa

        if( kz .ne. last_kz )                         recompute = .true.
        if(recompute) last_kz = kz
 
        if( kip .ne. last_kip )                       recompute = .true.
        if(recompute) last_kip = kip

       
      if(recompute) then
          tau(0)    = 0
          tau(1:nt) = xvec(1:nt) * 1d-9
          kzv       =        kz  * 1d10
          kipv      =        kip * 1d10
          K         = Boltzmannkonstante *  temp / d * kappa
          call romanov(n, B, K, d, gam, eta3, temp, rho, d_evanescent, nt, tau, nk, kzv, kipv, & 
                       skt, nq, qmax, nr, rmax, rcutfactor )

          write(6,*) " ============= romanov recomputation ============= "
          write(6,'(a,e16.7)')  "B:=     ",B
          write(6,'(a,e16.7)')  "K:=     ",K
          write(6,'(a,e16.7)')  "gam:=   ",gam
          write(6,'(a,e16.7)')  "T:=     ",temp
          write(6,'(a,e16.7)')  "eta3:=  ",eta3
          write(6,'(a,e16.7)')  "deva:=  ",d_evanescent
          write(6,'(a,e16.7)')  "qmax:=  ",qmax
          write(6,'(a,e16.7)')  "rmax:=  ",rmax
          write(6,'(a,e16.7)')  "rcf:=   ",rcutfactor
          write(6,'(a,i8)')     "n:=     ",n
          write(6,'(a,i8)')     "nq:=    ",nq
          write(6,'(a,i8)')     "nr:=    ",nr
          write(6,'(a,20F12.6)')"kz:=    ",kzv*1d-10
          write(6,'(a,20f12.6)')"kip:=   ",kipv*1d-10
          do it=0,nt
            write(6,'(f12.6,2x,10e14.7)') tau(it)*1d9,  (skt(ik,it)*1d20,ik=1,nk)
          enddo

      endif

dtab: do it=1 ,nt
        if( abs( x - xvec(it) ) <= 1d-5 * abs(x) ) then
          Sqt = skt(1,it)
          Sq  = skt(1,0)
          exit dtab
        endif
      enddo dtab


!  ---- calculate theory here -----


       th_romanov =  a0 * Sqt/Sq
  
       return
       end function th_romanov




!!> program romanov_v0  
!!> !
!!> ! Scattering form layer stack on an interface computed starting from the
!!> ! model exposed in  V. P. Romanov and S. V. Ul’yanov PHYSICAL REVIEW E 66, 061701 (2002)
!!> ! with added computation of the intermediate scattering function as can be observed in
!!> ! NEWSES
!!> !
!!> ! Author: Michael Monkenbusch
!!> !         Juelich Center for Neutron Science
!!> !         Forschungszentrum Juelich GmbH
!!> !         Juelich Germany
!!> ! Email:  m.monkenbusch@fz-juleich.de
!!> !
!!> !
!!> ! Use:    FC -o romanov_v0 romanov_v0.f90 -O2
!!> !         to compile an link, where FC stands for any modern (2003) fortran compiler, e.g. gfortran (V > 4.9)  
!!> !         romanov_v0 help
!!> !         to see possible call parameters
!!> !         romanov_v0 parnam <val> .....
!!> !         to run (takes a few minutes)
!!> !
!!> 
!!> 
!!>   implicit none
!!>   double precision, parameter :: Pi = 4*atan(1d0)
!!>   double precision, parameter :: Boltzmannkonstante  =     1.3806488d-23
!!> 
!!>   double precision   :: kip        = 0.01d10   !! experimental scattering vector components
!!>   double precision   :: kz         = 0.2d10    !! experimental scattering vector components
!!> 
!!> 
!!>   integer            :: n      = 10            !! number of layers
!!>   double precision   :: B      = 1.0d5         !! compression modulus in SI-units
!!>   double precision   :: gam    = 0.030d0       !! free surface tension in SI-units
!!>   double precision   :: eta3   = 0.002d0       !! layer viscosity in SI-units
!!>   double precision   :: d      = 60d-10        !! layer-layer distance in SI-units
!!>   double precision   :: rho    = 1000d0        !! density (of fluid)
!!>   double precision   :: kappa  = 10d0          !! bending modulus in kT units
!!>   double precision   :: temp   = 300d0         !! temperature in Kelvin
!!>   double precision   :: K      = 0d0           !! (bulk) bending modulus in SI units
!!> 
!!>   
!!>   double precision   :: d_evanescent = 500d-10  !! evanescent waveextension
!!> 
!!>   integer            :: nr = 200                ! number of in plane radius points 
!!>   integer            :: nq = 300                ! number of q integration points
!!> 
!!>   double precision   :: qmax       = 0.3d10     ! maximum for q integration
!!>   double precision   :: rmax       = 2000d-10   ! maximum for r-integration (this is the size of the patterson function carrier)
!!>   double precision   :: rcutfactor = 0.3d0      ! rmax rcutfactor = Gaussian width (plain) of Gaussian Patteron fucntion for Gnm
!!> 
!!> 
!!> 
!!>   integer            :: nk = 5
!!>   integer            :: nz = 1
!!>   integer            :: nt = 50
!!>   double precision   :: taumax     = 15d-9
!!>   double precision   :: taumin     = 0.05d-9
!!>   double precision, allocatable   :: tau(:)
!!>   double precision, allocatable   :: skt(:,:)
!!>   double precision, allocatable   :: kzv(:), kipv(:)
!!> 
!!>   Integer            :: it, ik
!!> 
!!>   integer            :: i, length, status
!!>   character(len=256) :: val
!!>   
!!>   do i=1,30
!!> !    write(6,*)'==========',i,'=========='
!!>     call get_command_argument(i,val,length, status)
!!>     if(status /= 0) exit
!!> !    write(6,*)'val:    ',trim(val)
!!> !    write(6,*)'len:    ',length
!!> !    write(6,*)'status: ',status
!!>     if(trim(val) .eq. 'n') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) n
!!>     endif
!!>     if(trim(val) .eq. 'B') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) B
!!>     endif
!!>     if(trim(val) .eq. 'kappa') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) kappa
!!>     endif
!!>     if(trim(val) .eq. 'gam') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) gam 
!!>     endif
!!>     if(trim(val) .eq. 'eta3') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) eta3
!!>     endif
!!>     if(trim(val) .eq. 'd') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) d
!!>     endif
!!>     if(trim(val) .eq. 'K') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) K
!!>     endif
!!>     if(trim(val) .eq. 'T') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) temp
!!>     endif
!!>     if(trim(val) .eq. 'rho') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) rho
!!>     endif
!!>     if(trim(val) .eq. 'kz') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) kz
!!>     endif
!!>     if(trim(val) .eq. 'kip') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) kip
!!>     endif
!!>     if(trim(val) .eq. 'nt') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) nt
!!>     endif
!!>     if(trim(val) .eq. 'nk') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) nk
!!>     endif
!!>     if(trim(val) .eq. 'nz') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) nz
!!>     endif
!!>     if(trim(val) .eq. 'nr') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) nr
!!>     endif
!!>     if(trim(val) .eq. 'rcf') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) rcutfactor
!!>     endif
!!>     if(trim(val) .eq. 'nq') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) nq
!!>     endif
!!>     if(trim(val) .eq. 'qmax') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) qmax
!!>     endif
!!>     if(trim(val) .eq. 'rmax') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) rmax
!!>     endif
!!>     if(trim(val) .eq. 'tmax') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) taumax
!!>     endif
!!>     if(trim(val) .eq. 'tmin') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) taumin
!!>     endif
!!>     if(trim(val) .eq. 'deva') then
!!>        call get_command_argument(i+1,val,length, status)
!!>        read(val,*) d_evanescent
!!>     endif
!!>     if(trim(val) .eq. 'help' .or. trim(val) .eq. '--help') then
!!>        write(6,*)"============================================================================"
!!>        write(6,*)"= romanov dynamic function of layer stack in evanescent field              ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=    romanov_xx nam <val> ...                                              ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=    nam           value                                                   ="
!!>        write(6,*)"=    ---           -----                                                   ="
!!>        write(6,*)"=    B             layer compression modulus (SI)                          ="
!!>        write(6,*)"=    K             bulk bending modulus (SI) if not given use kappa        ="
!!>        write(6,*)"=    kappa         layer bending modulus in kT units                       ="
!!>        write(6,*)"=    gam           surface tension (SI)                                    ="
!!>        write(6,*)"=    eta3          layer viscosity in (SI)                                 ="
!!>        write(6,*)"=    d             layer distance in (SI) (m)                              ="
!!>        write(6,*)"=    n             number of layers                                        ="
!!>        write(6,*)"=    deva          decay lengh of evanescent wave (SI)                     ="
!!>        write(6,*)"=    T             temperature in Kelvin                                   ="
!!>        write(6,*)"=    rho           (fluid) density (SI)                                    ="
!!>        write(6,*)"=    nt            number of log spaced tau points                         ="
!!>        write(6,*)"=    tmin          min tau log scale (if=0 -> lin scale)  in (SI) (s)      ="
!!>        write(6,*)"=    tmax          max tau in (SI) (s)                                     ="
!!>        write(6,*)"=    kz            z-component of q-vector (SI) (m**-1)                    ="
!!>        write(6,*)"=    kip           in plane q-vector component in (SI) (m**-1)             ="
!!>        write(6,*)"=    nk            number of in-plane k-values betwenn 0..kip              ="
!!>        write(6,*)"=    nz            number of kz-values betwenn 0..kz                       ="
!!>        write(6,*)"=    nq            nuber of q-integration points (internal) 200+           ="
!!>        write(6,*)"=    qmax          max internal q-integration range (SI) (1d10)            ="
!!>        write(6,*)"=    nr            number of r integration points (internal) (200)         ="
!!>        write(6,*)"=    rmax          r-range of patterson function (SI) (2000d-10)           ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=    result is written to skt.dgli                                         ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=                                                                          ="
!!>        write(6,*)"=    Author: M. Monkenbusch, JCNS, FZJ                                     ="
!!>        write(6,*)"============================================================================"
!!>        stop
!!>     endif
!!>   enddo
!!>   
!!> 
!!>   allocate(tau(0:nt))
!!>   allocate(skt(max(nk,nz),0:nt))
!!>   allocate(kipv(max(nk,nz)))
!!>   allocate(kzv(max(nk,nz)))
!!> 
!!>   if(nz == 1) then
!!>     kzv(1:nk)   = kz  
!!>     kipv(1:nk)  = kip * [(i,i=1,nk)]/dble(nk)
!!>   elseif( nk == 1) then
!!>     nk = nz
!!>     kzv(1:nk)   = kz  * [(i,i=1,nk)]/dble(nk)
!!>     kipv(1:nk)  = kip 
!!>   else
!!>     stop "either nz or nk must be one..."
!!>   endif  
!!> 
!!> 
!!>   if(K == 0d0) then
!!>       K          = Boltzmannkonstante *  temp / d * kappa
!!>   endif
!!> 
!!>   kappa          = K / (Boltzmannkonstante *  temp / d)
!!> 
!!> !! -----------------------------------------------------------------
!!>   write(6,'(a                )')" ------------------------------------------------------------------------------ "
!!>   write(6,'(a                )')"REF:  V. P. Romanov and S. V. Ul’yanov PHYSICAL REVIEW E 66, 061701 (2002)"
!!>   write(6,'(a                )')"      scattering without contrats and formfactors in 5.2 but with dynamics "
!!>   write(6,'(a                )')" ------------------------------------------------------------------------------ "
!!>   write(6,'(a                )')"Parameters:                                                               "
!!>   write(6,'(a                )')"name  ............................ value ....... SI units to be entered .... conversions.."
!!>   write(6,'(a,e14.7,a,e14.7,a)')"B     (layer compression modulus) ",B    ," [N/m**2] = ",B*1d5/1d4,  " [dyn/cm**2]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"K     (bulk  bending modulus)     ",K    ," [N]      = ",K*1d5,      " [dyn]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"kappa (layer bending modulus/kT)  ",kappa," [kT]     = ",kappa*Boltzmannkonstante*temp*1d7," [erg]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"gam   (surface tension)           ",gam  ," [N/m]    = ",gam*1d5/1d2," [dyn/cm]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"eta3  (layer viscosity)           ",eta3 ," [kg/ms]  = ",eta3*1d3/1d2," [g/cms]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"rho   (fluid density)             ",rho  ," [kg/m**3]= ",rho/1d3," [g/cm**3]"
!!>   write(6,'(a,e14.7,a        )')"temp  (temperature    )           ",temp ," [K]  "
!!>   write(6,'(a,e14.7,a,e14.7,a)')"d     (layer layer distance)      ",d    ," [m]      = ",d*1d10," [Angstroem]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"deva  (evanescent field depth)    ",d_evanescent ," [m]      = ",d_evanescent*1d10," [Angstroem]"
!!>   write(6,'(a,i8             )')"n     (number of layers)          ",n
!!>   write(6,'(a                )')" ------------------------------------------------------------------------------ "
!!>   write(6,'(a                )')"Parameters controlling the internal computation: "
!!>   write(6,'(a,e14.7,a,e14.7,a)')"qmax  (q-integration width,eq 3.9)",qmax ," [m]      = ",qmax*1d-10," [1/Angstroem]"
!!>   write(6,'(a,i8             )')"nq    (number of q intpoints)     ",nq
!!>   write(6,'(a,e14.7,a,e14.7,a)')"rmax  (r-integration width,eq 5.3)",rmax ," [m]      = ",rmax*1d10," [Angstroem]"
!!>   write(6,'(a,i8             )')"nr    (number of q intpoints)     ",nr
!!>   write(6,'(a,e14.7          )')"rcf   (Gaussian Patterson width)  ",rcutfactor 
!!>   write(6,'(a                )')" ------------------------------------------------------------------------------ "
!!>   write(6,'(a                )')"Parameters controlling the amount of calculation: "
!!>   write(6,'(a,e14.7,a,e14.7,a)')"tmin  (min nozero tau )           ",taumin ," [s]      = ",taumin*1d9," [ns]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"tmax  (max        tau )           ",taumax ," [s]      = ",taumax*1d9," [ns]"
!!>   write(6,'(a,i8             )')"nt    (number of tau points )     ",nt
!!>   write(6,'(a,e14.7,a,e14.7,a)')"kz    (z-component of Q )         ",kz ,    " [1/m]    = ",kz*1d-10," [1/Angstroem]"
!!>   write(6,'(a,e14.7,a,e14.7,a)')"kip   (largest plane comp of Q )  ",kip ,   " [1/m]    = ",kip*1d-10," [1/Angstroem]"
!!>   write(6,'(a,i8             )')"nk    (number of in plane Qs)     ",nk
!!>   write(6,'(a,i8             )')"nz    (number of Qs z-comp. vals) ",nz
!!>   write(6,'(a                )')" ------------------------------------------------------------------------------ "
!!> 
!!> 
!!> 
!!> 
!!> 
!!>   if(taumin  > 0d0) then
!!>     tau(1:nt)  = taumin * (taumax/taumin)**([(i,i=0,nt-1)]/dble(nt-1))  ! here we need an extra entry with t=0 in element 0
!!>   else
!!>     tau(0:nt)  = [(i,i=0,nt)] * (taumax/nt)                             ! linear tau scaling
!!>   endif
!!>   tau(0)     = 0
!!> 
!!> 
!!> 
!!>   call romanov(n, B, K, d, gam, eta3, temp, rho, d_evanescent, nt, tau, nk, kzv, kipv, skt, nq, qmax, nr, rmax, rcutfactor )
!!>   
!!>   write(6,*)" DONE ... "  
!!>      
!!> !!! output of results !!!
!!> !!! if gnuplot is available a plot is created automatically
!!> !!! in addition the Sqt data are writen to skt.dgli and skt.dtr 
!!>      open(16,file='skt.gp')
!!>      open(17,file='skt.dgli')
!!>      write(6,*)" Skt: written to skt.dgli"
!!>      write(17,'(a,e16.7)')"B:=     ",B
!!>      write(17,'(a,e16.7)')"K:=     ",K
!!>      write(17,'(a,e16.7)')"gam:=   ",gam
!!>      write(17,'(a,e16.7)')"T:=     ",temp
!!>      write(17,'(a,e16.7)')"eta3:=  ",eta3
!!>      write(17,'(a,e16.7)')"deva:=  ",d_evanescent
!!>      write(17,'(a,e16.7)')"qmax:=  ",qmax
!!>      write(17,'(a,e16.7)')"rmax:=  ",rmax
!!>      write(17,'(a,e16.7)')"rcf:=   ",rcutfactor
!!>      write(17,'(a,i8)')   "n:=     ",n
!!>      write(17,'(a,i8)')   "nq:=    ",nq
!!>      write(17,'(a,i8)')   "nr:=    ",nr
!!>      write(17,'(a,20F12.6)')"kz:=",kzv*1d-10
!!>      write(17,'(a,20f12.6)')"kip:=",kipv*1d-10
!!>      do it=0,nt
!!>        write(6,'(21f14.3)') tau(it)*1d9,  (skt(ik,it)*1d20,ik=1,nk)
!!>        write(16,'(21f14.3)') tau(it)*1d9, (skt(ik,it)/skt(ik,0),ik=1,nk)        
!!>        write(17,'(21f14.3)') tau(it)*1d9, (skt(ik,it)*1d20,ik=1,nk)
!!>      enddo
!!>      close(16)
!!>      close(17)
!!> !!!
!!> !!! and as alternative output in (jcns datreat) format
!!> !!!
!!>      open(18,file='skt.dtr')
!!>        do ik=1,nk
!!>          write(18,'(a         )')"romanov output"
!!>          write(18,'(a,i0,a,i8 )')"r",ik," sqt/sq  vs  t/ns ",ik
!!>          write(18,'(a         )')"parameter"
!!>          write(18,'(a,e16.7   )')"B       ",B
!!>          write(18,'(a,e16.7   )')"K       ",K
!!>          write(18,'(a,e16.7   )')"gam     ",gam
!!>          write(18,'(a,e16.7   )')"temp    ",temp
!!>          write(18,'(a,e16.7   )')"eta3    ",eta3
!!>          write(18,'(a,e16.7   )')"deva    ",d_evanescent
!!>          write(18,'(a,e16.7   )')"qmax    ",qmax
!!>          write(18,'(a,e16.7   )')"rmax    ",rmax
!!>          write(18,'(a,e16.7   )')"rcf     ",rcutfactor
!!>          write(18,'(a,i8      )')"n       ",n
!!>          write(18,'(a,i8      )')"nq      ",nq
!!>          write(18,'(a,i8      )')"nr      ",nr
!!>          write(18,'(a,f12.6   )')"kz      ",kzv(ik) *1d-10
!!>          write(18,'(a,f12.6   )')"kip     ",kipv(ik)*1d-10
!!>          write(18,*)
!!>          write(18,'(a         )')"values "
!!>          do it=0,nt
!!>             write(18,'(f12.6,3x,f12.6)') tau(it)*1d9, skt(ik,it)/skt(ik,0)
!!>          enddo
!!>          write(6,*)
!!>          if(ik < nk) then
!!>             write(18,'(a      )')"#nxt"
!!>          else
!!>             write(18,'(a      )')"#eod"
!!>          endif
!!>        enddo
!!>      close(18)
!!> 
!!> 
!!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!> !! the following section just writes a gnuplot macro to display the curves
!!> !! written into skt.gp and the execute it using the system command
!!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!>      
!!>    open(11,file='skplot.gp')
!!>      write(11,'(a        )') 'set xlabel "{/=20 tau/ns}" '
!!>      write(11,'(a        )') 'set ylabel "{/=20 S(k,t)/S(k,0)}" '
!!>      write(11,'(a,e14.7,a )') 'set label 1 "{/=10 B = ', &
!!>                               B,                          ' N/m**2}" at graph 0.6,0.4  left '
!!>      write(11,'(a,f8.2,a )') 'set label 2 "{/=10 kappa = ', &
!!>                               kappa,                       ' kT  }" at graph 0.6,0.375 left '
!!>      write(11,'(a,f8.5,a )') 'set label 3 "{/=10 eta3 = ', &
!!>                               eta3,' kg/ms  }" at graph 0.6,0.35  left '
!!>      write(11,'(a,f8.5,a )') 'set label 4 "{/=10  gam   = ', &
!!>                              gam,                   ' N/m       }" at graph 0.6,0.325 left '
!!>      write(11,'(a,i8,a )') 'set label 5 "{/=10  n   = ', &
!!>                              n,                   '      }" at graph 0.6,0.30 left '
!!>      write(val,'(a,i0,a)')"(a,",nk,"f9.4, a)"
!!>      write(11,trim(val)) 'set label 6 "{/=10  kz   = ', &
!!>                              kzv*1d-10,                   '      }" at graph 0.6,0.2750 left '
!!>      write(11,trim(val)) 'set label 7 "{/=10  kip   = ', &
!!>                              kipv*1d-10,                   '     }" at graph 0.6,0.250 left '
!!>      write(11,'(a,f8.2,a )') 'set label 8 "{/=10  deva   = ', &
!!>                              d_evanescent*1d10, '                  }" at graph 0.6,0.225 left '
!!> 
!!>      write(11,'(a         )',advance='no') "plot "
!!>      do ik=1,nk
!!>        write(11,'(a,i0,a,i0,a )',advance='no') " 'skt.gp' using 1:",ik+1," with lines lw 3 lc ",ik+1,"  "
!!>        if(ik < nk) then
!!>           write(11,'(a)',advance='no')", "
!!>        else
!!>           write(11,'(a)')" "
!!>        endif
!!>      enddo
!!>      write(11,'(a        )') "set terminal postscript eps enhanced color font 'Helvetica,10'" 
!!>      write(11,'(a        )') 'set output "lastgplot.eps" ' 
!!>      write(11,'(a         )',advance='no') "plot "
!!>      do ik=1,nk
!!>        write(11,'(a,i0,a,i0,a )',advance='no') " 'skt.gp' using 1:",ik+1," with lines lw 3 lc ",ik+1,"  "
!!>        if(ik < nk) then
!!>           write(11,'(a)',advance='no')", "
!!>        else
!!>           write(11,'(a)')" "
!!>        endif
!!>      enddo
!!>  
!!>    close(11)
!!> 
!!>    call system('gnuplot skplot.gp -p')
!!> 
!!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!> !! end of gnuplot macro generation
!!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!> 
!!> 
!!> 
!!> end program romanov_v0
!!> 






subroutine romanov &
  (n, B, K, d, gam, eta3, Temperatur, rho, d_evanescent, nt, tau, nk, kz, kip, skt, nq, qmax, nr, rmax, rcutfactor )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  LIT:
!!  V. P. Romanov and S. V. Ul’yanov PHYSICAL REVIEW E 66, 061701 􏰀2002􏰁
!!
!!  Implementation of the solution for dynamics of a solid supported layer stack
!!  up to the computation of the scattering function analog Eq 52.
!!  modifed such that the dynamics is included and the necessary modification
!!  in Eq. 5.3 to take into account that this is analog to a Patterson fucntion and
!!  for a finite support needs a self convolution of the support as
!!  modulation to give reasonable results (i.e. nonnegative scattering intensities)
!!  this is achieved by assuming a Gaussian modulation.
!!  The prefactor of Eq. 5.2 with the layer formfactor and contrast is not included here!
!!  The amplitudes at the layers are modulated by the evanescent wave amplitudes at that
!!  distance.
!!  Observe: the fixed surface it at index n+1 !
!!
!! Author: Michael Monkenbusch, JCNS-1, FZ-Juelich
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none

  integer, intent(in)            :: n              ! number of layers
  double precision, intent(in)   :: B              ! layer comperssion modulus in SI-units [N/m**2] (typ 1e6..1e7)
  double precision, intent(in)   :: K              ! bending modulus in SI-units [N] (typ Boltzmannkonstante *  Temperatur / d * 10d0)
  double precision, intent(in)   :: d              ! layer distance in SI-units [m] (typ 60e-10)
  double precision, intent(in)   :: gam            ! surface tension of free surface in SI-units [N/m] (typ 0.030)
  double precision, intent(in)   :: eta3           ! layer viscosity in SI-units [] (typ 0.001..0.01) ?
  double precision, intent(in)   :: Temperatur     ! temprature in SI-units [K] (typ 300)
  double precision, intent(in)   :: rho            ! density in SI-unis [kg/m**3] (typ 1000)
  double precision, intent(in)   :: d_evanescent   ! penetration depth of evanescent wave in SI-units [m] (typ 500e-10)
  integer         , intent(in)   :: nt             ! number of tau values to be computed
  double precision, intent(in)   :: tau(0:nt)      ! list of tau values in SI-Units [s]  /typ (0..20e-9), tau(0)=0
  integer         , intent(in)   :: nk             ! number of k-values
  double precision, intent(in)   :: kz(nk)         ! scattering vectors z-component in SI-units [1/m] (typ 0.15e10)
  double precision, intent(in)   :: kip(nk)        ! in plane component of scattering vector in SI-units [1/m] (typ 0.02e10)
  double precision, intent(out)  :: skt(nk,0:nt)   ! s(k,t) OUTPUT
  integer         , intent(in)   :: nq             ! number of q integration points (200+)
  double precision, intent(in)   :: qmax           ! q integration rangei in SI  (1.0d10 m**-1)
  integer         , intent(in)   :: nr             ! number of r integration points (200)
  double precision, intent(in)   :: rmax           ! maximum for r-integration (this is the size of the patterson function carrier) (2000d-10)
  double precision, intent(in)   :: rcutfactor     ! rmax rcutfactor = Gaussian width (plain) of Gaussian Patteron fucntion for Gnm (0.3)

  
  double precision, parameter :: Pi = 4*atan(1d0)
  double precision, parameter :: Boltzmannkonstante  =     1.3806488d-23


  double precision   :: ul(n,n), uul(n,n), faa(n)
  double precision   :: U(-1:n), x(n), alp

  
  integer            :: i, j, l, M, ik

  double precision   :: famax = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer            :: ir, it, iq, nn, mm
  double precision   :: ulq(n,n), faa1(n), faa2(n), fa1, fa2, ft1, ft2, aal, bbx
  complex(kind= 8)   :: omega(2,n)
  complex(kind= 8)   :: omega1(2,n)
  complex(kind= 8)   :: omega2(2,n)
  complex(kind= 8)   :: a_omega(2,n) 
  complex(kind= 8)   :: ci=(0d0,1d0), c1=(1d0,0d0), cs 
  complex(kind= 8)   :: om1, om2
  double precision   :: beta1, beta2
  double precision   :: uu_nmrt(n,n,0:nr,0:nt)
  double precision   :: gnmt(n,n,0:nt)
  double precision   :: rip(0:nr)
  double precision   :: qip(0:nq)
  


  double precision, external :: bessel_trapezint
  double precision, external :: bessel_box0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! the <u_n(r,t) u_m(0)> sequence store to uu_nmrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!
!! prepare vectors r in plane and tau
!!

   rip(0:nr)  = [(i+0.001d0,i=0,nr)] * (rmax/nr) 
   qip(0:nq)  = [(i,i=0,nq)] * (qmax/nq)

   uu_nmrt    = 0
         write(6,'(a,e14.7,a,i4,a)')"Preparing <u_n(r,t) u_m(0,0)> with q = 0..",qmax," m**-1 in ",nq," steps"
         write(6,'(a)',advance='no')"Progress(X = 10 q-values): "
dq:   do iq=0,nq-1
!         write(6,*)"iq = ",iq,"    qip= ",qip(iq)*1d-10
         if(mod(iq,10)==0) write(6,'(a1)',advance='no') "X"

         call romanov_layers_omsolv(n, qip(iq), B, K, gam, rho, eta3, d, omega, ul)
         call romanov_mode_energies(n, qip(iq), B, K, gam, d, ul, faa )
         ulq     = ul/2
         faa1    = faa *(2*Pi) / ( Boltzmannkonstante * Temperatur )
         omega1  = omega
         call romanov_layers_omsolv(n, qip(iq+1), B, K, gam, rho, eta3, d, omega, ul)
         call romanov_mode_energies(n, qip(iq+1), B, K, gam, d, ul, faa )
         ulq     = ulq+ul/2
         faa2    = faa *(2*Pi) / ( Boltzmannkonstante * Temperatur )
         omega2  = omega



dt:      do it=0,nt
!            write(6,*)"it = ",it,"  tau = ",tau(it)
dr:         do ir=0,nr
!               write(6,*)"ir = ",ir," r= ",rip(ir)
               bbx = bessel_box0(rip(ir),qip(iq),qip(iq+1))
dl:            do l=1,n

!! using expresions for the damped harmonic oscillator form 
!! Uhlenbeck and Ornstein, Rhys. Rev. 36 (1930), p 823-841, Eq. 33                  
!! for the time dependece of the displacement amplitude correlator 
                  beta1    =  - ( imag( omega1(1,l) + omega1(2,l) )  )   
                  om1      =    ( omega1(1,l) - omega1(2,l) ) / 2d0

                  beta2    =  - ( imag( omega2(1,l) + omega2(2,l) )  )     
                  om2      =    ( omega2(1,l) - omega2(2,l) ) / 2d0
!> ---- the following commented expressions have to be replaced by the 
!>      mathematical identical expressions that follow after the comment block !>
!>      since otherwise overflows in the exponentials (resp. cos--> cosh) will interfere 
!>      with the proper result (note om1,2 is complex, functions generic, ci=I)
!>
!>                  ft1  = exp(-(beta1/2d0)*tau(it)) * ( beta1/(2d0*(om1)) * sin((om1)*tau(it)) + cos((om1)*tau(it)) )
!>                  ft2  = exp(-(beta2/2d0)*tau(it)) * ( beta2/(2d0*(om2)) * sin((om2)*tau(it)) + cos((om2)*tau(it)) )
!>
             
                  ft1 = (-ci)*beta1/(2d0*om1) *    &
                              0.5d0*(exp((-(beta1/2)+om1*ci)*tau(it))-exp((-(beta1/2)-om1*ci)*tau(it))) &
                         + 0.5d0*(exp((-(beta1/2)+om1*ci)*tau(it))+exp((-(beta1/2)-om1*ci)*tau(it))) 
 
                  ft2 = (-ci)*beta2/(2d0*om2) *     &
                              0.5d0*(exp((-(beta2/2)+om2*ci)*tau(it))-exp((-(beta2/2)-om2*ci)*tau(it))) &
                         + 0.5d0*(exp((-(beta2/2)+om2*ci)*tau(it))+exp((-(beta2/2)-om2*ci)*tau(it))) 
 
                  
                  fa1  = ft1/faa1(l)
                  fa2  = ft2/faa2(l)

                  aal  = bbx * (fa1+fa2)/2       
              
                  do nn=1,n
                    do mm=1,n
                      uu_nmrt(nn,mm,ir,it) =  uu_nmrt(nn,mm,ir,it) + aal * ulq(nn,l)*ulq(mm,l)
                    enddo
                  enddo
             enddo dl 

             enddo dr
         enddo dt
      enddo dq
   

   write(6,*)" Ready "


!!! compute scattering function
!!!
!!! start with G(n,m,t)(kip,kz)  !!!!! ! here we need definitively a Gaussian cutoff to mimick the convolution involved in the
!!!                              !!!!! ! Patterson function involved here, a sharp cutoff causes spurious oscillations and negative
!!!                              !!!!! ! intensities since it ignores the drop of combined density for any finite summation range
!!!                              !!!!! ! Gaussian works fine since Gaussian convoluted with Gaussian = Gaussian
!!!                              !!!!! ! one might also use the real function disk (*) disk, which yields positive intensity but
!!!                              !!!!! ! still some sidelobes of the specular, which might be disturbing....
!!!                              !!!!! ! Gaussian multipliction !!
!!!                              !!!!! ! the plain disk Patterson function on carrier 0..2R is
!!! $ p(r) = {R}^{2}\pi -1/2\,R\sqrt {4-{\frac {{r}^{2}}{{R}^{2}}}}r-2\,{R}^{2}\arcsin \left( 1/2\,{\frac {r}{R}} \right) $


!!!! the scattering function !!!!
    skt = 0

sk: do ik=1,nk

      gnmt = 0
     gn:   do nn=1,n
     gm:     do mm=1,n
     gt:       do it=0,nt
     gr:         do ir=0,nr-1    
                      
                   gnmt(nn,mm,it) =  gnmt(nn,mm,it) & 
                                  +  bessel_box0(kip(ik),rip(ir),rip(ir+1))        &
                                  *  exp(-((rip(ir)+rip(ir+1))/(2*rcutfactor*rmax))**2) &  !! this is the Gaussian to account for Patterson function
                                  *  (   exp( (kz(ik)**2)*uu_nmrt(nn,mm,ir,it))    &
                                       + exp( (kz(ik)**2)*uu_nmrt(nn,mm,ir+1,it))  &
                                      ) / 2
         
              enddo gr
            enddo gt
          enddo gm
        enddo gn
        
st:  do it=0,nt
sn:    do nn=1,n
sm:       do mm=1,n
            skt(ik,it) = skt(ik,it) + cos( kz(ik) * d * (nn-mm) )             &
                              * exp( -d *(n+1-nn)/d_evanescent )    &
                              * exp( -d *(n+1-mm)/d_evanescent )    &
                              * exp( -(kz(ik)**2)*(uu_nmrt(nn,nn,0,0)+uu_nmrt(mm,mm,0,0))/2 ) &
                              * gnmt(nn,mm,it)
          enddo sm
       enddo sn
      enddo st
     enddo sk
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!! ---> procedures that are used by the romanov subroutine ---------------------------------------- !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     


contains


subroutine romanov_mode_energies(n, q_plane, B, K, gam, d, ul, faa )
  ! determine the mode contribution coefficients to the free energy
  !
  ! use SI UNITS
  !
  ! see eq 3.1 we give the mode dependent analog to M in faa
  !
  implicit none
  integer,         intent(in)  :: n          ! number of free layers (layer n+1 is the fixed boundary)
                                             ! layer 1 is the free surface
  double precision,intent(in)  :: q_plane    ! q-mode-vector in the layer plane  (SI!)
  double precision,intent(in)  :: B          ! compression modulus [N/m**2]
  double precision,intent(in)  :: K          ! bending modulus [N]
  double precision,intent(in)  :: gam        ! surface tension [N/m|
  double precision,intent(in)  :: d          ! layer spacing   [m]
  double precision,intent(in)  :: ul(n,n)    ! eigenvectors of the modes
  double precision,intent(out) :: faa(n)     ! the energy contribution coefficient of <a(l)a(l)>


  double precision             :: alpha, y
  double precision             :: mu(n)
  integer                      :: i, l

  alpha   =     ( d * gam * q_plane**2 ) / B
  y       =     -1d0 - d**2 * K * q_plane**4 / (2*B)

  do l=1,n

     mu(1)    =   (2*y+1-alpha) * ul(1,l) + ul(2,l)  
     if(n > 2) then
        do i=2,n-1
           mu(i) = ul(i-1,l) + 2*y * ul(i,l) + ul(i+1,l)
        enddo
     endif
     mu(n)    =  ul(n-1,l) + 2*y * ul(n,l) 

     faa(l)   = dot_product(ul(:,l),mu) * (-B/d)

  enddo
    
end subroutine romanov_mode_energies
   

subroutine romanov_layers_omsolv(n, q_plane, B, K, gam, rho, eta3, d, omega, ul)
  ! solve the eigenvalue proble 2.6 from ROMANOV & UL'YANOV PYSICAL REVIEW E 66, 061701 (2002)
  ! ...doing this numerically (here) allows further (future) extensions to layer-layer variations
  !    of the parameters (not yet implemented).
  !
  ! use SI UNITS 
  !
  implicit none
  integer,         intent(in)  :: n          ! number of free layers (layer n+1 is the fixed boundary)
                                             ! layer 1 is the free surface
  double precision,intent(in)  :: q_plane    ! q-mode-vector in the layer plane  (SI!)
  double precision,intent(in)  :: B          ! compression modulus [N/m**2]
  double precision,intent(in)  :: K          ! bending modulus [N]
  double precision,intent(in)  :: gam        ! surface tension [N/m|
  double precision,intent(in)  :: rho        ! density (fluid) [kg/m**3]
  double precision,intent(in)  :: eta3       ! layer viscosity [N..]
  double precision,intent(in)  :: d          ! layer spacing   [m]
  complex(kind= 8),intent(out) :: omega(2,n) ! frequencies of mode [1/s]
  double precision,intent(out) :: ul(n,n)    ! eigenvectors of the modes


  double precision             :: a_diag(n)
  double precision             :: a_off_diag(n)
  double precision             :: lambda(n)     ! eigenvalues --> omega

  double precision             :: x20, alpha
  complex(kind= 8)             :: cs
  complex(kind= 8),parameter   :: ci = (0.0d0, 1.0d0), c1 = (1d0,0d0)

  integer                      :: i, ierr, info


  ! here we formulate the proble such that the Eq. (2.6) becomes a tridiagonal
  ! Eigenvalue problem 

  x20     =   - ( K * d**2 * q_plane**4 + 2*B ) / B
  alpha   =     ( d * gam * q_plane**2 ) / B

  a_diag(2:n)     =  x20
  a_diag(1)       =  x20 + (1-alpha)
  a_off_diag(1:n) =  1

  call tridiag_diagonalize (n, n, a_diag, a_off_diag, lambda, ul, ierr)
  if(ierr .ne. 0) then
    write(6,*)"tridiag: ierr = ",ierr
    stop
  endif

  ! now we solve the lambda value for the two omegas
  
  do i=1,n
    cs         = sqrt( c1 * ( -d**2*eta3**2*q_plane**4 - 4*B*rho*lambda(i) ) )
    omega(1,i) = -1d0/(2*rho*d) * ( ci * eta3 * q_plane**2 * d + cs)
    omega(2,i) = -1d0/(2*rho*d) * ( ci * eta3 * q_plane**2 * d - cs)

!  write(6,'(i8,f12.6,2x,2e16.7,2x,2e16.7,2x,2e16.7)') i, lambda(i), cs , omega(1,i), omega(2,i)
  enddo
  

  
end subroutine romanov_layers_omsolv








SUBROUTINE tridiag_diagonalize (NM, N, Din, Ein, eigenvalues, Z, IERR)

! adaption from SLATEC (translated to double precision and modern fortran)
!
!***BEGIN PROLOGUE  IMTQL2
!***PURPOSE  Compute the eigenvalues and eigenvectors of a symmetric
!            tridiagonal matrix using the implicit QL method.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (IMTQL2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure IMTQL2,
!     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson,
!     as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a SYMMETRIC TRIDIAGONAL matrix by the implicit QL method.
!     The eigenvectors of a FULL SYMMETRIC matrix can also
!     be found if  TRED2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        Din contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional double precision  ::  array, dimensioned D(N).
!
!        Ein contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional double precision  ::  array, dimensioned
!          E(N).
!
!        Z contains the transformation matrix produced in the reduction
!          by  TRED2,  if performed.  This transformation matrix is
!          necessary if you want to obtain the eigenvectors of the full
!          symmetric matrix.  If the eigenvectors of the symmetric
!          tridiagonal matrix are desired, Z must no longer  contain the identity
!          matrix  Z is a two-dimensional double precision  ::  array, dimensioned
!          Z(NM,N).
!
!      On OUTPUT
!
!        D is kept
!
!        Ein is kept
!
!        Eigenvalues  contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are unknown
!
!
!        Z contains orthonormal eigenvectors of the full symmetric
!          or symmetric tridiagonal matrix, depending on what it
!          contained on input.  If an error exit is made,  Z contains
!          the eigenvectors associated with the stored eigenvalues.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues and eigenvectors should be correct
!                     for indices 1, 2, ..., IERR-1, but the eigenvalues
!                     are not ordered.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IMTQL2
!

      implicit none
   
      integer         , intent(in)   :: NM              ! leading dimension of Z
      integer         , intent(in)   :: N               ! dimension 
      double precision, intent(in)   :: Din(N)
      double precision, intent(in)   :: Ein(N)
      double precision, intent(out)  :: eigenvalues(N)
      double precision, intent(out)  :: Z(NM,N)         ! eigenvectors
      integer         , intent(out)  :: ierr            ! error code
      

      integer, parameter             :: maxit = 30

      INTEGER           ::  I,J,K,L,M,II,MML
      double precision  ::  D(N),E(N)
      double precision  ::  B,C,F,G,P,R,S,S1,S2
!      double precision  ::  PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  IMTQL2
      IERR = 0
      IF (N .EQ. 1) return
  
      D(1:N)   = Din(1:N)  

      E(1:N-1) = Ein(2:N)
      E(N)     = 0.0d0

      Z(1:N,1:N) = 0d0
      do i=1,N
        Z(i,i) = 1d0
      enddo

d0:   do L = 1, N
         J = 0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
diter:  do 
d1:      do m = l, n
            if(m == n) exit d1
            S1 = ABS(D(M)) + ABS(D(M+1))
            S2 = S1 + ABS(E(M))
            if( s2 == s1 ) exit d1
         enddo  d1

         P = D(L)
         IF (M == L)    cycle d0
         IF (J >= maxit)   then
            ierr = L
            return
         endif

         J = J + 1

!     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0d0 * E(L))
!         R = PYTHAG(G,1.0d0)
         r = sqrt(g*g+1)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = 1.0d0
         C = 1.0d0
         P = 0.0d0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
d2:      DO  II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (ABS(F) >= ABS(G)) then
               C = G / F
               R = SQRT(C*C+1.0d0)
               E(I+1) = F * R
               S = 1.0d0 / R
               C = C * S
            else
               S = F / G
               R = SQRT(S*S+1.0d0)
               E(I+1) = G * R
               C = 1.0d0 / R
               S = S * C
            endif
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0d0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
!     .......... FORM VECTOR ..........
 d4:        DO  K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
            enddo d4
!
         enddo d2
!
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0d0
         enddo diter
  enddo d0

 

!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
ds1:   DO II = 2, N
         I = II - 1
         K = I
         P = D(I)

         DO  J = II, N
            IF (D(J) .GE. P) cycle
            K = J
            P = D(J)
         enddo

         IF (K .EQ. I) cycle ds1
         D(K) = D(I)
         D(I) = P

         DO  J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
         enddo
      enddo ds1


      eigenvalues = D
!
      END subroutine tridiag_diagonalize



   subroutine ChebyshevU( x, n, u )
     implicit none
     double precision, intent(in) :: x
     integer         , intent(in) :: n
     double precision, intent(out):: u(-1:n)

     integer  :: m

     u(-1) = 0
     u(0)  = 1
     if( n == 0 ) return
     u(1) = 2 * x
     if( n == 1 ) return
   
     do m = 2,n
       u(m) = 2 * x * u(m-1) - u(m-2)
     enddo

   end subroutine ChebyshevU


  end subroutine romanov

function bessel_box0( r, a, b ) result (y)
    implicit none
    double precision, intent(in)  :: r
    double precision, intent(in)  :: a
    double precision, intent(in)  :: b

    double precision              :: y

     y = -(a*Bessel_J1(r*a)-b*Bessel_J1(r*b))/r

end function bessel_box0



