 FUNCTION th_ring2dm4(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen
!  S. Goossen et al., PRL 2014, 113, 168302
      use theory_description 
      implicit none 
      real    :: th_ring2dm4
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
     double precision :: diff       ! limiting diffusion (NMR value) in [A**2/ns]                                     
     double precision :: r021       ! reference mean squared displacement at first transition point to medium timed di
     double precision :: alpha      ! sub-diffusion  exponent of SHORT time diffusion                                 
     double precision :: r022       ! reference mean squared displacement at second transition point to D0 in [A**2]  
     double precision :: beta       ! sub-diffusion  exponent of medium time diffusion                                
     double precision :: a_cross    ! transition exponent between short and long time diffusion (sharper kink for larg
     integer          :: n          ! number of segments in one ring                                                  
     double precision :: l          ! effective segment length                                                        
     double precision :: nu0        ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
     double precision :: nue        ! chain statistics exponent (larger dist)                     
     double precision :: wl4        ! Rouse rate in [A**4/ns]                                                         
     double precision :: ntrans       ! transition mode number between simple ring-Rouse and large p modification       
     double precision :: nwidth     ! sharpness of transition                                                         
     double precision :: pmin       !                                                                                 
     double precision :: pwidth     ! 

     double precision :: pneglim       !                                                                                 
     double precision :: pposlim     ! 
     double precision :: nuan       ! 
     integer          :: mode       ! 

! the recin parameter representation 
     double precision :: q          ! q                                                                               
! the reout parameter representation 
     double precision :: Rg         ! predicted ring radius of gyration                                               
 
     double precision :: th
 
     double precision :: t
     double precision :: rr, nu, Sq, Sqt, Sqt0
     double precision :: msd, msd_im
     double precision :: tauinf_ren

     logical          :: parameter_is_q, parameter_is_t


!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ring2dm4'
       nparx =       20
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ring2dm4 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen with nue(nn) and Ap(ip)" //cr//parspace//&
       " compute S(q,t)/S(q) as function of t, if mode = 0 and record contains q as parameter and no t"  //cr//parspace//&
       "         S(q,t)      as function of q, if mode = 0 and record contains t as parameter and no q" //cr//parspace//&
       " compute MSD         as function of t, if mode = 1 and record contains q as parameter and no t"  //cr//parspace//&
       "         S(q,t)      as function of q, if mode = 2 and record contains t as parameter and no q" //cr//parspace//&
       " compute MSD_int     as function of t, if mode = 2 and record contains q as parameter and no t"  //cr//parspace//&
       "         S_int(q,t)  as function of q, if mode = 2 and record contains t as parameter and no q" //cr//parspace//&
       "         S_inc(q,t)  as function of t, if mode = 3 and record contains q as parameter and no t"       

       th_citation(idesc)     = " S. Goossen et al., PRL 2014, 113, 168302"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'diff    '  ! limiting diffusion (NMR value) in [A**2/ns]                                     
        parnam ( 3) = 'r021    '  ! reference mean squared displacement at first transition point to medium timed di
        parnam ( 4) = 'alpha   '  ! sub-diffusion  exponent of SHORT time diffusion                                 
        parnam ( 5) = 'r022    '  ! reference mean squared displacement at second transition point to D0 in [A**2]  
        parnam ( 6) = 'beta    '  ! sub-diffusion  exponent of medium time diffusion                                
        parnam ( 7) = 'a_cross '  ! transition exponent between short and long time diffusion (sharper kink for larg
        parnam ( 8) = 'n       '  ! number of segments in one ring                                                  
        parnam ( 9) = 'l       '  ! effective segment length                                                        
        parnam (10) = 'nu0     '  ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
        parnam (11) = 'nue     '  ! chain statistics exponent (nu=0.5 => random walk, Gaussian) for large dist
        parnam (12) = 'wl4     '  ! Rouse rate in [A**4/ns]                                                         
        parnam (13) = 'ntrans  '  ! transition mode number between simple ring-Rouse and large p modification       
        parnam (14) = 'nwidth  '  ! sharpness of transition                                                         
        parnam (15) = 'pmin    '  ! Rouse mode amplitude transition pmin value                                      
        parnam (16) = 'pwidth  '  ! Rouse mode amplitude transition width                                           
        parnam (17) = 'pneglim '  ! Rouse mode amplitude transition value at p = -infty                             
        parnam (18) = 'pposlim '  ! Rouse mode amplitude transition value at p = +infty                             
        parnam (19) = 'nuan    '  ! Nue Modification N**(2*nuan) in A                         
        parnam (20) = 'mode    '  ! what to compute                      
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "     prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "ABS: limiting diffusion (NMR value) in [A**2/ns]" !//cr//parspace//&
        th_param_desc( 3,idesc) = "ABS: reference mean squared displacement at first transition point to medium timed diff." !//cr//parspace//&
        th_param_desc( 4,idesc) = "ABS: sub-diffusion  exponent of SHORT time diffusion" !//cr//parspace//&
        th_param_desc( 5,idesc) = "ABS: reference mean squared displacement at second transition point to D0 in [A**2]" !//cr//parspace//&
        th_param_desc( 6,idesc) = "ABS: sub-diffusion  exponent of medium time diffusion" !//cr//parspace//&
        th_param_desc( 7,idesc) = "ABS: transition exponent between short and long time diffusion (sharper kink for larger a)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "INT: number of segments in one ring" !//cr//parspace//&
        th_param_desc( 9,idesc) = "ABS: effective segment length" !//cr//parspace//&
        th_param_desc(10,idesc) = "ABS: chain statistics exponent (nu=0.5 => random walk, Gaussian)" !//cr//parspace//&
        th_param_desc(11,idesc) = "ABS: chain statistics exponent (nu=0.5 => random walk, larger dist)" !//cr//parspace//&
        th_param_desc(12,idesc) = "ABS: Rouse rate in [A**4/ns]" !//cr//parspace//&
        th_param_desc(13,idesc) = "ABS: transition mode number between simple ring-Rouse and large p modification" !//cr//parspace//&
        th_param_desc(14,idesc) = "ABS: sharpness of transition" !//cr//parspace//&
        th_param_desc(15,idesc) = "     Fermi function pmin for Rouse mode amplitude       " !//cr//parspace//&
        th_param_desc(16,idesc) = "     Fermi function pwidth for Rouse mode amplitude " !//cr//parspace//&
        th_param_desc(17,idesc) = "     Fermi function amplitude for Rouse mode at p=-infty       " !//cr//parspace//&
        th_param_desc(18,idesc) = "     Fermi function amplitude for Rouse mode at p=+infty " !//cr//parspace//&
        th_param_desc(19,idesc) = "     nue exponent for N in A ~ N**(2*nuan) " !//cr//parspace//&
        th_param_desc(20,idesc) = "     0: Sqt, 1: msd+diff 2:msd-int 3:sqt_inc (*ampli) " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "Rg       > predicted ring radius of gyration"
! 
        th_ring2dm4 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =       (  pa( 1) )
      diff     =    abs(  pa( 2) )
      r021     =    abs(  pa( 3) )
      alpha    =    abs(  pa( 4) )
      r022     =    abs(  pa( 5) )
      beta     =    abs(  pa( 6) )
      a_cross  =    abs(  pa( 7) )
      n        =   nint(  pa( 8) )
      l        =    abs(  pa( 9) )
      nu0      =    abs(  pa(10) )
      nue      =    abs(  pa(11) )
      wl4      =    abs(  pa(12) )
      ntrans   =    abs(  pa(13) )
      nwidth   =    abs(  pa(14) )
      pmin     =       (  pa(15) )
      pwidth   =       (  pa(16) )
      pneglim  =       (  pa(17) )
      pposlim  =       (  pa(18) )
      nuan     =       (  pa(19) )
      mode     =   nint(  pa(20) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q
      xh =   0.1
      call parget('q       ',xh,iadda,ier)
      parameter_is_q = (ier == 0)
      q        = xh
! >>> extract: t
      xh =   0.1
      call parget('t       ',xh,iadda,ier)
      parameter_is_t = (ier == 0)
      t        = xh

 
    if(parameter_is_q .and. parameter_is_t) then
      write(*,*)"th_ring2am4: record parameters of ",iadda," may only contain q or t "
      th_ring2dm4 = 0.0
      return
    endif
    if(.not.parameter_is_q .and. .not.parameter_is_t) then
      write(*,*)"th_ring2am4: record parameters of ",iadda," must contain either q or t "
      th_ring2dm4 = 0.0
      return
    endif
!
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 

     if(parameter_is_q) then
       t   = x
     else if(parameter_is_t) then
       q   = x
     endif


     call Nrouse_ring2dm4(q,t,N, l, nu0, nue, wl4, ntrans, nwidth, pmin, pwidth, pneglim, pposlim, nuan, &
                                diff, r021, alpha, r022, beta, a_cross, Rg, Sq, Sqt, msd, msd_im, Sqt0)

     if(mode == 1) then
       th =  msd
       if(parameter_is_t) th = ampli * Sqt
     else if (mode == 2) then
       th = msd_im
       if(parameter_is_t) th = ampli * Sqt0
     else if (mode == 3 ) then
       th = ampli*exp(-msd * q*q / 6d0) 
     else
       th =  ampli * Sqt/Sq
       if(parameter_is_t) th = ampli * Sqt
     endif

     th_ring2dm4 = th

 
! ---- writing computed parameters to the record >>>  
      call parset('Rg      ',sngl(Rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
       subroutine Nrouse_ring2dm4(q,t,N, l, nu0, nue, wl4, ntrans, nwidth, pmin, pwidth, pneglim, pposlim, nuan, &
                                diff, r021, alpha, r022, beta, a_cross, Rg, Sq, Sqt, msd, msd_im, Sqt0)
!      ================================================================================================
!
! Rouse expression for a ring polymer
! with Nb segments each, the ideal Ree of the linear version of the polymer is given by R.
! amod contains mode modifiers
! Input parameters: (polymer and model descriptor are contained shared with main )
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    ... 
!    ...   ----> see cmments below
! Output parameters:
!    Rg    <--- radius of gyration
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)



       double precision, intent(in) ::  q      ! q-value
       double precision, intent(in) ::  t      ! time
       integer         , intent(in) ::  N      ! number of segments
       double precision, intent(in) ::  l      ! segment length
       double precision, intent(in) ::  nu0    !  below ntrans
       double precision, intent(in) ::  nue    !  exponent
       double precision, intent(in) ::  wl4    ! Rouse rate
       double precision, intent(in) :: ntrans    ! transition mode number between 
       double precision, intent(in) :: nwidth  ! sharpness of transition  

       double precision, intent(in) :: pmin    ! transition mode number between 
       double precision, intent(in) :: pwidth  ! sharpness of transition  

       double precision, intent(in) :: pneglim ! transition mode number between 
       double precision, intent(in) :: pposlim ! sharpness of transition  
       double precision, intent(in) :: nuan     ! in Ap  N==> N**(2*nuan)  
     
       double precision, intent(in) :: diff     ! limiting diffusion (NMR value) in [A**2/ns]
       double precision, intent(in) :: r021     ! reference mean squared displacement at first transition point to medium timed di
       double precision, intent(in) :: alpha    ! sub-diffusion  exponent of SHORT time diffusion                                 
       double precision, intent(in) :: r022     ! reference mean squared displacement at second transition point to D0 in [A**2]  
       double precision, intent(in) :: beta     ! sub-diffusion  exponent of medium time diffusion                                
       double precision, intent(in) :: a_cross  ! transition exponent between short and long time diffusion (sharper kink for             
       double precision, intent(out):: Rg      ! radius of gyration                                              
       double precision, intent(out):: Sq      ! S(q)                                              
       double precision, intent(out):: Sqt     ! S(q,t)                                              

       double precision, intent(out):: msd     ! msd incl diffusion and internal modes                                           
       double precision, intent(out):: msd_im  ! msd internal modes only                                             
       double precision, intent(out):: Sqt0    ! Sqt without centre-of-mass diffusion                                             
       integer                          ::  nn,mm,p

       double precision :: ff2, nu(0:N), rr
       double precision :: rate, tau_R, W
       double precision :: nuofp

       double precision :: cosarray(0:N,N/2), ewfac(N/2)
       double precision :: traf(0:N), tram(0:N), Ap(0:N)



       integer          :: ip, N2


       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif


!       nu    = 2*nue
! ---- and the Rousetime ----
       W     = Wl4 / l**4
       tau_R = N**2/(W * Pi**2)

! ---- init sums ----
       N2 = N/2
! p(even) = 2*ip in the following...
!$OMP PARALLEL DO
  do nn=0,N
    tram(nn) = 1d0/(1.0d0+exp((nn-ntrans)/nwidth))
!                                 -----   ------
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
  do ip=1,N2
    traf(ip) = tram(nint(dble(N)/dble(2*ip)))
  enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do nn=0,N
    nu(nn) = 2*(nu0*tram(nn) + (1-tram(nn))*nue)
  enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(nuofp)
  do ip=1,N2
    nuofp   =   nue * 1d0/(1.0d0+exp((2*ip-pmin)/pwidth))       &  !! new
             +  nu0 * (1d0-1d0/(1.0d0+exp((2*ip-pmin)/pwidth)))    !! new
    nuofp   = nuofp * nuan
    Ap(ip) = ( pneglim * 1d0/(1.0d0+exp((2*ip-pmin)/pwidth))   &
             + pposlim * (1d0-1d0/(1.0d0+exp((2*ip-pmin)/pwidth)))   &
             ) / (2*ip)**2  &
           * (-2d0* dble(N)**(2*nuofp) *(l*q)**2/(3d0*pi**2))  !! HIER KEIN N**2nu OK ??
  enddo
!$OMP END PARALLEL DO



!$OMP PARALLEL DO 
       do nn=0,N 
        do ip=1,N2
         
         cosarray(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))               &
                         * Ap(ip)
        enddo
       enddo
!$OMP END PARALLEL DO

       msd_im = 0
!$OMP PARALLEL DO PRIVATE(rate) REDUCTION(+:msd_im)
       do ip=1,N/2  
         rate      = (2*ip)**2 / tau_R
         ewfac(ip) = 1.0d0-exp(-t * rate )
         msd_im = msd_im -Ap(ip) * ewfac(ip)
       enddo
!$OMP END PARALLEL DO
       msd_im = 6 * msd_im / q**2
 

           Sq  = 0
           Sqt = 0
           Rg  = 0

!??           ff2  = -2d0* dble(N)**nu *(l*q)**2/(3d0*pi**2)
           ff2 = 1d0    ! wird schon in Zeile 262 anmultipliziert mit nu(n)  IST DAS KONSISTENT??

! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt,Rg)
           do nn = 1,N
            do mm = 1,N


              Sq  = Sq  + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) )
              Sqt = Sqt + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm))  &
                              + ff2*sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)) )

              Rg  = Rg  + l**2  * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm))

            enddo
           enddo
!$OMP END PARALLEL DO

           Sq  = Sq /N
           Sqt = Sqt/N

! and the com diffusion
     
           rr=((exp(((-alpha+beta)*log(r021/r022)+alpha*beta*(log(2d0)+log(3d0)-      &
           log(r022/diff)))/beta)*r022*t**alpha)**a_cross+(exp((log(2d0)+log(3d0)- &
           log(r022/diff))*beta)*r022*t**beta)**a_cross+(6d0*diff*t)**a_cross)**(1d0/a_cross)

           Sqt0 = sqt
           Sqt = sqt * exp(-rr * q*q / 6d0) 

           Rg  = sqrt(Rg/2d0/N**2)

           msd = rr + msd_im

       end subroutine Nrouse_ring2dm4

 
 end function th_ring2dm4
