 FUNCTION th_ring2am4(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen
!  S. Goossen et al., PRL 2014, 113, 168302
      use theory_description 
      implicit none 
      real    :: th_ring2am4
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
     double precision :: nue        ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
     double precision :: wl4        ! Rouse rate in [A**4/ns]                                                         
     double precision :: pmin       ! transition mode number between simple ring-Rouse and large p modification       
     double precision :: pwidth     ! sharpness of transition                                                         
     double precision :: f0         ! prefactor f(p) limit for small p values (default 1)                             
     double precision :: finf       ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
     double precision :: pexinf     ! large p tau(p) = tauinf/p**pexinf     
     integer          :: mode                                          
! the recin parameter representation 
     double precision :: q          ! q                                                                               
! the reout parameter representation 
     double precision :: Rg         ! predicted ring radius of gyration                                               
 
     double precision :: th
 
     double precision :: t
     double precision :: rr, nu, Sq, Sqt
     double precision :: msd, msd_im
     double precision :: tau_R, W, Ne0, taue

     logical          :: parameter_is_q, parameter_is_t

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ring2am4'
       nparx =       18
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ring2am4 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen"//cr//parspace//&
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
        parnam (10) = 'nue     '  ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
        parnam (11) = 'wl4     '  ! Rouse rate in [A**4/ns]                                                         
        parnam (12) = 'pmin    '  ! transition mode number between simple ring-Rouse and large p modification       
        parnam (13) = 'pwidth  '  ! sharpness of transition                                                         
        parnam (14) = 'f0      '  ! prefactor f(p) limit for small p values (default 1)                             
        parnam (15) = 'finf    '  ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
        parnam (16) = 'pexinf  '  ! large p tau(p) = tauinf/p**pexinf                                               
        parnam (17) = 'mode    '  ! mode sqt oder msd                                               
        parnam (18) = 'ne0     '  ! experimental ne0                                              
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
        th_param_desc(11,idesc) = "ABS: Rouse rate in [A**4/ns]" !//cr//parspace//&
        th_param_desc(12,idesc) = "ABS: transition mode number between simple ring-Rouse and large p modification" !//cr//parspace//&
        th_param_desc(13,idesc) = "ABS: sharpness of transition" !//cr//parspace//&
        th_param_desc(14,idesc) = "ABS: prefactor f(p) limit for small p values (default 1)" !//cr//parspace//&
        th_param_desc(15,idesc) = "ABS: prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwidth" !//cr//parspace//&
        th_param_desc(16,idesc) = "ABS: large p tau(p) = tauinf/p**pexinf" !//cr//parspace//&
        th_param_desc(17,idesc) = "0: sqt, 1: msd, 2:msd-intern" !//cr//parspace//&
        th_param_desc(18,idesc) = "0: Ne0 from N/pmin else Ne0" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "Rg       > predicted ring radius of gyration"
! 
        th_ring2am4 = 0.0
 
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
      nue      =    abs(  pa(10) )
      wl4      =    abs(  pa(11) )
      pmin     =    abs(  pa(12) )
      pwidth   =    abs(  pa(13) )
      f0       =    abs(  pa(14) )
      finf     =    abs(  pa(15) )
      pexinf   =    abs(  pa(16) )
      mode     =   nint(  pa(17) )
      ne0      =    abs(  pa(18) )
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
      th_ring2am4 = 0.0
      return
    endif
    if(.not.parameter_is_q .and. .not.parameter_is_t) then
      write(*,*)"th_ring2am4: record parameters of ",iadda," must contain either q or t "
      th_ring2am4 = 0.0
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

     rr=((exp(((-alpha+beta)*log(r021/r022)+alpha*beta*(log(2d0)+log(3d0)-      &
         log(r022/diff)))/beta)*r022*t**alpha)**a_cross+(exp((log(2d0)+log(3d0)- &
         log(r022/diff))*beta)*r022*t**beta)**a_cross+(6d0*diff*t)**a_cross)**(1d0/a_cross)

     nu  = 2*nue

     call Nrouse_ring2am4(q,t,n)
     th =  ampli *  exp(-rr * q*q / 6d0) * Sqt/Sq

     msd = msd_im + rr

     if(mode == 1) then
       th =  msd
       if(parameter_is_t) th = ampli *  exp(-rr * q*q / 6d0) * Sqt
     else if (mode == 2) then
       th = msd_im
       if(parameter_is_t) th = ampli * Sqt
     else if (mode == 3 ) then
       th = ampli*exp(-msd * q*q / 6d0) 
     else
       th =  ampli *  exp(-rr * q*q / 6d0) * Sqt/Sq
       if(parameter_is_t) th = ampli *  exp(-rr * q*q / 6d0) * Sqt
     endif
!
     th_ring2am4 = th
 
! ---- writing computed parameters to the record >>>  
      call parset('Rg      ',sngl(Rg)   ,iadda,ier)
      call parset('tau_R   ',sngl(tau_R), iadda, ier)
      call parset('taue    ',sngl(taue) , iadda, ier)
      call parset('W       ',sngl(W)    , iadda, ier)
      call parset('Ne0     ',sngl(Ne0)  , iadda, ier)
      call parset('tau_pmin',sngl(tau_R/pmin**2)  , iadda, ier)
      call parset('tau_px  ',sngl(1/( (W*pi**2) /Ne0**2)  )  , iadda, ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

! subroutines and functions entered here are private to this theory and share its variables

       subroutine Nrouse_ring2am4(q,t,N)
!      ================================================================
!
! Rouse expression for a ring polymer
! with Nb segments each, the ideal Ree of the linear version of the polymer is given by R.
! amod contains mode modifiers
! Input parameters: (polymer and model descriptor are contained shared with main )
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
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



       double precision, intent(in)     ::  q,t
       integer, intent(in)              ::  N

       integer                          ::  nn,mm,p

       double precision :: rate, traf, nun ! , tau_R, W, Ne0, taue


       double precision :: cosarray(0:N,N/2), ewfac(N/2), ff2(N/2)




       integer          :: ip, N2


       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- and the Rousetime ----
       W     = Wl4 / l**4
       tau_R = N**2/(W * Pi**2)
       if(Ne0 <= 0d0)  Ne0   = N/pmin
       taue  = Ne0**2 * l**4/(Wl4*pi**2)
!!       tauinf= taue * pmin**pexinf !! taue * (N/Ne0)**pexinf 
!!       rp=(Wl4*pi^2)*(p/pmin)**pexp/Ne0^2.

! ---- init sums ----
       N2 = N/2
! p(even) = 2*ip in the following...

!$OMP PARALLEL DO PRIVATE(traf,nun)
       do nn=0,N
        do ip=1,N2
         traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
         nun  = traf * nu + (1d0-traf)
         ff2(ip) =  -2d0* dble(N)**nun *(l*q)**2/(3d0*pi**2)
         cosarray(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))               &
                          /  dble(2*ip)**(2d0)                         &
                          *  (f0*(1d0-traf)+finf*(traf))
        enddo
       enddo
!$OMP END PARALLEL DO


!       ff2  = -2d0* dble(N)**nu *(l*q)**2/(3d0*pi**2)

       msd_im = 0
!$OMP PARALLEL DO PRIVATE(rate,traf) REDUCTION(+:msd_im)
!!$OMP PARALLEL DO PRIVATE(rate,traf)
       do ip=1,N/2
         traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))

         rate      = dble(2d0*ip)**(2d0)/tau_R*(1d0-traf) + &
                     (W*pi**2)* (2d0*ip/pmin)**pexinf/Ne0**2 * traf

         ewfac(ip) = 1.0d0-exp(-t * rate )


         msd_im = msd_im - ff2(ip)  /  dble(2*ip)**(1d0+1.d0)                      &
                          *  (f0*(1d0-traf)+finf*(traf))                       &
                          *  ewfac(ip)
       enddo
!$OMP END PARALLEL DO
       msd_im = 6 * msd_im / q**2


           Sq  = 0
           Sqt = 0
           Rg  = 0


! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt,Rg)
           do nn = 1,N
            do mm = 1,N


              Sq  = Sq  + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu )
              Sqt = Sqt + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu  &
                              + sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)*ff2(1:N2)) )

              Rg  = Rg  + l**2  * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu

            enddo
           enddo
!$OMP END PARALLEL DO

           Sq  = Sq /N
           Sqt = Sqt/N

           Rg  = sqrt(Rg/2d0/N**2)

       end subroutine Nrouse_ring2am4

 end function th_ring2am4
