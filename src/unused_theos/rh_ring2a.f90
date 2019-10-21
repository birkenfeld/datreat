 FUNCTION th_ring2a(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen
!  S. Goossen et al., PRL 2014, 113, 168302
      use theory_description 
      implicit none 
      real    :: th_ring2a
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
     double precision :: tauinf     ! large p tau(p) = tauinf/p**pexinf                                               
     double precision :: pexinf     ! large p tau(p) = tauinf/p**pexinf                                               
! the recin parameter representation 
     double precision :: q          ! q                                                                               
! the reout parameter representation 
     double precision :: Rg         ! predicted ring radius of gyration                                               
 
     double precision :: th
 
     double precision :: t
     double precision :: rr, nu, Sq, Sqt
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ring2a'
       nparx =       17
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ring2a = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Dynamic scattering function for a polymer ring melt as developed in the PRL of S. Goossen"
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
        parnam (16) = 'tauinf  '  ! large p tau(p) = tauinf/p**pexinf                                               
        parnam (17) = 'pexinf  '  ! large p tau(p) = tauinf/p**pexinf                                               
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
        th_param_desc(17,idesc) = "ABS: large p tau(p) = tauinf/p**pexinf" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "Rg       > predicted ring radius of gyration"
! 
        th_ring2a = 0.0
 
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
      tauinf   =    abs(  pa(16) )
      pexinf   =    abs(  pa(17) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q
      xh =   0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     rr=((exp(((-alpha+beta)*log(r021/r022)+alpha*beta*(log(2d0)+log(3d0)-      &
         log(r022/diff)))/beta)*r022*t**alpha)**a_cross+(exp((log(2d0)+log(3d0)- &
         log(r022/diff))*beta)*r022*t**beta)**a_cross+(6d0*diff*t)**a_cross)**(1d0/a_cross)

     nu  = 2*nue

     call Nrouse_ring2a(q,t,n)
     th =  ampli *  exp(-rr * q*q / 6d0) * Sqt/Sq
     th_ring2a = th
 
! ---- writing computed parameters to the record >>>  
      call parset('Rg      ',sngl(Rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

! subroutines and functions entered here are private to this theory and share its variables

       subroutine Nrouse_ring2a(q,t,N)
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

       double precision :: ff2
       double precision :: rate, tau_R, traf, W


       double precision :: cosarray(0:N,N/2), ewfac(N/2)




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

! ---- init sums ----
       N2 = N/2
! p(even) = 2*ip in the following...

!$OMP PARALLEL DO PRIVATE(traf)
       do nn=0,N
        do ip=1,N2
         traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
         cosarray(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))               &
                          /  dble(2*ip)**(1d0+nu)                      &
                          *  (f0*(1d0-traf)+finf*(traf))
        enddo
       enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(rate,traf)
       do ip=1,N/2
         traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
         rate      = dble(2d0*ip)**(2d0*nu)/tau_R*(1d0-traf) + (2d0*ip)**pexinf/tauinf * traf

         ewfac(ip) = 1.0d0-exp(-t * rate )
       enddo
!$OMP END PARALLEL DO

           Sq  = 0
           Sqt = 0
           Rg  = 0

           ff2  = -2d0* dble(N)**nu *(l*q)**2/(3d0*pi**2)

! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt,Rg)
           do nn = 1,N
            do mm = 1,N


              Sq  = Sq  + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu )
              Sqt = Sqt + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu  &
                              + ff2*sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)) )

              Rg  = Rg  + l**2  * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu

            enddo
           enddo
!$OMP END PARALLEL DO

           Sq  = Sq /N
           Sqt = Sqt/N

           Rg  = sqrt(Rg/2d0/N**2)

       end subroutine Nrouse_ring2a

 end function th_ring2a
