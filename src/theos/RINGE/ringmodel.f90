MODULE ring

 double precision :: Rg
 double precision :: tau_R
 double precision :: taue
 double precision :: W
 double precision :: Ne0
 double precision :: tauring
 

CONTAINS

function ring2am7(q, t, &
                  diff, r021, alpha, r022, beta, a_cross, &
                  n, l, nue, nu0, wl4, pmin ,pwidth, ntrans, nwidth, &
                  f0, finf, pexinf, ne_0)                       result(sq0t)
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
 !      parameter(kb=1.380662d-23)
       double precision, parameter   :: pi=4d0*atan(1d0)


       double precision, intent(in) :: q
       double precision, intent(in) :: t
!
       double precision, intent(in) :: diff     ! limiting diffusion (NMR value) in [A**2/ns]                                     
       double precision, intent(in) :: r021     ! reference mean squared displacement at first transition point to medium timed di
       double precision, intent(in) :: alpha    ! sub-diffusion  exponent of SHORT time diffusion                                 
       double precision, intent(in) :: r022     ! reference mean squared displacement at second transition point to D0 in [A**2]  
       double precision, intent(in) :: beta     ! sub-diffusion  exponent of medium time diffusion                                
       double precision, intent(in) :: a_cross  ! transition exponent between short and long time diffusion (sharper kink for larg
       integer         , intent(in) :: n        ! number of segments in one ring                                                  
       double precision, intent(in) :: l        ! effective segment length                                                        
       double precision, intent(in) :: nue      ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
       double precision, intent(in) :: nu0      ! chain statistics exponent (nu=0.5 => random walk, Gaussian)                     
       double precision, intent(in) :: wl4      ! Rouse rate in [A**4/ns]                                                         
       double precision, intent(in) :: pmin     ! transition mode number between simple ring-Rouse and large p modification       
       double precision, intent(in) :: pwidth   ! sharpness of transition     
       double precision, intent(in) :: ntrans   ! structure transitions from nu0 to nue
       double precision, intent(in) :: nwidth   ! structure transitions width
       double precision, intent(in) :: f0       ! prefactor f(p) limit for small p values (default 1)                             
       double precision, intent(in) :: finf     ! prefactor f(p) limit for large p values (default F=0.9??) transitin width is pwi
       double precision, intent(in) :: pexinf   ! large p tau(p) = tauinf/p**pexinf     
       double precision, intent(in) :: ne_0     ! experimental ne0                                              
   
       double precision                 :: sq, sqt, sq0t(2)
! -- local? vars
       double precision :: msd_im
       double precision :: rr


       integer          ::  nn,mm,p
       double precision :: rate, traf, nun ! , tau_R, W, Ne0, taue
       double precision :: cosarray(0:N,N/2), ewfac(N/2), ff2(N/2), taus(N/2)
       double precision :: nu(0:N)
       double precision :: tram(0:N)

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
       Ne0   = ne_0
       if(Ne0 <= 0d0)  Ne0   = N/pmin
       taue  = Ne0**2 * l**4/(Wl4*pi**2)
!!       tauinf= taue * pmin**pexinf !! taue * (N/Ne0)**pexinf 
!!       rp=(Wl4*pi^2)*(p/pmin)**pexp/Ne0^2.

! ---- init sums ----
       N2 = N/2
! p(even) = 2*ip in the following...a

!$OMP PARALLEL DO
  do nn=0,N
    tram(nn) = 1d0/(1.0d0+exp((nn-ntrans)/nwidth))
!                                 -----   ------
  enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do nn=0,N
    nu(nn) = 2*(nu0*tram(nn) + (1-tram(nn))*nue)
  enddo
!$OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(traf,nun)
       do nn=0,N
        do ip=1,N2
         traf = 1d0/(1.0d0+exp((2d0*ip-pmin)/pwidth))
         nun  = traf * (2*nue) + (1d0-traf) * (2*nu0)
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

         taus(ip)  = 1d0/rate

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


              Sq  = Sq  + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) )
              Sqt = Sqt + exp(-((q*l)**2)/6d0 * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) &
                              + sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)*ff2(1:N2)) )

              Rg  = Rg  + l**2  * ( dble(abs(nn-mm))*dble((N-abs(nn-mm)))/dble(N))**nu(abs(nn-mm)) 

            enddo
           enddo
!$OMP END PARALLEL DO

           Sq  = Sq /N
           Sqt = Sqt/N

          rr=((exp(((-alpha+beta)*log(r021/r022)+alpha*beta*(log(2d0)+log(3d0)-      &
             log(r022/diff)))/beta)*r022*t**alpha)**a_cross+(exp((log(2d0)+log(3d0)- &
             log(r022/diff))*beta)*r022*t**beta)**a_cross+(6d0*diff*t)**a_cross)**(1d0/a_cross)

           Sqt = Sqt * exp(-rr * q*q / 6d0) 
           sq0t = [sq,sqt]

           Rg  = sqrt(Rg/2d0/N**2)

           tauring = taus(1)

end function ring2am7

END MODULE RING




program test
 use ring
  implicit none


  double precision :: diff     
  double precision :: r021       
  double precision :: alpha      
  double precision :: r022       
  double precision :: beta       
  double precision :: a_cross    
  double precision :: n          
  double precision :: l          
  double precision :: nue        
  double precision :: wl4        
  double precision :: pmin       
  double precision :: pwidth     
  double precision :: f0         
  double precision :: finf       
  double precision :: pexinf     
  double precision :: mode       
!  double precision :: ne0        
  double precision :: nu0       
  double precision :: ntrans     
  double precision :: nwidth  



  double precision :: q, t, s(2)
  integer          :: i


diff    =   0.6200000E-01 
r021    =    1100.000     
alpha   =   0.4463293     
r022    =    2750.000     
beta    =   0.7500000     
a_cross =    8.000000     
n       =    190.0000     
l       =    10.70000     
nue     =   0.4300000     
wl4     =    14890.00     
pmin    =    22.00000     
pwidth  =    2.000000     
f0      =    1.000000     
finf    =    1.000000     
pexinf  =    2.302588     
mode    =    0.000000     
ne0     =    0.000000     
nu0     =   0.5000000     
ntrans  =    11.25000     
nwidth  =    1.000000      


q = 0.1d0

do i=0,50
  t = i
  s = ring2am7(q, t, &
                  diff, r021, alpha, r022, beta, a_cross, &
                  nint(n), l, nue, nu0, wl4, pmin ,pwidth, ntrans, nwidth, &
                  f0, finf, pexinf, ne0)  
  write(*,*) t, s, s(2)/s(1)

enddo

n= 30
t=0
do i=0,50
 q = 0.01*i
   s = ring2am7(q, t, &
                  diff, r021, alpha, r022, beta, a_cross, &
                  nint(n), l, nue, nu0, wl4, pmin ,pwidth, ntrans, nwidth, &
                  f0, finf, pexinf, ne0)  
   write(*,*) q, s(1)
 
enddo

 


end program test
