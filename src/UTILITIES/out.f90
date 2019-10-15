 FUNCTION th_nvilgis(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Scattering factor of a Rouse  chain in a harmonic potential. Using direct summation over an effective number of beads, however, using the infinite chain rms distance form Eq(24) of the publication.
!  T.A. Vilgis and F. Boue, Journal of Polymer Science Part B Polymer Physics 26, 2291-2301 (1988)
      use theory_description 
      implicit none 
      real    :: th_nvilgis
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
     double precision :: n          ! number of segments (do not fit, its integer)                                    
     double precision :: l          ! effective segment length Rg = l * N**nu /sqrt(6)                                
     double precision :: wl4        ! Rouse rate                                                                      
     double precision :: rmesh      ! effective potential parameter in terms mesh size, see paper                     
     double precision :: tmesh      ! lifetime of the rmesh constraint (experimental)                                 
! the recin parameter representation 
     double precision :: q          ! q-value    default value                                                        
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision :: th
 
   double precision :: q0
   integer          :: ni
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nvilgis'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nvilgis = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Scattering factor of a Rouse  chain in a harmonic potential. Using direct summation over an effective number of beads, however, using the infinite chain rms distance form Eq(24) of the publication."
       th_citation(idesc)     = " T.A. Vilgis and F. Boue, Journal of Polymer Science Part B Polymer Physics 26, 2291-2301 (1988)"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments (do not fit, its integer)                                    
        parnam ( 3) = 'l       '  ! effective segment length Rg = l * N**nu /sqrt(6)                                
        parnam ( 4) = 'wl4     '  ! Rouse rate                                                                      
        parnam ( 5) = 'rmesh   '  ! effective potential parameter in terms mesh size, see paper                     
        parnam ( 6) = 'tmesh   '  ! lifetime of the rmesh constraint (experimental)                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 4,idesc) = "Rouse rate" !//cr//parspace//&
        th_param_desc( 5,idesc) = "effective potential parameter in terms mesh size, see paper" !//cr//parspace//&
        th_param_desc( 6,idesc) = "lifetime of the rmesh constraint (experimental)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_nvilgis = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        =      pa( 2)
      l        =      pa( 3)
      wl4      =      pa( 4)
      rmesh    =      pa( 5)
      tmesh    =      pa( 6)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     ni  = nint(n)
     rg  =  l * dble(ni)**nu / sqrt(6d0)
     th  = ampli * sq


     th_nvilgis = th
 
! ---- writing computed parameters to the record >>>  
      call parset('rg      ',sngl(rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


       subroutine NrouseP(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    pmin  ----> minimum p
!    pmax  ----> maximum p
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4, pmin, pmax
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)

       integer :: ipmin, ipmax, i

!       integer iout

       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/N)

! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

!$OMP PARALLEL DO
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip
        enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
       do i=1,N
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*i)/dfloat(N)))*t))
       enddo
!$OMP END PARALLEL DO

       ipmin = max(1,nint(pmin))
       ipmax = min(N,nint(pmax))

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w

       return
       end





 end function th_nvilgis
