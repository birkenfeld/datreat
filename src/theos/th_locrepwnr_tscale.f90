 FUNCTION th_locrep2nrt(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrep2nrt
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
     double precision :: b          ! fluctuation intensity (relative)
     double precision :: a          ! length scale                    
     double precision :: tau        ! timescale                       
     double precision :: lz         ! total length                    
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: sqt, sqt0
     double precision :: teps = 1d-6
     double precision :: re   
     integer          :: ne   
     double precision :: temp  
     double precision :: wl4 , wlmod
     double precision :: tscale
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'lrep2nrt'
       nparx =        10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrep2nrt = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " but with finite summation of integrals and lenthscale," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters" //cr//parspace//&
                                " Multiplied with a ROUSE blob, with Wl4 connected to a**4 * tau "//cr//parspace//&
                                " ACCELERATED (OMP+IMPROVED ALG) VERSION!"

       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'b       '  ! fluctuation intensity (relative)                                                
        parnam ( 3) = 'a       '  ! length scale                                                                    
        parnam ( 4) = 'tau     '  ! timescale                                                                       
        parnam ( 5) = 'lz      '  ! total length                                                                    
        parnam ( 6) = 'teps    '  ! near zero for t=0 evaluation   
        parnam ( 7) = 're      '  ! blob re ca. dtube
        parnam ( 8) = 'n       '  ! effective number of segments in blob
        parnam ( 9) = 'wlmod   '  ! modification factor for wl4 in Rosueblob
        parnam (10) = 'tscale  '  ! modification factor for wl4 in Rosueblob
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fluctuation intensity (relative)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 5,idesc) = "total length" !//cr//parspace//&
        th_param_desc( 6,idesc) = "teps near zero for t=0 evaluation" !//cr//parspace//&
        th_param_desc( 7,idesc) = "effective blob Re (ca. tube)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "n effective segments in blob" !//cr//parspace//&
        th_param_desc( 9,idesc) = "wl4 modification factor in Rosue blob (def=1) (0==>1)" !//cr//parspace//&
        th_param_desc(10,idesc) = "time scaling factor" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(1,idesc)   = "a2sqt   = a**2/sqrt(tau) "
        th_out_param(2,idesc)   = "tau0_lr = tau/a**4 "
        th_out_param(3,idesc)   = "Wl4 = a2sqt**2 "
! 
        th_locrep2nrt = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      b        =   abs( pa( 2) )
      a        =   abs( pa( 3) )
      tau      =   abs( pa( 4) )
      lz       =   abs( pa( 5) )
      teps     =        pa( 6)
      re       =   abs( pa( 7) )
      ne       =  nint( pa( 8) )
      wlmod    =  abs ( pa( 9) )
      tscale   =  abs ( pa(10) )
      if (teps <= 0d0) teps = 1d-6
      if(wlmod == 0d0) wlmod = 1d0
      if(tscale == 0d0) tscale = 1d0

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      xh =     400d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh

      t        = x * tscale
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     sqt0 =  local_reptation2( q*a, teps,  b, lz)
     sqt  =  local_reptation2( q*a, t/tau, b, lz)
     th_locrep2nrt = ampli * sqt/sqt0

     wl4 = a*a*a*a/tau  * wlmod
     call  NrousePbl(q,t,temp,wl4,ne,re, sqt0,sqt)

!write(*,*)"Nrouse:",q,t,temp,wl4,ne,re,sqt0,sqt

     th_locrep2nrt = th_locrep2nrt * sqt/sqt0

    

!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda)
      call parset("tau0_lr ", sngl(tau/a**4),iadda)
      call parset("wl4     ", sngl(wl4),iadda)
                   
 CONTAINS

! subroutines and functions entered here are private to this theory and share its variables 
 
  function local_reptation2(q, t, B, L) result(val)
    implicit none
    double precision, intent(in)   :: q, t, B, L
    double precision, parameter    :: A = 1d0
    double precision, parameter    :: Pi = 3.141592653589793d0 
    double precision               :: val
    double precision, parameter    :: sqp = sqrt(4*atan(1d0))

    double precision               :: ec1, ec2, dec, edec, z


    z   = sqrt(t) * q**2

    ec1 = erfc((q ** 2 * t + 3 * L) * t ** (-0.5d0) / 6d0)
    ec2 = erfc(sqrt(t) * q ** 2 / 6d0)

    if( isnan(ec1) ) ec1 = 1d0 - erf((q ** 2 * t + 3 * L) * t ** (-0.5d0) / 6d0)
    if( isnan(ec2) ) ec2 = 1d0 - erf(sqrt(t) * q ** 2 / 6d0)

    dec  =   (-ec1 + ec2 )
    edec = exp(t * q ** 4 / 36d0) *  ( dec ) 

    if( isnan(edec) .or. z > 50d0 ) then
          edec = - ( -6/(sqrt(Pi)*z)+108/(sqrt(Pi)*z**3)-5832/(sqrt(Pi)*z**5) )
    endif

      val= 72d0 * (sqrt(t) * q ** 4 * exp((-2d0 * L * q ** 2 * t - 3 * L ** 2) / t / 12d0) / 36d0 + &
           sqrt(Pi) * (q ** 2 * t / 3 + L) * q ** 4 * edec / 72d0 - sqrt(t) * q ** 4 / 36d0) &
           * B / q ** 4 * Pi ** (-0.5d0) / L +  72d0 * (A * exp(-q ** 2 * L / 6d0) * sqrt(Pi) &
           + (A * L * q ** 2 / 6d0 - A) * sqrt(Pi)) * Pi ** (-0.5d0) / q ** 4 / L

  end function local_reptation2






       subroutine NrousePbl(q,t,temp,wl4,N,R, Sq,Sqt)
!      ==============================================
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
! Output parameters:
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in)  :: q
       double precision, intent(in)  :: t
       double precision, intent(in)  :: temp
       double precision, intent(in)  :: R
       integer         , intent(in)  :: N
       double precision, intent(out) :: Sq ,Sqt
       double precision, intent(in)  :: wl4

       integer :: nn,mm,ip

       double precision :: l, tau_p, kbt, arg1, arg2
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2
       double precision :: p
       double precision :: xi, W
       double precision :: Dr 

       double precision :: cosarray(N,N), ewfac(N)


      
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
       Dr  = 0d0
! ---- init sums ----

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip 
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO    
       do ip=1,N
         ewfac(ip) = (1d0-exp(-2*W*(1-cos((pi*ip)/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    




       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)


! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                ff2* sum(cosarray(nn,:) * cosarray(mm,:) *  ewfac(:) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!       write(6,'(1x,a,7E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w, l', q,t,Sq,Sqt, Sqt/Sq, w, l 

       return
    
 end subroutine NrousePbl


end function th_locrep2nrt

