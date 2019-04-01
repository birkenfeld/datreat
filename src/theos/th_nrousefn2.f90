 FUNCTION th_nrousefn2(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Rouse by discrete summation with mode restriction
!  Rouse, Doi_Edwards
      use theory_description 
      implicit none 
      real    :: th_nrousefn2
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
     double precision :: amplitu    ! prefactor                                                                       
     double precision :: wl4        ! rouse rate                                                                      
     integer          :: n_segmen   ! effective N for summation                                                       
     double precision :: re         ! end-to-end radius of Gaussian coil                                              
     double precision :: temp       ! temperature if not in parameters                                                
     double precision :: com_diff   ! center of mass diffusion (if 0, use Rouse default)                              
     double precision :: q_width    ! q_averaging width                                                  
     double precision :: re0        ! re at t=0                                                  
     double precision :: re_inf     ! re at t=inf                                                  
     double precision :: dtau       ! relaxation time re0-> re_inf                                                  
     double precision :: beta_re    ! streched exp                                                  

     integer          :: ifx        ! end fixed (1) or not (0)                                                  
     integer          :: nmin       ! summation start for mix of varoius lengths 
! the recin parameter representation 
     double precision :: q          ! momentum transfer                                                               
                                                              
! the reout parameter representation 
     double precision :: Dr         ! diffusion   
     double precision :: W          ! rouse     W
     double precision :: l          ! sqgment length    
                                                                    
     double precision :: sqt, sq, t
     double precision :: sums, sumnorm, sofq0
     double precision :: qzz, Rn, fn, dqw
     integer          :: Nn, nqw, i
     logical          :: sqof0
     real             :: tget, qget

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nrousefy'
       nparx =        12
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nrousefn2 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Rouse by discrete summation with mode restriction, fix ends and var length" & 
                //cr//parspace//" and expewrimental time dependent Re value " &
                //cr//parspace//" replaces nrousefn, ..fx, ACCELEARTED (OMP+IMPROVED ALG.) VERSION"

       th_citation(idesc)     = " Rouse, Doi_Edwards"
!       --------------> set the parameter names --->
        parnam ( 1) = 'intensit'  ! prefactor    
        parnam ( 2) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 3) = 'n_segmnt'  ! effective N for summation                                                       
        parnam ( 4) = 're      '  ! end-to-end radius of Gaussian coil                                              
        parnam ( 5) = 'temp    '  ! temperature if not in parameters                                                
        parnam ( 6) = 'com_diff'  ! center of mass diffusion (if 0, use Rouse default)                              
        parnam ( 7) = 'q_width '  ! minimum mode to be included                                                     
        parnam ( 8) = 'fixend  '  ! fix end or not                                                    
        parnam ( 9) = 'nmin    '  ! start of length summation                                                   
        parnam (10) = 're_inf  '  ! infinite time re value
        parnam (11) = 'dtau    '  ! timeconstant re relaxation                                                  
        parnam (12) = 'beta_re '  ! streched exp re relaxation                                                   
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective N for summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature if not in parameters" !//cr//parspace//&
        th_param_desc( 6,idesc) = "center of mass diffusion (if 0, use Rouse default)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "q-integration width (default=0)" !//cr//parspace//&
        th_param_desc( 8,idesc) = "fixed-end=1, free-chain=0" !//cr//parspace//&
        th_param_desc( 9,idesc) = "start of length summation" !//cr//parspace//&
        th_param_desc(10,idesc) = "infinite time re value   " !//cr//parspace//&
        th_param_desc(11,idesc) = "timeconstant re relaxation" !//cr//parspace//&
        th_param_desc(12,idesc) = "streched exp re relaxation" !//cr//parspace//&
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "diff     > diffusion in cm**2/s"
        th_out_param(  2,idesc) = "l        > effective segment length"
        th_out_param(  3,idesc) = "w        > rouse rate: W"
        th_out_param(  4,idesc) = "wl4      > inferred value of Wl4"
! 
        th_nrousefn2 = 0.0
 
        RETURN
     ENDIF
!
      sqof0 = .false.

! ---- transfer parameters -----
      amplitu  =      pa( 1)
      wl4      =  abs(pa( 2))
      n_segmen = nint(pa( 3))
      re0      =  abs(pa( 4))
      temp     =      pa( 5)
      com_diff =  abs(pa( 6))
      q_width  =      pa( 7)
      ifx      = nint(pa( 8))
      nmin     = nint(pa( 9))
      re_inf   = abs(pa(10))
      dtau     = abs(pa(11))
      beta_re  = abs(pa(12))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: momentum transfer
!      xh =     0.1d0 
!      call parget('q       ',xh,iadda,ier)
!      q        = xh
! >>> extract: temperature
!      xh =    300d0 
!      call parget('temp    ',xh,iadda,ier)
!      if( xh > 0d0) temp     = xh
 
     t   = x


       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       q   = qget
       if(ier.ne.0) then 
         write(6,*)'Assume Q as independend variable!'
         q   = x
         t   = 0.0d0
         call parget('tau     ',qget,iadda,ier)
         if(ier.eq.0) then
          t  = qget
         endif
         write(6,*)'tau = ',t 
         sqof0 = .true.
       endif
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         if( tget > 0d0) temp = tget
       endif

! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
!

      th_nrousefn2     = 0
      sums           = 0
      sumnorm        = 0
      nqw            = 15
      dqw            = 4*q_width/nqw
      if(q_width.eq.0) then
        nqw = 0
        q_width = 1.0d0
      endif 

      re = re_inf + (re0-re_inf) * exp(-(t/dtau)**beta_re)

      do i=-nqw,nqw
       qzz = q + i*dqw
       if(qzz.gt.0) then 
         fn   = dexp( -(i*dqw/q_width)**2)

! --- include center of mass diffusion ---
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
!          dr = diff
         Dr  = com_diff * 1d-9 / 1d-16  ! in A**2/ns
   
          do Nn=nmin, N_segmen      ! loop over partial chains

            Rn = Re * sqrt(dble(Nn)/dble(N_segmen))
            call NrousePX(qzz,t,temp,dr,wl4,Nn,Rn, W, l,ifx, Sq,Sqt)
            sums = sums + fn*Sqt       
            sumnorm = sumnorm + Sq*fn

          enddo

        endif
       enddo

       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            th_nrousefn2 = sums
          else
            th_nrousefn2 = sums/sumnorm
          endif
       else
          th_nrousefn2 = 0
       endif 
       th_nrousefn2 = th_nrousefn2*amplitu







 
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
      call        parset('w       ',sngl(W),iadda,ier)      !     "
      call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
  
      Dr        = Dr /( 1d-9 / 1d-16 ) ! in cm**2/s
      call parset('diff    ',sngl(Dr),iadda,ier)
 
 CONTAINS 
 

       subroutine NrousePX(q,t,temp,Dr,wl4,N,R, W, l,ifx , Sq,Sqt)
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
!    ifx   ----> fixend or not
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

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4
       integer N, nn,mm,ifx,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix

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

       if(ifx.eq.0) then
         p0fix = 0
       else
         p0fix = -0.5d0
       endif

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*(ip+p0fix)*nn)/dfloat(N)) / (ip+p0fix)
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO PRIVATE(p)    
       do i=1,N 
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*(i+p0fix))/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    

       ipmin = 1
       ipmax = N

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0


! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt,arg2)

       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)
    
          arg2 = 0
          arg20= 0

          do ip=ipmin, ipmax
            
!            ffc   = cosarray(nn,ip) * cosarray(mm,ip) 

!            arg2  = arg2 +  (1d0-exp(-2*W*(1-cos((pi*ip)/dfloat(N)))*t)) * ffc
            arg2  = arg2 +  ewfac(ip) * cosarray(nn,ip) * cosarray(mm,ip) 

                                
!            arg20 =  arg20 + ffc

          enddo 
          

!          arg2  = arg2  * ff2
!          arg20 = arg20 * ff2

!          aa1 = arg10

!          aa2 = arg1+ff2*arg2


          Sq  = Sq  + exp(arg10)
          Sqt = Sqt + exp(arg1+ff2*arg2)

        enddo

       enddo

!$OMP END PARALLEL DO



       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end subroutine NrousePX
 

 end function th_nrousefn2
