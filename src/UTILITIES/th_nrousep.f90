 FUNCTION th_nrousep(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Rouse by discrete summation with mode restriction
!  Rouse, Doi_Edwards
      use theory_description 
      implicit none 
      real    :: th_nrousep
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
     double precision :: pmin       ! minimum mode to be included                                                     
     double precision :: pmax       ! maximum mode to be included                                                     
! the recin parameter representation 
     double precision :: q          ! momentum transfer                                                               
                                                              
! the reout parameter representation 
     double precision :: Dr         ! diffusion   
     double precision :: W          ! rouse     W
     double precision :: l          ! sqgment length    
                                                                    
     double precision :: sqt, sq, t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nrousep'
       nparx =        8
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nrousep = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Rouse by discrete summation with mode restriction"
       th_citation(idesc)     = " Rouse, Doi_Edwards"
!       --------------> set the parameter names --->
        parnam ( 1) = 'amplitu '  ! prefactor                                                                       
        parnam ( 2) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 3) = 'n_segmen'  ! effective N for summation                                                       
        parnam ( 4) = 're      '  ! end-to-end radius of Gaussian coil                                              
        parnam ( 5) = 'temp    '  ! temperature if not in parameters                                                
        parnam ( 6) = 'com_diff'  ! center of mass diffusion (if 0, use Rouse default)                              
        parnam ( 7) = 'pmin    '  ! minimum mode to be included                                                     
        parnam ( 8) = 'pmax    '  ! maximum mode to be included                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective N for summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature if not in parameters" !//cr//parspace//&
        th_param_desc( 6,idesc) = "center of mass diffusion (if 0, use Rouse default)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "minimum mode to be included" !//cr//parspace//&
        th_param_desc( 8,idesc) = "maximum mode to be included" !//cr//parspace//&
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
        th_nrousep = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      amplitu  =      pa( 1)
      wl4      =  abs(pa( 2))
      n_segmen = nint(pa( 3))
      re       =  abs(pa( 4))
      temp     =      pa( 5)
      com_diff =  abs(pa( 6))
      pmin     =      pa( 7)
      pmax     =      pa( 8)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: momentum transfer
      xh =     0.1d0 
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh =    300d0 
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

     Dr  = com_diff * 1d-9 / 1d-16  ! in A**2/ns


     call NrouseP(q,t,temp,Dr,wl4,n_segmen,Re, W, l,pmin,pmax, Sq,Sqt)


     th_nrousep = amplitu * sqt/sq

 
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
      call        parset('w       ',sngl(W),iadda,ier)      !     "
      call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
  
      Dr        = Dr /( 1d-9 / 1d-16 ) ! in cm**2/s
      call parset('diff    ',sngl(Dr),iadda,ier)
 
 CONTAINS 
 

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

       integer :: ipmin, ipmax

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

    
       p0fix = 0
       pfac  = 1

       ipmin = max(1,nint(pmin))
       ipmax = min(N,nint(pmax))

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)
    
          arg2 = 0
          arg20= 0
          do ip=ipmin, ipmax
            p = ip+p0fix
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))
            a0    = -t*tau_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)
            ffc   = ffc*pfac

            arg2  = arg2  + ffc*e0
            arg20 = arg20 + ffc

          enddo   
          arg2  = arg2  * ff2
          arg20 = arg20 * ff2

          aa1 = arg10
          if(aa1.lt.-300.0d0) aa1 = -300.0d0
          if(aa1.gt. 300.0d0) aa1 =  300.0d0

          aa2 = arg1+arg2
          if(aa2.lt.-300.0d0) aa2 = -300.0d0
          if(aa2.gt. 300.0d0) aa2 =  300.0d0


          Sq  = Sq  + exp(aa1)
          Sqt = Sqt + exp(aa2)

        enddo
       enddo

       Sq  = Sq /N
       Sqt = Sqt/N

 !      if(iout().gt.1)write(6,'(1x,a,6E14.6)')
 !    *        'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

 end function th_nrousep
