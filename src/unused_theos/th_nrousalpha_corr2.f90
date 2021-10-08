 FUNCTION th_nrousalpha(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Rouse by discrete summation with mode restriction
!  Rouse, Doi_Edwards
      use theory_description 
      implicit none 
      real    :: th_nrousalpha
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
     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset

     integer          :: mode

 !  a = alpha0 * exp(-((1d0-log(t/talphamax))/talphawd)**2) + aoffset



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
       thnam = 'nrousalP'           ! fqq0 correct version of nrousealp
       nparx =        13
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nrousalpha = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Rouse by discrete summation with mode restriction"//cr//parspace//&
                                " with non-Gauss alpha (Guenza2014) "//cr//parspace//&
                                " alpha(1) =alpha0*exp(-(((log(t+eps)-log(talpmax))/talpwd)**2)/2)+alpoff"//cr//parspace//&
                                " (eps=1e-3),            ACCELERATED (OMP+IMPOVED ALG.) VERSION!"

       th_citation(idesc)     = " Rouse, Doi_Edwards + M-.Guenza PHYSICAL REVIEW E 89, 052603 (2014)"
!       --------------> set the parameter names --->
        parnam ( 1) = 'amplitu '  ! prefactor                                                                       
        parnam ( 2) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 3) = 'n_segmen'  ! effective N for summation                                                       
        parnam ( 4) = 're      '  ! end-to-end radius of Gaussian coil                                              
        parnam ( 5) = 'temp    '  ! temperature if not in parameters                                                
        parnam ( 6) = 'com_diff'  ! center of mass diffusion (if 0, use Rouse default)                              
        parnam ( 7) = 'pmin    '  ! minimum mode to be included                                                     
        parnam ( 8) = 'pmax    '  ! maximum mode to be included                                                     
        parnam ( 9) = 'alpha0  '  !                                                     
        parnam (10) = 'talpmax '  !                                                     
        parnam (11) = 'talpwd  '  !                                                    
        parnam (12) = 'alpoffs '  !                                                     
        parnam (13) = 'mode    '  !                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective N for summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature if not in parameters" !//cr//parspace//&
        th_param_desc( 6,idesc) = "center of mass diffusion (if 0, use Rouse default)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "minimum mode to be included" !//cr//parspace//&
        th_param_desc( 8,idesc) = "maximum mode to be included" !//cr//parspace//&
        th_param_desc( 9,idesc) = "prefactor of non Gaussianity alpha(t) expression" !//cr//parspace//&
        th_param_desc(10,idesc) = "t of max in log distribution alpha(t)" !//cr//parspace//&
        th_param_desc(11,idesc) = "width of Gaussion max. in log t of alpha(t)" !//cr//parspace//&
        th_param_desc(12,idesc) = "offset in alpha(t)" !//cr//parspace//&
        th_param_desc(13,idesc) = "mode=0: normal sqt, mode=1: alpha" !//cr//parspace//&
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
        th_nrousalpha = 0.0
 
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

      alpha0   =      pa( 9)
      talphamax=  abs(pa(10))
      talphawd =  abs(pa(11))
      aoffset  =      pa(12)
      mode     = nint(pa(13))



! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: momentum transfer
      xh =     0.1d0 
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh =    300d0 
      call parget('temp    ',xh,iadda,ier)
      if( xh > 0d0) temp     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

     Dr  = com_diff * 1d-9 / 1d-16  ! in A**2/ns

     select case(mode)
       case(0)
         call nrousalpha(q,t,temp,Dr,wl4,n_segmen,Re, W, l,pmin,pmax, Sq,Sqt)
         th_nrousalpha = amplitu * sqt/sq
       case(1)
         th_nrousalpha = alpha(t)
       case default 
         write(*,*)"th_nroualp: mode must be 0 or 1"
     end select

 
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
      call        parset('w       ',sngl(W),iadda,ier)      !     "
      call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
  
      Dr        = Dr /( 1d-9 / 1d-16 ) ! in cm**2/s
      call parset('diff    ',sngl(Dr),iadda,ier)
 
 CONTAINS 
 

       subroutine nrousalpha(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
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
       double precision :: rmm, fqq, fqq0

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

       rmm = 0
!$OMP PARALLEL DO REDUCTION(+:rmm)
       do mm = 1,N
             rmm = rmm + 4d0*N*l**2/(pi**2) * &
                   sum(cosarray(mm,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )
       enddo
!$OMP END PARALLEL DO
       rmm = rmm/N

       fqq  = 1d0 -q**2 * rmm/12d0 * alpha(t)       !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
!       fqq0 = 1d0 -q**2 * rmm/12d0 * alpha(0d0)     !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
       fqq0 = 1

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(- fqq0*(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(- fqq* (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

!          Sq  = Sq  + exp(- (q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
!          Sqt = Sqt + exp(- (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
!                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
   double precision, intent(in) :: t
   double precision             :: a

! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
!                                        are shared (common) with those of the th-function
! the model
   a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   
end function alpha






 end function th_nrousalpha
