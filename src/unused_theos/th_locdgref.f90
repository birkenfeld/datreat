 FUNCTION th_locdgref(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locdgref
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
     double precision :: a          ! length scale  = segemnt length l_seg         
     double precision :: tau        ! timescale           
     double precision :: wl4        ! rouse rate       
     integer          :: n_segmen   ! effective N for Rouse summation   
     double precision :: re         ! end-to-end radius of Gaussian coil                                              
     double precision :: blocrep
     double precision :: dtlocrep
                                                              
! the reout parameter representation 
     double precision :: W0, wr      ! rouse     W
     double precision :: l          ! sqgment length    
     double precision :: b
     double precision :: imod

     double precision :: ne, n
                                                                    

 
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: sqt, sqt0
     double precision :: sqdebt, sqdeb0

     double precision, parameter   :: teps = 1d-6
!     integer                       :: mode
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locdgref'
       nparx =         10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locdgref = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " with fixing of the sq-debye problem  " //cr//parspace//&
                                "  "
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor       
        parnam ( 2) = 'a       '  ! length scale = real segment length   
        parnam ( 3) = 'ne      '  ! entanglemen ne length  
        parnam ( 4) = 'n       '  ! n               
        parnam ( 5) = 'wl4     '  ! rouse rate       
        parnam ( 6) = 'n_segmen'  ! effective N for summation  (scales with rouse)
        parnam ( 7) = 're      '  ! end-to-end radius of Gaussian coil of Rouse blob 
        parnam ( 8) = 'blocrep '  !                                                                 
        parnam ( 9) = 'dtlocrep'  !                                                                 
        parnam (10) = 'imod    '  !                                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "length scale = l_segment" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "total n" !//cr//parspace//&
        th_param_desc( 5,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 6,idesc) = "effective N for Rouse summation" !//cr//parspace//&
        th_param_desc( 7,idesc) = "end-to-end radius of blob Gaussian coil" !//cr//parspace//&
        th_param_desc( 8,idesc) = "extra prefactor to locrep (=b)" !//cr//parspace//&
        th_param_desc( 9,idesc) = "tzero shift of locrep" !//cr//parspace//&
        th_param_desc(10,idesc) = "modulation between eq 43/46 for Lrep" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
!
                                             
! >>>>> describe parameters >>>>>>> 
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "l        > effective segment length (Rouse summation)"
        th_out_param(  2,idesc) = "w0       > rouse rate: W for WL4/a**4"
        th_out_param(  3,idesc) = "wl4      > inferred value of Wl4"
        th_out_param(  4,idesc) = "rga      > sqrt(ne) * a"
!
        th_locdgref = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =        pa( 1)
      a        =   abs( pa( 2) )
      ne       =   abs( pa( 3) )
      n        =   abs( pa( 4) )
      wl4      =   abs( pa( 5))
      n_segmen =  nint( pa( 6))
      re       =   abs( pa( 7))
      blocrep  =   abs( pa( 8))
      dtlocrep =   abs( pa( 9))
      imod     =   abs( pa(10))
 

      b   = blocrep
      W0  = wl4/a**4

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 

!    call NrouseP        (q,t,temp,Dr,wl4,n_segmen,Re, Wr, l,pmin,pmax, Sqdeb0,Sqdebt)
     call nrouseP0(q, t, wl4, N_segmen,Re, Wr, l, Sqdeb0,Sqdebt) 
   
     sqt0  =  local_reptationdr2( q, teps+dtlocrep,         a, W0, n, ne,  b, 0)
     sqt   =  local_reptationdr2( q, t+dtlocrep   ,         a, W0, n, ne,  b, 0)

      th_locdgref = ampli * sqt/sqt0  * Sqdebt/Sqdeb0

!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda)
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda)      ! in ns A units
      call        parset('w0      ',sngl(W0),iadda)      !     "
      call        parset('wl4     ',sngl(W0*l**4),iadda) !     "
      call        parset('rga     ',sngl(sqrt(ne)*a),iadda) !     "
  


CONTAINS


 
  function local_reptationdr2(q, t, lseg, W0, n, ne, b, modelr) result(val)   !! with series expansion
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: lseg       !  lsegment
    double precision, intent(in)   :: W0        ! (Rouse?) rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments/entanglement
    double precision, intent(in)   :: b        ! b-fluctuation intensity
    integer         , intent(in)   :: modelr   ! full deGennes 0, or only T2 (locrep) 1
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)

    double precision :: T1, T2, T3, ec1, ec2, dec
    double precision :: x, xa
    double precision, parameter :: xlimser = 15d0
    double precision :: W

    double precision :: Z, ae



!??!   W = w0*(1d0 - (1d0-ww0)/(1+exp((t-t0)/tww)))
    W = W0
    Z    = n/ne

    T1 = 72d0 * (exp(-q**2 * n * lseg**2 / 6d0) + q**2 * n * lseg**2 / 6d0 - 1d0) / q**4 / lseg**4
!! here we tweak the normalisation to the new interpretation
    T1 = T1/n**2   !! Normalisation to 1 at q=0 

     


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! now T2 the local reptation part
    x = sqrt(lseg**4*q**4*w*t/36)
    xa= n/(2*sqrt(w*t))   !! oder SO ??   Das entspricht Degennes DAS MUSS SO SEIN
    
    if( abs(x) < xlimser ) then
      dec = exp(x**2) *(erfc(x)-erfc(x+xa))
    else
      dec =  1d0/(sqrt(Pi)*x)-1d0/(2*sqrt(Pi)*x**3)+3d0/(4*sqrt(Pi)*x**5) &
          - (x**4-xa*x**3+(xa**2-1d0/2d0)*x**2+(-xa**3+3d0/2d0*xa)*x+xa**4-3*xa**2+3d0/4d0) &
           * exp(-xa**2)*exp(-2*xa*x)/(sqrt(Pi)*x**5)
    endif


    T2 = b*((exp(-1/6*lseg**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
        ne*(n/3+lseg**2*q**2*w*t/9) &
        * dec )

    T2 = T2*3d0/(n*ne*b)    !! Normalisation to T2(t=0) = 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ae  = exp(1d0)
     val = ( T1**ae  + (1d0/Z + (1d0-imod)*b/(3*Z)*T2)**ae ) ** (1d0/ae)  &
            + (imod* b/(3*Z)*T2)
!     val = ( T1**ae + (1d0/Z)**ae ) ** (1d0/ae) + b/(3*Z)*T2   ??
      


  end function local_reptationdr2




       subroutine NrouseP0(q, t, wl4, N, R, W, l, Sq,Sqt)
!      ===================================================
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

       double precision q,t,xi,R, W,Sq,Sqt, wl4
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, Sq0, arg1, arg2, dr
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
       W   = wl4/l**4 ! in 1/ns
       Dr  = 0d0


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

       ipmin = 1
       ipmax = N




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
       end  subroutine NrouseP0




 end function th_locdgref

 
