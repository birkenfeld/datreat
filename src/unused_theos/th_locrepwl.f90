 FUNCTION th_locrepwl(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrepwl
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

     double precision :: ampdebye    ! prefactor                                                                       
     double precision :: wl4        ! rouse rate                                                                      
     integer          :: n_segmen   ! effective N for summation                                                       
     double precision :: re         ! end-to-end radius of Gaussian coil                                              
     double precision :: temp       ! temperature if not in parameters                                                
     double precision :: com_diff   ! center of mass diffusion (if 0, use Rouse default) 
     double precision :: pmin       ! minimum mode to be included 
     double precision :: pmax       ! maximum mode to be included  
     double precision :: t0_locrep   ! locrep cutoff
     double precision :: tw_locrep    ! cutoff width
                                                              
! the reout parameter representation 
     double precision :: Dr         ! diffusion   
     double precision :: W0, wr      ! rouse     W
     double precision :: l          ! sqgment length    

     double precision :: ne, n
                                                                    
!     double precision :: sqdeb0, sqdebt



     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset


!
 
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: local_reptation, sqt, sqt0
     double precision :: sqdebt, sqdeb0
     double precision :: slope_l , tauexp, alphaex
     double precision :: ww0, t0, tww
     double precision :: dnse, deff, dblob
     double precision :: f_re, f_dnse, f_deff, wd_mod
     double precision :: sumlik, sumlik0, mue, st, alpha_lik
     

     double precision, parameter   :: teps = 1d-6
!     integer                       :: mode
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locrepwl'
       nparx =        12
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrepwl = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " fermi W modificatio and b,w,ne coupling," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters  locrep t=0 limit exp(-t) "//cr//parspace//&
                                " correlation in timescales by factors to Re, dNSE, deff "
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'b       '  ! fluctuation intensity 
        parnam ( 2) = 'a       '  ! length scale                                                                    
        parnam ( 3) = 'ne      '  ! entanglemen ne length        
        parnam ( 4) = 'n       '  ! n                                                                      
        parnam ( 5) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 6) = 're      '  ! end-to-end radius of Gaussian coil                                              
        parnam ( 7) = 'f_re    '  ! cutoff of locrep
        parnam ( 8) = 'f_dnse  '  ! width of fermi    
        parnam ( 9) = 'f_deff  '  !                                                     
        parnam (10) = 'wd_mod  '  !                                                     
        parnam (11) = 'alpha0  '  !                                                     
        parnam (12) = 'alpha_l '  !                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "fluctuation intensity" !//cr//parspace//&
        th_param_desc( 2,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "totla n" !//cr//parspace//&
        th_param_desc( 5,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 6,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 7,idesc) = "dblob factor for Re" !//cr//parspace//&
        th_param_desc( 8,idesc) = "dblob factor for dnse" !//cr//parspace//&
        th_param_desc( 9,idesc) = "dblob factor for deff" !//cr//parspace//&
        th_param_desc(10,idesc) = "width modification factro w(t), alpha(t)" !//cr//parspace//&
        th_param_desc(11,idesc) = "NG alpha strength" !//cr//parspace//&
        th_param_desc(12,idesc) = "CLF Likthmann alpha (0=off)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(1,idesc)   = "a2sqt = a**2/sqrt(tau) "
!
                                             
! >>>>> describe parameters >>>>>>> 
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "l        > effective segment length (Rouse summation)"
        th_out_param(  2,idesc) = "w0       > rouse rate: W for WL4/a**4"
        th_out_param(  3,idesc) = "wl4      > inferred value of Wl4"
        th_out_param(  4,idesc) = "deff     > sqrt((dnse**2+2*re**2)/3)"
        th_out_param(  5,idesc) = "dnse     > sqrt(b*ne) * a"
        th_out_param(  6,idesc) = "dblob    > f_re*re+f_dnse*dnse+f_deff*deff"
        th_out_param(  7,idesc) = "t0_locrep> t0 of w(t) (=t_blob)"
        th_out_param(  8,idesc) = "tw_locrep> width of w(t) fermi modification "
!


 
        th_locrepwl = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      b        =   abs( pa( 1) )
      a        =   abs( pa( 2) )
      ne       =   abs( pa( 3) )
      n        =   abs( pa( 4) )
      wl4      =   abs( pa( 5))
      re       =   abs( pa( 6))
      f_re     =   abs( pa( 7))
      f_dnse   =   abs( pa( 8))
      f_deff   =   abs( pa( 9))
      wd_mod   =   abs( pa(10))
      alpha0   =      ( pa(11))
      alpha_lik=      ( pa(12))
 

      ampli     = 1d0
      ampdebye  = 1d0

     
      W0   = wl4/a**4
      ww0 = 1d0/ne

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      xh =     400d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh

      dnse = sqrt(b*ne)*a
      deff = sqrt( (dnse**2 + 2*re**2)/3d0 )

      dblob = f_re * re + f_dnse * dnse + f_deff * deff

      t0_locrep = dblob**4 /(Pi**2 * wl4)  
      tw_locrep = t0_locrep * wd_mod

      talphamax = t0_locrep
      talphawd  = log( talphamax * wd_mod )

      n_segmen = ne
      com_diff = 1d-14
      pmin     = 1
      pmax     = n_segmen
      aoffset  = 0
      temp     = 400d0     ! ist irrelevant



      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 



     Dr  = com_diff * 1d-9 / 1d-16  ! in A**2/ns
!    call NrouseP        (q,t,temp,Dr,wl4,n_segmen,Re, Wr, l,pmin,pmax, Sqdeb0,Sqdebt)
     call nrousalpha_corr(q,t,temp,Dr,wl4,N_segmen,Re, Wr, l,pmin,pmax, Sqdeb0,Sqdebt) 

!!     t0  = t0_locrep * N_segmen**2/Wr/Pi**2
     t0  = t0_locrep
     tww = tw_locrep 


     sumlik = 1
     if(alpha_lik > 0d0) then 
        mue=q**2*n*a**2/12.d0   
!        st = 0.5d0*((1.5d0/num)*(t/taue)**0.25)
        st = 0.5d0*((ne*alpha_lik/n)*(t/t0)**0.25)

        sumlik = (n/(2*mue**2))*(2*mue+dexp(-2*mue)+2.d0-4*mue*st-4* &
            dexp(-2*mue*st)+dexp(-4*mue*st))


!write(*,*)"lik:",q,n,q,alpha_lik,q,t,t0,st,sumlik
!write(*,*)"lik:",q,t,t0,mue,st,sumlik
        sumlik0 = (n/(2*mue**2))*(2*mue+exp(-2*mue)-1.d0)
        sumlik  = sumlik/sumlik0
!write(*,*)"TST:",sumlik,sumlik0
     endif

!?     sqt0  =  local_reptationdr2( q, teps,      a, W0, n, ne,  b, 0)
     sqt0  = 1d0
     sqt   =  local_reptationdr2( q, t,         a, W0, n, ne,  b, 0)

   

!     Sqdeb0 = ampdebye * Sqdeb0 
!     Sqdebt = ampdebye * Sqdebt 


 ! write(*,*)"q,t=",q,t,"  sqt0, sqt=",sqt0,sqt, "debye=", Sqdeb0,Sqdebt
      th_locrepwl = ampli * sqt/sqt0  * Sqdebt/Sqdeb0 

      th_locrepwl = th_locrepwl * sumlik

!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda)
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda)      ! in ns A units
      call        parset('w0      ',sngl(W0),iadda)      !     "
      call        parset('wl4     ',sngl(W0*l**4),iadda) !     "
      call        parset('deff    ',sngl(deff),iadda) !     "
      call        parset('dnse    ',sngl(dnse),iadda) !     "
      call        parset('dblob   ',sngl(dblob),iadda) !     "
      call        parset('t0_locrep     ',sngl(t0),iadda) !     "
      call        parset('tw_locrep     ',sngl(tww),iadda) !     "
  
!      Dr        = Dr /( 1d-9 / 1d-16 ) ! in cm**2/s
!      call parset('diff    ',sngl(Dr),iadda,ier)
 

                   


CONTAINS


 
  function local_reptationdr2(q, t, a, W0, n, ne, b, modelr) result(val)   !! with series expansion
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: a        ! step (entanglement?) length
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

    W = w0*(1d0 - (1d0-ww0)/(1+exp((t-t0)/tww)))

    T1 = 72d0 * (exp(-q**2 * n * a**2 / 6d0) + q**2 * n * a**2 / 6d0 - 1d0) / q**4 / a**4
!    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
!        ne*(n/3+a**2*q**2*w*t/9)*exp(a**4*q**4*w*t/36)* &
!         (erfc(a**2*q**2*sqrt(w*t)/6)-erfc(n/(2*sqrt(w*t))+a**2*q**2*sqrt(w*t)/6)))
    T3 = 72d0/((q**4*a**4))*(exp(-(q**2*n*a**2)/6d0)+(q**2*n*a**2)/6d0-1d0)+b*ne*n/3d0
!
!    val = (T1+T2)/T3

    x = sqrt(a**4*q**4*w*t/36)
    xa= n/(2*sqrt(w*t))
    
    if( abs(x) < xlimser ) then
      dec = exp(x**2) *(erfc(x)-erfc(x+xa))
    else
      dec =  1d0/(sqrt(Pi)*x)-1d0/(2*sqrt(Pi)*x**3)+3d0/(4*sqrt(Pi)*x**5) &
          - (x**4-xa*x**3+(xa**2-1d0/2d0)*x**2+(-xa**3+3d0/2d0*xa)*x+xa**4-3*xa**2+3d0/4d0) &
           * exp(-xa**2)*exp(-2*xa*x)/(sqrt(Pi)*x**5)
    endif


    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
        ne*(n/3+a**2*q**2*w*t/9) &
        * dec )

    if(modelr == 1) then
      val = T2
    else
      val = (T1+T2)/T3
    endif


  end function local_reptationdr2





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
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*i)/dfloat(N)))*t)) &
         * (((erf(-pmin+i)-1d0)/2+1) * ((erf(pmax-i)-1d0)/2+1)) 
       enddo
!$OMP END PARALLEL DO    


       ipmin = max(1,nint(pmin)-4)
       ipmax = min(N,nint(pmax)+4)

!       ipmin = 1
!       ipmax = N




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
       end  subroutine NrouseP




 function fermi(t, t0, w) result(f)
  double precision, intent(in) :: t
  double precision, intent(in) :: t0
  double precision, intent(in) :: w
  double precision             :: f

  f = 1d0/(1d0+exp((t-t0)/w))

 end function fermi



 

       subroutine nrousalpha_corr(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
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
       fqq0 = 1d0                                   !! since rmm(0)=0


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






 end function th_locrepwl

 
