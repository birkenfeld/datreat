MODULE reptation2
implicit none
!=================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
!  and modified to account for the Gaussian chain condition (Debye)
!==================================================================================
!
! Author: Michael Monkenbusch, JCNS-1, Forschungszentrum Juelich
!         m.monkenbusch@fz-juelich.de
!
!==================================================================================
!  GNU GENERAL PUBLIC LICENSE 
!  Version 3, 29 June 2007
!==================================================================================
!
!  complile with 
!  gfortran -c reptation2_module -fopenmp -O2
!
!==================================================================================


  double precision :: b_fluctuation = 1
  double precision :: b_z_floor     = 1
  double precision :: b_mcdefault   = sqrt(2d0)-1d0

CONTAINS

  function reptation2_sqt(q,t, N, lseg, Ne, Re, wl4, alpha0, talphamax,talphawd) result(sqr0t)
    implicit none
    double precision, intent(in)  :: q     ! the Q-value 
    double precision, intent(in)  :: t     ! the time
    double precision, intent(in)  :: N     ! total number of segments in chain
    double precision, intent(in)  :: lseg  ! segment length
    double precision, intent(in)  :: Ne    ! average number of segments per entanglement strand
    double precision, intent(in)  :: Re    ! entanglement "blob" size (R-end-end)
    double precision, intent(in)  :: Wl4   ! Rouse rate
    double precision, intent(in)  :: alpha0     ! Non-Gaussianity correction amplitude
    double precision, intent(in)  :: talphamax  !     "  distribution max-location
    double precision, intent(in)  :: talphawd   !     "        "      width

    double precision              :: sqr0t(2)

    double precision, parameter   :: teps = 1d-6
    double precision :: W, Wr, lr
    double precision :: tmax, twidth
    double precision :: sq_rouse(2), sq_locrep(2)
    integer          :: Nrsum

    double precision :: b=1

    if(alpha0 == 0d0) then
       tmax   = 1d0
       twidth = 1d0
    else
       tmax   = talphamax
       twidth = talphawd
    endif
    
    W         = Wl4 / lseg**4
    Nrsum     = min(300,nint(Ne))

    sq_rouse = nrouse_ngcorr( q, t, Nrsum, Re, Wl4, alpha0, tmax, twidth) 

    sq_locrep(1)   =  local_reptationdr( q, teps          ,  lseg, W, n, ne)
    sq_locrep(2)   =  local_reptationdr( q, max(t,teps)   ,  lseg, W, n, ne)

    sqr0t = sq_rouse * sq_locrep

  end function reptation2_sqt






 
  function local_reptationdr(q, t, lseg, W, n, ne) result(val)   
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: lseg     ! lsegment
    double precision, intent(in)   :: W        ! Rouse rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments per entanglement strand
  
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)
    double precision, parameter    :: xlimser = 15d0  ! where to switch erf expression to series expansion

    double precision :: T1, T2, dec
    double precision :: x, xa
    double precision :: Z

   

    Z    = n/ne

    T1 = 72d0 * (exp(-q**2 * n * lseg**2 / 6d0) + q**2 * n * lseg**2 / 6d0 - 1d0) / q**4 / lseg**4
!!>> limit(T1(q=0)) is  n**2  
!!>> here we tweak the normalisation to the new interpretation
    T1 = T1 / n**2   !! Normalisation to 1 at q=0 
!! now T2 the local reptation part
    x  = sqrt(lseg**4*q**4*w*t/36)
    xa = n/(2*sqrt(w*t))
    
    if( abs(x) < xlimser ) then
      dec = exp(x**2) *(erfc(x)-erfc(x+xa))
    else
      dec =  1d0/(sqrt(Pi)*x)-1d0/(2*sqrt(Pi)*x**3)+3d0/(4*sqrt(Pi)*x**5) &
          - (x**4-xa*x**3+(xa**2-1d0/2d0)*x**2+(-xa**3+3d0/2d0*xa)*x+xa**4-3*xa**2+3d0/4d0) &
           * exp(-xa**2)*exp(-2*xa*x)/(sqrt(Pi)*x**5)
    endif

    T2  = ((exp(-1/6*lseg**2*n*q**2-n**2/(4*w*t))-1)*2/3*ne*sqrt(w*t/Pi)+ &
        ne*(n/3+lseg**2*q**2*w*t/9) &
        * dec )

    T2  = T2*3d0 / (n*ne)    !! Normalisation to T2(t=0) = 1

!    val = 1d0/(1+1/Z) *  (T1 + (2d0/3d0 + (1d0/3d0) * T2) / Z )
    val = 1d0/(1+1/Z*b_z_floor) *  (T1 + (1d0 + b_mcdefault*b_fluctuation * T2) / Z * b_z_floor )

  end function local_reptationdr



       function nrouse_ngcorr(q,t,N,R,Wl4,alpha0,talphamax,talphawd) result(sq0t)
!      =========================================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q             ----> momentum transfer in A**-1
!    t             ----> time in nano-sec
!    N             ----> effective number of segments
!    R             ----> Re of "blob", l will be computed from N and R
!    W             ----> Rouse rate
!    alpha0        ----> strength of Non-Gaussianity (NG) correction
!    talphamax     ----> peak of the NG-alpha function
!    talphawd      ----> width of the (log) distribution alpha(t)
!    sqt(1:2)      <--- [S(Q), S(Q,t)]
! ------------------------------------------------------------------------------
!
       implicit none
       double precision, intent(in)  :: q
       double precision, intent(in)  :: t
       integer         , intent(in)  :: N
       double precision, intent(in)  :: R
       double precision, intent(in)  :: Wl4
       double precision, intent(in)  :: alpha0
       double precision, intent(in)  :: talphamax
       double precision, intent(in)  :: talphawd

       double precision              :: sq0t(2)


       double precision, parameter    :: Pi = 4*atan(1d0)
     
       integer :: nn,mm,ip

       double precision :: Sq, Sqt
       double precision :: Dr, l, ff2, W

       double precision :: cosarray(N,N), ewfac(N)
       double precision :: rmm, fqq, fqq0

       integer :: ipmin, ipmax, i

!       integer iout
       
       if(N.le.0) then
         sq0t = 0
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       Dr  = 0d0                ! in this module: no com-diffusion of the "blobs"
       W = Wl4/l**4

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

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N
       Sq0t = [Sq, Sqt]

       
 contains
 

   function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
      double precision, intent(in) :: t
      double precision             :: a
   
   ! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
   ! the model
      a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) 
      
   end function alpha

end function nrouse_ngcorr

END MODULE reptation2

