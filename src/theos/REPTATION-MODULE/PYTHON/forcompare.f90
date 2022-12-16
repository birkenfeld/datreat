program forcompare
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
!  gfortran -c reptation_module -fopenmp -O2
!
!==================================================================================


    double precision :: q     ! the Q-value 
    double precision :: t     ! the time
    double precision :: N     ! total number of segments in chain
    double precision :: lseg  ! segment length
    double precision :: Ne    ! average number of segments per entanglement strand
    double precision :: Re    ! entanglement "blob" size (R-end-end)
    double precision :: Wl4   ! Rouse rate
    double precision :: alpha0     ! Non-Gaussianity correction amplitude
    double precision :: talphamax  !     "  distribution max-location
    double precision :: talphawd   !     "        "      width

    double precision              :: nrouse_ngcorr

    double precision, parameter   :: teps = 1d-6
    double precision :: W, Wr, lr
    double precision :: tmax, twidth
    double precision :: sq_rouse
    integer          :: Nrsum, i

    double precision :: b=1

    wl4    = 30000
    Re     = 50d0
    Nrsum  = 100  
    alpha0 = 0.10
    tmax   = 1
    twidth = 1
    q      = 0.15d0

    do i=0,50
      t = 1.1d0**i
      sq_rouse = nrouse_ngcorr( q, t, Nrsum, Re, Wl4, alpha0, tmax, twidth) 
      write(*,*) t, sq_rouse
    enddo

 contains

  end program forcompare



       function nrouse_ngcorr(q,t,N,R,Wl4,alpha0,talphamax,talphawd) result(sqt)
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
       

       
 contains
 

   function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
      double precision, intent(in) :: t
      double precision             :: a
   
   ! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
   ! the model
      a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) 
      
   end function alpha

end function nrouse_ngcorr


