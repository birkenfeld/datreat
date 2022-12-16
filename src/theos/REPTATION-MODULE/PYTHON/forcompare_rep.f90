program forcompare
use reptation
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


    double precision, parameter   :: teps = 1d-6
    double precision :: W, Wr, lr
    double precision :: tmax, twidth
    double precision :: sq_rouse(2), sqrep(2)
    integer          :: Nrsum, i

    double precision :: b=1

    wl4    = 30000
    Re     = 50d0
    Nrsum  = 100 
    n      = 5000
    ne     = 100
    lseg   = 5d0 
    alpha0 = 0.30
    tmax   = 1
    twidth = 1
    q      = 0.15d0

    w = wl4/lseg**4

    do i=0,50
      t = 0.01 * 1.3d0**i
      sq_rouse = nrouse_ngcorr( q, t, Nrsum, Re, Wl4, alpha0, tmax, twidth) 
      lr = local_reptationdr( q, t, lseg, w,  n ,ne )
      sqrep= reptation_sqt(q,t, N, lseg, Ne, Re, wl4, alpha0, tmax,twidth)
      write(*,*) t, sqrep(2)/sqrep(1),lr, sq_rouse
    enddo

 contains

end program forcompare

