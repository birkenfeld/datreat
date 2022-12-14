program testvilgi
 implicit none
 double precision :: l     = 20d0
 double precision :: wl4   = 100000d0
 double precision :: rmesh = 500
 integer          :: N     = 100

 double precision :: q, t, sqt, sqt0, sqti

 integer          :: m=10
 integer          :: i
 q = 0.05d0
 t = 0.1d0
 
 call NNvilgis(q,0d0,wl4,l,rmesh,N,sqt0)

 do i=0,m
   t = 2*i
   call NNvilgis(q,t,wl4,l,rmesh,N,sqt,sqti)
   write(*,*)q,t,sqt0,sqt/sqt0, sqti/sqt0
 enddo

CONTAINS


 subroutine NNvilgis(q,t,wl4,l,rmesh,N,sqt, sqtint)
!==================================================
!
! Implementation of Eq. 30 of 
! [Vilgis and Boue, J. Polymer Science Part B: Polymer Physics, Vol 26, 2291-2301 (1988)]
! the Fourier integral can be done analytically an yields terms as b*exp(-a*t/(p**2+p0**2)/(p**2+p0**2)
! the p-sum is then done explicitly (efficiency increse by observing that all factors can be written as cos*
! NOTE: there are several typos in the formulae of the paper !! 
! 
  implicit none

  double precision, intent(in)  :: q        ! momentum transfer
  double precision, intent(in)  :: t        ! time
  double precision, intent(in)  :: wl4      ! rouse rate
  double precision, intent(in)  :: l        ! segment length
  double precision, intent(in)  :: rmesh    ! "confinement" length due to the harmonic potential
  integer         , intent(in)  :: N        ! number of segments (summation)
  double precision, intent(out) :: sqt      ! the sqt-value
  double precision, intent(out), optional :: sqtint      ! the sqt-value internal modes only


  double precision, parameter :: Pi = 4d0*atan(1d0)

  double precision :: coser(1:N,1:N)
  double precision :: eser(1:N)
  double precision :: rrarr(N,N)  ! just for checking can be removed later

  double precision :: p0, q0, qq0, W, pidn, rr, tt

  integer          :: nn, mm, p


  W    = wl4/l**4
  pidn = Pi / N
  q0   = l**2/(3*rmesh**2)
  p0   = q0 * N / Pi
!$OMP PARALLEL DO 
  do p=1,N
    do nn=1,N
      coser(nn,p) = cos(Pi*p*nn/dfloat(N))  / sqrt((p**2 + p0**2))
    enddo
  enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO 
  do p=1,N
      eser(p) = exp(-(Pi**2/N**2)*W*(p**2 + p0**2)*t)
! write(*,*)"eser:",p,t,W,(Pi**2/N**2)*W*(p**2 + p0**2),eser(p) 
  enddo
!$OMP END PARALLEL DO 

  sqt = 0
!$OMP PARALLEL DO PRIVATE(rr) REDUCTION(+:sqt)
  do nn=1,N
    do mm=1,N
      rr = 0
!$OMP PARALLEL DO REDUCTION(+:rr)
      do p=1,N
        rr = rr + coser(nn,p)**2 + coser(mm,p)**2 - 2*coser(nn,p)*coser(mm,p) * eser(p)
!        rr = rr + cos(Pi*p*nn/dfloat(N))**2/p**2 + cos(Pi*p*mm/dfloat(N))**2/p**2 - &
!                2*cos(Pi*p*nn/dfloat(N))*cos(Pi*p*mm/dfloat(N))/p**2
      enddo
!$OMP END PARALLEL DO 
      rr  = rr * 4*N / (2*Pi) /Pi 
!                              === dieser Faktor is noetig damit rr(n,m,0)=abs(n-m) 
!                              === fehlt im Paper ?!
! write(*,*)nn,mm,rr
      sqt = sqt + exp(-(q*l)**2/6d0 * rr)
      rrarr(mm,nn) = rr
    enddo
  enddo
!$OMP END PARALLEL DO 
  sqt = sqt / (N*N) 
  if(present(sqtint)) sqtint = sqt
 
! Diffusion
  qq0 = q0**2
  if(abs(W*qq0*t) > 1d-6) then
     rr     = 2*l**2/N/qq0*(1-exp(-W*qq0*t))  ! a factor of 2 different from Vilgis ...
  else
     rr     = 2*l**2*W*t/N-l**2*W**2*t**2*qq0/N
  endif

  sqt    = sqt*exp(-(q)**2/6d0 * rr)

end subroutine NNvilgis

end program testvilgi 

