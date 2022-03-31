      FUNCTION th_roufric (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rouse <--------                                              
!                                                                       
!                                    
      implicit none                                   
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      real :: pa, qq, x, th_roufric, zpi
      integer :: npar, ini, nparx, ier, iert, ierq
      integer                     :: mbuf
      integer, intent(inout)      :: nopar                 ! Anzahl der Parameter data
      character(len=80), intent(inout) :: napar(mbuf)           ! name des parameters n
      real, intent(inout)         :: params(mbuf)          ! value des parameters n
      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) a0, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff, wl4, re_arm 
      REAL qget, tget, tauget 

      integer :: f_arm, n_arm
      double precision :: fric0, fricn, nfric, fricc0
!      integer :: onearm, fixcenter
      double precision :: slope,  sqt0
       
      logical :: x_is_tau = .true.                                                                 
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'roufric' 
         nparx = 8
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_roufric = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit'     ! general prefactor 
         parnam (2) = 'wl4'          ! Rouse Rate W l**4 in A**4/ns
         parnam (3) = 're_arm'       ! Re of one arm     in A
         parnam (4) = 'n_arm '       ! number of beads per arm
         parnam (5) = 'fric0'        ! friction modulation fmod= fric0+fricn*i**nfric
         parnam (6) = 'fricn'        !      " 
         parnam (7) = 'nfric'        !      " 
         parnam (8) = 'fricc0'       ! friction modifer center
                                                                        
!                                                                       
         th_roufric = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau   = x 
      a0    = pa (1) 
      wl4   = abs (pa (2) ) 
      re_arm    = abs (pa (3) ) 
      n_arm     = nint (pa (4) ) 
      fric0     = pa(5)
      fricn     = pa(6)
      nfric     = pa(7)
      fricc0    = pa(8)

                              ! in cm**2/sec                            
      diff  = 0d0 

      f_arm=1
                                                                          
      qget = 0.01 
      CALL getpar ('q       ', qget,nopar ,params,napar,mbuf, ierq)  
      if(ierq == 0) then 
         x_is_tau =.true.
         qz       = qget
         tau      = x 
      else
         tauget = 0
         x_is_tau = .false.
         CALL getpar ('tau     ', tauget,nopar ,params,napar,mbuf, iert)  
         tau = tauget
         qz  = x
      endif      
      IF (ierq.ne.0) write (6, * ) 'Warning q not found' 

      IF (temp.eq.0.0d0) then 
         tget = 300.0 
      CALL getpar ('temp    ', tget,nopar ,params,napar,mbuf, ier) 
         temp = tget 
      ENDIF 
        
      if(x_is_tau) then                                                                
        sqt0   =       pericoroufric_sqt(qz,0d0,f_arm,n_arm,Re_arm,Wl4, fric0, fricn, nfric)                                                       
        th_roufric=  a0 * pericoroufric_sqt(qz,tau,f_arm,n_arm,Re_arm,Wl4,fric0,fricn, nfric)/sqt0
      else
        th_roufric=  a0 * pericoroufric_sqt(qz,tau,f_arm,n_arm,Re_arm,Wl4, fric0, fricn, nfric)
      endif 

      if(diff .ne. 0d0) then
        th_roufric = th_roufric * exp(-qz*qz*diff*tau)
      endif 
                                                                       
     


CONTAINS
   
double precision function pericoroufric_sqt(q,t,f_arm,n_arm,Re_arm,Wl4, fric0, fricn, nfric)
!----------------------------------------------------------------------------------------------
!! star following the derivation of:
!!      M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!
!! this module written by: michael Monkenbusch, Jcns
!!
  implicit none

  double precision, intent(in)  :: q
  double precision, intent(in)  :: t
  integer,          intent(in)  :: f_arm
  integer,          intent(in)  :: n_arm 
  double precision, intent(in)  :: Re_arm
  double precision, intent(in)  :: Wl4
  double precision, intent(in)  :: fric0, fricn, nfric
  double precision, parameter   :: Pi=3.141592654d0
  integer                       :: N
  integer                       :: i, j, ip, l, i1, i2, job, info
  double precision              :: l0         ! effective segment length
  double precision              :: tauR       ! basic Rouse-time of a chain of 2 arm lengths

  double precision              :: p, s11, s12, dec


 logical, save                :: first_run = .true.
 integer, save                :: f_arm0, n_arm0
 double precision, save       :: fric00, fricn0, nfric0, fricc00 
  

 real                         :: U(1+n_arm*f_arm,1+n_arm*f_arm)
 real, allocatable,save       :: H(:,:)
 real                         :: X(1+n_arm*f_arm,1+n_arm*f_arm)
 real                         :: M(1+n_arm*f_arm,1+n_arm*f_arm)
 real, allocatable,save       :: A(:,:)  
 real                         :: HA(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: Ah(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: phi(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: work(2*(1+n_arm*f_arm))
 complex                      :: csum
 complex, allocatable, save   :: Eigenvalues(:)  
 complex, allocatable, save   :: Eigenvectors(:,:) 
 complex                      :: CX(1+n_arm*f_arm,1+n_arm*f_arm)
 complex, allocatable, save   :: CM(:,:)


 ! real                         :: std_diffusion
 ! real                         :: std_bead_friction
!  real                         :: effective_diffusion



  N      = f_arm * n_arm + 1
  l0     = Re_arm/sqrt(dble(n_arm))

!! compute Eigenvectors and Eigenvalues for the problem:
!! see: M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!
!! modification with increased friction after each n_repet beads
!! 


    N=1+n_arm*f_arm
   


eig:  if(  first_run                    &
     .or.  n_arm    .ne. n_arm0         &
     .or.  abs(fric0-fric00)   > 1e-8 &
     .or.  abs(fricc0-fricc00) > 1e-8 &
     .or.  abs(fricn-fricn0)   > 1e-8 &
     .or.  abs(nfric-nfric0)   > 1d-8 ) then
    
    first_run    = .false.
    n_arm0       = n_arm
    fric00       = fric0
    fricc00      = fricc0
    fricn0       = fricn
    nfric0       = nfric
   

    if(allocated(Eigenvalues)) then
      deallocate(A)
      deallocate(H)
      deallocate(Eigenvalues)
      deallocate(Eigenvectors)
      deallocate(CM)
    endif

    allocate(Eigenvalues(N))
    allocate(A(N,N))
    allocate(H(N,N))
    allocate(Eigenvectors(N,N))
    allocate(CM(N,N))

    U = 0            ! inverse of static bond correlation matrix <l_i*l_j>/l**2  (here  1)
    M = 0            ! bead-to-bond vector transformation matrix
    A = 0            ! A = Mt(U)M   (were the first col and row of U are 0)
    H = 0            ! H = hydrodynamic matrix (so far here 1) friction is lumped in the Wl4 parameter
   

 forall (i=2:N)   U(i,i) = 1.0
!? forall (i=1:N)   M(1,i) = 1.0/N
M(1,1)=1 !?
!?<<
 forall (i=2:N)   M(i,i) = 1.0
 forall (i=1:N)   H(i,i) = 1.0



 H(1,1) = H(1,1)/fricc0
! do j=1,f_arm
   do i=1,n_arm
!      i1 = 1+n_arm*(j-1)+i
      i1 = 1+i
      H(i1,i1) = H(i1,i1)/(fric0+fricn*(dble(i)/n_arm)**nfric)    !!! hier ??? !!!
   enddo
! enddo 

! check
! write(6,'(i3,f12.6)')(i,H(i,i),i=1,n)
! check
 
 do j=1,f_arm
   i1 = 2+n_arm*(j-1)
   i2 = j*n_arm
   forall (i=i1:i2)   M(i+1,i) = -1
 enddo
 do j=1,f_arm
   i = 2+(j-1)*n_arm
   M(i,1) = -1
 enddo

!! testing >>>>
!! write(*,*)"M:"
!! do i=1,N
!!   write(*,'(100i2)')nint(M(1:N,i))
!! enddo
!! testing <<<<
!
! compute A = M^t (U) M
!

do i=1,N
 do j=1,N
  X(i,j) = dot_product(U(i,1:N),M(1:N,j))
 enddo
enddo
do i=1,N
 do j=1,N
  A(i,j) = dot_product(M(1:N,i),X(1:N,j))
 enddo
enddo

! and now HA 
do i=1,N
 do j=1,N
  HA(i,j) = dot_product(H(i,1:N),A(1:N,j))
 enddo
enddo


 Ah = HA
 job = 1

 call SGEEV (Ah, N, N, Eigenvalues, Eigenvectors, N, WORK, JOB, INFO)

! C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
! C            of a real general matrix.
! C***LIBRARY   SLATEC
! C***CATEGORY  D4A2
! C***TYPE      SINGLE PRECISION (SGEEV-S, CGEEV-C)
! C***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
! C***AUTHOR  Kahaner, D. K., (NBS)
! C           Moler, C. B., (U. of New Mexico)
! C           Stewart, G. W., (U. of Maryland)

!! computation of mua !!
do i=1,N
 do j=1,N
  CX(i,j) = dot_product(A(i,1:N),Eigenvectors(1:N,j))
 enddo
enddo
do i=1,N
 do j=1,N
  CM(i,j) = dot_product((Eigenvectors(1:N,i)),CX(1:N,j))
 enddo
enddo
   

! ATTENTION: eigenvalues are not sorted !!!
 write(6,*) &
' #i     Eigenvalues(Re Im)                CM(ll)(Re Im)          |EV(1,i)|     |EV(i,1)|      H(i,i): '
  do i=1,N
    write(6,'(i5,2f12.6,5x,2f12.6,4x,f12.6,2x,f12.6,3x,f12.6)') i,eigenvalues(i), CM(i,i), &
                                    abs(Eigenvectors(1,i)), abs(Eigenvectors(i,1)), H(i,i)
  enddo
 write(*,*)"Number of near zeroes(1e-3) EV(1): ",count(abs(Eigenvectors(1,1:N))<1d-3/N)
 write(*,*)"Number of nera zeroes(1e-6) EV(1): ",count(abs(Eigenvectors(1,1:N))<1d-6/N)
!! OMIT (below) ALL ZEROES IF CENTER SHALL BE FIXED

 endif eig



! -- fill the (symmetric) phi(n,m) matrix --
phi = 0
!?! philp2: do i=1,N
!?!   do j=1,N
!?!    csum = 0
!?!    do l=1,N
!?!     if(fixcenter==1 .and. abs(Eigenvectors(1,l))>1d-3/N) cycle
!?!     if(fixcenter==1 .and. abs(Eigenvectors(1,j))>1d-3/N) cycle  !!??
!?!     if(abs(CM(l,l)) > 0.00001) then
!?!     csum = csum + &
!?!       (Eigenvectors(i,l))**2+(Eigenvectors(j,l))**2
!?! !    if(fixcenter==1 .and. abs(Eigenvectors(1,l))>1d-3/N) cycle
!?! !!    if(fixcenter==1 .and. abs(Eigenvectors(1,j))>1d-3/N) cycle 
!?!     csum = csum - &
!?!     2*(Eigenvectors(i,l))*Eigenvectors(j,l)*exp(-wl4/(l0**4)*Eigenvalues(l)*t)/CM(l,l)
!?!     endif
!?!    enddo
!?!    phi(i,j) = csum * l0**2
!?!   enddo
!?!  enddo philp2


philp2: do i=1,N
  do j=1,N
   csum = 0
   do l=1,N
    if(abs(CM(l,l)) > 0.00001) then
    csum = csum + &
      ((Eigenvectors(i,l))**2+(Eigenvectors(j,l))**2- &
    2*(Eigenvectors(i,l))*Eigenvectors(j,l)*exp(-wl4/(l0**4)*Eigenvalues(l)*t))/CM(l,l)
    endif
   enddo
   phi(i,j) = csum * l0**2
  enddo
 enddo philp2


!
! --> ok phi matrix filled 
!
  
!
! -- and now the scattering contibutions of the arms
!
        s11 = 0         ! the "self term"
!$OMP PARALLEL DO REDUCTION(+:S11)      
        do i=2,n_arm+1
!          if(fixcenter==1 .and. abs(Eigenvectors(1,i))>1d-3/N) cycle
           s11 = s11 + 2 * exp(-(q*q/6d0) * phi(1,i)) + exp(-(q*q/6d0)* phi(i,i)) 
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO REDUCTION(+:S11)  
        do i=2,n_arm
!         if(fixcenter==1 .and. abs(Eigenvectors(1,i))>1d-3/N) cycle
         do j=i+1,n_arm+1
!           if(fixcenter==1 .and. abs(Eigenvectors(1,j))>1d-3/N) cycle
           s11 = s11 + 2* exp(-(q*q/6d0)* phi(i,j))
         enddo
        enddo
!$OMP END PARALLEL DO


        pericoroufric_sqt = (1+f_arm*(s11+(f_arm-1)*s12))/N**2

! write(6,*) s11, s12, sum, sqrt(a*a*n_arm)




! tentatively add diffusion !
!     if(diff == 0d0) then
!        std_diffusion       =  wl4/(3*f_arm*n_arm*l0*l0)
!        std_bead_friction   =  1/(n*std_diffusion)
!        effective_diffusion =  1/(sum([(1.0/H(i,i),i=1,n)])*std_bead_friction)
!        pericoroufric_sqt = pericoroufric_sqt * exp( -q*q* effective_diffusion * t)
!     endif
!


end function pericoroufric_sqt

END FUNCTION th_roufric

