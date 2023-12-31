      FUNCTION th_star (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rouse <--------                                              
!                                                                       
!                                    
      implicit none                                   
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      real :: pa, qq, x, th_star, zpi
      integer :: npar, ini, nparx, ier, iert, ierq
      integer                     :: mbuf
      integer, intent(inout)      :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)           ! name des parameters n
      real, intent(inout)         :: params(mbuf)          ! value des parameters n
      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff, wl4, re_arm 
      REAL qget, tget, tauget 

      integer :: f_arm, n_arm
      double precision :: pericostar_sqt, sqt0
       
      logical :: x_is_tau = .true.                                                                 
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'star' 
         nparx = 6 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_star = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit' 
         parnam (2) = 'wl4' 
         parnam (3) = 're_arm' 
         parnam (4) = 'n_arm ' 
         parnam (5) = 'f_arm' 
         parnam (6) = 'com_diff' 
                                                                        
!                                                                       
         th_star = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau   = x 
      a0    = pa (1) 
      wl4   = abs (pa (2) ) 
      re_arm= abs (pa (3) ) 
      n_arm = nint (pa (4) ) 
      f_arm = nint (pa (5) )
                              ! in cm**2/sec                            
      diff = abs (pa (6) ) 
      diff = diff * 1d-9 / 1d-16 
                                                                        
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
        sqt0   =       pericostar_sqt(qz,0d0,f_arm,n_arm,Re_arm,Wl4,diff)                                                       
        th_star=  a0 * pericostar_sqt(qz,tau,f_arm,n_arm,Re_arm,Wl4,diff)/sqt0
      else
        th_star=  a0 * pericostar_sqt(qz,tau,f_arm,n_arm,Re_arm,Wl4,diff)
      endif 

      if(diff .ne. 0d0) then
        th_star = th_star * exp(-qz*qz*diff*tau)
      endif 
                                                                       
      RETURN 
      END FUNCTION th_star



   
double precision function pericostar_sqt(q,t,f_arm,n_arm,Re_arm,Wl4,diff)
!-------------------------------------------------------------------
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
  double precision, intent(in)  :: diff


  double precision, parameter   :: Pi=3.141592654d0
  integer                       :: N
  integer                       :: i, j, ip, l, i1, i2, job, info
  double precision              :: l0         ! effective segment length
  double precision              :: tauR       ! basic Rouse-time of a chain of 2 arm lengths

  double precision              :: p, sum, s11, s12, dec

  integer, parameter   :: real_kind = 4
 
  real(kind=real_kind)         :: U(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: H(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: X(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: M(1+n_arm*f_arm,1+n_arm*f_arm)
  real(kind=real_kind)         :: A(1+n_arm*f_arm,1+n_arm*f_arm)  
  real(kind=real_kind)         :: phi(1+n_arm*f_arm,1+n_arm*f_arm)  
  real(kind=real_kind)         :: work((1+n_arm*f_arm)*2)
  real(kind=real_kind), allocatable, save   :: Eigenvalues(:)  
  real(kind=real_kind), allocatable, save   :: Eigenvectors(:,:)  

  logical, save                :: first_run = .true.
  integer, save                :: f_arm0, n_arm0   




  N      = f_arm * n_arm + 1
  l0     = Re_arm/sqrt(dble(n_arm))      ! must be changed if a nontrivial U-Matrix will be used

!! compute Eigenvectors and Eigenvalues for the problem:
!! see: M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!

    N=1+n_arm*f_arm
!
! the first index pertains to the center of the star
!


 first_run = .true.    ! allocated are not save in all cases thus we need this as a workaround   


eig:  if(first_run .or. n_arm .ne. n_arm0 .or. f_arm .ne. f_arm0) then
    
    first_run = .false.
    f_arm0    = f_arm
    n_arm0    = n_arm

    if(allocated(Eigenvalues)) then
      deallocate(Eigenvalues)
      deallocate(Eigenvectors)
    endif

    allocate(Eigenvalues(N))
    allocate(Eigenvectors(N,N))

    U = 0            ! inverse of static bond correlation matrix <l_i*l_j>/l**2  (here  1)
    M = 0            ! bead-to-bond vector transformation matrix
    A = 0            ! A = Mt(U)M   (were the first col and row of U are 0)
    H = 0            ! H = hydrodynamic matrix (so far here 1) friction is lumped in the Wl4 parameter
   
    forall (i=2:N)   U(i,i) = 1.0
    forall (i=1:N)   M(1,i) = 1.0/N
    forall (i=2:N)   M(i,i) = 1.0
    forall (i=1:N)   H(i,i) = 1.0
    do j=1,f_arm
      i1 = 2+n_arm*(j-1)
      i2 = j*n_arm
      forall (i=i1:i2)   M(i+1,i) = -1_real_kind
    enddo
    do j=1,f_arm
      i = 2+(j-1)*n_arm
      M(i,1) = -1
    enddo
   
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
   
   ! and for completeness HA !
   do i=1,N
    do j=1,N
     X(i,j) = dot_product(H(i,1:N),A(1:N,j))
    enddo
   enddo

   ! check symmetry
   do i=1,N
    do j=i+1,N
     if(abs(X(i,j)-X(j,i)) > 1e-5) then
       write(6,*)'ERROR nonsymmetric matrix HA '
       write(6,*)'USE star2 to deal with cases where H is nontrivial'
       pericostar_sqt = 0
       return
     endif
    enddo
   enddo

   
   Eigenvectors = X
   job = 1            ! signal to ssiev that the matrix is to be replaced by the eigenvectors
   !
   ! Slatec:
    call SSIEV (Eigenvectors, N , N, Eigenvalues, WORK, JOB, INFO)
   
 endif eig



! -- fill the (symmetric) phi(n,m) matrix -- ( the L**2 * d_ij(t) values of the papaers) 
philp: do i=1,N
        do j=1,N   
          sum = 0d0
          do l = 2 ,N  ! sjkip the first eigenvalue = 0 
           p = Eigenvalues(l)    
           sum = sum + (  Eigenvectors(i,l)**2+Eigenvectors(j,l)**2- &
                        2*Eigenvectors(i,l)*Eigenvectors(j,l)*exp(-wl4/(l0**4)*p*t))/p
          enddo
          phi(j,i) = (l0**2) * sum 
        enddo
       enddo philp
!
! --> ok phi matrix filled (check p-summation role, fixed middle point --> p odd)
!
  
!
! -- and now the scattering contibutions of the arms
!
      
        s11 = 0         ! the "self term"
        do i=2,n_arm+1
           s11 = s11 + 2 * exp(-(q*q/6d0) * phi(1,i)) + exp(-(q*q/6d0)* phi(i,i)) 
        enddo
        do i=2,n_arm
         do j=i+1,n_arm+1
           s11 = s11 + 2* exp(-(q*q/6d0)* phi(i,j))
         enddo
        enddo

        s12 = 0         ! the "cross term"
        if(f_arm > 1) then
          do i=2,n_arm+1
             s12 = s12 + exp(-(q*q/6d0) * phi(i+n_arm,i))
          enddo
          do i=2,n_arm
             do j=i+1,n_arm+1
               s12 = s12 + 2 * exp(-(q*q/6d0) * phi(i+n_arm,j))
             enddo
          enddo
        endif

        sum = (1+f_arm*(s11+(f_arm-1)*s12))/N**2

     pericostar_sqt = sum


! tentatively add diffusion !
     if(diff == 0d0) then
        pericostar_sqt = pericostar_sqt * exp( -q*q* wl4/(3*f_arm*n_arm*l0*l0) * t)
     endif
!


end function pericostar_sqt
