      FUNCTION th_roufbead (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rouse star an with 'hevy' beads each nrep beads, witn extra friction and contrast <--------                      
!                                                                       
!                                    
      implicit none                                   
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      real :: pa, qq, x, th_roufbead, zpi
      integer :: npar, ini, nparx, ier, iert, ierq
      integer                     :: mbuf
      integer, intent(inout)      :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)           ! name des parameters n
      real, intent(inout)         :: params(mbuf)          ! value des parameters n
      REAL(8) temp, qz, tau, eta, yz, SQ_rouse, a, b, xi 
      REAL(8) a0, sum, sumnorm, q_width, dqw, qzz, fn 
      REAL(8) epsilon, diff, wl4, re_arm , extra_frict, extra_contrast_factor
      REAL qget, tget, tauget 

      integer :: f_arm, n_arm, n_repeat, contrast_select
      double precision :: plinear1_sqt, sqt0, sqt
       
      logical :: x_is_tau = .true.                                                                 
                                                                        
                                                                        
      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'roufbead' 
         nparx = 9 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_roufbead = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'intensit'     ! general prefactor 
         parnam (2) = 'wl4'          ! Rouse Rate W l**4 in A**4/ns
         parnam (3) = 're_arm'       ! Re of one arm     in A
         parnam (4) = 'n_arm '       ! number of beads per arm
         parnam (5) = 'contrast'     ! contrast selsct: 0=full coherent, 1:=-AB-, 2:blockwise pseudo inc.
         parnam (6) = 'extfrict'     ! factor for frictin of each nrepeat bead
         parnam (7) = 'nrepeat'      ! each nrepeat bead has assigned an extra friction
         parnam (8) = 'com_diff'     ! common diffusion in cm**2/s, if =0 use Rouse default
         parnam (9) = 'extcontr'     ! contrast factor for the beads with extrafriction 
                                                                        
!                                                                       
         th_roufbead = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      tau             = x 
      a0              = pa (1) 
      wl4             = abs (pa (2) ) 
      re_arm          = abs (pa (3) ) 
      n_arm           = nint (pa (4) ) 
      contrast_select = nint (pa (5) )
      extra_frict     = abs(pa(6))
      n_repeat        = nint(pa(7))

                              ! in cm**2/sec                            
      diff = abs (pa (8) ) 
      diff = diff * 1d-9 / 1d-16 

      extra_contrast_factor = pa(9)
                                                                        
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
        
!!!!----> observe sqt0(q=0) is not normalize to one but yields N**2 for a coherent scattering situation with
!!!!----> uniform b (contrast+select=0) (see implementation in plinear_v1)

      if(x_is_tau) then                                                                
        sqt0       =     plinear1_sqt(qz,0d0,n_arm,Re_arm,Wl4,extra_frict,n_repeat,contrast_select,diff)
        sqt        =     plinear1_sqt(qz,tau,n_arm,Re_arm,Wl4,extra_frict,n_repeat,contrast_select,diff) 
        th_roufbead=a0 * sqt/sqt0
      else
        th_roufbead=a0 * plinear1_sqt(qz,tau,n_arm,Re_arm,Wl4,extra_frict,n_repeat,contrast_select,diff)
      endif 

      if(diff .ne. 0d0) then
        th_roufbead = th_roufbead * exp(-qz*qz*diff*tau)
      endif 

      RETURN 
      END FUNCTION th_roufbead



   
   
double precision function plinear1_sqt(q,t,n_arm,Re_arm,Wl4,extra_frict,n_repeat,contrast_select,diff)
!-------------------------------------------------------------------------------------------------------------
!! star following the derivation of:
!!      M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!
!! this module written by: michael Monkenbusch, Jcns
!!
  implicit none

  double precision, intent(in)  :: q
  double precision, intent(in)  :: t
  integer, parameter            :: f_arm = 1
  integer,          intent(in)  :: n_arm 
  double precision, intent(in)  :: Re_arm
  double precision, intent(in)  :: Wl4
  double precision, intent(in)  :: extra_frict
  integer,          intent(in)  :: n_repeat
  integer,          intent(in)  :: contrast_select   ! 0=uniform coherent, 1=-AB- blockwise, 2=pseudo incoherent blockwise
  double precision, intent(in)  :: diff

  double precision, parameter   :: Pi=3.141592654d0
  integer                       :: N
  integer                       :: i, j, ip, l, i1, i2, job, info
  integer                       :: ib
  double precision              :: l0         ! effective segment length
  double precision              :: tauR       ! basic Rouse-time of a chain of 2 arm lengths

  double precision              :: p, s11, s12, dec

  integer, parameter   :: sort_x_and_perm_inc=2
  integer, parameter   :: sort_x_and_perm_dec=-2
  integer, parameter   :: sort_x_set_perm_inc=1
  integer, parameter   :: sort_x_set_perm_dec=-1



 logical, save                :: first_run = .true.
 logical, save                :: new_parameters
 logical, save                :: inc_summing 
 integer, save                :: f_arm0, n_arm0, n_repeat0, contrast_select0
 double precision, save       :: extra_frict0   

 real                         :: U(1+n_arm*f_arm,1+n_arm*f_arm)
 real, allocatable,save       :: H(:,:)
 real                         :: X(1+n_arm*f_arm,1+n_arm*f_arm)
 real                         :: M(1+n_arm*f_arm,1+n_arm*f_arm)
 real, allocatable,save       :: A(:,:)  
 real, allocatable,save       :: r(:)           ! random numbers
 real, allocatable,save       :: b(:)           ! scattering length contibutions
 real                         :: HA(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: Ah(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: phi(1+n_arm*f_arm,1+n_arm*f_arm)  
 real                         :: work(4*(1+n_arm*f_arm))
 complex                      :: csum
 complex, allocatable, save   :: Eigenvalues(:)  
 complex, allocatable, save   :: Eigenvectors(:,:) 
 complex                      :: CX(1+n_arm*f_arm,1+n_arm*f_arm)
 complex, allocatable, save   :: CM(:,:)


 real                         :: std_diffusion
 real                         :: std_bead_friction
 real                         :: effective_diffusion

 integer                      :: perm(1+n_arm*f_arm)
 real                         :: lambdas(1+n_arm*f_arm)
 integer                      :: ier, ilp
 character(len=20)            :: fname




  N      = f_arm * n_arm + 1
 
! write(6,*) N

! write(6,*)dble(n_arm), Re_arm

  l0     = Re_arm/sqrt(dble(n_arm))

!! compute Eigenvectors and Eigenvalues for the problem:
!! see: M. Guenza, M. Mormino and A. Perico, Macromolecules 1991, 24, 6166-6174
!! and  M. Guenza, A. Perico Macromolecules 1993, 26, 4196-4202
!!
!! modification with increased friction after each n_repet beads
!! 


! first_run = .true.

    N=1+n_arm*f_arm


    new_parameters =    n_arm           .ne. n_arm0                &
                  .or.  f_arm           .ne. f_arm0                &
                  .or.  n_repeat        .ne. n_repeat0             &
                  .or.  contrast_select .ne. contrast_select0      &
                  .or.  abs(extra_frict-extra_frict0) > 1e-6       

   

eig2: if(  first_run .or. new_parameters ) then
    

!eig2: if(t==0) then
    first_run        = .false.
    f_arm0           = f_arm
    n_arm0           = n_arm
    n_repeat0        = n_repeat
    extra_frict0     = extra_frict
    contrast_select0 = contrast_select

    if(allocated(Eigenvalues)) then
      deallocate(A)
      deallocate(H)
      deallocate(Eigenvalues)
      deallocate(Eigenvectors)
      deallocate(CM)  
      deallocate(r)
      deallocate(b)
    endif

    allocate(Eigenvalues(N))
    allocate(A(N,N))
    allocate(H(N,N))
    allocate(Eigenvectors(N,N))
    allocate(CM(N,N))
    allocate(r(N))
    allocate(b(N))

!   CALL init_random_seed()
!   CALL RANDOM_NUMBER(r)



    U = 0            ! inverse of static bond correlation matrix <l_i*l_j>/l**2  (here  1)
    M = 0            ! bead-to-bond vector transformation matrix
    A = 0            ! A = Mt(U)M   (were the first col and row of U are 0)
    H = 0            ! H = hydrodynamic matrix (so far here 1) friction is lumped in the Wl4 parameter
   
    b = 1.0          ! contrast factors  

 forall (i=2:N)   U(i,i) = 1.0
 forall (i=1:N)   M(1,i) = 1.0/N
 forall (i=2:N)   M(i,i) = 1.0
 forall (i=1:N)   H(i,i) = 1

 inc_summing = .false.

 select case(contrast_select) 
   case (0)
     b = 1
     if(new_parameters) write(6,*)'contrast is uniform, selector coherent'
   case (1)
     forall (i=1:N)  b(i) = b(i)*(-1)**(nint(0.5+(n_repeat/2.0+i-1)/n_repeat))
     if(new_parameters) write(6,*)'contrast is alternating left right block'
   case (2)
     b = 1
     inc_summing = .true.
     if(new_parameters) write(6,*)'contrast computation: blockwise incoherent'
   case default
    b = 1
 end select 

! loopit: do ilp=1,2

 H(1,1) = H(1,1)/extra_frict / 0.5
 do i=n_repeat,n_arm-n_repeat,n_repeat
      i1 = 1+i
      H(i1,i1) = H(i1,i1)/extra_frict 
 enddo

!?? new
 H(n,n) = 1.0/extra_frict / 0.5
!?? new

 forall (i=2:N)   U(i,i) = 1.0
 forall (i=1:N)   M(1,i) = 1.0/N
 forall (i=2:N)   M(i,i) = 1.0

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


!! just for esthetics: normalize eigenvectors

do l=1,n
 Eigenvectors(1:N,l) =  Eigenvectors(1:N,l)/sqrt(dot_product(conjg(Eigenvectors(1:N,l)),Eigenvectors(1:N,l)))
enddo


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
!write(6,*)'Eigenvalues, CM(ll): '
!do i=1,N
!   write(6,'(i5,2f12.6,5x,2f12.6)') i,eigenvalues(i), CM(i,i)
!enddo

! sort here
 
 lambdas(1:n) = Eigenvalues(1:n)
 perm(1:n)    = [(i,i=1,n)]
 call spsort(lambdas,n,perm,sort_x_set_perm_inc,ier)

if(new_parameters) then
 write(6,*)'Eigenvalues sorted: '
 write(6,'(i4,2f12.8,10x,2f12.8)') (i,Eigenvalues(perm(i)),cm(perm(i),perm(i)),i=1,n)
 open(10,file='ew.dgli')
   write(10,*)' f_arm        := ', f_arm
   write(10,*)' n_arm        := ', n_arm
   write(10,*)' extra_frict  := ', extra_frict
   write(10,*)' n_repeat     := ', n_repeat
   write(10,'(i4,2f12.8,10x,2f12.8)') (i,Eigenvalues(perm(i)),cm(perm(i),perm(i)),i=1,n)
 close(10)

 do i=1,n
  write(fname,'(a,i0.0,a)')'ev',i,'.dgli'
  open(11,file=fname)
   write(11,*)' f_arm        := ', f_arm
   write(11,*)' n_arm        := ', n_arm
   write(11,*)' extra_frict  := ', extra_frict
   write(11,*)' n_repeat     := ', n_repeat
   write(11,*)' eigenvector  := ',i
   write(11,*)' n            := ',n
   write(11,'(i4,2f12.8)') (j,Eigenvectors(j,perm(i)),j=1,n)
  close(11)
 enddo

 write(6,*)'Contrastfactors and frictions: '
 write(6,'(i4,2f12.8)') (i,b(i),1.0/H(i,i),i=1,n)
 open(12,file='bfric.dgli')
   write(12,'(i4,2f12.8)')(i,b(i),1.0/H(i,i),i=1,n)
   if(inc_summing)write(12,*)'blockwise incoherent summing!'
   if(inc_summing)write(6 ,*)'blockwise incoherent summing!'
 close(12)
endif
!* endif eig

endif eig2

! -- fill the (symmetric) phi(n,m) matrix --


philp2: do i=1,N
  do j=1,N
   csum = 0
   do l=1,N
    if(abs(CM(l,l)) > 0.0001) then
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
      
     if(.not. inc_summing) then
        s11 = b(1)*b(1)         ! the "self term"
        do i=2,n_arm+1
           s11 = s11 + 2 * exp(-(q*q/6d0) * phi(1,i))*b(1)*b(i) + exp(-(q*q/6d0)* phi(i,i))*b(i)*b(i) 
        enddo
        do i=2,n_arm
         do j=i+1,n_arm+1
           s11 = s11 + 2* exp(-(q*q/6d0)* phi(i,j))*b(i)*b(j)
         enddo
        enddo 
!        plinear1_sqt = (s11)/(N*N)
        plinear1_sqt = (s11)
     else

! --- independent summation of sections between "heavy" links --> pseudo incoherent approx 
       s11 = 0
       do ib=1,n-n_repeat,n_repeat
        do i=ib,ib+n_repeat-1
         do j=ib,ib+n_repeat-1
            s11 = s11 +  exp(-(q*q/6d0)* phi(i,j))*b(i)*b(j)
         enddo
        enddo
       enddo

! ..and the rest
       ib = j
       do i=ib,n
         do j=ib,n
            s11 = s11 +  exp(-(q*q/6d0)* phi(i,j))*b(i)*b(j)
         enddo
        enddo
!        plinear1_sqt = (s11)/(N*n_repeat)
        plinear1_sqt = (s11)
     endif


!        s12 = 0         ! the "cross term"
!        if(f_arm > 1) then
!        do i=2,n_arm+1
!           s12 = s12 + exp(-(q*q/6d0) * phi(i+n_arm,i))
!        enddo
!        do i=2,n_arm
!           do j=i+1,n_arm+1
!             s12 = s12 + 2 * exp(-(q*q/6d0) * phi(i+n_arm,j))
!           enddo
!        enddo
!        endif
!
!        plinear1_sqt = (1+f_arm*(s11+(f_arm-1)*s12))/N**2
!

! write(6,*) s11, s12, sum, sqrt(a*a*n_arm)

! tentatively add diffusion !
     if(diff == 0d0) then
        std_diffusion       =  wl4/(3*f_arm*n_arm*l0*l0)
        std_bead_friction   =  1/(n*std_diffusion)
        effective_diffusion =  1/(sum([(1.0/H(i,i),i=1,n)])*std_bead_friction)
        plinear1_sqt = plinear1_sqt * exp( -q*q* effective_diffusion * t)
     endif
!



!> do i=1,n-1
!>  do j=i+1,n
!>   H(i,j) = 0.25/sqrt(phi(i,j)/phi(1,2))
!>   H(j,i) = H(i,j)
!>  enddo
!> enddo
!>
!> enddo loopit

end function plinear1_sqt
