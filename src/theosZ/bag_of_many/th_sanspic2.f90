      FUNCTION thsanspic2 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!
! scattering of a platelety with a lattice of pores
!
      implicit none
 
      double precision, parameter :: Pi=3.141592654d0

      real    :: thsanspic2

      real    ::   x, pa
      integer ::   ini, npar, nparx


      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
!		DATA zpi / 6.283185 / 


      double precision :: gauss_rand
      double precision :: y, yd
      double precision, parameter :: radeg=57.29577951d0

      double precision  :: sigma(2), sigma_lam, psi, theta, lambda, k(3), ki(3), kf(3), Q(3)
      double precision  :: pore_scat1, L, R, Phi
      double precision  :: pore_scat2, a, r_domain, f_domain, psi_domain, sq_hexlat2
      double precision  :: s, sum, dtheta, qq, dq, qx, qy, qz


      integer :: i, Nsamp, j, Ntheta, m
      integer :: Narad
      double precision  :: sumr, da, sec1, sec2 


      integer :: iadda, ier
      common/thiadd/iadda
      real    :: xh
!
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'sanspic2' 
         nparx = 16
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           thsanspic2 = 0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1)  = 'intensit' 
         parnam (2)  = 'a'          ! distance of pores 
         parnam (3)  = 'r'          ! radius of pores
         parnam (4)  = 'l'          ! length of pores
         parnam (5)  = 'nsamp'      ! no mc-sampling points per resolution
         parnam (6)  = 'dq'         ! delta-q
         parnam (7)  = 'phi'        ! sample rotation angle
         parnam (8)  = 'lambda'     ! wavelength
         parnam (9)  = 'fwhmsel'    ! selector width
         parnam (10) = 'collim1'    ! collimation fwhm in degrees
         parnam (11) = 'collim2'    ! collimation fwhm in degrees
         parnam (12) = 'r_domain'   ! domain radius
         parnam (13) = 'f_domain'   ! filling factor of domain
         parnam (14) = 'psi_dom'    ! domain rotation angle
         parnam (15) = 'sec1'       ! mask-sector angle start
         parnam (16) = 'sec2'       ! mask-sector-angle end
      
  
!                                                                       
         thsanspic2 = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----
      a           = pa(2)
      R           = pa(3)
                  !=== ............Pore radius
      L           = pa(4)
                  !====== .........Pore length

     Nsamp       = nint(pa(5))
                  !=== ............ Number of random Q-samples in resolution

      dq          = pa(6)
                  !===== ...Q-width of computation
                  !====............Pores distance
      phi          = pa(7)/radeg
                  !===.............Sample rotation
      lambda       = pa(8)
                  !===.............Wavelength
      sigma_lam    = pa(9)*lambda
                  !======.............FWHM lambda spread
      sigma(1)     = pa(10)/radeg
      sigma(2)     = pa(11)/radeg
                  !===.............FWHM angular divergence

      r_domain    = pa(12)
                  !====............Domain width
      f_domain    = pa(13)
                  !=== ............Pore fraction
      psi_domain  = pa(14)/radeg

      sec1        = pa(15)/radeg
      sec2        = pa(16)/radeg


!
! --> get resolution params from paramsec if available
!
       call        parget('collim1 ',xh ,iadda,ier)
       if(ier.eq.0 .and. sigma(1).eq.0d0 )  sigma(1)  = xh/radeg
       call        parget('collim2 ',xh ,iadda,ier)
       if(ier.eq.0 .and. sigma(2).eq.0d0 )  sigma(2)  = xh/radeg
       call        parget('lambda  ',xh ,iadda,ier)
       if(ier.eq.0 .and. lambda  .eq.0d0 )  lambda    = xh
       call        parget('fwhmsel ',xh ,iadda,ier)
       if(ier.eq.0 .and. sigma_lam.eq.0d0 )  sigma_lam   = xh

       call        parget('sec1    ',xh ,iadda,ier)
       if(ier.eq.0 .and. sec1.eq.0d0 )   sec1      = xh/radeg
       call        parget('sec2    ',xh ,iadda,ier)
       if(ier.eq.0 .and. sec2.eq.0d0 )   sec2      = xh/radeg
      



      qq    = x                                     

!! radial average 
      Narad = 2*(sec2-sec1)*qq/dq
      if(Narad.le.0) Narad = 1
      da    = (sec2-sec1)/Narad
      sumr = 0    

      do i=0,Narad

       qy = qq*sin(i*da+sec1)
       qz = qq*cos(i*da+sec1)

       theta = -2*asin(sqrt(qy*qy+qz*qz)*lambda/4/Pi)
       psi   = atan2(qz,qy)+Pi
    
       sum = 0
       do m=1,Nsamp
!         call get_Q(lambda, sigma_lam, sigma(1) , theta, psi, ki, kf, Q )
         call get2_Q(lambda, sigma_lam, sigma , theta, psi, ki, kf, Q )
         s = pore_scat1( Q, L, R, phi)*sq_hexlat2(q, a, r_domain, f_domain,psi_domain) 
!         write(6,'(9f12.6)')ki,kf,Q
         sum = sum + s
       enddo
       sum = sum / Nsamp
!       write(6,*) qy,qz, sum, theta*radeg, psi*radeg
       sumr = sumr + sum
      enddo
      sumr = sumr/Narad

      thsanspic2 = pa(1)*sumr

                                                             
      RETURN 
      END FUNCTION thsanspic2



     double precision function sq_hexlat2(q, a, r_domain, f_domain, psi_domain)
!    --------------------------------------------------------------------------
! --> structure factor of a finite (r_domain)  fractionally filled 
!     triangular lattice   
!

     implicit none

     real rand
    
     double precision, parameter :: Pi=3.141592654d0

     double precision, intent(in)  :: q(3)       !  Q
     double precision, intent(in)  :: a          !  lattice constant of triangular
     double precision, intent(in)  :: r_domain   !  rdius of domain
     double precision, intent(in)  :: f_domain   !  occupancy fo domain
     double precision, intent(in)  :: psi_domain !  rotation angle of lattice

     double precision, save        :: a0          !  lattice constant of triangular
     double precision, save        :: r_domain0   !  rdius of domain
     double precision, save        :: f_domain0   !  occupancy fo domain
     double precision, save        :: psi_domain0 !  occupancy fo domain

     integer, parameter            :: lmax=10000

     double precision              :: b1(2), b2(2), c(2), dq, qa, central, den
     double precision              :: dbesj0, dbesj1, cp, sp
     double precision, save        :: lattice(2,lmax)
     integer, save                 :: nlat, nlat0
     integer                       :: i,j,k,n

     complex*16                    :: ci=(0d0,1d0)
     complex*16                    :: csum, cq

! fill cluster  (maybe done just once)
     if(a.ne.a0 .or. r_domain   .ne. r_domain0  & 
                .or. f_domain   .ne. f_domain0  & 
                .or. psi_domain .ne. psi_domain0) then
       a0           = a
       f_domain0    = f_domain
       r_domain0    = r_domain
       psi_domain0  = psi_domain

    
         b1(1) = a
         b1(2) = 0
         b2(1) = a/2
         b2(2) = sqrt(3d0)/2d0*a

! rotate it
         cp = cos(psi_domain)
         sp = sin(psi_domain)
         c(1) = cp*b1(1) + sp*b1(2)
         c(2) =-sp*b1(1) + cp*b1(2)
         b1   = c
         c(1) = cp*b2(1) + sp*b2(2)
         c(2) =-sp*b2(1) + cp*b2(2)
         b2   = c
         write(6,'(a,7f12.6)')'ps:', psi_domain,cp,sp,b1,b2
! rotate it
         
         nlat0 = 0
         nlat  = 0
         n = nint(r_domain/a)+1
         do i=-n,n
          do j=-n,n
            c = i*b1 + j*b2
            if(c(1)**2+c(2)**2 .lt. r_domain**2) then
             nlat0 = nlat0+1
             if(rand().lt.f_domain) then
               nlat = nlat + 1
               if(nlat.gt.lmax) then
                 write(6,*)'sq_hexagonal: too large domain! '
                 stop
               endif
               lattice(1:2,nlat) = c(1:2)
             endif
            endif 
          enddo
         enddo
         write(6,'(a,2I8)')'nlat: ',nlat0,nlat
         if(nlat0.le.6) then 
           Write(6,*)'NOT removing central peak..'
         else
           Write(6,*)'Removing central peak !'
         endif
    !
    endif
!
! determine sof q
     csum = 0
     do i=1,nlat
        cq = (lattice(1,i)*Q(2) + lattice(2,i)*Q(3))*ci
        csum = csum + exp(cq)
     enddo
!     csum = csum/nlat    
!     sq_hexlat2 = csum*conjg(csum)/nlat/nlat
      sq_hexlat2 = csum*conjg(csum)

! subtract the central peak by using a circular embedding
!
! the continuum density in the domain is
 
      if(nlat0.gt.6) then
    
       den = nlat / (Pi*r_domain**2) 
       qa  = sqrt(Q(2)**2+Q(3)**2)
       central = ( den*2*Pi*r_domain*dbesj1(qa*r_domain)/qa )**2

       sq_hexlat2 = sq_hexlat2 - central  ! here we subtract the shape scattering of the domain 
                                          ! since it is embedded into other domains
      endif

! normalize to the domain size
      sq_hexlat2 = sq_hexlat2 / nlat0 
      


     return   
     end function sq_hexlat2




     double precision function pore_scat1( Q, L, R, phi )
!    ----------------------------------------------------
! scattering crossection (in terms of Volume**2) of one cylindrical pore
! of length L, radius R axis along x and rotated by phi around z
!
     implicit none
     double precision, intent(in)  :: Q(3), L, R, phi
     double precision              :: qx, qr, fl, fr, qrot(3), cr, sr
     double precision              :: dbesj1
     double precision, parameter   :: Pi=3.141592654d0

     cr      = cos(phi)
     sr      = sin(phi)
     qrot(1) =  cr*Q(1) + sr*Q(2)
     qrot(2) = -sr*Q(1) + cr*Q(2)
     qrot(3) = Q(3)


     qx = qrot(1) 
     qr = sqrt(qrot(2)**2+qrot(3)**2) 
     
     fr = 2*dbesj1(qr*R)/(qr*R)
     if(qx.ne.0) then 
       fl = sin(qx*L/2)/(qx*L/2)
     else
       fl = 1d0
     endif

     pore_scat1 = (fl*fr)**2 
    

     return
     end function pore_scat1



     subroutine get2_Q( lambda, sigma_lambda, sigma_angle, theta, psi, ki, kf, Q)
!    ---------------------------------------------------------------------------
!    create randomized ki and kf and Q
!
     implicit none

     double precision, intent(in)  :: lambda, sigma_lambda, sigma_angle(2), theta, psi
     double precision, intent(out) :: ki(3), kf(3), Q(3)

     double precision              :: lambda_r, gauss_rand

     lambda_r = abs( lambda + gauss_rand(sigma_lambda) )
 
     call get_kvector( lambda_r, sigma_angle(1), 0d0, 0d0, ki)
     call get_kvector( lambda_r, sigma_angle(2), theta, psi, kf)

     Q = kf - ki

     return
     end subroutine get2_Q
     
     

     subroutine get2_Q2( lambda, sigma_lambda, sigma_angle, theta, psi, ki, kf, Q)
!    ----------------------------------------------------------------------------
!    create randomized ki and kf and Q
!
     implicit none

     double precision, intent(in)  :: lambda, sigma_lambda, sigma_angle(2), theta(2), psi(2)
     double precision, intent(out) :: ki(3), kf(3), Q(3)

     double precision              :: lambda_r, gauss_rand

     lambda_r = abs( lambda + gauss_rand(sigma_lambda) )
 
     call get_kvector( lambda_r, sigma_angle(1), theta(1), psi(1), ki)
     call get_kvector( lambda_r, sigma_angle(2), theta(2), psi(2), kf)

     Q = kf - ki

     return
     end subroutine get2_Q2
     
     


     subroutine get_Q( lambda, sigma_lambda, sigma_angle, theta, psi, ki, kf, Q)
!    ---------------------------------------------------------------------------
!    create randomized ki and kf and Q
!
     implicit none

     double precision, intent(in)  :: lambda, sigma_lambda, sigma_angle, theta, psi
     double precision, intent(out) :: ki(3), kf(3), Q(3)

     double precision              :: lambda_r, gauss_rand

     lambda_r = abs( lambda + gauss_rand(sigma_lambda) )
 
     call get_kvector( lambda_r, sigma_angle, 0d0, 0d0, ki)
     call get_kvector( lambda_r, sigma_angle, theta, psi, kf)

     Q = kf - ki

     return
     end subroutine get_Q
     
     

     subroutine get_Q2( lambda, sigma_lambda, sigma_angle, theta, psi, ki, kf, Q)
!    ----------------------------------------------------------------------------
!    create randomized ki and kf and Q
!
     implicit none

     double precision, intent(in)  :: lambda, sigma_lambda, sigma_angle, theta(2), psi(2)
     double precision, intent(out) :: ki(3), kf(3), Q(3)

     double precision              :: lambda_r, gauss_rand

     lambda_r = abs( lambda + gauss_rand(sigma_lambda) )
 
     call get_kvector( lambda_r, sigma_angle, theta(1), psi(1), ki)
     call get_kvector( lambda_r, sigma_angle, theta(2), psi(2), kf)

     Q = kf - ki

     return
     end subroutine get_Q2
     
     




     subroutine get_kvector( lambda, sigma_angle, theta, psi, k )
!    ------------------------------------------------------------
!    get a kvector with randomize divergence width sigma
!
     implicit none

     double precision, intent(in)  :: lambda, sigma_angle, theta, psi
     double precision, intent(out) :: k(3)
     double precision              :: gauss_rand
     double precision              :: eta1, eta2
     double precision              :: ct, st, cp, sp, c1, s1, c2, s2, h(3)

     double precision, parameter   :: Pi=3.141592654d0

! a k vector with direction x

     k(1) = 2*Pi/lambda
     k(2) = 0
     k(3) = 0

! randomize by rotations around z and y
   
!    write(6,*)'lambda= ',lambda
!    write(6,*)'sgma  = ',sigma_angle
!    write(6,*)'psi   = ',psi
!    write(6,*)'theta = ',theta
!    write(6,*)'k     = ',k


     eta1 = gauss_rand(sigma_angle)
     eta2 = gauss_rand(sigma_angle)
     c1   = cos(eta1)
     c2   = cos(eta2)
     s1   = sin(eta1)
     s2   = sin(eta2)



     h(1) = c1*k(1) +        s1*k(3)
     h(2) =           k(2)
     h(3) =-s1*k(1) +        c1*k(3)

     k = h

!     write(6,*)'k     = ',k

     h(1) =  c2*k(1) + s2*k(2)
     h(2) = -s2*k(1) + c2*k(2)
     h(3) =                     k(3)

     k = h

!     write(6,*)'k     = ',k

! and rotation to the scattering angle theta
     ct = cos(theta)
     st = -sin(theta)
     cp = cos(psi)
     sp = -sin(psi)


     h(1) =  ct*k(1) + st*k(2)
     h(2) = -st*k(1) + ct*k(2)
     h(3) =                     k(3)
     
     k = h

!     write(6,*)'k     = ',k

! and psi

     h(1) =  k(1)
     h(2) =         cp*k(2) + sp*k(3)
     h(3) =        -sp*k(2) + cp*k(3)


     k = h

!     write(6,*)'k     = ',k

     return
     end subroutine get_kvector

 
!/***************************************************************************C*/
!/* Methode: nach Numerical Recipies 7.2  (Box-Muller Methode)               C*/
!/*                                                                          C*/
!/***************************************************************************C*/
  double precision function gauss_rand(width)
! -------------------------------------------> FWHM width

  implicit none

  real rand  

  double precision, intent(in)  :: width 
  double precision              :: r1, r2, r, f, grand1, grand2
  double precision, parameter   :: Pi=3.141592654d0
 
    r = 2
    do while ((r >= 1d0) .or. (r == 0)) 
      r1 = 2d0*rand()-1d0
      r2 = 2d0*rand()-1d0
      r  = r1*r1 +r2*r2
    enddo
 
    f = sqrt(-2*log(r)/r)
    grand1 = r1*f
!   grand2 = r2*f ! verwerfe 2te um keine Korrelationen zu haben 
 
    gauss_rand  = grand1*width/sqrt(2*Pi)   ! should match with FWHM then
  
    return 
    end function gauss_rand

