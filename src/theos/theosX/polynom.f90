!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! 
!! RPA model: 3 (polymeric) components 
!! the S-matrix expressions are derived for 3 components without extra mutial interaction
!! the expressions come form Eqs. 10-18 of the Paper
!! Ackcasu and Tombakoglu, Macromoleculse 23, 607 (1990)
!! 
!! The method relies on numerical determination of the denominator and nuerator
!! coefficients of the polynomials in s that give the matrix elements Sii
!! the numerics is a bit touchy since we may habe higher order polynomials that
!! may obtain very large values if the support of the table for fitting is not well
!! adjusted. The spcing is controlled by dss, the number of (extra) support points
!! by npp_plus and a general scaling of the polynoimal values by sfak0
!! as long as those values are suitably set to avoid too extreme numerical acrobatics
!! the results should NOT furher depend on them!
!! I am still looking for a relaible automatic setting alogorithm for these parameters
!! until that suceeded, some careful "manual" adjustment may be needed.
!!
!! michael monkenbusch, JCNS-1, Forschungszentrum Juelich
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


module polynom
use rpa_laplace
 
 double precision    :: sfak0       = 1d-3        ! base scaling factor for polynom data
 double precision    :: dss         = 0.01d0      ! modification factor for 
 double precision    :: ss_offset   = 0d0         ! offset for polynom evaluation range
                                                  ! step width for polynomial evaluation to fit coeffs
 double precision    :: rwarn_limit = 1d-12       ! limit for accuracy warning of polynomial model
 integer             :: npp_plus    = 1           ! additional points in table for polynomila fits


 integer             :: is_sel                    ! naming ? Is it needed ??

contains


       subroutine polyfit(x, y, m, a, n, r)
! -----====================================----------------------------
!  polynomial fit to data in tables x(1..m), y(1..m)  (input)
!  coefficients are returned in a(0..n) (output)
!  max. degree is given by n (input)
!  residual error is given by r (rms) (output)
!
! use nag routine(s):
!     f04aef       real symmetric equations solution
! ---------------------------------------------------------------------
!
       implicit real*8 (a-h,o-z)
       parameter(maxdeg=100)
       parameter(mdim=maxdeg+1)

       dimension am(mdim,mdim),b(mdim,1),c(mdim,1),bb(mdim,1), &
                 wkspce(mdim)
       dimension al(mdim,mdim)

       dimension xx(2*mdim)
       dimension x(m),y(m),a(*)

       ifail = -1

       if(n.gt.maxdeg) then
         r = 9.99999d99
         write(6,*)'polyfit degree n=',n,' is larger than maxdeg'
         return
       endif

       if(n+1.gt.m) then
         n = m-1
         write(6,*)'polyfit: too few points   (',m,')'
         write(6,*)'         degree reduced to ',n
       endif

! --- x-scalingfactor ----
       xscale = 0.1d0 * x(m)
!cc    xscale =1.0d0

       do j=1,2*n+1
         xx(j) = 0.0d0
         b(j,1)= 0.0d0
       enddo

       do i=1,m
         xh = 1.0d0
         do j=1,2*n+1
           xx(j) = xx(j) + xh
           if(j.le.n+1) then
             b(j,1) = b(j,1) + xh * y(i)
           endif
           xh    = xh * x(i)/xscale
         enddo
       enddo

! --- prepare a-matrix ----
      do i=1,n+1
        do j=1,n+1
          am(j,i) = xx(1+(i-1)+(j-1))
        enddo
      enddo
!
! --- solve equation ----
      ier = 0
      call f04aef_wrap(am,mdim,b,mdim,n+1,1,c,mdim, &
                  wkspce)
!
! --- scale c and put it to output
      xh   = 1.0d0
      do i=1,n+1
        a(i) = c(i,1)*xh
        xh   = xh / xscale
      enddo
!!
      if(ier.ne.0) then     !! use simple linear interpolation
        slope = (y(m)-y(1))/(x(m)-x(1))
        a(1)  =  y(1) -slope*x(1)
        do i=2,n+1
          if(i.eq.2) then
            a(i) = slope
          else
            a(i) = 0.0d0
          endif
        enddo
        write(6,*)'polyfit due to ier=',ier,' --> linear interp.'
      endif
!!
!
! --- determine the quality of solution
      r    = 0.0d0
      do i=1,m
        sumv = 0.0d0
        xh   = 1.0d0
        do j=1,n+1
           sumv = sumv + xh * a(j)
           xh  = xh * x(i)
        enddo
        d = sumv-y(i)
!c-test
!c      write(6,*)i,d
!c-teste
        r = r + d*d
      enddo
      r = sqrt(r) / sum( abs( y(1:m) ) )

      if(ier.ne.0) then
        write(6,901) a(1),a(2),r
901     format(' polyfit: a0= ',e14.6,' a1=',e14.6,' r=',e14.6)
      endif


      return
      end

      function polyeval(x,a,n)
!     ========================
! evaluate polynom at x
! coefficients are in a(0:n)
! degree in n
! ==================================================================

      implicit none
      integer  n
      real*8   x, a(0:n), polyeval
      real*8   xx
      integer  i

      polyeval = a(0)
      xx       = x

      do i=1,n
        polyeval = polyeval + a(i)*xx
        xx = xx * x
      enddo

      return
      end


  subroutine coeff_multiplication_polynom(a1,n1,a2,n2, aa, nn)
! ============================================================
 implicit none
 double precision, intent(in)   :: a1(0:n1)    ! coefficient of first polynom
 integer         , intent(in)   :: n1          ! max degree of first polynom
 double precision, intent(in)   :: a2(0:n2)    ! coefficient of second polynom
 integer         , intent(in)   :: n2          ! max degree of second polynom
 double precision, intent(out)  :: aa(0:n1+n2) ! coefficients of resulting product polynom
 integer         , intent(out)  :: nn          ! max degree of resulting polynom

 double precision, parameter    :: epsilon = 1d-9
 double precision               :: epx
 integer :: i, j

 aa = 0
 do i=0,n1
  do j=0,n2
    aa(i+j) = aa(i+j) + a1(i)*a2(j)
  enddo
 enddo

 epx = maxval(abs(aa))*epsilon

 
 d1: do i=n1+n2,1,-1
       if(abs(aa(i)) > epx) then
         nn = i
         exit d1
        else
         aa(i) = 0d0
        endif
     enddo d1

end subroutine coeff_multiplication_polynom
      

function sproduct_polynom( s, r, n) result( p )
  implicit none
  double precision, intent(in)  :: s
  double precision, intent(in)  :: r(n)
  integer         , intent(in)  :: n
  double precision              :: p   ! = (s+r1)(s+r2)..(s+rn)

  integer :: i
  
  p=1d0
  do i=1,n
   p = p * (s+r(i))
  enddo

end function sproduct_polynom 


 subroutine coeff_sproduct_polynom( r, n, aa )
  implicit none
  double precision, intent(in)  :: r(n)     ! r-parameters of the sproduct polynom
  integer         , intent(in)  :: n        ! number pf products = degree of polynomial
  double precision, intent(out) :: aa(0:n)  ! resulting polynom coefficients 

  double precision              :: a3(0:n)
 
  integer :: i, nn, nx


  aa=0
  aa(0)=r(1)
  aa(1)=1d0
  nn=1

  if(n <= 1) return

  do i=2,n
    call coeff_multiplication_polynom([r(i),1d0],1,aa,nn, a3, nx)
    aa = a3
    nn = nx
  enddo
  
end subroutine coeff_sproduct_polynom 





 function Fr(s,n,a,r)
!--------------------
  implicit none
  real   (kind=XPREC)               :: Fr    ! Laplace transform of the n-exp timefunctio
  real   (kind=XPREC) , intent(in)  :: s     ! Laplace 'frequency' variable         
  integer             , intent(in)  :: n     ! number of single simple exponentials represented
  real   (kind=XPREC) , intent(in)  :: a(n)  ! amplitude coefficients  
  real   (kind=XPREC) , intent(in)  :: r(n)  ! rate coefficients

  integer :: i
  
  Fr = 0
  do i = 1, n
    Fr = Fr +  a(i) / ( s + r(i) ) 
  enddo

! eventuell geht mittlerweile auch:  Fn = sum( a/(s+r) )
             
 end function Fr



  function Ss_kernel_real(u) result(v)
 !----------------------------------
  implicit none 
 
  real   (kind=XPREC)              :: v           ! real part of integrand
  real   (kind=XPREC), intent(in)  :: u           ! integration path variable
 
  real   (kind=XPREC)  :: s           ! Laplace variable
  real   (kind=XPREC)  :: Fc          ! relaxation function "matrix" polymer (cc-component) 
  real   (kind=XPREC)  :: Fs1         ! relaxation function of polymer component 1
  real   (kind=XPREC)  :: Fs2         ! relaxation function of polymer component 2
  real   (kind=XPREC)  :: Ss11        ! S_11(q,s)
  real   (kind=XPREC)  :: Ss22        ! S_22(q,s)
  real   (kind=XPREC)  :: Ss12        ! S_12(q,s)
  
  real   (kind=XPREC)  :: Ssxx

! Reference Maple: RpaLaplace_corr_2D_3.mw  (23)-(27)

!! the following variables shoul be global for the module --> move to module section
!!  integer              :: is_sel      ! selector 1: S11, 2: S12, 3: S12, 4: S22 
!!  complex(kind=XPREC)  :: Ssxx        ! return variable
!!  real   (kind=XPREC)  :: phi1        ! volume fraction of polymer component 1
!!  real   (kind=XPREC)  :: phi2        ! volume fraction of polymer component 2
!!  real   (kind=XPREC)  :: Scc00       ! unperturbed structure factor S(Q) of "matrix" polymers
!!  real   (kind=XPREC)  :: S0011       ! unperturbed structure factor S(Q) of polymer 1
!!  real   (kind=XPREC)  :: S0022       ! unperturbed structure factor S(Q) of polymer 2
!!  integer              :: nexpcc      ! number of exp-functions to describe background
!!  integer              :: nexp1       ! number of exp-functions to describe component1
!!  integer              :: nexp2       ! number of exp-functions to describe component2
!!  integer, parameter   :: mex = 32    ! maximum number of exponentials (typ 3-6 !)
!!  real   (kind=XPREC)  :: aexp_cc(mex)! amplitude coeffs for laplace-exp representation of "matrix"
!!  real   (kind=XPREC)  :: rexp_cc(mex)! rate      coeffs for laplace-exp representation of "matrix"
!!  real   (kind=XPREC)  :: aexp_s1(mex)! amplitude coeffs for laplace-exp representation of polymer 1
!!  real   (kind=XPREC)  :: rexp_s1(mex)! rate      coeffs for laplace-exp representation of polymer 1
!!  real   (kind=XPREC)  :: aexp_s2(mex)! amplitude coeffs for laplace-exp representation of polymer 2
!!  real   (kind=XPREC)  :: rexp_s2(mex)! rate      coeffs for laplace-exp representation of polymer 2

!! note expternal parameter: t_param for the time in the resulting S(Q,t)
!! and the apodisation epsilon parameter epap
!! xil the parallel shift of the path
     
 

  s     = u

 
  Fc   =  Fr (s, nexpcc , aexp_cc, rexp_cc )
  Fs1  =  Fr (s, nexp1  , aexp_s1, rexp_s1 )
  Fs2  =  Fr (s, nexp2  , aexp_s2, rexp_s2 )


  select case(is_sel) 

  case(1)
   Ss11 = &
     phi1*S0011*(-Fs1*(-1+phi1+phi2)**2*(Fc*s-1)*Scc00**2+(-1+phi1+phi2)*(((-2+(Fc+Fs2)*s)* &
     S0022*phi2+Fc*s*phi1*S0011)*Fs1-Fc*phi1*S0011)*Scc00-phi2*((phi2*S0022*(Fs2*s-1)+Fs2* &
     S0011*phi1*s)*Fs1-Fs2*phi1*S0011)*S0022)/((-1+phi1+phi2)*Scc00*(Fc*s-1)-Fs1*S0011*phi1*s+ &
     (-Fs2*S0022*s+S0022)*phi2+phi1*S0011)/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
   Ssxx = Ss11
   

   case(2:3) 
    Ss12 = &
     phi1*S0011*phi2*S0022*(-(-1+phi1+phi2)*(Fs1*Fs2*s+Fc-Fs1-Fs2)*Scc00+Fs2*S0011*(Fs1*s-1)* &
     phi1+Fs1*phi2*S0022*(Fs2*s-1))/((-1+phi1+phi2)*Scc00*(Fc*s-1)+(-Fs1*S0011*s+S0011)*phi1- &
     phi2*S0022*(Fs2*s-1))/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
    Ssxx = Ss12
    

   case(4)
    Ss22 = &
     (-Fs2*(-1+phi1+phi2)**2*(Fc*s-1)*Scc00**2+(((-2+(Fc+Fs1)*s)*S0011*phi1+phi2*S0022*Fc*s)* &
     Fs2-Fc*phi2*S0022)*(-1+phi1+phi2)*Scc00-phi1*S0011*((S0011*(Fs1*s-1)*phi1+Fs1*S0022*phi2*s)* &
     Fs2-Fs1*phi2*S0022))*phi2*S0022/((-1+phi1+phi2)*Scc00*(Fc*s-1)-Fs2*S0022*phi2*s+ &
     (-Fs1*S0011*s+S0011)*phi1+phi2*S0022)/((1-phi1-phi2)*Scc00+phi1*S0011+phi2*S0022)
    Ssxx = Ss22
   
   
   case default
    write(6,*)"Error in Ss_kernel, invalid selection (outside 1:4), is_sel=", is_sel
    stop   ! this indicates a programming error, thus we stop here -> fix the bug !

  end select

    v = Ssxx
 

end function Ss_kernel_real


     function S_denominator(s) result(n)
       implicit none
       double precision, intent(in) :: s
       double precision             :: n

       double precision :: S1, S2, Sc, f1, f2, fc
   
       S1 =            phi1 * S0011
       S2 =            phi2 * S0022
       Sc = (1d0-phi1-phi2) * Scc00

      
       fc  =  Fr (s, nexpcc , aexp_cc, rexp_cc )
       f1  =  Fr (s, nexp1  , aexp_s1, rexp_s1 )
       f2  =  Fr (s, nexp2  , aexp_s2, rexp_s2 )


       N   = (S1+S2+Sc)*(S2*(f2*s-1d0)+(f1*s-1d0)*S1+Sc*(fc*s-1d0))

     end function S_denominator
 

     function S11_numerator(s) result(s11)
       implicit none
       double precision, intent(in) :: s
       double precision             :: s11

       double precision :: S1, S2, Sc, f1, f2, fc
   
       S1 =            phi1 * S0011
       S2 =            phi2 * S0022
       Sc = (1d0-phi1-phi2) * Scc00

      
       fc  =  Fr (s, nexpcc , aexp_cc, rexp_cc )
       f1  =  Fr (s, nexp1  , aexp_s1, rexp_s1 )
       f2  =  Fr (s, nexp2  , aexp_s2, rexp_s2 )

       s11=S1*(f1*(f2*s-1d0)*S2*S2+f2*(f1*s-1d0)*S2*S1+fc*(f1*s-1d0)*Sc*S1+&
               f1*(f2*s+fc*s-2d0)*Sc*S2+Sc*Sc*(fc*s-1d0)*f1)
     end function S11_numerator
 

     function S22_numerator(s) result(s22)
       implicit none
       double precision, intent(in) :: s
       double precision             :: s22

       double precision :: S1, S2, Sc, f1, f2, fc
   
       S1 =            phi1 * S0011
       S2 =            phi2 * S0022
       Sc = (1d0-phi1-phi2) * Scc00

      
       fc  =  Fr (s, nexpcc , aexp_cc, rexp_cc )
       f1  =  Fr (s, nexp1  , aexp_s1, rexp_s1 )
       f2  =  Fr (s, nexp2  , aexp_s2, rexp_s2 )

       s22=S2*(f2*(f1*s-1d0)*S1*S1+f2*(f1*s+fc*s-2d0)*S1*Sc+&
               f1*(f2*s-1d0)*S1*S2+fc*(f2*s-1d0)*S2*Sc+Sc*Sc*(fc*s-1d0)*f2)
             
     end function S22_numerator
 

     function S12_numerator(s) result(s12)
       implicit none
       double precision, intent(in) :: s
       double precision             :: s12

       double precision :: S1, S2, Sc, f1, f2, fc
   
       S1 =            phi1 * S0011
       S2 =            phi2 * S0022
       Sc = (1d0-phi1-phi2) * Scc00

      
       fc  =  Fr (s, nexpcc , aexp_cc, rexp_cc )
       f1  =  Fr (s, nexp1  , aexp_s1, rexp_s1 )
       f2  =  Fr (s, nexp2  , aexp_s2, rexp_s2 )

       s12 = -S1*S2*((f1*f2*s-f1-f2+fc)*Sc+f2*S1*(f1*s-1d0)+f1*S2*(f2*s-1d0))
              
             
     end function S12_numerator
 


     function denominator_polynomial(s) result (p)
       implicit none
       double precision, intent(in) :: s
       double precision             :: p

       integer                      :: i

       p = S_denominator(s) 

       do i=1, nexp1
        p = p * ( s + rexp_s1(i))
       enddo

       do i=1, nexp2
        p = p * ( s + rexp_s2(i))
       enddo

       do i=1, nexpcc
        p = p * ( s + rexp_cc(i))
       enddo

     end function denominator_polynomial




     function s11_polynomial(s) result (p)
       implicit none
       double precision, intent(in) :: s
       double precision             :: p

       integer                      :: i

       p = S11_numerator(s) 

       do i=1, nexp1
        p = p * ( s + rexp_s1(i))
       enddo

       do i=1, nexp2
        p = p * ( s + rexp_s2(i))
       enddo

       do i=1, nexpcc
        p = p * ( s + rexp_cc(i))
       enddo

     end function s11_polynomial



     function s22_polynomial(s) result (p)
       implicit none
       double precision, intent(in) :: s
       double precision             :: p

       integer                      :: i

       p = S22_numerator(s) 

       do i=1, nexp1
        p = p * ( s + rexp_s1(i))
       enddo

       do i=1, nexp2
        p = p * ( s + rexp_s2(i))
       enddo

       do i=1, nexpcc
        p = p * ( s + rexp_cc(i))
       enddo

     end function s22_polynomial



     function s12_polynomial(s) result (p)
       implicit none
       double precision, intent(in) :: s
       double precision             :: p

       integer                      :: i

       p = S12_numerator(s) 

       do i=1, nexp1
        p = p * ( s + rexp_s1(i))
       enddo

       do i=1, nexp2
        p = p * ( s + rexp_s2(i))
       enddo

       do i=1, nexpcc
        p = p * ( s + rexp_cc(i))
       enddo

     end function s12_polynomial



    subroutine get_spoints_vector( sx, n )
! ist so noch ungeeignet !
       implicit none
       double precision, intent(out)   :: sx(n)
       integer,          intent(inout) :: n

       double precision                :: sxs(n)
       double precision                :: smin, smax
       integer                         :: perm(n)
       integer   :: mp, ier, j, i
       double precision                :: epsilon = 1d-4

       mp = nexp1 + nexp2 + nexpcc 
       
       if( n < mp+2 ) then
        write(6,*)"too few spoints specified !"
        stop
       endif

       sxs(1             :  nexp1)              = rexp_s1(1:nexp1)
       sxs(nexp1+1       :  nexp1+nexp2)        = rexp_s2(1:nexp2)
       sxs(nexp1+nexp2+1 :  nexp1+nexp2+nexpcc) = rexp_cc(1:nexpcc)

       perm(1:n)    = [(i,i=1,mp)]
       call dpsort(sxs,mp,perm,2,ier)
      
       j=1
       do i=2,mp
         if(abs(sxs(j)-sxs(i)) > epsilon) then
           j=j+1
           sxs(j) = sxs(i)
         endif
       enddo
       mp = j


       smin = sxs(1)
       smax = sxs(mp)
       
       sx(1) = 0
       sx(2) = 0.5d0 * smin
       do i=1,mp-1
         sx(i+2) = (sxs(i)+sxs(i+1))/2
       enddo
       
        sx(mp+2) = 1.5d0 * smax

        n = mp+2
       
     end subroutine get_spoints_vector

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Inverse Laplace coeffcient determinetion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 subroutine get_invlaplace_coeffs( il_coeffs11, il_coeffs12, il_coeffs22 , il_alphas, nnsum) 
!-------------------------------------------------------------------------------------------
  implicit none
  complex(kind=XPREC), intent(out)  :: il_coeffs11(1:)
  complex(kind=XPREC), intent(out)  :: il_coeffs12(1:)
  complex(kind=XPREC), intent(out)  :: il_coeffs22(1:)
  complex(kind=XPREC), intent(out)  :: il_alphas(1:)
  integer            , intent(out)  :: nnsum


  integer, parameter :: mcoeff  = 50   ! --> hier noch besser formulieren
  integer, parameter :: npoints = 500  !      "    "    "     "
  double precision   :: x(npoints), y(npoints) 
  double precision   :: epsilon = 1d-6
  double precision   :: cpoly_den(0:mcoeff)   ! denominator polynomial coefficients
  double precision   :: cpoly_s11(0:mcoeff)   ! s11 polynomial coefficients
  double precision   :: cpoly_s12(0:mcoeff)   ! s12 polynomial coefficients
  double precision   :: cpoly_s22(0:mcoeff)   ! s22 polynomial coefficients
  double precision   :: a(0:mcoeff)

  complex(kind=XPREC) :: csum_d
  complex(kind=XPREC) :: ccoeffs(0:mcoeff)
  complex(kind=XPREC) :: alpha(1:mcoeff)
  complex(kind=XPREC) :: ccs11(mcoeff), ccs12(mcoeff), ccs22(mcoeff)

  integer             :: degree
  complex(kind=XPREC) :: csum_11, csum_12, csum_22

!!  double precision    :: dss = 1/7d0
  double precision    :: r, sfak, dssx, smax
 


!!  double precision    :: rwarn_limit = 1d-12
  logical             :: rwarning

  integer :: i,j,nn, nx, npp, npt

! ------------------  

  rwarning = .false.

  sfak = sfak0

  sfak = sfak / maxval(rexp_s1(1:nexp1))
  sfak = sfak / maxval(rexp_s2(1:nexp2))
  sfak = sfak / maxval(rexp_cc(1:nexpcc))

 

  npp =  nexp1+nexp2+nexpcc - 1-1 
  npt =  min(npp+3+npp_plus,npoints)  !! ?? optimieren !?
  npp =  min(npp,mcoeff)

 ! call get_spoints_vector( x, npt )  ????? 


 smax =  max( maxval(rexp_s1(1:nexp1)),  maxval(rexp_s2(1:nexp2)), maxval(rexp_cc(1:nexpcc)))
 smax =  max( maxval(rexp_s1(1:nexp1-1)),  maxval(rexp_s2(1:nexp2-1)), maxval(rexp_cc(1:nexpcc-1)))
 
 dssx = dss * smax / npt



 ! Write(6,*)"Denominator Polynomial "
 do j=1,npt
   x(j) =  (j-npt/2) * dssx - ss_offset

   y(j) =  denominator_polynomial(x(j)) * sfak
 enddo

 a   = 0
 call polyfit(x, y, npt, a,npp+1, r)
 cpoly_den = a

 if( r > rwarn_limit ) rwarning = .true.

!!pr!! write(6,*)"denpol r = ", r
!!pr!! do i=0,npp+1
!!pr!!   write(6,*) i, a(i), cpoly_den(i)
!!pr!! enddo
!!pr!! do j=1,npt
!!pr!!   write(6,*) j, x(j), y(j), polyeval(x(j),a,npp+1) 
!!pr!! enddo

 
!  Write(6,*)"S11 Polynomial "

 do j=1,npt
   x(j) =  (j-npt/2) * dssx - ss_offset

   y(j) =  S11_polynomial(x(j))  * sfak
 enddo

 a = 0
 call polyfit(x, y, npt, a,npp, r)
 cpoly_s11 = a

 if( r > rwarn_limit ) rwarning = .true.


!!pr!! write(6,*)"s11pol r = ", r
!!pr!! do i=0,npp
!!pr!!   write(6,*) i, a(i), cpoly_s11(i)
!!pr!! enddo
!!pr!! do j=1,npt
!!pr!!   write(6,*) j, x(j), y(j), polyeval(x(j),a,npp) 
!!pr!! enddo


 
 
! Write(6,*)"S22 Polynomial "

 do j=1,npt
   x(j) =  (j-npt/2) * dssx - ss_offset


   y(j) =  S22_polynomial(x(j)) * sfak
 enddo

 a = 0
 call polyfit(x, y, npt, a,npp, r)
 cpoly_s22 = a

 if( r > rwarn_limit ) rwarning = .true.

!!pr!! write(6,*)"s22pol r = ", r
!!pr!! do i=0,npp
!!pr!!   write(6,*) i, a(i), cpoly_s22(i)
!!pr!! enddo
!!pr!! do j=1,npt
!!pr!!   write(6,*) j, x(j), y(j), polyeval(x(j),a,npp) 
!!pr!! enddo



 
 
! Write(6,*)"S12 Polynomial "

 do j=1,npt
   x(j) =  (j-npt/2) * dssx - ss_offset


   y(j) =  S12_polynomial(x(j)) * sfak
 enddo

 a = 0
 call polyfit(x, y, npt, a,npp, r)
 cpoly_s12 = a

 if( r > rwarn_limit ) rwarning = .true.

 
!!pr!! write(6,*)"s12pol r = ", r
!!pr!! do i=0,npp
!!pr!!   write(6,*) i, a(i), cpoly_s12(i)
!!pr!! enddo
!!pr!! do j=1,npt
!!pr!!   write(6,*) j, x(j), y(j), polyeval(x(j),a,npp) 
!!pr!! enddo






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NEXT STEP
!!
!! komplexe Nullstellen alpha_i des Nennerpolynoms (coeffs cpoly_den) berechnen.
!!
 
  ccoeffs = 0 
  ccoeffs = cpoly_den
  degree = npp+1
dl:  do i=degree,1,-1
       if(abs(ccoeffs(i)) > epsilon * sum(abs(ccoeffs(0:degree)))) then
         degree = i
         exit dl
       endif
     enddo dl

  call cmplx_roots_gen(alpha, ccoeffs, degree, .true. , .false. )

! write(6,*)'roots: ',alpha

  call check_degeneracy(degree, alpha)

  ccs11(1:degree) = 0  
  ccs22(1:degree) = 0 
  ccs12(1:degree) = 0 

  do j=1,degree
   csum_11 = 0
   csum_12 = 0
   csum_22 = 0
   csum_d  = 0
   do i = 0,degree 
    csum_11 = csum_11 + cpoly_s11(i) * alpha(j)**i 
    csum_12 = csum_12 + cpoly_s12(i) * alpha(j)**i 
    csum_22 = csum_22 + cpoly_s22(i) * alpha(j)**i 
    csum_d  = csum_d  + i * cpoly_den(i) * alpha(j)**(i-1)
   enddo
   ccs11(j) = csum_11 / csum_d
   ccs12(j) = csum_12 / csum_d
   ccs22(j) = csum_22 / csum_d
  enddo


!!pr!! write(6,*)"cctest1a"
!!pr!! do i=1,degree
!!pr!! write(6,'(i5,8e16.7)')i,ccs11(i), ccs12(i), ccs22(i) , alpha(i) 
!!pr!! enddo


  nnsum = degree

  il_coeffs11(1:degree)  = ccs11(1:degree)
  il_coeffs12(1:degree)  = ccs12(1:degree)
  il_coeffs22(1:degree)  = ccs22(1:degree)
  il_alphas(1:degree)    = alpha(1:degree)



 end subroutine get_invlaplace_coeffs


  
  subroutine compute_invlaplace(t, il_coeffs11, il_coeffs12, il_coeffs22 , il_alphas, nnsum, s11, s12, s22) 
! ---------------------------------------------------------------------------------------------------------
     implicit none  
     double precision   , intent(in)  :: t
     complex(kind=XPREC), intent(in)  :: il_coeffs11(nnsum)
     complex(kind=XPREC), intent(in)  :: il_coeffs12(nnsum)
     complex(kind=XPREC), intent(in)  :: il_coeffs22(nnsum)
     complex(kind=XPREC), intent(in)  :: il_alphas(nnsum)
     integer            , intent(in)  :: nnsum
     double precision   , intent(out) :: s11
     double precision   , intent(out) :: s12
     double precision   , intent(out) :: s22

     complex(kind=XPREC) :: csum_11, csum_12, csum_22
     integer             :: j

     csum_11 = 0
     csum_12 = 0
     csum_22 = 0
     do j=1,nnsum
        csum_11 = csum_11 + il_coeffs11(j) * exp(il_alphas(j)*t)
        csum_12 = csum_12 + il_coeffs12(j) * exp(il_alphas(j)*t)
        csum_22 = csum_22 + il_coeffs22(j) * exp(il_alphas(j)*t)
     enddo

     s11 = RealPart(csum_11)
     s12 = RealPart(csum_12)
     s22 = RealPart(csum_22)
 
  end subroutine  compute_invlaplace

 



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine f04aef_wrap(A, IA, B, IB, N, M, C, IC, WKSPCE )
       implicit real*8 (A-H,O-Z)
       integer ind, itask
       integer iwksp(IA)
       dimension A(IA,IA),C(IA,1), B(IB,1)
       double precision :: WKSPCE(IA)
!      write(*,*)'nag routine f04aef replaced by slatec...'
       itask = 1
       ind = 1
       call dgefs(A,IA,N,B,itask,ind,WKSPCE,iwksp)
       do i=1,N
          C(i,1)=B(i,1)
       enddo

       return
       END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module polynom
