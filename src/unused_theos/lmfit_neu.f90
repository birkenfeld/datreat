


 MODULE lmfit_neu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! matching a relaxation curve given by the table xv(1..np), yv(1..np) with a sum of 
!! nexps simple exponentials with prefactors aexp(1..nexp) and decay rates rexp(1..nexp)
!! the quality of matching is indicated by a small value < 10e-4 of ssq
!! the main intention of this procedure is to get an analytical model for the Laplace
!! transform of any computed relaxation curves from the usual (or unusual) models e.g. in datreat
!! The analytic representation as F(s) = sum(i=1,nexps) aexp(i)/(s+rexp(i))
!! allows combination of curves within the RPA expressions and either "direct" or numerical
!! backtransformation (in the complex regime (s --> i omega)
!!
!! michael monkenbusch, JCNS-1, Forschungszentrum Juelich
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!
!! based on an encapsulation of the MINPACK routines:
!     Argonne National Laboratory. MINPACK Project. November 1996.      
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'            
!
! And using (adapted) SLATEC sources to make this self-contained
!


  abstract interface 
      double precision function f_model(t, px, npar)
         double precision, intent(in)              :: t
         double precision, intent(in)              :: px(npar)
         integer         , intent(in)              :: npar
      end function f_model
   end interface 


  double precision, public :: NEXP_DEVIATION_TOLERANCE = 1d-3
  double precision, public :: NEXP_MINIMUM_RATE        = 1d-6
  double precision, public :: NEXP_TABRANGE_EXTENDER   = 3d0


PRIVATE


public :: fit_simple
public :: f_model

public :: fm
public :: nexp_match
public :: nexp_match2
public :: match_exp
public :: prepare_ttable



CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! matching a relaxation curve given by the table xv(1..np), yv(1..np) with a sum of 
!! ne1 simple exponentials (maximum me) with prefactors aexp(1..ne1) and decay rates rexp(1..ne1)
!! the quality of matching is indicated by a small value of maxdev
!! the main intention of this procedure is to get an analytical model for the Laplace
!! transform of any computed relaxation curves from the usual (or unusual) models e.g. in datreat
!! The analytic representation as F(s) = sum(i=1,nexps) aexp(i)/(s+rexp(i))
!! allows combination of curves within the RPA expressions and either "direct" or numerical
!! backtransformation (in the complex regime (s --> i omega)
!! USE: PREPARE_TTABLE to compute a suitable xv table and the fill with some external S:  yv=S(Q,xv)
!!      then call MATCH_EXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine match_exp ( xv,yv,np,me,ne1,a,r,maxdev) 
 implicit none
 double precision, intent(in) :: xv(:)     ! table with (log-space t values of tabulated S(Q,t)
 double precision, intent(in) :: yv(:)     ! table of coresponding S(Q,t=xv) values
 integer,          intent(in) :: np        ! xv, yv table size 
 integer,          intent(in) :: me        ! max number of exponentials
 integer,          intent(out):: ne1       ! actual number of exps (outcome of this fitting proc)
 double precision, intent(out):: a(me)     ! amplitudes
 double precision, intent(out):: r(me)     ! rates
 double precision, intent(out):: maxdev    ! max deviation


 integer             :: i, ier
 double precision    :: ydiff(np), chisq

 if(me < 3) stop "match_exp: me must be at least 3"

!! TBD new >>>  
ne1 = 3   
a(1:ne1) = 1d0/3d0
r(2)     = 1d0/sum( (yv(1:np-1)-yv(np)) * ( xv(2:np)-xv(1:np-1) ) )
r(1)     = sqrt(r(2)/xv(1))
r(3)     = sqrt(r(2)/(xv(np)/NEXP_TABRANGE_EXTENDER))  

call nexp_match2(xv,yv,np,ne1,a,r,chisq)

! check
do i=1,np
  ydiff(i) = yv(i) - sum(a(1:ne1)*exp(-xv(i)*r(1:ne1)))
!  write(*,'(i5,f16.9,2x,f16.7,4x,f16.12)') i, xv(i), yv(i), ydiff(i)
enddo
maxdev = maxval(abs(ydiff(1:np)))

if(me <=4) return


DO while( 2*ne1 - 1   <= me .and. maxdev  > NEXP_DEVIATION_TOLERANCE)

do i=ne1,1,-1
  a(2*i-1) = a(i)
  r(2*i-1) = r(i)
enddo
ne1 = ne1+(ne1-1)
do i=2,ne1,2
  a(i) = 1d-3
  r(i) = sqrt(r(i-1)*r(i+1))
enddo
call nexp_match2(xv,yv,np,ne1,a,r,chisq)    ! TBD INCLUDE

! check
do i=1,np
  ydiff(i) = yv(i) - sum(a(1:ne1)*exp(-xv(i)*r(1:ne1)))
!  write(*,'(i5,f16.9,2x,f16.7,4x,f16.12)') i, xv(i), yv(i), ydiff(i)
enddo

maxdev = maxval(abs(ydiff(1:np)))


ENDDO


end subroutine match_exp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! USE: PREPARE_TTABLE to compute a suitable xv table and the fill with some external S:  yv=S(Q,xv)
!!      then call MATCH_EXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine prepare_ttable(t0, t1, np, xv)
 implicit none
 double precision, intent(in)   :: t0     ! smallest time
 double precision, intent(in)   :: t1     ! largest time (range)
 integer         , intent(in)   :: np     ! number of points
 double precision, intent(out)  :: xv(np) ! xtable containing taus
 
 ! double precision    :: t_table_spacing1 = NEXP_TABRANGE_EXTENDER 
 ! double precision    :: t_table_spacing2 = 1d0 
 integer             :: i
 do i=1,np
  xv(i)  =  exp(i*log(t1*NEXP_TABRANGE_EXTENDER/t0)/np) * t0    !! TBD new here is the table spacing
 enddo
     
end subroutine prepare_ttable





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! matching a relaxation curve given by the table xv(1..np), yv(1..np) with a sum of 
!! nexps simple exponentials with prefactors aexp(1..nexp) and decay rates rexp(1..nexp)
!! the quality of matching is indicated by a small value < 10e-4 of ssq
!! the main intention of this procedure is to get an analytical model for the Laplace
!! transform of any computed relaxation curves from the usual (or unusual) models e.g. in datreat
!! The analytic representation as F(s) = sum(i=1,nexps) aexp(i)/(s+rexp(i))
!! allows combination of curves within the RPA expressions and either "direct" or numerical
!! backtransformation (in the complex regime (s --> i omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     subroutine nexp_match2(xv, yv, np, nexps, aexp, rexp, ssq, iout) 
!    ---------------------------------------------------------------
     implicit none
     double precision, intent(in)     :: xv(np)   ! x-values (i.e. time) of cuve to be matched
     double precision, intent(in)     :: yv(np)   ! y-values (i.e. S(q,t)/S(q))
     integer         , intent(in)     :: np       ! number of values
     integer         , intent(inout)  :: nexps    ! number of exponentials to be fitted
     double precision, intent(inout)    :: aexp(nexps) ! exp prefactors (ampltiues)
     double precision, intent(inout)    :: rexp(nexps) ! exp rates:  y = sum(1...nexps) aexp(i)*exp(-t*rexp(i))
     double precision, intent(out)    :: ssq         ! quality of matching
     integer         , intent(in), optional :: iout  ! output levevl


     double precision    :: px(2*nexps) , pxscale(2*nexps), pxerr(2*nexps)
     integer             :: iperm(nexps)
     double precision    :: a0(nexps), r0(nexps)

     double precision    :: ye(np)
     double precision    :: tmax, t
     integer             :: i, j, ier, nxs

     logical             :: verbose = .false.

     if(present(iout)) then
       verbose = (iout > 0)
     endif


      ! make estimate
     
      px      = 0d0
      pxscale = 1d0
      ye      = 1d0
      ssq     = 0d0

      
      tmax = maxval(xv)
      do i=1,nexps
        if(aexp(i) .eq. 0d0) then
           px(2*i-1)       = 1d0/nexps
           pxscale(2*i-1)  = 0.1d0
           px(2*i)         = exp(-i*log(tmax)/nexps)
           pxscale(2*i)    = px(2*i)
        else
           px(2*i-1)       = aexp(i)
           pxscale(2*i-1)  = 0.1d0
           px(2*i)         = rexp(i)
           pxscale(2*i)    = px(2*i)
        endif
      enddo
      
      ssq = fit_simple(fm,2*nexps,px ,pxscale     , np, xv, yv, ye, px, pxerr)
! combine similar rates
  


  
!1      if(ssq > 1d-4) then
!1        write(6,*)"WARNING: bad matching of exp-model to time function!",ssq
!1      endif
     
      do i=1,nexps
         a0(i) =  px(2*i-1)
         r0(i) =  px(2*i)
!!         write(6,'(i8,2f18.7,6x,2f18.7)')i, a0(i), pxerr(2*i-1), r0(i), pxerr(2*i)
         if( r0(i) < NEXP_MINIMUM_RATE ) then
            if(verbose) write(6,*)"WARNING(nexp_match): rate = ",r0(i)," is set to ",NEXP_MINIMUM_RATE 
            r0(i)   = NEXP_MINIMUM_RATE 
            px(2*i) = r0(i)
         endif 
      enddo

!3      if(verbose) then 
!3        write(6,'(a,i4,a,2f12.6,a,2f12.6,a)')  &
!3           "# exp-fit n, range:",np,"(", xv(1), yv(1),")-->(", xv(np), yv(np),")" 
!3        write(6,'(a,i2,a,e9.2,a,20(a,f6.3," t=",f13.2,"|"))')&
!3           "# ",nexps," exp model(",ssq,"):|",("a=",px(2*i-1),1d0/px(2*i),i=1,nexps)
!3      endif

!2      do i=1,np
!2        ye(i) = yv(i)-fm(xv(i),px,2*nexps) 
!2      enddo
!2      if(maxval(ye(1:np)) > 0.01d0) then
!2         write(*,'("Maximum deviation (model-n-exp-rpresentation): ",f12.6," at t= ",f12.4)') &
!2                 maxval(ye(1:np)), xv(maxloc(ye(1:np)))
!2      endif

     
    ! sorting fastes rate first
      iperm = [(i,i=1,nexps)]
      call  DPSORT (r0, nexps, IPERM,  -1 , ier)
     
      aexp = a0(iperm(1:nexps))
      rexp = r0(iperm(1:nexps))
     
     
     end subroutine nexp_match2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! matching a relaxation curve given by the table xv(1..np), yv(1..np) with a sum of 
!! nexps simple exponentials with prefactors aexp(1..nexp) and decay rates rexp(1..nexp)
!! the quality of matching is indicated by a small value < 10e-4 of ssq
!! the main intention of this procedure is to get an analytical model for the Laplace
!! transform of any computed relaxation curves from the usual (or unusual) models e.g. in datreat
!! The analytic representation as F(s) = sum(i=1,nexps) aexp(i)/(s+rexp(i))
!! allows combination of curves within the RPA expressions and either "direct" or numerical
!! backtransformation (in the complex regime (s --> i omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     subroutine nexp_match(xv, yv, np, nexps, aexp, rexp, ssq, iout) 
!    ---------------------------------------------------------------
     implicit none
     double precision, intent(in)     :: xv(np)   ! x-values (i.e. time) of cuve to be matched
     double precision, intent(in)     :: yv(np)   ! y-values (i.e. S(q,t)/S(q))
     integer         , intent(in)     :: np       ! number of values
     integer         , intent(in)     :: nexps    ! number of exponentials to be fitted
     double precision, intent(out)    :: aexp(nexps) ! exp prefactors (ampltiues)
     double precision, intent(out)    :: rexp(nexps) ! exp rates:  y = sum(1...nexps) aexp(i)*exp(-t*rexp(i))
     double precision, intent(out)    :: ssq         ! quality of matching
     integer         , intent(in), optional :: iout  ! output levevl


     double precision    :: px(2*nexps) , pxscale(2*nexps), pxerr(2*nexps)
     integer             :: iperm(nexps)
     double precision    :: a0(nexps), r0(nexps)

     double precision    :: ye(np)
     double precision    :: tmax, t
     integer             :: i, ier, nxs

     logical             :: verbose = .false.

     if(present(iout)) then
       verbose = (iout > 0)
     endif


      ! make estimate
     
      px      = 0d0
      pxscale = 1d0
      ye      = 1d0
      ssq     = 0d0

      
      tmax = maxval(xv)
      do i=1,nexps
        px(2*i-1)       = 1d0/nexps
        pxscale(2*i-1)  = 0.1d0
        px(2*i)         = exp(-i*log(tmax)/nexps)
        pxscale(2*i)    = px(2*i)
      enddo
      
      ssq = fit_simple(fm,2*nexps,px ,pxscale     , np, xv, yv, ye, px, pxerr)
!      ssq = fit_simple(fm,2*nexps,px ,pxscale     , np, xv, yv, ye, px, pxerr)
  


  
      if(ssq > 1d-4) then
        write(6,*)"WARNING: bad matching of exp-model to time function!",ssq
      endif
     
      do i=1,nexps
         a0(i) =  px(2*i-1)
         r0(i) =  px(2*i)
!!         write(6,'(i8,2f18.7,6x,2f18.7)')i, a0(i), pxerr(2*i-1), r0(i), pxerr(2*i)
         if( r0(i) < NEXP_MINIMUM_RATE ) then
            if(verbose) write(6,*)"WARNING(nexp_match): rate = ",r0(i)," is set to ",NEXP_MINIMUM_RATE 
            r0(i)   = NEXP_MINIMUM_RATE 
            px(2*i) = r0(i)
         endif 
      enddo

      if(verbose) then 
        write(6,'(a,i4,a,2f12.6,a,2f12.6,a)')  &
           "# exp-fit n, range:",np,"(", xv(1), yv(1),")-->(", xv(np), yv(np),")" 
        write(6,'(a,i2,a,e9.2,a,20(a,f6.3," t=",f13.2,"|"))')&
           "# ",nexps," exp model(",ssq,"):|",("a=",px(2*i-1),1d0/px(2*i),i=1,nexps)
      endif

      do i=1,np
        ye(i) = yv(i)-fm(xv(i),px,2*nexps) 
      enddo
      if(maxval(ye(1:np)) > 0.01d0) then
         write(*,'("Maximum deviation (model-n-exp-rpresentation): ",f12.6," at t= ",f12.4)') &
                 maxval(ye(1:np)), xv(maxloc(ye(1:np)))
      endif

     
    ! sorting fastes rate first
      iperm = [(i,i=1,nexps)]
      call  DPSORT (r0, nexps, IPERM,  -1 , ier)
     
      aexp = a0(iperm)
      rexp = r0(iperm)
     
     
     end subroutine nexp_match


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!! being a necessary ingredient for the nexp_match function his may serve as a general means
!! for nonlinear curve fitting (even recursive use may be considered)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     recursive function fit_simple( fmodel, npar, parin, scalw, np, xval, yval, yerr, par, parerr) result(ssq)
!    --------------------------------------------------------------------------------------------------------
     implicit none
     procedure(f_model)               :: fmodel      ! model function (see abstact intefrace def above=
     integer, intent(in)              :: npar        ! number of model parameters 
     double precision, intent(in)     :: parin(npar) ! parameter vector (start) 
     double precision, intent(in)     :: scalw(npar) ! scale parameter (0=fixed=

     integer, intent(in)              :: np        ! number of data points      
     double precision, intent(in)     :: xval(np)  ! x values of data
     double precision, intent(in)     :: yval(np)  ! y values of data
     double precision, intent(in)     :: yerr(np)  ! error values of data
     double precision, intent(out)    :: par(npar) ! resulting optimized parameters
     double precision, intent(out)    :: parerr(npar) ! estimated parameter errors

     double precision                 :: ssq

     double precision                 :: xstart(npar)
     double precision                 :: pa_new(npar)
!     double precision                 :: fvec(np)
     double precision                 :: fvec(np)
     double precision                 :: diag(npar) ! = 1d0

     double precision                 :: xjac(np,npar)
     integer                          :: ipvt(npar)
     double precision                 :: qtf(npar)

     double precision                 :: ftol = 1d-6
     double precision                 :: xtol = 1d-6
     double precision                 :: gtol = 0
     integer                          :: maxfev = 3000
     double precision                 :: epsfcn = 1d-6
     integer                          :: mode   = 2
     double precision                 :: factor = 1d0
     integer                          :: nprint = 0
     integer                          :: info
     integer                          :: nfev

     integer                          :: i,j,l,ll, n

     double precision                 :: ermat(npar,np), sume
     double precision, allocatable    :: gmat(:,:), ginv(:,:), gmul(:,:)

     double precision                 :: wa1(npar), wa2(npar), wa3(npar)
     double precision                 :: wa4(np)
 

     
     wa1 = 0
     wa2 = 0
     wa3 = 0
     wa4 = 0



     !
     ! prepare startvector
     !
     par    = parin
     n      = 0
     xstart = 0
     do i=1,npar
       if(scalw(i) .ne. 0d0 ) then
         n = n+1
         xstart(n) = par(i) / scalw(i)
       endif
     enddo
     !
     if( n == 0 ) then
       write(6,*)"ERROR: fit_simple, no parmeters to fit:",scalw
     endif

     diag = 1d0

 !    write(6,*)"no parameters:",n
 !    write(6,*)"parin :",parin
 !    write(6,*)"xstart:",xstart
    
 !    write(6,*)"check 1"

     nprint = 1


     call lmdif(fcn,np,n,xstart,fvec,ftol,xtol,gtol,maxfev,epsfcn,     &
     &                 diag,mode,factor,nprint,info,nfev,xjac,np,   &
     &                 ipvt,qtf, wa1, wa2, wa3, wa4 )
  

 !    write(6,*)"check 2"

     call fdjac2(fcn,np,n,xstart,fvec,xjac,np,2,epsfcn,wa4)

     par = pa_new
     

!!!================ ERROR DETERMINATION ==================================
     if(allocated(gmat)) deallocate(gmat)
                         allocate(gmat(n,n))
     if(allocated(ginv)) deallocate(ginv)
                         allocate(ginv(n,n))
     if(allocated(gmul)) deallocate(gmul)
                         allocate(gmul(n,n))
     gmat = 0
     ginv = 0
     gmul = 0
! ----- error-determination ?? test ?? -----------------------------
!
!!         write(6,*)'error determination...'
         do i=1,n
          do j=1,n
           sume = 0
           do l=1,np
            sume = sume + xjac(l,i)*xjac(l,j)
           enddo
           gmat(i,j) = sume
          enddo
         enddo
! --- invertiere gmat ---
         call matrix_inverse(gmat,ginv,n)
 !!        write(6,*)'gmat'
 !!        do i=1,n ; write(6,'(4f12.6)') gmat(i,:) ; enddo
 !!        write(6,*)'ginv'
 !!        do i=1,n ; write(6,'(4f12.6)') ginv(i,:) ; enddo
 !!        write(6,*) 'mul'
 !!        gmul =  matmul(gmat,ginv)
 !!        do i=1,n ; write(6,'(4f12.6)') gmul(i,:) ; enddo

! --- bilde g**-1 f  (f=xjac) ---
         do i=1,np
          do l=1,n
           sume = 0
           do j=1,n
            sume = sume + ginv(l,j)*xjac(i,j)
           enddo
           ermat(l,i) = sume
          enddo
         enddo
! --- bilde fehlerquadratmatrix ---
         do l=1,n
          do ll=1,n
           sume = 0
           do i=1,np
            sume = sume + ermat(l,i)*ermat(ll,i)
           enddo
           gmat(l,ll) = sume
          enddo
         enddo

         j = 0
         do i=1,npar
           if(scalw(i) .ne. 0) then
             j = j+1
             parerr(i) = sqrt(gmat(j,j)) * scalw(i)
           else
             parerr(i) = 0
           endif
         enddo

!     write(6,'("ssq, par=",f12.6,3x,4(f16.7,"+-",f12.6,2x))') ssq, (par(i),parerr(i),i=1,npar)



     CONTAINS


       subroutine fcn(m,n,x,fvec)!  - unused (PAZ) ,iflag) 
         implicit none
         integer, intent(in)             :: m
         integer, intent(in)             :: n
         double precision, intent(in)    :: x(n)
         double precision, intent(out)   :: fvec(m)
         !double precision, intent(inout),optional :: iflag
         
!         double precision              :: fmodel
         integer                       :: i, j



 !       write(6,*)"fcn: ",m,n,x

       ! unscale parameters !
         pa_new = par
         j      = 0
         do i=1,npar
            if(scalw(i) .ne. 0d0 ) then
              j = j+1
              pa_new(i) = x(j) * scalw(i)
            endif
         enddo

    
       ! create scale errors
         do i=1,m
           if(yerr(i) > 0d0) then
              fvec(i) = (fmodel(xval(i),pa_new,npar) - yval(i)) / yerr(i)
           else
              fvec(i) = (fmodel(xval(i),pa_new,npar) - yval(i))
           endif  
         enddo
        
         ssq = 0
         do i=1,m
           ssq = ssq + fvec(i)**2
         enddo 
         ssq = ssq / m
 
       end subroutine fcn

   end function fit_simple







      recursive &                                  
      subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,     &
     &                 diag,mode,factor,nprint,info,nfev,fjac,ldfjac,   &
     &                 ipvt,qtf,wa1,wa2,wa3,wa4)                        
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac 
      integer ipvt(n) 
      double precision ftol,xtol,gtol,epsfcn,factor 
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),      &
     &                 wa1(n),wa2(n),wa3(n),wa4(m)                      
      external fcn 
!     **********                                                        
!                                                                       
!     subroutine lmdif                                                  
!                                                                       
!     the purpose of lmdif is to minimize the sum of the squares of     
!     m nonlinear functions in n variables by a modification of         
!     the levenberg-marquardt algorithm. the user must provide a        
!     subroutine which calculates the functions. the jacobian is        
!     then calculated by a forward-difference approximation.            
!                                                                       
!     the subroutine statement is                                       
!                                                                       
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,   
!                        diag,mode,factor,nprint,info,nfev,fjac,        
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)               
!                                                                       
!     where                                                             
!                                                                       
!       fcn is the name of the user-supplied subroutine which           
!         calculates the functions. fcn must be declared                
!         in an external statement in the user calling                  
!         program, and should be written as follows.                    
!                                                                       
!         subroutine fcn(m,n,x,fvec,iflag)                              
!         integer m,n,iflag                                             
!         double precision x(n),fvec(m)                                 
!         ----------                                                    
!         calculate the functions at x and                              
!         return this vector in fvec.                                   
!         ----------                                                    
!         return                                                        
!         end                                                           
!                                                                       
!         the value of iflag should not be changed by fcn unless        
!         the user wants to terminate execution of lmdif.               
!         in this case set iflag to a negative integer.                 
!                                                                       
!       m is a positive integer input variable set to the number        
!         of functions.                                                 
!                                                                       
!       n is a positive integer input variable set to the number        
!         of variables. n must not exceed m.                            
!                                                                       
!       x is an array of length n. on input x must contain              
!         an initial estimate of the solution vector. on output x       
!         contains the final estimate of the solution vector.           
!                                                                       
!       fvec is an output array of length m which contains              
!         the functions evaluated at the output x.                      
!                                                                       
!       ftol is a nonnegative input variable. termination               
!         occurs when both the actual and predicted relative            
!         reductions in the sum of squares are at most ftol.            
!         therefore, ftol measures the relative error desired           
!         in the sum of squares.                                        
!                                                                       
!       xtol is a nonnegative input variable. termination               
!         occurs when the relative error between two consecutive        
!         iterates is at most xtol. therefore, xtol measures the        
!         relative error desired in the approximate solution.           
!                                                                       
!       gtol is a nonnegative input variable. termination               
!         occurs when the cosine of the angle between fvec and          
!         any column of the jacobian is at most gtol in absolute        
!         value. therefore, gtol measures the orthogonality             
!         desired between the function vector and the columns           
!         of the jacobian.                                              
!                                                                       
!       maxfev is a positive integer input variable. termination        
!         occurs when the number of calls to fcn is at least            
!         maxfev by the end of an iteration.                            
!                                                                       
!       epsfcn is an input variable used in determining a suitable      
!         step length for the forward-difference approximation. this    
!         approximation assumes that the relative errors in the         
!         functions are of the order of epsfcn. if epsfcn is less       
!         than the machine precision, it is assumed that the relative   
!         errors in the functions are of the order of the machine       
!         precision.                                                    
!                                                                       
!       diag is an array of length n. if mode = 1 (see                  
!         below), diag is internally set. if mode = 2, diag             
!         must contain positive entries that serve as                   
!         multiplicative scale factors for the variables.               
!                                                                       
!       mode is an integer input variable. if mode = 1, the             
!         variables will be scaled internally. if mode = 2,             
!         the scaling is specified by the input diag. other             
!         values of mode are equivalent to mode = 1.                    
!                                                                       
!       factor is a positive input variable used in determining the     
!         initial step bound. this bound is set to the product of       
!         factor and the euclidean norm of diag*x if nonzero, or else   
!         to factor itself. in most cases factor should lie in the      
!         interval (.1,100.). 100. is a generally recommended value.    
!                                                                       
!       nprint is an integer input variable that enables controlled     
!         printing of iterates if it is positive. in this case,         
!         fcn is called with iflag = 0 at the beginning of the first    
!         iteration and every nprint iterations thereafter and          
!         immediately prior to return, with x and fvec available        
!         for printing. if nprint is not positive, no special calls     
!         of fcn with iflag = 0 are made.                               
!                                                                       
!       info is an integer output variable. if the user has             
!         terminated execution, info is set to the (negative)           
!         value of iflag. see description of fcn. otherwise,            
!         info is set as follows.                                       
!                                                                       
!         info = 0  improper input parameters.                          
!                                                                       
!         info = 1  both actual and predicted relative reductions       
!                   in the sum of squares are at most ftol.             
!                                                                       
!         info = 2  relative error between two consecutive iterates     
!                   is at most xtol.                                    
!                                                                       
!         info = 3  conditions for info = 1 and info = 2 both hold.     
!                                                                       
!         info = 4  the cosine of the angle between fvec and any        
!                   column of the jacobian is at most gtol in           
!                   absolute value.                                     
!                                                                       
!         info = 5  number of calls to fcn has reached or               
!                   exceeded maxfev.                                    
!                                                                       
!         info = 6  ftol is too small. no further reduction in          
!                   the sum of squares is possible.                     
!                                                                       
!         info = 7  xtol is too small. no further improvement in        
!                   the approximate solution x is possible.             
!                                                                       
!         info = 8  gtol is too small. fvec is orthogonal to the        
!                   columns of the jacobian to machine precision.       
!                                                                       
!       nfev is an integer output variable set to the number of         
!         calls to fcn.                                                 
!                                                                       
!       fjac is an output m by n array. the upper n by n submatrix      
!         of fjac contains an upper triangular matrix r with            
!         diagonal elements of nonincreasing magnitude such that        
!                                                                       
!                t     t           t                                    
!               p *(jac *jac)*p = r *r,                                 
!                                                                       
!         where p is a permutation matrix and jac is the final          
!         calculated jacobian. column j of p is column ipvt(j)          
!         (see below) of the identity matrix. the lower trapezoidal     
!         part of fjac contains information generated during            
!         the computation of r.                                         
!                                                                       
!       ldfjac is a positive integer input variable not less than m     
!         which specifies the leading dimension of the array fjac.      
!                                                                       
!       ipvt is an integer output array of length n. ipvt               
!         defines a permutation matrix p such that jac*p = q*r,         
!         where jac is the final calculated jacobian, q is              
!         orthogonal (not stored), and r is upper triangular            
!         with diagonal elements of nonincreasing magnitude.            
!         column j of p is column ipvt(j) of the identity matrix.       
!                                                                       
!       qtf is an output array of length n which contains               
!         the first n elements of the vector (q transpose)*fvec.        
!                                                                       
!       wa1, wa2, and wa3 are work arrays of length n.                  
!                                                                       
!       wa4 is a work array of length m.                                
!                                                                       
!     subprograms called                                                
!                                                                       
!       user-supplied ...... fcn                                        
!                                                                       
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac            
!                                                                       
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod                 
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i,iflag,iter,j,l 
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,   &
     &                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,  &
     &                 sum,temp,temp1,temp2,xnorm,zero                  
      double precision xdpmpar,xenorm 
      data one,p1,p5,p25,p75,p0001,zero                                 &
     &     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/             
!                                                                       
!     epsmch is the machine precision.                                  
!         
!!      write(6,*)"entered lmdif 1"
                                                               
      epsmch = dpmpar(1) 
!                                                                       
      info = 0 
      iflag = 0 
      nfev = 0 
!                                                                       
!     check the input parameters for errors.                            
!                                                                       
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m                     &
     &    .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero   &
     &    .or. maxfev .le. 0 .or. factor .le. zero) go to 300           
      if (mode .ne. 2) go to 20 
      do 10 j = 1, n 
         if (diag(j) .le. zero) go to 300 
   10    continue 
   20 continue 
!                                                                       
!     evaluate the function at the starting point                       
!     and calculate its norm.                                           
!                                                                       
      iflag = 1 

      call fcn(m,n,x,fvec,iflag) 


      nfev = 1 
      if (iflag .lt. 0) go to 300 
      fnorm = enorm(m,fvec) 
!                                                                       
!     initialize levenberg-marquardt parameter and iteration counter.   
!                                                                       
      par = zero 
      iter = 1 
!                                                                       
!     beginning of the outer loop.

                                       
!                                                                       
   30 continue 
!                                                                       
!        calculate the jacobian matrix.                                 
!                                                                       
         iflag = 2 
         call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4) 
         nfev = nfev + n 
         if (iflag .lt. 0) go to 300 
!                                                                       
!        if requested, call fcn to enable printing of iterates.         
!                                                                       
         if (nprint .le. 0) go to 40 
         iflag = 0 
         if (mod(iter-1,nprint) .eq. 0) call fcn(m,n,x,fvec,iflag) 
         if (iflag .lt. 0) go to 300 
   40    continue 
!                                                                       
!        compute the qr factorization of the jacobian.                  
!                                                                       
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3) 
!                                                                       
!        on the first iteration and if mode is 1, scale according       
!        to the norms of the columns of the initial jacobian.           
!                                                                       
         if (iter .ne. 1) go to 80 
         if (mode .eq. 2) go to 60 
         do 50 j = 1, n 
            diag(j) = wa2(j) 
            if (wa2(j) .eq. zero) diag(j) = one 
   50       continue 
   60    continue 
!                                                                       
!        on the first iteration, calculate the norm of the scaled x     
!        and initialize the step bound delta.                           
!                                                                       
         do 70 j = 1, n 
            wa3(j) = diag(j)*x(j) 
   70       continue 
         xnorm = enorm(n,wa3) 
         delta = factor*xnorm 
         if (delta .eq. zero) delta = factor 
   80    continue 
!                                                                       
!        form (q transpose)*fvec and store the first n components in    
!        qtf.                                                           
!                                                                       
         do 90 i = 1, m 
            wa4(i) = fvec(i) 
   90       continue 
         do 130 j = 1, n 
            if (fjac(j,j) .eq. zero) go to 120 
            sum = zero 
            do 100 i = j, m 
               sum = sum + fjac(i,j)*wa4(i) 
  100          continue 
            temp = -sum/fjac(j,j) 
            do 110 i = j, m 
               wa4(i) = wa4(i) + fjac(i,j)*temp 
  110          continue 
  120       continue 
            fjac(j,j) = wa1(j) 
            qtf(j) = wa4(j) 
  130       continue 
!                                                                       
!        compute the norm of the scaled gradient.                       
!                                                                       
         gnorm = zero 
         if (fnorm .eq. zero) go to 170 
         do 160 j = 1, n 
            l = ipvt(j) 
            if (wa2(l) .eq. zero) go to 150 
            sum = zero 
            do 140 i = 1, j 
               sum = sum + fjac(i,j)*(qtf(i)/fnorm) 
  140          continue 
            gnorm = dmax1(gnorm,dabs(sum/wa2(l))) 
  150       continue 
  160       continue 
  170    continue 
!                                                                       
!        test for convergence of the gradient norm.                     
!                                                                       
         if (gnorm .le. gtol) info = 4 
         if (info .ne. 0) go to 300 
!                                                                       
!        rescale if necessary.                                          
!                                                                       
         if (mode .eq. 2) go to 190 
         do 180 j = 1, n 
            diag(j) = dmax1(diag(j),wa2(j)) 
  180       continue 
  190    continue 
!                                                                       
!        beginning of the inner loop.                                   
!                                                                       
  200    continue 
!                                                                       
!           determine the levenberg-marquardt parameter.                
!                                                                       
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,   &
     &                 wa3,wa4)                                         
!                                                                       
!           store the direction p and x + p. calculate the norm of p.   
!                                                                       
            do 210 j = 1, n 
               wa1(j) = -wa1(j) 
               wa2(j) = x(j) + wa1(j) 
               wa3(j) = diag(j)*wa1(j) 
  210          continue 
            pnorm = enorm(n,wa3) 
!                                                                       
!           on the first iteration, adjust the initial step bound.      
!                                                                       
            if (iter .eq. 1) delta = dmin1(delta,pnorm) 
!                                                                       
!           evaluate the function at x + p and calculate its norm.      
!                                                                       
            iflag = 1 
            call fcn(m,n,wa2,wa4,iflag) 
            nfev = nfev + 1 
            if (iflag .lt. 0) go to 300 
            fnorm1 = enorm(m,wa4) 
!                                                                       
!           compute the scaled actual reduction.                        
!                                                                       
            actred = -one 
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2 
!                                                                       
!           compute the scaled predicted reduction and                  
!           the scaled directional derivative.                          
!                                                                       
            do 230 j = 1, n 
               wa3(j) = zero 
               l = ipvt(j) 
               temp = wa1(l) 
               do 220 i = 1, j 
                  wa3(i) = wa3(i) + fjac(i,j)*temp 
  220             continue 
  230          continue 
            temp1 = enorm(n,wa3)/fnorm 
            temp2 = (dsqrt(par)*pnorm)/fnorm 
            prered = temp1**2 + temp2**2/p5 
            dirder = -(temp1**2 + temp2**2) 
!                                                                       
!           compute the ratio of the actual to the predicted            
!           reduction.                                                  
!                                                                       
            ratio = zero 
            if (prered .ne. zero) ratio = actred/prered 
!                                                                       
!           update the step bound.                                      
!                                                                       
            if (ratio .gt. p25) go to 240 
               if (actred .ge. zero) temp = p5 
               if (actred .lt. zero)                                    &
     &            temp = p5*dirder/(dirder + p5*actred)                 
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1 
               delta = temp*dmin1(delta,pnorm/p1) 
               par = par/temp 
               go to 260 
  240       continue 
               if (par .ne. zero .and. ratio .lt. p75) go to 250 
               delta = pnorm/p5 
               par = p5*par 
  250          continue 
  260       continue 
!                                                                       
!           test for successful iteration.                              
!                                                                       
            if (ratio .lt. p0001) go to 290 
!                                                                       
!           successful iteration. update x, fvec, and their norms.      
!                                                                       
            do 270 j = 1, n 
               x(j) = wa2(j) 
               wa2(j) = diag(j)*x(j) 
  270          continue 
            do 280 i = 1, m 
               fvec(i) = wa4(i) 
  280          continue 
            xnorm = enorm(n,wa2) 
            fnorm = fnorm1 
            iter = iter + 1 
  290       continue 
!                                                                       
!           tests for convergence.                                      
!                                                                       
            if (dabs(actred) .le. ftol .and. prered .le. ftol           &
     &          .and. p5*ratio .le. one) info = 1                       
            if (delta .le. xtol*xnorm) info = 2 
            if (dabs(actred) .le. ftol .and. prered .le. ftol           &
     &          .and. p5*ratio .le. one .and. info .eq. 2) info = 3     
            if (info .ne. 0) go to 300 
!                                                                       
!           tests for termination and stringent tolerances.             
!                                                                       
            if (nfev .ge. maxfev) info = 5 
            if (dabs(actred) .le. epsmch .and. prered .le. epsmch       &
     &          .and. p5*ratio .le. one) info = 6                       
            if (delta .le. epsmch*xnorm) info = 7 
            if (gnorm .le. epsmch) info = 8 
            if (info .ne. 0) go to 300 
!                                                                       
!           end of the inner loop. repeat if iteration unsuccessful.    
!                                                                       
            if (ratio .lt. p0001) go to 200 
!                                                                       
!        end of the outer loop.                                         
!                                                                       
         go to 30 
  300 continue 
!                                                                       
!     termination, either normal or user imposed.                       
!                                                                       
      if (iflag .lt. 0) info = iflag 
      iflag = 0 
      if (nprint .gt. 0) call fcn(m,n,x,fvec,iflag) 
      return 
!                                                                       
!     last card of subroutine lmdif.                                    
!                                                                       
      END  

   


    
      double precision function dpmpar(i) 
      integer i 
!     **********                                                        
!                                                                       
!     Function dpmpar                                                   
!                                                                       
!     This function provides double precision machine parameters        
!     when the appropriate set of data statements is activated (by      
!     removing the c from column 1) and all other data statements are   
!     rendered inactive. Most of the parameter values were obtained     
!     from the corresponding Bell Laboratories Port Library function.   
!                                                                       
!     The function statement is                                         
!                                                                       
!       double precision function dpmpar(i)                             
!                                                                       
!     where                                                             
!                                                                       
!       i is an integer input variable set to 1, 2, or 3 which          
!         selects the desired machine parameter. If the machine has     
!         t base b digits and its smallest and largest exponents are    
!         emin and emax, respectively, then these parameters are        
!                                                                       
!         dpmpar(1) = b**(1 - t), the machine precision,                
!                                                                       
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,            
!                                                                       
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.     
!                                                                       
!     Argonne National Laboratory. MINPACK Project. November 1996.      
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'            
!                                                                       
!     **********                                                        
      integer mcheps(4) 
      integer minmag(4) 
      integer maxmag(4) 
      double precision dmach(3) 
      equivalence (dmach(1),mcheps(1)) 
      equivalence (dmach(2),minmag(1)) 
      equivalence (dmach(3),maxmag(1)) 
!                                                                       
      data dmach(1) /2.22044604926d-16/ 
      data dmach(2) /2.22507385852d-308/ 
      data dmach(3) /1.79769313485d+308/ 
!                                                                       
      dpmpar = dmach(i) 
      return 
!                                                                       
!     Last card of function dpmpar.                                     
!                                                                       
      END 


                                          
      double precision recursive function enorm(n,x) 
      integer n 
      double precision x(n) 
!     **********                                                        
!                                                                       
!     function enorm                                                    
!                                                                       
!     given an n-vector x, this function calculates the                 
!     euclidean norm of x.                                              
!                                                                       
!     the euclidean norm is computed by accumulating the sum of         
!     squares in three different sums. the sums of squares for the      
!     small and large components are scaled so that no overflows        
!     occur. non-destructive underflows are permitted. underflows       
!     and overflows do not occur in the computation of the unscaled     
!     sum of squares for the intermediate components.                   
!     the definitions of small, intermediate and large components       
!     depend on two constants, rdwarf and rgiant. the main              
!     restrictions on these constants are that rdwarf**2 not            
!     underflow and rgiant**2 not overflow. the constants               
!     given here are suitable for every known computer.                 
!                                                                       
!     the function statement is                                         
!                                                                       
!       double precision function enorm(n,x)                            
!                                                                       
!     where                                                             
!                                                                       
!       n is a positive integer input variable.                         
!                                                                       
!       x is an input array of length n.                                
!                                                                       
!     subprograms called                                                
!                                                                       
!       fortran-supplied ... dabs,dsqrt                                 
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i 
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,   &
     &                 x1max,x3max,zero                                 
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/ 
      s1 = zero 
      s2 = zero 
      s3 = zero 
      x1max = zero 
      x3max = zero 
      floatn = n 
      agiant = rgiant/floatn 
      do 90 i = 1, n 
         xabs = dabs(x(i)) 
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70 
            if (xabs .le. rdwarf) go to 30 
!                                                                       
!              sum for large components.                                
!                                                                       
               if (xabs .le. x1max) go to 10 
                  s1 = one + s1*(x1max/xabs)**2 
                  x1max = xabs 
                  go to 20 
   10          continue 
                  s1 = s1 + (xabs/x1max)**2 
   20          continue 
               go to 60 
   30       continue 
!                                                                       
!              sum for small components.                                
!                                                                       
               if (xabs .le. x3max) go to 40 
                  s3 = one + s3*(x3max/xabs)**2 
                  x3max = xabs 
                  go to 50 
   40          continue 
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2 
   50          continue 
   60       continue 
            go to 80 
   70    continue 
!                                                                       
!           sum for intermediate components.                            
!                                                                       
            s2 = s2 + xabs**2 
   80    continue 
   90    continue 
!                                                                       
!     calculation of norm.                                              
!                                                                       
      if (s1 .eq. zero) go to 100 
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max) 
         go to 130 
  100 continue 
         if (s2 .eq. zero) go to 110 
            if (s2 .ge. x3max)                                          &
     &         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))            
            if (s2 .lt. x3max)                                          &
     &         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))             
            go to 120 
  110    continue 
            enorm = x3max*dsqrt(s3) 
  120    continue 
  130 continue 
      return 
!                                                                       
!     last card of function enorm.                                      
!                                                                       
      END 


      recursive &                                    
      subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa) 
      integer m,n,ldfjac,iflag 
      double precision epsfcn 
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(m) 
!     **********                                                        
!                                                                       
!     subroutine fdjac2                                                 
!                                                                       
!     this subroutine computes a forward-difference approximation       
!     to the m by n jacobian matrix associated with a specified         
!     problem of m functions in n variables.                            
!                                                                       
!     the subroutine statement is                                       
!                                                                       
!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)   
!                                                                       
!     where                                                             
!                                                                       
!       fcn is the name of the user-supplied subroutine which           
!         calculates the functions. fcn must be declared                
!         in an external statement in the user calling                  
!         program, and should be written as follows.                    
!                                                                       
!         subroutine fcn(m,n,x,fvec,iflag)                              
!         integer m,n,iflag                                             
!         double precision x(n),fvec(m)                                 
!         ----------                                                    
!         calculate the functions at x and                              
!         return this vector in fvec.                                   
!         ----------                                                    
!         return                                                        
!         end                                                           
!                                                                       
!         the value of iflag should not be changed by fcn unless        
!         the user wants to terminate execution of fdjac2.              
!         in this case set iflag to a negative integer.                 
!                                                                       
!       m is a positive integer input variable set to the number        
!         of functions.                                                 
!                                                                       
!       n is a positive integer input variable set to the number        
!         of variables. n must not exceed m.                            
!                                                                       
!       x is an input array of length n.                                
!                                                                       
!       fvec is an input array of length m which must contain the       
!         functions evaluated at x.                                     
!                                                                       
!       fjac is an output m by n array which contains the               
!         approximation to the jacobian matrix evaluated at x.          
!                                                                       
!       ldfjac is a positive integer input variable not less than m     
!         which specifies the leading dimension of the array fjac.      
!                                                                       
!       iflag is an integer variable which can be used to terminate     
!         the execution of fdjac2. see description of fcn.              
!                                                                       
!       epsfcn is an input variable used in determining a suitable      
!         step length for the forward-difference approximation. this    
!         approximation assumes that the relative errors in the         
!         functions are of the order of epsfcn. if epsfcn is less       
!         than the machine precision, it is assumed that the relative   
!         errors in the functions are of the order of the machine       
!         precision.                                                    
!                                                                       
!       wa is a work array of length m.                                 
!                                                                       
!     subprograms called                                                
!                                                                       
!       user-supplied ...... fcn                                        
!                                                                       
!       minpack-supplied ... dpmpar                                     
!                                                                       
!       fortran-supplied ... dabs,dmax1,dsqrt                           
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i,j 
      double precision eps,epsmch,h,temp,zero 
      double precision xdpmpar 
      data zero /0.0d0/ 
!                                                                       
!     epsmch is the machine precision.                                  
!                                                                       
      epsmch = dpmpar(1) 
!                                                                       
      eps = dsqrt(dmax1(epsfcn,epsmch)) 
      do 20 j = 1, n 
         temp = x(j) 
         h = eps*dabs(temp) 
         if (h .eq. zero) h = eps 
         x(j) = temp + h 
         call fcn(m,n,x,wa,iflag) 
         if (iflag .lt. 0) go to 30 
         x(j) = temp 
         do 10 i = 1, m 
            fjac(i,j) = (wa(i) - fvec(i))/h 
   10       continue 
   20    continue 
   30 continue 
      return 
!                                                                       
!     last card of subroutine fdjac2.                                   
!                                                                       
      END   

      recursive &                                   
      subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,     &
     &                 wa2)                                             
      integer n,ldr 
      integer ipvt(n) 
      double precision delta,par 
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),    &
     &                 wa2(n)                                           
!     **********                                                        
!                                                                       
!     subroutine lmpar                                                  
!                                                                       
!     given an m by n matrix a, an n by n nonsingular diagonal          
!     matrix d, an m-vector b, and a positive number delta,             
!     the problem is to determine a value for the parameter             
!     par such that if x solves the system                              
!                                                                       
!           a*x = b ,     sqrt(par)*d*x = 0 ,                           
!                                                                       
!     in the least squares sense, and dxnorm is the euclidean           
!     norm of d*x, then either par is zero and                          
!                                                                       
!           (dxnorm-delta) .le. 0.1*delta ,                             
!                                                                       
!     or par is positive and                                            
!                                                                       
!           abs(dxnorm-delta) .le. 0.1*delta .                          
!                                                                       
!     this subroutine completes the solution of the problem             
!     if it is provided with the necessary information from the         
!     qr factorization, with column pivoting, of a. that is, if         
!     a*p = q*r, where p is a permutation matrix, q has orthogonal      
!     columns, and r is an upper triangular matrix with diagonal        
!     elements of nonincreasing magnitude, then lmpar expects           
!     the full upper triangle of r, the permutation matrix p,           
!     and the first n components of (q transpose)*b. on output          
!     lmpar also provides an upper triangular matrix s such that        
!                                                                       
!            t   t                   t                                  
!           p *(a *a + par*d*d)*p = s *s .                              
!                                                                       
!     s is employed within lmpar and may be of separate interest.       
!                                                                       
!     only a few iterations are generally needed for convergence        
!     of the algorithm. if, however, the limit of 10 iterations         
!     is reached, then the output par will contain the best             
!     value obtained so far.                                            
!                                                                       
!     the subroutine statement is                                       
!                                                                       
!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,       
!                        wa1,wa2)                                       
!                                                                       
!     where                                                             
!                                                                       
!       n is a positive integer input variable set to the order of r.   
!                                                                       
!       r is an n by n array. on input the full upper triangle          
!         must contain the full upper triangle of the matrix r.         
!         on output the full upper triangle is unaltered, and the       
!         strict lower triangle contains the strict upper triangle      
!         (transposed) of the upper triangular matrix s.                
!                                                                       
!       ldr is a positive integer input variable not less than n        
!         which specifies the leading dimension of the array r.         
!                                                                       
!       ipvt is an integer input array of length n which defines the    
!         permutation matrix p such that a*p = q*r. column j of p       
!         is column ipvt(j) of the identity matrix.                     
!                                                                       
!       diag is an input array of length n which must contain the       
!         diagonal elements of the matrix d.                            
!                                                                       
!       qtb is an input array of length n which must contain the first  
!         n elements of the vector (q transpose)*b.                     
!                                                                       
!       delta is a positive input variable which specifies an upper     
!         bound on the euclidean norm of d*x.                           
!                                                                       
!       par is a nonnegative variable. on input par contains an         
!         initial estimate of the levenberg-marquardt parameter.        
!         on output par contains the final estimate.                    
!                                                                       
!       x is an output array of length n which contains the least       
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,    
!         for the output par.                                           
!                                                                       
!       sdiag is an output array of length n which contains the         
!         diagonal elements of the upper triangular matrix s.           
!                                                                       
!       wa1 and wa2 are work arrays of length n.                        
!                                                                       
!     subprograms called                                                
!                                                                       
!       minpack-supplied ... dpmpar,enorm,qrsolv                        
!                                                                       
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt                     
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i,iter,j,jm1,jp1,k,l,nsing 
      double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,    &
     &                 sum,temp,zero                                    
      double precision xdpmpar,xenorm 
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/ 
!                                                                       
!     dwarf is the smallest positive magnitude.                         
!                                                                       
      dwarf = dpmpar(2) 
!                                                                       
!     compute and store in x the gauss-newton direction. if the         
!     jacobian is rank-deficient, obtain a least squares solution.      
!                                                                       
      nsing = n 
      do 10 j = 1, n 
         wa1(j) = qtb(j) 
         if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1 
         if (nsing .lt. n) wa1(j) = zero 
   10    continue 
      if (nsing .lt. 1) go to 50 
      do 40 k = 1, nsing 
         j = nsing - k + 1 
         wa1(j) = wa1(j)/r(j,j) 
         temp = wa1(j) 
         jm1 = j - 1 
         if (jm1 .lt. 1) go to 30 
         do 20 i = 1, jm1 
            wa1(i) = wa1(i) - r(i,j)*temp 
   20       continue 
   30    continue 
   40    continue 
   50 continue 
      do 60 j = 1, n 
         l = ipvt(j) 
         x(l) = wa1(j) 
   60    continue 
!                                                                       
!     initialize the iteration counter.                                 
!     evaluate the function at the origin, and test                     
!     for acceptance of the gauss-newton direction.                     
!                                                                       
      iter = 0 
      do 70 j = 1, n 
         wa2(j) = diag(j)*x(j) 
   70    continue 
      dxnorm = enorm(n,wa2) 
      fp = dxnorm - delta 
      if (fp .le. p1*delta) go to 220 
!                                                                       
!     if the jacobian is not rank deficient, the newton                 
!     step provides a lower bound, parl, for the zero of                
!     the function. otherwise set this bound to zero.                   
!                                                                       
      parl = zero 
      if (nsing .lt. n) go to 120 
      do 80 j = 1, n 
         l = ipvt(j) 
         wa1(j) = diag(l)*(wa2(l)/dxnorm) 
   80    continue 
      do 110 j = 1, n 
         sum = zero 
         jm1 = j - 1 
         if (jm1 .lt. 1) go to 100 
         do 90 i = 1, jm1 
            sum = sum + r(i,j)*wa1(i) 
   90       continue 
  100    continue 
         wa1(j) = (wa1(j) - sum)/r(j,j) 
  110    continue 
      temp = enorm(n,wa1) 
      parl = ((fp/delta)/temp)/temp 
  120 continue 
!                                                                       
!     calculate an upper bound, paru, for the zero of the function.     
!                                                                       
      do 140 j = 1, n 
         sum = zero 
         do 130 i = 1, j 
            sum = sum + r(i,j)*qtb(i) 
  130       continue 
         l = ipvt(j) 
         wa1(j) = sum/diag(l) 
  140    continue 
      gnorm = enorm(n,wa1) 
      paru = gnorm/delta 
      if (paru .eq. zero) paru = dwarf/dmin1(delta,p1) 
!                                                                       
!     if the input par lies outside of the interval (parl,paru),        
!     set par to the closer endpoint.                                   
!                                                                       
      par = dmax1(par,parl) 
      par = dmin1(par,paru) 
      if (par .eq. zero) par = gnorm/dxnorm 
!                                                                       
!     beginning of an iteration.                                        
!                                                                       
  150 continue 
         iter = iter + 1 
!                                                                       
!        evaluate the function at the current value of par.             
!                                                                       
         if (par .eq. zero) par = dmax1(dwarf,p001*paru) 
         temp = dsqrt(par) 
         do 160 j = 1, n 
            wa1(j) = temp*diag(j) 
  160       continue 
         call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2) 
         do 170 j = 1, n 
            wa2(j) = diag(j)*x(j) 
  170       continue 
         dxnorm = enorm(n,wa2) 
         temp = fp 
         fp = dxnorm - delta 
!                                                                       
!        if the function is small enough, accept the current value      
!        of par. also test for the exceptional cases where parl         
!        is zero or the number of iterations has reached 10.            
!                                                                       
         if (dabs(fp) .le. p1*delta                                     &
     &       .or. parl .eq. zero .and. fp .le. temp                     &
     &            .and. temp .lt. zero .or. iter .eq. 10) go to 220     
!                                                                       
!        compute the newton correction.                                 
!                                                                       
         do 180 j = 1, n 
            l = ipvt(j) 
            wa1(j) = diag(l)*(wa2(l)/dxnorm) 
  180       continue 
         do 210 j = 1, n 
            wa1(j) = wa1(j)/sdiag(j) 
            temp = wa1(j) 
            jp1 = j + 1 
            if (n .lt. jp1) go to 200 
            do 190 i = jp1, n 
               wa1(i) = wa1(i) - r(i,j)*temp 
  190          continue 
  200       continue 
  210       continue 
         temp = enorm(n,wa1) 
         parc = ((fp/delta)/temp)/temp 
!                                                                       
!        depending on the sign of the function, update parl or paru.    
!                                                                       
         if (fp .gt. zero) parl = dmax1(parl,par) 
         if (fp .lt. zero) paru = dmin1(paru,par) 
!                                                                       
!        compute an improved estimate for par.                          
!                                                                       
         par = dmax1(parl,par+parc) 
!                                                                       
!        end of an iteration.                                           
!                                                                       
         go to 150 
  220 continue 
!                                                                       
!     termination.                                                      
!                                                                       
      if (iter .eq. 0) par = zero 
      return 
!                                                                       
!     last card of subroutine lmpar.                                    
!                                                                       
      END 

      recursive &                                          
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) 
      integer m,n,lda,lipvt 
      integer ipvt(lipvt) 
      logical pivot 
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n) 
!     **********                                                        
!                                                                       
!     subroutine qrfac                                                  
!                                                                       
!     this subroutine uses householder transformations with column      
!     pivoting (optional) to compute a qr factorization of the          
!     m by n matrix a. that is, qrfac determines an orthogonal          
!     matrix q, a permutation matrix p, and an upper trapezoidal        
!     matrix r with diagonal elements of nonincreasing magnitude,       
!     such that a*p = q*r. the householder transformation for           
!     column k, k = 1,2,...,min(m,n), is of the form                    
!                                                                       
!                           t                                           
!           i - (1/u(k))*u*u                                            
!                                                                       
!     where u has zeros in the first k-1 positions. the form of         
!     this transformation and the method of pivoting first              
!     appeared in the corresponding linpack subroutine.                 
!                                                                       
!     the subroutine statement is                                       
!                                                                       
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)    
!                                                                       
!     where                                                             
!                                                                       
!       m is a positive integer input variable set to the number        
!         of rows of a.                                                 
!                                                                       
!       n is a positive integer input variable set to the number        
!         of columns of a.                                              
!                                                                       
!       a is an m by n array. on input a contains the matrix for        
!         which the qr factorization is to be computed. on output       
!         the strict upper trapezoidal part of a contains the strict    
!         upper trapezoidal part of r, and the lower trapezoidal        
!         part of a contains a factored form of q (the non-trivial      
!         elements of the u vectors described above).                   
!                                                                       
!       lda is a positive integer input variable not less than m        
!         which specifies the leading dimension of the array a.         
!                                                                       
!       pivot is a logical input variable. if pivot is set true,        
!         then column pivoting is enforced. if pivot is set false,      
!         then no column pivoting is done.                              
!                                                                       
!       ipvt is an integer output array of length lipvt. ipvt           
!         defines the permutation matrix p such that a*p = q*r.         
!         column j of p is column ipvt(j) of the identity matrix.       
!         if pivot is false, ipvt is not referenced.                    
!                                                                       
!       lipvt is a positive integer input variable. if pivot is false,  
!         then lipvt may be as small as 1. if pivot is true, then       
!         lipvt must be at least n.                                     
!                                                                       
!       rdiag is an output array of length n which contains the         
!         diagonal elements of r.                                       
!                                                                       
!       acnorm is an output array of length n which contains the        
!         norms of the corresponding columns of the input matrix a.     
!         if this information is not needed, then acnorm can coincide   
!         with rdiag.                                                   
!                                                                       
!       wa is a work array of length n. if pivot is false, then wa      
!         can coincide with rdiag.                                      
!                                                                       
!     subprograms called                                                
!                                                                       
!       minpack-supplied ... dpmpar,enorm                               
!                                                                       
!       fortran-supplied ... dmax1,dsqrt,min0                           
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i,j,jp1,k,kmax,minmn 
      double precision ajnorm,epsmch,one,p05,sum,temp,zero 
      double precision xdpmpar,xenorm 
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/ 
!                                                                       
!     epsmch is the machine precision.                                  
!                                                                       
      epsmch = dpmpar(1) 
!                                                                       
!     compute the initial column norms and initialize several arrays.   
!                                                                       
      do 10 j = 1, n 
         acnorm(j) = enorm(m,a(1,j)) 
         rdiag(j) = acnorm(j) 
         wa(j) = rdiag(j) 
         if (pivot) ipvt(j) = j 
   10    continue 
!                                                                       
!     reduce a to r with householder transformations.                   
!                                                                       
      minmn = min0(m,n) 
      do 110 j = 1, minmn 
         if (.not.pivot) go to 40 
!                                                                       
!        bring the column of largest norm into the pivot position.      
!                                                                       
         kmax = j 
         do 20 k = j, n 
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k 
   20       continue 
         if (kmax .eq. j) go to 40 
         do 30 i = 1, m 
            temp = a(i,j) 
            a(i,j) = a(i,kmax) 
            a(i,kmax) = temp 
   30       continue 
         rdiag(kmax) = rdiag(j) 
         wa(kmax) = wa(j) 
         k = ipvt(j) 
         ipvt(j) = ipvt(kmax) 
         ipvt(kmax) = k 
   40    continue 
!                                                                       
!        compute the householder transformation to reduce the           
!        j-th column of a to a multiple of the j-th unit vector.        
!                                                                       
         ajnorm = enorm(m-j+1,a(j,j)) 
         if (ajnorm .eq. zero) go to 100 
         if (a(j,j) .lt. zero) ajnorm = -ajnorm 
         do 50 i = j, m 
            a(i,j) = a(i,j)/ajnorm 
   50       continue 
         a(j,j) = a(j,j) + one 
!                                                                       
!        apply the transformation to the remaining columns              
!        and update the norms.                                          
!                                                                       
         jp1 = j + 1 
         if (n .lt. jp1) go to 100 
         do 90 k = jp1, n 
            sum = zero 
            do 60 i = j, m 
               sum = sum + a(i,j)*a(i,k) 
   60          continue 
            temp = sum/a(j,j) 
            do 70 i = j, m 
               a(i,k) = a(i,k) - temp*a(i,j) 
   70          continue 
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80 
            temp = a(j,k)/rdiag(k) 
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2)) 
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80 
            rdiag(k) = enorm(m-j,a(jp1,k)) 
            wa(k) = rdiag(k) 
   80       continue 
   90       continue 
  100    continue 
         rdiag(j) = -ajnorm 
  110    continue 
      return 
!                                                                       
!     last card of subroutine qrfac.                                    
!                                                                       
      END  

      recursive &                                         
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa) 
      integer n,ldr 
      integer ipvt(n) 
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n) 
!     **********                                                        
!                                                                       
!     subroutine qrsolv                                                 
!                                                                       
!     given an m by n matrix a, an n by n diagonal matrix d,            
!     and an m-vector b, the problem is to determine an x which         
!     solves the system                                                 
!                                                                       
!           a*x = b ,     d*x = 0 ,                                     
!                                                                       
!     in the least squares sense.                                       
!                                                                       
!     this subroutine completes the solution of the problem             
!     if it is provided with the necessary information from the         
!     qr factorization, with column pivoting, of a. that is, if         
!     a*p = q*r, where p is a permutation matrix, q has orthogonal      
!     columns, and r is an upper triangular matrix with diagonal        
!     elements of nonincreasing magnitude, then qrsolv expects          
!     the full upper triangle of r, the permutation matrix p,           
!     and the first n components of (q transpose)*b. the system         
!     a*x = b, d*x = 0, is then equivalent to                           
!                                                                       
!                  t       t                                            
!           r*z = q *b ,  p *d*p*z = 0 ,                                
!                                                                       
!     where x = p*z. if this system does not have full rank,            
!     then a least squares solution is obtained. on output qrsolv       
!     also provides an upper triangular matrix s such that              
!                                                                       
!            t   t               t                                      
!           p *(a *a + d*d)*p = s *s .                                  
!                                                                       
!     s is computed within qrsolv and may be of separate interest.      
!                                                                       
!     the subroutine statement is                                       
!                                                                       
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)             
!                                                                       
!     where                                                             
!                                                                       
!       n is a positive integer input variable set to the order of r.   
!                                                                       
!       r is an n by n array. on input the full upper triangle          
!         must contain the full upper triangle of the matrix r.         
!         on output the full upper triangle is unaltered, and the       
!         strict lower triangle contains the strict upper triangle      
!         (transposed) of the upper triangular matrix s.                
!                                                                       
!       ldr is a positive integer input variable not less than n        
!         which specifies the leading dimension of the array r.         
!                                                                       
!       ipvt is an integer input array of length n which defines the    
!         permutation matrix p such that a*p = q*r. column j of p       
!         is column ipvt(j) of the identity matrix.                     
!                                                                       
!       diag is an input array of length n which must contain the       
!         diagonal elements of the matrix d.                            
!                                                                       
!       qtb is an input array of length n which must contain the first  
!         n elements of the vector (q transpose)*b.                     
!                                                                       
!       x is an output array of length n which contains the least       
!         squares solution of the system a*x = b, d*x = 0.              
!                                                                       
!       sdiag is an output array of length n which contains the         
!         diagonal elements of the upper triangular matrix s.           
!                                                                       
!       wa is a work array of length n.                                 
!                                                                       
!     subprograms called                                                
!                                                                       
!       fortran-supplied ... dabs,dsqrt                                 
!                                                                       
!     argonne national laboratory. minpack project. march 1980.         
!     burton s. garbow, kenneth e. hillstrom, jorge j. more             
!                                                                       
!     **********                                                        
      integer i,j,jp1,k,kp1,l,nsing 
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero 
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/ 
!                                                                       
!     copy r and (q transpose)*b to preserve input and initialize s.    
!     in particular, save the diagonal elements of r in x.              
!                                                                       
      do 20 j = 1, n 
         do 10 i = j, n 
            r(i,j) = r(j,i) 
   10       continue 
         x(j) = r(j,j) 
         wa(j) = qtb(j) 
   20    continue 
!                                                                       
!     eliminate the diagonal matrix d using a givens rotation.          
!                                                                       
      do 100 j = 1, n 
!                                                                       
!        prepare the row of d to be eliminated, locating the            
!        diagonal element using p from the qr factorization.            
!                                                                       
         l = ipvt(j) 
         if (diag(l) .eq. zero) go to 90 
         do 30 k = j, n 
            sdiag(k) = zero 
   30       continue 
         sdiag(j) = diag(l) 
!                                                                       
!        the transformations to eliminate the row of d                  
!        modify only a single element of (q transpose)*b                
!        beyond the first n, which is initially zero.                   
!                                                                       
         qtbpj = zero 
         do 80 k = j, n 
!                                                                       
!           determine a givens rotation which eliminates the            
!           appropriate element in the current row of d.                
!                                                                       
            if (sdiag(k) .eq. zero) go to 70 
            if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40 
               cotan = r(k,k)/sdiag(k) 
               sin = p5/dsqrt(p25+p25*cotan**2) 
               cos = sin*cotan 
               go to 50 
   40       continue 
               tan = sdiag(k)/r(k,k) 
               cos = p5/dsqrt(p25+p25*tan**2) 
               sin = cos*tan 
   50       continue 
!                                                                       
!           compute the modified diagonal element of r and              
!           the modified element of ((q transpose)*b,0).                
!                                                                       
            r(k,k) = cos*r(k,k) + sin*sdiag(k) 
            temp = cos*wa(k) + sin*qtbpj 
            qtbpj = -sin*wa(k) + cos*qtbpj 
            wa(k) = temp 
!                                                                       
!           accumulate the tranformation in the row of s.               
!                                                                       
            kp1 = k + 1 
            if (n .lt. kp1) go to 70 
            do 60 i = kp1, n 
               temp = cos*r(i,k) + sin*sdiag(i) 
               sdiag(i) = -sin*r(i,k) + cos*sdiag(i) 
               r(i,k) = temp 
   60          continue 
   70       continue 
   80       continue 
   90    continue 
!                                                                       
!        store the diagonal element of s and restore                    
!        the corresponding diagonal element of r.                       
!                                                                       
         sdiag(j) = r(j,j) 
         r(j,j) = x(j) 
  100    continue 
!                                                                       
!     solve the triangular system for z. if the system is               
!     singular, then obtain a least squares solution.                   
!                                                                       
      nsing = n 
      do 110 j = 1, n 
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1 
         if (nsing .lt. n) wa(j) = zero 
  110    continue 
      if (nsing .lt. 1) go to 150 
      do 140 k = 1, nsing 
         j = nsing - k + 1 
         sum = zero 
         jp1 = j + 1 
         if (nsing .lt. jp1) go to 130 
         do 120 i = jp1, nsing 
            sum = sum + r(i,j)*wa(i) 
  120       continue 
  130    continue 
         wa(j) = (wa(j) - sum)/sdiag(j) 
  140    continue 
  150 continue 
!                                                                       
!     permute the components of z back to components of x.              
!                                                                       
      do 160 j = 1, n 
         l = ipvt(j) 
         x(l) = wa(j) 
  160    continue 
      return 
!                                                                       
!     last card of subroutine qrsolv.                                   
!                                                                       
      END                                           


     recursive &
     subroutine matrix_inverse( A , Ainverse, n)
       implicit none
       integer, intent(in)          :: n
       double precision, intent(in) :: A(n,n)
       double precision,intent(out) :: Ainverse(n,n)
 
       double precision :: Aw(n,n), wk(n), rcond
       integer          :: ipvt(n)
       double precision :: B
       integer :: j

       if(n==1) then
         Ainverse(1,1) = 1d0 / A(1,1)
         return
       endif

       Aw = A
       CALL DGECO(Aw,N,N,IPVT,RCOND,wk)
       B = 0
       do j=1,n
         Ainverse(j,j) = 1d0
       enddo
       do j=1,n
          CALL DGESL(Aw,N,N,IPVT,Ainverse(1,j),0)
       enddo

     end subroutine matrix_inverse


!!!! the generic n-exp function

       function fm(t, px, npar)
         implicit none
         double precision, intent(in)              :: t
         double precision, intent(in)              :: px(npar)
         integer         , intent(in)              :: npar

         double precision              :: fm 
         integer                       :: i


         fm = 0d0
     
         do i=1,npar/2
            fm = fm + px(2*i-1) * exp( -t * px(2*i) )       
         enddo
         

       end function fm







     
!DECK DASUM                                                             
      DOUBLE PRECISION FUNCTION DASUM (N, DX, INCX) 
!***BEGIN PROLOGUE  DASUM                                               
!***PURPOSE  Compute the sum of the magnitudes of the elements of a     
!            vector.                                                    
!***LIBRARY   SLATEC (BLAS)                                             
!***CATEGORY  D1A3A                                                     
!***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)             
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR       
!***AUTHOR  Lawson, C. L., (JPL)                                        
!           Hanson, R. J., (SNLA)                                       
!           Kincaid, D. R., (U. of Texas)                               
!           Krogh, F. T., (JPL)                                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  Subprogram                                    
!    Description of Parameters                                          
!                                                                       
!     --Input--                                                         
!        N  number of elements in input vector(s)                       
!       DX  double precision vector with N elements                     
!     INCX  storage spacing between elements of DX                      
!                                                                       
!     --Output--                                                        
!    DASUM  double precision result (zero if N .LE. 0)                  
!                                                                       
!     Returns sum of magnitudes of double precision DX.                 
!     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)),                  
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.              
!                                                                       
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.     
!                 Krogh, Basic linear algebra subprograms for Fortran   
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.          
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791001  DATE WRITTEN                                                
!   890531  Changed all specific intrinsics to generic.  (WRB)          
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900821  Modified to correct problem with a negative increment.      
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DASUM                                                 
      DOUBLE PRECISION DX(*) 
      INTEGER I, INCX, IX, M, MP1, N 
!***FIRST EXECUTABLE STATEMENT  DASUM                                   
      DASUM = 0.0D0 
      IF (N .LE. 0) RETURN 
!                                                                       
      IF (INCX .EQ. 1) GOTO 20 
!                                                                       
!     Code for increment not equal to 1.                                
!                                                                       
      IX = 1 
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1 
      DO 10 I = 1,N 
        DASUM = DASUM + ABS(DX(IX)) 
        IX = IX + INCX 
   10 END DO 
      RETURN 
!                                                                       
!     Code for increment equal to 1.                                    
!                                                                       
!     Clean-up loop so remaining vector length is a multiple of 6.      
!                                                                       
   20 M = MOD(N,6) 
      IF (M .EQ. 0) GOTO 40 
      DO 30 I = 1,M 
        DASUM = DASUM + ABS(DX(I)) 
   30 END DO 
      IF (N .LT. 6) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,6 
        DASUM = DASUM + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2)) +      &
     &          ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))              
   50 END DO 
      RETURN 
      END                                           
!DECK DAXPY                                                             
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY) 
      IMPLICIT NONE ! (PAZ)
!***BEGIN PROLOGUE  DAXPY                                               
!***PURPOSE  Compute a constant times a vector plus a vector.           
!***LIBRARY   SLATEC (BLAS)                                             
!***CATEGORY  D1A7                                                      
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)              
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR                       
!***AUTHOR  Lawson, C. L., (JPL)                                        
!           Hanson, R. J., (SNLA)                                       
!           Kincaid, D. R., (U. of Texas)                               
!           Krogh, F. T., (JPL)                                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  Subprogram                                    
!    Description of Parameters                                          
!                                                                       
!     --Input--                                                         
!        N  number of elements in input vector(s)                       
!       DA  double precision scalar multiplier                          
!       DX  double precision vector with N elements                     
!     INCX  storage spacing between elements of DX                      
!       DY  double precision vector with N elements                     
!     INCY  storage spacing between elements of DY                      
!                                                                       
!     --Output--                                                        
!       DY  double precision result (unchanged if N .LE. 0)             
!                                                                       
!     Overwrite double precision DY with double precision DA*DX + DY.   
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +  
!       DY(LY+I*INCY),                                                  
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is    
!     defined in a similar way using INCY.                              
!                                                                       
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.     
!                 Krogh, Basic linear algebra subprograms for Fortran   
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.          
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791001  DATE WRITTEN                                                
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)           
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DAXPY                                                 
      DOUBLE PRECISION DX(*), DY(*), DA 
      INTEGER N
      INTEGER INCX,INCY
!
      INTEGER IX, IY
      INTEGER I, M, MP1, NS
!***FIRST EXECUTABLE STATEMENT  DAXPY                                   
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN 
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60 
!                                                                       
!     Code for unequal or nonpositive increments.                       
!                                                                       
    5 IX = 1 
      IY = 1 
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1 
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        DY(IY) = DY(IY) + DA*DX(IX) 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!     Code for both increments equal to 1.                              
!                                                                       
!     Clean-up loop so remaining vector length is a multiple of 4.      
!                                                                       
   20 M = MOD(N,4) 
      IF (M .EQ. 0) GO TO 40 
      DO 30 I = 1,M 
        DY(I) = DY(I) + DA*DX(I) 
   30 END DO 
      IF (N .LT. 4) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,4 
        DY(I) = DY(I) + DA*DX(I) 
        DY(I+1) = DY(I+1) + DA*DX(I+1) 
        DY(I+2) = DY(I+2) + DA*DX(I+2) 
        DY(I+3) = DY(I+3) + DA*DX(I+3) 
   50 END DO 
      RETURN 
!                                                                       
!     Code for equal, positive, non-unit increments.                    
!                                                                       
   60 NS = N*INCX 
      DO 70 I = 1,NS,INCX 
        DY(I) = DA*DX(I) + DY(I) 
   70 END DO 
      RETURN 
      END                                           
!DECK DDOT                                                              
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY) 
      IMPLICIT NONE ! (PAZ)
!***BEGIN PROLOGUE  DDOT                                                
!***PURPOSE  Compute the inner product of two vectors.                  
!***LIBRARY   SLATEC (BLAS)                                             
!***CATEGORY  D1A4                                                      
!***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)                
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR               
!***AUTHOR  Lawson, C. L., (JPL)                                        
!           Hanson, R. J., (SNLA)                                       
!           Kincaid, D. R., (U. of Texas)                               
!           Krogh, F. T., (JPL)                                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  Subprogram                                    
!    Description of Parameters                                          
!                                                                       
!     --Input--                                                         
!        N  number of elements in input vector(s)                       
!       DX  double precision vector with N elements                     
!     INCX  storage spacing between elements of DX                      
!       DY  double precision vector with N elements                     
!     INCY  storage spacing between elements of DY                      
!                                                                       
!     --Output--                                                        
!     DDOT  double precision dot product (zero if N .LE. 0)             
!                                                                       
!     Returns the dot product of double precision DX and DY.            
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),    
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is    
!     defined in a similar way using INCY.                              
!                                                                       
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.     
!                 Krogh, Basic linear algebra subprograms for Fortran   
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.          
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791001  DATE WRITTEN                                                
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)           
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DDOT                                                  
      DOUBLE PRECISION DX(*), DY(*) 
      INTEGER N
      INTEGER INCX,INCY
!***
      INTEGER IX, IY
      INTEGER I, M, MP1, NS
!***FIRST EXECUTABLE STATEMENT  DDOT                                    
      DDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60 
!                                                                       
!     Code for unequal or nonpositive increments.                       
!                                                                       
    5 IX = 1 
      IY = 1 
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1 
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        DDOT = DDOT + DX(IX)*DY(IY) 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!     Code for both increments equal to 1.                              
!                                                                       
!     Clean-up loop so remaining vector length is a multiple of 5.      
!                                                                       
   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO 30 I = 1,M 
         DDOT = DDOT + DX(I)*DY(I) 
   30 END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,5 
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +   &
     &              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)                   
   50 END DO 
      RETURN 
!                                                                       
!     Code for equal, positive, non-unit increments.                    
!                                                                       
   60 NS = N*INCX 
      DO 70 I = 1,NS,INCX 
        DDOT = DDOT + DX(I)*DY(I) 
   70 END DO 
      RETURN 
      END                                           
!DECK DGECO 
      recursive &                                                            
      SUBROUTINE DGECO (A, LDA, N, IPVT, RCOND, Z) 
!***BEGIN PROLOGUE  DGECO                                               
!***PURPOSE  Factor a matrix using Gaussian elimination and estimate    
!            the condition number of the matrix.                        
!***LIBRARY   SLATEC (LINPACK)                                          
!***CATEGORY  D2A1                                                      
!***TYPE      DOUBLE PRECISION (SGECO-S, DGECO-D, CGECO-C)              
!***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION                                      
!***AUTHOR  Moler, C. B., (U. of New Mexico)                            
!***DESCRIPTION                                                         
!                                                                       
!     DGECO factors a double precision matrix by Gaussian elimination   
!     and estimates the condition of the matrix.                        
!                                                                       
!     If  RCOND  is not needed, DGEFA is slightly faster.               
!     To solve  A*X = B , follow DGECO by DGESL.                        
!     To compute  INVERSE(A)*C , follow DGECO by DGESL.                 
!     To compute  DETERMINANT(A) , follow DGECO by DGEDI.               
!     To compute  INVERSE(A) , follow DGECO by DGEDI.                   
!                                                                       
!     On Entry                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                the matrix to be factored.                             
!                                                                       
!        LDA     INTEGER                                                
!                the leading dimension of the array  A .                
!                                                                       
!        N       INTEGER                                                
!                the order of the matrix  A .                           
!                                                                       
!     On Return                                                         
!                                                                       
!        A       an upper triangular matrix and the multipliers         
!                which were used to obtain it.                          
!                The factorization can be written  A = L*U  where       
!                L  is a product of permutation and unit lower          
!                triangular matrices and  U  is upper triangular.       
!                                                                       
!        IPVT    INTEGER(N)                                             
!                an INTEGER vector of pivot indices.                    
!                                                                       
!        RCOND   DOUBLE PRECISION                                       
!                an estimate of the reciprocal condition of  A .        
!                For the system  A*X = B , relative perturbations       
!                in  A  and  B  of size  EPSILON  may cause             
!                relative perturbations in  X  of size  EPSILON/RCOND . 
!                If  RCOND  is so small that the logical expression     
!                           1.0 + RCOND .EQ. 1.0                        
!                is true, then  A  may be singular to working           
!                precision.  In particular,  RCOND  is zero  if         
!                exact singularity is detected or the estimate          
!                underflows.                                            
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                a work vector whose contents are usually unimportant.  
!                If  A  is close to a singular matrix, then  Z  is      
!                an approximate null vector in the sense that           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                                                                       
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.     
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.            
!***ROUTINES CALLED  DASUM, DAXPY, DDOT, DGEFA, DSCAL                   
!***REVISION HISTORY  (YYMMDD)                                          
!   780814  DATE WRITTEN                                                
!   890531  Changed all specific intrinsics to generic.  (WRB)          
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900326  Removed duplicate information from DESCRIPTION section.     
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DGECO                                                 
      INTEGER LDA,N,IPVT(*) 
      DOUBLE PRECISION A(LDA,*),Z(*) 
      DOUBLE PRECISION RCOND 
!                                                                       
      DOUBLE PRECISION xDDOT,EK,T,WK,WKM 
      DOUBLE PRECISION ANORM,S,xDASUM,SM,YNORM 
      INTEGER INFO,J,K,KB,KP1,L 
!                                                                       
!     COMPUTE 1-NORM OF A                                               
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DGECO                                   
      ANORM = 0.0D0 
      DO 10 J = 1, N 
         ANORM = MAX(ANORM,DASUM(N,A(1,J),1)) 
   10 END DO 
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DGEFA(A,LDA,N,IPVT,INFO) 
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    
!     OVERFLOW.                                                         
!                                                                       
!     SOLVE TRANS(U)*W = E                                              
!                                                                       
      EK = 1.0D0 
      DO 20 J = 1, N 
         Z(J) = 0.0D0 
   20 END DO 
      DO 100 K = 1, N 
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,-Z(K)) 
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30 
            S = ABS(A(K,K))/ABS(EK-Z(K)) 
            CALL DSCAL(N,S,Z,1) 
            EK = S*EK 
   30    CONTINUE 
         WK = EK - Z(K) 
         WKM = -EK - Z(K) 
         S = ABS(WK) 
         SM = ABS(WKM) 
         IF (A(K,K) .EQ. 0.0D0) GO TO 40 
            WK = WK/A(K,K) 
            WKM = WKM/A(K,K) 
         GO TO 50 
   40    CONTINUE 
            WK = 1.0D0 
            WKM = 1.0D0 
   50    CONTINUE 
         KP1 = K + 1 
         IF (KP1 .GT. N) GO TO 90 
            DO 60 J = KP1, N 
               SM = SM + ABS(Z(J)+WKM*A(K,J)) 
               Z(J) = Z(J) + WK*A(K,J) 
               S = S + ABS(Z(J)) 
   60       CONTINUE 
            IF (S .GE. SM) GO TO 80 
               T = WKM - WK 
               WK = WKM 
               DO 70 J = KP1, N 
                  Z(J) = Z(J) + T*A(K,J) 
   70          CONTINUE 
   80       CONTINUE 
   90    CONTINUE 
         Z(K) = WK 
  100 END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL(N,S,Z,1) 
!                                                                       
!     SOLVE TRANS(L)*Y = W                                              
!                                                                       
      DO 120 KB = 1, N 
         K = N + 1 - KB 
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1) 
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 110 
            S = 1.0D0/ABS(Z(K)) 
            CALL DSCAL(N,S,Z,1) 
  110    CONTINUE 
         L = IPVT(K) 
         T = Z(L) 
         Z(L) = Z(K) 
         Z(K) = T 
  120 END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL(N,S,Z,1) 
!                                                                       
      YNORM = 1.0D0 
!                                                                       
!     SOLVE L*V = Y                                                     
!                                                                       
      DO 140 K = 1, N 
         L = IPVT(K) 
         T = Z(L) 
         Z(L) = Z(K) 
         Z(K) = T 
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1) 
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 130 
            S = 1.0D0/ABS(Z(K)) 
            CALL DSCAL(N,S,Z,1) 
            YNORM = S*YNORM 
  130    CONTINUE 
  140 END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL(N,S,Z,1) 
      YNORM = S*YNORM 
!                                                                       
!     SOLVE  U*Z = V                                                    
!                                                                       
      DO 160 KB = 1, N 
         K = N + 1 - KB 
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150 
            S = ABS(A(K,K))/ABS(Z(K)) 
            CALL DSCAL(N,S,Z,1) 
            YNORM = S*YNORM 
  150    CONTINUE 
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K) 
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0 
         T = -Z(K) 
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1) 
  160 END DO 
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL(N,S,Z,1) 
      YNORM = S*YNORM 
!                                                                       
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM 
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0 
      RETURN 
      END                                           
!DECK DGEFA   
      recursive &                                                          
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO) 
!***BEGIN PROLOGUE  DGEFA                                               
!***PURPOSE  Factor a matrix using Gaussian elimination.                
!***LIBRARY   SLATEC (LINPACK)                                          
!***CATEGORY  D2A1                                                      
!***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)              
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,                  
!             MATRIX FACTORIZATION                                      
!***AUTHOR  Moler, C. B., (U. of New Mexico)                            
!***DESCRIPTION                                                         
!                                                                       
!     DGEFA factors a double precision matrix by Gaussian elimination.  
!                                                                       
!     DGEFA is usually called by DGECO, but it can be called            
!     directly with a saving in time if  RCOND  is not needed.          
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .                   
!                                                                       
!     On Entry                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                the matrix to be factored.                             
!                                                                       
!        LDA     INTEGER                                                
!                the leading dimension of the array  A .                
!                                                                       
!        N       INTEGER                                                
!                the order of the matrix  A .                           
!                                                                       
!     On Return                                                         
!                                                                       
!        A       an upper triangular matrix and the multipliers         
!                which were used to obtain it.                          
!                The factorization can be written  A = L*U  where       
!                L  is a product of permutation and unit lower          
!                triangular matrices and  U  is upper triangular.       
!                                                                       
!        IPVT    INTEGER(N)                                             
!                an integer vector of pivot indices.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  normal value.                                     
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error       
!                     condition for this subroutine, but it does        
!                     indicate that DGESL or DGEDI will divide by zero  
!                     if called.  Use  RCOND  in DGECO for a reliable   
!                     indication of singularity.                        
!                                                                       
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.     
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.            
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX                               
!***REVISION HISTORY  (YYMMDD)                                          
!   780814  DATE WRITTEN                                                
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900326  Removed duplicate information from DESCRIPTION section.     
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DGEFA                                                 
      INTEGER LDA,N,IPVT(*),INFO 
      DOUBLE PRECISION A(LDA,*) 
!                                                                       
      DOUBLE PRECISION T 
      INTEGER xIDAMAX,J,K,KP1,L,NM1 
!                                                                       
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
!                                                                       
!***FIRST EXECUTABLE STATEMENT  DGEFA                                   
      INFO = 0 
      NM1 = N - 1 
      IF (NM1 .LT. 1) GO TO 70 
      DO 60 K = 1, NM1 
         KP1 = K + 1 
!                                                                       
!        FIND L = PIVOT INDEX                                           
!                                                                       
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1 
         IPVT(K) = L 
!                                                                       
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
!                                                                       
         IF (A(L,K) .EQ. 0.0D0) GO TO 40 
!                                                                       
!           INTERCHANGE IF NECESSARY                                    
!                                                                       
            IF (L .EQ. K) GO TO 10 
               T = A(L,K) 
               A(L,K) = A(K,K) 
               A(K,K) = T 
   10       CONTINUE 
!                                                                       
!           COMPUTE MULTIPLIERS                                         
!                                                                       
            T = -1.0D0/A(K,K) 
            CALL DSCAL(N-K,T,A(K+1,K),1) 
!                                                                       
!           ROW ELIMINATION WITH COLUMN INDEXING                        
!                                                                       
            DO 30 J = KP1, N 
               T = A(L,J) 
               IF (L .EQ. K) GO TO 20 
                  A(L,J) = A(K,J) 
                  A(K,J) = T 
   20          CONTINUE 
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1) 
   30       CONTINUE 
         GO TO 50 
   40    CONTINUE 
            INFO = K 
   50    CONTINUE 
   60 END DO 
   70 CONTINUE 
      IPVT(N) = N 
      IF (A(N,N) .EQ. 0.0D0) INFO = N 
      RETURN 
      END                                           
!DECK DGESL     
      recursive &                                                        
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB) 
!***BEGIN PROLOGUE  DGESL                                               
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the      
!            factors computed by DGECO or DGEFA.                        
!***LIBRARY   SLATEC (LINPACK)                                          
!***CATEGORY  D2A1                                                      
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)              
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE                    
!***AUTHOR  Moler, C. B., (U. of New Mexico)                            
!***DESCRIPTION                                                         
!                                                                       
!     DGESL solves the double precision system                          
!     A * X = B  or  TRANS(A) * X = B                                   
!     using the factors computed by DGECO or DGEFA.                     
!                                                                       
!     On Entry                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                the output from DGECO or DGEFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                the leading dimension of the array  A .                
!                                                                       
!        N       INTEGER                                                
!                the order of the matrix  A .                           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                the pivot vector from DGECO or DGEFA.                  
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                the right hand side vector.                            
!                                                                       
!        JOB     INTEGER                                                
!                = 0         to solve  A*X = B ,                        
!                = nonzero   to solve  TRANS(A)*X = B  where            
!                            TRANS(A)  is the transpose.                
!                                                                       
!     On Return                                                         
!                                                                       
!        B       the solution vector  X .                               
!                                                                       
!     Error Condition                                                   
!                                                                       
!        A division by zero will occur if the input factor contains a   
!        zero on the diagonal.  Technically this indicates singularity  
!        but it is often caused by improper arguments or improper       
!        setting of LDA .  It will not occur if the subroutines are     
!        called correctly and if DGECO has set RCOND .GT. 0.0           
!        or DGEFA has set INFO .EQ. 0 .                                 
!                                                                       
!     To compute  INVERSE(A) * C  where  C  is a matrix                 
!     with  P  columns                                                  
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            
!           IF (RCOND is too small) GO TO ...                           
!           DO 10 J = 1, P                                              
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        
!        10 CONTINUE                                                    
!                                                                       
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.     
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.            
!***ROUTINES CALLED  DAXPY, DDOT                                        
!***REVISION HISTORY  (YYMMDD)                                          
!   780814  DATE WRITTEN                                                
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900326  Removed duplicate information from DESCRIPTION section.     
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DGESL                                                 
      INTEGER LDA,N,IPVT(*),JOB 
      DOUBLE PRECISION A(LDA,*),B(*) 
!                                                                       
      DOUBLE PRECISION xDDOT,T 
      INTEGER K,KB,L,NM1 
!***FIRST EXECUTABLE STATEMENT  DGESL                                   
      NM1 = N - 1 
      IF (JOB .NE. 0) GO TO 50 
!                                                                       
!        JOB = 0 , SOLVE  A * X = B                                     
!        FIRST SOLVE  L*Y = B                                           
!                                                                       
         IF (NM1 .LT. 1) GO TO 30 
         DO 20 K = 1, NM1 
            L = IPVT(K) 
            T = B(L) 
            IF (L .EQ. K) GO TO 10 
               B(L) = B(K) 
               B(K) = T 
   10       CONTINUE 
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1) 
   20    CONTINUE 
   30    CONTINUE 
!                                                                       
!        NOW SOLVE  U*X = Y                                             
!                                                                       
         DO 40 KB = 1, N 
            K = N + 1 - KB 
            B(K) = B(K)/A(K,K) 
            T = -B(K) 
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1) 
   40    CONTINUE 
      GO TO 100 
   50 CONTINUE 
!                                                                       
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
!        FIRST SOLVE  TRANS(U)*Y = B                                    
!                                                                       
         DO 60 K = 1, N 
            T = DDOT(K-1,A(1,K),1,B(1),1) 
            B(K) = (B(K) - T)/A(K,K) 
   60    CONTINUE 
!                                                                       
!        NOW SOLVE TRANS(L)*X = Y                                       
!                                                                       
         IF (NM1 .LT. 1) GO TO 90 
         DO 80 KB = 1, NM1 
            K = N - KB 
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1) 
            L = IPVT(K) 
            IF (L .EQ. K) GO TO 70 
               T = B(L) 
               B(L) = B(K) 
               B(K) = T 
   70       CONTINUE 
   80    CONTINUE 
   90    CONTINUE 
  100 CONTINUE 
      RETURN 
      END                                           
!DECK DSCAL                                                             
      SUBROUTINE DSCAL (N, DA, DX, INCX) 
!***BEGIN PROLOGUE  DSCAL                                               
!***PURPOSE  Multiply a vector by a constant.                           
!***LIBRARY   SLATEC (BLAS)                                             
!***CATEGORY  D1A6                                                      
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)              
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR                       
!***AUTHOR  Lawson, C. L., (JPL)                                        
!           Hanson, R. J., (SNLA)                                       
!           Kincaid, D. R., (U. of Texas)                               
!           Krogh, F. T., (JPL)                                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  Subprogram                                    
!    Description of Parameters                                          
!                                                                       
!     --Input--                                                         
!        N  number of elements in input vector(s)                       
!       DA  double precision scale factor                               
!       DX  double precision vector with N elements                     
!     INCX  storage spacing between elements of DX                      
!                                                                       
!     --Output--                                                        
!       DX  double precision result (unchanged if N.LE.0)               
!                                                                       
!     Replace double precision DX by double precision DA*DX.            
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX), 
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.              
!                                                                       
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.     
!                 Krogh, Basic linear algebra subprograms for Fortran   
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.          
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791001  DATE WRITTEN                                                
!   890831  Modified array declarations.  (WRB)                         
!   890831  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900821  Modified to correct problem with a negative increment.      
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  DSCAL                                                 
      DOUBLE PRECISION DA, DX(*) 
      INTEGER I, INCX, IX, M, MP1, N 
!***FIRST EXECUTABLE STATEMENT  DSCAL                                   
      IF (N .LE. 0) RETURN 
      IF (INCX .EQ. 1) GOTO 20 
!                                                                       
!     Code for increment not equal to 1.                                
!                                                                       
      IX = 1 
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1 
      DO 10 I = 1,N 
        DX(IX) = DA*DX(IX) 
        IX = IX + INCX 
   10 END DO 
      RETURN 
!                                                                       
!     Code for increment equal to 1.                                    
!                                                                       
!     Clean-up loop so remaining vector length is a multiple of 5.      
!                                                                       
   20 M = MOD(N,5) 
      IF (M .EQ. 0) GOTO 40 
      DO 30 I = 1,M 
        DX(I) = DA*DX(I) 
   30 END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,5 
        DX(I) = DA*DX(I) 
        DX(I+1) = DA*DX(I+1) 
        DX(I+2) = DA*DX(I+2) 
        DX(I+3) = DA*DX(I+3) 
        DX(I+4) = DA*DX(I+4) 
   50 END DO 
      RETURN 
      END                                           
!DECK IDAMAX                                                            
      INTEGER FUNCTION IDAMAX (N, DX, INCX) 
!***BEGIN PROLOGUE  IDAMAX                                              
!***PURPOSE  Find the smallest index of that component of a vector      
!            having the maximum magnitude.                              
!***LIBRARY   SLATEC (BLAS)                                             
!***CATEGORY  D1A2                                                      
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)           
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR           
!***AUTHOR  Lawson, C. L., (JPL)                                        
!           Hanson, R. J., (SNLA)                                       
!           Kincaid, D. R., (U. of Texas)                               
!           Krogh, F. T., (JPL)                                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  Subprogram                                    
!    Description of Parameters                                          
!                                                                       
!     --Input--                                                         
!        N  number of elements in input vector(s)                       
!       DX  double precision vector with N elements                     
!     INCX  storage spacing between elements of DX                      
!                                                                       
!     --Output--                                                        
!   IDAMAX  smallest index (zero if N .LE. 0)                           
!                                                                       
!     Find smallest index of maximum magnitude of double precision DX.  
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)), 
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.              
!                                                                       
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.     
!                 Krogh, Basic linear algebra subprograms for Fortran   
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.          
!***ROUTINES CALLED  (NONE)                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   791001  DATE WRITTEN                                                
!   890531  Changed all specific intrinsics to generic.  (WRB)          
!   890531  REVISION DATE from Version 3.2                              
!   891214  Prologue converted to Version 4.0 format.  (BAB)            
!   900821  Modified to correct problem with a negative increment.      
!           (WRB)                                                       
!   920501  Reformatted the REFERENCES section.  (WRB)                  
!***END PROLOGUE  IDAMAX                                                
      DOUBLE PRECISION DX(*), DMAX, XMAG 
      INTEGER I, INCX, IX, N 
!***FIRST EXECUTABLE STATEMENT  IDAMAX                                  
      IDAMAX = 0 
      IF (N .LE. 0) RETURN 
      IDAMAX = 1 
      IF (N .EQ. 1) RETURN 
!                                                                       
      IF (INCX .EQ. 1) GOTO 20 
!                                                                       
!     Code for increments not equal to 1.                               
!                                                                       
      IX = 1 
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1 
      DMAX = ABS(DX(IX)) 
      IX = IX + INCX 
      DO 10 I = 2,N 
        XMAG = ABS(DX(IX)) 
        IF (XMAG .GT. DMAX) THEN 
          IDAMAX = I 
          DMAX = XMAG 
        ENDIF 
        IX = IX + INCX 
   10 END DO 
      RETURN 
!                                                                       
!     Code for increments equal to 1.                                   
!                                                                       
   20 DMAX = ABS(DX(1)) 
      DO 30 I = 2,N 
        XMAG = ABS(DX(I)) 
        IF (XMAG .GT. DMAX) THEN 
          IDAMAX = I 
          DMAX = XMAG 
        ENDIF 
   30 END DO 
      RETURN 
      END                                           




!DECK DPSORT                                                            
      SUBROUTINE DPSORT (DX, N, IPERM, KFLAG, IER) 
!***BEGIN PROLOGUE  DPSORT                                              
!***PURPOSE  Return the permutation vector generated by sorting a given 
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order. 
!            A slightly modified quicksort algorithm is used.           
!***LIBRARY   SLATEC                                                    
!***CATEGORY  N6A1B, N6A2B                                              
!***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H) 
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)                                        
!           Rhoads, G. S., (NBS)                                        
!           Wisniewski, J. A., (SNLA)                                   
!***DESCRIPTION                                                         
!                                                                       
!   DPSORT returns the permutation vector IPERM generated by sorting    
!   the array DX and, optionally, rearranges the values in DX.  DX may  
!   be sorted in increasing or decreasing order.  A slightly modified   
!   quicksort algorithm is used.                                        
!                                                                       
!   IPERM is such that DX(IPERM(I)) is the Ith value in the             
!   rearrangement of DX.  IPERM may be applied to another array by      
!   calling IPPERM, SPPERM, DPPERM or HPPERM.                           
!                                                                       
!   The main difference between DPSORT and its active sorting equivalent
!   DSORT is that the data are referenced indirectly rather than        
!   directly.  Therefore, DPSORT should require approximately twice as  
!   long to execute as DSORT.  However, DPSORT is more general.         
!                                                                       
!   Description of Parameters                                           
!      DX - input/output -- double precision array of values to be      
!           sorted.  If ABS(KFLAG) = 2, then the values in DX will be   
!           rearranged on output; otherwise, they are unchanged.        
!      N  - input -- number of values in array DX to be sorted.         
!      IPERM - output -- permutation array such that IPERM(I) is the    
!              index of the value in the original order of the          
!              DX array that is in the Ith location in the sorted       
!              order.                                                   
!      KFLAG - input -- control parameter:                              
!            =  2  means return the permutation vector resulting from   
!                  sorting DX in increasing order and sort DX also.     
!            =  1  means return the permutation vector resulting from   
!                  sorting DX in increasing order and do not sort DX.   
!            = -1  means return the permutation vector resulting from   
!                  sorting DX in decreasing order and do not sort DX.   
!            = -2  means return the permutation vector resulting from   
!                  sorting DX in decreasing order and sort DX also.     
!      IER - output -- error indicator:                                 
!          =  0  if no error,                                           
!          =  1  if N is zero or negative,                              
!          =  2  if KFLAG is not 2, 1, -1, or -2.                       
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm  
!                 for sorting with minimal storage, Communications of   
!                 the ACM, 12, 3 (1969), pp. 185-187.                   
!***ROUTINES CALLED  XERMSG                                             
!***REVISION HISTORY  (YYMMDD)                                          
!   761101  DATE WRITTEN                                                
!   761118  Modified by John A. Wisniewski to use the Singleton         
!           quicksort algorithm.                                        
!   870423  Modified by Gregory S. Rhoads for passive sorting with the  
!           option for the rearrangement of the original data.          
!   890619  Double precision version of SPSORT created by D. W. Lozier. 
!   890620  Algorithm for rearranging the data vector corrected by R.   
!           Boisvert.                                                   
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.        
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.     
!   920507  Modified by M. McClain to revise prologue text.             
!   920818  Declarations section rebuilt and code restructured to use   
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)                             
!***END PROLOGUE  DPSORT                                                
!     .. Scalar Arguments ..                                            
      INTEGER IER, KFLAG, N 
!     .. Array Arguments ..                                             
      DOUBLE PRECISION DX(*) 
      INTEGER IPERM(*) 
!     .. Local Scalars ..                                               
      DOUBLE PRECISION R, TEMP 
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN 
!     .. Local Arrays ..                                                
      INTEGER IL(21), IU(21) 
!     .. External Subroutines ..                                        
      EXTERNAL XERMSG 
!     .. Intrinsic Functions ..                                         
      INTRINSIC ABS, INT 
!***FIRST EXECUTABLE STATEMENT  DPSORT                                  
      IER = 0 
      NN = N 
      IF (NN .LT. 1) THEN 
         IER = 1 
!         CALL XERMSG ('SLATEC', 'DPSORT',                               &
!     &    'The number of values to be sorted, N, is not positive.',     &
!     &    IER, 1) 
         write(6,*)"DPSORT: The number of values to be sorted, N, is not positive."                                                      
         RETURN 
      ENDIF 
!                                                                       
      KK = ABS(KFLAG) 
      IF (KK.NE.1 .AND. KK.NE.2) THEN 
         IER = 2 
!         CALL XERMSG ('SLATEC', 'DPSORT',                               &
!     &    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!     &    IER, 1)   
   
          write(6,*)"DPSORT: The sort control parameter, KFLAG, is not 2, 1, -1, or -2."                                                 
         RETURN 
      ENDIF 
!                                                                       
!     Initialize permutation vector                                     
!                                                                       
      DO 10 I=1,NN 
         IPERM(I) = I 
   10 END DO 
!                                                                       
!     Return if only one value is to be sorted                          
!                                                                       
      IF (NN .EQ. 1) RETURN 
!                                                                       
!     Alter array DX to get decreasing order if needed                  
!                                                                       
      IF (KFLAG .LE. -1) THEN 
         DO 20 I=1,NN 
            DX(I) = -DX(I) 
   20    CONTINUE 
      ENDIF 
!                                                                       
!     Sort DX only                                                      
!                                                                       
      M = 1 
      I = 1 
      J = NN 
      R = .375D0 
!                                                                       
   30 IF (I .EQ. J) GO TO 80 
      IF (R .LE. 0.5898437D0) THEN 
         R = R+3.90625D-2 
      ELSE 
         R = R-0.21875D0 
      ENDIF 
!                                                                       
   40 K = I 
!                                                                       
!     Select a central element of the array and save it in location L   
!                                                                       
      IJ = I + INT((J-I)*R) 
      LM = IPERM(IJ) 
!                                                                       
!     If first element of array is greater than LM, interchange with LM 
!                                                                       
      IF (DX(IPERM(I)) .GT. DX(LM)) THEN 
         IPERM(IJ) = IPERM(I) 
         IPERM(I) = LM 
         LM = IPERM(IJ) 
      ENDIF 
      L = J 
!                                                                       
!     If last element of array is less than LM, interchange with LM     
!                                                                       
      IF (DX(IPERM(J)) .LT. DX(LM)) THEN 
         IPERM(IJ) = IPERM(J) 
         IPERM(J) = LM 
         LM = IPERM(IJ) 
!                                                                       
!        If first element of array is greater than LM, interchange      
!        with LM                                                        
!                                                                       
         IF (DX(IPERM(I)) .GT. DX(LM)) THEN 
            IPERM(IJ) = IPERM(I) 
            IPERM(I) = LM 
            LM = IPERM(IJ) 
         ENDIF 
      ENDIF 
      GO TO 60 
   50 LMT = IPERM(L) 
      IPERM(L) = IPERM(K) 
      IPERM(K) = LMT 
!                                                                       
!     Find an element in the second half of the array which is smaller  
!     than LM                                                           
!                                                                       
   60 L = L-1 
      IF (DX(IPERM(L)) .GT. DX(LM)) GO TO 60 
!                                                                       
!     Find an element in the first half of the array which is greater   
!     than LM                                                           
!                                                                       
   70 K = K+1 
      IF (DX(IPERM(K)) .LT. DX(LM)) GO TO 70 
!                                                                       
!     Interchange these elements                                        
!                                                                       
      IF (K .LE. L) GO TO 50 
!                                                                       
!     Save upper and lower subscripts of the array yet to be sorted     
!                                                                       
      IF (L-I .GT. J-K) THEN 
         IL(M) = I 
         IU(M) = L 
         I = K 
         M = M+1 
      ELSE 
         IL(M) = K 
         IU(M) = J 
         J = L 
         M = M+1 
      ENDIF 
      GO TO 90 
!                                                                       
!     Begin again on another portion of the unsorted array              
!                                                                       
   80 M = M-1 
      IF (M .EQ. 0) GO TO 120 
      I = IL(M) 
      J = IU(M) 
!                                                                       
   90 IF (J-I .GE. 1) GO TO 40 
      IF (I .EQ. 1) GO TO 30 
      I = I-1 
!                                                                       
  100 I = I+1 
      IF (I .EQ. J) GO TO 80 
      LM = IPERM(I+1) 
      IF (DX(IPERM(I)) .LE. DX(LM)) GO TO 100 
      K = I 
!                                                                       
  110 IPERM(K+1) = IPERM(K) 
      K = K-1 
      IF (DX(LM) .LT. DX(IPERM(K))) GO TO 110 
      IPERM(K+1) = LM 
      GO TO 100 
!                                                                       
!     Clean up                                                          
!                                                                       
  120 IF (KFLAG .LE. -1) THEN 
         DO 130 I=1,NN 
            DX(I) = -DX(I) 
  130    CONTINUE 
      ENDIF 
!                                                                       
!     Rearrange the values of DX if desired                             
!                                                                       
      IF (KK .EQ. 2) THEN 
!                                                                       
!        Use the IPERM vector as a flag.                                
!        If IPERM(I) < 0, then the I-th value is in correct location    
!                                                                       
         DO 150 ISTRT=1,NN 
            IF (IPERM(ISTRT) .GE. 0) THEN 
               INDX = ISTRT 
               INDX0 = INDX 
               TEMP = DX(ISTRT) 
  140          IF (IPERM(INDX) .GT. 0) THEN 
                  DX(INDX) = DX(IPERM(INDX)) 
                  INDX0 = INDX 
                  IPERM(INDX) = -IPERM(INDX) 
                  INDX = ABS(IPERM(INDX)) 
                  GO TO 140 
               ENDIF 
               DX(INDX0) = TEMP 
            ENDIF 
  150    CONTINUE 
!                                                                       
!        Revert the signs of the IPERM values                           
!                                                                       
         DO 160 I=1,NN 
            IPERM(I) = -IPERM(I) 
  160    CONTINUE 
!                                                                       
      ENDIF 
!                                                                       
      RETURN 
      END                                           



end module lmfit_neu

!9 ! test program
!9 program testprog
!9 use lmfit
!9 implicit none
!9 
!9 integer, parameter  :: NP = 200
!9 integer, parameter  :: ME = 30
!9 double precision    :: t0 = 0.001d0       ! start of log list
!9 double precision    :: t1 = 1000d0        ! end   of log list
!9 double precision    :: xv(NP), yv(NP), ydiff(NP)
!9 double precision    :: a(ME) , r(ME)
!9 double precision    :: chisq
!9 integer             :: i, ier
!9 double precision    :: ts
!9 
!9 double precision    :: t_table_spacing1 = 3d0 
!9 double precision    :: t_table_spacing2 = 1d0 
!9 ! exp(i*log(tmax*300d0)/nxpoints)/100d0
!9 ! log t table, first point at exp(log(tmax*t_table_spacing1)/nxpoints)/t_table_spacing2 
!9 !              last  point at exp(log(tmax*t_table_spacing1))/t_table_spacing2 = tmax*(t_table_spacing1/t_table_spacing2) 
!9 
!9 integer :: ne = 8
!9 integer :: ne1
!9 
!9 
!9  NEXP_DEVIATION_TOLERANCE = 1d-3
!9 
!9  
!9 !7  CALL PREPARE_TTABLE(t0, t1, np, xv)
!9 !7  
!9 !7  write(*,*)"Initial table:"
!9 !7  do i=1,np
!9 !7  !  ts      =  exp(i*log(t1*t_table_spacing1/t0)/np)/(t_table_spacing2/t0)    !! TBD new here is the table spacing
!9 !7  !  xv(i)   =  ts
!9 !7    yv(i)   =  testfun(xv(i))
!9 !7  enddo
!9 !7       
!9 !7  
!9 !7  
!9 !7  
!9 !7  !! TBD new >>>  
!9 !7  ne1 = 3   
!9 !7  a(1:ne1) = 1d0/3d0
!9 !7  r(2)     = 1d0/sum( (yv(1:np-1)-yv(np)) * ( xv(2:np)-xv(1:np-1) ) )
!9 !7  r(1)     = sqrt(r(2)/xv(1))
!9 !7  r(3)     = sqrt(r(2)/(xv(np)/3))  
!9 !7  
!9 !7  write(*,*)"Initial0 params a:",a(1:ne1)
!9 !7  write(*,*)"Initial0 params r:",r(1:ne1)
!9 !7  
!9 !7  call nexp_match2(xv,yv,np,ne1,a,r,chisq)
!9 !7  
!9 !7  ! check
!9 !7  do i=1,np
!9 !7    ydiff(i) = yv(i) - sum(a(1:ne1)*exp(-xv(i)*r(1:ne1)))
!9 !7    write(*,'(i5,f16.9,2x,f16.7,4x,f16.12)') i, xv(i), yv(i), ydiff(i)
!9 !7  enddo
!9 !7  write(*,*)"Initial3 params a:",a(1:ne1)
!9 !7  write(*,*)"Initial3 params r:",r(1:ne1)
!9 !7  
!9 !7  
!9 !7  
!9 !7  
!9 !7  DO while( 2*ne1 - 1   <= ME .and. maxval(abs(ydiff(1:np))) > 1d-3)
!9 !7  
!9 !7  do i=ne1,1,-1
!9 !7    a(2*i-1) = a(i)
!9 !7    r(2*i-1) = r(i)
!9 !7  enddo
!9 !7  ne1 = ne1+(ne1-1)
!9 !7  do i=2,ne1,2
!9 !7    a(i) = 1d-3
!9 !7    r(i) = sqrt(r(i-1)*r(i+1))
!9 !7  enddo
!9 !7  call nexp_match2(xv,yv,np,ne1,a,r,chisq)    ! TBD INCLUDE
!9 !7  
!9 !7  
!9 !7  
!9 !7  ! check
!9 !7  do i=1,np
!9 !7    ydiff(i) = yv(i) - sum(a(1:ne1)*exp(-xv(i)*r(1:ne1)))
!9 !7    write(*,'(i5,f16.9,2x,f16.7,4x,f16.12)') i, xv(i), yv(i), ydiff(i)
!9 !7  enddo
!9 !7  
!9 !7  
!9 !7  write(*,*)"NEXP,      AEXP             REXP   "
!9 !7  write(*,'(i5,f16.9,2x,f16.7)')(i, a(i), r(i), i=1,ne1)
!9 !7  
!9 !7  write(*,*)
!9 !7  write(*,'("Nexp=",i3,"  SSQ=",f12.9,"  Delta= ",2f12.9," suma=",f12.9)') &
!9 !7            ne1, chisq, minval(ydiff), maxval(ydiff), sum(a(1:ne1))
!9 !7  
!9 !7  write(*,*)
!9 !7  write(*,*)"================================================================================"
!9 !7  write(*,*)
!9 !7  
!9 !7  ENDDO
!9 !7  
!9 !7  
!9 !7 
!9 
!9 
!9 
!9  CALL PREPARE_TTABLE(t0, t1, np, xv)
!9 
!9  do i=1,np
!9   yv(i)   =  testfun(xv(i))
!9   write(*,'(i5,f16.9,2x,f16.7)') i, xv(i), yv(i)
!9  enddo
!9   
!9 
!9 
!9 call match_exp(xv,yv,np,me,ne1,a,r,chisq)
!9 write(*,*)"ME: NEXP,      AEXP             REXP   "
!9 write(*,'(i5,f16.9,2x,f16.7)')(i, a(i), r(i), i=1,ne1)
!9 
!9 write(*,*)
!9 write(*,'("Nexp=",i3,"  Delta=",f12.9,"   suma=",f12.9)') &
!9           ne1, chisq, sum(a(1:ne1))
!9 
!9 write(*,*)
!9 write(*,*)"================================================================================"
!9 write(*,*)
!9 
!9 
!9 
!9 
!9 !! << TBD new end
!9 !3 enddo
!9 
!9 CONTAINS
!9 
!9  function testfun(t) result(y)
!9  implicit none
!9  double precision, intent(in) :: t
!9  double precision             :: y
!9  double precision             :: a1=0.3d0,  beta1=0.5d0, r1=1d0
!9  double precision             :: a2=0.5d0,  beta2=0.8d0, r2=20d0
!9 
!9  y = a1 * exp(-(t*r1)**beta1) +  &
!9      a2 * exp(-(t*r2)**beta2) +  &
!9      (1d0-a1-a2)
!9 
!9 end function testfun
!9 
!9 
!9 
!9 end program testprog