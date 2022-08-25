      function th37(x,pa,thnam,parnam,npar,ini)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ========================================
!
!                          -------> cloizeaux <--------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)


         CHARACTER*8     thnam,parnam(20)
         REAL            th37, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision w,d,q,n,ne,cloizeaux,t,a,wl4


 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'cloizeaux'
         nparx =  3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th37   = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitud'          
         parnam(2) = 'wl4     '
         parnam(3) = 'dtube   '
c
         th37  = 0
         return
       endif
c
c ---- calculate theory here -----
        a       = pa(1)   ! amplitude (should be 1)
        wl4     = pa(2)   ! Rouse rate W*l**4
        d       = pa(3)   ! Tube diameter

        call parget('q       ',xh ,iadda,ier)
        if(ier.eq.0) then
          q       = xh
        else
          q       = 1.0  
        endif
 
        t = x 
 
        th37 = cloizeaux(t,q,d,Wl4)

        th37 = th37 * a

        return
        end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formel (23) aus Cloizeaux, J. de Physique I, 3 (1993) p.1533 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function cloizeaux(t,q,d,Wl4)

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-2 , gwidth=4.0d0, maxit=1000)

      double precision FA_ker, FB_ker, FYZ

      double precision d, q, Wl4, t, cadapint, gamma
      double precision go,I1,I2,result,n1,p,dp
      integer          j        
      double precision QFIX(200)
      double precision QX(200)
      double precision WX(200)
      double precision alpha, beta 

      external FYZ    
      
      save     qx,wx,n

      data n/0/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(n.eq.0) then
       call erset (0,0,0)  ! supress imsl error messages erf underflow
       write(6,*)'Creating Gauss-Laguerre Coeffs' 
       alpha = 0.0d0
       beta  = 1.0d0

       nfix  = 0
       iweight = 6

       n = 6

       call DGQRUL(n,iweight, alpha, beta, nfix, qfix, qx, wx)
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       gamma = Wl4/9

       Zav = q*q*d*d/2
       Z   = Zav
       Y   = q*q*sqrt(gamma*t)

        i1 = 0.0d0
        n1 = 0.0d0
        dp = 0.2d0
        do j=1,n
          p = qx(j)
          i1 = i1 + FYZ(p)*wx(j)
        enddo
        
      result     =   i1 + log(1+Zav)/Zav
      cloizeaux = result
      write(6,*)t,result
      return

      end




      double precision function FYZ(p)

      implicit none 

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-2, gwidth=4.0d0, maxit=1000)

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      double precision p, badapint, gu, go, fa, fb
      double precision FA_ker, FB_ker, FAB_ker
      double precision upbound, lowbound      
      double precision result

      external FA_ker, FB_ker, FAB_ker      
      external upbound, lowbound


      data Faf/1.0d0/
      data Fbf/1.0d0/

      Z = p*Zav

      gu = 0.0d0
      go = 1.0d0
    
       call dtwodq(FAB_ker,gu,go,lowbound,upbound,
     *            epsilon, epsilon*10, 2, result, errac)
       Fyz = result * Z
!       write(6,*)' 2: ',Fyz, errac

      return
      end


      double precision function lowbound(x)
      implicit none
      double precision x
      
      lowbound = x
      return
      end

      double precision function upbound(x)
      implicit none
      double precision x
      
      upbound = 1.0d0
      return
      end



      double precision function FA_ker(a)
      
      implicit none
      double precision a, arg, derfc, ssap, ssa, ssaq, u1, arg2
      integer          p, plim

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-12, gwidth=4.0d0, maxit=1000)
      double precision edwp   
      parameter(edwp=0.564189583d0)        !! = 1/sqrt(pi)


      ssa = 0      
      plim = 0.5d0*gwidth*Y/Z + 2   

!      write(6,*)'Fa plim:',plim, ' a=',a
      do p=-plim,plim
        u1   = abs((a-2*p)*Z)
        arg2 = -((u1/Y)**2)
        if(arg2.lt.-150.0d0) arg2=-150.0d0
        ssap = edwp*exp(arg2)*Y - u1*derfc(u1/Y) 
        ssa  = ssa + ssap
!        write(6,*)p,ssa,ssap
      enddo

      arg = -Z*a - ssa
      FA_ker = exp(arg) * Faf

!      write(6,*)'FA_ker=', FA_ker      

      return
      end

      double precision function FB_ker(b)

      implicit none
      double precision b, arg, derfc, ssbp, ssb, ssbq, u1, arg2
      integer          p, plim

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-12, gwidth=4.0d0, maxit=1000)
      double precision edwp   
      parameter(edwp=0.564189583d0)        !! = 1/sqrt(pi)


      ssb = 0      
      plim = 0.5d0*gwidth*Y/Z + 2   


!      write(6,*)'Fb plim:',plim, ' b=',b

      do p=-plim,plim
        u1   = abs((b-2*p)*Z)
        arg2 = -((u1/Y)**2)
        if(arg2.lt.-150.0d0) arg2=-150.0d0
        ssbp = edwp*exp(arg2)*Y - u1*derfc(u1/Y) 
        ssb  = ssb + ssbp
!       write(6,*)p,ssb,ssbp
      enddo

      arg = ssb
      FB_ker = exp(arg) * Fbf


!      write(6,*)'FB_ker=', FB_ker      

      return
      end




      double precision function FAB_ker(a,b)
      
      implicit none
      double precision a, arg, derfc, ssap, ssa, ssaq, u1, arg2
      double precision b,u2,arg22, ssbp
      integer          p, plim, nn

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-12, gwidth=4.0d0, maxit=1000)
      double precision edwp   
      parameter(edwp=0.564189583d0)        !! = 1/sqrt(pi)


      ssa = 0      
      plim = 0.5d0*gwidth*Y/Z + 2   

!      write(6,*)'Fa plim:',plim, ' a=',a
      do p=-plim,plim
        u1   = abs((a-2*p)*Z)
        arg2 = -((u1/Y)**2)
        if(arg2.lt.-150.0d0) arg2=-150.0d0
        ssap = edwp*exp(arg2)*Y - u1*derfc(u1/Y) 
  
        u2   = abs((b-2*p)*Z)
        arg22 = -((u2/Y)**2)
        if(arg22.lt.-150.0d0) arg22=-150.0d0
        ssbp = edwp*exp(arg22)*Y - u2*derfc(u2/Y) 
        ssa  = ssa + ssap - ssbp
!        write(6,*)p,ssa,ssap
      enddo

      arg = -Z*a - ssa
      FAB_ker = exp(arg) 

!      write(6,*)'FAB_ker=', FAB_ker      

      return
      end




      double precision function icheck(u)
      double precision aa, bb , u
      common/cic/aa,bb


      icheck = exp(-(aa/u)**2)
    

      return
      end













       function bgaus8a (f,xu,xo)
c-----------------------------------------------------------------------
c      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
c-----------------------------------------------------------------------
       parameter ( ndim = 8 )
       parameter ( ndim2=ndim/2)
       implicit real*8 (a-h,o-z)
       dimension a(ndim2),x(ndim2)
       data a / 0.362683783378362d0,
     *          0.313706645877887d0,
     *          0.222381034453374d0,
     *          0.101228536290376d0/
       data x / 0.183434642495650d0,
     *          0.525532409916329d0,
     *          0.796666477413627d0,
     *          0.960289856497536d0/
 
       xave = (xo+xu)*0.5d0
       range= (xo-xu)*0.5d0
       sum = 0.d0
       do i=1,ndim2
         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) )
       enddo
       bgaus8a = sum * range
 
       return
       end


       function badapint(f,a,b,epsilon0,maxiter0,erroraccu)
c      ===================================================
c
c      lokal adaptives integrationsverfahren 2te ver. wg. rekur.
c
       implicit real*8 (a-h,o-z)
c
       parameter( maxstack = 30)
c                 --------------->  stacktiefe
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack)
c
       logical cray
       real*4 xxxx, yyyy, ptxf                          !! aix
       common/outlev/iot,ibild,ierrs,inka, iibuf, xxxx,yyyy, ptxf(20) !! aix
c
!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!
       character*1 answer
       common/cailff/cmue,cbeta,crres   !!!! test only !!!
!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!
 
       external f
c
c      -------------------------------------------------------------------
         ifrom            = erroraccu+0.1
         maxiter          = 10000
         epsilon = epsilon0
          erroraccu        = 0.d0
         iterationcounter = 0
         erg              = 0.0d0
         itop      = 0
         xa(itop)  = a
         xb(itop)  = b
         s(itop)   = bgaus8a( f, a, b )
 
         do i=1,maxiter
           if (itop.ge.(maxstack-2)) then
              erg = erg + s(itop)
              xbb = xb(itop)
              itm = itop-1
              do itop=itm,0,-1
                 if(xa(itop).eq.xbb) goto 1
              enddo
1             continue
             epsilon = 5*epsilon 
             write(6,*)'warning! badaptint stack overflow!'
           else
              iterationcounter = iterationcounter + 1
              itop             = itop +1
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0
              xb(itop) =  xb(itop-1)
              s(itop)  = bgaus8a( f, xa(itop), xb(itop) )
              itop             = itop +1
              xa(itop) =  xa(itop-2)
              xb(itop) =  xa(itop-1)
              s(itop)  = bgaus8a( f, xa(itop), xb(itop) )
              error   =  dabs((s(itop)+s(itop-1))-s(itop-2))
 
              if (iot.gt.2) then
                 write(6,'(1x,i3,i5,4e13.6)')itop,iterationcounter,
     *                 xa(itop),xb(itop),(s(itop)+s(itop-1)),error
              endif
 
              if (error.lt.epsilon) then
                 erg = erg + s(itop)+s(itop-1)
                 erroraccu = erroraccu + error
                 itop = itop-1
                 xbb = xb(itop)
                 itop = itop-1
                 itm = itop-1
                 do itop=itm,0,-1
                    if(xa(itop).eq.xbb) goto 2
                 enddo
2                continue
              endif
            endif
            if (itop.le.0) goto 3
         enddo
         write(6,*)                                                     '
     *    'badapint(',ifrom,') fatal error iterationnumber exceeded!'
         write(6,*) 'limits a=',a
         write(6,*) 'limits b=',b
         write(6,*) 'cmue    =',cmue
         write(6,*) 'cbeta   =',cbeta
         write(6,*) 'crres   =',crres
         write(6,*) 'W(eiter), F(ile), A(bbruch) ---> ?'
         read(5,'(A)') answer
         if(answer.eq.'a') stop
         if(answer.eq.'f') then
           open(22,file='test.out',status='UNKNOWN')
           nn = 100
           dx = (b-a)/nn
           do i=0,100
            xx = a + i*dx
            yy = f(xx)
            write(22,*) xx,yy
           enddo
           close(22)
           write(6,*)'integrand written to test.out'
         endif
         badapint = 9.99999d33
         return
3        continue
 
 
         badapint = erg
 
         return
 
       end

!! ------------------------------------------------------------------------
!! IMSL Name:  TWODQ/DTWODQ (Single/Double precision version)
!!  
!! Revised:    May 31, 1991
!!  
!! Purpose:    Compute a two-dimensional iterated integral.
!!  
!! Usage:      CALL TWODQ (F, A, B, G, H, ERRABS, ERRREL, IRULE, RESULT,
!!                         ERREST)
!!  
!! Arguments:
!!    F      - User-supplied FUNCTION to be integrated.  The form is
!!             F(X, Y), where
!!             X      - First argument of F.  (Input)
!!             Y      - Second argument of F.  (Input)
!!             F      - The function value.  (Output)
!!             F must be declared EXTERNAL in the calling program.
!!    A      - Lower limit of outer integral.  (Input)
!!    B      - Upper limit of outer integral.  (Input)
!!    G      - User-supplied FUNCTION to evaluate the lower limits of
!!             the inner integral.  The form is
!!             G(X), where
!!             X      - Only argument of G.  (Input)
!!             G      - The function value.  (Output)
!!             G must be declared EXTERNAL in the calling program.
!!    H      - User-supplied FUNCTION to evaluate the upper limits of
!!             the inner integral.  The form is
!!             H(X), where
!!             X      - Only argument of H.  (Input)
!!             H      - The function value.  (Output)
!!             H must be declared EXTERNAL in the calling program.
!!    ERRABS - Absolute accuracy desired.  (Input)
!!    ERRREL - Relative accuracy desired.  (Input)
!!    IRULE  - Choice of quadrature rule.  (Input)
!!             The Gauss-Kronrod rule is used with the following
!!             points:
!!              IRULE  Points
!!                1      7-15
!!                2     10-21
!!                3     15-31
!!                4     20-41
!!                5     25-51
!!                6     30-61
!!             If the function has a peak singularity, use IRULE = 1.
!!             If the function is oscillatory, use IRULE = 6.
!!    RESULT - Estimate of the integral from A to B of F.  (Output)
!!    ERREST - Estimate of the absolute value of the error.  (Output)
!!  
!! Remarks:
!! 1. Automatic workspace usage is
!!             TWODQ    2500 units, or
!!             DTWODQ   4500 units.
!!    Workspace may be explicitly provided, if desired, by use of
!!    T2ODQ/DT2ODQ.  The reference is
!!             CALL T2ODQ (F, A, B, G, H, ERRABS, ERRREL, IRULE,
!!                         RESULT, ERREST, MAXSUB, NEVAL, NSUBIN,
!!                         ALIST, BLIST, RLIST, ELIST, IORD, WK, IWK)
!!    The additional arguments are as follows:
!!    MAXSUB - Number of subintervals allowed.  (Input)
!!             A value of 250 is used by TWODQ.
!!    NEVAL  - Number of evaluations of F.  (Output)
!!    NSUBIN - Number of subintervals generated in the outer integral.
!!             (Output)
!!    ALIST  - Array of length MAXSUB containing a list of the NSUBIN
!!             left endpoints for the outer integral.  (Output)
!!    BLIST  - Array of length MAXSUB containing a list of the NSUBIN
!!             right endpoints for the outer integral.  (Output)
!!    RLIST  - Array of length MAXSUB containing approximations to the
!!             NSUBIN integrals over the intervals defined by ALIST,
!!             BLIST, pertaining only to the outer integral.  (Output)
!!    ELIST  - Array of length MAXSUB containing the error estimates
!!             of the NSUBIN values in RLIST.  (Output)
!!    IORD   - Array of length MAXSUB.  (Output)
!!             Let K be NSUBIN if NSUBIN .LE. (MAXSUB/2+2),
!!             MAXSUB+1-NSUBIN otherwise.  Then the first fortranK
!!             locations contain pointers to the error estimates over
!!             the corresponding subintervals, such that
!!             ELIST(IORD(1)), ..., ELIST(IORD(K)) form a decreasing
!!             sequence.
!!    WK     - Work array of length 4*MAXSUB, needed to evaluate the
!!             inner integral.
!!    IWK    - Work array of length MAXSUB, needed to evaluate the
!!             inner integral.
!!  
!! 2. Informational errors
!!    Type Code
!!      4   1  The maximum number of subintervals allowed has been
!!             reached.
!!      3   2  Roundoff error, preventing the requested tolerance from
!!             being achieved, has been detected.
!!      3   3  A degradation in precision has been detected.
!!  
!! 3. If EXACT is the exact value, TWODQ attempts to find RESULT such
!!    that ABS(EXACT-RESULT) .LE. MAX(ERRABS,ERRREL*ABS(EXACT)).  To
!!    specify only a relative error, set ERRABS to zero.  Similarly, to
!!    specify only an absolute error, set ERRREL to zero.
!!  
!! Keyword:    Quadrature
!!  
!! GAMS:       H2b1a1
!!  
!! Chapter:    MATH/LIBRARY Integration and Differentiation
!!  
!! Page No.:   MATH/LIBRARY User's Manual page 715
!!  
!! ------------------------------------------------------------------------
!! IMSL Name:  GQRUL/DGQRUL (Single/Double precision version)
!!  
!! Revised:    April 19, 1991
!!  
!! Purpose:    Compute a Gauss, Gauss-Radau, or Gauss-Lobatto
!!             quadrature rule with various classical weight
!!             functions.
!!  
!! Usage:      CALL GQRUL (N, IWEIGH, ALPHA, BETA, NFIX, QXFIX, QX, QW)
!!  
!! Arguments:
!!    N      - Number of quadrature points.  (Input)
!!    IWEIGH - Index of the weight function.  (Input)
!!             IWEIGH  WT(X)             Interval     Name
!!                1    1                 (-1,+1)      Legendre
!!                2    1/SQRT(1-X**2)    (-1,+1)      Chebyshev 1st kind
!!                3    SQRT(1-X**2)      (-1,+1)      Chebyshev 2nd kind
!!                4    EXP(-X**2)        (-inf,+inf)  Hermite
!!                5    (1-X)**ALPHA*
!!                     (1+X)**BETA       (-1,+1)      Jacobi
!!                6    EXP(-X)*X**ALPHA  (0,+inf)     Generalized
!!                                                    Laguerre
!!                7    1/COSH(X)         (-inf,+inf)  COSH
!!    ALPHA  - Parameter used in the weight function with some values
!!             of IWEIGH, otherwise it is ignored.  (Input)
!!    BETA   - Parameter used in the weight function with some values
!!             of IWEIGH, otherwise it is ignored.  (Input)
!!    NFIX   - Number of fixed quadrature points.  (Input)
!!             NFIX = 0, 1 or 2.  For the usual Gauss quadrature rules,
!!             NFIX = 0.
!!    QXFIX  - Array of length NFIX (ignored if NFIX = 0) containing
!!             the preset quadrature point(s).  (Input)
!!    QX     - Array of length N containing quadrature points.
!!             (Output)
!!    QW     - Array of length N containing quadrature weights.
!!             (Output)
!!  
!! Remarks:
!! 1. Automatic workspace usage is
!!             GQRUL    N units, or
!!             DGQRUL   2*N units.
!!    Workspace may be explicitly provided, if desired, by use of
!!    G2RUL/DG2RUL.  The reference is
!!             CALL G2RUL (N, IWEIGH, ALPHA, BETA, NFIX, QXFIX,
!!                         QX, QW, WK)
!!    The additional argument is
!!    WK     - Work array of length N.
!!  
!! 2. If IWEIGH specifies the weight WT(X) and the interval (a,b), then
!!    approximately
!!               b                  N
!!         INTEGRAL F(X)*WT(X) dX = SUM F(QX(I))*QW(I)
!!               a                 I=1
!!  
!! 3. Gaussian quadrature is always the method of choice when the
!!    function F(X) behaves like a polynomial.  Gaussian quadrature is
!!    also useful on infinite intervals (with appropriate weight
!!    functions), because other techniques often fail.
!!  
!! 4. The weight function 1/COSH(X) behaves like a polynomial near zero
!!    and like EXP(-ABS(X)) far from zero.
!!  
!! Keywords:   Univariate quadrature; Legendre; Chebyshev; Hermite;
!!             Jacobi; Laguerre; Numerical integration
!!  
!! GAMS:       H2c
!!  
!! Chapter:    MATH/LIBRARY Integration and Differentiation
!!  
!! Page No.:   MATH/LIBRARY User's Manual page 723
!!  
