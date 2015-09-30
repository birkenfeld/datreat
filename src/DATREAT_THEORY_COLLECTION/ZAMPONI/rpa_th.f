      function th15(x,pa,thnam,parnam,npar,ini)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                     -------> rpa <-------                            !
!                                                                      !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
!!
!! Umsetzung RPA siehe z.B. Vilgis et al. Macromolecules Vol. 24 p4481 (1991)
!!
!!

      implicit none
      
      double precision Pi
      Parameter       (Pi=3.141592654d0)
      
      
      CHARACTER*8     thnam,parnam(20)
      REAL            th15, x, pa, qq, zpi, xh, vol_frac         
      INTEGER         ini, npar, nparx
      DIMENSION       pa(20),qq(3)
      DATA            zpi/6.283185/
      
      double precision q,n1,n2,n3,rpa,leff,incoh1,a1,a2,incoh2
      real frac

c aus rpa.f
      integer maxd
      parameter (maxd=15)    ! maximum dimension

      double precision S0 (maxd,maxd)
      double precision S  (maxd,maxd)
      double precision V  (maxd,maxd)
      double precision S0i(maxd,maxd)
      double precision Si (maxd,maxd)
      double precision CHI(maxd,maxd) 
c      double precision B  (maxd)
      real B  (maxd)

c      integer nblock,ng(2,maxd),n
      integer nblock
      double precision ng(2,maxd),n
      integer ndim
      real b1,b2,b3,b4
      double precision Sab,intSab
      double precision l_eff
      double precision v0
      double precision Phi0
      double precision sum
      integer i,j,l
c      
             
      INTEGER iadda, ier, iout, iot                
      COMMON/thiadd/iadda       !! transfer of address for parget

c
c ----- initialisation -----
      if(ini.eq.0) then
         thnam = 'rpa '

         nparx =  8
c         nparx =  7

         if(npar.lt.nparx) then
            write(6,1)thnam,nparx,npar
 1          format(' theory: ',a8,' no of parametrs=',i8,
     *           ' exceeds current max. = ',i8)
            th15   = 0
            return
         endif
         npar = nparx


c     --------------> set the number of parameters
         parnam(1) = 'ampl1   '          
         parnam(2) = 'ampl2   '          
         parnam(3) = 'l       '
         parnam(4) = 'n1      '
         parnam(5) = 'n2      '
         parnam(6) = 'n3      '
         parnam(7) = 'inc_bk1 '
         parnam(8) = 'inc_bk2 '
c         parnam(6) = 'monofrac'

c     
          th15  = 0
         return
      endif
c     
c     ---- calculate theory here -----
      a1      = abs(pa(1))      ! amplitude
      a2      = abs(pa(2))      ! amplitude
      leff    = abs(pa(3))      ! segment length
      n1      = abs(pa(4))      ! first block  
      n2      = abs(pa(5))      ! second block
      n3      = abs(pa(6))      ! third block
      incoh1  = pa(7)      ! incoh. bkgr 
      incoh2  = pa(8)      ! incoh. bkgr 
c      frac    = pa(6)           ! fraction of homopolymer (matrix)

      
      call parget('frac    ',frac,iadda,ier)
      phi0 = frac 

      q = x 
      
      
! 0. Limit Parameter
      v0 = 1d5                          ! v0 --> infinity to realize excluded volume

      nblock = 3
      ndim   = nblock+1 

      ng(1,1) = 1                       ! Block limits in terms of segment numbers
      ng(2,1) = n1
      ng(1,2) = n1+1
      ng(2,2) = n1+n2
      ng(1,3) = n1+n2+1
      ng(2,3) = n1+n2+n3
      n       = ng(2,3)                 ! unconnected block of total length (matrix polymer)
c      ng(1,1) = 1                       ! Block limits in terms of segment numbers
c      ng(2,1) = n1
c      ng(1,2) = n1+1
c      ng(2,2) = n2
c      ng(1,3) = n2+1
c      ng(2,3) = n3
c      n       = ng(2,3)                 ! unconnected block of total length (matrix polymer)


      l_eff   = leff/sqrt(6.0d0)


      do i=1,4
         do j=1,4
            if(i.le.3 .and. j.le.3) then
c               S(i,j) = Sab(q,l_eff,nint(n),nint(ng(1,i)),nint(ng(2,i)),
c     #              nint(ng(1,j)),nint(ng(2,j)))*(1-Phi0)
              S(i,j) = intSab(q,l_eff,n,ng(1,i),ng(2,i),ng(1,j),ng(2,j))
     *              *(1-Phi0)
            else
               S(i,j) = 0
            endif
            V(i,j) = v0
         enddo
      enddo
c      S(4,4) = Sab(q,l_eff,nint(n),1,nint(n),1,nint(n))*Phi0 
      S(4,4) = intSab(q,l_eff,n,1.d0,n,1.d0,n)*Phi0 
      
! call IMSL inversion of S
 
      CALL DLINRG (ndim, S, maxd, Si, maxd)
 
! Add V
      do i=1,ndim
         do j=1,ndim
            S0i(i,j) = Si(i,j) + V(i,j)
         enddo
      enddo
   
! Add Chi (to be done if needed)


! and invert it to get S0
      CALL DLINRG (ndim, S0i, maxd, S0, maxd)
      
! and check the symmetry 
      B(1)  = 0
      B(2)  = 1
      B(3)  = 0
      B(4)  = 0 
c
      call parget ('b1      ',B(1),iadda,ier)
      call parget ('b2      ',B(2),iadda,ier)
      call parget ('b3      ',B(3),iadda,ier)
      call parget ('b4      ',B(4),iadda,ier)
c
      sum = 0
      do i=1,ndim
         do j=1,ndim
            sum = sum + B(i)*S0(i,j)*B(j)
         enddo
      enddo        

c      th15 = S0(2,2)
      th15 = sum
      
c      th15 = th15 * a1 + incoh1

      if(frac.gt.0.001) then
         th15 = th15 * a1 + incoh1
      else
         th15 = th15 * a2 + incoh2
      endif
           

 999  continue

      return
      end
      

c ==============================================================================

      double precision function Sab(q,l_eff,N,na1,na2,nb1,nb2)
!     --------------------------------------------------------
!
! see e.g. Sakurai et al. Macromolecules Vol. 25 p2679 (1992)
!
! partial structure functions of an multiblock-polymer
! where the a-block extends from segments na1...na2
! and   the b-block         from segments nb1...nb2
!
! the total number of segments is N
!
! the effective segment length of a is l_eff
!
! such that R_a(i,j) = (l_a**2)*|i-j|
!
       implicit none 

       double precision q,l_eff
       integer N,na1,na2,nb1,nb2
       
       double precision x,sum
       integer i,j


       sum = 0
       do i=na1,na2
        do j=nb1,nb2
          x = -abs(i-j)*(q*l_eff)**2
          sum = sum+exp(x)  
        enddo
       enddo
       
       Sab = sum/N

       return
       end


c ==============================================================================

      double precision function intSab(q,l_eff,N,na1,na2,nb1,nb2)

      implicit none 

      double precision q,l_eff
      double precision N,na1,na2,nb1,nb2,Heaviside
      double precision x11,x12,x21,x22
      double precision s1,s2,t0


      x11 = na1-nb1
      x12 = na1-nb2
      x21 = na2-nb1
      x22 = na2-nb2

c      t0=2.*Heaviside(x22)*(sinh((q*l_eff)**2.*x22)-(q*l_eff)**2.*x22)
c     #-2.*Heaviside(x12)*(sinh((q*l_eff)**2.*x12)-(q*l_eff)**2.*x12)
c     #-2.*Heaviside(x21)*(sinh((q*l_eff)**2.*x21)-(q*l_eff)**2.*x21)
c     #+2.*Heaviside(x11)*(sinh((q*l_eff)**2.*x11)-(q*l_eff)**2.*x11)    
c     #-(exp((q*l_eff)**2.*x11)+exp((q*l_eff)**2.*x22)-
c     #exp((q*l_eff)**2.*x21)-exp((q*l_eff)**2.*x12))

      t0=(Heaviside(x22)-1.d0)*exp((q*l_eff)**2*x22)-Heaviside(x22)*
     #exp(-(q*l_eff)**2*x22)-2.*Heaviside(x22)*(q*l_eff)**2*x22+
     #(Heaviside(x11)-1.d0)*exp((q*l_eff)**2*x11)-Heaviside(x11)*
     #exp(-(q*l_eff)**2*x11)-2.*Heaviside(x11)*(q*l_eff)**2*x11-
     #(Heaviside(x12)-1.d0)*exp((q*l_eff)**2*x12)+Heaviside(x12)*
     #exp(-(q*l_eff)**2*x12)+2.*Heaviside(x12)*(q*l_eff)**2*x12-
     #(Heaviside(x21)-1.d0)*exp((q*l_eff)**2*x21)+Heaviside(x21)*
     #exp(-(q*l_eff)**2*x21)+2.*Heaviside(x21)*(q*l_eff)**2*x21

      t0=t0/(q*l_eff)**4

      intSab = t0/N
      
      return
      end


c ==============================================================================

      double precision function Heaviside(x)

      implicit none

      double precision x
      
      if (x.lt.0.) then
         Heaviside = 0.
      else
         Heaviside = 1.
      endif

      return
      end

c ==============================================================================
