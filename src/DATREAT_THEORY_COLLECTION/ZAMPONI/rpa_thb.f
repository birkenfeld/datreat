      function th16(x,pa,thnam,parnam,npar,ini)
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
      REAL            th16, x, pa, qq, zpi, xh, vol_frac         
      INTEGER         ini, npar, nparx
      DIMENSION       pa(20),qq(3)
      DATA            zpi/6.283185/
      double precision n1a,n2a,n3a,n1b,n2b,n3b,n4b,n5b
      double precision q,rpa,leff,incoh1,a1,a2,incoh2,ratio
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
      double precision B  (maxd)

      integer nblka,nblkb
      double precision ng(2,maxd),n
      integer ndim
    
      double precision bSab,bintSab
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
         thnam = 'rpab '

C         nparx = 14
         nparx = 12

         if(npar.lt.nparx) then
            write(6,1)thnam,nparx,npar
 1          format(' theory: ',a8,' no of parametrs=',i8,
     *           ' exceeds current max. = ',i8)
            th16   = 0
            return
         endif
         npar = nparx


c     --------------> set the number of parameters
         parnam(1) = 'ampl1   '          
         parnam(2) = 'ampl2   '          
         parnam(3) = 'l       '
         parnam(4) = 'inc_bk1 '
         parnam(5) = 'inc_bk2 '
         parnam(6) = 'ratio   '
         parnam(7) = 'n1a     '
         parnam(8) = 'n2a     '
         parnam(9) = 'n3a     '
         parnam(10) = 'n1b     '
         parnam(11) = 'n2b     '
         parnam(12) = 'n3b     '
C         parnam(13) = 'n4b     '
C         parnam(14) = 'n5b     '
c     
          th16  = 0
         return
      endif
c     
c     ---- calculate theory here -----
      a1      = abs(pa(1))      ! amplitude
      a2      = abs(pa(2))      ! amplitude
      leff    = abs(pa(3))      ! segment length
      incoh1  = abs(pa(4))      ! incoh. bkgr 
      incoh2  = abs(pa(5))      ! incoh. bkgr 
      ratio    = abs(pa(6))           ! fraction of homopolymer (matrix)
      n1a      = abs(pa(7))      ! first block 
      n2a      = abs(pa(8))      ! second block
      n3a      = abs(pa(9))      ! third block
      n1b      = abs(pa(10))      ! fourth block
      n2b      = abs(pa(11))      ! fifth block
      n3b      = abs(pa(12))      ! fifth block
c      n4b      = abs(pa(13))      ! fifth block
c      n5b      = abs(pa(14))      ! fifth block
      
      call parget('frac    ',frac,iadda,ier)
      phi0 = frac 

      q = x 
      
      
! 0. Limit Parameter
      v0 = 1d5                          ! v0 --> infinity to realize excluded volume

      nblka = 3
c      nblkb = 8
      nblkb = 6
      ndim   = nblkb+1 

      ng(1,1) = 1                       ! Block limits in terms of segment numbers
      ng(2,1) = n1a
      ng(1,2) = n1a+1
      ng(2,2) = n1a+n2a
      ng(1,3) = n1a+n2a+1
      ng(2,3) = n1a+n2a+n3a
      ng(1,4) = 1
      ng(2,4) = n1b
      ng(1,5) = n1b+1
      ng(2,5) = n1b+n2b
      ng(1,6) = n1b+n2b+1
      ng(2,6) = n1b+n2b+n3b
c      ng(1,7) = n1b+n2b+n3b+1
c      ng(2,7) = n1b+n2b+n3b+n4b
c      ng(1,8) = n1b+n2b+n3b+n4b+1
c      ng(2,8) = n1b+n2b+n3b+n4b+n5b
      n       = ng(2,3)                 ! unconnected block of total length (matrix polymer)

      l_eff   = leff/sqrt(6.0d0)

      do i=1,ndim
         do j=1,ndim
            if(i.le.nblka .and. j.le.nblka) then
               S(i,j)=bintSab(q,l_eff,n,ng(1,i),ng(2,i),ng(1,j),ng(2,j))
     *              *(1-Phi0-ratio)
            elseif(i.gt.nblka.and.i.le.nblkb.and.j.gt.nblka.and.
     *              j.le.nblkb)then
               S(i,j)=bintSab(q,l_eff,n,ng(1,i),ng(2,i),ng(1,j),ng(2,j))
     *              *ratio
            else
               S(i,j) = 0
            endif
            V(i,j) = v0
         enddo
      enddo
      S(ndim,ndim) = bintSab(q,l_eff,n,1.d0,n,1.d0,n)*Phi0 

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
      B(5)  = 1 
      B(6)  = 0 
      B(7)  = 0 
c      B(8)  = 0 
c      B(9)  = 0 
      sum = 0
      do i=1,ndim
         do j=1,ndim
            sum = sum + B(i)*S0(i,j)*B(j)
         enddo
      enddo        

      th16 = sum
      
      if(frac.eq.0.6) then
         th16 = th16 * a1 + incoh1
      else
         th16 = th16 * a2 + incoh2
      endif
           

 999  continue

      return
      end
      

c ==============================================================================

      double precision function bSab(q,l_eff,N,na1,na2,nb1,nb2)
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
       
       bSab = sum/N

       return
       end


c ==============================================================================

      double precision function bintSab(q,l_eff,N,na1,na2,nb1,nb2)

      implicit none 

      double precision q,l_eff
      double precision N,na1,na2,nb1,nb2,Bhs
      double precision x11,x12,x21,x22
      double precision s1,s2,t0


      x11 = na1-nb1
      x12 = na1-nb2
      x21 = na2-nb1
      x22 = na2-nb2

      t0=(Bhs(x22)-1.d0)*exp((q*l_eff)**2*x22)-Bhs(x22)*
     #exp(-(q*l_eff)**2*x22)-2.*Bhs(x22)*(q*l_eff)**2*x22+
     #(Bhs(x11)-1.d0)*exp((q*l_eff)**2*x11)-Bhs(x11)*
     #exp(-(q*l_eff)**2*x11)-2.*Bhs(x11)*(q*l_eff)**2*x11-
     #(Bhs(x12)-1.d0)*exp((q*l_eff)**2*x12)+Bhs(x12)*
     #exp(-(q*l_eff)**2*x12)+2.*Bhs(x12)*(q*l_eff)**2*x12-
     #(Bhs(x21)-1.d0)*exp((q*l_eff)**2*x21)+Bhs(x21)*
     #exp(-(q*l_eff)**2*x21)+2.*Bhs(x21)*(q*l_eff)**2*x21

      t0=t0/(q*l_eff)**4

      bintSab = t0/N
      
      return
      end


c ==============================================================================

      double precision function Bhs(x)

      implicit none

      double precision x
      
      if (x.lt.0.) then
         Bhs = 0.
      else
         Bhs = 1.
      endif

      return
      end

c ==============================================================================
