module lcalc_Zomega_cc
!
! gfortran adaptive_integral_module.f90 lcalc_20_Zomega_cc.f90 lib/slatec/libslatec.a -O2
!
! Compute inductance matrix for an assembly of concentric cylindrical coils ic=1..ncoils
! each with current distributed between radii:  ri(ic) and ra(ic)
! over axis extension                        :  z(1,ic) .. z(2,ic)
! with number of windings N                  :  nturns(ic)
!    if N >  1 it is assumed that the current density is homogeneous within the winding zone
!    if N == 1 it is assumed that the current density is proportional to 1/r
!                 applies to winding forms etc.
!
!
!
!
! use adaptive_integral
 implicit none

 private
 public :: Z_llcoil,  Zomega


   abstract interface 
     function f_kernel(x)
     double precision  :: f_kernel
     double precision  :: x
     end function f_kernel
   end interface 


!   abstract interface 
!     function f_kernel_mp(x,i)
!     double precision  :: f_kernel_mp
!     double precision  :: x
!     integer           :: i
!     end function f_kernel_mp
!   end interface 




      integer, parameter, private           :: maxstack = 100000  ! --> stacktiefe

      integer, parameter, private           :: ndim = 8
      integer, Parameter, private           :: ndim2 = ndim / 2
      double precision, parameter, private  :: a (ndim2) = [0.362683783378362d0, 0.313706645877887d0,  &
                                                            0.222381034453374d0, 0.101228536290376d0 ] 
      double precision, parameter, private  :: x (ndim2) = [0.183434642495650d0, 0.525532409916329d0,  &
                                                            0.796666477413627d0, 0.960289856497536d0 ]  


                       
!! integrand
!! y=cos(phi)/sqrt((r1-r2*cos(phi))**2+(r1-r2*sin(phi))**2+(z1-z2)**2) * r1 * r2
!!
  double precision, parameter :: Pi =4d0*atan(1d0)
  double precision, parameter :: mu0=1.25663706212d-6

  double precision, parameter :: rho_Cu      = 1.7d-8         !! specific resistance copper 20C
  double precision, parameter :: rho_Cu_tc   = 0.004d0        !! temp coeff specific resistance copper
  double precision, parameter :: rho_AlMg    = 5.0d-8         !! specific resistance AlMg alloy
  double precision, parameter :: rho_AlMg_tc = 0.004d0        !! temp coeff 
  double precision, parameter :: rho_AlSC     = 5.0d-14      !! rho for Al winding at 6K, UNCERTAIN!
                                                              !! the residual resitance varies steeply
                                                              !! with purity over many orders of magnitude
                                                              !! check with manufacturer
  double precision, parameter :: density_Cu  = 8960.d0        !! copper density

  double precision            :: rhoAlcore


!===============================================================
! COMMON shared with integration kernels to perform the 3D integration
! sequence, by these variables the parameters are made known to the kernels
  double precision :: r1,r2
  double precision :: rj1,rj2
  double precision :: R1i, R1a, Z1,  Z2
  double precision :: R2i, R2a, ZZ1, ZZ2
!===============================================================

  double precision :: Rinnen, Raussen, Laenge

  double precision :: ZArea1, ZArea2
  double precision :: erg   

  integer          :: maxit = 50000
  double precision :: eps0  = 1d-6

  integer, parameter :: mcoil                = 100
  integer            :: ncoils               = 1
  logical            :: oneturn(mcoil)       = .false.
  integer            :: nturns(mcoil)        = 1
  integer            :: nparts(mcoil)        = 1
  
  double precision   :: ri(mcoil)
  double precision   :: ra(mcoil)
  double precision   :: z(2,mcoil)

  double precision   :: ll(mcoil,mcoil) = 0  !! the resulting induction matrix

  double precision   :: nt, dlength, wwidth, rsep, zshift
 
  integer            :: ic1, ic2
  integer            :: i, j, ioc, ios, irho

  integer, parameter :: mllcoil = 10
  double precision   :: omega, freq
  integer            :: iom
  complex, parameter :: ci = (0d0,1d0)
  complex            :: cll(mllcoil,mllcoil)        
  complex            :: ca(mllcoil,mllcoil)        
  complex            :: clw(mllcoil,mllcoil)        
  complex            :: cwl(mllcoil,mllcoil)        
  complex            :: cww(mllcoil,mllcoil)        
  complex            :: cwwinv(mllcoil,mllcoil) 
  complex            :: ciw(mllcoil) , cul(mllcoil), cuw(mllcoil) !! currents in windingcores, voltages on partial coils      
  complex            :: zwk(mcoil)        
  complex            :: cs    
  complex            :: cdet(2)    
  complex            :: Zom   
  double precision   :: Rw(mllcoil)  = 0.d0
  double precision   :: CC(mllcoil)  = 0.d0
  ! double precision   :: CperSQM      = 50d-9   !! ASSUMING 50nF per m**2
  double precision   :: CperSQM      = 1d-9   !! ASSUMING 50nF per m**2


  integer            :: ipvt(mcoil)
  real               :: rcond

  double precision   :: cparallel = 5d-9  ! overall parallel capacitance 
 

CONTAINS

complex function Z_llcoil(f, rhoAl, cpara, cpersq, dlen, wwid, rsep, Rmodif) result(z)
  double precision, intent(in) :: f
  double precision, intent(in) :: rhoAl
  double precision, intent(in) :: cpara
  double precision, intent(in) :: cpersq
  double precision, intent(in) :: dlen
  double precision, intent(in) :: wwid
  double precision, intent(in) :: rsep  
  double precision, intent(in) :: Rmodif(mllcoil)

  double precision, save ::  rhoAl0, cpara0, cpersq0, dlen0, wwid0, rsep0

  if((dlen0   .ne. dlen)  .or.  &
     (wwid0   .ne. wwid)  .or.  &
     (rsep0   .ne. rsep)      ) then
     call init_ll(dlen, wwid, rsep)
  endif 

  z = Zomega(f, rhoAl, cpara, cpersq, Rmodif) 

  rhoAl0  = rhoAl
  cpara0  = cpara
  cpersq0 = cpersq
  dlen0   = dlen
  wwid0   = wwid
  rsep0   = rsep  

end function Z_llcoil

 
subroutine init_ll(dlength, wwidth, rsep)
 double precision, intent(in) :: dlength
 double precision, intent(in) :: wwidth
 double precision, intent(in) :: rsep
  
  open(10, file="coiltables.plain")
!  open(10, file="/Users/Apple/Desktop/HOMEWORK/NIST/LCALC/coiltables.plain")



  ncoils = 10         !! these are (mist be) the 10 partial coils of solmain1
  do i=1,ncoils
    read(10,*,iostat=ios)ic1, z(1,i), z(2,i), ri(i), ra(i), nt
    nturns(i) = nint(nt)
! write(*,*)i,ic1, z(1,i), z(2,i), ri(i), ra(i), nturns(i)
  enddo
  close(10)
! zshift   = (minval(z(1,1:ncoils))+maxval(z(2,1:ncoils)))/2
  zshift   = z(1,1)

  do i=1,ncoils
    z(1,i) = z(1,i) -  zshift 
    z(2,i) = z(2,i) -  zshift 

    z(1,i+ncoils) =  z(1,i) - dlength       !! here add (next 10 coils) the winding cores
    z(2,i+ncoils) =  z(2,i) + dlength       !! using (estimated) thickness of core cyl.
    ra(i+ncoils)  =  ri(i)-rsep             !! end extra width
    ri(i+ncoils)  =  ri(i)-rsep-wwidth
    nturns(i+ncoils) = 1 
  enddo
  
   ncoils = 2*ncoils
  

  oneturn(:)    = (nturns(:) == 1)


 write(*,*)
 write(*,*)"========================== coil model ======================================== "
 write(*,*)"coil:  z1         z2          ri        ra            N                        "
 do ic1=1,ncoils
   write(*,'(i4,":",2f10.6,3x,2f10.6,i8)') ic1, z(:,ic1), ri(ic1), ra(ic1), nturns(ic1)
 enddo
 write(*,*)"============================================================================== "
 write(*,*)


  ll = 0

 dc1: do ic1=1,ncoils
  dc2: do ic2=1,ncoils 
         R1i = ri(ic1)
         R1a = ra(ic1)
         Z1  = z(1,ic1)
         Z2  = z(2,ic1)
         if(oneturn(ic1)) then
           ZArea1 =  1d0/(log(R1a/R1i)*(Z2-Z1))
         else
           ZArea1 =  1d0/((R1a-R1i)*(Z2-Z1))
         endif

         R2i = ri(ic2)
         R2a = ra(ic2)
         ZZ1  = z(1,ic2)
         ZZ2  = z(2,ic2)
         if(oneturn(ic2)) then
            ZArea2 =  1d0/(log(R2a/R2i)*(ZZ2-ZZ1))
         else
            ZArea2 =  1d0/((R2a-R2i)*(ZZ2-ZZ1))
         endif

write(*,*)"ic1,ic2=",ic1,ic2

         ll(ic1,ic2) =  a3dapint(fker3,R2i,R2a,eps0*10,maxit)  &
                     * mu0/4*Pi *(2*Pi) * nturns(ic1) * nturns(ic2) &
                     * ZArea1 * ZArea2 &
                      /Pi/Pi
 ! write(*,*)"Erg = ",erg, erg*Nturns**2  /Pi/Pi !!???? PI**2 factor ??  
 !                                       ======> dieser Faktor scheint nötig, aber wo
 !                                               in der Herleitung sind die Pi's her???
 !                                               Check!!  L=   erg*Nturns**2  /Pi/Pi  STIMMT!!
 !                                                                            ------> wieso der Faktor
 ! 
 !        write(*,'("L(",i0,",",i0,")=",f16.12)')ic1,ic2, ll(ic1,ic2)
write(*,*)ll(ic1,ic2)
   enddo dc2
 enddo dc1

 write(*,*)
 write(*,*)"========================== Inductance Matrix [Hy] ============================ "
 do ic1=1,ncoils
   do ic2=1,ic1
      write(*,'(2i8,2x,e15.8)')ic1,ic2,ll(ic1,ic2) 
   enddo
 enddo
 write(*,*)"============================================================================== "
 write(*,*)


 write(*,*)
 write(*,*)"========================== Inductance Matrix [Hy] ============================ "
 do ic1=1,ncoils
   write(*,'(20f10.6)')ll(ic1,1:ncoils) 
 enddo
 write(*,*)"============================================================================== "
 write(*,*)

 write(*,*)
 write(*,'("sun(Lij)[Hy]  =",f16.12)') sum(ll) 


!!  effective mutual inductance: sqrt((L12**2/R2 + L13**2/R3 + L14**2/R4)*RR) !!
!!  with RR=1/(1/R2+1/R3+1/R4)
!!  sqrt((L12**2*R3*R4 + L13**2*R2*R4 + L14**2*R2*R3)/(R2*R3 + R2*R4 + R3*R4))
!!  Resistance of single cyl: 2*Pi*rho/(ln(Ra/Ri) * L)

!  double precision   :: omega, freq
!  complex, parameter :: ci = (0d0,1d0)
!  complex            :: cll(mllcoil,mllcoil)        
!  complex            :: clw(mllcoil,mllcoil)        
!  complex            :: cwl(mllcoil,mllcoil)        
!  complex            :: cww(mllcoil,mllcoil)        
!  complex            :: cwwinv(mllcoil,mllcoil)        
!  complex            :: cs    
!  double precision   :: Rw(mllcoil)  


!!! complex impedance of the 10 coil system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(ncoils/2 .ne. mllcoil) STOP "ERROR: numbers of subcoils dont match!"

end subroutine init_ll


!! HIER NOCH RW BERECHNEN!!
  
function Zomega(freq, rhoAlcore, cparallel, CperSQM, Rmodifier) result(zom)
  double precision, intent(in) :: freq
  double precision, intent(in) :: rhoAlcore
  double precision, intent(in) :: cparallel
  double precision, intent(in) :: cperSQM
  double precision, intent(in) :: Rmodifier(mllcoil)
  complex                      :: zom
 


  do i = 1,mllcoil
    ic1 = mllcoil+i
    Rw(i) = 2*Pi*rhoAlcore/(log(ra(ic1)/ri(ic1)) * (z(2,ic1)-z(1,ic1)))  * Rmodifier(i)
    CC(i) = 2*Pi*ri(i)*(z(2,i)-z(1,i)) * CperSQM
  enddo



  omega  = 2*Pi*freq

  do i=1,mllcoil
   do j=1,mllcoil
    cll(i,j) = ll(i,j)                  * ci * omega                  
    clw(i,j) = ll(i,j+mllcoil)          * ci * omega
    cwl(i,j) = ll(i+mllcoil,j)          * ci * omega
    cww(i,j) = ll(i+mllcoil,j+mllcoil)  * ci * omega
   enddo
   cww(i,i)  = cww(i,i) + Rw(i)                    
  enddo

!! include capacitance to the coil-matrix Lc = 1/(1/L+i*omega*C)
  ca = cll
  CALL CGECO (ca, mllcoil, mllcoil, IPVT, RCOND, zwk)
  CALL CGEDI (ca, mllcoil, mllcoil, IPVT, cDET, zwk, 11)
  do i=1,mllcoil
    ca(i,i) = ca(i,i) + ci * omega * CC(i)
  enddo
  CALL CGECO (ca, mllcoil, mllcoil, IPVT, RCOND, zwk)
  CALL CGEDI (ca, mllcoil, mllcoil, IPVT, cDET, zwk, 11)
  cll = ca


  ca = cww
  CALL CGECO (ca, mllcoil, mllcoil, IPVT, RCOND, zwk)
  CALL CGEDI (ca, mllcoil, mllcoil, IPVT, cDET, zwk, 11)
  cwwinv = ca

!! voltages on windingcores
    do i=1,mllcoil
     cs = (0.0,0.0)
     do j=1,mllcoil
       cs = cs + cwl(i,j)  
     enddo
     cuw(i) = cs
    enddo
!! --> currents through winding cores
    do i=1,mllcoil
     cs = (0.0,0.0)
     do j=1,mllcoil
       cs = cs - cwwinv(i,j) * cuw(j) !!  -  !! 
     enddo
     ciw(i) = cs
    enddo
!! --> Voltage vector on partical coils
    do i=1,mllcoil
     cs = (0.0,0.0)
     do j=1,mllcoil
       cs = cs + cll(i,j) + clw(i,j) * ciw(j)
     enddo
     cul(i) = cs
    enddo
!! --> total clamp voltage yields the impedance (unit input current assumed)
    Zom = sum(cul(1:mllcoil))

!! --> and for the effect of a parallel capacitance
    Zom = 1d0/(1d0/Zom + ci*omega*cparallel)
    
end function zomega



function fker1(phi) result (fval) ! ??
  double precision              :: phi
  double precision              :: fval

  double precision :: t1 
  double precision :: t2 
  double precision :: t3 
  double precision :: t7 
  double precision :: t8 
  double precision :: t9 
  double precision :: t10
  double precision :: t12
  double precision :: t14
  double precision :: t16
  double precision :: t20
  double precision :: t22
  double precision :: t25
  double precision :: t28
  double precision :: t32
  double precision :: t34
  double precision :: t36
  double precision :: t41
  double precision :: t43
  double precision :: t46
  double precision :: t49
  double precision :: t51
  double precision :: t53

  t1 =cos(phi)
  t2 =t1*r1
  t3 =Z1**2
  t7 =2*r2*t2
  t8 =r1**2
  t9 =r2**2
  t10=ZZ1**2
  t12=sqrt(-2*ZZ1*Z1+t10+t3-t7+t8+t9)
  t14=log(-Z1+ZZ1+t12)
  t16=Z2**2
  t20=sqrt(-2*ZZ1*Z2+t10+t16-t7+t8+t9)
  t22=log(-Z2+ZZ1+t20)
  t25=log(Z1-ZZ1+t12)
  t28=log(Z2-ZZ1+t20)
  t32=ZZ2**2
  t34=sqrt(-2*ZZ2*Z1+t3+t32-t7+t8+t9)
  t36=log(-Z1+ZZ2+t34)
  t41=sqrt(-2*ZZ2*Z2+t16+t32-t7+t8+t9)
  t43=log(-Z2+ZZ2+t41)
  t46=log(Z1-ZZ2+t34)
  t49=log(Z2-ZZ2+t41)
  t51=Z1*t14-Z1*t36-Z2*t22+Z2*t43+ZZ1*t25-ZZ1*t28-ZZ2*t46+ZZ2*t49+t12-t20-t34+t41
  t53=t51*r2*t2
  fval = t53 * rj1 * rj2

end function fker1




function fker2(r) result (fval)
  double precision              :: r
  double precision              :: fval

  r1 = r
  if(oneturn(ic1)) then
    rj1   = 1d0/r 
  else
    rj1   = 1
  endif
  fval = adapint(fker1,0.d0,2*Pi, eps0, maxit)
 
end function fker2


function fker3(r) result (fval)
  double precision              :: r
  double precision              :: fval

  r2  =  r
  if(oneturn(ic2)) then
    rj2   = 1d0/r 
  else
    rj2   = 1
  endif

  fval = a2dapint(fker2,R1i,R1a, eps0, maxit)
 
end function fker3





      FUNCTION adapint (f, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                         
                                                                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  adapint
!               
      integer           :: iterationcounter , maxiter
      integer           :: i, itop, itm  
      double precision  :: erg, xbb, error, erraccu, eps
      double precision   ::  s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
                                                         
!                                                                        
!                                                                       
!      -----------------------------------------------------------------
      if(present(epso)) then
        eps = epso
      else
        eps = abs(f((a+b)/2d0)) * 1d-10
      endif  
                 
      if(present(maxitero)) then
        maxiter = maxitero
      else
        maxiter = 1000
      endif                   
                                                     
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      erraccu = 0
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = gaus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! adaptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
                                                                        
         IF (error.lt.eps) then 
            erg = erg + s (itop) + s (itop - 1) 
            erraccu = erraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'adapint fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      adapint = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION adapint        


      FUNCTION a2dapint (f, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                       
                                                                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  a2dapint
!               
      integer           :: iterationcounter , maxiter
      integer           :: i, itop, itm  
      double precision  :: erg, xbb, error, erraccu, eps
      double precision   ::  s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
                                                         
!                                                                       
!      -----------------------------------------------------------------
      if(present(epso)) then
        eps = epso
      else
        eps = abs(f((a+b)/2d0)) * 1d-10
      endif  
                 
      if(present(maxitero)) then
        maxiter = maxitero
      else
        maxiter = 1000
      endif                   
                                                     
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      erraccu = 0
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = g2aus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! a2daptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
                                                                        
         IF (error.lt.eps) then 
            erg = erg + s (itop) + s (itop - 1) 
            erraccu = erraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'a2dapint fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      a2dapint = erg 

                                                                        
      END FUNCTION a2dapint        



      FUNCTION a3dapint (f, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f
                                                                 
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  a3dapint
!               
      integer           :: iterationcounter , maxiter
      integer           :: i, itop, itm  
      double precision  :: erg, xbb, error, erraccu, eps
      double precision   ::  s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
                                                         
!                                                                       
!      -----------------------------------------------------------------
      if(present(epso)) then
        eps = epso
      else
        eps = abs(f((a+b)/2d0)) * 1d-10
      endif  
                 
      if(present(maxitero)) then
        maxiter = maxitero
      else
        maxiter = 1000
      endif                   
                                                     
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      erraccu = 0
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = g2aus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! a3daptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
                                                                        
         IF (error.lt.eps) then 
            erg = erg + s (itop) + s (itop - 1) 
            erraccu = erraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'a3dapint fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      a3dapint = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION a3dapint        


 
      FUNCTION gaus8a (f, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f

      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  gaus8a
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
      gsum = gsum + a (i) * (f (xave+grange * x (i) ) + f (xave-grange * x ( i) ) )                                                            
      enddo 
      gaus8a = gsum * grange 
                                                                        
      END FUNCTION gaus8a  


      FUNCTION g2aus8a (f, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f

      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  g2aus8a
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
      gsum = gsum + a (i) * (f (xave+grange * x (i) ) + f (xave-grange * x ( i) ) )                                                            
      enddo 
      g2aus8a = gsum * grange 
                                                                        
      END FUNCTION g2aus8a  


      FUNCTION g3aus8a (f, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel)            ::  f

      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  g3aus8a
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
      gsum = gsum + a (i) * (f (xave+grange * x (i) ) + f (xave-grange * x ( i) ) )                                                            
      enddo 
      g3aus8a = gsum * grange 
                                                                        
      END FUNCTION g3aus8a  




end module lcalc_Zomega_cc


!
!program testll
! use lcalc_20_Zomega_cc
! implicit none
!
! double precision :: frequenz
! double precision :: R0, zabs, zphase
! complex          :: zx
! integer          :: iomeg
!
!
!  dlength      = 0.01d0
!  wwidth       = 0.01d0  
!  rsep         = 0.0005d0
!  CperSQM      = 1d-9   
!  cparallel    = 5d-9 
!  rhoAlcore    = 2d-9
!  R0           = 0.0d0
!
!  
!diom: do iomeg=0,2500
!  
!  frequenz   = 0.0001d0 * 2d0**(iomeg*0.01d0)
!
!  zx =  Z_llcoil(frequenz, rhoAlcore, cparallel, cpersqm, dlength, wwidth, rsep) + R0
!
!  write(*,'("f=",f12.6,"   |Z|=",f18.4,"  Zr=",f18.4,"  Zi=",f18.6,"  phase=",F12.4)') &
!               frequenz,      abs(zx),    zx, atan2(aimag(zx),real(zx))*180d0/Pi
!
!end do diom
! 
! end program testll
