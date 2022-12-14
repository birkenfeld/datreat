MODULE ADAPTIVE_INTEGRAL_MP
!
! Adative integration with internal STACK, faster than the recursive function from the "integration" module
!
! use "adapint_mp"   for the top level integration
!     "a2dapint_mp"  for integrations in the kernel functions used by adapint_mp
!     "a3dapint_mp   for integrations in the kernel functions used b2 a2dapint_mp
!
! with the addition mp_intparam_select variable the parameters in kernel functions
! can be selecetd according to the summation index in a parllelized do loop
!


!   abstract interface 
!     function f_kernel(x)
!     double precision  :: f_kernel
!     double precision  :: x
!     end function f_kernel
!   end interface 
!

   abstract interface 
     function f_kernel_mp(x,mp_intparam_select)
     double precision  :: f_kernel_mp
     double precision  :: x
     integer           :: mp_intparam_select
     end function f_kernel_mp
   end interface 



      integer, parameter, private           :: maxstack = 2000  ! --> stacktiefe

      integer, parameter, private           :: ndim = 8
      integer, Parameter, private           :: ndim2 = ndim / 2
      double precision, parameter, private  :: a (ndim2) = [0.362683783378362d0, 0.313706645877887d0,  &
                                                            0.222381034453374d0, 0.101228536290376d0 ] 
      double precision, parameter, private  :: x (ndim2) = [0.183434642495650d0, 0.525532409916329d0,  &
                                                            0.796666477413627d0, 0.960289856497536d0 ]  
      private :: gaus8a_mp, g2aus8a_mp, g3aus8a_mp

CONTAINS
                       

      FUNCTION adapint_mp (f, ip, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                         
                                                                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
      integer         , intent(in)   ::  ip
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  adapint_mp
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
        eps = abs(f((a+b)/2d0, ip)) * 1d-10
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
      s (itop) = gaus8a_mp (f,ip, a, b) 
                                                                        
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
         s (itop) = gaus8a_mp (f,ip, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = gaus8a_mp (f,ip, xa (itop), xb (itop) ) 
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
      WRITE (6, * ) 'adapint_mp fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      adapint_mp = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION adapint_mp        


      FUNCTION a2dapint_mp (f,ip, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                       
                                                                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
      integer         , intent(in)   ::  ip
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  a2dapint_mp
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
        eps = abs(f((a+b)/2d0, ip)) * 1d-10
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
      s (itop) = g2aus8a_mp (f,ip, a, b) 
                                                                        
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
         s (itop) = g2aus8a_mp (f,ip, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = g2aus8a_mp (f,ip, xa (itop), xb (itop) ) 
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
      WRITE (6, * ) 'a2dapint_mp fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      a2dapint_mp = erg 

                                                                        
      END FUNCTION a2dapint_mp        



      FUNCTION a3dapint_mp (f,ip, a, b, epso, maxitero, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
      IMPLICIT none 
!                                                                       
                       
!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
                                                                 
      integer         , intent(in)   ::  ip
      double precision, intent(in)   ::  a
      double precision, intent(in)   ::  b
      double precision, intent(in), optional   ::  epso
      integer         , intent(in), optional   ::  maxitero
      double precision, intent(out),optional   ::  erroraccu
      double precision               ::  a3dapint_mp
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
        eps = abs(f((a+b)/2d0, ip)) * 1d-10
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
      s (itop) = g2aus8a_mp (f, ip, a, b) 
                                                                        
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
         s (itop) = g2aus8a_mp (f,ip, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = g2aus8a_mp (f,ip, xa (itop), xb (itop) ) 
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
      WRITE (6, * ) 'a3dapint_mp fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
      if(present(erroraccu)) erroraccu = erraccu      
                                                                  
      a3dapint_mp = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION a3dapint_mp        


 
      FUNCTION gaus8a_mp (f, mp_intparam_select, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
      integer         , intent(in)   ::  mp_intparam_select 
      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  gaus8a_mp
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
        gsum = gsum + a (i) * (f (xave+grange * x (i)  , mp_intparam_select ) &
                            +  f (xave-grange * x ( i) , mp_intparam_select  ) ) 
      enddo 
      gaus8a_mp = gsum * grange 
                                                                        
      END FUNCTION gaus8a_mp  


      FUNCTION g2aus8a_mp (f, mp_intparam_select, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
      integer         , intent(in)   ::  mp_intparam_select 

      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  g2aus8a_mp
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
        gsum = gsum + a (i) * (f (xave+grange * x (i)  , mp_intparam_select ) &
                            +  f (xave-grange * x ( i) , mp_intparam_select  ) )   
      enddo 
      g2aus8a_mp = gsum * grange 
                                                                        
      END FUNCTION g2aus8a_mp  


      FUNCTION g3aus8a_mp (f , mp_intparam_select,xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      implicit none

!     double precision, EXTERNAL     ::  f 
      procedure(f_kernel_mp)         ::  f
      integer         , intent(in)   ::  mp_intparam_select 

      double precision, intent(in)   ::  xu
      double precision, intent(in)   ::  xo
      double precision               ::  g3aus8a_mp
 
      integer          :: i
      double precision :: xave, grange, gsum
                                                                  
      xave = (xo + xu) * 0.5d0 
      grange = (xo - xu) * 0.5d0 
      gsum = 0.d0 
      DO i = 1, ndim2 
          gsum = gsum + a (i) * (f (xave+grange * x (i)  , mp_intparam_select ) &
                              +  f (xave-grange * x ( i) , mp_intparam_select  ) )  
      enddo 
      g3aus8a_mp = gsum * grange 
                                                                        
      END FUNCTION g3aus8a_mp  



END MODULE ADAPTIVE_INTEGRAL_MP
!! 
!! 
!! 
!!  double precision function ff(x,i)
!!    double precision, intent(in) :: x
!!    ff = 1d0/sqrt(x) * i**2
!!  end function ff
!!  
!! program testada
!!   use ADAPTIVE_INTEGRAL_MP
!! 
!!   implicit none
!!   integer :: mp_intparam_select
!!   double precision :: ea = 0
!!   double precision, external :: ff
!!   
!!   mp_intparam_select = 1
!!   write(*,*) adapint_mp(ff,mp_intparam_select,0d0,1d0,1d-15,1000,ea)
!!   write(*,*) a2dapint_mp(ff,mp_intparam_select,0d0,1d0)
!!   write(*,*) a3dapint_mp(ff,mp_intparam_select,0d0,1d0)
!! 
!!   
!!   mp_intparam_select = 3
!!   write(*,*) adapint_mp(ff,mp_intparam_select,0d0,1d0,1d-15,1000,ea)
!!   write(*,*) a2dapint_mp(ff,mp_intparam_select,0d0,1d0)
!!   write(*,*) a3dapint_mp(ff,mp_intparam_select,0d0,1d0)
!! 
!!   write(*,*)ea
!! 
!! end program testada
