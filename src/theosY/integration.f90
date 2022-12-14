

!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------
MODULE integration 
!----------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------


!!!! adaptiv integration to fortran !!!! 
!!!! iterationcounter should not be global to the module
 

   implicit none

   abstract interface 
     function f_kernel(x)
     double precision  :: f_kernel
     double precision  :: x
     end function f_kernel
   end interface 


   abstract interface 
     double precision recursive function a_integrator(fx, a, b)
     double precision, external   :: fx
     double precision, intent(in) :: a, b
     end function a_integrator
   end interface 



 

   double precision, dimension(8), parameter :: xgausslaguerre8 =     &
       [ 0.170279632305d0,                                            &
         0.903701776799d0,                                            &
         2.251086629866d0,                                            &
         4.266700170288d0,                                            &
         7.045905402393d0,                                            &
        10.758516010181d0,                                            &
        15.740678641278d0,                                            & 
        22.863131736889d0 ]
 
  double precision, dimension(8), parameter ::  wgausslaguerre8 =       & 
       [ 3.69188589342d-1,                                              & 
         4.18786780814d-1,                                              & 
         1.75794986637d-1,                                              & 
         3.33434922612d-2,                                              & 
         2.79453623523d-3,                                              & 
         9.07650877336d-5,                                              & 
         8.48574671627d-7,                                              & 
         1.04800117487d-9 ]
 

 double precision, dimension(8), parameter ::  xgauss16        =        &
       [ 0.095012509837637440185d0,                                       &
         0.281603550779258913230d0,                                       &
         0.458016777657227386342d0,                                       &
         0.617876244402643748447d0,                                       &
         0.755404408355003033895d0,                                       &
         0.865631202387831743880d0,                                       &
         0.944575023073232576078d0,                                       &
         0.989400934991649932596d0 ]
 
 double precision, dimension(8), parameter :: wgauss16        =        &
       [ 0.189450610455068496285d0,                                      &
         0.182603415044923588867d0,                                      &
         0.169156519395002538189d0,                                      &
         0.149595988816576732081d0,                                      &
         0.124628971255533872052d0,                                      &
         0.095158511682492784810d0,                                      &
         0.062253523938647892863d0,                                      &
         0.027152459411754094852d0 ]
 

double precision, dimension(4), parameter ::       xgauss8         =       &
       [ 0.183434642495650d0      ,                                        &
         0.525532409916329d0      ,                                        &
         0.796666477413627d0      ,                                        &
         0.960289856497536d0       ]                    
 
double precision, dimension(4), parameter ::       wgauss8         =      &
       [ 0.362683783378362d0      ,                                       &
         0.313706645877887d0      ,                                       &
         0.222381034453374d0      ,                                       &
         0.101228536290376d0       ]      
 
double precision, dimension(2), parameter ::       xgauss4         =      &
       [ 0.339981043584856d0      ,                                       &
         0.861136311594053d0       ]
double precision, dimension(2), parameter ::       wgauss4         =      &
       [ 0.652145154862546d0      ,                                       &
         0.347854845137454d0       ]
 
 
 
double precision, private ::      epsilon_i = 1.0d-10
double precision, private ::      minabdist = 1.0d-20     ! should be more local
double precision, private ::      erroraccu = 0.0d0
integer         , private ::      maxit     = 150
integer         , private ::      printit   = 0
integer         , private ::      methode      
logical         , private ::      accwarning= .false.
integer         , private ::      nesting   = 0


procedure(a_integrator), pointer, private  :: intmethod 
procedure(a_integrator), pointer, private  :: intadause 
 
integer                         , private :: unitmsg = 6

                                  private createwarnings
                                  private integratecheck
                                  private adaptiveintegral

CONTAINS 
 

SUBROUTINE initint()
   implicit none
   printit = 0
   nesting = 0
   unitmsg = 6
END SUBROUTINE initint

SUBROUTINE setwriteunit(unit)
   implicit none
   integer, intent(in) :: unit
   unitmsg = unit
END SUBROUTINE setwriteunit
  
SUBROUTINE setmaxiteration( maxiteration )
   implicit none
   integer, intent(in) :: maxiteration
   maxit   = maxiteration
   if(printit > 3) write(unitmsg,*)"integral.maxiteration set to ...............:",maxit
END SUBROUTINE setmaxiteration
 
 
SUBROUTINE integral_setaccuracy( accuracy ) 
   implicit none
   double precision :: accuracy
   epsilon_i = accuracy
   if(printit > 3) write(unitmsg,*)"\nintegral.accuracy set to ...................: %e",epsilon_i
END SUBROUTINE integral_setaccuracy
 
 
 
 
double precision FUNCTION integral_getaccuracy() 
       integral_getaccuracy = erroraccu
END FUNCTION integral_getaccuracy
 
 
 
 
SUBROUTINE integral_setmininterval( mininterval ) 
   implicit none
   double precision, intent(in) :: mininterval
   minabdist = mininterval
   if(printit > 3) write(unitmsg,*)"integral.minimal integrationinterval set to : %e",minabdist
END SUBROUTINE integral_setmininterval 
 
 
 
SUBROUTINE intprint(ip) 
   implicit none
   integer, intent(in) :: ip
   printit = ip
   if(printit > 3) write(unitmsg,*)"integral.printing set to ...................: on"
END SUBROUTINE intprint
 
 
 
double precision recursive FUNCTION gausslaguerreintegral( f ) 
   implicit none
   double precision, external :: f
 
   double precision ::    s
   integer          ::    l
 
    s = 0.0

    do l=1,size(xgausslaguerre8)
      s = s + f(xgausslaguerre8(l)) * wgausslaguerre8(l)
    enddo

    gausslaguerreintegral = s
END function gausslaguerreintegral 
 
 
 
 
double precision recursive FUNCTION gint4( f , a, b )
    implicit none
    double precision, external    :: f
   
    double precision, intent(in)  :: a, b
 
    integer                       ::     i 
    double  precision             :: hw, ce , s, x1, x2, w    ;
 
     hw = (b-a)*0.5d0
     ce = (b+a)*0.5d0
     s  = 0.0d0
     do i=1, SIZE(xgauss4)
       w  = wgauss4(i)
       x1 = hw * xgauss8(i)
       x2 = ce + x1
       x1 = ce - x1
       s  = s + w*( f(x1) + f(x2) )
     enddo
     gint4 = hw*s
END FUNCTION gint4 
 
 
 
double precision recursive FUNCTION gint8( f , a, b )
    implicit none
    double precision, external    :: f
    double precision, intent(in)  :: a, b
 
    integer                       ::     i 
    double  precision             :: hw, ce , s, x1, x2, w    ;
 
     hw = (b-a)*0.5d0
     ce = (b+a)*0.5d0
     s  = 0.0d0
     do i=1, SIZE(xgauss8)
       w  = wgauss8(i)
       x1 = hw * xgauss8(i)
       x2 = ce + x1
       x1 = ce - x1
       s  = s + w*( f(x1) + f(x2) )
     enddo
     gint8 = hw*s
END FUNCTION gint8 
 
 
 
double precision recursive FUNCTION gint16( f , a, b )
    implicit none
    double precision, external    :: f
    double precision, intent(in)  :: a, b
 
    integer                       ::     i 
    double  precision             :: hw, ce , s, x1, x2, w    ;
 
     hw = (b-a)*0.5d0
     ce = (b+a)*0.5d0
     s  = 0.0d0
     do i=1, SIZE(xgauss16)
       w  = wgauss16(i)
       x1 = hw * xgauss16(i)
       x2 = ce + x1
       x1 = ce - x1
       s  = s + w*( f(x1) + f(x2) )
     enddo
     gint16 = hw*s
END FUNCTION gint16 
 
 
SUBROUTINE createwarnings(iterationcounter)
   implicit none
   integer, intent(in) :: iterationcounter

   if (accwarning ) then
       write(unitmsg,*)"adaint: low accuracy erroraccu=",erroraccu
   endif
   if (iterationcounter > maxit) then
      write(unitmsg,*)"adaint: maxit < iterations=%i",iterationcounter
   endif
END SUBROUTINE createwarnings
 
 
 
double precision recursive FUNCTION integratecheck( intf,f,ax,bx,estimate, iterationcounter ) RESULT (int_chk)
  implicit none
  procedure(a_integrator),pointer  :: intf
  procedure(f_kernel)              :: f
  double precision, intent(in)     :: ax, bx
  double precision, intent(in)     :: estimate
  integer         , intent(inout)  :: iterationcounter

  double precision :: irechts, ilinks, ab, sum, error 
  logical          ::     toosmallinterval                
 
   ab               = 0.5d0 * (ax+bx)
   irechts          = intf(f,ax,ab)
   ilinks           = intf(f,ab,bx)
   sum              = irechts + ilinks
   error            = abs( sum - estimate )
   toosmallinterval = ( abs(bx-ax) < minabdist )
   accwarning       = accwarning .or. toosmallinterval
 
   if (printit > 0 ) then
       write(unitmsg,'(4e14.7)') ax,bx,sum,estimate
       if (toosmallinterval) write(unitmsg,'(a)')" ! "
   endif
 
   if (( error < epsilon_i ) .or. toosmallinterval .or. (iterationcounter > maxit)) then
      erroraccu = erroraccu + error
!      integratecheck =  (sum); return
      int_chk =  (sum); return
   else  
      iterationcounter =  iterationcounter + 1
 
      irechts = integratecheck(intf,f,ax,ab,irechts, iterationcounter)
      ilinks  = integratecheck(intf,f,ab,bx,ilinks, iterationcounter )
 
   endif
 
!   integratecheck = (irechts+ilinks)
   int_chk = (irechts+ilinks)
 
END FUNCTION  integratecheck 
 
 
 
double precision recursive FUNCTION adaptiveintegral( f, a0,  b0 )
   implicit none
   double precision, external    :: f
   double precision, intent(in)  :: a0, b0
   double precision              :: wert
   integer                       :: iterationcounter
 
!   /* init */
   accwarning       = .false.
   iterationcounter = 0
   erroraccu        = 0.0d0
 
   wert     = intadause( f, a0,b0 )
   wert     = integratecheck( intadause, f, a0,b0, wert, iterationcounter    )
 
   call createwarnings(iterationcounter)
 
   adaptiveintegral =  (wert)    
 
END FUNCTION adaptiveintegral
 
 
 
 
SUBROUTINE integral_setmethod(  method ) 
   implicit none
   integer, intent(in)  :: method
   character(len=255)   :: buffer
   methode = method;
   select case(methode) 
   case ( 1  )  
       intmethod => gint4 
       intadause => gint4 
   case ( 2  )   
       intmethod => gint8 
       intadause => gint8  
   case ( 3  )   
       intmethod => gint16 
       intadause => gint16
   case ( 4  )   
       intadause => gint4
       intmethod => adaptiveintegral
   case ( 5  )   
       intadause => gint8
       intmethod => adaptiveintegral
   case ( 6  )   
       intadause => gint16
       intmethod => adaptiveintegral
   case ( 7  )   
       intadause => gint4
       intmethod => adaptiveintegral
   case ( 8  )   
       intadause => gint8
       intmethod => adaptiveintegral
   case ( 9  )   
       intadause => gint16
       intmethod => adaptiveintegral
   case default
               write(unitmsg,'(a,i4)')"intsetmethod invalid method: ",methode
               intadause => gint4
               intmethod => adaptiveintegral
               methode   = 4

  end select
!  write(unitmsg,'(a,i4)')"integral.method   set to ...................: %d",methode
END SUBROUTINE integral_setmethod
 
 
 
double precision recursive FUNCTION integral( f, a0, b0, maxiter, aimed_error, erraccu ) 
  implicit none
  double precision, external   :: f
  double precision, intent(in) :: a0, b0 
  integer,          intent(in), optional   :: maxiter
  double precision, intent(in), optional   :: aimed_error
  double precision, intent(inout),optional :: erraccu        ! on call set to zero
  double precision             :: wert

  if(present(maxiter)) then
    maxit = maxiter
  endif

  if(present(aimed_error)) then
    call integral_setaccuracy(aimed_error)
    call integral_setmethod(5)
  endif

  wert     = intmethod( f, a0,b0 )
  integral = wert

  if(present(erraccu)) then
    erraccu = erraccu+integral_getaccuracy()
  endif

END FUNCTION integral
 
END MODULE integration


