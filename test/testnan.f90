module nan

 interface is_inf
   module procedure is_inf_real
   module procedure is_inf_dbl
 end interface 

 interface is_nan
   module procedure is_nan_real
   module procedure is_nan_dbl
 end interface 

 interface stay_numeric
   module procedure stay_numeric_real
   module procedure stay_numeric_dbl
   module procedure stay_numeric_dvec
   module procedure stay_numeric_rvec
 end interface


contains

  logical function is_nan_real( x )
    implicit none
    real, intent(in) :: x
    is_nan_real = ( x .ne. x )
  end function is_nan_real


  logical function is_nan_dbl( x )
    implicit none
    real(kind=kind(1d0)), intent(in) :: x
    is_nan_dbl = ( x .ne. x )
  end function is_nan_dbl

  
  logical function is_inf_real( x )
    implicit none
    real, intent(in) :: x
    is_inf_real = ( x + Huge(x) == x )
  end function is_inf_real

  logical function is_inf_dbl( x )
    implicit none
    real(kind=kind(1d0)), intent(in) :: x
    is_inf_dbl = ( x + Huge(x) == x )
  end function is_inf_dbl

  function stay_numeric_real( x ) result(y)
    implicit none
    real, parameter  :: myNumInf   = 3.3333e33  ! or Huge(x)
    real, parameter  :: NaNreplace = 3.3333e-33 ! ???
    real, intent(in) :: x
    real             :: y
    
! deal with infinies
    if( x+Huge(x) == x ) then
      if( x*x == x ) then
        y =  myNumInf
      else
        y = -myNumInf
      endif
      return
    endif    

! deal with Nans
    if( x .ne. x) then
      y = NaNreplace
      return
    endif

    y = x

  end function stay_numeric_real


  function stay_numeric_dbl( x) result(y)
    implicit none
    real(kind=kind(1d0)), parameter  :: myNumInf   = 3.3333d33  ! or Huge(x)                                                                    
    real(kind=kind(1d0)), parameter  :: NaNreplace = 3.3333d-33 ! ???                                                                           
    real(kind=kind(1d0)), intent(in) :: x
    real(kind=kind(1d0))             :: y

! deal with infinies                                                                                                            
    if( x+Huge(x) == x ) then
      if( x*x == x ) then
        y =  myNumInf
      else
        y = -myNumInf
      endif
      return
    endif

! deal with Nans                                                                                                                
    if(x .ne. x) then
       y= NaNreplace
      return
    endif

    y = x

  end function stay_numeric_dbl



  function stay_numeric_dvec( x, n) result(y)
    implicit none
    real(kind=kind(1d0)), parameter  :: myNumInf   = 3.3333d33  ! or Huge(x)                                                                    
    real(kind=kind(1d0)), parameter  :: NaNreplace = 3.3333d-33 ! ???                                                                           
    real(kind=kind(1d0)), intent(in) :: x(n)
    real(kind=kind(1d0))             :: y(n)
    integer, intent(in)              :: n

    integer :: i

    do i=1,n
! deal with infinies                                                                                                            
    if( x(i)+Huge(x(i)) == x(i) ) then
      if( x(i)*x(i) == x(i) ) then
        y(i) =  myNumInf
      else
        y(i) = -myNumInf
      endif
    elseif(x(i) .ne. x(i)) then ! deal with Nans                                                                                                    
       y(i)= NaNreplace
    else
     y(i) = x(i)
    endif
   enddo

  end function stay_numeric_dvec


  function stay_numeric_rvec( x, n) result(y)
    implicit none
    real(kind=kind(1.0)), parameter  :: myNumInf   = 3.3333e33  ! or Huge(x)                                                                    
    real(kind=kind(1.0)), parameter  :: NaNreplace = 3.3333e-33 ! ???                                                                           
    real(kind=kind(1.0)), intent(in) :: x(n)
    real(kind=kind(1.0))             :: y(n)
    integer, intent(in)              :: n

    integer :: i

    do i=1,n
! deal with infinies                                                                                                            
    if( x(i)+Huge(x(i)) == x(i) ) then
      if( x(i)*x(i) == x(i) ) then
        y(i) =  myNumInf
      else
        y(i) = -myNumInf
      endif
    elseif(x(i) .ne. x(i)) then ! deal with Nans                                                                                                    
       y(i)= NaNreplace
    else
     y(i) = x(i)
    endif
   enddo

  end function stay_numeric_rvec




end module nan


program testnan
use nan
implicit none
real(kind=kind(1.d0)) :: x, y(-3:3), z
integer             :: i


do i=-3,3
  x = i
  write(6,'(i8,3g18.8)')i,x,Huge(x),-Huge(x)
  y(i) = 1d0/x
  z = sqrt(x)
  write(6,*)i, x, y(i), z
!  if( y+Huge(y) == y                    ) write(6,≈π*)"is infinite"
!  if( y+Huge(x) == y .and. y*y  == y    ) write(6,*)"is positive infinite"
!  if( y+Huge(x) == y .and. y*y .ne. y   ) write(6,*)"is negative infinite"
!  if( z .ne. z                          ) write(6,*)"z is NaN"

  if( is_inf(y(i))) write(6,*)"is infinite", stay_numeric(y(i))
  if( is_nan(z)) write(6,*)"z is NaN"   , stay_numeric(z)

enddo

   write(6,*) y
   write(6,*) stay_numeric(y,size(y))

end program testnan
