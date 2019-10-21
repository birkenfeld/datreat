program cerftest
implicit none
complex(kind=8) :: z, cer, cder
integer         :: i

do i=0,30
  z = cmplx(1d-40,0.5d0*i)
  cer = cerf(z)
  write(*,*)"CERF  : ",Imag(z), cer, cer*exp(z**2)
  call cerror(z,cer)
  write(*,*)"CERROR: ",Imag(z), cer, cer*exp(z**2)
  
enddo


do i=0,30
  z = cmplx(0.5d0*i,0d0)
  cer = cerf(z)
  write(*,*)"CERF  : ",Imag(z), cer, erf(Real(z))
  call cerror(z,cer)
  write(*,*)"CERROR: ",Imag(z), cer, erf(Real(z))
  
enddo


CONTAINS

function cerf(z) result(cer)
 implicit none
 complex(kind=8), intent(in) :: z
 complex(kind=8)             :: cer
 complex(kind=8)             :: cder

 call cerf_orig (z, cer, cder)

end function cerf

end program cerftest



subroutine cerf_orig ( z, cer, cder )

!*****************************************************************************80
!
!! CERF computes the error function and derivative for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ), the argument.
!
!    Output, complex ( kind = 8 ) CER, CDER, the values of erf(z) and erf'(z).
!
  implicit none

  complex ( kind = 8 ) c0
  complex ( kind = 8 ) cder
  complex ( kind = 8 ) cer
  complex ( kind = 8 ) cs
  real ( kind = 8 ) ei1
  real ( kind = 8 ) ei2
  real ( kind = 8 ) eps
  real ( kind = 8 ) er
  real ( kind = 8 ) er0
  real ( kind = 8 ) er1
  real ( kind = 8 ) er2 
  real ( kind = 8 ) eri
  real ( kind = 8 ) err
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) ss
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  complex ( kind = 8 ) z

  eps = 1.0D-12
  pi = 3.141592653589793D+00
  x = real ( z, kind = 8 )
  y = imag ( z )
  x2 = x * x

  if ( x <= 3.5D+00 ) then

    er = 1.0D+00
    r = 1.0D+00
    do k = 1, 100
      r = r * x2 / ( k + 0.5D+00 )
      er = er + r
      if ( abs ( er - w ) <= eps * abs ( er ) ) then
        exit
      end if
      w = er
    end do

    c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
    er0 = c0 * er

  else

    er = 1.0D+00
    r = 1.0D+00
    do k = 1, 12
      r = - r * ( k - 0.5D+00 ) / x2
      er = er + r
    end do
    c0 = exp ( - x2 ) / ( x * sqrt ( pi ) )
    er0 = 1.0D+00 - c0 * er

  end if

  if ( y == 0.0D+00 ) then
    err = er0
    eri = 0.0D+00
  else
    cs = cos ( 2.0D+00 * x * y )
    ss = sin ( 2.0D+00 * x * y )
    er1 = exp ( - x2 ) * ( 1.0D+00 - cs ) / ( 2.0D+00 * pi * x )
    ei1 = exp ( - x2 ) * ss / ( 2.0D+00 * pi * x )
    er2 = 0.0D+00
    do n = 1, 100
      er2 = er2 + exp ( - 0.25D+00 * n * n ) &
        / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x &
        - 2.0D+00 * x * cosh ( n * y ) * cs &
        + n * sinh ( n * y ) * ss )
      if ( abs ( ( er2 - w1 ) / er2 ) < eps ) then
        exit
      end if
      w1 = er2
    end do

    c0 = 2.0D+00 * exp ( - x2 ) / pi
    err = er0 + er1 + c0 * er2
    ei2 = 0.0D+00
    do n = 1, 100
      ei2 = ei2 + exp ( - 0.25D+00 * n * n ) &
        / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x &
        * cosh ( n * y ) * ss + n * sinh ( n * y ) * cs )
      if ( abs ( ( ei2 - w2 ) / ei2 ) < eps ) then
        exit
      end if
      w2 = ei2
    end do

    eri = ei1 + c0 * ei2

  end if

  cer = cmplx ( err, eri, kind = 8 )
  cder = 2.0D+00 / sqrt ( pi ) * exp ( - z * z )

  return
end subroutine cerf_orig







subroutine cerror ( z, cer )

!*****************************************************************************80
!
!! CERROR computes the error function for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CER, the function value.
!
  implicit none

  real ( kind = 8 ) a0
  complex ( kind = 8 ) c0
  complex ( kind = 8 ) cer
  complex ( kind = 8 ) cl
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) cs
  integer ( kind = 4 ) k
  real ( kind = 8 ) pi
  complex ( kind = 8 ) z
  complex ( kind = 8 ) z1

  a0 = abs ( z )
  c0 = exp ( - z * z )
  pi = 3.141592653589793D+00
  z1 = z

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    z1 = - z
  end if

  if ( a0 <= 5.8D+00 ) then    

    cs = z1
    cr = z1
    do k = 1, 120
      cr = cr * z1 * z1 / ( k + 0.5D+00 )
      cs = cs + cr
      if ( abs ( cr / cs ) < 1.0D-15 ) then
        exit
      end if
    end do

    cer = 2.0D+00 * c0 * cs / sqrt ( pi )

  else

    cl = 1.0D+00 / z1              
    cr = cl
    do k = 1, 13
      cr = -cr * ( k - 0.5D+00 ) / ( z1 * z1 )
      cl = cl + cr
      if ( abs ( cr / cl ) < 1.0D-15 ) then
        exit
      end if
    end do

    cer = 1.0D+00 - c0 * cl / sqrt ( pi )

  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    cer = -cer
  end if

  return
end
 
