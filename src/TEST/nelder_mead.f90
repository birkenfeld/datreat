program match_scal
implicit none

integer, parameter   :: mpt=1000
double precision     :: x1(mpt), y1(mpt), y1err(mpt), y1d(mpt)
double precision     :: x2(mpt), y2(mpt), y2err(mpt), y2d(mpt)
double precision     :: mpar(2), errmin

character(len=128)   :: line
integer              :: i, n1, n2
integer              :: ios, unit1, unit2

open(newunit=unit1,file="/Users/monk/Desktop/SNS/KRUTYEVA/DENDRIMER/SANS/sHB2_1o3m.txt")
open(newunit=unit2,file="/Users/monk/Desktop/SNS/KRUTYEVA/DENDRIMER/SANS/sHB2_4o5m.txt")

read(unit1,'(a)',iostat=ios) line
write(*,*)"A",line
read(unit1,'(a)',iostat=ios) line
write(*,*)"B",line

i = 1
d1: do
     read(unit1,*,iostat=ios) x1(i), y1(i), y1err(i), y1d(i)
     if(ios .ne. 0) exit d1
     if(y1(i) .eq. 0) cycle d1
     i = i+1 
    enddo d1
    if(y1(i) .eq.0) n1 = i-1

read(unit2,'(a)',iostat=ios) line
write(*,*)"C",line
read(unit2,'(a)',iostat=ios) line
write(*,*)"D",line

i = 1
r2: do
     read(unit2,*,iostat=ios) x2(i), y2(i), y2err(i), y2d(i)
     if(ios .ne. 0) exit r2
     if(y2(i) .eq. 0) cycle r2
     i = i+1 
    enddo r2
    if(y2(i) .eq.0) n2 = i-1

write(*,*)"x1"
write(*,'(3e15.7)')(x1(i),y1(i),y1err(i),i=1,n1)

write(*,*)"x2"
write(*,'(3e15.7)')(x2(i),y2(i),y2err(i),i=1,n2)


! call of the test


mpar = [0d0, 0.1d0] * i
call smatch(x1,y1,y1err,n1, x2,y2,y2err,n2, mpar, errmin )

write(*,'("errmin, mpar=",3e14.7)')errmin, mpar

 
end program match_scal


subroutine smatch(x1,y1,y1err,n1, x2,y2,y2err,n2, mpar, errmin )
  implicit none
  double precision, intent(in ) :: x1(n1)      ! x-values of master record
  double precision, intent(in ) :: y1(n1)      ! y-values of master record
  double precision, intent(in ) :: y1err(n1)   ! y-error  of master record
  integer         , intent(in ) :: n1          ! nuber of values in master
  double precision, intent(in ) :: x2(n2)      ! x-values of slave record
  double precision, intent(in ) :: y2(n2)      ! y-values of slave record
  double precision, intent(in ) :: y2err(n2)   ! y-error  of slave record
  integer         , intent(in ) :: n2          ! nuber of values in slave
  double precision, intent(inout) :: mpar(2)     ! offset and scale ....
  double precision, intent(out) :: errmin      ! residual matching error
  
  double precision :: mpar_min(size(mpar))
  double precision :: step(size(mpar))
  double precision :: reqmin = 1d-7
  integer          :: konvge = 1
  integer          :: kcount = 10000
  integer          :: icount
  integer          :: numres
  integer          :: ifault


  step   = 0.05d0
  errmin = match_err(mpar, size(mpar))

  call nelmin ( match_err, size(mpar), mpar, mpar_min, errmin, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

  mpar = mpar_min


contains 

function match_err(sp, n ) result(val) 
  implicit none
  double precision, intent(in ) :: sp(n)
  integer         , intent(in ) :: n 

  double precision :: val
  double precision :: yinterp, yinterp_err, xt, yt, yte, p
  integer          :: i, ma, mb, ncompare

  !! go through the points of slave and compare with master
  val      =  0
  ncompare = 0
ds1: do i=1,n2
       xt  = x2(i)
       yt  = (y2(i) + sp(1))* sp(2)
       yte =  y2err(i)      * sp(2)
       ma = minloc(abs(x1(1:n1)-xt),dim=1)
       mb = min(ma + 1,n1)
       if(mb>n1) mb = ma-1
       p             = (xt - x1(ma))/(x1(mb) - x2(ma))
       if(abs(p) > 2d0) cycle ds1
       yinterp      = p * y1(ma) + (1d0-p) * y1(mb)
       yinterp_err  = sqrt(  (p*y1err(ma))**2  +((1d0-p)*y1err(mb))**2 )
       ncompare = ncompare + 1
       val = val + (yinterp-yt)**2 / (yinterp_err**2 + yte**2)
     enddo ds1
 
     if(ncompare > 0) then
       val = val / ncompare
     else
       val = Huge(val)
     endif
 
!!TP: write(*,'(3f12.6)') sp, val

end function match_err

end subroutine smatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Suroutines from external sources                                          !!
!!!! ASA047                                                                    !!
!!!! Nelder-Mead Minimization                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
  icount, numres, ifault )

!*****************************************************************************
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, n )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
  real ( kind = 8 ) del
  real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
  real ( kind = 8 ), parameter :: eps = 0.001D+00
  real ( kind = 8 ), external :: fn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcount
  integer ( kind = 4 ) kcount
  integer ( kind = 4 ) konvge
  integer ( kind = 4 ) l
  integer ( kind = 4 ) numres
  real ( kind = 8 ) p(n,n+1)
  real ( kind = 8 ) p2star(n)
  real ( kind = 8 ) pbar(n)
  real ( kind = 8 ) pstar(n)
  real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  real ( kind = 8 ) reqmin
  real ( kind = 8 ) rq
  real ( kind = 8 ) start(n)
  real ( kind = 8 ) step(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin(n)
  real ( kind = 8 ) y(n+1)
  real ( kind = 8 ) y2star
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ynewlo
  real ( kind = 8 ) ystar
  real ( kind = 8 ) z
!
!  Check the input parameters.
!
  if ( reqmin <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( n < 1 ) then
    ifault = 1
    return
  end if

  if ( konvge < 1 ) then
    ifault = 1
    return
  end if
!
!  Initialization.
!
  icount = 0
  numres = 0
  jcount = konvge
  del = 1.0D+00
  rq = reqmin * real ( n, kind = 8 )
!
!  Initial or restarted loop.
!
  do

    p(1:n,n+1) = start(1:n)
    y(n+1) = fn ( start, n )
    icount = icount + 1
!
!  Define the initial simplex.
!
    do j = 1, n
      x = start(j)
      start(j) = start(j) + step(j) * del
      p(1:n,j) = start(1:n)
      y(j) = fn ( start, n )
      icount = icount + 1
      start(j) = x
    end do
!
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
!
    ilo = minloc ( y(1:n+1), 1 )
    ylo = y(ilo)
!
!  Inner loop.
!
    do while ( icount < kcount )
!
!  YNEWLO is, of course, the HIGHEST value???
!
      ihi = maxloc ( y(1:n+1), 1 )
      ynewlo = y(ihi)
!
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
!
      do i = 1, n
        pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
      end do
!
!  Reflection through the centroid.
!
      pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
      ystar = fn ( pstar, n )
      icount = icount + 1
!
!  Successful reflection, so extension.
!
      if ( ystar < ylo ) then

        p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
        y2star = fn ( p2star, n )
        icount = icount + 1
!
!  Retain extension or contraction.
!
        if ( ystar < y2star ) then
          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
        else
          p(1:n,ihi) = p2star(1:n)
          y(ihi) = y2star
        end if
!
!  No extension.
!
      else

        l = 0
        do i = 1, n + 1
          if ( ystar < y(i) ) then
            l = l + 1
          end if
        end do

        if ( 1 < l ) then

          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
!
!  Contraction on the Y(IHI) side of the centroid.
!
        else if ( l == 0 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
          y2star = fn ( p2star, n )
          icount = icount + 1
!
!  Contract the whole simplex.
!
          if ( y(ihi) < y2star ) then

            do j = 1, n + 1
              p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
              xmin(1:n) = p(1:n,j)
              y(j) = fn ( xmin, n )
              icount = icount + 1
            end do

            ilo = minloc ( y(1:n+1), 1 )
            ylo = y(ilo)

            cycle
!
!  Retain contraction.
!
          else
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          end if
!
!  Contraction on the reflection side of the centroid.
!
        else if ( l == 1 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
          y2star = fn ( p2star, n )
          icount = icount + 1
!
!  Retain reflection?
!
          if ( y2star <= ystar ) then
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          else
            p(1:n,ihi) = pstar(1:n)
            y(ihi) = ystar
          end if

        end if

      end if
!
!  Check if YLO improved.
!
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      end if

      jcount = jcount - 1

      if ( 0 < jcount ) then
        cycle
      end if
!
!  Check to see if minimum reached.
!
      if ( icount <= kcount ) then

        jcount = konvge

        x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
        z = sum ( ( y(1:n+1) - x )**2 )

        if ( z <= rq ) then
          exit
        end if

      end if

    end do
!
!  Factorial tests to check that YNEWLO is a local minimum.
!
    xmin(1:n) = p(1:n,ilo)
    ynewlo = y(ilo)

    if ( kcount < icount ) then
      ifault = 2
      exit
    end if

    ifault = 0

    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn ( xmin, n )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) - del - del
      z = fn ( xmin, n )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) + del
    end do

    if ( ifault == 0 ) then
      exit
    end if
!
!  Restart the procedure.
!
    start(1:n) = xmin(1:n)
    del = eps
    numres = numres + 1

  end do

  return
end
