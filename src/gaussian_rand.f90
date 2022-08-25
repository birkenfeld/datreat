program testrg
implicit none

double precision :: gauss_rand
double precision :: rg, s1, s2
integer          :: i, n=300000

integer          :: m = 10, md=20
double precision :: x(1000), x0(1000)

call init_random_seed()
!! 
!! s1 = 0
!! s2 = 0
!! do i=1,n
!!  rg = gauss_rand()
!!  s1 = s1+rg
!!  s2 = s2+rg**2
!!  if(mod(i,1000)==0) then
!!    write(6,*) i, rg, s1/i, s2/i
!!  endif
!! enddo
!! 
call random_ndimensional_direction ( md, x0 )

do i=1,m
  call random_ndimensional_direction ( md, x )
  write(6,'(32f8.5)') x(1:md), dot_product(x(1:md),x(1:md)), dot_product(x(1:md),x0(1:md)) 
enddo
end program testrg


subroutine random_ndimensional_direction ( n, x )
  implicit none
  integer         , intent(in)  :: n
  double precision, intent(out) :: x(n)

  integer          :: i
  double precision :: xn
  double precision :: gauss_rand

  xn = 0
  do i=1,n
    x(i) = gauss_rand()
    xn   = xn + x(i)**2
  enddo
  x = x / sqrt(xn)
end subroutine random_ndimensional_direction 



!--------------------------------------
 double precision function gauss_rand() 
!--------------------------------------

  implicit none

  double precision ::  r1, r2, r, f, grand1, grand2, now, rnd(2)

  
    do
      CALL RANDOM_NUMBER(rnd)
      r1 = rnd(1)*2.d0-1.d0
      r2 = rnd(2)*2.d0-1.d0
      r  = r1*r1 +r2*r2
      if(( (r<1.d0) .and. (r .ne. 0.d0) )) exit
    enddo
 
    f = sqrt(-2.d0*log(r)/r)
    gauss_rand = r1*f
    
 
  end function gauss_rand 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! gnu example:

          subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un=99, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)

            write(6,*)'Random seed: ',seed

          end subroutine init_random_seed




