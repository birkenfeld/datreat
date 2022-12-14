program dynamic
implicit none

!        real xwerte(c_MWERT,c_MBUF)
!        real ywerte(c_MWERT,c_MBUF)
!        real yerror(c_MWERT,c_MBUF)
!            character*80 xname(c_MBUF)
!        character*80 yname(c_MBUF)
!        character*80 name(c_MBUF)
!        integer nwert(c_MBUF)
!        integer numor(c_MBUF)
!        integer :: nbuf = 0
!        integer :: iw_perm(c_MWERT)
! 
!        character*80 coment(c_MBUF)*80
!        real params(c_MPAR,c_MBUF)
!        character*80 napar(c_MPAR,c_MBUF)
!        integer nopar(c_MBUF)
!                integer :: params_display_level(c_MPAR,c_MBUF) = 0
!


type dat
  double precision, allocatable :: x(:)
  double precision, allocatable :: y(:)
end type dat

double precision, pointer :: xv(:)


!type(dat), target :: d(100)
type(dat), target, allocatable :: d(:)

integer :: i

allocate(d(100))

write(*,*) size(d)

allocate(d(1)%x(20))
allocate(d(1)%y(20))


allocate(d(2)%x(200))
allocate(d(2)%y(200))

d(1)%x = [(i,i=1,size(d(1)%x))]
d(1)%y = [(i,i=1,size(d(1)%y))]

xv => d(1)%x 

d(2)%x = [(i,i=1,size(d(2)%x))]
d(2)%y = [(i,i=1,size(d(2)%y))]

write(*,*) size(d), size(d(1)%x), size(d(2)%x),xv


end program dynamic 
