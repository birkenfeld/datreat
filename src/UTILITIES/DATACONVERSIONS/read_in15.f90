program readin15 
implicit none
integer, parameter :: mtau=100
integer, parameter :: mq  =3
integer            :: nt  =0
integer            :: nq  =0
integer            :: numor 
integer            :: ioin, ios, ioout 
character(len=80)  :: fname
character(len=256) :: line
integer            :: i, i1, i2, j

double precision   :: phirous, l, nrous, nlin
double precision   :: tau(mtau)=0, q(0:mq)=0, sqt(mtau,0:mq)=0, sqt_err(mtau,0:mq)=0

read(*,*) numor, phirous, l, nrous, nlin
write(fname,'(i0,a)')numor,".txt"

open(newunit=ioin,file=trim(fname))
   read(ioin,'(a)') line
 ! extract q-values from first line  
   do i = 1,mq
     i1 = index(line,"EQ_")
     if(i1==0) stop "missing  Q"
     i1 = i1+3
     line = line(i1:)
     i2 = index(line," ")
     read(line(:i2),*) q(i)
     line = line(i2:)
   enddo
   q(0) = sum(q(1:3))/3d0
! write(*,*)"q: ",q

dr: do nt=1,mtau 
     read(ioin,*,iostat=ios) tau(nt), (sqt(nt,i),sqt_err(nt,i),i=0,3)
     if(ios.ne.0) exit dr
    enddo dr
    nt = nt-1

close(ioin)

! write(*,*)"q: ",q
! write(*,*)"t: ",tau(1:nt)
! write(*,*)"sqt: ",sqt(1:nt,0)

! write datreat compatible format
ioout = 6
do i=0,mq
  write(ioout,'(a,i0,a)')"in15 data numor:  ",numor," " 
  write(ioout,'(a,i0,a)')"in15da   sqt vs t ",numor," " 
  write(ioout,'(a)')"parameter"
  write(ioout,'(a,i0)')   "qstrip       ",i
  write(ioout,'(a,f12.6)')"q            ",q(i)
  write(ioout,'(a,f12.6)')"nrous        ",nrous
  write(ioout,'(a,f12.6)')"phirous      ",phirous
  write(ioout,'(a,f12.6)')"l            ",l
  write(ioout,'(a,f12.6)')"nlin         ",nlin
  write(ioout,'(a,f12.6)')"nlin_cc      ",nlin
  write(ioout,'(a,f12.6)')"philin       ",0.01
  write(ioout,'(a,f12.6)')"alin         ",0.0
  write(ioout,'(a,f12.6)')"arous        ",1.0
  write(ioout,'(a,f12.6)')"analytic     ",2.0
  write(ioout,'(a)')"   "
  write(ioout,'(3f16.8)') (tau(j),sqt(j,i),sqt_err(j,i), j=2,nt)
  write(ioout,*)
  write(ioout,'(a)')"#nxt"

enddo


end program readin15
