program read_ino
integer, parameter :: mchannel = 1000
integer, parameter :: mblocks  = 300
character(len=128) :: filename = "239953_239956_sqw_5K_dPBcross.ino"

double precision   :: minimal_counts = 1d0
character(len=128) :: line
character(len=128) :: comment
integer            :: ios, iu
integer            :: nchannel, nblock
integer            :: last_block = 0 

double precision   :: x__val(mchannel,mblocks)
double precision   :: y__val(mchannel,mblocks)
double precision   :: errval(mchannel,mblocks)
integer            :: meta1(8,mblocks)
double precision   :: meta2(5,mblocks)
double precision   :: meta3(3,mblocks)

double precision, parameter :: Pi                  = 4d0*atan(1d0)
double precision, parameter :: Planckkonstante     =     6.62606957d-34
double precision, parameter :: Elektronenladung    =     1.602176565d-19
double precision, parameter :: conversion_meV_GHz  = 1d-3*Elektronenladung*2*Pi/Planckkonstante/1d9


integer            :: length, status, i, j
character(len=256) :: val

character(len=16)  :: nam

call get_command_argument(1,val,length, status)
if(status == 0) filename = val
write(*,*)"Reading from file: ",trim(filename)
 

open(newunit=iu,file=trim(filename),iostat=ios)
if(ios .ne. 0) stop "File not found"

nblock = 0
dbloc: do j=1,mblocks
  read(iu,'(a)',iostat=ios) line
!write(*,*)"t1",j,nblocj,trim(line)
  if(ios .ne. 0) exit dbloc
  nblock = nblock+1
  read(line,*) meta1(:,nblock)
  nchannel = meta1(1,nblock)-3
  read(iu,'(a)') comment
! write(*,*) nblock,"   ",trim(comment)
  read(iu,'(a)',iostat=ios) line
!write(*,*)"t2",j,nblocj,trim(line)
  read(line,*) meta2(:,nblock)
  read(iu,'(a)',iostat=ios) line
!write(*,*)"t3",j,nblocj,trim(line)
  read(line,*) meta3(:,nblock)

  dcha: do ichannel=1,nchannel
   read(iu,*) x__val(ichannel,nblock),y__val(ichannel,nblock),errval(ichannel,nblock)
  enddo dcha


  if(maxval(y__val(1:nchannel,nblock)) <= minimal_counts) nblock = nblock-1
  
  if(nblock < 1  .or. nblock == last_block ) cycle dbloc

  last_block = nblock

  write(*,'("Block: ",i8,"  total sum=",e14.7,2x,8f12.6,2x3f12.6)')  &
           nblock, sum(y__val(1:nchannel,nblock)), meta2(:,nblock), meta3(:,nblock)

enddo dbloc

  if(nblock > 0) then
    write(*,*) nblock,"   ",trim(comment)
  else
    write(*,*) "no file or file with valid data found !"
    write(*,*) "specify a filename following the command !"
  endif

close(iu)


!! -- now write datreat format files --
!! create a name
 j = index(filename,".ino")
 if(j > 0) then
   filename = filename(1:j)//"dtr"
 else
   filename = "datreat_file.dtr2"
 endif

 write(*,*)"writing datreat file:", "d"//trim(filename)


open(newunit=iu,file="d"//trim(filename))

i   = index(filename,"_sqw_")
nam = filename(i+5:j-1)

write(*,*) nam

do i=1,nblock
  write(iu,'(a,a,a)') trim(nam)," ",trim(comment)
  write(iu,'(a,a,i8)')trim(nam),"  S(q,om)  vs om/GHz  ",nblock+1000
  write(iu,'(a,f12.6)')"  q       ",meta2(1,i)
  do j=1,size(meta2(:,1))
    write(iu,'(a,i0,3x,e14.7)')"  meta2_",j,meta2(j,i)
  enddo
  do j=1,size(meta3(:,1))
    write(iu,'(a,i0,3x,e14.7)')  "meta3_",j,meta3(j,i)
  enddo
  write(iu,*)
  do j=1,nchannel
     if(errval(j,i) <= 0d0) cycle
     write(iu,'(3e18.7)')  x__val(j,i) *  conversion_meV_GHz   ,y__val(j,i), errval(j,i)
  enddo
  write(iu,*)
  write(iu,'(a)')"#nxt"

enddo

close(iu)


end program read_ino
