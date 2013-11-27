       program test_get_pair
         
        real :: x1(100),y1(100),e1(100)
        real :: x2(100),y2(100),e2(100)
        real :: x3(100),y3(100),e3(100)
        real :: y4(100),e4(100)

        integer :: n1, n2, nout, i


        n1 = 10
        n2 = 10
        
       do i=1,n1

         x1(i) = i
         x2(i) = i
         y1(i) = i
         y2(i) = i+0.5
         e1(i) = 0.1*i
         e2(i) = 0.3

       enddo


       call get_pair_of_points(x1,y1,e1,n1,x2,y2,e2,n2,x3,y3,e3,y4,e4,nout)

       write(6,*)'nout = ',nout
       do i=1,n1
         write(6,'(11f10.4)')x1(i),y1(i),e1(i),x2(i),y2(i),e2(i),x3(i),y3(i),e3(i),y4(i),e4(i)
       enddo


     write(6,*)'====================================================='

      do i=1,n1

         x1(i) = i
         x2(i) = i+0.1
         y1(i) = i
         y2(i) = i+0.5
         e1(i) = 0.1*i
         e2(i) = 0.3

       enddo


       call get_pair_of_points(x1,y1,e1,n1,x2,y2,e2,n2,x3,y3,e3,y4,e4,nout)

       write(6,*)'nout = ',nout
       do i=1,n1
         write(6,'(11f10.4)')x1(i),y1(i),e1(i),x2(i),y2(i),e2(i),x3(i),y3(i),e3(i),y4(i),e4(i)
       enddo



     write(6,*)'====================================================='

      do i=1,n1

         x1(i) = i
         x2(i) = i-0.1
         y1(i) = i
         y2(i) = i+0.5
         e1(i) = 0.1*i
         e2(i) = 0.3

       enddo


       call get_pair_of_points(x1,y1,e1,n1,x2,y2,e2,n2,x3,y3,e3,y4,e4,nout)

       write(6,*)'nout = ',nout
       do i=1,n1
         write(6,'(11f10.4)')x1(i),y1(i),e1(i),x2(i),y2(i),e2(i),x3(i),y3(i),e3(i),y4(i),e4(i)
       enddo


     write(6,*)'====================================================='

      do i=1,n1

         x1(i) = i
         x2(i) = i-0.2-0.005*i*i
         y1(i) = i*i
         y2(i) = i+0.5
         e1(i) = 0.1*i
         e2(i) = 0.3

       enddo


       call get_pair_of_points(x1,y1,e1,n1,x2,y2,e2,n2,x3,y3,e3,y4,e4,nout)

       write(6,*)'nout = ',nout
       do i=1,n1
         write(6,'(11f10.4)')x1(i),y1(i),e1(i),x2(i),y2(i),e2(i),x3(i),y3(i),e3(i),y4(i),e4(i)
       enddo



       end program test_get_pair




       subroutine get_pair_of_points(xin1,yin1,yer1,n1, xin2,yin2,yer2,n2, xout, yout1, yerout1, yout2, yerout2, nout)
!      ===============================================================================================================
!
       implicit none

       real, intent(in)     :: xin1(n1), yin1(n1), yer1(n1)
       real, intent(in)     :: xin2(n2), yin2(n1), yer2(n2)
       integer, intent(in)  :: n1, n2

       real, intent(out)    :: xout(nout), yout1(nout), yerout1(nout), yout2(nout), yerout2(nout) 
       integer              :: nout
         

       integer              :: i, j, jc, ii
       real                 :: x, xx, y1, y2, e1, e2, dist,p,y,ye
 
       nout = 0
! assuming nontrivial length check
       if(n1.lt.2 .or. n2.lt.2) then
         write(6,*)'at least one of the vectors has too few components:',n1,n2
         return
       endif
! assuming equal length check
       if(n1.ne.n2) then
         write(6,*)'vectors have different number of components:',n1,n2
         return
       endif
! assuming ordered vectors CHECK
       do i=1,n1-1
        if(xin1(i).ge.xin1(i+1)) then
          write(6,*)'vector 1 is not ordered '
          return
        endif
       enddo
       do i=1,n2-1
        if(xin2(i).ge.xin2(i+1)) then
          write(6,*)'vector 2 is not ordered '
          return
        endif
       enddo

! do the assignments !
       do i=1,n1
!      search for matching points       
!      1. closest mate
         dist = abs(xin1(i)-xin2(1))
         jc   = 1
         do j=2,n2
           xx = abs(xin1(i)-xin2(j))
           if(xx.lt.dist) then
              dist = xx
              jc   = j
           endif
         enddo


!        now our closest neighbour pair is i<-->jc
!        if they are equal we are done
         if(xin1(i).eq.xin2(jc)) then
           nout = nout + 1
           xout(nout)    = xin1(i)
           yout1(nout)   = yin1(i)
           yerout1(nout) = yer1(i)
           yout2(nout)   = yin2(jc) 
           yerout2(nout) = yer2(jc)
           cycle
        endif

        if(xin1(i).lt.xin2(jc)) then
          do ii=1,n1
            if(xin1(ii).gt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)  !! this would be good if full statistical independeence
                                                           !! is assumed, however interpolation creates correlation
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif

        if(xin1(i).gt.xin2(jc)) then
          do ii=n1,1,-1
            if(xin1(ii).lt.xin2(jc)) then
             exit
            endif
          enddo
          p  = (xin2(jc)-xin1(i))/(xin1(ii)-xin1(i))
          x  = xin2(jc)
          y  = yin1(ii)*p + yin1(i)*(1d0-p)
          ye = sqrt((yer1(ii)*p)**2+(yer1(i)*(1d0-p))**2)
          nout = nout+1
          xout(nout)    = x
          yout1(nout)   = y
          yerout1(nout) = ye
          yout2(nout)   = yin2(jc)
          yerout2(nout) = yer2(jc)
          cycle
        endif
         
       enddo
       

       return
       end subroutine get_pair_of_points



