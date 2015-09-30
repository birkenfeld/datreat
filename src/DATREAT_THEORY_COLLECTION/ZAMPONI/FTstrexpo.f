       implicit none
       real*8 a,b, domega, o0
       real*8 adapint, sum, result, result2, err, erraccu
       integer maxit,i
       real*8 strex_kernel
       external strex_kernel


       real*8 Omega, str_beta, str_delta
       common /cstrex1/Omega, str_beta, str_delta

       real*8 Pi
       parameter(Pi=3.141592654d0)

       str_beta  = 0.8
       str_delta = 0.8

       domega    = 0.01d0

       maxit = 2000
       err   = 1d-8

       a = 0
       b = 5.0d0/str_delta


       sum = 0  
       open(10,file='st.dgli',status='unknown')  
       do i=0,1000
        Omega  = domega*i
        result  = adapint(strex_kernel,a,b,err,maxit,erraccu)*2/Pi
        result2 = adapint(strex_kernel,a,2*b,err,maxit,erraccu)*2/Pi
!       write(6,*) Omega, result, result2
        write(10,*) Omega, result, result2
        sum = sum+domega*result2 
      enddo
       O0 = Omega  
       domega = 10*domega
       do i=1,2000
        Omega  = O0+domega*i
        result  = adapint(strex_kernel,a,b,err,maxit,erraccu)*2/Pi
        result2 = adapint(strex_kernel,a,2*b,err,maxit,erraccu)*2/Pi
!       write(6,*) Omega, result, result2
        write(10,*) Omega, result, result2
        sum = sum+domega*result2 
      enddo

       close(10)      

       write(6,*)'sum=',sum

       stop
       end



       function strex_kernel(t)
!      ------------------------
       implicit none

       real*8 strex_kernel, t
       real*8 Omega, str_beta, str_delta
       common /cstrex1/Omega, str_beta, str_delta

       strex_kernel=
     * exp(-t**str_beta-1d0/4d0*(str_delta*t)**2)*cos(t*Omega)    

       return
       end




 
       function gaus8a (f,xu,xo)
c-----------------------------------------------------------------------
c      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
c-----------------------------------------------------------------------
       parameter ( ndim = 8 )
       parameter ( ndim2=ndim/2)
       implicit real*8 (a-h,o-z)
       dimension a(ndim2),x(ndim2)
       data a / 0.362683783378362d0,
     *          0.313706645877887d0,
     *          0.222381034453374d0,
     *          0.101228536290376d0/
       data x / 0.183434642495650d0,
     *          0.525532409916329d0,
     *          0.796666477413627d0,
     *          0.960289856497536d0/
 
       xave = (xo+xu)*0.5d0
       range= (xo-xu)*0.5d0
       sum = 0.d0
       do i=1,ndim2
         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) )
       enddo
       gaus8a = sum * range
 
       return
       end
c*ds
c*ed
       function adapint(f,a,b,epsilon,maxiter,erroraccu)
c      =================================================                --
c
c      lokal adaptives integrationsverfahren 2te ver. wg. rekur.
c
       implicit real*8 (a-h,o-z)
c
       parameter( maxstack = 200 )
c                 --------------->  stacktiefe
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack)
c
       logical cray
!!     common/outlev/iot,ibild,ierrs,inka, cray
       real*4 xxxx,yyyy,ptxf                                         !! aix
       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)  !! aix
c
       external f
c
c      -------------------------------------------------------------------
 
 
         iterationcounter = 0
         erg              = 0.0d0
         itop      = 0
         xa(itop)  = a
         xb(itop)  = b
         s(itop)   = gaus8a( f, a, b )
 
         do i=1,maxiter
           if (itop.ge.(maxstack-2)) then
              erg = erg + s(itop)
              xbb = xb(itop)
              itm = itop-1
              do itop=itm,0,-1
                 if(xa(itop).eq.xbb) goto 1
              enddo
1             continue
              write(6,*)'warning! adaptint stack overflow!'
           else
              iterationcounter = iterationcounter + 1
              itop             = itop +1
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0
              xb(itop) =  xb(itop-1)
              s(itop)  = gaus8a( f, xa(itop), xb(itop) )
              itop             = itop +1
              xa(itop) =  xa(itop-2)
              xb(itop) =  xa(itop-1)
              s(itop)  = gaus8a( f, xa(itop), xb(itop) )
              error   =  dabs((s(itop)+s(itop-1))-s(itop-2))
 
              if (iot.gt.2) then
                 write(6,'(1x,i3,i5,4e13.6)')itop,iterationcounter,
     *                 xa(itop),xb(itop),(s(itop)+s(itop-1)),error
              endif
 
              if (error.lt.epsilon) then
                 erg = erg + s(itop)+s(itop-1)
                 erroraccu = erroraccu + error
                 itop = itop-1
                 xbb = xb(itop)
                 itop = itop-1
                 itm = itop-1
                 do itop=itm,0,-1
                    if(xa(itop).eq.xbb) goto 2
                 enddo
2                continue
              endif
            endif
            if (itop.le.0) goto 3
         enddo
         write(6,*)                                                     '
     *    'adapint fatal error iterationnumber exceeded!'
3        continue
 
 
         adapint = erg
 
         return
 
       end
 
