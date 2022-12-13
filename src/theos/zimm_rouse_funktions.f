!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ZIMM function subroutines                                           !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       real*8 function h_zimm(u,epsilon)

         implicit none
         real*8  dDsi, u, Aintsep, epsilon, erraccu, Pi
         real    Si
         real*8  Part1, Part2, x
         real*8  a2dapint
         real*8  uc
         common/hzkc/ uc

         real*8 axp, epsilon_c
         real*8 lgeps
         common /czimrou/ axp, epsilon_c

         integer maxit
         Parameter (Pi=3.141592654d0)

         external h_zimm_kernel 

         maxit   = 1300
         lgeps   = dlog(epsilon)
         if(lgeps.gt.0.d0) lgeps = -0.1d0 
         Aintsep = (-lgeps*sqrt(2.0d0))**axp
 
         if(u.lt.1d-7) u = 1d-7
         uc = u

         x = Aintsep * u
         Part2 = (dcos(x)+dDSi(x)*x)/Aintsep - Pi*u/2
!         Part2 = (dcos(x)+Si(sngl(x))*x)/Aintsep - Pi*u/2

         Part1 = a2dapint(h_zimm_kernel, 0.0d0, Aintsep,
     *                    epsilon, maxit, erraccu)

         h_zimm = (Part1+Part2)*2/Pi

       return
       end

      
       real*8 function h_zimm_kernel(x)

         implicit none
         real*8  uc
         common/hzkc/ uc
       
         real*8 x, axpi, sq2
         Parameter(Sq2=1.414213562d0)
         real*8 axp, epsilon_c
         common /czimrou/ axp, epsilon_c

         if(x.lt.1d-7) x=1.0d-7
         axpi = 1.0d0/axp

         h_zimm_kernel = dcos(x*uc)*(1.0d0-dexp(-(x**axpi)/Sq2))/(x*x)
 
       return
       end 

       subroutine create_zimm_table(epsilon2)
!!     =====================================

        implicit none
        real*8 t, Q, temp, eta,kb,Pi
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu

        real*8 SQ_zimm_kernel
        real*8 GamktAxp
        common/sqzkc/GamktAxp 
        real*8 axp, epsilon_c
        common /czimrou/ axp, epsilon_c
        integer maxit2 

        integer nzimtab
        Parameter(nzimtab=200) 
        real*8 dgy, zimm_tab(nzimtab)
        common/czimtab/dgy, zimm_tab 

        integer i

        external SQ_zimm_kernel

 
        epsilon_c = epsilon2
        Bint_end = -dlog(epsilon_c)
        maxit2   = 1300

        Axp = 2.0d0/3.0d0


!        t_SI = t * 1d-9
!        Q_SI = Q * 1d10

!        Gamkt = (kb * temp)/(6*Pi*eta) 
!        Gamkt = Gamkt*Q_SI*Q_SI*(Q_SI*t_SI)
!        GamktAxp = Gamkt**Axp


        dgy = 0.025d0
        zimm_tab(1) = 1.0d0

        do i=2,nzimtab
          GamktAxp = i*dgy

          zimm_tab(i) = adapint( SQ_zimm_kernel, 0d0, Bint_end,
     *                      epsilon2, maxit2, erraccu)
          
          write(6,*)i,GamktAxp,zimm_tab(i)
        enddo
     
       return
       end


       real*8 function zimm_table(y)
!!     =====================================

        implicit none

        integer nzimtab
        Parameter(nzimtab=200) 
        real*8 dgy, zimm_tab(nzimtab)
        common/czimtab/dgy, zimm_tab 

        real*8 ya,y,h, yout,arg
        integer i1,i2,i3
     
        ya = dabs(y)
        i1 = ya/dgy + 1
        i2 = i1+1
        i3 = i1+2

        if(i2.ge.nzimtab) then
          arg = -1.35*ya
          if(arg.lt.-50.0d0) arg = -50.0d0
          zimm_table = dexp(arg)
          return
        endif

        

        yout =  zimm_tab(i1)*(ya-i3*dgy)*(ya-i2*dgy)/((i1-i3)*(i1-i2))
     *         +zimm_tab(i2)*(ya-i1*dgy)*(ya-i3*dgy)/((i2-i1)*(i2-i3))
     *         +zimm_tab(i3)*(ya-i1*dgy)*(ya-i2*dgy)/((i3-i1)*(i3-i2))
        yout = yout/(dgy*dgy)

        zimm_table = yout

       return
       end

                

       real*8 function SQ_zimm(t,Q,temp,eta,epsilon2)
!!     ==============================================
!!  from scratch
!!

        implicit none
        real*8 t, Q, temp, eta,kb,Pi
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu

        real*8 SQ_zimm_kernel
        real*8 GamktAxp
        common/sqzkc/GamktAxp 
        real*8 axp, epsilon_c
        common /czimrou/ axp, epsilon_c
        integer maxit2 

        external SQ_zimm_kernel

 
        epsilon_c = epsilon2
        Bint_end = -dlog(epsilon_c)
        maxit2   = 1300

        Axp = 2.0d0/3.0d0


        t_SI = t * 1d-9
        Q_SI = Q * 1d10

        Gamkt = (kb * temp)/(6*Pi*eta) 
        Gamkt = Gamkt*Q_SI*Q_SI*(Q_SI*t_SI)
        GamktAxp = Gamkt**Axp

        SQ_zimm = adapint( SQ_zimm_kernel, 0d0, Bint_end,
     *                    epsilon2, maxit2, erraccu)
     
       return
       end

       real*8 function SQ_zimmT(t,Q,temp,eta,epsilon2)
!!     ===============================================
!!     from table
!!
        implicit none
        real*8 t, Q, temp, eta,kb,Pi
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu

        real*8 SQ_zimm_kernel
        real*8 GamktAxp
        common/sqzkc/GamktAxp 
        real*8 axp, epsilon_c
        common /czimrou/ axp, epsilon_c
        integer maxit2 
        real*8 zimm_table

 
        if(epsilon2.ne.epsilon_c) then
          call create_zimm_table(epsilon2)
        endif


        t_SI = t * 1d-9
        Q_SI = Q * 1d10

        Gamkt = (kb * temp)/(6*Pi*eta) 
        Gamkt = Gamkt*Q_SI*Q_SI*(Q_SI*t_SI)
        GamktAxp = Gamkt**Axp

        SQ_zimmT = zimm_table(GamktAxp)

     
       return
       end

                
       real*8 function  SQ_zimm_kernel(u)
!!     ----------------------------------

        implicit none
        real*8 u, arg
        real*8 h_Zimm
        real*8 epseff

        real*8 GamktAxp
        common/sqzkc/GamktAxp 
        real*8 axp, epsilon_c
        common /czimrou/ axp, epsilon_c

        if(GamktAxp.lt.1d-7) GamktAxp = 1d-7 
        epseff =0.1d0*epsilon_c*dexp(u)/GamktAxp       

        arg = u + GamktAxp * h_Zimm( u / GamktAxp, epseff)

        if(arg.gt.50.0d0) arg = 50.0d0

        SQ_zimm_kernel = dexp(-arg)

        return
        end


 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ROUSE function subroutines                                          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       real*8 function h_rouse(u,epsilon)

         implicit none
         real*8  dDsi, u, Aintsep, epsilon, erraccu, Pi
         real    Si
         real*8  Part1, Part2, x
         real*8  a2dapint
         real*8  uc
         common/hroukc/ uc

         real*8 axp, epsilon_c
         real*8 lgeps
         common /crouse/ axp, epsilon_c

         integer maxit
         Parameter (Pi=3.141592654d0)

         external h_rouse_kernel 

         maxit   = 1300
         lgeps   = dlog(epsilon)
         if(lgeps.gt.0.d0) lgeps = -0.1d0 
         Aintsep = (-lgeps)**axp
  
         if(u.lt.1d-7) u = 1d-7
         uc = u

         x = Aintsep * u
         Part2 = (dcos(x)+dDSi(x)*x)/Aintsep - Pi*u/2
!         Part2 = (dcos(x)+Si(sngl(x))*x)/Aintsep - Pi*u/2

         Part1 = a2dapint(h_rouse_kernel, 0.0d0, Aintsep,
     *                    epsilon, maxit, erraccu)

         h_rouse = (Part1+Part2)*2/Pi

       return
       end

      
       real*8 function h_rouse_kernel(x)

         implicit none
         real*8  uc
         common/hroukc/ uc
       
         real*8 x, axpi, sq2
         Parameter(Sq2=1.414213562d0)
         real*8 axp, epsilon_c
         common /crouse/ axp, epsilon_c

         if(x.lt.1d-7) x=1.0d-7
         axpi = 1.0d0/axp

         h_rouse_kernel = dcos(x*uc)*(1.0d0-dexp(-(x**axpi)))/(x*x)
 
       return
       end 

       subroutine create_rouse_table(epsilon2)
!!     =====================================

        implicit none
        real*8 t, Q, temp, eta,kb,Pi
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu

        real*8 SQ_rouse_kernel
        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c
        integer maxit2 

        integer nroutab
        Parameter(nroutab=1000) 
        real*8 dgy, rouse_tab(nroutab)
        common/croutab/dgy, rouse_tab 

        integer i

        external SQ_rouse_kernel

 
        epsilon_c = epsilon2
        Bint_end = -dlog(epsilon_c)
        maxit2   = 1300

        Axp = 1.0d0/2.0d0


!        t_SI = t * 1d-9
!        Q_SI = Q * 1d10

!        Gamkt = (kb * temp)/(6*Pi*eta) 
!        Gamkt = Gamkt*Q_SI*Q_SI*(Q_SI*t_SI)
!        GamktAxp = Gamkt**Axp


        dgy = 0.025d0
        rouse_tab(1) = 1.0d0

        do i=2,nroutab
          GamktAxp = i*dgy

          rouse_tab(i) = adapint( SQ_rouse_kernel, 0d0, Bint_end,
     *                      epsilon2, maxit2, erraccu)
          
          write(6,*)i,GamktAxp,rouse_tab(i)
        enddo
     
       return
       end


       real*8 function rouse_table(y)
!!     =====================================

        implicit none

        integer nroutab
        Parameter(nroutab=1000) 
        real*8 dgy, rouse_tab(nroutab)
        common/croutab/dgy, rouse_tab 

        real*8 ya,y,h, yout,arg
        integer i1,i2,i3
     
        ya = abs(y)
        i1 = ya/dgy + 1
        i2 = i1+1
        i3 = i1+2

!        write(6,*)'routab:', i1,i2,i3, y, dgy

        if(i3.ge.nroutab) then
          arg = -1.77*ya
          if(arg.lt.-50.0d0) arg = -50.0d0
          rouse_table = dexp(arg)
          return
        endif
        if(i1.lt.0) then
          rouse_table = 1d0
          return
        endif

        yout =  rouse_tab(i1)*(ya-i3*dgy)*(ya-i2*dgy)/((i1-i3)*(i1-i2))
     *         +rouse_tab(i2)*(ya-i1*dgy)*(ya-i3*dgy)/((i2-i1)*(i2-i3))
     *         +rouse_tab(i3)*(ya-i1*dgy)*(ya-i2*dgy)/((i3-i1)*(i3-i2))
        yout = yout/(dgy*dgy)

!        write(6,*)'routab::', i1,i2,i3, yout


        rouse_table = yout

       return
       end

                

       real*8 function SQ_rouse(t,Q,temp,xi,b,epsilon2)
!!     ================================================
!!  from scratch
!!

        implicit none
        real*8 t, Q, temp, eta,kb,Pi,Wl4
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu
        real*8 b,xi

        real*8 SQ_rouse_kernel
        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c
        integer maxit2 

        integer iadda
        common/thiadd/iadda
 

        external SQ_rouse_kernel

 
        epsilon_c = epsilon2
        Bint_end = -dlog(epsilon_c)
        maxit2   = 1300

        Axp = 1.0d0/2.0d0


        t_SI = t * 1d-9
        Q_SI = Q * 1d10

        Gamkt = (kb * temp)*b*b/(12*xi) 
        Wl4   = 36*Gamkt 
        call parset('wl4     ',sngl(wl4*1d31),iadda) ! in A**4ns
        Gamkt = (Gamkt*Q_SI)*Q_SI*(Q_SI*t_SI)*Q_SI
        GamktAxp = Gamkt**Axp

        SQ_rouse = adapint( SQ_rouse_kernel, 0d0, Bint_end,
     *                    epsilon2, maxit2, erraccu)
     
       return
       end

       real*8 function SQ_rouseT(t,Q,temp,xi,b,epsilon2)
!!     =================================================
!!     from table
!!
        implicit none
        real*8 t, Q, temp, eta,kb,Pi, Wl4
        Parameter (kb=1.380662d-23)
        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu
        real*8 b,xi

        real*8 SQ_rouse_kernel
        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c
        integer maxit2 
        real*8 rouse_table

        integer iadda
        common/thiadd/iadda
 
        if(epsilon2.ne.epsilon_c) then
          call create_rouse_table(epsilon2)
        endif


        t_SI = t * 1d-9
        Q_SI = Q * 1d10

        Gamkt = (kb * temp)*b*b/(12*xi)
        Wl4   = 36*Gamkt 
        call parset('wl4     ',sngl(wl4*1d31),iadda) ! in A**4 ns
        Gamkt = (Gamkt*Q_SI)*Q_SI*(Q_SI*t_SI)*Q_SI
        GamktAxp = Gamkt**Axp

        SQ_rouseT = rouse_table(GamktAxp)

     
       return
       end

       real*8 function SQ_rouseW(t,Q,Wl4)
!!     ==================================
!!  from scratch
!!

        implicit none
        double precision, intent(in) ::  t, Q, Wl4
!        Parameter (kb=1.380662d-23)
!        Parameter (Pi=3.141592654d0)
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu
        real*8 b,xi

        real*8 SQ_rouse_kernel
        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c
        integer maxit2 

        integer iadda
        common/thiadd/iadda
 

        external SQ_rouse_kernel

        epsilon2  = 1d-6
        epsilon_c = epsilon2
        Bint_end = -dlog(epsilon_c)
        maxit2   = 1300

        Axp = 1.0d0/2.0d0


  !      t_SI = t * 1d-9
  !      Q_SI = Q * 1d10
  !
  !      Gamkt = (kb * temp)*b*b/(12*xi) 
  !      Wl4   = 36*Gamkt 
  !      call parset('wl4     ',sngl(wl4*1d31),iadda) ! in A**4ns
  !      Gamkt = (Gamkt*Q_SI)*Q_SI*(Q_SI*t_SI)*Q_SI
        Gamkt = (Wl4/36) * Q**4 * t
        GamktAxp = Gamkt**Axp

        SQ_rouseW = adapint( SQ_rouse_kernel, 0d0, Bint_end,
     *                    epsilon2, maxit2, erraccu)
     
       return
       end

       real*8 function SQ_rouseTW(t,Q,Wl4)
!!     ===================================
!!     from table
!!
        implicit none
        double precision, intent(in) ::  t, Q, Wl4
 
        real*8 t_SI, Q_SI, Gamk
        real*8 Gamkt
        real*8 Bint_end
        real*8 adapint
        real*8 epsilon2, erraccu
        real*8 b,xi

        real*8 SQ_rouse_kernel
        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c
        integer maxit2 
        real*8 rouse_table

        integer iadda
        common/thiadd/iadda
 
        epsilon2 = 1d-6
        if(epsilon2.ne.epsilon_c) then
          call create_rouse_table(epsilon2)
        endif

        Axp = 1.0d0/2.0d0


        t_SI = t * 1d-9
        Q_SI = Q * 1d10

!        Gamkt = (kb * temp)*b*b/(12*xi)
!        Wl4   = 36*Gamkt 
!        call parset('wl4     ',sngl(wl4*1d31),iadda) ! in A**4 ns
        Gamkt = (wL4/36)*q**4*t
        GamktAxp = Gamkt**Axp

        SQ_rouseTW = rouse_table(GamktAxp)

     
       return
       end

                
       real*8 function  SQ_rouse_kernel(u)
!!     ----------------------------------

        implicit none
        real*8 u, arg
        real*8 h_Rouse
        real*8 epseff

        real*8 GamktAxp
        common/sqrokc/GamktAxp 
        real*8 axp, epsilon_c
        common /crouse/ axp, epsilon_c

        if(GamktAxp.lt.1d-7) GamktAxp = 1d-7 
        epseff =0.1d0*epsilon_c*dexp(u)/GamktAxp       

        arg = u + GamktAxp * h_Rouse( u / GamktAxp, epseff)

        if(arg.gt.50.0d0) arg = 50.0d0

        SQ_rouse_kernel = dexp(-arg)

        return
        end


 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ROUSE function nach Akcasu Benmouna Gln. 61                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      real*8 function rouse61(t, q, xi, b, temp, epsilon )
c     ----------------------------------------------------
       implicit none
       real*8 kb,Pi
       Parameter (kb=1.380662d-23)
       Parameter (Pi=3.141592654d0)
       real*8   t, q , xi, b, epsilon, epsilon2, temp
       real*8   sum, term0, W, Wt, alpha
       real*8   sum_norm , wint, wint0
       real*8   leps, a1, exp1, exp2
       real*8   t_SI, Q_SI
       real*8   adapint,erraccu 
       real*8   rouse61_kernel
       integer  rs_max, i     
       integer  maxiter
       Parameter (maxiter=300)

       integer m
       common /crouse61/ m
       integer iadda
       common/thiadd/iadda
       external rouse61_kernel
      
       leps = dlog(epsilon)
       if(leps.ge.-1.d0) leps = -1.0d0

       t_SI = t * 1d-9
       Q_SI = q * 1d10
       if(xi.lt.1d-50) xi = 1d-50

       alpha = ((Q_SI*b)**2)/6.0d0
       W     = (3.0d0/(b**2))*kb*temp/xi
       Wt    = 2*W*t_SI
       
       rs_max = -leps/alpha + 1

c       write(6,*)' alpha = ',alpha
c       write(6,*)' W     = ',W
c       write(6,*)' Wt    = ',Wt
c       write(6,*)' t     = ',t
c       write(6,*)' rs_max= ',rs_max      

        call parset('w       ',sngl(w),iadda)

       sum      = 0.0d0
       sum_norm = 0.0d0
       

       do i=1,rs_max
         m   = i
         a1  = alpha*i
         exp1 = dexp(-a1)
         epsilon2 = epsilon/(exp1*alpha)
         erraccu = 0
c       write(6,*)'          m = ',m,i
c       write(6,*)'       exp1 = ',exp1
c       write(6,*)'    epsilon2 = ',epsilon2
         wint=adapint(rouse61_kernel,0.0d0,Wt,epsilon2,maxiter,erraccu)
c       write(6,*)'        wint = ',wint
c       write(6,*)'     erraccu = ',erraccu

         exp2 = dexp(-alpha*wint)
c       write(6,*)'       exp2=', exp2

         sum       = sum + exp1 * exp2
         sum_norm  = sum_norm + exp1
       enddo 
      
       m        = 0
       exp1     = 1.0d0
       epsilon2 = epsilon/(exp1*alpha)
       erraccu  = 0
c       write(6,*)'          m = ',m,i
c       write(6,*)'       exp1 = ',exp1
c       write(6,*)'    epsilon2 = ',epsilon2
       wint=adapint(rouse61_kernel,0.0d0,Wt,epsilon2,maxiter,erraccu)
c       write(6,*)'        wint = ',wint
c       write(6,*)'     erraccu = ',erraccu
       exp2 = dexp(-alpha*wint)

       sum      = 2*sum + exp1*exp2
       sum_norm = 2*sum_norm + exp1

       rouse61 = sum/sum_norm

       return
       end 







       real*8 function rouse61_kernel(u)
c      ---------------------------------
       implicit none
       integer maxord
       parameter(maxord=300)
       real*8   u, xnu, x, bsv(0:maxord)
       real*8   uu
       complex*16 z,bzv(0:maxord)


       parameter(xnu=0.0d0) 
       real*8   DBSIES 
       integer  nz, ifail

       integer m
       common /crouse61/ m

c       if(m.lt.0 .or. m.gt.maxord) then
c         write(6,*)'Error rouse61_kernel Besselfunction order=',m
c         rouse61_kernel = 0.0d0
c         return
c       endif

       if(m.le.100) then
       
!!         call DBSIES(xnu, u, m+1, bsv)
!! nag!!
           ifail = 1
           z     = u
           call s18def(xnu,z,m+1,'S',bzv,nz,ifail)
           bsv(m) = bzv(m)
!!nag!!

         rouse61_kernel = bsv(m)
      
       else
         uu = u * (m/100.0d0)**2
!!         call DBSIES(xnu, uu, 101 , bsv)
!! nag!!
           ifail = 1
           z     = uu
           call s18def(xnu,z,101,'S',bzv,nz,ifail)
           bsv(100) = bzv(100)
!!nag!!
         rouse61_kernel = bsv(100)*(m/100.0d0)
       endif     

       return
       end



      double precision function ddsi (x)
c december 1980 edition, w. fullerton, bell labs.
      double precision x, sics(18), pi2, xsml, absx, f, g, cosx,
     1  dcsevl, d1mach, dcos, dsin, dsqrt
      external d1mach, dcos, dcsevl, dsin, dsqrt, initds
c
c series for si   on the interval  0.00000e+00 to  1.60000e+01
c                                        with weighted error   8.58e-32
c                                         log weighted error  31.07
c                               significant figures required  30.53
c                                    decimal places required  31.69
c
      data si  cs(  1) / -0.1315646598 1848419289 0427517300 0457d0/
      data si  cs(  2) / -0.2776578526 9736018920 4828766015 7299d0/
      data si  cs(  3) /  0.0354414054 8666591797 4913546471 0086d0/
      data si  cs(  4) / -0.0025631631 4479339776 5875278836 1530d0/
      data si  cs(  5) /  0.0001162365 3904970092 8126492148 2985d0/
      data si  cs(  6) / -0.0000035904 3272416060 4267000434 7148d0/
      data si  cs(  7) /  0.0000000802 3421237057 1016230865 2976d0/
      data si  cs(  8) / -0.0000000013 5629976925 4025064993 1846d0/
      data si  cs(  9) /  0.0000000000 1794407215 9973677556 7759d0/
      data si  cs( 10) / -0.0000000000 0019083873 4308714549 0737d0/
      data si  cs( 11) /  0.0000000000 0000166699 8958682433 0853d0/
      data si  cs( 12) / -0.0000000000 0000001217 3098836850 3042d0/
      data si  cs( 13) /  0.0000000000 0000000007 5418186699 3865d0/
      data si  cs( 14) / -0.0000000000 0000000000 0401417884 2446d0/
      data si  cs( 15) /  0.0000000000 0000000000 0001855369 0716d0/
      data si  cs( 16) / -0.0000000000 0000000000 0000007516 6966d0/
      data si  cs( 17) /  0.0000000000 0000000000 0000000026 9113d0/
      data si  cs( 18) / -0.0000000000 0000000000 0000000000 0858d0/
c
      data pi2 / 1.5707963267 9489661923 1321691639 75 d0 /
      data nsi, xsml /0, 0.0d0/
c
      if (nsi.ne.0) go to 10
      nsi = initds (sics, 18, 0.1*sngl(d1mach(3)))
      xsml = sqrt(d1mach(3))
c
 10   absx = abs(x)
      if (absx.gt.4.0d0) go to 20
      ddsi = x
      if (absx.lt.xsml) return
c
      ddsi = x*(0.75d0 + dcsevl ((x*x-8.0d0)*.125d0, sics, nsi))
      return
c
 20   call d9sifg (absx, f, g)
      cosx = cos (absx)
      call erroff
      ddsi = pi2 - f*cosx - g*sin(x)
      if (x.lt.0.0d0) ddsi = -ddsi
c
      return
      end

      subroutine d9sifg (x, f, g)
c december 1980 edition.  w. fullerton, bell labs.
      double precision x, f, g, f1cs(43), f2cs(99), g1cs(44),
     1  g2cs(44), g3cs(56), xbnd, xbndg, xbig, xmaxf, xmaxg,
     2  dcsevl, d1mach, dexp, dlog, dsqrt
      external d1mach, dcsevl, dexp, dlog, dsqrt, initds
c
c series for f1   on the interval  2.00000e-02 to  6.25000e-02
c                                        with weighted error   2.45e-32
c                                         log weighted error  31.61
c                               significant figures required  30.42
c                                    decimal places required  32.43
c
      data f1  cs(  1) / -0.1191081969 0513636103 4820196582 8918d0/
      data f1  cs(  2) / -0.0247823144 9962362475 9007415082 3133d0/
      data f1  cs(  3) /  0.0011910281 4533578212 6812036305 4457d0/
      data f1  cs(  4) / -0.0000927027 7143885617 4830860036 0706d0/
      data f1  cs(  5) /  0.0000093373 1415682709 9686820458 2766d0/
      data f1  cs(  6) / -0.0000011058 2878205571 4393897942 6306d0/
      data f1  cs(  7) /  0.0000001464 7720714601 6216933655 0799d0/
      data f1  cs(  8) / -0.0000000210 6944962876 8953260122 7548d0/
      data f1  cs(  9) /  0.0000000032 2934923668 4823638285 7374d0/
      data f1  cs( 10) / -0.0000000005 2065296175 2937582801 4986d0/
      data f1  cs( 11) /  0.0000000000 8748788845 7027875026 8316d0/
      data f1  cs( 12) / -0.0000000000 1521761870 5612366829 4574d0/
      data f1  cs( 13) /  0.0000000000 0272571924 0541957390 0583d0/
      data f1  cs( 14) / -0.0000000000 0050070530 7596855629 0255d0/
      data f1  cs( 15) /  0.0000000000 0009402409 0272606851 1779d0/
      data f1  cs( 16) / -0.0000000000 0001800144 4479180367 8336d0/
      data f1  cs( 17) /  0.0000000000 0000350626 2143274178 5826d0/
      data f1  cs( 18) / -0.0000000000 0000069352 8292676914 9709d0/
      data f1  cs( 19) /  0.0000000000 0000013909 2513645421 6568d0/
      data f1  cs( 20) / -0.0000000000 0000002824 8688507417 0585d0/
      data f1  cs( 21) /  0.0000000000 0000000580 3130569357 9081d0/
      data f1  cs( 22) / -0.0000000000 0000000120 4690157337 5820d0/
      data f1  cs( 23) /  0.0000000000 0000000025 2505244365 5940d0/
      data f1  cs( 24) / -0.0000000000 0000000005 3398026880 5594d0/
      data f1  cs( 25) /  0.0000000000 0000000001 1385578627 4122d0/
      data f1  cs( 26) / -0.0000000000 0000000000 2446286150 5259d0/
      data f1  cs( 27) /  0.0000000000 0000000000 0529365932 0439d0/
      data f1  cs( 28) / -0.0000000000 0000000000 0115318494 0277d0/
      data f1  cs( 29) /  0.0000000000 0000000000 0025278656 8318d0/
      data f1  cs( 30) / -0.0000000000 0000000000 0005573864 5378d0/
      data f1  cs( 31) /  0.0000000000 0000000000 0001235824 5621d0/
      data f1  cs( 32) / -0.0000000000 0000000000 0000275435 0842d0/
      data f1  cs( 33) /  0.0000000000 0000000000 0000061690 6808d0/
      data f1  cs( 34) / -0.0000000000 0000000000 0000013881 7443d0/
      data f1  cs( 35) /  0.0000000000 0000000000 0000003137 5329d0/
      data f1  cs( 36) / -0.0000000000 0000000000 0000000712 1249d0/
      data f1  cs( 37) /  0.0000000000 0000000000 0000000162 2778d0/
      data f1  cs( 38) / -0.0000000000 0000000000 0000000037 1206d0/
      data f1  cs( 39) /  0.0000000000 0000000000 0000000008 5221d0/
      data f1  cs( 40) / -0.0000000000 0000000000 0000000001 9633d0/
      data f1  cs( 41) /  0.0000000000 0000000000 0000000000 4538d0/
      data f1  cs( 42) / -0.0000000000 0000000000 0000000000 1052d0/
      data f1  cs( 43) /  0.0000000000 0000000000 0000000000 0245d0/
c
c series for f2   on the interval  0.00000e+00 to  2.00000e-02
c                                        with weighted error   2.38e-32
c                                         log weighted error  31.62
c                               significant figures required  30.01
c                                    decimal places required  32.62
c
      data f2  cs(  1) / -0.0348409253 8970132330 8360497337 45577d0/
      data f2  cs(  2) / -0.0166842205 6779596873 2467863122 78676d0/
      data f2  cs(  3) /  0.0006752901 2412377385 0452078592 39727d0/
      data f2  cs(  4) / -0.0000535066 6225447013 6287855775 57429d0/
      data f2  cs(  5) /  0.0000062693 4217790075 2670507594 31626d0/
      data f2  cs(  6) / -0.0000009526 6388019916 6806777904 14293d0/
      data f2  cs(  7) /  0.0000001745 6292242509 8804255044 27666d0/
      data f2  cs(  8) / -0.0000000368 7954030653 0933070976 46628d0/
      data f2  cs(  9) /  0.0000000087 2026777051 3952640758 16938d0/
      data f2  cs( 10) / -0.0000000022 6019703919 7387485304 23167d0/
      data f2  cs( 11) /  0.0000000006 3246249765 2506125204 44877d0/
      data f2  cs( 12) / -0.0000000001 8889118884 7178692409 11480d0/
      data f2  cs( 13) /  0.0000000000 5967746729 9978133726 20472d0/
      data f2  cs( 14) / -0.0000000000 1980443117 3722390111 96007d0/
      data f2  cs( 15) /  0.0000000000 0686413954 7721033837 13264d0/
      data f2  cs( 16) / -0.0000000000 0247310193 0701991060 74890d0/
      data f2  cs( 17) /  0.0000000000 0092263594 5499414041 96042d0/
      data f2  cs( 18) / -0.0000000000 0035523634 9992617844 97297d0/
      data f2  cs( 19) /  0.0000000000 0014076049 6253515914 61820d0/
      data f2  cs( 20) / -0.0000000000 0005726228 4997476527 94311d0/
      data f2  cs( 21) /  0.0000000000 0002386537 5454131718 10106d0/
      data f2  cs( 22) / -0.0000000000 0001017141 8907645971 42232d0/
      data f2  cs( 23) /  0.0000000000 0000442594 5310783644 24968d0/
      data f2  cs( 24) / -0.0000000000 0000196344 9330491897 61979d0/
      data f2  cs( 25) /  0.0000000000 0000088688 7483148104 61024d0/
      data f2  cs( 26) / -0.0000000000 0000040743 3450273115 46948d0/
      data f2  cs( 27) /  0.0000000000 0000019016 8372156753 39859d0/
      data f2  cs( 28) / -0.0000000000 0000009009 7072974780 42442d0/
      data f2  cs( 29) /  0.0000000000 0000004329 2112740956 68667d0/
      data f2  cs( 30) / -0.0000000000 0000002108 1444653224 79526d0/
      data f2  cs( 31) /  0.0000000000 0000001039 6379070264 52274d0/
      data f2  cs( 32) / -0.0000000000 0000000518 8910079489 31936d0/
      data f2  cs( 33) /  0.0000000000 0000000261 9553248698 99371d0/
      data f2  cs( 34) / -0.0000000000 0000000133 6903999513 01570d0/
      data f2  cs( 35) /  0.0000000000 0000000068 9410577029 31664d0/
      data f2  cs( 36) / -0.0000000000 0000000035 9053626104 37250d0/
      data f2  cs( 37) /  0.0000000000 0000000018 8780772557 91706d0/
      data f2  cs( 38) / -0.0000000000 0000000010 0161252655 94380d0/
      data f2  cs( 39) /  0.0000000000 0000000005 3607256915 78228d0/
      data f2  cs( 40) / -0.0000000000 0000000002 8931989749 44827d0/
      data f2  cs( 41) /  0.0000000000 0000000001 5740651002 02625d0/
      data f2  cs( 42) / -0.0000000000 0000000000 8630271064 31206d0/
      data f2  cs( 43) /  0.0000000000 0000000000 4767156028 62288d0/
      data f2  cs( 44) / -0.0000000000 0000000000 2652227399 98504d0/
      data f2  cs( 45) /  0.0000000000 0000000000 1485828650 63866d0/
      data f2  cs( 46) / -0.0000000000 0000000000 0837972359 23135d0/
      data f2  cs( 47) /  0.0000000000 0000000000 0475659164 22711d0/
      data f2  cs( 48) / -0.0000000000 0000000000 0271690733 53112d0/
      data f2  cs( 49) /  0.0000000000 0000000000 0156127388 81686d0/
      data f2  cs( 50) / -0.0000000000 0000000000 0090245550 78347d0/
      data f2  cs( 51) /  0.0000000000 0000000000 0052460970 49119d0/
      data f2  cs( 52) / -0.0000000000 0000000000 0030664508 18697d0/
      data f2  cs( 53) /  0.0000000000 0000000000 0018019962 50957d0/
      data f2  cs( 54) / -0.0000000000 0000000000 0010644430 50752d0/
      data f2  cs( 55) /  0.0000000000 0000000000 0006319421 58881d0/
      data f2  cs( 56) / -0.0000000000 0000000000 0003770138 12246d0/
      data f2  cs( 57) /  0.0000000000 0000000000 0002259975 42918d0/
      data f2  cs( 58) / -0.0000000000 0000000000 0001361008 44814d0/
      data f2  cs( 59) /  0.0000000000 0000000000 0000823332 32003d0/
      data f2  cs( 60) / -0.0000000000 0000000000 0000500259 86091d0/
      data f2  cs( 61) /  0.0000000000 0000000000 0000305262 45684d0/
      data f2  cs( 62) / -0.0000000000 0000000000 0000187051 64021d0/
      data f2  cs( 63) /  0.0000000000 0000000000 0000115084 04393d0/
      data f2  cs( 64) / -0.0000000000 0000000000 0000071087 14611d0/
      data f2  cs( 65) /  0.0000000000 0000000000 0000044080 65533d0/
      data f2  cs( 66) / -0.0000000000 0000000000 0000027437 60867d0/
      data f2  cs( 67) /  0.0000000000 0000000000 0000017141 44851d0/
      data f2  cs( 68) / -0.0000000000 0000000000 0000010747 68860d0/
      data f2  cs( 69) /  0.0000000000 0000000000 0000006762 59777d0/
      data f2  cs( 70) / -0.0000000000 0000000000 0000004269 81348d0/
      data f2  cs( 71) /  0.0000000000 0000000000 0000002705 00637d0/
      data f2  cs( 72) / -0.0000000000 0000000000 0000001719 33331d0/
      data f2  cs( 73) /  0.0000000000 0000000000 0000001096 36138d0/
      data f2  cs( 74) / -0.0000000000 0000000000 0000000701 32573d0/
      data f2  cs( 75) /  0.0000000000 0000000000 0000000450 01784d0/
      data f2  cs( 76) / -0.0000000000 0000000000 0000000289 63835d0/
      data f2  cs( 77) /  0.0000000000 0000000000 0000000186 97009d0/
      data f2  cs( 78) / -0.0000000000 0000000000 0000000121 04646d0/
      data f2  cs( 79) /  0.0000000000 0000000000 0000000078 59065d0/
      data f2  cs( 80) / -0.0000000000 0000000000 0000000051 16867d0/
      data f2  cs( 81) /  0.0000000000 0000000000 0000000033 40627d0/
      data f2  cs( 82) / -0.0000000000 0000000000 0000000021 86851d0/
      data f2  cs( 83) /  0.0000000000 0000000000 0000000014 35340d0/
      data f2  cs( 84) / -0.0000000000 0000000000 0000000009 44523d0/
      data f2  cs( 85) /  0.0000000000 0000000000 0000000006 23117d0/
      data f2  cs( 86) / -0.0000000000 0000000000 0000000004 12101d0/
      data f2  cs( 87) /  0.0000000000 0000000000 0000000002 73208d0/
      data f2  cs( 88) / -0.0000000000 0000000000 0000000001 81558d0/
      data f2  cs( 89) /  0.0000000000 0000000000 0000000001 20934d0/
      data f2  cs( 90) / -0.0000000000 0000000000 0000000000 80737d0/
      data f2  cs( 91) /  0.0000000000 0000000000 0000000000 54022d0/
      data f2  cs( 92) / -0.0000000000 0000000000 0000000000 36227d0/
      data f2  cs( 93) /  0.0000000000 0000000000 0000000000 24348d0/
      data f2  cs( 94) / -0.0000000000 0000000000 0000000000 16401d0/
      data f2  cs( 95) /  0.0000000000 0000000000 0000000000 11074d0/
      data f2  cs( 96) / -0.0000000000 0000000000 0000000000 07497d0/
      data f2  cs( 97) /  0.0000000000 0000000000 0000000000 05091d0/
      data f2  cs( 98) / -0.0000000000 0000000000 0000000000 03470d0/
      data f2  cs( 99) /  0.0000000000 0000000000 0000000000 02377d0/
c
c series for g1   on the interval  2.00000e-02 to  6.25000e-02
c                                        with weighted error   7.23e-32
c                                         log weighted error  31.14
c                               significant figures required  30.35
c                                    decimal places required  31.96
c
      data g1  cs(  1) / -0.3040578798 2534959544 9972668209 1083d0/
      data g1  cs(  2) / -0.0566890984 5971205877 3133915611 8269d0/
      data g1  cs(  3) /  0.0039046158 1732756439 1998407155 4082d0/
      data g1  cs(  4) / -0.0003746075 9592022606 1861933986 7489d0/
      data g1  cs(  5) /  0.0000435431 5565598436 7955222084 0065d0/
      data g1  cs(  6) / -0.0000057417 2944530250 4656197072 3475d0/
      data g1  cs(  7) /  0.0000008282 5521045026 2974193761 6492d0/
      data g1  cs(  8) / -0.0000001278 2458925946 4272788391 3223d0/
      data g1  cs(  9) /  0.0000000207 9783529486 8788443925 7529d0/
      data g1  cs( 10) / -0.0000000035 3132059219 9079804203 2682d0/
      data g1  cs( 11) /  0.0000000006 2108242363 0895106863 1449d0/
      data g1  cs( 12) / -0.0000000001 1252154744 4629264933 6987d0/
      data g1  cs( 13) /  0.0000000000 2090889176 8442160526 7019d0/
      data g1  cs( 14) / -0.0000000000 0397158317 3768172768 9158d0/
      data g1  cs( 15) /  0.0000000000 0076904313 1427208993 9005d0/
      data g1  cs( 16) / -0.0000000000 0015146967 4273161351 9826d0/
      data g1  cs( 17) /  0.0000000000 0003028921 4655235968 4119d0/
      data g1  cs( 18) / -0.0000000000 0000613997 0383470882 5400d0/
      data g1  cs( 19) /  0.0000000000 0000126006 0582951093 3553d0/
      data g1  cs( 20) / -0.0000000000 0000026150 2925093948 3683d0/
      data g1  cs( 21) /  0.0000000000 0000005482 7884489179 6821d0/
      data g1  cs( 22) / -0.0000000000 0000001160 3818212952 6571d0/
      data g1  cs( 23) /  0.0000000000 0000000247 7165410712 9795d0/
      data g1  cs( 24) / -0.0000000000 0000000053 3067275322 3389d0/
      data g1  cs( 25) /  0.0000000000 0000000011 5566607559 8465d0/
      data g1  cs( 26) / -0.0000000000 0000000002 5228054774 4957d0/
      data g1  cs( 27) /  0.0000000000 0000000000 5542903855 0786d0/
      data g1  cs( 28) / -0.0000000000 0000000000 1225220842 1297d0/
      data g1  cs( 29) /  0.0000000000 0000000000 0272366431 8684d0/
      data g1  cs( 30) / -0.0000000000 0000000000 0060870783 1422d0/
      data g1  cs( 31) /  0.0000000000 0000000000 0013672487 4476d0/
      data g1  cs( 32) / -0.0000000000 0000000000 0003085662 6806d0/
      data g1  cs( 33) /  0.0000000000 0000000000 0000699521 2319d0/
      data g1  cs( 34) / -0.0000000000 0000000000 0000159258 7569d0/
      data g1  cs( 35) /  0.0000000000 0000000000 0000036405 1056d0/
      data g1  cs( 36) / -0.0000000000 0000000000 0000008353 9465d0/
      data g1  cs( 37) /  0.0000000000 0000000000 0000001924 0303d0/
      data g1  cs( 38) / -0.0000000000 0000000000 0000000444 6816d0/
      data g1  cs( 39) /  0.0000000000 0000000000 0000000103 1182d0/
      data g1  cs( 40) / -0.0000000000 0000000000 0000000023 9887d0/
      data g1  cs( 41) /  0.0000000000 0000000000 0000000005 5976d0/
      data g1  cs( 42) / -0.0000000000 0000000000 0000000001 3100d0/
      data g1  cs( 43) /  0.0000000000 0000000000 0000000000 3074d0/
      data g1  cs( 44) / -0.0000000000 0000000000 0000000000 0723d0/
c
c series for g2   on the interval  5.00000e-03 to  2.00000e-02
c                                        with weighted error   3.25e-32
c                                         log weighted error  31.49
c                               significant figures required  30.32
c                                    decimal places required  32.31
c
      data g2  cs(  1) / -0.1211802894 7316462635 4183404685 8267d0/
      data g2  cs(  2) / -0.0316761386 3949502867 0140792350 5610d0/
      data g2  cs(  3) /  0.0013383199 7788626801 6381942949 2182d0/
      data g2  cs(  4) / -0.0000895511 0113922524 2553190506 9518d0/
      data g2  cs(  5) /  0.0000079155 5629617182 1311524946 7924d0/
      data g2  cs(  6) / -0.0000008438 7933222415 2018141898 2080d0/
      data g2  cs(  7) /  0.0000001029 9804256775 3014664722 7274d0/
      data g2  cs(  8) / -0.0000000139 2957506051 8383579583 4444d0/
      data g2  cs(  9) /  0.0000000020 4227039598 7598040067 7594d0/
      data g2  cs( 10) / -0.0000000003 1965346942 0642703543 4752d0/
      data g2  cs( 11) /  0.0000000000 5281478326 5726769861 5312d0/
      data g2  cs( 12) / -0.0000000000 0913395546 7267103373 5289d0/
      data g2  cs( 13) /  0.0000000000 0164262512 3896776044 4819d0/
      data g2  cs( 14) / -0.0000000000 0030558970 3932266000 2410d0/
      data g2  cs( 15) /  0.0000000000 0005856558 2578577971 7892d0/
      data g2  cs( 16) / -0.0000000000 0001152291 9773094012 0563d0/
      data g2  cs( 17) /  0.0000000000 0000232094 6911998853 7310d0/
      data g2  cs( 18) / -0.0000000000 0000047743 5583417753 5025d0/
      data g2  cs( 19) /  0.0000000000 0000010009 9676580018 0573d0/
      data g2  cs( 20) / -0.0000000000 0000002135 3377808225 6704d0/
      data g2  cs( 21) /  0.0000000000 0000000462 7719077736 7671d0/
      data g2  cs( 22) / -0.0000000000 0000000101 7580741022 7657d0/
      data g2  cs( 23) /  0.0000000000 0000000022 6765739988 4672d0/
      data g2  cs( 24) / -0.0000000000 0000000005 1163077607 6426d0/
      data g2  cs( 25) /  0.0000000000 0000000001 1676701491 3108d0/
      data g2  cs( 26) / -0.0000000000 0000000000 2693542767 2470d0/
      data g2  cs( 27) /  0.0000000000 0000000000 0627566584 1146d0/
      data g2  cs( 28) / -0.0000000000 0000000000 0147588055 7531d0/
      data g2  cs( 29) /  0.0000000000 0000000000 0035014531 4739d0/
      data g2  cs( 30) / -0.0000000000 0000000000 0008375773 2152d0/
      data g2  cs( 31) /  0.0000000000 0000000000 0002019181 5152d0/
      data g2  cs( 32) / -0.0000000000 0000000000 0000490356 7705d0/
      data g2  cs( 33) /  0.0000000000 0000000000 0000119912 3348d0/
      data g2  cs( 34) / -0.0000000000 0000000000 0000029517 0610d0/
      data g2  cs( 35) /  0.0000000000 0000000000 0000007311 3112d0/
      data g2  cs( 36) / -0.0000000000 0000000000 0000001821 7843d0/
      data g2  cs( 37) /  0.0000000000 0000000000 0000000456 5148d0/
      data g2  cs( 38) / -0.0000000000 0000000000 0000000115 0151d0/
      data g2  cs( 39) /  0.0000000000 0000000000 0000000029 1267d0/
      data g2  cs( 40) / -0.0000000000 0000000000 0000000007 4125d0/
      data g2  cs( 41) /  0.0000000000 0000000000 0000000001 8953d0/
      data g2  cs( 42) / -0.0000000000 0000000000 0000000000 4868d0/
      data g2  cs( 43) /  0.0000000000 0000000000 0000000000 1256d0/
      data g2  cs( 44) / -0.0000000000 0000000000 0000000000 0325d0/
c
c series for g3   on the interval  0.00000e+00 to  5.00000e-03
c                                        with weighted error   3.83e-32
c                                         log weighted error  31.42
c                               significant figures required  29.71
c                                    decimal places required  32.29
c
      data g3  cs(  1) / -0.0280574367 8094729284 0281526433 5299d0/
      data g3  cs(  2) / -0.0137271597 1622369754 0910050808 9556d0/
      data g3  cs(  3) /  0.0002894032 6387602960 2744894127 3751d0/
      data g3  cs(  4) / -0.0000114129 2393911971 4590874362 2517d0/
      data g3  cs(  5) /  0.0000006813 9655907262 4299772020 7302d0/
      data g3  cs(  6) / -0.0000000547 9522896046 5236366905 8052d0/
      data g3  cs(  7) /  0.0000000055 2074299182 1252910940 6521d0/
      data g3  cs(  8) / -0.0000000006 6414641993 2292002249 1428d0/
      data g3  cs(  9) /  0.0000000000 9223736634 8704110856 4960d0/
      data g3  cs( 10) / -0.0000000000 1442990888 8668286261 1718d0/
      data g3  cs( 11) /  0.0000000000 0249639048 9203071024 8705d0/
      data g3  cs( 12) / -0.0000000000 0047082406 7587524472 2971d0/
      data g3  cs( 13) /  0.0000000000 0009572176 5921675998 8140d0/
      data g3  cs( 14) / -0.0000000000 0002078899 6609580903 0537d0/
      data g3  cs( 15) /  0.0000000000 0000478750 9997087743 1627d0/
      data g3  cs( 16) / -0.0000000000 0000116190 7058337717 3759d0/
      data g3  cs( 17) /  0.0000000000 0000029565 0896926783 6974d0/
      data g3  cs( 18) / -0.0000000000 0000007852 9498825649 2025d0/
      data g3  cs( 19) /  0.0000000000 0000002169 2226436825 6612d0/
      data g3  cs( 20) / -0.0000000000 0000000621 1351583167 6342d0/
      data g3  cs( 21) /  0.0000000000 0000000183 8456883845 0977d0/
      data g3  cs( 22) / -0.0000000000 0000000056 1088748213 7276d0/
      data g3  cs( 23) /  0.0000000000 0000000017 6186280528 0062d0/
      data g3  cs( 24) / -0.0000000000 0000000005 6811105054 1451d0/
      data g3  cs( 25) /  0.0000000000 0000000001 8778627958 2313d0/
      data g3  cs( 26) / -0.0000000000 0000000000 6353169415 1124d0/
      data g3  cs( 27) /  0.0000000000 0000000000 2196880236 8238d0/
      data g3  cs( 28) / -0.0000000000 0000000000 0775466655 0395d0/
      data g3  cs( 29) /  0.0000000000 0000000000 0279101835 6581d0/
      data g3  cs( 30) / -0.0000000000 0000000000 0102317852 5247d0/
      data g3  cs( 31) /  0.0000000000 0000000000 0038169340 3919d0/
      data g3  cs( 32) / -0.0000000000 0000000000 0014476789 5606d0/
      data g3  cs( 33) /  0.0000000000 0000000000 0005577951 2634d0/
      data g3  cs( 34) / -0.0000000000 0000000000 0002181723 9071d0/
      data g3  cs( 35) /  0.0000000000 0000000000 0000865664 6309d0/
      data g3  cs( 36) / -0.0000000000 0000000000 0000348215 7895d0/
      data g3  cs( 37) /  0.0000000000 0000000000 0000141918 8130d0/
      data g3  cs( 38) / -0.0000000000 0000000000 0000058571 4314d0/
      data g3  cs( 39) /  0.0000000000 0000000000 0000024466 0482d0/
      data g3  cs( 40) / -0.0000000000 0000000000 0000010338 7099d0/
      data g3  cs( 41) /  0.0000000000 0000000000 0000004417 7299d0/
      data g3  cs( 42) / -0.0000000000 0000000000 0000001908 0079d0/
      data g3  cs( 43) /  0.0000000000 0000000000 0000000832 6038d0/
      data g3  cs( 44) / -0.0000000000 0000000000 0000000366 9553d0/
      data g3  cs( 45) /  0.0000000000 0000000000 0000000163 2875d0/
      data g3  cs( 46) / -0.0000000000 0000000000 0000000073 3357d0/
      data g3  cs( 47) /  0.0000000000 0000000000 0000000033 2327d0/
      data g3  cs( 48) / -0.0000000000 0000000000 0000000015 1906d0/
      data g3  cs( 49) /  0.0000000000 0000000000 0000000007 0020d0/
      data g3  cs( 50) / -0.0000000000 0000000000 0000000003 2539d0/
      data g3  cs( 51) /  0.0000000000 0000000000 0000000001 5240d0/
      data g3  cs( 52) / -0.0000000000 0000000000 0000000000 7193d0/
      data g3  cs( 53) /  0.0000000000 0000000000 0000000000 3420d0/
      data g3  cs( 54) / -0.0000000000 0000000000 0000000000 1638d0/
      data g3  cs( 55) /  0.0000000000 0000000000 0000000000 0790d0/
      data g3  cs( 56) / -0.0000000000 0000000000 0000000000 0383d0/
c
      data nf1, nf2, ng1, ng2, ng3 / 5*0 /
      data xbnd, xbndg, xbig, xmaxf, xmaxg / 5*0.0d0 /
c
      if (nf1.ne.0) go to 10
      tol = 0.1d0*d1mach(3)
      nf1 = initds (f1cs, 43, tol)
      nf2 = initds (f2cs, 99, tol)
      ng1 = initds (g1cs, 44, tol)
      ng2 = initds (g2cs, 44, tol)
      ng3 = initds (g3cs, 56, tol)
c
      xbig = sqrt(1.0d0/d1mach(3))
      xmaxf = exp (dmin1 (-log(d1mach(1)), log(d1mach(2))) - .01d0)
      xmaxg = 1.0d0/sqrt(d1mach(1))
      xbnd = sqrt (50.0d0)
      xbndg = sqrt (200.d0)
c
 10   if (x.lt.4.0d0) write(*,*)" 34hd9sifg  pproxs invalid for x lt 4"
c
      if (x.gt.xbnd) go to 20
      f = (1.0d0 + dcsevl ((1.d0/x**2-0.04125d0)/.02125d0, f1cs, nf1))/x
      g = (1.0d0 + dcsevl((1.d0/x**2-.04125d0)/.02125d0, g1cs,ng1))/x**2
      return
c
 20   if (x.gt.xbig) go to 30
      f = (1.0d0 + dcsevl (100.d0/x**2-1.d0, f2cs, nf2))/x
      if (x.le.xbndg) g = (1.0d0 + dcsevl ((10000.d0/x**2-125.d0)/75.d0,
     1  g2cs, ng2))/x**2
      if (x.gt.xbndg) g = (1.0d0 + dcsevl (400.d0/x**2-1.d0, g3cs,
     1  ng3))/x**2
      return
c
 30   f = 0.d0
      if (x.lt.xmaxf) f = 1.0d0/x
      g = 0.d0
      if (x.lt.xmaxg) g = 1.0d0/x**2
      return
c
      end

      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end

      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
c       call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
c        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   continue
c     call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end

      SUBROUTINE ERROFF
C
C  TURNS OFF THE ERROR STATE OFF BY SETTING LERROR=0.
C
      I=I8SAVE(1,0,.TRUE.)
      RETURN
C
      END



      SUBROUTINE E9RINT(MESSG,NW,NERR,SAVE)
C
C  THIS ROUTINE STORES THE CURRENT ERROR MESSAGE OR PRINTS THE OLD ONE,
C  IF ANY, DEPENDING ON WHETHER OR NOT SAVE = .TRUE. .
C
C  CHANGED, BY P.FOX, MAY 18, 1983, FROM THE ORIGINAL VERSION IN ORDER
C  TO GET RID OF THE FORTRAN CARRIAGE CONTROL LINE OVERWRITE
C  CHARACTER +, WHICH HAS ALWAYS CAUSED TROUBLE.
C  FOR THE RECORD, THE PREVIOUS VERSION HAD THE FOLLOWING ARRAY
C  AND CALLS -   (WHERE CCPLUS WAS DECLARED OF TYPE INTEGER)
C
C      DATA CCPLUS  / 1H+ /
C
C      DATA FMT( 1) / 1H( /
C      DATA FMT( 2) / 1HA /
C      DATA FMT( 3) / 1H1 /
C      DATA FMT( 4) / 1H, /
C      DATA FMT( 5) / 1H1 /
C      DATA FMT( 6) / 1H4 /
C      DATA FMT( 7) / 1HX /
C      DATA FMT( 8) / 1H, /
C      DATA FMT( 9) / 1H7 /
C      DATA FMT(10) / 1H2 /
C      DATA FMT(11) / 1HA /
C      DATA FMT(12) / 1HX /
C      DATA FMT(13) / 1HX /
C      DATA FMT(14) / 1H) /
C
C        CALL S88FMT(2,I1MACH(6),FMT(12))
C        WRITE(IWUNIT,FMT) CCPLUS,(MESSGP(I),I=1,NWP)
C
       CHARACTER*1 MESSG(NW)

      LOGICAL SAVE
C
C  MESSGP STORES AT LEAST THE FIRST 72 CHARACTERS OF THE PREVIOUS
C  MESSAGE. ITS LENGTH IS MACHINE DEPENDENT AND MUST BE AT LEAST
C
C       1 + 71/(THE NUMBER OF CHARACTERS STORED PER INTEGER WORD).
C
      CHARACTER*1 MESSGP(72),FMT(10)
      CHARACTER*10 FMT10
      EQUIVALENCE (FMT(1),FMT10)
C
C  START WITH NO PREVIOUS MESSAGE.
C

      DATA MESSGP(1)/'1'/, NWP/0/, NERRP/0/

C
C  SET UP THE FORMAT FOR PRINTING THE ERROR MESSAGE.
C  THE FORMAT IS SIMPLY (A1,14X,72AXX) WHERE XX=I1MACH(6) IS THE
C  NUMBER OF CHARACTERS STORED PER INTEGER WORD.
C
      DATA FMT( 1) / '(' /
      DATA FMT( 2) / '3' /
      DATA FMT( 3) / 'X' /
      DATA FMT( 4) / ',' /
      DATA FMT( 5) / '7' /
      DATA FMT( 6) / '2' /
      DATA FMT( 7) / 'A' /
      DATA FMT( 8) / 'X' /
      DATA FMT( 9) / 'X' /
      DATA FMT(10) / ')' /

C
      IF (.NOT.SAVE) GO TO 20
C
C  SAVE THE MESSAGE.
C
        NWP=NW
        NERRP=NERR
        DO 10 I=1,NW
 10     MESSGP(I)=MESSG(I)
C
        GO TO 30
C
 20   IF (I8SAVE(1,0,.FALSE.).EQ.0) GO TO 30
C
C  PRINT THE MESSAGE.
C
        IWUNIT=I1MACH(4)
        WRITE(IWUNIT,9000) NERRP
 9000   FORMAT(7H ERROR ,I4,4H IN )
c        CALL S88FMT(2, 1, FMT(8))

c        WRITE(IWUNIT,FMT10) (MESSGP(I),I=1,NWP)
C
 30   RETURN
C
      END


      SUBROUTINE EPRINT
C
C  THIS SUBROUTINE PRINTS THE LAST ERROR MESSAGE, IF ANY.
C
      CHARACTER*1 MESSG(1)

C
      CALL E9RINT(MESSG,1,1,.FALSE.)
      RETURN
C
      END

      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
