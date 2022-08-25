!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  ZIMM function subroutines                                           !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       real*8 function h_zimm(u,epsilon)

         implicit none
         real*8  Dsi, u, Aintsep, epsilon, erraccu, Pi
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
         Aintsep = (-lgeps*dsqrt(2.0d0))**axp
 
         if(u.lt.1d-7) u = 1d-7
         uc = u

         x = Aintsep * u
         Part2 = (dcos(x)+DSi(x)*x)/Aintsep - Pi*u/2

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

          zimm_tab(i) = adapint( SQ_zimm_kernel, 0, Bint_end,
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

        SQ_zimm = adapint( SQ_zimm_kernel, 0, Bint_end,
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
         real*8  Dsi, u, Aintsep, epsilon, erraccu, Pi
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
         Part2 = (dcos(x)+DSi(x)*x)/Aintsep - Pi*u/2

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
        Parameter(nroutab=200) 
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

          rouse_tab(i) = adapint( SQ_rouse_kernel, 0, Bint_end,
     *                      epsilon2, maxit2, erraccu)
          
          write(6,*)i,GamktAxp,rouse_tab(i)
        enddo
     
       return
       end


       real*8 function rouse_table(y)
!!     =====================================

        implicit none

        integer nroutab
        Parameter(nroutab=200) 
        real*8 dgy, rouse_tab(nroutab)
        common/croutab/dgy, rouse_tab 

        real*8 ya,y,h, yout,arg
        integer i1,i2,i3
     
        ya = dabs(y)
        i1 = ya/dgy + 1
        i2 = i1+1
        i3 = i1+2

        if(i2.ge.nroutab) then
          arg = -1.77*ya
          if(arg.lt.-50.0d0) arg = -50.0d0
          rouse_table = dexp(arg)
          return
        endif

        

        yout =  rouse_tab(i1)*(ya-i3*dgy)*(ya-i2*dgy)/((i1-i3)*(i1-i2))
     *         +rouse_tab(i2)*(ya-i1*dgy)*(ya-i3*dgy)/((i2-i1)*(i2-i3))
     *         +rouse_tab(i3)*(ya-i1*dgy)*(ya-i2*dgy)/((i3-i1)*(i3-i2))
        yout = yout/(dgy*dgy)

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

        SQ_rouse = adapint( SQ_rouse_kernel, 0, Bint_end,
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
