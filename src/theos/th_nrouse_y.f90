      FUNCTION th_nrouse_y (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> nrouse_y <--------                                              
!                                                                       
!                
      implicit none  
      double precision, parameter :: Pi                  =     3.141592653589793238462643d0    
                                                 
      CHARACTER(len=8), intent(inout) :: thnam, parnam (20) 
      real, intent(in)                :: x,  pa (20)
      integer, intent(in)             :: ini, mbuf
      integer, intent(inout)          :: npar
      integer, intent(inout)          :: nopar            ! Anzahl der Parameter data
      character*80, intent(inout)     :: napar(mbuf)      ! name des parameters n
      real, intent(inout)             :: params(mbuf)     ! value des parameters n
      real                            :: th_nrouse_y

      integer          :: ier, nparx

      double precision :: temp, tau, q, wl4, l 
      double precision :: a0 , qz, R, dr
      double precision :: p_trans, w_scale, a_scale
      double precision :: diff 
      double precision :: sq, sqt
      REAL             :: qget 
      integer          :: n
      logical          :: sqof0

      integer          :: iadda 
      common/thiadd/iadda

      double precision :: a,b,c,d,beta0,betae,Z 



 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'nrouse_y'
         nparx = 14
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_nrouse_y = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '
         parnam(3) = 'n_segmnt'   ! : if <0 return non-normalized s(q,t)
         parnam(4) = 're      '   ! : re**2 = N * b**2
         parnam(5) = 'com_diff'
         parnam(6) = 'p_trans '   ! transition mode number
         parnam(7) = 'w_scale '   ! rate scale for low modes
         parnam(8) = 'a_scale '   ! amplitude scale for low modes
         parnam(9) = 'bet_a   '   ! beta model a parameter
         parnam(10)= 'bet_c   '   ! beta model c parameter
         parnam(11)= 'bet_d   '   ! beta model d parameter
         parnam(12)= 'beta0   '   ! beta model beta0 parameter
         parnam(13)= 'betae   '   ! beta model min beta parameter
         parnam(14)= 'z       '   ! beta model Z parameter
 
         th_nrouse_y = 0
         return
       endif
 
       sqof0 = .false.
!  ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       Wl4      = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       diff     = abs(pa(5))  ! in cm**2/sec
       p_trans  = pa(6)
       w_scale  = abs(pa(7))
       a_scale  = abs(pa(8))
       a        = abs(pa(9))
       c        = pa(10)
       d        = abs(pa(11))
       beta0    = abs(pa(12))
       betae    = abs(pa(13))
       z        = abs(pa(14))

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) then 
         write(6,*)'Assume Q as independend variable!'
         qz   = x
         tau = 0.0d0
         call parget('tau     ',qget,iadda,ier)
         if(ier.eq.0) then
          tau = qget
         endif
         write(6,*)'tau = ',tau 
         sqof0 = .true.
       endif

       dr = diff
             
       call Nrouse3(qz,tau,dr,Wl4,N,R, p_trans, w_scale, a_scale,a,b,c,d,beta0,betae,Z, l,Sq,Sqt)
 
       if(sqof0) then
          th_nrouse_y =  a0 * Sqt
       else
          th_nrouse_y =  a0 * Sqt/Sq
       endif

 


       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
 
       return
       end




       subroutine Nrouse3(q,t,Dr,Wl4,N,R,p_trans, w_scale, a_scale, a,b,c,d, beta0,betae,Z, l,Sq,Sqt)
!      ==============================================================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    Wl4   ----> Rouse rate in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
! p_trans  ----> mode number for transition from scaled to pure Rouse
! w_trans  ----> rate scale
! a_trans  ----> amplitiude scale
! a,b,c,beta0,betae,Z --> mode dependent beta model params
! Output parameters:
!    l     <--- segment length
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in)     ::  q,t,Wl4,R
       double precision, intent(in)     ::  p_trans, w_scale, a_scale
       double precision, intent(in)     ::  a, c, d, beta0, betae, z
       double precision, intent(inout)  ::  Dr
       double precision, intent(out)    ::  l, Sq, Sqt
       integer,          intent(in)     ::  N

       integer                          ::  nn,mm,p

       double precision :: W, rate_p, kbt, Sq0, arg1, arg2
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2
       double precision, save :: b3old

       double precision :: wscal(n), ascal(n), ppsin(n), cospn(n,n)
       double precision :: beta(n),  tteff(n)
       integer          :: pv(n)

       double precision :: b, b0
       integer          :: pmax

       integer          :: iout, itrans
        
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

       itrans = nint(p_trans)
       wscal  = 1
       ascal  = 1
       if(itrans > 0) then
         wscal(1:min(n,itrans)) = w_scale
         ascal(1:min(n,itrans)) = a_scale
       endif


! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----

       W   = Wl4 / l**4

! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = Wl4 / ( 3*N *l**2) 
       endif

! to increase efficiency prepare p**2  and the matrix of cos factors 
       ppsin(1:N) = 4 * (N/Pi)**2 * sin( (pi/(2*N)) * [(p,p=1,N)] )**2
       do nn =1,N
          do p=1,N
            cospn(p,nn) = cos((pi*p*nn)/dfloat(N))
          enddo
       enddo

! and for variation of beta with mode number
! default beta = b0 is obtained for
!      
       if(beta0 .ne. 0d0) then
          b0         = beta0
          b          = abs(tan(betae*Pi/2))
          pv(1:n)    = [(p,p=1,N)]
          beta(1:n)  = 2*atan(a*(log(dfloat(N)/pv)-log(Z))**2+b)/Pi*(1-c+c*erf(d*(log(dfloat(N)/pv)-log(Z))))*b0
          tteff(1:n) = Gamma(1d0/beta)/beta
          if(beta(n/3) .ne. b3old) then
             b3old = beta(n/3)
             open(22,file='nrouse_y_betas.dgli')
                write(22,'(i8,3e16.7)') (p, dfloat(n)/dfloat(p), log(dfloat(n)/dfloat(p)), beta(p), p=1,N)
             close(22)
             open(22,file='nry.gp')
                write(22,'(a        )') 'set term x11 nopersist font Courier '
                write(22,'(a        )') 'set xlabel "{/=20 log N/p}" '
                write(22,'(a        )') 'set ylabel "{/=20 beta}" '
                write(22,'(a        )') 'set yrange [0:1.1] '
                write(22,'(a,f8.4,a )') 'set label 1 "{/=20 Z  = ', Z ,'}" at graph 0.6,0.9  left '
                write(22,'(a,f8.4,a )') 'set label 2 "{/=20 a  = ', a ,'}" at graph 0.6,0.85 left '
                write(22,'(a,f8.4,a )') 'set label 3 "{/=20 c  = ', c ,'}" at graph 0.6,0.8  left '
                write(22,'(a,f8.4,a )') 'set label 4 "{/=20 d  = ', d ,'}" at graph 0.6,0.75 left '
                write(22,'(a,f8.4,a )') 'set label 5 "{/=20 b0 = ', beta0 ,'}" at graph 0.6,0.7 left '
                write(22,'(a,f8.4,a )') 'set label 6 "{/=20 de = ', betae ,'}" at graph 0.6,0.65 left '
                write(22,'(a        )') "plot 'nrouse_y_betas.dgli' using 3:4 with lines lw 3 lc 2" 
                write(22,'(a        )') 'set output "lastgplot.eps" ' 
                write(22,'(a        )') 'set term postscript eps  ' 
                write(22,'(a        )') "plot 'nrouse_y_betas.dgli' using 3:4 with lines lw 3 lc 2" 
                write(22,'(a        )') "pause 2" 
             close(22)
             call system("gnuplot nry.gp &")
           endif
       else
          beta       = 1d0
          tteff      = 1d0
       endif


! estimate pmax for the present q  TO BE DONE
       pmax = N
       
       
! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)
    
          arg2 = 0
          arg20= 0
          do p = 1,pmax
!            rate_p = 2*W*(1-cos((pi*p)/dfloat(N))) * wscal(p)
            rate_p = W * (pi/N)**2 * ppsin(p)    * wscal(p) * tteff(p)
            a0    = - (t*rate_p)**beta(p)
!            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
!            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = cospn(p,nn) * cospn(p,mm)
!            ffc   = ffc / (p**2)
            ffc   = ffc / ppsin(p)  * ascal(p)

            arg2  = arg2  + ffc*e0
            arg20 = arg20 + ffc

          enddo   
          arg2  = arg2  * ff2
          arg20 = arg20 * ff2

          aa1 = arg10
 !         if(aa1.lt.-300.0d0) aa1 = -300.0d0
 !         if(aa1.gt. 300.0d0) aa1 =  300.0d0

          aa2 = arg1+arg2
 !         if(aa2.lt.-300.0d0) aa2 = -300.0d0
 !         if(aa2.gt. 300.0d0) aa2 =  300.0d0


          Sq  = Sq  + exp(aa1)
          Sqt = Sqt + exp(aa2)

        enddo
       enddo

       Sq  = Sq /N
       Sqt = Sqt/N

       if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end

