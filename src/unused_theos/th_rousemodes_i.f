      function thrmi(x,pa,thnam,parnam,npar,ini,nopar,params,napar,mbuf)
c     ===================================================================
c
c -------> rousemod_inc <---------
c
c
       implicit none
       character*8 thnam,parnam(20)
       real     thrmi, x, pa(20),qq(3)
       integer  mbuf,npar,ini,nparx,ier,i
       integer nopar ! Anzahl der Parameter data
       character*80 napar(mbuf)
       real*8 params(mbuf)
       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, l, Sq0, Sq, Sqt
       real*8    a0, sum, sumnorm, q_width, dqw, qzz, fn
       real*8    sumi, sumnormi, rat_ci, sqi, sqti
       real*8    epsilon, diff, dr
       real*4    qget, tget, an
       integer   n, nqw
       logical   sqof0

       double precision Rp, Rw, f0, wl4

       integer iadda
       common/thiadd/iadda

       double precision zpi
       data zpi/6.283185d0/
c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Documentation section                                                      !!
!! see: ~/DATREAT05/datreat1107/ndatreat/doc/THEORIE_DESCRIPTION/th_template.tex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!:
!:
!:
!: \thname{rousemoi}{Rouse model by discrete summation with mode rate modification function and icoherent contribution}
!: 
!:\thdescript{
!:Normalized intermediate scattering function of the Rouse model as obtained by explicit summation.
!:In addition the mode relaxation rates may be modified by a modulation function that reduces the
!:rate for low modes down to $f_0$ and asymptotically reaches 1.   
!:The summation assumes $N$ freely joined (virtual) segments. Their length $l$ is automatically 
!:choosen such that the value of $R$ is preserved. The temperatur may be given, but if it is
!:there, it will be taken from the parameter section in the data file. The diffusion constant is 
!:computed according to Rouse if the parameter $D$ is set to zero. The mode modification is controlled
!:by the parameters $R_p$, $R_w$ and $f_0$. This leads to a transition mode number of $p_t=(R/R_p)^2$
!:and a width of the transition $p_w=R^2/(2 \,R_w \, R_p)$. 
!:}
!:\parameters
!: \parnam{intensit }{ $A_0$      }{amplitude                   }{ -                }{yes}  
!: \parnam{wl4      }{ $Wl^4$     }{Rouse-constant              }{ $\Aa/ns$         }{yes}
!: \parnam{n\_segmnt}{ $N$        }{virtual number of segments  }{ -                }{NO}
!: \parnam{re       }{ $R_e$      }{end-to-end distance         }{ $\Aa$            }{yes}
!: \parnam{temp     }{ $T$        }{temperature                 }{ $K$              }{X}
!: \parnam{com\_diff}{ $D$        }{centre-of-mass diffusion    }{ $cm^2/s$         }{yes}
!: \parnam{q\_width }{ $\Delta q$ }{q-resolution (def=0)        }{ $\rm{\AA}^{-1}$  }{X}
!: \parnam{rmode    }{ $R_m$      }{mode restriction radius     }{ $\rm{\AA}$       }{yes}
!: \parnam{rmode\_wd}{ $R_w$      }{width for mode restriction  }{ $\rm{\AA}$       }{yes}
!: \parnam{flowlim  }{ $f_0$      }{mode scaling factor         }{ -                }{yes}
!:\eparameters


c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'rousemoi'
         nparx = 11
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           thrmi = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '
         parnam(3) = 'n_segmnt'   ! : if <0 return non-normalized s(q,t)
         parnam(4) = 're      '   ! : re**2 = N * b**2
         parnam(5) = 'temp    '
         parnam(6) = 'com_diff'
         parnam(7) = 'q_width '
         parnam(8) = 'rmode   '   ! Uebergangslaenge
         parnam(9) = 'rmode_wd'   ! Schaerfe des Uebergangs
         parnam(10)= 'flowlim '   ! Factor (< 1) fuer langwellige Moden (f=1--> keine Modenanalyse)
         parnam(11)= 'rat_scsi'   ! Ratio: sigma_coherent pro effektivem Segment : sigma_incoherent pro eff. Segment
                                  ! Achtung auf den 4Pi Faktor zwischen ds/do und sigma !
                                  ! 
                                  ! ins Verhaeltnis zu setzen ist:
                                  ! 4* Pi (N_monomer/n_segment*V_monomer*rho_kontrast)**2 /
                                  !        N_monomer/n_segment*(Summe_monomer[sigma_inc])


c
         thrmi = 0
         return
       endif
c
       sqof0 = .false.
c ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       wl4      = abs(pa(2))
       an       = pa(3)
       N        = nint(abs(an))
       R        = abs(pa(4))
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       q_width  = pa(7)
       Rp       = pa(8)
       Rw       = pa(9)
       f0       = pa(10)
       rat_ci   = pa(11)

       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

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
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif
!:
!:\expectedfileparameters
!: \efileparnam{q       }{$q$       }{scattering wavevector length}{ $\Aa^{-1}$       }{XOR tau}
!: \efileparnam{tau     }{$t$       }{Fourier-time                }{ ns               }{XOR q  }
!: \efileparnam{temp    }{$T$       }{temperature                 }{ K                }{overrides thparam}
!:\eexpectedfileparameters
!:


      thrmi   = 0
      sum     = 0
      sumnorm = 0
      sumi    = 0
      sumnormi= 0
      nqw     = 15
      dqw     = 4*q_width/nqw
      if(q_width.eq.0) then
        nqw = 0
        q_width = 1.0d0
      endif 

      do i=-nqw,nqw
       qzz = qz + i*dqw
       if(qzz.gt.0) then 
         fn   = dexp( -(i*dqw/q_width)**2)

c --- include center of mass diffusion ---
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
          dr = diff
       
          call Mrouse_i(qzz,tau,temp,Dr,Wl4,N,R,Rp,Rw,f0,W,l,xi,Sq,Sqt,
     *                                                        Sqi,Sqti)

          sum      = sum + fn*Sqt       
          sumnorm  = sumnorm + Sq*fn

          sumi     = sumi + fn*Sqti       
          sumnormi = sumnormi + Sqi*fn

        endif
       enddo



       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            thrmi = rat_ci*sum -0.333333*sumi 
          else
            if(an.gt.0) then
              thrmi = (rat_ci*sum -0.333333*sumi)/
     *                (rat_ci*sumnorm-0.333333*sumnormi)
            else
              thrmi = rat_ci*sum -0.333333*sumi
            endif
          endif
       else
          thrmi = 0
       endif 
       thrmi = thrmi*a0

!: \fileparameters
!: \fileparnam{l    }{$l $ }{effective virtual segment length}{$\Aa$}
!: \fileparnam{w    }{$W $ }{Rouse rate in                   }{$1/ns$}
!: \fileparnam{xi   }{$\xi$}{friction parameter in SI units  }{SI}
!: \fileparnam{diff }{$D$  }{centre-of-mass diffusion in     }{$cm^2/s$}
!: \efileparameters
!: 
!: \uvariables
!: \uvar{-}{-}{-}{-}
!: \euvariables

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('xi      ',sngl(xi),iadda,ier)     !    SI
       dr        = dr /( 1d-9 / 1d-16 )                      ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
c
       return
       end




       Subroutine Mrouse_i(q,t,temp,Dr,Wl4,N,R,Rp,Rw,f0,W,l,xi,Sq,Sqt,
     *                                                       Sqi,Sqti)
!      ===============================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Wl4   ----> Rouse parameter in nsA units
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule   in A
!    Rp    ----> transition length for mode restriction        in A
!    Rpw   ----> width of transition                           in A
!    f0    ----> long mode rate-factor
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    xi    ----> friction coefficient (Si: Ns/m=kg/s, hier: kg/ns)
!    Sq    <--- S(Q)_coh
!    Sqt   <--- S(Q,t)_coh
!    Sqi   <--- S(Q)_inc
!    Sqt_i <--- S(Q,t)_inc
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt,Sqi,Sqti
       double precision Wl4, Rp,Rw,f0
       double precision pt, pw, f_mode_i

       integer N, nn,mm,p

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2

       integer iout
       
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif


!:  
!: \algorithm{
!:        W   &=& Wl4/l^4                        \cr  
!:        xi  &=& 3*kbt/(W*(l^2))                \cr
!:    pt      &=& (R/Rp)^2                       \cr
!:    pw      &=& (R^2)/(2*Rw*Rp)                \cr
!:    f_{mode}&=& (1-f_0)/(1+\exp((pt-p)/pw))+f0 \cr
!:         a1 &=& -(q^2)*(Dr*t + |n-m|*(l^2)/6)                 \cr
!:         f2 &=& -2*N*(l*q)**2/(3*\pi^2)                           \cr
!:         a0 &=& -t* 2*W*(1-cos(\pi*p/N)) * f_{mode}(p,f0,pt,pw)   \cr
!:         fc &=& cos(\pi*p*n)/N) * cos(\pi*p*m)/N)/p^2             \cr
!:     S(q,t) &=& \sum_{n,m,p} { fc (1-\exp(a0))}                   \cr
!:     th     &=& A_0 \; S(q,t)/S(q,0)           
!: }
!: 


! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2

!      W   = 3*kbt/(xi*(l**2))  ! in 1/ns
       W   = Wl4/l**4           ! in 1/ns    
       xi  = 3*kbt/(W*(l**2))   ! in SI


!      transition mode numbers
       if(Rp.gt.0d0 .and. Rw.gt.0d0) then
         pt = (R/Rp)**2
         pw = (R**2)/(2*Rw*Rp)
       else
         pt = 0d0
         pw = 1d0
       endif
!       write(6,'(a,2f12.6)')' pt,pw= ',pt,pw


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       Sqi = 0
       Sqti= 0

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)
    
          arg2 = 0
          arg20= 0
          do p = 1,N
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))        ! this is the rate !!!
            tau_p = tau_p * f_mode_i(dfloat(p),f0,pt,pw) ! mode rate modificator
            a0    = -t*tau_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)
            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)

            arg2  = arg2  + ffc*e0
            arg20 = arg20 + ffc

          enddo   
          arg2  = arg2  * ff2
          arg20 = arg20 * ff2

          aa1 = arg10
          if(aa1.lt.-300.0d0) aa1 = -300.0d0
          if(aa1.gt. 300.0d0) aa1 =  300.0d0

          aa2 = arg1+arg2
          if(aa2.lt.-300.0d0) aa2 = -300.0d0
          if(aa2.gt. 300.0d0) aa2 =  300.0d0


          Sq  = Sq  + exp(aa1)
          Sqt = Sqt + exp(aa2)

          if(nn.eq.mm) then
           Sqi = Sqi + exp(aa1)
           Sqti= Sqti+ exp(aa2)
          endif


        enddo
       enddo


       if(iout().gt.0)write(6,'(1x,7E14.6)')q,t,Sq,Sqt,Sqi,Sqti, w 

       return
       end



       double precision function f_mode_i(p,f0,pt,pw)
!      ----------------------------------------------
!      Uebergangsfunktion 
!
       implicit none

       double precision f0, pt, pw, p
 
       if(pw.ne.0) then
        f_mode_i=(1d0-f0)/(1d0+exp((pt-p)/pw))+f0
       else
        if(p.lt.pt) f_mode_i = f0
        if(p.eq.pt) f_mode_i = f0+0.5*(1-f0)
        if(p.gt.pt) f_mode_i = 1d0
       endif
 
       return
       end
       
