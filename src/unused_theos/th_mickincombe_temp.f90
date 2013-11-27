       function th_mickincomx(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     =====================================================      
!
! -------> micellkinetics: mickin  <--------
!          combined model for <P> and scattering                          
!                                                                       
!                                                                       
       implicit none                     
                                                    
       character*8 thnam,parnam(20) 

       integer     , intent(inout) :: mbuf
       integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
       character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
                                                                       
       real*4    x, pa, qq, th_mickincomx
       dimension pa(20),qq(3) 
       integer   npar, ini, nparx,idum 
                                                                        
       real*4    pa0 
       logical   ipa0 
       dimension pa0(20), ipa0(20) 
       save pa0, ipa0 
                                                                        
       double precision ts 
                                                                        
       integer iadda
       common/thiadd/iadda

!!! specifics of the micellar model !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parameters of the micelles via common block                           
                      
                                ! Molecularweigts kg
       double precision Mw_core    ,  Mw_shell 
                                                    ! Densities       kg
       double precision den_coremat,  den_shellmat 
                                                    ! Volumefractions   
       double precision phi_core,     phi_shell 
                                                    ! Scatteringlength d
       double precision sigma_core,   sigma_shell 
                                                    ! interface width / 
       double precision smear_core,   smear_shell,  exp_shell 
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell,  exp_shell                      
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 

       double precision intensity, sum, SQ_csN, FQ_csN, FQ_csNue                                                           
       double precision eps_mc, eps_mr

       double precision sf1,sf2,sx,ff,Px,Rs,Rc,den,eps,nu,Vc,Vs
       double precision s_of_q, peryev2, blob1, cte_blob, xi
       double precision epsilon_r, epsilon_d, nue
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                                                        
       integer nsmax, nq 
       double precision sval, dq, q 
       parameter (nsmax=500) 
       common/cscatdat/sval(nsmax),dq,nq 
                                                                        
       integer i, ier, j 
       logical newcalc 
       double precision tx, rho0 , v0
       real*4           pget 
                                                                        
       INTEGER N, IW, IFAIL 
       parameter(N=500,iw=(12+N)*N+150) 
       real*8 Xs, XEND, Y(N), TOL, G, W(IW),D02EJW 
       double precision gamma, beta, alpha, f1, coscal, kd_p
       double precision tau_t1, tau_t2, tauscal, gamma0, temp_t
       double precision tau_off, amp_T
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax 
       common/cpar2/ epsilon_f
       double precision epsilon_f

       common/cpar3/rexpon
       double precision rexpon          ! exponent of rate dependence on rho(1)
       common/cpar4/tx, tau_t1, tau_t2, tauscal
       common/cpar5/gamma0, temp_t, tau_off, amp_T                                                          
       save  rho0, xs, xend, tol 
       save  Y 
       save  nparx
                                                                        
       integer NTmax, Ntt 
       parameter(NTmax=500) 
       double precision P(NTmax), PP(NTmax), Wp(NTmax,N), tt(NTmax)  ! Tabulation of Results
                                                                     ! P average mic
                                                                     ! Wp("time",aggrer#) distribution


      
       integer sel 

       common/cresmic/P, PP, tt, Wp, Ntt 
                                                                        
!      P(i)   = average micellesize NEU: gemaess SANS Mittelung==>2,.Moment.1.Moment                        
!      PP(i)  = average squared micellesize                             
!      tt(i)  = scale time                                              
                                                                        
       CHARACTER*1 RELABS 
       EXTERNAL FCNC, D02EJY, OUTPUCT, D02EJW 
                                                                        
                                                           
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'mickincx' 
         nparx = 19 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_mickincomx = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
                                         ! table start                  
         parnam(1) = 't_start' 
         ipa0  (1) = .TRUE. 
                                         ! table end time               
         parnam(2) = 't_end' 
         ipa0  (2) = .TRUE. 
                                         ! dimension (number of equation
         parnam(3) = 'p_max'             
         ipa0  (3) = .TRUE. 
                                         
         parnam(4) = 'gamma0'             ! reduced surface tension / kT         
         ipa0  (4) = .TRUE. 
                                         
         parnam(5) = 'beta'              ! [--> beta=factor in front of p**(3/2) Fm]
         ipa0  (5) = .TRUE. 
                                        
         parnam(6) = 'coscal'            ! [--> concentration scale}
         ipa0  (6) = .TRUE. 
                                         ! time scale                   
         parnam(7) = 'tauscal' 
         ipa0  (7) = .FALSE. 
                                         ! accuracy parameter tol       
         parnam(8) = 'tol' 
         ipa0  (8) = .TRUE. 

         parnam(9) = 'tau_t1'             ! time constant temp decay
         ipa0  (9) = .FALSE. 

         parnam(10) = 'alpha'            ! prefactor for term in GpN default=1
         ipa0  (10) = .TRUE. 

         parnam(11) = 'tau_t2'         ! scaling of Mw_core
         ipa0  (11) = .FALSE. 
                                                                        
                                         ! rho0   ACHTUNG aus Datenparameterfeld                      
         parnam(12) = 'tau_off'         ! scaling of Mw_shell
         ipa0  (12) = .FALSE.
                                  
         parnam(13) = 'amp_T'              ! exponent for rho**(-nu) 
         ipa0  (13) = .FALSE.
                                  
         parnam(14) = 'phicore'          ! 1.0d0 volume fraction polymerchain segments in core
         ipa0  (14) = .FALSE.
          
         parnam(15) = 'phishell'         ! 0.3d0 
         ipa0  (15) = .FALSE.

         parnam(16) = 'rexpon'           ! = 2/3 
         ipa0  (16) = .TRUE.

         parnam(17) = 'expshell'         ! = 2/3 
         ipa0  (17) = .FALSE.

         parnam(18) = 'intensity'        ! scattering intensityscale
         ipa0  (18) = .FALSE.

         parnam(19) = 'xi_blob'          ! blob_xi
         ipa0  (19) = .FALSE.

        
                                                                        
!                                                                       
         th_mickincomx = 0       
           return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
       cte_blob = 10.0
       epsilon_f  = 0
!        epsil_mc =1 
!       epsil_ms =1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! scattering  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!             data should be given as paramteres in the data files   !!
!!             reading to be implememted                              !!
!!             this are the defaults                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         den_coremat    = 1000.0         ! kg/m**3                             
         den_shellmat   = 1000.0         ! 
         sigma_core     = -8.846d13      ! Scatteringlength densities m/m**3=m**-2
         sigma_shell    = 13.83d13       ! differences compared to 90:10 DMF:water
         Mw_core        = 1.4d0          ! kg/mol
         Mw_shell       = 23.9d0         ! kg/mol

         nue            = 1.3333d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                                                        
       pget = rho0 
    
       call    getpar('rho0    ',pget, nopar ,params,napar,mbuf, ier) 
       rho0 = pget 
       pa(nparx+1) = rho0
       ipa0(nparx+1) = .TRUE.
                  
       pget = 1000 ! "arbitrary default"
       call    getpar('time    ',pget, nopar ,params,napar,mbuf, ier) 
       if(ier.eq.0) then                                                
           tx       = pget / pa(7)
           Q        = x*1e10
           sel      = 1             ! compute S(Q)
       else
           tx       = x / pa(7)
           sel      = 0             ! compute <P>
       endif              


!!neu
!!       tx = tx * rho0**2
!!
                                                    
                                                                        
       xs           = pa(1) 
       xend         = pa(2) 
       kmax         = NINT(abs(pa(3))) 
       gamma0        = pa(4) 
       beta         = pa(5) 
       coscal       = pa(6) 
       tol          = pa(8)
       tau_t1        = pa(9)
       alpha        = pa(10)
                     epsilon_r   = 1 
                     epsilon_d   = 1 
       tau_t2       = pa(11) 
       tau_off       = pa(12) 
       amp_t          = pa(13)
       phi_core     = pa(14) 
       phi_shell    = pa(15)
       rexpon       = pa(16) 
                      smear_core   = 5d-10
       exp_shell    = pa(17)
                      smear_shell  = 5d-10
       intensity    = pa(18)
       xi           = pa(19)
                                                                  
!       temp_t = (25.0 + 14*exp(- x/tau_t)) + 273
!       gamma = gamma0*298/temp_t    
                                                                        
       newcalc = .FALSE. 
       do i=1,nparx+1
        if(pa(i).ne.pa0(i) .and. ipa0(i) ) newcalc = .TRUE. 
       enddo 
                                                                        
                                                                        
       if(newcalc) then 
                                                                        
         do i=1,nparx+1
            pa0(i) = pa(i) 
         enddo 
                                                                        
         if(kmax.gt.N) then 
           write(6,*)'Limit kmax = ',kmax,'  to ',N 
           kmax = N 
         endif 
                                                                        
                                                                        
         ifail = 0 
         TOL   = 1D-6 
         RELABS='M' 
                                                                        
         f1            = 0 
                                                                        
                                                                        
         write(6,*)'xs   .......... ',xs 
         write(6,*)'xend .......... ',xend 
         write(6,*)'rho0 .......... ',rho0 
         write(6,*)'gamma.......... ',gamma 
         write(6,*)'beta  .......... ',beta 
         write(6,*)'alpha.......... ',alpha 
         write(6,*)'epsilon_f...... ',epsilon_f 
         write(6,*)'coscal ........ ',coscal 
         write(6,*)'tol  .......... ',tol 
         write(6,*)'kmax .......... ',kmax 
                                                                        
         open(10,file='t.dgli',status='unknown') 
         open(11,file='t3d.dgli',status='unknown') 
         open(12,file='s3d.dgli',status='unknown') 
                                                                        
         do i=10,12 
          write(i,*)'xs   .......... ',xs 
          write(i,*)'xend .......... ',xend 
          write(i,*)'rho0 .......... ',rho0 
          write(i,*)'gamma.......... ',gamma 
          write(i,*)'beta .......... ',beta 
          write(i,*)'alpha.......... ',alpha 
          write(i,*)'epsilon_f...... ',epsilon_f 
          write(i,*)'coscal......... ',coscal 
          write(i,*)'tol  .......... ',tol 
          write(i,*)'kmax .......... ',kmax 
         enddo 
                                                                        
                                                                        
!!!!!!!!!!! NOT YET USED FOR FIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<>!                                                                         
!<>!                                  ! kg/mol                              
!<>!         Mw_core        = 1.4d0 
!<>!                                  ! kg/mol                              
!<>!         Mw_shell       = 23.9d0 
!<>!                                  ! kg/m**3                             
!<>!         den_coremat    = 1000.0 
!<>!                                  ! kg/m**3                             
!<>!         den_shellmat   = 1000.0 
!<>!         phi_core       = 1.0d0 
!<>!         phi_shell      = 0.3d0 
!<>!         sigma_core     = 4d13 
!<>!                                  ! Scatteringlength densities m/m**3=m*
!<>!         sigma_shell    = 4d13 
!<>!                                                                        
!<>!         smear_core     = 5d-10 
!<>!         smear_shell    =20d-10 
!<>! 


                                                                       
         nq = 100 
         dq = 2e9/nq 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
         do i=2,kmax 
          y(i) = 0 
         enddo 
         y(1)  = rho0 
                                                                        
       Ntt = 0 
       ifail = -1
                                                                        
       call D02EJF(XS, XEND, kmax, Y, FCNC, D02EJY, TOL, RELABS, OUTPUCT,  &
     &             D02EJW, W, IW, IFAIL)                                  
                                                                        
        write(6,*)'ifail = ',ifail 
        close(10) 
        close(11) 
        close(12) 
                                                                        
       endif 
                                                                        
!!! ---> extract value from tabulated results !!!                       
       if(sel.eq.0) then
         if(tx.lt.tt(1)) then 
           write(6,*)'t out of range low tx=',tx 
           th_mickincomx = P(1) 
           return 
         endif 
                                                                        
         if(tx.ge.tt(Ntt)) then 
           write(6,*)'t out of range high tx=',tx 
           th_mickincomx = P(Ntt) 
           return 
         endif 
                                                                        
         do i=1,Ntt-1 
           if(tx.ge.tt(i) .and. tx .lt. tt(i+1) ) then 
            th_mickincomx =                                                   &
     &      log((-exp(P(i))*(tx-tt(i+1))+exp(P(i+1))*(tx-tt(i)))/         &
     &      (tt(i+1)-tt(i)))                                              
           return 
           endif 
         enddo 
       endif

       if(sel.eq.1) then
         if(tx.lt.tt(1)) then 
           write(6,*)'t out of range low tx=',tx 
           th_mickincomx = 0 
           return 
         endif 
                                                                        
         if(tx.ge.tt(Ntt)) then 
           write(6,*)'t out of range high tx=',tx 
           th_mickincomx = 0 
           return 
         endif 
                                                                        
         do i=1,Ntt-1 
           if(tx.ge.tt(i) .and. tx .lt. tt(i+1) ) then 
            th_mickincomx =                                                    &
     &      log((-exp(P(i))*(tx-tt(i+1))+exp(P(i+1))*(tx-tt(i)))/              &
     &      (tt(i+1)-tt(i)))                                              
          
            do j=1,kmax
             y(j) =                                                             &
     &      log((-exp(WP(i,j))*(tx-tt(i+1))+exp(WP(i+1,j))*(tx-tt(i)))/         &
     &      (tt(i+1)-tt(i)))                                              
           enddo
           endif 
         enddo 
         
!!<         sum = 0
!!<         do j=1,kmax
!!<           sum = sum + SQ_csN(q,j)*y(j)
!!<         enddo

!!! new: following M.Kotlarchyk and SH Chen, J.Chem.Phys.79,2461(1983)

         Mw_core        = Mw_core  * eps_mc   !**> neu
         Mw_shell       = Mw_shell * eps_mr   !**> neu

!         write(6,*)'Mws=',Mw_core,Mw_shell

         nu  = nue

         sf1 = 0
         sf2 = 0
         sum = 0
         do j=1,kmax
          if((y(j)*j*j/rho0).gt.1d-11) then    ! Auslassen der "leeren" Aggregatzahlen
            ff  = FQ_csNue(q,j,nu)
            sf1 = sf1 + ff    *y(j)
            sf2 = sf2 + ff*ff *y(j)
            sum = sum + y(j)
          endif
         enddo
        
         if(sum.gt.0d0) then
           sf1 = sf1/sum
           sf2 = sf2/sum
         endif

         sf1 = sf1 * sf1 * Na * sum
         sf2 = sf2       * Na * sum

         

         eps = 1d-12  
         den = 0      
         sx  = 0                         
         do j=1,kmax                           ! NO (k=2: unimere ausgelassen !?
          den = den + y(j)
          sx  = sx + j*y(j)
         enddo

         Px   = sx/den
         
         den= den/sx*Na*rho0*(den_coremat*Mw_core+den_shellmat*Mw_shell) &
     &        /(Mw_core+Mw_shell)**2
         
!!         Rc=(3/fourpi*Px*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
!!         Rs=(                                                             &
!!     &       (3/fourpi*Px*Mw_shell/(Na*den_shellmat*phi_shell))-Rc**3    &
!!     &      )**(1d0/3d0) 

         Vc= Px*Mw_core/(Na*den_coremat*phi_core)                                                                
         Rc=(3/fourpi*Vc)**(1d0/3d0) 
         Vs= Px*Mw_shell/(Na*den_shellmat*phi_shell)  * Rc**(-nu)         ! wegen Skalierung rho(r)~(r/Rc)**(-nu)
         Rs=(Vs/fourpi*(3d0-nu)+Rc**(3d0-nu))**(1d0/(3d0-nu))                                                
 
!!??> Versuche alternative den Definition       
!!??>         den = rho0/(fourpi/3*Rs**3)

         Rs  = Rs  * epsilon_r
         den = den * epsilon_d
        


!        s_of_q  = peryev2(q,Rs,den,eps)
         s_of_q  = 1

!dbg         write(6,'(e13.4,f12.6,e13.4,e13.4,f12.4,e13.4,e13.4)')q,s_of_q,Rs,den,Px,sf1,sf2


         sum  = sf2 - sf1 + sf1 * s_of_q
           
         sum  = sum+blob1(q, xi, 1d-9)*sigma_shell**2*xi**3 *cte_blob





!         th_mickincomx = sum * intensity * rho0 * 1d4 
         v0      =  Mw_core/(den_coremat) + Mw_shell/(den_shellmat)
         th_mickincomx = sum * intensity * 1d4 / v0

         call parset('nu      ',sngl(nu),iadda) 
 
         return  
       endif

       th_mickincomx = 0
       write(6,*)'th_mickincomx: should not be here'
                                                               
      END                                           
                                                                        
                                                                        
       subroutine FCNC(x,y,f) 
!      ---------------------                                            
       implicit none 
       real*8 x, y(*), f(*) 
       integer i 
       double precision NaN
       double precision JkN, sum, gamma0, temp_t0
       double precision gamma, beta, alpha, f1, coscal, kd_p 
       double precision tau_t1, tau_t2, tau_off, temp_t
       double precision tx, tauscal, amp_T
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax
       common/cpar4/ tau_t1, tau_t2, tauscal
       common/cpar5/gamma0, temp_t, tau_off, amp_T 
       temp_t = amp_T*(1 - exp(-(x -tau_off)/tau_t2))**4.2*                &
     &   (exp(-(x-tau_off)/tau_t1)) + 298
   
!       temp_t = (25.0 + 14*exp(-x/tau_t)) + 273
!       if(temp_t.lt.297.0) then
!          temp_t =298.0
!         endif
! 
!       if(temp_t.eq.NaN) then
!          temp_t =298.0
!         endif

       gamma = gamma0*298.0/temp_t 
!      write(6,*)'beta   .......... ',beta
!      write(6,*)'coscal   .......... ',coscal
!      write(6,*)'gamma   .......... ',gamma
!      write(6,*)'temp_t   .......... ',temp_t

       sum = 0 
       do i=2,kmax 
        f(i) = JkN(i-1,y)-JkN(i,y) 
        sum = sum + f(i)*i 
       enddo 
                                                                        
       f(1) = -sum 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
       subroutine OUTPUCT(x,y) 
!      ----------------------                                           
! increment x and compute derived values from the solution              
                                                                        
       real*8 x, y(*) 
                                                                        
       double precision sum, sum0, sum2, s 
       integer i,j 
       double precision gamma, beta, alpha, f1, coscal, kd_p
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax 
!       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax
        common/cpar4/ tau_t1, tau_t2, tauscal
       common/cpar5/gamma0, temp_t, tau_off, amp_T                                                                   
       integer nsmax, nq 
       double precision sval,dq 
       parameter (nsmax=500) 
       common/cscatdat/sval(nsmax),dq,nq 
                                                                        
       double precision q, SQ_csN 
                           
       integer N
       parameter (N=500)                                             
       integer NTmax, Ntt 
       parameter(NTmax=500) 
       double precision P(NTmax), PP(NTmax), Wp(NTmax,N), tt(NTmax) 
       common/cresmic/P, PP, tt, Wp, Ntt
!       common/cpar4/tx, tau_t1, tau_t2, tauscal
!       common/cpar5/gamma0, temp_t, tau_off, amp_T
                                                                 
       sum  = 0 
       sum0 = 0 
       sum2 = 0 
       do i=1, kmax 
        sum0= sum0+ y(i) 
        sum = sum + i*y(i) 
        sum2= sum2+ i*i*y(i) 
       enddo 
       s = sum 
                                                                        
       if(sum0.gt.0) then 
         sum  = sum/sum0 
         sum2 = sum2/sum0 
       endif 
                                                                        
!       write(6,'(11e10.3)')x,sum,sum2,s,(y(i),i=1,7) 
       write(10,'(11e10.3)')x,sum,sum2,s,(y(i),i=1,7) 
       write(11,'(A,3e13.6)')'X <k> S: ',x,sum,s 
                                                                        
       do i=1, kmax 
         write(11,*)y(i) 
       enddo 
                                                                        
       if(Ntt.lt.NTmax) then 
         Ntt = Ntt + 1 
         tt(Ntt) = x 
!         P(Ntt)  = sum 
         P(Ntt)  = sum2/sum 
         PP(Ntt) = sum2 
         do i=1, kmax 
           Wp(Ntt,i) = y(i) 
         enddo 
       endif 
                                                                        
!--------------------!                                                  
! compute scattering !                                                  
!--------------------!                                                  
                                                                        
!!<>!       do j=1,nq                                                   
!!<>!        q = dq*(j-1)                                               
!!<>!        sum  = 0                                                   
!!<>!        do i=1, kmax                                               
!!<>!         sum = sum + SQ_csN(q,i)*y(i)                               
!!<>!        enddo                                                      
!!<>!        sval(j) = sum                                              
!!<>! !       write(6,*)j,q,sum                                         
!!<>!        write(12,*)sval(j)                                         
!!<>!       enddo                                                       
!!<>!                                                                   
!!<>!                                                                   
                                                                        
!       x = x*1.05d0                                                    
!       x = x*1.25d0 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        

                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Micellar Modelling Functions                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
!!  Thinking about the Neu paper it seems that we should find the minima
!!  free energy at equilibrium using the translational entropy term  -  
!!  we can use the free energy also in the expression for  J(p) but must
!!  careful with a ln(phi factor) (note the eqs. 2.12-2.13) and the fact
!!  that they have defined the energy the opposite than normal (eps_P = 
!!  I suggest to put:                                                   
!!                                                                      
!!  GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0))-p*f - 
!!                                                                      
!!  f = gamma*((5d0/3d0)*beta**(-1d0/3d0) )  - ln(coscal)         
!!                                                                      
!!                                                                      
!!  and JkN=ka_pN*rho(k)*rho(1) - kd_p*rho(k+1)                           
!!  Detailed  balance gives using the free energy and rho(p) = exp(- GpN)
!!                                                                      
!!  ka_pN**rho(1)=kd_p*exp(-(GpN(k+1)-GpN(k)))                             
!!                                                                      
!!  or                                                                  
!!   JkN = kd_p*(exp(-(GpN(k+1)-GpN(k)))rho(k)-rho_plus)                   
!!  There must be a minus sign now in front because of how GpN is        
!!  defined (the energy we also define opposite than Neu).              
!!                                                                      
!!  I think this will give phi(1)/coscal factor in front of the exp which n
!!  If for example coscal is of the order of 10e-5- it should give a better
!!  value for rho0. Could you try it?                                   
!!                                                                      
!!!  Aggregation number dependent Free Energy (omitting translational en
!! Dear Michael, thank you for all.                                     
!! I saw and immediate error which may very well lead to unlimited growt
!! the term of the coronal free energy is not right (not growing with P.
!!                                                                      
!! GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0))        
!!     *      -p*f                                                      
!!     *      -p*log(rho(1))                                            
!!                                                                      
!! should be replaced by:                                               
!!                                                                      
!! GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0)*p**(3.d0
!!     *      -p*f                                                      
!!     *      -p*log(rho(1))                                            
!!                                                                      
!! I would guess that would change things?  In the mean time I will keep
!! 
!! Dear Michael, 
!! I have a small correction. I have ignored a term (P-1)* and wrote P* 
!! that gives an excess log(phi(1)) which any way cancel because the GpN
!! (P+1) -GpN
!! (P) is calculated. This however may give some corrections? (also the 
!! ground state): 
!! 
!! 
!! GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0)*
!!      *       p**(3.d0/2.d0))
!!      *      -p*f
!!      *      -(p-1)*log(rho(1))
!! 
!! (last term changed from ( -p*log(rho(1)))
!! 
!! which gives a slightly changed F1:
!! 
!! f = gamma*((5d0/3d0)*beta**(-1d0/3d0) ) - ((beta - 
!! 1.0d0)/beta)*log(coscal)
!! 




                                                                        
                                                                        
       double precision function GpN(k,rho) 
!      -----------------------------------                              
       implicit none 
       integer k 
       double precision p,f 
       double precision rho(*)
       double precision tx, tau_t1, tau_t2, tau_off, tauscal
 
       double precision gamma, beta, alpha, f1, coscal, kd_p 
       double precision gamma0
       double precision temp_t, amp_T
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax
       common/cpar2/ epsilon_f
       double precision epsilon_f
       common/cpar4/ tau_t1, tau_t2, tauscal
       common/cpar5/gamma0, temp_t, tau_off, amp_T                                                         
!        GpN = k-1                                                       
!          temp_t = (25.0 + 14*exp(- x/tau_t)) + 273
!          gamma = gamma0*298/temp_t
          beta = dabs(beta)
!            write(6,*)'temp_t   .......... ',temp_t 
!           write(6,*)'gamma   .......... ',gamma 
! Version 7.12.
! 
! f = gamma*((5d0/3d0)*beta**(-1d0/3d0) )  - log(cmc)
! 
! 
! GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0)*     
!      *       p**(3.d0/2.d0))
!      *      -p*f
!      *      -(p-1)*log(rho(1)) - log(cmc)


                                                                        
       p  = k 
!>!       if(f1.ne.0d0)then 
!>!         f = f1 
!>!       else                  
!>!!!         f = gamma*((5d0/3d0)*beta**(-1d0/3d0) )  -            &
!>!!!     &    ((beta - 1.0d0)/beta)*log(cmc) 
!>!         f = gamma*((5d0/3d0)*beta**(-1d0/3d0) ) * (1+epsilon_f) &
!>!       &      - log(cmc)
!>!       endif 
!>!                                                                        
!>!!!        GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0)*   &
!>!!!     &       p**(3.d0/2.d0))                                            &
!>!!!     &      -p*f                                                        &
!>!!!     &      -(p-1)*log(rho(1))                                              
!>!          GpN = gamma*(p**(2d0/3d0)+(2d0/3d0)*beta**(-5d0/6d0)*   &     
!>!       &       p**(3.d0/2.d0))                                            &
!>!       &      -p*f                                                        &
!>!       &      -(p-1)*log(rho(1)) - log(cmc)
!>!

!          GpN = beta*p**(3.0d0/2.0d0)+gamma*p**(2.0d0/3.0d0)   & ! Star-like
!      &       - p*(beta+gamma)                                 &
!      &        - (p-1)*log(rho(1)/coscal) 

                                                                        
         GpN = beta*p**(5.0d0/3.0d0)+gamma*p**(2.0d0/3.0d0)   & ! SURFACTANTS
      &       - p*(beta+gamma)                                &
      &        - (p-1)*log(rho(1)/coscal) 

                                                                                                        
       return 
      END                                           
                                                                       
                                                                        
       double precision function DkN(k) 
!      -------------------------------                                  
       implicit none 
       integer k 
                                                                        
       DkN = 1d0 
                                                                        
       return 
      END                                           
                                                                        
!    REIDAR: ka_pN(k) From Halperin & Alexander 22, 2403 (1989)/ !       
!    Dormidontova MM 32 7630 (1999)                                     
                                                                        
       double precision function ka_pN(k) 
!      -------------------------------                                  
       implicit none 
       double precision tau0
       
!       double precision tau_t1, tau_t2, tau_off, temp_t
!       double precision tx, tauscal, amp_T
!       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax
!       common/cpar4/ tau_t1, tau_t2, tauscal
!       common/cpar5/gamma0, temp_t, tau_off, amp_T 
        
       integer k 
       double precision gamma, beta, alpha, f1, coscal, kd_p 
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, coscal, kd_p, kmax
                                                                       
!!! set the timescale !!!!!!!!!!!!!!!!!!!!!!!!!!!                       
                     ! set the timescale                                
       tau0 = 1.0d-6 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
                                                                        
!       ka_pN = 1.d0/tau0*exp(-beta*k**(0.5d0)) 
 
        ka_pN = 1.d0/tau0*4*3.14159265*k**(0.333333333333333333)                                                                
 !!   TO BE ACTIVATED
 !           temp_t = amp_T*(1 - exp(-(x -tau_off)/tau_t2))**4.2*                &
  !   &   (exp(-(x-tau_off)/tau_t1)) + 298
  !      ka_pN =  ka_pN*temp_t/298  !! Insertion constant scales with viscosity - proportional to 1/T







       return 
      END                                           
                                                                        
                                                                           
       double precision function JkN(k,rho) 
!      -----------------------------------                              
       implicit none 
       integer k 
       double precision rho(*) 
       double precision rho_plus 
!       double precision DkN, GpN, ka_pN                                   
       double precision  GpN, ka_pN 
       double precision gamma, beta, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, beta, alpha, f1, cmc, kd_p, kmax 
       common/cpar3/rexpon
       double precision rexpon          ! exponent of rate dependence on rho(1)

                                                                       
       if(k.ge.kmax) then 
         rho_plus = 0 
       else 
         rho_plus = rho(k+1) 
       endif 
                                                                        
                                                                        
!       JkN = exp(GpN(k+1)-GpN(k))*rho(1)*rho(k)-rho_plus                  
!       JkN = kd_p*(exp(-(GpN(k+1,rho)-GpN(k,rho)))*rho(k)-rho_plus)       
!       JkN = JkN*DkN(k)                                                   
                                                                        
!  REIDAR: Alternatively expressing via insertion coeffecient, ka_pN:    
        JkN = rho(1)*(rho(k)- rho_plus*exp((GpN(k+1,rho)-GpN(k,rho)))) 
        JkN = ka_pN(k)*JkN 

        JkN = JkN/(k**alpha)
!!neu dazu
!        JkN = JkN/(rho(1)**(-rexpon))
                                                                        
 !  We now have a clean expression for ka_pN from Halperin & Alexander- s
 !                                                                      
       return 
      END                                           
                                                                        
                                                                        
                                                                        
!!mi>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mi>!! simple scattering function core-shell                                
!!mi>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mi>!                                                                       
!!mi>! Basic F(Q) Factor for a sphere                                        
!!mi>!                                                                       
!!mi>       double precision function FQ_sphere(Q,R) 
!!mi>!      ----------------------------------------                         
!!mi>       implicit none 
!!mi>                                 ! Q-value                              
!!mi>       double precision Q 
!!mi>                                 ! Radius fo the sphere                 
!!mi>       double precision R 
!!mi>                                                                        
!!mi>       double precision qr 
!!mi>       double precision fourpi 
!!mi>       parameter (fourpi=12.56637061d0) 
!!mi>                                                                        
!!mi>       if(Q.eq.0.0) then 
!!mi>         FQ_sphere = (fourpi/3)*(R**3) 
!!mi>       else 
!!mi>         FQ_sphere = -fourpi*(Q*R*cos(Q*R)-sin(Q*R))/Q**3 
!!mi>       endif 
!!mi>                                                                        
!!mi>       return 
!!mi>      END                                           
!!mi>                                                                        
                                                                        
       double precision function SQ_csN(Q,P) 
!      ------------------------------------                             
! core shell scattering factor as function of aggregation number P      
! rest of parameters via common                                         
       implicit none 
                                                                        
                             ! Q-value                                  
       double precision Q 
                             ! Aggregation number                       
       integer          P 
                                                                        
! parameters of the micelles via common block                           
                                                    ! Molecularweigts kg
       double precision Mw_core    ,  Mw_shell 
                                                    ! Densities       kg
       double precision den_coremat,  den_shellmat 
                                                    ! Volumefractions   
       double precision phi_core,     phi_shell 
                                                    ! Scatteringlength d
       double precision sigma_core,   sigma_shell 
                                                    ! interface width / 
       double precision smear_core,   smear_shell ,  exp_shell
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell,  exp_shell                       
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 
                                                                        
                                                                        
       double precision Vc, Vs, Rc, Rs, F, Fc, Fs 
       double precision FQ_sphere 
                                                                        
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                                                        
                                                                        
       Rc=(3/fourpi*P*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
       Rs=(                                                             &
     &       (3/fourpi*P*Mw_shell/(Na*den_shellmat*phi_shell))+Rc**3    &
     &    )**(1d0/3d0)                                                  
                                                                        
       F =   (sigma_core *phi_core *FQ_sphere(Q,Rc) +                   &
     &        sigma_shell*phi_shell*(FQ_sphere(Q,Rs)-FQ_sphere(Q,Rc)))  
                                                                        
! smear interfaces                                                      
       F = F * exp(-(Q*smear_core)**2) 
! other smearing for the shell                                          
       F = F+sigma_shell*phi_shell*(                                    &
     &     FQ_sphere(Q,Rs)*                                             &
     &     (exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2)))          
                                                                        
                                                                        
       SQ_csN = Na * F*F 
                                                                        
!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_csN                 
                                                                        
       return 
      END                                           


                                                                        
       double precision function FQ_csN(Q,P) 
!      ------------------------------------                             
! core shell scattering factor as function of aggregation number P      
! rest of parameters via common                                         
       implicit none 
                                                                        
                             ! Q-value                                  
       double precision Q 
                             ! Aggregation number                       
       integer          P 
                                                                        
! parameters of the micelles via common block                           
                                                    ! Molecularweigts kg
       double precision Mw_core    ,  Mw_shell 
                                                    ! Densities       kg
       double precision den_coremat,  den_shellmat 
                                                    ! Volumefractions   
       double precision phi_core,     phi_shell 
                                                    ! Scatteringlength d
       double precision sigma_core,   sigma_shell 
                                                    ! interface width / 
       double precision smear_core,   smear_shell,  exp_shell 
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell,  exp_shell                       
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 
                                                                        
                                                                        
       double precision Vc, Vs, Rc, Rs, F, Fc, Fs, nu 
       double precision FQ_sphere
       double precision FQ_shell, Fss
                                                                        
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                                                        
                                                                        
!!       Rc=(3/fourpi*P*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
!!       Rs=(                                                             &
!!     &       (3/fourpi*P*Mw_shell/(Na*den_shellmat*phi_shell))+Rc**3    &
!!     &    )**(1d0/3d0)                                                  


       nu = 0

       Vc= P*Mw_core/(Na*den_coremat*phi_core)                                                                
       Rc=(3/fourpi*Vc)**(1d0/3d0) 
       Vs= P*Mw_shell/(Na*den_shellmat*phi_shell)   * Rc**(-nu)
       Rs=(Vs/fourpi*(3d0-nu)+Rc**(3d0-nu))**(1d0/(3d0-nu))                                                
                                                                        
       Fc =   sigma_core *phi_core *FQ_sphere(Q,Rc) 
       Fs =   sigma_shell*phi_shell*(FQ_sphere(Q,Rs)-FQ_sphere(Q,Rc))
!!       Fss=   sigma_shell*phi_shell*(FQ_shell(Q,Rc,Rs,nu))
       F  =   Fc+Fs  
                                                          
!!       write(6,'(e13.4,2x,2e13.4,2x,2e13.4,f10.3)')q,Fs,Fss,Rc,Rs,nu
              
! smear interfaces                                                      
       F = F * exp(-(Q*smear_core)**2) 
! other smearing for the shell                                          
       F = F+sigma_shell*phi_shell*(                                    &
     &     FQ_sphere(Q,Rs)*                                             &
     &     (exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2)))          
                                                                        
                                                                        
       FQ_csN = F 
                                                                        
!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_csN                 
                                                                        
       return 
      END                                           


                                                                        
       double precision function FQ_csNue(Q,P,nu) 
!      ------------------------------------------                             
! core shell scattering factor as function of aggregation number P
! with shell with r**-nu      
! rest of parameters via common                                         
       implicit none 
                                                                        
                             ! Q-value                                  
       double precision Q 
                             ! Aggregation number                       
       integer          P 
                                                                        
! parameters of the micelles via common block                           
                                                    ! Molecularweigts kg
       double precision Mw_core    ,  Mw_shell 
                                                    ! Densities       kg
       double precision den_coremat,  den_shellmat 
                                                    ! Volumefractions   
       double precision phi_core,     phi_shell 
                                                    ! Scatteringlength d
       double precision sigma_core,   sigma_shell 
                                                    ! interface width / 
       double precision smear_core,   smear_shell ,  exp_shell
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell,  exp_shell                       
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 
                                                                        
                                                                        
       double precision Vc, Vs, Rc, Rs, F, Fc, Fs 
       double precision FQ_sphere, FQ_shell
       double precision nu 
       double precision den_sh1
                                                                        
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                     

       den_sh1 = den_shellmat * (P/20d0)**(exp_shell)
!      ----------------------------------------------
!                                  ----

                                   
       Vc= P*Mw_core/(Na*den_coremat*phi_core)                                                                
       Rc=(3/fourpi*Vc)**(1d0/3d0) 
       Vs= P*Mw_shell/(Na*den_sh1*phi_shell)   * Rc**(-nu)
       Rs=(Vs/fourpi*(3d0-nu)+Rc**(3d0-nu))**(1d0/(3d0-nu))                                                
                                                                        
       Fc =   sigma_core *phi_core *FQ_sphere(Q,Rc) 
       Fs =   sigma_shell*phi_shell*(FQ_shell(Q,Rc,Rs,nu))
       F  =   Fc+Fs  
                                                                        
! smear interfaces                                                      
       F = F * exp(-(Q*smear_core)**2) 
! other smearing for the shell                                          
       F = F+Fs*(exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2))          
                                                                        
                                                                        
       FQ_csNue = F 
                                                                        
!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_csNue                 
                                                                        
       return 
      END                                           




!c  percus-yevick s(q)
!c -------------------
!c lit. m.s. wertheim, prl 10 (1963) 321
!c -------------------------------------
!c
       function peryev2(q,rr,den,eps)
!c      -----------------------------
!c      q     = impulsuebertrag
!c      rr    = teilchenradius
!c      den   = teilchendichte
!c      eps   = realteil ---> 0
!c
       parameter(pi=3.141592654)
       implicit complex*16 (c)
       complex*16           dcmplx, cdexp
       real*8   eta, dimag, deps, dr, dq, q, r, rr, den, eps
       real*8   peryev
!c
       r   = 2*rr
!c --- r ist der durchmesser ! der teilchen !
       dr  = r
       dq  = q
       deps= eps
       eta = pi*den*r*r*r/6.d0
!c
       ct  = dcmplx ( deps , -dq*dr )
       cst = (((1.d0-eta)**2*ct+6.d0*eta*(1.d0-eta))*ct  &
     &      +18.d0*eta*eta)*ct                           &
     &      -12.d0*eta*(1.d0+2.d0*eta)
       clt = 12.d0*eta*((1.d0+0.5d0*eta)*ct+(1.d0+2.d0*eta))
!c
       cgt = (ct*clt) / ( 12.d0*eta*(clt+cst*cdexp(ct)) )
!c
       peryev2 =  1.d0 + den * (4.d0*pi/dq) * dr*dr * dimag( cgt )
!c
       return
       end



       

       double precision function blob1(q, xi, b)
!      -----------------------------------------
!
! Ref: C.M. Marques et al, Eur.Phys.J. B3, 353-358(1998) Eq.11
!
       implicit none
       double precision q, xi, b

       double precision x, c

       x = q*xi
       c = (xi/b)**(1d0/3d0)

       blob1 = c * sqrt(5d0/3d0) * sin((2d0/3d0)*atan(sqrt(27d0/20d0)*x)) &
      &                          / (x**3+27d0/20d0*x**5)**(1d0/3d0)



       return
       end








!!        double precision rhocore = 8.1354d10
!!        double precision rhocorona = 1.0403d11
!!        double precision rho0 = 9.02d+10 ! Scattering Length Density of 90:10 DMF/water

!!        double precision VPEP = 1459.0d0
!!        double precision VPEO = 19478.0d0
!!        double precision phi = 1.0d0



         double precision function FQ_shell(Q,Ri,Ra,nu) 
!        ----------------------------------------------
!        Scattering of a shell with rho ~ r**-nu 
!                        
         implicit none 
                                 ! Q-value                              
         double precision Q 
                                 ! Radii of the shell                 
         double precision Ri, Ra
         double precision nu     ! Exponent 
                                                                        
         double precision eps, erracc, adapint
         integer maxit 
         double precision fourpi 
         parameter (fourpi=12.56637061d0) 

         double precision q_fqs, nu_fqs, ri_fqs
         common /cfqs/ q_fqs, nu_fqs, ri_fqs

         double precision fqs_kernel
         external fqs_kernel

         eps   = (1d-6)*(Ra)**(3-nu)
         maxit = 3000
 
         q_fqs  = q
         nu_fqs = nu
         ri_fqs = Ri

                                                                       
         if(Q.eq.0.0) then 
           FQ_shell = (fourpi/(3d0-nu))*(Ra**(3d0-nu)-Ri**(3d0-nu))*Ri**nu 
         else 
           FQ_shell = fourpi*adapint(fqs_kernel,Ri,Ra,eps,maxit,erracc) 
         endif 
                                                                        
         return 
         END                                           

      
         double precision function fqs_kernel(r)
!        ---------------------------------------
         implicit none
         double precision r
         double precision q_fqs, nu_fqs, ri_fqs
         common /cfqs/ q_fqs, nu_fqs, ri_fqs


         fqs_kernel = sin(q_fqs*r)/(q_fqs*r)*(r/ri_fqs)**(-nu_fqs)*r*r

         return
         end
