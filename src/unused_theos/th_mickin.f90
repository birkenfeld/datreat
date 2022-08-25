       function th_mickin(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
                                            
!      micellkinetics   
!                                                                       
!                                                                       
       implicit none                     
                                                    
       character*8 thnam,parnam(20) 

  		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
                                                                       
       real*4    x, pa, qq, th_mickin 
       dimension pa(20),qq(3) 
       integer   npar, ini, nparx,idum 
                                                                        
       real*4    pa0 
       logical   ipa0 
       dimension pa0(20), ipa0(20) 
       save pa0, ipa0 
                                                                        
       double precision ts 
                                                                        
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
       double precision smear_core,   smear_shell 
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell                       
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 
                                                                        
                                                                        
       integer nsmax, nq 
       double precision sval, dq 
       parameter (nsmax=500) 
       common/cscatdat/sval(nsmax),dq,nq 
                                                                        
       integer i, iadda, ier 
       logical newcalc 
       double precision tx, rho0 
       real*4           pget 
                                                                        
       INTEGER N, IW, IFAIL 
       parameter(N=500,iw=(12+N)*N+150) 
       real*8 Xs, XEND, Y(N), TOL, G, W(IW),D02EJW 
       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax 
       common/cpar2/ epsilon_f
       double precision epsilon_f
                                                                       
       save  rho0, xs, xend, tol 
       save  Y 
                                                                        
       integer NTmax, Ntt 
       parameter(NTmax=500) 
       double precision P(NTmax), PP(NTmax), tt(NTmax) 
       common/cresmic/P, PP, tt, Ntt 
                                                                        
!      P(i)   = average micellesize                                     
!      PP(i)  = average squared micellesize                             
!      tt(i)  = scale time                                              
                                                                        
       CHARACTER*1 RELABS 
       EXTERNAL FCN, D02EJY, OUTPUT, D02EJW 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'mickin' 
         nparx = 11 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_mickin = 0 
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
                                         ! surface tension / kT         
         parnam(4) = 'gamma' 
         ipa0  (4) = .TRUE. 
                                         ! "equilibrium aggregation numb
         parnam(5) = 'p_equil' 
         ipa0  (5) = .TRUE. 
                                         ! critical micelle concentratio
         parnam(6) = 'cmc' 
         ipa0  (6) = .TRUE. 
                                         ! time scale                   
         parnam(7) = 'tauscal' 
         ipa0  (7) = .FALSE. 
                                         ! accuracy parameter tol       
         parnam(8) = 'tol' 
         ipa0  (8) = .TRUE. 

         parnam(9) = 'rhoscal'           ! sacling of input concentrations   
         ipa0  (9) = .TRUE. 

         parnam(10) = 'alpha'            ! prefactor for term in Gp default=1
         ipa0  (10) = .TRUE. 

         parnam(11) = 'epsilonf'         ! prefactor-1 for term F1 in Gp default=0
         ipa0  (11) = .TRUE. 
                                                                        
                                         ! rho0   ACHTUNG aus Datenparameterfeld                      
         ipa0  (12) = .TRUE. 
                                                                        
!                                                                       
         th_mickin = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
       pget = rho0 
       call        getpar('rho0    ',pget, nopar ,params,napar,mbuf, ier) 
       rho0 = pget 
       pa(12) = pget*pa(9) 
                                                                        
                                                                        
       tx           = x / pa(7) 
                                                                        
                                                                        
       xs           = pa(1) 
       xend         = pa(2) 
       kmax         = NINT(abs(pa(3))) 
       gamma        = pa(4) 
       p_equilibrium= abs(pa(5)) 
       cmc          = pa(6) 
       tol          = pa(8)
       alpha        = pa(10)
       epsilon_f    = pa(11) 
       rho0         = pa(12) 
                                                                        
       newcalc = .FALSE. 
       do i=1,12 
        if(pa(i).ne.pa0(i) .and. ipa0(i) ) newcalc = .TRUE. 
       enddo 
                                                                        
                                                                        
       if(newcalc) then 
                                                                        
         do i=1,12
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
         Write(6,*)'rho0 .......... ',rho0 
         write(6,*)'gamma.......... ',gamma 
         write(6,*)'p_eq .......... ',p_equilibrium 
         write(6,*)'alpha.......... ',alpha 
         write(6,*)'epsilon_f...... ',epsilon_f 
         write(6,*)'cmc  .......... ',cmc 
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
          write(i,*)'p_eq .......... ',p_equilibrium 
          write(i,*)'alpha.......... ',alpha 
          write(i,*)'epsilon_f...... ',epsilon_f 
          write(i,*)'cmc  .......... ',cmc 
          write(i,*)'tol  .......... ',tol 
          write(i,*)'kmax .......... ',kmax 
         enddo 
                                                                        
                                                                        
!!!!!!!!!!! NOT YET USED FOR FIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
                                  ! kg/mol                              
         Mw_core        = 1.4d0 
                                  ! kg/mol                              
         Mw_shell       = 23.9d0 
                                  ! kg/m**3                             
         den_coremat    = 1000.0 
                                  ! kg/m**3                             
         den_shellmat   = 1000.0 
         phi_core       = 1.0d0 
         phi_shell      = 0.3d0 
         sigma_core     = 4d13 
                                  ! Scatteringlength densities m/m**3=m*
         sigma_shell    = 4d13 
                                                                        
         smear_core     = 5d-10 
         smear_shell    =20d-10 
                                                                        
         nq = 100 
         dq = 2e9/nq 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
         do i=2,kmax 
          y(i) = 0 
         enddo 
         y(1)  = rho0 
                                                                        
       Ntt = 0 
       ifail = -1
                                                                        
!!       call D02EJF(XS, XEND, kmax, Y, FCN, D02EJY, TOL, RELABS, OUTPUT,D02EJW, W, IW, IFAIL)                                
         call DDRIV1 (kmax, XS, Y, FCN, XEND,-1 , TOL, W, IW,IFAIL)
        write(6,*)'ifail = ',ifail 
        close(10) 
        close(11) 
        close(12) 
                                                                        
       endif 
                                                                        
!!! ---> extract value from tabulated results !!!                       
       if(tx.lt.tt(1)) then 
         write(6,*)'t out of range low tx=',tx 
         th_mickin = P(1) 
         return 
       endif 
                                                                        
       if(tx.ge.tt(Ntt)) then 
         write(6,*)'t out of range high tx=',tx 
         th_mickin = P(Ntt) 
         return 
       endif 
                                                                        
       do i=1,Ntt-1 
         if(tx.ge.tt(i) .and. tx .lt. tt(i+1) ) then 
          th_mickin =                                                   &
     &    log((-exp(P(i))*(tx-tt(i+1))+exp(P(i+1))*(tx-tt(i)))/         &
     &    (tt(i+1)-tt(i)))                                              
         return 
         endif 
       enddo 
                                                                        
      END                                           
                                                                        
                                                                        
       subroutine FCN(x,y,f) 
!      ---------------------                                            
       implicit none 
       real*8 x, y(*), f(*) 
       integer i 
       double precision Jk, sum 
       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax 
                                                                        
       sum = 0 
       do i=2,kmax 
        f(i) = Jk(i-1,y)-Jk(i,y) 
        sum = sum + f(i)*i 
       enddo 
                                                                        
       f(1) = -sum 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
       subroutine OUTPUT(x,y) 
!      ----------------------                                           
! increment x and compute derived values from the solution              
                                                                        
       real*8 x, y(*) 
                                                                        
       double precision sum, sum0, sum2, s 
       integer i,j 
       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax 
                                                                        
       integer nsmax, nq 
       double precision sval,dq 
       parameter (nsmax=500) 
       common/cscatdat/sval(nsmax),dq,nq 
                                                                        
       double precision q, SQ_cs 
                                                                        
       integer NTmax, Ntt 
       parameter(NTmax=500) 
       double precision P(NTmax), PP(NTmax), tt(NTmax) 
       common/cresmic/P, PP, tt, Ntt 
                                                                        
                                                                        
                                                                        
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
                                                                        
       write(6,'(11e10.3)')x,sum,sum2,s,(y(i),i=1,7) 
       write(10,'(11e10.3)')x,sum,sum2,s,(y(i),i=1,7) 
       write(11,'(A,3e13.6)')'X <k> S: ',x,sum,s 
                                                                        
       do i=1, kmax 
         write(11,*)y(i) 
       enddo 
                                                                        
       if(Ntt.lt.NTmax) then 
         Ntt = Ntt + 1 
         tt(Ntt) = x 
         P(Ntt)  = sum 
         PP(Ntt) = sum2 
       endif 
                                                                        
!--------------------!                                                  
! compute scattering !                                                  
!--------------------!                                                  
                                                                        
!!<>!       do j=1,nq                                                   
!!<>!        q = dq*(j-1)                                               
!!<>!        sum  = 0                                                   
!!<>!        do i=1, kmax                                               
!!<>!         sum = sum + SQ_cs(q,i)*y(i)                               
!!<>!        enddo                                                      
!!<>!        sval(j) = sum                                              
!!<>! !       write(6,*)j,q,sum                                         
!!<>!        write(12,*)sval(j)                                         
!!<>!       enddo                                                       
!!<>!                                                                   
!!<>!                                                                   
                                                                        
!       x = x*1.05d0                                                    
       x = x*1.25d0 
                                                                        
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
!!  Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0))-p*f - 
!!                                                                      
!!  f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  - ln(cmc)         
!!                                                                      
!!                                                                      
!!  and Jk=ka_p*rho(k)*rho(1) - kd_p*rho(k+1)                           
!!  Detailed  balance gives using the free energy and rho(p) = exp(- Gp)
!!                                                                      
!!  ka_p**rho(1)=kd_p*exp(-(Gp(k+1)-Gp(k)))                             
!!                                                                      
!!  or                                                                  
!!   Jk = kd_p*(exp(-(Gp(k+1)-Gp(k)))rho(k)-rho_plus)                   
!!  There must be a minus sign now in front because of how Gp is        
!!  defined (the energy we also define opposite than Neu).              
!!                                                                      
!!  I think this will give phi(1)/cmc factor in front of the exp which n
!!  If for example cmc is of the order of 10e-5- it should give a better
!!  value for rho0. Could you try it?                                   
!!                                                                      
!!!  Aggregation number dependent Free Energy (omitting translational en
!! Dear Michael, thank you for all.                                     
!! I saw and immediate error which may very well lead to unlimited growt
!! the term of the coronal free energy is not right (not growing with P.
!!                                                                      
!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0))        
!!     *      -p*f                                                      
!!     *      -p*log(rho(1))                                            
!!                                                                      
!! should be replaced by:                                               
!!                                                                      
!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*p**(3.d0
!!     *      -p*f                                                      
!!     *      -p*log(rho(1))                                            
!!                                                                      
!! I would guess that would change things?  In the mean time I will keep
!! 
!! Dear Michael, 
!! I have a small correction. I have ignored a term (P-1)* and wrote P* 
!! that gives an excess log(phi(1)) which any way cancel because the Gp
!! (P+1) -Gp
!! (P) is calculated. This however may give some corrections? (also the 
!! ground state): 
!! 
!! 
!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*
!!      *       p**(3.d0/2.d0))
!!      *      -p*f
!!      *      -(p-1)*log(rho(1))
!! 
!! (last term changed from ( -p*log(rho(1)))
!! 
!! which gives a slightly changed F1:
!! 
!! f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) ) - ((p_equilibrium - 
!! 1.0d0)/p_equilibrium)*log(cmc)
!! 




                                                                        
                                                                        
       double precision function Gp(k,rho) 
!      -----------------------------------                              
       implicit none 
       integer k 
       double precision p,f 
       double precision rho(*) 
       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax
       common/cpar2/ epsilon_f
       double precision epsilon_f
                                                                        
!        Gp = k-1                                                       


! Version 7.12.
! 
! f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  - log(cmc)
! 
! 
! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*     
!      *       p**(3.d0/2.d0))
!      *      -p*f
!      *      -(p-1)*log(rho(1)) - log(cmc)


                                                                        
       p  = k 
       if(f1.ne.0d0)then 
         f = f1 
       else                  
!!         f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  -            &
!!     &    ((p_equilibrium - 1.0d0)/p_equilibrium)*log(cmc) 
         f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) ) * (1+epsilon_f) &
       &      - log(cmc)
       endif 
                                                                        
!!        Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*   &
!!     &       p**(3.d0/2.d0))                                            &
!!     &      -p*f                                                        &
!!     &      -(p-1)*log(rho(1))                                              
          Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*   &     
       &       p**(3.d0/2.d0))                                            &
       &      -p*f                                                        &
       &      -(p-1)*log(rho(1)) - log(cmc)
                                                                        
                                                                        
       return 
      END                                           
                                                                       
                                                                        
       double precision function Dk(k) 
!      -------------------------------                                  
       implicit none 
       integer k 
                                                                        
       Dk = 1d0 
                                                                        
       return 
      END                                           
                                                                        
!    REIDAR: ka_p(k) From Halperin & Alexander 22, 2403 (1989)/ !       
!    Dormidontova MM 32 7630 (1999)                                     
                                                                        
       double precision function ka_p(k) 
!      -------------------------------                                  
       implicit none 
       double precision tau0 
       integer k 
                                                                        
!!! set the timescale !!!!!!!!!!!!!!!!!!!!!!!!!!!                       
                     ! set the timescale                                
       tau0 = 1.0d-6 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
                                                                        
       ka_p = 1.d0/tau0*exp(-k**(0.5d0)) 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function Jk(k,rho) 
!      -----------------------------------                              
       implicit none 
       integer k 
       double precision rho(*) 
       double precision rho_plus 
!       double precision Dk, Gp, ka_p                                   
       double precision  Gp, ka_p 
       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
       integer kmax 
       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax 
                                                                        
       if(k.ge.kmax) then 
         rho_plus = 0 
       else 
         rho_plus = rho(k+1) 
       endif 
                                                                        
                                                                        
!       Jk = exp(Gp(k+1)-Gp(k))*rho(1)*rho(k)-rho_plus                  
!       Jk = kd_p*(exp(-(Gp(k+1,rho)-Gp(k,rho)))*rho(k)-rho_plus)       
!       Jk = Jk*Dk(k)                                                   
                                                                        
!  REIDAR: Alternatively expressing via insertion coeffecient, ka_p:    
        Jk = rho(1)*(rho(k)- rho_plus*exp((Gp(k+1,rho)-Gp(k,rho)))) 
        Jk = ka_p(k)*Jk 

        Jk = Jk/(k**alpha)
   
                                                                        
 !  We now have a clean expression for ka_p from Halperin & Alexander- s
 !                                                                      
       return 
      END                                           
                                                                        
                                                                        
!!mis>!! find these in th_mickinscat                                                                        
!!mis>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mis>!! simple scattering function core-shell                                
!!mis>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mis>!                                                                       
!!mis>! Basic F(Q) Factor for a sphere                                        
!!mis>!                                                                       
!!mis>       double precision function FQ_sphere(Q,R) 
!!mis>!      ----------------------------------------                         
!!mis>       implicit none 
!!mis>                                 ! Q-value                              
!!mis>       double precision Q 
!!mis>                                 ! Radius fo the sphere                 
!!mis>       double precision R 
!!mis>                                                                        
!!mis>       double precision qr 
!!mis>       double precision fourpi 
!!mis>       parameter (fourpi=12.56637061d0) 
!!mis>                                                                        
!!mis>       if(Q.eq.0.0) then 
!!mis>         FQ_sphere = (fourpi/3)*(R**3) 
!!mis>       else 
!!mis>         FQ_sphere = -fourpi*(Q*R*cos(Q*R)-sin(Q*R))/Q**3 
!!mis>       endif 
!!mis>                                                                        
!!mis>       return 
!!mis>      END                                           
!!mis>                                                                        
!!mis>                                                                        
!!mis>       double precision function SQ_cs(Q,P) 
!!mis>!      ------------------------------------                             
!!mis>! core shell scattering factor as function of aggregation number P      
!!mis>! rest of parameters via common                                         
!!mis>       implicit none 
!!mis>                                                                        
!!mis>                             ! Q-value                                  
!!mis>       double precision Q 
!!mis>                             ! Aggregation number                       
!!mis>       integer          P 
!!mis>                                                                        
!!mis>! parameters of the micelles via common block                           
!!mis>                                                    ! Molecularweigts kg
!!mis>       double precision Mw_core    ,  Mw_shell 
!!mis>                                                    ! Densities       kg
!!mis>       double precision den_coremat,  den_shellmat 
!!mis>                                                    ! Volumefractions   
!!mis>       double precision phi_core,     phi_shell 
!!mis>                                                    ! Scatteringlength d
!!mis>       double precision sigma_core,   sigma_shell 
!!mis>                                                    ! interface width / 
!!mis>       double precision smear_core,   smear_shell 
!!mis>                                                                        
!!mis>                                                    ! Molecularweigts kg
!!mis>       common/csatmod/  Mw_core    ,  Mw_shell,                         &
!!mis>     &                  den_coremat,  den_shellmat,                     &
!!mis>     &                  phi_core,     phi_shell,                        &
!!mis>     &                  sigma_core,   sigma_shell,                      &
!!mis>     &                  smear_core,   smear_shell                       
!!mis>                                                    ! Densities       kg
!!mis>                                                    ! Volumefractions   
!!mis>                                                    ! Scatteringlength d
!!mis>                                                    ! interface width / 
!!mis>                                                                        
!!mis>                                                                        
!!mis>       double precision Vc, Vs, Rc, Rs, F, Fc, Fs 
!!mis>       double precision FQ_sphere 
!!mis>                                                                        
!!mis>       double precision fourpi, Na 
!!mis>       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
!!mis>                                                                        
!!mis>                                                                        
!!mis>       Rc=(3/fourpi*P*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
!!mis>       Rs=(                                                             &
!!mis>     &       (3/fourpi*P*Mw_shell/(Na*den_shellmat*phi_shell))-Rc**3    &
!!mis>     &    )**(1d0/3d0)                                                  
!!mis>                                                                        
!!mis>       F =   (sigma_core *phi_core *FQ_sphere(Q,Rc) +                   &
!!mis>     &        sigma_shell*phi_shell*(FQ_sphere(Q,Rs)-FQ_sphere(Q,Rc)))  
!!mis>                                                                        
!!mis>! smear interfaces                                                      
!!mis>       F = F * exp(-(Q*smear_core)**2) 
!!mis>! other smearing for the shell                                          
!!mis>       F = F+sigma_shell*phi_shell*(                                    &
!!mis>     &     FQ_sphere(Q,Rs)*                                             &
!!mis>     &     (exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2)))          
!!mis>                                                                        
!!mis>                                                                        
!!mis>       SQ_cs = Na * F*F 
!!mis>                                                                        
!!mis>!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_cs                 
!!mis>                                                                        
!!mis>       return 
!!mis>      END                                           
