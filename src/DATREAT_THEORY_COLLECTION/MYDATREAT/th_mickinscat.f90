       function th_mickinscat(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     =====================================================      
!      FUNCTION thgauss (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)                                                                        
! -------> micellkinetics: mickin  <--------                            
!                                                                       
!                                                                       
       implicit none                     
                                                    
       character*8 thnam,parnam(20) 

  		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
                                                                       
       real*4    x, pa, qq, th_mickinscat 
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

       double precision intensity, sum, SQ_cs                                                                        
                                                                        
       integer nsmax, nq 
       double precision sval, dq, q 
       parameter (nsmax=500) 
       common/cscatdat/sval(nsmax),dq,nq 
                                                                        
       integer i, iadda, ier, j 
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
       double precision P(NTmax), PP(NTmax), Wp(NTmax,N), tt(NTmax)  ! Tabulation of Results
                                                                     ! P average mic
                                                                     ! Wp("time",aggrer#) distribution
       common/cresmic/P, PP, tt, Wp, Ntt 
                                                                        
!      P(i)   = average micellesize                                     
!      PP(i)  = average squared micellesize                             
!      tt(i)  = scale time                                              
                                                                        
       CHARACTER*1 RELABS 
       EXTERNAL FCNS, D02EJY, OUTPUST, D02EJW 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'mickins' 
         nparx = 18 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_mickinscat = 0 
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
         parnam(12) = 'mwcore'           ! 1.4    kg/mol                              
         ipa0  (12) = .FALSE.
                                  
         parnam(13) = 'mwshell'          ! 23.9  kg/m**3 
         ipa0  (13) = .FALSE.
                                  
         parnam(14) = 'phicore'          ! 1.0d0 volume fraction polymerchain segments in core
         ipa0  (14) = .FALSE.
          
         parnam(15) = 'phishell'         ! 0.3d0 
         ipa0  (15) = .FALSE.

         parnam(16) = 'smearcore'        !     = 5d-10 
         ipa0  (16) = .FALSE.
         parnam(17) = 'smearshel'        !    =20d-10 
         ipa0  (17) = .FALSE.

         parnam(18) = 'intensity'        ! scattering intensityscale
         ipa0  (18) = .FALSE.

         ipa0  (19) = .TRUE. 
                                                                        
!                                                                       
         th_mickinscat = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
       pget = rho0 
       call        getpar('rho0    ',pget, nopar ,params,napar,mbuf, ier) 
       rho0 = pget 
       pa(19) = pget*pa(9) 
                  
       pget = 1000 ! "arbitrary default"
       call        getpar('time    ',pget, nopar ,params,napar,mbuf, ier) 
                                                                        
       tx           = pget / pa(7) 
       
       Q            = x*1e10                                                                  
                                                                        
       xs           = pa(1) 
       xend         = pa(2) 
       kmax         = NINT(abs(pa(3))) 
       gamma        = pa(4) 
       p_equilibrium= abs(pa(5)) 
       cmc          = pa(6) 
       tol          = pa(8)
       alpha        = pa(10)
       epsilon_f    = pa(11) 
       mw_core      = pa(12) 
       mw_shell     = pa(13)
       phi_core     = pa(14) 
       phi_shell    = pa(15)
       smear_core   = pa(16) 
       smear_shell  = pa(17)
       intensity    = pa(18)
       rho0         = pa(19) 
                                                                        
       newcalc = .FALSE. 
       do i=1,19 
        if(pa(i).ne.pa0(i) .and. ipa0(i) ) newcalc = .TRUE. 
       enddo 
                                                                        
                                                                        
       if(newcalc) then 
                                                                        
         do i=1,19
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
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! scattering  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         den_coremat    = 1000.0         ! kg/m**3                             
         den_shellmat   = 1000.0         ! 
         sigma_core     = 4d13           ! Scatteringlength densities m/m**3=m*
         sigma_shell    = 4d13 

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
                                                                        
       call D02EJF(XS, XEND, kmax, Y, FCNS, D02EJY, TOL, RELABS, OUTPUST,  &
     &             D02EJW, W, IW, IFAIL)                                
                                                                        
        write(6,*)'ifail = ',ifail 
        close(10) 
        close(11) 
        close(12) 
                                                                        
       endif 
                                                                        
!!! ---> extract value from tabulated results !!!                       
       if(tx.lt.tt(1)) then 
         write(6,*)'t out of range low tx=',tx 
         th_mickinscat = 0 
         return 
       endif 
                                                                        
       if(tx.ge.tt(Ntt)) then 
         write(6,*)'t out of range high tx=',tx 
         th_mickinscat = 0 
         return 
       endif 
                                                                        
       do i=1,Ntt-1 
         if(tx.ge.tt(i) .and. tx .lt. tt(i+1) ) then 
          th_mickinscat =                                                    &
     &    log((-exp(P(i))*(tx-tt(i+1))+exp(P(i+1))*(tx-tt(i)))/              &
     &    (tt(i+1)-tt(i)))                                              
          
          do j=1,kmax
           y(j) =                                                             &
     &    log((-exp(WP(i,j))*(tx-tt(i+1))+exp(WP(i+1,j))*(tx-tt(i)))/         &
     &    (tt(i+1)-tt(i)))                                              
         enddo
         endif 
       enddo 
         
       sum = 0
       do j=1,kmax
         sum = sum + SQ_cs(q,j)*y(j)
       enddo

       th_mickinscat = sum * intensity 

                                                               
      END                                           
                                                                        
                                                                        
       subroutine FCNS(x,y,f) 
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
                                                                        
                                                                        
                                                                        
                                                                        
       subroutine OUTPUST(x,y) 
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
                           
       integer N
       parameter (N=500)                                             
       integer NTmax, Ntt 
       parameter(NTmax=500) 
       double precision P(NTmax), PP(NTmax), Wp(NTmax,N), tt(NTmax) 
       common/cresmic/P, PP, tt, Wp, Ntt 
                                                                        
                                                                        
                                                                        
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
                                                                        
                                                                        
                                                                        
                                                                        
!!mi> Find these in th_mickin !!                                                                        
!!mi>                                                                        
!!mi>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mi>!!!!!!! Micellar Modelling Functions                                    
!!mi>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!mi>                                                                        
!!mi>!!  Thinking about the Neu paper it seems that we should find the minima
!!mi>!!  free energy at equilibrium using the translational entropy term  -  
!!mi>!!  we can use the free energy also in the expression for  J(p) but must
!!mi>!!  careful with a ln(phi factor) (note the eqs. 2.12-2.13) and the fact
!!mi>!!  that they have defined the energy the opposite than normal (eps_P = 
!!mi>!!  I suggest to put:                                                   
!!mi>!!                                                                      
!!mi>!!  Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0))-p*f - 
!!mi>!!                                                                      
!!mi>!!  f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  - ln(cmc)         
!!mi>!!                                                                      
!!mi>!!                                                                      
!!mi>!!  and Jk=ka_p*rho(k)*rho(1) - kd_p*rho(k+1)                           
!!mi>!!  Detailed  balance gives using the free energy and rho(p) = exp(- Gp)
!!mi>!!                                                                      
!!mi>!!  ka_p**rho(1)=kd_p*exp(-(Gp(k+1)-Gp(k)))                             
!!mi>!!                                                                      
!!mi>!!  or                                                                  
!!mi>!!   Jk = kd_p*(exp(-(Gp(k+1)-Gp(k)))rho(k)-rho_plus)                   
!!mi>!!  There must be a minus sign now in front because of how Gp is        
!!mi>!!  defined (the energy we also define opposite than Neu).              
!!mi>!!                                                                      
!!mi>!!  I think this will give phi(1)/cmc factor in front of the exp which n
!!mi>!!  If for example cmc is of the order of 10e-5- it should give a better
!!mi>!!  value for rho0. Could you try it?                                   
!!mi>!!                                                                      
!!mi>!!!  Aggregation number dependent Free Energy (omitting translational en
!!mi>!! Dear Michael, thank you for all.                                     
!!mi>!! I saw and immediate error which may very well lead to unlimited growt
!!mi>!! the term of the coronal free energy is not right (not growing with P.
!!mi>!!                                                                      
!!mi>!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0))        
!!mi>!!     *      -p*f                                                      
!!mi>!!     *      -p*log(rho(1))                                            
!!mi>!!                                                                      
!!mi>!! should be replaced by:                                               
!!mi>!!                                                                      
!!mi>!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*p**(3.d0
!!mi>!!     *      -p*f                                                      
!!mi>!!     *      -p*log(rho(1))                                            
!!mi>!!                                                                      
!!mi>!! I would guess that would change things?  In the mean time I will keep
!!mi>!! 
!!mi>!! Dear Michael, 
!!mi>!! I have a small correction. I have ignored a term (P-1)* and wrote P* 
!!mi>!! that gives an excess log(phi(1)) which any way cancel because the Gp
!!mi>!! (P+1) -Gp
!!mi>!! (P) is calculated. This however may give some corrections? (also the 
!!mi>!! ground state): 
!!mi>!! 
!!mi>!! 
!!mi>!! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*
!!mi>!!      *       p**(3.d0/2.d0))
!!mi>!!      *      -p*f
!!mi>!!      *      -(p-1)*log(rho(1))
!!mi>!! 
!!mi>!! (last term changed from ( -p*log(rho(1)))
!!mi>!! 
!!mi>!! which gives a slightly changed F1:
!!mi>!! 
!!mi>!! f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) ) - ((p_equilibrium - 
!!mi>!! 1.0d0)/p_equilibrium)*log(cmc)
!!mi>!! 
!!mi>
!!mi>
!!mi>
!!mi>
!!mi>                                                                        
!!mi>                                                                        
!!mi>       double precision function Gp(k,rho) 
!!mi>!      -----------------------------------                              
!!mi>       implicit none 
!!mi>       integer k 
!!mi>       double precision p,f 
!!mi>       double precision rho(*) 
!!mi>       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
!!mi>       integer kmax 
!!mi>       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax
!!mi>       common/cpar2/ epsilon_f
!!mi>       double precision epsilon_f
!!mi>                                                                        
!!mi>!        Gp = k-1                                                       
!!mi>
!!mi>
!!mi>! Version 7.12.
!!mi>! 
!!mi>! f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  - log(cmc)
!!mi>! 
!!mi>! 
!!mi>! Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*     
!!mi>!      *       p**(3.d0/2.d0))
!!mi>!      *      -p*f
!!mi>!      *      -(p-1)*log(rho(1)) - log(cmc)
!!mi>
!!mi>
!!mi>                                                                        
!!mi>       p  = k 
!!mi>       if(f1.ne.0d0)then 
!!mi>         f = f1 
!!mi>       else                  
!!mi>!!         f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) )  -            &
!!mi>!!     &    ((p_equilibrium - 1.0d0)/p_equilibrium)*log(cmc) 
!!mi>         f = gamma*((5d0/3d0)*p_equilibrium**(-1d0/3d0) ) * (1+epsilon_f) &
!!mi>       &      - log(cmc)
!!mi>       endif 
!!mi>                                                                        
!!mi>!!        Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*   &
!!mi>!!     &       p**(3.d0/2.d0))                                            &
!!mi>!!     &      -p*f                                                        &
!!mi>!!     &      -(p-1)*log(rho(1))                                              
!!mi>          Gp = gamma*(p**(2d0/3d0)+(2d0/3d0)*p_equilibrium**(-5d0/6d0)*   &     
!!mi>       &       p**(3.d0/2.d0))                                            &
!!mi>       &      -p*f                                                        &
!!mi>       &      -(p-1)*log(rho(1)) - log(cmc)
!!mi>                                                                        
!!mi>                                                                        
!!mi>       return 
!!mi>      END                                           
!!mi>                                                                       
!!mi>                                                                        
!!mi>       double precision function Dk(k) 
!!mi>!      -------------------------------                                  
!!mi>       implicit none 
!!mi>       integer k 
!!mi>                                                                        
!!mi>       Dk = 1d0 
!!mi>                                                                        
!!mi>       return 
!!mi>      END                                           
!!mi>                                                                        
!!mi>!    REIDAR: ka_p(k) From Halperin & Alexander 22, 2403 (1989)/ !       
!!mi>!    Dormidontova MM 32 7630 (1999)                                     
!!mi>                                                                        
!!mi>       double precision function ka_p(k) 
!!mi>!      -------------------------------                                  
!!mi>       implicit none 
!!mi>       double precision tau0 
!!mi>       integer k 
!!mi>                                                                        
!!mi>!!! set the timescale !!!!!!!!!!!!!!!!!!!!!!!!!!!                       
!!mi>                     ! set the timescale                                
!!mi>       tau0 = 1.0d-6 
!!mi>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
!!mi>                                                                        
!!mi>       ka_p = 1.d0/tau0*exp(-k**(0.5d0)) 
!!mi>                                                                        
!!mi>       return 
!!mi>      END                                           
!!mi>                                                                        
!!mi>                                                                        
!!mi>       double precision function Jk(k,rho) 
!!mi>!      -----------------------------------                              
!!mi>       implicit none 
!!mi>       integer k 
!!mi>       double precision rho(*) 
!!mi>       double precision rho_plus 
!!mi>!       double precision Dk, Gp, ka_p                                   
!!mi>       double precision  Gp, ka_p 
!!mi>       double precision gamma, p_equilibrium, alpha, f1, cmc, kd_p 
!!mi>       integer kmax 
!!mi>       common/cpar/ gamma, p_equilibrium, alpha, f1, cmc, kd_p, kmax 
!!mi>                                                                        
!!mi>       if(k.ge.kmax) then 
!!mi>         rho_plus = 0 
!!mi>       else 
!!mi>         rho_plus = rho(k+1) 
!!mi>       endif 
!!mi>                                                                        
!!mi>                                                                        
!!mi>!       Jk = exp(Gp(k+1)-Gp(k))*rho(1)*rho(k)-rho_plus                  
!!mi>!       Jk = kd_p*(exp(-(Gp(k+1,rho)-Gp(k,rho)))*rho(k)-rho_plus)       
!!mi>!       Jk = Jk*Dk(k)                                                   
!!mi>                                                                        
!!mi>!  REIDAR: Alternatively expressing via insertion coeffecient, ka_p:    
!!mi>        Jk = rho(1)*(rho(k)- rho_plus*exp((Gp(k+1,rho)-Gp(k,rho)))) 
!!mi>        Jk = ka_p(k)*Jk 
!!mi>
!!mi>        Jk = Jk/(k**alpha)
!!mi>   
!!mi>                                                                        
!!mi> !  We now have a clean expression for ka_p from Halperin & Alexander- s
!!mi> !                                                                      
!!mi>       return 
!!mi>      END                                           
!!mi>                                                                        
                                                                        
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! simple scattering function core-shell                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
! Basic F(Q) Factor for a sphere                                        
!                                                                       
       double precision function FQ_sphere(Q,R) 
!      ----------------------------------------                         
       implicit none 
                                 ! Q-value                              
       double precision Q 
                                 ! Radius fo the sphere                 
       double precision R 
                                                                        
       double precision qr 
       double precision fourpi 
       parameter (fourpi=12.56637061d0) 
                                                                        
       if(Q.eq.0.0) then 
         FQ_sphere = (fourpi/3)*(R**3) 
       else 
         FQ_sphere = -fourpi*(Q*R*cos(Q*R)-sin(Q*R))/Q**3 
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function SQ_cs(Q,P) 
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
                                                                        
                                                                        
       double precision Vc, Vs, Rc, Rs, F, Fc, Fs 
       double precision FQ_sphere 
                                                                        
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                                                        
                                                                        
       Rc=(3/fourpi*P*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
       Rs=(                                                             &
     &       (3/fourpi*P*Mw_shell/(Na*den_shellmat*phi_shell))-Rc**3    &
     &    )**(1d0/3d0)                                                  
                                                                        
       F =   (sigma_core *phi_core *FQ_sphere(Q,Rc) +                   &
     &        sigma_shell*phi_shell*(FQ_sphere(Q,Rs)-FQ_sphere(Q,Rc)))  
                                                                        
! smear interfaces                                                      
       F = F * exp(-(Q*smear_core)**2) 
! other smearing for the shell                                          
       F = F+sigma_shell*phi_shell*(                                    &
     &     FQ_sphere(Q,Rs)*                                             &
     &     (exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2)))          
                                                                        
                                                                        
       SQ_cs = Na * F*F 
                                                                        
!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_cs                 
                                                                        
       return 
      END                                           
