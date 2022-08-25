         function riggch (x,pa,thnam,parnam,npar,ini, nopar ,params,napar,mbuf) 
!                                                                     ! 
!        rigid gaussian chain:   riggch                               ! 
!                                                                     ! 
!        according to:                                                ! 
!        C.M. Marques and G.H. Fredrickson                            ! 
!        J.Phys. II France 7 (1997) 1805-1816                         ! 
!        Eqn. 11 and 14                                               ! 
!        Attention: observe error in eqn. 11 m<n, (m-n)->(n-m)        ! 
!        Parameters:                                                  ! 
!        ===========                                                  ! 
!        volfrac        -> Phi-Polymer                                ! 
!        rho            -> density of polymer in g/cm**3              ! 
!        molwgt         -> molecular weight   in g/mol                ! 
!        b_chain        -> scattering length density in cm**-2        ! 
!        nseg           -> no. of segments                            ! 
!        rg             -> radius of gyration                         ! 
!        lpers          -> presistence length                         ! 
!        select         -> selects: 0=rigid gaussin chain (rgc)       ! 
!                       +> selects: 1=approximant to rgc Eq 21        ! 
!                       +> selects: 2=debye function                  ! 
!                                                                     ! 
!        Required data set parameters:                                ! 
!        =============================                                ! 
!        bsolv          -> scattering length density of solvent       ! 
!                                                                     ! 
                                                                        
         implicit none 
                                                                        
         DOUBLE PRECISION Pi 
         Parameter       (Pi=3.141592654d0) 
                                                                        
                                                                        
         CHARACTER*8     thnam,parnam(20) 
 		 integer :: mbuf
		 integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
         character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		 real, intent(inout) :: params(mbuf)             ! value des parameters n
		 REAL            riggch, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx 
         DIMENSION       pa(20),qq(3) 
         DATA            zpi/6.283185/ 
         real            bsolv_s 
         double precision Navog 
         parameter (Navog=6.022045d23) 
                                                                        
         double precision volfrac, lpers, rg, lp , bseg, Y, q 
         double precision vchain 
         double precision sq_rgc, sq_rgc_approx21, qrgq 
         double precision b_chain, bsolv, drho, rho , molwgt 
         integer          nseg, select 
                                                                        
         INTEGER  ier, iot 
!                                                                       
!                                                                       
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'riggch  ' 
         nparx = 8 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           riggch   = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
                                   ! -> Phi-Polymer                     
         parnam(1) = 'volfrac ' 
                                   ! -> density of polymer in g/cm**3   
         parnam(2) = 'rho     ' 
                                   ! -> molecular weight   in g/mol     
         parnam(3) = 'molwgt  ' 
                                   ! -> scattering length density in cm*
         parnam(4) = 'b_chain ' 
                                   ! -> no. of segments                 
         parnam(5) = 'nseg    ' 
                                   ! -> radius of gyration              
         parnam(6) = 'rg      ' 
                                   ! -> presistence length              
         parnam(7) = 'lpers   ' 
                                   ! -> selects: 0=rigid gaussin chain (
         parnam(8) = 'select  ' 
                                   ! -> selects: 1=approximant to rgc   
                                   ! -> selects: 2=debye function       
!                                                                       
         riggch  = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
                                                                        
        q       = x 
                                                                        
        volfrac = pa(1) 
        rho     = pa(2) 
        molwgt  = pa(3) 
        b_chain = pa(4) 
        nseg    = NINT(pa(5)) 
        rg      = pa(6) 
        lpers   = pa(7) 
        select  = NINT(pa(8)) 
                                                                        
                                                                        
        call getpar('bsolv   ',bsolv_s,nopar ,params,napar,mbuf, ier) 
                                                                        
        drho = b_chain-bsolv_s 
                                                                        
! ---> volume of one chain in cm**3 <---                                
        vchain = molwgt/rho/Navog 
        call setpar('vchain   ',sngl(vchain),nopar ,params,napar,mbuf, ier) 
                                                                        
! ---> Parameters for rigid gaussian chain from Rg and lpers:           
                                                 ! "chemical" persistenc
        lp   = 0.33333333d0*nseg*(lpers/Rg)**2 
                                                 ! segment length       
        bseg = lpers/lp 
                                                                        
        call  set_gr_chain(nseg,bseg,lp) 
                                                                        
        if(select.eq.0) then 
         call setpar('lp_rgc    ',sngl(lp),nopar ,params,napar,mbuf, ier) 
         call setpar('bseg_rgc  ',sngl(bseg),nopar ,params,napar,mbuf, ier) 
         Y = sq_rgc(q) / nseg 
        else if(select.eq.1) then 
         call setpar('lp_rgca   ',sngl(lp),nopar ,params,napar,mbuf, ier) 
         call setpar('bseg_rgca ',sngl(bseg),nopar ,params,napar,mbuf, ier) 
         Y = sq_rgc_approx21(q) / nseg 
        else if(select.eq.2) then 
         call setpar('rg_debye  ',sngl(Rg),nopar ,params,napar,mbuf, ier) 
         qrgq = (Rg*q)**2 
         if(qrgq.lt.1d-8) qrgq = 1d-8 
         Y = (2/qrgq**2)*(exp(-qrgq)-1+qrgq) 
        else 
         Write(6,*)'select not within allowed range 0..2, :',select 
         Y = 0 
        endif 
                                                                        
                                                                        
        riggch = volfrac * (vchain *drho**2) * Y 
                                                                        
        return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!       program rigid_gaussian_chain                                    
!!      ----------------------------                                    
!                                                                       
!       double precision q, y1, y2                                      
!       double precision sq_grc, sq_rgc_approx21                        
!                                                                       
!       double precision bseg, lpers, Rg, qrgq, deb                     
!       integer          nseg                                           
!                                                                       
!                                                                       
! 1     continue                                                        
!       write(6,*)'Enter: Nseg, bseg, lpers >'                          
!       read(5,*,end=999) nseg, bseg, lpers                             
!                                                                       
!       call  set_gr_chain(nseg,bseg,lpers)                             
!                                                                       
!                                                                       
!! compute Rg:                                                          
!                                                                       
!       Rg = sqrt(2.0d0*nseg*lpers*bseg*bseg/6.0d0)                     
!                                                                       
!       open(10,file='tmp.dgli',status='unknown')                       
!       write(10,'(a,i5)')      'Nseg=',nseg                            
!       write(10,'(a,F8.2)')    'bseg=',bseg                            
!       write(10,'(a,F8.2)')    'lper=',lpers                           
!       write(10,'(a,F8.2)')    'Rg  =',Rg                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!       do i=0,200                                                      
!        q  = i*0.01d0                                                  
!                                                                       
!        qrgq = (Rg*q)**2                                               
!        if(qrgq.lt.1d-8) qrgq = 1d-8                                   
!        deb = (2/qrgq**2)*(exp(-qrgq)-1+qrgq)*nseg                     
!                                                                       
!        y1 = sq_rgc(q)                                                 
!        y2 = sq_rgc_approx21(q)                                        
!        write(10,'(4F18.8)')q,y1,y2,deb                                
!       enddo                                                           
!                                                                       
!       close(10)                                                       
!       goto 1                                                          
! 999   continue                                                        
!                                                                       
!       end                                                             
!                                                                       
                                                                        
       double precision function Rmn_squared(mm,nn) 
!      ------------------------------------------                       
!                                                                       
!--> nach: C.M. Marques and G.H. Fredrickson, J.Phys II France 7 (1997) 
!          "Rigid Gaussian Chains I: the Scattering Function"           
!                                                                       
!          Equation 11                                                  
!                                                                       
! !!!! --> this equation contains an error (m-n) is to be replaced by (n
!          it holds for n>m and NOT for m>n   !!!!                      
!                                                                       
                                                                        
       implicit none 
                                                                        
                                   ! segment length                     
       double precision    b 
                                   ! presistence length                 
       double precision    lp 
                                    ! number of segments per molecule   
       integer             Ns 
                                                                        
       common/cgrc/ b,lp,Ns 
                                                                        
                                                                        
       integer m, n , mm, nn 
       double precision    f1, f2, x1, x2 
                                                                        
                                                                        
                                                                        
! --> the expression only holds for m > n <--                           
                                                                        
      if(mm.eq.nn) then 
       Rmn_squared = 0.0d0 
       return 
      endif 
                                                                        
      if(nn.gt.mm) then 
       m = mm 
       n = nn 
      else 
       m = nn 
       n = mm 
      endif 
                                                                        
      if(m.gt.Ns) then 
       write(6,*)'m > Ns',m,Ns 
      endif 
                                                                        
! ---------------------------------------------                         
                                                                        
       f1 = b**2*Ns*sinh(1.0d0/lp)/(1.0d0+Ns/tanh(Ns/lp)*sinh(1.0d0/lp)) 
       f2 = 1.0d0/(4*sinh(Ns/lp)*sinh(0.5d0/lp)**2) 
                                                                        
       x1 = (n-m)/tanh(0.5d0/lp) 
       x2 = cosh((Ns-2*m)/lp)+cosh((Ns-2*n)/lp)                         &
     &    + 2*cosh((Ns+m-n)/lp)-2*cosh((Ns-m-n)/lp)-2*cosh(Ns/lp)       
                                                                        
                                                                        
                                                                        
                                                                        
       Rmn_squared = f1*(x1+f2*x2) 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
       double precision function sq_rgc(q) 
!      -----------------------------------                              
!                                                                       
! S(Q) for the rigid gaussian chain using gaussian quadratic form proper
!                                                                       
       implicit none 
                                                                        
       double precision q 
       double precision Rmn_squared 
       double precision sum, arg 
                                                                        
       integer n,m 
                                                                        
                                    ! segment length                    
       double precision    b 
                                    ! presistence length                
       double precision    lp 
                                    ! number of segments per molecule   
       integer             Ns 
                                                                        
       common/cgrc/ b,lp,Ns 
                                                                        
                                                                        
       sum = 0.0d0 
                                                                        
       do m=1,Ns 
        do n=1,m-1 
         arg = -q*q * Rmn_squared(m,n) / 6.0d0 
         if(arg.lt.-100d0) arg=-100.d0 
         sum = sum + exp(arg) 
        enddo 
       enddo 
                                                                        
       sum = 2*sum+Ns 
                                                                        
       sq_rgc  = sum/Ns 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function sq_rgc_approx21(q) 
!      --------------------------------------------                     
!                                                                       
! S(Q) Approximant according to Eqn 21 in J.Physique II 7 (1997) 1805   
!                                                                       
                                                                        
       implicit none 
                                                                        
       double precision Pi 
       parameter (Pi= 3.141592654d0) 
                                                                        
                                                                        
       double precision q, xh 
                                                                        
                                    ! segment length                    
       double precision    b 
                                    ! presistence length                
       double precision    lp 
                                    ! number of segments per molecule   
       integer             Ns 
       common/cgrc/ b,lp,Ns 
                                                                        
       double precision a 
                                                                        
       a = b 
!      a = sqrt(2.14)*b                                                 
                                                                        
                                                                        
!                                                                       
! first alterantive                                                     
!                                                                       
!        xh = q*q*Ns*2*lp*a*a/(12*(1+q*lp*a*sqrt(Pi/6)))                
!                                                                       
! second alternative                                                    
!                                                                       
       xh = q*q*Ns*2*lp*a*a/(12*(1+q*lp*a*(Pi/6))) 
                                                                        
                                                                        
       sq_rgc_approx21 = Ns/(1.0d0+xh) 
                                                                        
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
       subroutine set_gr_chain(N,bs,l) 
!      ------------------------------                                   
                                                                        
       implicit none 
                                                                        
                              ! no segments/chain                       
       integer N 
                              ! segment length, presistence length      
       double precision bs, l 
                                                                        
                                    ! segment length                    
       double precision    b 
                                    ! presistence length                
       double precision    lp 
                                    ! number of segments per molecule   
       integer             Ns 
       common/cgrc/ b,lp,Ns 
                                                                        
                                                                        
       b  = bs 
       lp = l 
       Ns = N 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
