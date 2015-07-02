       function th_zgbicbic(x,pa,thnam,parnam,npar, ini, nopar ,params,napar,mbuf) 
!      ===================================================              
!                                                                       
!       zilman-granek for bicontinuous phase
!
         ! cite:
         ! M. Monkenbusch et al., CONDENSED MATTER J. Phys.: Condens. Matter 17 (2005) S2903–S2909
         ! doi:10.1088/0953-8984/17/31/017
         ! and the Original Zilman Granek PRL 
         
!                                                                       

!        implicit none                                                  
        CHARACTER(8) thnam,parnam (20) 
        DIMENSION pa (20), qq (3) 
        integer :: mbuf
		integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
        character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n
		INTEGER i, ifail, lim, npts, nxi, na
        REAL(8) ti, ws 
        REAL(8) a_miu, a_ro, miu, ro, qu, sum_u, sum_ro 
        REAL(8) akappa, eta, omega, aa, qh, axi, amp 
        REAL(8) a1, a2, sigma, xi, pi 
        REAL(8) f_miu, s17aef, d01ahf, eps, errrel
        REAL(8)  o_na, o_nxi 
        REAL(8) o_q, o_xi, o_aa, o_akappa, o_eta, o_axi, o_amp 
                                                                        
        common/com_pa/qh, akappa, eta, xi, aa, axi, amp, ti, ws, nxi, na 
        common/miuro/miu, ro 
                                                                        
        save o_q, o_xi, o_nxi, o_aa, o_na, o_akappa, o_eta, o_axi, o_amp 
        logical changed 
                                                                        
        data zpi/6.283185/ 
                                                                        
        external f_miu 
!                                                                       
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'zgbic' 
         nparx = 8 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_zgbic  = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam(1) = 'xi' 
         parnam(2) = 'nxi' 
         parnam(3) = 'aa' 
         parnam(4) = 'na' 
         parnam(5) = 'akappa' 
         parnam(6) = 'eta' 
         parnam(7) = 'axi' 
         parnam(8) = 'amp' 
                                                                        
         th_zgbic = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
                          !!--> correlation length                      
       xi  = pa(1) 
                          !!--> additional scaling factor for r_max limi
       nxi = pa(2) 
                          !!--> molecular size [A]                      
       aa = pa(3) 
                          !!--> na=0                                    
       na = pa(4) 
                          !!--> bending modulus [k_BT]                  
       akappa = pa(5) 
                          !!--> solvent viscosity [k_BT A^3/ns]         
       eta = pa(6) 
                          !!--> scaling factor of the limit over the mod
       axi = pa(7) 
                                !! vector (k) spectrum: k_min=pi/(axi xi
                          !!--> scaling factor amp*S(q,t)/S(q)          
       amp = pa(8) 
                                                                        
                                                                        
       call  getpar('q       ',q,nopar ,params,napar,mbuf, ier) 
       q=abs(q) 
                                                                        
                                                                        
       if (o_q.NE.q.OR.o_xi.NE.xi.OR.o_nxi.NE.nxi.OR.o_aa.NE.aa.OR.     &
     &  o_na.NE.na.OR.o_akappa.NE.akappa.OR.o_eta.NE.eta.OR.            &
     &  o_axi.NE.axi.OR.o_amp.NE.amp) then                           
                                                                        
         changed = .true. 
         o_q = q 
         o_xi = xi 
         o_nxi = nxi 
         o_aa = aa 
         o_na = na 
         o_akappa = akappa 
         o_eta = eta 
         o_axi = axi 
         o_amp = amp 
       else 
          changed = .false. 
       endif 
                                                                        
       qh = q 
       pi = 4.*atan(1.) 
       a1 = 0.0d0 
       a2 = 1.0d0 
                                                                        
       eps = 1.d-4 
       lim = 5.d3 
       ifail = -1
       errrel = 0.1
       npts = 1000 
       ti = x 
       write(*,*) 'q=', qh 
       write(*,*) 't=', ti 
                                                                        
                                                                        
       if (changed) then 
         th_zgbic_0 = d01ahf(a1, a2, eps, npts, errrel, f_miu, lim, ifail) 
         call  setpar('th_zgbic_0 ',th_zgbic_0,nopar ,params,napar,mbuf, ier) 
        write(*,*) 'th0=', th_zgbic_0 
       endif 
                                                                        
       call  getpar('th_zgbic_0 ',th_zgbic_0,nopar ,params,napar,mbuf, ier) 
       th_zgbic = d01ahf(a1, a2, eps, npts, errrel, f_miu, lim, ifail) 
       th_zgbic = amp*th_zgbic/th_zgbic_0 
       write(*,*)'th_zgbic, th_zgbic_0, amp=',th_zgbic, th_zgbic_0, amp
!                                                                       
       return 
       END                                           
                                                                        
         real*8 function f_miu(a_miu) 
!====================================================================   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        implicit none 
        integer i, ifail, maxitro, nxi, na 
        real*8 a_miu, miu, ro 
        real*8 qh, akappa, eta, xi, aa, axi, amp 
        real*8 ti, ws 
        real*8 romin, romax, sigma 
        real*8 pi, epsro, erracro 
        real*8 f_ro1, radapint, s17aef, sum_ro 
        common/com_pa/qh, akappa, eta, xi, aa, axi, amp, ti, ws, nxi, na 
        common/miuro/miu, ro 
                                                                        
        external f_ro1 
                                                                        
        miu = a_miu 
        pi = 4.*atan(1.) 
                                                                        
        epsro = 1.d-8 
        maxitro = 1.d8 
        erracro = 0.0d0 
                                                                        
        romin = na*aa 
        romax = nxi*axi*xi 
                                                                        
        sum_ro = radapint(f_ro1, romin, romax, epsro, maxitro, erracro) 
                                                                        
        f_miu = sum_ro 
                                                                        
        return 
       END                                           
                                                                        
                                                                        
        real*8 function f_ro1(a_ro) 
!====================================================================   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        implicit none 
        integer i, ifail, maxitu, nxi, na 
        real*8 a_ro, miu, ro 
        real*8 qh, akappa, eta, xi, aa, axi, amp 
        real*8 ti, ws 
        real*8 umin, umax, sigma 
        real*8 pi, epsu, erracu 
        real*8 f_u1, uadapin1, s17aef, sum_u 
        common/com_pa/qh, akappa, eta, xi, aa, axi, amp, ti, ws, nxi, na 
        common/miuro/miu, ro 
                                                                        
        external f_u1 
        ro = a_ro 
        pi = 4.*atan(1.) 
                                                                        
        epsu = 1.d-6 
        maxitu = 1.d8 
        erracu = 0.0d0 
                                                                        
        sigma = axi*xi 
                                                                        
        umin = pi/(axi*xi) 
        umax = pi/aa 
                                                                        
        ifail=-1 
                                                                        
        sum_u = uadapin1(f_u1, umin, umax, epsu, maxitu, erracu) 
                                                                        
        f_ro1 =a_ro*s17aef(qh*a_ro*sqrt(1-miu**2),ifail)*               &
     &        exp(-1./(2.*pi*akappa)*miu**2*qh**2*sum_u                 &
     &        -a_ro**2/(2.*sigma**2))                                   
                                                                        
                                                                        
        return 
       END                                           
                                                                        
                                                                        
        real*8 function f_u1(u) 
!====================================================================   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        implicit none 
        integer i, ifail, nxi, na 
        real*8 pi, ro, u, umin, umax, sigma, miu 
        real*8 qh, akappa, eta, xi, omega, aa, axi, amp 
        real*8 ti, sigmau, ws 
        real*8 s17aef 
        common/com_pa/qh, akappa, eta, xi, aa, axi, amp, ti, ws, nxi, na 
        common/miuro/miu, ro 
                                                                        
        pi = 4.*atan(1.) 
                                                                        
        ifail=-1 
                                                                        
        omega = akappa/(4.0d0*eta) 
!        write(*,*) 'omega=', omega                                     
                                                                        
        f_u1 = (1.d0-exp(-omega*ti*u**3)*s17aef(u*ro, ifail))*          &
     &        1.d0/u**3                                                 
                                                                        
                                                                        
        return 
       END                                           
                                                                        
                                                                        
       function rgaus8a (f,xu,xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
!-----------------------------------------------------------------------
       parameter ( ndim = 8 ) 
       parameter ( ndim2=ndim/2) 
       implicit real*8 (a-h,o-z) 
       dimension a(ndim2),x(ndim2) 
       data a / 0.362683783378362d0,                                    &
     &          0.313706645877887d0,                                    &
     &          0.222381034453374d0,                                    &
     &          0.101228536290376d0/                                    
       data x / 0.183434642495650d0,                                    &
     &          0.525532409916329d0,                                    &
     &          0.796666477413627d0,                                    &
     &          0.960289856497536d0/                                    
                                                                        
       xave = (xo+xu)*0.5d0 
       range= (xo-xu)*0.5d0 
       sum = 0.d0 
       do i=1,ndim2 
         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) ) 
       enddo 
       rgaus8a = sum * range 
                                                                        
       return 
       END                                           
!*ds                                                                    
!*ed                                                                    
       function radapint(f,a,b,epsilon,maxiter,erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
       use outlev
       implicit real*8 (a-h,o-z) 
!                                                                       
       parameter( maxstack = 200 ) 
!                 --------------->  stacktiefe                          
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack) 
!                                                                       
!       logical cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
!       real*4 xxxx,yyyy,ptxf 
                                                                     !! 
!       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20) 
!                                                                       
       external f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
         iterationcounter = 0 
         erg              = 0.0d0 
         itop      = 0 
         xa(itop)  = a 
         xb(itop)  = b 
         s(itop)   = rgaus8a( f, a, b ) 
                                                                        
         do i=1,maxiter 
           if (itop.ge.(maxstack-2)) then 
              erg = erg + s(itop) 
              xbb = xb(itop) 
              itm = itop-1 
              do itop=itm,0,-1 
                 if(xa(itop).eq.xbb) goto 1 
              enddo 
    1         continue 
              write(6,*)'warning! adaptint stack overflow!' 
           else 
              iterationcounter = iterationcounter + 1 
              itop             = itop +1 
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0 
              xb(itop) =  xb(itop-1) 
              s(itop)  = rgaus8a( f, xa(itop), xb(itop) ) 
              itop             = itop +1 
              xa(itop) =  xa(itop-2) 
              xb(itop) =  xa(itop-1) 
              s(itop)  = rgaus8a( f, xa(itop), xb(itop) ) 
              error   =  dabs((s(itop)+s(itop-1))-s(itop-2)) 
                                                                        
              if (iot.gt.2) then 
                 write(6,'(1x,i3,i5,4e13.6)')itop,iterationcounter,     &
     &                 xa(itop),xb(itop),(s(itop)+s(itop-1)),error      
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
    2            continue 
              endif 
            endif 
            if (itop.le.0) goto 3 
         enddo 
         write(6,*)                                                     &
     &   'radapint fatal error iterationnumber exceeded!'               
    3    continue 
                                                                        
                                                                        
         radapint = erg 
                                                                        
         return 
                                                                        
       END                                           
                                                                        
       function ugaus8a1 (f,xu,xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
!-----------------------------------------------------------------------
       parameter ( ndim = 8 ) 
       parameter ( ndim2=ndim/2) 
       implicit real*8 (a-h,o-z) 
       dimension a(ndim2),x(ndim2) 
       data a / 0.362683783378362d0,                                    &
     &          0.313706645877887d0,                                    &
     &          0.222381034453374d0,                                    &
     &          0.101228536290376d0/                                    
       data x / 0.183434642495650d0,                                    &
     &          0.525532409916329d0,                                    &
     &          0.796666477413627d0,                                    &
     &          0.960289856497536d0/                                    
                                                                        
       xave = (xo+xu)*0.5d0 
       range= (xo-xu)*0.5d0 
       sum = 0.d0 
       do i=1,ndim2 
         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) ) 
       enddo 
       ugaus8a1 = sum * range 
                                                                        
       return 
       END                                           
!*ds                                                                    
!*ed                                                                    
       function uadapin1(f,a,b,epsilon,maxiter,erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!      
       use outlev                                                                 
       implicit real*8 (a-h,o-z) 
!                                                                       
       parameter( maxstack = 200 ) 
!                 --------------->  stacktiefe                          
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack) 
!                                                                       
!       logical cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
!       real*4 xxxx,yyyy,ptxf 
                                                                     !! 
!       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20) 
!                                                                       
       external f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
         iterationcounter = 0 
         erg              = 0.0d0 
         itop      = 0 
         xa(itop)  = a 
         xb(itop)  = b 
         s(itop)   = ugaus8a1( f, a, b ) 
                                                                        
         do i=1,maxiter 
           if (itop.ge.(maxstack-2)) then 
              erg = erg + s(itop) 
              xbb = xb(itop) 
              itm = itop-1 
              do itop=itm,0,-1 
                 if(xa(itop).eq.xbb) goto 1 
              enddo 
    1         continue 
              write(6,*)'warning! adaptint stack overflow!' 
           else 
              iterationcounter = iterationcounter + 1 
              itop             = itop +1 
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0 
              xb(itop) =  xb(itop-1) 
              s(itop)  = ugaus8a1( f, xa(itop), xb(itop) ) 
              itop             = itop +1 
              xa(itop) =  xa(itop-2) 
              xb(itop) =  xa(itop-1) 
              s(itop)  = ugaus8a1( f, xa(itop), xb(itop) ) 
              error   =  dabs((s(itop)+s(itop-1))-s(itop-2)) 
                                                                        
              if (iot.gt.2) then 
                 write(6,'(1x,i3,i5,4e13.6)')itop,iterationcounter,     &
     &                 xa(itop),xb(itop),(s(itop)+s(itop-1)),error      
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
    2            continue 
              endif 
            endif 
            if (itop.le.0) goto 3 
         enddo 
         write(6,*)                                                     &
     &   'uadapin1 fatal error iterationnumber exceeded!'               
    3    continue 
                                                                        
                                                                        
         uadapin1 = erg 
                                                                        
         return 
                                                                        
       END                                           
