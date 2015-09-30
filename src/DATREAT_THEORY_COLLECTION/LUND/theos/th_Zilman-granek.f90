        function th_zglam(x,pa,thnam,parnam,npar,ini, nopar ,params,napar,mbuf) 
!       ===================================================             
!                                                                       
! -------> zilman-granek for lammelar phase  <--------                  
!                                                                       
!                                                                       
!        implicit none                                                  
        character*8 thnam,parnam(20) 
        dimension pa(20),qq(3)
				integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n

        integer i, ifail, maxitro 
        integer na 
        real*8 epsro, erracro, romin, romax 
        real*8 a_miu, a_ro, miu, ro, qu, sum_u, sum_ro 
        real*8 qh, akappa, eta, xi, aa, ti, axi, sint2 
        real*8 sigma, pi, amp 
        real*8 f_ro, s17aef, zadapint 
        real*8 o_q, o_xi, o_axi, o_aa, o_na, o_akappa, o_eta, o_amp 
        common/com_pa/qh, akappa, eta, xi, aa, ti, sint2, axi, na 
        common/miuro/ro 
                                                                        
        save o_q, o_xi, o_axi, o_aa, o_na, o_akappa, o_eta, o_amp 
        logical changed 
                                                                        
        data zpi/6.283185/ 
                                                                        
        external f_ro 
!                                                                       
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'zglam' 
         nparx = 7 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_zglam  = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam(1) = 'xi' 
         parnam(2) = 'axi' 
         parnam(3) = 'aa' 
         parnam(4) = 'na' 
         parnam(5) = 'akappa' 
         parnam(6) = 'eta' 
         parnam(7) = 'amp' 
                                                                        
         th_zglam = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here -----                                      
                              !!--> correlation length                  
       xi  = abs(pa(1)) 
                              !!--> scaling factor of the limit over the
       axi = pa(2) 
                                       !! vector (k) spectrum k_min=pi/(
       aa  = abs(pa(3)) 
                              !!--> short wavelength limit k_max=pi/(na 
       na = pa(4) 
                              !!--> bending modulus (k_BT)              
       akappa = abs(pa(5)) 
                              !!--> solvent viscosity [k_BT A^3/ns]     
       eta = abs(pa(6)) 
                            !!--> should be ~ 1                         
       amp = abs(pa(7)) 
                                                                        
       q = 1.0 
                                                                        
       call  getpar('q       ',q,nopar ,params,napar,mbuf, ier) 
       q = abs(q) 
                                                                        
                                                                        
       if (o_q.NE.q.OR.o_xi.NE.xi.OR.o_axi.NE.axi.OR.o_aa.NE.aa.OR.     &
     &    o_na.NE.na.OR.o_akappa.NE.akappa.OR.o_eta.NE.eta.OR.          &
     &    o_amp.NE.amp) then                                            
                                                                        
         changed = .true. 
         o_q = q 
         o_xi = xi 
         o_axi = axi 
         o_aa = aa 
         o_na = na 
         o_akappa = akappa 
         o_eta = eta 
         o_amp = amp 
       else 
          changed = .false. 
       endif 
                                                                        
        qh = q 
                                                                        
        alambda = 8.d0 
        pi = 4.*atan(1.) 
                            !! -->the angle between the normal to the   
        sint2 = 1d-3 
                            !! lamellar planes and the q-vector is ~0   
                            !! in our set-up. (No angle average required
                                                                        
        ti = x 
        write(*,*) 't=', ti 
                                                                        
        epsro = 1.d-4 
        maxitro = 1.d8 
        erracro = 0.0d0 
                                                                        
        romin = na*aa 
        romax = axi*xi 
                                                                        
                                                                        
       if (changed) then 
        th_zglam_0 = zadapint(f_ro, romin, romax, epsro, maxitro, erracro) 
        write(*,*) 'th0=', th_zglam_0 
       endif 
                                                                        
       th_zglam = zadapint(f_ro, romin, romax, epsro, maxitro, erracro) 
       th_zglam = amp*th_zglam/th_zglam_0 
!                                                                       
       return 
      END                                           
                                                                        
                                                                        
        real*8 function f_ro(a_ro) 
!====================================================================   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        implicit none 
        integer i, ifail, maxitu 
        integer  na 
        real*8 a_ro, miu, ro 
        real*8 qh, akappa, eta, xi, aa, ti, axi, sint2 
        real*8 umin, umax, sigma 
        real*8 pi, epsu, erracu 
        real*8 f_u, uadapint, s17aef, sum_u 
        common/com_pa/qh, akappa, eta, xi, aa, ti, sint2, axi, na 
        common/miuro/ro 
                                                                        
        external f_u 
        ro = a_ro 
        pi = 4.*atan(1.) 
                                                                        
        epsu = 1.d-8 
        maxitu = 1.d8 
        erracu = 0.0d0 
                                                                        
!        sigma = axi*xi                                                 
                           !! set the limits of integration k_min and k_
        umin = pi/(axi*xi) 
                           !! over the undulation mode spectrum         
        umax = pi/aa 
                                                                        
        ifail=-1 
                                                                        
!       write(*,*)'cost2=', dsqrt(1.d0-sint2**2)                        
                                                                        
        sum_u = uadapint(f_u, umin, umax, epsu, maxitu, erracu) 
                                                                        
        f_ro = a_ro*s17aef(qh*a_ro*sint2,ifail)*                        &
     &         exp(-1./(2.*pi*akappa)*(1.d0-sint2**2)*qh**2*sum_u)      
!     *         -a_ro**2/(2.*sigma**2)) !! the smooth cutoff using a gau
                                        !! no more necessary!           
                                                                        
        return 
      END                                           
                                                                        
                                                                        
        real*8 function f_u(u) 
!====================================================================   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        implicit none 
        integer i, ifail 
        integer  na 
        real*8 pi, ro, u, umin, umax, sigma, miu 
        real*8 qh, akappa, eta, xi, aa, ti, axi, sint2 
        real*8 sigmau, omega 
        real*8 s17aef 
        common/com_pa/qh, akappa, eta, xi, aa, ti, sint2, axi, na 
        common/miuro/ro 
                                                                        
        pi = 4.*atan(1.) 
                                                                        
        ifail=-1 
                                                                        
        omega = akappa/(4.0d0*eta) 
                                                                        
        f_u = (1.d0-exp(-omega*ti*u**3)*s17aef(u*ro, ifail))*           &
     &        1.d0/u**3                                                 
!     *(1.-exp(-u**2/(2.*sigmau**2)))                                   
                                                                        
                                                                        
        return 
      END                                           
                                                                        
                                                                        
       function zgaus8a (f,xu,xo) 
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
       zgaus8a = sum * range 
                                                                        
       return 
      END                                           
!*ds                                                                    
!*ed                                                                    
       function zadapint(f,a,b,epsilon,maxiter,erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
       implicit real*8 (a-h,o-z) 
!                                                                       
       parameter( maxstack = 200 ) 
!                 --------------->  stacktiefe                          
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack) 
!                                                                       
       logical cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
       real*4 xxxx,yyyy,ptxf 
                                                                     !! 
       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20) 
!                                                                       
       external f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
         iterationcounter = 0 
         erg              = 0.0d0 
         itop      = 0 
         xa(itop)  = a 
         xb(itop)  = b 
         s(itop)   = zgaus8a( f, a, b ) 
                                                                        
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
              s(itop)  = zgaus8a( f, xa(itop), xb(itop) ) 
              itop             = itop +1 
              xa(itop) =  xa(itop-2) 
              xb(itop) =  xa(itop-1) 
              s(itop)  = zgaus8a( f, xa(itop), xb(itop) ) 
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
     &   'zadapint fatal error iterationnumber exceeded!'               
    3    continue 
                                                                        
                                                                        
         zadapint = erg 
                                                                        
         return 
                                                                        
      END                                           
                                                                        
       function ugaus8a (f,xu,xo) 
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
       ugaus8a = sum * range 
                                                                        
       return 
      END                                           
!*ds                                                                    
!*ed                                                                    
       function uadapint(f,a,b,epsilon,maxiter,erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
       implicit real*8 (a-h,o-z) 
!                                                                       
       parameter( maxstack = 200 ) 
!                 --------------->  stacktiefe                          
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack) 
!                                                                       
       logical cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
       real*4 xxxx,yyyy,ptxf 
                                                                     !! 
       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20) 
!                                                                       
       external f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
         iterationcounter = 0 
         erg              = 0.0d0 
         itop      = 0 
         xa(itop)  = a 
         xb(itop)  = b 
         s(itop)   = ugaus8a( f, a, b ) 
                                                                        
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
              s(itop)  = ugaus8a( f, xa(itop), xb(itop) ) 
              itop             = itop +1 
              xa(itop) =  xa(itop-2) 
              xb(itop) =  xa(itop-1) 
              s(itop)  = ugaus8a( f, xa(itop), xb(itop) ) 
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
     &  'uadapint fatal error iterationnumber exceeded!'                
    3    continue 
                                                                        
                                                                        
         uadapint = erg 
                                                                        
         return 
                                                                        
      END                                           
