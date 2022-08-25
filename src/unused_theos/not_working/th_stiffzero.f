      function th23(x,pa,thnam,parnam,npar,ini,nopar,params,napar,mbuf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ============================================             !!  
!                                                              !!  
!           -------> stiffpol <--------                        !!  
!                                                              !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                             !!
!! Rouse dynamik fuer steife Polymere nach:                    !!
!!                                                             !!
!! L.Harnau,R.Winkler,P.Reineker, J.Chem.Phys 106 p2469 (1997) !!
!! ----------------------------------------------------------- !!
!!                                                             !!
!! Normaler Output:   S(q,t)/S(q)                              !!
!! ----------------                                            !!
!! dabei ist            x --> t                                !!
!! Datenparameter       q        in A**-1                      !!
!!                      temp     in K                          !!
!!                                                             !!
!!                                                             !!
!! Parameteruebersetzung:                                      !!
!!   gamma ----> N*zeta/L            (gamx)                    !!
!!   p     ----> 1/(l0*Cinfinity)                              !!
!!                                                             !!
!! Einheiten hier: g A ns                                      !!
!!           dh.   typische Werte von zeta --> 1e-17           !!
!!                                                             !!
!!                                                             !!
!! Sonderfunktion:    dSigma(q)/dOmega                         !!
!! ---------------                                             !!
!! aktiv wenn Datenparameter  q NICHT vorhanden                !!
!!                      x --> q                                !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         include 'stiffzero_inc1.inc'
   
         CHARACTER*8     thnam,parnam(20)
         REAL            th23, x, pa, qq, zpi, xh, vol_frac, paold(20)
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         DOUBLE PRECISION ampli, tau, Qex, diff, sq0, qold
         double precision p_old, L_old, Rg, sqtnag
         integer          maxm_old
         REAL             qh, bh

         integer mbuf
         integer nopar ! Anzahl der Parameter data
         character*80 napar(mbuf)
         real*8 params(mbuf)

         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

         double precision bsolv, bchain, molwgt, vchain, drho, rho
         double precision sigma_factor
         double precision Navog
         parameter (Navog=6.022045d23)

         SAVE    paold, sq0
         LOGICAL change
        
         integer typ
         SAVE typ
         Data typ/0/
        


c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'stiffpol'
         nparx = 8
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th23   = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'amplitud'          
         parnam(2) = 'zeta    '
         parnam(3) = 'l0      '
         parnam(4) = 'nsegment'
         parnam(5) = 'cinfty  '
         parnam(6) = 'modes   '
         parnam(7) = 'temp    '
         parnam(8) = 'diffscal'
c
         th23  = 0
         return
       endif
c
c ---- calculate theory here -----
        call allparset

        ampli   = pa(1)  
        zeta    = abs(pa(2))
        l0      = abs(pa(3))
        nsegment= NINT(abs(pa(4)))
        cinf    = abs(pa(5))
        maxmode = NINT(pa(6))
        temp    = pa(7)
        diffsca = abs(pa(8))
 

        ier = -1
        call parget('q       ',qh ,iadda,ier)
        if(ier.eq.0) then  !
          qex       = qh
          tau       = x
          if(typ.ne.0) Write(6,*)' Computing S(Q,t)/S(Q) ..'
          typ       = 0
        else               ! assume that S(Q) shall be computed
          qex       = x
          tau       = 0
          if(typ.eq.0) Write(6,*)' Computing S(Q)=S(Q,t=0)..'
          typ       = 1
          bh        = 0
          call parget('bsolv   ',bh,iadda,ier)
          bsolv     = bh
          bh        =1e11 
          call parget('bchain  ',bh,iadda,ier)
          bchain    = bh
          bh        = 1e4
          call parget('molwgt  ',bh,iadda,ier)
          molwgt    = bh
          bh        = 1.0
          call parget('rho     ',bh,iadda,ier)
          rho       = bh
          drho = bchain-bsolv
          ! ---> volume of one chain in cm**3 <--- 
          vchain = molwgt/rho/Navog
          call parset('vchain   ',sngl(vchain),iadda)
          sigma_factor = vchain*drho**2
        endif

        call parget('temp    ',qh ,iadda,ier)
        if(ier.eq.0) then
          temp      = qh
        endif


        L = nsegment * l0
        p = 1.0d0/(l0*Cinf)

        Rg= sqrt(abs(L/(6*p)-1/(4*p**2)+1/(4*p**3*L)-
     *               (1-exp(-2*p*L))/(8*p**4*L**2)))

        
        gamx = nsegment *zeta / L



        diff = kb*Temp/nsegment/zeta*diffsca
        call        parset('diff    ',sngl(diff),iadda,ier)
        call        parset('rg      ',sngl(Rg)  ,iadda,ier)
        call        parset('lpers   ',sngl(0.5/p),iadda,ier)


        change = .false.     
        if(L.ne.L_old .or. p.ne.p_old .or. maxmode.ne.maxm_old) then
          call create_zeros
          change   = .true.
          L_old    = L
          p_old    = p
          maxm_old = maxmode
        endif
        
        if(change .or. qold.ne.qex .and. typ.ne.1) then
         sq0  = sqtnag(qex,0.0d0)
         write(6,*)'S(q,t=0) ',qex,sq0
         qold = qex         
        endif 
   
        if(typ.eq.0) then
         th23 = ampli * sqtnag(qex,tau)/sq0 
        else
         th23 = ampli * sqtnag(qex,0.0d0) * sigma_factor
        endif
       

        end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Optional cross-checks                                              !!
!!   not needed for evaluation purposes                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine check_stiffzero
!     --------------------------
      include 'stiffzero_inc1.inc'
      

      double precision gu, go, epsilon, erraccu, adapint
      integer          maxit
      double precision ye, yo, yoe, qq, tt, yn, ds
      integer i,j
      integer IERSVR, IPACT,  ISACT

      double precision phi_sst, s1, s2, t, phi_sst_rouse,sqt,sqts
      double precision sqtii, sqtnag

      double precision q_vec(10),yv(10)

      double precision dot_kernel_ev, dot_kernel_odd, dot_kernel_oe,
     *                 dot_kernel_ll
      external         dot_kernel_ev, dot_kernel_odd, dot_kernel_oe,
     *                 dot_kernel_ll

      q_vec(1) = 0.03
      q_vec(2) = 0.04
      q_vec(3) = 0.08
      q_vec(4) = 0.10
      q_vec(5) = 0.15
      q_vec(6) = 0.20
      q_vec(7) = 0.25
      q_vec(8) = 0.30
      q_vec(9) = 0.35
      q_vec(10) = 0.40





      maxmode =50
      l0   = 1.54
      Cinf = 6.73
      p = 1.0d0/(l0*Cinf)
      L = 212
      tim = 1
      q   = 0.1
      diffsca = 1

      zeta = 0.43d-17
      Temp = 470

!      write(6,*)'enter p, maxmode'
!      read(5,*,end=999)p, maxmode
    
      write(6,*)'p = ',p,'   maxmode = ',maxmode
    
      call create_zeros


      maxit    = 1000
      epsilon  = 1d-8
      gu       = -L/2
      go       =  L/2

      do ip=1,10
       do iq=1,10
!       ye  = adapint(dot_kernel_ev,gu,go,epsilon,maxit,erraccu)
!       yo  = adapint(dot_kernel_odd,gu,go,epsilon,maxit,erraccu)
!       yoe = adapint(dot_kernel_oe,gu,go,epsilon,maxit,erraccu)
       yoe = adapint(dot_kernel_ll,gu,go,epsilon,maxit,erraccu)
!       write(6,'(2I4,3F18.6)')ip,iq,ye,yo,yoe
       write(6,'(2I4,F18.6)')ip,iq,yoe
       enddo
      enddo




! display the zeros !
      write(6,*)'odd zeros deviation'
      do i=1,iz1
       write(6,*)i,exzeros1(i)*L/Pi-(1+(i-1)*2)
      enddo
      write(6,*)'even zeros deviation'
      do i=1,iz2
       write(6,*)i,exzeros2(i)*L/Pi-(i*2)
      enddo

     
!       goto 2222
! Write S(Q,t)
       open(10,file='sqtv2.dgli',status='unknown')
       write(10,*)'L    = ',L
       write(10,*)'p    = ',p
       write(10,*)'zeta = ',zeta
       write(10,*)'temp = ',temp
       write(10,*)'diffs= ',diffsca
       write(10,*)'maxm = ',maxmode
       write(10,'(a,a,10F12.7)')'tau     ','  q=',q_vec

       tt = 0
       do j=1,10
        qq = q_vec(j)
        yv(j) = sqtnag(qq,tt)
       enddo
       write(10,'(11F12.7)')tt,yv

       do i=1,100
        tt = 0.01d0*2**(14.0d-2*(i-1))
        do j=1,10
         qq = q_vec(j)
         yv(j) = sqtnag(qq,tt)
        enddo
        write(10,'(11F12.7)')tt,yv
       enddo
!       do i=1,9
!        tt = i*1.0d0
!        do j=1,10
!         qq = q_vec(j)
!         yv(j) = sqtnag(qq,tt)
!        enddo
!        write(10,'(11F12.7)')tt,yv
!       enddo
!       do i=1,10
!        tt = i*10.0d0
!        do j=1,10
!         qq = q_vec(j)
!         yv(j) = sqtnag(qq,tt)
!        enddo
!        write(10,'(11F12.7)')tt,yv
!       enddo
       close(10)
 2222  continue

! 222  continue
!      write(6,*)'S(Q,t) enter t:'
!      read(5,*,end=888)tt
!      open(10,file='sqt.dgli',status='unknown')
!      write(10,*)'t= ',t
!      write(10,*)'L= ',L
!      write(10,*)'p= ',p
!      ye = 0
!      yo = 0
!      yoe= 0
!      do i=1,50
!        qq = (i-1)*0.01d0
!!        ye = sqts(qq,tt)
!!        yo = sqtii(qq,tt)
!        yoe= sqtnag(qq,tt)
!        write(10,'(4F15.7)')qq,ye,yo,yoe
!      enddo
!      close(10)
!      goto 222
! 888  continue

 !111  continue
      open(10,file='phiss.dgli',status='unknown')
      write(6,*)'ds '
      read(5,*,end=999)ds
      write(10,*)'ds = ',ds
      diffsca = 0
      do i=1,200
       t  = 0.01d0*(2**((i-1)/12))
       ye = 0
       yo = 0
       yn = 0
     
       do j = 1, L/l0
        s1  = j*l0-L/2
        s2  = s1+ds
        if(s2.lt.L/2) then 
         ye = ye + phi_sst(s1,s2,t)
         yo = yo + phi_sst_rouse(s1,s2,t)
         yn = yn+1
        endif
       enddo
     
       ye = ye/yn
       yo = yo/yn
       write(10,*)t,ye,yo
      enddo
      close(10)
!      goto 111
 999  continue

      return
      end


      subroutine create_zeros
!     -----------------------

      include 'stiffzero_inc1.inc'

      double precision z1, z2, alpha, zero, dalpha
      double precision z1v(10000), z2v(10000), a(10000)
      integer i, iout
      
      double precision errabs,errrel,eps,eta
      integer          nroot,itmax, maxm

      double precision gu, go, epsilon, erraccu, adapint,yi
      integer          maxit

      double precision normin_odd, normin_even

      double precision dot_kernel_ev, dot_kernel_odd, dot_kernel_oe
      external         dot_kernel_ev, dot_kernel_odd, dot_kernel_oe
     
      external z1,z2
      
      maxm = maxmode/2
      if(maxm.gt.maxn) maxm=maxn
      write(6,*)'maxm = ',maxm

      iz1 = 0
      iz2 = 0
      dalpha = 0.1d0*pi/L
      do i=1,10000
        if(iz1.ge.maxm .or. iz2.ge.maxm) goto 10
        alpha = dalpha*i
        a(i)    = alpha
        z1v(i)  = z1(alpha)
        z2v(i)  = z2(alpha)
        if(i.gt.1) then
         if(z1v(i)*z1v(i-1).le.0.0d0) then
           zero = a(i)-(a(i)-a(i-1))/(z1v(i)-z1v(i-1))*z1v(i)
           iz1   = iz1+1
           zeros1(iz1) = zero
           if(iout().gt.0) then
            write(6,*)'zero1 at ',zero*L/pi,' (',z1(zero),')'
           endif
         endif
         if(z2v(i)*z2v(i-1).le.0.0d0) then
           zero = a(i)-(a(i)-a(i-1))/(z2v(i)-z2v(i-1))*z2v(i)
           iz2  = iz2+1
           zeros2(iz2) = zero
           if(iout().gt.0) then
            write(6,*)'zero2 at ',zero*L/pi,' (',z2(zero),')'
           endif
         endif
        endif
      enddo

 10   continue
      write(6,*)'Zeros found, odd:',iz1,'  even:',iz2

      errabs = 1d-8
      errrel = 1d-8
      eps    = dalpha
      eta    = dalpha
      nroot  = iz1
      itmax  = 50

      call dzreal(z1,errabs,errrel,eps,eta,nroot,itmax,
     *            zeros1,exzeros1,infozer1)
      if(iout().gt.0) then
      write(6,'(a,f14.8,2x,f14.8,2x,i8,2x,f14.10)')
     *  (' odd : ',exzeros1(i),exzeros1(i)*L/pi,infozer1(i),
     *    z1(exzeros1(i)),i=1,iz1)
      endif

      nroot  = iz2
      call dzreal(z2,errabs,errrel,eps,eta,nroot,itmax,
     *            zeros2,exzeros2,infozer2)

      if(iout().gt.0) then
      write(6,'(a,f14.8,2x,f14.8,2x,i8,2x,f14.10)')
     *  (' even: ',exzeros2(i),exzeros2(i)*L/pi,infozer2(i),
     *    z2(exzeros2(i)),i=1,iz2)
      endif


      maxit    = 1000
      epsilon  = 1d-6
      gu       = -L/2
      go       =  L/2

! establish normalization constants
      do ip = 1,iz1   ! odd
       iq = ip
       norm1(ip) = 1.0d0
!       yi  = adapint(dot_kernel_odd,gu,go,epsilon,maxit,erraccu)
       yi = 1/(normin_odd(ip)**2)
       norm1(ip) = sqrt(1.0d0/yi)
       if(iout().gt.0) then
        write(6,*)'norm odd (',ip,')=',norm1(ip)
       endif
      enddo
      do ip = 1,iz2   ! even
       iq = ip
       norm2(ip) = 1.0d0
!       yi  = adapint(dot_kernel_ev,gu,go,epsilon,maxit,erraccu)
       yi = 1/(normin_even(ip)**2)
       norm2(ip) = sqrt(1.0d0/yi)
       if(iout().gt.0) then
        write(6,*)'norm even(',ip,')=',norm2(ip)
       endif
      enddo

      return
      end





      double precision function z1(alpha)
!     -----------------------------------> odd zeros

      include 'stiffzero_inc1.inc'
      double precision alpha


      z1 = alpha**3*sin(alpha*L/2)-sqrt(alpha**2+4*p**2)**3*cos(alpha*L/
     #2)*tanh(sqrt(alpha**2+4*p**2)*L/2)-2*p*(2*alpha**2+4*p**2)*cos(alp
     #ha*L/2)

      return
      end



      double precision function z2(alpha)
!     -----------------------------------> even zeros

      include 'stiffzero_inc1.inc'
      double precision alpha

      z2 = alpha**3*cos(alpha*L/2)+sqrt(alpha**2+4*p**2)**3*sin(alpha*L/
     #2)/tanh(sqrt(alpha**2+4*p**2)*L/2)+2*p*(2*alpha**2+4*p**2)*sin(alp
     #ha*L/2)

      return
      end


     
      double precision function psi_odd(m,s)
!     --------------------------------------
      include 'stiffzero_inc1.inc'
      integer  m
      double precision a,b, s
 
      if(m.lt.1 .or. m.gt. iz1) then
       write(6,*)'psi_odd m out of range:',m,iz1
      endif
      a = exzeros1(m)
      b = sqrt(a**2+4*p**2)

      psi_odd = a*sin(a*s)/cos(a*L/2)+b*sinh(b*s)/cosh(b*L/2)
      psi_odd = psi_odd*norm1(m)

      return
      end    

     
      double precision function psi_even(m,s)
!     --------------------------------------
      include 'stiffzero_inc1.inc'
      integer  m
      double precision a,b, s
 
      if(m.lt.1 .or. m.gt. iz2) then
       write(6,*)'psi_odd m out of range:',m,iz2
      endif
      a = exzeros2(m)
      b = sqrt(a**2+4*p**2)

      psi_even = -a*cos(a*s)/sin(a*L/2)+b*cosh(b*s)/sinh(b*L/2)
      psi_even = psi_even*norm2(m)

      return
      end    


      double precision functiondot_kernel_ev(s)
!     -----------------------------------------
      include 'stiffzero_inc1.inc'
      double precision psi_even, s

      dot_kernel_ev = psi_even(ip,s)*psi_even(iq,s)

      return
      end

      double precision function dot_kernel_odd(s)
!     -------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision psi_odd, s

      dot_kernel_odd = psi_odd(ip,s)*psi_odd(iq,s)

      return
      end

      double precision function dot_kernel_oe(s)
!     ------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision psi_even, psi_odd, s

      dot_kernel_oe = psi_odd(ip,s)*psi_even(iq,s)

      return
      end


      double precision function dot_kernel_ll(s)
!     ------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision psi_l,  s

      dot_kernel_ll = psi_l(ip,s)*psi_l(iq,s)

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Funktionen fuer die Normierung !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ACHTUNG: overflow Problem fuer hohe Modennummern
!          dat sinh*cosh sehr gross wird
!          z.B. in  sinhsinh = -(-2*cosh(L*b/2)*sinh(L*b/2)+L*b)/b/2


      double precision function sinsin(a)
!     -----------------------------------
      include 'stiffzero_inc1.inc'
      double precision a

      sinsin = (-2*cos(L*a/2)*sin(L*a/2)+L*a)/a/2

      return 
      end

      double precision function sinhsinh(b)
!     -------------------------------------
      include 'stiffzero_inc1.inc'
      double precision b

      sinhsinh = -(-2*cosh(L*b/2)*sinh(L*b/2)+L*b)/b/2

      return 
      end

      double precision function sinsinh(a,b)
!     --------------------------------------
      include 'stiffzero_inc1.inc'
      double precision a,b,t0

      t0 = (b*sin(L*a/2)*exp(-L*b)-a*cos(L*a/2)+a*cos(L*a/2)*exp(-L*b)+b
     #*sin(L*a/2))*exp(L*b/2)/(b**2+a**2)

      sinsinh = t0       

      return 
      end

      double precision function normin_odd(m)
!     -------------------------------------
      include 'stiffzero_inc1.inc'
      double precision a,b,am1,am2
      double precision sinsin, sinsinh, sinhsinh
      integer m

      a = exzeros1(m)
      b = sqrt(a**2+4*p**2)
      am1 = a/cos(a*L/2)
      am2 = b/cosh(b*L/2)

      normin_odd = am1*am1*sinsin(a) +
     *             2*am1*am2*sinsinh(a,b) +
     *             am2*am2*sinhsinh(b)

      normin_odd = 1.0d0/sqrt(normin_odd)

      return
      end


      double precision function coscos(a)
!     -----------------------------------
      include 'stiffzero_inc1.inc'
      double precision a

      coscos = (2*cos(L*a/2)*sin(L*a/2)+L*a)/a/2

      return 
      end

      double precision function coshcosh(b)
!     -------------------------------------
      include 'stiffzero_inc1.inc'
      double precision b

      coshcosh = (2*cosh(L*b/2)*sinh(L*b/2)+L*b)/b/2

      return 
      end

      double precision function coscosh(a,b)
!     --------------------------------------
      include 'stiffzero_inc1.inc'
      double precision a,b,t0

      t0 = (-b*cos(L*a/2)*exp(-L*b)+a*sin(L*a/2)+a*sin(L*a/2)*exp(-L*b)+
     #b*cos(L*a/2))*exp(L*b/2)/(b**2+a**2)

      coscosh = t0

      return 
      end

      double precision function normin_even(m)
!     ----------------------------------------
      include 'stiffzero_inc1.inc'
      double precision a,b,am1,am2
      double precision coscos, coscosh, coshcosh
      integer m

      a = exzeros2(m)
      b = sqrt(a**2+4*p**2)
      am1 = -a/sin(a*L/2)
      am2 = b/sinh(b*L/2)

      normin_even = am1*am1*coscos(a) +
     *             2*am1*am2*coscosh(a,b) +
     *             am2*am2*coshcosh(b)

      normin_even = 1.0d0/sqrt(normin_even)

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Relaxationsspektrum                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function tau_l(n)
!     ----------------------------------
      include 'stiffzero_inc1.inc'
      integer n, i
      double precision a

      i = n/2

      if(2*i.eq.n) then      ! --> EVEN
        a = exzeros2(i)
!!      tau_l = (zeta/l0)/(3*kb*Temp*((a**4)/(4*p)+p*(a**2)))
        tau_l =    (gamx)/(3*kb*Temp*((a**4)/(4*p)+p*(a**2)))
        return
      else                   ! --> ODD
        a = exzeros1(i+1)
!!      tau_l = (zeta/l0)/(3*kb*Temp*((a**4)/(4*p)+p*(a**2)))
        tau_l =    (gamx)/(3*kb*Temp*((a**4)/(4*p)+p*(a**2)))
        return
      endif 
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Durchnummerierte Eigenfunktionen                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function psi_l(n,s)
!     ------------------------------------
      include 'stiffzero_inc1.inc'
      integer n, i
      double precision s, psi_even, psi_odd

      i = n/2

      if(2*i.eq.n) then      ! --> EVEN
        psi_l = psi_even(i,s)
        return
      else                   ! --> ODD
        psi_l = psi_odd(i+1,s)
        return
      endif 
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Phi(s,s*,t)                                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function phi_sst(s1,s2,t)
!     ------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision sum, s1, s2, t 
      double precision psi_l, tau_l, tau
      double precision p1, p2
      integer i, n

      n = 2*(iz1-1)
      
      sum = 0
      
      do i=1,n
       tau = tau_l(i)
       p1  = psi_l(i,s1)
       p2  = psi_l(i,s2)
       sum = sum + tau*(p1**2+p2**2-2*p1*p2*exp(-t/tau))
      enddo

!!    phi_sst = 3*kb*Temp/((zeta/l0)) * sum
!! 
!!    phi_sst = phi_sst+diffsca*6*kb*Temp/(L/l0)/zeta*t  ! ??

      phi_sst = 3*kb*Temp/((gamx)) * sum
   
      phi_sst = phi_sst+diffsca*6*kb*Temp/nsegment/zeta*t  

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Phi(s,s*,t) Rouse                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function phi_sst_rouse(s1,s2,t)
!     -----------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision sum, s1, s2, t 
      double precision psi_l, tau_l, tau
      double precision p1, p2, al0
      double precision etau
      integer i, n, nn

      nn  = NINT(L*p)
      al0 = 1/p

      n = 2*(iz1-1) ! ---> max modenumber
    
!      if(n.gt.nn) n = nn
      
      sum = 0
      
      do i=1,n 
       etau = 3*kb*Temp/(zeta/l0)*((i*pi/L)**2/(Cinf*l0)+
     *                             (i*pi/L)**4*Cinf*l0/4)
       tau = 1.0d0/etau
       p1  = cos(i*pi*((s1)/al0)/nn)
       p2  = cos(i*pi*((s2)/al0)/nn)
       sum = sum + tau*(p1**2+p2**2-2*p1*p2*exp(-t/tau))
      enddo

      phi_sst_rouse = 3*kb*Temp/((zeta/l0)) * sum * 4.0d0/(nn**2)

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Twod-integration Kernel for S(q,t)                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function phi_sst_kernel(m,s)
!     ---------------------------------------------
      include 'stiffzero_inc1.inc'
      integer m
      double precision s(m), pss
      double precision phi_sst
      
      pss            = phi_sst(s(1),s(2),tim)       
      phi_sst_kernel = exp(-q**2/6*pss)
!      write(6,1)'pk:',q,tim,s,pss,phi_sst_kernel
! 1    format(a,6F16.8)  
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! inner Kernel for two step integration S(q,t)              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function phi_sst_kernel1(s1)
!     ---------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision s1, s2, phi_sst, pss
      common /cisst/s2
      
      pss            = phi_sst(s1,s2,tim)       
      phi_sst_kernel1 = exp(-q**2/6*pss)

      return
      end

      double precision function phi_sst_kernel2(s20)
!     ----------------------------------------------
      include 'stiffzero_inc1.inc'
      double precision s1, s2, s20
      double precision gu , go , eps, errac
      integer maxit
      common /cisst/s2
      double precision phi_sst_kernel1, adapint
      external         phi_sst_kernel1
   
      s2 = s20
      
      eps    = 1d-5
      maxit  = 200
      gu = -L/2
      go =  L/2   

      phi_sst_kernel2 = 
     *   adapint(phi_sst_kernel1,gu,go,eps,maxit,errac)/L

      return
      end


! And S(Q,t) by integration (adapint)

      double precision function sqtii(qq,t)
!     -------------------------------------
      include 'stiffzero_inc1.inc'
      double precision s1, s2, s20, qq, t
      double precision gu , go , eps, errac
      integer maxit
      double precision phi_sst_kernel2, bdapint
      external         phi_sst_kernel2

      q   = qq
      tim = t

      eps    = 1d-4
      maxit  = 200
      gu = -L/2
      go =  L/2   

      sqtii = 
     * bdapint(phi_sst_kernel2,gu,go,eps,maxit,errac)/L

      return
      end




! And S(Q,t) by integration (nag) 

      double precision function sqtnag(qq,t)
!     --------------------------------------
      include 'stiffzero_inc1.inc'
      double precision qq, t
      double precision a(2), b(2), errabs, errrel, result
      double precision eps, acc
      integer          minpts, maxpts, lenwrk, ifail
      integer          iout
      parameter(maxpts=20000, lenwrk=10000)
      double precision wrk(lenwrk), finval
      double precision phi_sst_kernel
      external         phi_sst_kernel



      q   = qq
      tim = t

      a(1) = -L/2
      a(2) = -L/2
      b(1) =  L/2
      b(2) =  L/2

      eps     = 1d-3
      minpts  = 200
      ifail   = -1

      call D01FCF(2,a,b,MINPTS,MAXPTS,phi_sst_kernel, 
     *            EPS, ACC, LENWRK, WRK, FINVAL, 
     *                      IFAIL) 


      sqtnag = finval/(L**2)

      if(iout().gt.0)then 
        write(6,'(a,e12.4,a,I5,a,I2,a,f12.7,a,f12.7,a,f12.7)')
     * ' acc=',acc,' n=',minpts,' ifail=',ifail,
     * ' satnag(',q,',',t,')=',sqtnag
      endif

      return
      end



! And S(Q,t) by summation

      double precision function sqts(qq,t)
!     ------------------------------------
      include 'stiffzero_inc1.inc'
      double precision qq, t
      double precision s(2), sum, ds
      integer i, j, nsum
      double precision phi_sst_kernel

      q   = qq
      tim = t
      
      nsum = NINT(L*p) 
      ds  = L/nsum
      sum  = 0

      do i = 1,nsum
       s(1) = i*ds-L/2-ds/2 
       do j = 1,nsum      
         s(2) = j*ds-L/2-ds/2
         sum = sum + phi_sst_kernel(2,s)
       enddo
      enddo

      sqts = sum/nsum**2

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
       function gaus8b (f,xu,xo)
c-----------------------------------------------------------------------
c      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
c-----------------------------------------------------------------------
       parameter ( ndim = 8 )
       parameter ( ndim2=ndim/2)
       implicit real*8 (a-h,o-z)
       dimension a(ndim2),x(ndim2)
       data a / 0.362683783378362d0,
     *          0.313706645877887d0,
     *          0.222381034453374d0,
     *          0.101228536290376d0/
       data x / 0.183434642495650d0,
     *          0.525532409916329d0,
     *          0.796666477413627d0,
     *          0.960289856497536d0/
 
       xave = (xo+xu)*0.5d0
       range= (xo-xu)*0.5d0
       sum = 0.d0
       do i=1,ndim2
         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) )
       enddo
       gaus8b = sum * range
 
       return
       end
c*ds
c*ed
       function bdapint(f,a,b,epsilon,maxiter,erroraccu)
c      =================================================                --
c
c      lokal adaptives integrationsverfahren 2te ver. wg. rekur.
c
       implicit real*8 (a-h,o-z)
c
       parameter( maxstack = 200 )
c                 --------------->  stacktiefe
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack)
c
       logical cray
!!     common/outlev/iot,ibild,ierrs,inka, cray
       real*4 xxxx,yyyy,ptxf                                         !! aix
       common/outlev/iot ,ibild,ierrs,inka,iibuf,xxxx,yyyy,ptxf(20)  !! aix
c
       external f
c
c      -------------------------------------------------------------------
 
 
         iterationcounter = 0
         erg              = 0.0d0
         itop      = 0
         xa(itop)  = a
         xb(itop)  = b
         s(itop)   = gaus8b( f, a, b )
 
         do i=1,maxiter
           if (itop.ge.(maxstack-2)) then
              erg = erg + s(itop)
              xbb = xb(itop)
              itm = itop-1
              do itop=itm,0,-1
                 if(xa(itop).eq.xbb) goto 1
              enddo
1             continue
              write(6,*)'warning! adaptint stack overflow!'
           else
              iterationcounter = iterationcounter + 1
              itop             = itop +1
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0
              xb(itop) =  xb(itop-1)
              s(itop)  = gaus8b( f, xa(itop), xb(itop) )
              itop             = itop +1
              xa(itop) =  xa(itop-2)
              xb(itop) =  xa(itop-1)
              s(itop)  = gaus8b( f, xa(itop), xb(itop) )
              error   =  dabs((s(itop)+s(itop-1))-s(itop-2))
 
              if (iot.gt.2) then
                 write(6,'(1x,i3,i5,4e13.6)')itop,iterationcounter,
     *                 xa(itop),xb(itop),(s(itop)+s(itop-1)),error
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
2                continue
              endif
            endif
            if (itop.le.0) goto 3
         enddo
         write(6,*)                                                     '
     *    'bdapint fatal error iterationnumber exceeded!'
3        continue
 
 
         bdapint = erg
 
         return
 
       end
