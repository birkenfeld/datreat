         function th22(x,pa,thnam,parnam,npar,ini)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ========================================
!
!                          -------> Debye <--------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         DOUBLE PRECISION Pi, Navog
         Parameter       (Pi=3.141592654d0)
         Parameter       (Navog=6.022045d23)


         CHARACTER*8     thnam,parnam(20)
         REAL            th22, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision v, nphi

         DOUBLE PRECISION ampli, q, rg, qrg, gamma, sq, sqi
         common/cdebik/gamma,qrg

         DOUBLE PRECISION mw, bpoly, drho, conc, rhopoly,fac
         real bsolv 


         DOUBLE PRECISION dbik, adapint, erra
         external dbik

 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget
c
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'debye2  '
         nparx = 8
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           th22   = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters
         parnam(1) = 'intensit'          
         parnam(2) = 'rg      '
         parnam(3) = 'gamma   '
         parnam(4) = 'molwght '
         parnam(5) = 'density '
         parnam(6) = 'bpolym  '
         parnam(7) = 'v       '
         parnam(8) = 'volfrac '
c
         th22  = 0
         return
       endif
c
c ---- calculate theory here -----
        ampli   = abs(pa(1))  
        rg      = abs(pa(2))
        gamma   = abs(pa(3))
        mw      = abs(pa(4))
        rhopoly = abs(pa(5))
        bpoly   = pa(6)
        v       = pa(7)
        conc    = abs(pa(8)) 



        call parget('bsolv   ',bsolv,iadda,ier)
        drho = bpoly-bsolv

        fac = conc * mw / rhopoly / Navog 

        q       = x
        qrg     = q*q*rg*rg

        if(gamma.eq.0.0d0) then
          sq =(2/(qrg**2))*(qrg-1+exp(-qrg))
        else
          sq = adapint(dbik,0.0d0,1.0d0,1d-8,1000,erra)
        endif

! Zimm formula for interacting chains
        sqi = 1.0d0/(fac*sq) + v*conc*conc
        if(abs(sqi).lt.1d-20) sqi = 1d-20

        th22 = ampli*drho**2/sqi

        return
        end 




