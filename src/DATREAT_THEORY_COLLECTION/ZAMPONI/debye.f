         function th22(x,pa,thnam,parnam,npar,ini)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ========================================
!
!          -------> Debye <--------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         DOUBLE PRECISION Pi, Navog
         Parameter       (Pi=3.141592654d0)
         Parameter       (Navog=6.022045d23)


         CHARACTER*8     thnam,parnam(20)
         REAL            th22, x, pa, qq, zpi, xh, vol_frac,bk
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision v, nphi

         DOUBLE PRECISION ampli, q, rg, qrg, gamma, sq
         common/cdebik/gamma,qrg

         DOUBLE PRECISION dbik, adapint, erra      
         external dbik

         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget
c
c
c ----- initialisation -----
       if(ini.eq.0) then

         thnam = 'debye   '

         nparx = 3

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
         parnam(3) = 'bk      '
c
         th22  = 0
         return
       endif
c
c ---- calculate theory here -----
        ampli   = abs(pa(1))  
        rg      = abs(pa(2))
        bk      = abs(pa(3))

        q       = x

        qrg     = q*q*rg*rg

        sq =(2/(qrg**2))*(qrg-1+exp(-qrg))

        th22 = ampli*sq + bk

        return
        end 
