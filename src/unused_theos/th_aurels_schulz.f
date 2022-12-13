           function check(x, pa, thnam, parnam, npar,idum,ini)


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)


         CHARACTER*8     thnam,parnam(20)
         REAL            check, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx,idum
         DIMENSION       pa(20),qq(3)
        DATA            zpi/6.283185/

 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

         DOUBLE PRECISION pepep_1, Q,ds
     *                    rhosolv,  rhopep, rhope, 
     *                    d, sigc, al, sigb, R1, R2, c,
     *                    a, dl,i0, contrast,  
     *                    alpha, xi, anu,
     *                    porod
     *                    aint,conc


         DOUBLE PRECISION r0, rsigma, r_min,ds,cconc
         INTEGER          n 

         INTEGER          ns
         DOUBLE PRECISION xunit,b_w,b_s
         Parameter       (xunit = 1d-8)      !! cm pro Angstroem 
         DOUBLE PRECISION FFcs,FFss
         DOUBLE PRECISION con_cs,ampl_cs,con_ss,ampl_ss
         DOUBLE PRECISION d_C12E5,SLD_D20,SLD_C12E5,r_total
         DOUBLE PRECISION cQ, crhosolv,  crhopep, crhope,cd, csigc 
         DOUBLE PRECISION cal, csigb, cc,ca, cdl,calpha, cxi, canu 
         DOUBLE PRECISION  ss,fc,xx,fs
         DOUBLE PRECISION cSLD_D20, cSLD_C12E5,cd_C12E5
         DOUBLE PRECISION cporod, cr0, crsigma, ci0, ccontrast
         INTEGER          cns, cn
         COMMON/ker10r/ cQ, crhosolv,  crhopep, crhope,cd, csigc,
     *                 cal, csigb, cc,ca, cdl,calpha, cxi, canu,
     *                 cporod, cr0, crsigma, cns, cn, ci0, ccontrast,
     *                 cSLD_D20, cSLD_C12E5,cd_C12E5,cconc
        
         DOUBLE PRECISION adapint, epsilon, erraccu, gu , go
         INTEGER          maxit 

         DOUBLE PRECISION  pepep_1_kernel_r,norm,wsize_kernel_r
         EXTERNAL          pepep_1_kernel_r,wsize_kernel_r
c
c
c ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'check'
         nparx = 7
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,
     *      ' exceeds current max. = ',i8)
           check   = 0
           return
         endif
         npar = nparx
c        --------------> set the number of parameters                 

         parnam(1)='d_C12E5' 
         parnam(2)='conc' 
         parnam(3)='SLD_D20'
         parnam(4)='SLD_C12E5'
         parnam(5)='z_schulz'
         parnam(6)='r_rod'
         parnam(7)='r_min'





         check  = 0
         return
       endif
c
c ---- calculate theory here -----
 
        d_C12E5 = pa(1)
           conc = pa(2) 
        SLD_D20 = pa(3) 
      SLD_C12E5 = pa(4)
         rsigma = pa(5)
             r0 = pa(6)  
          r_min = pa(7)

 

        Q = x / xunit

        cQ        = Q
        

        
        cr0       = r0 
        crsigma   = rsigma
        cSLD_D20  = SLD_D20
        cSLD_C12E5 = SLD_C12E5
        cd_C12E5 = d_C12E5                
        cconc=conc      
        


        epsilon = 1d-7
        maxit   = 100
        erraccu = 0.0
       
        gu  = r_min
        go  = 30*r0
        if(gu.lt.0.0d0) gu = 0.0d0 
!-------------------------------------------------------------------------------


           Q = x / xunit



         check = adapint(wsize_kernel_r,gu,go,epsilon,maxit,erraccu)

         norm = adapint(wsize_kernel_r,gu,go,epsilon,maxit,erraccu)
        
          check = check/norm


!              check=pepep_1(Q,SLD_D20,SLD_C12E5,d_C12E5,conc,r0)
        
        
        return
        end function check
!----------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         function   pepep_1(qi,                        !! momentum transfer 
     *                      b_w,                   !! scattering length densities in cm**-2
     *                      b_s,ds,cc,R)
                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         implicit none
         
         DOUBLE PRECISION pepep_1
         DOUBLE PRECISION Pi,cc
         Parameter       (Pi=3.141592654d0)
         DOUBLE PRECISION 
     *                    rhosolv,  rhopep, rhope, 
     *                    ds, sigc, L, sigb, R, c,
     *                    a, dl, contrast, i0, 
     *                    alpha, xi, anu,
     *                    porod, dbesj1

         INTEGER          ns

         DOUBLE PRECISION qi,xx,xunit
         INTEGER          js
         DOUBLE PRECISION blob, alpha0, qxi, qne
         DOUBLE PRECISION DDAWS, r1, r2
         DOUBLE PRECISION r_total, x, Q,cs,ss
         DOUBLE PRECISION Fc,Fs,b_w,b_s
         Parameter       (xunit = 1d-8) 
          
          
          Q=qi                                      

         if (abs(qi) .le. 1.d-7) then
          write(6,*)'Warning qi=',qi,' to small, force to 1e-7' 
          qi = 1.d-7
         endif        


! core shell model:

              r_total = r + ds
        
              xx=Q*r_total   
              
             ss = (b_s - b_w)**2

          
           Fc=(r**3)*(3*(sin(Q*r)-Q*r*cos(Q*r))/((Q*r)**3))

            

            Fs=(r_total**3)*(3*(sin(xx)-xx*cos(xx))/(xx**3))
             
              
              
             pepep_1= cc*ss*((4*pi/3)**2)*(Fc + Fs)**2



         return
         end


         Function pepep_1_kernel_r(r)
!!       ----------------------------
!!
         implicit none
         
         DOUBLE PRECISION pepep_1_kernel_r, r, pepep_1
         DOUBLE PRECISION size_distribution2,c,pschulz
         DOUBLE PRECISION SLD_D20, SLD_C12E5,d_C12E5
         
         DOUBLE PRECISION cSLD_D20, cSLD_C12E5,cd_C12E5
         DOUBLE PRECISION cconc
         DOUBLE PRECISION cQ, crhosolv,  crhopep, crhope,cd, csigc 
         DOUBLE PRECISION cal, csigb, cc,ca, cdl,calpha, cxi, canu 
         DOUBLE PRECISION cporod, cr0, crsigma, ci0, ccontrast
         INTEGER          cns, cn
         COMMON/ker10r/ cQ, crhosolv,  crhopep, crhope,cd, csigc,
     *                 cal, csigb, cc,ca, cdl,calpha, cxi, canu,
     *                 cporod, cr0, crsigma, cns, cn, ci0, ccontrast,
     *                 cSLD_D20, cSLD_C12E5,cd_C12E5
     *                  

        pepep_1_kernel_r =  pepep_1(cQ,cSLD_D20,cSLD_C12E5,
     *                       cd_c12E5,cconc,r)                   


        pepep_1_kernel_r =  pepep_1_kernel_r * 
     *                      size_distribution2(r, cr0, crsigma)
!                              
         
        return
        end


         Function wsize_kernel_r(r)
!!       ----------------------------
!!
         implicit none
         
         DOUBLE PRECISION wsize_kernel_r, r
         DOUBLE PRECISION size_distribution2,pschulz 

 
         DOUBLE PRECISION cd_C12E5,cSLD_D20,cSLD_C12E5
         DOUBLE PRECISION cQ, crhosolv,  crhopep, crhope,cd, csigc 
         DOUBLE PRECISION cal, csigb, cc,ca, cdl,calpha, cxi, canu 
         DOUBLE PRECISION cporod, cr0, crsigma, ccontrast, ci0
         INTEGER          cns, cn
         COMMON/ker10r/ cQ, crhosolv,  crhopep, crhope,cd, csigc,
     *                 cal, csigb, cc,ca, cdl,calpha, cxi, canu,
     *                 cporod, cr0, crsigma, cns, cn, ci0, ccontrast,
     *                 cSLD_D20, cSLD_C12E5,cd_C12E5
        
         wsize_kernel_r = size_distribution2(r, cr0, crsigma)

        

              return
        end


             Function size_distribution2(r, r0, sigma)
!!      -------------------------------------------
!!
!! non-normalized !!!!
!! here: Schulz-Flory, r0==rquer, sigma==Z
!!
        implicit none

        DOUBLE PRECISION size_distribution2
        DOUBLE PRECISION arg, r, r0, sigma, dgamma
        DOUBLE PRECISION Z, rquer
        INTEGER          n
    
        Z = abs(sigma)
        if(Z.eq.0.0d0)   Z = 1d-4  

        rquer = abs(r0)
        if(r0.eq.0.0d0)  r0 = 1d-4

        arg = - (Z+1)*r/rquer
        if(arg.lt.-50.0d0) arg = -50.0d0

 
        size_distribution2=(r**Z)*exp(arg)*((Z+1)/r0)**(z+1)/dgamma(z+1) 

        return
        end
