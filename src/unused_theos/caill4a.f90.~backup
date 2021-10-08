        subroutine setcai(ceps,ceps2,cxmuswit,cstrechf,crellim,         &
     & csplitf4,crres,ie1sw,ipolpo)                                     
! ----===================---------------------------------------------- 
!  global parameters that determine the performance of the caille       
!  system may be set from outside                                       
! ----------------------------------------------------------------------
!                                                                       
       implicit real*8(a-h,o-z) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
                                                                        
       if(ceps    .ne. 0.0d0 ) eps       = ceps 
       if(ceps2   .ne. 0.0d0 ) eps2      = ceps2 
       if(cxmuswit.ne. 0.0d0 ) xmuswit   = cxmuswit 
       if(cstrechf.ne. 0.0d0 ) strechf   = cstrechf 
       if(crellim .ne. 0.0d0 ) rellim    = crellim 
       if(csplitf4.ne. 0.0d0 ) splitf4   = csplitf4 
       if(crres   .ne. 0.0d0 ) rres      = crres 
       if(ie1sw   .ge. 0     ) ie1swit   = ie1sw 
       if(ipolpo  .ne. 0     ) nppow     = ipolpo 
                                                                        
       if(iout().gt.1) then 
         write(6,*)'---------- setcai -------------------------' 
         write(6,*)'eps       = ',eps 
         write(6,*)'eps2      = ',eps2 
         write(6,*)'xmuswit   = ',xmuswit 
         write(6,*)'xupper    = ',xupper 
         write(6,*)'rres      = ',rres 
         write(6,*)'strechf   = ',strechf 
         write(6,*)'rellim    = ',rellim 
         write(6,*)'wrnlim    = ',wrnlim 
         write(6,*)'rlim4     = ',rlim4 
         write(6,*)'splitf4   = ',splitf4 
         write(6,*)'maxit     = ',maxit 
         write(6,*)'maxit2    = ',maxit2 
         write(6,*)'ie1swit   = ',ie1swit 
         write(6,*)'nppow     = ',nppow 
         write(6,*)'inf4      = ',inf4 
         write(6,*)'nstart4   = ',nstart4 
         write(6,*)'-------------------------------------------' 
       endif 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       block data caillbd 
! ----===================---------------------------------------------- 
!  global parameters that determine the performance of the caille       
!  system                                                               
! ----------------------------------------------------------------------
!                                                                       
       implicit real*8(a-h,o-z) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
! ---->                                                                 
!      eps ...... adaptive integration accuracy for caille1 (x<xmuswit) 
!      eps2 ..... series term limit for caille1 (x>xmuswit, eps of xnuej
!      xmuswit .. xmue value separating the gamma-fktn+adapint approx. a
!                 the series (xnuej0) approx. in caille2, must be less t
!      xupper ... upper integral limit for the series expansion in caill
!      rres ..... real space resolution parameter for the caillrx type r
!      strechf .. strech factor to create subsequent intervals of       
!                 caille2 approximation in caill4a                      
!      rellim ... relative limit size of caille2 values in table generat
!                 of caille4a                                           
!x     abslim ... absolute limit size of caille2 values in table generat
!                 of caille3                                            
!      wrnlim ... accuracy waring limit in caille4a                     
!x     rlim4  ... accuracy limit for intervalsplit in caille4           
!      splitf4 .. interval split-factor for the first interval in caill4
!      maxit .... maximum adapint iterations in caille1/2 (xmu<xmuswit) 
!      maxit2 ... maximum series terms in caille1/2       (xmu>xmuswit) 
!      ie1swit... switch for the e1(q0*z*sqrt(alpha)) term in the dwf (c
!      nppow .... polynomial approximation power in caill4a             
!      inf4 ..... caille4a summation to infinity ?                      
!                                                                       
!x <-- items with x are obsolete                                        
!                                                                       
       data eps          /1.0d-6/ 
       data eps2         /1.0d-15/ 
       data xmuswit      /1.5d0/ 
       data xupper       /32.2013246993d0/ 
       data rres         /10.0d0/ 
       data strechf      /6.0d0/ 
       data rellim       /2.0d-4/ 
       data abslim       /0.0d0/ 
       data wrnlim       /5.0d-2/ 
       data rlim4        /0.0d0/ 
       data splitf4      /1.0d0/ 
       data maxit        /50/ 
       data maxit2       /50/ 
       data ie1swit      /0/ 
       data nppow        /5/ 
       data inf4         /1/ 
       data nstart4      /1/ 
                                                                        
      END                                           
                                                                        
                                                                        
                                                                        
       function ccheck2(qortho, qz, a, alpha, z, ztest) 
!      ------------------------------------------------                 
       implicit real*8 (a-h,o-z) 
       parameter    (q0=300.0d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
       save aa, am 
                                                                        
       if(z.lt.1.0d0) then 
         z1 = caille2(qortho, qz, q0, a, alpha, ztest) 
         z2 = caille2(qortho, qz, q0, a, alpha, ztest+1.0d0) 
         aa = dlog(z2)-dlog(z1) 
         am = z1/dexp(aa*ztest) 
         write(6,*)'amplitude =',am 
         write(6,*)'decay     =',aa 
       endif 
       zz = caille2(qortho,qz,q0,a,alpha,z) 
       ccheck2 = zz/(am*dexp(aa*z)) 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function caill4a (qortho, qz, q0, d , a, alpha) 
!------==============================================------------------ 
!                                                                       
! z-fouriersum of the caille theory by tabulation, polynomial*exp       
! expansion and summation of the polynomial terms by explicit           
! formulae                                                              
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       parameter    (mppow=5) 
       parameter    (mtable=1000) 
       parameter    (minterv=20) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
       dimension pcoeff(0:mppow) 
       dimension ztable(mtable), ctable(mtable), ratio(mtable) 
                                                                        
       if(nppow.gt.mppow) then 
         nppow = mppow 
         write(6,*)'caille4a: warning nppow reduced to mppow!' 
       endif 
                                                                        
! -- explicit summation of the first terms up to nstart - 1 --          
       zeps    = 0.02d0*d 
       sum     = caille2(qortho,qz,q0,a,alpha,zeps) 
                                                                        
! --- determine required z-range ---                                    
          z1 = d 
          z2 = z1+d 
          c1 = caille2(qortho, qz, q0, a, alpha, z1) 
          c2 = caille2(qortho, qz, q0, a, alpha, z2) 
          if(c1.le.0.d0) then 
             caill4a = sum 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caill4a( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)' ((c1)) )       =', caill4a 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
             return 
          endif 
          if(c2.le.0.d0) then 
             an  = d * qz 
             c   = dcos(an) 
             sum = sum + 2.0d0 * c1 * c 
             caill4a = sum 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caill4a( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)' ((c2)) )       =', caill4a 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
             return 
          endif 
          dz = z2-z1 
          decay0 = (dlog(c2)-dlog(c1))/dz 
          if(iout().gt.2) write(6,*)'decay-0 =',decay0 
! --- estimation of range (0) ---                                       
          zmax = dlog(rellim)/decay0 
          if(iout().gt.2) write(6,*)'range first estimate=',zmax 
! --- iterate range estimate ---                                        
          z1 = zmax/2 
          z2 = z1+d 
          c1 = caille2(qortho, qz, q0, a, alpha, z1) 
          c2 = caille2(qortho, qz, q0, a, alpha, z2) 
          dz = z2-z1 
          if(c1.le.0.0d0 .or. c2.le.0.0d0) then 
            decay1 = decay0 
          else 
            decay1= (dlog(c2)-dlog(c1))/dz 
          endif 
          if(iout().gt.2) write(6,*)'decay-1 =',decay1 
! --- estimation of range (1) ---                                       
          zmax = dlog(rellim)/decay1 
          if(iout().gt.2) write(6,*)'range 2nd   estimate=',zmax 
!                                                                       
                                                                        
! -- decide whether direct summation is more effective --               
       nsum = zmax/d 
       nsum = nsum*1.5d0 
       if(nsum.lt.40) then 
         if(iout().gt.2) write(6,*)'direct sum ',nsum, zmax 
         do n=1,nsum 
            z   = n * d 
            an  = z * qz 
            c   = dcos(an) 
            sum = sum + 2.0d0 * caille2(qortho,qz,q0,a,alpha,z) * c 
         enddo 
         caill4a = sum 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caill4a( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)' ((dsum))       =', caill4a 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
         return 
       endif 
                                                                        
       if(nstart4.gt.1) then 
         do n=1,nstart4-1 
            z   = n * d 
            an  = z * qz 
            c   = dcos(an) 
            sum = sum + 2.0d0 * caille2(qortho,qz,q0,a,alpha,z) * c 
         enddo 
       endif 
                                                                        
                                                                        
! --- determine splitting of the table ---                              
! ---  first interval:                                                  
       zint1 = (-1.0d0/decay0) * splitf4 
       if(iout().gt.2) write(6,*)'zint1 = ',zint1 
                                                                        
       rerr   = 0.0d0 
! --- number of points/interval ---                                     
!cx    ilength=   nppow * 3                                             
       ilength=   nppow * 2 
                                                                        
       n = 0 
       zi1 = d * nstart4 
       n1  = nstart4 
       zi2 = zint1 
       n2  = zint1/d 
       if(n1.ge.n2) then 
          write(6,*)'error: n1.ge.n2' ,n1,n2 
          caill4a = sum 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caill4a( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)' ((n1)) )       =', caill4a 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
          return 
       endif 
                                                                        
       do interv=1,minterv 
         if(iout().gt.2) then 
           write(6,*)'interval:',interv,' --> ',n1,n2 
           write(6,*)'interval:',interv,' --> ',zi1,zi2 
         endif 
         dzi = (zi2-zi1)/(ilength-1) 
                                                                        
         do j=1,ilength 
           z = zi1 + (j-1)*dzi 
           ca = caille2(qortho, qz, q0, a, alpha, z ) 
           ctable(j) = ca 
           ztable(j) = z 
           if(iout().gt.2) write(6,*)'ctable:',j,z,ca 
         enddo 
                                                                        
          c1 = ctable(1  ) 
          c2 = ctable(2  ) 
          z1 = ztable(1  ) 
          z2 = ztable(2  ) 
          dz = z2-z1 
          decay = (dlog(c2)-dlog(c1))/dz 
          if(iout().gt.2) write(6,*)'decay(',interv,') =',decay 
                                                                        
! --- prepare dataset for polyfit ---                                   
          do i=1,ilength 
            ratio(i) = ctable(i)*dexp(-decay*ztable(i))/ctable(1) 
            if(iout().gt.2) write(6,*) 'ratio =',i,ratio(i) 
          enddo 
!                                                                       
! --- and do the polyfit ...                                            
          call polyfit(ztable,ratio,ilength,pcoeff,nppow,r) 
          do i=0,nppow 
            pcoeff(i) = pcoeff(i) * ctable(1) 
            if(iout().gt.2) write(6,*)'pcoeff(',i,')= ',pcoeff(i) 
          enddo 
          if(r.gt.rerr) rerr = r 
          if(iout().gt.2) write(6,*)'r = ',r 
                                                                        
          dpow   = 1.0d0 
          do j = 0,nppow 
            sum = sum + pcoeff(j)*  eisum(qz,d,-decay,j,n1,n2) * dpow 
            dpow = dpow * d 
          enddo 
                                                                        
! --- prepare the next interval ---                                     
         if(zi2.ge.zmax) goto 2101 
         dz = (zi2 - zi1)*strechf 
         n1 = n2+1 
         zi1= d*n1 
         zi2= zi1+dz 
         if(zi2.ge.zmax) zi2 = zmax 
         if(zi2.le.zi1 ) then 
           zi2 = zi1 
           goto 2101 
         endif 
         n2 = zi2/d 
         zi2= n2*d 
                                                                        
        enddo 
        write(6,*)'too few intervals' 
 2101   continue 
!                                                                       
                                                                        
          if(inf4.ne.0) then 
!                       -------> sum to infinity                        
            n1  = n2 + 1 
            n2  = -1 
            dpow   = 1.0d0 
            do j = 0,nppow 
              sum = sum + pcoeff(j)*  eisum(qz,d,-decay,j,n1,n2) * dpow 
              dpow = dpow * d 
            enddo 
          endif 
                                                                        
                                                                        
       caill4a = sum 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caill4a( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'        )       =', caill4a 
        write(6,*)'-----------------==============================' 
        write(6,*)'   rerr = ',rerr 
        write(6,*)'-----------------------------------------------' 
       endif 
                                                                        
       if(rerr.gt.wrnlim) then 
         write(6,*)'caille4a:warning wrnlim<rerr=',rerr 
         if(iout().gt.2) then 
           write(6,*)'qortho   =',qortho 
           write(6,*)'qz       =',qz 
           write(6,*)'q0       =',q0 
           write(6,*)'d        =',d 
           write(6,*)'a        =',a 
           write(6,*)'alpha    =',alpha 
           write(6,*)'caille4a =',sum 
         endif 
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
       subroutine caillea(qortho,qz,q0,a,alpha,decay,pcoeff,zmin,zmax,r) 
!------===============================================================--
!                                                                       
! approximation of caille2 values by exp(-decay*z)*polynom(z)           
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       parameter    (mtable=200) 
       parameter    (mppow=5) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
       dimension ztable(mtable), ctable(mtable), ratio(mtable) 
       dimension pcoeff(0:mppow) 
                                                                        
       if(nppow.gt.mppow) then 
         nppow = mppow 
         write(6,*)'caillea: warning nppow reduced to mppow!' 
       endif 
                                                                        
! --- create table ----                                                 
       z = zmin 
       do i=1,mtable 
          ca = caille2(qortho, qz, q0, a, alpha, z) 
          ztable(i) = z 
          ctable(i) = ca 
          if(iout().gt.3) write(6,*)'table(',i,'): ',z,ca 
          if(i.gt.nppow+1) then 
            if(z            .gt. zmax  ) goto 1 
            if(ca/ctable(1) .lt. rellim) goto 1 
            if(ca           .lt. abslim) goto 1 
          endif 
          z  = z * strechf 
       enddo 
                                                                        
       write(6,*)'caillea: mtable=',mtable,                             &
     &           ' is not sufficient with strechf=',strechf             
       write(6,*)'qortho   =',qortho 
       write(6,*)'qz       =',qz 
       write(6,*)'q0       =',q0 
       write(6,*)'a        =',a 
       write(6,*)'alpha    =',alpha 
                                                                        
                                                                        
    1  continue 
       ntable = i 
       zmax   = ztable(ntable) 
!                                                                       
! --- determine the exp ---                                             
!cc      c1 = ctable(1)                                                 
!cc      c2 = ctable(2)                                                 
!cc      z1 = ztable(1)                                                 
!cc      z2 = ztable(2)                                                 
         c1 = ctable(ntable/2) 
         c2 = ctable(ntable/2+1) 
         z1 = ztable(ntable/2) 
         z2 = ztable(ntable/2+1) 
!ccc     c1 = ctable(ntable/2)                                          
!ccc     c2 = ctable(ntable)                                            
!ccc     z1 = ztable(ntable/2)                                          
!ccc     z2 = ztable(ntable)                                            
         dz = z2-z1 
         decay = (dlog(c2)-dlog(c1))/dz 
         if(iout().gt.2) write(6,*)'decay     =',decay 
                                                                        
! --- prepare dataset for polyfit ---                                   
       do i=1,ntable 
         ratio(i) = ctable(i)*dexp(-decay*ztable(i))/ctable(1) 
         if(iout().gt.2) write(6,*) 'ratio =',i,ratio(i) 
       enddo 
!                                                                       
! --- and do the polyfit ...                                            
       call polyfit(ztable,ratio,ntable,pcoeff,nppow,r) 
                                                                        
       do i=0,nppow 
        pcoeff(i) = pcoeff(i) * ctable(1) 
       enddo 
                                                                        
!c     write(6,*)'caillea: r=',r                                        
                                                                        
       return 
      END                                           
                                                                        
       function caille4(qortho, qz, q0, d , a, alpha) 
!------==============================================------------------ 
!                                                                       
! z-fouriersum of the caille theory by tabulation, polynomial*exp       
! expansion and summation of the polynomial terms by explicit           
! formulae                                                              
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       parameter    (mppow=5) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
       dimension pcoeff(0:mppow) 
                                                                        
       sum     = 0.0d0 
       zm      = 9.9999d33 
       zmin    = d * nstart4 
                                                                        
       if(iout().gt.2) then 
         write(6,*)'caille4 called with parameters:' 
         write(6,*)'qortho   =',qortho 
         write(6,*)'qz       =',qz 
         write(6,*)'q0       =',q0 
         write(6,*)'d        =',d 
         write(6,*)'a        =',a 
         write(6,*)'alpha    =',alpha 
       endif 
                                                                        
       call caillea(qortho,qz,q0,a,alpha,deca ,pcoeff,zmin,zm,r) 
       zm1 = zm 
       rerr = r 
                                                                        
       if(iout().gt.2) then 
           write(6,*)'---caille4 --------- 1 ----------------------' 
           write(6,*)' decay = ',deca 
           write(6,*)' zmax  = ',zm 
           write(6,*)' r     = ',r 
           do i=0,nppow 
             write(6,*)' pcoeff(',i,')= ',pcoeff(i) 
           enddo 
           write(6,*)'---------------------------------------------' 
       endif 
                                                                        
       if( r .gt. rlim4 ) then 
          zmi = zmin 
          zma = zm1 * splitf4 
          call caillea(qortho,qz,q0,a,alpha,deca ,pcoeff,zmi,zma,r) 
          rerr = r 
          if(iout().gt.2) then 
             write(6,*)'---caille4---------- 2a ---------------------' 
             write(6,*)' decay = ',deca 
             write(6,*)' zmax  = ',zma 
             write(6,*)' r     = ',r 
             do i=0,nppow 
               write(6,*)' pcoeff(',i,')= ',pcoeff(i) 
             enddo 
             write(6,*)'---------------------------------------------' 
          endif 
                                                                        
          n1     = nstart4 
          n2     = zma/d 
          zma    = d*n2 
          dpow   = 1.0d0 
          do n = 0,nppow 
            sum = sum + pcoeff(n) * eisum(qz,d,-deca,n,n1,n2) * dpow 
            dpow = dpow * d 
          enddo 
                                                                        
          zmi = zma 
          zma = 9.99999d33 
          call caillea(qortho,qz,q0,a,alpha,deca ,pcoeff,zmi,zma,r) 
          if(r.gt.rerr) rerr = r 
                                                                        
          if(iout().gt.2) then 
             write(6,*)'---caille4---------- 2b ---------------------' 
             write(6,*)' decay = ',deca 
             write(6,*)' zmax  = ',zma 
             write(6,*)' r     = ',r 
             do i=0,nppow 
               write(6,*)' pcoeff(',i,')= ',pcoeff(i) 
             enddo 
             write(6,*)'---------------------------------------------' 
          endif 
                                                                        
          n1     = n2 + 1 
          n2     = zma/d 
          zma    = d*n2 
          if(inf4.ne.0) n2 = -1 
!                       -------> sum to infinity                        
          dpow   = 1.0d0 
          do n = 0,nppow 
            sum = sum + pcoeff(n) * eisum(qz,d,-deca,n,n1,n2) * dpow 
            dpow = dpow * d 
          enddo 
                                                                        
       else 
                                                                        
          n1     = nstart4 
          n2     = zm1/d 
          zma    = d*n2 
          if(inf4.ne.0) n2 = -1 
!                       -------> sum to infinity                        
          dpow   = 1.0d0 
          do n = 0,nppow 
            sum = sum + pcoeff(n)*  eisum(qz,d,-deca,n,n1,n2) * dpow 
            dpow = dpow * d 
          enddo 
                                                                        
       endif 
                                                                        
                                                                        
! --- add the rest (beginning of sum.. )                                
                                                                        
       if(nstart4.gt.1) then 
       do n=1,nstart4-1 
            z   = n * d 
            an  = z * qz 
            c   = dcos(an) 
            ac  =  caille2(qortho,qz,q0,a,alpha,z) 
            sum = sum + 2.0d0 * caille2(qortho,qz,q0,a,alpha,z) * c 
         enddo 
       endif 
                                                                        
                                                                        
       zeps    = 0.02d0*d 
       sum     = sum + caille2(qortho,qz,q0,a,alpha,zeps) 
                                                                        
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caille4( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'        )       =', caille4 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
       caille4 = sum 
                                                                        
       if(rerr.gt.wrnlim) then 
         write(6,*)'caille4:warning wrnlim<rerr=',rerr 
         if(iout().gt.2) then 
           write(6,*)'qortho  =',qortho 
           write(6,*)'qz      =',qz 
           write(6,*)'q0      =',q0 
           write(6,*)'d       =',d 
           write(6,*)'a       =',a 
           write(6,*)'alpha   =',alpha 
           write(6,*)'caille4 =',sum 
         endif 
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function caille3(qortho, qz, q0, d , a, alpha, nmax) 
!------====================================================------------ 
!                                                                       
! z-fouriersum of the caille theory                                     
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       zeps = 0.01d0 
                                                                        
       sum = 0.0d0 
                                                                        
       do n=1,nmax 
          z   = n * d 
          an  = z * qz 
          c   = dcos(an) 
          ac  =  caille2(qortho,qz,q0,a,alpha,z) 
          sum = sum + ac * c 
          if(ac.lt.eps) goto 100 
       enddo 
  100  continue 
                                                                        
       sum = 2.0d0*sum + caille2(qortho,qz,q0,a,alpha,zeps) 
                                                                        
       caille3 = sum 
       if(iout().gt.0) then 
        write(6,*)'caille3( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'        )       =', caille3 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
!c-test-                                                                
!c     write(6,*)'n=',n,' ac=',ac,'caille3=',sum                        
!c-teste-                                                               
       return 
      END                                           
                                                                        
                                                                        
       function caillr3(qortho, qz, q0, d , a, alpha, nmax) 
!------====================================================------------ 
!                                                                       
! z-fouriersum of the caille theory  (with resolution)                  
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       parameter    (sqrtpi=1.772453850906d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       zeps = 0.01d0 
                                                                        
       sum = 0.0d0 
       sum2= 0.0d0 
                                                                        
       do n=1,nmax 
          z   = n * d 
          an  = z * qz 
          c   = dcos(an) 
          ac  =  caillr2(qortho,qz,q0,a,alpha,z) 
!!          ac2  =  caille2(qortho,qz,q0,a,alpha,z)                     
!!      write(6,*) ac , ac2                                             
          ac  = ac * dexp( - (z/rres)**2 ) 
          sum = sum + ac * c 
!!          ac2 = ac2* dexp( - (z/rres)**2 )                            
!!          sum2 = sum2 + ac2*c                                         
          if(ac.lt.eps) goto 100 
       enddo 
       write(6,*)'caillr3: nmax not sufficient!' 
  100  continue 
                                                                        
       sum = 2.0d0*sum + caillr2(qortho,qz,q0,a,alpha,zeps) 
!!       sum2 = 2.0d0*sum2 + caille2(qortho,qz,q0,a,alpha,zeps)         
                                                                        
       caillr3 = sum 
!!???  caillr3 = sum /( rres *sqrtpi )   !!????                         
!!!!       write(6,*)qz,qortho,sum                                      
                                                                        
       if(iout().gt.1) then 
        write(6,*)'caillr3( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'         rres   =', rres 
        write(6,*)'        )       =', cailler3 
        write(6,*)'-----------------==============================' 
       endif 
!c-test-                                                                
!c     write(6,*)'n=',n,' ac=',ac,'caille3=',sum                        
!c-teste-                                                               
       return 
      END                                           
                                                                        
                                                                        
       function cail2r3(qortho, qz, q0, d , a, alpha, nmax) 
!------====================================================------------ 
!                                                                       
! z-fouriersum of the caille theory  (with resolution)                  
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       parameter    (sqrtpi=1.772453850906d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       zeps = 0.01d0 
                                                                        
       sum = 0.0d0 
       sum2= 0.0d0 
                                                                        
       do n=1,nmax 
          z   = n * d 
          an  = z * qz 
          c   = dcos(an) 
!!          ac  =  caillr2(qortho,qz,q0,a,alpha,z)                      
          ac2  =  caille2(qortho,qz,q0,a,alpha,z) 
!!      write(6,*) ac , ac2                                             
!!          ac  = ac * dexp( - (z/rres)**2 )                            
!!          sum = sum + ac * c                                          
          ac2 = ac2* dexp( - (z/rres)**2 ) 
          sum2 = sum2 + ac2*c 
          if(ac2.lt.eps) goto 100 
       enddo 
       write(6,*)'cail2r3: nmax not sufficient!' 
  100  continue 
                                                                        
!!       sum = 2.0d0*sum + caillr2(qortho,qz,q0,a,alpha,zeps)           
       sum2 = 2.0d0*sum2 + caille2(qortho,qz,q0,a,alpha,zeps) 
                                                                        
       cail2r3 = sum2 
!!???  caillr3 = sum /( rres *sqrtpi )   !!????                         
!!       write(6,*)qz,qortho,sum2                                       
                                                                        
       if(iout().gt.1) then 
        write(6,*)'cail2r3( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'         rres   =', rres 
        write(6,*)'        )       =', cailler3 
        write(6,*)'-----------------==============================' 
       endif 
!c-test-                                                                
!c     write(6,*)'n=',n,' ac=',ac,'cail2r3=',sum                        
!c-teste-                                                               
       return 
      END                                           
                                                                        
                                                                        
       function caills3(qortho, qz, q0, d , a, alpha, peakf,dwf) 
!------====================================================------------ 
!                                                                       
! series expansion of the caille theory                                 
!                                                                       
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       parameter    (pi=3.1415926535898d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       earg = -a*qz*qz*egamma*dwf 
       qzq  =  qz*qz 
       qo4  = qortho**4 
       prefak  = dexp(earg) * a * dsqrt(alpha) / (2 *pi*pi) 
       caills3 =           qzq / (qzq + alpha*qo4) *(2*pi)**2 
! --- add the bragg peak --                                             
       qbragg  = 2*pi/d 
       qzb     = (qz-qbragg)**2 
       caills3 = caills3 + qzq / (qzb+ alpha*qo4)*peakf 
       caills3 = caills3 * prefak 
                                                                        
       if(iout().gt.0) then 
        write(6,*)'caills3( qz     =', qz 
        write(6,*)'         qortho =', qortho 
        write(6,*)'         a      =', a 
        write(6,*)'         alpha  =', alpha 
        write(6,*)'         d      =', d 
        write(6,*)'         q0     =', q0 
        write(6,*)'         peakf  =', peakf 
        write(6,*)'        )       =', cailles3 
        write(6,*)'-----------------==============================' 
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function caille2(qortho, qz, q0, a, alpha, z) 
!------=============================================------------------- 
!                                                                       
!                                                                       
!    exp(-mue/2*(2*gamma+e1(q0*z*sqrt(alpha)) * 1/qortho**2 *           
!                                                                       
!    (q0/qortho)^(-mue) * caille1(mue,beta)                             
!                                                                       
!                                                                       
!  with caille1 =                                                       
!                                                                       
!     inf                                                               
!      /                                                                
!      | j0(x) x^(1-mue) exp(-1/2 mue e1(x^2/beta)) dx                  
!      /                                                                
!     0                                                                 
!                                                                       
!  mue    = a * qz^2                                                    
!  beta   = 4*qortho^2*sqrt(alpha)*z                                    
!                                                                       
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! uses the following nag-routines:                                      
!      s13aaf(x,ier)   --->  e1(x)             exponential-integral     
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       ifail = -1 
                                                                        
                                                                        
! ---- prepare the parameters ---                                       
       xmue   =  a * qz*qz 
!old   sqa    =  dsqrt(alpha)                                           
       sqa    =  dsqrt(dabs(alpha)) 
       qoq    =  qortho * qortho 
       beta   =  4.d0*qoq*sqa*z 
                                                                        
! ---- argument of exponential prefactor (debey-waller-type) ---        
       if(ie1swit.ne.0) then 
         ier    = ifail 
         arg    = -0.5d0*xmue*(2.0d0*egamma+s13aaf(q0*z*sqa,ier)) 
       else 
         arg    = -0.5d0*xmue*(2.0d0*egamma) 
       endif 
       dwf    = dexp(arg) 
                                                                        
       caille2 = caille1(xmue,beta) * (dwf/qoq) *                       &
     &           (q0*0.5d0/qortho)**(-xmue)                             
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
       function caillr2(qortho, qz, q0, a, alpha, z) 
!------=============================================------------------- 
!                                                                       
!                                                                       
!    exp(-mue/2*(2*gamma+e1(q0*z*sqrt(alpha)) * 1/qortho**2 *           
!                                                                       
!    (q0/qortho)^(-mue) * caillr1(mue,beta)                             
!                                                                       
!                                                                       
!  with caille1 =                                                       
!                                                                       
!     inf                                                               
!      /                                                                
!      | j0(x) x^(1-mue) exp(-1/2 mue e1(x^2/beta)) dx                  
!      /                                                                
!     0                                                                 
!                                                                       
!  mue    = a * qz^2                                                    
!  beta   = 4*qortho^2*sqrt(alpha)*z                                    
!                                                                       
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! uses the following nag-routines:                                      
!      s13aaf(x,ier)   --->  e1(x)             exponential-integral     
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       common/cailff/cmue,cbeta,crres 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       ifail = -1 
                                                                        
                                                                        
! ---- prepare the parameters ---                                       
       xmue   =  a * qz*qz 
!old   sqa    =  dsqrt(alpha)                                           
       sqa    =  dsqrt(dabs(alpha)) 
       qoq    =  qortho * qortho 
       beta   =  4.d0*qoq*sqa*z 
       crres  =  rres * qortho 
!!!!!--> neu 10.11.98 --> !!!!!                                         
       if(crres.lt.80.0d0) crres = 80.0d0 
!!!!!<-----------------------------------                               
                                                                        
! ---- argument of exponential prefactor (debey-waller-type) ---        
       if(ie1swit.ne.0) then 
         ier    = ifail 
         arg    = -0.5d0*xmue*(2.0d0*egamma+s13aaf(q0*z*sqa,ier)) 
       else 
         arg    = -0.5d0*xmue*(2.0d0*egamma) 
       endif 
       dwf    = dexp(arg) 
                                                                        
      q0 = abs(q0) 
                                                                        
      caillr2 = caillr1(xmue,beta) * (dwf/qoq) *                        &
     &           (q0*0.5d0/qortho)**(-xmue)                             
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function caille1(xmue, beta) 
!------============================------------------------------------ 
!                                                                       
!     inf                                                               
!      /                                                                
!      | j0(x) x^(1-mue) exp(-1/2 mue e1(x^2/beta)) dx                  
!      /                                                                
!     0                                                                 
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! uses the following nag-routines:                                      
!      s17aef(x,ier)   --->  j0(x)             bessel-function j0       
!      s13aaf(x,ier)   --->  e1(x)             exponential-integral     
!      s14aaf(x,ier)   --->  gamma(x)          gamma-function           
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       external cailf1 
       external cailf2 
       common/cailff/cmue,cbeta,crres 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       ifail = -1 
       maxi  = maxit 
                                                                        
       if(    xmue.lt.0.0d0      .or.                                   &
     &        beta.lt.0.0d0      .or.                                   &
     &        eps .le.0.0d0        ) then                               
         write(6,*)'caille1: wrong arguments !' 
         write(6,*)'xmue   =',xmue 
         write(6,*)'beta   =',beta 
         write(6,*)'eps    =',eps 
         caille1 = 0.0d0 
         return 
       endif 
                                                                        
                                                                        
! ---- determine upper limit for numerical integration ....             
       xm1   = dsqrt(beta) 
       aleb  = -beta * dlog(2.0d0*eps/(xmue+1.0d-6)) 
       if(aleb.le.0) then 
         xm = xm1 
       else 
         xm2 = dsqrt(aleb) 
         if(xm2 .gt. xm1) then 
            xm = xm2 
         else 
            xm = xm1 
         endif 
       endif 
!c-test-                                                                
!c     write(6,*)'xm=    ',xm                                           
!c-teste                                                                
                                                                        
       cmue  = xmue 
       cbeta = beta 
                                                                        
!c - xm-check                                                           
       if(iout().gt.3) then 
         ier = ifail 
         e1  = s13aaf(xm*xm/cbeta,ier) 
         e11 = dexp(-0.5d0*cmue*e1) - 1.0d0 
         write(6,*)'caille1 check e11 = ',e11 
       endif 
!c - xm-checkend                                                        
                                                                        
! ---- section i  for mue < 2 ----------------------------------------- 
       if(xmue.lt.xmuswit) then 
!                                                       gamma(1-mue/2)  
! ---- part1: int(0..inf)[j0(x)*x^(1-mue) dx]=2^(1-mue)*--------------  
!                                                       gamma( mue/2)   
!                                                                       
         xmue2 = 0.5d0 * xmue 
         ier   = ifail 
         g1    = s14aaf(1.0d0-xmue2,ier) 
         ier   = ifail 
         g2    = s14aaf(xmue2,ier) 
         part1 = (g1/g2) * 2.0d0**(1.0d0-xmue) 
!c-test-                                                                
!c     write(6,*)'part11 ',part1                                        
!c-teste                                                                
                                                                        
! ---- part2: int(0..xmax)[j0(x)*x^(1-mue)*(exp(-mue/2*e1(x^2/beta))-1) 
         a     = 0.0d0 
         b     = xm 
         epsi  = 0.1d0*eps 
                                  ! flag fuer Fehlermeldung             
         erracc= 111 
         part2 = cadapint(cailf1, a, b,  epsi   ,maxi   ,erracc   ) 
         part2 = part2 - (b**(2.0d0-xmue))/(2.0d0-xmue) 
!c-test-                                                                
!c     write(6,*)'part21 ',part2                                        
!c-teste                                                                
! -----                                                                 
         caille1 = part1 + part2 
                                                                        
       else 
! ---- section ii for mue > 2 ----------------------------------------- 
!                                                                       
! ---- part1: int(0..xm)[j0(x)*x^(1-mue)*exp(-mue/2*e1(x^2/beta))*dx]   
!                                                                       
         a = 0.0d0 
         b = xm 
         epsi  = 0.1d0*eps 
                                  ! flag fuer Fehlermeldung             
         erracc= 222 
         part1 = cadapint(cailf2, a, b,  epsi   ,maxi   ,erracc   ) 
!c-test-                                                                
!c     write(6,*)'part12 ',part1                                        
!c-teste                                                                
!                                                                       
! ---- part2: int(xm..inf)[j0(x)*x^(1-mue)*dx]  = series-expansion      
!                                                 (see xnuej0)          
!                                                                       
! ---- for optimum accuracy the upper limit must be between 25..30      
!      it should be approximately equal to (n+1/4)*pi                   
!      more accuraetly xupper is the solution of:                       
!      cos(z+pi/4)*z - sin(z+pi/4) = 0  which is near (n+1/4)*pi        
!                                                                       
!c       xupper = 25.87954d0                                            
!c       xupper = 25.93                                                 
!c       eps2   = 1.0d-15                                               
         max2   = maxit2 
         xnu    = 1.0d0-xmue 
         if(xm.lt.xupper) then 
            part2  = xnuej0(xupper,xnu,eps2,max2) -                     &
     &               xnuej0(xm    ,xnu,eps2,max2)                       
         else 
           write(6,*)'caille1 warning: xupper < xm =',xm 
            part2  = xnuej0(xupper,xnu,eps2,max2) -                     &
     &               xnuej0(xm    ,xnu,eps2,max2)                       
           write(6,*)'            part2 = 0',part2 
         endif 
!c-test-                                                                
!c     write(6,*)'part22 ',part2                                        
!c-teste                                                                
! -----                                                                 
         caille1 = part1 + part2 
                                                                        
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function caillr1(xmue, beta) 
!------============================------------------------------------ 
!                                                                       
!     inf                                                               
!      /                                                                
!      | j0(x) x^(1-mue) exp(-1/2 mue e1(x^2/beta)) * exp(-(x/crres)^2) 
!      /                                                                
!     0                                                                 
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! uses the following nag-routines:                                      
!      s17aef(x,ier)   --->  j0(x)             bessel-function j0       
!      s13aaf(x,ier)   --->  e1(x)             exponential-integral     
!      s14aaf(x,ier)   --->  gamma(x)          gamma-function           
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       external cailfr 
       common/cailff/cmue,cbeta,crres 
       common/caisys/ eps, eps2, xmuswit, xupper, rres,                 &
     &                strechf, rellim, abslim, wrnlim,rlim4,splitf4,    &
     &                maxit, maxit2, ie1swit, nppow, inf4, nstart4      
!                                                                       
                                                                        
       ifail = -1 
       maxi  = maxit 
                                                                        
       if(    xmue.lt.0.0d0      .or.                                   &
     &        beta.lt.0.0d0      .or.                                   &
     &        eps .le.0.0d0        ) then                               
         write(6,*)'caillr1: wrong arguments !' 
         write(6,*)'xnue   =',xnue 
         write(6,*)'beta   =',beta 
         write(6,*)'eps    =',eps 
         caillr1 = 0.0d0 
         return 
       endif 
                                                                        
                                                                        
! ---- determine upper limit for numerical integration ....             
       aleb  = - dlog(eps) 
       if(aleb.le.0) then 
         xm = 1.0d0 
       else 
         xm  = dsqrt(aleb) * crres 
       endif 
!c-test-                                                                
!c     write(6,*)'xm=    ',xm                                           
!c-teste                                                                
                                                                        
       cmue  = xmue 
       cbeta = beta 
                                                                        
!c - xm-check                                                           
       if(iout().gt.3) then 
         e11 = cailfr(xm) 
         write(6,*)'caille1 check e11 = ',e11 
       endif 
!c - xm-checkend                                                        
                                                                        
! --- integration ---                                                   
         a     = 0.0d0 
         b     = xm 
         epsi  =  eps 
                                  ! flag fuer Fehlermeldung             
         erracc= 333 
         part1 = cadapint(cailfr, a, b,  epsi   ,maxi   ,erracc   ) 
                                                                        
         caillr1 = part1 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       function xnuej0(x,xnue,eps,maxit) 
! -----=================================------------------------------- 
!  series expansion for:                                                
!                        maxit|term<eps                                 
!    x                      ----                                        
!    /                      \    (-)^n * x^(nue+1+2n)                   
!    | x**xnue * j0(x) dx =      -----------------------                
!    /                      /    2^(2n)*n!*n!*(nue+1+2n)                
!    ...                    ----                                        
!                           n=0                                         
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
                                                                        
         xnfak = 1.0d0 
         two2n = 1.0d0 
         xnu   = x**(xnue+1.0d0) 
         xx    = x*x 
         x2n   = 1.0d0 
         sum   = 1.0d0/(xnue+1.0d0) 
         do n=1,maxit 
           xnfak = xnfak * n / x 
           two2n =-two2n*4.0d0 
           add   = 1.0d0/(xnfak*xnfak*two2n*(2.0d0*n+xnue+1.0d0)) 
           sum   = sum + add 
           if(dabs(add).lt.eps) then 
             xnuej0 = sum * xnu 
             return 
           endif 
         enddo 
         write(6,*)'xnuej0 did not converge !' 
         write(6,*)'x      =',x 
         write(6,*)'eps    =',eps 
         write(6,*)'add    =',add 
         write(6,*)'maxit  =',maxit 
         write(6,*)'sum    =',sum 
         xnuej0 = sum * xnu 
                                                                        
       return 
      END                                           
                                                                        
       function cailf1(x) 
!------==================---------------------------------------------- 
! integrand function for caille1 :                                      
!                                                                       
!   j0(x) * x^(1-mue) * ( exp(-mue/2*e1(x^2/beta)) - 1)                 
!                                                                       
!   the parameters mue (cmue) and beta (cbeta) are transferred via      
!   common /cailff/                                                     
!                                                                       
! uses nag routines:                                                    
!   s13aaf(x,ier) ---> e1(x)                 exponential-integral       
!   s17aef(x,ier) ---> j0(x)                 bessel-function            
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       common/cailff/cmue,cbeta,crres 
                                                                        
       ifail =  -1 
       ier   = ifail 
       e1  = s13aaf(x*x/cbeta,ier) 
       e11 = dexp(-0.5d0*cmue*e1) - 1.0d0 
       ier   =  ifail 
       cailf1 = (x**(1.0d0-cmue))*( s17aef(x,ier)*e11 + 1.0d0 ) 
       return 
      END                                           
                                                                        
       function cailf2(x) 
!------==================---------------------------------------------- 
! integrand function for caille1 :                                      
!                                                                       
!   j0(x) * x^(1-mue) *  exp(-mue/2*e1(x^2/beta))                       
!                                                                       
!   the parameters mue (cmue) and beta (cbeta) are transferred via      
!   common /cailff/                                                     
!                                                                       
! uses nag routines:                                                    
!   s13aaf(x,ier) ---> e1(x)                 exponential-integral       
!   s17aef(x,ier) ---> j0(x)                 bessel-function            
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       common/cailff/cmue,cbeta,crres 
                                                                        
       ifail =  -1 
       ier   = ifail 
       e1  = s13aaf(x*x/cbeta,ier) 
       e11 = dexp(-0.5d0*cmue*e1) 
       ier   =  ifail 
       cailf2 = (x**(1.0d0-cmue))*s17aef(x,ier)*e11 
       return 
      END                                           
                                                                        
       function cailfr(x) 
!------==================---------------------------------------------- 
! integrand function for caillr1 :                                      
!                                                                       
!   j0(x) * x^(1-mue) *  exp(-mue/2*e1(x^2/beta)) * exp(-(x/crres)^2)   
!                                                                       
!   the parameters mue (cmue) and beta (cbeta) and (crres)              
!   are transferred via                                                 
!   common /cailff/                                                     
!                                                                       
! uses nag routines:                                                    
!   s13aaf(x,ier) ---> e1(x)                 exponential-integral       
!   s17aef(x,ier) ---> j0(x)                 bessel-function            
!---------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z) 
       common/cailff/cmue,cbeta,crres 
                                                                        
       ifail =  -1 
       ier   = ifail 
       x = abs(x) 
       e1  = s13aaf(x*x/cbeta,ier) 
       e11 = dexp(-0.5d0*cmue*e1-(x/crres)**2) 
       ier   =  ifail 
       cailfr = (x**(1.0d0-cmue))*s17aef(x,ier)*e11 
       return 
      END                                           
                                                                        
                                                                        
                                                                        
       subroutine polyfit(x, y, m, a, n, r) 
! -----====================================---------------------------- 
!  polynomial fit to data in tables x(1..m), y(1..m)  (input)           
!  coefficients are returned in a(0..n) (output)                        
!  max. degree is given by n (input)                                    
!  residual error is given by r (rms) (output)                          
!                                                                       
! use nag routine(s):                                                   
!     f04aef       real symmetric equations solution                    
! --------------------------------------------------------------------- 
!                                                                       
       implicit real*8 (a-h,o-z) 
       parameter(maxdeg=30) 
       parameter(mdim=maxdeg+1) 
                                                                        
       dimension am(mdim,mdim),b(mdim,1),c(mdim,1),bb(mdim,1),          &
     &           wkspce(mdim)                                           
       dimension al(mdim,mdim) 
                                                                        
       dimension xx(2*mdim) 
       dimension x(m),y(m),a(*) 
                                                                        
       ifail = -1 
                                                                        
       if(n.gt.maxdeg) then 
         r = 9.99999d33 
         write(6,*)'polyfit degree n=',n,' is larger than maxdeg' 
         return 
       endif 
                                                                        
       if(n+1.gt.m) then 
         n = m-1 
         write(6,*)'ployfit: too few points   (',m,')' 
         write(6,*)'         degree reduced to ',n 
       endif 
                                                                        
! --- x-scalingfactor ----                                              
       xscale = 0.1d0 * x(m) 
!cc    xscale =1.0d0                                                    
                                                                        
       do j=1,2*n+1 
         xx(j) = 0.0d0 
         b(j,1)= 0.0d0 
       enddo 
                                                                        
       do i=1,m 
         xh = 1.0d0 
         do j=1,2*n+1 
           xx(j) = xx(j) + xh 
           if(j.le.n+1) then 
             b(j,1) = b(j,1) + xh * y(i) 
           endif 
           xh    = xh * x(i)/xscale 
         enddo 
       enddo 
                                                                        
! --- prepare a-matrix ----                                             
      do i=1,n+1 
        do j=1,n+1 
          am(j,i) = xx(1+(i-1)+(j-1)) 
        enddo 
      enddo 
!                                                                       
! --- solve equation ----                                               
      ier = ifail 
      call f04aef(am,mdim,b,mdim,n+1,1,c,mdim,                          &
     &            wkspce,al,mdim,bb,mdim,ier)                           
!                                                                       
! --- scale c and put it to output                                      
      xh   = 1.0d0 
      do i=1,n+1 
        a(i) = c(i,1)*xh 
        xh   = xh / xscale 
      enddo 
!                                                                       
! --- determine the quality of solution                                 
      r    = 0.0d0 
      do i=1,m 
        sum = 0.0d0 
        xh   = 1.0d0 
        do j=1,n+1 
           sum = sum + xh * a(j) 
           xh  = xh * x(i) 
        enddo 
        d = sum-y(i) 
!c-test                                                                 
!c      write(6,*)i,d                                                   
!c-teste                                                                
        r = r + d*d 
      enddo 
      r = dsqrt(r) 
!c-test                                                                 
!c    write(6,*)'r=',r                                                  
!c-teste                                                                
      return 
      END                                           
                                                                        
                                                                        
                                                                        
      function eisum(q,d,a,nppow,n1,n2) 
! --- ================================= --------------------------------
!       n2                                                              
!     -----                                                             
!      \                                  nppow                         
!      |  2  cos(n*q*d) * exp(-a*n*d) * n^                              
!      /                                                                
!     -----                                                             
!     n=n1                                                              
!                                                                       
!                                                                       
!  if n2 is negative: then formula for n2--> inf is used.               
!                                                                       
! ----------------------------------------------------------------------
                                                                        
      implicit real*8 (a-h,o-s,u-z) 
      implicit complex*16 (t) 
      complex*16 cdexp, dcmplx 
                                                                        
      parameter(mppow=5) 
                                                                        
      if(iout().gt.3) then 
        write(6,*)'eisum called with parameters:' 
        write(6,*)'q     = ',q 
        write(6,*)'d     = ',d 
        write(6,*)'a     = ',a 
        write(6,*)'nppow = ',nppow 
        write(6,*)'n1    = ',n1 
        write(6,*)'n2    = ',n2 
      endif 
                                                                        
                                                                        
                                                                        
      if(nppow.gt.mppow.or.nppow.lt.0) then 
        write(6,*)'eisum: nppow is out of range:',nppow 
        eisum = 0.0d0 
        return 
      endif 
                                                                        
      if(n2.lt.0) goto 1000 
                                                                        
      iselect = nppow + 1 
      goto(10,20,30,40,50,60),iselect 
                                                                        
! -se.f  ---------------------------------------------------------------
   10 continue 
      t2 =  dexp(d*a) 
      t7 = 1/(cdexp(dcmplx(0.d0,1.d0)*d*q)-t2)*t2 
      t9 = (n2+1)*d 
      t17 = n1*d 
      t25 = t7*cdexp(dcmplx(0.d0,1.d0)*t9*q)/cdexp(t9*a)-               &
     & t7*cdexp(dcmplx(0.d0,1.d0)*t17*q)/cdexp(t17*a)                   
      eisum = dreal(t25) * 2.0d0 
      return 
! -se1.f ---------------------------------------------------------------
   20 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t4**2 
      t7 = t2**2 
      t9 = 1/(t6-2*t5+t7) 
      t10 = n2+1 
      t11 = d*t10 
      t13 = cdexp(dcmplx(0.d0,1.d0)*t11*q) 
      t17 = 1/cdexp(t11*a) 
      t22 = t2/(t4-t2) 
      t26 = d*n1 
      t28 = cdexp(dcmplx(0.d0,1.d0)*t26*q) 
      t32 = 1/cdexp(t26*a) 
      t38 = -t5*t9*t13*t17+t22*t13*t17*t10+t5*t9*t28*t32-t22*t28*t32*n1 
      eisum = dreal(t38) * 2.0d0 
      return 
! -se2.f ---------------------------------------------------------------
   30 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t7 = t5*(t4+t2) 
      t8 = t4**2 
      t11 = t2**2 
      t15 = 1/(t8*t4-3*t8*t2+3*t4*t11-t11*t2) 
      t16 = n2+1 
      t17 = d*t16 
      t19 = cdexp(dcmplx(0.d0,1.d0)*t17*q) 
      t23 = 1/cdexp(t17*a) 
      t28 = t5/(t8-2*t5+t11) 
      t29 = t19*t23 
      t34 = t2/(t4-t2) 
      t35 = t16**2 
      t38 = d*n1 
      t40 = cdexp(dcmplx(0.d0,1.d0)*t38*q) 
      t44 = 1/cdexp(t38*a) 
      t47 = t40*t44 
      t50 = n1**2 
      t53 = t7*t15*t19*t23-2*t28*t29*t16+t34*t29*t35-t7*t15*t40*t44+    &
     & 2*t28*t47*n1-t34*t47*t50                                         
      eisum = dreal(t53) * 2.0d0 
      return 
! -se3.f ---------------------------------------------------------------
   40 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t4**2 
      t9 = t5*(t6+4*t5+t7) 
      t10 = t6**2 
      t11 = t6*t2 
      t14 = t7*t4 
      t16 = t7**2 
      t18 = 1/(t10-4*t11*t4+6*t6*t7-4*t2*t14+t16) 
      t19 = n2+1 
      t20 = d*t19 
      t22 = cdexp(dcmplx(0.d0,1.d0)*t20*q) 
      t26 = 1/cdexp(t20*a) 
      t30 = t5*(t4+t2) 
      t34 = 1/(t11-3*t4*t6+3*t7*t2-t14) 
      t41 = t2/(-t4+t2) 
      t42 = t22*t26 
      t43 = t19**2 
      t49 = t5/(t7-2*t5+t6) 
      t52 = d*n1 
      t54 = cdexp(dcmplx(0.d0,1.d0)*t52*q) 
      t58 = 1/cdexp(t52*a) 
      t65 = t54*t58 
      t66 = n1**2 
      t72 = -t9*t18*t22*t26-3*t30*t34*t22*t26*t19-t41*t42*t43*t19-      &
     & 3*t49*t42*t43+t9*t18*t54*t58+3*t30*t34*t54*t58*n1+               &
     & t41*t65*t66*n1+3*t49*t65*t66                                     
      eisum = dreal(t72) * 2.0d0 
      return 
! -se4.f ---------------------------------------------------------------
   50 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t6*t2 
      t8 = t4*t6 
      t9 = t4**2 
      t10 = t9*t2 
      t11 = t9*t4 
      t13 = t5*(t7+11*t8+11*t10+t11) 
      t14 = t6**2 
      t19 = t9**2 
      t23 = 1/(t14*t2-5*t14*t4+10*t7*t9-10*t6*t11+5*t2*t19-t19*t4) 
      t24 = n2+1 
      t25 = d*t24 
      t27 = cdexp(dcmplx(0.d0,1.d0)*t25*q) 
      t31 = 1/cdexp(t25*a) 
      t35 = t5*(t6+4*t5+t9) 
      t40 = 1/(t14-4*t7*t4+6*t6*t9-4*t2*t11+t19) 
      t47 = t5/(t9-2*t5+t6) 
      t48 = t27*t31 
      t49 = t24**2 
      t54 = t5*(t4+t2) 
      t56 = 1/(t7-3*t8+3*t10-t11) 
      t63 = t2/(-t4+t2) 
      t64 = t49**2 
      t67 = d*n1 
      t69 = cdexp(dcmplx(0.d0,1.d0)*t67*q) 
      t73 = 1/cdexp(t67*a) 
      t80 = t69*t73 
      t81 = n1**2 
      t89 = t81**2 
      t92 = -t13*t23*t27*t31-4*t35*t40*t27*t31*t24-4*t47*t48*t49*t24-   &
     & 6*t54*t56*t27*t31*t49-t63*t48*t64+t13*t23*t69*t73+               &
     & 4*t35*t40*t69*t73*n1+4*t47*t80*t81*n1+                           &
     & 6*t54*t56*t69*t73*t81+t63*t80*t89                                
      eisum = dreal(t92) * 2.0d0 
      return 
! -se5.f ---------------------------------------------------------------
   60 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t6**2 
      t8 = t6*t2 
      t9 = t8*t4 
      t10 = t4**2 
      t11 = t6*t10 
      t12 = t10*t4 
      t13 = t2*t12 
      t14 = t10**2 
      t16 = t5*(t7+26*t9+66*t11+26*t13+t14) 
      t18 = t7*t2 
      t23 = t14*t4 
      t27 = 1/(t7*t6-6*t18*t4+15*t7*t10-20*t8*t12+                      &
     & 15*t6*t14-6*t2*t23+t14*t10)                                      
      t28 = n2+1 
      t29 = d*t28 
      t31 = cdexp(dcmplx(0.d0,1.d0)*t29*q) 
      t35 = 1/cdexp(t29*a) 
      t38 = t4*t6 
      t39 = t10*t2 
      t41 = t5*(t8+11*t38+11*t39+t12) 
      t47 = 1/(t18-5*t7*t4+10*t8*t10-10*t6*t12+5*t2*t14-t23) 
      t54 = t2/(-t4+t2) 
      t55 = t31*t35 
      t56 = t28**2 
      t57 = t56**2 
      t62 = t5*(t4+t2) 
      t64 = 1/(t8-3*t38+3*t39-t12) 
      t71 = t5*(t6+4*t5+t10) 
      t73 = 1/(t7-4*t9+6*t11-4*t13+t14) 
      t80 = t5/(t10-2*t5+t6) 
      t83 = d*n1 
      t85 = cdexp(dcmplx(0.d0,1.d0)*t83*q) 
      t89 = 1/cdexp(t83*a) 
      t96 = t85*t89 
      t97 = n1**2 
      t98 = t97**2 
      t113 = -t16*t27*t31*t35-5*t41*t47*t31*t35*t28-t54*t55*t57*t28-    &
     & 10*t62*t64*t31*t35*t56*t28-10*t71*t73*t31*t35*t56-               &
     & 5*t80*t55*t57+t16*t27*t85*t89+5*t41*t47*t85*t89*n1+              &
     & t54*t96*t98*n1+10*t62*t64*t85*t89*t97*n1+                        &
     & 10*t71*t73*t85*t89*t97+5*t80*t96*t98                             
      eisum = dreal(t113) * 2.0d0 
      return 
!                                                                       
! ----- infinit upper bound of sum ....                                 
 1000 continue 
!                                                                       
                                                                        
      iselect = nppow + 1 
      goto(110,120,130,140,150,160),iselect 
                                                                        
! -sei.f ---------------------------------------------------------------
  110 continue 
      t2 =  dexp(d*a) 
      t8 = d*n1 
      t16 = -t2/(cdexp(dcmplx(0.d0,1.d0)*d*q)-t2)*                      &
     & cdexp(dcmplx(0.d0,1.d0)*t8*q)/cdexp(t8*a)                        
      eisum = dreal(t16) * 2.0d0 
      return 
! -sei1.f --------------------------------------------------------------
  120 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t4**2 
      t7 = t2**2 
      t10 = d*n1 
      t12 = cdexp(dcmplx(0.d0,1.d0)*t10*q) 
      t16 = 1/cdexp(t10*a) 
      t25 = t5/(t6-2*t5+t7)*t12*t16-t2/(t4-t2)*t12*t16*n1 
      eisum = dreal(t25) * 2.0d0 
      return 
! -sei2.f --------------------------------------------------------------
  130 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t8 = t2**2 
      t11 = t4**2 
      t16 = d*n1 
      t18 = cdexp(dcmplx(0.d0,1.d0)*t16*q) 
      t22 = 1/cdexp(t16*a) 
      t28 = t18*t22 
      t34 = n1**2 
      t37 = t5*(t4+t2)/(t8*t2-3*t4*t8+3*t11*t2-t11*t4)*t18*t22+         &
     & 2*t5/(t11-2*t5+t8)*t28*n1+t2/(-t4+t2)*t28*t34                    
      eisum = dreal(t37) * 2.0d0 
      return 
! -sei3.f --------------------------------------------------------------
  140 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t4**2 
      t10 = t6**2 
      t11 = t6*t2 
      t14 = t7*t4 
      t16 = t7**2 
      t19 = d*n1 
      t21 = cdexp(dcmplx(0.d0,1.d0)*t19*q) 
      t25 = 1/cdexp(t19*a) 
      t41 = t21*t25 
      t42 = n1**2 
      t51 = t5*(t6+4*t5+t7)/(t10-4*t11*t4+6*t6*t7-4*t2*t14+t16)*t21*t25+&
     & 3*t5*(t4+t2)/(t11-3*t4*t6+3*t7*t2-t14)*t21*t25*n1+               &
     & t2/(-t4+t2)*t41*t42*n1+3*t5/(t7-2*t5+t6)*t41*t42                 
      eisum = dreal(t51) * 2.0d0 
      return 
! -sei4.f --------------------------------------------------------------
  150 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t6*t2 
      t8 = t4*t6 
      t9 = t4**2 
      t10 = t9*t2 
      t11 = t9*t4 
      t14 = t6**2 
      t19 = t9**2 
      t24 = d*n1 
      t26 = cdexp(dcmplx(0.d0,1.d0)*t24*q) 
      t30 = 1/cdexp(t24*a) 
      t47 = t26*t30 
      t48 = n1**2 
      t63 = t48**2 
      t66 = t5*(t7+11*t8+11*t10+t11)/(t14*t2-5*t14*t4+10*t7*t9-         &
     & 10*t6*t11+5*t2*t19-t19*t4)*t26*t30+4*t5*(t6+4*t5+t9)/            &
     & (t14-4*t7*t4+6*t6*t9-4*t2*t11+t19)*t26*t30*n1+4*t5/              &
     & (t9-2*t5+t6)*t47*t48*n1+6*t5*(t4+t2)/                            &
     & (t7-3*t8+3*t10-t11)*t26*t30*t48+t2/                              &
     & (-t4+t2)*t47*t63                                                 
      eisum = dreal(t66) * 2.0d0 
      return 
! -sei5.f --------------------------------------------------------------
  160 continue 
      t2 =  dexp(d*a) 
      t4 = cdexp(dcmplx(0.d0,1.d0)*d*q) 
      t5 = t2*t4 
      t6 = t2**2 
      t7 = t6**2 
      t8 = t6*t2 
      t9 = t8*t4 
      t10 = t4**2 
      t11 = t6*t10 
      t12 = t10*t4 
      t13 = t2*t12 
      t14 = t10**2 
      t18 = t7*t2 
      t23 = t14*t4 
      t28 = d*n1 
      t30 = cdexp(dcmplx(0.d0,1.d0)*t28*q) 
      t34 = 1/cdexp(t28*a) 
      t37 = t4*t6 
      t38 = t10*t2 
      t54 = t30*t34 
      t55 = n1**2 
      t56 = t55**2 
      t82 = t5*(t7+26*t9+66*t11+26*t13+t14)/(t7*t6-6*t18*t4+15*t7*t10-20&
     &*t8*t12+15*t6*t14-6*t2*t23+t14*t10)*t30*t34-5*t5*(t8+11*t37+11*t38&
     &+t12)/(t23-5*t2*t14+10*t6*t12-10*t8*t10+5*t7*t4-t18)*t30*t34*n1-t2&
     &/(t4-t2)*t54*t56*n1-10*t5*(t4+t2)/(t12-3*t38+3*t37-t8)*t30*t34*t55&
     &*n1+10*t5*(t6+4*t5+t10)/(t7-4*t9+6*t11-4*t13+t14)*t30*t34*t55+5*t5&
     &/(t10-2*t5+t6)*t54*t56                                            
      eisum = dreal(t82) * 2.0d0 
      return 
                                                                        
! -------------- end eisum ---------------------------------------------
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
       function fbsj0(x) 
       implicit real*8 (a-h,o-z) 
       common /cfbsj0/qortho,xexpo,decay 
                                                                        
       fbsj0 = (x**(1.0d0-xexpo))*s17aef(x,ier)*dexp(-decay*x) 
       return 
      END                                           
                                                                        
       function f3(x) 
       implicit real*8 (a-h,o-z) 
       common /cf3/xmu3,beta3,decay3 
                                                                        
       ier   =  0 
       e1  = s13aaf(x*x/beta3,ier) 
       e11 = dexp(-0.5d0*xmu3*e1) - 1.0d0 
       ier   =  0 
       f3    = (x**(1.0d0-xmu3))*s17aef(x,ier)*e11 
       return 
      END                                           
                                                                        
       function f4(x) 
       implicit real*8 (a-h,o-z) 
       common /cf3/xmu3,beta3,decay3 
                                                                        
       ier   =  0 
       e1  = s13aaf(x*x/beta3,ier) 
       e11 = dexp(-0.5d0*xmu3*e1) 
       ier   =  0 
       f4    = (x**(1.0d0-xmu3))*s17aef(x,ier)*dexp(-decay3*x)*e11 
       return 
      END                                           
                                                                        
       function f1(x) 
       implicit real*8 (a-h,o-z) 
       common /cf1/ alpha,rho,zz 
       ier = 0 
!old   f1 =  (1.d0-s17aef(x*rho,ier)*dexp(-dsqrt(alpha)*x*x*zz))/x      
       f1 =  (1.d0-s17aef(x*rho,ier)*                                   &
     &        dexp(-dsqrt(dabs(alpha))*x*x*zz))/x                       
       return 
      END                                           
                                                                        
       function f1i(rho,z,alpha,q0) 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       ier = 0 
       e11 = s13aaf(q0*q0*dsqrt(alpha)*z,ier) 
       ier = 0 
       e12 = s13aaf(rho*rho/(dsqrt(alpha)*4.0d0*z),ier) 
       f1i=  0.5d0 * ( 2.0d0*egamma + dlog(q0*q0*rho*rho/4.d0)          &
     &                              + e11 + e12 )                       
       return 
      END                                           
                                                                        
       function f2i(xmu) 
       implicit real*8 (a-h,o-z) 
       parameter    (egamma=0.5772156649d0) 
       ier = 0 
       g11 = s14aaf(0.5d0*(1.0d0-xmu)+0.5d0,ier) 
       ier = 0 
       g12 = s14aaf(0.5d0-0.5d0*(1.0d0-xmu),ier) 
       f2i=  g11/g12 * 2.0d0**(1.0d0-xmu) 
       return 
      END                                           
                                                                        
       function cgaus8a (f,xu,xo) 
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
       cgaus8a = sum * range 
                                                                        
       return 
      END                                           
!*ds                                                                    
!*ed                                                                    
       function cadapint(f,a,b,epsilon0,maxiter0,erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
       implicit real*8 (a-h,o-z) 
!                                                                       
       parameter( maxstack = 30) 
!                 --------------->  stacktiefe                          
       dimension s(0:maxstack), xa(0:maxstack), xb(0:maxstack) 
!                                                                       
       logical cray 
                                                        !! aix          
       real*4 xxxx, yyyy, ptxf 
                                                                      !!
       common/outlev/iot,ibild,ierrs,inka, iibuf, xxxx,yyyy, ptxf(20) 
!                                                                       
!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!                                
       character*1 answer 
                                        !!!! test only !!!              
       common/cailff/cmue,cbeta,crres 
!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!                                
                                                                        
       external f 
!                                                                       
!      -----------------------------------------------------------------
         ifrom            = erroraccu+0.1 
         maxiter          = 10000 
         epsilon = epsilon0 
          erroraccu        = 0.d0 
         iterationcounter = 0 
         erg              = 0.0d0 
         itop      = 0 
         xa(itop)  = a 
         xb(itop)  = b 
         s(itop)   = cgaus8a( f, a, b ) 
                                                                        
         do i=1,maxiter 
           if (itop.ge.(maxstack-2)) then 
              erg = erg + s(itop) 
              xbb = xb(itop) 
              itm = itop-1 
              do itop=itm,0,-1 
                 if(xa(itop).eq.xbb) goto 1 
              enddo 
    1         continue 
             epsilon = 5*epsilon 
             write(6,*)'warning! cadaptint stack overflow!' 
           else 
              iterationcounter = iterationcounter + 1 
              itop             = itop +1 
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0 
              xb(itop) =  xb(itop-1) 
              s(itop)  = cgaus8a( f, xa(itop), xb(itop) ) 
              itop             = itop +1 
              xa(itop) =  xa(itop-2) 
              xb(itop) =  xa(itop-1) 
              s(itop)  = cgaus8a( f, xa(itop), xb(itop) ) 
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
     &    'cadapint(',ifrom,') fatal error iterationnumber exceeded!'   
         write(6,*) 'limits a=',a 
         write(6,*) 'limits b=',b 
         write(6,*) 'cmue    =',cmue 
         write(6,*) 'cbeta   =',cbeta 
         write(6,*) 'crres   =',crres 
         write(6,*) 'W(eiter), F(ile), A(bbruch) ---> ?' 
         read(5,'(A)') answer 
         if(answer.eq.'a') stop 
         if(answer.eq.'f') then 
           open(22,file='test.out',status='UNKNOWN') 
           nn = 100 
           dx = (b-a)/nn 
           do i=0,100 
            xx = a + i*dx 
            yy = f(xx) 
            write(22,*) xx,yy 
           enddo 
           close(22) 
           write(6,*)'integrand written to test.out' 
         endif 
         cadapint = 9.99999d33 
         return 
    3    continue 
                                                                        
                                                                        
         cadapint = erg 
                                                                        
         return 
                                                                        
      END                                           
