!*ed                                                                    
! *********** tdpur   fortran ******************************************
!                                                                      *
! ---> version for general data fitting without IMSL or NAG calls       
!                                                                       
!***********************************************************************
!*                                                                     *
!*    all calls to IMSL routines are redifined with calls              *
!*    to minpack for fitting, fftpack5 for fourier transforms and      *
!*    slatec for other functions                                       *
!*                                                                     *
!*    Olaf Holderer, IFF, FZ-Juelich, 2005                             *
!*                                                                     *
!***********************************************************************
!                                                                       
!*ds                                                                    
!*ds                                                                    



subroutine unlsf(func,m,n,xguess,xscale,fscale,iparam,rparam,    &
     &                   x,f,xjac,ixjac)                                
!      Fitting with minpack instead of imsl...                          
                                                                        
       use cincom
       use cincoc
       use xoutxx
       use xroxxx
       use cdata
       use outlev
       use theory
       use constants
       
	character*8 ci 
	real*8 getval 
!                                                                       
!                                                                       
       common/therrc/therro(mtpar,mtcal) 
!                                                                       

       common/selist/isels(mbuf),ifits(mbuf),nsel,numpls 
 
!                                                                       
      dimension iparam(6),rparam(7),x(mfit),f(msmpl),xjac(msmpl,mfit),  &
     &          xguess(mfit),xscale(mfit),fscale(msmpl)                 
      dimension numv(minc),numn(minc) 
      dimension ermat(mfit,msmpl), gmat(mfit,mfit), ginv(mfit,mfit) 
!                                                                       
      dimension   xcenter(mfit), xstepsc(mfit) 
      character*8 xpname(mfit) 
!                                                                       
      external func 
!                                                                       
      logical sqwght,sqwbuf 
      logical autox1,autox2 
      logical found, folgt 
      logical final_thc, lbuffer 
      common/cfunc/iprt,sqwght,x1,x2,autox1,autox2,ferror(msmpl) 
                                                                        
       real    fcssq 
       logical lerrel, lwrtfitdat 
       common/cfunce/ lerrel, lwrtfitdat, fcssq 
                                                                        
                                                                        
! ---- OH: additional parameters for minpack and lapack                 
       double precision xdoub1 
       real xreal1 
       integer icode 
                                                                        
       integer nfev, mode, iflag 
       integer iwa(mfit) 
       real qtf(mfit),wa1(mfit) 
       real wa2(mfit),wa3(mfit),wa4(msmpl) 
       real epsfcn, zero, ftol, gtol, xtol 
       data zero/0.0d0/ 
                                                                        
       data nsig/6/,eps/1.d-5/,delta/1.d-5/,maxfn/50/,iopt/0/ 
                                                                        
       ftol = eps 
       xtol = eps 
       gtol = zero 
       epsfcn = 2.d-7 
       mode = 2 
       iflag = 2 
       write(6,*)' so weit, so gut  *****' 
       call lmdif(func,m,n,x,f,ftol,xtol,gtol,maxfn,                    &
     &            epsfcn,xscale,mode,1.0,iprint,icode,nfev,xjac,ixjac,  &
     &            iwa, qtf,wa1,wa2, wa3,wa4)                            
       call fdjac2(func,m,n,x,f,xjac,ixjac,iflag,epsfcn,wa4) 
                                                                        
       write(6,*)' und nu???' 
        if(icode.eq.0) write(6,*)' **-0- : improper input params  *****' 
        if(icode.eq.1) write(6,*)' **-1- : ftol reached           *****' 
        if(icode.eq.2) write(6,*)' **-2- : xtol reached           *****' 
        if(icode.eq.3) write(6,*)' **-3- : ftol and xtol reached  *****' 
        if(icode.eq.4) write(6,*)' **-4- : fvec orthogonal to jac *****' 
        if(icode.eq.5) write(6,*)' **-5- : max no. fncalls reached ****' 
        if(icode.eq.6) write(6,*)' **-6- : ftol too small          ****' 
        if(icode.eq.7) write(6,*)' **-7- : xtol too small          ****' 
        if(icode.eq.8) write(6,*)' **-8- : gtol too small          ****' 
                                                                        
        return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       function iercd() 
!      returns a dummy error code of the fitting routine.               
!      The real error code is                                           
!      plotted by the modified fitting routine itself.                  
       iercd = 0 
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       subroutine linrg(n, gmat, mfit, ginv, mfit1) 
!      call lapack library instead of IMSL:                             
       dimension gmat(mfit,mfit), ginv(mfit,mfit) 
       integer inf111, mfit, mfit1, n, lwork 
       integer iwa(mfit) 
                                                                        
       lwork=n 
       call sgetrf(n,n,gmat, mfit, iwa, inf111) 
       call sgetri(n, gmat, mfit, iwa, ginv, lwork, inf111) 
       do i=1,n 
         do l=1,n 
           ginv(i,l) = gmat(i,l) 
         enddo 
       enddo 
                                                                        
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       function bsj0(x) 
!      imsl bessel fn                                                   
       real x 
       bsj0 = besj0(x) 
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       function s17aef(x,ifail) 
!      nag bessel fn                                                    
       real*8 s17aef, x , dbesj0
       integer ifail 
       external dbesj0 
                                                                        
       ifail = 0 
       s17aef = dbesj0(x) 
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function s14aaf(x,ifail) 
!      nag gamma fn                                                    
       real*8 s14aaf, x 
       integer ifail 
       external dgamma 
                                                                        
       ifail = 0 
       s14aaf = dgamma(x) 
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function s13aaf(x,ifail) 
!      nag exponential integral fn                                                
       real*8 s13aaf, x 
       integer ifail 
       external de1 
                                                                        
       ifail = 0 
       s13aaf = de1(x) 
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function s15adf(x,ifail) 
!      nag compl. err. fn erfc(x)                                                    
       real*8 s15adf, x 
       integer ifail 
       external derfc 
                                                                        
       ifail = 0 
       s15adf = derfc(x) 
       return 
      END                                           
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       subroutine s18def(fnu, z, n, scale, cy, nz, ifail) 
!      nag modified bessel function
                                                    
       real*8 fnu 
       integer ifail, n, nz, kode 
       complex z, cy
       character*1 scale
       dimension cy(n)
       external cbesi  
                                                                        
       ifail = 0
       kode=2
       if ( scale == 'U' ) kode=1
       call cbesi( z, fnu, kode, n, cy, nz, ifail ) 
       return 
       END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       subroutine s17def(fnu, z, n, scale, cy, nz, ifail) 
!      nag bessel function
                                                    
       real*8 fnu 
       integer ifail, n, nz, kode 
       complex z, cy
       character*1 scale
       dimension cy(n)
       external cbesi  
                                                                        
       ifail = 0
       kode=2
       if (( scale == 'U' ).or.(scale == 'u')) kode=1
       call cbesj( z, fnu, kode, n, cy, nz, ifail ) 
       return 
       END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       function dbsj1(x) 
!      IMSL bessel fn                                                   
       real*8 dbsj1, x 
       external dbesj1 
                                                                        
       dbsj1 = dbesj0(x) 
       return 
      END                                           
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       subroutine fftcf(N, CA, CB) 
!      call fftpack5 fft instead of IMSL:                               
       parameter(LENC=1024, LENSAV=5000, LENWRK=5000) 
       integer IER,INC 
       real    WSAVE(LENSAV), WORK(LENWRK) 
       complex C(LENC), CA(LENC), CB(LENC) 
                                                                        
       INC=1 
       C = CA 
                                                                        
       call cfft1i(N, WSAVE, LENSAV, IER) 
       call cfft1f(N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK, IER) 
                                                                        
!      CB = C                                                           
!      The following line ensures the same behaviour as the IMSL routine
       CB = C*N 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       subroutine fftcb(N, CA, CB) 
!      call fftpack5 fft instead of IMSL:                               
       parameter(LENC=1024, LENSAV=5000, LENWRK=5000) 
       integer IER,INC 
       real    WSAVE(LENSAV), WORK(LENWRK) 
       complex C(LENC), CA(LENC), CB(LENC) 
                                                                        
       INC=1 
       C = CA 
                                                                        
       call cfft1i(N, WSAVE, LENSAV, IER) 
       call cfft1b(N, INC, C, LENC, WSAVE, LENSAV,WORK, LENWRK, IER) 
                                                                        
       CB = C 
                                                                        
                                                                        
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       subroutine fft2d(NRO, NCO, CA, LDA, CB, LDCOEF ) 
!      call fftpack5 fft instead of IMSL:                               
       parameter(LENC=1025, LENSAV=5000, LENWRK=5000) 
       integer IER 
       real    WSAVE(LENSAV), WORK(LENWRK) 
       complex C(LENC,LENC), CA(LENC,LENC), CB(LENC,LENC) 
                                                                        
       C = CA 
                                                                        
       call cfft2i(NRO, NCO, WSAVE, LENSAV, IER) 
       call cfft2f(LDA,NRO,NCO,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER) 
                                                                        
       write(*,*)'Normalization of 2d-fft not yet tested properly...' 
                                                                        
       CB = C 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       subroutine fft2b(NRO, NCO, CA, LDA, CB, LDCOEF ) 
!      call fftpack5 fft instead of IMSL:                               
       parameter(LENC=1025, LENSAV=5000, LENWRK=5000) 
       integer IER 
       real    WSAVE(LENSAV), WORK(LENWRK) 
       complex C(LENC,LENC), CA(LENC,LENC), CB(LENC,LENC) 
                                                                        
       C = CA 
                                                                        
       call cfft2i(NRO, NCO, WSAVE, LENSAV, IER) 
       call cfft2b(LDA,NRO,NCO,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER) 
                                                                        
       CB = C 
                                                                        
                                                                        
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       subroutine fcost(N, CA, CB) 
!      call fftpack5 fft instead of IMSL:                               
       parameter(LENC=1024, LENSAV=5000, LENWRK=5000) 
       integer IER,INC 
       real    WSAVE(LENSAV), WORK(LENWRK) 
       real R(LENC), CA(LENC), CB(LENC) 
                                                                        
       INC=1 
       R = CA 
                                                                        
       call cost1i(N, WSAVE, LENSAV, IER) 
       call cost1f(N, INC, R, LENC, WSAVE, LENSAV,WORK, LENWRK, IER) 
                                                                        
       write(*,*)'Normalization of 2d-fft not yet tested properly...' 
       CB = R 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function d01ahf(gunten,goben,errabs,neval,errel,fdes,limit,ifail) 
!      parameter(lenw=20000) 
      integer ier, neval, last, limit, iwork(limit), ifail, lenw 
      real*8 fdes
      external fdes 
      real*8 gunten, goben, errabs, errel,errret, work(limit*4), result 
      real*8 d01ahf
      integer lim1
                                                                         
!      errel = 0.1
!      neval = 1000
!      lim1 = 5000
!      errabs = 10
      lenw = limit * 4
      call dqags(fdes,gunten,goben,errabs,errel,result,errret,neval,    &
     &          ier, limit, lenw, last , iwork, work ) 
      d01ahf = result
      return 
      END                                           
                                                                       
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
! qdng  => qng  from quadpack                                           
! qdags => qags from quadpack                                           
                                                                        
      subroutine qdng(fdes,gunten,goben,errabs,errel,result,errret) 
      integer ier, neval 
      external fdes 
      real gunten, goben, errabs, errel,errret 
                                                                        
      write(*,*)'Test the subroutine qng!!!!' 
      call qng(fdes,gunten,goben,errabs,errel,result,errret,neval,ier) 
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      subroutine qdags(fdes,gunten,goben,errabs,errel,result,errret) 
      parameter(lenw=2000, limit=500) 
      integer ier, neval, last, iwork(limit) 
      external fdes 
      real gunten, goben, errabs, errel,errret, work(lenw) 
                                                                        
      write(*,*)'Test the subroutine qags!!!!' 
      call qags(fdes,gunten,goben,errabs,errel,result,errret,neval,ier, &
     &          limit, lenw, last , iwork, work )                       
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      subroutine qdag(fdes,gunten,goben,errabs,errel,                   &
     &                irule, result,errret)                             
      parameter(lenw=2000, limit=500) 
      integer ier, neval, last, iwork(limit), irule 
      external fdes 
      real gunten, goben, errabs, errel,errret, work(lenw) 
                                                                        
      write(*,*)'Test the subroutine qag!!!!' 
      call qag(fdes,gunten,goben,errabs,errrel, key,result,errret,neval,&
     &          limit, lenw, last , iwork, work )                       
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       subroutine qand(f,c,alim,blim,erraec,errrec,                     &
     &                 maxint,result, errest)                           
!      calculates the echo form                                         
       integer minpts, maxpts, key,restar, ifail, neval, c 
       parameter (numfun=1, nw=5000) 
       integer abserr(numfun) 
       real work(nw), resu(numfun) 
       common/partran/j1echo,j2echo,j0delta,cdelta 
       common/wlntran/alam0, dalam 
       common/sqtran/tau 
       external f, funsub 
                                                                        
                                                                        
       minpts=0 
       maxpts=maxint 
       key=0 
       restar=0 
       write(*,*)'ONLY for echo production!!!!' 
       call scuhre(c,numfun,alim,blim,minpts, maxpts, funsub,           &
     &            erraec,errrec,key,                                    &
     &            nw, restar, resu, abserr , neval , ifail ,work)       
       result=resu(1) 
       errest=abserr(1) 
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine funsub(c, x, numfun, funvls) 
       integer c, numfun 
       dimension x(c), funvls(numfun) 
       external fecho 
       common/partran/j1echo,j2echo,j0delta,cdelta 
       common/wlntran/alam0, dalam 
       common/sqtran/tau 
                                                                        
       funvls=fecho(c,x) 
!       funvls=1.0                                                      
!       write(6,*)'funsub: ',funvls                                     
       return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       function fecho(n,x) 
!      ---------------                                                  
! --- integrand ---                                                     
       implicit real*4 (a-h,o-z) 
       real*4 j1echo,j2echo,j0delta 
       common/partran/j1echo,j2echo,j0delta,cdelta 
       dimension x(n) 
!                                                                       
!      x(1) = lambda                                                    
!      x(2) = omega                                                     
!                                                                       
       data gamma/18303.33e0/ 
! ---- larmorkonstante in rad/s/gauss                                   
       data xmh  /2.50607e-4/ 
! ---- neutronenmasse durch h in 2/m/angstroem                          
       data zpi  /6.283185307e0/ 
!       write(*,*)j1echo,j2echo,j0delta,cdelta                          
       a  = gamma*xmh*x(1)/(1.d0+xmh*x(1)*x(1)*x(2)*1d-10/zpi) 
       b  = gamma*xmh*x(1) 
       del= j0delta**2 + ( (j1echo+j2echo)*0.5e0*cdelta )**2 
       fecho  = ( 1.d0 +  exp(-(a**2+b**2)*del*0.25e0 )*                &
     &               cos( a*j1echo - b*j2echo )                         &
     &      ) * sofqom1(x(2)) * wlambda1(x(1))                          
       return 
      END                                           
                                                                        
       function wlambda1( alam ) 
!      ------------------------                                         
! --- repraesentiert die wellenlaengenverteilung                        
       implicit real*4 (a-h,o-z) 
       common/wlntran/alam0, dalam 
          arg     = ( (alam-alam0)/dalam )**2 
          if(arg.lt.50.e0) then 
            wlambda1 =  exp( -arg ) 
          else 
            wlambda1 = 0.e0 
          endif 
       return 
      END                                           
       function sofqom1( omega ) 
!      ------------------------                                         
! --- repraesentiert die streufunktion, omega in s**-1                  
       implicit real*4 (a-h,o-z) 
       common/sqtran/tau 
!      tau in sekunden                                                  
                                                                        
       x      = omega * tau 
       sofqom1 = tau /(1.e0+x*x) 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      subroutine erset(a,b,c) 
      write(*,*)'subroutine erset not needed in this version...' 
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      subroutine rnset(iseed) 
                                                                        
      write(*,*)'rnset: not extensively tested...' 
      CALL RANDOM_SEED 
                                                                        
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function rnunf() 
      real R 
      CALL RANDOM_NUMBER(R) 
      rnunf=R 
                                                                        
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       function erfi(X) 
!      calculates the inverse error function from tabulated values...   
       REAL X,b 
       integer i 
       dimension xer(100), yer(100) 
!      Tabulated inverse error function:                                
       data xer / 0.0E0,0.01E0,0.02E0,0.03E0,0.04E0,0.05E0,0.06E0,      &
     &     0.07E0,0.08E0,0.09E0,0.1E0,0.11E0,0.12E0,0.13E0,             &
     &     0.14E0,0.15E0,0.16E0,0.17E0,0.18E0,0.19E0,0.2E0,             &
     &     0.21E0, 0.22E0, 0.23E0, 0.24E0, 0.25E0,  0.26E0,  0.27E0,    &
     &     0.28E0, 0.29E0, 0.3E0, 0.31E0, 0.32E0,  0.33E0,  0.34E0,     &
     &     0.35E0, 0.36E0, 0.37E0, 0.38E0, 0.39E0,  0.4E0,  0.41E0,     &
     &     0.42E0, 0.43E0, 0.44E0, 0.45E0, 0.46E0,  0.47E0,  0.48E0,    &
     &     0.49E0, 0.5E0, 0.51E0, 0.52E0, 0.53E0,  0.54E0,  0.55E0,     &
     &     0.56E0, 0.57E0, 0.58E0, 0.59E0, 0.6E0,  0.61E0,  0.62E0,     &
     &     0.63E0, 0.64E0, 0.65E0, 0.66E0, 0.67E0,  0.68E0,  0.69E0,    &
     &     0.7E0, 0.71E0, 0.72E0, 0.73E0, 0.74E0,  0.75E0,  0.76E0,     &
     &     0.77E0, 0.78E0, 0.79E0, 0.8E0, 0.81E0,  0.82E0,  0.83E0,     &
     &     0.84E0, 0.85E0, 0.86E0, 0.87E0, 0.88E0,  0.89E0,  0.9E0,     &
     &     0.91E0, 0.92E0, 0.93E0, 0.94E0, 0.95E0,  0.96E0,  0.97E0,    &
     &     0.98E0, 0.99E0/                                              
       data  yer / 0.0, 0.0088625, 0.0177264,                           &
     &     0.02659308,0.03546394,0.04434039,0.05322383,                 &
     &     0.06211568,  0.07101736,  0.07993032,  0.08885599,           &
     &        0.09779584,  0.10675136,  0.11572403,  0.12471537,        &
     &     0.13372692,  0.14276025,  0.15181693,  0.1608986 ,           &
     &     0.17000688,  0.17914345,  0.18831004,  0.19750838,           &
     &     0.20674027,  0.21600754,  0.22531206,  0.23465575,           &
     &     0.2440406,  0.25346864,  0.26294196,  0.27246271,            &
     &     0.28203312,  0.29165548,  0.30133215,  0.31106558,           &
     &     0.32085832,  0.33071299,  0.34063232,  0.35061914,           &
     &     0.3606764,  0.37080716,  0.38101461,  0.39130209,            &
     &     0.40167307,  0.41213118,  0.42268024,  0.43332422,           &
     &     0.44406731,  0.4549139,  0.4658686,  0.47693628,             &
     &     0.48812205,  0.49943133,  0.51086984,  0.52244362,           &
     &     0.53415909,  0.54602306,  0.55804277,  0.57022593,           &
     &     0.58258077,  0.59511608,  0.60784127,  0.62076642,           &
     &     0.63390239,  0.64726086,  0.66085443,  0.67469672,           &
     &     0.68880253,  0.70318791,  0.71787038,  0.73286908,           &
     &     0.74820497,  0.76390113,  0.77998301,  0.79647881,           &
     &     0.81341985,  0.83084113,  0.84878189,  0.86728635,           &
     &     0.88640462,  0.9061938,  0.92671938,  0.94805698,            &
     &     0.97029462,  0.99353563,  1.01790246,  1.04354184,           &
     &     1.07063171,  1.09939095,  1.13009321,  1.16308715,           &
     &     1.19882722,  1.23792199,  1.28121432,  1.32992191,           &
     &    1.38590382,  1.45221978,  1.53448562,  1.64497636,            &
     &    1.82138637/                                                   
                                                                        
       if (X.ge.(1.0)) then 
          erfi=5.0 
          write(6,*),'Inverse Error Function out of bounds!!!!!!!' 
       endif 
       loop:do i=1,100 
          if(xer(i).ge.X) then 
            b=(yer(i)-yer(i-1))/(xer(i)-xer(i-1))*(X-xer(i-1)) 
            erfi=yer(i-1)+b 
            exit 
          endif 
       enddo loop 
       return 
      END                                           
                                                                        
! Tabulated inverse error function from python: scipy.special.erfinv(x):
                                                                        
!x= array([ 0.  ,  0.01,  0.02,  0.03,  0.04,  0.05,  0.06,  0.07,  0.08
!             0.11,  0.12,  0.13,  0.14,  0.15,  0.16,  0.17,  0.18,  0.
!             0.21,  0.22,  0.23,  0.24,  0.25,  0.26,  0.27,  0.28,  0.
!             0.31,  0.32,  0.33,  0.34,  0.35,  0.36,  0.37,  0.38,  0.
!             0.41,  0.42,  0.43,  0.44,  0.45,  0.46,  0.47,  0.48,  0.
!             0.51,  0.52,  0.53,  0.54,  0.55,  0.56,  0.57,  0.58,  0.
!             0.61,  0.62,  0.63,  0.64,  0.65,  0.66,  0.67,  0.68,  0.
!             0.71,  0.72,  0.73,  0.74,  0.75,  0.76,  0.77,  0.78,  0.
!             0.81,  0.82,  0.83,  0.84,  0.85,  0.86,  0.87,  0.88,  0.
!             0.91,  0.92,  0.93,  0.94,  0.95,  0.96,  0.97,  0.98,  0.
                                                                        
!In [35]: y                                                             
!Out[35]:                                                               
!array([ 0.        ,  0.0088625 ,  0.0177264 ,  0.02659308,  0.03546394,
!             0.05322383,  0.06211568,  0.07101736,  0.07993032,  0.0888
!             0.09779584,  0.10675136,  0.11572403,  0.12471537,  0.1337
!             0.14276025,  0.15181693,  0.1608986 ,  0.17000688,  0.1791
!             0.18831004,  0.19750838,  0.20674027,  0.21600754,  0.2253
!             0.23465575,  0.2440406 ,  0.25346864,  0.26294196,  0.2724
!             0.28203312,  0.29165548,  0.30133215,  0.31106558,  0.3208
!             0.33071299,  0.34063232,  0.35061914,  0.3606764 ,  0.3708
!             0.38101461,  0.39130209,  0.40167307,  0.41213118,  0.4226
!             0.43332422,  0.44406731,  0.4549139 ,  0.4658686 ,  0.4769
!             0.48812205,  0.49943133,  0.51086984,  0.52244362,  0.5341
!             0.54602306,  0.55804277,  0.57022593,  0.58258077,  0.5951
!             0.60784127,  0.62076642,  0.63390239,  0.64726086,  0.6608
!             0.67469672,  0.68880253,  0.70318791,  0.71787038,  0.7328
!             0.74820497,  0.76390113,  0.77998301,  0.79647881,  0.8134
!             0.83084113,  0.84878189,  0.86728635,  0.88640462,  0.9061
!             0.92671938,  0.94805698,  0.97029462,  0.99353563,  1.0179
!             1.04354184,  1.07063171,  1.09939095,  1.13009321,  1.1630
!            1.19882722,  1.23792199,  1.28121432,  1.32992191,  1.38590
!            1.45221978,  1.53448562,  1.64497636,  1.82138637])        
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!       subroutine csscv(nn,x,y,iequal,break,cscoef)                    
!       integer md, val, nc, ier, m, k                                  
!       real wk(nn*13+6), c(nn*4,1)                                     
!       dimension weight(nn), wy(1), x(nn),y(nn), break(nn), cscoef(4,nn
!       data wy/1D0/                                                    
!                                                                       
!       md=2                                                            
!       val=1                                                           
!       nc=nn*4                                                         
!       m=2                                                             
!       k=1                                                             
!       do i=1,nn                                                       
!           weight(i)=1.0                                               
!       enddo                                                           
!       write(*,*)'cubic spline smoothing,'                             
!       write(*,*)'estimated smoothing parameters...'                   
!       write(*,*)'not properly tested without nag or IMSL...'          
!       write(*,*)'equivalent nag routine: G10ACF'                      
!       call GCVSPL (x,y,nn,weight,wy,m,nn,k,md,val,c, nc, wk, ier)     
! c      The calculation of the smoothed fn is done by datreat...       
! c      value at point i:                                              
! c      Q(i) = SPLDER ( 0, m, nn, x(i), x, c, J, wk )                  
! c      here: break=x (?)                                              
!        do i=1,nn                                                      
!            break(i)=x(i)                                              
!        enddo                                                          
!        k=0                                                            
!        do i=1,4                                                       
!            do j=1,nn                                                  
!                k=k+1                                                  
!                cscoef(i,j)=c(k)                                       
!            enddo                                                      
!        enddo                                                          
!                                                                       
!       return                                                          
!       end                                                             
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine csscv(nn,x,y,iequal,brk,cscoef) 
      integer ier, m, k, cn 
      real wk(nn*13+6), c0(nn), c(nn-1,3) 
      dimension weight(nn), wy(1),x(nn),y(nn), brk(nn), cscoef(4,nn) 
                                                                        
      cn=nn-1 
      write(*,*)'cubic spline smoothing,' 
      write(*,*)'estimated smoothing parameters...' 
      write(*,*)'not properly tested without nag or IMSL...' 
      write(*,*)'equivalent nag routine: G10ACF' 
      job=0 
      var=-1 
!      call GCVSPL (x,y,nn,weight,wy,m,nn,k,md,val,c, nc, wk, ier)      
      call CUBGCV (dble(x),dble(y),dble(weight),nn,dble(c0),dble(c),    &
     &cn,dble(var),job,dble(se),dble(wk),ier) 
		 do i=1,nn 
           brk(i)=x(i) 
           cscoef(1,i)=c0(i) 
       enddo 
       k=0 
       do i=1,3 
           do j=1,nn-1 
               k=k+1 
               cscoef(i+1,j)=c(j,i) 
           enddo 
       enddo 
                                                                        
      return 
      END                                           
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      subroutine cssmh(nn,x,y,weight,smpar,brk,cscoef) 
      integer ier, m, k, cn 
      real wk(nn*13+6), c0(nn), c(nn-1,3) 
      dimension weight(nn), wy(1),x(nn),y(nn), brk(nn), cscoef(4,nn) 
                                                                        
      cn=nn-1 
      write(*,*)'cubic spline smoothing,' 
      write(*,*)'ATTENTION: is weight parameter used correctly' 
      write(*,*)'in current routine (smoothing seems to work  ' 
      write(*,*)'properly...)' 
      write(*,*)'equivalent nag routine: G10ABF' 
      job=0 
      var=-1 
!      call GCVSPL (x,y,nn,weight,wy,m,nn,k,md,val,c, nc, wk, ier)      
      call CUBGCV (x,y,weight,nn,c0,c,cn,var,job,se,wk,ier) 
       do i=1,nn 
           brk(i)=x(i) 
           cscoef(1,i)=c0(i) 
       enddo 
       k=0 
       do i=1,3 
           do j=1,nn-1 
               k=k+1 
               cscoef(i+1,j)=c(j,i) 
           enddo 
       enddo 
                                                                        
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      function BSI0E(x) 
                                                                        
      write(*,*)'IMSL function BSI0E replaced by slatec/besi0e.f...' 
      BSI0E=besi0e(x) 
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
      subroutine u4lsf(a,b) 
                                                                        
      write(*,*)'IMSL routine u4lsf called,' 
      write(*,*)'not used with minpack...' 
                                                                        
      return 
      END                                           
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
      subroutine dsvrgp(a,b,c,d) 
                                                                        
      write(*,*)'IMSL routine dsvrgp called,' 
      write(*,*)'not replaced yet...' 
      write(*,*)'ERROR!!!!!!!!' 
                                                                        
                                                                        
      return 
      END                                           
                                                                        
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
      subroutine dpermu(a,b,c,d,e) 
                                                                        
      write(*,*)'IMSL routine dpermu called,' 
      write(*,*)'not replaced yet...' 
      write(*,*)'ERROR!!!!!!!!' 
                                                                        
      return 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
      function Dsi(x) 
! where is the function in use????                                      
                                                                        
      write(*,*)'Function Dsi called...' 
      write(*,*)'not yet replaced!!!!!!!!' 
      write(*,*)'ERROR!!!!!!!!' 
      Dsi=0 
      return 
      END          
                                 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      function dlngam(x)
!  Slatec-Funktion gleichen Namens eingefgt
!      write(*,*)'Function dlngam called...'
!      write(*,*)'not yet replaced!!!!!!!!'
!      write(*,*)'ERROR!!!!!!!!'
!      dlngam=alngam
!      return
!      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

