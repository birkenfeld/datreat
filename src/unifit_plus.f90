      PROGRAM unift 
                                                                        
!     Universal Fourier Transformer                                     
!     Version 3.2                                                       
!                                                                       
!     This program Fourier transforms TOF and BS data with INX format.  
!     Constant-Q-ing should be performed prior to this program by INTERP
!     or SIGRID. Special features:                                      
!                                                                       
!       1. non-equidistant E scale (no FFT used)                        
!       2. "mirroring" of negative Es                                   
!       3. optional energy range restriction                            
!       4. symmetrization of S(Q,E) (classical)                         
!       5. time step adjusted to E range individually                   
!       6. error propagation (v. 2.x)                                   
!       7. maximum number of t-values -> 3000                           
!       8. optimal time step can be overridden                          
!       9. linear interpolation of S(Q,E) data (v. 3.x)                 
!      10. averaging of spectra (only for identical E scales)           
!                                                                       
!                                                                       
!     Item 2 may be realized more correctly.                            
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z) 
                                                                        
      PARAMETER (NE=2048,NQ=128) 
!     Commons for CRIN use:                                             
      REAL*4 SPECPAR,ANGLES,NUMSPE 
      REAL*4 X,Y,YE 
      COMMON/SAMPL1/SPECPAR(30),ANGLES,NUMSPE 
      COMMON/SAMPL2/X(NE),Y(NE),YE(NE) 
      CHARACTER*80 FILIN 
      CHARACTER*40 TITLIN 
      PARAMETER (mt=3000) 
!     Arrays and variables for UNIFT use:                               
      DIMENSION es(ne),ss(ne),ds(ne) 
      DIMENSION ebegs,jbegs,eends,jends 
      DIMENSION es1(ne),ss1(ne),ds1(ne),ns1 
      DIMENSION er(ne),sr(ne),dr(ne) 
      DIMENSION ebegr,jbegr,eendr,jendr 
      DIMENSION er1(ne),sr1(ne),dr1(ne),nr1 
      DIMENSION ndats,ndatr 
      DIMENSION fts(mt),ftr(mt) 
      DIMENSION dfts(mt),dftr(mt),dsqt(mt) 
      CHARACTER asym,ans 
      CHARACTER*80 parse,inp1,ofnam,inp2 
                                                                        
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,*) '*** UNIFT Version 3.2 ***' 
      WRITE(6,*) 'Fourier transform INX files' 
      WRITE(6,*) 
      WRITE(6,*) 
                                                                        
                                                                        
!     Sort to ascending order:                                          
                                                                        
      DO i=1,nqs 
        IF (eends(i).gt.ebegs(i)) THEN 
          DO j=1,ndats(i) 
            es(j,i)=x(j+jbegs(i)-1,i) 
            ss(j,i)=y(j+jbegs(i)-1,i) 
            ds(j,i)=ye(j+jbegs(i)-1,i) 
          END DO 
        ELSE 
          DO j=1,ndats(i) 
            es(j,i)=x(jends(i)-j+1,i) 
            ss(j,i)=y(jends(i)-j+1,i) 
            ds(j,i)=ye(jends(i)-j+1,i) 
          END DO 
        END IF 
      END DO 
      do i=1,nqs 
        do j=1,ndats(i) 
          write(16,'(f9.5,2e13.5)') es(j,i),ss(j,i),ds(j,i) 
        end do 
        write(16,'(a)') '-----------------------------------' 
      end do 
                                                                        
!     Read the resolution file:                                         
                                                                        
                                                                        
!     Sort to ascending order:                                          
                                                                        
      DO i=1,nqr 
        IF (eendr(i).gt.ebegr(i)) THEN 
          DO j=1,ndatr(i) 
            er(j,i)=x(j+jbegr(i)-1,i) 
            sr(j,i)=y(j+jbegr(i)-1,i) 
            dr(j,i)=ye(j+jbegr(i)-1,i) 
          END DO 
        ELSE 
          DO j=1,ndatr(i) 
            er(j,i)=x(jendr(i)-j+1,i) 
            sr(j,i)=y(jendr(i)-j+1,i) 
            dr(j,i)=ye(jendr(i)-j+1,i) 
          END DO 
        END IF 
      END DO 
                                                                        
!     Check sample/resolution consistency:                              
                                                                        
      WRITE(6,*) 
      IF (nqs.NE.nqr) THEN 
        WRITE(6,*) 'Warning: Different numbers of Q values,' 
        WRITE(6,*) '         select corresponding detectors manually!' 
      ELSE 
        DO i=1,nqs 
          IF (ebegs(i).NE.ebegr(i)) GOTO 20 
          IF (eends(i).NE.eendr(i)) GOTO 20 
        END DO 
        GOTO 21 
   20   WRITE(6,*) 'Warning: Different energy ranges for some Qs,' 
        WRITE(6,*) '         restrict energy range manually!' 
   21 END IF 
                                                                        
!     Start treating a Q value:                                         
                                                                        
                                                                        
!     Ask energy limit for Fourier transform:                           
                                                                        
   50 WRITE(6,1011) 
 1011 FORMAT('Limit for E range? (>0, =0: no limit) ',$) 
      READ(5,*) emir 
      IF (emir.eq.0) emir=1.0e30 
                                                                        
!     Calculate range of sample data effectively used:                  
                                                                        
      DO j=1,ndats(iqs) 
        IF (abs(es(j,iqs)).le.emir) GOTO 45 
      END DO 
   45 jbegse=j 
      DO j=jbegse,ndats(iqs) 
        IF (abs(es(j,iqs)).gt.emir) GOTO 46 
      END DO 
   46 jendse=j-1 
      write(6,1021) jbegse,jendse,es(jbegse,iqs),es(jendse,iqs) 
 1021 FORMAT('Effective limits for sample:    ',                        &
     &       i4,'...',i4,f10.3,'...',f8.3)                              &
                                                                        
!     Calculate range of resolution data effectively used:              
                                                                        
      DO j=1,ndatr(iqr) 
        IF (abs(er(j,iqr)).le.emir) GOTO 40 
      END DO 
   40 jbegre=j 
      DO j=jbegre,ndatr(iqr) 
        IF (abs(er(j,iqr)).gt.emir) GOTO 41 
      END DO 
   41 jendre=j-1 
      write(6,1020) jbegre,jendre,er(jbegre,iqr),er(jendre,iqr) 
 1020 FORMAT('Effective limits for resolution:',                        &
     &       i4,'...',i4,f10.3,'...',f8.3)                              &
                                                                        
      emaxr=max(abs(er(jbegre,iqr)),abs(er(jendre,iqr))) 
      emaxs=max(abs(es(jbegse,iqs)),abs(es(jendse,iqs))) 
      dt=2.067834/((emaxs+emaxr)*0.5)  !!!>>>>!!! [ h/eunit/2 / tunit] /  ((emaxs+emaxr)*0.5) !! here eunit=meV, tunit=ps   
      WRITE(6,1014) dt 
 1014 FORMAT('Timestep calculated from energy ranges:',f10.3,' ps') 
                                                                        
                                                                        
!     Ask limit for mirroring:                                          
!     (The convention is: emir positive, but limit is -emir.)           
                                                                        
   55 WRITE(6,1012) 
 1012 FORMAT('Limit for E mirroring? (0: no mirroring) ',$) 
      READ(5,*) emir 
      emir=abs(emir) 
                                                                        
!     Ask for symmetrization temperature:                               
                                                                        
      WRITE(6,1013) 
 1013 FORMAT('Temperature for symmetrization? (0: no symm.) ',$) 
      READ(5,*) temp 
      IF (temp.EQ.0.0) THEN 
        fdbf=0.0 
      ELSE 
        fdbf=5.80222/temp       !!!>>>>!!! [eunit/2/k_B] / temp  (here eunit=meV)
      END IF 
                                                                        
!     Mirroring, Symmetrisation, and conversion E->omega:               
                                                                        
      call mirror(iqs,es,ss,ds,jbegse,jendse,fdbf,emir,                 &
     &            es1,ss1,ds1,ns1)                                      &
      call mirror(iqr,er,sr,dr,jbegre,jendre,0.0d0,emir,                &
     &            er1,sr1,dr1,nr1)                                      &
                                                                        
!     Ask for time limit:                                               
                                                                        
      WRITE(6,1029) 
 1029 FORMAT('Maximum time to calculate (ps)?  ',$) 
      READ(5,*) tmax 
      nt=int(tmax/dt) 
                                                                        
!     Ask for output file name:                                         
                                                                        
      WRITE(6,1030) 
 1030 FORMAT('Output file name? (w/o extension) ',$) 
      READ(5,'(a20)') ofnam 
      DO ich=1,20 
        IF (ofnam(ich:ich).EQ.' ') GOTO 60 
      END DO 
   60 nchofn=ich-1 
                                                                        
!     Outer loop over times (sample):                                   
                                                                        
      WRITE(6,*) 
      WRITE(6,*) 'Transforming sample data:' 
      WRITE(6,*)                                                        &
     & '       t [ps]      real      imag      ',                       &
     & ' abs       arg      norm      rerr'                             &
      DO it=1,nt+1 
        t=(it-1)*dt 
                                                                        
!     Fourier transform the sample data:                                
                                                                        
        call fourier(es1,ss1,ds1,iqs,ns1(iqs),t,                        &
     &               ft1,ft2,fts(it),dfts(it))                          &
                                                                        
        WRITE(6,1040) it,t,ft1,ft2,fts(it),ft2/ft1,                     &
     &                fts(it)/fts(1),dfts(it)/fts(it)                   &
 1040   FORMAT(i4,7f10.3) 
                                                                        
      END DO 
                                                                        
!     Write Fourier transformed SAMPLE DATA to file:                    
                                                                        
      OPEN(unit=41,file=ofnam(1:nchofn)//'.sam',status='unknown') 
      DO it=1,nt+1 
        t=(it-1)*dt 
        WRITE(41,1045) t,fts(it),dfts(it) 
 1045   FORMAT(3f12.4) 
      END DO 
      CLOSE(41) 
                                                                        
!     Outer loop over times (resolution):                               
                                                                        
      WRITE(6,*) 
      WRITE(6,*) 'Transforming resolution data:' 
      WRITE(6,*)                                                        &
     & '       t [ps]      real      imag      ',                       &
     & ' abs       arg      norm      rerr'                             &
      DO it=1,nt+1 
        t=(it-1)*dt 
                                                                        
!     Fourier transform the resolution data:                            
                                                                        
        call fourier(er1,sr1,dr1,iqr,nr1(iqr),t,                        &
     &               ft1,ft2,ftr(it),dftr(it))                          &
                                                                        
        WRITE(6,1040) it,t,ft1,ft2,ftr(it),ft2/ft1,                     &
     &                ftr(it)/ftr(1),dftr(it)/ftr(it)                   &
                                                                        
      END DO 
                                                                        
!     Write Fourier transformed RESOLUTION DATA to file:                
                                                                        
      OPEN(unit=42,file=ofnam(1:nchofn)//'.res',status='unknown') 
      DO it=1,nt+1 
        t=(it-1)*dt 
        WRITE(42,1045) t,ftr(it),dftr(it) 
      END DO 
      CLOSE(42) 
                                                                        
!     Write DECONVOLUTED S(Q,t) DATA to file:                           
                                                                        
      OPEN(unit=43,file=ofnam(1:nchofn)//'.sqt',status='unknown') 
      DO it=1,nt+1 
        t=(it-1)*dt 
        sqt=fts(it)/ftr(it) 
        sqterr=sqt*sqrt((dftr(it)/ftr(it))**2+(dfts(it)/fts(it))**2) 
        WRITE(43,1045) t,sqt,sqterr 
      END DO 
      CLOSE(43) 
                                                                        
!     Play it again, Sam...                                             
                                                                        
      WRITE(6,*) 
      GOTO 30 
                                                                        
!     This is the end, my only friend, the end...                       
                                                                        
  999 END                                           
                                                                        
                                                                        
      SUBROUTINE mirror(iq,es,ss,ds,jbeg,jend,fdbf,emir,                &
     &                  es1,ss1,ds1,ns1)                                &
                                                                        
!     This subroutine mirrors the high energy part of the spectrum      
!     to complement the energy loss side. At the same time it transforms
!     from S(Q,E) to S(Q,omega).                                        
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z) 
                                                                        
      PARAMETER (NE=2048,NQ=128) 
      DIMENSION es(ne),ss(ne),ds(ne) 
      DIMENSION es1(ne),ss1(ne),ds1(ne)) 
                                                                        
      if (es(jbeg).gt.-emir) then 
        if (es(jbeg).gt.0.0) then 
          write(6,'(a,i3)') 'Fatal: No elastic line at #',iq 
          ns1(iq)=0 
          goto 900 
        else 
          write(6,'(a,i3)') 'Warning: Cut elastic line at #',iq 
        end if 
      end if 
                                                                        
      j1=1 
                                                                        
      if (emir.le.0.0) then 
        jlim=jbeg 
        goto 120 
      end if 
                                                                        
      do je=jend,jbeg,-1 
        ee=es(je) 
        if (ee.le.emir) goto 100 
        es1(j1)=-1.519257*ee                !!! >>> !!! -[2*pi*eunit/h * tunit] *ee 
        dbf=exp(fdbf*ee) 
        ss1(j1)=dbf*ss(je)*0.658216      !!! >>> !!! dbf*ss(je)*[1/(2*pi*eunit/h * tunit)]
        ds1(j1)=dbf*ds(je)*0.658216 
        j1=j1+1 
      end do 
      stop 
                                                                        
  100 do je=jbeg,jend 
        if (es(je).ge.-emir) goto 110 
      end do 
      stop 
                                                                        
  110 jlim=je 
  120 do je=jlim,jend 
        ee=es(je) 
        es1(j1)=1.519257*ee 
        dbf=exp(fdbf*ee) 
        ss1(j1)=dbf*ss(je)*0.658216 
        ds1(j1)=dbf*ds(je)*0.658216 
        j1=j1+1 
      end do 
                                                                        
      ns1(iq)=j1-1 
                                                                        
  900 return 
      END                                           
                                                                        
                                                                        
      SUBROUTINE fourier(w,s,ds,iq,n,t,ft1,ft2,ft,dft) 
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z) 
                                                                        
      PARAMETER (NE=2048,NQ=128) 
      DIMENSION w(ne),s(ne),ds(ne) 
      DIMENSION a1(ne),a2(ne),swt(ne),cwt(ne) 
                                                                        
      if (t.ne.0.0) then 
                                                                        
        do i=1,n 
          swt(i)=sin(w(i)*t) 
          cwt(i)=cos(w(i)*t) 
        end do 
        t2=1.0/(t*t) 
                                                                        
        a1(1)=-swt(1)/t+(cwt(1)-cwt(2))/(w(2)-w(1))*t2 
        a2(1)=cwt(1)/t+(swt(1)-swt(2))/(w(2)-w(1))*t2 
        do i=2,n-1 
          a1(i)=((cwt(i)-cwt(i+1))/(w(i+1)-w(i))                  &
     &          +(cwt(i-1)-cwt(i))/(w(i-1)-w(i)))*t2              &
          a2(i)=((swt(i)-swt(i+1))/(w(i+1)-w(i))                  &
     &          +(swt(i-1)-swt(i))/(w(i-1)-w(i)))*t2              &
        end do 
        a1(n)=swt(n)/t+(cwt(n-1)-cwt(n))/(w(n-1)-w(n))*t2 
        a2(n)=-cwt(n)/t+(swt(n-1)-swt(n))/(w(n-1)-w(n))*t2 
                                                                        
      else 
                                                                        
        a1(1)=(w(2)-w(1))*0.5 
        a2(1)=0.0 
        do i=2,n 
          a1(i)=(w(i+1)-w(i-1))*0.5 
          a2(i)=0.0 
        end do 
        a1(n)=(w(n)-w(n-1))*0.5 
        a2(n)=0.0 
                                                                        
      end if 
                                                                        
      ft1=0.0 
      ft2=0.0 
      do i=1,n 
        ft1=ft1+a1(i)*s(i) 
        ft2=ft2+a2(i)*s(i) 
      end do 
      ft=sqrt(ft1*ft1+ft2*ft2) 
                                                                        
      dft2=0.0 
      do i=1,n 
        dft2=dft2+((a1(i)*ft1+a2(i)*ft2)/ft*ds(i))**2 
      end do 
      dft=sqrt(dft2) 
                                                                        
      return 
      END                                           
                                                                        

*DECK DPSORT
!      SUBROUTINE DPSORT (DX, N, IPERM, KFLAG, IER)
C***BEGIN PROLOGUE  DPSORT
C***PURPOSE  Return the permutation vector generated by sorting a given
C            array and, optionally, rearrange the elements of the array.
C            The array may be sorted in increasing or decreasing order.
C            A slightly modified quicksort algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A1B, N6A2B
C***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
C***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
C***AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DPSORT returns the permutation vector IPERM generated by sorting
C   the array DX and, optionally, rearranges the values in DX.  DX may
C   be sorted in increasing or decreasing order.  A slightly modified
C   quicksort algorithm is used.
C
C   IPERM is such that DX(IPERM(I)) is the Ith value in the
C   rearrangement of DX.  IPERM may be applied to another array by
C   calling IPPERM, SPPERM, DPPERM or HPPERM.
C
C   The main difference between DPSORT and its active sorting equivalent
C   DSORT is that the data are referenced indirectly rather than
C   directly.  Therefore, DPSORT should require approximately twice as
C   long to execute as DSORT.  However, DPSORT is more general.
C
C   Description of Parameters
C      DX - input/output -- double precision array of values to be
C           sorted.  If ABS(KFLAG) = 2, then the values in DX will be
C           rearranged on output; otherwise, they are unchanged.
C      N  - input -- number of values in array DX to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the value in the original order of the
C              DX array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting DX in increasing order and sort DX also.
C            =  1  means return the permutation vector resulting from
C                  sorting DX in increasing order and do not sort DX.
C            = -1  means return the permutation vector resulting from
C                  sorting DX in decreasing order and do not sort DX.
C            = -2  means return the permutation vector resulting from
C                  sorting DX in decreasing order and sort DX also.
C      IER - output -- error indicator:
C          =  0  if no error,
C          =  1  if N is zero or negative,
C          =  2  if KFLAG is not 2, 1, -1, or -2.
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.





Module PhysConstants

!**** Physikalische Konstanten ******************************/
! Larmorkonstante       2913.06598e+4  /* Hz/Tesla    */
! Mue0                  1.25663706e-6  /* Tesla/(A/m) */
! Neutronenmasse        1.6747e-27     /* kg          */
! Planckkonstante       6.626176e-34   /* Js          */
! Erdbeschleunigung     9.80665        /* m/s^2       */
! Boltzmannkonstante    1.380662e-23   /* J/K         */
! RaDeg                 57.29577951    /* Konversion rad in grad*/

    double precision, save :: Pi                  =     3.141592654d0
    
    
    double precision, save :: Larmorkonstante     =     2913.06598d+4
    double precision, save :: Mue0                =     1.25663706d-6
    double precision, save :: Neutronenmasse      =     1.6747d-27 
    double precision, save :: Planckkonstante     =     6.626176d-34
    double precision, save :: Erdbeschleunigung   =     9.80665d0       
    double precision, save :: Boltzmannkonstante  =     1.380662d-23 
    double precision, save :: Elektronenladung    =     1.6021892d-19
    double precision, save :: Avogadrozahl        =     6.022045d+23
    double precision, save :: Molvolumen_Gas      =     0.02241383d0
    double precision, save :: Elektronenmasse     =     9.109534d-31
    double precision, save :: amu                 =     1.6605655d-27
    double precision, save :: Epsilon0            =     8.854187818d-12  ! Vakuumdielktrizitaetskonstante
    double precision, save :: Lichtgeschwindigkeit=     299792458.0d0
    double precision, save :: Gaskonstante        =     8.31441d0
    
    
    double precision, save :: unit_Angstroem      =     1d-10
    double precision, save :: unit_cm             =     1d-2
    double precision, save :: unit_eV             =     1.6021892d-19
    double precision, save :: unit_barn           =     1d-28
    
    
    double precision, save :: Elimit_classical    =     1d-17            ! Limit Neutron energy for relativistic






CONTAINS





double precision function NeutronWavelength_fromE(Energy)
!--------------------------------------------------------
    implicit none

    double precision, intent(in) :: Energy 

    if(Energy < Elimit_classical) then

      ! in non-relativistic approximation

      NeutronWavelength_fromE =  Planckkonstante /sqrt(2*Energy*Neutronenmasse)
    
!      write(6,*)'wlc =', NeutronWavelength_fromE

    else

    ! relativistic expression

       NeutronWavelength_fromE =  Planckkonstante * Lichtgeschwindigkeit / &
       sqrt((Energy+Neutronenmasse*Lichtgeschwindigkeit**2)**2-(Neutronenmasse*Lichtgeschwindigkeit**2)**2)


!       write(6,*)'wlr =', NeutronWavelength_fromE
    endif
    
end function NeutronWavelength_fromE




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




double precision function NeutronEnergy_fromLambda(Lambda)
!---------------------------------------------------------
    implicit none

    double precision, intent(in) :: Lambda
    double precision             :: E0, hc

     ! in non-relativistic approximation
   
     if(Lambda >  Planckkonstante /sqrt(2*Elimit_classical*Neutronenmasse)) then

       NeutronEnergy_fromLambda =  Planckkonstante**2 /(2d0*Lambda**2*Neutronenmasse)

!       write(6,*)'Ecl = ',NeutronEnergy_fromLambda

    else


    ! relativistic expression

     hc = Planckkonstante*Lichtgeschwindigkeit
     E0 = Neutronenmasse * Lichtgeschwindigkeit**2     

     NeutronEnergy_fromLambda = -E0+sqrt(E0**2+(hc/Lambda)**2)  

!     write(6,*)'Erl = ',NeutronEnergy_fromLambda

    endif    


end function NeutronEnergy_fromLambda



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



double precision function NeutronVelocity_fromLambda(Lambda)
!-----------------------------------------------------------
    implicit none

    double precision, intent(in) :: Lambda
    double precision             :: E0, hc, E

! in non-relativistic approximation
!   
     if(Lambda >  Planckkonstante /sqrt(2*Elimit_classical*Neutronenmasse)) then

       E =  Planckkonstante**2 /(2d0*Lambda**2*Neutronenmasse)
       NeutronVelocity_fromLambda = sqrt(2*E/Neutronenmasse)

!       write(6,*)'vcl = ',NeutronVelocity_fromLambda

    else
!

! relativistic expression
!------------------------
     hc = Planckkonstante*Lichtgeschwindigkeit
     E0 = Neutronenmasse * Lichtgeschwindigkeit**2     
     E  = -E0+sqrt(E0**2+(hc/Lambda)**2)
     NeutronVelocity_fromLambda = Lichtgeschwindigkeit*sqrt(1d0-(E0/(E+E0))**2)  
!
!     
!     write(6,*)'vrl = ',NeutronVelocity_fromLambda

    endif    


end function NeutronVelocity_fromLambda






end Module PhysConstants
