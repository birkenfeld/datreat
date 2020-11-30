      PROGRAM unift

C     Universal Fourier Transformer
C     Version 3.2
C
C     This program Fourier transforms TOF and BS data with INX format.
C     Constant-Q-ing should be performed prior to this program by INTERPOL
C     or SIGRID. Special features:
C
C       1. non-equidistant E scale (no FFT used)
C       2. "mirroring" of negative Es
C       3. optional energy range restriction
C       4. symmetrization of S(Q,E) (classical)
C       5. time step adjusted to E range individually
C       6. error propagation (v. 2.x)
C       7. maximum number of t-values -> 3000
C       8. optimal time step can be overridden
C       9. linear interpolation of S(Q,E) data (v. 3.x)
C      10. averaging of spectra (only for identical E scales)
C
C
C     Item 2 may be realized more correctly.

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NE=2048,NQ=128)
C     Commons for CRIN use:
      REAL*4 SPECPAR,ANGLES,NUMSPE
      REAL*4 X,Y,YE
      COMMON/SAMPL1/SPECPAR(30),ANGLES(NQ),NUMSPE(NQ)
      COMMON/SAMPL2/X(NE,NQ),Y(NE,NQ),YE(NE,NQ)
      CHARACTER*80 FILIN
      CHARACTER*40 TITLIN
      PARAMETER (mt=3000)
C     Arrays and variables for UNIFT use:
      DIMENSION es(ne,nq),ss(ne,nq),ds(ne,nq)
      DIMENSION ebegs(nq),jbegs(nq),eends(nq),jends(nq)
      DIMENSION es1(ne,nq),ss1(ne,nq),ds1(ne,nq),ns1(nq)
      DIMENSION er(ne,nq),sr(ne,nq),dr(ne,nq)
      DIMENSION ebegr(nq),jbegr(nq),eendr(nq),jendr(nq)
      DIMENSION er1(ne,nq),sr1(ne,nq),dr1(ne,nq),nr1(nq)
      DIMENSION ndats(nq),ndatr(nq)
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

C     Read the sample file:

 1    WRITE(6,1001)
 1001 FORMAT('Sample file name [zero treatment]? ',$)
      READ(5,'(a80)') inp1
      filin=parse(inp1,1)
      itrz=int(0.5d0+val(parse(inp1,2)))
      CALL crin(filin,titlin,ifail)
      IF (ifail.EQ.1) GOTO 1
      nqs=int(specpar(4))
      nes=int(specpar(19))

C     Skip 'zero channels' and write sample file statistics:

      WRITE(6,*) ' #  angle/Q   first channel   last channel'
      DO i=1,nqs
        DO j=1,nes
          IF ((x(j,i).eq.0.0).and.(y(j,i).eq.0.0)) GOTO 13
          if ((y(j,i).eq.0.0).and.(ye(j,i).eq.0.0)) then
            if (itrz.eq.1) goto 13
            if (itrz.eq.0) then 
              write(6,*) 'Read from file:'
              write(6,1111) x(j,i),y(j,i),ye(j,i)
 1111         format(5X,F10.5,E13.5,E12.4)
              write(6,1112)
 1112         format('Treat this as (1) missing value',
     &               ' or (2) zero value? ',$)
              read(5,'(i1)') itrz
              if (itrz.eq.1) goto 13
            end if
          end if
          goto 11
 13     END DO
 11     ebegs(i)=x(j,i)
        jbegs(i)=j
        DO j=jbegs(i),nes
          IF ((x(j,i).EQ.0.0).AND.(y(j,i).EQ.0.0)) GOTO 12
          if ((y(j,i).eq.0.0).and.(ye(j,i).eq.0.0)) then
            if (itrz.eq.1) goto 12
            if (itrz.eq.0) then 
              write(6,*) 'Read from file:'
              write(6,1111) x(j,i),y(j,i),ye(j,i)
              write(6,1112)
              read(5,'(i1)') itrz
              if (itrz.eq.1) goto 12
            end if
          end if
        END DO
 12     eends(i)=x(j-1,i)
        jends(i)=j-1
        ndats(i)=jends(i)-jbegs(i)+1
        WRITE(6,1003) i,angles(i),jbegs(i),ebegs(i),jends(i),eends(i)
 1003   FORMAT(i3,f9.3,i7,f9.3,i6,f9.3)
      END DO

C     Sort to ascending order:

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

C     Read the resolution file:

      WRITE(6,*)
 2    WRITE(6,1002)
 1002 FORMAT('Resolution file name? ',$)
      READ(5,'(a80)') filin
      CALL crin(filin,titlin,ifail)
      IF (ifail.ne.0) GOTO 2
      nqr=int(specpar(4))
      ner=int(specpar(19))

C     Skip 'zero channels' and write resolution file statistics:

      WRITE(6,*) ' #  angle/Q   first channel   last channel'
      DO i=1,nqr
        DO j=1,ner
          IF ((x(j,i).eq.0.0).and.(y(j,i).eq.0.0)) GOTO 19
          if ((y(j,i).eq.0.0).and.(ye(j,i).eq.0.0)) then
            if (itrz.eq.1) goto 19
            if (itrz.eq.0) then 
              write(6,*) 'Read from file:'
              write(6,1111) x(j,i),y(j,i),ye(j,i)
              write(6,1112)
              read(5,'(i1)') itrz
              if (itrz.eq.1) goto 19
            end if
          end if
          goto 18
 19     END DO
 18     ebegr(i)=x(j,i)
        jbegr(i)=j
        DO j=jbegr(i),nes
          IF ((x(j,i).EQ.0.0).AND.(y(j,i).EQ.0.0)) GOTO 14
          if ((y(j,i).eq.0.0).and.(ye(j,i).eq.0.0)) then
            if (itrz.eq.1) goto 14
            if (itrz.eq.0) then 
              write(6,*) 'Read from file:'
              write(6,1111) x(j,i),y(j,i),ye(j,i)
              write(6,1112)
              read(5,'(i1)') itrz
              if (itrz.eq.1) goto 14
            end if
          end if
        END DO
 14     eendr(i)=x(j-1,i)
        jendr(i)=j-1
        ndatr(i)=jendr(i)-jbegr(i)+1
        WRITE(6,1003) i,angles(i),jbegr(i),ebegr(i),jendr(i),eendr(i)
      END DO

C     Sort to ascending order:

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

C     Check sample/resolution consistency:

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
 20     WRITE(6,*) 'Warning: Different energy ranges for some Qs,'
        WRITE(6,*) '         restrict energy range manually!'
 21   END IF

C     Start treating a Q value:

      WRITE(6,*)
 30   WRITE(6,1010)
 1010 FORMAT('Spectra: sample (from-to)',
     &       ' resolution (from-to)? (0: end) ',$)
      read(5,'(a80)') inp2
      iqs=int(0.5+val(parse(inp2,1)))
      iqr=int(0.5+val(parse(inp2,2)))
      iqse=int(0.5+val(parse(inp2,3)))
      iqre=int(0.5+val(parse(inp2,4)))
      IF (iqs.EQ.0) GOTO 999
      
C     Q-average data if desired:

      if ((iqse.ne.0).and.(iqse.ne.iqs)) then
        nave=iqse-iqs+1
        do j=1,ndats(iqs)
          do iq=iqs+1,iqse
            es(j,iqs)=es(j,iqs)+es(j,iq)
            ss(j,iqs)=ss(j,iqs)+ss(j,iq)
            ds(j,iqs)=sqrt(ds(j,iqs)**2+ds(j,iq)**2)
          end do
          es(j,iqs)=es(j,iqs)/nave
          ss(j,iqs)=ss(j,iqs)/nave
          ds(j,iqs)=ds(j,iqs)/nave
        end do
      end if
      if ((iqre.ne.0).and.(iqre.ne.iqr)) then
        nave=iqre-iqr+1
        do j=1,ndatr(iqr)
          do iq=iqr+1,iqre
            er(j,iqr)=er(j,iqr)+er(j,iq)
            sr(j,iqr)=sr(j,iqr)+sr(j,iq)
            dr(j,iqr)=sqrt(dr(j,iqr)**2+dr(j,iq)**2)
          end do
          er(j,iqr)=er(j,iqr)/nave
          sr(j,iqr)=sr(j,iqr)/nave
          dr(j,iqr)=dr(j,iqr)/nave
        end do
      end if

C     Ask energy limit for Fourier transform:

 50   WRITE(6,1011)
 1011 FORMAT('Limit for E range? (>0, =0: no limit) ',$)
      READ(5,*) emir
      IF (emir.eq.0) emir=1.0e30

C     Calculate range of sample data effectively used:

      DO j=1,ndats(iqs)
        IF (abs(es(j,iqs)).le.emir) GOTO 45
      END DO
 45   jbegse=j
      DO j=jbegse,ndats(iqs)
        IF (abs(es(j,iqs)).gt.emir) GOTO 46
      END DO
 46   jendse=j-1
      write(6,1021) jbegse,jendse,es(jbegse,iqs),es(jendse,iqs)
 1021 FORMAT('Effective limits for sample:    ',
     1       i4,'...',i4,f10.3,'...',f8.3)

C     Calculate range of resolution data effectively used:

      DO j=1,ndatr(iqr)
        IF (abs(er(j,iqr)).le.emir) GOTO 40
      END DO
 40   jbegre=j
      DO j=jbegre,ndatr(iqr)
        IF (abs(er(j,iqr)).gt.emir) GOTO 41
      END DO
 41   jendre=j-1
      write(6,1020) jbegre,jendre,er(jbegre,iqr),er(jendre,iqr)
 1020 FORMAT('Effective limits for resolution:',
     1       i4,'...',i4,f10.3,'...',f8.3)

      emaxr=max(abs(er(jbegre,iqr)),abs(er(jendre,iqr)))
      emaxs=max(abs(es(jbegse,iqs)),abs(es(jendse,iqs)))
      dt=2.067834/((emaxs+emaxr)*0.5)
      WRITE(6,1014) dt
 1014 FORMAT('Timestep calculated from energy ranges:',f10.3,' ps')

      WRITE(6,1017)
 1017 FORMAT('Accept these values? (y/n) ',$)
      READ(5,'(a)') ans
      IF (.NOT.((ans.EQ.'n').OR.(ans.EQ.'N'))) GOTO 55
      WRITE(6,1018)
 1018 FORMAT('Forced timestep? (0: recalculate with new E_max) ',$)
      READ(5,*) dtf
      if (dtf.eq.0.0) goto 50
      dt=dtf

C     Ask limit for mirroring:
C     (The convention is: emir positive, but limit is -emir.)

 55   WRITE(6,1012)
 1012 FORMAT('Limit for E mirroring? (0: no mirroring) ',$)
      READ(5,*) emir
      emir=abs(emir)

C     Ask for symmetrization temperature:

      WRITE(6,1013)
 1013 FORMAT('Temperature for symmetrization? (0: no symm.) ',$)
      READ(5,*) temp
      IF (temp.EQ.0.0) THEN
        fdbf=0.0
      ELSE
        fdbf=5.80222/temp
      END IF

C     Mirroring, Symmetrisation, and conversion E->omega:

      call mirror(iqs,es,ss,ds,jbegse,jendse,fdbf,emir,
     &            es1,ss1,ds1,ns1)
      call mirror(iqr,er,sr,dr,jbegre,jendre,0.0d0,emir,
     &            er1,sr1,dr1,nr1)
      do j=1,ns1(iqs)
        write(17,'(f9.5,2e13.5)') es1(j,iqs),ss1(j,iqs),ds1(j,iqs)
      end do
      do j=1,nr1(iqr)
        write(18,'(f9.5,2e13.5)') er1(j,iqr),sr1(j,iqr),dr1(j,iqr)
      end do

C     Ask for time limit:

      WRITE(6,1029)
 1029 FORMAT('Maximum time to calculate (ps)?  ',$)
      READ(5,*) tmax
      nt=int(tmax/dt)

C     Ask for output file name:

      WRITE(6,1030)
 1030 FORMAT('Output file name? (w/o extension) ',$)
      READ(5,'(a20)') ofnam
      DO ich=1,20
        IF (ofnam(ich:ich).EQ.' ') GOTO 60
      END DO
 60   nchofn=ich-1

C     Outer loop over times (sample):

      WRITE(6,*)
      WRITE(6,*) 'Transforming sample data:'
      WRITE(6,*)
     1 '       t [ps]      real      imag      ',
     2 ' abs       arg      norm      rerr'
      DO it=1,nt+1
        t=(it-1)*dt

C     Fourier transform the sample data:

        call fourier(es1,ss1,ds1,iqs,ns1(iqs),t,
     1               ft1,ft2,fts(it),dfts(it))
     
        WRITE(6,1040) it,t,ft1,ft2,fts(it),ft2/ft1,
     1                fts(it)/fts(1),dfts(it)/fts(it)
 1040   FORMAT(i4,7f10.3)

      END DO

C     Write Fourier transformed SAMPLE DATA to file:

      OPEN(unit=41,file=ofnam(1:nchofn)//'.sam',status='unknown')
      DO it=1,nt+1
        t=(it-1)*dt
        WRITE(41,1045) t,fts(it),dfts(it)
 1045   FORMAT(3f12.4)
      END DO
      CLOSE(41)

C     Outer loop over times (resolution):

      WRITE(6,*)
      WRITE(6,*) 'Transforming resolution data:'
      WRITE(6,*)
     1 '       t [ps]      real      imag      ',
     2 ' abs       arg      norm      rerr'
      DO it=1,nt+1
        t=(it-1)*dt

C     Fourier transform the resolution data:

        call fourier(er1,sr1,dr1,iqr,nr1(iqr),t,
     1               ft1,ft2,ftr(it),dftr(it))

        WRITE(6,1040) it,t,ft1,ft2,ftr(it),ft2/ft1,
     1                ftr(it)/ftr(1),dftr(it)/ftr(it)

      END DO

C     Write Fourier transformed RESOLUTION DATA to file:

      OPEN(unit=42,file=ofnam(1:nchofn)//'.res',status='unknown')
      DO it=1,nt+1
        t=(it-1)*dt
        WRITE(42,1045) t,ftr(it),dftr(it)
      END DO
      CLOSE(42)

C     Write DECONVOLUTED S(Q,t) DATA to file:

      OPEN(unit=43,file=ofnam(1:nchofn)//'.sqt',status='unknown')
      DO it=1,nt+1
        t=(it-1)*dt
        sqt=fts(it)/ftr(it)
        sqterr=sqt*sqrt((dftr(it)/ftr(it))**2+(dfts(it)/fts(it))**2)
        WRITE(43,1045) t,sqt,sqterr
      END DO
      CLOSE(43)

C     Play it again, Sam...

      WRITE(6,*)
      GOTO 30

C     This is the end, my only friend, the end...

 999  END


C     ILL routine CRIN:

          SUBROUTINE CRIN(FILIN,TITLIN,IFAIL)
 
C......... Data input of CROSSX files
 
           PARAMETER (NE=2048,NQ=128)
           COMMON/SAMPL1/SPECPAR(30),ANGLES(NQ),NUMSPE(NQ)
           COMMON/SAMPL2/X(NE,NQ),Y(NE,NQ),YE(NE,NQ)
           CHARACTER*80 FILIN
           CHARACTER*40 TITLIN
           DIMENSION NZONE(6),NCS(NQ)
 
4     FORMAT(8I5)
5     FORMAT(A40)
6     FORMAT(1X,F6.2,F8.3,F8.4,F9.3,F6.1,I2)
7     FORMAT(16X,3F8.4)
8     FORMAT(1X,A1)
9     FORMAT(5X,F10.5,E13.5,E12.4)


       IR=5
       IW=6
       IFAIL=0
       IFILE2=21
       OPEN (UNIT=IFILE2,FILE=FILIN,STATUS='OLD',ERR=95)

C      Test input replacement:
C      nspect=1
C      j=1
C220   read(ifile2,*,end=210) x(j,1),y(j,1),ye(j,1)
C      j=j+1
C      goto 220
C210   ncmin=j
C      einc=25
C      angles(1)=90
C      numspe(1)=1
C      goto 230
 
       NSPECT=0
       NCMIN=NE
       DO I=1,NQ
        READ (IFILE2,4,END=90,err=90) NLINES,NZONE,NCHANS
        NSPECT=NSPECT+1
        NUMSPE(I)=I!This is normally to identify the detector-Cannot know where
c       the file comes from => just take the index as detector number but
c       beware of programs testing this number
        NCS(I)=NCHANS!Different spectra can have different number of channels
        IF (NCHANS.LT.NCMIN) NCMIN=NCHANS !..but we don't like that so take min
        DO J=1,6
         IF (NZONE(J).LE.0) GOTO 100
         IF (J.EQ.1) THEN
         READ (IFILE2,5) TITLIN
         ELSE
         IF (J.EQ.2) THEN
         READ (IFILE2,6) ANGLES(I),EINC,QINC,TEMP,AMASS,ISYM
         READ (IFILE2,7)DELTAEN,DELTATAU,DELTAK
         ELSE
         IF (J.EQ.3) THEN
         READ (IFILE2,8) DATYPE
         END IF
         END IF
         END IF
100     END DO
 
       DO J=1,NCHANS
        READ(IFILE2,9) X(J,I),Y(J,I),YE(J,I)
       END DO
       END DO
90     CONTINUE
      PI=3.141592653
230   SPECPAR(4)=NSPECT
      SPECPAR(19)=NCMIN  !N channels: we take the min to be sure it's filled
c                          with data
      SPECPAR(9)=TEMP
      SPECPAR(18)=DELTATAU
      IF (EINC.EQ.0.) THEN
         WRITE (IW,33)
33       FORMAT(1X,' No incident energy data found:Give it (meV)->',$)
         READ (IR,*) EINC
      END IF
      SPECPAR(21)=SQRT(81.799/EINC)
 
       CLOSE (UNIT=IFILE2)
       RETURN
95     IFAIL=1
       RETURN
       END


      SUBROUTINE mirror(iq,es,ss,ds,jbeg,jend,fdbf,emir,
     1                  es1,ss1,ds1,ns1)

C     This subroutine mirrors the high energy part of the spectrum
C     to complement the energy loss side. At the same time it transforms
C     from S(Q,E) to S(Q,omega).
     
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NE=2048,NQ=128)
      DIMENSION es(ne,nq),ss(ne,nq),ds(ne,nq)
      DIMENSION es1(ne,nq),ss1(ne,nq),ds1(ne,nq),ns1(nq)

      if (es(jbeg,iq).gt.-emir) then
        if (es(jbeg,iq).gt.0.0) then
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
        ee=es(je,iq)
        if (ee.le.emir) goto 100
        es1(j1,iq)=-1.519257*ee
        dbf=exp(fdbf*ee)
        ss1(j1,iq)=dbf*ss(je,iq)*0.658216
        ds1(j1,iq)=dbf*ds(je,iq)*0.658216
        j1=j1+1
      end do
      stop
        
 100  do je=jbeg,jend
        if (es(je,iq).ge.-emir) goto 110
      end do
      stop

 110  jlim=je
 120  do je=jlim,jend
        ee=es(je,iq)
        es1(j1,iq)=1.519257*ee
        dbf=exp(fdbf*ee)
        ss1(j1,iq)=dbf*ss(je,iq)*0.658216
        ds1(j1,iq)=dbf*ds(je,iq)*0.658216
        j1=j1+1
      end do
        
      ns1(iq)=j1-1
      
 900  return
      end

      
      SUBROUTINE fourier(w,s,ds,iq,n,t,ft1,ft2,ft,dft)
     
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NE=2048,NQ=128)
      DIMENSION w(ne,nq),s(ne,nq),ds(ne,nq)
      DIMENSION a1(ne),a2(ne),swt(ne),cwt(ne)
      
      if (t.ne.0.0) then
      
        do i=1,n
          swt(i)=sin(w(i,iq)*t)
          cwt(i)=cos(w(i,iq)*t)
        end do
        t2=1.0/(t*t)
        
        a1(1)=-swt(1)/t+(cwt(1)-cwt(2))/(w(2,iq)-w(1,iq))*t2
        a2(1)=cwt(1)/t+(swt(1)-swt(2))/(w(2,iq)-w(1,iq))*t2
        do i=2,n-1
          a1(i)=((cwt(i)-cwt(i+1))/(w(i+1,iq)-w(i,iq))
     1          +(cwt(i-1)-cwt(i))/(w(i-1,iq)-w(i,iq)))*t2
          a2(i)=((swt(i)-swt(i+1))/(w(i+1,iq)-w(i,iq))
     1          +(swt(i-1)-swt(i))/(w(i-1,iq)-w(i,iq)))*t2
        end do
        a1(n)=swt(n)/t+(cwt(n-1)-cwt(n))/(w(n-1,iq)-w(n,iq))*t2
        a2(n)=-cwt(n)/t+(swt(n-1)-swt(n))/(w(n-1,iq)-w(n,iq))*t2
        
      else
      
        a1(1)=(w(2,iq)-w(1,iq))*0.5
        a2(1)=0.0
        do i=2,n
          a1(i)=(w(i+1,iq)-w(i-1,iq))*0.5
          a2(i)=0.0
        end do
        a1(n)=(w(n,iq)-w(n-1,iq))*0.5
        a2(n)=0.0
        
      end if
      
      ft1=0.0
      ft2=0.0
      do i=1,n
        ft1=ft1+a1(i)*s(i,iq)
        ft2=ft2+a2(i)*s(i,iq)
      end do
      ft=sqrt(ft1*ft1+ft2*ft2)
      
      dft2=0.0
      do i=1,n
        dft2=dft2+((a1(i)*ft1+a2(i)*ft2)/ft*ds(i,iq))**2
      end do
      dft=sqrt(dft2)
      
      return
      end
          

      character*80 function parse(line,nparq)

C     returns the nparq-th component of line (separated by blank space)

      character*80 line
      character*80 param
      character och,TAB
C     parameter(TAB=char(9))

      TAB=char(9)
      ich0=1
      npar=1

C     skip blanks and tabs

 10   do ich=ich0,80
        if ((line(ich:ich).ne.' ').and.(line(ich:ich).ne.TAB)) goto 20
      end do
      param=' '
      goto 99
 20   ich0=ich

C     read parameter

      jch=1
      param=' '
      do ich=ich0,80
        if ((line(ich:ich).ne.' ').and.(line(ich:ich).ne.TAB)) then
          param(jch:jch)=line(ich:ich)
          jch=jch+1
        else
          goto 30
        end if
      end do
      goto 99
 30   if (npar.eq.nparq) goto 99
      npar=npar+1
      ich0=ich
      goto 10

 99   parse=param

      return
      end


      real*8 function val(str)

      character*80 str

C     write(*,*) '!',str,'!'
      if (str.eq.'') goto 10
      read(str,*,err=10) val
      goto 20

 10   val=0.0

 20   return
      end


