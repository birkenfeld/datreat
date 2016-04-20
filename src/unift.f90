   Module unift

                                                                        
!     Universal Fourier Transformer                                     
!     Version 3.2                                                       
!                                                                       
!     This program Fourier transforms TOF and BS data with INX format.  
!     Constant-Q-ing should be performed prior to this program by INTERP
!     or SIGRID. Special features:                                      
!         
! ---> based on routines of Reiner Zorn unift_3_2 <----
!
!
   use PhysicalConstantsPlus
   use cdata
   use selist
   use constants

   contains

      SUBROUTINE uni_ft(ipoint_da, ipoint_re, ipoint_fts) 
!     ---------------------------------------------------

      implicit none

      integer                :: ipoint_da             !! Addresse der "Daten" 
      integer                :: ipoint_re             !! Addresse der "Referenz" 
      integer                :: ipoint_fts            !! StartAddresse der Ergebnisse

      
      double precision       :: es(mwert) , ss(mwert) , ds(mwert)    !! Daten Energie, Werte, Fehler
      double precision       :: es1(mwert), ss1(mwert), ds1(mwert)
      double precision       :: er(mwert) , sr(mwert) , dr(mwert)    !! Referenz
      double precision       :: er1(mwert), sr1(mwert), dr1(mwert)
      integer                :: ndats, ndats1                        !! number of valid data-points
      integer                :: ndatr, ndatr1                        !! number of valid ref-points
      integer                :: ns1, nr1                             !! number of points of mirrored vec's

      double precision       :: tv(mwert)                            !! Timevektor for FT's
      double precision       :: fts(mwert), dfts(mwert)              !! FT-Daten      roh
      double precision       :: ftr(mwert), dftr(mwert)              !! FT-Aufloesung roh
      double precision       :: sqt(mwert), sqterr(mwert)            !! FT-Aufloesung roh
  
      double precision       :: t, ft1, ft2
      double precision       :: dt                                   !! timestep
      double precision       :: emaxs, emaxr, emir                   !! "energy" range specifiers
      double precision       :: fdbf                                 !! detailed balance exp-argument

      double precision       :: domega, tmax

      double precision       :: eunit=Elektronenladung*1d-6          !! energy scale unit in J
      double precision       :: tunit=1d-9                           !! desired time unit for result in s
      character*8            :: eunit_name='micro-eV'
      character*8            :: tunit_name='ns      '
       
      real                   :: temp

      integer                :: i, it, ier
      integer                :: ipoint_out
      integer                :: iout
      integer                :: ida1 = 2
      integer                :: ielastic, nft
!     ------------------------------------------------------------

!---  get the data -->

      if(ipoint_da.gt.0) then

        call DataGet(ipoint_da,es,ss,ds,ndats)
        write(6,*)'uni_ft: got ',ndats,' data points from record:',ipoint_da,' :  ', name(ipoint_da)

        emaxs = es(ndats)
        emir  = emaxs
        write(6,*)'E-Range = ',emaxs

        if     (index(xname(ipoint_da),'eV')  .gt.0 .and. index(xname(ipoint_da),'mic') .gt. 0  ) then
          eunit=Elektronenladung*1d-6 
          eunit_name='micro-eV'
        else if (index(xname(ipoint_da),'eV')  .gt.0 .and. index(xname(ipoint_da),'milli') .gt. 0) then
          eunit=Elektronenladung*1d-3
          eunit_name='meV     '
        else if (index(xname(ipoint_da),'meV') .gt.0                                             ) then
          eunit=Elektronenladung*1d-3
        else if (index(xname(ipoint_da),'GHz') .gt.0                                             ) then
          !! assuming omega not nue 
          eunit=1d9*Planckkonstante/2/Pi
          eunit_name='GHz     '
        else if (index(xname(ipoint_da),'omega') .gt.0                                           ) then
          !! assuming omega not nue is in GigaRad/s !
          eunit=1d9*Planckkonstante/2/Pi
          eunit_name='GHz     '
        else
          Write(6,*)'uni_ft: eunit not reckognizable from xaxis, assuming defaults: '
          write(6,*)'      : ',eunit_name,' ==> ',eunit,' J'
        endif
      
        temp = 300.0  
        call parget('temp    ',temp ,ipoint_da, ier) 
       
      else
        write(6,*)'uni_ft: no data selected ! NO ACTION '
        return
      endif

!---  get the reference -->
      if(ipoint_re.gt.0) then
        call DataGet(ipoint_re,er,sr,dr,ndatr)
        write(6,*)'uni_ft: got ',ndatr,' ref  points from record:',ipoint_re,' :  ', name(ipoint_re) 
        emaxs = er(ndatr)
        write(6,*)'E-Range = ',emaxr
        if(xname(ipoint_da).ne.xname(ipoint_re)) then
         write(6,*)'uni_ft: e-units of data and refernce are different! No action'
        endif
      else
        write(6,*)'uni_ft: no refernce selected: just transforming data only... '
      endif

!---- optimum time step -->
  
      dt = Planckkonstante/eunit/2/tunit / emaxs       !!! [ h/eunit/2 / tunit] /  ((emaxs+emaxr)*0.5) !! here eunit=meV, tunit=ps        write(6,*)'uni_ft: timestep = ',dt,tunit_name

      fdbf  = eunit/2/Boltzmannkonstante/temp
! Marshall Lovesey Convention omega = energy loss of the neutron
      fdbf  = -fdbf


      write(6,'(a)        ') "#####################  uni_ft ###############################################"
      write(6,'(a,e14.7,a)') "#  will use detailed balance symmetrisation: exp(",fdbf," * omega) "
      write(6,'(a,f7.2,a) ') "#  assumed temperature is T =",temp," K"
      write(6,'(a,e14.7,a)') "#  a negative coeff corresponds to Marshall-Lovesey definition: "
      write(6,'(a,e14.7,a)') "#  positive omega is energy LOSS of the NEUTRON "
      write(6,'(a)        ') "###### >>>>>>>>>>>> "

!     conversion E->omega:               
                                                                        
      ns1 = mwert     ! will be updated by uni_mirror
      nr1 = mwert     !  "        "           "

      
      ipoint_out = ipoint_fts
      call uni_mirror(es,ss,ds,ida1,ndats,fdbf ,emir, es1,ss1,ds1,ns1, eunit, tunit)
      call DataCopy(ipoint_da,ipoint_out)             !! copy data here aonly to transfer auxiliary params...
      call DataPut(ipoint_out,es1,ss1,ds1,ns1)
      numor(ipoint_out)     = numor(ipoint_out) + 100000
      name(ipoint_out)(8:8) = 'M'
      write(coment(ipoint_out),'(a,a)') 'uni_ft M:',coment(ipoint_out)(1:71)    !! check ob das so geht, besser hilfskft.
      xname(ipoint_out) ='1/'//tunit_name(1:6) 
      call parset('tunit   ', sngl(tunit), ipoint_out)   


      if(ipoint_re.gt.0) then
       ipoint_out=ipoint_out+1
       call uni_mirror(er,sr,dr,ida1,ndatr,0.0d0,emir, er1,sr1,dr1,nr1, eunit, tunit)
       call DataCopy(ipoint_re,ipoint_out)             !! copy data here only to transfer auxiliary params..
       call DataPut(ipoint_out,er1,sr1,dr1,nr1)      
       numor(ipoint_out)     = numor(ipoint_out) + 100000
       name(ipoint_out)(8:8) = 'M'
       write(coment(ipoint_out),'(a,a)') 'uni_ft M:',coment(ipoint_out)(1:71)    !! check ob das so geht, besser hilfskft.
       xname(ipoint_out) ='1/'//tunit_name(1:6) 
       call parset('tunit   ', sngl(tunit), ipoint_out)   
      endif

      ipoint_out=ipoint_out+1

!!! Fouriertransform of data !!

!!>> neu:
!! determine maximum useful time of ft
       ielastic = max( minloc(abs(es1(1:ns1)),1) ,2 )
       domega   = abs(es1(ielastic)-es1(ielastic-1))
       tmax     = 2*Pi/domega    
       nft      = min( nint(tmax/dt+1)  , mwert) 
       write(6,'(a,i8)'   )"elastic channel nr. = ",ielastic
       write(6,'(a,e14.7)')"d_omega             = ",domega
       write(6,'(a,e14.7)')"tmax                = ",tmax
       write(6,'(a,e14.7)')"dt                  = ",dt
       write(6,'(a,i8)'   )"nft                 = ",nft
!!<<
      DO it=1,nft 
        t=(it-1)*dt 
        tv(it) = t
        call uni_fourier(es1,ss1,ds1,ns1,t,ft1,ft2,fts(it),dfts(it))
        
        if(iout().gt.0) then          
          if(it.eq.1) then
           write(6,*) ' FT of data: '
           WRITE(6,*) '       t           real      imag       abs       arg      norm      rerr'
          endif                                                  
          WRITE(6,'(i4,7f10.3)') it,t,ft1,ft2,fts(it),ft2/ft1,fts(it)/fts(1),dfts(it)/fts(it)
        endif                                                                
      END DO 

      call DataCopy(ipoint_da,ipoint_out)             !! copy data here only to transfer auxiliary params...
      call DataPut(ipoint_out,tv,fts,dfts,nft)
      numor(ipoint_out)     = numor(ipoint_out) + 200000
      name(ipoint_out)(8:8) = 'F'
      write(coment(ipoint_out),'(a,a)') 'uni_ft F:',coment(ipoint_out)(1:71)    !! check ob das so geht, besser hilfskft.
      xname(ipoint_out) = tunit_name 
      call parset('tunit   ', sngl(tunit), ipoint_out)   

      

!!! Fouriertransform of res !!

       if(ipoint_re.ne.0) then
        ipoint_out = ipoint_out+1
        DO it=1,nft 
          t=(it-1)*dt 
          tv(it) = t
          call uni_fourier(er1,sr1,dr1,nr1,t,ft1,ft2,ftr(it),dftr(it))
          
          if(iout().gt.0) then          
            if(it.eq.1) then
             write(6,*) ' FT of resolution: '
             WRITE(6,*) '       t           real      imag       abs       arg      norm      rerr'
            endif                                                  
            WRITE(6,'(i4,7f10.3)') it,t,ft1,ft2,fts(it),ft2/ft1,fts(it)/fts(1),dfts(it)/fts(it)
          endif                                                                
        END DO 
  
        call DataCopy(ipoint_re,ipoint_out)             !! copy data here only to transfer auxiliary params...
        call DataPut(ipoint_out,tv,ftr,dftr,nft)
        numor(ipoint_out)     = numor(ipoint_out) + 200000
        name(ipoint_out)(8:8) = 'F'
        write(coment(ipoint_out),'(a,a)') 'uni_ft F:',coment(ipoint_out)(1:71)    !! check ob das so geht, besser hilfskft.
        xname(ipoint_out) = tunit_name 
        call parset('tunit   ', sngl(tunit), ipoint_out)   

        ipoint_out = ipoint_out+1
     


!!!
!     Write DECONVOLUTED S(Q,t) DATA to file:                           
                                                                        
        DO it=1,nft
          t=(it-1)*dt
          tv(it) = t 
          sqt(it)=fts(it)/ftr(it) 
          sqterr(it)=sqt(it)*sqrt((dftr(it)/ftr(it))**2+(dfts(it)/fts(it))**2) 
                   if(iout().gt.0) then          
              if(it.eq.1) then
               write(6,*) ' deconvouted data: '
               WRITE(6,*) '       t           S(Q,t)              sqterr '
              endif                                                  
              write(6,'(f12.6,e16.8,e12.4)') tv(it),sqt(it),sqterr(it) 
            endif                                                                
        END DO 

        call DataCopy(ipoint_da,ipoint_out)             !! copy data here only to transfer auxiliary params...
        call DataPut(ipoint_out,tv,sqt,sqterr,nft)
        numor(ipoint_out)     = numor(ipoint_out) + 300000
        name(ipoint_out)(8:8) = 'S'
        write(coment(ipoint_out),'(a,a)') 'uni_ft S:',coment(ipoint_out)(1:71)    !! check ob das so geht, besser hilfskft.
        xname(ipoint_out) = tunit_name 
        call parset('tunit   ', sngl(tunit), ipoint_out)   

      endif

      write(6,*)'uni_ft: finished: --> dir  to identify the final and intermediate results '                                                                       

      end subroutine uni_ft

                                                            
                                                                        
      SUBROUTINE uni_mirror(es,ss,ds,jbeg,jend,fdbf,emir,es1,ss1,ds1,ns1,eunit,tunit)
                                                                        
!     This subroutine mirrors the high energy part of the spectrum      
!     to complement the energy loss side. At the same time it transforms
!     from S(Q,E) to S(Q,omega).                                        
                                                                        
      implicit none
                                                                 
      integer, intent(in)           :: jbeg, jend
      integer, intent(inout)        :: ns1      
      double precision, intent(in)  :: es(jend),ss(jend),ds(jend) 
      double precision, intent(out) :: es1(ns1),ss1(ns1),ds1(ns1) 
      double precision, intent(in)  :: fdbf, emir, eunit, tunit

      integer                       :: j1, je, jlim
      double precision              :: ee, dbf
      double precision              :: ef_conversion


      ef_conversion = 2*Pi*eunit/Planckkonstante*tunit


                                                                        
      if (es(jbeg).gt.-emir) then 
        if (es(jbeg).gt.0.0) then 
          write(6,'(a,i3)') 'Fatal: No elastic line ' 
          ns1=0 
          return 
        else 
          write(6,'(a,i3)') 'Warning: Cut elastic line ' 
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
        es1(j1)=-ef_conversion*ee             !!! >>> !!! -[2*pi*eunit/h * tunit] *ee 
        dbf=exp(fdbf*ee)                      !!! detailed balance factor
        ss1(j1)=dbf*ss(je)/ef_conversion      !!! >>> !!! dbf*ss(je)*[1/(2*pi*eunit/h * tunit)]
        ds1(j1)=dbf*ds(je)/ef_conversion
        j1=j1+1
        if(j1.gt.ns1) then
          write(6,*)'uni_mirror: dimension of output vector too low (1) !'
          return
        endif 
      end do 
      stop 
                                                                        
  100 do je=jbeg,jend 
        if (es(je).ge.-emir) goto 110 
      end do 
      stop 
                                                                        
  110 jlim=je 
  120 do je=jlim,jend 
        ee=es(je) 
        es1(j1)=ef_conversion*ee 
        dbf=exp(fdbf*ee) 
        ss1(j1)=dbf*ss(je)/ef_conversion
        ds1(j1)=dbf*ds(je)/ef_conversion
        j1=j1+1 
        if(j1.gt.ns1) then
          write(6,*)'uni_mirror: dimension of output vector too low (2) !'
          return
        endif 
      end do 
                                                                        
      ns1=j1-1 
                                                                        
  900 return 
      END  SUBROUTINE uni_mirror                                         
                                                                        




                                                                        
      SUBROUTINE uni_fourier(w,s,ds,n,t,ft1,ft2,ft,dft) 
!     -------------------------------------------------
!
! Fourier integration over histogram data s(i) vs w(i)
! including error determinateion if data error is in ds(i)
! 
! spacing may be non-unigorm
!
                                                                       
      implicit none

      integer,          intent(in)    :: n                   ! length of data vectors
      double precision, intent(in)    :: w(n), s(n), ds(n)   ! omega, data, data-error
      double precision, intent(in)    :: t                   ! time parameter (result)
      double precision, intent(out)   :: ft1,ft2,ft,dft      ! real, imaginary, modulus, err(modulus) FT-values
 
      double precision                :: a1(n),a2(n),swt(n),cwt(n)
      double precision                :: t2, dft2  
      integer                         :: i
                                                                       
                                                                        
      if (t.ne.0.0) then 
                                                                        
        do i=1,n 
          swt(i)=sin(w(i)*t) 
          cwt(i)=cos(w(i)*t) 
        end do 
        t2=1.0d0/(t*t) 
                                                                        
        a1(1)=-swt(1)/t+(cwt(1)-cwt(2))/(w(2)-w(1))*t2 
        a2(1)=cwt(1)/t+(swt(1)-swt(2))/(w(2)-w(1))*t2 
        do i=2,n-1 
          a1(i)=((cwt(i)-cwt(i+1))/(w(i+1)-w(i)) + (cwt(i-1)-cwt(i))/(w(i-1)-w(i)))*t2              
          a2(i)=((swt(i)-swt(i+1))/(w(i+1)-w(i)) + (swt(i-1)-swt(i))/(w(i-1)-w(i)))*t2              
        end do 
        a1(n)=swt(n)/t+(cwt(n-1)-cwt(n))/(w(n-1)-w(n))*t2 
        a2(n)=-cwt(n)/t+(swt(n-1)-swt(n))/(w(n-1)-w(n))*t2 
                                                                        
      else 
                                                                        
        a1(1)=(w(2)-w(1))*0.5d0 
        a2(1)=0.0d0 
        do i=2,n 
          a1(i)=(w(i+1)-w(i-1))*0.5d0 
          a2(i)=0.0d0 
        end do 
        a1(n)=(w(n)-w(n-1))*0.5d0 
        a2(n)=0.0d0 
                                                                        
      end if 
                                                                        
      ft1=0.0d0 
      ft2=0.0d0 
      do i=1,n 
        ft1=ft1+a1(i)*s(i) 
        ft2=ft2+a2(i)*s(i) 
      end do 
      ft=sqrt(ft1*ft1+ft2*ft2) 
                                                                        
      dft2=0.0d0 
      do i=1,n 
        dft2=dft2+((a1(i)*ft1+a2(i)*ft2)/ft*ds(i))**2 
      end do 
      dft=sqrt(dft2) 
                                                                        
      return
 
      END  SUBROUTINE uni_fourier                                        
                                                                        
                                                                        
   end module unift
       
