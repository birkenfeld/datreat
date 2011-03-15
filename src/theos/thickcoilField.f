C*DS
C*ED
       Subroutine Thick_Coil_Field(r,xa,xe,rii,raa,anii, B)
C      ----------------------------------------------------
C      Feld einer dicken zylindersymmetrischen Spule
C      Eingabedaten:
C      r(1..3)     = Aufpunkt
C      xa(1..3)    = Koordinaten des Achsanfangs der Spule
C      xe(1..3)    = Koordinaten des Achsendes   der Spule
C      rii         = Innenradius der Spule
C                    falls negativ ==> symmetrisch konische Stromv.
C      raa         = Aussenradius der Spule
C      anii        = Amperewindungsprodukt der Spule
C      Ausgabedaten:
C      B(1..3)     = Magnetfeld(Tesla) am Ort r(1..3)
C
C Alle Angaben in SI-Einheiten ( Laenge-->m,Strom-->A,Mag.Feld-->Tesla )
C-----------------------------------------------------------------------
C
       implicit real*8 (a-h,o-z)
       dimension xa(3),xe(3),r(3),b(3)
       dimension xn(3),xr(3)
       parameter( Pihalb=1.57079632679489 )
       external baxi,bradf,baxic,bradfc
       common/cbfitr/ani,xlen,ri,ra,rr,zz
       common/tceps/ Epsilon, MaxItA
C
C ---- Koordinatenumrechnung auf Spule in STD-Koordinaten
C
C Achslaenge und Richtungsvektor
       do i=1,3
         xn(i) = xe(i)-xa(i)
       enddo
       xlen = dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
       do i=1,3
         xn(i) = xn(i)/xlen
       enddo
C Achsabstand von Xa
       zz=    xn(1)*r(1) +xn(2)*r(2) +xn(3)*r(3) -
     *     (  xn(1)*xa(1)+xn(2)*xa(2)+xn(3)*xa(3) )
C Radialrichtung und Vektor
       do i=1,3
         xr(i) = r(i) - xa(i) - zz*xn(i)
       enddo
       rr = dsqrt( xr(1)**2+xr(2)**2+xr(3)**2 )
       if(rr.gt.0.d0) then
         do i=1,3
           xr(i) = xr(i)/rr
         enddo
       endif
       if(rii.lt.0) zz = zz-xlen/2
 
C
C --- Feldberechnung im STD- Koordinatensystem ---
C
       ani = anii
       ra  = raa
       ri  = rii
 
       if(rr.eq.0.d0 .and. zz.eq.0.d0) zz = Epsilon
       if(rr.eq.0.d0 .and. zz.eq.xlen) zz = zz + Epsilon
 
       if(ri.ge.0.d0) then
C -- Axialfeld:
cccc     ba=  2 * Gauss8(baxi ,pihalb,3*pihalb)
         ba=  2 * AdaptivIntegral2
     *         (baxi,pihalb,3*pihalb,Epsilon,Maxita ,ErrorAccu)
C -- Radialfeld
cccc     br=  2 * Gauss8(bradf,pihalb,3*pihalb)
         br=  2 * AdaptivIntegral2
     *         (bradf,pihalb,3*pihalb,Epsilon,Maxita ,ErrorAccu)
       else
C -- Axialfeld:
cccc     ba=  2 * Gauss8(baxic,pihalb,3*pihalb)
         ba=  2 * AdaptivIntegral2
     *         (baxic,pihalb,3*pihalb,Epsilon,Maxita ,ErrorAccu)
C -- Radialfeld
cccc     br=  2 * Gauss8(bradfc,pihalb,3*pihalb)
         br=  2 * AdaptivIntegral2
     *         (bradfc,pihalb,3*pihalb,Epsilon,Maxita ,ErrorAccu)
       endif
C
C --- Ruecktransformation auf aktuelle Koordinaten ---
C
      do i=1,3
        b(i) = ( ba*xn(i) + br*xr(i) ) * ani
      enddo
 
      return
      end
 
 
C*DS
C*ED
 
       Function bradf(phi)
c-----------------------------------------------------------------------
c      Integrand function of the phi integration for the radial
c      flux density component
c      All parameters must be supplied by the common block /cbfitr/
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cbfitr/ani,xlen,ri,ra,r ,z
       Parameter ( Pi  =3.14159265358979312d0)
 
       area  = xlen*(Ra-Ri)
       rq    = r*r
       a1q   = ri*ri
       a2q   = ra*ra
       z1q   = z*z
       z2q   =(z-xlen)**2
       s     = dsin(phi)
       rs    = r*s
       a11   = dsqrt(rq+a1q+z1q-2*ri*rs)
       a12   = dsqrt(rq+a2q+z1q-2*ra*rs)
       a21   = dsqrt(rq+a1q+z2q-2*ri*rs)
       a22   = dsqrt(rq+a2q+z2q-2*ra*rs)
       erg   = s* ( (a11-a12)+(a22-a21) +
     *              rs * dlog( (a11+ri-rs)*(a22+ra-rs) /
     *                        ((a12+ra-rs)*(a21+ri-rs)) ) )
       bradf  = erg*1d-7/area
 
       return
       end
 
 
C*DS
C*ED
       Function baxi(phi)
c-----------------------------------------------------------------------
c      Integrand function of the phi integration for the axial
c      flux density component
c      All parameters must be supplied by the common block /cbfitr/
c-----------------------------------------------------------------------
       implicit real*8  (a-h,o-s,u-z)
       common/cbfitr/ani,xlen,ri,ra,r ,z
       complex*16 bc,rarsbc,rirsbc,carg1,carg2,carg11,carg22,cll,
     *            cdsqrt,dcmplx,cdlog
       Parameter ( Pi  =3.14159265358979312d0)
 
       area  = xlen*(Ra-Ri)
       rq    = r*r
       s     = dsin(phi)
       c     = dsqrt(1.d0-s*s)
       rc    = r*c
       rs    = r*s
       bc    = dcmplx(0.d0,rc)
       rarsbc= rs+bc-ra
       rirsbc= rs+bc-ri
       ars1  = 2*ri*rs
       aq1   = ri*ri
       ars2  = 2*ra*rs
       aq2   = ra*ra
       z01   = 0
       zz01  = z01-z
       zz0q1 = zz01**2
       z02   = xlen
       zz02  = z02-z
       zz0q2 = zz02**2
       aaa11 = dsqrt(rq+aq1-ars1+zz0q1)
       aaa12 = dsqrt(rq+aq1-ars1+zz0q2)
       aaa21 = dsqrt(rq+aq2-ars2+zz0q1)
       aaa22 = dsqrt(rq+aq2-ars2+zz0q2)
       carg1 = (aaa22*zz02+bc*(rarsbc)-zz0q2)
     *        /(aaa21*zz01+bc*(rarsbc)-zz0q1)
       carg2 = (aaa11*zz01+bc*(rirsbc)-zz0q1)
     *        /(aaa12*zz02+bc*(rirsbc)-zz0q2)
       carg11= rarsbc / rarsbc
       carg22= rirsbc / rirsbc
       rars  = ra-rs
       rirs  = ri-rs
       aa2212= (aaa22+rars)/(aaa12+rirs)
       aa2111= (aaa21+rars)/(aaa11+rirs)
       t1    =( dlog( aa2212 )*z02  - dlog( aa2111 )*z01 ) -
     *          dlog( aa2212 / aa2111 ) * z
c ---- 4 complex logarithms are required to be sure to ----------------
c      stay in the same Riemman sheet, if the field
c      inside the current carrying zone is NOT to be computed
c      cll may be evaluated as: cdlog(carg1/carg11)+cdlog(carg2/carg22)
 
        cll   = (cdlog(carg1) - cdlog(carg11)) +
     *          (cdlog(carg2) - cdlog(carg22))
c ---------------------------------------------------------------------
        t1    = t1 + dreal(cll)*rs - dimag(cll)*rc
        baxi  = t1*1d-7/area
 
       return
       end
 
C*DS
C*ED
       Function Gauss8( F,xu,xo)
C-----------------------------------------------------------------------
C      8-Punkte Gauss-Integration  : Int(f,xu..xo)
C-----------------------------------------------------------------------
       Parameter ( ndim = 8 )
       Parameter ( ndim2=ndim/2)
       implicit real*8 (a-h,o-z)
       dimension A(ndim2),X(ndim2)
       Data A / 0.362683783378362d0,
     *          0.313706645877887d0,
     *          0.222381034453374d0,
     *          0.101228536290376d0/
       Data X / 0.183434642495650d0,
     *          0.525532409916329d0,
     *          0.796666477413627d0,
     *          0.960289856497536d0/
 
       xave = (xo+xu)*0.5d0
       range= (xo-xu)*0.5d0
       sum = 0.d0
       do i=1,ndim2
         sum = sum + A(i)*( F(xave+range*X(i))+ F(xave-range*X(i)) )
       enddo
       Gauss8 = sum * range
 
       return
       end
C*DS
C*ED
       Function AdaptivIntegral2(F,A,B,Epsilon,Maxiter,ErrorAccu)
c      ==========================================================       --
c
c      lokal adaptives Integrationsverfahren
c
       Implicit real*8 (a-h,o-z)
c
       Parameter( MaxStack = 200 )
c                 --------------->  Stacktiefe
       Dimension s(0:MaxStack), xa(0:MaxStack), xb(0:MaxStack)
c
       logical Cray
       real*4 xxxx,yyyy,ptxf                        !! aix
       common/outlev/iot,ibild,ierrs,inka, iibuf, xxxx,yyyy,ptxf(20)  !! aix
c
       External F
c
c      -------------------------------------------------------------------
 
 
         IterationCounter = 0
         erg              = 0.0d0
         itop      = 0
         xa(itop)  = A
         xb(itop)  = B
         s(itop)   = Gauss8( F, A, B )
 
         Do i=1,MaxIter
           IF (itop.ge.(MaxStack-2)) THEN
              erg = erg + s(itop)
              xbb = xb(itop)
              itm = itop-1
              do itop=itm,0,-1
                 if(xa(itop).eq.xbb) goto 1
              enddo
1             continue
              Write(6,*)'WARNING! AdaptiveIntegral2 Stack Overflow!'
           ELSE
              IterationCounter = IterationCounter + 1
              itop             = itop +1
              xa(itop) = (xa(itop-1)+xb(itop-1))*0.5d0
              xb(itop) =  xb(itop-1)
              s(itop)  = Gauss8( F, xa(itop), xb(itop) )
              itop             = itop +1
              xa(itop) =  xa(itop-2)
              xb(itop) =  xa(itop-1)
              s(itop)  = Gauss8( F, xa(itop), xb(itop) )
              Error   =  DABS((s(itop)+s(itop-1))-s(itop-2))
 
              IF (iot.gt.2) THEN
                 Write(6,'(1x,i3,i5,4E13.6)')itop,IterationCounter,
     *                 xa(itop),xb(itop),(s(itop)+s(itop-1)),Error
              ENDIF
 
              IF (Error.lt.Epsilon) THEN
                 erg = erg + s(itop)+s(itop-1)
                 ErrorAccu = ErrorAccu + Error
                 itop = itop-1
                 xbb = xb(itop)
                 itop = itop-1
                 itm = itop-1
                 do itop=itm,0,-1
                    if(xa(itop).eq.xbb) goto 2
                 enddo
2                continue
              ENDIF
            ENDIF
            IF (itop.le.0) goto 3
         enddo
         Write(6,*)                                                     '
     *    'Adaptivintegral fatal Error Iterationnumber exceeded!'
3        continue
 
 
         AdaptivIntegral2 = erg
 
         RETURN
 
       END
 
C*DS
C*ED
       Function bradfc(phi)
C-----------------------------------------------------------------------
C      Integrand-Funktion der phi-Integration fuer den radialen
C      Feldanteil. Es werden Stromdichteskalierung und mue0
C      als Faktoren angebracht.
C      Alle Parameter werden aus Common-Block /cbfitr/
C      entnommen.
C-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cbfitr/ani,xlen,ri,ra,rr,zz
       Parameter       ( amu0=1.256637061d-6 )
       Parameter       ( Pi  =3.141592654d0  )
 
       Flaeche = Xlen*Ra/2
       erg     = bradikc(xlen,ra,ri,rr,zz,phi)
       bradfc  = erg*amu0/(4*Pi*Flaeche)
 
       return
       end
 
 
C*DS
C*ED
       Function baxic(phi)
C-----------------------------------------------------------------------
C      Integrand-Funktion der phi-Integration fuer den radialen
C      Feldanteil. Es werden Stromdichteskalierung und mue0
C      als Faktoren angebracht.
C      Alle Parameter werden aus Common-Block /cbfitr/
C      entnommen.
C-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       common/cbfitr/ani,xlen,ri,ra,rr,zz
       Parameter       ( amu0=1.256637061d-6 )
       Parameter       ( Pi  =3.141592654d0  )
 
       Flaeche = Xlen*Ra/2
       erg  = baxiakc(xlen,ra,ri,rr,zz,phi)
       baxic= erg*amu0/(4*Pi*Flaeche)
 
       return
       end
 
 
C*DS
C*ED
       Function Baxiakc(Zlen,Ra,Ri,rr,zz,phi)
C-----------------------------------------------------------------------
C  Axialer Feldanteil des phi-Sektors der Spule im STD-System.
C  D.h. Spulenanfang bei  (0,  0,  0   )
C       Spulenende   bei  (0,  0,  Zlen)
C       Innenradius   =   ri
C       Aussenradius  =   ra
C       Aufpunkt     bei  (rr, 0,  zz  )  (o.d radialsymm. Aeuquiv.)
C
C       Baxiak =
C       int(int((a**2-a*r*s)/(r**2+a**2+(z-z0)**2-2*r*a*s)**(3/2),
C              z0=-t*a..t*a)
C           a=0..ra )
C       mit s=sin(phi), z02=Zlen, a1=ri, a2=ra
C       Dieses Integral laesst sich mit Maple in der oben
C       geschriebnen Form loesen und in die unten enthaltene
C       Fortransequenz uebersetzen.
C
C       Um Fallunterscheidungen etc. beim Auftreten von negativen
C       sqrt und alog - Argumenten zu vermeiden wird Complex
C       gerechnet.
C
C       Um das Feld zu erhalten muss die phi-Integration noch
C       numerisch durchgefuehrt werden (s.o.)
C
C-----------------------------------------------------------------------
C
       implicit real*8 (a-h,o-s,u-z)
       implicit complex*16 (t)
       complex*16 cdsqrt,cdlog,dcmplx
       real*8 t
       Parameter (Flmin=1d-10)
 
       z02  = Zlen
       a1   = Ri
       a2   = Ra
       z    = zz
       if(dabs(z).lt.flmin) z = flmin
       r    = rr+1d-9
C                ----> vermeide clog(0,0) bei r=0
       s    = dsin(phi)
       rs   = r*s
       t    = Zlen*0.5d0/Ra
 
c --- Stromverteilung = symmetrischer Diabolo ( Ri=0 (kommt nicht vor))
      t1 = -rs
      t4 =  dcmplx(0.d0,dsqrt(r-rs))
      t6 =  dsqrt(rs+r)
      t7 = t4*t6
      t9 =cdlog(ra+t1+t7)
      t10 = rs**2
      t11 = r**2
      t14 = cdsqrt(t10-t11)
      t15 = t**2
      t16 = 1+t15
      t17 = cdsqrt(t16)
      t18 = t17*t16
      t20 = 2*t15*t10
      t21 = t15*rs*t4*t6
      t22 = 2*t21
      t24 = -t11*t15
      t25 = z*t*rs
      t26 = -2*t25
      t27 = z*t*t4*t6
      t28 = -2*t27
      t29 = z**2
      t31 = cdsqrt(t20+t22+t24+t26+t28+t29)
      t32 = 2*t25
      t33 = 2*t27
      t35 = cdsqrt(t20+t22+t24+t32+t33+t29)
      t36 = -2*t21
      t38 = cdsqrt(t20+t36+t24+t32+t28+t29)
      t41 = ra**2
      t42 = t15*t41
      t43 = z*t*ra
      t46 = -2*rs*ra
      t48 = cdsqrt(t42+t41-2*t43+t46+t29+t11)
      t50 = -t27
      t51 = ra*t4*t6
      t52 = rs*t4*t6
      t53 = -t52
      t54 = -t10
      t55 = t15*rs*ra
      t56 = -t25
      t57 = t15*ra*t4*t6
      t58 = -t43
      t60 =cdlog(t31*t48+t50+t29+t51+t53+t54+t11+t55+t56+t57+t58)
      t62 = cdsqrt(t20+t36+t24+t26+t33+t29)
      t67 =cdlog(-ra+rs+t7)
      t68 = t*t15
      t72 = -t51
      t73 = -t57
      t75 =cdlog(t62*t48+t27+t29+t72+t52+t54+t11+t55+t56+t73+t58)
      t84 = cdsqrt(t42+t41+2*t43+t46+t29+t11)
      t87 =cdlog(t35*t84+t27+t29+t51+t53+t54+t11+t55+t25+t57+t43)
      t93 =cdlog(t38*t84+t50+t29+t72+t52+t54+t11+t55+t25+t73+t43)
      t194 = t67*t14*z*t18*t31*t62*t38-t87*t14*z*t18*t31*t62*t38-t9*t11*
     +t68*t18*t31*t62*t35+t93*t11*t68*t18*t31*t62*t35-t9*t11*t*t18*t31*t
     +62*t35+2*t75*t14*t*rs*t18*t31*t35*t38+2*t9*t10*t*t18*t31*t35*t38+t
     +87*t11*t*t18*t31*t62*t38+t93*t11*t*t18*t31*t62*t35-t9*t14*z*t18*t1
     +5*t31*t62*t35+t93*t14*z*t18*t15*t31*t62*t35-2*t9*t14*t68*rs*t18*t3
     +1*t35*t38-2*t75*t10*t*t18*t31*t35*t38+2*t9*t10*t68*t18*t31*t62*t35
     ++2*t75*t14*t68*rs*t18*t31*t35*t38-2*t93*t10*t68*t18*t31*t62*t35+t6
     +7*t14*z*t18*t15*t31*t62*t38-t87*t14*z*t18*t15*t31*t62*t38-2*t9*t14
     +*t*rs*t18*t31*t62*t35+2*t67*t14*t68*rs*t18*t31*t62*t38-2*t9*t14*t6
     +8*rs*t18*t31*t62*t35+2*t93*t14*t68*rs*t18*t31*t62*t35+2*t93*t14*t*
     +rs*t18*t31*t62*t35+t9*t14*z*t18*t15*t31*t35*t38-2*t87*t10*t*t18*t3
     +1*t62*t38+2*t67*t10*t*t18*t31*t62*t38+2*t67*t10*t68*t18*t62*t35*t3
     +8+t60*z*rs*t18*t15*t62*t35*t38-t9*t11*t68*t18*t31*t35*t38-2*t60*t1
     +4*t68*rs*t18*t62*t35*t38-2*t87*t10*t68*t18*t31*t62*t38-2*t60*t10*t
     +68*t18*t62*t35*t38+2*t67*t10*t68*t18*t31*t62*t38+t75*z*rs*t18*t31*
     +t35*t38-t67*t14*z*t18*t62*t35*t38-t9*z*rs*t18*t31*t35*t38+t75*t11*
     +t68*t18*t31*t35*t38-t67*t11*t68*t18*t62*t35*t38+2*t67*t14*t*rs*t18
     +*t62*t35*t38+t60*t14*z*t18*t15*t62*t35*t38-t87*z*rs*t18*t31*t62*t3
     +8+t67*z*rs*t18*t31*t62*t38+t60*t14*z*t18*t62*t35*t38-2*t60*t14*t*r
     +s*t18*t62*t35*t38-t67*z*rs*t18*t62*t35*t38-t9*z*rs*t18*t15*t31*t35
     +*t38-t87*z*rs*t18*t15*t31*t62*t38+t67*z*rs*t18*t15*t31*t62*t38+t75
     +*z*rs*t18*t15*t31*t35*t38-t67*z*rs*t18*t15*t62*t35*t38+2*t67*t14*t
     +68*rs*t18*t62*t35*t38-t93*z*rs*t18*t31*t62*t35-t75*t14*z*t18*t15*t
     +31*t35*t38+t93*t14*z*t18*t31*t62*t35+t60*t11*t68*t18*t62*t35*t38+t
     +9*z*rs*t18*t31*t62*t35+2*t9*t10*t*t18*t31*t62*t35-2*t93*t10*t*t18*
     +t31*t62*t35-2*t9*t14*t*rs*t18*t31*t35*t38-t67*t11*t*t18*t31*t62*t3
     +8-t9*t14*z*t18*t31*t62*t35-2*t60*t10*t*t18*t62*t35*t38+t9*z*rs*t18
     +*t15*t31*t62*t35-t67*t14*z*t18*t15*t62*t35*t38-t93*z*rs*t18*t15*t3
     +1*t62*t35+2*t67*t10*t*t18*t62*t35*t38
      t218 = dlog(2.d0)
      t219 = t218*t14*z*t18*t15*t62*t35*t38
      t221 = 2*t218*t14*t*rs*t18*t31*t35*t38
      t223 = -t218*z*rs*t18*t15*t31*t62*t38
      t224 = t218*t14*z*t18*t15*t31*t62*t35
      t225 = t218*z*rs*t18*t15*t62*t35*t38
      t227 = -t218*z*rs*t18*t15*t31*t62*t35
      t229 = 2*t218*t14*t68*rs*t18*t31*t35*t38
      t231 = 2*t218*t14*t68*rs*t18*t31*t62*t35
      t233 = -2*t218*t10*t68*t18*t31*t62*t35
      t235 = -2*t218*t10*t68*t18*t62*t35*t38
      t236 = t218*t11*t*t18*t62*t35*t38
      t238 = -2*t218*t10*t*t18*t31*t35*t38
      t240 = -t218*t14*z*t18*t15*t31*t62*t38
      t241 = t218*t11*t68*t18*t62*t35*t38
      t242 = t218*t14*z*t18*t62*t35*t38
      t244 = -t218*t14*z*t18*t31*t35*t38
      t246 = -t218*t14*z*t18*t15*t31*t35*t38
      t247 = t218*t11*t68*t18*t31*t35*t38
      t248 = t218*z*rs*t18*t62*t35*t38
      t251 = -2*t218*t10*t*t18*t62*t35*t38
      t253 = -2*t218*t10*t68*t18*t31*t35*t38
      t254 = t218*t11*t*t18*t31*t35*t38
      t255 = t218*z*rs*t18*t15*t31*t35*t38
      t257 = 2*t218*t14*t*rs*t18*t31*t62*t35
      t259 = -t218*t14*z*t18*t31*t62*t38
      t261 = -2*t218*t10*t68*t18*t31*t62*t38
      t262 = t218*t11*t*t18*t31*t62*t38
      t264 = -2*t218*t10*t*t18*t31*t62*t35
      t265 = t218*t11*t68*t18*t31*t62*t38
      t267 = -2*t218*t14*t*rs*t18*t62*t35*t38
      t269 = -2*t218*t10*t*t18*t31*t62*t38
      t271 = -2*t218*t14*t*rs*t18*t31*t62*t38
      t273 = -2*t218*t14*t68*rs*t18*t31*t62*t38
      t274 = t218*t11*t*t18*t31*t62*t35
      t276 = -t218*z*rs*t18*t31*t62*t35
      t277 = t218*t11*t68*t18*t31*t62*t35
      t278 = t218*z*rs*t18*t31*t35*t38
      t279 = t218*t14*z*t18*t31*t62*t35
      t281 = -2*t218*t14*t68*rs*t18*t62*t35*t38
      t283 = -t218*z*rs*t18*t31*t62*t38
      t287 = 2*t15*ra
      t288 = z*t
      t289 = 2*t288
      t290 = 2*ra
      t291 = -2*rs
      t293 =cdlog(2*t17*t84+t287+t289+t290+t291)
      t294 = t15**2
      t295 = t*t294
      t302 = -2*t288
      t304 =cdlog(2*t17*t48+t287+t302+t290+t291)
      t326 = t278+t277+t264+t259+t261+t254+t276+t255+t257+6*rs*t304*t68*
     +t31*t62*t35*t38+t244+t235+t231+t233+t229+t251+t224+t225+t227+t219+
     +t221+t223+t75*t11*t*t18*t31*t35*t38+t87*t11*t68*t18*t31*t62*t38+2*
     +t9*t10*t68*t18*t31*t35*t38+2*z*t15*t293*t31*t62*t35*t38+t60*z*rs*t
     +18*t62*t35*t38+6*rs*t293*t68*t31*t62*t35*t38-2*t304*z*t31*t62*t35*
     +t38+t279+t246+t9*t14*z*t18*t31*t35*t38-t75*t14*z*t18*t31*t35*t38+t
     +281+t273+t274+t271+t269+t265+t267+t262+t253+t247+t248+t283+t241+t2
     +42+t236+t238+t240-2*t75*t10*t68*t18*t31*t35*t38-t9*t11*t*t18*t31*t
     +35*t38-t67*t11*t68*t18*t31*t62*t38-2*t87*t14*t*rs*t18*t31*t62*t38-
     +t67*t11*t*t18*t62*t35*t38-2*t87*t14*t68*rs*t18*t31*t62*t38+t60*t11
     +*t*t18*t62*t35*t38+2*t48*t18*t*t31*t62*t35*t38+2*t84*t18*t*t31*t62
     +*t35*t38-2*z*t15*t304*t31*t62*t35*t38+4*t304*t*rs*t31*t62*t35*t38+
     +2*t67*t14*t*rs*t18*t31*t62*t38+2*t293*z*t31*t62*t35*t38+2*rs*t304*
     +t295*t31*t62*t35*t38+4*t293*t*rs*t31*t62*t35*t38+2*rs*t293*t295*t3
     +1*t62*t35*t38
      t327 = t194+t326
      t328 = t16**2
      t331 = t17/t16/t328
      t332 = 1/t31
      t333 = 1/t62
      t334 = 1/t35
      t335 = 1/t38
      t339 = cdsqrt(t29+t11)
      t342 =cdlog(t35*t339+t11+t27+t53+t54+t25+t29)
      t347 =cdlog(t31*t339+t11+t50+t53+t54+t56+t29)
      t350 =cdlog(t1+t7)
      t357 =cdlog(t38*t339+t11+t50+t52+t54+t25+t29)
      t363 =cdlog(rs+t7)
      t380 =cdlog(t62*t339+t11+t27+t52+t54+t56+t29)
      t430 = -2*t350*t14*t68*rs*t18*t31*t35*t38+t350*t14*z*t18*t31*t35*t
     +38+t254+2*t380*t14*t68*rs*t18*t31*t35*t38+t347*t11*t68*t18*t62*t35
     +*t38+t255-t350*t11*t*t18*t31*t35*t38+2*t350*t10*t68*t18*t31*t35*t3
     +8+t244+2*t363*t10*t*t18*t62*t35*t38-2*t347*t10*t*t18*t62*t35*t38-2
     +*t350*t14*t*rs*t18*t31*t35*t38+t380*t11*t68*t18*t31*t35*t38+t235+t
     +231+t233+t229+t251-t380*t14*z*t18*t31*t35*t38+t224+2*t363*t10*t68*
     +t18*t62*t35*t38+t225+t227-t363*t11*t68*t18*t62*t35*t38-2*t380*t10*
     +t*t18*t31*t35*t38+t219+t221+t223-t363*t14*z*t18*t15*t62*t35*t38+t3
     +47*t11*t*t18*t62*t35*t38-2*t347*t10*t68*t18*t62*t35*t38-2*t357*t10
     +*t68*t18*t31*t62*t35-2*t350*t14*t68*rs*t18*t31*t62*t35-t350*t14*z*
     +t18*t15*t31*t62*t35-t350*z*rs*t18*t15*t31*t35*t38+t246-2*t380*t10*
     +t68*t18*t31*t35*t38+t380*t11*t*t18*t31*t35*t38-t363*z*rs*t18*t62*t
     +35*t38-t363*t14*z*t18*t62*t35*t38+t347*t14*z*t18*t62*t35*t38+t347*
     +z*rs*t18*t62*t35*t38-t342*t14*z*t18*t15*t31*t62*t38-t350*t11*t68*t
     +18*t31*t35*t38+t363*t14*z*t18*t15*t31*t62*t38+t253+t247+t248-t342*
     +z*rs*t18*t15*t31*t62*t38+2*t350*t10*t*t18*t31*t35*t38-t363*t11*t*t
     +18*t62*t35*t38+t363*z*rs*t18*t15*t31*t62*t38+t241+t242+t350*t14*z*
     +t18*t15*t31*t35*t38+t236+t238+t240+2*t350*t10*t68*t18*t31*t62*t35+
     +2*t357*t14*t68*rs*t18*t31*t62*t35+t380*z*rs*t18*t15*t31*t35*t38-2*
     +t342*t14*t68*rs*t18*t31*t62*t38+t347*t14*z*t18*t15*t62*t35*t38-t38
     +0*t14*z*t18*t15*t31*t35*t38+2*t380*t14*t*rs*t18*t31*t35*t38
      t497 = 2*t17*t339
      t499 =cdlog(t497+t291+t302)
      t509 =cdlog(t497+t291+t289)
      t525 = 6*rs*t499*t68*t31*t62*t35*t38+t278+t277-t350*t14*z*t18*t31*
     +t62*t35+t357*t14*z*t18*t31*t62*t35-t350*z*rs*t18*t31*t35*t38+t380*
     +z*rs*t18*t31*t35*t38-t350*t11*t68*t18*t31*t62*t35+2*t509*z*t31*t62
     +*t35*t38+t264-t357*z*rs*t18*t31*t62*t35+2*t363*t14*t*rs*t18*t31*t6
     +2*t38-2*t342*t14*t*rs*t18*t31*t62*t38-2*t342*t10*t*t18*t31*t62*t38
     ++t357*t14*z*t18*t15*t31*t62*t35+t259+t261+t342*t11*t68*t18*t31*t62
     +*t38-2*t357*t10*t*t18*t31*t62*t35-t363*t11*t*t18*t31*t62*t38+t276+
     +2*t363*t10*t68*t18*t31*t62*t38-2*t342*t10*t68*t18*t31*t62*t38+2*z*
     +t15*t509*t31*t62*t35*t38+t257-t363*z*rs*t18*t15*t62*t35*t38+2*t363
     +*t14*t68*rs*t18*t62*t35*t38+2*rs*t499*t295*t31*t62*t35*t38+2*t363*
     +t14*t*rs*t18*t62*t35*t38+2*rs*t509*t295*t31*t62*t35*t38+6*rs*t509*
     +t68*t31*t62*t35*t38+4*t339*t18*t*t31*t62*t35*t38+t347*z*rs*t18*t15
     +*t62*t35*t38-2*t499*z*t31*t62*t35*t38-2*z*t15*t499*t31*t62*t35*t38
     ++t363*z*rs*t18*t31*t62*t38-t342*z*rs*t18*t31*t62*t38+t350*z*rs*t18
     +*t31*t62*t35+t357*t11*t68*t18*t31*t62*t35-t350*t11*t*t18*t31*t62*t
     +35+2*t363*t14*t68*rs*t18*t31*t62*t38-t363*t11*t68*t18*t31*t62*t38-
     +2*t347*t14*t*rs*t18*t62*t35*t38+2*t350*t10*t*t18*t31*t62*t35+t363*
     +t14*z*t18*t31*t62*t38+2*t357*t14*t*rs*t18*t31*t62*t35+t342*t11*t*t
     +18*t31*t62*t38+4*t499*t*rs*t31*t62*t35*t38+2*t363*t10*t*t18*t31*t6
     +2*t38+t279+t357*t11*t*t18*t31*t62*t35+t281-2*t350*t14*t*rs*t18*t31
     +*t62*t35+t273+t274-2*t347*t14*t68*rs*t18*t62*t35*t38+t271+t269+t26
     +5+t267+t262+t283+t350*z*rs*t18*t15*t31*t62*t35+4*t509*t*rs*t31*t62
     +*t35*t38-t342*t14*z*t18*t31*t62*t38-t357*z*rs*t18*t15*t31*t62*t35
      t526 = t430+t525
      t529 = t327*t331*t332*t333*t334*t335/2-t526*t331*t332*t333*t334*t3
     +35/2
 
       baxiakc= t529
 
       return
       end
 
 
 
C*DS
C*ED
       Function Bradikc(Zlen,Ra,Ri,rr,zz,phi)
C-----------------------------------------------------------------------
C  Radialer Feldanteil des phi-Sektors der Spule im STD-System.
C  D.h. Spulenanfang bei  (0,  0,  0   )
C       Spulenende   bei  (0,  0,  Zlen)
C       Innenradius   =   ri
C       Aussenradius  =   ra
C       Aufpunkt     bei  (rr, 0,  zz  )  (o.d radialsymm. Aeuquiv.)
C
C       Baxiak =
C       int(int((a*(z-z0)*s)/(r**2+a**2+(z-z0)**2-2*r*a*s)**(3/2),
C              z0=-t*a..t*a),
C           a=ri..ra)
C       mit s=sin(phi), z02=Zlen, a1=ri, a2=ra
C       Dieses Integral laesst sich mit Maple in der oben
C       geschriebenen Form loesen und in die unten enthaltene
C       Fortransequenz uebersetzen.
C
C       Um das Feld zu erhalten muss die phi-Integration noch
C       numerisch durchgefuehrt werden (s.o.)
C
C-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
 
       z02  = Zlen
       a1   = Ri
       a2   = Ra
       z    = zz
       r    = rr+1d-9
C                ----> vermeide clog(0,0) bei r=0
       s    = dsin(phi)
       rs   = r*s
       t    = Zlen*0.5d0/ra
 
c --- Stromverteilung = symmetrischer Diabolo ( Ri=0 (kommt nicht vor))
 
      t1 = t**2
      t2 = ra**2
      t3 = t1*t2
      t4 = z*t*ra
      t7 = -2*r*s*ra
      t8 = z**2
      t9 = r**2
      t11 = dsqrt(t3+t2-2*t4+t7+t8+t9)
      t12 = 1+t1
      t13 = dsqrt(t12)
      t14 = t13*t12
      t17 = t*t1
      t20 = dsqrt(t3+t2+2*t4+t7+t8+t9)
      t24 = 2*t1*ra
      t25 = z*t
      t26 = 2*t25
      t27 = 2*ra
      t29 = -2*r*s
      t31 = dlog(2*t13*t20+t24+t26+t27+t29)
      t36 = -2*t25
      t38 = dlog(2*t13*t11+t24+t36+t27+t29)
      t53 = t12**2
      t56 = t13/t12/t53
      t59 = ri**2
      t60 = t1*t59
      t61 = z*t*ri
      t64 = -2*r*s*ri
      t66 = dsqrt(t60+t59-2*t61+t64+t8+t9)
      t71 = dsqrt(t60+t59+2*t61+t64+t8+t9)
      t75 = 2*t1*ri
      t76 = 2*ri
      t78 = dlog(2*t13*t71+t75+t26+t76+t29)
      t84 = dlog(2*t13*t66+t75+t36+t76+t29)
      t100 = -s*(-t11*t14-z*t17*t31-z*t*t38+r*s*t31-z*t17*t38+r*s*t31*t1
     +-r*s*t38-r*s*t38*t1+t20*t14-z*t*t31)*t56+s*(-t66*t14-z*t17*t78-z*t
     +*t84+r*s*t78-z*t17*t84+r*s*t78*t1-r*s*t84-r*s*t84*t1+t71*t14-z*t*t
     +78)*t56
 
       bradikc= t100
 
       return
       end
