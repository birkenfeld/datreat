      FUNCTION th_cloizeaux (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formel (23) aus Cloizeaux, J. de Physique I, 3 (1993) p.1533 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 / 


         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)


         REAL            th_cloizeaux, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx
         DIMENSION       qq(3)
         DATA            zpi/6.283185/

         double precision w,d,q,n,ne,cloizeaux,t,a,wl4


 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

!
! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'cloizeaux'
         nparx =  3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_cloizeaux   = 0
           return
         endif
         npar = nparx
!        --------------> set the number of parameters
         parnam(1) = 'amplitud'          
         parnam(2) = 'wl4     '
         parnam(3) = 'dtube   '
!
         th_cloizeaux  = 0
         return
       endif
!
! ---- calculate theory here -----
        a       = pa(1)   ! amplitude (should be 1)
        wl4     = pa(2)   ! Rouse rate W*l**4
        d       = pa(3)   ! Tube diameter

        call parget('q       ',xh ,iadda,ier)
        if(ier.eq.0) then
          q       = xh
        else
          q       = 1.0  
        endif
 
        t = x 
 
        th_cloizeaux = a * cloizeaux(t,q,d,Wl4) 

         return
        end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formel (23) aus Cloizeaux, J. de Physique I, 3 (1993) p.1533 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function cloizeaux(t,q,d,Wl4)

  implicit none

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

!      integer maxit
!      double precision epsilon, errac, gwidth
!      parameter (epsilon=1d-2 , gwidth=4.0d0, maxit=1000)

      double precision FA_ker, FB_ker, FYZ

      double precision d, q, Wl4, t, cadapint_r, gamma
      double precision go,I1,I2,res,n1,p
      integer          j        
     
 
! for Gauss Laguerre Intergration !      
  integer, parameter :: ngl64 = 64

  double precision, dimension(ngl64), parameter :: xgl64 =  &
 [ 2.241587414670528002281188481904d-02 , &
   1.181225120967704797974664367103d-01 , &
   2.903657440180364839991300663853d-01 , &
   5.392862212279790393181449478122d-01 , &
   8.650370046481139446199550747097d-01 , &
   1.267814040775241398115708877694d+00 , &
   1.747859626059436252829963951288d+00 , &
   2.305463739307508718548070543894d+00 , &
   2.940965156725251840679468152108d+00 , &
   3.654752650207290527035397912093d+00 , &
   4.447266343313094356742550980161d+00 , &
   5.318999254496390343522109859188d+00 , &
   6.270499046923653912911064646334d+00 , &
   7.302370002587395747223498409519d+00 , &
   8.415275239483024194495218591205d+00 , &
   9.609939192796108035762882049548d+00 , &
   1.088715038388637214259455042018d+01 , &
   1.224776450424430161816236929072d+01 , &
   1.369270784554750515272993257464d+01 , &
   1.522298111152472884800826878343d+01 , &
   1.683966365264873721052883803916d+01 , &
   1.854391817085919052361962597108d+01 , &
   2.033699594873023550114981580640d+01 , &
   2.222024266595087653992215434709d+01 , &
   2.419510487593325398988644388024d+01 , &
   2.626313722711848578512602395478d+01 , &
   2.842601052750102729949977152675d+01 , &
   3.068552076752597177104858239845d+01 , &
   3.304359923643782912552021068053d+01 , &
   3.550232389114120958697877853505d+01 , &
   3.806393216564646826035731791502d+01 , &
   4.073083544445862636573186951324d+01 , &
   4.350563546642152985270318493173d+01 , &
   4.639114297861619207360539994238d+01 , &
   4.939039902562468667923580082269d+01 , &
   5.250669934134630165017927698046d+01 , &
   5.574362241327838046333571129115d+01 , &
   5.910506191901710660884874209180d+01 , &
   6.259526440015139559605501790119d+01 , &
   6.621887325124756438221376267101d+01 , &
   6.998098037714682922853465797223d+01 , &
   7.388718723248296321095740311351d+01 , &
   7.794367743446312031368797587057d+01 , &
   8.215730377831930429519586834224d+01 , &
   8.653569334945651821021627837527d+01 , &
   9.108737561313309014563671534928d+01 , &
   9.582194001552073209476721543650d+01 , &
   1.007502319695139796292592614512d+02 , &
   1.058845994687999493563604278509d+02 , &
   1.112392075244395820634847366376d+02 , &
   1.168304450513064984633866690774d+02 , &
   1.226774602685385765774196905650d+02 , &
   1.288028787692376725127536230543d+02 , &
   1.352337879495258278339804988789d+02 , &
   1.420031214899315190251400382911d+02 , &
   1.491516659000493885872934629324d+02 , &
   1.567310751326711612336169608139d+02 , &
   1.648086026551505229931901090252d+02 , &
   1.734749468364242745221528448669d+02 , &
   1.828582046914314636463427945096d+02 , &
   1.931511360370729114793855274172d+02 , &
   2.046720284850594559490644333426d+02 , &
   2.180318519353285163324523844483d+02 , &
   2.348095791713261647130555297254d+02   ]
   
  double precision, dimension(ngl64), parameter :: wgl64 = &
 [ 5.625284233902984574102185450625d-02 , &
   1.190239873124260278149035058888d-01 , &
   1.574964038621445238201964347055d-01 , &
   1.675470504157739478809044116588d-01 , &
   1.533528557792366180854547925640d-01 , &
   1.242210536093297445126137821926d-01 , &
   9.034230098648505773897410920156d-02 , &
   5.947775576835502421224699743972d-02 , &
   3.562751890403607185416573533689d-02 , &
   1.948041043116640604333738027148d-02 , &
   9.743594899382002240107960271379d-03 , &
   4.464310364166275292364824190544d-03 , &
   1.875359581323114826750128222516d-03 , &
   7.226469815750051227191087066195d-04 , &
   2.554875328334967097144484021795d-04 , &
   8.287143534396942179063224039877d-05 , &
   2.465686396788558745973370221716d-05 , &
   6.726713878829668527612545536274d-06 , &
   1.681785369964088897821201086506d-06 , &
   3.850812981546684414827886110704d-07 , &
   8.068728040990499790415110915505d-08 , &
   1.545723706757688828003703926158d-08 , &
   2.704480147617481409988607446575d-09 , &
   4.316775475427200912314164217177d-10 , &
   6.277752541761452201652955763263d-11 , &
   8.306317376288958063878932234380d-12 , &
   9.984031787220164055897291499895d-13 , &
   1.088353887116662685326123546986d-13 , &
   1.074017403441590186482804810228d-14 , &
   9.575737231574442105585024583281d-16 , &
   7.697028023648586098862970088678d-17 , &
   5.564881137454025366525461526815d-18 , &
   3.609756409010446498298646897469d-19 , &
   2.095095369548946234767748078464d-20 , &
   1.084793301097549361202631562742d-21 , &
   4.994699486363804115792443627063d-23 , &
   2.037836974598822310658063498146d-24 , &
   7.339537564278837039106367650444d-26 , &
   2.323783082198694261304776157500d-27 , &
   6.438234706908762420384341347468d-29 , &
   1.553121095788275270621507732395d-30 , &
   3.244250092019537314465982686763d-32 , &
   5.832386267836201501282207418418d-34 , &
   8.963254833102854061314281996410d-36 , &
   1.168703989550736241198572917330d-37 , &
   1.282055984359980381497645585305d-39 , &
   1.172094937405002291828122995135d-41 , &
   8.835339672328604981296921087574d-44 , &
   5.424955590306186594338147130486d-46 , &
   2.675542666678893828949793699428d-48 , &
   1.042917031411367078110598371701d-50 , &
   3.152902351957772623653127987116d-53 , &
   7.229541910647522339706610677135d-56 , &
   1.224235301230082264459594302229d-58 , &
   1.482168504901910411780618275588d-61 , &
   1.232519348814518808064365504012d-64 , &
   6.691499004571269526814069374687d-68 , &
   2.220465941850448995507392045782d-71 , &
   4.120946094738876249979125791383d-75 , &
   3.774399061896489170417616059552d-79 , &
   1.414115052917619417463147014780d-83 , &
   1.591833064041367917861031849654d-88 , &
   2.989484348860634307741312694817d-94 , &
   2.089063508436952770828154254400d-101 ]
  ! --- end GaussLaguerre coeffs -------------
  
      double precision alpha, beta 

      external FYZ    
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      
       gamma = Wl4/9

       Zav = q*q*d*d/2
       Z   = Zav
       Y   = q*q*sqrt(gamma*t)

        i1 = 0.0d0
        n1 = 0.0d0
   
        do j=1,ngl64
          p = xgl64(j)
          i1 = i1 + FYZ(p)*wgl64(j)
        enddo
        
      res     =   i1 + log(1+Zav)/Zav
      cloizeaux = res
!      write(6,*)t,res
      return

      end function cloizeaux




      double precision function FYZ(p)

      implicit none 

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-6, gwidth=4.0d0, maxit=1000)

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      double precision p, badapint_r, gu, go, fa, fb
      double precision FA_ker, FB_ker, FAB_ker
      double precision upbound, lowbound      
      double precision res, adapint_r

      external FA_ker, FB_ker, FAB_ker      
      external upbound, lowbound


      data Faf/1.0d0/
      data Fbf/1.0d0/

      Z = p*Zav

      gu = 0.0d0
      go = 1.0d0
    
!       call dtwodq(FAB_ker,gu,go,lowbound,upbound,
      !     *            epsilon, epsilon*10, 2, result, errac)

      res = adapint_r(FA_ker, gu, go, epsilon, maxit, errac ) 
      
       Fyz = res * Z
!       write(6,*)' 2: ',Fyz, errac

      return
      end


      double precision function FA_ker(a)
        implicit none
         double precision, intent(in) :: a
         double precision             :: a_cloik2
         double precision             :: adapint_r
         double precision, parameter  :: epsilon = 1d-9
         double precision             :: erracc
         integer, parameter           :: maxit = 300 
         common/cloik2/ a_cloik2
         double precision, external   :: FB_ker

         a_cloik2 = a

         FA_ker = adapint_r(FB_ker , a, 1d0, epsilon, maxit, erracc)
            
      end function FA_ker  

      

      double precision function FB_ker(b)
        implicit none
        double precision, intent(in) :: b
        double precision :: FAB_ker

        double precision :: a_cloik2
        common/cloik2/ a_cloik2

        FB_ker = FAB_ker(a_cloik2,b)
      
      end function FB_ker

      



      double precision function FAB_ker(a,b)
      
      implicit none
      double precision a, arg, derfc, ssap, ssa, ssaq, u1, arg2
      double precision b,u2,arg22, ssbp
      integer          p, plim, nn

      double precision Z,Y,Zav, FaF, Fbf
      common/cloik1/Z,Y,Zav, FaF, Fbf

      integer maxit
      double precision epsilon, errac, gwidth
      parameter (epsilon=1d-12, gwidth=4.0d0, maxit=1000)
      double precision edwp   
      parameter(edwp=0.564189583d0)        !! = 1/sqrt(pi)


      ssa = 0      
      plim = 0.5d0*gwidth*Y/Z + 2   

!      write(6,*)'Fa plim:',plim, ' a=',a
      do p=-plim,plim
        u1   = abs((a-2*p)*Z)
        arg2 = -((u1/Y)**2)
        if(arg2.lt.-150.0d0) arg2=-150.0d0
        ssap = edwp*exp(arg2)*Y - u1*derfc(u1/Y) 
  
        u2   = abs((b-2*p)*Z)
        arg22 = -((u2/Y)**2)
        if(arg22.lt.-150.0d0) arg22=-150.0d0
        ssbp = edwp*exp(arg22)*Y - u2*derfc(u2/Y) 
        ssa  = ssa + ssap - ssbp
!        write(6,*)p,ssa,ssap
      enddo

      arg = -Z*a - ssa
      FAB_ker = exp(arg) 

!      write(6,*)'FAB_ker=', FAB_ker      

      return
      end function FAB_ker





      RECURSIVE FUNCTION adapint_r (f, a, b, epsilon, maxiter, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
!#!      use outlev
      IMPLICIT real (8)(a - h, o - z) 
!                                                                       
      PARAMETER (maxstack = 200) 
!                 --------------->  stacktiefe                          
      DIMENSION s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
!                                                                       
!      LOGICAL cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
!      REAL xxxx, yyyy, ptxf 
                                                                     !! 
!      COMMON / outlev / iot, ibild, ierrs, inka, iibuf, xxxx, yyyy,     &
!      ptxf (20)                                                         
!                                                                       
      EXTERNAL f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = gaus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! adaptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
!         IF (iot.gt.2) then 
!            WRITE (6, '(1x,i3,i5,4e13.6)') itop, iterationcounter, xa ( &
!            itop) , xb (itop) , (s (itop) + s (itop - 1) ) , error      
!         ENDIF 
                                                                        
         IF (error.lt.epsilon) then 
            erg = erg + s (itop) + s (itop - 1) 
            erroraccu = erroraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'adapint_r fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
                                                                        
      adapint_r = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION adapint_r                          

!!#>       function gaus8a (f,xu,xo)
!!#> !-----------------------------------------------------------------------
!!#> !      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
!!#> !-----------------------------------------------------------------------
!!#>       parameter ( ndim = 8 )
!!#>       parameter ( ndim2=ndim/2)
!!#>       implicit real*8 (a-h,o-z)
!!#>       dimension a(ndim2),x(ndim2)
!!#>       data a / 0.362683783378362d0,&
!!#>                0.313706645877887d0,&
!!#>                0.222381034453374d0,&
!!#>                0.101228536290376d0/
!!#>       data x / 0.183434642495650d0,&
!!#>                0.525532409916329d0,&
!!#>                0.796666477413627d0,&
!!#>                0.960289856497536d0/
!!#> 
!!#>       xave = (xo+xu)*0.5d0
!!#>       range= (xo-xu)*0.5d0
!!#>       sum = 0.d0
!!#>       do i=1,ndim2
!!#>         sum = sum + a(i)*( f(xave+range*x(i))+ f(xave-range*x(i)) )
!!#>       enddo
!!#>       gaus8a = sum * range
!!#> 
!!#>       return
!!#>       end function gaus8a
