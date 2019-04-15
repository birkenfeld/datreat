 FUNCTION th_rpastar(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  applies rpa mixture tor the scattering function for a star-linear polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))
! 
      use theory_description 
      implicit none 
      real    :: th_rpastar
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor                                                                       
     double precision :: wl4star    ! rouse rate star                                                                     
     integer          :: narmeff    ! effective segment star arm                                                      
     double precision :: diffstar   ! centre-of-mass diffusino star                                                   
     double precision :: locr2_b    ! locrep2 parameter b                                                             
     double precision :: locr2_a    ! locrep2 parameter a                                                             
     double precision :: locr2_ta   ! locrep2 parameter tau                                                           
     double precision :: locr2_lz   ! locrep2 parameter lz                                                            
     double precision :: locr2_te   ! locrep2 parameter teps                                                          
     integer          :: nro_me     ! effective segments for Me-Modelling                                             
     double precision :: re_me      ! Re of Me (entangle strand Re)                                                   
     double precision :: wl4        ! rouse rate                                                                      
! the recin parameter representation 
     double precision :: q          ! scattering wavevector                                                           
     double precision :: temp       ! temperature                                                                     
     integer          :: narm       ! number of segments per star arm                                                 
     double precision :: l          ! segment length                                                                  
     integer          :: f          ! number of arms per star                                                         
     double precision :: phistar    ! volume fraction of star                                                         
     integer          :: nlin       ! number of segments of linear polymer                                            
     double precision :: philin     ! volume fraction of linera polymer                                               
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)                                                       
 
     double precision   :: t
    
     double precision   :: sqt, sqt0
     double precision   :: locrep2
     double precision   :: pstar0, pstar    
     double precision   :: plin0, plin
     double precision   :: dr, wx, lx
     double precision   :: Re_arm
     double precision   :: Sinv0, Sinv
     double precision   :: local_reptation2
     double precision   :: pericostar_sqt

     integer            :: ifix

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'rpastar'
       nparx =       12
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rpastar = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " applies rpa mixture tor the scattering function for a star-linear polymer micture and extends this into the dynamic regime (i.e. S(Q) --> S(Q,t))"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'wl4star '  ! rouse rate star                                                                 
        parnam ( 3) = 'narmeff '  ! effective segment star arm                                                      
        parnam ( 4) = 'diffstar'  ! centre-of-mass diffusino star                                                   
        parnam ( 5) = 'locr2_b '  ! locrep2 parameter b                                                             
        parnam ( 6) = 'locr2_a '  ! locrep2 parameter a                                                             
        parnam ( 7) = 'locr2_ta'  ! locrep2 parameter tau                                                           
        parnam ( 8) = 'locr2_lz'  ! locrep2 parameter lz                                                            
        parnam ( 9) = 'locr2_te'  ! locrep2 parameter teps                                                          
        parnam (10) = 'nro_me  '  ! effective segments for Me-Modelling                                             
        parnam (11) = 're_me   '  ! Re of Me (entangled strand Re)                                                
        parnam (12) = 'wl4     '  ! Wl4 (local entanglement strand Re)                                                   
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate star" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment star arm" !//cr//parspace//&
        th_param_desc( 4,idesc) = "centre-of-mass diffusino star" !//cr//parspace//&
        th_param_desc( 5,idesc) = "locrep2 parameter b" !//cr//parspace//&
        th_param_desc( 6,idesc) = "locrep2 parameter a" !//cr//parspace//&
        th_param_desc( 7,idesc) = "locrep2 parameter tau" !//cr//parspace//&
        th_param_desc( 8,idesc) = "locrep2 parameter lz" !//cr//parspace//&
        th_param_desc( 9,idesc) = "locrep2 parameter teps" !//cr//parspace//&
        th_param_desc(10,idesc) = "effective segments for Me-Modelling" !//cr//parspace//&
        th_param_desc(11,idesc) = "Re of Me (entangle strand Re)" !//cr//parspace//&
        th_param_desc(12,idesc) = "Wl4 entaglement strand" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  3,idesc) = "narm     > number of segments per star arm"
        th_file_param(  4,idesc) = "l        > segment length"
        th_file_param(  5,idesc) = "f        > number of arms per star"
        th_file_param(  6,idesc) = "phistar  > volume fraction of star"
        th_file_param(  7,idesc) = "nlin     > number of segments of linear polymer"
        th_file_param(  8,idesc) = "philin   > volume fraction of linera polymer"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_rpastar = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4star  = abs( pa( 2))
      narmeff  = nint(pa( 3))
      diffstar = abs( pa( 4))
      locr2_b  = abs( pa( 5))
      locr2_a  = abs( pa( 6))
      locr2_ta = abs( pa( 7))
      locr2_lz = abs( pa( 8))
      locr2_te = abs( pa( 9))
      nro_me   = nint(pa(10))
      re_me    = abs( pa(11))
      wl4      = abs( pa(12))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: scattering wavevector
      xh = 0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh = 400d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: number of segments per star arm
      xh = 1
      call parget('narm    ',xh,iadda,ier)
      narm     = xh
! >>> extract: segment length
      xh = 3d0
      call parget('l       ',xh,iadda,ier)
      l        = xh
! >>> extract: number of arms per star
      xh = 1
      call parget('f       ',xh,iadda,ier)
      f        = nint(xh)
! >>> extract: volume fraction of star
      xh = 0.01d0
      call parget('phistar ',xh,iadda,ier)
      phistar  = xh
! >>> extract: number of segments of linear polymer
      xh = 10
      call parget('nlin    ',xh,iadda,ier)
      nlin     = nint(xh)
! >>> extract: volume fraction of linera polymer
      xh = 0.9d0
      call parget('philin  ',xh,iadda,ier)
      philin   = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

! linear chain function according to the reptation interpolation model locrep2

     sqt0    =  local_reptation2( q*locr2_a, locr2_te  , locr2_b, locr2_lz)
     sqt     =  local_reptation2( q*locr2_a, t/locr2_ta, locr2_b, locr2_lz)
     locrep2 =  sqt/sqt0
     dr      = 1d-20
     ifix    = 0
     call  NrouseY(q,t,temp,Dr,wl4,nro_me,re_me, Wx, lx,ifix, sqt0,sqt)
     plin0   = nlin  * Debye_qnl(q, nlin, l) 
     plin    = plin0 * locrep2 * sqt / sqt0

 ! start       
              
     Re_arm = sqrt(narm * l**2)                                              
     sqt0   =  pericostar_sqt(q,0d0,f,narmeff,Re_arm,Wl4star,diffstar)    ! automatic diffusion if diffstar==0     
     sqt    =  pericostar_sqt(q,t  ,f,narmeff,Re_arm,Wl4star,diffstar)  
 ! fix normalization
     pstar0 = SQstar(q,narm,f,l)
     pstar  = sqt/sqt0 * pstar0
     pstar  =  pstar * exp( -diffstar * q*q *t )                       !! add diffusion expilicitly

! Rpa (assuming contrast only between h-lin polymer and the rest

!     Sqt0 = Plin0*(phistar+philin-1d0)*(Plin0*philin+Pstar0*phistar)/((Plin0-Pstar0)*phistar-Plin0)
!     Sqt  = Plin *(phistar+philin-1d0)*(Plin *philin+Pstar *phistar)/((Plin -Pstar )*phistar-Plin )


      Sqt0   = philin*((phistar+philin-1)*Plin0-Pstar0*phistar)*Plin0/(Plin0*(phistar-1)-Pstar0*phistar)
      Sqt    = philin*((phistar+philin-1)*Plin -Pstar *phistar)*Plin /(Plin *(phistar-1)-Pstar *phistar)

     th_rpastar = Sqt/Sqt0

 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
      function Debye_qnl( q, nlin, l) result(p0)
        implicit none
        double precision, intent(in)   :: q
        integer         , intent(in)   :: nlin
        double precision, intent(in)   :: l
        double precision               :: p0

        double precision               :: x

        x   = (q*l)**2 * (nlin / 6d0)
        p0  = 2d0/x**2 *((exp(-x)-1d0)+x)

      end function  Debye_qnl

      function SQstar( q, narm, f, l) result(p0)
        implicit none
        double precision, intent(in)   :: q
        integer         , intent(in)   :: narm
        integer         , intent(in)   :: f
        double precision, intent(in)   :: l

        double precision               :: p0

        double precision               :: x

        x   = (q*l)**2 * (narm / 6d0)
        p0  = narm * (2d0/x**2)
        p0  = p0 * ( x - (1d0-exp(-x) ) + (f-1)/2d0 * (1d0-exp(-x))**2) 

      end function  SQstar





 end function th_rpastar
