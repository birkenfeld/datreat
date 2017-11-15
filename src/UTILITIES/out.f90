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
     double precision :: wl4        ! rouse rate                                                                      
     double precision :: narmeff    ! effective segment star arm                                                      
     double precision :: diffstar   ! centre-of-mass diffusino star                                                   
     double precision :: locr2_b    ! locrep2 parameter b                                                             
     double precision :: locr2_a    ! locrep2 parameter a                                                             
     double precision :: locr2_ta   ! locrep2 parameter tau                                                           
     double precision :: locr2_lz   ! locrep2 parameter lz                                                            
     double precision :: locr2_te   ! locrep2 parameter teps                                                          
     double precision :: nro_me     ! effective segments for Me-Modelling                                             
     double precision :: re_me      ! Re of Me (entangle strand Re)                                                   
! the recin parameter representation 
     double precision :: q          ! scattering wavevector                                                           
     double precision :: temp       ! temperature                                                                     
     double precision :: narm       ! number of segments per star arm                                                 
     double precision :: l          ! segment length                                                                  
     double precision :: f          ! number of arms per star                                                         
     double precision :: phistar    ! volume fraction of star                                                         
     double precision :: nlin       ! number of segments of linear polymer                                            
     double precision :: philin     ! volume fraction of linera polymer                                               
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)                                                       
 
     double precision   :: Navogadro = 6.022045d23
     double precision   :: q
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'rpastar'
       nparx =       11
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
        parnam ( 2) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 3) = 'narmeff '  ! effective segment star arm                                                      
        parnam ( 4) = 'diffstar'  ! centre-of-mass diffusino star                                                   
        parnam ( 5) = 'locr2_b '  ! locrep2 parameter b                                                             
        parnam ( 6) = 'locr2_a '  ! locrep2 parameter a                                                             
        parnam ( 7) = 'locr2_ta'  ! locrep2 parameter tau                                                           
        parnam ( 8) = 'locr2_lz'  ! locrep2 parameter lz                                                            
        parnam ( 9) = 'locr2_te'  ! locrep2 parameter teps                                                          
        parnam (10) = 'nro_me  '  ! effective segments for Me-Modelling                                             
        parnam (11) = 're_me   '  ! Re of Me (entangle strand Re)                                                   
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment star arm" !//cr//parspace//&
        th_param_desc( 4,idesc) = "centre-of-mass diffusino star" !//cr//parspace//&
        th_param_desc( 5,idesc) = "locrep2 parameter b" !//cr//parspace//&
        th_param_desc( 6,idesc) = "locrep2 parameter a" !//cr//parspace//&
        th_param_desc( 7,idesc) = "locrep2 parameter tau" !//cr//parspace//&
        th_param_desc( 8,idesc) = "locrep2 parameter lz" !//cr//parspace//&
        th_param_desc( 9,idesc) = "locrep2 parameter teps" !//cr//parspace//&
        th_param_desc(10,idesc) = "effective segments for Me-Modelling" !//cr//parspace//&
        th_param_desc(11,idesc) = "Re of Me (entangle strand Re)" !//cr//parspace//&
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
        th_out_param(  1,idesc) = "intens0  > intensity factor (cm**-1)"
! 
        th_rpastar = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4      =      pa( 2)
      narmeff  =      pa( 3)
      diffstar =      pa( 4)
      locr2_b  =      pa( 5)
      locr2_a  =      pa( 6)
      locr2_ta =      pa( 7)
      locr2_lz =      pa( 8)
      locr2_te =      pa( 9)
      nro_me   =      pa(10)
      re_me    =      pa(11)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: scattering wavevector
      xh = 
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh = 
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: number of segments per star arm
      xh = 
      call parget('narm    ',xh,iadda,ier)
      narm     = xh
! >>> extract: segment length
      xh = 
      call parget('l       ',xh,iadda,ier)
      l        = xh
! >>> extract: number of arms per star
      xh = 
      call parget('f       ',xh,iadda,ier)
      f        = xh
! >>> extract: volume fraction of star
      xh = 
      call parget('phistar ',xh,iadda,ier)
      phistar  = xh
! >>> extract: number of segments of linear polymer
      xh = 
      call parget('nlin    ',xh,iadda,ier)
      nlin     = xh
! >>> extract: volume fraction of linera polymer
      xh = 
      call parget('philin  ',xh,iadda,ier)
      philin   = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

     rpastar  =



     th_rpastar = ! INSERT RESULTING VALUE OF TH EVALUATION HERE
 
! ---- writing computed parameters to the record >>>  
      call parset('intens0 ',sngl(intens0),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

 end function th_rpastar
