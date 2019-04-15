 FUNCTION th_sansfac(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  yields the intensity prefactor for polymer sans spectra. Polymer data: from the parameter section of data rec. Use as multiply theory in combination with an S(q) type form factor...
! 
      use theory_description 
      implicit none 
      real    :: th_sansfac
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
! the recin parameter representation 
     double precision :: Mw         ! molecular weight of objects  (g/mol)                                            
     double precision :: rho        ! scattering length density polymer (cm**-1)                                      
     double precision :: rhosolv    ! scattering length density background(solvent) (cm**-1)                          
     double precision :: volfrac    ! volume fraction of polymer                                                      
     double precision :: density    ! density of polymer (g/cm*+3)                                                    
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)                                                       
 
     double precision   :: Navogadro = 6.022045d23
     double precision   :: q
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sansfac'
       nparx =        1
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sansfac = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " yields the intensity prefactor for polymer sans spectra. "//cr//parspace//&
                                " Polymer data: from the parameter section of data rec. "//cr//parspace//&
                                " Use as multiply theory in combination with an S(q) type form factor..."
       th_citation(idesc)     = " "
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "Mw       > molecular weight of objects  (g/mol)"
        th_file_param(  2,idesc) = "rho      > scattering length density polymer (cm**-1)"
        th_file_param(  3,idesc) = "rhosolv  > scattering length density background(solvent) (cm**-1)"
        th_file_param(  4,idesc) = "volfrac  > volume fraction of polymer"
        th_file_param(  5,idesc) = "density  > density of polymer (g/cm*+3)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(  1,idesc) = "intens0  > intensity factor (cm**-1)"
! 
        th_sansfac = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: molecular weight of objects  (g/mol)
      xh = 1
      call parget('Mw      ',xh,iadda,ier)
      Mw       = xh
! >>> extract: scattering length density polymer (cm**-1)
      xh = 1
      call parget('rho     ',xh,iadda,ier)
      rho      = xh
! >>> extract: scattering length density background(solvent) (cm**-1)
      xh = 1
      call parget('rhosolv ',xh,iadda,ier)
      rhosolv  = xh
! >>> extract: volume fraction of polymer
      xh = 0.5
      call parget('volfrac ',xh,iadda,ier)
      volfrac  = xh
! >>> extract: density of polymer (g/cm*+3)
      xh = 1
      call parget('density ',xh,iadda,ier)
      density  = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x

     intens0  =   (rho-rhosolv)**2           &
                * volfrac*(1d0-volfrac)      &
                * Mw / ( density * Navogadro )



     th_sansfac = ampli * intens0
 
! ---- writing computed parameters to the record >>>  
      call parset('intens0 ',sngl(intens0),iadda,ier)
 
! CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

 end function th_sansfac
