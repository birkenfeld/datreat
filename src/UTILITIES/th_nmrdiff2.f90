 FUNCTION th_nmrdiff2(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  for for pfg nmr diffsuion curves Zeitbereichdiffusion (Kärger), two-state model
!  J. Kaerger, Annalen der Physik 7 Band 27(1971) 107-109
      use theory_description 
      implicit none 
      real    :: th_nmrdiff2
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
     double precision :: diff1      ! diffusion constant 1                                                            
     double precision :: diff2      ! diffusion constant 2                                                            
     double precision :: tau1       ! lifetime 1                                                                      
     double precision :: tau2       ! lifetime 2                                                                      
! the recin parameter representation 
     double precision :: delta      ! delta time of signal detection                                                  
! the reout parameter representation 
 
     double precision   :: qq
     double precision   :: a1, a2, p1p, p2p, p1, p2
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nmrdiff2'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nmrdiff2 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " for for pfg nmr diffsuion curves Zeitbereichdiffusion (Kärger), two-state model"
       th_citation(idesc)     = " J. Kaerger, Annalen der Physik 7 Band 27(1971) 107-109"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'diff1   '  ! diffusion constant 1                                                            
        parnam ( 3) = 'diff2   '  ! diffusion constant 2                                                            
        parnam ( 4) = 'tau1    '  ! lifetime 1                                                                      
        parnam ( 5) = 'tau2    '  ! lifetime 2                                                                      
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "diffusion constant 1" !//cr//parspace//&
        th_param_desc( 3,idesc) = "diffusion constant 2" !//cr//parspace//&
        th_param_desc( 4,idesc) = "lifetime 1" !//cr//parspace//&
        th_param_desc( 5,idesc) = "lifetime 2" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "delta    > delta time of signal detection"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_nmrdiff2 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      diff1    =      pa( 2)
      diff2    =      pa( 3)
      tau1     =      pa( 4)
      tau2     =      pa( 5)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: delta time of signal detection
      xh = 
      call parget('delta   ',xh,iadda,ier)
      delta    = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     qq  = x

     a1  = 0.5d0*(qq*(diff1+diff2) + 1d0/tau1 + 1d0/tau2)    &
         -  sqrt( (qq*(diff2-diff1) + 1d0/tau2 - 1d0/tau1)**2 + 4d0/(tau1*tau2) )

     a2  = 0.5d0*(qq*(diff1+diff2) + 1d0/tau1 + 1d0/tau2)    &
         +  sqrt( (qq*(diff2-diff1) + 1d0/tau2 - 1d0/tau1)**2 + 4d0/(tau1*tau2) )

     p1  = tau1/(tau1+tau2)
     p2  = tau2/(tau1+tau2)
     p2p = 1d0/(a1-a2) *(p1 * qq * diff1 + p2 * qq * diff2 - a1)
     p1p = 1-p2p

     th_nmrdiff2 = ampli * ( p1p*exp(-a1*delta) + p2p*exp(-a2*delta) )

 
 end function th_nmrdiff2
