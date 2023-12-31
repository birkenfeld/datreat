 FUNCTION th_hyperbra(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Simple approach (Dolgushev) to the dynamics of hyperbranched denrimers: exp(-q**2 * D * t) * exp(-(gam_q*q**q_exp)*t)**beta)
! 
      use theory_description 
      implicit none 
      real    :: th_hyperbra
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
     double precision :: diff       ! center of mass diffusion                                                        
     double precision :: gam_q      ! rate prefactor                                                                  
     double precision :: q_exp      ! Gamma q exponent                                                                
     double precision :: beta       ! streched exp                                                                    
! the recin parameter representation 
     double precision :: q          ! q                                                                               
! the reout parameter representation 
 
     double precision :: th
 
     double precision :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'hyperbra'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_hyperbra = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Simple approach (Dolgushev) to the dynamics of hyperbranched denrimers: exp(-q**2 * D * t) * exp(-(gam_q*q**q_exp)*t)**beta)"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'diff    '  ! center of mass diffusion                                                        
        parnam ( 3) = 'gam_q   '  ! rate prefactor                                                                  
        parnam ( 4) = 'q_exp   '  ! Gamma q exponent                                                                
        parnam ( 5) = 'beta    '  ! streched exp                                                                    
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "center of mass diffusion" !//cr//parspace//&
        th_param_desc( 3,idesc) = "rate prefactor" !//cr//parspace//&
        th_param_desc( 4,idesc) = "Gamma q exponent" !//cr//parspace//&
        th_param_desc( 5,idesc) = "streched exp" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_hyperbra = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      diff     =  abs(pa( 2))
      gam_q    =  abs(pa( 3))
      q_exp    =      pa( 4)
      beta     =      pa( 5)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q
      xh = 0.0d0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     th  = ampli * exp( -q*q*diff*t -((gam_q*q**q_exp)*t)**beta)
     th_hyperbra = th
 
! ---- writing computed parameters to the record >>>  

 end function th_hyperbra
