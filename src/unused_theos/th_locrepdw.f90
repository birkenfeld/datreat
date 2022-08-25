 FUNCTION th_locrepdw(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, but with finite summation of integrals and lenthscale, timescale and fluctuation ratio as parameters DR Version Nov2020
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrepdw
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
     double precision :: b          ! fluctuation intensity (relative)                                                
     double precision :: a          ! length scale                                                                    
     double precision :: W          ! (Rouse)rate                                                                     
     double precision :: n          ! No-segments                                                                     
     double precision :: ne         ! No-segments/entanglement                                                        
! the recin parameter representation 
     double precision :: q          ! q-value    default value
     double precision :: t0         ! w cross over time
     double precision :: tw         ! w transition width
     double precision :: ww0        ! transition depth
! the reout parameter representation 
 
     double precision :: th
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locrepdw'
       nparx =        9
       IF (npar.lt.nparx) then
           WRITE (9,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrepdw= 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, but with finite summation of integrals and lenthscale, timescale and fluctuation ratio as parameters DR Version Nov2020 with an additional Fermi cross over in w"
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'b       '  ! fluctuation intensity (relative)                                                
        parnam ( 3) = 'a       '  ! length scale                                                                    
        parnam ( 4) = 'W       '  ! (Rouse)rate                                                                     
        parnam ( 5) = 'n       '  ! No-segments                                                                     
        parnam ( 6) = 'ne      '  ! No-segments/entanglement
        parnam ( 7) = 't0      '  ! w cross over time t0
        parnam ( 8) = 'tw      '  ! width of w cross over
        parnam ( 9) = 'ww0     '  ! transition depth
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fluctuation intensity (relative)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "(Rouse)rate" !//cr//parspace//&
        th_param_desc( 5,idesc) = "No-segments" !//cr//parspace//&
        th_param_desc( 6,idesc) = "No-segments/entanglement" !//cr//parspace//&
        th_param_desc( 7,idesc) = "cross over time of w" !//cr//parspace//&
        th_param_desc( 8,idesc) = "cross over width of w" !//cr//parspace//&
        th_param_desc( 9,idesc) = "transition depth ww0" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_param_desc( 1,idesc) = "asqrtne     >  a*sqrt(ne)" !//cr//parspace//& 
        th_locrepdw = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      b        =      pa( 2)
      a        =  abs(pa( 3) )
      W        =  abs(pa( 4) )
      n        =  abs(pa( 5) )
      ne       =  abs(pa( 6) )
      t0       =  abs(pa( 7) )
      tw       =  abs(pa( 8) )
      ww0      =  abs(pa( 9) )
      
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t  = x              ! since we prefer to call the independent variable t, x must be copied to t
     th = ampli * local_reptationdr(q, t, a, W, n, ne, b,t0,tw,ww0)

     th_locrepdw = th

     call        parset('asqrtne   ',sngl(a*sqrt(ne)),iadda,ier)      !     "
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function local_reptationdr(q, t, a, W, n, ne, b,t0,tw,ww0) result(val)
    implicit none
    double precision, intent(in)   :: q, t
    double precision, intent(in)   :: a        ! step (entanglement?) length
    double precision, intent(in)   :: W        ! (Rouse?) rate
    double precision, intent(in)   :: n        ! No Segments
    double precision, intent(in)   :: ne       ! No Segments/entanglement
    double precision, intent(in)   :: b        ! b-fluctuation intensity
    double precision, intent(in)   :: t0       ! w cross over time
    double precision, intent(in)   :: tw       ! w cross-over width
    double precision, intent(in)   :: ww0      ! transition depth
    double precision               :: wt
    double precision               :: val
    double precision, parameter    :: Pi = 4*atan(1d0)

    double precision :: T1, T2, T3
    wt = w*(1 - (1-ww0)/(1+exp((t-t0)/tw)))

    T1 = 72d0 * (exp(-q**2 * n * a**2 / 6d0) + q**2 * n * a**2 / 6d0 - 1d0) / q**4 / a**4
    T2 = b*((exp(-1/6*a**2*n*q**2-n**2/(4*wt*t))-1)*2/3*ne*sqrt(wt*t/Pi)+ &
        ne*(n/3+a**2*q**2*wt*t/9)*exp(a**4*q**4*wt*t/36)* &
         (erfc(a**2*q**2*sqrt(wt*t)/6)-erfc(n/(2*sqrt(wt*t))+a**2*q**2*sqrt(wt*t)/6)))
    T3 = 72d0/((q**4*a**4))*(exp(-(q**2*n*a**2)/6d0)+(q**2*n*a**2)/6d0-1d0)+b*ne*n/3d0

    val = (T1+T2)/T3


  end function local_reptationdr
 end function th_locrepdw
