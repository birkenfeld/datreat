 FUNCTION th_rrofpm(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  rr**2 of Rouse by discrete summation with mode restriction
!  Rouse, Doi_Edwards
      use theory_description 
      implicit none 
      real    :: th_rrofpm
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
     double precision :: Re         ! undisturbed re of chain                                                         
     double precision :: N          ! number of segments                                                              
     double precision :: pmin       ! pmin                                                                            
     double precision :: pwidth     ! width                                                                           
     double precision :: mix        ! mix mix*fermi+(1-mix)*erf                                                       
! the recin parameter representation 
     double precision :: q          ! momentum transfer                                                               
     double precision :: temp       ! temperature                                                                     
! the reout parameter representation 
     double precision :: diff=0d0       ! diffusion                                                                       
 
     double precision :: th
 
     double precision   :: nseg
     double precision   :: amodif(11), a
     character(len=8)   :: parx
     integer            :: i
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'rrofpm'
       nparx =        16
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_rrofpm = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " rr**2 of Rouse by discrete summation with mode restriction"
       th_citation(idesc)     = " Rouse, Doi_Edwards"
!       --------------> set the parameter names --->
        parnam ( 1) = 'Re      '  ! undisturbed re of chain                                                         
        parnam ( 2) = 'N       '  ! number of segments                                                              
        parnam ( 3) = 'pmin    '  ! pmin                                                                            
        parnam ( 4) = 'pwidth  '  ! width                                                                           
        parnam ( 5) = 'mix     '  ! mix mix*fermi+(1-mix)*erf                                                       
        parnam ( 6) = 'amod1   '  ! extra mode modification factor   
        parnam ( 7) = 'amod2   '  ! extra mode modification factor   
        parnam ( 8) = 'amod3   '  ! extra mode modification factor   
        parnam ( 9) = 'amod4   '  ! extra mode modification factor   
        parnam (10) = 'amod5   '  ! extra mode modification factor   
        parnam (11) = 'amod6   '  ! extra mode modification factor   
        parnam (12) = 'amod7   '  ! extra mode modification factor   
        parnam (13) = 'amod8   '  ! extra mode modification factor   
        parnam (14) = 'amod9   '  ! extra mode modification factor   
        parnam (15) = 'amod10  '  ! extra mode modification factor   
        parnam (16) = 'amod11  '  ! extra mode modification factor   
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "undisturbed re of chain" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments" !//cr//parspace//&
        th_param_desc( 3,idesc) = "pmin" !//cr//parspace//&
        th_param_desc( 4,idesc) = "width" !//cr//parspace//&
        th_param_desc( 5,idesc) = "mix mix*fermi+(1-mix)*erf" !//cr//parspace//&
        th_param_desc(6:16,idesc) = "extra amplitude modifiers (on top of fermi/erf) " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "first 11 total mode amplitudes amod*fermierf"
! 
        th_rrofpm = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      Re       =      pa( 1)
      N        =      pa( 2)
      pmin     =      pa( 3)
      pwidth   =      pa( 4)
      mix      =      pa( 5)
      amodif   =      pa(6:16)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: momentum transfer
      xh =  0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh = 300d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     nseg   = x



     th = Re**2 * rr(nint(nseg),nint(N),pmin,pwidth,mix)


     th_rrofpm = th
 
! ---- writing computed parameters to the record >>>  
  call parset('diff    ',sngl(diff),iadda,ier)
     
     do i = 1,11
       a =  modecut(dble(i), pmin, pwidth, mix) * amodif(i)
       write(parx,'("amod",i0)')i
       call parset(parx,sngl(a),iadda,ier)
     enddo


     CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

function rr(nseg,N,pmin,pwidth,mix) result (rx)  ! <r**2> in terms of Re**2
  implicit none
  integer         , intent(in) :: nseg
  integer         , intent(in) :: N
  double precision, intent(in) :: pmin
  double precision, intent(in) :: pwidth
  double precision, intent(in) :: mix
  double precision :: rx

  integer :: p

  rx = 0d0
  do p = 1,N
    rx = rx + cos( (dble(p)*Pi*nseg)/dble(N) )**2 / p**2   * modecut(dble(p), pmin, pwidth, mix) * amodif(min(p,11))
  enddo

  rx = 2d0/(3*Pi**2) * rx

end function rr

function modecut(p,pmin,pwidth,mix) result(fp)
  double precision, intent(in) :: p,pmin,pwidth,mix
  double precision :: fp


  fp = (1d0-mix)*(1d0 +  erf((p-pmin)/pwidth))/2  +  mix/(1d0+exp((pmin-p)/pwidth))

end function modecut



 end function th_rrofpm
