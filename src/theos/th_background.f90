!*ds                                                                    
!*ed                                                                    
! *********** thdatr1 fortran ******************************************
!                                                                      *
! ---> version for general data fitting !!!!                            
!                                                                      *
!***********************************************************************
!*                                                                     *
!*    provides theories for triax program system                      * 
!*    please use the convention explained by the comments in the       *
!*    sample theories                                                  *
!*                                                                     *
!*    s(q,ome) is expected to be in units of sec                       *
!*             the scaling factor may be interpreted as                *
!*             n*sigma-inc in cm**2                                    *
!*             --> omit the 1/hquer factor of the usual def. of s(q,o) *
!*                                                                     *
!*    on input q (via qq(3) ) is given in miller-indizes               *
!*    om       is given in thz      (neg. = energy gain of neutron)    *
!*                                                                     *
!*                                                                     *
!***********************************************************************
!                                                                       
!*ds                                                                    
                                                                   
      FUNCTION bkgr (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!  bkgr 
!  linear background bkgr = pa (1) + pa (2) * x   ---------------
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character(len=80), intent(inout) :: napar(mbuf) ! name des parameters n
      real, intent(inout) :: params(mbuf)       ! value des parameters n
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
      thnam = 'bkgr    ' 
         nparx = 2 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th2 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'level ' 
         parnam (2) = 'slope ' 
!                                                                       
         th2 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
!                                                                       
      bkgr = pa (1) + pa (2) * x 
!      ------------------------                                         
!                                                                       
      RETURN 
      END FUNCTION bkgr                         
!                                                                       
