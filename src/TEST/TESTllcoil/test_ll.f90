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
                             
 subroutine test_llcoil
!================================================================================
!  compute the frequency dependent impedance of the sc-maincoi
      use lcalc_Zomega_cc
    
! the internal parameter representation 
     double precision :: dlengtha                                     
     double precision :: wwidtha                                      
     double precision :: rsepa                                        
     double precision :: cper                                        
     double precision :: cpa                                         
     double precision :: rhoal                                       
     double precision :: r0                                          
! the recin parameter representation
     double precision :: abspha                                      
! the reout parameter representation 
     double precision :: intens0    ! intensity factor (cm**-1)
     double precision :: Rwc_mod(10)
 
     double precision   :: frequenz
     complex            :: zx
!
! ---- transfer parameters -----
  
      dlengtha  =     0.01d0
      wwidtha   =     0.01d0
      rsepa     =     0.0005d0
      cper      =     1d-9
      cpa       =     1d-9
      rhoal     =     1d-5
      r0        =     1
      Rwc_mod   =     1
 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     x = 0.1d0
     frequenz   = x

write(*,*) "freq=",x
     zx         =  Z_llcoil(frequenz, rhoAl, cpa, cper, dlengtha, wwidtha, rsepa, Rwc_mod) + R0
write(*,*) "zx=",zx

 end subroutine test_llcoil
