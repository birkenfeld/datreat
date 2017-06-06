!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! simple scattering function core-shell                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
! Basic F(Q) Factor for a sphere                                        
!                                                                       
       double precision function FQ_sphere(Q,R) 
!      ----------------------------------------                         
       implicit none 
                                 ! Q-value                              
       double precision Q 
                                 ! Radius fo the sphere                 
       double precision R 
                                                                        
       double precision qr 
       double precision fourpi 
       parameter (fourpi=12.56637061d0) 
                                                                        
       if(Q.eq.0.0) then 
         FQ_sphere = (fourpi/3)*(R**3) 
       else 
         FQ_sphere = fourpi*(Q*R*cos(Q*R)-sin(Q*R))/Q**3 
       endif 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function SQ_cs(Q,P) 
!      ------------------------------------                             
! core shell scattering factor as function of aggregation number P      
! rest of parameters via common                                         
       implicit none 
                                                                        
                             ! Q-value                                  
       double precision Q 
                             ! Aggregation number                       
       integer          P 
                                                                        
! parameters of the micelles via common block                           
                                                    ! Molecularweigts kg
       double precision Mw_core    ,  Mw_shell 
                                                    ! Densities       kg
       double precision den_coremat,  den_shellmat 
                                                    ! Volumefractions   
       double precision phi_core,     phi_shell 
                                                    ! Scatteringlength d
       double precision sigma_core,   sigma_shell 
                                                    ! interface width / 
       double precision smear_core,   smear_shell 
                                                                        
                                                    ! Molecularweigts kg
       common/csatmod/  Mw_core    ,  Mw_shell,                         &
     &                  den_coremat,  den_shellmat,                     &
     &                  phi_core,     phi_shell,                        &
     &                  sigma_core,   sigma_shell,                      &
     &                  smear_core,   smear_shell                       
                                                    ! Densities       kg
                                                    ! Volumefractions   
                                                    ! Scatteringlength d
                                                    ! interface width / 
                                                                        
                                                                        
       double precision Vc, Vs, Rc, Rs, F, Fc, Fs 
       double precision FQ_sphere 
                                                                        
       double precision fourpi, Na 
       parameter (fourpi=12.56637061d0, Na=6.022045d23) 
                                                                        
                                                                        
       Rc=(3/fourpi*P*Mw_core/(Na*den_coremat*phi_core))**(1d0/3d0) 
       Rs=(                                                             &
     &       (3/fourpi*P*Mw_shell/(Na*den_shellmat*phi_shell))-Rc**3    &
     &    )**(1d0/3d0)                                                  
                                                                        
       F =   (sigma_core *phi_core *FQ_sphere(Q,Rc) +                   &
     &        sigma_shell*phi_shell*(FQ_sphere(Q,Rs)-FQ_sphere(Q,Rc)))  
                                                                        
! smear interfaces                                                      
       F = F * exp(-(Q*smear_core)**2) 
! other smearing for the shell                                          
       F = F+sigma_shell*phi_shell*(                                    &
     &     FQ_sphere(Q,Rs)*                                             &
     &     (exp(-(Q*smear_shell)**2)-exp(-(Q*smear_core)**2)))          
                                                                        
                                                                        
       SQ_cs = Na * F*F 
                                                                        
!       write(6,'(a,i3,5e13.6)')'sq ',p,q,Rc,Rs,F,SQ_cs                 
                                                                        
       return 
      END                                           
