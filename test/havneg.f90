program havneg
implicit none
complex(kind=8) :: y, z
real(kind=8)    :: om, alpha, beta, tau


 om = 1
 tau = 1
 
 do while (tau > 0d0) 
   write(*,*)" tau, alpha, beta "
   read(5,*)tau, alpha, beta
   
   y = (1d0,0d0) / ( (1d0,0d0) + ((0d0,1d0)*om*tau)**alpha)**beta

   write(*,*) y

 enddo

 
end program havneg 
