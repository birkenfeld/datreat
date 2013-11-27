
                                                                        
      SUBROUTINE uni_fourier(w,s,ds,n,t,ft1,ft2,ft,dft) 
!     -------------------------------------------------
!
! Fourier integration over histogram data s(i) vs w(i)
! including error determinateion if data error is in ds(i)
! 
! spacing may be non-unigorm
!
                                                                       
      implicit none

      integer,          intent(in)    :: n                   ! length of data vectors
      double precision, intent(in)    :: w(n), s(n), ds(n)   ! omega, data, data-error
      double precision, intent(in)    :: t                   ! time parameter (result)
      double precision, intent(out)   :: ft1,ft2,ft,dft      ! real, imaginary, modulus, err(modulus) FT-values
 
      double precision                :: a1(n),a2(n),swt(n),cwt(n)
      double precision                :: t2, dft2  
      integer                         :: i
                                                                       
                                                                        
      if (t.ne.0.0) then 
                                                                        
        do i=1,n 
          swt(i)=sin(w(i)*t) 
          cwt(i)=cos(w(i)*t) 
        end do 
        t2=1.0d0/(t*t) 
                                                                        
        a1(1)=-swt(1)/t+(cwt(1)-cwt(2))/(w(2)-w(1))*t2 
        a2(1)=cwt(1)/t+(swt(1)-swt(2))/(w(2)-w(1))*t2 
        do i=2,n-1 
          a1(i)=((cwt(i)-cwt(i+1))/(w(i+1)-w(i))                  &
     &          +(cwt(i-1)-cwt(i))/(w(i-1)-w(i)))*t2              
          a2(i)=((swt(i)-swt(i+1))/(w(i+1)-w(i))                  &
     &          +(swt(i-1)-swt(i))/(w(i-1)-w(i)))*t2              

        end do 
        a1(n) = swt(n)/t+(cwt(n-1)-cwt(n))/(w(n-1)-w(n))*t2 
        a2(n) =-cwt(n)/t+(swt(n-1)-swt(n))/(w(n-1)-w(n))*t2 
                                                                        
      else 
                                                                        
        a1(1)=(w(2)-w(1))*0.5d0 
        a2(1)=0.0d0 
        do i=2,n 
          a1(i)=(w(i+1)-w(i-1))*0.5d0 
          a2(i)=0.0d0 
        end do 
        a1(n)=(w(n)-w(n-1))*0.5d0 
        a2(n)=0.0d0 
                                                                        
      end if 
                                                                        
      ft1=0.0d0 
      ft2=0.0d0 
      do i=1,n 
        ft1=ft1+a1(i)*s(i) 
        ft2=ft2+a2(i)*s(i) 
      end do 
      ft=sqrt(ft1*ft1+ft2*ft2) 
                                                                        
      dft2=0.0d0 
      do i=1,n 
        dft2=dft2+((a1(i)*ft1+a2(i)*ft2)/ft*ds(i))**2 
      end do 
      dft=sqrt(dft2) 
                                                                        
      return
 
      END  SUBROUTINE uni_fourier                                        
                                                                        
