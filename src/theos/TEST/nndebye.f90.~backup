program nndebye_test
implicit none
double precision :: qq=0.1d0, ll=3d0, nu=0.3d0, y
integer          :: nn=1000, i
double precision :: ql2, mu, Rg, x, deb

   mu  = 2d0 * nu

   


do i=0,100
  qq = 0.001d0*i
   ql2 = (1d0/6d0) * (qq*ll)**2
   Rg  = sqrt((ll*ll*dble(nn)**mu) / 6d0)   
  x  = (qq*Rg)**2
  deb= (2d0/x**2)*(exp(-x)-1+x)
  y  = nndebye(qq, ll, nu,nn )
  write(*,*) i, qq, y, deb 
enddo


CONTAINS


function nndebye(q, l, nue,n ) result(val) 
implicit none
double precision ,intent(in) :: q           ! Q-value
double precision ,intent(in) :: l           ! effective segment length
double precision ,intent(in) :: nue         ! exponent
integer          ,intent(in) :: n           ! number of segments
double precision             :: val

integer          :: i, j
double precision :: eterms(0:n-1)
double precision :: ql2, mu, Rg

   ql2 = (1d0/6d0) * (q*l)**2
   mu  = 2d0 * nue
!  Rg  = sqrt((l*l*dble(n)**mu) / 6d0) )
!$OMP PARALLEL DO    
   do i=0,n-1
     eterms(i) = exp(-ql2*dble(i)**mu)
   enddo
!$OMP END PARALLEL DO 
   val = 0
!$OMP PARALLEL DO REDUCTION(+:val)
   do i=1,n
    do j=1,n
      val = val + eterms(abs(i-j))
    enddo
   enddo
!$OMP END PARALLEL DO 

   val = val / (n*n)
 
end function nndebye


end program nndebye_test


