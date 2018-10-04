       program pytest
         implicit none
         double precision :: q, rr, den, eps
         double precision :: peryev

         integer :: i
         
         rr   = 1
         den  = 0.01
         eps  = 1d-7

        do i=1,200
         q = i * 0.02
         write(6,*)i,q,peryev(q,rr,den,eps)
        enddo


       end program pytest

 
c  percus-yevick s(q)
c -------------------
c lit. m.s. wertheim, prl 10 (1963) 321
c -------------------------------------
c
       function peryev(q,rr,den,eps)
c      -----------------------------
c      q     = impulsuebertrag
c      rr    = teilchenradius
c      den   = teilchendichte
c      eps   = realteil ---> 0
c
       parameter(pi=3.141592654)
       implicit complex*16 (c)
       complex*16           dcmplx, cdexp
       real*8   eta, dimag, deps, dr, dq, q, r, rr, den, eps
       real*8   peryev
c
       r   = 2*rr
c --- r ist der durchmesser ! der teilchen !
       dr  = r
       dq  = q
       deps= eps
       eta = pi*den*r*r*r/6.d0
c
       ct  = dcmplx ( deps , -dq*dr )
       cst = (((1.d0-eta)**2*ct+6.d0*eta*(1.d0-eta))*ct
     *      +18.d0*eta*eta)*ct
     *      -12.d0*eta*(1.d0+2.d0*eta)
       clt = 12.d0*eta*((1.d0+0.5d0*eta)*ct+(1.d0+2.d0*eta))
c
       cgt = (ct*clt) / ( 12.d0*eta*(clt+cst*cdexp(ct)) )
c
       peryev =  1.d0 + den * (4.d0*pi/dq) * dr*dr * dimag( cgt )
c
       return
       end
