       function th_hohdiff(x,pa,thnam,parnam,npar, ini, nopar ,params,napar,mbuf) 
!=================================================================
!            Water dynamics for TOF vibrational+rotational+translational motion 
!            RANDOM JUMP - DIFFUSION SEARS
!		 J. Teixeira
!		Experimental determination of the nature of diffusive motions of water molecules at low temperatures
 !      Phys. Rev. A 31, 1913–1917 (1985)
!
implicit none
real th_hohdiff
character*8  thnam,parnam(20)    ! theory and parameternames
real pa(20)
real x
integer      npar,ini,ier
integer      mbuf
integer      nopar               ! Anzahl der Parameter data     , intent(inout) ::
character*80 napar(mbuf)         ! name des parameters n   , intent(inout) ::
real         params(mbuf)        ! value des parameters n  , intent(inout) ::
real       result
integer      nparx
real  amplitu
real  u_quadra
real  self_dif
real  tau0_dif
real  a_rot
real  d_rot
real  shiftt
real  elas_dw
real  linbgr
real  q
logical found_q
double precision pi,h,w,elastic,Dw,gamma,Transdif,j0_a_q,j1_a_q,Rotdiff,I_q_omega

real e_mev  ! meV in data ! internal val. name of independent varable

if(ini.eq.0) then     ! initialization of theories
  thnam = 'hohdiff'        ! name in datreat
  nparx = 7                 ! number of parameters
  if(npar.lt.nparx) then     
    write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
    th_hohdiff = 0
    return
  endif
  npar = nparx
  parnam(1) = 'amplitu'      ! name of theory parameters
  parnam(2) = 'u_quadra'
  parnam(3) = 'self_dif'
  parnam(4) = 'tau0_dif'
  parnam(5) = 'a_rot'
  parnam(6) = 'd_rot'
  parnam(7) = 'shiftt'
  th_hohdiff = 0
  return
endif
!                     Please give the Units and meaning here
e_mev    = x     ! meV      X
amplitu  = pa(1) !          overall amplitude
u_quadra = pa(2) ! A^2      Debye-Waller-Factor <u^2>
self_dif = pa(3) ! A^2/ns   Self Diffusion coefficient Random Jump diffusion
tau0_dif = pa(4) ! ns       residence time             Random Jump diffusion
a_rot    = pa(5) ! A        rotation radius (distance OH in h2O   = 0.98A)
d_rot    = pa(6) ! ns^-1    water rotational diffusion
shiftt   = pa(7) ! as pa(1) shift  elatisc bgr to 0    e_mev+shiftt*q
 q = 0.1         ! A^-1     Q wavevector
 !       get parameters from the data
 call  getpar ('q       ',q,nopar,params,napar,mbuf,ier)

 pi=3.1415926535897                             
 !h=4.13567E-15  eVs= 4.13567E-6  eVns  
 h = 4.13567E-6  
 ! convert meV to ns^-1
 w=    ((e_mev+shiftt*q)*0.001)/(h/2*pi)
 if (w.eq.0) then
 		w=1e-20
 endif
 ! Debye Waller factor
 DW=exp(-1*u_quadra*q**2)
!Translation Diffusion 
 gamma=self_dif*q**2/(1+self_dif*q**2*tau0_dif)
 Transdif=(1/pi)*(gamma)/(gamma**2+w**2)
! Rotational fiffusion
! sperical Bessel j0^2
 j0_a_q=(sin(q*a_rot)/(q*a_rot))**2    
 ! sperical Bessel j0^2
 j1_a_q=( sin(q*a_rot)/(q*a_rot)**2 - cos(q*a_rot)/(q*a_rot) )**2 

 Rotdiff=j0_a_q + (1/pi)*(2*1+1)*j1_a_q*(2*d_rot)/(4*d_rot**2+w**2)
 
 
 I_q_omega=DW*Transdif*Rotdiff
 result = I_q_omega*amplitu
!      ------
 th_hohdiff     =sngl(result)
! write some parameters calculated to the data 
       call  setpar ('u_quadra', u_quadra, nopar ,params,napar,mbuf, ier)
       call  setpar ('self_dif', self_dif, nopar ,params,napar,mbuf, ier)
       call  setpar ('tau0_dif', tau0_dif, nopar ,params,napar,mbuf, ier)
       call  setpar ('a_rot   ', a_rot   , nopar ,params,napar,mbuf, ier)
       call  setpar ('d_rot   ', d_rot   , nopar ,params,napar,mbuf, ier)
       call  setpar ('shiftt  ', shiftt  , nopar ,params,napar,mbuf, ier)
   return
 end 





