program dlocrep_Snm
implicit none

double precision :: b=1d0
double precision :: wt
integer          :: nm, i,j
double precision :: Snm, Snm_yy, Snm_s2, Snm_s20

write(*,'(//a//)')"==1==========================================================================="

do i=0,50
 wt = 0.1d0*i
 write(*,'(i3,f12.6,2x,11f12.6)')i,wt,(Snm(wt,b,nm),nm=0,10)
enddo

write(*,'(//a//)')"==2==========================================================================="

do i=0,50
 wt = 0.1d0*i
 write(*,'(i3,f12.6,2x,11f12.6)')i,wt,(Snm_yy(wt,b,nm),nm=0,10)
enddo


write(*,'(//a//)')"==3==========================================================================="

do i=0,50
 wt = 0.1d0*i
 write(*,'(i3,f12.6,2x,11f12.6)')i,wt,(Snm_s2(wt,b,nm),nm=0,10)
enddo

write(*,'(//a//)')"==2-1========================================================================="

do i=0,50
 wt = 0.1d0*i
 write(*,'(i3,f12.6,2x,11f12.6)')i,wt,(Snm_yy(wt,b,nm)-Snm(wt,b,nm),nm=0,10)
enddo


write(*,'(//a//)')"==1-1a========================================================================="

do i=0,50
 wt = 0.1d0*i
 write(*,'(i3,f12.6,2x,11f12.6)')i,wt,(Snm(wt,b,nm)-Snm_s20(wt,b,nm),nm=0,10)
enddo



end program dlocrep_Snm


function Snm(wt0, b, nm) result(t49)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|


! y:=1+beta * exp(-(nm+u-v)**2/(4*wt));
! Snm:=int(y,u=0..1,v=0..1)
! with
! beta=b/sqrt(36*Pi*wt)

! Maple expression:
! Snm:=
! sqrt(wt)*b*sqrt(exp(nm/wt))/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4)*(exp(1/wt))^(1/4))+
! sqrt(wt)*b/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4)*sqrt(exp(nm/wt))*(exp(1/wt))^(1/4))-
! 2*sqrt(wt)*b/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4))-(1/6)*b*erf(-1/(2*sqrt(wt))+
! nm/(2*sqrt(wt)))+(1/6)*b*erf(-1/(2*sqrt(wt))+nm/(2*sqrt(wt)))*nm+
! (1/6)*b*erf(1/(2*sqrt(wt))+nm/(2*sqrt(wt)))+(1/6)*b*erf(1/(2*sqrt(wt))+
! nm/(2*sqrt(wt)))*nm-(1/3)*erf(nm/(2*sqrt(wt)))*b*nm+1
! 


 double precision :: t49
 real(kind=16)    :: wt, t1,t2,t4,t5,t6,t7,t9,t10,t11,t13,t14,t16,t17,t18,t30,t31,t33,t34,t39,t40,t45

      if(wt0 < 1d-8) then
        wt = 1d-8
      else
        wt = wt0
      endif 

      t1 = sqrt(wt)
      t2 = sqrt(0.3141592654D1)
      t4 = 0.1D1 / t2 * t1
      t5 = b * t4
      t6 = 0.1D1 / wt
      t7 = nm ** 2
      t9 = exp(t7 * t6)
      t10 = t9 ** (0.1D1 / 0.4D1)
      t11 = 0.1D1 / t10
      t13 = exp(nm * t6)
      t14 = sqrt(t13)
      t16 = exp(t6)
      t17 = t16 ** (0.1D1 / 0.4D1)
      t18 = 0.1D1 / t17
      t30 = 0.1D1 / t1
      t31 = nm * t30
      t33 = erf(-t30 / 0.2D1 + t31 / 0.2D1)
      t34 = t33 * b
      t39 = erf(t30 / 0.2D1 + t31 / 0.2D1)
      t40 = t39 * b
      t45 = erf(t31 / 0.2D1)
      t49 = t18 * t14 * t11 * t5 / 0.3D1 + t18 / t14 * t11 * t5 / 0.3D1 &
          - 0.2D1 / 0.3D1 * t11 * b * t4 - t34 / 0.6D1 + nm * t34 / 0.6D1 + &
            t40 / 0.6D1 + nm * t40 / 0.6D1 - nm * b * t45 / 0.3D1 + 0.1D1
      if(isnan(t49)) t49=1

end function Snm


function Snm_yy(wt0, b, nm) result(yyy)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: yyy
 double precision :: eps=1d0/4d0, wt
 double precision, parameter :: Pi=4d0*atan(1d0)

 wt=wt0
 if(wt<1d-6) wt=1d-6

 yyy = 1d0 + b*exp(-(eps+nm)**2/(4*wt0))/(6*sqrt(Pi*wt))

end function Snm_yy


function Snm_s2(wt0, b, nm) result(s2)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: s2
 double precision :: eps=1d0/4d0, wt
 double precision, parameter :: Pi=4d0*atan(1d0)

 wt=wt0
 if(wt<1d-6) wt=1d-6

 s2 = 1d0 + b*(erf((2d0*nm+1d0)/(4d0*sqrt(wt)))-erf((2d0*nm-1d0)/(4d0*sqrt(wt))))/6d0

end function Snm_s2


function Snm_s20(wt0, b, nm) result(s2)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: s2
 double precision :: eps=1d0/4d0, wt, swt, acc
 double precision, parameter :: Pi=4d0*atan(1d0)
 double precision, parameter :: sPi=sqrt(4d0*atan(1d0))
 double precision            :: nmm, nmp

 wt=wt0
 if(wt<1d-6) wt=1d-6
 swt = sqrt(wt)

 nmm =nm-1
 nmp =nm+1

 acc =        wt*exp(-(nmm)**2/(4*wt))
 acc = acc +  wt*exp(-(nmp)**2/(4*wt))
 acc = acc -2*wt*exp(-(dble(nm)**2)/(4*wt))
 acc = acc +  sPi*swt*(-nmm)*erf((-nmm)/(2*swt))/2
 acc = acc +  sPi*swt*( nmp)*erf(( nmp)/(2*swt))/2
 acc = acc -  sPi*swt*2*erf(nm/(2*swt))*nm/2

 s2 = 1d0 + b * acc / (3*sPi*swt)


end function Snm_s20
