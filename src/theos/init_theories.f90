!################ 
subroutine init_theories(thenam,thparn,nthpar,thrapar,thramin,thramax,mth,mtpar)
    
    real x,par(20)       
    integer mth,mtpar,iadda
    integer nthpar(mth)
    character*8 , intent(inout):: thenam(mth)
    character*8 , intent(inout):: thparn(mtpar,mth) 
    character*8 , intent(inout):: thrapar(mth)
    real*4      , intent(inout):: thramin(mth),thramax(mth)
!  thenam(i)  = name of i-th theory                                     
!  thparn(j,i)= name of j-th parameter of i-th theory                   
!  nthpar(i)  = no of parameters required for i-th theory               
!  thparx(j,l)= parameters for the l-th activated theory                
!  thpsca(j,l)= corresponding fit scales   "                            
!  nthtab(l)  = no. of l-th activated theory                            
!  ntheos     = total no. of activated theories                         
!  multflg    = flag indicating a multiplicative theory 

!   integer  mbuf             ! maximum length of data arrays in calling main program   
!   integer  nopar            ! number of Parameter in data
!   character*80  napar(mbuf) ! name of  parameters(n)
!   real params(mbuf)         ! value of parameters(n)
        
               
    do i=1,mth                    
          thenam(i)  = ' '         
      thrapar(i) = '        '   
      thramin(i) = -1e33        
      thramax(i) = 1e33         
    enddo                        
! ---- initialize theory tables ----
! the include file contains this to initialise all availible theories
!    nthpar(1) = mtpar 
!    dummy = th_zglam(x,par,thenam(1),thparn(:,1),nthpar(1),ini,idum ,idum,idum,idum)
      INCLUDE "ini_theories.i90" 
        return
    
    
end 
!   ###################################################
    
real function theory_x(x,ntheos,nthtab,multflg,thparx,thrapar,thramin,thramax,nopar,params,napar,mbuf,mth,mtpar,mtcal) 
!      ----------------------                              from datreat nopar(iadda),params(:,iadda),napar(:,iadda)
! ---- compute the value of the active set of theories at value x       
!                                                                       
       character*8 thenam,thparn 
        integer nthtab(mtcal),multflg(mtcal),ntheos
        real thparx(mtpar,mtcal)
            !  thenam(i)  = name of i-th theory                                     
            !  thparn(j,i)= name of j-th parameter of i-th theory                   
            !  nthpar(i)  = no of parameters required for i-th theory               
            !  thparx(j,l)= parameters for the l-th activated theory                
            !  nthtab(l)  = no. of l-th activated theory                            
            !  ntheos     = total no. of activated theories                         
            !  multflg    = flag indicating a multiplicative theory      
             
            character*8 thrapar(mth) 
            real*4     thramin(mth), thramax(mth)
        !                            here only the params naparnopar of the current dataset is known 
        !                            (1 dim array like params(:,iadda)  in call of this function) 
            integer nopar,ier,mbuf
            real params(mbuf)
            character*80 napar(mbuf)
            LOGICAL :: imu, iins 
            real  pr_val 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                       
                         
            sum = 0 
! ----> sum over selected theories ----                                 
            
            do  it=1,ntheos 
                ith = nthtab(it) 
                imu  = (multflg(it).eq.1)
                iins = (multflg(it).eq.2)
                if(thrapar(it).ne.'        ') then 
                        call getpar(thrapar(it),pr_val ,nopar ,params,napar,mbuf,ier)
                        if(ier.ne.0)                    cycle 
                        if(pr_val.lt. thramin(it))      cycle 
                        if(pr_val.gt. thramax(it))      cycle 
                endif 
! ------ include routine zum Summieren und multiplizieren von Theorien S
        INCLUDE "sum_theories_up.i90"
            enddo 
       theory_x = sum 
       return 
      END  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
               
! ------ include routine zum Summieren und multiplizieren von Theorien S
!        include has for all theories:
!!1        continue
!!           if(imu) then
!!              sum=sum*th_zglam(x ,thparx(1,it),dum,dum,idum,1,nopar,params,napar,mbuf)
!!           else
!!               sum=sum+th_zglam(x ,thparx(1,it),dum,dum,idum,1,nopar,params,napar,mbuf)
!!           endif
!!           goto 250        
