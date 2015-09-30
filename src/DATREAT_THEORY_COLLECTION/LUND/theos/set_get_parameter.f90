
subroutine setpar (pname,pvalue,nopar ,params,napar,mbuf, ier) 
!      ==============================from datreat nopar(iadda),params(:,iadda),napar(:,iadda)=========                        
	character*8   pname            !parametername
	real          pvalue           ! parametervalue
	integer       mbuf             ! maximum lenght of following arrays  ; maximum of nopar
	integer       nopar,ier        ! number of used Parameter, error 
	character*80  napar(mbuf)      ! name of parameters n
	real          params(mbuf)     ! value of parameters n
!                                                                       
! ----------- look for the parameter -----------------------------------
	ier=1
	do i=1,nopar
		if ( trim(napar(i)).eq.trim(pname) ) then
			params(i) = pvalue 
			ier = 0 
			return
		endif
	enddo 
	! parameter not found so appand new one
	if (nopar+1.le.mbuf) then
		nopar = nopar + 1 
		params(nopar) = pvalue 
		napar(nopar) = pname 
		ier=0
	else
		ier=2
	endif
	return   
END  subroutine  

                     
subroutine getpar (pname,pvalue,nopar ,params,napar,mbuf, ier) 
!                 from datreat  nopar(iadda),params(:,iadda),napar(:,iadda)                                 
	character*8     pname                              !parametername
	real            pvalue               ! parametervalue
	integer         mbuf                 ! maximum lenght of following arrays  ; maximum of nopar
	integer         nopar,ier            ! number of used Parameter, error 
	character*80    napar(mbuf)          ! name of parameters n
	real            params(mbuf)         ! value of parameters n
		
!                                                                       
!     get the parameters from the data 
!     searched by the name of the parameters and assigned to variables
!     old was call to  parget from data
	ier=1
	do ii=1 , nopar 
		if (trim(pname) .eq. trim(napar(ii))) then
			pvalue = params(ii)
			ier=0
			return
		endif
	enddo
	! ier=2 says the pnam was not found in the napar list
	ier=2  
	return
END subroutine
