      FUNCTION HN (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> rousei <--------                                              
!                                                                       
!                
      implicit none                                                       
      CHARACTER(len=8), intent(inout) :: thnam, parnam (20) 
      real, intent(in)                :: x,  pa (20)
      integer, intent(in)             :: ini, mbuf
      integer, intent(inout)          :: npar
      integer, intent(inout)          :: nopar            ! Anzahl der Parameter data
      character*80, intent(inout)     :: napar(mbuf)      ! name des parameters n
      real, intent(inout)             :: params(mbuf)     ! value des parameters n
      real                            :: HN

      double precision :: q1, a1, b1, k1
      double precision :: hn1, dfmax1, dtmax1, deps1
      double precision :: da1, db1, dtau1, dk1

      integer           :: ier, nparx

    double precision, parameter :: pi                  =     3.141592653589793238462643d0
    

                                                                        
    integer iadda
    common/thiadd/iadda
                                                                   
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'HN' 
         nparx = 4 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            HN = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters   units: micro ns A , macro g, cm, mol                
         parnam (1)  = 'broadn' 
         parnam (2)  = 'sym' 
         parnam (3)  = 'tau' 
         parnam (4)  = 'int' 

!                                                                       
         HN = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      q1       = x 
      da1      = pa (1) 
      db1      = pa (2) 
      dtau1    = pa (3) 
      dk1      = pa (4)
  

            hn1=dk1*((1+(2*((q1*2*pi*dtau1)**(a1))*&
      (dcos(pi*a1/2)))+((q1*2*pi*dtau1)**(2*a1)))**&
      (-b1/2))*(dsin(b1*(datan((dsin(pi*a1/2))/&
      (((q1*2*pi*dtau1)**(-a1))+(dcos(pi*a1/2)))))))



           dfmax1=(((dsin((pi*a1)/(2*&
      (b1+1)))/dsin((pi*a1*b1)/(2*&
      (b1+1))))**(1/a1))*(dtau1**(-1)))/(2*pi)


           dtmax1=(dfmax1**(-1))

           deps1=dk1*((1+(2*((dfmax1*2*pi*dtau1)**(a1))*&
      (dcos(pi*a1/2)))+((dfmax1*2*pi*dtau1)**(2*a1)))**&
      (-b1/2))*(dsin(b1*(datan((dsin(pi*a1/2))/&
      (((dfmax1*2*pi*dtau1)**(-a1))+(dcos(pi*a1/2)))))))


 
      HN =  hn1
!      

      
      call parset('dfmax1     ',sngl(dfmax1),iadda) 
      call parset('dtmax1     ',sngl(dtmax1),iadda) 
      call parset('deps1      ',sngl(deps1),iadda) 
                                                                 
      RETURN 
      END
