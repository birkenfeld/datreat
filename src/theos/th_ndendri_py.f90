 FUNCTION th_ndendri_py(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  dendrimer random chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments chain statistics is expressed ny nu (nu=0.5 ==> Gaussian chain) in addition a chi parameter my be set: S(Q) = 1/(1/(N*P(Q) + 2*chi)
!  mm
      use theory_description 
      implicit none 
      real    :: th_ndendri_py
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor                                                                       
     double precision :: n          ! number of segments per subbranch (do not fit, its integer)                      
     double precision :: gen        ! number of dendrimewr generations                                                
     double precision :: branch     ! number of branching                                                             
     double precision :: l          ! effective segment length Rg = l * N**nu /sqrt(6)                                
     double precision :: nu         ! expansion parameter (Gaussian random walk =0.5)                                 
     double precision :: chi        ! chi parameter                                                                   
     double precision :: withsq     ! if ==1,  multiply with sq_molecule_dendrimer    

     double precision :: rr         ! Percus-Yevick hard-core radius

     double precision :: den        ! Percus-Yevick density in 1/A**3
                                
! the recin parameter representation 
!     double precision ::            !                                                                                 
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision :: th
 
     double precision :: q, pq(4), sq
     double precision :: sq_molecule_dendrimer
     double precision :: py, eps
     double precision :: peryev
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ndendpy'
       nparx =        9
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ndendri_py = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " dendrimer random chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments chain statistics is expressed ny nu (nu=0.5 ==> Gaussian chain) " //cr//parspace//&
                                " the Guinier part only is multiplied by Percus-Yevick S(Q)"
       th_citation(idesc)     = " mm"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments per subbranch (do not fit, its integer)                      
        parnam ( 3) = 'gen     '  ! number of dendrimewr generations                                                
        parnam ( 4) = 'branch  '  ! number of branching                                                             
        parnam ( 5) = 'l       '  ! effective segment length Rg = l * N**nu /sqrt(6)                                
        parnam ( 6) = 'nu      '  ! expansion parameter (Gaussian random walk =0.5)                                 
        parnam ( 7) = 'r_py    '  ! Percus-Yevick hard core radius                                
        parnam ( 8) = 'den_py  '  ! Percus-Yevick density n 1/A**3                               
        parnam ( 9) = 'withsq  '  ! if ==1,  multiply with sq_molecule_dendrimer                                    
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments per subbranch (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "number of dendrimewr generations" !//cr//parspace//&
        th_param_desc( 4,idesc) = "number of branching" !//cr//parspace//&
        th_param_desc( 5,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "expansion parameter (Gaussian random walk =0.5)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "Percus-Yevick hard core radius " !//cr//parspace//&
        th_param_desc( 8,idesc) = "Percus-Yevick density n 1/A**3" !//cr//parspace//&
        th_param_desc( 9,idesc) = "if ==1,  multiply with sq_molecule_dendrimer" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_ndendri_py = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        =      pa( 2)
      gen      =      pa( 3)
      branch   =      pa( 4)
      l        =      pa( 5)
      nu       =  abs(pa( 6))
      rr       =  abs(pa( 7))
      den      =  abs(pa( 8))
      withsq   =      pa( 9)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
!      xh = 
!      call parget('        ',xh,iadda,ier)
!               = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     pq  = dendrimer_sq(q,l,nu, nint(gen), nint(branch), nint(n))
     rg  = pq(4)

     eps = 1d-11
     py  = peryev (q, rr, den, eps)
 
!     sq  = 1d0 / ( 1d0/(pq(2)) + 2d0 * chi )  * pq(1)/pq(2)

     sq = pq(3) * py  + (pq(1)-pq(3))  ! Sq_Guinier * Sq_PY + Sq-Sq_Guinier

     th  = ampli * sq
     if(nint(withsq) == 1) then
        th = th * sq_molecule_dendrimer(q)
     endif



     th_ndendri_py = th
 
! ---- writing computed parameters to the record >>>  
      call parset('rg      ',sngl(pq(4)),iadda)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

function dendrimer_sq(q,l,nue,generations,multiplicity,branchlength) result(sqv)

    implicit none
    double precision, intent(in)       :: l
    double precision, intent(in)       :: q
    double precision, intent(in)       :: nue
    integer         , intent(in)       :: generations
    integer         , intent(in)       :: multiplicity
    integer         , intent(in)       :: branchlength
    double precision                   :: sqv(4) ! [sq normalized, sq, sq-Guinier, rg]

    double precision                   :: sq

    integer, save                      :: segments, nodes
    integer                            :: i, g, iroot
    integer                            :: inode, inode2
    integer                            :: ip1, g1, g2
    integer                            :: N, s1, s2, n1, n2

    integer, save                      :: ipar(3) = 0

    integer, allocatable,    save      :: root(:), gen(:), path(:,:)
    integer, allocatable,    save      :: gdistance(:,:)
    integer, allocatable,    save      :: sdistance(:,:)
    double precision, allocatable      :: eterms(:)

    integer                            :: maxdist
    character(len=80)                  :: fmt


    double precision                   :: ql2, mu
    double precision, save             :: llast 
    double precision, save             :: Rg 

    nodes    = (multiplicity**(generations+1) - multiplicity)/(multiplicity-1)
    segments = nodes * branchlength

!
! DETERMINE TOPOLOGICAL DISTANCE MATRIX
! once for a given topolgy as defined by: generations, multiplicity, branchlength
!
! first the distance beween nodes (implict center = node 0, not expilict in lists)
! --> gdistances
!
! the consider single segment distances, start od segment counting at node towards
! inner part till the previous generation node
! --> sdistanceds
!

ii:    if(maxval(abs(ipar-[generations,multiplicity,branchlength]))>0) then


    write(*,*)"Reestablish distance matrix.."

    llast = 0

    if(allocated(root))      deallocate(root)
    if(allocated(gen))       deallocate(gen)
    if(allocated(path))      deallocate(path)
    if(allocated(gdistance)) deallocate(gdistance)
    if(allocated(sdistance)) deallocate(sdistance)
    allocate(root(nodes+1))
    allocate(gen (nodes+1))
    allocate(path(nodes,0:generations))
    allocate(gdistance(nodes,nodes))
    allocate(sdistance(segments,segments))

    path      = 0
    gdistance = 0
    iroot     = 0


    do inode=1, nodes
      root(inode) = iroot
      if(mod(inode,multiplicity) == 0) iroot = iroot + 1
    !  write(*,*) inode, root(inode)
    enddo


    ! det generation
    do inode=1,nodes
         g = 0
         iroot = inode
         path(inode,0) = iroot
    dl:  do g=0,generations
          iroot = root(iroot)
          path(inode,g+1) = iroot
          if(iroot == 0) exit dl
         enddo dl
         gen(inode) = g
    !     write(*,*) inode, root(inode), gen(inode), 'P:', path(inode,:)
    enddo

    ! get connecting node
    do inode=1,nodes
      do inode2=1,nodes
    dg1:   do g1 = 0, gen(inode)
             ip1 = path(inode,g1)
             do g2 = 0, gen(inode2)
               if(ip1 == path(inode2,g2)) exit dg1
             enddo
           enddo dg1
     !      if(path(inode,g1) > 0) then
     !        write(*,*)"pair ",inode, inode2, gen(inode), gen(inode2)
     !        write(*,*) path(inode ,:)
     !        write(*,*) path(inode2,:)
     !        write(*,*) "g     : ", g1,g2
     !        write(*,*) "g-path: ", g1+g2
     !         write(*,*) "pair: ",inode, inode2,"   g-path: ", g1+g2
               gdistance(inode,inode2) = g1+g2
    !          if(g1>0 .and. g2>0) then
    !           write(*,*) "connect at:", path(inode,g1), path(inode2,g2)
    !         endif
    !       endif
      enddo
    enddo


    !!
    if(nodes < 80) then
       write(*,*)"Topological distance between dendrimer nodes (center is 0)"

       write(*,'(4x    ,*(i3))')(inode,inode=1,nodes)
       write(*,'(3x,"|",*(a3))')("---",inode=1,nodes)
       do inode=1,nodes
        write(*,'(i3,"|",*(i3))') inode, gdistance(inode,:)
       !  write(fmt,'(a,i0,a)') '(i3,"|",',inode,'("..."),*(i3))'
       !  write(*,fmt) inode, gdistance(inode,inode+1:)
       enddo
    endif

    !!!!! now establis segment segment distances

!$OMP PARALLEL DO PRIVATE(inode,inode2,N,n1,n2,s1,s2)

    do inode2=1,nodes
      do inode=1,nodes
         N = gdistance(inode,inode2) * branchlength
         do n1 = 0, branchlength-1
           do n2 = 0, branchlength-1
             s1 = (inode -1)*branchlength + n1 + 1
             s2 = (inode2-1)*branchlength + n2 + 1
             if(N==0) then
                sdistance(s1,s2) = abs(n1 -n2)
             else
                sdistance(s1,s2) = abs(N-n1 -n2)
             endif
           enddo
         enddo
      enddo
    enddo
 !$OMP END PARALLEL DO

    !!
    if(segments < 100) then
       write(*,*)"Topological distance between segments (center is 0)"

       write(*,'(4x    ,*(i3))')(s1,s1=1,segments)
       write(*,'(3x,"|",*(a3))')("---",s1=1,segments)
       do s1=1,segments
        write(*,'(i3,"|",*(i3))') s1, sdistance(s1,:)
       !  write(fmt,'(a,i0,a)') '(i3,"|",',inode,'("..."),*(i3))'
       !  write(*,fmt) inode, gdistance(inode,inode+1:)
       enddo

       write(*,*)"Topological symmetry check"

       write(*,'(4x    ,*(i3))')(s1,s1=1,segments)
       write(*,'(3x,"|",*(a3))')("---",s1=1,segments)
       do s1=1,segments
        write(*,'(i3,"|",*(i3))') s1, (sdistance(s1,:)-sdistance(:,s1))
       !  write(fmt,'(a,i0,a)') '(i3,"|",',inode,'("..."),*(i3))'
       !  write(*,fmt) inode, gdistance(inode,inode+1:)
       enddo
    endif

    ipar = [generations,multiplicity,branchlength]

endif ii

    maxdist = branchlength*generations*2
    ! maxdist = maxval(sdistance)

    ! write(*,*)"Maxdist=",maxdist, branchlength*generations*2

!! DETERMINATION OF SQ
!! by summation of the Gaussian (possibly nu-modified) over all
!! tpological distances (i.e. |n1-n2| in simple cases)


    if(allocated(eterms)) deallocate(eterms)
    allocate(eterms(0:maxdist))

      ql2 = (1d0/6d0) * (q*l)**2
      mu  = 2d0 * nue

!$OMP PARALLEL DO
       do i=0,maxdist
         eterms(i) = exp(-ql2*dble(i)**mu)
       enddo
!$OMP END PARALLEL DO

       sq = 0
!$OMP PARALLEL DO REDUCTION(+:sq)
       do s2=1,segments
        do s1=1,segments
          sq = sq + eterms(sdistance(s1,s2))
        enddo
       enddo
!$OMP END PARALLEL DO

      sq = sq/(segments**2)
      sqv(1) = sq
      sqv(2) = sq * segments

!! and the radius of gyration
      if( l .ne. llast ) then
!$OMP PARALLEL DO
         do i=0,maxdist
           eterms(i) = dble(i)**mu
         enddo
!$OMP END PARALLEL DO
         Rg = 0
!$OMP PARALLEL DO REDUCTION(+:Rg)
         do s2=1,segments
          do s1=1,segments
            Rg = Rg + eterms(sdistance(s1,s2))
          enddo
         enddo
!$OMP END PARALLEL DO
        Rg    = sqrt( (0.5d0 * l**2) * Rg/(segments**2) )
        llast = l
      endif
       
      sqv(3) = exp(-0.33333333333d0* (Rg*q)**2)
      sqv(4) = Rg


end function dendrimer_sq

 end function th_ndendri_py
