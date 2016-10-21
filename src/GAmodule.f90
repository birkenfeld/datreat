

   Module GAmodule

     integer               :: GAtrace = 1                ! Output tracing

     integer , parameter   :: maxindiv = 1000            ! maximum poulation size
     integer , parameter   :: maxparen = 2               ! maximum number of parents
     integer , parameter   :: maxchild = 2               ! maximum number of childs
     integer , parameter   :: maxchrom = 1000            ! maximum length of chromosome

     integer               :: number_of_parents    = 2
     integer               :: number_of_childs     = 2
     real                  :: mutation_probability = 0.002
     
     integer               :: bitcoding = 10             ! length of representation of bit-coding

     real            , private :: x_center(maxchrom)
     real            , private :: x_range(maxchrom)
     integer, private          :: n                      ! length of x vectors

     type animal
          integer             :: chromosome(maxchrom)       
          integer             :: nchrom
          real                :: x(maxchrom)             ! value vector after decoding of chromosome
          integer             :: nx                      ! length of value vector
          real                :: fitness
          integer             :: parents(maxparen)
          integer             :: cross_site
          integer             :: childs(maxchild)
          integer             :: nparents
          integer             :: nchilds
     end type animal


     type population
          type(animal)        :: individual(maxindiv)
          integer             :: size
          integer             :: fittest
          real                :: best_fitness
          real                :: worst_fitness
          real                :: average_fitness
          real                :: variance_fitness
     end type population


     type (population)        :: new_generation 
     type (population)        :: old_generation


   contains

   function gaevent( probability )
!  -------------------------------

!  creates true randomly with probability

    implicit none
    logical gaevent
    
    real             , intent(in)   :: probability
    real                            :: rand_number

    call random_number(rand_number)
    
    gaevent = (rand_number .lt. probability)

   end function gaevent     


   subroutine create_generation( generation, size, nchrom )
!  --------------------------------------------------------
!
! creates random population of size
!
   implicit none
   type(population), intent(inout) :: generation
   integer         , intent(inout) :: size
   integer         , intent(inout) :: nchrom

   ! logical                         :: gaevent
   integer                         :: i, j

   if(size.gt.maxindiv) then
     write(6,*)'create_generation: error: size=',size,' exceeds max=',maxindiv
     size = 0
     return
   endif

   if(nchrom.gt.maxchrom) then
     write(6,*)'create_generation: error: nchrom=',nchrom,' exceeds max=',maxchrom
     nchrom = 0
     return
   endif

   do i=1,size
     generation%individual(i)%fitness   = 0
     generation%individual(i)%parents   = 0
     generation%individual(i)%childs    = 0
     generation%individual(i)%nparents  = 0
     generation%individual(i)%nchilds   = 0
     generation%individual(i)%nchrom    = nchrom
     do j=1,nchrom
       if(gaevent(0.5)) then
         generation%individual(i)%chromosome(j) =  1
       else
         generation%individual(i)%chromosome(j) =  0
       endif 
     enddo
   enddo

   generation%size                    = size
   generation%fittest                 = 0
   generation%best_fitness            = 0
   generation%worst_fitness           = 0
   generation%average_fitness         = 0
   generation%variance_fitness        = 0

    
   if(GAtrace.gt.0) write(6,*)'GAtrace: generation created with size: ', size

   return

   end subroutine create_generation


   subroutine ga_xpress( chromosome, nchrom, x )
!  ---------------------------------------------
!
! expresses the contents of the chromosome into the real
! number vector x of length n
!
  implicit none
  integer, intent(in)           :: nchrom               ! length of chromosome 
  integer, intent(in)           :: chromosome(nchrom)   ! contents
  real            , intent(out) :: x(n)                 ! output x-vector


  integer          :: nbits                             ! number of bits
  real             :: f, fnorm
  
  integer          :: i, j

  ! write(6,*)'ga_xpress##: ',nchrom, n

  if(modulo(nchrom,n) .ne. 0) then
   write(6,*)'ga_xpress warning nchrom is not integer multiple of n ', nchrom,n
  endif

  nbits = nchrom / n

  if(nbits.lt.2) then
    write(6,*)'ga_xpress warning nchrom is too small for n ', nchrom,n
  endif

  fnorm = 2**nbits
  do i=1,n
    f = 0
    do j=1,nbits
      f = f + chromosome(nbits*(i-1)+j) * 2**(j-1)
    enddo
    x(i) = x_center(i)-x_range(i)*(f/fnorm-0.5d0)
  enddo


  if(GAtrace.gt.1) then
    write(6,*)'GAtrace: ga_xpress chromosome: '
    write(6,'(80I1)')chromosome(1:nchrom)
    write(6,*)'x='
    write(6,*)x(1:n)
  endif


  return
  end subroutine ga_xpress
  

  function select_individual(generation, sum_fitness)
! ---------------------------------------------------
  implicit none

  integer                          :: select_individual
  type(population) , intent(in)    :: generation
  real             , intent(inout) :: sum_fitness

  integer                          :: i
  real                             :: rannum, ranlevel, accumulated_fitness

  if(sum_fitness.eq.0d0) then
    do i=1, generation%size
     sum_fitness = sum_fitness + generation%individual(i)%fitness
    enddo
  endif

  call random_number(rannum)
  ranlevel = rannum * sum_fitness

  accumulated_fitness = 0
  do i=1, generation%size
     accumulated_fitness = accumulated_fitness + generation%individual(i)%fitness
     if(accumulated_fitness.ge.ranlevel) then
       select_individual = i
       return                                  ! ------------>
     endif
  enddo
  
  select_individual = generation%size          ! ------------>
  return
 
  end function select_individual


  function ga_breeding(parents,generation, next_generation)
! ---------------------------------------------------------
  implicit none
  integer                         :: ga_breeding      ! on exit size of nex_generation
  integer, intent(in)             :: parents(maxparen)
  type(population), intent(in)    :: generation
  type(population), intent(inout) :: next_generation

  integer                         :: i,k,l, i1,i2
  real                            :: rannum
  integer                         :: cross_site

  ! determine cross-site
  k = next_generation%size
  if(k+2.gt.maxindiv) then
    write(6,*) 'ga_breed: next population size becomes too large ',k, maxindiv
    return
  endif

  call random_number(rannum)
  cross_site = nint(generation%individual(parents(1))%nchrom * rannum + 0.5d0)

  next_generation%individual(k+1) = generation%individual(parents(1))
  next_generation%individual(k+2) = generation%individual(parents(2))

  do l=1,cross_site
    next_generation%individual(k+1)%chromosome(l) =  generation%individual(parents(1))%chromosome(l) 
    next_generation%individual(k+2)%chromosome(l) =  generation%individual(parents(2))%chromosome(l) 
  enddo
  do l=cross_site,generation%individual(parents(1))%nchrom
    next_generation%individual(k+1)%chromosome(l) =  generation%individual(parents(2))%chromosome(l) 
    next_generation%individual(k+2)%chromosome(l) =  generation%individual(parents(1))%chromosome(l) 
  enddo

  ! mutation
  do l=1,generation%individual(parents(1))%nchrom
    if(next_generation%individual(k+1)%chromosome(l).eq.0) then
      i1 = 1
    else
      i1 = 0
    endif 
    if(next_generation%individual(k+2)%chromosome(l).eq.0) then
      i2 = 1
    else
      i2 = 0
    endif 
    if(gaevent(mutation_probability)) then 
      next_generation%individual(k+1)%chromosome(l) = i1
      if(GAtrace.gt.0) write(6,*)'mutation: ',k+1,l,i1
    endif
    if(gaevent(mutation_probability)) then 
      next_generation%individual(k+2)%chromosome(l) = i2
      if(GAtrace.gt.0) write(6,*)'mutation: ',k+2,l,i2
    endif
 
  enddo 
  
  next_generation%individual(k+1)%parents = parents
  next_generation%individual(k+2)%parents = parents
  
  next_generation%individual(k+1)%cross_site = cross_site
  next_generation%individual(k+2)%cross_site = cross_site
   

  next_generation%size = k+2

  if(GAtrace.gt.1) then
   write(6,'(a,2i5,a,2i5)')'GAtrace: ga breeding parents= ',parents, ' --> ',k+1,k+2
  endif


  ga_breeding = next_generation%size

  return 
  end function ga_breeding 

  

   function ga_contest(generation)
!  -------------------------------
   implicit none
   integer                         :: ga_contest    ! returns the address of fittest individual
   type(population), intent(inout) :: generation

   integer                         :: i
   real                            :: fitness

 !   write(6,*)'###############0:',generation%size    


   generation%fittest                 = 1
   generation%best_fitness            = ga_fitness(generation%individual(1))
   generation%worst_fitness           = generation%best_fitness
   generation%average_fitness         = generation%best_fitness
  
  
   do i=2,generation%size
     fitness = ga_fitness(generation%individual(i))
     if(fitness .gt. generation%best_fitness) then
       generation%fittest      = i
       generation%best_fitness = fitness
     endif

     if(fitness .lt. generation%worst_fitness) then
       generation%worst_fitness = fitness
     endif

     generation%average_fitness         = generation%average_fitness + fitness

   enddo                    

 !  write(6,*)'###############1:',generation%size    

   generation%average_fitness =   generation%average_fitness / generation%size
   
   generation%variance_fitness   = 0
   do i=1,generation%size
     generation%variance_fitness = generation%variance_fitness + (generation%individual(i)%fitness-generation%average_fitness)**2
   enddo
 
   generation%variance_fitness = sqrt(generation%variance_fitness / generation%size)


   ga_contest = generation%fittest


   if(GAtrace.gt.0) then
     write(6,*)'GAtrace: ga_contest population=', generation%size
     write(6,*)'Fittest           = ', generation%fittest
     write(6,*)'Average Fitness   = ', generation%average_fitness
     write(6,*)'Best    Fitness   = ', generation%best_fitness
     write(6,*)'Worst   Fitness   = ', generation%worst_fitness
     write(6,*)'VarianceFitness   = ', generation%variance_fitness
   endif

   return
   end function ga_contest

   function ga_fitness(individual)
!  -------------------------------
   implicit none
   real                            :: ga_fitness
   type(animal), intent(inout)     :: individual

   real                            :: x(maxchrom)   ! differnt transmission needed

   real            , external      :: ga_func

 !  x_center(1:maxchrom) = 0
 !  x_range(1:maxchrom)  = 1d0

 !  write(6,*)'##### ga_fitness nchrom,bitcoding :',individual%nchrom, bitcoding

   n = individual%nchrom/bitcoding  !  oder global

   call ga_xpress( individual%chromosome, individual%nchrom, x)


   ga_fitness = ga_func(x,n)   ! < HIER INTERFACE

   individual%fitness = ga_fitness
   individual%x       = x
   individual%nx      = n

   return
   end function ga_fitness


   subroutine ga_evolve( xestimate ,xspan, nx, bitresolution, mutation_rate, npop, ngen, best_x, best_fitness)
!  -----------------------------------------------------------------------------------------------------------
   implicit none
   real             , intent(in)    :: xestimate(nx)      ! center of evaluation range
   real             , intent(in)    :: xspan(nx)          ! span of evaluation range
   integer          , intent(in)    :: nx                 ! dimension of x
   integer          , intent(in)    :: bitresolution      !
   real             , intent(in)    :: mutation_rate      !
   integer          , intent(inout) :: npop               ! size of population in a generation
   integer          , intent(in)    :: ngen               ! number of generations
   real             , intent(out)   :: best_x(nx)         ! final x
   real             , intent(inout) :: best_fitness       ! input criterion, output best    
   
   integer                          :: i, j, ibest, nchrom
   integer                          :: new_size
   integer                          :: next_parents(maxparen)
   type(population)                 :: actual_generation
   type(population)                 :: next_generation
   type(population)                 :: super_generation
   type(animal)                     :: best_individual
 
   real                             :: sum_fitness

  ! integer                          :: ga_contest

   ! xvalues to global (private)   SUBROUTINes hier uaf gloal xcenter ... umstellen
   ! n's, mutation rate ... setzen

   write(6,*)'========== ga_evolve =============='
   write(6,*)' population size : ', npop
   write(6,*)' bitresolution   : ', bitresolution
   write(6,*)' mutation rate   : ', mutation_rate
   write(6,*)' generations     : ', ngen
   write(6,*)' variables       : ', nx
   write(6,*)'==================================='

   if(nx <= 0) then
     call errsig(999,"ERROR GA evolve: number of variables is zero!$")
     return
   endif

   mutation_probability     = mutation_rate
   bitcoding                = bitresolution

   actual_generation%size   = 0
   next_generation%size     = 0
   super_generation%size    = 0

   best_individual%fitness  = 0

   x_center(1:nx) = xestimate(1:nx)
   x_range(1:nx)  = xspan(1:nx)
   n        = nx

   nchrom = nx * bitresolution    

   call create_generation(actual_generation, npop, nchrom )

   new_size = actual_generation%size

   do i=1,ngen

 !    write(6,*)'#0# gen = ',i, 'new_size =',new_size,actual_generation%size, n

      ibest           =    ga_contest(actual_generation)
 !     write(6,*)'#ibest= ',ibest

      sum_fitness     =    actual_generation%size * actual_generation%average_fitness
      if(best_individual%fitness .lt. actual_generation%individual(ibest)%fitness) then 
         best_individual =    actual_generation%individual(ibest)  
         if(GAtrace.gt.-1) write(6,*)'fitness(',i,'): ',best_individual%fitness,best_individual%x(1:nx) 
      else
!         actual_generation%individual(npop/2) = best_individual
! hier sind auch andere Ersetzungen denkbar: random, worst ....  !!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! verhindert aber offenbar das weitere explorieren aus einem lokalen Minimum, also so erstmal weglassen!  <
!
      endif
      
      do while (next_generation%size .lt. npop)
        next_parents(1) =    select_individual(actual_generation, sum_fitness)
        next_parents(2) =    select_individual(actual_generation, sum_fitness)  ! extra code to avoid inbreeding ?
        new_size        =    ga_breeding(next_parents, actual_generation, next_generation)      
      enddo

 !     write(6,*)'### gen = ',i, 'new_size =',new_size,actual_generation%size, next_generation%size, n

      actual_generation    = next_generation   !! ineffizient bei grosser dimesion der populationsize!!! <<<<<<<<<<<<<
      next_generation%size = 0
 !    write(6,*)'#2# gen = ',i, 'new_size =',new_size,actual_generation%size, next_generation%size, n

   enddo

   best_x(1:nx)       = best_individual%x(1:nx)
   best_fitness       = best_individual%fitness

   return
   end subroutine ga_evolve


   end Module GAmodule

