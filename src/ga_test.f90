    program gatest

    use GAmodule

    implicit none
    integer, parameter        :: m = 30
    integer                   :: n = 6
    integer                   :: np
    integer                   :: ng
    real                      :: xest(m), xsp(m), mutrate, bestx(m), best_fit
    integer                   :: bitres = 10


    xest(1:n) = 0
    xsp(1:n)  = 2
    mutrate = 0.003
    np = 100
    ng = 100

    write(6,*)'enter: npop, ngen, mutrat, bitres, trace > '
    read(5,*) np,ng,mutrate,bitres,GAtrace
 
    call ga_evolve( xest ,xsp, n, bitres, mutrate, np, ng, bestx, best_fit)

    write(6,*) 'best_fitness= ', best_fit
    write(6,*) 'best_x: ',bestx(1:n)

    end program gatest

    function ga_func(x,n)
!   =====================

    implicit none
    real                          :: ga_func
    real             , intent(in) :: x(n)
    integer          , intent(in) :: n

    integer          :: i
    real             :: y

    y = 1.0
    do i=1,n
      y = y + (x(i)-1.0)**2
    enddo
    
    ga_func = 1.0/y 

    end function ga_func
