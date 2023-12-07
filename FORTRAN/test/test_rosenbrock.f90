!> test_rosenbrock
!> Test of the implementation of the Rosenbrock function
!>
!> The function has a global minimum at x = 1
!>
!> ------------------------------------------------------
!> 09/01/23, J.B., Initial implementation

program test_rosenbrock

    use kind_parameter, only: dp
    use rosenbrock

    implicit none

    integer     :: i
    real(dp)    :: x1(10)
    real(dp)    :: x2(10)
    real(dp)    :: x3(10)

    real(dp)    :: f
    real(dp)    :: g(10)

    !> 
    !> Defining a set of vectors to test the function, and 
    !> subsequentially printing all results
    !>
    x1 = [(10, i = 1, 10)]
    x2 = [(0, i = 1, 10)]
    x3 = [(1, i = 1, 10)]

    call rosenf(x1,f)
    call roseng(x1,g)
    print *, 'f(x1)'
    print *, f
    print *, 'g(x1)'
    print *, g
    print *, ''

    call rosenf(x2,f)
    call roseng(x2,g)
    print *, 'f(x2)'
    print *, f
    print *, 'g(x2)'
    print *, g

    call rosenf(x3,f)
    call roseng(x3,g)
    print *, 'f(x3)'
    print *, f
    print *, 'g(x3)'
    print *, g

end program