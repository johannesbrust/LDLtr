!> rosenbrock.f90 is an implementation of the Rosenbrock function
!> 
!> This function is for even number of variables n
!>
!> f = sum(from 1 to n/2 : (100(x2im1^2 - x2i)^2 + (x2im1-1)^2)
!>
!> where x2imi = x_{2i-1} and x2i = x_{2i}
!>
!> The gradient is computed element-wise for even and odd variables
!>
!> df / dx2i    = -200 (x2im1^2 - x2i)
!> df / dx2im1  = 400 x2im1 (x2im1^2 - x2i) + 2 (x2im1-1)
!>
!>-----------------------------------------------------------------
!> 09/01/23, J.B., Initial version

module rosenbrock
    use kind_parameter, only : dp

    implicit none

    contains

    !> 
    !> Implementation of the rosenbrock objective
    !>
    subroutine rosenf(x,f)

        real(dp), intent(in)    :: x(:)
        real(dp), intent(out)   :: f
        
        integer                 :: n

        n       = size(x)
        f       = 0

        !> 
        !> Function assumes that n is even
        !>

        f = sum(100*(x(1:n-1:2)**2-x(2:n:2))**2+(x(1:n-1:2)-1)**2)

        ! do i = 1, n/2

        !     f = f + 100*(x(2*i-1)**2 - x(2*i))**2 + (x(2*i-1)-1)**2

        ! end do

    end subroutine

    !> 
    !> Implementation of the rosenbrock objective
    !>
    subroutine roseng(x,g)

        real(dp), intent(in)    :: x(:)
        real(dp), intent(out)   :: g(:)
        
        integer                 :: n

        n       = size(x)        

        !> 
        !> Function assumes that n is even
        !>
    
        g(2:n:2)    = -200*(x(1:n-1:2)**2-x(2:n:2))
        g(1:n-1:2)  = 400*x(1:n-1:2)*( x(1:n-1:2)**2 - x(2:n:2)) + &
                        2*(x(1:n-1:2)-1)

        ! do i = 1, n/2
        !     g(i) = 400*x(2*i-1)*(x(2*i-1)**2 - x(2*i)) + 2*x(2*i-1)
        ! end do

    end subroutine

end module