!> test_ldlqn1
!> Calling the ldlqn method
!>
!> This test program uses the rosenbrock function to try calling 
!> the algorithm
!> 
!>-----------------------------------------------------------------
!> 09/11/23, J.B., Initial implementation
!> 12/08/23, J.B., Preparation for release

program test_ldlqn1

    use kind_parameter, only : dp
    use rosenbrock
    use ldlqn
    !use feval
    !use geval

    implicit none

    integer         :: n, itn, nf, status, i

    real(dp)                :: f
    real(dp), allocatable   :: x(:)
    
    !> Initialize starting point
    n = 1000
    allocate(x(n),source=0.0_dp)

    !> Initialize starting point
    !x = [(0, i = 1,n)]
    
    !>
    !> Call the main algorithm
    !> 
    call ldlqn_bfgswtr(x, itn, nf, status)

    print *,''
    print '(a8,i2)', 'Status: ', status

end program

!>
!> Interface to the objective function and its gradient
!>

!interface
subroutine feval(n,x,f)

    use kind_parameter, only : dp
    use rosenbrock

    implicit none

    integer, intent(in)     :: n
    real(dp), intent(in)    :: x(n)    
    real(dp), intent(out)   :: f

    call rosenf(x,f)

end subroutine

! interface
subroutine geval(n,x,g)

    use kind_parameter, only : dp
    use rosenbrock

    implicit none

    integer, intent(in)     :: n
    real(dp), intent(in)    :: x(n)    
    real(dp), intent(out)   :: g(n)

    call roseng(x,g)

end subroutine
!end interface



