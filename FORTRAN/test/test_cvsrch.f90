!> test_cvsrch
!> Calling the cvsrch strong Wolfe line-search
!>
!> This test program uses the rosenbrock function to try accessing
!> the F77 line-search
!> 
!>-----------------------------------------------------------------
!> 09/01/23, J.B., Initial implementation
!> 09/06/23, J.B., Use of explicitly defined arrays to call 
!>                  subroutine cvsrch.f

program test_cvsrch

    use kind_parameter, only : dp
    use rosenbrock
    !use feval
    !use geval

    implicit none

    !>
    !> Interface for external subroutines
    !>
    ! interface 

    !     subroutine feval(n,x,f)
        
    !         use kind_parameter, only : dp

    !         integer, intent(in)     :: n
    !         real(dp), intent(in)    :: x(n)            
    !         real(dp), intent(out)   :: f
        
    !     end subroutine

    !     subroutine geval(n,x,g)
        
    !         use kind_parameter, only : dp

    !         integer, intent(in)     :: n
    !         real(dp), intent(in)    :: x(n)            
    !         real(dp), intent(out)   :: g(n)
                
    !     end subroutine


    ! end interface

    integer         :: n, info, nfev, i, maxfev
    real(dp)        :: f, x(10), g(10), s(10), wa(10), &
                        stp, ftol, gtol, xtol, stpmin, &
                        stpmax
  
    n = 10

    x = [(0, i = 1,10)]
    wa(:) = x(:)
    !call rosenf(x,f)
    !call roseng(x,g)

    call feval(n,x,f)
    call geval(n,x,g)

    !>
    !> For debugging
    !>
    !>----------------
    ! print *, 'Rosen g'
    ! print *, g
    ! print *, ''
    ! print *, 'Rosen ng'
    ! print *, sqrt(sum(g*g))
    ! print *, ''
    !>----------------

    s = -g / sqrt(sum(g*g))
    
    !>
    !> Parameters for cvsrch for ldlqn
    !>

    stp     = 1.0_dp
    gtol    = 0.9_dp
    ftol    = 1.0e-4_dp
    xtol    = 1e-12_dp
    stpmin  = 1e-14_dp
    stpmax  = 20_dp
    maxfev  = 200     
    nfev    = 1

    !>
    !> Additional printing for debugging
    !>
    ! print *, f
    ! print *, g
    ! print *, ''

    !> 
    !> Print values before and after calling cvsrch
    !>
    print '(a5,f6.3)', 'Curv: ', dot_product(g,s)
    print '(a,a,a)', '    it   ', 'fx    ', 'ngx   '
    print '(i5,f6.3,f6.3)', 0, f, sqrt(sum(g*g))      

    call cvsrch(n,x,f,g,s,stp,ftol,gtol,xtol, &
                      stpmin,stpmax,maxfev,info,nfev,wa)

    !print '(a,a)', 'fx', 'ngx'
    print '(i5, f6.3,f6.3)', 1, f, sqrt(sum(g*g)) 
    print *, '' 
    print '(a5,i6)', 'Info: ', info
    print '(a5, f6.3)', 'Step:', stp
    
    
    ! Parameters
    print  '(a5,f6.3,f6.3,f6.3,f6.3,f6.3,f6.2,i5,i5)', 'Pars:', 1.0_dp, &
     gtol, ftol, xtol, stpmin, stpmax, maxfev, nfev
    print *, 'x:'
    print *, x
    ! print *, g
    ! print *, f

    !contains

    !> 
    !> Subroutines for cvsrch
    !> Both function and gradient evaluations
    !>

end program test_cvsrch

!interface
subroutine feval(n,x,f)

    use kind_parameter, only : dp
    use rosenbrock

    implicit none

    integer, intent(in)     :: n
    real(dp), intent(in)    :: x(n)    
    real(dp), intent(out)   :: f

    ! print *, 'Rosen f (b):'
    ! print *, f
    ! print *, 'x (b)'
    ! print *, x

    call rosenf(x,f)

    !f = sum(100*(x(1:n-1:2)**2-x(2:n:2))**2+(x(1:n-1:2)-1)**2)

    !print '(a6,f3.6)', 'Rosen:', f
    ! print *, 'Rosen f (a):'
    ! print *, f
    ! print *, 'x (a)'
    ! print *, x
    !print *, n

end subroutine
!end interface

!interface
subroutine geval(n,x,g)

    use kind_parameter, only : dp
    use rosenbrock

    implicit none

    integer, intent(in)     :: n
    real(dp), intent(in)    :: x(n)    
    real(dp), intent(out)   :: g(n)

    call roseng(x,g)

    !g(2:n:2)    = -200*(x(1:n-1:2)**2-x(2:n:2))
    !g(1:n-1:2)  = 400*x(1:n-1:2)*( x(1:n-1:2)**2 - x(2:n:2)) + &
    !                   2*(x(1:n-1:2)-1)

end subroutine
!end interface



