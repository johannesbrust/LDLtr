!>
!> test_invLDLCG
!>
!> Test the inverse LDLT CG iteration on a small problem
!> This is a Fortran translation of the function
!> inv_ldl_CG1.m
!>
!>------------------------------------------------------------
!> 09/11/23, J.B., Initial version

program test_invLDLCG
    
    use kind_parameter, only: dp
    use ldlqn_invLDLCG ! module depends on ldlqn_triang

    implicit none

    !>
    !> Initialize test data
    !>

    integer     :: n, i    
    real(dp)    :: R(5,5)
    real(dp)    :: s(5)
    real(dp)    :: p(5)
    real(dp)    :: res(5)
    real(dp)    :: g(5)
    real(dp)    :: D(5)

    integer     :: maxit
    real(dp)    :: sig
    real(dp)    :: tol

    n = size(R,dim=1)
    R = reshape([1,0,0,0,0,& 
                 2,3,0,0,0,&
                 1,2,3,0,0,&
                 1,1,2,2,0,&
                 1,2,3,4,5]&
                 ,shape(R))    

    D   = [(2, i = 1,5)]
    g   = [(1, i = 1,5)]

    sig = 10.0
    tol = 1e-7
    maxit = 15

    !>
    !> Call to the main iteration
    !>

    print *, "Inverse LDLT CG iteration"
    print *, "s (b):"
    print *, s
    call invLDLCG(s, p, res, sig, g, R, D, tol, maxit )
    
    print *, "s (a):"
    print *, s

    print *, "norm(res):"
    print *, norm2(res)

end program