!> test_chol
!>
!> Test program to test the column oriented Cholesky factorization
!>
!>-----------------------------------------------------------------
!> 09/01/23, J.B., initial implementation

program test_chol

    use kind_parameter, only : dp
    use ldlqn_chol, only : cchol

    implicit none

    ! Initialization
    integer     :: n    
    real(dp)    :: R(3,3)
    real(dp)    :: A(3,3)
    
    n = size(A,dim=1)
    A = reshape([1,0,0,2,3,0,1,2,3],shape(A))
    A = matmul(transpose(A),A)

    print *, 'Test Cholesky'
    print *, 'A'
    print *, A
    print *, ''
    print *, 'R'
    call cchol(A,R)
    print *, R
    print *, ''

end program test_chol