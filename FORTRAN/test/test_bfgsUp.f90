!>
!> test_bfgsUp
!>
!> Test of updating the inverse BFGS LDLT factorization
!>
!>----------------------------------------------------------
!> 09/06/23, J.B., Initial implementation

program test_bfgsUp

    use kind_parameter, only: dp
    use ldlqn_bfgsRDRup

    implicit none

    !>
    !> Initialize test data
    !>

    integer     :: n, i    
    real(dp)    :: R(3,3)
    real(dp)    :: A(3,3)
    real(dp)    :: sk(3)
    real(dp)    :: yk(3)
    real(dp)    :: Hy(3)
    
    real(dp)    :: DI(3)
    real(dp)    :: GI(3)

    real(dp)    :: nsk
    real(dp)    :: nyk

    n = size(A,dim=1)
    R = reshape([1,0,0,2,3,0,1,2,3],shape(A))
    !A = matmul(transpose(A),A)

    DI  = [(1, i = 1,3)]
    GI  = [(2, i = 1,3)]
    sk  = [1.0, 2.0, 3.0]
    yk  = [0.0, 2.0, 1.0]
    Hy  = [-1.0,0.0,1.0]
    nsk = norm2(sk)
    nyk = norm2(yk)

    print *, 'Test LDLT Update'
    print *, 'R (o)'
    print *, R
    print *, ''
    print *, 'R (up)'
    call bfgsRDRupEig( R, DI, GI, sk, yk, Hy, nsk, nyk )
    print *, R
    print *, ''
    print *, 'D (up)'
    print *, DI



end program
