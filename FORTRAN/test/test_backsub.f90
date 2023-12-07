! Test the triangular library
!
! Test program for testing the implementations of 
! triangular matrices
!
! ------------------------------------------------------------
! 08/30/23, J.B., initial version

program test_backsub
    !use, intrinsic :: iso_fortran_env, only: dp=>real64
    use kind_parameter, only: dp
    use ldlqn_triang

    ! Initialization
    integer     :: n
    real(dp)    :: b(3)
    real(dp)    :: R(3,3)
    real(dp)    :: x(3)

    n           = 3            
    b = [1_dp, 1_dp, 1_dp]    
    R = reshape([1,0,0,2,3,0,1,2,3],shape(R))

    print *, 'Test back subsitution'
    call backsub(R,b,x)
    print *, x
    print *, ''

    print *, 'Test forward substitution'
    call forwardsub(transpose(R),b,x)
    print *, x
    print *, ''

    print *, 'Test upper triang. mult'
    call uppermult(R,b,x)
    print *, x
    print *, ''

    print *, 'Test lower triang. mult'
    call lowermult(transpose(R),b,x)
    print *, x
    print *, ''



end program