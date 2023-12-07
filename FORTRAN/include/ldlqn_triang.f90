! ldlqn_triang implements four operations with triangular matrices
!
! 1. Back substitution
! 2. Forward substitution
! 3. Mat-vec multiply with upper triangular
! 4. Mat-vec multiply with lower triangular
!
! -----------------------------------------------------------------
! 08/30/23, J.B., initial implementation
! 09/11/23, J.B., Update to make arguments inout in 
!                   matrix-vector multiplication

module ldlqn_triang
    !use ldlqn_types, only: dp
    !use, intrinsic :: iso_fortran_env, only: dp=>real64
    use kind_parameter, only: dp

    implicit none
    private
    public backsub, forwardsub, uppermult, lowermult

    contains

    ! Column-oriented back substitution (upper triangular)
    subroutine backsub(R,b,x)

        real(dp), intent(in) :: R(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) ::x(:)
        integer              :: n
        integer              :: j
        
        n = size(R,dim=2)        
        x(1:n) = b(1:n)

        do j = n, 2, -1

            x(j)        = x(j) / R(j,j)
            x(1:(j-1))  = x(1:(j-1)) - x(j) * R(1:(j-1),j) 

        end do

        x(1) = x(1)/ R(1,1)

    end

    ! Column-oriented forward substitution (lower triangular)
    subroutine forwardsub(R,b,x)

        real(dp), intent(in) :: R(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) ::x(:)
        integer              :: n
        integer              :: j
        
        n = size(R,dim=2)        
        x(1:n) = b(1:n)

        do j = 1, n-1

            x(j)        = x(j) / R(j,j)
            x((j+1):n)  = x((j+1):n) - x(j) * R((j+1):n,j) 

        end do

        x(n) = x(n)/ R(n,n)

    end

    ! Column-oriented upper-triangular times vector
    subroutine uppermult(R,x,b)

        real(dp), intent(in) :: R(:,:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(inout) ::b(:)
        integer              :: n
        integer              :: j
        
        n = size(R,dim=2)        
        b(1:n) = 0

        do j = n, 1, -1
                
            b(1:j)  = b(1:j) + x(j) * R(1:j,j) 

        end do

    end

    ! Column-oriented lower-triangular times vector
    subroutine lowermult(R,x,b)

        real(dp), intent(in) :: R(:,:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(inout) ::b(:)
        integer              :: n
        integer              :: j
        
        n = size(R,dim=2)        
        b(1:n) = 0

        do j = 1, n
            
            b(j:n)  = b(j:n) + x(j) * R(j:n,j)
            
            !print *, b

        end do

    end


end module