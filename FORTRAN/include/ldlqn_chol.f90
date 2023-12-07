!> ldlqn_chol implements a column oriented cholesky factorization
!>
!> This implementation uses the full matrix A, and returns a full
!> upper triangular matrix R
!>
!>---------------------------------------------------------------
!> 01/09/23, J.B., Initial implementation

module ldlqn_chol
    use kind_parameter, only : dp

    contains

    !> cchol perform a column oriented cholesky factorization
    !>
    !>    | r11            |   | r11 r21 . rn1 |
    !>    | r21 r22        | * |     r22 . rn2 | = A
    !>    |  .  .   .      |   |         .  .  |
    !>    | rn1 rn2  . rnn |   |           rnn |
    !>  
    !>      R' * R = A
    !> 
    !> ARGUMENTS:
    !> A: Symmetric positive definite 2D Array (INPUT)
    !> R: Upper triangular (OUTPUT)
    !>
    !>-----------------------------------------------------------
    subroutine cchol(A,R)

        !>
        !> Varialbe Initializations
        !>

        real(dp), intent(in) :: A(:,:)
        real(dp), intent(out) :: R(:,:)
        
        real(dp)    :: rjj
        integer     :: n
        integer     :: j

        n = size(A,dim=1)

        !>
        !> Main loop of factorization
        !>

        do j = 1, n 

            R(j:n,j) = A(j:n,j)

            if (j > 1) then

                R(j:n,j) = R(j:n,j) - &
                    matmul(R(j:n,1:(j-1)),R(j,1:(j-1))) 

            end if

            rjj = sqrt(R(j,j))

            R(j:n,j) = R(j:n,j) / rjj 


        end do

    end 

end module 
