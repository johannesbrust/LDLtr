!>
!> ldlqn_invLDLCG.f90
!>
!> Translation of a Matlab inverse LDLT CG iteration 
!>
!>--------------------------------------------------------------
!> 09/11/23, J.B., Initial implementation

module ldlqn_invLDLCG

    use kind_parameter, only : dp
    use ldlqn_triang

    implicit none

    contains

    !>
    !> Inverse iteration using CG
    !> (Comment from the original Matlab implementation:
    !>      inv_ldl_CG1.m)
    ! %inv_ldl_CG1 Approximate solve with a LDL factorization
    ! %
    ! % The system to solve is
    ! %
    ! %       (sig.inv(L)inv(L')+D)x = g
    ! %
    ! %-----------------------------------------------------------------%
    ! % 05/09/2023, J.B., Initial version
    
    subroutine invLDLCG(x, p, res, sig, g, R, D, tol, maxit )
        
        real(dp), intent(out)   :: x(:)
        real(dp), intent(inout) :: p(:)
        real(dp), intent(inout) :: res(:)
        real(dp), intent(in)    :: sig
        real(dp), intent(in)    :: g(:)
        real(dp), intent(in)    :: R(:,:)
        real(dp), intent(in)    :: D(:)
        real(dp), intent(in)    :: tol
        integer, intent(in)     :: maxit

        !> 
        !> Further local variables
        !>
        
        integer                 :: n, icg
        real(dp)                :: nr0
        real(dp)                :: rr
        real(dp)                :: alpc
        real(dp)                :: nr
        real(dp)                :: betc
        real(dp), allocatable   :: ap(:)
        real(dp), allocatable   :: ap0(:)

        !>
        !> Starting assignment
        !>
        
        n = size(x)
        allocate(ap0(n))
        allocate(ap(n))

        x(:)    = 0
        res(:)  = g(:)
        p(:)    = res(:)
        nr0     = norm2(res)

        !>
        !> Main CG loop
        !>
        ! print '(a7, a7)', 'iter   ', 'nr     '
        ! print '(a14)', '--------------'
        ! print '(i5, f6.3)', 0, nr0

        do icg = 1, maxit

            !>
            !> Upper triangular matrix-vector product
            !>

            call uppermult(R,p,ap0)
            call lowermult(transpose(R),ap0,ap)

            ap      = sig * ap + D*p
            rr      = dot_product(res,res)
            alpc    = rr / dot_product(p,ap)
            
            x       = x + alpc * p
            res     = res - alpc * ap
            nr      = norm2(res)

            !>
            !> For testing. Print-out of residuals
            !>
            ! print '(i5, f6.3)', icg, nr

            if (nr/nr0 < tol ) then
                exit
            end if
            
            betc    = dot_product(res,res) / rr
            p       = res + betc * p

        end do

    end subroutine

end module