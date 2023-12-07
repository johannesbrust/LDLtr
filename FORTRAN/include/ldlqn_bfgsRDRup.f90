!>
!> ldlqn_bfgsRDRup 
!>
!>
!> Module for updating the inverse LDLT factorization of a 
!> BFGS quasi-Newton matrix. 
!>
!>

module ldlqn_bfgsRDRup

    use kind_parameter, only: dp

    implicit none

    contains

    !>
    !> Performing the update using an eigendecomposition of 
    !> a rank-2 update. This function is a translation of 
    !> the Matlab function bfgsRDRUpdate_E1.m 
    !>
    subroutine bfgsRDRupEig(R, DI, GI, sk, yk, Hy, nsk, nyk )
    ! % bfgsRDRUpdate_E Updating the inverse LDLT factorization of the BFGS matrix
    ! % and its inverse using eigenvalues of the rank-2 update
    ! %   bfgsRDRUpdate( R, DI, GI, sk, yk, Hy, nsk, nyk, sy )
    ! %
    ! %   Updates the inverse of a LDL' factorization by a rank-2 matrix
    ! %
    ! %   L1 D1 L1' = L D L' - Bs*(Bs/sk'*Bs)' + yk*(yk/sk'*yk)'
    ! %
    ! %   where R = inv(L1), and DI = inv(D1)
    ! %
    ! %   On output D and R are overwritten by DI and R1.
    ! %   The vector GI contains the column norms of R
    ! %
    ! %   INPUTS:
    ! %   R : Upper triangular
    ! %   DI: Diagonal (stored as vector)
    ! %   GI : Column lengths of R (stored as vector)
    ! %   sk : Vector for the BFGS update
    ! %   yk : Vector for the BFGS update
    ! %   Hy : Vector for the inverse BFGS update
    ! %   nsk : norm(sk) (scalar)
    ! %   nyk : norm(yk) (scalar)
    ! %   sy : sk'*yk (scalar)
    ! %
    ! %   OUTPUTS:
    ! %   R : Updated upper triangular
    ! %   DI: Updated diagonal (stored as vector)
    ! %   GI : Diagonals of R'*R
    ! %
    ! %--------------------------------------------------------------------------
    ! % 01/26/23, J.B., Initial implementation
    ! % 03/02/23, J.B., Inverse updating
    ! % 03/10/23, J.B., Eigenvalue update
    ! % 05/23/23, J.B., Storing of diagonals

    real(dp), intent(out)   :: R(:,:)
    real(dp), intent(out)   :: DI(:)
    real(dp), intent(out)   :: GI(:)
    real(dp), intent(in)    :: sk(:)
    real(dp), intent(in)    :: yk(:)
    real(dp), intent(in)    :: Hy(:)
    real(dp), intent(in)    :: nsk
    real(dp), intent(in)    :: nyk
    
    !> 
    !> Temporary/ internal variables
    !>
    real(dp), allocatable   :: w(:)
    real(dp)                :: nHy
    real(dp)                :: nsy
    real(dp)                :: l1
    real(dp)                :: l2
    real(dp)                :: v1
    real(dp)                :: v2
    real(dp)                :: djI
    real(dp)                :: alp
    real(dp)                :: bet
    real(dp)                :: buff

    !>
    !> For testing
    !>
    real(dp), allocatable   :: w1(:), w2(:)

    integer                 :: n
    integer                 :: j, jI

    n = size(sk)
    allocate(w(n))

    !>
    !> Computing eigenvalues and eigenvectors of a small (rank-2)
    !> eigendecomposition. Eigenvalues are l1, and l2, however
    !> these variables are also used for temporary storage
    !> 
    
    nHy     = norm2(Hy)
    w(:)    = yk/nyk
    nsy     = dot_product(sk/nsk,w)
    
    !>
    !> Intermediate variables
    !> 
    w1      = sk/nsk
    w2      = Hy/nHy


    l1      = (nsk*(nsy/nyk) + dot_product(w,(Hy/nyk)))/nsy**2;
    l2      = -(nHy)/(nsy*nyk)  
    v1      = 1.0
    v2      = l2                                ! temporary store l2
    l1      = 0.5*(l1+sqrt(l1**2 + 4*l2**2))
    l2      = -l2**2 / l1
    v2      = v2/l1                         
    buff    = sqrt(v1**2 + v2**2)               ! temporary variable  
    v1      = v1 / buff
    v2      = v2 / buff

    !>
    !> Update the working vector w and initialize
    !> alp, bet for the main loops
    !>

    !w(:)    = (v1/nsk)*sk + (v2/nHy)*Hy
    w(:)    = (v1)*w1 + (v2)*w2

    alp     = l1

    !>
    !> First update loop
    !>
    do j = 1, n

        jI  = n-j+1 

        djI = DI(jI) + alp*w(jI)**2
        bet = (w(jI)*alp)/djI
        alp = (DI(jI)*alp)/djI
        
        DI(jI) = djI

        w(1:(n-j))      = w(1:(n-j)) - w(jI)*R(1:(n-j),jI);
        R(1:(n-j),jI)   = R(1:(n-j),jI) + bet*w(1:(n-j));

    end do

    !>
    !> Update the working vector w and initialize
    !> alp, bet for the main loops
    !>

    !w(:)    = (-v2/nsk)*sk + (v1/nHy)*Hy
    w(:)    = (-v2)*w1 + (v1)*w2
    alp     = l2

    !>
    !> Second update loop
    !>
    do j = 1, n

        jI  = n-j+1 

        djI = DI(jI) + alp*w(jI)**2
        bet = (w(jI)*alp)/djI
        alp = (DI(jI)*alp)/djI
        
        DI(jI) = djI

        w(1:(n-j))      = w(1:(n-j)) - w(jI)*R(1:(n-j),jI);
        R(1:(n-j),jI)   = R(1:(n-j),jI) + bet*w(1:(n-j));

        GI(jI)  = dot_product(R(1:jI,jI),R(1:jI,jI))

    end do

    !>
    !> Make matrix positive definite, if round-off errors
    !> resulted in an indefinite matrix
    !>

    if (minval(DI) < 0) then
        where (DI < 0) 
            DI = abs(DI)
        end where
        !DI(DI < 0) = abs(DI(DI < 0))
    end if

    end subroutine



end module ldlqn_bfgsRDRup