!>
!> ldlqn_pars
!>
!> Parameters for the ldlqn algorithm
!> This file is a translation of bfgs_parms.m. Some options from the  
!> original file are not implemented 
!>
!>-------------------------------------------------------------------------
!> 09/11/23, J.B, initial version

module ldlqn_pars

    use kind_parameter, only:dp

    implicit none

    !integer    :: maxIterations    = 1500        ! maximum iterates allowed
    !integer    :: maxIterations    = 3000        ! maximum iterates allowed
    integer     :: maxIterations    = 6000        ! maximum iterates allowed

    real(dp)    :: tiny             = 1.0e-12     ! tiny number
    real(dp)    :: tolStat          = 1.0e-4      ! stationary tolerance
    
    real(dp)    :: wolfeTols1       = 0.9         ! Wolfe gradient tolerance for quasi-Newton
    real(dp)    :: wolfeTols2       = 1.0e-4      ! Wolfe function tolerance
    integer     :: jfmax            = 20          ! max functions per line search
    real(dp)    :: gammaC           = 0.5         ! contraction factor in backtracking
    real(dp)    :: condMax          = 1.0e+16     ! maximum condition estimate
    real(dp)    :: dxMax            = 100         ! maximum  change in x    
    real(dp)    :: dxInf            = 1.0e+10     ! infinite change in x
    
    real(dp)    :: updateTol        = 0.99        ! Tolerance for rejecting BFGS update    
    integer     :: lines            = 15          ! lines between printing header    
    real(dp)    :: infinity         = 1.0e+15     ! user-defined +infinity    
    integer     :: count_near       = 1           ! count near optimal solutions
    integer     :: print_run_det    = 1           ! option to print run details
    
    real(dp)    :: unBoundedf       = -1.0e+9     ! definition of an unbounded objective value

    real(dp)    :: stp              = 1.0_dp
    real(dp)    :: gtol             = 0.9_dp
    real(dp)    :: ftol             = 1.0e-4_dp
    real(dp)    :: xtol             = 1e-12_dp
    real(dp)    :: stpmin           = 1e-14_dp
    real(dp)    :: stpmax           = 20_dp
    integer     :: maxfev           = 200 
    
    !>
    !> Trust-region parameters
    !>
    real(dp)    :: c1               = 1e-5        ! step acceptance
    real(dp)    :: c2               = 0.75       
    real(dp)    :: c3               = 0.8!0.8        
    real(dp)    :: c4               = 2           ! radius increase
    real(dp)    :: c5               = 0.1        
    real(dp)    :: c6               = 0.75       
    real(dp)    :: c7               = 0.5         ! radius decrease
    real(dp)    :: TRrelTol         = 1e-6        ! 1e-6 trust-region relative tolerance.
    real(dp)    :: TRabsTol         = 1e-3           ! trust-region absolute tolerance.
    integer     :: TRlog            = 0           ! trust-region log off/on
    integer     :: TRitnMax         = 30          ! trust-region iteration limit

    real(dp)    :: trTol            = 1e-7 ! 1e-7 ! tolerance for the TR subproblem
    integer     :: trMax            = 30          ! maximum iterations for the TR subproblem

    !>
    !> Inverse conjugate gradient parameters
    !>
    real(dp)    :: tolICG           = 1e-7        ! CG iteration complexity
    integer     :: maxitICG         = 35          ! Maximum CG iterations
    real(dp)    :: shrink           = 0.5         ! Factor to shrink the shift 
    integer     :: maxcg            = 5           ! Maximum iterations for Phase 2 (prev. default = 3)

    !>
    !> Problem size
    !> 
    integer     :: nmax             = 100         ! Size of MS algorithm
    
    !>
    !> Not-implemented parameters 
    !>

    !parms.write_file    = true       ! write final run statistics to file
    !parms.out_file      = false       ! write one line per iteration to file
    !parms.searchType    = 'Armijo'    ! Armijo           line search
    !parms.searchType    = 'Wolfe'     ! Wolfe            line search

    ! parms.updateTol     = 1.0e-8      ! Tolerance for rejecting BFGS update
    
    !parms.probsize          = 50000      ! Only run problems below n=probsize
    !parms.ResetFreq     = NaN         ! Reset frequency
    !parms.Assert        = 0           ! (0/1) check consistency
    !parms.trials        = 1           ! Repeat each problem # times metrics averaged
    !parms.show_warnings = false       ! (0/1) print warnings
    ! parms.dxMax         = 1.0e+10     ! maximum  change in x

end module