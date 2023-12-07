! Module that defines double precision variables
!
! Template from https://fortran-lang.org/learn/rosetta_stone/
! In order to define double precision variables in 
! a separate module one includes 
!   use ldlqn_types, only: 
! at the beginning of the module
!
! To declare double precision variables use e.g., 
! real(dp) :: f  
module ldlqn_types
    implicit none
    private
    public dp, hp
    integer, parameter :: dp=kind(0.d0), &          ! double precision
                            hp=selected_real_kind(15) ! high precision
end module