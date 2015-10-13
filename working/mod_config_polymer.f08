module mod_config_polymer
  ! Module with types and basic procedures for the simulation of a polymerisation.
  use numbers
!=================================================================================!
!--------------------------------------TYPES--------------------------------------!
!=================================================================================!
  type monomer
    ! Definition of a monomer.
    character(:), allocatable :: name(:) ! Monomer name.
    real(i16)                 :: amount  ! Amount of monomer.
    real(dp)                  :: mass    ! Monomer mass.
    real(dp), allocatable     :: k_ij(:) ! Kinetic reaction coefficients.
    real(dp), allocatable     :: r_i(:)  ! Ratio of kinetic reaction coefficients. Must have half the extent of k_ij.
    real(dp), allocatable     :: p(:)    ! Probability of reaction between monomers. Let n := # of monomers, and i := initiators, then extent(p) = n**2 + i
  end type monomer
!=================================================================================!
  type termination
    ! Definition of a termination.
    character(:), allocatable :: name(:) ! Termination name.
    real(dp)                  :: p       ! Termination probability.
    integer(i16)              :: l       ! Kinetic chain length.
  end type termination
!=================================================================================!
!------------------------------------INTERFACE------------------------------------!
!=================================================================================!
  interface operator(+)
    ! Defining addition of monomers and terminations.
    module procedure sum_monomer, sum_termination
  end interface
!=================================================================================!
  interface operator(/)
    ! Defining division of monomers and terminations.
    module procedure normalise_monomer, normalise_termination
  end interface
!=================================================================================!
!-----------------------------------PROCEDURE-------------------------------------!
!=================================================================================!
contains
!=================================================================================!
!-----------------------------------FUNCTIONS-------------------------------------!
!=================================================================================!
  function sum_monomer(a,b)
    ! Summing monomers.
    type(monomer)             :: sum_mon
    type(monomer), intent(in) :: a, b
    ! Adding masses.
    sum_mon % mass = a % mass + b % mass
    ! Adding amounts.
    sum_mon % amount = a % amount + b % amount
  end function sum_monomer
!=================================================================================!
  function normalise_monomer(a,b)
    ! Normalising reaction coefficients.
    type(monomer)             :: norm_mon
    type(monomer), intent(in) :: a, b
    real(dp)                  :: norm
    ! Adding corresponding reaction coefficients.
    ! K_I = K_IA + K_IB, K_A = K_AA + K_AB, or K_B = K_BA + K_BB
    norm = sum( b % k_ij )
    ! Normalising reaction coefficients.
    norm_mon % k_ij = a % k_ij / norm
  end function normalise_monomer
!=================================================================================!
  function sum_termination(a,b)
    ! Summing termination probabilities.
    type(termination)             :: sum_term
    type(termination), intent(in) :: a, b
    sum_term % p = a % p + b % p
  end function sum_termination
!=================================================================================!
  function normalise_termination(a,b)
    type(termination)             :: norm_term
    type(termination), intent(in) :: a, b
    norm_term % p = a % p / b % p
  end function normalise_termination
!=================================================================================!
!----------------------------------SUBROUTINES------------------------------------!
!=================================================================================!
  subroutine p_mon(poly)
    implicit none
    type(monomer), intent(inout) :: poly
  end subroutine p_mon
!=================================================================================!
end module mod_config_polymer
