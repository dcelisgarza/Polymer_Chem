!=================================================================================!
!----------------------------------CONFIG MODULE----------------------------------!
!=================================================================================!
module mod_config_polymer
  ! Module with types and basic procedures for the simulation of a polymerisation.
  use numbers
  implicit none

!=================================================================================!
!--------------------------------------TYPES--------------------------------------!
!=================================================================================!
  type monomer
    ! Definition of a monomer.
    character(:), allocatable :: name ! Monomer name.
    real(i16)                 :: amount  ! Amount of monomer.
    real(dp)                  :: mass    ! Monomer mass.
    real(dp), allocatable     :: k(:)    ! Kinetic reaction coefficients with other monomers.
    real(dp), allocatable     :: p(:)    ! Reaction probability. Let n := # of monomers, then extent(p) = n.
  end type monomer
!=================================================================================!
  type termination
    ! Definition of a termination.
    character(:), allocatable :: name ! Termination name.
    real(dp), allocatable     :: p(:)    ! Termination probability.
  end type termination
!=================================================================================!
  type chains
    integer(i16)              :: kl        ! Kinetic chain length.
    character(:), allocatable :: store(:)  ! Store chains.
    integer(i16)              :: index     ! Index of a chain.
    integer(i16)              :: rem       ! Number of chains removed.
    integer(i16), allocatable :: length(:) ! Chain lengths.
  end type chains

!=================================================================================!
!------------------------------------INTERFACE------------------------------------!
!=================================================================================!
  interface operator(+)
    ! Defining addition of monomers and terminations.
    module procedure sum_monomer
  end interface
!=================================================================================!
  interface operator(/)
    ! Defining division of monomers and terminations.
    module procedure normalise_termination
  end interface

!=================================================================================!
!-----------------------------------PROCEDURE-------------------------------------!
!=================================================================================!
contains
!=================================================================================!
!-----------------------------------FUNCTIONS-------------------------------------!
!=================================================================================!
  function sum_monomer(a,b)
    implicit none
    ! Summing monomers.
    type(monomer)             :: sum_monomer
    type(monomer), intent(in) :: a, b
    ! Adding masses.
    sum_monomer % mass = a % mass + b % mass
    ! Adding amounts.
    sum_monomer % amount = a % amount + b % amount
  end function sum_monomer
!=================================================================================!
  function normalise_termination(a,b)
    implicit none
    type(termination) :: normalise_termination
    type(termination), intent(in) :: a, b
    real(dp) :: norm
    norm = sum( b % p )
    normalise_termination % p = normalise_termination % p / norm
  end function normalise_termination
!=================================================================================!
end module mod_config_polymer
