!=================================================================================!
!----------------------------------CONFIG MODULE----------------------------------!
!=================================================================================!
module mod_config_polymer
  ! Module with types and basic procedures for the simulation of a polymerisation.
  ! Shared Nomenclature:
  ! name := Entity name.
  ! p(:) := Event probability
  use numbers
  implicit none
!=================================================================================!
!--------------------------------------TYPES--------------------------------------!
!=================================================================================!
  type monomer
    ! Definition of a monomer.
    character(:), allocatable :: name
    integer              :: amount  ! Amount of monomer.
    real(dp)                  :: mass    ! Monomer mass.
    real(dp), allocatable     :: k(:)    ! Kinetic reaction coefficients with other monomers.
    real(dp), allocatable     :: p(:)    ! Let n := # of monomers (not initiators), then extent(p) = n.
    integer(i1)               :: reacted ! If reacted = 1, the monomer reacted last time, so a new reaction probability has to be calculated. If reacted = 0 the monomer did not react last time.
  end type monomer
!=================================================================================!
  type termination
    ! Definition of a termination.
    character(:), allocatable :: name
    integer, allocatable :: kl(:) ! Kinetic chain length. extent(p) = n.
    real(dp), allocatable     :: p(:)  ! Let n := # of terminations, then
  end type termination
!=================================================================================!
  type chain
    character(:), allocatable :: store(:)  ! Chain storage.
    integer, allocatable :: length(:) ! Chain lengths.
    integer              :: index     ! Index of a chain. Let n := # of chains, then index = n + 1. Chains will be stored in the range [1,index-1], this prepares the array for the next chain.
    integer              :: rem       ! Number of chains removed.
  end type chain

!=================================================================================!
!------------------------------------INTERFACE------------------------------------!
!=================================================================================!
!  interface operator(+)
!    module procedure sum_termination, sum_monomer
!  end interface

!  interface operator(/)
!    module procedure normalise_termination, normalise_monomer
!  end interface

!=================================================================================!
!-----------------------------------PROCEDURE-------------------------------------!
!=================================================================================!
!contains
!=================================================================================!
!-----------------------------------FUNCTIONS-------------------------------------!
!=================================================================================!
!  function sum_termination(a,b)
!    type(termination) :: sum_termination
!    type(termination), intent(in) :: a, b

!    sum_termination % p = a % p + b % p
!  end function sum_termination

!  pure function normalise_termination(a,b)
!    type(termination) :: normalise_termination
!    type(termination), intent(in) :: a, b

!    normalise_termination % p = a % p / b % p
!  end function normalise_termination

!  function sum_monomer(a,b)
!    type(monomer) :: sum_monomer
!    type(monomer), intent(in) :: a, b

!    sum_monomer % k   = a % k + b % k
!    sum_monomer % amount = a % amount + b % amount
!  end function sum_monomer

!  pure function normalise_monomer(a,b)
!    type(monomer) :: normalise_monomer
!    type(monomer), intent(in) :: a, b
!    real(dp) :: norm

!    norm = sum( abs( b % k ) )
!    normalise_monomer % k = a % k / norm
!  end function normalise_monomer
end module mod_config_polymer
