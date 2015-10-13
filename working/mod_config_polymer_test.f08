module mod_config_polymer
  ! Module with types and basic procedures for the simulation of a polymerisation.
  use numbers

  type monomer
    ! Definition of a monomer.
    character(:), allocatable :: name(:) ! Monomer name.
    real(i16)                 :: amount  ! Amount of monomer.
    real(dp)                  :: mass    ! Monomer mass.
    real(dp), allocatable     :: k_ij(:) ! Kinetic reaction coefficients.
  end type monomer

  type termination
    ! Definition of a termination.
    character(:), allocatable :: name(:) ! Termination name.
    real(dp)                  :: p       ! Termination probability.
    integer(i16)              :: l       ! Kinetic chain length.
  end type termination

  type terminations
    type(termination) :: disp, tran, reco
  end type terminations

  interface operator(+)
    module procedure sum_termination, sum_monomer
  end interface

  interface operator(/)
    module procedure normalise_termination, normalise_monomer
  end interface

contains

  function sum_termination(a,b)
    type(termination) :: sum_term
    type(termination), intent(in) :: a, b

    sum_term % p = a % p + b % p
  end function sum_termination

  function normalise_termination(a,b)
    type(termination) :: norm_term
    type(termination), intent(in) :: a, b

    norm_term % p = a % p / b % p
  end function normalise_termination

  function sum_monomer(a,b)
    type(monomer) :: sum_mon
    type(monomer), intent(in) :: a, b

    sum_mon % k_ij   = a % k_ij + b % k_ij
    sum_mon % amount = a % amount + b % amount
  end function sum_monomer

  function normalise_monomer(a,b)
    type(monomer) :: norm_mon
    type(monomer), intent(in) :: a, b
    real(dp) :: norm

    norm = sum( abs( b % k_ij ) )
    norm_mon % k_ij = a % k_ij / norm
  end function normalise_monomer

end module mod_config_polymer
