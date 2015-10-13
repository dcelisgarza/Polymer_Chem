!=================================================================================!
!--------------------------TWO MONOMER DATA DECLARATION---------------------------!
!=================================================================================!
module two_monomer_data_declaration
  ! Module for a two polymer copolymerisation.
  use mod_config_polymer
  ! n_mon := # of monomers, n_init = # of initiators, n_tot := # of monomers + # of initiators, n_k := # of reaction coefficients, n_ki = # of reaction coefficients for a monomer with an initiator, n_term = # terminations.
  integer(i1), parameter :: n_mon = 2, n_init = 1, n_tot = n_mon + n_init, &
                            n_k = n_mon, n_ki = n_mon*n_init, n_prob = ( n_mon )**2 + n_init, &
                            n_term = 3
  ! We'll be simulating a dimer.
  type(monomer)          :: dimer(n_tot)
  ! We're having three terminations.
  type(termination)      :: term(n_term)
contains
  subroutine allocation
    ! Allocating monomer names.
    allocate ( character(len=1) :: dimer(1) % name, dimer(2) % name, dimer(3) % name )
    ! Allocating reaction coefficients for all monomers.
    allocate ( dimer(1) % k(n_ki), dimer(2) % k(n_k), dimer(3) % k(n_k) )
    ! Allocating reaction probabilities for all monomers.
    allocate ( dimer(1) % p(n_ki), dimer(2) % p(n_k), dimer(3) % p(n_k) )
    ! Allocating termination names.
    allocate ( character(len=18) :: term(1) % name(1) )
    allocate ( character(len=8)  :: term(2) % name(1) )
    allocate ( character(len=13) :: term(3) % name(1) )
  end subroutine allocation
end module two_monomer_data_declaration
!=================================================================================!
!------------------------------POLYMERISATION MODULE------------------------------!
!=================================================================================!
module polymerisation
  use mod_config_polymer
  use two_monomer_data_declaration
  implicit none
contains

  subroutine in_out
  end subroutine in_out

  subroutine chain_initiate
  end subroutine chain_initiate

  subroutine chain_grow
  end subroutine chain_grow

  subroutine chain_terminate
  end subroutine chain_terminate
end module polymerisation
