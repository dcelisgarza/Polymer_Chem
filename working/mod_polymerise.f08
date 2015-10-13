module two_monomer_data
  ! Module for a two polymer copolymerisation.
  use mod_config_polymer
  ! n_mon := # of monomers, n_k_ij := # of reaction coefficients, n_term = # terminations
  integer(i1), parameter :: n_mon = 3, n_k_ij = n_mon - 1, n_term = 3, n_prob = ( n_mon - 1 )**2 + n_mon
  ! We'll be simulating a dimer.
  type(monomer)          :: dimer(n_mon)
  ! We're having three terminations.
  type(termination)      :: term(n_term), norm_term
end module two_monomer_data

module polymerisation
  use mod_config_polymer
  implicit none
contains
  subroutine two_monomer_copolymerisation
  end subroutine two_monomer_copolymerisation

  subroutine in_out_direct
  end subroutine in_out_direct

  subroutine in_out_txt
  end subroutine in_out_txt

  subroutine chain_initiate
  end subroutine chain_initiate

  subroutine chain_grow
  end subroutine chain_grow

  subroutine chain_terminate
  end subroutine chain_terminate

end module polymerisation
