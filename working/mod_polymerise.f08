!=================================================================================!
!--------------------------TWO MONOMER DATA DECLARATION---------------------------!
!=================================================================================!
module two_monomer_data_declaration
  ! Module for a two polymer copolymerisation.
  use mod_config_polymer
  implicit none
  ! n_mon := # of monomers, n_init = # of initiators, n_tot := # of monomers + # of initiators, n_k := # of reaction coefficients, n_ki = # of reaction coefficients for a monomer with an initiator, n_term = # terminations.
  integer(i1), parameter    :: n_mon = 2, n_init = 1, n_tot = n_mon + n_init, &
                               n_k = n_mon, n_ki = n_mon*n_init, n_prob = ( n_mon )**2 + n_init, &
                               n_term = 3
  ! We'll be simulating a dimer.
  type(monomer)             :: dimer(n_tot)
  ! We're having three terminations.
  type(termination)         :: term(n_term)
  type(chains)              :: o_chain(n_term)
  character(:), allocatable :: c_chain
contains
  subroutine allocation
    implicit none
    ! Allocating monomer names.
    allocate( character :: dimer(1) % name(1), &
                           dimer(2) % name(1), &
                           dimer(3) % name(1) )
    ! Allocating reaction coefficients (k) and probabilities (p) for all monomers.
    allocate( dimer(1) % k(n_mon), dimer(1) % p(n_mon), &
              dimer(2) % k(n_mon), dimer(2) % p(n_mon), &
              dimer(3) % k(n_mon), dimer(3) % p(n_mon) )
    ! Allocating termination names.
    allocate ( character :: term(1) % name(1), &
                            term(2) % name(1), &
                            term(3) % name(1) )
    ! Allocating termination probabilities.
    allocate ( term(1) % p(1), &
               term(2) % p(1), &
               term(3) % p(1) )
    ! Allocating old chains, chain storage and chain lengths.
    allocate ( character :: o_chain(1) % chain(1), o_chain(1) % store(1), &
                            o_chain(2) % chain(1), o_chain(2) % store(1), &
                            o_chain(3) % chain(1), o_chain(3) % store(1) )
    allocate ( o_chain(1) % lengths(1), &
               o_chain(2) % lengths(1), &
               o_chain(3) % lengths(1) )
    ! Allocating current chain
    allocate ( character :: c_chain )
  end subroutine allocation
end module two_monomer_data_declaration
!=================================================================================!
!------------------------------POLYMERISATION MODULE------------------------------!
!=================================================================================!
module polymerisation
  use two_monomer_data_declaration
  implicit none
contains

  subroutine in_out
  end subroutine in_out

  subroutine chain_initiate
  end subroutine chain_initiate

  subroutine chain_grow(chain,mon,name_len)
    implicit none
    type(monomer), intent(in)                :: mon
    character(:), allocatable, intent(inout) :: chain
    !character(len = len(old) + 1_i2) :: new
    character(:), allocatable                :: work
    integer(i4)                              :: name_len
    !new = old // mon%name
    allocate(character( len=len( chain ) + len( mon % name ) ) :: work )
    work = chain // mon % name(1)(1:name_len)
    deallocate( chain )
    allocate( character( len=len( work ) ) :: chain )
    chain = work
  end subroutine chain_grow

  subroutine chain_terminate
  end subroutine chain_terminate

  subroutine chain_store(arr,chain)
    implicit none
    character(:), allocatable, intent(in)      :: chain
    character(len=*), allocatable, intent(out) :: arr(:)
  !  character(len=*), allocatable,
  end subroutine chain_store
end module polymerisation
