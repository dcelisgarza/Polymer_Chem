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
    allocate ( character :: o_chain(1) % store(1), &
                            o_chain(2) % store(1), &
                            o_chain(3) % store(1) )
    allocate ( o_chain(1) % length(1), &
               o_chain(2) % length(1), &
               o_chain(3) % length(1) )
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

  subroutine chain_grow(chain,name)
    implicit none
    ! Make sure the name has already been trimmed before it enters the subroutine.
    character(len=*), intent(in)             :: name
    integer(i4)                              :: name_len
    character(:), allocatable, intent(inout) :: chain
    character(:), allocatable                :: work
    name_len = len(name)
    allocate( character :: work )
    work = chain // name
    deallocate( chain )
    allocate( character :: chain )
    chain = work
  end subroutine chain_grow

  subroutine chain_terminate
  end subroutine chain_terminate

  !subroutine chain_store(o_chain,c_chain)
    !implicit none
    !character(:), allocatable, intent(inout)  :: o_chain(:) ! Store chains.
    !character(len=*), allocatable, intent(in) :: c_chain
    !type(chains), allocatable                 :: work
    !integer(i16)                              :: o_index, n_index
    !o_index = o_chain % index
    !n_index = o_index + 1_i16
    !allocate ( character :: work % store(o_index) )
    !work % store = o_chain % store
  !end subroutine chain_store

  subroutine chain_store(o_chain,c_chain)
    implicit none
    type(chains), allocatable, intent(inout)  :: o_chain
    character(len=*), allocatable, intent(in) :: c_chain
    type(chains), allocatable                 :: work
    integer(i16)                              :: o_index, n_index, i, c_chain_length

    ! Save the old index.
    o_index = o_chain % index
    ! The new index is one more than the old one.
    n_index = o_index + 1
    ! Obtain the current chain length.
    c_chain_length = len(c_chain)

    ! If this is not the first chain.
    if (o_index > 1) then
      ! Allocate the array where the data will be saved.
      allocate( character :: work % store(o_index) )
      allocate( work % length(o_index) )

      ! Save all the old chains and their lengths dummy arrays.
      do concurrent (i = 1: o_index)
        work % length(i) = o_chain % length(i)
        work % store(i)(1:work % length(i)) = o_chain % store(i)(1:o_chain % length(i))
      end do

      ! Deallocate the old chain and length storage and allocate space for one more entry (n_index = o_index + 1) in each.
      deallocate( o_chain % store, o_chain % length )
      allocate( character :: o_chain % store(n_index) )
      allocate( o_chain % length(n_index) )

      ! Move the saved chains and their lengths into the expanded storage array.
      do concurrent (i = 1: o_index)
        o_chain % length(i) = work % length(i)
        o_chain % store(i)(1:o_chain % length(i)) = work % store(i)(1:work % length(i))
      end do

      ! Add the new entries to the expanded length and chain storage arrays.
      o_chain % length(n_index) = c_chain_length
      o_chain % store(n_index)(1:c_chain_length) = c_chain
    else
      ! If this is the first entry of the array then there is nothing to save.
      ! So we only need to save the chain and its length.
      o_chain % length(o_index) = c_chain_length
      o_chain % store(o_index)(1:c_chain_length) = c_chain
    end if

    ! Update the chain index.
    o_chain % index = n_index
  end subroutine chain_store
end module polymerisation
