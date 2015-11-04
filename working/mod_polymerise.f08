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
    allocate( character :: dimer(1) % name, &
                           dimer(2) % name, &
                           dimer(3) % name )

    ! Allocating reaction coefficients (k) and probabilities (p) for all monomers.
    allocate( dimer(1) % k(n_mon), dimer(1) % p(n_mon), &
              dimer(2) % k(n_mon), dimer(2) % p(n_mon), &
              dimer(3) % k(n_mon), dimer(3) % p(n_mon) )

    ! Allocating termination names.
    allocate ( character :: term(1) % name, &
                            term(2) % name, &
                            term(3) % name )

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
    o_chain % index = 1
    o_chain % rem   = 0

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

  subroutine refresh_chain_storage(o_chain)
    implicit none

    type(chains), intent(inout) :: o_chain

    ! Reallocating memory.
    deallocate( o_chain % store, o_chain % length )
    allocate( character :: o_chain % store(1) )
    allocate( o_chain % length(1) )

    ! We must add an index if we want to refresh this for chain storage.
    o_chain % index = 1
    o_chain % rem = 0

  end subroutine refresh_chain_storage

  subroutine chain_grow(c_chain, name)
    implicit none

    ! Make sure the name has already been trimmed before it enters the subroutine.
    character(len=*), intent(in)             :: name
    character(:), allocatable, intent(inout) :: c_chain

    ! Concantenating the name onto the chain.
    c_chain = c_chain // name

  end subroutine chain_grow

  subroutine chain_terminate
  end subroutine chain_terminate

  subroutine chain_store(o_chain,c_chain)
    implicit none

    type(chains), intent(inout)               :: o_chain
    character(len=*), intent(in)              :: c_chain
    type(chains)                              :: work
    integer(i16)                              :: o_index, c_index, n_index, c_chain_length, i

    ! Current index
    c_index = o_chain % index
    ! Save the old index.
    o_index = c_index - 1
    ! The new index is one more than the old one.
    n_index = c_index + 1
    ! Obtain the current chain length.
    c_chain_length = len(c_chain)

    ! If this is not the first chain.
    if (c_index > 1) then
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
      o_chain % length(c_index) = c_chain_length
      o_chain % store(c_index)(1:c_chain_length) = c_chain
    else
      ! If this is the first entry of the array then there is nothing to save.
      ! So we only need to save the chain and its length.
      o_chain % length(c_index) = c_chain_length
      o_chain % store(c_index)(1:c_chain_length) = c_chain
    end if

    ! Update the chain index.
    o_chain % index = n_index

  end subroutine chain_store

  subroutine remove_chain(ol_chain, t_index)
    implicit none

    type(chains), intent(inout)              :: ol_chain ! Old chains, ol := lifted chain.
    integer(i16), intent(in)                 :: t_index ! Lifted index.
    ! Flagging the chain we just lifted as removed.
    ! Making the removed chain's length equal to zero.
    ol_chain % length(t_index) = 0
    ! If we don't do this Fortran places whatever character we set this to at the start of all other chains.
    ! By doing this we avoid such undesirable phenomenon. Works with (1:n) where n <= 0.
    ol_chain % store(t_index)(1:0) = ''
    ! Add a counter to the number of removed chains.
    ol_chain % rem = ol_chain % rem + 1

  end subroutine remove_chain

  subroutine transfer(ol_chain, ot_chain, c_chain, t_index)
    implicit none

    type(chains), intent(inout)              :: ol_chain, ot_chain ! Old chains, ol := lifted chain, os := old chains by transfer.
    character(:), allocatable, intent(inout) :: c_chain ! Current chain comes in, reactivated chain comes out.
    integer(i16), intent(in)                 :: t_index ! Lifted index.
    ! Storing the chain that just ended.
    call chain_store(ot_chain, c_chain)
    ! Swapping the current chain to the chain we just transfered the active site to.
    c_chain = ol_chain % store(t_index)(1:ol_chain % length(t_index))
    ! Removing chain.
    call remove_chain(ol_chain, t_index)

  end subroutine transfer

!  subroutine recombination(o_chain, c_chain, index)
!    implicit none
!    type(chains), intent(inout) :: os_chain
!  end subroutine recombination
end module polymerisation
