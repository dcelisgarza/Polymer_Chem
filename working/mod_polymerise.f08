!=================================================================================!
!------------------------------POLYMERISATION MODULE------------------------------!
!=================================================================================!
module polymerisation
  use mod_config_polymer
  implicit none
  real(kind(1.0d0)) :: ZBQLU01
  external ZBQLINI, ZBQLU01
contains
!=================================================================================!
  subroutine refresh_chain_storage(o_chain)
    ! Subroutine for refreshing chain storage.
    implicit none
    type(chains), intent(inout) :: o_chain

    ! Reallocating memory.
    deallocate( o_chain % store, o_chain % length )
    allocate( character :: o_chain % store(1) )
    allocate( o_chain % length(1) )
    ! Refreshing indexes, removed chains and kinetic length.
    o_chain % index = 1
    o_chain % rem = 0
    o_chain % kl = 0

  end subroutine refresh_chain_storage
!=================================================================================!
  subroutine refresh_chain(c_chain)
    implicit none
    character(:), allocatable, intent(inout) :: c_chain

    ! Deallocating c_chain.
    deallocate(c_chain)
    ! Allocating c_chain.
    allocate(character :: c_chain)
  end subroutine refresh_chain
!=================================================================================!
  subroutine chain_grow(left, right)
    implicit none
    ! Make sure the name has already been trimmed before it enters the subroutine.
    character(len=*), intent(in)             :: right
    character(:), allocatable, intent(inout) :: left

    ! Concantenating the name onto the chain.
    left = left // right
  end subroutine chain_grow
!=================================================================================!
  subroutine chain_terminate
  end subroutine chain_terminate
!=================================================================================!
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
!=================================================================================!
  subroutine chain_reverse(c_chain)
    implicit none
    ! Make sure c_chain has already been trimmed.
    character(len=*), intent(inout) :: c_chain ! Current chain.
    integer(i4)                     :: i, c_length

    c_length = len(c_chain)
    forall (i=1:c_length) c_chain(i:i) = c_chain(c_length-i+1:c_length-i+1)
  end subroutine chain_reverse
!=================================================================================!
  subroutine remove_chain(o_chain, r_index)
    implicit none
    type(chains), intent(inout)              :: o_chain ! Old chains, ol := lifted chain.
    integer(i16), intent(in)                 :: r_index
    integer(i16)                             :: o_index, n_index, i ! Lifted index.
    type(chains)                             :: work

    ! Calculating indices.
    o_index = o_chain % index
    ! Removing one chain so the index is one less.
    n_index = o_index - 1
    ! Allocating a work array for data transfer.
    allocate( character :: work % store(n_index) )
    ! Removing the old chain length, dynamically reallocating space. This shifts all indexes after r_index down by one.
    o_chain % length(r_index) = 0
    o_chain % length = pack(o_chain % length, o_chain % length /= 0)
    ! Moving the non-removed chains into the work array.
    do concurrent (i = 1: n_index)
      if (i < r_index) then
        work % store(i)(1:o_chain % length(i)) = o_chain % store(i)(1:o_chain % length(i))
      else
        work % store(i)(1:o_chain % length(i)) = o_chain % store(i+1)(1:o_chain % length(i))
      end if
    end do
    ! Deallocating the old chain storage.
    deallocate( o_chain % store )
    allocate( character :: o_chain % store(n_index)  )
    ! Moving the saved chains back.
    do concurrent (i = 1: n_index)
      o_chain % store(i)(1:o_chain % length(i)) = work % store(i)(1:o_chain % length(i))
    end do
    o_chain % index = n_index

    ! Flagging the chain we just lifted as removed.
    ! Making the removed chain's length equal to zero.
    ! If we don't do this Fortran places whatever character we set this to at the start of all other chains.
    ! By making a zero size array we avoid such an undesirable phenomenon. Works with (a:b) where b < a.
    !o_chain % store(r_index)(1:0) = ''
    ! Add a counter to the number of removed chains.
    o_chain % rem = o_chain % rem + 1
  end subroutine remove_chain
!=================================================================================!
  subroutine recombination(ol_chain, or_chain, c_chain, r_index)
    implicit none
    type(chains), intent(inout)              :: ol_chain, or_chain ! Old chains, ol := lifted chain, or := old chains by recombination.
    character(:), allocatable, intent(inout) :: c_chain ! Current chain.
    integer(i16), intent(in)                 :: r_index ! Index of the chain onto which c_chain is added.
    character(:), allocatable                :: work    ! Chain which will be recombined with c_chain.

    ! Reversing the current chain (c_chain).
    call chain_reverse(c_chain)
    ! Grabbing the old chain, onto which we will append the reversed c_chain.
    work = ol_chain % store(r_index)(1:ol_chain % length(r_index))
    ! Removing the lifted chan.
    call remove_chain(ol_chain, r_index)
    ! Appending the reversed chain onto the lifted work chain.
    c_chain = work // c_chain
    ! Storing recombined chain into the recombination array.
    call chain_store(or_chain, c_chain)
  end subroutine recombination
!=================================================================================!
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
end module polymerisation
!=================================================================================!
!--------------------------TWO MONOMER DATA DECLARATION---------------------------!
!=================================================================================!
module two_monomer_data_declaration
  use polymerisation
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
