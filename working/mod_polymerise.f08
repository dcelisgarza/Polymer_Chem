!=================================================================================!
!------------------------------POLYMERISATION MODULE------------------------------!
!=================================================================================!
module polymerisation
  use mod_config_polymer
  implicit none
  real(dp) :: ZBQLU01
  external ZBQLINI, ZBQLU01
contains
!=================================================================================!
  subroutine refresh_chain_storage(o_chain)
    ! Subroutine for restoring chain storage to its default values.
    implicit none
    type(chains), intent(inout) :: o_chain ! Chain storage.

    ! Reallocating memory.
    deallocate(o_chain % store, o_chain % length)
    allocate(character :: o_chain % store(1))
    allocate(o_chain % length(1))
    ! Refreshing indexes, removed chains and kinetic length.
    o_chain % index = 1
    o_chain % rem   = 0
  end subroutine refresh_chain_storage
!=================================================================================!
  subroutine refresh_chain(c_chain)
    ! Clear a chain.
    implicit none
    character(:), allocatable, intent(inout) :: c_chain ! Current chain.

    ! Deallocating c_chain.
    deallocate(c_chain)
    ! Allocating c_chain.
    allocate(character :: c_chain)
  end subroutine refresh_chain
!=================================================================================!
  subroutine term_prob_from_kinetic_chain_len(term)
    ! Calculate un-normalised termination probabilities from kinetic chain lengths.
    type(termination), intent(inout) :: term(:) ! Termination array.
    integer                          :: n_term, i  ! Number of terminations and counter variable.

    ! Number of terminations.
    n_term = size(term)
    ! Calculating un-normalised Termination Probabilites (TP).
    TP: do concurrent (i = 1: n_term)
      term(i) % p(1) = (term(i) % kl(1)) ** (-1.)
    end do TP
    ! Normalise terminaton probabilities.
    call norm_term_prob(term)
  end subroutine term_prob_from_kinetic_chain_len
!=================================================================================!
  subroutine norm_term_prob(term)
    ! Normalise termination probabilities.
    implicit none
    type(termination), intent(inout) :: term(:)      ! Termination array.
    real(dp), allocatable            :: norm_term(:) ! Normalisation constant.
    integer                          :: i, n_tot     ! Counter variables.

    ! Setting counter limit to the size of the termination array.
    n_tot = size(term)
    ! Allocating memory.
    allocate (norm_term(n_tot))
    ! Initialising normalisation constant.
    norm_term  = 0.
    ! Calculating the Normalisation Constant of the Terminations (NCT).
    NCT: do i = 1, n_tot
      norm_term = norm_term + term(i) % p
    end do NCT
    ! Normalised Terminations (NT).
    NT: do concurrent (i = 1: n_tot)
      term(i) % p = term(i) % p / norm_term
    end do NT
  end subroutine norm_term_prob
!=================================================================================!
  subroutine norm_rtn_coeff(mon)
    ! Normalise reaction coefficients.
    implicit none
    type(monomer), intent(inout) :: mon(:)        ! Monomer array.
    real(dp), allocatable        :: norm_coeff(:) ! Normalisation constant.
    integer                      :: i, j, n_tot   ! Counter variables.

    ! Setting counter limit to the size of the monomer array.
    n_tot = size(mon)
    ! Allocating memory.
    allocate(norm_coeff(n_tot))
    ! Initialising normalisation constant.
    norm_coeff = 0.
    ! Calculating the Normalisation Constant of the Reaction Coefficients (NCRC).
    ! We set i in K_ij.
    NCRC1: do i = 1, n_tot
      ! Then we go through all j's.
      NCRC2: do j = 1, n_tot
        norm_coeff(i) = norm_coeff(i) + mon(i) % k(j)
      end do NCRC2
    end do NCRC1
    ! Normalising Reaction Coefficients (NRC).
    NRC: do concurrent (i = 1: n_tot)
      mon(i) % k = mon(i) % k / norm_coeff(i)
    end do NRC
  end subroutine norm_rtn_coeff
!=================================================================================!
  pure subroutine rtn_prob(a_mon, o_mon, n_mon)
    ! Calculate reaction probabilities.
    implicit none
    type(monomer), intent(inout) :: a_mon(:)     ! Array of monomers.
    integer, intent(in)          :: o_mon, n_mon ! o_mon := monomer in a chain (or initiator). n_mon := free monomer.
    real(dp)                     :: num, denom   ! Numerator and denominator of the reaction probability.
    integer                      :: i, nbr_mon   ! Counter, number of monomers.

    ! Checking if a Reaction is Possible (RP)---only possible when k_ij > 0, or the amount of the free monomer > 0.
    RP: if ( a_mon(o_mon) % k(n_mon) == 0 .or. a_mon(n_mon) % amount == 0) then
      ! If reaction is impossible set the probability of reaction to zero
      a_mon(o_mon) % p(n_mon) = 0
      ! and skip to the end of the subroutine.
      return
    end if RP
    ! Setting counter limit to the size of the array of monomers.
    nbr_mon = size(a_mon)
    ! Preparing the numerator and denominator.
    denom = 0.
    num = 0.
    ! a_mon(o_mon) := Old monomer in a chain (or initiator)
    ! a_mon(n_mon) := Monomer which is potentially reacting with the old monomer.
    ! Calculating the numerator.
    ! Let "i" represent the monomer in a chain (or initiator) and "j" represent the free monomer:
    ! Reaction speed = k_ij * [n]
    num = a_mon(o_mon) % k(n_mon) * dble(a_mon(n_mon) % amount)
    ! Calculating the Denominator (D).
    ! Sum of all reaction speeds \sum_{j=1}^{N} ( k_ij * [j] )
    D: do i = 1, nbr_mon
      denom = denom + a_mon(o_mon) % k(i) * dble(a_mon(i) % amount)
    end do D
    ! Calculating the probability.
    a_mon(o_mon) % p(n_mon) = num / denom
end subroutine rtn_prob
!=================================================================================!
  subroutine prob_limits(limit, prob, nbr_limits)
    ! Calculate the limits of the ranges used in deciding what to do.
    implicit none
    real(dp), intent(out) :: limit(:)   ! Limits.
    real(dp), intent(in)  :: prob(:)    ! Array of probabilities.
    integer, intent(in)   :: nbr_limits ! Number of limits to calculate.
    integer               :: i

    ! Calculating Limits (CL). Our probabilities are normalised so they add up to one, no need to calculate the top limit as it's 1.
    CL: do concurrent (i = 1: nbr_limits)
      ICL: if (i > 1) then
        limit(i) = prob(i) + limit(i-1)
      else ICL
        limit(i) = prob(i)
      end if ICL
    end do CL
  end subroutine prob_limits
!=================================================================================!
  subroutine rtn_mon(mon, i_choice, limit)
    ! Choose monomer for reaction.
    implicit none
    type(monomer), intent(inout) :: mon(:)
    integer, intent(out) :: i_choice ! Choosing the monomer to react.
    real(dp), intent(in) :: limit(:)
  end subroutine rtn_mon
!=================================================================================!
  subroutine rtn(c_chain, mon, limit)
    implicit none
    character(:), allocatable, intent(inout) :: c_chain      ! Current chain.
    type(monomer), intent(inout)             :: mon(:)       ! Monomers.
    real(dp), allocatable, intent(out)       :: limit(:)
    character(:), allocatable                :: try_last_mon ! Trying to find the last monomer in c_chain.
    integer                                  :: i, n_mon, c_length, last_mon_idx, i_choice ! Counter variable, number of monomers, length of c_chain, index of the last monomer in mon(:), free monomer which will react
    real(dp), allocatable                    :: work(:)      ! Limit and work arrays for probabilities.

    n_mon = size(mon)
    ! Allocating work array.
    allocate( limit(n_mon), work(n_mon) )
    c_length = len(c_chain)
    ! Find the Last Monomer in the Chain (LMC).
    LMC: do i = 1, n_mon
      ! We're taking the last len( mon(i) % name ) characters in the chain, which are then checked again checked whether they are congruent with a given monomer's name.
      try_last_mon = c_chain(c_length - len( mon(i) % name ) + 1: c_length)
      ! Checking try_last_mon against mon(i) % name
      ILMC: if (try_last_mon == mon(i) % name) then
        ! If it's true then the last monomer is mon(i)
        last_mon_idx = i
        ! We exit the loop.
        exit LMC
      end if ILMC
    end do LMC
    ! We only need to calculate reaction probabilities between mon(last_mon_idx), or Last monomer in the chain and all Other monomers (LO).
    LO: do concurrent (i = 1: n_mon)
      call rtn_prob(mon, last_mon_idx, i)
    end do LO
    ! Passing reaction Probabilities onto the Work (PW) array.
    PW: do concurrent (i = 1: n_mon)
      work(i) = mon(last_mon_idx) % p(i)
    end do PW
    ! Calculate the limits of our ranges to help decide which monomer will react..
    call prob_limits(limit, work, n_mon)
  end subroutine rtn
!=================================================================================!
  subroutine grow_chain(c_chain, mon)
    implicit none
    ! Make sure the name has already been trimmed before it enters the subroutine.
    character(:), allocatable, intent(inout) :: c_chain ! Current chain.
    character(len=*), intent(in)             :: mon     ! Monomer to be added to the chain.

    ! Concantenating the name onto the chain.
    c_chain = c_chain // mon
  end subroutine grow_chain
!=================================================================================!
  subroutine reverse_chain(c_chain)
    implicit none
    ! Make sure c_chain has already been trimmed.
    character(len=*), intent(inout) :: c_chain     ! Current chain.
    integer(i16)                    :: i, c_length ! Counter variables.

    ! Setting the length to be equal to the chain length.
    c_length = len(c_chain)
    ! Reversing characters.
    forall (i=1:c_length) c_chain(i:i) = c_chain(c_length-i+1:c_length-i+1)
  end subroutine reverse_chain
!=================================================================================!
  subroutine store_chain(o_chain, c_chain)
    implicit none
    type(chains), intent(inout)  :: o_chain ! Old chains.
    character(len=*), intent(in) :: c_chain ! Current chain.
    type(chains)                 :: work    ! Work array for moving data.
    integer(i16)                 :: o_index, c_index, n_index, & ! Old, current and next index.
                                    c_chain_length, & ! Current chain length.
                                    i ! Counter variable.

    ! Current index
    c_index = o_chain % index
    ! Save the old index.
    o_index = c_index - 1
    ! The new index is one more than the old one.
    n_index = c_index + 1
    ! Obtain the current chain length.
    c_chain_length = len(c_chain)
    ! If this is not the First Chain (FC).
    FC: if (c_index > 1) then
      ! Allocate the array where the data will be saved.
      allocate( character :: work % store(o_index) )
      allocate( work % length(o_index) )
      ! Save all the Old Chains (SOC) and their lengths dummy arrays.
      SOC: do concurrent (i = 1: o_index)
        work % length(i) = o_chain % length(i)
        work % store(i)(1:work % length(i)) = o_chain % store(i)(1:o_chain % length(i))
      end do SOC
      ! Deallocate the old chain and length storage and allocate space for one more entry (n_index = o_index + 1) in each.
      deallocate( o_chain % store, o_chain % length )
      allocate( character :: o_chain % store(n_index) )
      allocate( o_chain % length(n_index) )
      ! Move the Saved Chains Back (MSCB) into the expanded storage array.
      MSCB: do concurrent (i = 1: o_index)
        o_chain % length(i) = work % length(i)
        o_chain % store(i)(1:o_chain % length(i)) = work % store(i)(1:work % length(i))
      end do MSCB
      ! Add the new entries to the expanded length and chain storage arrays.
      o_chain % length(c_index) = c_chain_length
      o_chain % store(c_index)(1:c_chain_length) = c_chain
    ! If this is the first entry of the array then there are no old chains to save.
    else FC
      ! We only need to save the current chain and its length.
      o_chain % length(c_index) = c_chain_length
      o_chain % store(c_index)(1:c_chain_length) = c_chain
    end if FC
    ! Update the chain index.
    o_chain % index = n_index
  end subroutine store_chain
!=================================================================================!
  subroutine remove_chain(o_chain, r_index)
    implicit none
    type(chains), intent(inout) :: o_chain             ! Old chains.
    integer(i16), intent(in)    :: r_index             ! Index to be removed.
    integer(i16)                :: o_index, n_index, i ! Old and new index.
    type(chains)                :: work                ! Work array for moving data.

    ! Calculating indices.
    o_index = o_chain % index
    ! Removing one chain so the index is one less than before.
    n_index = o_index - 1
    ! Allocating a work array for data transfer.
    allocate( character :: work % store(n_index) )
    ! Removing the old chain length, dynamically reallocating space. This shifts all indexes starting from r_index down by one.
    o_chain % length(r_index) = 0
    o_chain % length = pack(o_chain % length, o_chain % length /= 0)
    ! Saving the Remaining Chains (SRC).
    SRC: do concurrent (i = 1: n_index)
      ! Moving data has to be done differently depending on whether or not we're Behind the Removed Index (BRI).
      BRI: if (i < r_index) then
        work % store(i)(1:o_chain % length(i)) = o_chain % store(i)(1:o_chain % length(i))
      else BRI
        ! We shift i to i+1 in o_chain % store(i+1), because we're skipping the chain we're removing from the old array.
        ! We're still using i in o_chain % length(i) because we have previously removed the r_index'th entry with pack().
        work % store(i)(1:o_chain % length(i)) = o_chain % store(i+1)(1:o_chain % length(i))
      end if BRI
    end do SRC
    ! Deallocating the old chain storage.
    deallocate( o_chain % store )
    allocate( character :: o_chain % store(n_index)  )
    ! Moving the Saved Chains Back (MSCB) into the contracted storage array.
    MSCB: do concurrent (i = 1: n_index)
      o_chain % store(i)(1:o_chain % length(i)) = work % store(i)(1:o_chain % length(i))
    end do MSCB
    ! Setting the index to the new index.
    o_chain % index = n_index
    ! Registering the fact that we have removed one chain.
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
    call reverse_chain(c_chain)
    ! Grabbing the old chain, onto which we will append the reversed c_chain.
    work = ol_chain % store(r_index)(1:ol_chain % length(r_index))
    ! Removing the lifted chain from storage.
    call remove_chain(ol_chain, r_index)
    ! Appending the reversed chain onto the lifted work chain.
    c_chain = work // c_chain
    ! Storing recombined chain into the recombination array.
    call store_chain(or_chain, c_chain)
  end subroutine recombination
!=================================================================================!
  subroutine transfer(ol_chain, ot_chain, c_chain, t_index)
    implicit none
    type(chains), intent(inout)              :: ol_chain, ot_chain ! Old chains, ol := lifted chain, ot := old chains by transfer.
    character(:), allocatable, intent(inout) :: c_chain ! Current chain comes in, reactivated chain comes out.
    integer(i16), intent(in)                 :: t_index ! Index of the chain reactivated by transfer.

    ! Storing the chain that just ended.
    call store_chain(ot_chain, c_chain)
    ! Swapping the current chain to the chain we just transfered the active site to.
    c_chain = ol_chain % store(t_index)(1:ol_chain % length(t_index))
    ! Removing chain reactivated by transfer.
    call remove_chain(ol_chain, t_index)
  end subroutine transfer
!=================================================================================!
subroutine chain_availability(o_chain, c_chain, term_flag)
  implicit none
  type(chains), intent(inout)              :: o_chain(:) ! Old chain.
  character(:), allocatable, intent(inout) :: c_chain    ! Current chain.
  integer, intent(in)                      :: term_flag  ! Termination flag. term_flag = 2 -> recombination. term_flag = 3 -> transfer
  integer                                  :: i_choice   ! Choose the old chain array we perturb.
  integer(i16)                             :: index      ! Index of the chain we take.

  ! Check if there are Chains Available for Recombination (CAR).
  CAR: if (o_chain(1) % index > 1 .or. o_chain(3) % index > 1) then
    ! Narrow down the Choice for Lifted Chains (CLC). At least one of the following if conditions will be true.
    ! If there are no chains by disproportiation.
    CLC: if (o_chain(1) % index < 2) then
      ! Set choice to the only available alternative: transfer.
      i_choice = 3
    ! If there are no chains by transfer.
    else if (o_chain(3) % index < 2) then CLC
      ! Set choice to the only available alternative: disproportiation.
      i_choice = 1
    ! Otherwise randomise the choice.
    else CLC
      ! Random I_Choice (RIC)
      RIC: if (ZBQLU01(0) <= 0.5) then
        i_choice = 1
      else RIC
        i_choice = 3
      end if RIC
    end if CLC
    ! Deciding the lifted index.
    ! Remember that the last index of o_chain is only there in preparation for next step, so we must sample only the ones with stored chains, hence (o_chain(i_choice) % index - 1) as opposed to (o_chain(i_choice) % index).
    index = floor((o_chain(i_choice) % index - 1) * ZBQLU01(0)) + 1
    ! Decide whether Recombine Or Transfer (ROT)
    ROT: if (term_flag == 2) then
      ! Call the recombination subroutine with appropriate parameters.
      call recombination(o_chain(i_choice), o_chain(term_flag), c_chain, index)
    else if (term_flag == 3) then ROT
      call transfer(o_chain(i_choice), o_chain(term_flag), c_chain, index)
    end if ROT
  ! If there are no Chains Available for Recombination, end differently.
  else CAR
    ! If there is are no chains in either the disproportiation or transfer arrays, then the chain can only end by disproportiation.
    call store_chain(o_chain(1), c_chain)
  end if CAR
end subroutine chain_availability
!=================================================================================!
  subroutine decide_termination(o_chain, c_chain, term)
    type(chains), intent(inout)              :: o_chain(:) ! Old chains.
    character(:), allocatable, intent(inout) :: c_chain    ! Current chain.
    type(termination), intent(in)            :: term(:)
    integer                                  :: n_term, i, i_choice, term_flag  ! Number of terminations - 1, counter variable, integer choice and termination flag.
    real(dp)                                 :: choice ! Choose termination,
    real(dp), allocatable                    :: limit(:), work(:) ! Limits array, work array.

    ! Allocating size of limit to the number of terminations, this is unecessary but made for future-proofing the code should other terminations be added/removed. Makes updating the code easier.
    n_term = size(term) - 1
    ! Allocate array size.
    allocate( limit(n_term), work(n_term) )
    ! Moving the termination probabilities to a work array.
    do concurrent (i = 1: n_term)
      work(i) = term(i) % p(1)
    end do
    ! Calculating limits.
    call prob_limits(limit, work, n_term)
    ! Decide which Termination (DT) to use.
    DT: do
      ! Selecting which Termination (T) to use.
      choice = ZBQLU01(0)
      ! Disproportiation. Check if choice is in (0, P(disp)].
      T: if (choice <= limit(1)) then
        ! Store chain by disproportiation.
        call store_chain(o_chain(1), c_chain)
      ! Recombination. Check if choice is in (P(disp), P(disp) + P(reco)]
      else if (limit(1) < choice .and. choice <= limit(2)) then T
        ! Check if there are Chains Available for Recombination (CAR).
        term_flag = 2
        call chain_availability(o_chain, c_chain, term_flag)
      ! Transfer. Check if choice is in (P(disp) + P(reco), P(disp) + P(reco) + P(trns) = 1 ]
      else if (limit(2) < choice) then T
        term_flag = 3
        call chain_availability(o_chain, c_chain, term_flag)
      end if T
    end do DT
  end subroutine decide_termination
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
    allocate( dimer(1) % k(n_tot), dimer(1) % p(n_tot), &
              dimer(2) % k(n_tot), dimer(2) % p(n_tot), &
              dimer(3) % k(n_tot), dimer(3) % p(n_tot) )
    ! Allocating termination names.
    allocate ( character :: term(1) % name, &
                            term(2) % name, &
                            term(3) % name )
    ! Allocating termination probabilities.
    allocate ( term(1) % p(1), term(1) % kl(1), &
               term(2) % p(1), term(2) % kl(1), &
               term(3) % p(1), term(3) % kl(1) )
    ! Allocating old chains, chain storage, chain lengths and chain.
    allocate ( character :: o_chain(1) % store(1), &
                            o_chain(2) % store(1), &
                            o_chain(3) % store(1) )
    allocate ( o_chain(1) % length(1), &
               o_chain(2) % length(1), &
               o_chain(3) % length(1) )
    o_chain % index = 1
    o_chain % rem   = 0
    ! Allocating current chain
    allocate (character :: c_chain)
  end subroutine allocation
end module two_monomer_data_declaration
