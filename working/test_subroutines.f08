program test_subroutines
  !use two_monomer_data_declaration
  use polymerisation
  implicit none
  character(:), allocatable :: test_chain
  integer(i1) :: counter
  real(dp) :: start, end
  ! For two-step allocation of character arrays (tsca) i.e. jagged arrays with character entries.
  type tscat
    character(:), allocatable :: tsc
  end type tscat
  type(tscat), allocatable :: tsca(:)
  ! For derived-type one-step allocation of character arrays (dsca).
  type oscat
    character(:), allocatable :: oscac(:)
  end type oscat
  type(oscat) :: dsca
  ! For one-step allocation of character arrays (osca).
  character(:), allocatable :: osca(:)
  integer(i16) :: remove_entry
  integer :: i
  character(1) :: str
  double precision, parameter :: dummy = 0.0
  real(dp) :: ZBQLU01
  external ZBQLINI,ZBQLU01

  call cpu_time(start)
  print*, '---------------------------------'
  print*, 'Testing two step allocation of character arrays.'
  ! Allocate array size.
  allocate(tsca(2))
  ! Allocate arbitrary character size to entries.
  allocate(character :: tsca(1)%tsc, tsca(2)%tsc)
  tsca(1)%tsc(1:4) = 'word'
  tsca(2)%tsc(1:9) = 'two words'
  print*, 'tsca(1) = ', tsca(1)%tsc(1:4)
  print*, 'tsca(2) = ', tsca(2)%tsc(1:9)
  print*, '---------------------------------'

  print*, ''
  print*, 'Testing one-step allocation of derived-type character arrays.'
  ! Simultaneaous allocation of derived-type array size and arbitrary character size to entries.
  allocate(character :: dsca % oscac(2))
  dsca % oscac(1)(1:4) = 'word'
  dsca % oscac(2)(1:9) = 'two words'
  print*, 'dsca % oscac(1) = ', dsca % oscac(1)(1:4)
  print*, 'dsca % oscac(2) = ', dsca % oscac(2)(1:9)
  print*, '---------------------------------'

  print*, ''
  print*, 'Testing one-step allocation of character arrays.'
  ! Simultaneaous allocation of array size and arbitrary character size to entries.
  allocate(character :: osca(2))
  osca(1)(1:4) = 'word'
  osca(2)(1:9) = 'two words'
  print*, 'osca(1) = ', osca(1)(1:4)
  print*, 'osca(2) = ',osca(2)(1:9)
  print*, '---------------------------------'

  print*, 'Testing allocation of parameters to a monomer.'
  call allocation
  dimer(1) = monomer('IAB',50_i16,78.5_dp,[0.1_dp,0.2_dp],[0.5_dp,0.5_dp])
  print*, 'K = ', dimer(1) % k, ' || ', ' P = ', dimer(1) % p
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing chain elongation and storage for a chain accessible through use association (chain declared in a module).'
  c_chain = dimer(1) % name(1:1)
  print*, 'Before elongation: chain length = ', len(c_chain),' and chain = ', c_chain
  do counter = 1, 5
    call chain_store(o_chain(1),c_chain)
    call chain_grow(c_chain,trim(dimer(1) % name(1:1)))
  end do
  call chain_store(o_chain(1),c_chain)
  print*, 'After ', counter-1_i2, ' additions of initiator: chain length = ', len(c_chain),' and chain = ', c_chain
  print*, 'Printing the chain at every step'
  do counter = 1, 6
    print*, o_chain(1) % store(counter)(1:o_chain(1) % length(counter))
  end do
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing chain elongation and storage through host association (chain declared in the main program).'
  call refresh_chain_storage(o_chain(1))
  allocate(character :: test_chain)
  test_chain = trim(dimer(1) % name(1:3))
  print*, 'Before elongation: chain length = ', len(test_chain),' and chain = ', test_chain
  do counter = 1, 5
    call chain_store(o_chain(1),test_chain)
    call chain_grow(test_chain,trim(dimer(1) % name(1:3)))
  end do
    call chain_store(o_chain(1),test_chain)
  print*, 'After ', counter-1_i2, ' additions of initiator: chain length = ', len(test_chain),' and chain = ', test_chain
  print*, '-----------------------------------------'
  print*, 'Printing the chain at every step'
  do counter = 1, 6
    print*, o_chain(1) % store(counter)(1:o_chain(1) % length(counter))
  end do
  print*, '-----------------------------------------'

  print*, ''
  c_chain = 'TransferTest'
  !write(*, '(A)', advance = "no"), ' Which entry do you want to remove (keep it in [1,6])? '
  !read(*,*), remove_entry
  remove_entry = 6
  print*, "Testing transfer termination and storage. We're using the", remove_entry,'th entry of o_chain(1)'
  call transfer(o_chain(1), o_chain(2), c_chain, remove_entry)
  do counter = 1, 6
    if (o_chain(1) % length(counter) == 0) then
      print*, 'Old chain #', counter, ' has become the new current chain.'
      ! We go to the very end of this do iteration (we skip everything between here and end do).
      cycle
    endif
    print*, 'Old chain #', counter, ' = ' , o_chain(1) % store(counter)(1:o_chain(1) % length(counter))
  end do
  print*, ''
  print*, 'Chain ended by transfer = ', o_chain(2) % store(1)(1:o_chain(2) % length(1))
  print*, ''
  print*, 'New current chain = ', c_chain
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing internal write to string and appending it to another one.'
  deallocate (o_chain(2) % store)
  allocate (character(5) :: o_chain(2) % store(5))
  do counter = 1, 4
    write (Unit=str, FMT="(I1)") counter
    test_chain = 'abcd' // str
    o_chain(2) % store(counter)(1:5) = test_chain
    o_chain(2) % length(counter) = len(o_chain(2) % store(counter)(1:5))
    print*, counter, o_chain(2) % store(counter)(1:5), o_chain(2) % length(counter)
  end do
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing string reversal.'
  forall (i=1:len(test_chain)) test_chain(i:i) = test_chain(len(test_chain)-i+1:len(test_chain)-i+1)
  print*, test_chain
  print*, '-----------------------------------------'
  do counter = 1, 5
    !print*, counter, o_chain(1) % store(counter)(1:o_chain(1)%length(counter))
    test_chain = o_chain(1) % store(counter)(1:o_chain(1)%length(counter))
    print*, test_chain
    forall (i=1:len(test_chain)) test_chain(i:i) = test_chain(len(test_chain)-i+1:len(test_chain)-i+1)
    print*, test_chain
  end do
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing string reversal subroutine.'
  print*, 'Before reversal: ', test_chain
  call chain_reverse(test_chain)
  print*, 'After reversal:  ', test_chain
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing recombination subroutine.'
  if (remove_entry <= 6 .and. remove_entry > 2) then
    remove_entry = remove_entry - 2
  else
    remove_entry = remove_entry + 1
  end if
  print*, "Testing recombination termination and storage. We're using the", remove_entry,'th entry of o_chain(1)'
    print*, 'Old chain to be recombined: ', o_chain(1) % store(remove_entry)(1:o_chain(1) % length(remove_entry))
    print*, 'Current chain to be inverted and recombined: ', test_chain
  call recombination(o_chain(1), o_chain(3), test_chain, remove_entry)
  do counter = 1, 6
    if (o_chain(1) % length(counter) == 0) then
      print*, 'Old chain #', counter, ' has been removed.'
      ! We go to the very end of this do iteration (we skip everything between here and end do).
      cycle
    endif
    print*, 'Old chain #', counter, ' = ' , o_chain(1) % store(counter)(1:o_chain(1) % length(counter))
  end do
  print*, ''
  print*, 'Chain ended by recombination = ', test_chain
  print*, '-----------------------------------------'

  print*, ''
  print*, 'Testing random number generation.'
  call ZBQLINI(0)
  do counter = 1, 10
    print*, ZBQLU01(dummy)
  end do

  call cpu_time(end)
  print*, ''
  print*, 'Total execution time = ', end-start, 'seconds.'
  print*, '-----------------------------------------'
end program test_subroutines
