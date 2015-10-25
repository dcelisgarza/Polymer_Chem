program two_mon_copoly
  !use two_monomer_data_declaration
  use polymerisation
  implicit none
  character(:), allocatable :: test_chain
  integer(i1) :: counter
  real(dp) :: start, end

  call cpu_time(start)
  call allocation
  dimer(1) % name(1)(1:3) = 'IAB'
  dimer(1) % amount = 50_i16
  dimer(1) % mass = 78.5_dp
  dimer(1) % k = [0.1_dp,0.2_dp]
  dimer(1) % p = [0.5_dp,0.5_dp]
  o_chain % index = 1
  ! Correct assignation causes core dump :(
  ! dimer(1) = monomer('I',50_i16,78.5_dp,[0.1_dp,0.2_dp],[0.5_dp,0.5_dp])
  print*, 'K = ', dimer(1) % k, ' || ', ' P = ', dimer(1) % p
  print*, '-----------------------------------------'

  allocate(character :: test_chain)
  test_chain = trim(dimer(1) % name(1)(1:3))
  print*, 'Before elongation: chain length = ', len(test_chain),' and chain = ', test_chain
  do counter = 1, 5
    call chain_store(o_chain(1),test_chain)
    call chain_grow(test_chain,trim(dimer(1) % name(1)(1:3)))
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
  print*, 'Testing module allocation function for current chain.'
  c_chain = dimer(1) % name(1)(1:1)
  print*, 'Before elongation: chain length = ', len(c_chain),' and chain = ', c_chain
  do counter = 1, 5
    call chain_store(o_chain(1),c_chain)
    call chain_grow(c_chain,trim(dimer(1) % name(1)(1:1)))
  end do
  call chain_store(o_chain(1),c_chain)
  print*, 'After ', counter-1_i2, ' additions of initiator: chain length = ', len(c_chain),' and chain = ', c_chain
  call cpu_time(end)
  print*, '-----------------------------------------'
  print*, 'Total execution time = ', end-start, 'seconds.'
  print*, '-----------------------------------------'
end program two_mon_copoly
