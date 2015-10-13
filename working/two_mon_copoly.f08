program two_mon_copoly
  use two_monomer_data_declaration
  call allocation


  dimer(1) % name = 'I'
  dimer(1) % amount = 50_i16
  dimer(1) % mass = 78.5_dp
  dimer(1) % k = [0.1_dp,0.2_dp]
  dimer(1) % p = [0.5_dp,0.5_dp]
  print*, dimer(1) % k, dimer(1) % p

  ! Correct assignation causes core dump :(
  ! dimer(1) = monomer('I',50_i16,78.5_dp,[0.1_dp,0.2_dp],[0.5_dp,0.5_dp])
end program two_mon_copoly
