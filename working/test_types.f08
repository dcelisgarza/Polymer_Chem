program test
  use numbers
  use mod_config_polymer_test

  type(monomer)     :: I, A, B
  real(dp)          :: dum_sum
  type(termination) :: term_norm, disp, tran, reco
  type(monomer)     :: dimer(3)
  type(termination) :: term(4)

  ! Allocate monomer names and reaction coefficients.
  allocate ( character(len=1) :: I % name, A % name, B % name )
  allocate ( I % k_ij(2), A % k_ij(2), B % k_ij(2) )

  ! Allocate termination names.
  allocate ( character(len=18) :: disp % name )
  allocate ( character(len=8)  :: tran % name )
  allocate ( character(len=13) :: reco % name )

  ! Testing the termination normalisation function using termination.
  disp % p = 1._dp
  tran % p = 0.6_dp
  reco % p = 0.4_dp
  term_norm = disp + tran + reco
  disp = disp/term_norm
  tran = tran/term_norm
  reco = reco/term_norm
  print*, term_norm % p, disp % p, tran % p, reco % p
  ! It works.

  print*, '-------------------'

  ! Testing the reaction coefficient normalisation function of monomers.
  I % k_ij(1) = 1._dp
  I % k_ij(2) = 0.5_dp
  I = I/I

  A % k_ij(1) = 1._dp
  A % k_ij(2) = 1._dp
  A = A/A

  B % k_ij(1) = 3._dp
  B % k_ij(2) = 1._dp
  B = B/B
  print*, I % k_ij, A % k_ij, B % k_ij
  ! It works.

  print*, '-------------------'
  print*, '-------------------'

  ! Allocate polymer names and reaction coefficients.
  allocate ( character(len=1) :: dimer(1) % name, dimer(2) % name, dimer(3) % name )
  allocate ( dimer(1) % k_ij(2), dimer(2) % k_ij(2), dimer(3) % k_ij(2) )

  ! Allocate termination array names.
  allocate ( character(len=18) :: term(1) % name )
  allocate ( character(len=8)  :: term(2) % name )
  allocate ( character(len=13) :: term(3) % name )

  ! Testing the termination normalisation function using termination array.
  term(1) % p = 1._dp
  term(2) % p = 0.6_dp
  term(3) % p = 0.4_dp
  term(4) = term(1) + term(2) + term(3)
  term(1) = term(1)/term(4)
  term(2) = term(2)/term(4)
  term(3) = term(3)/term(4)
  print*, term(4) % p, term(1) % p, term(2) % p, term(3) % p
  ! It works.

  print*, '-------------------'

  ! Testing the reaction coefficient normalisation function of polymer.
  dimer(1) % k_ij(1) = 1._dp
  dimer(1) % k_ij(2) = 0.5_dp
  dimer(1) = dimer(1)/dimer(1)

  dimer(2) % k_ij(1) = 1._dp
  dimer(2) % k_ij(2) = 1._dp
  dimer(2) = dimer(2)/dimer(2)

  dimer(3) % k_ij(1) = 3._dp
  dimer(3) % k_ij(2) = 1._dp
  dimer(3) = dimer(3)/dimer(3)
  print*, dimer(1) % k_ij, dimer(2) % k_ij, dimer(3) % k_ij
  ! It works.

end program test
