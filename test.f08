program test
  use numbers
  use polymer

  type(monomer)     :: I, A, B
  real(dp) :: dum_sum
  type(termination) :: term_norm, disp, tran, reco



  allocate ( character(len=1) :: I % name(1), A % name(1), B % name(1) )
  allocate ( I % k_ij(2), A % k_ij(2), B % k_ij(2) )

  allocate ( character(len=18) :: disp % name(1) )
  allocate ( character(len=8)  :: tran % name(1) )
  allocate ( character(len=13) :: reco % name(1) )

  ! Testing the termination normalisation function using termination.
  disp%p = 1._dp
  tran%p = 0.6_dp
  reco%p = 0.4_dp
  term_norm = disp + tran + reco
  disp = disp/term_norm
  tran = tran/term_norm
  reco = reco/term_norm
  print*, term_norm%p, disp%p, tran%p, reco%p
  ! It works.

  print*, '_______________'

  ! Testing the reaction coefficient normalisation function.
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

end program test
