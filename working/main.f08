program main
  use two_monomer_data_declaration
  implicit none
  integer :: i, j
  ! Initial memory allocation.
  call allocation

  ! Initialising monomers.
  ! dimer(1) = monomer('I', amount, mass, [kii, kia ,kib], [pii, pia, pib], reacted)
  dimer(1) = monomer('I', 1, 1, [0., 1. ,1.], [0., 0., 0.], 0)
  ! dimer(2) = monomer('A', amount, mass, [kai, kaa ,kab], [pai, paa, pab], reacted)
  dimer(2) = monomer('A', 20, 1, [0., 2. ,1.], [0., 0., 0.], 0)
  ! dimer(3) = monomer('B', amount, mass, [kbi, kba ,kbb], [pbi, pba, pbb], reacted)
  dimer(3) = monomer('B', 80, 1, [0., 4. ,1.], [0., 0., 0.], 0)

  ! Initialising terminations.
  term(1) = termination('Disproportiation',[2.])
  ! Initialising terminations.
  term(2) = termination('Transfer',[1.])
  ! Initialising terminations.
  term(3) = termination('Recombination',[1.])


  call norm_termination(term)
  call norm_monomer(dimer)
  do i = 1, n_tot
    do j = 1, n_tot
      call reaction_probability(i, j, dimer)
    end do
  end do

  print*, ''
  print*, 'Normalised termination probabilities:'
  do i = 1, n_tot
    write(*,'(A, f0.5)') ' '//term(i) % name//' = ', term(i) % p(1)
  end do
  print*, ''

  print*, 'Normalised reaction constants:'
  print*, '    ',dimer(1) % name, '      ',dimer(2) % name, '      ', dimer(3) % name
  do i = 1, n_tot
    write(*,'(A, f0.5, A, f0.5, A, f0.5)') dimer(i) % name//' ', dimer(i) % k(1), ' ', dimer(i) % k(2), ' ', dimer(i) % k(3)
  end do
  print*, ''

  print*, 'Normalised reaction probabilites:'
  print*, '    ',dimer(1) % name, '      ',dimer(2) % name, '      ', dimer(3) % name
  do i = 1, n_tot
    write(*,'(A, f0.5, A, f0.5, A, f0.5)') dimer(i) % name//' ', dimer(i) % p(1), ' ', dimer(i) % p(2), ' ', dimer(i) % p(3)
  end do
end program main
