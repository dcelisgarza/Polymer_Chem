program copolymerisation
  call polymerise
end program copolymerisation

subroutine polymerise
  use polymer_constituents
  implicit none
  type(dimer)       :: polymer
  type(termination) :: term
  ! Allocate the reaction coefficients.
  ! The first element is reaction with A.
  ! The second is reaction with B.
  allocate (polymer%I%k_ij(2),&
            polymer%A%k_ij(2),&
            polymer%B%k_ij(2))

  polymer%I%name = 'I'
  polymer%A%name = 'A'
  polymer%B%name = 'B'

! Concentration.
  write(*,'(A)',advance='no') '[I] = '
  read(*,*) polymer%I%amount
  write(*,'(A)',advance='no') '[A] = '
  read(*,*) polymer%A%amount
  write(*,'(A)',advance='no') '[B] = '
  read(*,*) polymer%B%amount

! Mass.
  write(*,'(A)',advance='no') 'MM_I = '
  read(*,*) polymer%I%mass
  write(*,'(A)',advance='no') 'MM_A = '
  read(*,*) polymer%A%mass
  write(*,'(A)',advance='no') 'MM_B = '
  read(*,*) polymer%B%mass

! Reaction probabilities.
  write(*,'(A)',advance='no') 'K_IA = '
  read(*,*) polymer%I%k_ij(1)
  write(*,'(A)',advance='no') 'K_IB = '
  read(*,*) polymer%I%k_ij(2)
  write(*,'(A)',advance='no') 'K_AA = '
  read(*,*) polymer%A%k_ij(1)
  write(*,'(A)',advance='no') 'K_AB = '
  read(*,*) polymer%A%k_ij(2)
  write(*,'(A)',advance='no') 'K_BA = '
  read(*,*) polymer%B%k_ij(1)
  write(*,'(A)',advance='no') 'K_BB = '
  read(*,*) polymer%B%k_ij(2)
! Normalising reaction coefficients so they add up to 1.
  polymer = polymer/polymer

! Termination probabilities
  write(*,'(A)',advance='no') 'P(Disproportionation) = '
  read(*,*) term%disp
  write(*,'(A)',advance='no') 'P(Transfer) = '
  read(*,*) term%tran
  write(*,'(A)',advance='no') 'P(recombination) = '
  read(*,*) term%reco
! Normalising termination probabilities so they add up to 1.
  term = term/term
end subroutine polymerise
