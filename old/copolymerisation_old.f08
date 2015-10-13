program copolymerisation
  call polymerise
end program copolymerisation

subroutine polymerise
  ! Use the polymer constituents module. It defines the monomer, dimer and termination data types. It also defines the normalisation operations on the dimer and termination types.
  use polymer_constituents
  implicit none
  ! Defining the polymer as a two-monomer copolymer.
  type(dimer)       :: polymer
  ! Defining our termination variable.
  type(termination) :: term
  ! Dummy name variables.
  type(name) :: dname1, dname2, dname3, dname4
  allocate (character(1)  :: dname1%name(3))
  allocate (character(18) :: dname2%name(1))
  allocate (character(8)  :: dname3%name(1))
  allocate (character(13) :: dname4%name(1))
  ! Single character names.
  allocate (character(1) :: &
                            polymer % I % name % name(1), &
                            polymer % A % name % name(1), &
                            polymer % B % name % name(1))
  ! Allocate the reaction coefficients for each constituent of our two-monomer copolymer.
  ! The first element is of the k_ij array is the reaction with A,
  ! the second is reaction with B.
  allocate (polymer % I % k_ij(2), & ! K_IA and K_IB.
            polymer % A % k_ij(2), & ! K_AA and K_AB.
            polymer % B % k_ij(2))   ! K_BA and K_BB.

  ! Allocate different terminations.
  ! t: 1) Disproportionation, 2) Transfer, 3) Recombination.
  ! l: kinetic chain length, only 1 for this example.
  allocate (term % t(3), term % l(1))
  ! Termination name.
  allocate(character(18) :: term % name % name(3))

! What we're calling our monomers and initiator.
  !polymer % I % name % n(1) = 'I'
  !print*, polymer % I % name % n(1)
  !polymer % A % name % n(1) = 'A'
  !print*, polymer % A % name % n(1)
  !polymer % B % name (1) = 'B'
  !print*, polymer % I % name (1)

! Asking for the amount of each constituent (given as an integer quantity).
  write(*,'(A)',advance='no') '[I] = ' ! How much initiator?
  read(*,*) polymer % I % amount
  write(*,'(A)',advance='no') '[A] = ' ! How much A?
  read(*,*) polymer % A % amount
  write(*,'(A)',advance='no') '[B] = ' ! How much B?
  read(*,*) polymer % B % amount

! Asking for the mass of each constituent (given as a real max precision number).
  write(*,'(A)',advance='no') 'MM_I = ' ! Mass of the initiator?
  read(*,*) polymer % I % mass
  write(*,'(A)',advance='no') 'MM_A = ' ! Mass of A?
  read(*,*) polymer % A % mass
  write(*,'(A)',advance='no') 'MM_B = ' ! Mass of B?
  read(*,*) polymer % B % mass

! Asking for the reaction coefficients.
  write(*,'(A)',advance='no') 'K_IA = ' ! K_IA?
  read(*,*) polymer % I % k_ij(1)
  write(*,'(A)',advance='no') 'K_IB = ' ! K_IB?
  read(*,*) polymer % I % k_ij(2)
  write(*,'(A)',advance='no') 'K_AA = ' ! K_AA?
  read(*,*) polymer % A % k_ij(1)
  write(*,'(A)',advance='no') 'K_AB = ' ! K_AB?
  read(*,*) polymer % A % k_ij(2)
  write(*,'(A)',advance='no') 'K_BA = ' ! K_BA?
  read(*,*) polymer % B % k_ij(1)
  write(*,'(A)',advance='no') 'K_BB = ' ! K_BB?
  read(*,*) polymer % B % k_ij(2)
! Normalising reaction coefficients so they add up to 1.
  polymer = polymer/polymer

! Termination probabilities
  write(*,'(A)',advance='no') 'P(Disproportionation) = ' ! P(disproportionation) ?
  read(*,*) term % t(1)
  write(*,'(A)',advance='no') 'P(Transfer) = ' ! P(transfer) ?
  read(*,*) term % t(2)
  write(*,'(A)',advance='no') 'P(Recombination) = ' ! P(recombination) ?
  read(*,*) term % t(3)
! Normalising termination probabilities so they add up to 1.
  term = term/term
end subroutine polymerise
