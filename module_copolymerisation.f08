module polymer_constituents
  use numbers
  implicit none

  type monomer
    character(1) :: name
    integer(i16)  :: amount
    real(dp)  :: mass
    real(dp), dimension(:), allocatable :: k_ij
  end type monomer

  type dimer
    type(monomer) :: I, A, B
  end type dimer

  type termination
    ! Chain disproportionation
    real(dp) :: disp
    ! Chain transfer
    real(dp) :: tran
    ! Chain recombination
    real(dp) :: reco
  end type termination

  interface operator(/)
    module procedure norm_term, norm_dimer_kij
  end interface

contains
  function norm_term(a,b)
    type(termination) :: norm_term
    type(termination), intent(in) :: a
    type(termination), intent(in) :: b
    real(dp) :: norm

    norm = abs(b%disp) + abs(b%tran) + abs(b%reco)
    norm_term%disp = a%disp/norm
    norm_term%tran = a%tran/norm
    norm_term%reco = a%reco/norm
  end function norm_term

  function norm_dimer_kij(a,b)
    type(dimer) :: norm_dimer_kij
    type(dimer), intent(in) :: a
    type(dimer), intent(in) :: b
    real(dp) :: norm_i, norm_ab

    norm_i  = sum(abs(b%I%k_ij))
    norm_ab = sum(abs(b%A%k_ij)) + sum(abs(b%B%k_ij))

    norm_dimer_kij%I%k_ij = a%I%k_ij/norm_i
    norm_dimer_kij%A%k_ij = a%A%k_ij/norm_ab
    norm_dimer_kij%B%k_ij = a%B%k_ij/norm_ab
  end function norm_dimer_kij

end module polymer_constituents
