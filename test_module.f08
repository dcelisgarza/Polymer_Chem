module polymer
  use numbers

  type monomer
    character(:), allocatable :: name(:) ! Monomer name.
    real(i16) :: amount ! Amount of monomer.
    real(dp)  :: mass   ! Monomer mass.
    real(dp), dimension(:), allocatable :: k_ij ! Kinetic reaction coefficients.
  end type monomer

  type termination
    character(:), allocatable :: name(:) ! Termination name.
    real(dp)      :: p ! Termination probability.
    integer(i16)  :: l ! Kinetic chain length.
  end type termination

  type terminations
    type(termination) :: disp, tran, reco
  end type terminations

  interface operator(+)
    module procedure sum_term, sum_mon
  end interface

  interface operator(/)
    module procedure norm_term, norm_mon
  end interface

contains

  function sum_term(a,b)
    type(termination) :: sum_term
    type(termination), intent(in) :: a, b

    sum_term%p = a%p + b%p
  end function sum_term

  function norm_term(a,b)
    type(termination) :: norm_term
    type(termination), intent(in) :: a, b

    norm_term%p = a%p / b%p
  end function norm_term

  function sum_mon(a,b)
    type(monomer) :: sum_mon
    type(monomer), intent(in) :: a, b

    sum_mon%k_ij = a%k_ij + b%k_ij
    sum_mon%amount = a%amount + b%amount
  end function sum_mon

  function norm_mon(a,b)
    type(monomer) :: norm_mon
    type(monomer), intent(in) :: a, b
    real(dp) :: norm

    norm = sum(abs(b%k_ij))
    norm_mon%k_ij = a%k_ij/norm
  end function norm_mon


end module polymer
