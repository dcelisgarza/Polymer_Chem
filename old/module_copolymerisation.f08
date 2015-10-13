module polymer_constituents
  ! Use the numbers module which defines the precision and length of various numbers and constants.
  use numbers
  implicit none

  type name
    character(:), allocatable :: name(:)!, dimension(:), allocatable :: name(:)
  end type name

  type, extends(name) :: monomer
    ! Monomer type.
    !type(name)    :: name
    integer(i16)  :: amount ! Monomer amount.
    real(dp)      :: mass   ! Monomer mass.
    real(dp), dimension(:), allocatable :: k_ij ! Monomer reaction coefficients. The 'allocatable' attribute makes it possible to declare the extents of the defined ranks at compilation time.
  end type monomer

  type dimer
    ! Dimer type has three constituents: Initiator, Monomer A and Monomer B.
    !type(name_term) :: name
    type(monomer) :: I, A, B
  end type dimer

  type, extends(name) :: termination
    ! Termination type. Allows for different terminations and termination lengths.
    !type(name) :: name ! Termination name.
    real(dp)   :: prob ! Termination probability.
    real(dp)   :: klen ! Kinetic chain length.
  end type termination

  type terminations
    ! 1) Disproportiation, 2) Transfer, 3) Recombination
    type(termination) :: disp, tran, reco
  end type terminations

  interface operator(/)
    ! Defines an operator which calls subroutines where operations are defined for derived type structures. In this case, for normalisation operations on parts of data structures with type 'dimer' or 'termination'.
    module procedure norm_term, norm_dimer_kij
  end interface

  interface operator(==)
    ! Defines an operator which calls subroutines where operations are defined for derived type structures. In this case, for normalisation operations on parts of data structures with type 'dimer' or 'termination'.
    module procedure asgn_dimer
  end interface

  contains
    function norm_term(a,b)
      ! Normalises the termination probabilities.
      type(terminations) :: norm_term      ! Result.
      real(dp)           :: norm           ! Sum of all probabilities.
      type(terminations), intent(in) :: a,b

      norm = a % disp % prob + a % tran % prob + a % reco % prob ! Sum all probabilities.
      norm_term % disp % prob = a % disp % prob/norm  ! Normalise disp.
      norm_term % tran % prob = a % tran % prob/norm  ! Normalise tran.
      norm_term % reco % prob = a % reco % prob/norm  ! Normalise reco.
    end function norm_term

    function norm_dimer_kij(a,b)
      ! Normalise reaction coefficients.
      type(dimer) :: norm_dimer_kij           ! Result.
      real(dp)    :: norm_i, norm_a, norm_b   ! Sum of coefficients for possible  reactions.
      type(dimer), intent(in) :: a, b

      norm_i = sum(abs(a % I % k_ij)) ! Sum coefficients for I_ reactions.
      norm_a = sum(abs(a % A % k_ij)) ! Sum coefficients for A_ reactions.
      norm_b = sum(abs(a % B % k_ij)) ! Sum coefficients for B_ reactions.

      norm_dimer_kij % I % k_ij = a % I % k_ij/norm_i ! Normalise I_ reactions.
      norm_dimer_kij % A % k_ij = a % A % k_ij/norm_a ! Normalise A_ reactions.
      norm_dimer_kij % B % k_ij = a % B % k_ij/norm_b ! Normalise B_ reactions.
    end function norm_dimer_kij

    function asgn_dimer(a,b)
      type(dimer) :: asgn_dimer
      type(dimer), intent(in) :: a, b

      asgn_dimer % I % name   = a % I % name
      asgn_dimer % I % amount = a % I % amount
      asgn_dimer % I % mass   = a % I % mass
      asgn_dimer % I % k_ij   = a % I % k_ij

      asgn_dimer % A % name   = a % A % name
      asgn_dimer % A % amount = a % A % amount
      asgn_dimer % A % mass   = a % A % mass
      asgn_dimer % A % k_ij   = a % A % k_ij

      asgn_dimer % B % name   = a % B % name
      asgn_dimer % B % amount = a % B % amount
      asgn_dimer % B % mass   = a % B % mass
      asgn_dimer % B % k_ij   = a % B % k_ij
    end function asgn_dimer

end module polymer_constituents
