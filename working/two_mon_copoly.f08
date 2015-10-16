program two_mon_copoly
  !use two_monomer_data_declaration
  use polymerisation
  character(:), allocatable :: chain
  integer(i2) :: counter
  real(dp) :: start, end

  call cpu_time(start)
  call allocation
  dimer(1) % name = 'I'
  dimer(1) % amount = 50_i16
  dimer(1) % mass = 78.5_dp
  dimer(1) % k = [0.1_dp,0.2_dp]
  dimer(1) % p = [0.5_dp,0.5_dp]
  ! Correct assignation causes core dump :(
  ! dimer(1) = monomer('I',50_i16,78.5_dp,[0.1_dp,0.2_dp],[0.5_dp,0.5_dp])
  print*, 'K = ', dimer(1) % k, ' || ', ' P = ', dimer(1) % p
  print*, '-----------------------------------------'

  allocate(character(len=1) :: chain)
  chain = dimer(1) % name
  print*, 'Before elongation: chain length = ', len(chain),' and chain = ', chain
  do counter = 1, 5
    call chain_grow(chain,dimer(1))
  end do
  call cpu_time(end)
  print*, 'After ', counter-1_i2, ' additions of initiator: chain length = ', len(chain),' and chain = ', chain
  print*, 'Execution time = ', end-start, 'seconds.'

  !subroutine chain_grow(chain,mon)
  !  implicit none
  !  type(monomer), intent(in)        :: mon
  !  character(:), allocatable, intent(inout) :: chain
  !  !character(len = len(old) + 1_i2) :: new
  !  character(:), allocatable        :: work
  !  !new = old // mon%name
  !  allocate(character(len=len(chain) + 1_i2) :: work)
  !  work = chain // mon % name
  !  deallocate(chain)
  !  allocate(character(len=len(work)) :: chain)
  !  chain = work
  !end subroutine chain_grow
end program two_mon_copoly
