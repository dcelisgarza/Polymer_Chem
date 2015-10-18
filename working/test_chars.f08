program test_char_ar
  use polymerisation
  implicit none
  character(30) :: page(4)
  character(len=30) :: line(4)
  integer(i1) :: i
  character(len=30), allocatable :: dimd(:)
  character(:), allocatable :: dimc
  character(:), allocatable :: dimi(:)
  character(:), allocatable :: dimn(:)

  line(1) = '1'
  line(2) = '12345'
  line(3) = '1234567890'
  line(4) = '1234567890ABCDEFGHIJabcdefghij'

  print*, 'Testing with fixed dimension and length.'
  do i = 1, 4
    page(i)(1:len(trim(line(i)))) = line(i)
    print*, page(i)(1:len(trim(line(i))))
  end do
  print*, '---------------------------------'

  print*, 'Testing with allocatable dimension and fixed length.'
  allocate( dimd(4) )
  do i = 1, 4
    dimd(i)(1:len(trim(line(i)))) = line(i)
    print*, dimd(i)(1:len(trim(line(i))))
  end do
  print*, '---------------------------------'

  print*, 'Testing with fixed dimension and allocatable length (storing individual chains).'
  allocate( character :: dimc )
  do i = 1, 4
    dimc = line(i)
    print*, trim(dimc), len(trim(dimc))
  end do
  print*, '---------------------------------'

  print*, 'Testing with allocatable dimension and length (storing individual chains in an array).'
  allocate( character :: dimi(4))
  do i = 1, 4
    dimi(i)(1:len(trim(line(i)))) = line(i)
    print*, dimi(i)(1:len(trim(line(i)))), len(dimi(i)(1:len(trim(line(i)))))
  end do
  print*, '---------------------------------'

  print*, 'Testing information transfer from one character array to another with greater extent. &
         & This will help in dynamically growing the storage arrays.'
  allocate( character :: dimn(5))
  do i = 1, 4
    dimn(i)(1:len(trim(line(i)))) = dimi(i)(1:len(trim(line(i))))
    print*, dimn(i)(1:len(trim(line(i))))
  end do
  print*, '---------------------------------'

  print*, 'Testing going back to the initial array storage and adding a new entry.'
  deallocate(dimi)
  allocate( character :: dimi(5))
  dimi(5)(1:4) = 'WAKA'
  do i = 1, 4
    dimi(i)(1:len(trim(line(i)))) = dimn(i)(1:len(trim(line(i))))
    print*, dimi(i)(1:len(trim(line(i))))
  end do
  print*, dimi(5)(1:len('waka'))

end program test_char_ar
