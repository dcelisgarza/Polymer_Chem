program main
  use two_monomer_data_declaration
  implicit none

  ! Initial memory allocation.
  call allocation

  ! Setting initial monomer parameters.
  DI: do concurrent (i = 1: 3)
    dimer(i) % k(1) = 0.
    dimer(i) % reacted = 0
  end do DI
  dimer(1) % name = 'I'
  dimer(2) % name = 'A'
  dimer(3) % name = 'B'
  ! Read parameters from user input.
  write(*, '(A)', advance = "no") ' [I] = '; read(*,*) dimer(1) % amount
  write(*, '(A)', advance = "no") ' [A] = '; read(*,*) dimer(2) % amount
  write(*, '(A)', advance = "no") ' [B] = '; read(*,*) dimer(3) % amount
  write(*, '(A)', advance = "no") ' k_ia = '; read(*,*) dimer(1) % k(2)
  write(*, '(A)', advance = "no") ' k_ib = '; read(*,*) dimer(1) % k(3)
  write(*, '(A)', advance = "no") ' k_aa = '; read(*,*) dimer(2) % k(2)
  write(*, '(A)', advance = "no") ' k_ab = '; read(*,*) dimer(2) % k(3)
  write(*, '(A)', advance = "no") ' k_ba = '; read(*,*) dimer(3) % k(2)
  write(*, '(A)', advance = "no") ' k_bb = '; read(*,*) dimer(3) % k(3)

  ! Setting initial termination parameters.
  term(1) % name = 'Disproportiation'
  term(2) % name = 'Recombination'
  term(3) % name = 'Transfer'
  ! Reading parameters from user input.
  write(*, '(A)', advance = "no") ' P('//term(1) % name(1:4)//') = '; read(*,*) term(1) % p(1)
  write(*, '(A)', advance = "no") ' P('//term(2) % name(1:4)//') = '; read(*,*) term(2) % p(1)
  write(*, '(A)', advance = "no") ' P('//term(3) % name(1:4)//') = '; read(*,*) term(3) % p(1)
  write(*, '(A)', advance = "no") ' Lambda = '; read(*,*) term(1) % kl(1)
  KL: do concurrent (i = 2:3)
    term(i) % kl(1) = term(1) % kl(1)
  end do KL

  ! Initialising monomers.
  ! dimer(1) = monomer('I', amount, mass, [kii, kia ,kib], [pii, pia, pib], reacted)
  !dimer(1) = monomer('I', 1000, 1, [0., 1. ,1.], [0., 0., 0.], 0)
  ! dimer(2) = monomer('A', amount, mass, [kai, kaa ,kab], [pai, paa, pab], reacted)
  !dimer(2) = monomer('A', 40000, 1, [0., 2. ,1.], [0., 0., 0.], 0)
  ! dimer(3) = monomer('B', amount, mass, [kbi, kba ,kbb], [pbi, pba, pbb], reacted)
  !dimer(3) = monomer('B', 50000, 1, [0., 4. ,1.], [0., 0., 0.], 0)

  ! Initialising terminations.
  ! term(x) = termination(name, kinetic chain length, termination probability)
  !term(1) = termination('Disproportiation',[2000],[0.])
  !term(2) = termination('Recombination',[3000],[0.])
  !term(3) = termination('Transfer',[5000],[0.])

  !write(*,*) ''
  !call term_prob_from_kinetic_chain_len(term)
  !write(*,*) 'Termination probabilities from kinetic chain lengths:'
  !do i = 1, n_tot
  !  write(*,'(A, f0.10)') ' '//term(i) % name//' = ', term(i) % p
  !end do
  !write(*,*) ''

  call norm_rtn_coeff(dimer)

  IRP1: do concurrent (i = 1: n_tot)
    IRP2: do concurrent (j = 1: n_tot)
      call rtn_prob(dimer, i, j)
    end do IRP2
  end do IRP1

  call norm_term_prob(term)

  write(*,*) ''
  write(*,*) 'Normalised reaction constants:'
  write(*,*) '    ',dimer(1) % name, '      ',dimer(2) % name, '      ', dimer(3) % name
  NK: do i = 1, n_tot
    write(*,'(A, f0.5, A, f0.5, A, f0.5)') dimer(i) % name//' ', dimer(i) % k(1), ' ', dimer(i) % k(2), ' ', dimer(i) % k(3)
  end do NK

  write(*,*) ''
  write(*,*) 'Normalised reaction probabilites:'
  write(*,*) '    ',dimer(1) % name, '      ',dimer(2) % name, '      ', dimer(3) % name
  NP: do i = 1, n_tot
    write(*,'(A, f0.5, A, f0.5, A, f0.5)') dimer(i) % name//' ', dimer(i) % p(1), ' ', dimer(i) % p(2), ' ', dimer(i) % p(3)
  end do NP

  write(*,*) ''
  write(*,*) 'Normalised termination probabilities'
  NTP: do i = 1, n_tot
    write(*,'(A, f0.10)') ' P('//term(i) % name(1:4)//') = ', term(i) % p(1)
  end do NTP

  ! Initialise RNG seed.
  call ZBQLINI(0)

  ! Polymerisation Loop.
  ! Enters keeps going while there's material to polymerise.
  PL: do while (dimer(1) % amount > 0 .or. dimer(2) % amount > 0 .or. dimer(3) % amount > 0)

    ! Checking if the current Chain is Initial and we have Monomers. If we don't we exit.
    CCM: if ( c_chain == '' .and. dimer(1) % amount == 0 ) then
      ! If we need to start a new chain but can't for lack of initiators then we exit.
      exit PL
    ! If we need to start a new chain but we have initiators then we keep going.
    else if ( c_chain == '' .and. dimer(1) % amount > 0 ) then CCM
      ! We remove an initator because it's now part of the chain.
      c_chain = dimer(1) % name
      dimer(1) % amount = dimer(1) % amount - 1
      ! Calculate reaction Probabilities for IA and IB.
      RIAIB: do concurrent (i = 1: n_tot)
        call rtn_prob(dimer, 1, i)
        work(i) = dimer(1) % p(i)
      end do RIAIB
      choice = ZBQLU01(0)
      call prob_limits(limit, work, n_tot)
      ! Check if choice is in (0, PXI].
      M: if (choice <= limit(1)) then
        ! React with I.
        i_choice = 1
        call grow_chain(c_chain, dimer(i_choice) % name)
        dimer(i_choice) % amount = dimer(i_choice) % amount - 1
      ! Check if choice is in (PXI, PXI + PXA]
      else if (limit(1) < choice .and. choice <= limit(2)) then M
        ! React with A.
        i_choice = 2
        call grow_chain(c_chain, dimer(i_choice) % name)
        dimer(i_choice) % amount = dimer(i_choice) % amount - 1
      ! Check if choice is in (PXI + PXA, PXI + PXA + PXB = 1 ]
      else if (limit(2) < choice) then M
        ! React with B.
        i_choice = 3
        call grow_chain(c_chain, dimer(i_choice) % name)
        dimer(i_choice) % amount = dimer(i_choice) % amount - 1
      end if M
      ! Test polymer initiation.
      !write(*,*) choice, i_choice, c_chain, dimer(1) % amount, dimer(2) % amount, dimer(3) % amount
    end if CCM

    ! Exiting if there is no more initiator and chains can't continue growing via transfer or recombination.
    !ITR: if ( dimer(1) % amount == 0 .and. term(2) % p(1) == 0. .and. term(3) % p(1) == 0. ) then
    !  exit PL
    !end if ITR

  end do PL

end program main
