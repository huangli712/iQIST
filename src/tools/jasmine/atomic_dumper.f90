!=========================================================================!
! project : rambutan
! program : atomic_dumper
! history : May 07, 2011
! authors : duliang (duleung@gmail.com)
! purpose : this file is just for test purpose
! comment : 
!=========================================================================!
!>>> dump eigen-values of atomic Hamiltonian to file
  subroutine atomic_dump_heigs(ncfgs, heigs)
     use constants

     implicit none

! external arguments
! number of configurations
     integer :: ncfgs

! eigenvalues of atomic Hamiltonian
     real(dp) :: heigs(ncfgs)

! local variables
! loop index over basis
     integer :: ibas

! open output file
     open(mytmp, file="atom.eval.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         write(mytmp, "(I8, F17.10)") ibas, heigs(ibas)
     enddo ! over ivec={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_heigs

!>>> dump eigen-states of atomic Hamiltonian to file
  subroutine atomic_dump_heigv(norbs, ncfgs, heigv, invcd)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary code 
     integer, intent(in) :: invcd(norbs, ncfgs)

! eigenvalues of "ntot" subspace
     complex(dp), intent(in) :: heigv(ncfgs, ncfgs)

! local variables
! loop index over basis
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.eigs.in", status="unknown")

! dumper data into output file
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(heigv(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F17.10, 4X, 14I1)") ibas, jbas, heigv(ibas, jbas), invcd(:, ibas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_heigv

  subroutine atomic_dump_hmtrx(ncfgs, hmtrx)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! atomic Hamiltonian matrix
     complex(dp), intent(in) :: hmtrx(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.hmat.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (abs(hmtrx(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F17.10)") ibas, jbas, hmtrx(ibas, jbas)
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_hmtrx

  subroutine atomic_dump_amtrx(ncfgs, amtrx)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! transformation matrix to make local Fock terms absent
     complex(dp), intent(in) :: amtrx(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.amat1.in", status="unknown")

! dumper data into output file
     do jbas=1,ncfgs
         do ibas=1,ncfgs
             if (abs(amtrx(ibas, jbas)) .lt. eps6) cycle
             write(mytmp, "(2I8, 2F17.10)") ibas, jbas, amtrx(ibas, jbas)
         enddo ! over ibas={1,ncfgs} loop
     enddo ! over jbas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_amtrx

  subroutine atomic_dump_eloc(norbs, eloc)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: norbs

! local onsite energy levels (should be diagonal)
     real(dp), intent(in) :: eloc(norbs)

! local variables
! loop index over configurations
     integer :: ibas

! open output file
     open(mytmp, file="atom.eloc.in", status="unknown")

! dumper data into output file
     do ibas=1,norbs
         write(mytmp, "(I8, F17.10)") ibas, eloc(ibas)
     enddo ! over ibas={1,ncfgs} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_eloc

  subroutine atomic_dump_basis(norbs, ncfgs, basis, invsn, invcd)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! decimal representation of Fock basis
     integer, intent(in) :: basis(ncfgs)

! serial number of decimal represented Fock basis
     integer, intent(in) :: invsn(0:2**norbs-1)

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas

! open output file
     open(mytmp, file="atom.basis.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         write(mytmp, "(3I8,3X,14I2)") ibas, basis(ibas), &
         invsn(basis(ibas)), invcd(1:norbs, ibas)
     enddo ! over ibas={1,nconf(ntot)} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_basis

  subroutine atomic_dump_vpmsy(ncfgs, vpmsy)
     use constants

     implicit none

! number of configurations
     integer, intent(in) :: ncfgs

! symmetry of variational parameters
     integer, intent(in) :: vpmsy(ncfgs, ncfgs)

! local variables
! loop index over configurations
     integer :: ibas
     integer :: jbas

! open output file
     open(mytmp, file="atom.vpmsy.in", status="unknown")

! dumper data into output file
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (vpmsy(ibas, jbas) .eq. 0) cycle
             write(mytmp, "(3I8)") ibas, jbas, vpmsy(ibas, jbas)
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,nconf(ntot)} loop

! close output file
     close(mytmp)

     return
  end subroutine atomic_dump_vpmsy
