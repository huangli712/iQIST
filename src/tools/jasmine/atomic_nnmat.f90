
  subroutine atomic_make_nmat(norbs, ncfgs, invcd, nmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! nn operator matrix in many body basis
     complex(dp), intent(out) :: nmat(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     nmat = czero

     do ibas=1,ncfgs
         do iorb=1,norbs
             if (invcd(iorb, ibas) .eq. 1) then
                 nmat(ibas, ibas) = nmat(ibas, ibas) + done
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_make_nmat

  subroutine atomic_dump_upara(Uc, Uv, Jz, Js, Jp, lamb)
     use constants
     implicit none

     real(dp), intent(in) :: Uc
     real(dp), intent(in) :: Uv
     real(dp), intent(in) :: Jz
     real(dp), intent(in) :: Js
     real(dp), intent(in) :: Jp
     real(dp), intent(in) :: lamb

     open(mytmp, file="atom.Upara.in")
     write(mytmp, "(F17.10)") Uc
     write(mytmp, "(F17.10)") Uv
     write(mytmp, "(F17.10)") Jz
     write(mytmp, "(F17.10)") Js
     write(mytmp, "(F17.10)") Jp
     write(mytmp, "(F17.10)") lamb
     close(mytmp)

     return
  end subroutine atomic_dump_upara

 subroutine atomic_unit_matrix(ncfgs, zmtrx, heigv, deigs)
     use constants

     implicit none

! external arguments
! number of configurations
     integer, intent(in) :: ncfgs

! the input operator matrix in many body basis
     complex(dp), intent(in) :: zmtrx(ncfgs, ncfgs)

! eigenstate of atomic Hamiltonian
     complex(dp), intent(in) :: heigv(ncfgs, ncfgs)

! eigenvalue of the input operator matrix
     real(dp), intent(out) :: deigs(ncfgs)

! local variables
! loop index over orbits
     integer :: ibas
     integer :: jbas

! auxiliary complex(dp) matrix
     complex(dp) :: zmata(ncfgs, ncfgs)
     complex(dp) :: zmatb(ncfgs, ncfgs)

! unitary transformation of operator
     call zmat_zgemm0(ncfgs, zmtrx, heigv, zmata)
     call zmat_zgemm1(ncfgs, heigv, zmata, zmatb)

! make sure zmatb is a diagonal matrix
     do ibas=1,ncfgs
         do jbas=1,ncfgs
             if (ibas .eq. jbas) cycle
             if (abs(zmatb(ibas, jbas)) .gt. eps4) then
                 stop "severe error happened in atomic_unit_matrix"
             endif
         enddo ! over jbas={1,ncfgs} loop
     enddo ! over ibas={1,ncfgs} loop

! diagonal element is the eigenvalues of zmtrx
     do ibas=1,ncfgs
         deigs(ibas) = dble(zmatb(ibas, ibas))
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_unit_matrix

  subroutine atomic_build_nnmat(norbs, ncfgs, invcd, nnmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! nn operator matrix in many body basis
     complex(dp), intent(out) :: nnmat(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     nnmat = czero

     do ibas=1,ncfgs
         do iorb=1,norbs
             if (invcd(iorb, ibas) .eq. 1) then
                 nnmat(ibas, ibas) = nnmat(ibas, ibas) + done
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_build_nnmat

