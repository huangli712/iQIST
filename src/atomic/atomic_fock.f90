!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_ffmat
!!!           atomic_make_foccu
!!!           atomic_make_fspin
!!!           atomic_make_fhmat
!!! source  : atomic_full.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : computational subroutines for the calculations of occupation
!!!           number, spin moment, F-matrix, and atomic Hamiltonian in the
!!!           full Fock space.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_make_ffmat: make F-matrix for annihilation operators in
!!>>> full Hilbert space case
  subroutine atomic_make_ffmat()
     use control, only : ictqmc
     use control, only : norbs, ncfgs
     use m_fock, only : dec_basis, ind_basis
     use m_fock, only : fmat, evec

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! left Fock state
     integer :: left

! right Fock state
     integer :: right

! the sign change
     integer :: isgn

! evaluate F-matrix in the Fock basis
     do i=1,norbs
         do j=1,ncfgs
             right = dec_basis(j)
             if ( btest(right,i-1) .eqv. .true. ) then
                call atomic_make_c(i, right, left, isgn)
                k = ind_basis(left)
                fmat(k,j,i) = dble(isgn)
             endif ! back if ( btest(right,i-1) .eqv. .true. ) block
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

! rotate fmat from Fock basis to the atomic eigenvector basis
     if ( ictqmc == 1 ) then
         do i=1,norbs
             call atomic_tran_repr_real(ncfgs, fmat(:,:,i), evec)
         enddo ! over i={1,norbs} loop
     endif ! back if ( ictqmc == 1 ) block

     return
  end subroutine atomic_make_ffmat

!!>>> atomic_make_foccu: make occupancy for atomic eigenstates in the full
!!>>> Hilbert space case
  subroutine atomic_make_foccu()
     use constants, only : zero, one

     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis
     use m_fock, only : occu, evec

     implicit none

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! evaluate occupancy in the Fock basis
     occu = zero
     do ibas=1,ncfgs
         do iorb=1,norbs
             if ( bin_basis(iorb,ibas) == 1 ) then
                 occu(ibas,ibas) = occu(ibas,ibas) + one
             endif ! back if ( bin_basis(iorb,ibas ) == 1) block
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

! transform the occupancy from Fock basis to atomic eigenbasis
     call atomic_tran_repr_real(ncfgs, occu, evec)

     return
  end subroutine atomic_make_foccu

!!>>> atomic_make_fspin: make net Sz for atomic eigenstates in the full
!!>>> Hilbert space case
  subroutine atomic_make_fspin()
     use constants, only: zero, half

     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis
     use m_fock, only : spin, evec

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! evaluate spin moment in the Fock basis
     spin = zero
     do ibas=1,ncfgs
         do iorb=1,norbs
             if ( bin_basis(iorb,ibas) == 1 ) then
                 if ( mod(iorb,2) /= 0 ) then ! spin up
                     spin(ibas,ibas) = spin(ibas,ibas) + half
                 else                         ! spin down
                     spin(ibas,ibas) = spin(ibas,ibas) - half
                 endif ! back if ( mod(iorb,2) /= 0 ) block
             endif ! back if ( bin_basis(iorb,ibas ) == 1) block
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

! transform the net Sz from Fock basis to atomic eigenbasis
     call atomic_tran_repr_real(ncfgs, spin, evec)

     return
  end subroutine atomic_make_fspin

!!>>> atomic_make_fhmat: make atomic Hamiltonian in the full Hilbert space
  subroutine atomic_make_fhmat()
     use constants, only : one, czero, epst

     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis, dec_basis, ind_basis
     use m_fock, only : hmat
     use m_spmat, only : emat, umat

     implicit none

! local variables
! loop index
     integer :: i

! loop index for basis
     integer :: ibas, jbas

! loop index for orbital
     integer :: alpha, betta
     integer :: delta, gamma

! sign change due to fermion anti-commute relation
     integer :: sgn

! new atomic state after fermion operators act
     integer :: knew

! binary code form of an atomic state
     integer :: code(norbs)

! start to make Hamiltonian
! initialize hmat
     hmat = czero

! first, two fermion operator terms
     do jbas=1,ncfgs
         alploop: do alpha=1,norbs
             betloop: do betta=1,norbs

                 sgn = 0
                 knew = dec_basis(jbas)
                 code(1:norbs) = bin_basis(1:norbs,jbas)

! impurity level is too small
                 if ( abs( emat(alpha,betta) ) < epst ) CYCLE

! simulate one annihilation operator
                 if ( code(betta) == 1 ) then
                     do i=1,betta-1
                         if ( code(i) == 1 ) sgn = sgn + 1
                     enddo ! over i={1,betta-1} loop
                     code(betta) = 0

! simulate one creation operator
                     if ( code(alpha) == 0 ) then
                         do i=1,alpha-1
                             if ( code(i) == 1 ) sgn = sgn + 1
                         enddo ! over i={1,alpha-1} loop
                         code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                         knew = knew - 2**(betta-1)
                         knew = knew + 2**(alpha-1)
                         sgn  = mod(sgn,2)
! now ibas means the index for the new state
                         ibas = ind_basis(knew)
                         if ( ibas == 0 ) then
                             call s_print_error('atomic_make_fhmat','error while determining new state!')
                         endif ! back if ( ibas == 0 ) block
                         hmat(ibas,jbas) = hmat(ibas,jbas) + emat(alpha,betta) * (-one)**sgn
                     endif ! back if (code(alpha) == 0) block
                 endif ! back if (code(betta) == 1) block

             enddo betloop ! over betta={1,norbs} loop
         enddo alploop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop

! second, four fermion operator terms (coulomb interaction)
     do jbas=1,ncfgs
         alphaloop : do alpha=1,norbs
             bettaloop : do betta=1,norbs
                 gammaloop : do gamma=1,norbs
                     deltaloop : do delta=1,norbs
                         sgn  = 0
                         knew = dec_basis(jbas)
                         code(1:norbs) = bin_basis(1:norbs,jbas)

! applying Pauli principle
                         if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE

! U-matrix element is too small
                         if ( abs( umat(alpha,betta,delta,gamma) ) < epst ) CYCLE

! simulate two annihilation operators
                         if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) then
                             do i=1,gamma-1
                                 if ( code(i) == 1 ) sgn = sgn + 1
                             enddo ! over i={1,gamma-1} loop
                             code(gamma) = 0
                             do i=1,delta-1
                                 if ( code(i) == 1 ) sgn = sgn + 1
                             enddo ! over i={1,delta-1} loop
                             code(delta) = 0

! simulate two creation operator
                             if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) then
                                 do i=1,betta-1
                                     if ( code(i) == 1 ) sgn = sgn + 1
                                 enddo ! over i={1,betta-1} loop
                                 code(betta) = 1
                                 do i=1,alpha-1
                                     if ( code(i) == 1 ) sgn = sgn + 1
                                 enddo ! over i={1,alpha-1} loop
                                 code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                                 knew = knew - 2**(gamma-1) - 2**(delta-1)
                                 knew = knew + 2**(betta-1) + 2**(alpha-1)
                                 sgn  = mod(sgn,2)
! now ibas means the index for the new state
                                 ibas = ind_basis(knew)
                                 if ( ibas == 0 ) then
                                     call s_print_error('atomic_make_fhmat','error while determining new state')
                                 endif ! back if ( ibas == 0 ) block
                                 hmat(ibas,jbas) = hmat(ibas,jbas) + umat(alpha,betta,delta,gamma) * (-one)**sgn
                             endif ! back if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) block
                         endif ! back if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) block

                     enddo deltaloop ! over delta={1,norbs} loop
                 enddo gammaloop ! over gamma={1,norbs} loop
             enddo bettaloop ! over betta={1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop

     return
  end subroutine atomic_make_fhmat

!!
!! @sub atomic_diag_fhmat
!!
!! diagonalize atomic Hamiltonian directly by a lapack call
!!
  subroutine atomic_diag_fhmat()
     use constants, only : mystd

     use control, only : ncfgs

     use m_fock, only : hmat
     use m_fock, only : eval, evec

     implicit none

!! [body

     call s_eig_sy(ncfgs, ncfgs, real(hmat), eval, evec)
     write(mystd,'(4X,a)') 'eigenvalues and eigenvectors are built'

!! body]

     return
  end subroutine atomic_diag_fhmat

!!
!! @sub atomic_check_fhmat
!!
!! verify whether the atomic Hamiltonian is real
!!
  subroutine atomic_check_fhmat()
     use constants, only : eps6
     use constants, only : mystd

     use m_fock, only : hmat

     implicit none

!! [body

     if ( any( abs( aimag(hmat) ) > eps6 ) ) then
         call s_print_error('atomic_check_fhmat', &
             & 'atomic Hamiltonian is not real!')
     else
         write(mystd,'(4X,a)') 'atomic Hamiltonian is valid'
     endif ! back if ( any( abs( aimag(hmat) ) > eps6 ) ) block

!! body]

     return
  end subroutine atomic_check_fhmat
