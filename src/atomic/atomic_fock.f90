!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_make_ffmat
!!!           atomic_make_foccu
!!!           atomic_make_fspin
!!!           atomic_make_fhmat
!!! source  : atomic_fock.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/23/2024 by li huang (last modified)
!!! purpose : core subroutines for solving atomic eigenvalue problem in
!!!           the Fock space.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub atomic_make_ffmat
!!
!! construct annihilation operator matrix in the Fock space, and then
!! rotate it to the atomic eigenbasis
!!
  subroutine atomic_make_ffmat()
     use constants, only : mystd

     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis
     use m_fock, only : dec_basis
     use m_fock, only : ind_basis
     use m_fock, only : evec
     use m_fock, only : fmat

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

     ! left Fock state, < bra |
     integer :: left

     ! right Fock state, | ket >
     integer :: right

     ! the sign change
     integer :: isgn

!! [body

     ! evaluate annihilation operator matrix in the Fock basis
     do i=1,norbs
         write(mystd,'(4X,a,i2,a)') 'build f(alpha =', i, ') in Fock basis'
         do j=1,ncfgs
             ! get | ket >
             right = dec_basis(j)
             if ( btest(right,i-1) .eqv. .true. ) then
                 ! get < bra |
                 call atomic_make_c(i, right, left, isgn)
                 k = ind_basis(left)
                 !
                 ! evaluate < bra | c | ket >
                 fmat(k,j,i) = dble(isgn)
                 !
                 ! write the Fock states and the matrix elements
                 write(mystd,'(4X,a)', advance = 'no') '< bra | = '
                 write(mystd,'(*(i1))', advance = 'no') bin_basis(:,k)
                 write(mystd,'(2X,a)', advance = 'no') '| ket > = '
                 write(mystd,'(*(i1))', advance = 'no') bin_basis(:,j)
                 write(mystd,'(2X,a,i2)') 'value = ', isgn
             endif ! back if ( btest(right,i-1) .eqv. .true. ) block
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

     ! rotate annihilation operator matrix from the Fock basis to
     ! the atomic eigenbasis
     do i=1,norbs
         write(mystd,'(4X,a,i2,a)') 'rotate c(alpha =', i, ') to atomic eigenbasis'
         call atomic_tran_repr_real(ncfgs, fmat(:,:,i), evec)
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_make_ffmat

!!
!! @sub atomic_make_foccu
!!
!! construct density matrix (occupancy matrix) in the Fock space, and
!! then rotate it to the atomic eigenbasis
!!
  subroutine atomic_make_foccu()
     use constants, only : mystd
     use constants, only : zero, one

     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis
     use m_fock, only : evec
     use m_fock, only : occu

     implicit none

!! local variables
     ! loop index over orbits
     integer :: iorb

     ! loop index over Fock states
     integer :: ibas

!! [body

     ! evaluate density matrix in the Fock basis
     occu = zero
     !
     do ibas=1,ncfgs
         do iorb=1,norbs
             if ( bin_basis(iorb,ibas) == 1 ) then
                 occu(ibas,ibas) = occu(ibas,ibas) + one
             endif ! back if ( bin_basis(iorb,ibas ) == 1) block
         enddo ! over iorb={1,norbs} loop
         write(mystd,'(4X,a)', advance = 'no') '| ket > = '
         write(mystd,'(*(i1))', advance = 'no') bin_basis(:,ibas)
         write(mystd,'(2X,a,f5.2)') 'N = ', occu(ibas,ibas)
     enddo ! over ibas={1,ncfgs} loop

     ! try to transform the density matrix from the Fock basis
     ! to the atomic eigenbasis
     write(mystd,'(4X,a)') 'rotate density matrix to atomic eigenbasis'
     call atomic_tran_repr_real(ncfgs, occu, evec)

!! body]

     return
  end subroutine atomic_make_foccu

!!
!! @sub atomic_make_fspin
!!
!! construct magnetic moment (Sz) in the Fock space, and
!! then rotate it to the atomic eigenbasis
!!
  subroutine atomic_make_fspin()
     use constants, only : mystd
     use constants, only: zero, half

     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis
     use m_fock, only : evec
     use m_fock, only : spin

!! local variables
     ! loop index over orbits
     integer :: iorb

     ! loop index over Fock states
     integer :: ibas

!! [body

     ! evaluate magnetic moment in the Fock basis
     spin = zero
     !
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
         write(mystd,'(4X,a)', advance = 'no') '| ket > = '
         write(mystd,'(*(i1))', advance = 'no') bin_basis(:,ibas)
         write(mystd,'(2X,a,f5.2)') 'Sz = ', spin(ibas,ibas)
     enddo ! over ibas={1,ncfgs} loop

     ! try to transform the magnetic moment from the Fock basis
     ! to the atomic eigenbasis
     write(mystd,'(4X,a)') 'rotate magnetic moment to atomic eigenbasis'
     call atomic_tran_repr_real(ncfgs, spin, evec)

!! body]

     return
  end subroutine atomic_make_fspin

!!
!! @sub atomic_make_fhmat
!!
!! assemble atomic Hamiltonian in the Fock space
!!
  subroutine atomic_make_fhmat()
     use constants, only : dp
     use constants, only : mystd
     use constants, only : one, czero
     use constants, only : epst

     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis
     use m_fock, only : dec_basis
     use m_fock, only : ind_basis
     use m_fock, only : hmat

     use m_spmat, only : emat
     use m_spmat, only : umat

     implicit none

!! local variables
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

     ! matrix element of the atomic Hamiltonian
     real(dp) :: val

!! [body

     ! start to make Hamiltonian
     ! initialize hmat
     hmat = czero

     ! compute two fermion operators term
     ! it is f^{\dagger}_{\alpha} f_{\beta}
     write(mystd,'(4X,a)') 'compute two fermion operators term'
     !
     do jbas=1,ncfgs
         alploop: do alpha=1,norbs ! loop over creation operators
         betloop: do betta=1,norbs ! loop over annihilation operators

         ! retrieve the Fock state |jbas>
         sgn = 0
         knew = dec_basis(jbas)
         code = bin_basis(:,jbas)

         ! impurity level is too small
         if ( abs( emat(alpha,betta) ) < epst ) CYCLE

         ! simulate one annihilation operator, f_{\beta}
         if ( code(betta) == 1 ) then
             do i=1,betta-1
                 if ( code(i) == 1 ) sgn = sgn + 1
             enddo ! over i={1,betta-1} loop
             code(betta) = 0

             ! simulate one creation operator, f^{\dagger}_{\alpha}
             if ( code(alpha) == 0 ) then
                 do i=1,alpha-1
                     if ( code(i) == 1 ) sgn = sgn + 1
                 enddo ! over i={1,alpha-1} loop
                 code(alpha) = 1

                 ! determine the new Fock state, <ibas|
                 ! now ibas means the index for the new Fock state
                 knew = knew - 2**(betta-1) + 2**(alpha-1)
                 ibas = ind_basis(knew)
                 if ( ibas == 0 ) then
                     call s_print_error('atomic_make_fhmat', &
                         & 'error while determining new state!')
                 endif ! back if ( ibas == 0 ) block
                 !
                 ! determine the matrix element between the two Fock
                 ! states, i.e., <ibas| and |jbas>
                 sgn = mod(sgn,2)
                 val = emat(alpha,betta) * (-one)**sgn
                 !
                 ! setup the two fermion operators term
                 hmat(ibas,jbas) = hmat(ibas,jbas) + val
                 !
                 ! write the Fock states and the operators
                 write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha = ', alpha, ')'
                 write(mystd,'(2X,a,i2,a)') 'f(beta = ', betta, ')'
             endif ! back if (code(alpha) == 0) block
         endif ! back if (code(betta) == 1) block

         enddo betloop ! over betta={1,norbs} loop
         enddo alploop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop

     ! compute four fermion operators term (coulomb interaction)
     ! it is f^{\dagger}_{\alpha} f^{\dagger}_{\beta} f_{\delta} f_{\gamma}
     write(mystd,'(4X,a)') 'compute four fermion operators term'
     !
     do jbas=1,ncfgs
         alphaloop : do alpha=1,norbs ! loop over creation operators
         bettaloop : do betta=1,norbs ! loop over creation operators
         gammaloop : do gamma=1,norbs ! loop over annihilation operators
         deltaloop : do delta=1,norbs ! loop over annihilation operators

         ! retrieve the Fock state |jbas>
         sgn  = 0
         knew = dec_basis(jbas)
         code = bin_basis(:,jbas)

         ! applying Pauli principle
         if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE

         ! U-matrix element is too small
         if ( abs( umat(alpha,betta,delta,gamma) ) < epst ) CYCLE

         ! simulate two annihilation operators
         ! they are f_{\delta} f_{\gamma}
         if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) then
             do i=1,gamma-1
                 if ( code(i) == 1 ) sgn = sgn + 1
             enddo ! over i={1,gamma-1} loop
             code(gamma) = 0
             !
             do i=1,delta-1
                 if ( code(i) == 1 ) sgn = sgn + 1
             enddo ! over i={1,delta-1} loop
             code(delta) = 0

             ! simulate two creation operators
             ! they are f^{\dagger}_{\alpha} f^{\dagger}_{\beta}
             if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) then
                 do i=1,betta-1
                     if ( code(i) == 1 ) sgn = sgn + 1
                 enddo ! over i={1,betta-1} loop
                 code(betta) = 1
                 !
                 do i=1,alpha-1
                     if ( code(i) == 1 ) sgn = sgn + 1
                 enddo ! over i={1,alpha-1} loop
                 code(alpha) = 1

                 ! determine the new Fock state, <ibas|
                 ! now ibas means the index for the new Fock state
                 knew = knew - 2**(gamma-1) - 2**(delta-1)
                 knew = knew + 2**(betta-1) + 2**(alpha-1)
                 ibas = ind_basis(knew)
                 if ( ibas == 0 ) then
                     call s_print_error('atomic_make_fhmat', &
                         & 'error while determining new state')
                 endif ! back if ( ibas == 0 ) block
                 !
                 ! determine the matrix element between the two Fock
                 ! states, i.e., <ibas| and |jbas>
                 sgn = mod(sgn,2)
                 val = umat(alpha,betta,delta,gamma) * (-one)**sgn
                 !
                 ! setup the four fermion operators term
                 hmat(ibas,jbas) = hmat(ibas,jbas) + val
                 !
                 ! write the Fock states and the operators
                 write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha = ', alpha, ')'
                 write(mystd,'(2X,a,i2,a)', advance = 'no') 'f^+(beta = ', betta, ')'
                 write(mystd,'(2X,a,i2,a)', advance = 'no') 'f(delta = ', delta, ')'
                 write(mystd,'(2X,a,i2,a)') 'f(gamma = ', gamma, ')'
             endif ! back if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) block
         endif ! back if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) block

         enddo deltaloop ! over delta={1,norbs} loop
         enddo gammaloop ! over gamma={1,norbs} loop
         enddo bettaloop ! over betta={1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop

!! body]

     return
  end subroutine atomic_make_fhmat
