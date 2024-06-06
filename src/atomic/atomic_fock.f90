!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_make_ffmat
!!!           atomic_make_foccu
!!!           atomic_make_fspin
!!!           atomic_make_fhmat
!!!           atomic_diag_fhmat
!!!           atomic_check_fhmat
!!! source  : atomic_fock.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/31/2024 by li huang (last modified)
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
     ! actually, we want to get < bra | f | ket >
     do i=1,norbs
         write(mystd,'(4X,a,i2,a)') 'build f(alpha =', i, ') in Fock basis'
         do j=1,ncfgs
             ! get | ket >
             right = dec_basis(j)
             !
             ! simulate an annihilation operator
             if ( btest(right,i-1) .eqv. .true. ) then
                 ! get < bra |
                 call atomic_make_c(i, right, left, isgn)
                 k = ind_basis(left)
                 !
                 ! evaluate < bra | f | ket >
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
     ! the atomic eigenbasis. here evec are eigenvectors of hmat
     do i=1,norbs
         write(mystd,'(4X,a,i2,a)') 'rotate f(alpha =', i, ') to atomic eigenbasis'
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
     use constants, only : zero, one
     use constants, only : mystd

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
     ! note that it is diagonal
     occu = zero
     !
     do ibas=1,ncfgs
         do iorb=1,norbs
             if ( bin_basis(iorb,ibas) == 1 ) then
                 occu(ibas,ibas) = occu(ibas,ibas) + one
             endif ! back if ( bin_basis(iorb,ibas ) == 1) block
         enddo ! over iorb={1,norbs} loop
         !
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
     use constants, only: zero, half
     use constants, only : mystd

     use control, only : nband, norbs, ncfgs

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
     ! note that it is diagonal
     spin = zero
     !
     do ibas=1,ncfgs
         do iorb=1,norbs
             if ( bin_basis(iorb,ibas) == 1 ) then
                 if ( iorb <= nband ) then ! spin up
                     spin(ibas,ibas) = spin(ibas,ibas) + half
                 else                      ! spin down
                     spin(ibas,ibas) = spin(ibas,ibas) - half
                 endif ! back if ( iorb <= nband ) block
             endif ! back if ( bin_basis(iorb,ibas ) == 1) block
         enddo ! over iorb={1,norbs} loop
         !
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
     use constants, only : one
     use constants, only : czero
     use constants, only : epst
     use constants, only : mystd

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

     ! loop index for Fock state
     integer :: ibas, jbas

     ! loop index for orbital
     integer :: alpha, betta
     integer :: delta, gamma

     ! sign change due to fermion anti-commutation relation
     integer :: isgn

     ! new Fock state after a fermion operator acts
     integer :: knew

     ! binary form of a Fock state
     integer :: code(norbs)

     ! matrix element of the atomic Hamiltonian
     real(dp) :: val

!! [body

     ! start to make atomic Hamiltonian
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
         isgn = 0
         knew = dec_basis(jbas)
         code = bin_basis(:,jbas)

         ! impurity level is too small
         if ( abs( emat(alpha,betta) ) < epst ) CYCLE

         ! simulate one annihilation operator, f_{\beta}
         if ( code(betta) == 1 ) then
             do i=1,betta-1
                 if ( code(i) == 1 ) isgn = isgn + 1
             enddo ! over i={1,betta-1} loop
             code(betta) = 0

             ! simulate one creation operator, f^{\dagger}_{\alpha}
             if ( code(alpha) == 0 ) then
                 do i=1,alpha-1
                     if ( code(i) == 1 ) isgn = isgn + 1
                 enddo ! over i={1,alpha-1} loop
                 code(alpha) = 1

                 ! determine the new Fock state, <ibas|
                 ! now ibas means the index for the new Fock state
                 knew = knew - 2**(betta-1) + 2**(alpha-1)
                 ibas = ind_basis(knew)
                 if ( ibas == 0 ) then
                     call s_print_error('atomic_make_fhmat', &
                         & 'error while determining new Fock state!')
                 endif ! back if ( ibas == 0 ) block
                 !
                 ! determine the matrix element between the two Fock
                 ! states, i.e., <ibas| and |jbas>
                 isgn = mod(isgn,2)
                 val = emat(alpha,betta) * (-one)**isgn
                 !
                 ! setup the two fermion operators term
                 hmat(ibas,jbas) = hmat(ibas,jbas) + val
                 !
                 ! write the combination of operators
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
         isgn  = 0
         knew = dec_basis(jbas)
         code = bin_basis(:,jbas)

         ! applying Pauli principle
         if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE

         ! U matrix element is too small
         if ( abs( umat(alpha,betta,delta,gamma) ) < epst ) CYCLE

         ! simulate two annihilation operators
         ! they are f_{\delta} f_{\gamma}
         if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) then
             do i=1,gamma-1
                 if ( code(i) == 1 ) isgn = isgn + 1
             enddo ! over i={1,gamma-1} loop
             code(gamma) = 0
             !
             do i=1,delta-1
                 if ( code(i) == 1 ) isgn = isgn + 1
             enddo ! over i={1,delta-1} loop
             code(delta) = 0

             ! simulate two creation operators
             ! they are f^{\dagger}_{\alpha} f^{\dagger}_{\beta}
             if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) then
                 do i=1,betta-1
                     if ( code(i) == 1 ) isgn = isgn + 1
                 enddo ! over i={1,betta-1} loop
                 code(betta) = 1
                 !
                 do i=1,alpha-1
                     if ( code(i) == 1 ) isgn = isgn + 1
                 enddo ! over i={1,alpha-1} loop
                 code(alpha) = 1

                 ! determine the new Fock state, <ibas|
                 ! now ibas means the index for the new Fock state
                 knew = knew - 2**(gamma-1) - 2**(delta-1)
                 knew = knew + 2**(betta-1) + 2**(alpha-1)
                 ibas = ind_basis(knew)
                 if ( ibas == 0 ) then
                     call s_print_error('atomic_make_fhmat', &
                         & 'error while determining new Fock state')
                 endif ! back if ( ibas == 0 ) block
                 !
                 ! determine the matrix element between the two Fock
                 ! states, i.e., <ibas| and |jbas>
                 isgn = mod(isgn,2)
                 val = umat(alpha,betta,delta,gamma) * (-one)**isgn
                 !
                 ! setup the four fermion operators term
                 hmat(ibas,jbas) = hmat(ibas,jbas) + val
                 !
                 ! write the combination of the operators
                 write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha = ', alpha, ')'
                 write(mystd,'(2X,a,i2,a)', advance = 'no') 'f^+(beta = ', betta, ')'
                 write(mystd,'(2X,a,i2,a)', advance = 'no') 'f(delta = ', delta, ')'
                 write(mystd,'(2X,a,i2,a)') 'f(gamma = ', gamma, ')'
             endif ! back if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) block
         endif ! back if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) block

         enddo deltaloop ! over delta={1,norbs} loop
         enddo gammaloop ! over gamma={1,norbs} loop
         enddo bettaloop ! over betta={1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop

     !call atomic_test_ad()

!! body]

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



  subroutine atomic_find_subspace(sind, find, sector_size, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: sector_size(ncfgs)
     integer, intent(in) :: sector_basis(ncfgs,ncfgs)

     integer :: m
     integer :: n

     sind = 0

     m1loop: do m=1,ncfgs
         do n=1,sector_size(m)
             if ( sector_basis(m,n) == find ) then
                 sind = m
                 EXIT m1loop
             endif
         enddo
     enddo m1loop

     call s_assert(sind /= 0)

     return
  end subroutine atomic_find_subspace

  subroutine atomic_test_ad()
     use constants, only : dp
     use constants, only : mystd

     use control, only : nband, norbs, ncfgs

     use m_fock, only : hmat
     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

     integer, parameter :: max_mapping = 200000

     integer :: i
     integer :: j
     integer :: m
     integer :: jold, jnew, iorb
     integer :: ia, ib, isgn
     integer :: iup, idn
     integer :: HA, HL, HU

     integer, allocatable :: sector_size(:), sector_size_(:)
     integer, allocatable :: sector_basis(:,:), sector_basis_(:,:)
     integer, allocatable :: Mup(:,:)
     integer, allocatable :: Mdn(:,:)

     allocate(sector_size(ncfgs))
     allocate(sector_size_(ncfgs))
     allocate(sector_basis(ncfgs,ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))
     allocate(Mup(max_mapping,2))
     allocate(Mdn(max_mapping,2))

     ! initialization
     sector_size = 1
     sector_basis = 0
     do i=1,ncfgs
         sector_basis(i,1) = i
     enddo

     ! phase 1
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(hmat(i,j)) > 0.0_dp ) then

                 call atomic_find_subspace(ib, i, sector_size, sector_basis)
                 call atomic_find_subspace(ia, j, sector_size, sector_basis)

                 if ( ia /= ib ) then
                     call s_assert(sector_size(ib) >= 1)
                     call s_assert(sector_size(ia) >= 1)
                     !
                     do m=1,sector_size(ia)
                         sector_basis(ib, sector_size(ib) + m) = sector_basis(ia,m)
                     enddo
                     !
                     sector_size(ib) = sector_size(ib) + sector_size(ia)
                     sector_size(ia) = 0
                     sector_basis(ia,:) = 0
                 endif

             endif
         enddo
     enddo

     ! print subspaces
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a,*(i6))') 'basis :', sector_basis(i,1:sector_size(i))
             print *
         endif
     enddo

     sector_size_ = 0
     sector_basis_ = 0
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             sector_size_(m) = sector_size(i)
             sector_basis_(m,:) = sector_basis(i,:)
         endif
     enddo
     sector_size = sector_size_
     sector_basis = sector_basis_

     ! create mapping
     iup = 0
     idn = 0
     Mup = 0
     Mdn = 0

     DO iorb = 1,2

     do i=1,ncfgs
         call atomic_find_subspace(ia, i, sector_size, sector_basis)

         if ( bin_basis(iorb,i) == 0 ) then
             jold = dec_basis(i)
             call atomic_make_cdagger(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call atomic_find_subspace(ib, j, sector_size, sector_basis)

             iup = iup + 1
             call s_assert(iup <= max_mapping)
             Mup(iup,1) = ia
             Mup(iup,2) = ib
         endif

         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call atomic_find_subspace(ib, j, sector_size, sector_basis)

             idn = idn + 1
             call s_assert(idn <= max_mapping)
             Mdn(idn,1) = ia
             Mdn(idn,2) = ib
         endif
     enddo

     ENDDO

     print *, 'number of Mup:', iup
     print *, 'number of Mdn:', idn
     !
     do i=1,iup
         print *, Mup(i,:), Mdn(i,:)
     enddo
     STOP


     do i=1,max_mapping
         HA = Mup(i,1)
         HL = Mup(i,1)
         HU = Mup(i,2)

         if ( HL /= 0 .and. HU /= 0 ) then
         call zigzag(1, HA, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
         endif
     enddo

     ! print subspaces
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a,*(i6))') 'basis :', sector_basis(i,1:sector_size(i))
             print *
         endif
     enddo
     !STOP

     return
  end subroutine atomic_test_ad

recursive &
  subroutine zigzag(up_or_down, HA, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: up_or_down
     integer, intent(in) :: HA
     integer, intent(in) :: HL
     integer, intent(in) :: HU
     integer, intent(in) :: max_mapping
     integer, intent(inout) :: sector_size(ncfgs)
     integer, intent(inout) :: sector_basis(ncfgs,ncfgs)
     integer, intent(inout) :: Mup(max_mapping,2)
     integer, intent(inout) :: Mdn(max_mapping,2)

     integer :: i, j
     integer :: HB

     if ( up_or_down == 1 ) then
         do i=1,max_mapping
             if ( Mup(i,1) == HA ) then
                 HB = Mup(i,2)
                 Mup(i,:) = 0

                 if ( HB /= HU ) then
                     print *, 'merge up:', HA, HB, HU
                     do j=1,sector_size(HB)
                         sector_basis(HU, sector_size(HU) + j) = sector_basis(HB,j)
                     enddo
                     sector_size(HU) = sector_size(HU) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(2, HB, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
             endif
         enddo 
     else
         do i=1,max_mapping
             if ( Mdn(i,1) == HA ) then
                 HB = Mdn(i,2)
                 Mdn(i,:) = 0

                 if ( HB /= HL ) then
                     print *, 'merge dn:', HA, HB, HL
                     do j=1,sector_size(HB)
                         sector_basis(HL, sector_size(HL) + j) = sector_basis(HB,j)
                     enddo
                     sector_size(HL) = sector_size(HL) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(1, HB, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
             endif
         enddo 
     endif

     return
  end subroutine zigzag
