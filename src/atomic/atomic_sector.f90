!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_sfmat
!!!           atomic_make_shmat
!!!           atomic_diag_shmat
!!!           atomic_make_sectors
!!! source  : atomic_sector.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : implement the sector algorithm, calculate the F-matrix and
!!!           build the atomic Hamiltonian sector-by-sector
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_make_sfmat: build F-matrix (fmat) for good quantum numbers
!!>>> (GQNs) algorithm
  subroutine atomic_make_sfmat()
     use constants, only : zero

     use control, only : norbs
     use m_fock, only : dec_basis, ind_basis
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_fmat

     implicit none

! local variables
! loop index for orbital
     integer :: iorb

! loop index for annihilation and creation operators
     integer :: ityp

! loop index for sector
     integer :: isec
     integer :: jsec

! loop index for basis
     integer :: ibas
     integer :: jbas

! sign change due to commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold
     integer :: jnew

! loop over all the sectors
     do isec=1,nsectors
! loop over all the orbitals
         do iorb=1,norbs
! loop over the creation and annihilation fermion operators
             do ityp=0,1

! get the next sector: jsec
                 jsec = sectors(isec)%next(iorb,ityp)
                 if ( jsec == -1 ) CYCLE

! allocate memory for fmat and then initialize it
                 sectors(isec)%fmat(iorb,ityp)%n = sectors(jsec)%ndim
                 sectors(isec)%fmat(iorb,ityp)%m = sectors(isec)%ndim
                 call cat_alloc_fmat(sectors(isec)%fmat(iorb,ityp))
                 sectors(isec)%fmat(iorb,ityp)%val = zero

! loop over the basis for the isec-th sector
                 do jbas=1,sectors(isec)%ndim
                     jold = dec_basis(sectors(isec)%basis(jbas))

! apply creation fermion operator
                     if ( ityp == 1 .and. ( btest(jold, iorb-1) .eqv. .false. ) ) then
                         call atomic_make_cdagger(iorb, jold, jnew, isgn)
! apply annihilation fermion operator
                     else if ( ityp == 0 .and. ( btest(jold, iorb-1) .eqv. .true. ) ) then
                         call atomic_make_c(iorb, jold, jnew, isgn)
                     else
                         CYCLE
                     endif ! back if ( ityp == 1 .and. ( btest(jold, iorb-1) .eqv. .false. ) ) block

! loop over the basis for the jsec-th sector
                     do ibas=1,sectors(jsec)%ndim
                         if ( sectors(jsec)%basis(ibas) == ind_basis(jnew) ) then
! build matrix element for F-matrix
                             sectors(isec)%fmat(iorb,ityp)%val(ibas,jbas) = dble(isgn)
                             EXIT
                         endif ! back if ( sectors(jsec)%basis(ibas) == ind_basis(jnew) ) block
                     enddo ! over ibas={1,sectors(jsec)%ndim} loop
                 enddo ! over jbas={1,sectors(isec)%ndim} loop

! roate fmat to atomic eigenstates basis
                 call atomic_tran_fmat(sectors(jsec)%ndim, &
                                       sectors(isec)%ndim, &
                                       sectors(jsec)%evec, &
                                       sectors(isec)%fmat(iorb,ityp)%val, &
                                       sectors(isec)%evec)

             enddo ! over ityp={0,1} loop
         enddo ! over iorb={1,norbs} loop
     enddo ! over isec={1,nsectors} loop

     return
  end subroutine atomic_make_sfmat

!!>>> atomic_make_shmat: make Hamiltonian for each sector one by one
  subroutine atomic_make_shmat()
     use constants, only : one, epst, czero

     use control, only : norbs
     use m_fock, only : bin_basis, dec_basis, ind_basis
     use m_spmat, only : emat, umat
     use m_sector, only : nsectors, sectors

     implicit none

! local variables
! loop index
     integer :: i

! loop index for sector
     integer :: isec

! loop index for basis
     integer :: ibas
     integer :: jbas

! loop index for orbital
     integer :: alpha, betta
     integer :: delta, gamma

! sign change due to fermion anti-commute relation
     integer :: isgn

! new basis state after four fermion operation
     integer :: knew

! binary form of a Fock state
     integer :: code(norbs)

! loop over all sectors
     do isec=1,nsectors
         sectors(isec)%hmat = czero

! two fermions term
!-------------------------------------------------------------------------
         do jbas=1,sectors(isec)%ndim
             alploop: do alpha=1,norbs
                 betloop: do betta=1,norbs

                     isgn = 0
                     knew = dec_basis(sectors(isec)%basis(jbas))
                     code(1:norbs) = bin_basis(1:norbs,sectors(isec)%basis(jbas))

! impurity level is too small
                     if ( abs(emat(alpha,betta)) < epst ) CYCLE

! simulate one annihilation operator
                     if ( code(betta) == 1 ) then
                         do i=1,betta-1
                             if ( code(i) == 1 ) isgn = isgn + 1
                         enddo ! over i={1,betta-1} loop
                         code(betta) = 0

! simulate one creation operator
                         if ( code(alpha) == 0 ) then
                             do i=1,alpha-1
                                 if ( code(i) == 1 ) isgn = isgn + 1
                             enddo ! over i={1,alpha-1} loop
                             code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                             knew = knew - 2**(betta-1)
                             knew = knew + 2**(alpha-1)
                             isgn = mod(isgn,2)

! now ind_basis(knew) means the index of new Fock state
                             if ( ind_basis(knew) == 0 ) then
                                 call s_print_error('atomic_make_shmat','error while determining new state!')
                             endif ! back if ( ind_basis(knew) == 0 ) block
                             do ibas=1,sectors(isec)%ndim
                                 if ( sectors(isec)%basis(ibas) == ind_basis(knew) ) then
                                     sectors(isec)%hmat(ibas,jbas) = &
                                     sectors(isec)%hmat(ibas,jbas) + &
                                     emat(alpha,betta) * (-one)**isgn
                                 endif ! back if ( sectors(isec)%basis(ibas) == ind_basis(knew) ) block
                             enddo ! over ibas={1,sectors(isec)%ndim} loop
                         endif ! back if ( code(alpha) == 0 ) block
                     endif ! back if ( code(betta_ == 1 ) block

                 enddo betloop ! over betta={1,norbs} loop
             enddo alploop ! over alpha={1,norbs} loop
         enddo ! over jbas={1,sectors(isect)%ndim} loop

! four fermions term
!-------------------------------------------------------------------------
         do jbas=1,sectors(isec)%ndim
             alphaloop: do alpha=1,norbs
                 bettaloop: do betta=1,norbs
                     gammaloop: do gamma=1,norbs
                         deltaloop: do delta=1,norbs

                             isgn = 0
                             knew = dec_basis(sectors(isec)%basis(jbas))
                             code(1:norbs) = bin_basis(1:norbs,sectors(isec)%basis(jbas))

! applying Pauli principle
                             if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE

! U-matrix element is too small
                             if ( abs(umat(alpha,betta,delta,gamma)) < epst ) CYCLE

! simulate two annihilation operators
                             if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) then
                                 do i=1,gamma-1
                                     if ( code(i) == 1 ) isgn = isgn + 1
                                 enddo ! over i={1,gamma-1} loop
                                 code(gamma) = 0
                                 do i=1,delta-1
                                     if ( code(i) == 1 ) isgn = isgn + 1
                                 enddo ! over i={1,delta-1} loop
                                 code(delta) = 0

! simulate two creation operators
                                 if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) then
                                     do i=1,betta-1
                                         if ( code(i) == 1 ) isgn = isgn + 1
                                     enddo ! over i={1,betta-1} loop
                                     code(betta) = 1
                                     do i=1,alpha-1
                                         if ( code(i) == 1 ) isgn = isgn + 1
                                     enddo ! over i={1,alpha-1} loop
                                     code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                                     knew = knew - 2**(gamma-1) - 2**(delta-1)
                                     knew = knew + 2**(betta-1) + 2**(alpha-1)
                                     isgn = mod(isgn,2)

! now ind_basis(knew) means the index of new Fock state
                                     if ( ind_basis(knew) == 0 ) then
                                         call s_print_error('atomic_make_shmat','error while determining new state!')
                                     endif ! back if ( ind_basis(knew) == 0 ) block
                                     do ibas=1,sectors(isec)%ndim
                                         if ( sectors(isec)%basis(ibas) == ind_basis(knew) ) then
                                             sectors(isec)%hmat(ibas,jbas) = &
                                             sectors(isec)%hmat(ibas,jbas) + &
                                             umat(alpha,betta,delta,gamma) * (-one)**isgn
                                         endif ! back if ( sectors(isec)%basis(ibas) == ind_basis(knew) ) block
                                     enddo ! over ibas={1,sectors(isec)%ndim} loop
                                 endif ! back if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) block
                             endif ! back if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) block

                         enddo deltaloop ! over delta={1,norbs} loop
                     enddo gammaloop ! over gamma={1,norbs} loop
                 enddo bettaloop ! over betta={1,norbs} loop
             enddo alphaloop ! over alpha={1,norbs} loop
         enddo ! over jbas={1,sectors(isect)%ndim} loop

     enddo ! over isec={1,nsectors} loop

     return
  end subroutine atomic_make_shmat

!!>>> atomic_diag_shmat: diagonalize the Hamiltonian for each sector
  subroutine atomic_diag_shmat()
     use constants, only : dp

     use m_sector, only : nsectors, sectors

     implicit none

! local variables
! loop index
     integer :: i

! dummy array
     real(dp), allocatable :: hmat(:,:)

     do i=1,nsectors
         allocate(hmat(sectors(i)%ndim,sectors(i)%ndim))
         hmat = real( sectors(i)%hmat )
         call s_eig_sy( sectors(i)%ndim, &
                        sectors(i)%ndim, &
                        hmat,            &
                        sectors(i)%eval, &
                        sectors(i)%evec )
         deallocate(hmat)
     enddo ! over i={1,nsectors} loop

     return
  end subroutine atomic_diag_shmat

!!>>> atomic_make_sectors: determine all the sectors with the good quantum
!!>>> numbers algorithm
!!>>> a sector consists of some many-particle Fock states labeled by
!!>>> good quantum numbers
  subroutine atomic_make_sectors()
     use constants, only : zero

     use control, only : ictqmc
     use control, only : nband, norbs, ncfgs
     use control, only : nmini, nmaxi
     use m_fock, only : dim_sub_n, bin_basis
     use m_sector, only : max_dim_sect, ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_sector
     use m_sector, only : cat_alloc_sectors

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k
     integer :: l

! total electrons
     integer :: my_ntot

! Sz value
     integer :: my_sz

! Jz value
     integer :: my_jz

! PS value
     integer :: my_ps

! a counter
     integer :: counter

! index of Fock basis
     integer :: ibasis

! number of sectors
     integer :: nsect

! which sector point to
     integer :: which_sect

! can point to next sector
     logical :: can

! the Sz and Jz values for each orbital
     integer :: orb_good_sz(norbs)
     integer :: orb_good_jz(norbs)

! good quantum number N, Sz, Jz, and PS for each Fock state
     integer :: fock_good_ntot(ncfgs)
     integer :: fock_good_sz(ncfgs)
     integer :: fock_good_jz(ncfgs)
     integer :: fock_good_ps(ncfgs)

! good quantum number N, Sz, Jz, and PS for each sector
     integer :: sect_good_ntot(ncfgs)
     integer :: sect_good_sz(ncfgs)
     integer :: sect_good_jz(ncfgs)
     integer :: sect_good_ps(ncfgs)

! dimension of each sector
     integer :: ndims(ncfgs)

! a temp binary form of Fock basis
     integer :: code(norbs)

! sector basis index
! the first index: dimension size of the sector
! the second index: the index of sector
     integer, allocatable :: sector_basis(:,:)

! dummy variable
     integer :: sum_dim

! initialize some variables
     sect_good_ntot = 0
     sect_good_sz = 0
     sect_good_jz = 0
     sect_good_ps = 0

! allocate memory
     allocate(sector_basis(ncfgs,ncfgs))

! make orb_good_sz and orb_good_jz
!-------------------------------------------------------------------------
     orb_good_sz = 0
     call atomic_make_gsz(orb_good_sz)

! jz only valid for nband==3, 5, 7
     orb_good_jz = 0
     if ( nband == 3 .or. nband == 5 .or. nband == 7 ) then
         call atomic_make_gjz(orb_good_jz)
     endif ! back if ( nband == 3 .or. nband == 5 .or. nband == 7 ) block

! build good quantum numbers for each Fock state
!-------------------------------------------------------------------------
     counter = 0
     fock_good_ntot = 0
     fock_good_sz = 0
     fock_good_jz = 0
     fock_good_ps = 0
! loop over all number of total electrons
     do i=0,norbs
! loop over each state
         do j=1,dim_sub_n(i)
! here counter denotes the index of Fock state
             counter = counter + 1
! build N
             fock_good_ntot(counter) = i
! build Sz
             my_sz = 0
             do k=1,norbs
                 my_sz = my_sz + orb_good_sz(k) * bin_basis(k,counter)
             enddo ! over k={1,norbs} loop
             fock_good_sz(counter) = my_sz
! build Jz
              my_jz = 0
              do k=1,norbs
                  my_jz = my_jz + orb_good_jz(k) * bin_basis(k,counter)
              enddo ! over k={1,norbs} loop
              fock_good_jz(counter) = my_jz
! build PS number
             do k=1,nband
                 fock_good_ps(counter) = &
                 fock_good_ps(counter) + (2**k) * &
                     (bin_basis(2*k-1,counter) - bin_basis(2*k,counter))**2
             enddo ! over k={1,nband} loop
         enddo ! over j={1,dim_sub_n(i)} loop
     enddo ! over i={0,norbs} loop

! loop over all the Fock states to determine sectors
!-------------------------------------------------------------------------
     nsect = 0
     ndims = 0
     sector_basis = 0
     do i=1,ncfgs
         my_ntot = fock_good_ntot(i)

! truncate the occupancy according to nmini and nmaxi
         if ( my_ntot < nmini  .or. my_ntot > nmaxi ) CYCLE

         if ( ictqmc == 3 .or. ictqmc == 4 ) then
             my_sz = fock_good_sz(i)
         endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
         if ( ictqmc == 4 ) then
             my_ps = fock_good_ps(i)
         endif ! back if ( ictqmc == 4 ) block
         if ( ictqmc == 5 ) then
             my_jz = fock_good_jz(i)
         endif ! back if ( ictqmc == 5 ) block

! determine the first sector
         if ( nsect == 0 ) then
             sect_good_ntot(1) = my_ntot
             if ( ictqmc == 3 .or. ictqmc == 4 ) then
                 sect_good_sz(1) = my_sz
             endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
             if (ictqmc == 4) then
                 sect_good_ps(1) = my_ps
             endif ! back if ( ictqmc == 4 ) block
             if ( ictqmc == 5 ) then
                 sect_good_jz(1) = my_jz
             endif ! back if ( ictqmc == 5 ) block

             nsect = nsect + 1
             ndims(1) = ndims(1) + 1
             sector_basis(ndims(1),1) = i
         else
! loop over the exists sectors
             which_sect = -1
             do j=1,nsect
! compare the current state with existing sectors
                 select case (ictqmc)
                     case (2)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             which_sect = j; EXIT
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block

                     case (3)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_sz(j) == my_sz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_good_sz(j) == my_sz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block

                     case (4)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_sz(j) == my_sz ) then
                                 if ( sect_good_ps(j) == my_ps) then
                                     which_sect = j; EXIT
                                 endif ! back if ( sect_good_ps(j) == my_ps) block
                             endif ! back if ( sect_good_sz(j) == my_sz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block

                     case (5)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_jz(j) == my_jz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_good_jz(j) == my_jz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block

                 end select
             enddo ! over j={1,nsect} loop

! we can not assign the current state into any existing sectors, so we
! have to define a new sector
             if ( which_sect == -1 ) then
                 nsect = nsect + 1
                 sect_good_ntot(nsect) = my_ntot
                 if ( ictqmc == 3 .or. ictqmc == 4 ) then
                     sect_good_sz(nsect) = my_sz
                 endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
                 if ( ictqmc == 4 ) then
                     sect_good_ps(nsect) = my_ps
                 endif ! back if ( ictqmc == 4 ) block
                 if ( ictqmc == 5 ) then
                     sect_good_jz(nsect) = my_jz
                 endif ! back if ( ictqmc == 5 ) block
                 ndims(nsect) = ndims(nsect) + 1
                 sector_basis(ndims(nsect),nsect) = i
! we assign the current state to one of the old sectors
             else
                 ndims(which_sect) = ndims(which_sect) + 1
                 sector_basis(ndims(which_sect),which_sect) = i
             endif ! back if ( which_sect == -1 ) block
         endif ! back if ( nsect == 0 ) block
     enddo ! over i={1,ncfgs} loop

! after we know how many sectors and the dimension of each sector are there,
! we can allocate memory for global variables for sectors
!-------------------------------------------------------------------------
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect
     call cat_alloc_sectors()
! now we will build each sector
     counter = 1
     do i=1,nsect
         sectors(i)%ndim = ndims(i)
         sectors(i)%nele = sect_good_ntot(i)
         sectors(i)%sz   = sect_good_sz(i)
         sectors(i)%jz   = sect_good_jz(i)
         sectors(i)%ps   = sect_good_ps(i)
         sectors(i)%nops = norbs
         sectors(i)%istart = counter
         counter = counter + ndims(i)
! allocate memory for each sector
         call cat_alloc_sector( sectors(i) )
! set basis for each sector
         do j=1,ndims(i)
             sectors(i)%basis(j) = sector_basis(j,i)
         enddo ! over j={1,ndims(i)} loop
     enddo ! over i={1,nsect} loop

! make next sector index
!-------------------------------------------------------------------------
! loop over all the sectors
     do i=1,nsectors
! loop over all the orbtials
         do j=1,norbs
! loop over creation and annihilation fermion operators
             do k=0,1
                 which_sect = -1
! we should check each state in this sector
                 can = .false.
                 do l=1,sectors(i)%ndim
                     ibasis = sectors(i)%basis(l)
! for creation fermion operator
                     if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) then
                         code = bin_basis(:,ibasis)
                         can = .true.
                         EXIT
! for annihilation fermion operator
                     else if ( k==0 .and. bin_basis(j, ibasis) == 1 ) then
                         code = bin_basis(:,ibasis)
                         can = .true.
                         EXIT
                     endif ! back if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) block
                 enddo ! over l={1,sectors(i)%ndim} loop

                 if ( can .eqv. .true. ) then
                     select case (ictqmc)
                         case (2)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     which_sect = l; EXIT
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (3)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_sz = sect_good_sz(i) + orb_good_sz(j)
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_sz = sect_good_sz(i) - orb_good_sz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_sz(l) == my_sz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_good_sz(l) == my_sz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (4)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_sz = sect_good_sz(i) + orb_good_sz(j)
                                 code(j) = 1
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_sz   = sect_good_sz(i) - orb_good_sz(j)
                                 code(j) = 0
                             endif ! back if ( k == 1 ) block
! calculate new PS number
                             my_ps = 0
                             do l=1,nband
                                 my_ps = my_ps + (2**l) * ( code(2*l-1) - code(2*l) )**2
                             enddo ! over l={1,nband} loop
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_sz(l) == my_sz ) then
                                         if ( sect_good_ps(l) == my_ps) then
                                             which_sect = l; EXIT
                                         endif ! back if ( sect_good_ps(l) == my_ps) block
                                     endif ! back if ( sect_good_sz(l) == my_sz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (5)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_jz = sect_good_jz(i) + orb_good_jz(j)
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_jz = sect_good_jz(i) - orb_good_jz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_jz(l) == my_jz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_good_jz(l) == my_jz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                     end select ! back select case (ictqmc) block
                 endif  ! back if ( can == .true. ) block
                 sectors(i)%next(j,k) = which_sect
             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,nsectors} loop

! calculate the maximum and average dimensions of sectors
!-------------------------------------------------------------------------
     max_dim_sect = maxval(ndims)
     sum_dim = 0
     do i=1,nsectors
         sum_dim = sum_dim + sectors(i)%ndim
     enddo
     ave_dim_sect = real(sum_dim) / real(nsectors)

! dump sector information for reference
!-------------------------------------------------------------------------
     call atomic_dump_sector(sect_good_ntot, sect_good_sz, sect_good_ps, sect_good_jz)

! deallocate memory
     deallocate(sector_basis)

     return
  end subroutine atomic_make_sectors

!!
!! @sub atomic_check_shmat
!!
!! verify whether the atomic Hamiltonian is real
!!
  subroutine atomic_check_shmat()
     use constants, only : eps6
     use constants, only : mystd

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! local variables
     ! loop index
     integer :: i

!! [body

     ! we should go through every subspace
     do i=1,nsectors
         if ( any( abs( aimag(sectors(i)%hmat) ) > eps6 ) ) then
             call s_print_error('atomic_check_shmat', &
                 & 'atomic Hamiltonian is not real!')
         else
             write(mystd,'(4X,a,i4,2X,a)') 'subspace: ', i, 'is valid'
         endif ! back if ( any( abs( aimag(sectors(i)%hmat) ) > eps6 ) ) block
     enddo ! over i={1,nsectors} loop

!! body]

     return
  end subroutine atomic_check_shmat
