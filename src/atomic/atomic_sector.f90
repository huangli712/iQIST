!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_make_sectors
!!!           atomic_make_sfmat
!!!           atomic_make_shmat
!!!           atomic_diag_shmat
!!!           atomic_check_shmat
!!! source  : atomic_sector.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/31/2024 by li huang (last modified)
!!! purpose : try to implement the subspace diagonalization algorithm.
!!!           it contains some subroutines to construct the atomic
!!!           Hamiltonian subspace by subspace, diagonalize it, and then
!!!           calculate the annihilation operator matrix.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub atomic_make_sectors
!!
!! determine all the subspaces with the good quantum numbers algorithm.
!! the subspace should consist of some many-particle Fock states labeled
!! by good quantum numbers
!!
  subroutine atomic_make_sectors()
     use constants, only : zero
     use constants, only : mystd

     use control, only : ictqmc
     use control, only : nband, norbs, ncfgs
     use control, only : nmini, nmaxi

     use m_fock, only : dim_sub_n
     use m_fock, only : bin_basis

     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_sector
     use m_sector, only : cat_alloc_sectors

     implicit none

!! local parameters
     ! maximum number of subspaces 
     integer, parameter :: max_num_sect = 400

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k
     integer :: l

     ! good quantum number: total electrons N
     integer :: my_ntot

     ! good quantum number: Sz
     integer :: my_sz

     ! good quantum number: Jz
     integer :: my_jz

     ! good quantum number: PS
     integer :: my_ps

     ! index of Fock state
     integer :: ibasis

     ! number of subspaces (sectors)
     integer :: nsect

     ! index of selected subspace
     integer :: which_sect

     ! dummy variable
     integer :: val

     ! can point to next subspace (sector)
     logical :: can

     ! Sz, Jz, and PS for all orbitals
     integer :: orb_sz(norbs)
     integer :: orb_jz(norbs)
     integer :: orb_ps(nband)

     ! N, Sz, Jz, and PS for all Fock states
     integer :: fock_ntot(ncfgs)
     integer :: fock_sz(ncfgs)
     integer :: fock_jz(ncfgs)
     integer :: fock_ps(ncfgs)

     ! N, Sz, Jz, and PS for all subspaces
     integer :: sect_ntot(ncfgs)
     integer :: sect_sz(ncfgs)
     integer :: sect_jz(ncfgs)
     integer :: sect_ps(ncfgs)

     ! dimension for subspaces
     integer :: ndims(ncfgs)

     ! binary form of Fock state
     integer :: code(norbs)

     ! global indices of Fock states of subspaces
     !
     ! the first index: local index of Fock state in the given subspace
     ! the second index: index of subspace
     !
     ! for example, sector_basis(j,i) gives the global index of Fock
     ! state for the j-th basis function in the i-th subspace
     integer, allocatable :: sector_basis(:,:)

!! [body

     ! allocate memory
     if (ncfgs > max_num_sect) then
         nsect = max_num_sect
     else
         nsect = ncfgs
     endif ! back if (ncfgs > max_num_sect) block
     !
     write(mystd,'(4X,a,i4)') 'maximum number of subspaces is', nsect
     allocate(sector_basis(ncfgs,nsect))

     ! build good quantum numbers for each orbital
     !--------------------------------------------------------------------
     ! calculate orb_sz
     orb_sz = 0
     call atomic_make_gsz(orb_sz)
     write(mystd,'(4X,a)') 'compute [Sz] for orbitals'
     !
     ! calculate orb_jz
     orb_jz = 0
     if ( nband == 3 .or. nband == 5 .or. nband == 7 ) then
         ! jz only valid for nband == 3, 5, 7
         call atomic_make_gjz(orb_jz)
     endif ! back if ( nband == 3 .or. nband == 5 .or. nband == 7 ) block
     write(mystd,'(4X,a)') 'compute [Jz] for orbitals'
     !
     ! calculate orb_ps
     orb_ps = 0
     call atomic_make_gps(orb_ps)
     write(mystd,'(4X,a)') 'compute [PS] for orbitals'

! build good quantum numbers for each Fock state
!-------------------------------------------------------------------------
     counter = 0
     fock_ntot = 0
     fock_sz = 0
     fock_jz = 0
     fock_ps = 0
! loop over all number of total electrons
     do i=0,norbs
! loop over each state
         do j=1,dim_sub_n(i)
! here counter denotes the index of Fock state
             counter = counter + 1
! build N
             fock_ntot(counter) = i
! build Sz
             my_sz = 0
             do k=1,norbs
                 my_sz = my_sz + orb_sz(k) * bin_basis(k,counter)
             enddo ! over k={1,norbs} loop
             fock_sz(counter) = my_sz
! build Jz
              my_jz = 0
              do k=1,norbs
                  my_jz = my_jz + orb_jz(k) * bin_basis(k,counter)
              enddo ! over k={1,norbs} loop
              fock_jz(counter) = my_jz
! build PS number
             do k=1,nband
                 fock_ps(counter) = &
                 fock_ps(counter) + (2**k) * &
                     (bin_basis(2*k-1,counter) - bin_basis(2*k,counter))**2
             enddo ! over k={1,nband} loop
         enddo ! over j={1,dim_sub_n(i)} loop
     enddo ! over i={0,norbs} loop

! loop over all the Fock states to determine sectors
!-------------------------------------------------------------------------
     nsect = 0
     ndims = 0
! initialize some variables
     sect_ntot = 0
     sect_sz = 0
     sect_jz = 0
     sect_ps = 0
     sector_basis = 0
     do i=1,ncfgs
         my_ntot = fock_ntot(i)

! truncate the occupancy according to nmini and nmaxi
         if ( my_ntot < nmini  .or. my_ntot > nmaxi ) CYCLE

         if ( ictqmc == 3 .or. ictqmc == 4 ) then
             my_sz = fock_sz(i)
         endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
         if ( ictqmc == 4 ) then
             my_ps = fock_ps(i)
         endif ! back if ( ictqmc == 4 ) block
         if ( ictqmc == 5 ) then
             my_jz = fock_jz(i)
         endif ! back if ( ictqmc == 5 ) block

! determine the first sector
         if ( nsect == 0 ) then
             sect_ntot(1) = my_ntot
             if ( ictqmc == 3 .or. ictqmc == 4 ) then
                 sect_sz(1) = my_sz
             endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
             if (ictqmc == 4) then
                 sect_ps(1) = my_ps
             endif ! back if ( ictqmc == 4 ) block
             if ( ictqmc == 5 ) then
                 sect_jz(1) = my_jz
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
                         if ( sect_ntot(j) == my_ntot ) then
                             which_sect = j; EXIT
                         endif ! back if ( sect_ntot(j) == my_ntot ) block

                     case (3)
                         if ( sect_ntot(j) == my_ntot ) then
                             if ( sect_sz(j) == my_sz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_sz(j) == my_sz ) block
                         endif ! back if ( sect_ntot(j) == my_ntot ) block

                     case (4)
                         if ( sect_ntot(j) == my_ntot ) then
                             if ( sect_sz(j) == my_sz ) then
                                 if ( sect_ps(j) == my_ps) then
                                     which_sect = j; EXIT
                                 endif ! back if ( sect_ps(j) == my_ps) block
                             endif ! back if ( sect_sz(j) == my_sz ) block
                         endif ! back if ( sect_ntot(j) == my_ntot ) block

                     case (5)
                         if ( sect_ntot(j) == my_ntot ) then
                             if ( sect_jz(j) == my_jz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_jz(j) == my_jz ) block
                         endif ! back if ( sect_ntot(j) == my_ntot ) block

                 end select
             enddo ! over j={1,nsect} loop

! we can not assign the current state into any existing sectors, so we
! have to define a new sector
             if ( which_sect == -1 ) then
                 nsect = nsect + 1
                 sect_ntot(nsect) = my_ntot
                 if ( ictqmc == 3 .or. ictqmc == 4 ) then
                     sect_sz(nsect) = my_sz
                 endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
                 if ( ictqmc == 4 ) then
                     sect_ps(nsect) = my_ps
                 endif ! back if ( ictqmc == 4 ) block
                 if ( ictqmc == 5 ) then
                     sect_jz(nsect) = my_jz
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
         sectors(i)%nele = sect_ntot(i)
         sectors(i)%sz   = sect_sz(i)
         sectors(i)%jz   = sect_jz(i)
         sectors(i)%ps   = sect_ps(i)
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
                                 my_ntot = sect_ntot(i) + 1
                             else
                                 my_ntot = sect_ntot(i) - 1
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_ntot(l) == my_ntot ) then
                                     which_sect = l; EXIT
                                 endif ! back if ( sect_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (3)
                             if ( k == 1 ) then
                                 my_ntot = sect_ntot(i) + 1
                                 my_sz = sect_sz(i) + orb_sz(j)
                             else
                                 my_ntot = sect_ntot(i) - 1
                                 my_sz = sect_sz(i) - orb_sz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_ntot(l) == my_ntot ) then
                                     if ( sect_sz(l) == my_sz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_sz(l) == my_sz ) block
                                 endif ! back if ( sect_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (4)
                             if ( k == 1 ) then
                                 my_ntot = sect_ntot(i) + 1
                                 my_sz = sect_sz(i) + orb_sz(j)
                                 code(j) = 1
                             else
                                 my_ntot = sect_ntot(i) - 1
                                 my_sz   = sect_sz(i) - orb_sz(j)
                                 code(j) = 0
                             endif ! back if ( k == 1 ) block
! calculate new PS number
                             my_ps = 0
                             do l=1,nband
                                 my_ps = my_ps + (2**l) * ( code(2*l-1) - code(2*l) )**2
                             enddo ! over l={1,nband} loop
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_ntot(l) == my_ntot ) then
                                     if ( sect_sz(l) == my_sz ) then
                                         if ( sect_ps(l) == my_ps) then
                                             which_sect = l; EXIT
                                         endif ! back if ( sect_ps(l) == my_ps) block
                                     endif ! back if ( sect_sz(l) == my_sz ) block
                                 endif ! back if ( sect_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (5)
                             if ( k == 1 ) then
                                 my_ntot = sect_ntot(i) + 1
                                 my_jz = sect_jz(i) + orb_jz(j)
                             else
                                 my_ntot = sect_ntot(i) - 1
                                 my_jz = sect_jz(i) - orb_jz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_ntot(l) == my_ntot ) then
                                     if ( sect_jz(l) == my_jz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_jz(l) == my_jz ) block
                                 endif ! back if ( sect_ntot(l) == my_ntot ) block
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
     call atomic_dump_sector(sect_ntot, sect_sz, sect_ps, sect_jz)

! deallocate memory
     deallocate(sector_basis)

     return
  end subroutine atomic_make_sectors

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
     use constants, only : dp
     use constants, only : one
     use constants, only : czero
     use constants, only : epst
     use constants, only : mystd

     use control, only : norbs

     use m_fock, only : bin_basis
     use m_fock, only : dec_basis
     use m_fock, only : ind_basis

     use m_spmat, only : emat
     use m_spmat, only : umat

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! loop index for subspace (sector)
     integer :: isec

     ! loop index for Fock state
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

!! [body

     ! loop over all subspaces
     do isec=1,nsectors

         write(mystd,'(4X,a,i4)') 'subspace: ', isec

         ! start to make atomic Hamiltonian
         ! we should initialize hmat at first
         sectors(isec)%hmat = czero

         ! compute two fermion operators term (onsite impurity energy)
         ! it is f^{\dagger}_{\alpha} f_{\beta}
         !----------------------------------------------------------------
         write(mystd,'(4X,a)') 'compute two fermion operators term'
         !
         do jbas=1,sectors(isec)%ndim
             alploop: do alpha=1,norbs
             betloop: do betta=1,norbs

             ! retrieve the Fock state |jbas>
             isgn = 0
             knew = dec_basis(sectors(isec)%basis(jbas))
             code = bin_basis(:,sectors(isec)%basis(jbas))

             ! impurity level is too small
             if ( abs(emat(alpha,betta)) < epst ) CYCLE

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
                             EXIT
                         endif ! back if block
                     enddo ! over ibas={1,sectors(isec)%ndim} loop
                     !
                     ! write the Fock states and the operators
                     write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha = ', alpha, ')'
                     write(mystd,'(2X,a,i2,a)') 'f(beta = ', betta, ')'
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

             ! retrieve the Fock state |jbas>
             isgn = 0
             knew = dec_basis(sectors(isec)%basis(jbas))
             code = bin_basis(:,sectors(isec)%basis(jbas))

             ! applying Pauli principle
             if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE

             ! U matrix element is too small
             if ( abs(umat(alpha,betta,delta,gamma)) < epst ) CYCLE

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
                     !
                     ! write the Fock states and the operators
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
         enddo ! over jbas={1,sectors(isect)%ndim} loop

     enddo ! over isec={1,nsectors} loop

!! body]

     return
  end subroutine atomic_make_shmat

!!
!! @sub atomic_diag_shmat
!!
!! diagonalize the atomic Hamiltonian subspace by subspace
!!
  subroutine atomic_diag_shmat()
     use constants, only : dp
     use constants, only : mystd

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! dummy array
     real(dp), allocatable :: hmat(:,:)

!! [body

     do i=1,nsectors

         ! we will not destroy the raw Hamiltonian data,
         ! so we usually make a copy of it
         allocate(hmat(sectors(i)%ndim,sectors(i)%ndim))
         hmat = real( sectors(i)%hmat )

         ! diagonalize it, eval and evec will be updated
         call s_eig_sy( sectors(i)%ndim, &
                        sectors(i)%ndim, &
                        hmat,            &
                        sectors(i)%eval, &
                        sectors(i)%evec )

         ! deallocate memory
         deallocate(hmat)

         write(mystd,'(4X,a,i4,2X,a)') 'subspace: ', i, 'done'

     enddo ! over i={1,nsectors} loop

!! body]

     return
  end subroutine atomic_diag_shmat

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
