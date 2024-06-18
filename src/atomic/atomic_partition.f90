!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : automatic_partition
!!!           sector_create
!!!           sector_refine
!!!           sector_filter
!!!           sector_copyto
!!!           sector_locate
!!!           sector_lookup
!!!           map_create
!!!           map_remove
!!!           map_verify
!!!           get_sector_ntot
!!!           get_sector_sz
!!!           get_sector_jz
!!!           get_sector_ap
!!! source  : atomic_partition.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 06/11/2024 by li huang (created)
!!!           06/18/2024 by li huang (last modified)
!!! purpose : implement the automatic partition algorithm to divide the
!!!           atomic Hamiltonian into many blocks.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub automatic_partition
!!
!! this subroutine implements the automatic partition algorithm.
!!
!! we wish to find a permutation of the basis vectors such that (1) the
!! local Hamiltonian is block-diagonal, and (2) all f^{+} and f operators
!! are block matrices with at most one non-zero block in each row and
!! column. such a permutation would group the basis states belonging to
!! the same subspace together.
!!
!! the automatic partitioning algorithm employs no additional a priori
!! information about the symmetry of Hamiltonian H_{loc}. the only input
!! data is the full set of basis Fock states and the Hamiltonian itself.
!! the algorithm consists of two sequential phases. in the first phase,
!! the finest possible partition, which satisfies condition (1) alone is
!! constructed. and then in the second phase, this partition is adjusted
!! to additionally satisfy condition (2) further.
!!
!! see: Computer Physics Communications 200, 274–284 (2016)
!!
  subroutine automatic_partition()
     use constants, only : zero
     use constants, only : mystd

     use control, only : ictqmc
     use control, only : norbs, ncfgs
     use control, only : nmini, nmaxi

     use m_fock, only : bin_basis, dec_basis, ind_basis

     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_sector
     use m_sector, only : cat_alloc_sectors

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k
     integer :: l

     ! any good quantum number
     integer :: q

     ! index of Fock state
     integer :: sib

     ! number of subspaces (sectors)
     integer :: nsect
     integer :: nsect_

     ! index of selected subspace
     integer :: which_sect

     ! can point to next subspace (sector)
     logical :: can

     ! N, Sz, Jz, and AP for all subspaces
     integer, allocatable :: sect_ntot(:)
     integer, allocatable :: sect_sz(:)
     integer, allocatable :: sect_jz(:)
     integer, allocatable :: sect_ap(:)

     ! dimension for subspaces
     integer, allocatable :: ndims(:)
     integer, allocatable :: ndims_(:)

     ! global indices of Fock states of subspaces
     !
     ! the first index: local index of Fock state in the given subspace
     ! the second index: index of subspace
     !
     ! for example, sector_basis(j,i) gives the global index of Fock
     ! state for the j-th basis function in the i-th subspace
     integer, allocatable :: sector_basis(:,:)
     integer, allocatable :: sector_basis_(:,:)

!! [body

     ! initialization
     allocate(ndims_(ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))

     ! phase 1
     call sector_create(ncfgs, ndims_, sector_basis_)

     nsect_ = count(ndims_ > 0)

     ! filter sectors
     allocate(ndims(nsect_))
     allocate(sector_basis(ncfgs,nsect_))
     !
     !j = 0
     !do i=1,ncfgs
     !    if ( ndims_(i) > 0 ) then
     !        j = j + 1
     !        ndims(j) = ndims_(i)
     !        sector_basis(:,j) = sector_basis_(:,i)
     !    endif
     !enddo
     !
     !call s_assert(j == nsect_)
     call sector_filter(ncfgs, ndims_, sector_basis_, &
                        & nsect_, ndims, sector_basis)
     deallocate(ndims_)
     deallocate(sector_basis_)

     ! phase 2
     call sector_refine(nsect_, ndims, sector_basis)

     nsect = count(ndims > 0)
     print *, 'number of sectors: ', nsect

     ! setup subspaces
     allocate(sect_ntot(nsect_))
     allocate(sect_sz(nsect_))
     allocate(sect_jz(nsect_))
     allocate(sect_ap(nsect_))
     sect_ntot = 0
     sect_sz = 0
     sect_jz = 0
     sect_ap = 0

     k = 0
     do i=1,nsect_
         if ( ndims(i) > 0 ) then
             k = k + 1
             call get_sector_ntot(q, ndims(i), sector_basis(:,i))
             sect_ntot(i) = q
             !
             call get_sector_sz(q, ndims(i), sector_basis(:,i))
             sect_sz(i) = q
             !
             call get_sector_jz(q, ndims(i), sector_basis(:,i))
             sect_jz(i) = q
             !
             call get_sector_ap(q, i, nsect_, sect_ntot, sect_sz, sect_jz)
             sect_ap(i) = q

             write(mystd,'(a,i6)') 'subspace -> ', k
             write(mystd,'(a,i6)') 'size :', ndims(i)
             write(mystd,'(a)') 'basis :'
             do j=1,ndims(i)
                 write(mystd,'(i,2X,14i1)') j, bin_basis(:,sector_basis(j,i))
             enddo
             write(mystd, '(a, i3)') 'N :', sect_ntot(i)
             write(mystd, '(a, i3)') 'Sz:', sect_sz(i)
             write(mystd, '(a, i3)') 'Jz:', sect_jz(i)
             write(mystd, '(a, i3)') 'AP:', sect_ap(i)
             write(mystd, *)
         endif
     enddo
     call s_assert(k == nsect)

     write(mystd,'(4X,a)') 'allocate memory for subspaces'
     !
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect ! do not forget to setup nsectors
     !
     call cat_alloc_sectors()
     !
     ! next we will build every subspace one by one
     sib = 1
     k = 0
     do i=1,nsect_
         if ( ndims(i) > 0 ) then
             k = k + 1
             sectors(k)%istart = sib
             sectors(k)%ndim = ndims(i)
             sectors(k)%nops = norbs
             !
             sectors(k)%nele = sect_ntot(i)
             sectors(k)%sz   = sect_sz(i)
             sectors(k)%jz   = sect_jz(i)
             sectors(k)%ps   = sect_ap(i)
             !
             sib = sib + ndims(i)

             ! allocate memory for the subspace
             call cat_alloc_sector( sectors(k) )

             ! setup basis for the subspace
             do j=1,ndims(i)
                 sectors(k)%basis(j) = sector_basis(j,i)
             enddo ! over j={1,ndims(i)} loop

             write(mystd,'(4X,a,i4)', advance = 'no') 'subspace:', k
             write(mystd,'(2X,a,i4)', advance = 'no') 'size:', ndims(i)
             write(mystd,'(2X,a,i4)') 'start:', sectors(k)%istart
         endif
     enddo ! over i={1,nsect_} loop

     ! make index for next subspace
     !--------------------------------------------------------------------
     write(mystd,'(4X,a)') 'simulate fermion operator acts on subspaces'
     do i=1,nsectors  ! loop over all the subspaces
         do j=1,norbs ! loop over all the orbtials
             do k=0,1 ! loop over creation and annihilation fermion operators

                 which_sect = -1

                 ! we should check each Fock state in this subspace
                 can = .false.
                 do l=1,sectors(i)%ndim
                     sib = sectors(i)%basis(l)

                     ! test creation fermion operator
                     if ( k == 1 .and. bin_basis(j,sib) == 0 ) then
                         can = .true.
                         EXIT
                     !
                     ! test annihilation fermion operator
                     else if ( k == 0 .and. bin_basis(j, sib) == 1 ) then
                         can = .true.
                         EXIT
                     !
                     endif ! back if ( k == 1 .and. bin_basis(j,sib) == 0 ) block
                 enddo ! over l={1,sectors(i)%ndim} loop

                 ! if can == .true., it means that the fermion operator
                 ! can act on the given subspace. next, we would like to
                 ! figure out the resulting subspace.
                 if ( can .eqv. .true. ) then

                     if ( k == 1 ) then
                         q = dec_basis(sib) + 2**(j-1)
                         call sector_lookup(which_sect, ind_basis(q))
                     endif

                     if ( k == 0 ) then
                         q = dec_basis(sib) - 2**(j-1)
                         call sector_lookup(which_sect, ind_basis(q))
                     endif

                 endif  ! back if ( can == .true. ) block

                 ! setup the next array
                 sectors(i)%next(j,k) = which_sect
                 call map_verify(i, j, k, which_sect)

                 if (k == 1) then
                     write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha =', j, ')'
                     write(mystd,'(2X,a,i4)', advance = 'no') '|subspace>_i:', i
                     write(mystd,'(2X,a,i4)') '|subspace>_f:', which_sect
                 else
                     write(mystd,'(4X,a,i2,a)', advance = 'no') 'f  (alpha =', j, ')'
                     write(mystd,'(2X,a,i4)', advance = 'no') '|subspace>_i:', i
                     write(mystd,'(2X,a,i4)') '|subspace>_f:', which_sect
                 endif ! back if (k == 1) block
                 !print *

             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,nsectors} loop

     ! calculate the maximum and average dimensions of subspaces
     !--------------------------------------------------------------------
     max_dim_sect = maxval(ndims)
     ave_dim_sect = sum(ndims) / real(nsectors)
     !
     write(mystd,'(4X,a,i4)') 'maximum dimension of subspaces:', max_dim_sect
     write(mystd,'(4X,a,f6.2)') 'averaged dimension of subspaces:', ave_dim_sect

     ! dump subspace information for reference
     !--------------------------------------------------------------------
     !call atomic_dump_sector(sect_ntot, sect_sz, sect_ap, sect_jz)

     ! deallocate memory
     deallocate(sector_basis)

!! body]

     return
  end subroutine automatic_partition

!!
!! @sub sector_create
!!
!! phase 1 of the automatic partition algorithm. we just permute the Fock
!! states such that the resulting local Hamiltonian is block-diagonal.
!!
!! see: Computer Physics Communications 200, 274–284 (2016)
!!
  subroutine sector_create(nsect, ndims, sector_basis)
     use constants, only : zero

     use control, only : ncfgs

     use m_fock, only : hmat

     implicit none

!! external arguments
     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(inout) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

!! local variables
     ! loop index for Fock states
     integer :: i
     integer :: j

     ! index for subspaces
     integer :: ia
     integer :: ib

!! [body

     ! to start, one creates a data structure, which stores information
     ! about the way N basis states (Fock states) are partitioned into
     ! a number of subsets. initially each basis state resides alone in
     ! its own subset.
     ndims = 1
     sector_basis = 0
     do i=1,nsect
         sector_basis(1,i) = i
     enddo

     ! in the main loop of the algorithm, the Hamiltonian is sequentially
     ! applied to each basis state (initial state). this typically gives
     ! a linear combination of the basis states with only a few non-zero
     ! coefficients, since H_{loc} is usually sparse.
     !
     ! the algorithm iterates in an inner loop over all basis vectors
     ! with non-zero coefficients (final states).
     !
     ! if the initial and final states reside in different subsets, these
     ! subsets are merged together. once the main loop is over, the
     ! partition of the basis is done.
     do i=1,ncfgs     ! for initial states
         do j=1,ncfgs ! for final states
             if ( abs(hmat(i,j)) > zero ) then
                 ! get related subspaces for the initial and final states
                 call sector_locate(ia, i, nsect, ndims, sector_basis)
                 call sector_locate(ib, j, nsect, ndims, sector_basis)

                 ! merge the two subspaces if they are not the same
                 if ( ia /= ib ) then
                     call sector_copyto(ia, ib, nsect, ndims, sector_basis)
                 endif ! back if ( ia /= ib ) block
             endif
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,ncfgs} loop

!! body]

     return
  end subroutine sector_create

!!
!! @sub sector_refine
!!
!! phase 2 of the automatic partition algorithm. in this stage, we should
!! make sure all f^+ and f operators are block matrices with at most one
!! non-zero block in each row and column.
!!
!! see: Computer Physics Communications 200, 274–284 (2016)
!!
  subroutine sector_refine(nsect, ndims, sector_basis)
     use control, only : norbs, ncfgs

     implicit none

!! external arguments
     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(inout) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

!! local variables
     ! loop index for orbitals
     integer :: i

     ! loop index for subspace-to-subspace connections
     integer :: j

     ! index for subspaces
     integer :: HA, HL, HU

     ! upward subspace-to-subspace connections
     integer, allocatable :: Mup(:,:)

     ! downward subspace-to-subspace connections
     integer, allocatable :: Mdn(:,:)

!! [body

     ! allocate memory for subspace-to-subspace connections
     allocate(Mup(ncfgs/2,2))
     allocate(Mdn(ncfgs/2,2))

     do i=1,norbs
         ! build upward and downward subspace-to-subspace connections
         call map_create(i, nsect, ndims, sector_basis, Mup, Mdn)

         ! scan all the connections, try to remove them and permute
         ! the Fock states in different subspaces
         do j=1,ncfgs/2
             ! get subspace-to-subspace connection
             HA = Mup(j,1)
             HL = Mup(j,1)
             HU = Mup(j,2)
             !
             ! well, meet empty connection, just skip it
             if ( HL == 0 ) CYCLE
             if ( HU == 0 ) CYCLE
             !
             ! remove connections and group the Fock states belonging
             ! to the same subspace
             call map_remove(1, HA, HL, HU, &
                             & nsect, ndims, sector_basis, Mup, Mdn)
         enddo ! over j={1,ncfgs/2} loop
     enddo ! over i={1,norbs} loop

     ! deallocate memory
     deallocate(Mup)
     deallocate(Mdn)

!! body]

     return
  end subroutine sector_refine

!!
!! @sub sector_filter
!!
!! compress the original subspaces (the empty subspaces are discarded)
!!
  subroutine sector_filter(nsect_, ndims_, sector_basis_, nsect, ndims, sector_basis)
     use control, only : ncfgs

     implicit none

!! external arguments
     ! number of subspaces (old)
     integer, intent(in) :: nsect_

     ! dimension for subspaces (old)
     integer, intent(in) :: ndims_(nsect_)

     ! global indices of Fock states of subspaces (old)
     integer, intent(in) :: sector_basis_(ncfgs,nsect_)

     ! number of subspaces (new)
     integer, intent(inout) :: nsect

     ! dimension for subspaces (new)
     integer, intent(inout) :: ndims(nsect)

     ! global indices of Fock states of subspaces (new)
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

!! local variables
     ! loop index
     integer :: i

     ! counter
     integer :: c

!! [body

     c = 0
     !
     do i=1,nsect_
         if ( ndims_(i) > 0 ) then
             c = c + 1
             ndims(c) = ndims_(i)
             sector_basis(:,c) = sector_basis_(:,i)
         endif
     enddo
     !
     call s_assert(c == nsect)

!! body]

     return
  end subroutine sector_filter

!!
!! @sub sector_copyto
!!
!! merge two subspaces (subspace src is at first merged with subspace
!! dst, and then src is cleaned.)
!!
  subroutine sector_copyto(dst, src, nsect, ndims, sector_basis)
     use control, only : ncfgs

     implicit none

!! external arguments
     ! index for the destination subspace
     integer, intent(in) :: dst

     ! index for the source subspace
     integer, intent(in) :: src

     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(inout) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

!! local variables
     ! loop index
     integer :: m

!! [body

     ! check the two subspaces (dst and src). they should not be empty.
     call s_assert(ndims(dst) >= 1)
     call s_assert(ndims(src) >= 1)

     ! copy Fock states from src to dst
     do m=1,ndims(src)
           sector_basis(ndims(dst) + m, dst) = sector_basis(m, src)
     enddo ! over m={1,ndims(src)} loop
     ndims(dst) = ndims(dst) + ndims(src)
     !
     ! clear Fock states in src
     sector_basis(:,src) = 0
     ndims(src) = 0

!! body]

     return
  end subroutine sector_copyto

!!
!! @sub sector_locate
!!
!! figure out the subspace that contains the given Fock state
!!
  subroutine sector_locate(sind, find, nsect, ndims, sector_basis)
     use control, only : ncfgs

     implicit none

!! external arguments
     ! index of subspace that contains the given Fock state
     integer, intent(out) :: sind

     ! index of Fock state
     integer, intent(in) :: find

     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(in) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(in) :: sector_basis(ncfgs,nsect)

!! local variables
     ! loop index
     integer :: m
     integer :: n

!! [body

     sind = 0
     !
     ! loop over subspaces
     SECTOR_LOOP: do m=1,nsect
         ! loop over Fock states in m-th subspace
         do n=1,ndims(m)
             ! find the given Fock state
             if ( sector_basis(n,m) == find ) then
                 ! record index of the current subspace
                 sind = m
                 EXIT SECTOR_LOOP
             endif
         enddo
     enddo SECTOR_LOOP
     !
     ! we can always find the required subspace
     call s_assert(sind /= 0)

!! body]

     return
  end subroutine sector_locate

!!
!! @sub sector_lookup
!!
!! figure out the subspace that contains the given Fock state
!!
  subroutine sector_lookup(sind, find)
     use m_fock, only : bin_basis

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! external arguments
     ! index of subspace that contains the given Fock state
     integer, intent(out) :: sind

     ! index of Fock state
     integer, intent(in) :: find

!! local variables
     ! loop index
     integer :: m
     integer :: n

     ! occupancy of Fock state
     integer :: ntot

!! [body

     ! get occupancy for the given Fock state
     ntot = sum(bin_basis(:,find))

     sind = 0
     !
     ! loop over subspaces
     SECTOR_LOOP: do m=1,nsectors
         !
         ! check occupancy of subspace
         if ( sectors(m)%nele /= ntot ) CYCLE
         !
         ! loop over Fock states in m-th subspace
         do n=1,sectors(m)%ndim
             ! find the given Fock state
             if ( sectors(m)%basis(n) == find ) then
                 ! record index of the current subspace
                 sind = m
                 EXIT SECTOR_LOOP
             endif
         enddo
         !
     enddo SECTOR_LOOP
     !
     ! we can always find the required subspace
     call s_assert(sind /= 0)

!! body]

     return
  end subroutine sector_lookup

!!
!! @sub map_create
!!
!! try to create upward and downward subspace-to-subspace connections.
!!
!! see: Computer Physics Communications 200, 274–284 (2016)
!!
  subroutine map_create(iorb, nsect, ndims, sector_basis, Mup, Mdn)
     use control, only : ncfgs

     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

!! external arguments
     integer, intent(in) :: iorb

     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(in) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(in) :: sector_basis(ncfgs,nsect)

     ! upward subspace-to-subspace connections
     integer, intent(out) :: Mup(ncfgs/2,2)

     ! downward subspace-to-subspace connections
     integer, intent(out) :: Mdn(ncfgs/2,2)

!! local variables
     ! loop index for Fock state
     integer :: i

     ! decimal form of the resulting Fock state
     integer :: k

     ! index for subspace
     ! the connection is from ia to ib
     integer :: ia
     integer :: ib

     ! counter for upward and downward subspace-to-subspace connection
     integer :: cup
     integer :: cdn

!! [body

     ! init the counters
     cup = 0
     cdn = 0

     ! init the connections
     Mup = 0
     Mdn = 0

     do i=1,ncfgs
         ! find the subspace ia that contains i-th Fock state
         call sector_locate(ia, i, nsect, ndims, sector_basis)

         ! build connection for f^{+} operator
         if ( bin_basis(iorb,i) == 0 ) then
             ! get new Fock state k
             k = dec_basis(i) + 2**(iorb-1)

             ! find the subspace ib that contains the new Fock state k
             call sector_locate(ib, ind_basis(k), nsect, ndims, sector_basis)

             ! store the connection
             cup = cup + 1
             Mup(cup,1) = ia
             Mup(cup,2) = ib
         endif ! back if ( bin_basis(iorb,i) == 0 ) block

         ! build connection for f operator
         if ( bin_basis(iorb,i) == 1 ) then
             ! get new Fock state k
             k = dec_basis(i) - 2**(iorb-1)

             ! find the subspace ib that contains the new Fock state k
             call sector_locate(ib, ind_basis(k), nsect, ndims, sector_basis)

             ! store the connection
             cdn = cdn + 1
             Mdn(cdn,1) = ia
             Mdn(cdn,2) = ib
         endif ! back if ( bin_basis(iorb,i) == 1 ) block
     enddo ! over i={1,ncfgs} loop

     ! verify the number of subspace-to-subspace connections
     call s_assert(cup == ncfgs / 2)
     call s_assert(cdn == ncfgs / 2)

!! body]

     return
  end subroutine map_create

!!
!! @sub map_remove
!!
!! remove the upward and downward subspace-to-subspace connections by
!! using the zigzag algorithm.
!!
!! see: Computer Physics Communications 200, 274–284 (2016)
!!
  recursive &
  subroutine map_remove(dir, HA, HL, HU, nsect, ndims, sector_basis, Mup, Mdn)
     use control, only : ncfgs

     implicit none

!! external arguments
     ! direction of the zigzag algorithm
     ! dir = 1 means upward, and dir = 2 means downward
     integer, intent(in) :: dir

     ! indices for subspaces
     !
     ! here L and U stand for 'lower' and 'upper' respectively
     ! f^{+} |H_L> = |H_U> and f |H_U> = |H_L>
     integer, intent(in) :: HA
     integer, intent(in) :: HL
     integer, intent(in) :: HU

     ! number of subspaces
     integer, intent(in) :: nsect

     ! dimension for subspaces
     integer, intent(inout) :: ndims(nsect)

     ! global indices of Fock states in subspaces
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

     ! upward subspace-to-subspace connections
     integer, intent(inout) :: Mup(ncfgs/2,2)

     ! downward subspace-to-subspace connections
     integer, intent(inout) :: Mdn(ncfgs/2,2)

!! local variables
     ! loop index for connections
     integer :: i

     ! index for subspace
     integer :: HB

!! [body

     ! for upward connection
     !
     ! search all upward connections starting at H_A.
     ! each such connection H_A -> H_B is removed from Mup.
     ! H_B is merged with H_U.
     ! call map_remove(2, ...)
     if ( dir == 1 ) then
         do i=1,ncfgs/2
             if ( Mup(i,1) /= HA ) CYCLE
             HB = Mup(i,2)
             call s_assert(HB /= 0)
             !
             Mup(i,:) = 0
             !
             if ( HB /= HU .and. ndims(HB) > 0 ) then
                 call sector_copyto(HU, HB, nsect, ndims, sector_basis)
             endif
             !
             call map_remove(2, HB, HL, HU, &
                                 & nsect, ndims, sector_basis, Mup, Mdn)
         enddo ! over i={1,ncfgs/2} loop

     ! for downward connection
     !
     ! search all downward connections starting at H_A.
     ! each such connection H_A -> H_B is removed from Mdn.
     ! H_B is merged with H_L.
     ! call map_remove(1, ...)
     else
         do i=1,ncfgs/2
             if ( Mdn(i,1) /= HA ) CYCLE
             HB = Mdn(i,2)
             call s_assert(HB /= 0)
             !
             Mdn(i,:) = 0
             !
             if ( HB /= HL .and. ndims(HB) > 0 ) then
                 call sector_copyto(HL, HB, nsect, ndims, sector_basis)
             endif
             !
             call map_remove(1, HB, HL, HU, &
                                 & nsect, ndims, sector_basis, Mup, Mdn)
         enddo ! over i={1,ncfgs/2} loop

     endif ! back if ( dir == 1 ) block

!! body]

     return
  end subroutine map_remove

!!
!! @sub map_verify
!!
!! an operator f^{+}_{j} or f_{j} is acted on a given subspace i. this
!! subroutine will try to figure out the resulting subspace, and judge
!! whether it coincides with which_sect. actually, this subroutine is
!! used to check the mapping between two subspaces.
!!
  subroutine map_verify(i, j, k, which_sect)
     use m_fock, only : bin_basis, dec_basis, ind_basis

     use m_sector, only : sectors

     implicit none

!! external arguments
     ! index of initial subspace
     integer, intent(in) :: i

     ! index of orbital
     integer, intent(in) :: j

     ! type of operator, k = 1 for creator and k = 0 for destroyer
     integer, intent(in) :: k

     ! index of final subspace
     integer, intent(in) :: which_sect

!! local variables
     ! loop index for Fock state
     integer :: n

     ! index for resulting Fock state
     integer :: m

     ! index for the n-th Fock state in the i-th subspace
     integer :: sib

     ! decimal form of the resulting Fock state
     integer :: knew

!! [body

     ! for f^{+}_j operator
     if ( k == 1 ) then
         do n=1,sectors(i)%ndim
             ! go through each Fock state in the i-th subspace
             sib = sectors(i)%basis(n)
             !
             ! evaluate the resulting Fock state and then make sure it is
             ! in a subspace that is specified by which_sect
             if ( bin_basis(j,sib) == 0 ) then
                 knew = dec_basis(sib) + 2**(j-1)
                 m = ind_basis(knew)
                 call s_assert( count(sectors(which_sect)%basis == m) == 1 )
             endif
         enddo ! over n={1,sectors(i)%ndim} loop
     endif ! back if ( k == 1 ) block

     ! for f_j operator
     if ( k == 0 ) then
         do n=1,sectors(i)%ndim
             ! go through each Fock state in the i-th subspace
             sib = sectors(i)%basis(n)
             !
             ! evaluate the resulting Fock state and then make sure it is
             ! in a subspace that is specified by which_sect
             if ( bin_basis(j,sib) == 1 ) then
                 knew = dec_basis(sib) - 2**(j-1)
                 m = ind_basis(knew)
                 call s_assert( count(sectors(which_sect)%basis == m) == 1 )
             endif
         enddo ! over n={1,sectors(i)%ndim} loop
     endif ! back if ( k == 1 ) block

!! body]

     return
  end subroutine map_verify

!!
!! @sub get_sector_ntot
!!
!! return good quantum number N for the given subspace (sector)
!!
  subroutine get_sector_ntot(N, nsize, sector_basis)
     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis

     implicit none

!! external arguments
     ! good quantum number N for the given subspace
     integer, intent(out) :: N

     ! capacity for the given subspace
     integer, intent(in) :: nsize

     ! basis for the given subspace
     integer, intent(in) :: sector_basis(ncfgs)

!! local variables
     ! loop index
     integer :: i

     ! good quantum number N
     integer :: N_

     ! Fock state in the given subspace
     integer :: code(norbs)

!! [body

     N = 0
     !
     do i=1,nsize
         ! visit each Fock state in the subspace
         code = bin_basis(:,sector_basis(i))

         ! get N for the current Fock state
         N_ = sum(code)

         ! record N for the first Fock state
         if ( i == 1 ) then
             N = N_
         ! all Fock states in the subspace should share the same N
         else
             if ( N /= N_ ) then
                 call s_print_error('get_sector_ntot','wrong good &
                     & quantum number N for this subspace!')
             endif ! back if ( N /= N_ ) block
         endif ! back if ( i == 1 ) block
     enddo ! over i={1,nsize} loop

!! body]

     return
  end subroutine get_sector_ntot

!!
!! @sub get_sector_sz
!!
!! return good quantum number Sz for the given subspace (sector)
!!
  subroutine get_sector_sz(Sz, nsize, sector_basis)
     use control, only : isoc
     use control, only : nband, norbs, ncfgs

     use m_fock, only : bin_basis

     implicit none

!! external arguments
     ! good quantum number Sz for the given subspace
     integer, intent(out) :: Sz

     ! capacity for the given subspace
     integer, intent(in) :: nsize

     ! basis for the given subspace
     integer, intent(in) :: sector_basis(ncfgs)

!! local variables
     ! loop index
     integer :: i

     ! good quantum number Sz
     integer :: Sz_

     ! Fock state in the given subspace
     integer :: code(norbs)

!! [body

     Sz = 0
     !
     if ( isoc == 1 ) return
     !
     do i=1,nsize
         ! visit each Fock state in the subspace
         code = bin_basis(:,sector_basis(i))

         ! get Sz for the current Fock state
         Sz_ = sum(code(1:nband)) - sum(code(nband+1:norbs))

         ! record Sz for the first Fock state
         if ( i == 1 ) then
             Sz = Sz_
         ! all Fock states in the subspace should share the same Sz
         else
             if ( Sz /= Sz_ ) then
                 call s_print_error('get_sector_sz','wrong good &
                     & quantum number Sz for this subspace!')
             endif ! back if ( Sz /= Sz_ ) block
         endif ! back if ( i == 1 ) block
     enddo ! over i={1,nsize} loop

!! body]

     return
  end subroutine get_sector_sz

!!
!! @sub get_sector_jz
!!
!! return good quantum number Jz for the given subspace (sector)
!!
  subroutine get_sector_jz(Jz, nsize, sector_basis)
     use control, only : icf, isoc
     use control, only : norbs, ncfgs

     use m_fock, only : bin_basis

     implicit none

!! external arguments
     ! good quantum number Jz for the given subspace
     integer, intent(out) :: Jz

     ! capacity for the given subspace
     integer, intent(in) :: nsize

     ! basis for the given subspace
     integer, intent(in) :: sector_basis(ncfgs)

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! good quantum number Jz
     integer :: Jz_

     ! Fock state in the given subspace
     integer :: code(norbs)

     ! Jz quantum numbers for every orbitals
     integer :: good_jz(norbs)

!! [body

     Jz = 0
     !
     if ( isoc == 0 ) return
     if ( isoc == 1 .and. icf > 0 ) return
     !
     call atomic_make_gjz(good_jz)
     !
     do i=1,nsize
         ! visit each Fock state in the subspace
         code = bin_basis(:,sector_basis(i))

         ! get Jz for the current Fock state
         Jz_ = 0
         do k=1,norbs
             Jz_ = Jz_ + good_jz(k) * code(k)
         enddo ! over k={1,norbs} loop

         ! record Jz for the first Fock state
         if ( i == 1 ) then
             Jz = Jz_
         ! all Fock states in the subspace should share the same Jz
         else
             if ( Jz /= Jz_ ) then
                 call s_print_error('get_sector_jz','wrong good &
                     & quantum number Jz for this subspace!')
             endif ! back if ( Jz /= Jz_ ) block
         endif ! back if ( i == 1 ) block
     enddo ! over i={1,nsize} loop

!! body]

     return
  end subroutine get_sector_jz

!!
!! @sub get_sector_ap
!!
!! return good quantum number AP for the given subspace (sector)
!!
  subroutine get_sector_ap(AP, ind, nsect, sect_ntot, sect_sz, sect_jz)
     use control, only : icf, isoc

     implicit none

!! external arguments
     ! good quantum number AP for the given subspace
     integer, intent(out) :: AP

     ! index for the given subspace
     integer, intent(in) :: ind

     ! number of subspaces
     integer, intent(in) :: nsect

     ! GQN N from the 1st to ind-th subspaces
     integer, intent(in) :: sect_ntot(nsect)

     ! GQN Sz from the 1st to ind-th subspaces
     integer, intent(in) :: sect_sz(nsect)

     ! GQN Jz from the 1st to ind-th subspaces
     integer, intent(in) :: sect_jz(nsect)

!! local variables
     ! loop index
     integer :: j

     ! GQN N for the ind-th subspace
     integer :: N

     ! GQN Sz for the ind-th subspace
     integer :: Sz

     ! GQN Jz for the ind-th subspace
     integer :: Jz

!! [body

     AP = 0

     ! get N, Sz, and Jz for the ind-th subspace
     N  = sect_ntot(ind)
     Sz = sect_sz(ind)
     Jz = sect_jz(ind)

     ! spin-orbit coupling is disabled
     ! crystal field splitting is disabled or enabled
     ! good quantum numbers are N, Sz, and AP
     if ( isoc == 0 ) then
         AP = 1
         do j=1,ind-1
             if ( ( sect_ntot(j) == N ) .and. ( sect_sz(j) == Sz ) ) then
                 AP = AP + 1
             endif ! back if block
         enddo ! over j={1,ind-1} loop
     endif ! back if ( isoc == 0 ) block

     ! spin-orbit coupling is enabled
     ! crystal field splitting is disabled
     ! good quantum numbers are N, Jz, and AP
     if ( isoc == 1 .and. icf == 0 ) then
         AP = 1
         do j=1,ind-1
             if ( ( sect_ntot(j) == N ) .and. ( sect_jz(j) == Jz ) ) then
                 AP = AP + 1
             endif ! back if block
         enddo ! over j={1,ind-1} loop
     endif ! back if ( isoc == 1 .and. icf == 0 ) block

     ! spin-orbit coupling is enabled
     ! crystal field splitting is enabled
     ! good quantum numbers are N and AP
     if ( isoc == 1 .and. icf == 1 ) then
         AP = 1
         do j=1,ind-1
             if ( sect_ntot(j) == N ) then
                 AP = AP + 1
             endif ! back if ( sect_ntot(j) == N ) block
         enddo ! over j={1,ind-1} loop
     endif ! back if ( isoc == 1 .and. icf == 1 ) block

!! body]

     return
  end subroutine get_sector_ap
