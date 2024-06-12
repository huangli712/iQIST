!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : automatic_partition
!!!           sector_create
!!!           sector_refine
!!!           sector_locate
!!!           sector_lookup
!!!           sector_merge
!!!           sector_print
!!!           zigzag
!!!           try_next
!!!           get_sector_ntot
!!!           get_sector_sz
!!!           get_sector_jz
!!!           get_sector_ap
!!! source  : atomic_partition.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 06/11/2024 by li huang (created)
!!!           06/12/2024 by li huang (last modified)
!!! purpose : implement the automatic partition algorithm to divide the
!!!           atomic Hamiltonian into many blocks.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub automatic_partition
!!
!!
  subroutine automatic_partition()
     use constants, only : dp, mystd
     use constants, only : zero

     use control, only : ncfgs
     use control, only : norbs

     use m_fock, only : bin_basis, dec_basis, ind_basis

     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_sector
     use m_sector, only : cat_alloc_sectors

     implicit none

     integer :: i
     integer :: j
     integer :: k
     integer :: l
     integer :: m
     integer :: nsect, nsize, nsect_, nsize_
     integer :: N, Sz, Jz, Ap
     integer :: jold, jnew, isgn

     ! index of selected subspace
     integer :: which_sect

     ! index of Fock state
     integer :: ibasis

     ! can point to next subspace (sector)
     logical :: can

     integer, allocatable :: ndims(:)
     integer, allocatable :: ndims_(:)
     integer, allocatable :: sector_basis(:,:)
     integer, allocatable :: sector_basis_(:,:)

     integer, allocatable :: sect_ntot(:)
     integer, allocatable :: sect_sz(:)
     integer, allocatable :: sect_jz(:)
     integer, allocatable :: sect_ap(:)

     ! initialization
     allocate(ndims_(ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))
     call sector_create(ncfgs, ndims_, sector_basis_)

     nsect_ = count(ndims_ > 0)
     nsize_ = maxval(ndims_)
     print *, 'number of sectors: ', nsect_
     print *, 'maximum size of sectors: ', nsize_

     ! filter sectors
     allocate(ndims(nsect_))
     allocate(sector_basis(ncfgs,nsect_))
     !
     j = 0
     do i=1,ncfgs
         if ( ndims_(i) > 0 ) then
             j = j + 1
             ndims(j) = ndims_(i)
             sector_basis(:,j) = sector_basis_(:,i)
         endif
     enddo
     !
     call s_assert(j == nsect_)
     deallocate(ndims_)
     deallocate(sector_basis_)

     ! phase 2
     do i=1,norbs
         call sector_refine(i, nsect_, ndims, sector_basis)
     enddo

     nsect = count(ndims > 0)
     nsize = maxval(ndims)
     print *, 'number of sectors: ', nsect
     print *, 'maximum size of sectors: ', nsize

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
             call get_sector_ntot(N, ndims(i), sector_basis(:,i))
             call get_sector_sz(Sz, ndims(i), sector_basis(:,i))
             call get_sector_jz(Jz, ndims(i), sector_basis(:,i))
             call get_sector_ap(Ap, i, N, Sz, Jz, nsect_, sect_ntot, sect_sz, sect_jz)

             sect_ntot(i) = N
             sect_sz(i) = Sz
             sect_jz(i) = Jz
             sect_ap(i) = Ap
         endif
     enddo
     call s_assert(k == nsect)

     k = 0
     do i=1,nsect_
         if ( ndims(i) > 0 ) then
             k = k + 1
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

     write(mystd,'(4X,a)') 'allocate memory for subspaces'
     !
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect ! do not forget to setup nsectors
     !
     call cat_alloc_sectors()
     !
     ! next we will build every subspace one by one
     ibasis = 1
     k = 0
     do i=1,nsect_
         if ( ndims(i) > 0 ) then
             k = k + 1
             sectors(k)%istart = ibasis
             sectors(k)%ndim = ndims(i)
             sectors(k)%nops = norbs
             !
             sectors(k)%nele = sect_ntot(i)
             sectors(k)%sz   = sect_sz(i)
             sectors(k)%jz   = sect_jz(i)
             sectors(k)%ps   = sect_ap(i)
             !
             ibasis = ibasis + ndims(i)

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
                     ibasis = sectors(i)%basis(l)

                     ! test creation fermion operator
                     if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) then
                         can = .true.
                         EXIT
                     !
                     ! test annihilation fermion operator
                     else if ( k == 0 .and. bin_basis(j, ibasis) == 1 ) then
                         can = .true.
                         EXIT
                     !
                     endif ! back if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) block
                 enddo ! over l={1,sectors(i)%ndim} loop

                 ! if can == .true., it means that the fermion operator
                 ! can act on the given subspace. next, we would like to
                 ! figure out the resulting subspace.
                 if ( can .eqv. .true. ) then

                     if ( k == 1 ) then
                         jold = dec_basis(ibasis)
                         call atomic_make_cdagger(j, jold, jnew, isgn)
                         m = ind_basis(jnew)
                         call sector_lookup(which_sect, m)
                     endif

                     if ( k == 0 ) then
                         jold = dec_basis(ibasis)
                         call atomic_make_c(j, jold, jnew, isgn)
                         m = ind_basis(jnew)
                         call sector_lookup(which_sect, m)
                     endif

                 endif  ! back if ( can == .true. ) block

                 ! setup the next array
                 sectors(i)%next(j,k) = which_sect
                 call try_next(i, j, k, which_sect)

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

  subroutine sector_create(nsect, ndims_, sector_basis_)
     use constants, only : zero

     use control, only : ncfgs

     use m_fock, only : hmat

     implicit none

     integer, intent(in) :: nsect
     integer, intent(inout) :: ndims_(nsect)
     integer, intent(inout) :: sector_basis_(ncfgs,nsect)

     integer :: i, j
     integer :: ia, ib

     ndims_ = 1
     sector_basis_ = 0
     do i=1,nsect
         sector_basis_(1,i) = i
     enddo

     ! phase 1
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(hmat(i,j)) > zero ) then
                 call sector_locate(ia, i, nsect, ndims_, sector_basis_)
                 call sector_locate(ib, j, nsect, ndims_, sector_basis_)

                 if ( ia /= ib ) then
                     call sector_merge(ia, ib, nsect, ndims_, sector_basis_)
                 endif
             endif
         enddo
     enddo

     return
  end subroutine sector_create

  subroutine sector_refine(iorb, nsect, ndims, sector_basis)
     use control, only : ncfgs

     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

     integer, intent(in) :: iorb
     integer, intent(in) :: nsect
     integer, intent(inout) :: ndims(nsect)
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

     integer :: i
     integer :: j
     integer :: iup
     integer :: idn
     integer :: jnew, jold, isgn
     integer :: ia, ib
     integer :: HA, HL, HU

     integer, allocatable :: Mup(:,:)
     integer, allocatable :: Mdn(:,:)
     allocate(Mup(ncfgs/2,2))
     allocate(Mdn(ncfgs/2,2))

     iup = 0
     idn = 0

     do i=1,ncfgs
         call sector_locate(ia, i, nsect, ndims, sector_basis)

         ! c^+
         if ( bin_basis(iorb,i) == 0 ) then
             jold = dec_basis(i)
             call atomic_make_cdagger(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call sector_locate(ib, j, nsect, ndims, sector_basis)

             iup = iup + 1
             Mup(iup,1) = ia
             Mup(iup,2) = ib
         endif

         ! c
         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call sector_locate(ib, j, nsect, ndims, sector_basis)

             idn = idn + 1
             Mdn(idn,1) = ia
             Mdn(idn,2) = ib
         endif
     enddo

     call s_assert(iup == ncfgs / 2)
     call s_assert(idn == ncfgs / 2)

     print *, '# orb: ', iorb
     print *, 'number of Mup:', iup
     print *, 'number of Mdn:', idn


     do i=1,ncfgs/2
         HA = Mup(i,1)
         HL = Mup(i,1)
         HU = Mup(i,2)

         if ( HL /= 0 .and. HU /= 0 ) then
             call zigzag(1, HA, HL, HU, nsect, ndims, sector_basis, Mup, Mdn)
         endif
     enddo

     deallocate(Mup)
     deallocate(Mdn)

     return
  end subroutine sector_refine

  subroutine sector_locate(sind, find, nsect, ndims, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: nsect
     integer, intent(in) :: ndims(nsect)
     integer, intent(in) :: sector_basis(ncfgs,nsect)

     integer :: m
     integer :: n

     sind = 0

     SECTOR: do m=1,nsect
         do n=1,ndims(m)
             if ( sector_basis(n,m) == find ) then
                 sind = m
                 EXIT SECTOR
             endif
         enddo
     enddo SECTOR

     call s_assert(sind /= 0)

     return
  end subroutine sector_locate

  subroutine sector_lookup(sind, find)
     use m_sector, only : nsectors, sectors

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find

     integer :: m
     integer :: n

     sind = 0

     SECTOR: do m=1,nsectors
         do n=1,sectors(m)%ndim
             if ( sectors(m)%basis(n) == find ) then
                 sind = m
                 EXIT SECTOR
             endif
         enddo
     enddo SECTOR

     call s_assert(sind /= 0)

     return
  end subroutine sector_lookup

  subroutine sector_merge(ia, ib, nsect, ndims, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: ia
     integer, intent(in) :: ib
     integer, intent(in) :: nsect
     integer, intent(inout) :: ndims(nsect)
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

     integer :: m

     call s_assert(ndims(ia) >= 1)
     call s_assert(ndims(ib) >= 1)

     do m=1,ndims(ib)
           sector_basis(ndims(ia) + m, ia) = sector_basis(m, ib)
     enddo
     !
     ndims(ia) = ndims(ia) + ndims(ib)
     ndims(ib) = 0
     sector_basis(:,ib) = 0

     return
  end subroutine sector_merge

  subroutine print_sector(nsect, ndims, sector_basis)
     use constants, only : mystd

     use control, only : ncfgs

     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: ndims(nsect)
     integer, intent(in) :: sector_basis(ncfgs,nsect)

     integer :: i
     integer :: j
     integer :: m

     m = 0
     do i=1,nsect
         if ( ndims(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', ndims(i)
             write(mystd,'(a)') 'basis :'
             do j=1,ndims(i)
                 write(mystd,'(14i1)') bin_basis(:,sector_basis(j,i))
             enddo
             write(mystd,*)
         endif
     enddo

     return
  end subroutine print_sector

recursive &
  subroutine zigzag(up_or_down, HA, HL, HU, nsect, ndims, sector_basis, Mup, Mdn)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: up_or_down
     integer, intent(in) :: HA
     integer, intent(in) :: HL
     integer, intent(in) :: HU
     integer, intent(in) :: nsect
     integer, intent(inout) :: ndims(nsect)
     integer, intent(inout) :: sector_basis(ncfgs,nsect)
     integer, intent(inout) :: Mup(ncfgs/2,2)
     integer, intent(inout) :: Mdn(ncfgs/2,2)

     integer :: i, j
     integer :: HB

     if ( up_or_down == 1 ) then
         do i=1,ncfgs/2
             if ( Mup(i,1) == HA ) then
                 HB = Mup(i,2)
                 Mup(i,:) = 0

                 if ( HB /= HU ) then
                     print *, 'merge up:', HA, HB, HU
                     do j=1,ndims(HB)
                         sector_basis(ndims(HU) + j, HU) = sector_basis(j,HB)
                     enddo
                     ndims(HU) = ndims(HU) + ndims(HB)
                     ndims(HB) = 0
                 endif

                 call zigzag(2, HB, HL, HU, nsect, ndims, sector_basis, Mup, Mdn)
             endif
         enddo
     else
         do i=1,ncfgs/2
             if ( Mdn(i,1) == HA ) then
                 HB = Mdn(i,2)
                 Mdn(i,:) = 0

                 if ( HB /= HL ) then
                     print *, 'merge dn:', HA, HB, HL
                     do j=1,ndims(HB)
                         sector_basis(ndims(HL) + j, HL) = sector_basis(j,HB)
                     enddo
                     ndims(HL) = ndims(HL) + ndims(HB)
                     ndims(HB) = 0
                 endif

                 call zigzag(1, HB, HL, HU, nsect, ndims, sector_basis, Mup, Mdn)
             endif
         enddo
     endif

     return
  end subroutine zigzag

  subroutine try_next(i, j, k, which_sect)
     use m_fock, only : dec_basis, ind_basis, bin_basis
     use m_sector, only : sectors

     implicit none

     integer, intent(in) :: i
     integer, intent(in) :: j
     integer, intent(in) :: k
     integer, intent(in) :: which_sect

     integer :: m
     integer :: n
     integer :: jold, jnew, isgn
     integer :: ibasis

     if ( k == 1 ) then
         do n=1,sectors(i)%ndim
             ibasis = sectors(i)%basis(n)
             if ( bin_basis(j,ibasis) == 0 ) then
                 jold = dec_basis(ibasis)
                 call atomic_make_cdagger(j, jold, jnew, isgn)
                 m = ind_basis(jnew)
                 call s_assert( count(sectors(which_sect)%basis == m) == 1 )
                 !print *, n, m
             endif
         enddo
         !print *, 'which_sect:', which_sect, sectors(which_sect)%basis
     endif

     if ( k == 0 ) then
         do n=1,sectors(i)%ndim
             ibasis = sectors(i)%basis(n)
             if ( bin_basis(j,ibasis) == 1 ) then
                 jold = dec_basis(ibasis)
                 call atomic_make_c(j, jold, jnew, isgn)
                 m = ind_basis(jnew)
                 call s_assert( count(sectors(which_sect)%basis == m) == 1 )
                 !print *, n, m
             endif
         enddo
         !print *, 'which_sect:', which_sect, sectors(which_sect)%basis
     endif

     return
  end subroutine try_next

!!
!! @sub get_ntot
!!
!!
  subroutine get_sector_ntot(GQN_N, ndims, sector_basis)
     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: ndims
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_N

     integer :: i
     integer :: basis(norbs)
     integer :: N

     GQN_N = 999
     do i=1,ndims
         basis = bin_basis(:,sector_basis(i))
         N = sum(basis)
         if ( i == 1 ) then
             GQN_N = N
         else
             if ( N /= GQN_N ) then
                 STOP "wrong in GQN(N)"
             endif
         endif
     enddo

     return
  end subroutine get_sector_ntot

!!
!! @sub get_sz
!!
!!
  subroutine get_sector_sz(GQN_Sz, ndims, sector_basis)
     use control, only : isoc
     use control, only : nband, norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: ndims
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Sz

     integer :: i
     integer :: basis(norbs)
     integer :: Sz

     GQN_Sz = 999
     if ( isoc == 1 ) return
     do i=1,ndims
         basis = bin_basis(:,sector_basis(i))
         Sz = sum(basis(1:nband)) - sum(basis(nband+1:norbs))
         if ( i == 1 ) then
             GQN_Sz = Sz
         else
             if ( Sz /= GQN_Sz ) then
                 STOP "wrong in GQN(Sz)"
             endif
         endif
     enddo

     return
  end subroutine get_sector_sz

!!
!! @sub get_jz
!!
!!
  subroutine get_sector_jz(GQN_Jz, ndims, sector_basis)
     use control, only : isoc
     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: ndims
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Jz

     integer :: i, k
     integer :: basis(norbs)
     integer :: Jz
     integer :: good_jz(norbs)

     call atomic_make_gjz(good_jz)

     GQN_Jz = 999
     if ( isoc == 0 ) return

     do i=1,ndims
         basis = bin_basis(:,sector_basis(i))

         Jz = 0
         do k=1,norbs
             Jz = Jz + good_jz(k) * basis(k)
         enddo ! over k={1,norbs} loop

         if ( i == 1 ) then
             GQN_Jz = Jz
         else
             if ( Jz /= GQN_Jz ) then
                 STOP "wrong in GQN(Jz)"
             endif
         endif
     enddo

     return
  end subroutine get_sector_jz

!!
!! @sub get_ap
!!
!!
  subroutine get_sector_ap(Ap, i, N, Sz, Jz, nsect, sect_ntot, sect_sz, sect_jz)
     use control, only : isoc, icf

     implicit none

     integer, intent(out) :: Ap
     integer, intent(in) :: i
     integer, intent(in) :: N
     integer, intent(in) :: Sz
     integer, intent(in) :: Jz
     integer, intent(in) :: nsect
     integer, intent(in) :: sect_ntot(nsect)
     integer, intent(in) :: sect_sz(nsect)
     integer, intent(in) :: sect_jz(nsect)

     integer :: j

     Ap = 0

     if ( isoc == 0 ) then
         Ap = 1
         do j=1,i-1
             if ( ( sect_ntot(j) == N ) .and. ( sect_sz(j) == Sz ) ) then
                 Ap = Ap + 1
             endif
         enddo
     endif

     if ( isoc == 1 .and. icf == 0 ) then
         Ap = 1
         do j=1,i-1
             if ( ( sect_ntot(j) == N ) .and. ( sect_jz(j) == Jz ) ) then
                 Ap = Ap + 1
             endif
         enddo
     endif

     if ( isoc == 1 .and. icf == 1 ) then
         Ap = 1
         do j=1,i-1
             if ( sect_ntot(j) == N ) then
                 Ap = Ap + 1
             endif
         enddo
     endif

     return
  end subroutine get_sector_ap
