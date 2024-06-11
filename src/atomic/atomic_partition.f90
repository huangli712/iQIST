  subroutine automatic_partition()
     use constants, only : dp, mystd
     use constants, only : zero

     use control, only : ncfgs
     use control, only : norbs

     use m_fock, only : hmat, bin_basis

     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : cat_alloc_sector
     use m_sector, only : cat_alloc_sectors

     implicit none

     integer :: i
     integer :: j
     integer :: k
     integer :: ia, ib
     integer :: iorb
     integer :: nsect, nsize, nsect_, nsize_
     integer :: N, Sz, Jz, Ap

     ! index of Fock state
     integer :: ibasis

     integer, external :: get_nsect
     integer, external :: get_nsize

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
     ndims_ = 1
     sector_basis_ = 0
     do i=1,ncfgs
         sector_basis_(1,i) = i
     enddo

     ! phase 1
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(hmat(i,j)) > zero ) then
                 call locate_sector(ia, i, ncfgs, ndims_, sector_basis_)
                 call locate_sector(ib, j, ncfgs, ndims_, sector_basis_)

                 if ( ia /= ib ) then
                     call merge_sector(ia, ib, ncfgs, ndims_, sector_basis_) 
                 endif
             endif
         enddo
     enddo

     nsect_ = get_nsect(ncfgs, ndims_)
     nsize_ = get_nsize(ncfgs, ndims_)
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
     do iorb=1,norbs
         call refine_sector(iorb, nsect_, ndims, sector_basis)
     enddo

     nsect = get_nsect(nsect_, ndims)
     nsize = get_nsize(nsect_, ndims)
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
             call atomic_sector_N(N, ndims(i), sector_basis(:,i))
             call atomic_sector_Sz(Sz, ndims(i), sector_basis(:,i))
             call atomic_sector_Jz(Jz, ndims(i), sector_basis(:,i))

             Ap = 1
             do j=1,i-1
                 if ( ( sect_ntot(j) == N ) .and. ( sect_sz(j) == Sz ) ) then
                     Ap = Ap + 1
                 endif
             enddo

             sect_ntot(i) = N
             sect_sz(i) = Sz
             sect_jz(i) = Jz
             sect_ap(i) = Ap
         endif
     enddo
     call s_assert(k == nsect)

     !k = 0
     !do i=1,nsect_
     !    if ( ndims(i) > 0 ) then
     !        k = k + 1
     !        write(mystd,'(a,i6)') 'subspace -> ', k
     !        write(mystd,'(a,i6)') 'size :', ndims(i)
     !        write(mystd,'(a)') 'basis :'
     !        do j=1,ndims(i)
     !            write(mystd,'(i,2X,14i1)') j, bin_basis(:,sector_basis(j,i))
     !        enddo
     !        write(mystd, '(a, i3)') 'N :', sect_ntot(i)
     !        write(mystd, '(a, i3)') 'Sz:', sect_sz(i)
     !        write(mystd, '(a, i3)') 'Jz:', sect_jz(i)
     !        write(mystd, '(a, i3)') 'AP:', sect_ap(i)
     !        write(mystd, *)
     !    endif 
     !enddo

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
                         code = bin_basis(:,ibasis)
                         can = .true.
                         EXIT
                     !
                     ! test annihilation fermion operator
                     else if ( k == 0 .and. bin_basis(j, ibasis) == 1 ) then
                         code = bin_basis(:,ibasis)
                         can = .true.
                         EXIT
                     !
                     endif ! back if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) block
                 enddo ! over l={1,sectors(i)%ndim} loop

                 ! if can == .true., it means that the fermion operator
                 ! can act on the given subspace. next, we would like to
                 ! figure out the resulting subspace.
                 if ( can .eqv. .true. ) then


                 endif  ! back if ( can == .true. ) block

                 ! setup the next array
                 sectors(i)%next(j,k) = which_sect

                 if (k == 1) then
                     write(mystd,'(4X,a,i2,a)', advance = 'no') 'f^+(alpha =', j, ')'
                     write(mystd,'(2X,a,i4)', advance = 'no') '|subspace>_i:', i
                     write(mystd,'(2X,a,i4)') '|subspace>_f:', which_sect
                 else
                     write(mystd,'(4X,a,i2,a)', advance = 'no') 'f  (alpha =', j, ')'
                     write(mystd,'(2X,a,i4)', advance = 'no') '|subspace>_i:', i
                     write(mystd,'(2X,a,i4)') '|subspace>_f:', which_sect
                 endif ! back if (k == 1) block

             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,nsectors} loop




     STOP

     return
  end subroutine automatic_partition

  subroutine refine_sector(iorb, nsect, ndims, sector_basis)
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
         call locate_sector(ia, i, nsect, ndims, sector_basis)

         ! c^+
         if ( bin_basis(iorb,i) == 0 ) then
             jold = dec_basis(i)
             call atomic_make_cdagger(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call locate_sector(ib, j, nsect, ndims, sector_basis)

             iup = iup + 1
             Mup(iup,1) = ia
             Mup(iup,2) = ib
         endif

         ! c
         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call locate_sector(ib, j, nsect, ndims, sector_basis)

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
  end subroutine refine_sector

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

  subroutine locate_sector(sind, find, nsect, ndims, sector_basis)
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
  end subroutine locate_sector

  subroutine merge_sector(ia, ib, nsect, ndims, sector_basis)
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
  end subroutine merge_sector

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

  function get_nsect(nsect, ndims) result(val)
     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: ndims(nsect)

     integer :: val

     integer :: i

     val = 0
     do i=1,nsect
         if ( ndims(i) > 0 ) then
             val = val + 1
         endif
     enddo

     return
  end function get_nsect

  function get_nsize(nsect, ndims) result(val)
     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: ndims(nsect)

     integer :: val

     val = maxval(ndims)

     return
  end function get_nsize

  subroutine atomic_sector_N(GQN_N, ndims, sector_basis)
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
  end subroutine atomic_sector_N

  subroutine atomic_sector_Sz(GQN_Sz, ndims, sector_basis)
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
  end subroutine atomic_sector_Sz

  subroutine atomic_sector_Jz(GQN_Jz, ndims, sector_basis)
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
  end subroutine atomic_sector_Jz
