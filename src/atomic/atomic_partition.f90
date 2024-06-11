  subroutine automatic_partition()
     use constants, only : dp
     use constants, only : zero

     use control, only : ncfgs
     use control, only : norbs

     use m_fock, only : hmat

     implicit none

     integer :: i
     integer :: j
     integer :: ia, ib
     integer :: iorb
     integer :: nsect, nsize, nsect_, nsize_
     integer :: N, Sz, Jz, Ap

     integer, external :: get_nsect
     integer, external :: get_nsize

     integer, allocatable :: sector_size(:)
     integer, allocatable :: sector_size_(:)
     integer, allocatable :: sector_basis(:,:)
     integer, allocatable :: sector_basis_(:,:)

     integer, allocatable :: sect_ntot(:)
     integer, allocatable :: sect_sz(:)
     integer, allocatable :: sect_jz(:)
     integer, allocatable :: sect_ap(:)

     ! initialization
     allocate(sector_size_(ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))
     sector_size_ = 1
     sector_basis_ = 0
     do i=1,ncfgs
         sector_basis_(1,i) = i
     enddo

     ! phase 1
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(hmat(i,j)) > zero ) then
                 call locate_sector(ia, i, ncfgs, sector_size_, sector_basis_)
                 call locate_sector(ib, j, ncfgs, sector_size_, sector_basis_)

                 if ( ia /= ib ) then
                     call merge_sector(ia, ib, ncfgs, sector_size_, sector_basis_) 
                 endif
             endif
         enddo
     enddo

     nsect_ = get_nsect(ncfgs, sector_size_)
     nsize_ = get_nsize(ncfgs, sector_size_)
     print *, 'number of sectors: ', nsect_
     print *, 'maximum size of sectors: ', nsize_

     ! filter sectors
     allocate(sector_size(nsect_))
     allocate(sector_basis(ncfgs,nsect_))
     !
     j = 0
     do i=1,ncfgs
         if ( sector_size_(i) > 0 ) then
             j = j + 1
             sector_size(j) = sector_size_(i)
             sector_basis(:,j) = sector_basis_(:,i)
         endif
     enddo
     !
     call s_assert(j == nsect_)
     deallocate(sector_size_)
     deallocate(sector_basis_)

     ! phase 2
     do iorb=1,norbs
         call refine_sector(iorb, nsect_, sector_size, sector_basis)
     enddo

     nsect = get_nsect(nsect_, sector_size)
     nsize = get_nsize(nsect_, sector_size)
     print *, 'number of sectors: ', nsect
     print *, 'maximum size of sectors: ', nsize

     allocate(sect_ntot(nsect_))
     allocate(sect_sz(nsect_))
     allocate(sect_jz(nsect_))
     allocate(sect_ap(nsect_))
     sect_ntot = 0
     sect_sz = 0
     sect_jz = 0
     sect_ap = 0

     do i=1,nsect_
         if ( sector_size(i) > 0 ) then
             call atomic_sector_N(N, sector_size(i), sector_basis(:,i))
             call atomic_sector_Sz(Sz, sector_size(i), sector_basis(:,i))
             call atomic_sector_Jz(Jz, sector_size(i), sector_basis(:,i))

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
 
     STOP

     return
  end subroutine automatic_partition

  subroutine refine_sector(iorb, nsect, sector_size, sector_basis)
     use control, only : ncfgs

     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

     integer, intent(in) :: iorb
     integer, intent(in) :: nsect
     integer, intent(inout) :: sector_size(nsect)
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
         call locate_sector(ia, i, nsect, sector_size, sector_basis)

         ! c^+
         if ( bin_basis(iorb,i) == 0 ) then
             jold = dec_basis(i)
             call atomic_make_cdagger(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call locate_sector(ib, j, nsect, sector_size, sector_basis)

             iup = iup + 1
             Mup(iup,1) = ia
             Mup(iup,2) = ib
         endif

         ! c
         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call locate_sector(ib, j, nsect, sector_size, sector_basis)

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
             call zigzag(1, HA, HL, HU, nsect, sector_size, sector_basis, Mup, Mdn)
         endif
     enddo

     deallocate(Mup)
     deallocate(Mdn)

     return
  end subroutine refine_sector

recursive &
  subroutine zigzag(up_or_down, HA, HL, HU, nsect, sector_size, sector_basis, Mup, Mdn)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: up_or_down
     integer, intent(in) :: HA
     integer, intent(in) :: HL
     integer, intent(in) :: HU
     integer, intent(in) :: nsect
     integer, intent(inout) :: sector_size(nsect)
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
                     do j=1,sector_size(HB)
                         sector_basis(sector_size(HU) + j, HU) = sector_basis(j,HB)
                     enddo
                     sector_size(HU) = sector_size(HU) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(2, HB, HL, HU, nsect, sector_size, sector_basis, Mup, Mdn)
             endif
         enddo 
     else
         do i=1,ncfgs/2
             if ( Mdn(i,1) == HA ) then
                 HB = Mdn(i,2)
                 Mdn(i,:) = 0

                 if ( HB /= HL ) then
                     print *, 'merge dn:', HA, HB, HL
                     do j=1,sector_size(HB)
                         sector_basis(sector_size(HL) + j, HL) = sector_basis(j,HB)
                     enddo
                     sector_size(HL) = sector_size(HL) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(1, HB, HL, HU, nsect, sector_size, sector_basis, Mup, Mdn)
             endif
         enddo 
     endif

     return
  end subroutine zigzag

  subroutine locate_sector(sind, find, nsect, sector_size, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)
     integer, intent(in) :: sector_basis(ncfgs,nsect)

     integer :: m
     integer :: n

     sind = 0

     SECTOR: do m=1,nsect
         do n=1,sector_size(m)
             if ( sector_basis(n,m) == find ) then
                 sind = m
                 EXIT SECTOR
             endif
         enddo
     enddo SECTOR

     call s_assert(sind /= 0)

     return
  end subroutine locate_sector

  subroutine merge_sector(ia, ib, nsect, sector_size, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: ia
     integer, intent(in) :: ib
     integer, intent(in) :: nsect
     integer, intent(inout) :: sector_size(nsect)
     integer, intent(inout) :: sector_basis(ncfgs,nsect)

     integer :: m

     call s_assert(sector_size(ia) >= 1)
     call s_assert(sector_size(ib) >= 1)

     do m=1,sector_size(ib)
           sector_basis(sector_size(ia) + m, ia) = sector_basis(m, ib)
     enddo
     !
     sector_size(ia) = sector_size(ia) + sector_size(ib)
     sector_size(ib) = 0
     sector_basis(:,ib) = 0

     return
  end subroutine merge_sector

  subroutine print_sector(nsect, sector_size, sector_basis)
     use constants, only : mystd

     use control, only : ncfgs

     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)
     integer, intent(in) :: sector_basis(nsect,ncfgs)

     integer :: i
     integer :: j
     integer :: m

     m = 0
     do i=1,nsect
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a)') 'basis :'
             do j=1,sector_size(i)
                 write(mystd,'(14i1)') bin_basis(:,sector_basis(i,j))
             enddo
             write(mystd,*)
         endif
     enddo

     return
  end subroutine print_sector

  subroutine print_sector_new(nsect, sector_size, sector_basis)
     use constants, only : mystd

     use control, only : ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)
     integer, intent(in) :: sector_basis(nsect, ncfgs)

     integer :: m
     integer :: i
     integer :: j

     integer :: N
     integer :: Sz
     integer :: Ap
     integer :: sector_N(nsect)
     integer :: sector_Sz(nsect)
     integer :: sector_Ap(nsect)

     sector_N = 0
     sector_Sz = 0
     sector_Ap = 0

     ! print subspaces
     m = 0
     do i=1,nsect
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a)') 'basis :'
             do j=1,sector_size(i)
                 write(mystd,'(i,2X,14i1)') j, bin_basis(:,sector_basis(i,j))
             enddo
             call atomic_sector_N(N, sector_size(i), sector_basis(i,:))
             call atomic_sector_Sz(Sz, sector_size(i), sector_basis(i,:))
             !call atomic_sector_Jz(Jz, sector_size(i), sector_basis(i,:))
             Ap = 1
             do j=1,m-1
                 if ( ( sector_N(j) == N ) .and. ( sector_Sz(j) == Sz ) ) then
                     Ap = Ap + 1
                 endif
             enddo
             write(mystd, '(a, i3)') 'N :', N
             write(mystd, '(a, i3)') 'Sz:', Sz
             !write(mystd, '(a, i3)') 'Jz:', Jz
             write(mystd, '(a, i3)') 'AP:', Ap
             sector_N(m) = N
             sector_Sz(m) = Sz
             !sector_Jz(m) = Jz
             sector_Ap(m) = Ap
             print *
         endif
     enddo
 
     return
  end subroutine print_sector_new

  function get_nsect(nsect, sector_size) result(val)
     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)

     integer :: val

     integer :: i

     val = 0
     do i=1,nsect
         if ( sector_size(i) > 0 ) then
             val = val + 1
         endif
     enddo

     return
  end function get_nsect

  function get_nsize(nsect, sector_size) result(val)
     implicit none

     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)

     integer :: val

     val = maxval(sector_size)

     return
  end function get_nsize

  subroutine atomic_sector_N(GQN_N, sector_size, sector_basis)
     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_N

     integer :: i
     integer :: basis(norbs)
     integer :: N

     GQN_N = 999
     do i=1,sector_size
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

  subroutine atomic_sector_Sz(GQN_Sz, sector_size, sector_basis)
     use control, only : isoc
     use control, only : nband, norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Sz

     integer :: i
     integer :: basis(norbs)
     integer :: Sz

     GQN_Sz = 999
     if ( isoc == 1 ) return
     do i=1,sector_size
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

  subroutine atomic_sector_Jz(GQN_Jz, sector_size, sector_basis)
     use control, only : isoc
     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Jz

     integer :: i, k
     integer :: basis(norbs)
     integer :: Jz
     integer :: good_jz(norbs)

     call atomic_make_gjz(good_jz)

     GQN_Jz = 999
     if ( isoc == 0 ) return

     do i=1,sector_size
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
