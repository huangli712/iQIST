  subroutine automatic_partition()
     use constants, only : dp
     use constants, only : zero

     use control, only : ncfgs

     use m_fock, only : hmat
     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

     integer :: i
     integer :: j
     integer :: ia, ib
     integer :: nsect, nsize
     integer :: iorb, iup, idn
     integer :: jnew, jold, isgn
     integer, external :: get_nsect
     integer, external :: get_nsize

     integer, allocatable :: sector_size(:)
     integer, allocatable :: sector_size_(:)
     integer, allocatable :: sector_basis(:,:)
     integer, allocatable :: sector_basis_(:,:)

     allocate(sector_size_(ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))
     sector_size_ = 1
     sector_basis_ = 0
     do i=1,ncfgs
         sector_basis_(i,1) = i
     enddo

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

     call print_sector(ncfgs, sector_size_, sector_basis_)

     nsect = get_nsect(ncfgs, sector_size_)
     nsize = get_nsize(ncfgs, sector_size_)

     print *, 'number of sectors: ', nsect
     print *, 'maximum size of sectors: ', nsize

     allocate(sector_size(nsect))
     allocate(sector_basis(nsect,ncfgs))
     !
     j = 0
     do i=1,ncfgs
         if ( sector_size_(i) > 0 ) then
             j = j + 1
             sector_size(j) = sector_size_(i)
             sector_basis(j,:) = sector_basis_(i,:)
         endif
     enddo
     !
     call s_assert(j == nsect)

     STOP

     return
  end subroutine automatic_partition

  subroutine refine_sector()
     implicit none

     iorb = 2
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
             !call s_assert(iup <= max_mapping)
             !Mup(iup,1) = ia
             !Mup(iup,2) = ib
         endif

         ! c
         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call locate_sector(ib, j, nsect, sector_size, sector_basis)

             idn = idn + 1
             !call s_assert(idn <= max_mapping)
             !Mdn(idn,1) = ia
             !Mdn(idn,2) = ib
         endif
     enddo
     print *, '# orb: ', iorb
     print *, 'number of Mup:', iup
     print *, 'number of Mdn:', idn

     return
  end subroutine refine_sector

  subroutine locate_sector(sind, find, nsect, sector_size, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: nsect
     integer, intent(in) :: sector_size(nsect)
     integer, intent(in) :: sector_basis(nsect,ncfgs)

     integer :: m
     integer :: n

     sind = 0

     SECTOR: do m=1,nsect
         do n=1,sector_size(m)
             if ( sector_basis(m,n) == find ) then
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
     integer, intent(inout) :: sector_basis(nsect,ncfgs)

     integer :: m

     call s_assert(sector_size(ia) >= 1)
     call s_assert(sector_size(ib) >= 1)

     do m=1,sector_size(ib)
           sector_basis(ia, sector_size(ia) + m) = sector_basis(ib,m)
     enddo
     !
     sector_size(ia) = sector_size(ia) + sector_size(ib)
     sector_size(ib) = 0
     sector_basis(ib,:) = 0

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
