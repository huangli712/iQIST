  subroutine automatic_partition()
     use constants, only : dp

     use control, only : ncfgs

     use m_fock, only : hmat

     implicit none

     integer :: i
     integer :: j
     integer :: ia, ib

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
             if ( abs(hmat(i,j)) > 0.0_dp ) then
                 call locate_sector(ia, i, ncfgs, ncfgs, sector_size_, sector_basis_)
                 call locate_sector(ib, j, ncfgs, ncfgs, sector_size_, sector_basis_)

                 if ( ia /= ib ) then
                     call merge_sector(ia, ib, ncfgs, ncfgs, sector_size_, sector_basis_) 
                 endif
             endif
         enddo
     enddo

     call print_sector(ncfgs, ncfgs, sector_size_, sector_basis_)

     STOP

     return
  end subroutine automatic_partition

  subroutine locate_sector(sind, find, max_num, max_size, sector_size, sector_basis)
     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: max_num
     integer, intent(in) :: max_size
     integer, intent(in) :: sector_size(max_num)
     integer, intent(in) :: sector_basis(max_num,max_size)

     integer :: m
     integer :: n

     sind = 0

     SECTOR: do m=1,max_num
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

  subroutine merge_sector(ia, ib, max_num, max_size, sector_size, sector_basis)
     implicit none

     integer, intent(in) :: ia
     integer, intent(in) :: ib
     integer, intent(in) :: max_num
     integer, intent(in) :: max_size
     integer, intent(inout) :: sector_size(max_num)
     integer, intent(inout) :: sector_basis(max_num,max_size)

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

  subroutine print_sector(max_num, max_size, sector_size, sector_basis)
     use constants, only : mystd

     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: max_num
     integer, intent(in) :: max_size
     integer, intent(in) :: sector_size(max_num)
     integer, intent(in) :: sector_basis(max_num,max_size)

     integer :: i
     integer :: j
     integer :: m

     m = 0
     do i=1,max_num
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

   function get_max_num(max_num, sector_size) result(val)
     implicit none

     integer, intent(in) :: max_num
     integer, intent(in) :: sector_size(max_num)

     integer :: val

     integer :: i

     val = 0
     do i=1,max_num
         if ( sector_size(i) > 0 ) then
             val = val + 1
         endif
     enddo

     return
   end function get_max_num

   function get_max_size(max_num, sector_size) result(val)
     implicit none

     integer, intent(in) :: max_num
     integer, intent(in) :: sector_size(max_num)

     integer :: val

     val = maxval(sector_size)

     return
   end function get_max_size
