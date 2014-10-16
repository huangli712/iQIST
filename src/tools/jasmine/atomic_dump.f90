!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_dump_basis
!!!           atomic_dump_natural
!!!           atomic_dump_eimp
!!!           atomic_dump_feigval
!!!           atomic_dump_feigvec
!!!           atomic_dump_fcix
!!!           atomic_dump_seigval
!!!           atomic_dump_seigvec
!!!           atomic_dump_scix
!!! source  : atomic_dump.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : write output files
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_dump_basis: write basis of fullspace to file 'atom.basis.dat'
  subroutine atomic_dump_basis()
     use constants, only : mytmp
     use control, only : ncfgs
     use m_full, only : dec_basis, index_basis, bin_basis
  
     implicit none
  
! local variables
     integer :: i
  
! open file 'atom.basis.dat' to write
     open(mytmp, file='atom.basis.dat')
! write the header
     write(mytmp, '(a)') '#      i |  decimal(i) |    index(i) |      binary(i) |'
     do i=1, ncfgs
         write(mytmp, "(I10,4X,I10,4X,I10,8X,14I1)") i, dec_basis(i), &
                               index_basis(dec_basis(i)), bin_basis(:,i)   
     enddo 
     close(mytmp)
  
     return
  end subroutine atomic_dump_basis

!!>>> atomic_write_natural: write the transformation matrix 
!!>>> from the original basis to natural basis
  subroutine atomic_dump_natural(info)
     use constants, only : dp, mytmp
     use control, only : norbs
     use m_spmat, only : tmat
  
! external variables
     character(len=*), intent(in) :: info
  
! local variables
     integer :: i,j
  
     open(mytmp, file='atom.natural.dat')
     write(mytmp,'(a)') info
     write(mytmp,'(a)') '#      i |       j |    umat_real(i,j) |    umat_imag(i,j) |'
     do i=1, norbs
         do j=1, norbs
             write(mytmp, '(2I10,2F20.10)') j, i, tmat(j,i)
         enddo
     enddo 
     close(mytmp)
  
     return
  end subroutine atomic_dump_natural

!!>>> atomic_dump_eimp: write on-site impurity energy on natural basis 
  subroutine atomic_dump_eimp()
     use constants, only : mytmp
     use control, only : nband, norbs, isoc

     use m_spmat, only : eimpmat
  
     implicit none
  
! local variables
     integer :: i
     integer :: s_order
 
     open(mytmp, file='atom.eimp.dat') 
     do i=1, norbs
         if (isoc==0) then
             if (i <= nband) then
                 s_order = 2*i-1  
             else
                 s_order = 2*(i-nband)
             endif
         else 
             s_order = i
         endif 
         write(mytmp, '(I10,F20.10)') i, real(eimpmat(s_order,s_order))
     enddo
     close(mytmp)
  
     return
  end subroutine atomic_dump_eimp

!!>>> atomic_dump_feigval: write eigenvalue of fullspace 
!!>>> to file 'atom.eigval.dat'
  subroutine atomic_dump_feigval()
     use constants, only : mytmp
     use control, only : ncfgs
     use m_full, only : hmat_eigval, occu_mat
  
     implicit none
  
! local variables
     integer :: i
     
! open file 'atom.eigval.dat' to write
     open(mytmp, file='atom.eigval.dat')
     write(mytmp, '(a)') '#       i |     eigenvalue(i) |      occupancy(i) |'
     do i=1, ncfgs
         write(mytmp, "(I10, 2F20.10)") i, hmat_eigval(i), occu_mat(i,i)
     enddo
     close(mytmp)
  
     return
  end subroutine atomic_dump_feigval

!!>>> atomic_dump_feigvec: write eigenvector of 
!!>>> fullspace to file 'atom.eigvec.dat'
  subroutine atomic_dump_feigvec()
     use constants, only : mytmp, eps6
     use control, only : ncfgs

     use m_full, only : bin_basis
     use m_full, only : hmat_eigvec
  
     implicit none
  
! local variables
     integer :: i, j
   
! open file 'atom.eigvec.dat' to write
     open(mytmp, file='atom.eigvec.dat')
     write(mytmp, '(a)') '#      i |       j |  eigenvector(i,j) |       fock_basis(i) |'
     do i=1, ncfgs
         do j=1, ncfgs
             if ( abs(hmat_eigvec(j,i)) < eps6 ) cycle
             write(mytmp, "(2I10, F20.10, 8X, 14I1)") j, i, hmat_eigvec(j,i), bin_basis(:,j) 
         enddo
     enddo 
     close(mytmp)
  
     return
  end subroutine atomic_dump_feigvec

!!>>> atomic_dump_fcix: write atom.cix for CTQMC input
!!>>> for ictqmc == 1 case.
  subroutine atomic_dump_fcix()
     use constants, only : mytmp, zero
     use control, only : nband, norbs, ncfgs, isoc

     use m_full, only : hmat_eigval, occu_mat, anni_fmat
  
     implicit none
  
! local variables
     integer :: i, j, k
     integer :: s_order
  
! open file 'atom.cix' to write
     open(mytmp, file='atom.cix')
  
! write eigenvalues
     write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(I10,3F20.10)') i, hmat_eigval(i), occu_mat(i,i), zero
     enddo 
  
! write fmat
     write(mytmp,'(a)') '# f matrix element: alpha | beta | orbital | fmat'
     do i=1,norbs
! for non-soc case, the spin order of CTQMC is like up, up, up, dn, dn, dn
! but the spin order of this program is up, dn, up, dn, up, dn
! we adjust it here
! for soc case, it doesn't matter
         if (isoc==0) then
             if (i <= nband) then
                 s_order = 2*i-1  
             else
                 s_order = 2*(i-nband)
             endif
         else 
             s_order = i
         endif 
  
         do j=1,ncfgs
             do k=1,ncfgs
                 write(mytmp,'(3I10, F20.10)') k, j, i, anni_fmat(k, j, s_order)
             enddo 
         enddo 
     enddo 
     close(mytmp)
  
     return
  end subroutine atomic_dump_fcix

!!>>> atomic_dump_seigval: write eigenvalue of sectors 
!!>>> to file 'atom.eigval.dat'
  subroutine atomic_dump_seigval()
     use constants, only : mytmp

     use m_sector, only : nsectors, sectors
  
     implicit none
  
! local variables
     integer :: i, j
     integer :: counter
  
! open file 'atom.eigval.dat' to write
     open(mytmp, file='atom.eigval.dat')
     counter = 0
     write(mytmp, '(a)') "#      i |     sect(i) | electron(i) |           j |   eigenvalue(j,i) |"
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             counter = counter + 1
             write(mytmp, "(I10,4X,I10,4X,I10,4X,I10, F20.10)") counter, i, &
                               sectors(i)%nelectron, j, sectors(i)%myeigval(j)
         enddo
     enddo
     close(mytmp)
  
     return
  end subroutine atomic_dump_seigval

!!>>> atomic_dump_seigvec: write eigenvector of sectors 
!!>>> to file 'atom.eigval.dat'
  subroutine atomic_dump_seigvec()
     use constants, only : mytmp, eps6

     use m_full, only : bin_basis
     use m_sector, only : nsectors, sectors
  
     implicit none
  
! local variables
     integer :: i, j, k
     integer :: counter
  
     open(mytmp, file="atom.eigvec.dat")
     counter = 0
     write(mytmp, '(a)') '#      i |       j |       k |       eigvec(j,k) |       fock_basis(j) |'
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             do k=1, sectors(i)%ndim
                 if ( abs(sectors(i)%myeigvec(k,j)) < eps6 ) cycle
                     write(mytmp, "(3I10, F20.10, 8X, 14I1)") i, k+counter, j+counter, &
                     sectors(i)%myeigvec(k,j), bin_basis(:,sectors(i)%mybasis(k))
             enddo
         enddo
         counter = counter + sectors(i)%ndim
     enddo
     close(mytmp)
  
     return
  end subroutine atomic_dump_seigvec

!!>>> atomic_dump_scix: write atom.cix for CTQMC input, 
!!>>> for good quantum numbers (GQNs) algorithm
  subroutine atomic_dump_scix()
     use constants, only : mytmp
     use control, only : Uc, Uv, Jz, Js, Jp, Ud, JH, F0, F2, F4, F6, icu, isoc

     use m_sector, only : nsectors, sectors, max_dim_sect, ave_dim_sect
  
     implicit none
  
! local variables
     integer :: i, j, k, ii
     integer :: s_order
   
! open 'atom.cix' to write
     open(mytmp, file='atom.cix')
! write header
     write(mytmp, "(a)") "# This file is generated by JASMINE, can be &
                             ONLY used by PANSY and MANJUSHAKA CTQMC codes."
     write(mytmp, "(a)") "# WARNING: please DON'T change this file manually !"
     if (icu==1) then
         write(mytmp, "(a)") "# Kanamori-parameterized Coulomb interaction: "
         write(mytmp, "(5(a,F9.5,2X))") "# Uc:", Uc, "Uv:", Uv, "Jz:", Jz, "Js:", Js, "Jp:", Jp
     else
         write(mytmp, "(a)") "# Slater-parameterized Coulomb interaction: "
         write(mytmp, "(6(a,F9.5,2X))") "# Ud:", Ud, "JH:", JH, "F0:", F0, "F2:", F2, "F4:", F4, "F6:", F6
     endif
     write(mytmp, "(a)") "#================================================================================"

! write number of sectors
     write(mytmp, "(a)") "# NUMBER OF SECTORS | MAXIMUM DIMENSION OF SECTORS | AVERAGE DIMENSION OF SECTORS"
     write(mytmp, "(I10,15X,I10, 15X, F20.10)") nsectors, max_dim_sect, ave_dim_sect 
  
! write dimension, total electrons, next_sector, eigenvalue of each sector
     do i=1, nsectors
         write(mytmp, "(a)") "# SECT_INFO: INDEX  |  NDIM  |  NELEC  |   NOPS  |  ISTART"  
         write(mytmp, "(6X,5I10)") i, sectors(i)%ndim, sectors(i)%nelectron, sectors(i)%nops, sectors(i)%istart
  
! write next_sector
         write(mytmp, "(4X,a)") "# NEXT_SECTOR    F     F^{\dagger}"
         do j=1, sectors(i)%nops
! adjust the orbital order for CTQMC, up, up, up, dn, dn, dn
             if (isoc==0) then
                 if (j <= sectors(i)%nops / 2) then
                     s_order = 2*j-1  
                 else
                     s_order = 2*(j - sectors(i)%nops / 2)
                 endif
             else 
                 s_order = j
             endif 
             write(mytmp, "(2X, 3I10)") j, sectors(i)%next_sector(s_order,0), sectors(i)%next_sector(s_order,1)  
         enddo
  
! write eigeanvalue
         write(mytmp, "(4X,a)") "# EIGENVALUES"
         do j=1, sectors(i)%ndim
             write(mytmp, "(2X,I10, F20.10)") j, sectors(i)%myeigval(j) 
         enddo
     enddo
     close(mytmp)
  
! write fmat
     open(mytmp, file='atom.fmat', form='unformatted')
     do i=1, nsectors
         do j=1, sectors(i)%nops
             if (isoc==0) then
                 if (j <= sectors(i)%nops / 2) then
                     s_order = 2*j-1  
                 else
                     s_order = 2*(j-sectors(i)%nops / 2)
                 endif
             else 
                 s_order = j
             endif 
             do k=0,1
                 ii = sectors(i)%next_sector(s_order,k)
                 if (ii == -1) cycle 
                 write(mytmp)  sectors(i)%myfmat(s_order,k)%item(:,:)
             enddo  ! over k={0,1} loop
         enddo ! over j={1, sectors(i)%nops} loop
     enddo  ! over i={1, nsect} loop
     close(mytmp)
  
     return
  end subroutine atomic_dump_scix
