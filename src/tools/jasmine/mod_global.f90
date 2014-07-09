!-------------------------------------------------------------------------
! project : jasmine
! program : m_basis_fullspace
!         : m_spmat 
!         : m_glob_fullspace
!         : m_glob_sectors
! source  : mod_global.f90
! type    : modules 
! authors : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014
! purpose : global variables
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> Fock basis of full Hilbert space
module m_basis_fullspace
    use control, only: norbs, ncfgs

    implicit none

    ! dimension of total electron N's subspace 
    integer, public, allocatable, save :: dim_sub_n(:)

    ! binary form of Fock basis
    integer, public, allocatable, save :: bin_basis(:,:)

    ! decimal form of Fock basis
    integer, public, allocatable, save :: dec_basis(:)

    ! index of Fock basis, given their decimal number
    integer, public, allocatable, save :: index_basis(:)

    contains

    !>>> allocate memory for these matrices
    subroutine alloc_m_basis_fullspace()
        implicit none

        ! allocate them
        allocate(dim_sub_n(0:norbs))
        allocate(bin_basis(norbs, ncfgs))
        allocate(dec_basis(ncfgs))
        allocate(index_basis(0:ncfgs-1))

        ! init them
        dim_sub_n = 0
        bin_basis = 0
        dec_basis = 0
        index_basis = 0

        return
    end subroutine alloc_m_basis_fullspace

    ! deallocate memory for these matrices
    subroutine dealloc_m_basis_fullspace()
        implicit none

        ! deallocate them
        if (allocated(dim_sub_n))   deallocate(dim_sub_n)
        if (allocated(bin_basis))   deallocate(bin_basis)
        if (allocated(dec_basis))   deallocate(dec_basis) 
        if (allocated(index_basis)) deallocate(index_basis) 

        return
    end subroutine dealloc_m_basis_fullspace

end module m_basis_fullspace

!>>> single particle related matrices, including: 
! crystal field, spin-orbital coupling, Coulomb interaction U tensor
module m_spmat
    use constants, only: dp, czero
    use control,   only: norbs
    implicit none

    ! crystal field (CF)
    complex(dp), public, allocatable, save :: cfmat(:,:) 

    ! spin-orbital coupling (SOC)
    complex(dp), public, allocatable, save :: socmat(:,:)

    ! on-site energy (CF+SOC) of impurity
    complex(dp), public, allocatable, save :: eimpmat(:,:)

    ! Coulomb interaction U tensor
    complex(dp), public, allocatable, save :: cumat(:,:,:,:)

    ! the transformation matrix from origional basis to natural basis 
    complex(dp), public, allocatable, save :: tran_umat(:,:)

    contains

    !>>> allocate memory for these matrix
    subroutine alloc_m_spmat()
        implicit none

        ! allocate them
        allocate(cfmat(norbs, norbs))
        allocate(socmat(norbs, norbs))
        allocate(eimpmat(norbs, norbs))
        allocate(cumat(norbs, norbs, norbs, norbs))
        allocate(tran_umat(norbs, norbs))

        ! init them
        cfmat    = czero
        socmat   = czero
        eimpmat  = czero
        cumat    = czero
        tran_umat= czero

        return
    end subroutine alloc_m_spmat

    !>>> deallocate memory for these matrix
    subroutine dealloc_m_spmat()
        implicit none

        ! deallocate them
        if (allocated(cfmat))      deallocate(cfmat)
        if (allocated(socmat))     deallocate(socmat)
        if (allocated(eimpmat))    deallocate(eimpmat)
        if (allocated(cumat))      deallocate(cumat)
        if (allocated(tran_umat))  deallocate(tran_umat)

        return
    end subroutine dealloc_m_spmat

end module m_spmat

!>>> global variables for full space algorithm case
module m_glob_fullspace
    use constants,  only: dp, zero, czero
    use control,    only: norbs, ncfgs

    implicit none

    ! atomic Hamiltonian (CF + SOC + CU)
    complex(dp), public, allocatable, save :: hmat(:,:)

    ! eigen value of hmat
    real(dp), public, allocatable, save :: hmat_eigval(:)

    ! eigen vector of hmat
    real(dp), public, allocatable, save :: hmat_eigvec(:, :)

    ! fmat for annihilation fermion operators
    real(dp), public, allocatable, save :: anni_fmat(:,:,:)

    ! occupany number for atomic eigenstates
    real(dp), public, allocatable, save :: occu_mat(:,:)

    contains

    subroutine alloc_m_glob_fullspace()
        implicit none

        allocate(hmat(ncfgs, ncfgs))
        allocate(hmat_eigval(ncfgs))
        allocate(hmat_eigvec(ncfgs, ncfgs))
        allocate(anni_fmat(ncfgs, ncfgs, norbs))
        allocate(occu_mat(ncfgs, ncfgs))

        ! init them
        hmat = czero
        hmat_eigval = zero
        hmat_eigvec = zero
        anni_fmat = zero
        occu_mat = zero

        return
    end subroutine alloc_m_glob_fullspace

    subroutine dealloc_m_glob_fullspace()
        implicit none

        if(allocated(hmat))        deallocate(hmat)
        if(allocated(hmat_eigval)) deallocate(hmat_eigval)
        if(allocated(hmat_eigvec)) deallocate(hmat_eigvec)
        if(allocated(anni_fmat))   deallocate(anni_fmat)
        if(allocated(occu_mat))    deallocate(occu_mat)

        return
    end subroutine dealloc_m_glob_fullspace

end module m_glob_fullspace

!>>> global variables for sectors algorithm case
module m_glob_sectors
    use m_sector,  only: t_sector, nullify_one_sector, dealloc_one_sector

    implicit none

    ! number of sectors
    integer, public, save :: nsectors

    ! all the sectors
    type(t_sector), public, allocatable, save :: sectors(:)

    contains

    subroutine alloc_m_glob_sectors()
        implicit none

        ! local variables
        integer :: i

        allocate(sectors(nsectors)) 

        ! nullify each sector
        do i=1, nsectors
            call nullify_one_sector(sectors(i))
        enddo          

        return
    end subroutine alloc_m_glob_sectors

    subroutine dealloc_m_glob_sectors()
        implicit none
        
        integer :: i
        ! we will deallocate memory for pointers in type_subsapce 
        ! before deallocating subspaces to avoid memory leak 
   
        do i=1, nsectors
            call dealloc_one_sector(sectors(i)) 
        enddo 

        if (allocated(sectors))  deallocate(sectors) 
       
        return
    end subroutine dealloc_m_glob_sectors

end module m_glob_sectors
