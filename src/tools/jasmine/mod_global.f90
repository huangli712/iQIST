!=========================================================================!
! project : clematis
! program : context module
! history : Apr 26, 2011
! authors : xidai and duliang {duleung@gmail.com}
! purpose : important arrays defined for main program 
! comment :
!=========================================================================!

!-------------------------------------------------------------------------!
!>>> basis states related matrix
!-------------------------------------------------------------------------!
module mod_glob_basis
    implicit none

    ! number of dimensions of each (isub) subspace
    integer, public, allocatable, save :: nstat(:)

    ! decimal representation of basis state in Fock space
    integer, public, allocatable, save :: basis(:)

    ! the jz number for single particle basis for spin-orbit coupling case 
    integer, public, allocatable, save :: good(:)

    ! serial number of a decimal representated configuration
    integer, public, allocatable, save :: invsn(:)

    ! binary representation of a decimal represented basis state
    integer, public, allocatable, save :: invcd(:, :)

end module mod_glob_basis

!-------------------------------------------------------------------------!
!>>> atomic Hamiltonian martix related variables
!-------------------------------------------------------------------------!
module mod_glob_hmtrx
    use constants, only: dp
    implicit none

    ! eigenvalues of atomic single particle Hamiltonian
    real(dp), public, allocatable, save :: eloc(:)

    ! eigenvalues of atomic many particle Hamiltonian
    real(dp), public, allocatable, save :: eigs(:)

    ! eigenvector of atomic many particle Hamiltonian
    complex(dp), public, allocatable, save :: eigv(:, :)

    ! impurity energy level
    complex(dp), public, allocatable, save :: cemat(:, :)

    ! atomic many particle Hamiltonian matrix (u + soc)
    complex(dp), public, allocatable, save :: hmat(:, :)

    ! transformation matrix from orginal basis to natural basis
    complex(dp), public, allocatable, save :: amtrx(:, :)

    ! spin orbit coupling matrix in {px, sz} single particle basis
    complex(dp), public, allocatable, save :: somat(:, :)

    ! auxiliary complex(dp) matrix for temperary use
    complex(dp), public, allocatable, save :: zauxs(:, :)

    ! coefficents matrix for generalized interaction U in real-orbital basis
    complex(dp), public, allocatable, save :: cumat(:, :, :, :)

    ! coefficents matrix for generalized interaction U in natural basis
    complex(dp), public, allocatable, save :: cumat_t(:, :, :, :)

end module mod_glob_hmtrx

module mod_glob_subspaces
    use control, only: norbs
    use mod_subspace
    implicit none

    ! number of subspaces of $H_{\text{loc}}$
    integer, public, save :: nsubs
    ! data structure for all of the subspaces
    type(type_subspace), public, allocatable, save :: subspaces(:)
    ! the index from one subspace to another subspace 
    ! for create operator $d_{j^2,j_z}^{\dagger}\Ket{N,J_z} = \Ket{N+1,J_z+j_z}$
    integer, public, allocatable, save :: c_towhich(:,:)
    ! if truncated by $N$
    integer, public, allocatable, save :: c_towhich_trunk(:,:)
    ! for destory operator $d_{j^2,j_z}\Ket{N,J_z} = \Ket{N-1,J_z-j_z}$
    integer, public, allocatable, save :: d_towhich(:,:)
    ! if truncated by $N$
    integer, public, allocatable, save :: d_towhich_trunk(:,:)

    contains

    subroutine alloc_glob_subspaces()
        implicit none

        allocate(subspaces(nsubs)) 
        allocate(c_towhich(nsubs, norbs)) 
        allocate(c_towhich_trunk(nsubs, norbs)) 
        allocate(d_towhich(nsubs, norbs)) 
        allocate(d_towhich_trunk(nsubs, norbs)) 

        return
    end subroutine alloc_glob_subspaces

    subroutine dealloc_glob_subspaces()
        implicit none
        
        integer :: i
        ! we will deallocate memory for pointers in type_subsapce 
        ! before deallocating subspaces to avoid memory leak 
   
        do i=1, nsubs
            call dealloc_one_subspace(subspaces(i)) 
        enddo 

        if (allocated(subspaces))  deallocate(subspaces) 
        if (allocated(c_towhich))  deallocate(c_towhich) 
        if (allocated(c_towhich_trunk))  deallocate(c_towhich_trunk) 
        if (allocated(d_towhich))  deallocate(d_towhich) 
        if (allocated(d_towhich_trunk))  deallocate(d_towhich_trunk) 
       
        return
    end subroutine dealloc_glob_subspaces

end module mod_glob_subspaces

module mod_glob_fmat
    use constants, only: dp
    use control, only: norbs
    use mod_glob_subspaces
    use mod_fmat
    implicit none

    type(type_fmat), public, save, allocatable :: c_fmat(:,:)

    contains

    subroutine alloc_glob_fmat()
        implicit none

        allocate(c_fmat(nsubs,norbs))

        return
    end subroutine alloc_glob_fmat

    subroutine dealloc_glob_fmat()
        implicit none

        ! local variables
        ! loop index
        integer :: i,j 
        
        ! first, deallocate the memory of each c_fmat(i,j)
        do i=1, norbs
            do j=1, nsubs
                call dealloc_one_fmat(c_fmat(j, i))
            enddo
        enddo

        ! then, deallocate c_fmat
        if (allocated(c_fmat)) deallocate(c_fmat)

        return
    end subroutine dealloc_glob_fmat

end module mod_glob_fmat
!-------------------------------------------------------------------------!
!>>> memory managment subroutines (allocate and deallocate memory)
!-------------------------------------------------------------------------!
module mod_global
    use constants
    use control
    use mod_glob_basis
    use mod_glob_hmtrx
    use mod_glob_subspaces
    use mod_glob_fmat
    implicit none

    ! status flag
    integer, private :: istat

    contains

    !>>> allocate memory atomic_basis
    subroutine atomic_allocate_memory_basis()

        implicit none

        allocate(nstat(0:norbs), stat=istat); nstat = 0
        allocate(basis(1:ncfgs), stat=istat); basis = 0
        allocate(good(1:norbs), stat=istat); good = 0
        allocate(invsn(0:2**norbs-1), stat=istat); invsn = 0
        allocate(invcd(1:norbs, 1:ncfgs), stat=istat); invcd = 0

        if (istat /= 0) then
            stop "error happened in atomic_allocate_memory_basis"
        endif ! back if (istat /= 0) block

        return
    end subroutine atomic_allocate_memory_basis

    !>>> allocate memory atomic_Hmtrx
    subroutine atomic_allocate_memory_hmtrx()
        implicit none

        allocate(eloc(norbs), stat=istat); eloc = 0.0D0
        allocate(eigs(ncfgs), stat=istat); eigs = 0.0D0
        allocate(eigv(ncfgs, ncfgs), stat=istat); eigv = dcmplx(0.0D0, 0.0D0)
        allocate(amtrx(norbs, norbs), stat=istat); amtrx = dcmplx(0.0D0, 0.0D0)
        allocate(somat(norbs, norbs), stat=istat); somat = dcmplx(0.0D0, 0.0D0)
        allocate(cemat(norbs, norbs), stat=istat); cemat = dcmplx(0.0D0, 0.0D0)
        allocate(hmat(ncfgs, ncfgs), stat=istat); hmat = dcmplx(0.0D0, 0.0D0)
        allocate(cumat(norbs, norbs, norbs, norbs), stat=istat); cumat = dcmplx(0.0D0, 0.0D0)
        allocate(cumat_t(norbs, norbs, norbs, norbs), stat=istat); cumat_t = dcmplx(0.0D0, 0.0D0)

        if (istat /= 0) then
            stop "error happened in atomic_allocate_memory_hmtrx"
        endif ! back if (istat /= 0) block

        return
    end subroutine atomic_allocate_memory_hmtrx

    !>>> deallocate memory atomic_basis
    subroutine atomic_deallocate_memory_basis()
        implicit none

        if (allocated(nstat)) deallocate(nstat)
        if (allocated(invsn)) deallocate(invsn)
        if (allocated(invcd)) deallocate(invcd)
        if (allocated(basis)) deallocate(basis)
        if (allocated(good))  deallocate(good)

        return
    end subroutine atomic_deallocate_memory_basis

    !>>> deallocate memory atomic_hmtrx
    subroutine atomic_deallocate_memory_hmtrx()
        implicit none

        if (allocated(cemat)) deallocate(cemat)
        if (allocated(hmat)) deallocate(hmat)

        if (allocated(eigs)) deallocate(eigs)
        if (allocated(eigv)) deallocate(eigv)

        if (allocated(cumat)) deallocate(cumat)
        if (allocated(somat)) deallocate(somat)
        if (allocated(cumat_t)) deallocate(cumat_t)
        
        return
    end subroutine atomic_deallocate_memory_hmtrx

    subroutine atomic_make_good()
        use mod_mkgood
        implicit none
        
        select case(norbs)
        case(6)
            call make_good_3band(good)
        case(10)
            call make_good_5band(good)
        case(14)
            call make_good_7band(good)           
        end select
         
        return
    end subroutine atomic_make_good

end module mod_global
