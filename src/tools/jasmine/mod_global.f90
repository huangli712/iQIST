!=========================================================================!
! project : jasmine
! program : context module
! history : Apr 26, 2011
! authors : xidai and duliang {duleung@gmail.com}
! purpose : important arrays defined for main program 
! comment :
!=========================================================================!

!>>> Fock basis related matrix, fullspace case
  module m_basis_fullspace
     use control
     implicit none

! DIMension of each SUBspace labeled by total electron number N
     integer, public, allocatable, save :: dim_sub_N(:)

! BINary form of a BASIS
     integer, public, allocatable, save :: bin_basis(:,:)

! DECimal form of BASIS
     integer, public, allocatable, save :: dec_basis(:)

! INDEX of BASIS, given their decimal number
     integer, public, allocatable, save :: index_basis(:)

     contains

! allocate memory for these matrices
     subroutine alloc_m_basis_fullspace()
        implicit none

! allocate them
        allocate(dim_sub_N(0:norbs))
        allocate(bin_basis(norbs, ncfgs))
        allocate(dec_basis(ncfgs))
        allocate(index_basis(ncfgs))

! init them
        dim_sub_N = 0
        bin_basis = 0
        dec_bsis = 0
        index_basis = 0

        return
     end subroutine alloc_m_basis_fullspace

! deallocate memory for these matrices
     subroutine dealloc_m_basis_fullspace()
        implicit none

! deallocate them
        if (allocated(dim_sub_N))   deallocate(dim_sub_N)
        if (allocated(bin_basis))   deallocate(bin_basis)
        if (allocated(dec_basis))   deallocate(dec_basis) 
        if (allocated(index_basis)) deallocate(index_basis) 

        return
     end subroutine dealloc_m_basis_fullspace

  end module m_basis_fullspace

!>>> Single Particle MATrix
  module m_sp_mat
     use constants, only dp
     use control
     implicit none

! Single Particle Crystal Field MATrix
     complex(dp), public, allocatable, save :: sp_cf_mat(:,:) 

! Single Particle Spin-Orbital Coupling MATrix
     complex(dp), public, allocatable, save :: sp_soc_mat(:,:)

! Single Particle Coulomb U MATrix
     complex(dp), public, allocatable, save :: sp_cu_mat(:,:,:,:)

     contains

! allocate memory for these matrix
     subroutine alloc_m_sp_mat()
        implicit none

! allocate them
        allocate(sp_cf_mat(norbs, norbs))
        allocate(sp_soc_mat(norbs, norbs)
        allocate(sp_cu_mat(norbs, norbs, norbs, norbs))

! init them
        sp_cf_mat  = czero
        sp_soc_mat = czero
        sp_cu_mat  = czero

        return
     end subroutine alloc_m_sp_mat

! deallocate memory for these matrix
     subroutine dealloc_m_sp_mat()
        implicit none

! allocate them
        if (allocated(sp_cf_mat))  deallocate(sp_cf_mat)
        if (allocated(sp_soc_mat)) deallocate(sp_soc_mat)
        if (allocated(sp_cu_mat))  deallocate(sp_cu_mat)

        return
     end subroutine alloc_m_sp_mat
 
  end module m_sp_mat

!>>> atomic Hamiltonian martix related variables
module m_hmat_fullspace
    use constants, only: dp
    implicit none

    ! Many Particle Hamiltonian MATrix (CF + SOC + CU)
    complex(dp), public, allocatable, save :: mp_hmat(:, :)

    ! Many Particle Hamiltonian MATrix's EIGen VALues
    real(dp),    public, allocatable, save :: mp_hmat_eigval(:)

    ! Many Particle Hamiltonian MATrix's EIGen VECtors
    complex(dp), public, allocatable, save :: mp_hmat_eigvec(:, :)

end module m_hmat_fullspace

module m_subspaces
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
