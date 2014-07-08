!>>> module for defining good quantum numbers
module m_good_quantum

    ! good quantum number for N
    type :: t_good_n
        ! total number of electrons
        integer :: ntot  
    end type t_good_n

    ! good quantum number for N, Sz, PS
    type :: t_good_nszps
        ! total number of electrons
        integer :: ntot
        ! total Sz
        integer :: sz
        ! PS number
        integer :: ps
    end type t_good_nszps

    ! good quantum numbers for N, Jz
    type :: t_good_njz
        ! the total number of electrons
        integer :: ntot
        ! the total jz 
        integer :: jz  
    end type t_good_njz

    contains

    subroutine make_good_jz(good_jz)
        use control

        implicit none

        ! external variables
        integer, intent(out) :: good_jz(norbs)

        if (nband == 3) then
            ! j=1/2
            good_jz(1) = -1
            good_jz(2) =  1
            ! j=3/2
            good_jz(3) = -3
            good_jz(4) = -1
            good_jz(5) =  1
            good_jz(6) =  3
        elseif (nband == 5) then
            ! j=3/2
            good_jz(1) = -3
            good_jz(2) = -1
            good_jz(3) =  1
            good_jz(4) =  3
            ! j=5/2
            good_jz(5) = -5
            good_jz(6) = -3
            good_jz(7) = -1
            good_jz(8) =  1
            good_jz(9) =  3
            good_jz(10)=  5
        elseif (nband == 7) then
            ! j=5/2
            good_jz(1) = -5
            good_jz(2) = -3
            good_jz(3) = -1
            good_jz(4) =  1
            good_jz(5) =  3
            good_jz(6) =  5
            ! j=7/2
            good_jz(7) = -7
            good_jz(8) = -5
            good_jz(9) = -3
            good_jz(10)= -1
            good_jz(11)=  1
            good_jz(12)=  3
            good_jz(13)=  5
            good_jz(14)=  7
        endif

        return
    end subroutine make_good_jz

end module m_good_quantum

!>>> data structure of sector
module m_sector
    use constants, only: dp
    implicit none

    ! the fmat between any two subspaces, it is just a matrix
    type t_fmat
        ! the dimension
        integer :: n, m
        ! the items of the matrix
        real(dp), pointer :: item(:,:)
    end type t_fmat

    ! one subspace
    type :: t_sector 
        ! the dimension of this subspace
        integer :: ndim
        ! total number of electrons n
        integer :: nelectron 
        ! number of operators
        integer :: nops
        ! the start index of this sector
        integer :: istart
        ! the Fock basis index of this subspace
        integer, pointer :: mybasis(:)
        ! the Hamiltonian of this subspace
        complex(dp), pointer :: myham(:,:)
        ! the eigenvalues
        real(dp), pointer :: myeigval(:) 
        ! the eigenvectors, Hamiltonian should be real
        real(dp), pointer :: myeigvec(:,:) 
        ! the next sector it points to when a fermi operator acts on this sector
        ! 0 for outside of the space, otherwise, it is the index of sector
        ! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
        integer, pointer :: next_sector(:,:)
        ! the fmat between this sector and all other sectors
        ! if this sector doesn't point to some other sectors, the pointer is null
        ! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
        type(t_fmat), pointer :: myfmat(:,:)
    end type t_sector
    
    contains

    !>>> nullify one fmat
    subroutine nullify_one_fmat(one_fmat)
        implicit none

        ! external variables
        type(t_fmat), intent(inout) :: one_fmat

        nullify(one_fmat%item)

        return
    end subroutine nullify_one_fmat

    !>>> allocate one fmat
    subroutine alloc_one_fmat(one_fmat)
        implicit none

        ! external variables
        type(t_fmat), intent(inout) :: one_fmat

        allocate(one_fmat%item(one_fmat%n, one_fmat%m))

        ! initialize it
        one_fmat%item = zero

        return
    end subroutine alloc_one_fmat

    !>>> deallocate one fmat
    subroutine dealloc_one_fmat(one_fmat)
        implicit none

        ! external variables
        type(t_fmat), intent(inout) :: one_fmat

        if (associated(one_fmat%item)) deallocate(one_fmat%item)

        return
    end subroutine dealloc_one_fmat

    !>>> nullify one sector
    subroutine nullify_one_sector(one_sector)
        implicit none

        ! external variables
        type(t_sector), intent(inout) :: one_sector

        nullify( one_sector%mybasis )
        nullify( one_sector%myham )
        nullify( one_sector%myeigval )
        nullify( one_sectore%myeigvec )
        nullify( one_sectore%next_sector )
        nullify( one_sectore%myfmat )

        return
    end subroutine nullify_one_sector

    !>>> allocate memory for one sector
    subroutine alloc_one_sector(one_sector)
        implicit none

        ! external variables
        type(t_sector), intent(inout) :: one_sector

        ! local variables
        integer :: i, j

        allocate(one_sector%mybasis(one_sector%ndim)) 
        allocate(one_sector%myham(one_sector%ndim, one_sector%ndim)) 
        allocate(one_sector%myeigval(one_sector%ndim))
        allocate(one_sector%myeigvec(one_sector%ndim, one_sector%ndim)) 
        allocate(one_sector%next_sector(one_sector%nops,0:1))
        allocate(one_sector%myfmat(one_sector%nops,0:1))

        ! init them
        one_sector%mybasis = 0
        one_sector%myham = czero
        one_sector%myeigval = zero
        one_sector%myeigvec = zero
        one_sector%next_sector = 0

        ! init myfmat one by one
        do i=1, nops 
           do j=0, 1
               one_sector%myfmat(i,j)%n = 0
               one_sector%myfmat(i,j)%m = 0
               call nullify_one_fmat(one_sector%myfmat(i,j))
           enddo
        enddo

        return
    end subroutine alloc_one_sector

    ! deallocate memory for onespace
    subroutine dealloc_one_sector(one_sector)
        implicit none

        ! external variables
        type(t_sector), intent(inout) :: one_sector 

        ! local variables  
        integer :: i, j

        if (associated(one_sector%mybasis))      deallocate(one_sector%mybasis)
        if (associated(one_sector%myham))        deallocate(one_sector%myham)
        if (associated(one_sector%myeigval))     deallocate(one_sector%myeigval)
        if (associated(one_sector%myeigvec))     deallocate(one_sector%myeigvec)
        if (associated(one_sector%next_sector))  deallocate(one_sector%next_sector)

        ! deallocate myfmat one by one
        do i=1, nops
            do j=0,1
                call dealloc_one_fmat(one_sector%myfmat(i,j))
            enddo
        enddo 

        return
    end subroutine dealloc_one_sector
end module m_sector
