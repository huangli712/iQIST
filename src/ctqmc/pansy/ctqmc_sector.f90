!-------------------------------------------------------------------------
! project : pansy
! program : m_sector
! source  : mod_control.f90
! type    : modules
! authors : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014
! purpose : define data structure for good quantum number algorithm
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> data structure for good quantum number algorithm
module m_sector
    use constants,  only: dp, zero

    implicit none

    ! the fmat between any two sectors, it is just a matrix
    type :: t_fmat
        ! the dimension
        integer :: n, m
        ! the items of the matrix
        real(dp), pointer :: item(:,:)
    end type t_fmat

    ! one sector
    type :: t_sector 
        ! the dimension of this sector
        integer :: ndim
        ! total number of electrons n
        integer :: nelectron 
        ! number of fermion operators
        integer :: nops
        ! the start index of this sector
        integer :: istart
        ! the eigenvalues
        real(dp), pointer :: myeigval(:) 
        ! the next sector it points to when a fermion operator acts on this sector
        ! -1: outside of the Hilbert space, otherwise, it is the index of next sector
        ! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
        integer, pointer :: next_sector(:,:)
        ! for trunk of Hilbert space
        integer, pointer :: next_sector_trunk(:,:)
        ! the fmat between this sector and all other sectors
        ! if this sector doesn't point to some other sectors, the pointer is null
        ! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
        type(t_fmat), pointer :: myfmat(:,:)
        ! the final product matrices, which will be used to calculate the nmat
        real(dp), pointer :: final_product(:,:)
        ! matrices of occupancy operator c^{\dagger}c 
        real(dp), pointer :: occu(:,:,:)
        ! matrices of double occupancy operator c^{\dagger}cc^{\dagger}c
        real(dp), pointer :: double_occu(:,:,:,:)
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

        if ( associated(one_fmat%item) ) deallocate(one_fmat%item)

        return
    end subroutine dealloc_one_fmat

    !>>> nullify one sector
    subroutine nullify_one_sector(one_sector)
        implicit none

        ! external variables
        type(t_sector), intent(inout) :: one_sector

        nullify( one_sector%myeigval )
        nullify( one_sector%next_sector )
        nullify( one_sector%next_sector_trunk )
        nullify( one_sector%myfmat )

        return
    end subroutine nullify_one_sector

    !>>> allocate memory for one sector
    subroutine alloc_one_sector(one_sector)
        implicit none

        ! external variables
        type(t_sector), intent(inout) :: one_sector

        ! local variables
        integer :: i, j

        allocate(one_sector%myeigval(one_sector%ndim))
        allocate(one_sector%next_sector(one_sector%nops,0:1))
        allocate(one_sector%next_sector_trunk(one_sector%nops,0:1))
        allocate(one_sector%myfmat(one_sector%nops,0:1))
        allocate(final_product(one_sector%ndim, one_sector%ndim))
        allocate(occu(one_sector%ndim, one_sector%ndim, one_sector%nops))
        allocate(double_occu(one_sector%ndim, one_sector%ndim, one_sector%nops, one_sector%nops))

        ! init them
        one_sector%myeigval = zero
        one_sector%next_sector = 0
        one_sector%next_sector_trunk = 0
        final_product = zero
        occu = zero
        double_occu = zero

        ! init myfmat one by one
        do i=1, one_sector%nops 
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

        if (associated(one_sector%myeigval))           deallocate(one_sector%myeigval)
        if (associated(one_sector%next_sector))        deallocate(one_sector%next_sector)
        if (associated(one_sector%next_sector_trunk))  deallocate(one_sector%next_sector_trunk)
        if (associated(one_sector%final_product))      deallocate(one_sector%final_product)
        if (associated(one_sector%occu))               deallocate(one_sector%occu)
        if (associated(one_sector%double_occu))        deallocate(one_sector%double_occu)

        ! deallocate myfmat one by one
        do i=1, one_sector%nops
            do j=0,1
                call dealloc_one_fmat(one_sector%myfmat(i,j))
            enddo
        enddo 

        return
    end subroutine dealloc_one_sector
end module m_sector
