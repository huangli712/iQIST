module mod_fmat
    use constants, only: dp
    implicit none

    ! struct of F matrix $<y|f^{\dagger}|x>$
    type :: type_fmat
        ! which orbital
        integer :: iorb
        ! from x subspace to y subspace
        integer :: x
        integer :: y
        ! dimension of x and y subspaces
        integer :: ndimx
        integer :: ndimy
        ! matrix elements of this fmat
        real(dp), pointer :: elem(:,:)
    end type type_fmat

    ! utility subroutines for constructing fmat
    contains
    
    subroutine alloc_one_fmat(this_fmat)
        implicit none

        ! external variables
        type(type_fmat), intent(inout) :: this_fmat
        allocate(this_fmat%elem(this_fmat%ndimy, this_fmat%ndimx))

        return
    end subroutine alloc_one_fmat

    subroutine dealloc_one_fmat(this_fmat)
        implicit none

        ! external variables
        type(type_fmat), intent(inout) :: this_fmat

        if (associated(this_fmat%elem)) deallocate(this_fmat%elem)

        return
    end subroutine dealloc_one_fmat

    subroutine nullify_one_fmat(this_fmat)
        implicit none

        ! external variables
        type(type_fmat), intent(inout) :: this_fmat

        nullify(this_fmat%elem)

        return
    end subroutine nullify_one_fmat

end module mod_fmat
