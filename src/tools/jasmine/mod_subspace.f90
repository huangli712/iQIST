!====================================================================
! Module mod_subspace defines the type of good quantum numbers and
! the type of subspace.
!====================================================================
module mod_subspace
    use constants, only: dp
    implicit none

    ! good quantum numbers
    type :: type_gqn
        ! the total number of electrons
        integer :: ntot
        ! the total jz 
        integer :: jz  
    end type type_gqn

    type :: type_subspace 
        ! the index of this subspace
        integer :: indx 
        ! the dimension of this subspace
        integer :: ndim
        ! the good quantum number which labels this subspace
        type(type_gqn) :: mygqn
        ! the Fock basis index of this subspace
        integer, pointer :: mybasis(:)
        ! the Hamiltonian of this subspace
        complex(dp), pointer :: myham(:,:)
        ! the eigenvalues
        real(dp), pointer :: myeigval(:) 
        ! the eigenvectors, Hamiltonian should be real
        real(dp), pointer :: myeigvec(:,:) 
    end type type_subspace
    
    contains

    ! compare two good quantum numbers
    subroutine compare_gqn(obj1, obj2, equal)
        implicit none

        ! external variables
        type(type_gqn), intent(in)  :: obj1
        type(type_gqn), intent(in)  :: obj2
        logical, intent(out)        :: equal

        if ( obj1%ntot == obj2%ntot .and. obj1%jz == obj2%jz ) then
            equal = .true.
        else
            equal = .false. 
        endif

        return
    end subroutine compare_gqn

    ! allocate memory for onespace
    subroutine alloc_one_subspace(myspace)
        implicit none

        ! external variables
        type(type_subspace), intent(inout) :: myspace  

        allocate( myspace%mybasis(myspace%ndim) )
        allocate( myspace%myham(myspace%ndim, myspace%ndim) )
        allocate( myspace%myeigval(myspace%ndim) )
        allocate( myspace%myeigvec(myspace%ndim, myspace%ndim) )

        return
    end subroutine alloc_one_subspace

    ! deallocate memory for onespace
    subroutine dealloc_one_subspace(myspace)
        implicit none

        ! external variables
        type(type_subspace), intent(inout) :: myspace  

        if (associated(myspace%mybasis))    deallocate( myspace%mybasis )
        if (associated(myspace%myham))      deallocate( myspace%myham )
        if (associated(myspace%myeigval))   deallocate( myspace%myeigval )
        if (associated(myspace%myeigvec))   deallocate( myspace%myeigvec )

        return
    end subroutine dealloc_one_subspace

    ! nullify onespace
    subroutine nullify_one_subspace(myspace)
        implicit none

        ! external variables
        type(type_subspace), intent(inout) :: myspace  

        nullify( myspace%mybasis )
        nullify( myspace%myham )
        nullify( myspace%myeigval )
        nullify( myspace%myeigvec )

        return
    end subroutine nullify_one_subspace

end module mod_subspace


