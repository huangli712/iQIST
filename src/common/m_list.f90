!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : list
!!! source  : m_list.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!! purpose : this purpose of this module is to implement a typical and
!!!           useful data structure --- linked list.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module list
     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

! self-defined data structure which will be stored in the linked list
     type T_data
     end type T_data

! node for the linked list, it contains a self-pointer
     type T_node
         type (T_node), pointer :: next
         type (T_data)          :: data
     end type T_node

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     private :: T_node

     public  :: T_data

     public  :: list_create
     public  :: list_destroy
     public  :: list_count
     public  :: list_next
     public  :: list_insert
     public  :: list_insert_head
     public  :: list_delete_element
     public  :: list_get_data
     public  :: list_set_data

  contains ! encapsulated functionality

!     ! Arguments:
!     list       Pointer to new linked list
!     data       The data for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use list_destroy first to
!     destroy up an old list.
!
!!>>> list_create: create and initialise a list
  subroutine list_create( list, data )
     implicit none

     type(LINKED_LIST), pointer  :: list
     type(LIST_DATA), intent(in) :: data

     allocate( list )
     list%next => null()
     list%data =  data

     return
  end subroutine list_create

! list_destroy --
!     Destroy an entire list
! Arguments:
!     list       Pointer to the list to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!
subroutine list_destroy( list )
    type(LINKED_LIST), pointer  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    current => list
    do while ( associated(current%next) )
        next => current%next
        deallocate( current )
        current => next
    enddo
end subroutine list_destroy

! list_count --
!     Count the number of items in the list
! Arguments:
!     list       Pointer to the list
!
integer function list_count( list )
    type(LINKED_LIST), pointer  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    if ( associated(list) ) then
        list_count = 1
        current => list
        do while ( associated(current%next) )
            current => current%next
            list_count = list_count + 1
        enddo
    else
        list_count = 0
    endif
end function list_count

! list_next
!     Return the next element (if any)
! Arguments:
!     elem       Element in the linked list
! Result:
!
function list_next( elem ) result(next)
    type(LINKED_LIST), pointer :: elem
    type(LINKED_LIST), pointer :: next

    next => elem%next

end function list_next

! list_insert
!     Insert a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!
subroutine list_insert( elem, data )
    type(LINKED_LIST), pointer  :: elem
    type(LIST_DATA), intent(in) :: data

    type(LINKED_LIST), pointer :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data
end subroutine list_insert

! list_insert_head
!     Insert a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!
subroutine list_insert_head( list, data )
    type(LINKED_LIST), pointer  :: list
    type(LIST_DATA), intent(in) :: data

    type(LINKED_LIST), pointer :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem
end subroutine list_insert_head

! list_delete_element
!     Delete an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be
!                removed
!
subroutine list_delete_element( list, elem )
    type(LINKED_LIST), pointer  :: list
    type(LINKED_LIST), pointer  :: elem

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: prev

    if ( associated(list,elem) ) then
        list => elem%next
        deallocate( elem )
    else
        current => list
        prev    => list
        do while ( associated(current) )
            if ( associated(current,elem) ) then
                prev%next => current%next
                deallocate( current ) ! Is also "elem"
                exit
            endif
            prev    => current
            current => current%next
        enddo
    endif
!    allocate(next)
!
!    next%next => elem%next
!    elem%next => next
!    next%data =  data
end subroutine list_delete_element

! list_get_data
!     Get the data stored with a list element
! Arguments:
!     elem       Element in the linked list
!
function list_get_data( elem ) result(data)
    type(LINKED_LIST), pointer :: elem

    type(LIST_DATA)            :: data

    data = elem%data
end function list_get_data

! list_put_data
!     Store new data with a list element
! Arguments:
!     elem       Element in the linked list
!     data       The data to be stored
!
subroutine list_put_data( elem, data )
    type(LINKED_LIST), pointer  :: elem
    type(LIST_DATA), intent(in) :: data

    elem%data = data
end subroutine list_put_data

  end module list
