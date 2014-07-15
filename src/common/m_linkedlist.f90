!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : linkedlist
!!!           linkedlist@T_node
!!!           linkedlist@T_data
!!!           linkedlist@list_create
!!!           linkedlist@list_destroy
!!!           linkedlist@list_insert
!!!           linkedlist@list_insert_head
!!!           linkedlist@list_delete
!!!           linkedlist@list_delete_head
!!!           linkedlist@list_get_data
!!!           linkedlist@list_set_data
!!!           linkedlist@list_next
!!!           linkedlist@list_count
!!! source  : m_linkedlist.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!! purpose : this purpose of this module is to implement a typical and
!!!           useful data structure --- linked list.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!!>>> VERY IMPORTANT
!! in order to use the linked list (linkedlist) data structure, the user
!! should define his/her own data type at first, namely, type T_data. Note
!! that there are no pointers within it.
!!
!!

  module linkedlist
     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

! self-defined data structure which will be stored in the linked list
     type T_data
         character(len = 32) :: str_key
         character(len = 32) :: str_value
     end type T_data

! node for the linked list, it contains a self-pointer
     type T_node
         type (T_node), pointer :: next
         type (T_data)          :: data
     end type T_node

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: T_node
     public :: T_data

     public :: list_create
     public :: list_destroy
     public :: list_insert
     public :: list_insert_head
     public :: list_delete
     public :: list_delete_head
     public :: list_get_data
     public :: list_set_data
     public :: list_next
     public :: list_count

  contains ! encapsulated functionality

!!>>> list_create: create and initialise a list
  subroutine list_create( list, data )
     implicit none

! external arguments
! pointer to new linked list
     type(T_node), pointer    :: list

! the data for the first element
     type(T_data), intent(in) :: data

     allocate( list )
     list%next => null()
     list%data =  data

     return
  end subroutine list_create

!!>>> list_destroy: destroy an entire list
  subroutine list_destroy( list )
     implicit none

! external arguments
! pointer to the list to be destroyed
     type(T_node), pointer  :: list

! local variables
     type(T_node), pointer  :: current
     type(T_node), pointer  :: next

     current => list
     do while ( associated(current%next) )
         next => current%next
         deallocate( current )
         current => next
     enddo ! over do while loop

     return
  end subroutine list_destroy

!!>>> list_insert: insert a new element
  subroutine list_insert( elem, data )
     implicit none

! external arguments
! element in the linked list after which the new element should be inserted
     type(T_node), pointer    :: elem

! the data for the new element
     type(T_data), intent(in) :: data

! local variables
     type(T_node), pointer :: next

     allocate(next)
     next%data =  data
     next%next => elem%next
     elem%next => next

     return
  end subroutine list_insert

!!>>> list_insert_head: insert a new element before the first element
  subroutine list_insert_head( list, data )
     implicit none

! external arguments
! start of the list
     type(T_node), pointer    :: list

! the data for the new element
     type(T_data), intent(in) :: data

! local variables
     type(T_node), pointer :: elem

     allocate(elem)
     elem%data =  data
     elem%next => list
     list      => elem

     return
  end subroutine list_insert_head

!!>>> list_delete: delete an element from the list
  subroutine list_delete( list, elem )
     implicit none

! external arguments
! header of the list
     type(T_node), pointer  :: list

! element in the linked list to be removed
     type(T_node), pointer  :: elem

! local variables
     type(T_node), pointer  :: current
     type(T_node), pointer  :: prev

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
                 EXIT
             endif
             prev    => current
             current => current%next
         enddo ! over do while loop
     endif ! back if block

     return
  end subroutine list_delete

!!>>> list_delete_head: delete an element from the list
  subroutine list_delete_head( list )
     implicit none

! external arguments
! header of the list
     type(T_node), pointer  :: list

! local variables
     type(T_node), pointer  :: current

! more than 1 node
     if ( associated(list%next) ) then
         current => list
         list => list%next
         deallocate( current )
! only 1 node
     else
         deallocate( list )
     endif ! back if block

     return
  end subroutine list_delete_head

!!>>> list_get_data: get the data stored with a list element
  function list_get_data( elem ) result(data)
     implicit none

! external arguments
! element in the linked list
     type(T_node), pointer :: elem

! return value
     type(T_data)          :: data

     data = elem%data

     return
  end function list_get_data

!!>>> list_put_data: store new data with a list element
  subroutine list_set_data( elem, data )
     implicit none

! external arguments
! element in the linked list
     type(T_node), pointer    :: elem

! the data to be stored
     type(T_data), intent(in) :: data

     elem%data = data

     return
  end subroutine list_set_data

!!>>> list_next: return the next element (if any)
  function list_next( elem ) result(next)
     implicit none

! external arguments
! element in the linked list
     type(T_node), pointer :: elem
     type(T_node), pointer :: next

     next => elem%next

     return
  end function list_next

!!>>> list_count: count the number of items in the list
  integer &
  function list_count( list )
     implicit none

! external arguments
! pointer to the list
     type(T_node), pointer :: list

! local variables
     type(T_node), pointer :: current

     if ( associated(list) ) then
         list_count = 1
         current => list
         do while ( associated(current%next) )
             current => current%next
             list_count = list_count + 1
         enddo ! over do while loop
     else
         list_count = 0
     endif ! back if block

     return
  end function list_count

  end module linkedlist
