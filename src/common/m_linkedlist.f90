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
!!!           linkedlist@list_get
!!!           linkedlist@list_set
!!!           linkedlist@list_next
!!!           linkedlist@list_count
!!!           linkedlist@list_navigator
!!!           linkedlist@list_display
!!! source  : m_linkedlist.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/18/2014 by li huang
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
         logical             :: is_valid  ! used to judge whether it is fine
         character(len = 32) :: str_key   ! string for key
         character(len = 32) :: str_value ! string for value
     end type T_data

! node for the linked list, it contains a pointer pointing to next node,
! and a T_data structure which is used to store data
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
     public :: list_get
     public :: list_set
     public :: list_next
     public :: list_count
     public :: list_display

  contains ! encapsulated functionality

!!>>> list_create: create and initialise a list
  subroutine list_create(self, data)
     implicit none

! external arguments
! pointer to new linked list
     type(T_node), pointer    :: self

! the data for the first element
     type(T_data), intent(in) :: data

     allocate(self)
     self%next => null()
     self%data =  data

     return
  end subroutine list_create

!!>>> list_destroy: destroy an entire list
  subroutine list_destroy(self)
     implicit none

! external arguments
! pointer to the list to be destroyed
     type(T_node), pointer :: self

! local variables
     type(T_node), pointer :: curr
     type(T_node), pointer :: next

     curr => self
     do while ( associated(curr%next) )
         next => curr%next
         deallocate(curr)
         nullify(curr)
         curr => next
     enddo ! over do while loop

     return
  end subroutine list_destroy

!!>>> list_insert: insert a new element
  subroutine list_insert(self, data)
     implicit none

! external arguments
! element in the linked list after which the new element should be inserted
     type(T_node), pointer    :: self

! the data for the new element
     type(T_data), intent(in) :: data

! local variables
     type(T_node), pointer :: next

     allocate(next)
     next%data =  data
     next%next => self%next
     self%next => next

     return
  end subroutine list_insert

!!>>> list_insert_head: insert a new element before the first element
  subroutine list_insert_head(self, data)
     implicit none

! external arguments
! start of the list
     type(T_node), pointer    :: self

! the data for the new element
     type(T_data), intent(in) :: data

! local variables
     type(T_node), pointer :: elem

     allocate(elem)
     elem%data =  data
     elem%next => self
     self      => elem

     return
  end subroutine list_insert_head

!!>>> list_delete: delete an element from the list
  subroutine list_delete(self, elem)
     implicit none

! external arguments
! header of the list
     type(T_node), pointer :: self

! element in the linked list to be removed
     type(T_node), pointer :: elem

! local variables
     type(T_node), pointer :: curr
     type(T_node), pointer :: prev

     if ( associated(self,elem) ) then
         self => elem%next
         deallocate(elem)
     else
         curr => self
         prev => self
         do while ( associated(curr) )
             if ( associated(curr,elem) ) then
                 prev%next => curr%next
                 deallocate(curr) ! Is also "elem"
                 EXIT
             endif
             prev => curr
             curr => curr%next
         enddo ! over do while loop
     endif ! back if block

     return
  end subroutine list_delete

!!>>> list_delete_head: delete an element from the list
  subroutine list_delete_head(self)
     implicit none

! external arguments
! header of the list
     type(T_node), pointer :: self

! local variables
     type(T_node), pointer :: curr

! more than 1 node
     if ( associated(self%next) ) then
         curr => self
         self => self%next
         deallocate(curr)
! only 1 node
     else
         deallocate(self)
     endif ! back if block

     return
  end subroutine list_delete_head

!!>>> list_get_data: get the data stored with a list element
  function list_get(elem) result(data)
     implicit none

! external arguments
! element in the linked list
     type(T_node), pointer :: elem

! return value
     type(T_data)          :: data

     data = elem%data

     return
  end function list_get

!!>>> list_put_data: store new data with a list element
  subroutine list_set(elem, data)
     implicit none

! external arguments
! element in the linked list
     type(T_node), pointer    :: elem

! the data to be stored
     type(T_data), intent(in) :: data

     elem%data = data

     return
  end subroutine list_set

!!>>> list_next: return the next element (if any)
  function list_next(elem) result(next)
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
  function list_count(self)
     implicit none

! external arguments
! pointer to the list
     type(T_node), pointer :: self

! local variables
     type(T_node), pointer :: curr

     if ( associated(self) ) then
         list_count = 1
         curr => self
         do while ( associated(curr%next) )
             list_count = list_count + 1
             curr => curr%next
         enddo ! over do while loop
     else
         list_count = 0
     endif ! back if block

     return
  end function list_count

!!>>> list_display: display all of the nodes in the list one by one
  subroutine list_navigator(self)
     implicit none

! external arguments
! pointer to the list
     type(T_node), pointer :: self

! local variables
! pointer to the current node
     type(T_node), pointer :: curr

! node counter
     integer :: cntr

     if ( associated(self) ) then
         cntr = 1
         curr => self
         call list_display(cntr, curr)
         do while ( associated(curr%next) )
             cntr = cntr + 1
             curr => curr%next
             call list_display(cntr, curr)
         enddo ! over do while loop
     endif ! back if block

     return
  end subroutine list_navigator

!!>>> list_display: display the data contained in the node
  subroutine list_display(cntr, curr)
     implicit none

! external arguments
! pointer to the current node
     type(T_node), pointer  :: curr

! node counter
     integer, intent(in)    :: cntr

! local variables
! data stored in the current node
     type(T_data) :: data

     data = list_get(curr)
     write(*,*) cntr, data%is_valid, data%str_key, data%str_value

     return
  end subroutine list_display

  end module linkedlist
