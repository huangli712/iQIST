!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : linkedlist
!!!           linkedlist@list_d
!!!           linkedlist@list_t
!!!           linkedlist@list_init
!!!           linkedlist@list_free
!!!           linkedlist@list_insert
!!!           linkedlist@list_put
!!!           linkedlist@list_get
!!!           linkedlist@list_next
!!!           linkedlist@list_count
!!! source  : m_linkedlist.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/23/2014 by li huang
!!! purpose : this purpose of this module is to implement a typical and
!!!           useful data structure --- linked list. it is a generic
!!!           linked list, capable of storing arbitrary data.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module linkedlist
     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

! a public variable used as a mold for transfer() subroutine
     integer, dimension(:), allocatable :: list_d

! node for the linked list, it contains a pointer pointing to next node,
! and the integer pointer array is used to store data
     type list_t
         private
         integer, dimension(:), pointer :: data => null()
         type(list_t), pointer :: next => null()
     end type list_t

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: list_d
     public :: list_t

     public :: list_init
     public :: list_free

     public :: list_insert
     public :: list_put
     public :: list_get
     public :: list_next
     public :: list_count

  contains ! encapsulated functionality

!!>>> list_init: initialize a head node self and optionally store the provided data
  subroutine list_init(self, data)
     implicit none

! external arguments
! pointer to new linked list
     type(list_t), pointer :: self

! the data for the first element
     integer, dimension(:), intent(in), optional :: data

! allocate memory for linked list
     allocate(self)
     nullify(self%next)

! check whether we should make an empty node
     if ( present(data) ) then
         allocate( self%data( size(data) ) )
         self%data = data
     else
         nullify(self%data)
     endif ! back if block

     return
  end subroutine list_init

!!>>> list_free: free the entire list and all data, beginning at self
  subroutine list_free(self)
     implicit none

! external arguments
! pointer to the list to be destroyed
     type(list_t), pointer :: self

! local variables
! pointer to the current node
     type(list_t), pointer :: curr

! pointer to the next node
     type(list_t), pointer :: next

! go through the whole linked list
     curr => self
     do while ( associated(curr) )
         next => curr%next
! release memory for the internal data
         if ( associated(curr%data) ) then
             deallocate(curr%data)
             nullify(curr%data)
         endif
! release memory for the node itself
         deallocate(curr)
         nullify(curr)
         curr => next
     enddo ! over do while loop

     return
  end subroutine list_free

!!>>> list_insert: insert a node after self containing data (optional)
  subroutine list_insert(self, data)
     implicit none

! external arguments
! element in the linked list after which the new element should be inserted
     type(list_t), pointer :: self

! the data for the new element
     integer, dimension(:), intent(in), optional :: data

! local variables
! pointer to new node
     type(list_t), pointer :: next

! allocate memory for new node
     allocate(next)

! whether we should build an empty node
     if ( present(data) ) then
         allocate( next%data( size(data) ) )
         next%data = data
     else
         nullify(next%data)
     endif ! back if block

! update the linked list
     next%next => self%next
     self%next => next

     return
  end subroutine list_insert

!!>>> list_put: store the encoded data in list node self
  subroutine list_put(self, data)
     implicit none

! external arguments
! element in the linked list
     type(list_t), pointer :: self

! the data to be stored
     integer, dimension(:), intent(in) :: data

! release old memory at first
     if ( associated(self%data) ) then
         deallocate(self%data)
         nullify(self%data)
     endif ! back if block

! allocate new memory
     allocate( self%data( size(data) ) )

! save the data
     self%data = data

     return
  end subroutine list_put

!!>>> list_get: return the data stored in the node self
  function list_get(self) result(data)
     implicit none

! external arguments
! element in the linked list
     type(list_t), pointer :: self

! function value
     integer, dimension(:), pointer :: data

     data => self%data

     return
  end function list_get

!!>>> list_next: return the next node after self
  function list_next(self) result(next)
     implicit none

! external arguments
! pointer to the list
     type(list_t), pointer :: self

! function value
     type(list_t), pointer :: next

     next => self%next

     return
  end function list_next

!!>>> list_count: count the number of items in the list
  function list_count(self) result(counter)
     implicit none

! external arguments
! pointer to the list
     type(list_t), pointer :: self

! function value
     integer :: counter

! local variables
! pointer to current node
     type(list_t), pointer :: curr

     if ( associated(self) ) then
         counter = 1
         curr => self
         do while ( associated(curr%next) )
             counter = counter + 1
             curr => curr%next
         enddo ! over do while loop
     else
         counter = 0
     endif ! back if block

     return
  end function list_count

  end module linkedlist
