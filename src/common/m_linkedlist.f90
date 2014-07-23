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

  module linkedlist
     implicit none

     private
     public :: list_node_t
     public :: list_data
     public :: list_init
     public :: list_free
     public :: list_insert
     public :: list_put
     public :: list_get
     public :: list_next

  ! A public variable used as a MOLD for transfer()
     integer, dimension(:), allocatable :: list_data

  ! Linked list node
     type :: list_node_t
         private
         integer, dimension(:), pointer :: data => null()
         type(list_node_t), pointer :: next => null()
     end type list_node_t

  contains

!!>>> Initialize a head node SELF and optionally store the provided DATA.
  subroutine list_init(self, data)
     type(list_node_t), pointer :: self
     integer, dimension(:), intent(in), optional :: data

     allocate(self)
     nullify(self%next)

     if (present(data)) then
         allocate(self%data(size(data)))
         self%data = data
     else
         nullify(self%data)
     end if
  end subroutine list_init

!!>>> Free the entire list and all data, beginning at SELF
  subroutine list_free(self)
     type(list_node_t), pointer :: self
     type(list_node_t), pointer :: current
     type(list_node_t), pointer :: next

     current => self
     do while (associated(current))
         next => current%next
         if (associated(current%data)) then
             deallocate(current%data)
             nullify(current%data)
         end if
         deallocate(current)
         nullify(current)
         current => next
     end do
  end subroutine list_free

!!>>> Insert a list node after SELF containing DATA (optional)
  subroutine list_insert(self, data)
     type(list_node_t), pointer :: self
     integer, dimension(:), intent(in), optional :: data
     type(list_node_t), pointer :: next

     allocate(next)

     if (present(data)) then
         allocate(next%data(size(data)))
         next%data = data
     else
         nullify(next%data)
     end if

     next%next => self%next
     self%next => next
  end subroutine list_insert

!!>>> Store the encoded DATA in list node SELF
  subroutine list_put(self, data)
     type(list_node_t), pointer :: self
     integer, dimension(:), intent(in) :: data

     if (associated(self%data)) then
         deallocate(self%data)
         nullify(self%data)
     end if
     self%data = data
  end subroutine list_put

!!>>> Return the DATA stored in the node SELF
  function list_get(self) result(data)
     type(list_node_t), pointer :: self
     integer, dimension(:), pointer :: data
     data => self%data
  end function list_get

!!>>> Return the next node after SELF
  function list_next(self)
     type(list_node_t), pointer :: self
     type(list_node_t), pointer :: list_next
     list_next => self%next
  end function list_next

end module linkedlist

program test_list
  use linkedlist
  implicit none

  ! Data is stored in data_t
  type :: data_t
     real :: x
  end type data_t

  ! A container for storing data_t pointers
  type :: data_ptr
     type(data_t), pointer :: p
  end type data_ptr


  type(list_node_t), pointer :: list => null()
  type(data_ptr) :: ptr

  ! Allocate a new data element
  allocate(ptr%p)
  ptr%p%x = 2.7183

  ! Initialize the list with the first data element
  call list_init(list, transfer(ptr, list_data))
  print *, 'Initializing list with data:', ptr%p

  ! Allocate a second data element
  allocate(ptr%p)
  ptr%p%x = 0.5772

  ! Insert the second into the list
  call list_insert(list, transfer(ptr, list_data))
  print *, 'Inserting node with data:', ptr%p

  ! Retrieve data from the second node and free memory
  ptr = transfer(list_get(list_next(list)), ptr)
  print *, 'Second node data:', ptr%p
  deallocate(ptr%p)

  ! Retrieve data from the head node and free memory
  ptr = transfer(list_get(list), ptr)
  print *, 'Head node data:', ptr%p
  deallocate(ptr%p)

  ! Free the list
  call list_free(list)
end program test_list
