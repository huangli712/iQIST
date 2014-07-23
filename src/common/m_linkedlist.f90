!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : linkedlist
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
     public :: list_t
     public :: list_d

     public :: list_init
     public :: list_free

     public :: list_insert
     public :: list_put
     public :: list_get
     public :: list_next
     public :: list_count

! a public variable used as a MOLD for transfer()
     integer, dimension(:), allocatable :: list_d

! linked list node
     type list_t
         private
         integer, dimension(:), pointer :: data => null()
         type(list_t), pointer :: next => null()
     end type list_t

  contains

!!>>> list_init: initialize a head node SELF and optionally store the provided DATA.
  subroutine list_init(self, data)
     implicit none

! external arguments
     type(list_t), pointer :: self
     integer, dimension(:), intent(in), optional :: data

     allocate(self)
     nullify(self%next)

     if ( present(data) ) then
         allocate( self%data( size(data) ) )
         self%data = data
     else
         nullify(self%data)
     endif

     return
  end subroutine list_init

!!>>> list_free: free the entire list and all data, beginning at SELF
  subroutine list_free(self)
     implicit none

! external arguments
     type(list_t), pointer :: self

! local variables
     type(list_t), pointer :: current
     type(list_t), pointer :: next

     current => self
     do while ( associated(current) )
         next => current%next
         if ( associated(current%data) ) then
             deallocate(current%data)
             nullify(current%data)
         endif
         deallocate(current)
         nullify(current)
         current => next
     enddo

     return
  end subroutine list_free

!!>>> list_insert: insert a list node after SELF containing DATA (optional)
  subroutine list_insert(self, data)
     implicit none

! external arguments
     type(list_t), pointer :: self
     integer, dimension(:), intent(in), optional :: data

! local variables
     type(list_t), pointer :: next

     allocate(next)

     if ( present(data) ) then
         allocate( next%data( size(data) ) )
         next%data = data
     else
         nullify(next%data)
     endif

     next%next => self%next
     self%next => next

     return
  end subroutine list_insert

!!>>> list_put: store the encoded DATA in list node SELF
  subroutine list_put(self, data)
     implicit none

! external arguments
     type(list_t), pointer :: self
     integer, dimension(:), intent(in) :: data

     if ( associated(self%data) ) then
         deallocate(self%data)
         nullify(self%data)
     endif
     self%data = data

     return
  end subroutine list_put

!!>>> list_get: return the DATA stored in the node SELF
  function list_get(self) result(data)
     implicit none

! external arguments
     type(list_t), pointer :: self
     integer, dimension(:), pointer :: data
     data => self%data

     return
  end function list_get

!!>>> list_next: return the next node after SELF
  function list_next(self)
     implicit none

! external arguments
     type(list_t), pointer :: self
     type(list_t), pointer :: list_next
     list_next => self%next

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

  end module linkedlist
