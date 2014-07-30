!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : stack
!!!           stack@istack_create
!!!           stack@istack_clean
!!!           stack@istack_destroy
!!!           stack@istack_copyer
!!!           stack@istack_setter
!!!           stack@istack_getter
!!!           stack@istack_push
!!!           stack@istack_pop
!!!           stack@istack_display
!!!           stack@istack_gettop
!!!           stack@istack_getrest
!!!           stack@istack_getsize
!!!           stack@istack_isfull
!!!           stack@istack_isempty
!!! source  : m_stack.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/14/2009 by li huang
!!!           07/09/2014 by li huang
!!!           07/30/2014 by li huang
!!! purpose : the purpose of this module is to define a stack-type (LIFO)
!!!           data structure in fortran version
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!

  module stack
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

! stack size limit, default value
     integer, private, parameter :: limit = 1024

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! define integer type stack
     type istack

! top position of stack
         integer :: top

! size of allocatable array
         integer :: nsize

! allocatable array, which is used to store elements in stack
         integer, allocatable :: item(:)

     end type istack

! define generic type stack
     type gstack

! top position of stack
         integer :: top

! size of allocatable array
         integer :: nsize

! allocatable array, which is used to store elements in stack
! note: it is an unlimited polymorphic object
         class(*), allocatable :: item(:)

     end type gstack

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: istack_create
     public :: istack_clean
     public :: istack_destroy

     public :: gstack_create
     public :: gstack_clean
     public :: gstack_destroy

     public :: istack_copyer
     public :: istack_setter
     public :: istack_getter

     public :: gstack_copyer
     public :: gstack_setter
     public :: gstack_getter

     public :: istack_push
     public :: istack_pop

     public :: gstack_push
     public :: gstack_pop

     public :: istack_display

     public :: gstack_display

     public :: istack_gettop
     public :: istack_getrest
     public :: istack_getsize

     public :: gstack_gettop
     public :: gstack_getrest
     public :: gstack_getsize

     public :: istack_isfull
     public :: istack_isempty

     public :: gstack_isfull
     public :: gstack_isempty

  contains ! encapsulated functionality

!!>>> create and initialize a integer type stack
  type (istack) &
  function istack_create(n) result (s)
     implicit none

! external arguments
! size of stack
     integer, optional, intent(in) :: n

! local variables
! status flag
     integer :: istat

! determine the capacity of stack
     if ( present (n) ) then
         s%nsize = n
     else
         s%nsize = limit
     endif

! setup the top position
     s%top = 0

! allocate memory for item array
     allocate(s%item(s%nsize), stat=istat)

     return
  end function istack_create

!!>>> reset the integer type stack, clean all its elements
  subroutine istack_clean(s)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! reset top position
     s%top = 0

     return
  end subroutine istack_clean

!!>>> destroy and finalize a integer type stack
  subroutine istack_destroy(s)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! deallocate memory
     if ( allocated(s%item) ) deallocate(s%item)

! reset top position
     s%top = 0

     return
  end subroutine istack_destroy

!!>>> copy an istack object to another
  subroutine istack_copyer(sa, sb)
     implicit none

! external arguments
! integer type stack, input
     type (istack), intent(in) :: sa

! integer type stack, output
     type (istack), intent(inout) :: sb

! check nsize at first
     if ( sa%nsize /= sb%nsize ) then
         write(mystd,'(a)') 'istack: the sizes of two stacks are not equal'
         STOP
     endif

! sync the data
     sb%top = sa%top
     sb%item = sa%item

     return
  end subroutine istack_copyer

!!>>> update the item's value of istack at special position
  subroutine istack_setter(s, item, pos)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! elements to be setted
     integer, intent(in) :: item

! position of the element to be updated
     integer, intent(in) :: pos

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         s%item(pos) = item
     endif

     return
  end subroutine istack_setter

!!>>> return the item's value of istack at special position
  integer &
  function istack_getter(s, pos) result (item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

! position of the element
     integer, intent(in) :: pos

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         item = s%item(pos)
     endif

     return
  end function istack_getter

!!>>> push item on top of stack
  subroutine istack_push(s, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! elements to be pushed in the stack
     integer, intent(in) :: item

     if ( s%top == s%nsize ) then
         write(mystd,'(a)') 'istack: the stack is full, can not push any item on it'
         STOP
     else
         s%top = s%top + 1
         s%item(s%top) = item
     endif

     return
  end subroutine istack_push

!!>>> pop off item from the top of stack
  integer &
  function istack_pop(s) result (item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'istack: the stack is empty, can not pop off any item from it'
         STOP
     else
         item = s%item(s%top)
         s%top = s%top - 1
     endif

     return
  end function istack_pop

!!>>> istack_display: display the top item in the stack without pop it off
  subroutine istack_display(s, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

! the top item in the stack
     integer, intent(out)      :: item

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'istack: the stack is empty, can not return the top item of it'
         STOP
     else
         item = s%item(s%top)
     endif ! back if ( s%top == 0 ) block

     return
  end subroutine istack_display

!!>>> gstack_display: display the top item in the stack without pop it off
  subroutine gstack_display(s, item)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

! the top item in the stack
     class(*), intent(out)     :: item

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'gstack: the stack is empty, can not return the top item of it'
         STOP
     else
         select type (v => s%item)
             type is (integer)
                 select type (item)
                     type is (integer)
                         item = v(s%top)
                 end select
             type is (logical)
                 select type (item)
                     type is (logical)
                         item = v(s%top)
                 end select
             type is (real(dp))
                 select type (item)
                     type is (real(dp))
                         item = v(s%top)
                 end select
             type is (complex(dp))
                 select type (item)
                     type is (complex(dp))
                         item = v(s%top)
                 end select
             type is (character(len = *))
         end select
     endif ! back if ( s%top == 0 ) block

     return
  end subroutine gstack_display

!!>>> istack_gettop: return the top position of the stack, i.e, the number
!!>>> of items stored in the stack currently
  integer &
  function istack_gettop(s) result (t)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     t = s%top

     return
  end function istack_gettop

!!>>> gstack_gettop: return the top position of the stack, i.e, the number
!!>>> of items stored in the stack currently
  integer &
  function gstack_gettop(s) result (t)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

     t = s%top

     return
  end function gstack_gettop

!!>>> istack_getrest: return the number of empty sites of the stack
  integer &
  function istack_getrest(s) result (r)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     r = s%nsize - s%top

     return
  end function istack_getrest

!!>>> gstack_getrest: return the number of empty sites of the stack
  integer &
  function gstack_getrest(s) result (r)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

     r = s%nsize - s%top

     return
  end function gstack_getrest

!!>>> istack_getsize: return the actual capacity of the stack
  integer &
  function istack_getsize(s) result (n)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     n = s%nsize

     return
  end function istack_getsize

!!>>> gstack_getsize: return the actual capacity of the stack
  integer &
  function gstack_getsize(s) result (n)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

     n = s%nsize

     return
  end function gstack_getsize

!!>>> istack_isfull: check whether the stack is full of items
  logical &
  function istack_isfull(s) result (b)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     b = ( s%top == s%nsize )

     return
  end function istack_isfull

!!>>> gstack_isfull: check whether the stack is full of items
  logical &
  function gstack_isfull(s) result (b)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

     b = ( s%top == s%nsize )

     return
  end function gstack_isfull

!!>>> istack_isempty: check whether the stack is empty
  logical &
  function istack_isempty(s) result (b)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     b = ( s%top == 0 )

     return
  end function istack_isempty

!!>>> gstack_isempty: check whether the stack is empty
  logical &
  function gstack_isempty(s) result (b)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

     b = ( s%top == 0 )

     return
  end function gstack_isempty

  end module stack






  module AAA
     type nstack
         class(*), allocatable :: item(:)
     end type nstack

  contains

  subroutine nstack_init(s, n)
     implicit none

     type (nstack) :: s
     integer :: n, i

     allocate(s%item(n), source=0)
     select type (v => s%item)
         type is (integer)
             do i = 1,n
                 v(i) = i
             enddo
             print *, v
     end select
  end subroutine nstack_init

  subroutine nstack_pop(s, n, item)
     type (nstack) :: s
     integer :: n
     class(*), intent(out) :: item

     select type (item)
         type is (integer)
             select type (p => s%item)
                 type is (integer)
                     item = p(n)
             end select
     end select
  end subroutine nstack_pop

  end module AAA

  program test
     use AAA

     implicit none

     type (nstack) :: s
     integer :: i = 2
     call nstack_init(s, 10)
     call nstack_pop(s, 5, i)
     print *, i 
  end program test
