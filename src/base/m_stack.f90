!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : stack
!!! source  : m_stack.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/14/2009 by li huang (created)
!!!           04/19/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a stack-type (LIFO)
!!!           data structure in fortran version.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! In this module, we implement two types of stack, istack and gstack. The
!! istack type was designed to deal with integer numbers only. However,
!! gstack is a generic type stack. More specifically, it supports the
!! following four data types:
!!     integer,
!!     logical,
!!     real(dp),
!!     complex(dp).
!! The usages, subroutine parameters for these two stack types are almost
!! identical. To implement gstack type, we generally use the unlimited
!! polymorphic features in fortran 2003/2008 standard. Noted that not all
!! fortran compilers can support these features. This module was tested
!! using intel fortran compiler only. We do not guarantee it can work/be
!! compiled correctly for using the other fortran compilers. So please
!! use it carefully. In the iqist project, so far we only use the istack
!! type. However, in the future, we will turn to the gstack type.
!!
!! Usage
!! =====
!!
!! 1. import stack support
!! -----------------------
!!
!! use stack
!!
!! 2. declare stack
!! ----------------
!!
!! type (istack) :: is
!! type (gstack) :: gs
!!
!! 3. create stack
!! ---------------
!!
!! call istack_create(is, 1024)
!! call gstack_create(gs, 1, 1024) ! create stack to support integer
!! call gstack_create(gs, .true., 1024) ! create stack to support logical
!! call gstack_create(gs, 1.0_dp, 1024) ! create stack to support real(dp)
!! call gstack_create(gs, (1.0_dp, 1.0_dp), 1024) ! create stack to support complex(dp)
!!
!! Note: In istack_create(), the second parameter is the capacity of the
!! stack. However, in gstack_create(), the second parameter means the data
!! type that gstack will manipulate, and the third parameter will be used
!! to determine the capacity. It is an optional parameter.
!!
!! 4. push element
!! ---------------
!!
!! call istack_push(is, 1)
!! call gstack_push(gs, 1)
!! call gstack_push(gs, .true.)
!! call gstack_push(gs, 2.0_dp)
!! call gstack_push(gs, (1.0_dp, 1.0_dp))
!!
!! 5. pop element
!! --------------
!!
!! call istack_pop(is, i) ! i is an integer
!! call gstack_pop(is, j) ! j can be integer, logical, real(dp), and complex(dp)
!!
!! 6. check the status of stack
!! ----------------------------
!!
!! print *, istack_isfull(is)
!! print *, istack_isempty(is)
!! print *, istack_getsize(is)
!! print *, istack_getrest(is)
!!
!! print *, gstack_isfull(gs)
!! print *, gstack_isempty(gs)
!! print *, gstack_getsize(gs)
!! print *, gstack_getrest(gs)
!!
!! The above three function calls will tell you whether the stack is full,
!! whether it is empty, and its capacity.
!!
!! 7. clean the stack
!! ------------------
!!
!! call istack_clean(is)
!! call gstack_clean(gs)
!!
!! Note: This operation will reset the top position of the stack, instead
!! of releasing the memory of it. So you can still use the stack after that.
!!
!! 8. destroy the stack
!! --------------------
!!
!! call istack_destroy(is)
!! call gstack_destroy(gs)
!!
!! Note: When the stack was destroyed, you can not use it any more.
!!
!!

  module stack
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp    = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

! stack size limit, default value
     integer, private, parameter :: limit = 1024

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

!!
!! @sub istack_create
!!
!! create and initialize an integer type stack
!!
  subroutine istack_create(s, n)
     implicit none

! external arguments
! size of stack
     integer, optional, intent(in) :: n

! integer type stack
     type (istack), intent(out)    :: s

! determine the capacity of stack
     if ( present (n) ) then
         s%nsize = n
     else
         s%nsize = limit
     endif ! back if ( present (n) ) block

! setup the top position
     s%top = 0

! allocate memory for item array
     allocate(s%item(s%nsize))

     return
  end subroutine istack_create

!!>>> gstack_create: create and initialize a generic type stack
  subroutine gstack_create(s, t, n)
     implicit none

! external arguments
! size of stack
     integer, optional, intent(in) :: n

! mold for the elements in the stack
     class(*), intent(in)          :: t

! generic type stack
     type (gstack), intent(out)    :: s

! determine the capacity of stack
     if ( present (n) ) then
         s%nsize = n
     else
         s%nsize = limit
     endif ! back if ( present (n) ) block

! setup the top position
     s%top = 0

! allocate memory for item array
     select type (t)
         type is (integer)
             allocate(s%item(s%nsize), source = 0)

         type is (logical)
             allocate(s%item(s%nsize), source = .true.)

         type is (real(dp))
             allocate(s%item(s%nsize), source = 0.0_dp)

         type is (complex(dp))
             allocate(s%item(s%nsize), source = (0.0_dp, 0.0_dp))
     end select

     return
  end subroutine gstack_create

!!>>> istack_clean: reset the integer type stack, clean all its elements
  subroutine istack_clean(s)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! reset top position
     s%top = 0

     return
  end subroutine istack_clean

!!>>> gstack_clean: reset the generic type stack, clean all its elements
  subroutine gstack_clean(s)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(inout) :: s

! reset top position
     s%top = 0

     return
  end subroutine gstack_clean

!!>>> istack_destroy: destroy and finalize an integer type stack
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

!!>>> gstack_destroy: destroy and finalize a generic type stack
  subroutine gstack_destroy(s)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(inout) :: s

! deallocate memory
     if ( allocated(s%item) ) deallocate(s%item)

! reset top position
     s%top = 0

     return
  end subroutine gstack_destroy

!!>>> istack_copyer: copy an istack object to another
  subroutine istack_copyer(sa, sb)
     implicit none

! external arguments
! integer type stack, input
     type (istack), intent(in)    :: sa

! integer type stack, output
     type (istack), intent(inout) :: sb

! check nsize at first
     if ( sa%nsize /= sb%nsize ) then
         write(mystd,'(a)') 'istack: the sizes of two stacks are not equal'
         STOP
     endif ! back if ( sa%nsize /= sb%nsize ) block

! sync the data
     sb%top = sa%top
     sb%item = sa%item

     return
  end subroutine istack_copyer

!!>>> gstack_copyer: copy an gstack object to another
  subroutine gstack_copyer(sa, sb)
     implicit none

! external arguments
! generic type stack, input
     type (gstack), intent(in)    :: sa

! generic type stack, output
     type (gstack), intent(inout) :: sb

! check nsize at first
     if ( sa%nsize /= sb%nsize ) then
         write(mystd,'(a)') 'gstack: the sizes of two stacks are not equal'
         STOP
     endif ! back if ( sa%nsize /= sb%nsize ) block

! sync the data
     sb%top = sa%top
     select type (v => sb%item)
         type is (integer)
             select type (u => sa%item)
                 type is (integer)
                     v = u
             end select

         type is (logical)
             select type (u => sa%item)
                 type is (logical)
                     v = u
             end select

         type is (real(dp))
             select type (u => sa%item)
                 type is (real(dp))
                     v = u
             end select

         type is (complex(dp))
             select type (u => sa%item)
                 type is (complex(dp))
                     v = u
             end select
     end select

     return
  end subroutine gstack_copyer

!!>>> istack_setter: update the item's value of stack at specified position
  subroutine istack_setter(s, pos, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! position of the element to be updated
     integer, intent(in)          :: pos

! elements to be setted
     integer, intent(in)          :: item

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         s%item(pos) = item
     endif ! back if ( pos < 1 .or. pos > s%nsize ) block

     return
  end subroutine istack_setter

!!>>> gstack_setter: update the item's value of stack at specified position
  subroutine gstack_setter(s, pos, item)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(inout) :: s

! position of the element to be updated
     integer, intent(in)          :: pos

! elements to be setted
     class(*), intent(in)         :: item

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'gstack: the position is not correct'
         STOP
     else
         select type (v => s%item)
             type is (integer)
                 select type (item)
                     type is (integer)
                         v(pos) = item
                 end select

             type is (logical)
                 select type (item)
                     type is (logical)
                         v(pos) = item
                 end select

             type is (real(dp))
                 select type (item)
                     type is (real(dp))
                         v(pos) = item
                 end select

             type is (complex(dp))
                 select type (item)
                     type is (complex(dp))
                         v(pos) = item
                 end select
         end select
     endif ! back if ( pos < 1 .or. pos > s%nsize ) block

     return
  end subroutine gstack_setter

!!>>> istack_getter: return the item's value of stack at specified position
  subroutine istack_getter(s, pos, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

! position of the item
     integer, intent(in)       :: pos

! the item's value
     integer, intent(out)      :: item

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         item = s%item(pos)
     endif ! back if ( pos < 1 .or. pos > s%nsize ) block

     return
  end subroutine istack_getter

!!>>> gstack_getter: return the item's value of stack at specified position
  subroutine gstack_getter(s, pos, item)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(in) :: s

! position of the item
     integer, intent(in)       :: pos

! the item's value
     class(*), intent(out)     :: item

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'gstack: the position is not correct'
         STOP
     else
         select type (v => s%item)
             type is (integer)
                 select type (item)
                     type is (integer)
                         item = v(pos)
                 end select

             type is (logical)
                 select type (item)
                     type is (logical)
                         item = v(pos)
                 end select

             type is (real(dp))
                 select type (item)
                     type is (real(dp))
                         item = v(pos)
                 end select

             type is (complex(dp))
                 select type (item)
                     type is (complex(dp))
                         item = v(pos)
                 end select
         end select
     endif ! back if ( pos < 1 .or. pos > s%nsize ) block

     return
  end subroutine gstack_getter

!!>>> istack_push: push item on top of stack
  subroutine istack_push(s, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! elements to be pushed in the stack
     integer, intent(in)          :: item

     if ( s%top == s%nsize ) then
         write(mystd,'(a)') 'istack: the stack is full, can not push any item on it'
         STOP
     else
         s%top = s%top + 1
         s%item(s%top) = item
     endif ! back if ( s%top == s%nsize ) block

     return
  end subroutine istack_push

!!>>> gstack_push: push item on top of stack
  subroutine gstack_push(s, item)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(inout) :: s

! elements to be pushed in the stack
     class(*), intent(in)         :: item

     if ( s%top == s%nsize ) then
         write(mystd,'(a)') 'gstack: the stack is full, can not push any item on it'
         STOP
     else
         s%top = s%top + 1
         select type (v => s%item)
             type is (integer)
                 select type (item)
                     type is (integer)
                         v(s%top) = item
                 end select

             type is (logical)
                 select type (item)
                     type is (logical)
                         v(s%top) = item
                 end select

             type is (real(dp))
                 select type (item)
                     type is (real(dp))
                         v(s%top) = item
                 end select

             type is (complex(dp))
                 select type (item)
                     type is (complex(dp))
                         v(s%top) = item
                 end select
         end select
     endif ! back if ( s%top == s%nsize ) block

     return
  end subroutine gstack_push

!!>>> istack_pop: pop off item from the top of stack
  subroutine istack_pop(s, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! the top item in the stack
     integer, intent(out)         :: item

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'istack: the stack is empty, can not pop off any item from it'
         STOP
     else
         item = s%item(s%top)
         s%top = s%top - 1
     endif ! back if ( s%top == 0 ) block

     return
  end subroutine istack_pop

!!>>> gstack_pop: pop off item from the top of stack
  subroutine gstack_pop(s, item)
     implicit none

! external arguments
! generic type stack
     type (gstack), intent(inout) :: s

! the top item in the stack
     class(*), intent(out)        :: item

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'gstack: the stack is empty, can not pop off any item from it'
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
         end select
         s%top = s%top - 1
     endif ! back if ( s%top == 0 ) block

     return
  end subroutine gstack_pop

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
