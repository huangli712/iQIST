!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : parser
!!!           parser@p_create
!!!           parser@p_destroy
!!!           parser@p_parse
!!!           parser@p_get
!!! source  : m_parser.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/24/2014 by li huang
!!! purpose : this purpose of this module is to implement a generic and
!!!           flexible config/input file reader and analyzer.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module parser
     use constants
     use linkedlist 

     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

     type data_t
         logical             :: is_valid
         character(len = 32) :: str_key
         character(len = 32) :: str_value
     end type data_t

     type (data_t), pointer  :: data_ptr => null()
     type (list_t), pointer  :: list_ptr => null()

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     private :: data_ptr
     private :: list_ptr

     public  :: p_create
     public  :: p_destroy
     public  :: p_parse
     public  :: p_get

  contains

!!>>> p_create:
  subroutine p_create()
     implicit none

     allocate(data_ptr)
     data_ptr%is_valid = .false.
     data_ptr%str_key = "i_am_key"
     data_ptr%str_value = "i_am_value"
     call list_init(list_ptr, transfer(data_ptr, list_d))

     return
  end subroutine p_create

!!>>> p_destroy:
  subroutine p_destroy()
     implicit none

     call list_free(list_ptr)

     return
  end subroutine p_destroy

!!>>> p_parse:
  subroutine p_parse(in_file)
     use, intrinsic :: iso_fortran_env

     implicit none

! external arguments
     character(len = *), intent(in) :: in_file

! local variables
     character(len=100) :: string
     character(len=100) :: str_key
     character(len=100) :: str_value
     integer :: istat
     integer :: p, q

     type(list_t), pointer :: curr => null()

     open(mytmp, file = trim(in_file), form = 'formatted', status = 'unknown')

     do
! read one line from the input file
         read(mytmp, '(a100)', iostat = istat) string
! meet the end-of-file, exit main loop
         if ( istat == iostat_end ) then
             EXIT
         else
! get rid of the empty and tab in the string
             call s_str_compress(string)

! is it an empty string?
             if ( len_trim(string) == 0   ) then
                 CYCLE
             endif

! is it a comment line starting with '#'?
             if ( index(string, '#') == 1 ) then
                 CYCLE
             endif

! is it a comment line starting with '!'?
             if ( index(string, '!') == 1 ) then
                 CYCLE
             endif

! get rid of the comment part starting with '#'
             p = index(string, '#')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif

! get rid of the comment part starting with '!'
             p = index(string, '!')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif

! extract the key and value pair from the input string
! note that the user can use both ':' and '=' symbols to separate the key
! and value pair. however, the ':' and '=' symbols can not appear at the
! same time for the same string.
             p = index(string, ':')
             q = index(string, '=')
             if ( p == 0 .and. q == 0 ) then
                 call s_print_error('p_parse', 'wrong file format for '//trim(in_file))
             endif
             if ( p >  0 .and. q >  0 ) then
                 call s_print_error('p_parse', 'wrong file format for '//trim(in_file))
             endif
             if ( p > 0 ) then
                 str_key = string(0:p-1)
                 str_value = string(p+1:len(string))
             endif
             if ( q > 0 ) then
                 str_key = string(0:q-1)
                 str_value = string(q+1:len(string))
             endif

! convert str_key and str_value to lowercase
             call s_str_lowcase(str_key)
             call s_str_lowcase(str_value)

! store the key-value pair in the linked list structure
             allocate(data_ptr)
             data_ptr%is_valid = .true.
             data_ptr%str_key = trim(str_key)
             data_ptr%str_value = trim(str_value)
             call list_insert(list_ptr, transfer(data_ptr, list_d))
         endif
     enddo ! over do loop

     close(mytmp)

!!     curr => list_ptr
!!     do p=1,list_count(list_ptr)-1
!!         curr => list_next(curr)
!!         data_ptr = transfer(list_get(curr), data_ptr)
!!         print *, data_ptr%is_valid, data_ptr%str_key, data_ptr%str_value
!!     enddo

     return
  end subroutine p_parse

  subroutine p_get(in_key, out_value)
     implicit none

     character(len = *), intent(in) :: in_key
     class(*), intent(inout) :: out_value

     character(len = 32) :: str_key
     character(len = 32) :: str_value
     character(len = 32) :: str_key_cmp
     type(list_t), pointer :: curr => null()

     integer :: p, q

     str_key = in_key
     call s_str_compress(str_key)
     call s_str_lowcase(str_key)
     !!print *, str_key

     curr => list_ptr
     do p=1,list_count(list_ptr)-1
         curr => list_next(curr)
         data_ptr  = transfer(list_get(curr), data_ptr)
         str_key_cmp = data_ptr%str_key
         call s_str_compress(str_key_cmp)
         if ( trim(str_key) .eq. trim(str_key_cmp) ) then
             str_value = data_ptr%str_value
             call s_str_compress(str_value)
             EXIT
         endif
     enddo

     select type (out_value)
         type is (integer)
             print *, 'integer'
             !!read(out_value, *) dsdf
             print *, out_value
         type is (real(dp))
             print *, 'real(dp)'
         type is (logical)
             print *, 'logical'
     end select

     return
  end subroutine p_get

  end module parser
