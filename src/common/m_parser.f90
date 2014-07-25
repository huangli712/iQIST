!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : parser
!!!           parser@p_create
!!!           parser@p_destroy
!!!           parser@p_parse
!!!           parser@p_get
!!!           parser@p_get_vec
!!! source  : m_parser.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/26/2014 by li huang
!!! purpose : this purpose of this module is to implement a generic and
!!!           flexible config/input file reader and analyzer.
!!! status  : unstable
!!! comment : this module depends on constants and linkedlist modules
!!!-----------------------------------------------------------------------

  module parser
     use constants
     use linkedlist 

     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

! data structure for the items in config/input files: key-value pair
     type data_t
         logical             :: is_valid
         character(len = 32) :: str_key
         character(len = 32) :: str_value
     end type data_t

! pointer for the data_t structure
     type (data_t), pointer  :: data_ptr => null()

! pointer for the list_t structure
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
     public  :: p_get_vec

  contains ! encapsulated functionality

!!>>> p_create: init the linked list structure and insert a fake element
  subroutine p_create()
     implicit none

! allocate memory for data_ptr
     allocate(data_ptr)

! setup a fake element
     data_ptr%is_valid = .false.
     data_ptr%str_key = "i_am_key"
     data_ptr%str_value = "i_am_value"

! init the linked list structure and then insert the fake element
     call list_init(list_ptr, transfer(data_ptr, list_d))

     return
  end subroutine p_create

!!>>> p_destroy: free the linked list data structure
  subroutine p_destroy()
     implicit none

     call list_free(list_ptr)

     return
  end subroutine p_destroy

!!>>> p_parse: parse the input/config file, and store the key-value pairs
!!>>> to the linked list data structure
  subroutine p_parse(in_file)
     use, intrinsic :: iso_fortran_env, only : iostat_end

     implicit none

! external arguments
! filename for the input/config file
     character(len = *), intent(in) :: in_file

! local variables
! loop index
     integer :: p
     integer :: q

! file status flag
     integer :: istat

! original string for line reading
     character(len=100) :: string

! string representation for key of key-value pair
     character(len=100) :: str_key

! string representation for value of key-value pair
     character(len=100) :: str_value

! pointer for the linked list structure, used to access it
     type(list_t), pointer :: curr => null()

! open input/config file, here we do not judge whether the file exists
     open(mytmp, file = trim(in_file), form = 'formatted', status = 'unknown')

     FILE_READING: do
! read one line from the input file until we meet the end-of-file (EOF)
         read(mytmp, '(a100)', iostat = istat) string
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
             call s_str_compress(str_key)
             call s_str_lowcase(str_value)
             call s_str_compress(str_value)

! store the key-value pair in the linked list structure
             allocate(data_ptr)
             data_ptr%is_valid = .true.
             data_ptr%str_key = trim(str_key)
             data_ptr%str_value = trim(str_value)
             call list_insert(list_ptr, transfer(data_ptr, list_d))
         endif
     enddo FILE_READING ! over do loop

     close(mytmp)

     !!curr => list_ptr
     !!do p=1,list_count(list_ptr)-1
     !!    curr => list_next(curr)
     !!    data_ptr = transfer(list_get(curr), data_ptr)
     !!    print *, data_ptr%is_valid, data_ptr%str_key(1:3), data_ptr%str_value(1:3)
     !!enddo

     return
  end subroutine p_parse

  subroutine p_get(in_key, out_value)
     implicit none

     character(len = *), intent(in) :: in_key
     class(*), intent(inout) :: out_value

     character(len = 32) :: str_key
     character(len = 32) :: str_value
     type(list_t), pointer :: curr => null()

     integer :: p, q

     str_key = in_key
     call s_str_compress(str_key)
     call s_str_lowcase(str_key)

     curr => list_ptr
     do p=1,list_count(list_ptr)-1
         curr => list_next(curr)
         data_ptr  = transfer(list_get(curr), data_ptr)
         if ( trim(str_key) .eq. trim(data_ptr%str_key) ) then
             str_value = data_ptr%str_value
             call s_str_lowcase(str_value)
             call s_str_compress(str_value)
             EXIT
         endif
     enddo

     select type (out_value)
         type is (integer)
             !!print *, 'integer', str_value
             read (str_value,'(I10)') out_value
         type is (real(dp))
             !!print *, 'real(dp)', str_value
             read (str_value,'(F16.8)') out_value
         type is (logical)
             !!print *, 'logical', str_value
             read (str_value,'(L4)') out_value
         class default
             call s_print_error('p_get', 'unrecognize data type')
     end select

     return
  end subroutine p_get

  subroutine p_get_vec(in_key, out_value, nsize)
     implicit none

     integer, intent(in) :: nsize
     character(len = *), intent(in) :: in_key
     class(*), intent(inout) :: out_value(nsize)

     character(len = 32) :: str_key
     character(len = 32) :: str_value
     type(list_t), pointer :: curr => null()

     integer :: p, q, offset, int_aux
     real(dp) :: real_aux
     logical  :: bool_aux

     str_key = in_key
     call s_str_compress(str_key)
     call s_str_lowcase(str_key)

     curr => list_ptr
     do p=1,list_count(list_ptr)-1
         curr => list_next(curr)
         data_ptr  = transfer(list_get(curr), data_ptr)
         if ( trim(str_key) .eq. trim(data_ptr%str_key) ) then
             str_value = data_ptr%str_value
             call s_str_lowcase(str_value)
             call s_str_compress(str_value)
             EXIT
         endif
     enddo

     select type (out_value)
         type is (integer)
             print *, 'integer', str_value
             q=0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(I10)') int_aux
                 out_value(p) = int_aux
                 q = q + offset
             enddo
             read(str_value(q+1:), '(I10)') int_aux
             out_value(p) = int_aux
         type is (real(dp))
             print *, 'real(dp)', str_value
             q=0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(F16.8)') real_aux
                 out_value(p) = real_aux
                 q = q + offset
             enddo
             read(str_value(q+1:), '(F16.8)') real_aux
             out_value(p) = real_aux
         type is (logical)
             print *, 'logical', str_value
             q=0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(L4)') bool_aux
                 out_value(p) = bool_aux
                 q = q + offset
             enddo
             read(str_value(q+1:), '(L4)') bool_aux
             out_value(p) = bool_aux
         class default
             call s_print_error('p_get', 'unrecognize data type')
     end select

     return
  end subroutine p_get_vec

  end module parser
