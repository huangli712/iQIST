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

! open input/config file, here we do not judge whether the file exists
     open(mytmp, file = trim(in_file), form = 'formatted', status = 'unknown')

     FILE_READING: do ! start reading the file

! read one line from the input file until we meet the end-of-file (EOF)
! the line content is stored in string
         read(mytmp, '(a100)', iostat = istat) string
         if ( istat == iostat_end ) then
             EXIT
         else ! it is not the end
! get rid of the empty, tab, and null in the string
             call s_str_compress(string)

! is it an empty string? if yes, then we turn to the next line
             if ( len_trim(string) == 0   ) then
                 CYCLE
             endif ! back if ( len_trim(string) == 0   ) block

! is it a comment line starting with '#'? if yes, then we skip it
             if ( index(string, '#') == 1 ) then
                 CYCLE
             endif ! back if ( index(string, '#') == 1 ) block

! is it a comment line starting with '!'? if yes, then we skip it
             if ( index(string, '!') == 1 ) then
                 CYCLE
             endif ! back if ( index(string, '!') == 1 ) block

! get rid of the comment part starting with '#'
             p = index(string, '#')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif ! back if ( p > 0 ) block

! get rid of the comment part starting with '!'
             p = index(string, '!')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif ! back if ( p > 0 ) block

! extract the key and value pair from the input string
! note that the user can use both ':' and '=' symbols to separate the key
! and value pair. however, the ':' and '=' symbols can not appear at the
! same time for the same string.
             p = index(string, ':')
             q = index(string, '=')
             if ( p == 0 .and. q == 0 ) then
                 call s_print_error('p_parse', 'wrong file format for '//trim(in_file))
             endif ! back if ( p == 0 .and. q == 0 ) block
             if ( p >  0 .and. q >  0 ) then
                 call s_print_error('p_parse', 'wrong file format for '//trim(in_file))
             endif ! back if ( p >  0 .and. q >  0 ) block
             if ( p > 0 ) then
                 str_key = string(0:p-1)
                 str_value = string(p+1:len(string))
             endif ! back if ( p > 0 ) block
             if ( q > 0 ) then
                 str_key = string(0:q-1)
                 str_value = string(q+1:len(string))
             endif ! back if ( q > 0 ) block

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
         endif ! back if ( istat == iostat_end ) block

     enddo FILE_READING ! over do loop

! close the input/config file
     close(mytmp)

     return
  end subroutine p_parse

!!>>> p_get: retrieve the key-value pair from the linked list data structure
  subroutine p_get(in_key, out_value)
     implicit none

! external arguments
! string representation for the key of key-value pair
     character(len = *), intent(in) :: in_key

! polymorphic object for the value of key-value pair
     class(*), intent(inout)        :: out_value

! local variables
! loop index
     integer :: p

! string representation for the key
     character(len = 32) :: str_key

! string representation for the value
     character(len = 32) :: str_value

! pointer for the linked list data structure
     type(list_t), pointer :: curr => null()

! copy in_key to str_key and then postprocess it
     str_key = in_key
     call s_str_compress(str_key)
     call s_str_lowcase(str_key)

! visit the linked list and try to find out the required key-value pair
! whose key is the same with str_key
     curr => list_ptr
     do p=1,list_count(list_ptr)-1
! note that we skip the first element since it is invalid
         curr => list_next(curr)
         data_ptr  = transfer(list_get(curr), data_ptr)
! the required key-value pair is found, extract the value to str_value
         if ( trim(str_key) .eq. trim(data_ptr%str_key) ) then
             str_value = data_ptr%str_value
             call s_str_lowcase(str_value)
             call s_str_compress(str_value)
             EXIT
         endif ! back if block
     enddo ! over do loop
     curr => null()

! convert str_value to out_value, here we only support the following
! four cases: 1. integer; 2. logical; 3. real(dp); 4. character(len=*)
     select type (out_value)
         type is (integer)          ! for integer
             read (str_value,'(I10)')   out_value
         type is (logical)          ! for logical
             read (str_value,'(L4)')    out_value
         type is (real(dp))         ! for double precision number
             read (str_value,'(F16.8)') out_value
         type is (character(len=*)) ! for character
             out_value = str_value
         class default
             call s_print_error('p_get', 'unrecognize data type')
     end select

     return
  end subroutine p_get

!!>>> p_get_vec: retrieve the key-value pair from the linked list data structure
!!>>> here value is a array instead of single object
  subroutine p_get_vec(in_key, out_value, nsize)
     implicit none

! external arguments
! size of out_value
     integer, intent(in)            :: nsize

! string representation for the key of key-value pair
     character(len = *), intent(in) :: in_key

! polymorphic object for the value of key-value pair
     class(*), intent(inout)        :: out_value(nsize)

! local variables
! loop index, and used to tokenize the string
     integer  :: p
     integer  :: q
     integer  :: offset

! auxiliary variables used in the converting process 
     integer  :: int_aux
     logical  :: bool_aux
     real(dp) :: real_aux

! string representation for the key
     character(len = 32) :: str_key

! string representation for the value
     character(len = 32) :: str_value

! pointer for the linked list data structure
     type(list_t), pointer :: curr => null()

! copy in_key to str_key and then postprocess it
     str_key = in_key
     call s_str_compress(str_key)
     call s_str_lowcase(str_key)

! visit the linked list and try to find out the required key-value pair
! whose key is the same with str_key
     curr => list_ptr
     do p=1,list_count(list_ptr)-1
! note that we skip the first element since it is invalid
         curr => list_next(curr)
         data_ptr  = transfer(list_get(curr), data_ptr)
! the required key-value pair is found, extract the value to str_value
         if ( trim(str_key) .eq. trim(data_ptr%str_key) ) then
             str_value = data_ptr%str_value
             call s_str_lowcase(str_value)
             call s_str_compress(str_value)
             EXIT
         endif ! back if block
     enddo ! over do loop
     curr => null()

! convert str_value to out_value, here we only support the following
! four cases: 1. integer; 2. logical; 3. real(dp); 4. character(len=*)
     select type (out_value)
         type is (integer)          ! for integer
             q = 0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(I10)') int_aux
                 out_value(p) = int_aux
                 q = q + offset
             enddo ! over p={1,nsize-1} loop
             read(str_value(q+1:), '(I10)') int_aux
             out_value(nsize) = int_aux
         type is (logical)          ! for logical
             q = 0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(L4)') bool_aux
                 out_value(p) = bool_aux
                 q = q + offset
             enddo ! over p={1,nsize-1} loop
             read(str_value(q+1:), '(L4)') bool_aux
             out_value(nsize) = bool_aux
         type is (real(dp))         ! for double precision number
             q = 0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 read (str_value(q+1:q+offset-1), '(F16.8)') real_aux
                 out_value(p) = real_aux
                 q = q + offset
             enddo ! over p={1,nsize-1} loop
             read(str_value(q+1:), '(F16.8)') real_aux
             out_value(nsize) = real_aux
         type is (character(len=*)) ! for character
             q = 0
             do p=1,nsize-1
                 offset = index(str_value(q+1:), ',')
                 out_value(p) = str_value(q+1:q+offset-1)
                 q = q + offset
             enddo ! over p={1,nsize-1} loop
             out_value(nsize) = str_value(q+1:)
         class default
             call s_print_error('p_get_vec', 'unrecognize data type')
     end select

     return
  end subroutine p_get_vec

  end module parser
