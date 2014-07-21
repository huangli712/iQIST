!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : parser
!!!           parser@p_create
!!!           parser@p_destroy
!!!           parser@p_parse
!!!           parser@p_get_integer
!!!           parser@p_get_real
!!!           parser@p_get_logical
!!!           parser@p_get_character
!!! source  : m_parser.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!! purpose : this purpose of this module is to implement a generic and
!!!           flexible config/input file reader and analyzer.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module parser
     use constants
     use linkedlist 

     implicit none

     type (T_data) :: v_data
     type (T_node), pointer :: v_list 
     character(len = 32) :: v_file

     private :: v_data
     private :: v_list
     private :: v_file

     public  :: p_create
     public  :: p_destroy
     public  :: p_parse
     public  :: p_get_integer
     public  :: p_get_real
     public  :: p_get_logical
     public  :: p_get_character

  contains

  subroutine p_create(in_file)
     implicit none

! external arguments
     character(len=*), intent(in) :: in_file

     v_file = in_file
     v_data%str_key = "i_am_key"
     v_data%str_value = "i_am_value"
     call list_create(v_list, v_data)

     !! test the linkedlist
     !call list_navigator(v_list)
     !call list_insert_head(v_list, v_data)
     !call list_navigator(v_list)
     !call list_destroy(v_list)
     !call list_navigator(v_list)
 
     return
  end subroutine p_create

  subroutine p_destroy()
     implicit none
  
     return
  end subroutine p_destroy

  subroutine p_parse()
     use, intrinsic :: iso_fortran_env

     implicit none

     character(len=100) :: string
     character(len=100) :: str_key
     character(len=100) :: str_value
     character(len=100), external :: s_str_compress
     integer :: istat
     integer :: p, q

     open(mytmp, file = v_file, form = 'formatted', status = 'unknown')

     do
         read(mytmp, '(a100)', iostat = istat) string
         if ( istat == iostat_end ) then
             EXIT
         else
             string = s_str_compress(string)
             if ( len_trim(string) == 0   ) CYCLE
             if ( index(string, '#') == 1 ) CYCLE
             if ( index(string, '!') == 1 ) CYCLE
             p = index(string, '#')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif
             p = index(string, '!')
             if ( p > 0 ) then
                 string = string(0:p-1)
             endif

             p = index(string, ':')
             q = index(string, '=')
             if ( p == 0 .and. q == 0 ) STOP
             if ( p >  0 .and. q >  0 ) STOP
             if ( p > 0 ) then
                 str_key = string(0:p-1)
                 str_value = string(p+1:len(string))
             endif
             if ( q > 0 ) then
                 str_key = string(0:q-1)
                 str_value = string(q+1:len(string))
             endif
             v_data%str_key = trim(str_key)
             v_data%str_value = trim(str_value)
             call list_insert_head(v_list, v_data)
         endif
     enddo

     close(mytmp)

     call list_navigator(v_list)

     return
  end subroutine p_parse

  subroutine p_get_integer()
     implicit none
  
     return
  end subroutine p_get_integer

  subroutine p_get_real()
     implicit none
  
     return
  end subroutine p_get_real

  subroutine p_get_logical()
     implicit none
  
     return
  end subroutine p_get_logical

  subroutine p_get_character()
     implicit none
  
     return
  end subroutine p_get_character

  end module parser
