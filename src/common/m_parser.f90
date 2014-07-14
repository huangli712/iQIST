!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : parser
!!!           p_create
!!!           p_destroy
!!!           p_parse
!!!           p_get_integer
!!!           p_get_real
!!!           p_get_logical
!!!           p_get_character
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
     use mlist 

     implicit none

     type (T_data) :: m_data
     type (T_node), pointer :: m_list 
     character(len = 36) :: m_file

     private :: m_data
     private :: m_list
     private :: m_file

     public :: p_create
     public :: p_destroy
     public :: p_parse
     public :: p_get_integer
     public :: p_get_real
     public :: p_get_logical
     public :: p_get_character

  contains

  subroutine p_create()
     implicit none
 
     m_list => null()
 
     return
  end subroutine p_create

  subroutine p_destroy()
     implicit none
  
     return
  end subroutine p_destroy

  subroutine p_parse()
     implicit none
  
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
