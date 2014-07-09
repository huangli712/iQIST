!---------------------------------------------------------------
! project : maxent
! program : maxent_main
! source  : maxent_main.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 06/11/2013 by yilin wang
! purpose : 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

  program main
      use constants
      use control
      use context

      implicit none

! local variables
! print headers
      call maxent_print_header()

! set control parameters
      call maxent_config()

! print summary of control parameters
      call maxent_print_summary()

! setup key variables
      call maxent_setup_array()

! initialize the program
      call maxent_databins_init()

! call the driver
      call maxent_driver()

! finish the task, deallocate memory
      call maxent_final_array() 

! print footer
      call maxent_print_footer()

  end
