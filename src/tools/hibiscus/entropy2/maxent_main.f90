!!!=========+=========+=========+=========+=========+=========+=========+!
!!! ENTROPY2 @ iQIST                                                     !
!!!                                                                      !
!!! author  : Yilin Wang (at IOP/CAS)                                    !
!!! version : v2014.10.13T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment :                                                            !
!!!=========+=========+=========+=========+=========+=========+=========+!

  program main

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
