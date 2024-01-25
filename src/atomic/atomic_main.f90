!!!=========+=========+=========+=========+=========+=========+=========+!
!!! iQIST @ JASMINE                                                      !
!!!                                                                      !
!!! An atomic eigenvalue problem solver which is used to generate input  !
!!! file (atom.cix) for the hybridization expansion version continuous   !
!!! time quantum Monte Carlo (CT-HYB) quantum impurity solver            !
!!!                                                                      !
!!! author  : Yilin Wang (University of Science and Technology of China) !
!!!           Li Huang (China Academy of Engineering Physics)            !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : the atomic solver is based on Dr. Liang Du's rambutan code !
!!!           any question, please contact with huangli@caep.cn          !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!========================================================================
  PROGRAM ATOMIC_MAIN !                                                <<<
!!========================================================================

     implicit none

!! [body

     ! print the runtime information (start)
     call atomic_print_header()

     ! setup the control parameters
     call atomic_setup_param()

     ! verify the control parameters
     call atomic_check_param()

     ! allocate memories for global arrays
     call atomic_alloc_array()

     ! print the control parameters to terminal
     call atomic_print_summary()

     ! call the dispatcher to launch different tasks
     call atomic_dispatcher()

     ! deallocate memories for global arrays
     call atomic_final_array()

     ! print the runtime information (end)
     call atomic_print_footer()

!! body]

!!========================================================================
  END PROGRAM ATOMIC_MAIN !                                            <<<
!!========================================================================
