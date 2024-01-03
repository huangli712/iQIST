!!!=========+=========+=========+=========+=========+=========+=========+!
!!! iQIST @ JASMINE                                                      !
!!!                                                                      !
!!! An atomic eigenvalue problem solver which is used to generate input  !
!!! file (atom.cix) for the hybridization expansion version continuous   !
!!! time quantum Monte Carlo (CTQMC) quantum impurity solver             !
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

     ! print the running header
     call atomic_print_header()

     ! setup the parameters
     call atomic_setup_param()

     ! check validity of control parameters
     call atomic_check_param()

     ! allocate memory spaces
     call atomic_alloc_array()

     ! print the summary of control parameters
     call atomic_print_summary()

     ! call the drivers to perform different tasks
     call atomic_diapatcher()

     ! deallocate memory spaces
     call atomic_final_array()

     ! print footer
     call atomic_print_footer()

!! body]

!!========================================================================
  END PROGRAM ATOMIC_MAIN !                                            <<<
!!========================================================================
