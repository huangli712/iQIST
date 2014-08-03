!---------------------------------------------------------------
! project : maxent
! program : maxent_print_header
!         : maxent_print_footer
!         : maxent_print_summary
! source  : maxent_print.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 06/11/2013 by yilin wang
! purpose : print some information 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

!==============================================================
! subroutine: maxent_print_header
! purpose   : print welcome information
!==============================================================
  subroutine maxent_print_header()
      use constants

      implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call maxent_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'MAXENT'
     write(mystd,'(2X,a)') '>>> An Analytic Continuation Program Based On Maximum Entropy Method'
     write(mystd,*)

     write(mystd,'(2X,a)') 'version: 2013.06.20T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'develop: by yilin wang, IOP'
     write(mystd,'(2X,a)') 'support: qhwyl2006@126.com'
     write(mystd,'(2X,a)') 'license: GPL2 and later versions'
     write(mystd,*)

     write(mystd,'(2X,a)') 'MAXENT >>> start running at '//date_time_string
     write(mystd,*)

     return
  end subroutine maxent_print_header

!===============================================================
! subroutine: maxent_print_footer
! purpose   : print goodbye information
!===============================================================
  subroutine maxent_print_footer()
      use constants
  
      implicit none
      character (len = 20) :: date_time_string

! used to record the time usage information
      real(dp) :: tot_time

! obtain time usage information
      call cpu_time(tot_time)

! obtain current date and time
      call maxent_time_builder(date_time_string)
      write(mystd,'(2X,a)') 'MAXENT >>> Calculation Completed.'

      write(mystd,'(2X,a,f10.2,a)') 'MAXENT >>> total time spent:', tot_time, 's'
      write(mystd,*)

      write(mystd,'(4X,a)')  "Suggested references for the acknowledgment of this MAXENT program"
      write(mystd,'(4X,a)')  "This maxent program is based on the algorithm proposed by M. Jarrell et al. and R. K. Bryan."
      write(mystd,'(4X,a)')  "Parts of the codes are modified from the original codes written by M. Jarrell."
      write(mystd,'(4X,a)')  "Original code link: http://www.phys.lsu.edu/~jarrell/CODES/MEM/"
      write(mystd, *)
      write(mystd,'(4X,a)')  "We suggest to cite the following references: "
      write(mystd,'(4X,a)')  "[1] Bayesian Inference and the Analytic Continuation of Imaginary-Time Quantum Monte Carlo Data, " 
      write(mystd,'(4X,a)')  "M. Jarrell, and J.E. Gubernatis, Physics Reports Vol. 269 #3, pp133-195, (1996)"
      write(mystd,'(4X,a)')  "[2] Maximum entropy analysis of oversampled data problems,"
      write(mystd,'(4X,a)')  "R. K. Bryan, Eur. Biophys. J. 18. 165 (1990)"
      write(mystd,*)

      write(mystd,'(2X,a)') 'MAXENT >>> I am tired and want to go to bed. Bye!'
      write(mystd,'(2X,a)') 'MAXENT >>> happy ending at '//date_time_string
      write(mystd,*)

      return
  end subroutine maxent_print_footer

!================================================================
! subroutine: maxent_print_summary
! purpose   : print control parameters
!================================================================
  subroutine maxent_print_summary()
      use constants
      use control

      implicit none

      write(mystd,"(2X,a)") "MAXENT >>> parameters lists"
      write(mystd,"(4X,a,i5)") "imode : ", imode
      write(mystd,"(4X,a,i5)") "icov  : ", icov
      write(mystd,"(4X,a,i5)") "ntime : ", ntime
      write(mystd,"(4X,a,i5)") "nwhf  : ", nwhf
      write(mystd,"(4X,a,i5)") "nw    : ", nw
      write(mystd,"(4X,a,i5)") "nalpha: ", nalpha
      write(mystd,"(4X,a,i5)") "nbins : ", nbins

      write(mystd,*)
      write(mystd,"(4X,a,f10.5)") "beta     : ", beta
      write(mystd,"(4X,a,f10.5)") "step     : ", step
      write(mystd,"(4X,a,f10.5)") "sigma    : ", sigma
      write(mystd,"(4X,a,f10.5)") "max_alpha: ", max_alpha
      write(mystd,"(4X,a,f10.5)") "min_alpha: ", min_alpha
      write(mystd,*)

      return
  end subroutine maxent_print_summary

!================================================================
! subroutine: maxent_print_error
! purpose   : print error information and stop the program
!================================================================
  subroutine maxent_print_error( sub, msg )
      use constants

      implicit none
     
! external variables
      character(len=*), intent(in) :: sub
      character(len=*), intent(in) :: msg

      write(mystd,"(2X,4a)") "fatal error occurred in ", sub, ": ", msg 

! stop the program
      STOP

      return
  end subroutine maxent_print_error

!================================================================
! subroutine: maxent_print_exception
! purpose   : print exception information and continue 
!================================================================
  subroutine maxent_print_exception( sub, msg )
      use constants

      implicit none

! external variables
      character(len=*), intent(in) :: sub
      character(len=*), intent(in) :: msg

      write(mystd,"(2X,4a)") "exception occurred in ", sub, ": ", msg 

! just ignore the exception, and continue the program
      CONTINUE

      return
  end subroutine maxent_print_exception

!=================================================================
! subroutine: maxent_print_advice
! purpose   : print advice when diagonalize kmat and amat error
!=================================================================
  subroutine maxent_print_advice(msg)
      use constants

      implicit none

! external variables
      character(len=*), intent(in) :: msg

      write(mystd,"(2X,a,a)") "MAXENT >>> Fatal Error Occured When Diagonalize ", msg
      write(mystd,"(4X,a)") "This error is caused by the bad eigenvalues of covariance matrix."
      write(mystd,"(4X,a)") "There are some solutions which may fix the problem:"
      write(mystd,"(6X,a)") "1. Do not diagonalize the covariance matrix, set icov=0."
      write(mystd,"(6X,a)") "2. Adjust the max_alpha and the min_alpha."
      write(mystd,"(6X,a)") "3. Change the default model and have a try agagin."
      write(mystd,"(6X,a)") "4. Set the frequency range smaller and have a try agagin."
      write(mystd,"(4X,a)") "WARNING: all the above solutions may fix your problem, however, "
      write(mystd,"(4X,a)") "we cannot make sure that your result is correct. "
      write(mystd,"(4X,a)") "FINAL SOLUTIONS: "
      write(mystd,"(6X,a)") "1. Check your Monte Carlo data carefully, make sure they are not oversampled."
      write(mystd,"(6X,a)") "2. Increase your number of data bins so that Nbins > 2L." 
      write(mystd,"(2X,a)") "MAXENT >>> We Will Stop The Program Here."

      STOP

      return
  end subroutine maxent_print_advice
