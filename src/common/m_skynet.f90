!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : skynet
!!! source  : m_skynet.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 04/30/2010 by li huang
!!!           07/09/2014 by li huang
!!! purpose : this module is used to overload/define the basic arithmetic
!!!           operations for very large exponent numbers
!!! status  : unstable
!!! comment : this module is not activated until now
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! In this module, we define such a data structure (sunu) which is used
!! to represent very large float number. We overload the basic arithmetic
!! operations such +/-/*/\/ for it. Due to the performance issue, we do
!! not use it now. But it may be useful in the other fields.
!!
!! Usage
!! =====
!!
!! 1. import skynet support
!! ------------------------
!!
!! use skynet
!!
!! 2. declare sunu data type
!! -------------------------
!!
!! type (sunu) :: A
!!
!! 3. assignment operations for sunu structure
!! -------------------------------------------
!!
!! type (sunu) :: A
!! real(dp) :: B
!! real(dp) :: C
!! type (sunu) :: D
!!
!! B = 1.0_dp
!! A = 2.0_dp * B
!! C = A
!! D = A
!!
!! 4. basic arithmetic operations for sunu structure
!! -------------------------------------------------
!!
!! type (sunu) :: A
!! type (sunu) :: B
!! type (sunu) :: C
!! real(dp) :: d
!!
!! A = 1.0_dp
!! B = 3.0_dp
!!
!! C = A + B
!! d = C
!! print *, d
!!
!! C = A / B
!! d = C
!! print *, d
!!
!! C = A * B
!! d = C
!! print *, d
!!
!!

  module skynet
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp    = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare data structure                                           <<<
!!========================================================================

! SUper-large-NUmber structure
     public :: sunu
     type sunu
         real(dp) :: mant ! mantisa
         real(dp) :: expt ! exponent
     end type sunu

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! simulate assignment operations for sunu structure
     private :: sunu_to_sunu
     private :: sunu_to_real
     private :: real_to_sunu

! simulate add operations for sunu structure
     private :: sunu_add_sunu
     private :: sunu_add_real
     private :: real_add_sunu

! simulate substract operations for sunu structure
     private :: sunu_sub_sunu
     private :: sunu_sub_real
     private :: real_sub_sunu

! simulate multiply operations for sunu structure
     private :: sunu_mul_sunu
     private :: sunu_mul_real
     private :: real_mul_sunu

! simulate divide operations for sunu structure
     private :: sunu_div_sunu
     private :: sunu_div_real
     private :: real_div_sunu

! utility subroutines
     public :: skynet_make_sunu
     public :: skynet_done_sunu
     public :: skynet_test_sunu
     public :: skynet_zero_sunu

!!========================================================================
!!>>> declare interface and module procedure                           <<<
!!========================================================================

     interface assignment (=)
         module procedure sunu_to_sunu
         module procedure sunu_to_real
         module procedure real_to_sunu
     end interface assignment (=)

     interface operator (+)
         module procedure sunu_add_sunu
         module procedure sunu_add_real
         module procedure real_add_sunu
     end interface operator (+)

     interface operator (-)
         module procedure sunu_sub_sunu
         module procedure sunu_sub_real
         module procedure real_sub_sunu
     end interface operator (-)

     interface operator (*)
         module procedure sunu_mul_sunu
         module procedure sunu_mul_real
         module procedure real_mul_sunu
     end interface operator (*)

     interface operator (/)
         module procedure sunu_div_sunu
         module procedure sunu_div_real
         module procedure real_div_sunu
     end interface operator (/)

  contains ! encapsulated functionality

!!>>> sunu_to_sunu: assign sunu structure to sunu structure
  subroutine sunu_to_sunu(s1, s2)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(out) :: s1

! sunu structure
     type (sunu), intent(in)  :: s2

     s1%mant = s2%mant
     s1%expt = s2%expt

     return
  end subroutine sunu_to_sunu

!!>>> sunu_to_real: assign sunu structure to real number
  subroutine sunu_to_real(r1, s1)
     implicit none

! external arguments
! real number
     real(dp), intent(out)   :: r1

! sunu structure
     type (sunu), intent(in) :: s1

     r1 = s1%mant * exp( s1%expt )

! check the status of r1
     if ( isnan( r1 * r1 ) ) then
         write(mystd,'(a)') 'skynet: sunu_to_real, NaN error'
         STOP
     endif ! back if ( isnan( r1 * r1 ) ) block

     return
  end subroutine sunu_to_real

!!>>> real_to_sunu: assign real number to sunu structure
  subroutine real_to_sunu(s1, r1)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(out) :: s1

! real number
     real(dp), intent(in)     :: r1

     s1%mant = r1
     s1%expt = 0.0_dp

     return
  end subroutine real_to_sunu

!!>>> sunu_add_sunu: add sunu structure to sunu structure
  type (sunu) &
  function sunu_add_sunu(s1, s2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! sunu structure
     type (sunu), intent(in) :: s2

     if ( s1%expt > s2%expt ) then
         s%mant = s1%mant + s2%mant * exp( s2%expt - s1%expt )
         s%expt = s1%expt
     else
         s%mant = s2%mant + s1%mant * exp( s1%expt - s2%expt )
         s%expt = s2%expt
     endif ! back if ( s1%expt > s2%expt ) block

     return
  end function sunu_add_sunu

!!>>> sunu_add_real: add sunu structure to real number
  type (sunu) &
  function sunu_add_real(s1, r2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! real number
     real(dp), intent(in)    :: r2

! local variables
! sunu structure
     type (sunu) :: s2

     s2%mant = r2
     s2%expt = 0.0_dp

     s = sunu_add_sunu(s1, s2)

     return
  end function sunu_add_real

!!>>> real_add_sunu: add real number to sunu structure
  type (sunu) &
  function real_add_sunu(r1, s2) result(s)
     implicit none

! external arguments
! real number
     real(dp), intent(in)    :: r1

! sunu structure
     type (sunu), intent(in) :: s2

! local variables
! sunu structure
     type (sunu) :: s1

     s1%mant = r1
     s1%expt = 0.0_dp

     s = sunu_add_sunu(s1, s2)

     return
  end function real_add_sunu

!!>>> sunu_sub_sunu: substract sunu structure from sunu structure
  type (sunu) &
  function sunu_sub_sunu(s1, s2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! sunu structure
     type (sunu), intent(in) :: s2

     if ( s1%expt > s2%expt ) then
         s%mant = s1%mant - s2%mant * exp( s2%expt - s1%expt )
         s%expt = s1%expt
     else
         s%mant =-s2%mant + s1%mant * exp( s1%expt - s2%expt )
         s%expt = s2%expt
     endif ! back if ( s1%expt > s2%expt ) block

     return
  end function sunu_sub_sunu

!!>>> sunu_sub_real: substract real number from sunu structure
  type (sunu) &
  function sunu_sub_real(s1, r2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! real number
     real(dp), intent(in)    :: r2

! local variables
! sunu structure
     type (sunu) :: s2

     s2%mant = r2
     s2%expt = 0.0_dp

     s = sunu_sub_sunu(s1, s2)

     return
  end function sunu_sub_real

!!>>> real_sub_sunu: substract sunu structure from real number
  type (sunu) &
  function real_sub_sunu(r1, s2) result(s)
     implicit none

! external arguments
! real number
     real(dp), intent(in)    :: r1

! sunu structure
     type (sunu), intent(in) :: s2

! local variables
! sunu structure
     type (sunu) :: s1

     s1%mant = r1
     s1%expt = 0.0_dp

     s = sunu_sub_sunu(s1, s2)

     return
  end function real_sub_sunu

!!>>> sunu_mul_sunu: multiply sunu structure with sunu structure
  type (sunu) &
  function sunu_mul_sunu(s1, s2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! sunu structure
     type (sunu), intent(in) :: s2

     s%mant = s1%mant * s2%mant
     s%expt = s1%expt + s2%expt

     return
  end function sunu_mul_sunu

!!>>> sunu_mul_real: multiply sunu structure with real number
  type (sunu) &
  function sunu_mul_real(s1, r2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! real number
     real(dp), intent(in)    :: r2

     s%mant = s1%mant * r2
     s%expt = s1%expt

     return
  end function sunu_mul_real

!!>>> real_mul_sunu: multiply real number with sunu structure
  type (sunu) &
  function real_mul_sunu(r1, s2) result(s)
     implicit none

! external arguments
! real number
     real(dp), intent(in)    :: r1

! sunu structure
     type (sunu), intent(in) :: s2

     s%mant = r1 * s2%mant
     s%expt = s2%expt

     return
  end function real_mul_sunu

!!>>> sunu_div_sunu: divide sunu structure by sunu structure
  type (sunu) &
  function sunu_div_sunu(s1, s2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! sunu structure
     type (sunu), intent(in) :: s2

     s%mant = s1%mant / s2%mant
     s%expt = s1%expt - s2%expt

     return
  end function sunu_div_sunu

!!>>> sunu_div_real: divide sunu structure by real number
  type (sunu) &
  function sunu_div_real(s1, r2) result(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s1

! real number
     real(dp), intent(in)    :: r2

     s%mant = s1%mant / r2
     s%expt = s1%expt

     return
  end function sunu_div_real

!!>>> real_div_sunu: divide real number by sunu structure
  type (sunu) &
  function real_div_sunu(r1, s2) result(s)
     implicit none

! external arguments
! real number
     real(dp), intent(in)    :: r1

! sunu structure
     type (sunu), intent(in) :: s2

     s%mant = r1 / s2%mant
     s%expt = -s2%expt

     return
  end function real_div_sunu

!!>>> skynet_make_sunu: build sunu structure from a pair real number
  subroutine skynet_make_sunu(mant, expt, s)
     implicit none

! external arguments
! mantisa part
     real(dp), intent(in)     :: mant

! exponent part
     real(dp), intent(in)     :: expt

! sunu structure
     type (sunu), intent(out) :: s

     s%mant = mant
     s%expt = expt

     return
  end subroutine skynet_make_sunu

!!>>> skynet_done_sunu: make balance for the sunu structure
  subroutine skynet_done_sunu(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(inout) :: s

     if ( s%mant == 0.0_dp ) then
         s%expt = 0.0_dp
     endif ! back if ( s%mant == 0.0_dp ) block

     s%expt = s%expt + log( abs( s%mant ) )
     if ( s%mant > 0.0_dp ) then
         s%mant = 1.0_dp
     else
         s%mant =-1.0_dp
     endif ! back if ( s%mant > 0.0_dp ) block

     return
  end subroutine skynet_done_sunu

!!>>> skynet_test_sunu: check the status of sunu structure
  subroutine skynet_test_sunu(s)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s

     if ( isnan( s%mant * s%mant ) ) then
         write(mystd,'(a)') 'skynet: skynet_test_sunu, mantisa part is infinity'
         STOP
     endif ! back if ( isnan( s%mant * s%mant ) ) block

     if ( isnan( s%expt * s%expt ) ) then
         write(mystd,'(a)') 'skynet: skynet_test_sunu, exponent part is infinity'
         STOP
     endif ! back if ( isnan( s%expt * s%expt ) ) block

     return
  end subroutine skynet_test_sunu

!!>>> skynet_zero_sunu: test if sunu structure is zero in real form
  logical &
  function skynet_zero_sunu(s) result(z)
     implicit none

! external arguments
! sunu structure
     type (sunu), intent(in) :: s

     if ( s%mant == 0.0_dp ) then
         z = .true.
     else
         z = .false.
     endif ! back if ( s%mant == 0.0_dp ) block

     return
  end function skynet_zero_sunu

  end module skynet
