!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : spring
!!! source  : m_spring.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/03/2008 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : the purpose of this module is to define several modern, fast
!!!           highly reliable, ease-to-use, and state-of-the-art pseudo-
!!!           random number generators, which are essential in massively
!!!           parallel quantum Monte Carlo simulations.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! the following two random number generators (generates a random number
!! between 0 and 1, real double precision) are supported by now.
!!
!! (A) MT19937
!! Mersenne Twister pseudorandom number generator
!! by Makoto Matsumoto and Takuji Nishimura
!! << Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator >>
!! ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, January pp.3-30 (1998)
!!
!! (B) SFMT
!! SIMD-oriented Fast Mersenne Twister
!! by Mutsuo Saito and Makoto Matsumoto
!! << SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number Generator >>
!! Monte Carlo and Quasi-Monte Carlo Methods 2006, Springer, 2008, pp. 607-622
!!
!! Usage
!! =====
!!
!! 1. use MT19937
!! --------------
!!
!! use spring
!! call spring_mt_init(seed)
!! r = spring_mt_stream()
!! r = spring_mt_string()
!!
!! 2. use SFMT
!! -----------
!!
!! use spring
!! call spring_sfmt_init(seed)
!! r = spring_sfmt_stream()
!! r = spring_sfmt_string()
!!
!! note: since SFMT has a better performance, it is preferable
!!
!!

  module spring
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! kind types for 32-bit signed integers
     integer, private, parameter :: int32  = selected_int_kind( 9)

! kind types for 64-bit signed integers
     integer, private, parameter :: int64  = selected_int_kind(18)

! kind types for ieee 754/IEC 60559 single precision reals
     integer, private, parameter :: ieee32 = selected_real_kind( 6,  37)

! kind types for ieee 754/IEC 60559 double precision reals
     integer, private, parameter :: ieee64 = selected_real_kind(15, 307)

!!========================================================================
!!>>> declare common parameters for MT19937 generator                  <<<
!!========================================================================

! period parameters
     integer(int32), private, parameter :: N = 624_int32

! period parameters
     integer(int32), private, parameter :: M = 397_int32

!!========================================================================
!!>>> declare common parameters for SFMT generator                     <<<
!!========================================================================

! Mersenne Exponent. the period of the sequence is a multiple of 2^{ME}-1
     integer(int32), private, parameter :: ME  = 19937

! SFMT has an internal state array of 128-bit integers, and NS is its size.
     integer(int32), private, parameter :: NS  = ME / 128 + 1

! N32 is the size of internal state array when regarded as an array of 32-bit integers
     integer(int32), private, parameter :: N32 = NS * 4

! N64 is the size of internal state array when regarded as an array of 64-bit integers
     integer(int32), private, parameter :: N64 = NS * 2

! a parity check vector which certificate the period of 2^{ME}
     integer(int32), private, parameter :: parity(0:3) = (/1, 0, 0, 331998852/)

!!========================================================================
!!>>> declare common variables for MT19937 generator                   <<<
!!========================================================================

! states pointer, if mti == -1, itialization is necessary
     integer(int32), private, save :: mti = -1

! states array for twist generator
     integer(int64), private, save :: mt(0:N-1)

!!========================================================================
!!>>> declare common variables for SFMT generator                      <<<
!!========================================================================

! index counter to the 32-bit internal state array
     integer(int32), private, save :: idx = -1

! the 32-bit integer pointer to the 128-bit internal state array
     integer(int32), private, save :: pt32(0:N32-1)

! the 64-bit integer pointer to the 128-bit internal state array
     integer(int64), private, save :: pt64(0:N64-1)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public  :: spring_mt_init     ! MT19937 generator public  subroutines
     public  :: spring_mt_stream   !
     public  :: spring_mt_string   !
                                   !
     private :: spring_mt_source   ! MT19937 generator private subroutines
                                   !
     public  :: spring_sfmt_init   ! SFMT    generator public  subroutines
     public  :: spring_sfmt_stream !
     public  :: spring_sfmt_string !
                                   !
     private :: spring_sfmt_source ! SFMT    generator private subroutines
     private :: spring_sfmt_kernel !
     private :: spring_sfmt_core   !

  contains ! encapsulated functionality

!!========================================================================
!!>>> MT19937 random number generator subroutines                      <<<
!!========================================================================

!!>>> spring_mt_init: initializes the MT19937 generator with "seed"
  subroutine spring_mt_init(seed)
     implicit none

! external arguments
! seed for random number generator
     integer(int32), intent(in) :: seed

! local variables
! loop index
     integer :: i

! save seed
     mt(0) = seed

! set the seed using values suggested by Matsumoto & Nishimura, using
! a generator by Knuth. See original source for details
     do i=1,N-1
         mt(i) = iand(4294967295_int64,1812433253_int64*ieor(mt(i-1),ishft(mt(i-1),-30_int64))+i)
     enddo ! over i={1,N-1} loop

     mti = N

     return
  end subroutine spring_mt_init

!!>>> spring_mt_stream: obtain a psuedo random real number in the range
!!>>> (0,1), i.e., a number greater than 0 and less than 1
  function spring_mt_stream() result(r)
     implicit none

! local parameters
! pre-calculated to avoid division below
     real(ieee64), parameter :: factor = 1.0_ieee64 / 4294967296.0_ieee64

! local variables
! return type
     real(ieee64) :: r

! compute it
     r = (real(spring_mt_source(),ieee64) + 0.5_ieee64) * factor

     return
  end function spring_mt_stream

!!>>> spring_mt_string: obtain a psuedo random real number in the range
!!>>> [0,1], i.e., a number greater than or equal to 0 and less than or
!!>>> equal to 1
  function spring_mt_string() result(r)
     implicit none

! local parameters
! pre-calculated to avoid division below
     real(ieee64), parameter :: factor = 1.0_ieee64 / 4294967295.0_ieee64

! local variables
! return type
     real(ieee64) :: r

! compute it
     r = real(spring_mt_source(),ieee64) * factor

     return
  end function spring_mt_string

!!>>> spring_mt_source: obtain the next 64-bit integer in the psuedo
!!>>> random sequence
  function spring_mt_source() result(r)
     implicit none

! local parameters
     integer(int64), parameter :: MAGIC(0:1) = (/ 0_int64, -1727483681_int64 /)

! period parameters
     integer(int64), parameter :: UPPER_MASK =  2147483648_int64
     integer(int64), parameter :: LOWER_MASK =  2147483647_int64

! tempering parameters
     integer(int64), parameter :: TEMPERING_B = -1658038656_int64
     integer(int64), parameter :: TEMPERING_C =  -272236544_int64

! local variables
! return type
     integer(int64) :: r

! loop index
     integer(int32) :: k

! generate N words at a time
     if ( mti >= N ) then
! the value -1 acts as a flag saying that the seed has not been set.
         if ( mti == -1 ) call spring_mt_init(4357_int32)

! fill the mt array
         do k=0,N-M-1
             r = ior(iand(mt(k),UPPER_MASK),iand(mt(k+1),LOWER_MASK))
             mt(k) = ieor(ieor(mt(k + M),ishft(r,-1_int64)),MAGIC(iand(r,1_int64)))
         enddo ! over k={0,N-M-1} loop

         do k=N-M,N-2
             r = ior(iand(mt(k),UPPER_MASK),iand(mt(k+1),LOWER_MASK))
             mt(k) = ieor(ieor(mt(k + (M - N)),ishft(r,-1_int64)),MAGIC(iand(r,1_int64)))
         enddo ! over k={N-M,N-2} loop

         r = ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
         mt(N-1) = ieor(ieor(mt(M-1),ishft(r,-1)),MAGIC(iand(r,1_int64)))

! start using the array from first element
         mti = 0
     endif ! back if ( mti >= N ) block

! here is where we calculate the number with a series of transformations
     r = mt(mti)
     mti = mti + 1

     r = ieor(r,ishft(r,-11))
     r = iand(4294967295_int64,ieor(r,iand(ishft(r, 7),TEMPERING_B)))
     r = iand(4294967295_int64,ieor(r,iand(ishft(r,15),TEMPERING_C)))
     r = ieor(r,ishft(r,-18))

     return
  end function spring_mt_source

!!========================================================================
!!>>> SFMT random number generator subroutines                         <<<
!!========================================================================

!!>>> spring_sfmt_init: this function initializes the internal state array
!!>>> with a 32-bit integer seed
  subroutine spring_sfmt_init(seed)
     implicit none

! external arguments
! seed for random number generator
     integer(int32), intent(in) :: seed

! local variables
! loop index
     integer :: i
     integer :: j

! status flag
     integer :: inner = 0

! dummy variables
     integer(int32) :: work
     integer(int64) :: temp

! setup idx
     idx = N32

! setup pt32 array
     pt32(0) = seed
     do i=1,N32-1
         temp = ieor( pt32(i-1), ishft( pt32(i-1), -30 ) )
         temp = 1812433253_int64 * temp + i
         pt32(i) = int( ibits( temp, 0, 32 ), kind = int32 )
     enddo ! over i={1,N32-1} loop

! check period of random number generator
     do i=0,3
         work = iand(pt32(i), parity(i))
         do j=0,31
             inner = ieor(inner, iand(work, 1))
             work = ishft(work,-1)
         enddo ! over j={0,31} loop
     enddo ! over i={0,3} loop

! the period is OK? re-adjust pt32 array
     if ( inner /= 1 ) then
         adjust_period_loop: do i=0,3
             work = 1
             do j=0,31
                 if ( iand(work, parity(i)) /= 0 ) then
                     pt32(i) = ieor(pt32(i), work)
                     EXIT adjust_period_loop
                 endif ! back if ( iand(work, parity(i)) /= 0 ) block
                 work = ishft(work,1)
             enddo ! over j={0,31} loop
         enddo adjust_period_loop ! over i={0,3} loop
     endif ! back if ( inner /= 1 ) block

! setup pt64 array
     do i=0,N64-1
         pt64(i) = pt32(i*2 + 1)
         temp = pt32(i*2)
         temp = iand(temp, 4294967295_int64)
         pt64(i) = ior(ishft(pt64(i), 32), temp)
     enddo ! over i={0,N64-1} loop

     return
  end subroutine spring_sfmt_init

!!>>> spring_sfmt_stream: generates a pseudo random number on (0,1)
  function spring_sfmt_stream() result(r)
     implicit none

! local parameters
! pre-calculated to avoid division below
     real(ieee64), parameter :: factor = 1.0_ieee64 / 4294967296.0_ieee64

! local variables
! return type
     real(ieee64) :: r

! compute it
     r = (real(spring_sfmt_source(),ieee64) + 0.5_ieee64) * factor

     return
  end function spring_sfmt_stream

!!>>> spring_sfmt_string: generates a pseudo random number on [0,1]
  function spring_sfmt_string() result(r)
     implicit none

! local parameters
! pre-calculated to avoid division below
     real(ieee64), parameter :: factor = 1.0_ieee64 / 4294967295.0_ieee64

! local variables
! return type
     real(ieee64) :: r

! compute it
     r = real(spring_sfmt_source(),ieee64) * factor

     return
  end function spring_sfmt_string

!!>>> spring_sfmt_source: this function generates and returns 32-bit
!!>>> pseudo random number
  function spring_sfmt_source() result(r)
     implicit none

! local variables
! return type
     integer(int64) :: r

     if ( idx >= N32 ) then
         call spring_sfmt_kernel()
         idx = 0
     endif ! back if ( idx >= N32 ) block

     if ( btest(idx, 0) ) then
         r = ibits(pt64(idx/2), 32, 32)
     else
         r = ibits(pt64(idx/2), 0 , 32)
     endif ! back if ( btest(idx, 0) ) block

     idx = idx + 1

     return
  end function spring_sfmt_source

!!>>> spring_sfmt_kernel: this function fills the internal state array with
!!>>> pseudo random integers
  subroutine spring_sfmt_kernel()
     implicit none

! local variables
! loop index
     integer :: i

     integer(int64) :: r1Top
     integer(int64) :: r1Btm

     integer(int64) :: r2Top
     integer(int64) :: r2Btm

     r1Btm = pt64(N64 - 4)
     r1Top = pt64(N64 - 3)
     r2Btm = pt64(N64 - 2)
     r2Top = pt64(N64 - 1)

     do i=0,NS-122-1
         call spring_sfmt_core( pt64(i*2 + 1), pt64(i*2),   &
                                pt64(i*2 + 1), pt64(i*2),   &
                                pt64((i + 122)*2+1),        &
                                pt64((i + 122)*2),          &
                                r1Top, r1Btm, r2Top, r2Btm )
         r1Btm = r2Btm
         r1Top = r2Top
         r2Btm = pt64(i*2)
         r2Top = pt64(i*2 + 1)
     enddo ! over i={0,NS-122-1} loop

     do i=NS-122,NS-1
         call spring_sfmt_core( pt64(i*2 + 1), pt64(i*2),   &
                                pt64(i*2 + 1), pt64(i*2),   &
                                pt64((i + 122 - NS)*2 + 1), &
                                pt64((i + 122 - NS)*2),     &
                                r1Top, r1Btm, r2Top, r2Btm )
         r1Btm = r2Btm
         r1Top = r2Top
         r2Btm = pt64(i*2)
         r2Top = pt64(i*2 + 1)
     enddo ! over i={NS-122,NS-1} loop

     return
  end subroutine spring_sfmt_kernel

!!>>> spring_sfmt_core: this function represents the recursion formula
  subroutine spring_sfmt_core(rTop, rBtm, aTop, aBtm, bTop, bBtm, cTop, cBtm, dTop, dBtm)
     implicit none

! external arguments
     integer(int64), intent(in) :: aTop
     integer(int64), intent(in) :: aBtm

     integer(int64), intent(in) :: bTop
     integer(int64), intent(in) :: bBtm

     integer(int64), intent(in) :: cTop
     integer(int64), intent(in) :: cBtm

     integer(int64), intent(in) :: dTop
     integer(int64), intent(in) :: dBtm

     integer(int64), intent(out) :: rTop
     integer(int64), intent(out) :: rBtm

! local variables
     integer(int64) :: xTop
     integer(int64) :: xBtm

     integer(int64) :: yTop
     integer(int64) :: yBtm

     xTop = ior(ishft(aTop, 8), ishft(aBtm, -56))
     xBtm = ishft(aBtm, 8)

     yTop = ishft(cTop, -8)
     yBtm = ior(ishft(cTop, 56), ishft(cBtm, -8))

     rBtm = ieor(aBtm, xBtm)
     rTop = ieor(aTop, xTop)

     rBtm = ieor(rBtm, iand(ishft(bBtm, -11), 8667995624701935_int64))
     rTop = ieor(rTop, iand(ishft(bTop, -11), 9007156306837503_int64))

     rBtm = ieor(rBtm, yBtm)
     rTop = ieor(rTop, yTop)

     rBtm = ieor(rBtm, ishft(iand(dBtm, 70364449226751_int64), 18))
     rTop = ieor(rTop, ishft(iand(dTop, 70364449226751_int64), 18))

     return
  end subroutine spring_sfmt_core

  end module spring
