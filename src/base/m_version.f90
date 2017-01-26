!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : stack
!!! source  : m_stack.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/14/2009 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : the purpose of this module is to define a stack-type (LIFO)
!!!           data structure in fortran version.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! It is a common module which defines some common used numerical/physical
!! constants. We always need it.
!!
!! Usage
!! =====
!!
!! 1. import completely
!! --------------------
!!
!! use constants
!! implicit none
!!
!! real(dp) :: A
!! A = one
!!
!! 2. import partially
!! -------------------
!!
!! use constants, only : dp, one
!! implicit none
!!
!! real(dp) :: A
!! A = one
!!
!!

  module version
     implicit none

     character (len=20), public, parameter :: __version__ = 'v0.6.8 @ 2017.01.26D'
  end module version
