!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : version
!!! source  : m_version.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/26/2017 by li huang (created)
!!!           01/26/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a version string.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! It is a common module which defines the current version of iQIST.
!!
!! Usage
!! =====
!!
!! use version, only : __version__
!! implicit none
!!
!! print *, __version__
!!
!!

  module version
     implicit none

! version string
     character (len=20), public, parameter :: __version__ = 'v0.6.8 @ 2017.01.26D'

  end module version
