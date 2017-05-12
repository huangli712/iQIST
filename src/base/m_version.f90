!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : version
!!! source  : m_version.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/26/2017 by li huang (created)
!!!           05/12/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define version strings.
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
!! use version, only : FULL_VER
!! implicit none
!!
!! print *, FULL_VER
!!
!!

  module version
     implicit none

!!
!! @var FULL_VER
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: FULL_VER = 'v0.7.0 @ 2017.01.31D'

!!
!! @var CURR_VER
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: CURR_VER = 'v0.7.0'

!!
!! @var DATE_VER
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: DATE_VER = '2017.01.31'

!!
!! @var STAT_VER
!!
!! version string, only status info., D means devel, T testing, R released.
!!
     character(len=01), public, parameter :: STAT_VER = 'D'

!!
!! @var AUTH_VER
!!
!! version string, author info.
!!
     character(len=46), public, parameter :: AUTH_VER = 'by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'

!!
!! @var MAIL_VER
!!
!! version string, email info.
!!
     character(len=22), public, parameter :: MAIL_VER = 'lihuang.dmft@gmail.com'

!!
!! @var GPL3_VER
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: GPL3_VER = 'GNU General Public License version 3'

  end module version
