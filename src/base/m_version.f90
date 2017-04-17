!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : version
!!! source  : m_version.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/26/2017 by li huang (created)
!!!           04/17/2017 by li huang (last modified)
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
!! use version, only : FULL_VER
!! implicit none
!!
!! print *, FULL_VER
!!
!!

  module version
     implicit none

! version string: version number + date info. + status info.
     character(len=20), public, parameter :: FULL_VER = 'v0.7.0 @ 2017.01.31D'

! version string: only version number
     character(len=06), public, parameter :: CURR_VER = 'v0.7.0'

! version string: only date info.
     character(len=11), public, parameter :: DATE_VER = '2017.01.31'

! version string: only status info., D means devel, T testing, R released.
     character(len=01), public, parameter :: STAT_VER = 'D'

! version string: author info.
     character(len=46), public, parameter :: AUTH_VER = 'by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'

! version string: email info.
     character(len=22), public, parameter :: MAIL_VER = 'lihuang.dmft@gmail.com'

! version string: license info.
     character(len=36), public, parameter :: GPL3_VER = 'GNU General Public License version 3'

  end module version
