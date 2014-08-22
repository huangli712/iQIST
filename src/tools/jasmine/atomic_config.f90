!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_config
!!! source  : atomic_config.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : set control parameters 
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_config: read config parameters from file 'atom.config.in'
  subroutine atomic_config()
     use constants, only : dp, mytmp
     use control

     use parser, only : p_create, p_destroy, p_parse, p_get
  
     implicit none
  
! local variables
! file status
     logical :: exists
     
!----------------------------------------------------------------
     itask  = 1           ! type of task
     ictqmc = 1           ! type of CTQMC algorithm
     icf    = 0           ! type of crystal field
     isoc   = 0           ! type of spin-orbital coupling (SOC)
     icu    = 1           ! type of Coulomb interaction
     
!---------------------------------------------------------------- 
     nband = 1            ! number of bands
     nspin = 2            ! number of spins
     norbs = nband*nspin  ! number of orbits
     ncfgs = 2**norbs     ! number of many-body configurations

!----------------------------------------------------------------
     Uc = 2.00_dp         ! intraorbital Coulomb interaction
     Uv = 2.00_dp         ! interorbital Coulomb interaction
     Jz = 0.00_dp         ! Hund's exchange interaction
     Js = 0.00_dp         ! spin-flip interaction
     Jp = 0.00_dp         ! pair-hopping interaction
  
!----------------------------------------------------------------
     Ud = 2.00_dp         ! Ud
     JH = 0.00_dp         ! JH
     F0 = 0.00_dp         ! F0
     F2 = 0.00_dp         ! F2
     F4 = 0.00_dp         ! F4
     F6 = 0.00_dp         ! F6
  
!----------------------------------------------------------------
     lambda = 0.00_dp     ! spin-orbit coupling parameter
     mune   = 0.00_dp     ! chemical potential
  
!----------------------------------------------------------------
! file status
     exists = .false.
  
! inquire the input file status
     inquire( file="atom.config.in", exist=exists )
  
! read parameters from atom.config.in
     if ( exists .eqv. .true. ) then
!----------------------------------------------------------------
         call p_create()
         call p_parse('atom.config.in')
!----------------------------------------------------------------
         call p_get('nband', nband)
!----------------------------------------------------------------
         call p_get('itask',  itask)
         call p_get('ictqmc', ictqmc)
         call p_get('icf',    icf)
         call p_get('isoc',   isoc)
         call p_get('icu',    icu)
!----------------------------------------------------------------
         call p_get('Uc',     Uc) 
         call p_get('Uv',     Uv) 
         call p_get('Jz',     Jz) 
         call p_get('Js',     Js) 
         call p_get('Jp',     Jp) 
!----------------------------------------------------------------
         call p_get('Ud',     Ud) 
         call p_get('JH',     JH) 
!----------------------------------------------------------------
         call p_get('lambda', lambda)
         call p_get('mune',   mune)
!----------------------------------------------------------------
         call p_destroy()
!----------------------------------------------------------------
! calculate the norbs and ncfgs
         norbs = nband * nspin
         ncfgs = 2 ** norbs 

! calculate F0, F2, F4, F6 here
         F0 = Ud
         if (nband == 5) then
             F2 = JH * 14.0_dp / 1.625_dp 
             F4 = 0.625_dp * F2
         elseif(nband == 7) then
             F2 = JH * 6435.0_dp / (286.0_dp + (195.0_dp * 451.0_dp / 675.0_dp) &
                                            + (250.0_dp * 1001.0_dp / 2025.0_dp))
             F4 = 451.0_dp / 675.0_dp * F2
             F6 = 1001.0_dp / 2025.0_dp * F2
         endif
!----------------------------------------------------------------
     else
         call s_print_error('atomic_config', 'no file atom.config.in !')
     endif
  
     return
  end subroutine atomic_config
