! make single particle related matrices, including:
! crystal field, spin-orbital coupling, Coulomb inteartion U
subroutine atomic_make_spmat()
    use m_spmat

    implicit none

    ! first, allocate memory for these single particle matrices
    call alloc_m_spmat()

    ! second, make crystal field and spin-orbital coupling 
    if (itask == 1) then ! model calculation
        if (icf > 0) then
            call atomic_read_cf()
        else
            cfmat = czero
        endif

        if (isoc > 0) then
            call atomic_make_soc()
        else
            socmat = czero
        endif
    else ! material calculation
        ! read the eimp (CF+SOC) matrices on natural basis
        ! this matrix should be a diagonal matrix, and we just
        ! need its diagonal elements
        call atomic_read_eimp()
        ! read the transformation matrices used to transfer eimp 
        ! from the standard real orbitals basis to natural basis 
        call atomic_read_umat()
    endif

    ! third, make Coulomb interaction U
    if (icu == 1) then
        ! Kanamori parameters type
        call atomic_make_cumat_kanamori()
    else
        ! Slater-Cordon parameters type
        call atomic_make_cumat_slater()
    endif

    return
end subroutine atomic_make_spmat 

!>>> make spin-orbital coupling 
subroutine atomic_make_soc()
    use control
    use m_spmat

    implicit none

    if (nband == 3) then
        call atomic_make_soc_3band()
    elseif(nband == 5) then
        call atomic_make_soc_5band()
    elseif(nband == 7) then
        call atomic_make_soc_7band()
    else
        call atomic_print_error('atomic_make_soc', 'not implementd!')
    endif

    return 
end subroutine atomic_make_soc

!>>> make spin-orbital coupling matrix for 3 band
subroutine atomic_make_soc_3band()
    use constants
    use control
    use m_spmat

    implicit none

    ! we use the same default orbital order as in WANNIER90 package 
    ! the t2g orbital order of d(l=2) is 
    ! |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, |dxy,up>, |dxy,dn>
    ! the coressponding p(l=1) orbital order is
    ! |py,up>,  |py,dn>,  |px,up>,  |px,dn>,  |pz,up>,  |pz,dn>
    socmat = czero
    socmat(1,3) =  czi;     socmat(1,6) =  -czi
    socmat(2,4) = -czi;     socmat(2,5) =  -czi
    socmat(3,1) = -czi;     socmat(3,6) =  cone
    socmat(4,2) =  czi;     socmat(4,5) = -cone
    socmat(5,2) =  czi;     socmat(5,4) = -cone
    socmat(6,1) =  czi;     socmat(6,3) =  cone
    ! please note: minus sign for spin-orbital coupling strength         
    socmat = -socmat * lambda / two 

end subroutine atomic_make_soc_3band

!>>> make spin-orbital coupling matrix for 5 band
subroutine atomic_make_soc_5band()
    use constants
    use control
    use m_sp_mat

    implicit none

    ! local variables
    ! sqrt(3)
    real(dp) :: sqrt3

    sqrt3 = sqrt(3.0_dp)

    ! we use the same default orbital order as in WANNIER90 package 
    ! the orbital order is: 
    ! |dz2,up>, |dz2,dn>, |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, |dx2-y2,up>, |dx2-y2,dn>, |dxy,up>, |dxy,dn>
    socmat = czero
    socmat(1, 4) = -sqrt3*cone;      socmat(1, 6) = sqrt3*czi
    socmat(2, 3) =  sqrt3*cone;      socmat(2, 5) = sqrt3*czi
    socmat(3, 2) =  sqrt3*cone;      socmat(3, 5) = -czi
    socmat(3, 8) =       -cone;      socmat(3, 10)= czi
    socmat(4, 1) = -sqrt3*cone;      socmat(4, 6) = czi
    socmat(4, 7) =        cone;      socmat(4, 9) = czi
    socmat(5, 2) =  -sqrt3*czi;      socmat(5, 3) = czi
    socmat(5, 8) =        -czi;      socmat(5, 10)=-cone
    socmat(6, 1) =  -sqrt3*czi;      socmat(6, 4) = -czi
    socmat(6, 7) =        -czi;      socmat(6, 9) = cone
    socmat(7, 4) =        cone;      socmat(7, 6) = czi
    socmat(7, 9) =     -two*czi
    socmat(8, 3) =       -cone;      socmat(8, 5) = czi
    socmat(8,10) =      two*czi
    socmat(9, 4) =        -czi;      socmat(9, 6) = cone
    socmat(9, 7) =     two*czi;    
    socmat(10,3) =        -czi;      socmat(10,5) = -cone
    socmat(10,8) =    -two*czi; 

    ! scale the SOC strength lambda
    socmat = socmat * lambda / two

    return 
end subroutine atomic_make_soc_5band

!>>> make spin-orbital coupling matrix for 7 band
subroutine atomic_make_soc_7band()
    use constants
    use control
    use m_sp_mat

    implicit none
    
    write(mystd,*) "not implemented now!"

    return
end subroutine atomic_make_soc_7band

!>>> make Coulomb interaction U, this subroutine is taken from 
!>>> duliang's atomic software
subroutine atomic_make_cumat_kanamori()
    use constants
    use control
    use m_spmat

    implicit none

    ! local varibales
    ! loop index over orbits
    integer :: alpha, betta
    integer :: delta, gamma
    ! band index and spin index
    integer :: aband, bband 
    integer :: dband, gband 
    integer :: aspin, bspin 
    integer :: dspin, gspin 
    ! auxiliary variables
    real(dp) :: dtmp

    ! initialize cumat to zero
    cumat = czero

    ! loop for creation operators
    alphaloop: do alpha=1,norbs-1
    bettaloop: do betta=alpha+1,norbs

       ! loop for annihilation operators
       gammaloop: do gamma=1,norbs-1
       deltaloop: do delta=gamma+1,norbs
           aband = (alpha+1)/2; aspin = mod(alpha,2)
           bband = (betta+1)/2; bspin = mod(betta,2)
           gband = (gamma+1)/2; gspin = mod(gamma,2)
           dband = (delta+1)/2; dspin = mod(delta,2)

           ! here we use "res" due to overlap between "Uv and Jz"
           dtmp = zero

           ! intraorbital Coulomb interaction
           if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
               if ((aband.eq.bband) .and. (aspin.ne.bspin)) then
                   dtmp = dtmp + Uc
               endif
           endif

           ! interorbital Coulomb interaction
           if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
               if (aband .ne. bband) then
                   dtmp = dtmp + Uv
               endif
           endif

           ! Hund's exchange interaction 
           if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
               if ((aband.ne.bband) .and. (aspin.eq.bspin)) then
                   dtmp = dtmp - Jz
               endif
           endif
          
           ! spin flip term
           if ((aband.eq.gband) .and. (bband.eq.dband)) then
               if ((aspin.ne.gspin) .and. (bspin.ne.dspin) .and. (aspin.ne.bspin)) then
                   dtmp = dtmp - Js
               endif
           endif
         
           ! pair hopping term
           if ((aband.eq.bband) .and. (dband.eq.gband) .and. (aband.ne.dband)) then
               if ((aspin.ne.bspin) .and. (dspin.ne.gspin) .and. (aspin.eq.gspin)) then
                   dtmp = dtmp + Jp
               endif
           endif
                
           cumat(alpha, betta, delta, gamma) = dtmp

       enddo deltaloop ! over delta={gamma+1,norbs} loop
       enddo gammaloop ! over gamma={1,norbs-1} loop
    enddo bettaloop ! over betta={alpha+1,norbs} loop
    enddo alphaloop ! over alpha={1,norbs-1} loop

    return
end subroutine atomic_make_cumat_kanamori

!>>> make Coulomb interation U, Slater Integral version
subroutine atomic_make_cumat_slater()
    use constants
    use control
    use m_spmat

    implicit none

    ! local variables
    ! Slater-Cordon parameters
    real(dp), allocatable :: slater_cordon(:)
    ! gaunt coefficients
    real(dp), allocatable :: gaunt(:,:,:)

    ! orbital momentum quantum number
    integer :: l
    ! loop index
    integer :: k
    integer :: alpha, betta
    integer :: delta, gamma
    integer :: aband, aspin
    integer :: bband, bspin
    integer :: dband, dspin
    integer :: gband, gspin
    ! real(dp) auxiliary variables
    real(dp) :: res


    ! allocate memory for slater_cordon and gaunt and then build them
    if (nband == 5) then
        l = 2
        allocate(slater_cordon(0:2*l))     
        slater_cordon = zero
        slater_cordon(0) = F0
        slater_cordon(2) = F2
        slater_cordon(4) = F4
        allocate(gaunt(-l:l, -l:l, 0:2*l)
        call atomic_gaunt_5band(gaunt) 
    elseif(nband == 7) then
        l = 3
        allocate(slater_cordon(0:2*l))     
        slater_cordon = zero
        slater_cordon(0) = F0
        slater_cordon(2) = F2
        slater_cordon(4) = F4
        slater_cordon(6) = F6
        allocate(gaunt(-l:l, -l:l, 0:2*l)
        call atomic_gaunt_7band(gaunt)
    else
        call atomic_print_error('atomic_make_cumat_slater', 'nband is wrong!')
    endif

    ! make Coulomb interaction U matrix
    do alpha=1,norbs
    do betta=1,norbs
        aband = (alpha-1)/2-l
        bband = (betta-1)/2-l
        aspin = mod(alpha, 2)
        bspin = mod(betta, 2)

        do delta=1,norbs
        do gamma=1,norbs
            dband = (delta-1)/2-l
            gband = (gamma-1)/2-l
            dspin = mod(delta, 2)
            gspin = mod(gamma, 2)

            if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle

            if ((aband + bband) .ne. (dband + gband)) cycle
            if ((aspin .ne. gspin) .or. (bspin .ne. dspin)) cycle

            ! fixme  wrong in rpp 69, 2061 (2006), ref: prb 75, 155113 (2007)
            ! exchange of dband and bband
            res = zero
            do k=0, 2*l, 2
                res = res + gaunt(aband, gband, k) * gaunt(dband, bband, k) * slater_cordon(k)
            enddo
            cumat(alpha, betta, delta, gamma) = res
        enddo ! over gamma={1,norbs} loop
        enddo ! over delta={1,norbs} loop

    enddo ! over betta={1,norbs} loop
    enddo ! over alpha={1,norbs} loop

    cumat = half * cumat

    ! deallocate memory
    if (allocated(slater_cordon)) deallocate(slater_cordon) 
    if (allocated(gaunt))         deallocate(gaunt)

    return
end subroutine atomic_make_cumat_slater
