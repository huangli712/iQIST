!    This program analyze_lavender_soc reads in 
!    single-particle Green's funciton from
!        solver.grn.dat
!        solver.sgm.dat
!    and
!    two-particle functions from
!        analyze.twop_pp_singlet.dat
!	 analyze.twop_pp_triplet.dat
!        analyze.twop_ph_charge.dat
!        analyze.twop_ph_magnet.dat

!    Other required input files:
!        solver.ctqmc.in
!        parameterset
     program analyze_lavender_soc
     use constants
     use control
     implicit none


!   1. Read in configuration file "solver.ctqmc.in":
     call analyze_config()
!      Read in converged chemical potential from "parameterset",
!      and redefine some variables (ktime,...etc) for initial test
     call analyze_parameterset()
!      Print out the details of input on screen:
     call analyze_print_summary()

!   2. Allocate memory and initialize
     call analyze_setup_array()

!      Read in 1- and 2-ptle functions from files
!        solver.grn.dat
!        solver.sgm.dat
!      and
!        analyze.twop_pp_singlet.dat
!        analyze.twop_ph_charge.dat
!        analyze.twop_ph_magnet.dat
     write(6,*) "2. analyze_readin_array"
     call analyze_readin_array()

!   3. Solve Parquet eqs to obtain fully-irreducible vertex, lambda(w,w'):
!    Fix q=0 and nu=0, and
!    Gamma_s(w,w') = lambda(w,w')+(1/2)(Phid1 + Phid2)-(3/2)(Phim1 + Phim2),
!    where
!    Phi(d/m)1(q=0,nu=0)_p,p' = Phi(d/m)1(w'-w)_-w',w
!    Phi(d/m)2(q=0,nu=0)_p,p' = Phi(d/m)2(w'+w)_-w',-w

!   3.1: Prepare Gamma_s(w,w')
     write(6,*) "3. Solver Parquet eqs for fully-irreducible vertex lambda"
     write(6,*) "3.1 analyze_symmetrize_pairing_Gamma"
     call analyze_symmetrize_pairing_Gamma()
!   3.2: Prepare Gamma_c and Gamma_m
     write(6,*) "3.2 analyze_calculate_Gamma_local for c and m"
     call analyze_calculate_Gamma_local()
!   3.3: Calculate ladder 
!        Phid(w,w') = Sum  Gamma_c * chi_c * Gamma_c
!        Phim(w,w') = Sum  Gamma_m * chi_m * Gamma_m
     write(6,*) "3.3 analyze_calculate_ladder_local for c and m"
     call analyze_calculate_ladder_local()

!   3.4: Calculate lambda(w,w')
     write(6,*) "3.4 analyze_calculate_lambda"
     call analyze_calculate_lambda()
     write(6,*) "Obtain lambda!"

!   4. Solve Parquet eqs second time to obtain momentum dependent Gamma_pp_s:
!   4.1. Calculate local p-p reducible vertex Fc(iw,iw'):
     write(6,*) "4. Solver Parquet eqs pairing singlet irreducible vertex Gamma_pp_s"
!!     write(6,*) "4.1 analyze_calculate_Fc for 3 channels"
!!     call analyze_calculate_Fc()

!!     write(6,*) "4.2 analyze_symmetrize_pairing_Fc for pairing channel"
!!     call analyze_symmetrize_pairing_Fc()

     write(6,*) "4.3 analyze_calculate_Gl_chi0tilde"
     call analyze_calculate_Gl_chi0tilde()

     write(6,*) "4.4 analyze_calculate_vertexladder for c and m"
     call analyze_calculate_vertexladder()
!!    4.5: Deallocate the useless arrays
     write(6,*) "4.5 analyze_final_array_two_mat_local"
     call analyze_final_array_two_mat_local()
     
!   5. Construct irreducible vertex Gamma_pp_s from Parquet equation
     write(6,*) "allocate two mat arrays non local"
     call analyze_setup_array_two_mat_non_local()
     write(6,*) "5. analyze_calculate_Gamma for 1st order"
     call analyze_calculate_Gamma_1st_order()
     
     write(6,*) "5. analyze_calculate_Gamma for 2nd order"
     call analyze_calculate_vertexladder_2nd_order()

!   6. Calculate pairing matrix
     write(6,*) "6. analyze_pairing_matrix"
     call analyze_pairing_matrix()

! deallocate memory
     call analyze_final_array()

     stop
     end program analyze_lavender_soc
