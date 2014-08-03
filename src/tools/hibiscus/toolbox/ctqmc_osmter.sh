#!/bin/bash

#-------------------------------------------------------------------------
# ctqmc_osmter script is used to generate the necessary input files for
# orbital-selective mott transition (OSMT) calculations.
# it is written in Bash shell.
#
# author  : li huang (email:huangli712@yahoo.com.cn)
# version : 0.0.1
# status  : unstable
#-------------------------------------------------------------------------

# loop over Coulomb interaction U
for i in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0
do

# calculate Uv, Jz, Js, and Jp
Uv=`echo "scale=4; $i / 2.0" | bc`
Jz=`echo "scale=4; $i / 4.0" | bc`
Js=`echo "scale=4; $i / 4.0" | bc`
Jp=`echo "scale=4; $i / 4.0" | bc`

# setup eimp and sign
eimp=0.0
sign=1

# print parameters to the screen
echo "U= $i  Uv= $Uv  Jz= $Jz  Js= $Js  Jp= $Jp  Eimp=$eimp"

# dump the input file (solver.ctqmc.in) for ctqmc quantum impurity solver
cat > solver.ctqmc.in <<!
==========================================================================
    setup continuous time quantum Monte Carlo quantum impurity solver
==========================================================================
2            ! non-self-consistent (1) or self-consistent mode (2)
2            ! without symmetry    (1) or with symmetry   mode (2)
1            ! spin projection, PM (1) or AFM             mode (2)
1            ! without binning     (1) or with binning    mode (2)
--------------------------------------------------------------------------
3            ! number of bands
2            ! number of spin projection
6            ! number of orbitals (= nband * nspin)
64           ! number of atomic states
1024         ! maximum number of non-zero elements in sparse matrix style
40           ! maximum number of DMFT + CTQMC self-consistent iterations
--------------------------------------------------------------------------
$i           ! U : average Coulomb interaction
$i           ! Uc: intraorbital Coulomb interaction
$Uv          ! Uv: interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
$Jz          ! Jz: Hund's exchange interaction in z axis (Jz = Js = Jp = J)
$Js          ! Js: spin-flip term
$Jp          ! Jp: pair-hopping term
--------------------------------------------------------------------------
8.00         ! chemical potential or fermi level
100.         ! inversion of temperature
0.50         ! coupling parameter t for Hubbard model
0.70         ! mixing parameter for self-consistent engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1024         ! maximum perturbation expansions order
8193         ! maximum number of matsubara frequency
--------------------------------------------------------------------------
128          ! maximum number of matsubara frequency sampling by quantum impurity solver
1024         ! number of time slice
4            ! number of parts that the imaginary time axis is split
20000        ! flip period for spin up and spin down states
200000       ! maximum number of thermalization steps
40000000     ! maximum number of quantum Monte Carlo sampling steps
4000000      ! output period
100000       ! clean update period
100          ! how often to sampling the gmat and nmat
100          ! how often to sampling the gtau and prob
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!

# dump the input file (solver.eimp.in) for ctqmc quantum impurity solver
cat > solver.eimp.in <<!
1  $eimp 1
2  0.0   2
3  0.0   2
4  $eimp 1
5  0.0   2
6  0.0   2
!

# dump the input file (atom.in) for atomic eigenvalues problem solver
cat > atom.in <<!
--------------------------------------------------------------------------
>>>     atom.in: configuration file for atomic eigenvalues problem     <<<
--------------------------------------------------------------------------
1        ! imode: running mode, 1 = eigen basis mode; 2 = occupation number basis mode
1        ! ifock: source of fock space, 1 = internal mode; 2 = external mode
-1       ! isoce: with or without SOC, -1 = without; 0 = s; 1 = p; 2 = d; 3 = f; 99 = debug mode
--------------------------------------------------------------------------
3        ! nband: number of correlated bands
2        ! nspin: number of spin projection
6        ! norbs: number of correlated orbitals (= nband * nspin)
64       ! ncfgs: number of atomic states (= 2**norbs)
--------------------------------------------------------------------------
$i       ! U    : average Coulomb interaction
$i       ! Uc   : intraorbital Coulomb interaction
$Uv      ! Uv   : interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
$Jz      ! Jz   : Hund's exchange interaction in z axis (Jz = Js = Jp = J)
$Js      ! Js   : spin-flip term
$Jp      ! Jp   : pair-hopping term
0.00     ! lsoc : spin orbital coupling strength
--------------------------------------------------------------------------
$eimp    ! eimp(01) orbital, impurity level
 0.00000 ! eimp(02) orbital
 0.00000 ! eimp(03) orbital
$eimp    ! eimp(04) orbital
 0.00000 ! eimp(05) orbital
 0.00000 ! eimp(06) orbital
 0.00000 ! eimp(07) orbital
 0.00000 ! eimp(08) orbital
 0.00000 ! eimp(09) orbital
 0.00000 ! eimp(10) orbital
 0.00000 ! eimp(11) orbital
 0.00000 ! eimp(12) orbital
 0.00000 ! eimp(13) orbital
 0.00000 ! eimp(14) orbital
--------------------------------------------------------------------------
!

# execute atomic solver to generate atom.cix, which is necessary for ctqmc
# quantum impurity solver
./atom > atom.out

# generate correct directory, and prepare necessary files
if [ $sign -ge 0 ]
then
    mkdir p_d$eimp\_u$i
    mv solver.ctqmc.in p_d$eimp\_u$i
    mv solver.eimp.in  p_d$eimp\_u$i
    mv atom.cix        p_d$eimp\_u$i
    mv atom.in         p_d$eimp\_u$i
    mv atom.out        p_d$eimp\_u$i
else
    eimp=`echo ${eimp#-}`
    mkdir m_d$eimp\_u$i
    mv solver.ctqmc.in m_d$eimp\_u$i
    mv solver.eimp.in  m_d$eimp\_u$i
    mv atom.cix        m_d$eimp\_u$i
    mv atom.in         m_d$eimp\_u$i
    mv atom.out        m_d$eimp\_u$i
fi

# remove the temporary files in working directory
rm -f *.dat

done
