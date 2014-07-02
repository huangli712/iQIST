#!/usr/bin/env python3

# define the standard parameters set by a dictionary data structure
# item structure
# 0: full definition, string
# 1: short definition, string
# 2: default value, number
# 3: data type, string
para={}
para["_isscf "] = [
"""
   ! control flag: running mode
   ! if isscf == 1, one-shot non-self-consistent scheme, used in local density
   ! approximation plus dynamical mean field theory case
   ! if isscf == 2, self-consistent scheme, used in normal model hamiltonian
   ! plus dynamical mean field theory case
""",
"""
   ! non-self-consistent (1) or self-consistent mode (2)
""",
   2, "integer" ]

para["_issun "] = [
"""
   ! control flag: symmetry of bands
   ! if issun == 1, the bands are not symmetrized
   ! if issun == 2, the bands are symmetrized according to symmetry matrix
""",
"""
   ! without symmetry    (1) or with symmetry   mode (2)
""",
   2, "integer" ]

para["_isspn "] = [
"""
   ! control flag: symmetry of spin orientation
   ! if isspn == 1, enforce spin up = spin down
   ! if isspn == 2, let spin up and spin down states evolve independently
""",
"""
   ! spin projection, PM (1) or AFM             mode (2)
""",
   1, "integer" ]

para["_isbin "] = [
"""
   ! control flag: impurity green's function binning mode
   ! if isbin == 1, without binning mode
   ! if isbin == 2, with binning mode
""",
"""
   ! without binning     (1) or with binning    mode (2)
""",
   2, "integer" ]

para["_isort "] = [
"""
   ! control flag: apply orthogonal polynomial representation to perform measurement
   ! if isort == 1, use normal representation to measure G(\tau)
   ! if isort == 2, use legendre polynomial to measure G(\tau)
   ! if isort == 3, use chebyshev polynomial (the second kind) to measure G(\tau)
   ! if isort == 4, use normal representation to measure G(\tau) and F(\tau)
   ! if isort == 5, use legendre polynomial to measure G(\tau) and F(\tau)
   ! if isort == 6, use chebyshev polynomial (the second kind) to measure G(\tau) and F(\tau)
   ! note: if isort \in [1,3], we use ctqmc_make_hub1() to calculate the self
   ! energy function, or else we use ctqmc_make_hub2().
   ! note: as for the kernel polynomial representation, the default dirichlet
   ! kernel is applied automatically. if you want to choose the other kernel,
   ! please check the ctqmc_make_gtau() subroutine.
""",
"""
   ! normal measurement  (1) or legendre polynomial  (2) or chebyshev polynomial (3)
""",
   1, "integer" ]

para["_isvrt "] = [
"""
   ! control flag: whether we measure the high order correlation function
   ! if isvrt == 1, do nothing
   ! if isvrt == 2, calculate spin-spin correlation function
   ! if isvrt == 3, calculate orbital-orbital correlation function
   ! if isvrt == 4, calculate both two-particle green's function and vertex function
   ! if isvrt == 5, calculate both two-particle green's function and vertex function
   ! note: when isvrt == 4 and isvrt == 5, both the two-particle green's and
   ! vertex functions are computed by using two different algorithms.
""",
"""
   ! without vertex      (1) or with vertex function (2)
""",
   1, "integer" ]

para["_isscr "] = [
"""
   ! control flag: model need to be solved
   ! if isscr == 1, normal model
   ! if isscr == 2, holstein-hubbard model
   ! if isscr == 3, dynamic screening, palsmon pole model
   ! if isscr == 4, dynamic screening, ohmic model
   ! if isscr ==99, dynamic screening, realistic materials
   ! note: when isscr == 1, lc and wc are ignored. when isscr == 2, lc means
   ! the electron-phonon coupling constant \lambda, and wc phonon vibration
   ! frequency \omega_{0}. when isscr == 3, lc and wc just mean the control
   ! parameters \lambda and \omega^{'}, respectively. when isscr == 4, lc
   ! and wc just mean the control parameters \alpha and \omega_{c}, respectively.
   ! when isscr == 99, wc is ignored and lc means the shift for interaction
   ! matrix and chemical potential.
""",
"""
   ! normal (1) or holstein-hubbard (2) or plasmon pole (3) or ohmic model (4)
""",
   1, "integer" ]

para["_nband "] = [
"""
   ! number of correlated bands
""",
"""
   ! number of correlated bands
""",
   1, "integer" ]

para["_nspin "] = [
"""
   ! number of spin projection
""",
"""
   ! number of spin projection
""",
   2, "integer" ]

para["_norbs "] = [
"""
   ! number of correlated orbitals (= nband * nspin)
""",
"""
   ! number of correlated orbitals (= nband * nspin)
""",
   2, "integer" ]

para["_ncfgs "] = [
"""
   ! number of atomic states (= 2**norbs)
""",
"""
   ! number of atomic states
""",
   4, "integer" ]

para["_nzero "] = [
"""
   ! maximum allowed number of non-zero elements in F-matrix
""",
"""
   ! maximum number of non-zero elements in sparse matrix style
""",
   128, "integer" ]

para["_nvect "] = [
"""
   ! number of selected eigenvectors (maximum value is ncfgs, minimum value is 1)
""",
"""
   ! number of selected eigenvectors
""",
   4, "integer" ]

para["_nhmat "] = [
"""
   ! maximum allowed number of non-zero elements in H-matrix
""",
"""
   ! maximum number of non-zero elements of H in sparse matrix style
""",
   128, "integer" ]

para["_nfmat "] = [
"""
   ! maximum allowed number of non-zero elements in F-matrix
""",
"""
   ! maximum number of non-zero elements of F in sparse matrix style
""",
   128, "integer" ]

para["_niter "] = [
"""
   ! maximum number of continuous time quantum Monte Carlo quantum impurity
   ! solver plus dynamical mean field theory self-consistent iterations
""",
"""
   ! maximum number of DMFT + CTQMC self-consistent iterations
""",
   20, "integer" ]

para["_U     "] = [
"""
   ! average Coulomb interaction
""",
"""
   ! U : average Coulomb interaction
""",
   4.00, "float" ]

para["_Uc    "] = [
"""
   ! intraorbital Coulomb interaction
""",
"""
   ! Uc: intraorbital Coulomb interaction
""",
   4.00, "float" ]

para["_Uv    "] = [
"""
   ! interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
""",
"""
   ! Uv: interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
""",
   4.00, "float" ]

para["_Jz    "] = [
"""
   ! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
""",
"""
   ! Jz: Hund's exchange interaction in z axis (Jz = Js = Jp = J)
""",
   0.00, "float" ]

para["_Js    "] = [
"""
   ! spin-flip term
""",
"""
   ! Js: spin-flip term
""",
   0.00, "float" ]

para["_Jp    "] = [
"""
   ! pair-hopping term
""",
"""
   ! Jp: pair-hopping term
""",
   0.00, "float" ]

para["_lc    "] = [
"""
   ! strength of dynamical screening effect ( or electron-phonon coupling )
""",
"""
   ! lc: screening strength
""",
   1.00, "float" ]

para["_wc    "] = [
"""
   ! screening frequency ( or frequency for einstein phonons )
""",
"""
   ! wc: screening frequency
""",
   1.00, "float" ]

para["_mune  "] = [
"""
   ! chemical potential or fermi level
   ! note: it should be replaced with eimp
""",
"""
   ! chemical potential or fermi level
""",
   2.00, "float" ]

para["_beta  "] = [
"""
   ! inversion of temperature
""",
"""
   ! inversion of temperature
""",
   8.00, "float" ]

para["_part  "] = [
"""
   ! coupling parameter t for Hubbard model
""",
"""
   ! coupling parameter t for Hubbard model
""",
   0.50, "float" ]

para["_alpha "] = [
"""
   ! mixing parameter for dynamical mean field theory self-consistent engine
""",
"""
   ! mixing parameter for self-consistent engine
""",
   0.70, "float" ]

para["_lemax "] = [
"""
   ! maximum order for legendre polynomial
""",
"""
   ! maximum order for legendre polynomial
""",
   32, "integer" ]

para["_legrd "] = [
"""
   ! number of mesh points for legendre polynomial in [-1,1] range
""",
"""
   ! number of mesh points for legendre polynomial
""",
   20001, "integer" ]

para["_chmax "] = [
"""
   ! maximum order for chebyshev polynomial
""",
"""
   ! maximum order for chebyshev polynomial
""",
   32, "integer" ]

para["_chgrd "] = [
"""
   ! number of mesh points for chebyshev polynomial in [-1,1] range
""",
"""
   ! number of mesh points for chebyshev polynomial
""",
   20001, "integer" ]

para["_mkink "] = [
"""
   ! maximum perturbation expansion order
""",
"""
   ! maximum perturbation expansions order
""",
   1024, "integer" ]

para["_mfreq "] = [
"""
   ! maximum number of matsubara frequency point
""",
"""
   ! maximum number of matsubara frequency
""",
   8193, "integer" ]

para["_nffrq "] = [
"""
   ! number of matsubara frequency for the two-particle green's function
""",
"""
   ! number of matsubara frequency for the two-particle green's function
""",
   32, "integer" ]

para["_nbfrq "] = [
"""
   ! number of bosonic frequncy for the two-particle green's function
""",
"""
   ! number of bosonic frequncy for the two-particle green's function
""",
   8, "integer" ]

para["_nfreq "] = [
"""
   ! number of matsubara frequency sampling by continuous time quantum Monte
   ! Carlo quantum impurity solver
   ! note: the rest (mfreq - nfreq + 1 points) values are evaluated by using
   ! Hubbard-I approximation
""",
"""
   ! maximum number of matsubara frequency sampling by quantum impurity solver
""",
   128, "integer" ]

para["_ntime "] = [
"""
   ! number of imaginary time slice sampling by continuous time quantum Monte
   ! Carlo quantum impurity solver
""",
"""
   ! number of time slice
""",
   1024, "integer" ]

para["_nleja "] = [
"""
   ! maximum number of real leja points. the real leja points are used in the
   ! newton interpolation to evaluate exp( -\tau H ) |v> efficiently
""",
"""
   ! maximum number of real leja points
""",
   64, "integer" ]

para["_npart "] = [
"""
   ! number of parts that the imaginary time axis is split
   ! note: all operators in the imaginary time axis are grouped into npart
   ! parts according to their time values, in each Monte Carlo steps, only
   ! those changed parts are carefully dealt with, not all the parts.
   ! note: 2\sqrt{3 <k> nband} ~ 4\sqrt{3 <k> nband} may be the optimal value
   ! for npart to achieve maximum performance
""",
"""
   ! number of parts that the imaginary time axis is split
""",
   16, "integer" ]

para["_nflip "] = [
"""
   ! flip period for spin up and spin down states
   ! note: care must be taken to prevent the system from being trapped in a
   ! state which breaks a symmetry of local hamiltonian when it should not
   ! be. to avoid unphysical trapping, we introduce "flip" moves, which
   ! exchange the operators corresponding, for example, to up and down spins
   ! in a given orbital.
   ! note: in this code, nowadays the following flip schemes are supported
   !     if cflip = 1, flip inter-orbital spins randomly;
   !     if cflip = 2, flip intra-orbital spins one by one;
   !     if cflip = 3, flip intra-orbital spins globally.
   ! note: we use the sign of nflip to control flip schemes
   !     if nflip = 0, means infinite long period to do flip
   !     if nflip > 0, combine cflip = 2 (80%) and cflip = 3 (20%)
   !     if nflip < 0, combine cflip = 1 (80%) and cflip = 3 (20%)
   ! note: if nflip /= 0, the absolute value of nflip is the flip period
   ! note: when cflip = 1, the symmetry of all orbitals must be taken into
   ! consideration, otherwise the code may be trapped by a deadlock.
""",
"""
   ! flip period for spin up and spin down states
""",
   20000, "integer" ]

para["_ntherm"] = [
"""
   ! maximum number of thermalization steps
""",
"""
   ! maximum number of thermalization steps
""",
   200000, "integer" ]

para["_nsweep"] = [
"""
   ! maximum number of quantum Monte Carlo sampling steps
""",
"""
   ! maximum number of quantum Monte Carlo sampling steps
""",
   20000000, "integer" ]

para["_nwrite"] = [
"""
   ! output period for quantum impurity solver
""",
"""
   ! output period
""",
   2000000, "integer" ]

para["_nclean"] = [
"""
   ! clean update period for quantum impurity solver
""",
"""
   ! clean update period
""",
   100000, "integer" ]

para["_nmonte"] = [
"""
   ! how often to sampling the gmat and paux (nmat and nnmat)
   ! note: the measure periods for schi, sschi, ochi, oochi, g2_re, g2_im,
   ! h2_re, and h2_im are also controlled by nmonte parameter.
""",
"""
   ! how often to sampling the gmat and nmat
""",
   10, "integer" ]

para["_ncarlo"] = [
"""
   ! how often to sampling the gtau and prob
   ! note: the measure period for ftau is also controlled by ncarlo parameter.
""",
"""
   ! how often to sampling the gtau and prob
""",
   10, "integer" ]

# define the required parameters set for azalea code
azalea=[
   "=======",
   "CCCCCCC",
   "=======",
   "_isscf ",
   "_issun ",
   "_isspn ",
   "_isbin ",
   "-------",
   "_nband ",
   "_nspin ",
   "_norbs ",
   "_ncfgs ",
   "_niter ",
   "-------",
   "_U     ",
   "_Uc    ",
   "_Uv    ",
   "_Jz    ",
   "_Js    ",
   "_Jp    ",
   "-------",
   "_mune  ",
   "_beta  ",
   "_part  ",
   "_alpha ",
   "^^^^^^^",
   "_mkink ",
   "_mfreq ",
   "-------",
   "_nfreq ",
   "_ntime ",
   "_nflip ",
   "_ntherm",
   "_nsweep",
   "_nwrite",
   "_nclean",
   "_nmonte",
   "_ncarlo",
   "^^^^^^^"]

# define the required parameters set for gardenia code
gardenia=[
   "=======",
   "CCCCCCC",
   "=======",
   "_isscf ",
   "_issun ",
   "_isspn ",
   "_isbin ",
   "_isort ",
   "_isvrt ",
   "-------",
   "_nband ",
   "_nspin ",
   "_norbs ",
   "_ncfgs ",
   "_niter ",
   "-------",
   "_U     ",
   "_Uc    ",
   "_Uv    ",
   "_Jz    ",
   "_Js    ",
   "_Jp    ",
   "-------",
   "_mune  ",
   "_beta  ",
   "_part  ",
   "_alpha ",
   "^^^^^^^",
   "_lemax ",
   "_legrd ",
   "_chmax ",
   "_chgrd ",
   "-------",
   "_mkink ",
   "_mfreq ",
   "-------",
   "_nffrq ",
   "_nbfrq ",
   "_nfreq ",
   "_ntime ",
   "_nflip ",
   "_ntherm",
   "_nsweep",
   "_nwrite",
   "_nclean",
   "_nmonte",
   "_ncarlo",
   "^^^^^^^"]

# define the required parameters set for narcissus code
narcissus=[
   "=======",
   "CCCCCCC",
   "=======",
   "_isscf ",
   "_issun ",
   "_isspn ",
   "_isbin ",
   "_isort ",
   "_isvrt ",
   "_isscr ",
   "-------",
   "_nband ",
   "_nspin ",
   "_norbs ",
   "_ncfgs ",
   "_niter ",
   "-------",
   "_U     ",
   "_Uc    ",
   "_Uv    ",
   "_Jz    ",
   "_Js    ",
   "_Jp    ",
   "_lc    ",
   "_wc    ",
   "-------",
   "_mune  ",
   "_beta  ",
   "_part  ",
   "_alpha ",
   "^^^^^^^",
   "_lemax ",
   "_legrd ",
   "_chmax ",
   "_chgrd ",
   "-------",
   "_mkink ",
   "_mfreq ",
   "-------",
   "_nffrq ",
   "_nbfrq ",
   "_nfreq ",
   "_ntime ",
   "_nflip ",
   "_ntherm",
   "_nsweep",
   "_nwrite",
   "_nclean",
   "_nmonte",
   "_ncarlo",
   "^^^^^^^"]

# define escape color characters
BLACK='\033[30;1m'
RED='\033[31;1m'
GREEN='\033[32;1m'
YELLOW='\033[33;1m'
BLUE='\033[34;1m'
MAGENTA='\033[35;1m'
CYAN='\033[36;1m'
WHITE='\033[37;1m'
RESET='\033[0m'

# print the header
print(GREEN + "  pysci.py v0.01" + RESET)
print(GREEN + "  >>> A Flexible And Powerful SOLVER-CTQMC-IN Builder" + RESET)

print()
print()
print(RED + "  To build your own solver.ctqmc.in file, you need to follow the instructions and comments" + RESET)
print(RED + "  below, and answer a few questions carefully. The solver.ctqmc.in file generated by this " + RESET)
print(RED + "  script is compatible with iQIST v0.10 only. Good luck to you!" + RESET)

# Step A: select suitable quantum impurity solver
print()
print(GREEN + "  Step A:" + RESET)
print(GREEN + "  Please choose your favourite quantum impurity solver. Only the following three ctqmc" + RESET)
print(GREEN + "  impurity solvers are supported so far. Input the first letter of ctqmc program name " + RESET)
print(GREEN + "  in uppercase and then press the ENTER key." + RESET)
print("     A: azalea")
print("     G: gardenia")
print("     N: narcissus")

my_ctqmc = 'A'
user_ctqmc = input("  >>> enter your choice <<< ") or my_ctqmc

exe_ctqmc = []
str_ctqmc = ""
if   user_ctqmc == 'A' or user_ctqmc == 'a':
   exe_ctqmc = azalea
   str_ctqmc = "AZALEA: "
   print(RED + "  status: " + RESET + GREEN + "AZALEA is chosen." + RESET)
elif user_ctqmc == 'G' or user_ctqmc == 'g':
   exe_ctqmc = gardenia
   str_ctqmc = "GARDENIA: "
   print(RED + "  status: " + RESET + GREEN + "GARDENIA is chosen." + RESET)
elif user_ctqmc == 'N' or user_ctqmc == 'n':
   exe_ctqmc = narcissus
   str_ctqmc = "NARCISSUS: "
   print(RED + "  status: " + RESET + GREEN + "NARCISSUS is chosen." + RESET)
else:
   exe_ctqmc = azalea
   str_ctqmc = "AZALEA: "
   print(RED + "  status: " + RESET + GREEN + "AZALEA is chosen." + RESET)

# Step B: setup configuration parameters
print()
print(GREEN + "  Step B:" + RESET)
print(GREEN + "  Please setup the parameters for selected quantum impurity solver. You need to read" + RESET)
print(GREEN + "  the datatype, default value, and comment for the specific parameter carefully, and" + RESET)
print(GREEN + "  input your choose. If you don't have any idea, you can just press the ENTER key.  " + RESET)
print()

for p in exe_ctqmc:
   if p.startswith('_'):
      str_name = RED + "  >>> " + RESET + GREEN + p + RESET
      str_type = RED + "  type: " + RESET + GREEN + para[p][3] + RESET
      str_dflt = RED + "  default: " + RESET + GREEN + str(para[p][2]) + RESET
      print(str_name + str_type + str_dflt) # output the name, type, and default value
      print("  " + "-"*75)
      print(para[p][0]) # output the full definition
      my_input = para[p][2] # setup the default value
      user_input = input("  >>> enter your choice <<< ") or my_input # get the user input
      para[p][2] = user_input # overwrite the default value by using the user input
      print(RED + "  value: " + RESET + GREEN + str(para[p][2]) + RESET)
      print()
   else:
      continue

# Step C: output the solver.ctqmc.in file
print()
print(GREEN + "  Step C:" + RESET)
print(GREEN + "  Build the solver.ctqmc.in file for selected quantum impurity solver" + RESET)
print()

sci = open('solver.ctqmc.in','w') # open solver.ctqmc.in file
for p in exe_ctqmc:
   if p.startswith('='):
      sci.write("="*75 +"\n")
   if p.startswith('-'):
      sci.write("-"*75 +"\n")
   if p.startswith('^'):
      sci.write("^"*75 +"\n")
   if p.startswith('C'):
      sci.write(str_ctqmc + "continuous time quantum Monte Carlo quantum impurity solver" +"\n")
   if p.startswith('_'):
      value = str(para[p][2])
      empty = (13 - len(value))*" "
      sci.write(value + empty + (para[p][1]).strip() + "\n")
sci.close() # close solver.ctqmc.in file

# print the footer
print()
print("  Everything is OK! Then enjoy your ctqmc code.")
