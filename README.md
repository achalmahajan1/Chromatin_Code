# Chromatin Code
This code performs Brownian dynamics simulations of bead-spring polymer (chromatin) in the presence of a rigid confinement (Sphere/Elliposid at present). The input parameters needs to be in the Variables.f90 file with the following set of parameters
N: Number of beads/chain needed to be simulated
Nc: Number of chains in the system
Ntr: Number of brownian tracer particles (can be either the)

Rs = effective radius of the confinement
aoc, boc = ratio of a/c and b/c (= 1 for a sphere)
k_on/k_off = on and off rate of the molecular motors

# Prerequisites
gfortran; OpemMPI, and STKFMM as a dependency: https://github.com/wenyan4work/STKFMM. Install STKFMM and compile test codes. The instructions are available on the github page (feel free to reach out if you have any questions). I am also planning to add a simple version of the code to evaluate surface traction forces on an arbitraty domain.

# Compilation
Compiling the code is a straigt forward process: make


