# Chromatin Code
To simulate the genome organization at a resolution of ~100 kbp, along with phase separation of active and silent regions of the genome using first principles. We can model the chromatin chain (nucleosomes with nuclear enzymes) as a simplistic microscopic mechanical description in the form of a polymer chain with beads connected entropic springs. To differentiate, structurally and functionally heterochromatin and euchromatin domains of the genome, we can add to model the specific binding interaction of heterochromatin protein 1 (HP1) with H3K9me3 as crosslinks betweek heterochromatin domains. This code performs Brownian dynamics simulations of bead-spring polymer (chromatin) in the presence of a rigid confinement (Sphere/Elliposid at present). The input parameters needs to be in the Variables.f90 file with the following set of parameters (refer to the paper for more details on the model: https://www.pnas.org/content/115/45/11442.short).
N: Number of beads/chain needed to be simulated
Nc: Number of chains in the system
Ntr: Number of brownian tracer particles (can be either the)

Rs = effective radius of the confinement
aoc, boc = ratio of a/c and b/c (= 1 for a sphere)
k_on/k_off = on and off rate of the molecular motors

# Prerequisites
gfortran; OpemMPI, and STKFMM as a dependency: https://github.com/wenyan4work/STKFMM. Install STKFMM and compile test codes. The instructions are available on the github page (feel free to reach out if you have any questions). I am also planning to add a simple version of the code to evaluate surface traction forces on an arbitraty domain.

# Function/Subroutines (source files)

Each subroutine/functions performs specific functions. Here are the details of it:
1. ActiveForces.f90: This function evaluates active force dipoles on individual beads (+ve for extensile dipoles and -ve for contractile). The dipoles are stochastic in nature, and their binding/unbinding are dependens on k_{off} and k_{on} (parameters in Variables.f90). 
2. Cholesky.f90: This function evaluates the inverse of a symmetric positive-definite matrix. At present, it's not required in the calculations, but could be used for problems where the Brownian forces drives the flow and are not local.
3. Ext_mod2a.f90: This is main source file (still don't remember the origin of the name). This is where all the functions are used, and are tied in. All the local variables are defined here.
4. Geo_Val.f90: For calculating the surface normal on each surface grid. Again, this is not required for this problem specifically as we are dealing with a single layer potential, and the surface is rigid, but for problems where the surface deforms, we might need to evaluate the surface normal/tangent in order to dynamically move the surface.
5. Green_fs.f90: This function evaluates the surface integral of any vector field on the surface. 
6. InitConfig/Initdip.f90: Initilaize the chain inside the bounding domain. Here, I have initialized the chain as a random walk configuration on the surface of the sphere. I am also planning to add the part where the chain can be initialized as a fractal globule.
7. Interp.f90: This function interpolates the vector field on the right triangle in (eta-zeta) plane where the surface integrals are performed.
8. Leonard_jones.f90: This function evaluates the repulsive forces used to avoid crossing of the chain. Again the name is misleading because it's a generic third power repulsive potential and has nothing with Leonard Jones potential. 
9. gasdev.90: Random number generator for Brownian forces.
10. sgf_3d_fs.f90: This function calculate the greens functions (Stokeslet and Stresslet as well for the pressure).
# Compilation
Compiling the code is a straigt forward process: make


