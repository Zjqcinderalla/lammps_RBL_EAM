# lammps_RBL_EAM
Random batch method for molecular dynamics with EAM potential

# Usage Tutorial:
**LAMMPS Version**: lammps-7Aug19 or others version 

1. Download "USER-RBL-EAM" to the `src` folder in the LAMMPS source code:  
   ```bash
   cd src
   make yes-USER-RBL-EAM
   make mpi

2. The Examples-RBL folder provides test files. The RBL-EAM program part can be implemented using the following LAMMPS commands:
   ```bash
   # pair_style types cut_core_off  batch_size
   pair_style          eam/fs/rbl   2.6   10
   pair_coeff          *   *      Fe_mm.eam.fs  Fe


3. The Examples-lattice-constants folder provides a program for calculating the lattice constants of metal Fe based on the RBL-EAM algorithm.
