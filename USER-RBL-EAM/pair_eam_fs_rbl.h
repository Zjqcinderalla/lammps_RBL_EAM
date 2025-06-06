/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors: Zhang Jieqiong(NWU)
------------------------------------------------------------------------- */


#ifdef PAIR_CLASS

PairStyle(eam/fs/rbl,PairEAMFSRbl)

#else

#ifndef LMP_PAIR_EAM_FS_RBL_H
#define LMP_PAIR_EAM_FS_RBL_H

#include "pair_eam_rbl.h"

namespace LAMMPS_NS {

// need virtual public b/c of how eam/fs/opt inherits from it

class PairEAMFSRbl : virtual public PairEAMRbl {
 public:
  PairEAMFSRbl(class LAMMPS *);
  virtual ~PairEAMFSRbl() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: No matching element in EAM potential file

The EAM potential file does not contain elements that match the
requested elements.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect element names in EAM potential file

The element names in the EAM file do not match those requested.

E: Invalid EAM potential file

UNDOCUMENTED

*/
