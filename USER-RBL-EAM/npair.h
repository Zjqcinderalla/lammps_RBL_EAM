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

#ifndef LMP_NPAIR_H
#define LMP_NPAIR_H

#include "pointers.h"

namespace LAMMPS_NS {

class NPair : protected Pointers {
 public:
  int istyle;                   // 1-N index into pairnames
  class NBin *nb;               // ptr to NBin instance I depend on
  class NStencil *ns;           // ptr to NStencil instance I depend on
  bigint last_build;            // last timestep build performed

  double cutoff_custom;         // cutoff set by requestor

  NPair(class LAMMPS *);
  virtual ~NPair();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  void build_setup();
  virtual void build(class NeighList *) = 0;

 protected:
  double **mycutneighsq;         // per-type cutoffs when user specified

  // data from Neighbor class

  int includegroup;
  int exclude;
  double skin;
  double **cutneighsq;
  // 
  double **cutcoresq;              // zjqneighbor cutcoresq for each type pair
  // 
  double **cutneighghostsq;
  double cut_inner_sq;
  double cut_middle_sq;
  double cut_middle_inside_sq;
  double *bboxlo,*bboxhi;

  // exclusion data from Neighbor class

  int nex_type;                    // # of entries in type exclusion list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  int nex_mol;                     // # of entries in molecule exclusion list
  int *ex_mol_bit;                 // molecule group bits to exclude
  int *ex_mol_group;               // molecule group #'s to exclude
  int *ex_mol_intra;               // 0 = exclude if in 2 molecules (inter)
                                   // 1 = exclude if in same molecule (intra)

  // special data from Neighbor class

  int *special_flag;

  // data from NBin class

  int nbinx,nbiny,nbinz;
  int mbins;
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  double bininvx,bininvy,bininvz;
  int *atom2bin,*bins;
  int *binhead;

  // data from NStencil class

  int nstencil;
  int *stencil;
  int **stencilxyz;
  int *nstencil_multi;
  int **stencil_multi;
  double **distsq_multi;

  // data common to all NPair variants

  int molecular;

  // methods for all NPair variants

  virtual void copy_bin_info();
  virtual void copy_stencil_info();

  int exclusion(int, int, int,
                int, int *, tagint *) const;   // test for pair exclusion
  int coord2bin(double *);                     // mapping atom coord to a bin
  int coord2bin(double *, int &, int &, int&); // ditto

  // find_special: determine if atom j is in special list of atom i
  // if it is not, return 0
  // if it is and special flag is 0 (both coeffs are 0.0), return -1
  // if it is and special flag is 1 (both coeffs are 1.0), return 0
  // if it is and special flag is 2 (otherwise), return 1,2,3
  //   for which level of neighbor it is (and which coeff it maps to)

  inline int find_special(const tagint *list, const int *nspecial,
                          const tagint tag) const {
    const int n1 = nspecial[0];
    const int n2 = nspecial[1];
    const int n3 = nspecial[2];

    for (int i = 0; i < n3; i++) {
      if (list[i] == tag) {
        if (i < n1) {
          if (special_flag[1] == 0) return -1;
          else if (special_flag[1] == 1) return 0;
          else return 1;
        } else if (i < n2) {
          if (special_flag[2] == 0) return -1;
          else if (special_flag[2] == 1) return 0;
          else return 2;
        } else {
          if (special_flag[3] == 0) return -1;
          else if (special_flag[3] == 1) return 0;
          else return 3;
        }
      }
    }
    return 0;
  };

  int copymode;
  ExecutionSpace execution_space;
};

}

#endif

/* ERROR/WARNING messages:

E: Non-numeric positions - simulation unstable

UNDOCUMENTED

*/
