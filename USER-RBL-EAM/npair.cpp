/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair.h"
#include <cmath>
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "nstencil.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPair::NPair(LAMMPS *lmp)
  : Pointers(lmp), nb(NULL), ns(NULL), bins(NULL), stencil(NULL)
{
  last_build = -1;
  mycutneighsq = NULL;
  molecular = atom->molecular;
  copymode = 0;
  execution_space = Host;
}

/* ---------------------------------------------------------------------- */

NPair::~NPair()
{
  if (copymode) return;

  memory->destroy(mycutneighsq);
}

/* ---------------------------------------------------------------------- */

void NPair::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   done once per run
------------------------------------------------------------------------- */

void NPair::copy_neighbor_info()
{
  // general params

  includegroup = neighbor->includegroup;
  exclude = neighbor->exclude;
  skin = neighbor->skin;
  // 
  cutneighsq = neighbor->cutneighsq;
  cutcoresq = neighbor->cutcoresq;
  // 
  cutneighghostsq = neighbor->cutneighghostsq;
  cut_inner_sq = neighbor->cut_inner_sq;
  cut_middle_sq = neighbor->cut_middle_sq;
  cut_middle_inside_sq = neighbor->cut_middle_inside_sq;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // exclusion info

  nex_type = neighbor->nex_type;
  ex1_type = neighbor->ex1_type;
  ex2_type = neighbor->ex2_type;
  ex_type = neighbor->ex_type;

  nex_group = neighbor->nex_group;
  ex1_group = neighbor->ex1_group;
  ex2_group = neighbor->ex2_group;
  ex1_bit = neighbor->ex1_bit;
  ex2_bit = neighbor->ex2_bit;

  nex_mol = neighbor->nex_mol;
  ex_mol_group = neighbor->ex_mol_group;
  ex_mol_bit = neighbor->ex_mol_bit;
  ex_mol_intra = neighbor->ex_mol_intra;

  // special info

  special_flag = neighbor->special_flag;

  // overwrite per-type Neighbor cutoffs with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {
    memory->destroy(mycutneighsq);
    int n = atom->ntypes;
    memory->create(mycutneighsq,n+1,n+1,"npair:cutneighsq");
    int i,j;
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        mycutneighsq[i][j] = cutoff_custom * cutoff_custom;
    cutneighsq = mycutneighsq;
  }
}

/* ----------------------------------------------------------------------
   copy info from NBin class to this build class
------------------------------------------------------------------------- */

void NPair::copy_bin_info()
{
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;

  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;

  atom2bin = nb->atom2bin;
  bins = nb->bins;
  binhead = nb->binhead;
}

/* ----------------------------------------------------------------------
   copy info from NStencil class to this build class
------------------------------------------------------------------------- */

void NPair::copy_stencil_info()
{
  nstencil = ns->nstencil;
  stencil = ns->stencil;
  stencilxyz = ns->stencilxyz;
  nstencil_multi = ns->nstencil_multi;
  stencil_multi = ns->stencil_multi;
  distsq_multi = ns->distsq_multi;
}

/* ----------------------------------------------------------------------
   copy info from NBin and NStencil classes to this build class
------------------------------------------------------------------------- */

void NPair::build_setup()
{
  if (nb) copy_bin_info();
  if (ns) copy_stencil_info();

  // set here, since build_setup() always called before build()

  last_build = update->ntimestep;
}

/* ----------------------------------------------------------------------
   test if atom pair i,j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
------------------------------------------------------------------------- */

int NPair::exclusion(int i, int j, int itype, int jtype,
                     int *mask, tagint *molecule) const {
  int m;

  if (nex_type && ex_type[itype][jtype]) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask[i] & ex1_bit[m] && mask[j] & ex2_bit[m]) return 1;
      if (mask[i] & ex2_bit[m] && mask[j] & ex1_bit[m]) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)

      // intra-chain: exclude i-j pair if in same molecule
      // inter-chain: exclude i-j pair if in different molecules

      if (ex_mol_intra[m]) {
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] == molecule[j]) return 1;
      } else {
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] != molecule[j]) return 1;
      }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   same as coord2bin in Nbin, but also return ix,iy,iz offsets in each dim
   used by some of the ghost neighbor lists
------------------------------------------------------------------------- */

int NPair::coord2bin(double *x, int &ix, int &iy, int &iz)
{
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  return iz*mbiny*mbinx + iy*mbinx + ix;
}

