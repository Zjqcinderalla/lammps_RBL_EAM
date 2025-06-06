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

#include "npair_half_bin_atomonly_newton.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS; 
/* ---------------------------------------------------------------------- */

NPairHalfBinAtomonlyNewton::NPairHalfBinAtomonlyNewton(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfBinAtomonlyNewton::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  //-------------------zjq--------------------//
  int *neighptr_core, *neighptr_shell;
  int n1,n2,rbl_local;
  
  rbl_local = list->neigh_rbl;   // read the rbl value and start random batch calculation if rbl =1
  
  int *numneigh_core, **firstneigh_core, *numneigh_shell, **firstneigh_shell;
  MyPage<int> *ipage_core, *ipage_shell;
  // if the random batch list is needed, i.e., rbl = 1;
  if (rbl_local){
    // printf("npair_half_bin_atomonly_newton-rbl:%d \n", rbl_local);
    numneigh_core = list->numneigh_core;
    firstneigh_core = list->firstneigh_core;
    ipage_core = list->ipage_core;
    // 
    numneigh_shell = list->numneigh_shell;
    firstneigh_shell = list->firstneigh_shell;
    ipage_shell = list->ipage_shell;
  }  
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

  int inum = 0;
  ipage->reset();
  //-------------------//
  if (rbl_local){
    ipage_core->reset();
    ipage_shell->reset();
  }
  //^^^^^^^^^^^^^^^^^^^//
 
  for (i = 0; i < nlocal; i++) {
    n = 0; 
    neighptr = ipage->vget();
    //-------------------//
    if (rbl_local){
      n1 = n2 = 0;
      neighptr_core = ipage_core->vget();
      neighptr_shell = ipage_shell->vget();      
    }
    //^^^^^^^^^^^^^^^^^^^//
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i

    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }

      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
      //--------------zjq--------------------------------//
      if (rbl_local){
        if (rsq <= cutcoresq[itype][jtype]) {
          neighptr_core[n1++] = j; 
        } else if (rsq <= cutneighsq[itype][jtype]){
          neighptr_shell[n2++] = j;
        }
      }      
      //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^// 
    }
   
    // loop over all atoms in other bins in stencil, store every pair

    ibin = atom2bin[i];

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
        //--------------zjq--------------------------------//
        if (rbl_local){
          if (rsq <= cutcoresq[itype][jtype]) {
            neighptr_core[n1++] = j; 
          } else if (rsq <= cutneighsq[itype][jtype]){
            neighptr_shell[n2++] = j;
          }
        }      
      //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^// 
      }
    }
    
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n; 
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    //--------------------------------zjq--------------------------------//
    if (rbl_local){
      firstneigh_core[i] = neighptr_core;
      numneigh_core[i] = n1;
      ipage_core->vgot(n1);
      if (ipage_core->status())
        error->one(FLERR,"Neighbor_core list overflow, boost neigh_modify one");
      // 
      firstneigh_shell[i] = neighptr_shell;
      numneigh_shell[i] = n2;
      ipage_shell->vgot(n2);
      if (ipage_shell->status())
        error->one(FLERR,"Neighbor_shell list overflow, boost neigh_modify one");
    }
    
    //-------------------------------------------------------------------//
  }

  list->inum = inum;
}
