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

#ifdef PAIR_CLASS

PairStyle(eam,PairEAM)

#else

#ifndef LMP_PAIR_EAM_H
#define LMP_PAIR_EAM_H

#include "pair.h"
extern double times_rho;
extern double times_phi;

namespace LAMMPS_NS {


class PairEAM : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;   // Alex Stukowski option

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho,nr;
  int nfrho,nrhor,nz2r;
  double **frho,**rhor,**z2r;
  int *type2frho,**type2rhor,**type2z2r;

  // potentials in spline form used for force computation

  double dr,rdr,drho,rdrho,rhomax;
  double ***rhor_spline,***frho_spline,***z2r_spline;

  PairEAM(class LAMMPS *);
  virtual ~PairEAM();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void swap_eam(double *, double **);

 protected:
  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;

  // per-atom arrays

  double *rho,*fp;

  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct Funcfl {
    char *file;
    int nrho,nr;
    double drho,dr,cut,mass;
    double *frho,*rhor,*zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,**rhor,***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,***rhor,***z2r;
  };
  Fs *fs;

  // // ---------------------------------//
  // struct Atom_info{
  //   int id_atom;
  //   double rsq_atom;
  // };
  // Atom_info *info_atom;

  // struct Atom_str {
  //   int num_i;
  //   Atom_info   info_atom_i[200];
  // };
  // Atom_str * current_atom_str;


  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);

  virtual void read_file(char *);
  virtual void file2array();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Invalid EAM potential file

UNDOCUMENTED

*/
