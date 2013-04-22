#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  g_gamma = 5./3.;

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = 1.0/g_gamma;
  us[TRC] = 0.0;

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  *R;
  real    pjet, dnjet, vjet;
  real    scrh;

  scrh = 1.0/(g_gamma - 1.0);

  R     = grid[IDIR].xgc;
  pjet  = g_inputParam[PR_RATIO]*pow(2.0/(g_gamma + 1.0),g_gamma*scrh)/g_gamma;
  dnjet = g_inputParam[DN_RATIO]*pow(2.0/(g_gamma + 1.0),scrh);
  vjet  = sqrt(g_gamma*pjet/dnjet);

  if (side == X2_BEG){

    X2_BEG_LOOP(k,j,i){

      if (R[i] <= 1.) {
        d->Vc[RHO][k][j][i] = dnjet;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = vjet;
        d->Vc[PRS][k][j][i] = pjet;
      } else {
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][2*JBEG - j - 1][i];
      }
    }
  } 
}

