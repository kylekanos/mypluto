#include "pluto.h"

static double Profile(double r, int nv);
static void JETVAL (double *vjet);


/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int nv;
  double vjet[NVAR], scrh;

  g_gamma = 5.0/3.0;
  
  us[RHO] = g_inputParam[RHO_OUT];
  us[VX1] = 0.0;
  us[VX2] = 0.0;

  us[PRS] = g_inputParam[PRESS_IN];

  #if USE_FOUR_VELOCITY == YES
   scrh = 1.0/sqrt(1.0 - us[VX1]*us[VX1] - us[VX2]*us[VX2]);
   us[VX1] *= scrh;
   us[VX2] *= scrh;
  #endif

  if (x2 < -1.0){
    JETVAL(vjet);
    for (nv = 0; nv < NVAR; nv++)
      us[nv] = us[nv] +(vjet[nv] - us[nv])*Profile(x1,nv);
  }   

  g_smallPressure = us[PRS]/500.0;
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
  int   nv, i, j, k;
  double  r, vjet[NVAR], vout[NVAR];
 
  if (side == X2_BEG){

    JETVAL(vjet);
    X2_BEG_LOOP(k,j,i){
 
      r = grid[IDIR].x[i];
      for (nv = 0; nv < NVAR; nv++) 
        vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i];

      vout[VX2] *= -1.0;
      for (nv = 0; nv < NVAR; nv++) 
        d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);
    }
  }
}

/* ************************************************ */
void JETVAL (double *vjet)
/*
 *
 *
 *
 *
 *
 *
 ************************************************** */
{
  vjet[RHO] = g_inputParam[RHO_IN];
  vjet[VX1] = 0.0;
  vjet[VX2] = g_inputParam[BETA];
  #if USE_FOUR_VELOCITY == YES
   vjet[VX2] /= sqrt(1.0 - g_inputParam[BETA]*g_inputParam[BETA]);
  #endif
  vjet[PRS] = g_inputParam[PRESS_IN];

}
 
/* ************************************************ */
double Profile(double r, int nv)
/* 
 *
 *
 *
 ************************************************** */
{
  int xn = 14;
  real r0 = 1.0;

  if (nv == DN) r0 = 1.1;

  return 1.0/cosh(pow(r/r0,xn));
}
