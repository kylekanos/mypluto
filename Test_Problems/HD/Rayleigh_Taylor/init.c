/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Hydrodynamical Rayleigh-Taylor Instability.


  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double x,y, eps;

  x = x1;
  y = x2;

  eps = 5.e-3*(cos(2.0*CONST_PI*x1) + cos(-2.0*CONST_PI*x1));
  if (x2 < 2.0*(1.0 - eps) ){
    us[RHO] = 1.0;
    us[PRS] = 8.0 - x2;
  }else{
    us[RHO] = 2.0;
    us[PRS] = 10.0 - 2.0*x2;
  }
  us[VX1] = us[VX2] = 0.0;

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
  int  i, j, k;
  double x1, x2, x3;
  double rho, pr;
  
  if (side == X2_BEG){
    X2_BEG_LOOP(k,j,i){  
      x2 = grid[JDIR].x[j];  
      rho = 1.0;
      pr  = 8.0 - x2; 

      d->Vc[RHO][k][j][i] = rho;
      d->Vc[VX1][k][j][i] = d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[PRS][k][j][i] = pr;
    }
  } else if (side == X2_END) {
    X2_END_LOOP(k,j,i){  
      x2 = grid[JDIR].x[j];  

      rho = 2.0;
      pr  = 10.0 - 2.0*x2;

      d->Vc[RHO][k][j][i] = rho;
      d->Vc[VX1][k][j][i] = d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[PRS][k][j][i] = pr;
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] =  0.0;
  g[JDIR] = -1.0;
  g[KDIR] =  0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  return x2;
}
#endif
