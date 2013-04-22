/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Double Mach Reflection Problem.
 
  Sets the initial condition for a planar shock front making an angle of
  pi/3 with a reflecting wall.
  The wedge is represented by a reflecting boundary starting at x=1/6 
  along the lower y-boundary.
  As the shock reflects off the lower wall, a jet of denser gas forms.
  \b References
    - "The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shocks"\n
       Woodward, P.R., & Colella, P., JCP (1984) 54, 115.
    - "Comparison of some FLux Corrected Transport and Total variatiion
       diminishing numerical schemes for hydrodynamics and magnetohydrodynamic
       problems" \n
       Toth, G., Odstrcil, D., JCP (1996) 128, 82.

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
 *********************************************************************** */
{
  double alpha, xs;

  g_gamma = 1.4;
  alpha   = 1.0/3.0*CONST_PI;
  xs      = 1.0/6.0 + x2/tan(alpha);

  if (x1 > xs){
    us[RHO] = 1.4;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
    us[PRS] = 1.0;
  }else{
    us[RHO] = 8.0;
    us[VX1] =  8.25*sin(alpha);
    us[VX2] = -8.25*cos(alpha);
    us[PRS] = 116.5;
  }      
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
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
  real  x1, x2, x3;
  real  alpha,xs;
  
  alpha = 60./180.*CONST_PI;

  if (side == X1_BEG){

    X1_BEG_LOOP(k,j,i) {  

      d->Vc[RHO][k][j][i] = 8.0;
      d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
      d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
      d->Vc[PRS][k][j][i] = 116.5;
      
    }

  } else if (side == X2_BEG){

    X2_BEG_LOOP(k,j,i) {  

      x1 = grid[IDIR].x[i];
      if (x1 < 1.0/6.0){
        d->Vc[RHO][k][j][i] = 8.0;
        d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
        d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
        d->Vc[PRS][k][j][i] = 116.5;
      }else{                                   /* reflective boundary */

        d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =   d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = - d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][2*JBEG - j - 1][i];

      }                    
    }

  } else if (side == X2_END) {

    X2_END_LOOP(k,j,i){

      x1 = grid[IDIR].x[i];
      xs = 10.0*g_time/sin(alpha) + 1.0/6.0 + 1.0/tan(alpha);
      if (x1 < xs){
        d->Vc[RHO][k][j][i] = 8.0;
        d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
        d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
        d->Vc[PRS][k][j][i] = 116.5;
      }else{                                   /* reflective boundary */
        d->Vc[RHO][k][j][i] = 1.4;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = 0.;
        d->Vc[PRS][k][j][i] = 1.;
      }                    
    }
  }
}

