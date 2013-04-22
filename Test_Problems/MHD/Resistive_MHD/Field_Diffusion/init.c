#include "pluto.h"

void BoundValues (double *v, double x1, double x2, double x3, double t);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  BoundValues(us, x1, x2, x3, 1.0);
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *x, *y, *z;
  double *xp, *yp, *zp;
  double t, dt, vb[NVAR];

  if (side == 0){  
    #if GEOMETRY == CARTESIAN
    TOT_LOOP(k,j,i){
      d->Vc[VX1][k][j][i] = 0.0;
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[VX3][k][j][i] = 0.0;
    }
    #endif
  }

  t = g_time + 1.0;

  x = grid[IDIR].x; xp = grid[IDIR].xr;
  y = grid[JDIR].x; yp = grid[JDIR].xr;
  z = grid[KDIR].x; zp = grid[KDIR].xr;

  if (side == X1_BEG || side == X2_BEG || side == X3_BEG ||
      side == X1_END || side == X2_END || side == X3_END){  /* -- All Boundaries -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        BoundValues (vb, x[i], y[j], z[k], t);
        for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = vb[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i){
         BoundValues (vb, xp[i], y[j], z[k], t);
         d->Vs[BX1s][k][j][i] = vb[BX1];
       }
      #endif
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i){
         BoundValues (vb, x[i], yp[j], z[k], t);
         d->Vs[BX2s][k][j][i] = vb[BX2];
       }
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) {
         BoundValues (vb, x[i], y[j], zp[k], t);
         d->Vs[BX3s][k][j][i] = vb[BX3];
       }
      #endif
    }
  }
}

/* ******************************************************************* */
void BoundValues (double *v, double x1, double x2, double x3, double t)
/*
 *
 * Time-dependent exact solution
 * 
 ********************************************************************* */
{
  double Bx, By, Bz;
  double x, y, z;

/* -- find Cartesian coordinates from (x1,x2,x3) -- */

  #if GEOMETRY == CARTESIAN
   x = x1; y = x2; z = x3;
  #elif GEOMETRY == POLAR
   x = x1*cos(x2) - 5.0; 
   y = x1*sin(x2) - 5.0; 
   z = x3 - 5.0;
  #elif GEOMETRY == SPHERICAL
   x = x1*cos(x3)*sin(x2) - 5.0;
   y = x1*sin(x3)*sin(x2) - 5.0;
   z = x1*cos(x2) - 5.0;
  #else
   #error geometry not valid
  #endif

  v[RHO] = 1.0e9;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if EOS != ISOTHERMAL
   v[PRS] = 1.0;
  #endif
  #ifdef GLM_MHD
   v[PSI_GLM] = 0.0;
  #endif

/* -- find Cartesian components of magnetic field -- */

  Bx = exp(-0.25*(y*y)/g_inputParam[ETAZ]/t)*exp(-0.25*(z*z)/g_inputParam[ETAY]/t)/t;
  By = exp(-0.25*(x*x)/g_inputParam[ETAZ]/t)*exp(-0.25*(z*z)/g_inputParam[ETAX]/t)/t;
  Bz = exp(-0.25*(x*x)/g_inputParam[ETAY]/t)*exp(-0.25*(y*y)/g_inputParam[ETAX]/t)/t;

/* -----------------------------------------------------------
     Transform back to original coordinate system
   ----------------------------------------------------------- */
  
  #if GEOMETRY == CARTESIAN
   v[BX1] = Bx;
   v[BX2] = By;
   v[BX3] = Bz;
  #elif GEOMETRY == POLAR
   v[BX1] =  cos(x2)*Bx + sin(x2)*By;
   v[BX2] = -sin(x2)*Bx + cos(x2)*By;
   v[BX3] =  Bz;
  #elif GEOMETRY == SPHERICAL
   v[BX1] =  cos(x3)*sin(x2)*Bx + sin(x3)*sin(x2)*By + cos(x2)*Bz;
   v[BX2] =  cos(x3)*cos(x2)*Bx + sin(x3)*cos(x2)*By - sin(x2)*Bz;
   v[BX3] = -sin(x3)*Bx + cos(x3)*By;
  #else
   print1 ("! BoundValues: GEOMETRY not defined\n");
   QUIT_PLUTO(1);
  #endif  

}

