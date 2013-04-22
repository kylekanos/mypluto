#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* -------------------------------------------------------------
    Set inital condition for the fast rotor problem in 
    cartesian and polar geometries;
 
    Reference: 
 
     "On the divergence-free condition in Godunov-type schemes 
      for ideal MHD: the upwind constrained transport method"
  
     P. Londrillo, L. Del Zanna 
     JCP (2004), 195, 17
   ------------------------------------------------------------- */
  double r, r0, r1, Bx, f, omega;

  g_gamma   = 1.4;
  omega = 20.0;
  Bx    = 5.0/sqrt(4.0*CONST_PI);
  r0    = 0.1;
  r1    = 0.115;

  #if GEOMETRY == CARTESIAN

   r = sqrt(x1*x1 + x2*x2);

   us[PRS] = 1.0;
   us[BX1] = Bx;
   us[BX2] = 0.0;

   f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     us[RHO] = 10.0;
     us[VX1] = -omega*x2;
     us[VX2] =  omega*x1;
   }else if (r < r1) {
     us[RHO] = 1.0 + 9.0*f;
     us[VX1] = -f*omega*x2*r0/r;
     us[VX2] =  f*omega*x1*r0/r;
   } else {
     us[RHO] = 1.0;
     us[VX1] = 0.0;
     us[VX2] = 0.0;
   }

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = us[BX1]*x2;

  #elif GEOMETRY == POLAR

   r = x1;

   us[BX1] =  Bx*cos(x2);
   us[BX2] = -Bx*sin(x2);
   us[VX1] = 0.0;
   us[PRS] = 1.0;

  f = (r1 - r)/(r1 - r0);
   if (r <= r0) {
     us[RHO] = 10.0;
     us[VX2] = omega*r;
   }else if (r < r1) {
     us[RHO] = 1.0 + 9.0*f;
     us[VX2] = f*omega*r0;
   } else {
     us[RHO] = 1.0;
     us[VX2] = 0.0;
   }
  
   us[AX1] = us[AX2] = 0.0;
   us[AX3] = - r*sin(x2)*Bx;

  #endif

  #if BACKGROUND_FIELD == YES
   us[BX1] = us[BX2] = us[BX3] =
   us[AX1] = us[AX2] = us[AX3] = 0.0;
  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

#if PHYSICS == MHD
/* ************************************************************** */
void BackgroundField (real x1, real x2, real x3, real *B0)
/* 
 *
 *
 *
 *
 **************************************************************** */
{
  real Bx;
  Bx = 5.0/sqrt(4.0*CONST_PI);

  #if GEOMETRY == CARTESIAN
   B0[0] = Bx;
   B0[1] = 0.0;
   B0[2] = 0.0;
  #elif GEOMETRY == POLAR
   B0[0] =   Bx*cos(x2);
   B0[1] = - Bx*sin(x2);
   B0[2] = 0.0;
  #endif
}
#endif

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
  int   i, j, k, nv;
  double  *r, slp;

  if (side == X1_BEG){

    r = grid[IDIR].x;
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        slp = r[i]/r[IBEG];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[VX1][k][j][i] = slp*d->Vc[VX1][k][j][IBEG];
        d->Vc[VX2][k][j][i] = slp*d->Vc[VX2][k][j][IBEG]; 
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
      #endif
    }
  }
}

