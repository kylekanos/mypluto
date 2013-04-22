#include "pluto.h"
/*
real sb_A = -0.75;
real sb_vy, sb_Lx, sb_Ly;
*/
double sb_Omega = 1.0;
double sb_q     = 1.5;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double b0, kk, cosgam, singam, v0;
  double x,y,z;
  double eps = 0.01, r2;
  #if EOS == ISOTHERMAL
   double g_gamma = 1.0;
  #endif

  x = x1; y = x2; z = x3;

  #if DIMENSIONS == 3

   kk = 24.0*CONST_PI;
   b0 = g_inputParam[CS]*sqrt(2.0/g_gamma/g_inputParam[BETA]);

   us[RHO] = 1.0; 
   r2 = 0.1*g_inputParam[CS]*sin(2.e4*(x+2.0)*(y+12)*(z+33)); /* almost random pert. */

   us[VX1]  = 0.;
   us[VX2]  = 2.0*sb_A*x + r2;
   us[VX3]  = 0.;
   
   #if EOS == IDEAL
    us[PRS] = g_inputParam[CS]*g_inputParam[CS]/g_gamma;
   #elif EOS == ISOTHERMAL
    g_isoSoundSpeed  = g_inputParam[CS];
   #endif

   us[TR] = 0.0;

   #if PHYSICS == MHD || PHYSICS == RMHD
    /* 
    us[BX1] =  0.5*eps*b0*singam*(cos(kk*z)+cos(-kk*z));
    us[BX2] = -0.5*eps*b0*cosgam*(cos(kk*z)+cos(-kk*z)); */
    us[BX1] = 0.;
    us[BX2] = 0.;
    us[BX3] = b0;

   /* us[AX1] = - 0.5*eps*b0*cosgam*(sin(kk*z)-sin(-kk*z))/kk;
    us[AX2] = - 0.5*eps*b0*singam*(sin(kk*z)-sin(-kk*z))/kk + b0*x; */

    us[AX1] = 0.;
    us[AX2] = b0*x;
    us[AX3] = 0.0;
   #endif

  #elif DIMENSIONS == 2

   b0 = g_inputParam[CS]*sqrt(2.0/g_gamma/g_inputParam[BETA]);
   kk = 2.0*CONST_PI;

   us[RHO] = 1.0;
   us[VX1] = 0.0;
   us[VX2] = 2.0*sb_A*x;
   us[VX3] = 0.0;


   #if EOS == IDEAL
    us[PRS] = g_inputParam[CS]*g_inputParam[CS]/g_gamma;
   #elif EOS == ISOTHERMAL
    g_isoSoundSpeed  = g_inputParam[CS];
   #endif

   us[TR] = 0.0;

   #if PHYSICS == MHD || PHYSICS == RMHD
    us[BX1] =  b0*sin(kk*y);
    us[BX2] = 0.0;
    us[BX3] = 0.0;

    us[AX1] = 0.0;
    us[AX2] = 0.0;
    us[AX3] = -b0*cos(kk*y)/kk;
   #endif
  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */

{
  int i,j,k;
  real totrho;

  return;

  totrho = 0.0;
  for (k = KBEG; k <= KEND; k++){
  for (j = JBEG; j <= JEND; j++){
  for (i = IBEG; i <= IEND; i++){
  }}}

}

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
  int  i, j, k;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1.e-2) d->Vc[RHO][k][j][i] = 1.e-2;
    }
    return;
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
  #ifdef FARGO
   g[IDIR] = 0.0;
   g[JDIR] = sb_q*sb_Omega*v[VX1];
  #else
   g[IDIR] = 2.0*sb_Omega*sb_Omega*sb_q*x1;
   g[JDIR] = 0.0;
  #endif

/* -- vertical gravity ? -- */

  g[KDIR] = 0.0;
/*  g[KDIR] = -sb_Omega*sb_Omega*x3;    */
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
  return 0.0;
}
#endif
