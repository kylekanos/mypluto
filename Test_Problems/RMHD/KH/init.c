#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double x,y, arg, eps, kx;
  double alpha, beta;
  static double pr0;
 
  if (first_call == 1){
    arg = g_inputParam[MACH]/g_inputParam[VEL0];
    arg = arg*arg;
    #if EOS == IDEAL
     g_gamma = 4./3.;
     pr0 = (g_gamma - 1.0)/((g_gamma - 1.0)*arg - 1.0)/g_gamma;
    #else
    { 
      double a,b,c;
      a = 15.0 - 6.0*arg; b = 24.0 - 10.0*arg; c = 9.0;
      arg = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;
      pr0 = 2.0/3.0*sqrt(arg*arg/(1.0 - arg*arg));
    }
    #endif
    first_call = 0;
  }

  x = x1;
  y = x2;

  kx    = 2.0*CONST_PI*x1;
  alpha = 1.0/100.0;
  beta  = 1.0/10.0;

  arg = y*y/(beta*beta);
  eps = 0.01*g_inputParam[VEL0];

  us[RHO] = 1.0;
  us[VX1] = -g_inputParam[VEL0]*tanh(y/alpha);
  us[VX2] =  eps*0.5*(sin(kx) - sin(-kx))*exp(-arg);
  us[VX3] = 0.0;

/* --------------------------------------------
    find pressure values given the Mach number
   -------------------------------------------- */

  us[PRS] = pr0;

/* us[PRS] = 20.0;  */

  us[TR] = (y < 0.0 ? 1.0:-1.0);

  #if PHYSICS == MHD || PHYSICS == RMHD

   us[BX1] = sqrt(2.0*us[PRS]*g_inputParam[SIGMA_POL]);
   us[BX2] = 0.0;
   us[BX3] = sqrt(2.0*us[PRS]*g_inputParam[SIGMA_TOR]);

   us[AX1] = 0.0;
   us[AX2] = 0.0;
   us[AX3] = y*us[BX1];

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
void BACKGROUND_FIELD (real x1, real x2, real x3, real *B0)
/* 
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background 
 *   magnetic field.
 *
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    vector component of the background field.
 *
 *
 **************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
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
{ }

