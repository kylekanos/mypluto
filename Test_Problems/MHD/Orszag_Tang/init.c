#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */

{

/* ------------------------------------------------------
    Orszag Tang MHD Vortex Problem

    Reference:
    "Comparison of some Flux Corrected Transport
     and TVD numerical schemes for hydrodynamic
     and magnetohydrodynamic problems
   ------------------------------------------------------ */

  double x,y,z;

  x = x1;
  y = x2;
  z = x3;

  us[VX1] = - sin(y);
  us[VX2] = sin(x);
  us[VX3] = 0.0;
  us[BX1] = - sin(y);
  us[BX2] = sin(2.0*x);
  us[BX3] = 0.0;
  us[RHO] = 25./9.;
  #if EOS != ISOTHERMAL && EOS != BAROTROPIC
   us[PRS] = 5.0/3.0;
  #endif
  us[TR] = (y>CONST_PI)*1.0;

  us[AX1] = 0.0;
  us[AX2] = 0.0;
  us[AX3] = cos(y) + 0.5*cos(2.0*x);


 #if DIMENSIONS == 3 
 {
   double c0 = 0.8;
/*
   us[VX1] = - sin(y);
   us[VX2] =   sin(x);
   us[VX3] =   0.0;
   
   us[BX1] = c0*( - 2.0*sin(2.0*y) + sin(z));
   us[BX2] = c0*(       sin(z) + sin(x));
   us[BX3] = c0*(       sin(x) + sin(y));
   us[RHO] = 25./9.;
   us[PRS] = 5.0/3.0;

   us[AX1] = c0*( cos(y) - cos(z));
   us[AX2] = c0*(-cos(x) + cos(z));
   us[AX3] = c0*( cos(x) + cos(2.0*y));
*/
   us[VX2] = - sin(z);
   us[VX3] =   sin(y);
   us[VX1] =   0.0;
  
   us[BX2] = c0*( - 2.0*sin(2.0*z) + sin(x));
   us[BX3] = c0*(       sin(x) + sin(y));
   us[BX1] = c0*(       sin(y) + sin(z));
   us[RHO] = 25./9.;
   us[PRS] = 5.0/3.0;

   us[AX2] = c0*( cos(z) - cos(x));
   us[AX3] = c0*(-cos(y) + cos(x));
   us[AX1] = c0*( cos(y) + cos(2.0*z));

  }
  #endif

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  printf ("I'm calling analysis (%d)\n",prank);

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
{ }

