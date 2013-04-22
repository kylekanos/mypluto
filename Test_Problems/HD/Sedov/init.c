/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Sedov-Taylor Blast Wave.

  Set the initial condition for a Sedov-Taylor blast wave problem.
  The input parameters read from pluto.ini are labeled as:
  - DNST0: initial density
  - ENRG0: initial energy (NOT energy density)
  - GAMMA: specific heat ratio.


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
  double dr, vol, r;

  g_gamma = g_inputParam[GAMMA];

/* --------------------------------------------------
    dr is the size of the initial energy deposition 
    region: 2 ghost zones.
   -------------------------------------------------- */

  #if DIMENSIONS == 1
   dr = 2.0*(g_domEnd[IDIR]-g_domBeg[IDIR])/(double)NX1;
  #else
   dr = 3.5*(g_domEnd[IDIR]-g_domBeg[IDIR])/(double)NX1;
  #endif

/* ---------------------------------------
     compute region volume 
   --------------------------------------- */

  #if (GEOMETRY == CARTESIAN) && (DIMENSIONS == 1)
   vol = 2.0*dr;
  #elif (GEOMETRY == CYLINDRICAL && DIMENSIONS == 1)|| \
        (GEOMETRY == CARTESIAN   && DIMENSIONS == 2)
   vol = CONST_PI*dr*dr;
  #elif (GEOMETRY == SPHERICAL   && DIMENSIONS == 1)|| \
        (GEOMETRY == CYLINDRICAL && DIMENSIONS == 2)|| \
        (GEOMETRY == CARTESIAN   && DIMENSIONS == 3)
   vol = 4.0/3.0*CONST_PI*dr*dr*dr;
  #else
   print1 ("! Init: geometrical configuration not allowed\n");
   QUIT_PLUTO(1);
  #endif

  r = EXPAND(x1*x1, + x2*x2, +x3*x3);
  r = sqrt(r);

  us[RHO] = g_inputParam[DNST0];
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;

  if (r <= dr) {
    us[PRS] =  (g_gamma - 1.0)*g_inputParam[ENRG0]/vol;
  }else{
    us[PRS] = 1.e-5;
  }

  us[TRC]  = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

   us[BX1] = 0.0;
   us[BX2] = 0.0;
   us[BX3] = 0.0;

   #ifdef STAGGERED_MHD

    us[AX] = 0.0;
    us[AY] = 0.0;
    us[AZ] = 0.0;

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
{ }

