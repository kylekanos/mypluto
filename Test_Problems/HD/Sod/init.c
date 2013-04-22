/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Sod shock tube problem.

  The Sod shock tube problem is perhaps one of the most used 
  benchmark for shock-capturing schemes. 
  It is a one-dimensional problem with initial condition given 
  by a discontinuity separating two constant states:
  \f$ (\rho, p) = (1,1) \f$ for \f$ x < 1/2 \f$ and 
  \f$ (\rho, p) = (1/8,1/10) \f$ for \f$ x > 1/2 \f$.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*!
 *  Sets initial condition for the Sod shock-tube test problem.
 *
 *
 *********************************************************************** */
{
  g_gamma = 1.4;
  if (fabs(x1) < 0.5) {
    us[RHO] = 1.0;
    us[PRS] = 1.0;
  }else{
    us[RHO] = 0.125;
    us[PRS] = 0.1;
  }

  us[VX1] = 0.0;
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


