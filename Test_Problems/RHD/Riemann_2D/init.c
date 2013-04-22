#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double x, y, scrh;

  x = x1;
  y = x2;

  us[VX1] = us[VX2] = 0.0;
  if (x > 0.0 && y > 0.0){
    us[RHO] = 5.477875e-3;
    us[PRS] = 2.762987e-3;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
  }else if(x < 0.0 && y > 0.0){
    us[RHO] = 0.1;
    us[PRS] = 1.0;
    us[VX1] = 0.99;
    us[VX2] = 0.0;
  }else if(x < 0.0 && y < 0.0){
    us[RHO] = 0.5;
    us[PRS] = 1.0;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
  }else if(x > 0.0 && y < 0.0){
    us[RHO] = 0.1;
    us[PRS] = 1.0;
    us[VX1] = 0.0;
    us[VX2] = 0.99;
  }   

  #if USE_FOUR_VELOCITY == YES
   scrh = 1.0/sqrt(1.0 - us[VX1]*us[VX1] - us[VX2]*us[VX2]);
   us[VX1] *= scrh;
   us[VX2] *= scrh;
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

