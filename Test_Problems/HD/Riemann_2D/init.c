/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Two-Dimensional Riemann problem.

  Sets the initial condition for a 2D Riemann problem in terms of
  4 four constant states (PP, PM, MP, MM) defined by 4 flow quantities:
  density (DN), pressure (PR), x- and y- components of velocity (VX and VY).
  The 16 parameters are read from pluto.ini.

  \authors A. Mignone (mignone@ph.unito.it)
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
  double x,y;

  x = x1;
  y = x2;
  #if EOS == IDEAL
   g_gamma = 1.4;
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = 1.0;
  #endif

  if (x > 0.0 && y > 0.0){
    us[RHO] = g_inputParam[DN_PP];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_PP];
    #endif
    us[VX1] = g_inputParam[VX_PP];
    us[VX2] = g_inputParam[VY_PP];
  }else if(x < 0.0 && y > 0.0){
    us[RHO] = g_inputParam[DN_MP];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_MP];
    #endif
    us[VX1] = g_inputParam[VX_MP];
    us[VX2] = g_inputParam[VY_MP];
  }else if(x < 0.0 && y < 0.0){
    us[RHO] = g_inputParam[DN_MM];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_MM];
    #endif
    us[VX1] = g_inputParam[VX_MM];
    us[VX2] = g_inputParam[VY_MM];
  }else if(x > 0.0 && y < 0.0){
    us[RHO] = g_inputParam[DN_PM];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_PM];
    #endif
    us[VX1] = g_inputParam[VX_PM];
    us[VX2] = g_inputParam[VY_PM];
  }   
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

