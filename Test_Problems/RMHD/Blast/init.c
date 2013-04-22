#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *  Set the initial condition for the relativistic
 *  blast wave problem in either 2D or 3D.
 *
 *  The input parameters are:
 * 
 *   g_inputParam[P_IN]:    pressure inside the initial circular (2D) or
 *                 spherical (3D) region.
 *   g_inputParam[P_OUT]:   ambient pressure
 *   g_inputParam[RHO_OUT]: ambient density
 *   g_inputParam[BMAG]:    magnetic field intensity
 *   g_inputParam[THETA]:   angle between mag. field and z-axis
 *   g_inputParam[PHI]:     angle between xy mag. field and x-axis
 *   g_inputParam[RADIUS]:  radius of the initial over-pressurized 
 *                 region.
 *
 *  Configurations #1 and #2 are taken from 
 *
 *    Del Zanna et al, A&A (2003) 400,397
 *
 *  Configuration #3 is taken from 
 *
 *    Mignone et al, JCP (2007), 170, 228
 *    (rescaling {rho,p} -> 100*{rho,p}, B -> 10*B) 
 *
 *  Configuration #4 & #5 are taken from 
 *  
 *     Beckwith & Stone, ApJS (2011), 193, 6
 *     Strongly magnetized case, sec. 4.6 (Fig. 14)
 *
 *********************************************************************** */
{
  double r, rc, theta, phi;
  double dc, pc, de, pe;

  g_gamma = 4./3.;
                        
  #if DIMENSIONS == 2
   r = sqrt(x1*x1 + x2*x2);
  #elif DIMENSIONS == 3
   r = sqrt(x1*x1 + x2*x2 + x3*x3);
  #endif                                                                                                                                  
  rc = g_inputParam[RADIUS];

  dc = g_inputParam[RHO_IN];
  pc = g_inputParam[PR_IN];
  de = g_inputParam[RHO_OUT];
  pe = g_inputParam[PR_OUT];

  if (r <= rc) {
    us[RHO] = dc;
    us[PRS] = pc;
  #if DIMENSIONS == 3
   }else if (r > rc && r < 1.0){
     us[RHO] = de*(r - rc)/(1.0 - rc) + dc*(r - 1.0)/(rc - 1.0);
     us[PRS] = pe*(r - rc)/(1.0 - rc) + pc*(r - 1.0)/(rc - 1.0);
  #endif
  }else{
    us[RHO] = de;
    us[PRS] = pe;
  }

  us[VX1] = us[VX2] = us[VX3] = 0.0;
  us[AX1] = us[AX2] = us[AX3] = 0.0;

  theta = g_inputParam[THETA]*CONST_PI/180.0;
  phi   =   g_inputParam[PHI]*CONST_PI/180.0;

  us[BX1]  = g_inputParam[BMAG]*sin(theta)*cos(phi);
  us[BX2]  = g_inputParam[BMAG]*sin(theta)*sin(phi);
  us[BX3]  = g_inputParam[BMAG]*cos(theta);

  #if GEOMETRY == CARTESIAN
   us[AX1] = 0.0;
   us[AX2] = us[BX3]*x1;
   us[AX3] = -us[BX2]*x1 + us[BX1]*x2;
  #elif GEOMETRY == CYLINDRICAL
   us[AX1] = us[AX2] = 0.0;
   us[AX3] = 0.5*us[BX2]*x1;
  #endif

  g_smallPressure = 1.e-6;  
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

