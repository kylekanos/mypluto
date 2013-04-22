#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 * Set a 2D Blast wave problem with thermal conduction.
 *
 *  Input Parameters are:
 *
 *  T_IN, T_OUT: temperature inside and outside the circle (in K).
 *  RHO_IN, RHO_OUT: density inside and outside (in dimensionless units);
 *  BMAG:            Magnetic field strength (in Gauss)
 *  THETA:           Orientation of the field (in degrees).
 *
 *  Configuration 1-4 have the same initial condition and are done
 *  either with an explicit time stepping or STS, HD and MHD and 
 *  do not show evidence for any numerical artifact.
 *  Configuration 5-6, on the other hand, show that STS suffers from
 *  some kind of unstable behavior due to the flux limiter switching
 *  from classical to saturated regimes. Only small CFL (0.1 or less) 
 *  or larger values of STS_nu (e.g 0.05) mitigate the problem.
 *  Future improvement (RKC ?) should address this issue.
 *
 *********************************************************************** */
{
  static int first_call=1;
  double r, r0, mu, T, Unit_Pressure, prof;

  mu  = 1.26;
  g_gamma = 5.0/3.0;        

  g_unitDensity  = mu*CONST_mH;
  g_unitLength   = CONST_pc;
  g_unitVelocity = sqrt(g_gamma*CONST_kB*1.e6/(mu*CONST_mH)); /* cs=1 at T = 1.e6 K */
  Unit_Pressure = g_unitDensity*g_unitVelocity*g_unitVelocity;  
  
/* ----------------------------------------------
             Use c.g.s units 
   ---------------------------------------------- */

  r    = sqrt(EXPAND(x1*x1, + x2*x2, + x3*x3));
  r0   = 1.0;   /* -- cloud radius -- */
  prof = 1.0/cosh(10.0*pow(r/r0,10));

  us[VX1] = us[VX2] = 0.0;
  T      = g_inputParam[T_OUT]   + (g_inputParam[T_IN]   - g_inputParam[T_OUT])*(r <= r0);
  us[RHO] = g_inputParam[RHO_OUT] + (g_inputParam[RHO_IN] - g_inputParam[RHO_OUT])*prof;

  us[PRS] = T*us[RHO]/KELVIN;

  #if PHYSICS == MHD
   us[BX1] = g_inputParam[BMAG]*cos(g_inputParam[THETA]*CONST_PI/180.0);
   us[BX2] = g_inputParam[BMAG]*sin(g_inputParam[THETA]*CONST_PI/180.0);
   us[BX3] = 0.0;
   Unit_Pressure = g_unitDensity*g_unitVelocity*g_unitVelocity;  
   us[BX1] /= sqrt(Unit_Pressure*4.0*CONST_PI);
   us[BX2] /= sqrt(Unit_Pressure*4.0*CONST_PI);
   us[BX3] /= sqrt(Unit_Pressure*4.0*CONST_PI);

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = x2*us[BX1] - x1*us[BX2];
  #endif

  #ifdef GLM_MHD 
   us[PSI_GLM] = 0.0;
  #endif

}   
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{ }

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
  static int  first_call = 1;
  int   i, j, k, nv;
  static double vin[256];

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (first_call){
      Init(vin, 0.0, -10.0, 0.0);
      first_call = 0;
    }
    if (box->vpos == CENTER){
      for (nv = 0; nv < NVAR; nv++ ) BOX_LOOP(box,k,j,i){
        d->Vc[nv][k][j][i] = vin[nv];
      }
    }
  }
}

