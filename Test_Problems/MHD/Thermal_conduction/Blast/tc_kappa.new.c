#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, 
              double *kpar, double *knor, double *phi)
{
  double mu, nH, sqT, T;
  double ctts_mu = 1.265060;

  mu   = ctts_mu/2.0;
  T    = v[PRS]/v[RHO]*KELVIN*mu;
  sqT  = sqrt(T);

  *kpar = 5.6e-7*T*T*sqT;

  #if PHYSICS == MHD
   nH     = us[DN]*g_unitDensity/(1.26*CONST_mp);
   Bmag2  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]) + 1.e-12;
   Bmag2 *= 4.0*CONST_PI*g_unitDensity*g_unitVelocity*g_UnitVelocity;
   *knor  = 3.3e-16*nH*nH/sqT/Bmag2;
  #else
   *knor = 0.0;
  #endif

  *kpar *= CONST_amu/(g_unitDensity*g_unitVelocity*g_unitLength*CONST_kB);
  *knor *= CONST_amu/(g_unitDensity*g_unitVelocity*g_unitLength*CONST_kB);
  *phi = 0.3;
 
}

