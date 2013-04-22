#include "pluto.h"

const double Kpar_cgs = 5.6E-7;
const double Knor_cgs = 3.3E-16;

void TC_kappa(double *v, double x1, double x2, double x3, 
              double *kappar, double *kapnor, double *phi)
{
  double P0, T0, sqT0;
  double Kpar0, Knor0, Kpar, Knor;
  double mu, T5_2, scrh, sqT, T;
  double ctts_mu = 1.265060;

  mu   = ctts_mu/2.0;
  P0   = g_unitDensity*g_unitVelocity*g_unitVelocity;
  T0   = mu*CONST_mp*P0/(CONST_kB*g_unitDensity);
  sqT0 = sqrt(T0);

  Kpar0  = P0*g_unitVelocity*g_unitLength/(T0*T0*T0*sqT0);
  Knor0  = 4.0*CONST_PI*CONST_kB*CONST_kB*g_unitVelocity*g_unitLength;
  Knor0 *= T0*sqT0;

  Kpar = Kpar_cgs/Kpar0;
  Knor = Knor_cgs/Knor0;

  T    = v[PRS]/v[RHO];
  sqT  = sqrt(T);
  T5_2 = T*T*sqT;

/*  *Q = 5.0*phi*v[PRS]*sqT;   */
  *phi = 0.3;

  *kappar = Kpar*T5_2;  /* Coeff. di conduzione parallelo */
  #if PHYSICS == MHD
   scrh = v[RHO]*v[RHO]/sqT/(v[BX1]*v[BX1] + v[BX2]*v[BX2] + 1.e-12);
   *kapnor = Knor*scrh;  /* Coeff. di conduzione ortogonale */
  #else
   *kapnor = 0.0;
  #endif
/*
*kappar = 1.e-30;
*kapnor = 1.e-30;
*/
}

