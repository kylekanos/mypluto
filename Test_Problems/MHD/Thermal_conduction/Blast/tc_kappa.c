#include "pluto.h"

void TC_kappa(double *v, double x1, double x2, double x3, 
              double *kpar, double *knor, double *phi)
{
  double mu, nH, sqT, T, B2_cgs;
  double ctts_mu = 1.265060;

  mu   = ctts_mu/2.0;
  T    = v[PRS]/v[RHO]*(mu*CONST_mp*g_unitVelocity*g_unitVelocity/CONST_kB);
  sqT  = sqrt(T);

  *kpar = 5.6e-7*T*T*sqT;

  #if PHYSICS == MHD
   nH      = v[DN]*g_unitDensity/(mu*CONST_mp);
   B2_cgs  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]) + 1.e-12;
   B2_cgs *= 4.0*CONST_PI*g_unitDensity*g_unitVelocity*g_unitVelocity;
   *knor   = 3.3e-16*nH*nH/sqT/B2_cgs;
  #else
   *knor = 0.0;
  #endif

  *kpar *= CONST_mp*mu/(g_unitDensity*g_unitVelocity*g_unitLength*CONST_kB);
  *knor *= CONST_mp*mu/(g_unitDensity*g_unitVelocity*g_unitLength*CONST_kB);

  *phi = 0.3;




{
  double P0, T0, sqT0;
  double Kpar0, Knor0, Kpar, Knor;
  double mu, T5_2, scrh, sqT, T;
  double ctts_mu = 1.265060;
  double kpar_old, knor_old;

  mu   = ctts_mu/2.0;
  P0   = g_unitDensity*g_unitVelocity*g_unitVelocity;
  T0   = mu*CONST_mp*P0/(CONST_kB*g_unitDensity);
  sqT0 = sqrt(T0);

  Kpar0  = P0*g_unitVelocity*g_unitLength/(T0*T0*T0*sqT0);
  Knor0  = 4.0*CONST_PI*CONST_kB*CONST_kB*g_unitVelocity*g_unitLength;
  Knor0 *= T0*sqT0;

/*
if (fabs(4.0*CONST_PI*sqrt(CONST_kB*mu*CONST_mp)*mu*CONST_mp*pow(g_unitVelocity,4.0)*g_unitLength/Knor0-1.0) >
    1.e-12){
 print ("Wrong constant\n");
 QUIT_PLUTO(1);
}
*/
  Kpar = 5.6e-7/Kpar0;
  Knor = 3.3e-16/Knor0;

  T    = v[PRS]/v[RHO];
  sqT  = sqrt(T);
  T5_2 = T*T*sqT;


  kpar_old = Kpar*T5_2;  /* Coeff. di conduzione parallelo */
  #if PHYSICS == MHD
   scrh = v[RHO]*v[RHO]/sqT/(v[BX1]*v[BX1] + v[BX2]*v[BX2] + 1.e-12);
   knor_old = Knor*scrh;  /* Coeff. di conduzione ortogonale */
  #else
   knor_old = 0.0;
  #endif

if (fabs(*kpar/kpar_old - 1.0) > 1.e-12){
 print ("kpar = %12.6e, kpar_Old = %12.6e\n", *kpar, kpar_old);
 QUIT_PLUTO(1);
}

if (fabs(*knor/knor_old - 1.0) > 1.e-12){
 print ("knor = %12.6e, knor_Old = %12.6e\n", *knor, knor_old);
 QUIT_PLUTO(1);
}

}

}

