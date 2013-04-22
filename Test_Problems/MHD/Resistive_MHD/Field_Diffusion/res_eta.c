#include "pluto.h"

void ETA_Func(real *v, real x1, real x2, real x3, real *eta)
{
  eta[IDIR] = g_inputParam[ETAX];
  eta[JDIR] = g_inputParam[ETAY];
  eta[KDIR] = g_inputParam[ETAZ];
}
