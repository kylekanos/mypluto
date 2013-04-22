#include "pluto.h"
#include <stdio.h>
/*---------------------------------------------------------------------------*/
/*---- Specification of explicit first and second viscosity coefficients ----*/
/*---------------------------------------------------------------------------*/

void eta_visc_func(real *v, real x1, real x2, real x3, 
                   real *eta1_visc, real *eta2_visc )
{
double Reynolds = g_inputParam[Re]; 
	
 *eta1_visc = g_inputParam[OMEGA]*g_inputParam[R_INT]*(g_inputParam[R_EXT]-g_inputParam[R_INT])*v[RHO]/Reynolds;
 *eta2_visc = 0.0;
}
