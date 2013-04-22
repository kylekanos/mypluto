#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
/*
 *
 * PURPOSE
 *
 *  - Solve Riemann problem for the Euler equations 
 *    using th HLLE solver;
 * 
 *     Reference:    "On Godunov-Type Method near Low Densities"
 *                    B. Einfeldt, C.D. Munz, P.L. Roe,
 *                    JCP 92, 273-295 (1991)
 *   
 *
 *  - On input, it takes left and right primitive state
 *     vectors state->vL and state->vR at zone edge i+1/2;
 *     On output, return flux and pressure vectors at the
 *     same interface.
 *
 *  - Also, compute maximum wave propagation speed (cmax) 
 *    for  explicit time step computation
 *  
 * LAST_MODIFIED
 *
 *   April 11th 2008, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 **************************************************************************** */
{
  int    nv, i;
  double   scrh;
  static double *pL, *pR, *SL, *SR, *a2L, *a2R;
  static double **fL, **fR;
  double *uR, *uL;
  double bmax, bmin, *vL, *vR, aL, aR;

/* -- Allocate memory -- */

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);
    SR = ARRAY_1D(NMAX_POINT, double);
    SL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

  HLL_Speed (state->vL, state->vR, a2L, a2R, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i]  = scrh;

/*    
    bmin = MIN(0.0, SL[i]);
    bmax = MAX(0.0, SR[i]);
    scrh = 1.0/(bmax - bmin);    
    for (nv = NFLX; nv--;   ){
      state->flux[i][nv]  = bmin*bmax*(uR[i][nv] - uL[i][nv]) 
                        +   bmax*fL[i][nv] - bmin*fR[i][nv];
      state->flux[i][nv] *= scrh;
    }
    state->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh;
vL = state->vL[i];
vR = state->vR[i];
aL = sqrt(g_gamma*vL[PRS]/vL[RHO]);
aR = sqrt(g_gamma*vR[PRS]/vR[RHO]);
    *cmax  = MAX(*cmax, MAX(-bmin,bmax));
    scrh = (fabs(vL[VXn]) + fabs(vR[VXn]))/(aL + aR);
    g_maxMach = MAX(scrh, g_maxMach);
*/

    if (SL[i] > 0.0){
    
      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = pL[i];

    }else if (SR[i] < 0.0){

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = pR[i];

    }else{

      uR = state->uR[i];
      uL = state->uL[i];

      scrh = 1.0 / (SR[i] - SL[i]);
      for (nv = NFLX; nv--;  ) {
        state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;

    }

  } /* end loops on points */
}
