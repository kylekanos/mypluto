#include"pluto.h"

#if EOS == IDEAL
/* ***************************************************************************** */
void HLLC_Solver (const State_1D *state, int beg, int end, 
                  real *cmax, Grid *grid)
/* 
 *
 *
 * NAME
 *
 *   HLLC_SOLVER
 *
 *
 * PURPOSE
 *
 *  - Solve riemann problem for the MHD equations using the 
 *    two-state HLLC Riemann solver
 * 
 *     Reference:   "An HLLC RIemann Solver for MHD" 
 *                  S. Li, 2005, JCP 203, 344
 *
 *     However, our differ from Li's solver in the way
 *     transverse momenta are computed.
 * 
 *   - On input, it takes left and right primitive state
 *     vectors state->vL and state->vR at zone edge i+1/2;
 *     On output, return flux and pressure vectors at the
 *     same interface.
 *
 *   - Also, compute maximum wave propagation speed (cmax) 
 *     for  explicit time step computation
 *  
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 ******************************************************************************** */
{
  int   nv, i;
  real  scrh;
  real  pl, pr;
  real  vBl, usl[NFLX];
  real  vBr, usr[NFLX];

  real  Bxs, Bys, Bzs, ps, vBs;
  real  vxl, vxr, vxs, vys, vzs;
  real  Fhll[NFLX], alpha_l, alpha_r;
  real  **bgf, *vL, *vR, *uL, *uR, *SL, *SR;
  static real **fL, **fR, **Uhll;
  static double **VL, **VR, **UL, **UR;
  static double *pL, *pR, *a2L, *a2R;

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);
    Uhll = ARRAY_2D(NMAX_POINT, NFLX, double);
    pL  = ARRAY_1D(NMAX_POINT, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }
  
  #if BACKGROUND_FIELD == YES
   print ("! Background field splitting not allowed with HLLC solver\n");
   QUIT_PLUTO(1);
  #endif
  
  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, pL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, pR, beg, end);

 /* ----------------------------------------
       get max and min signal velocities
     ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, bgf, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {
    
    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

/* --------------------------------------------
              compute fluxes 
   -------------------------------------------- */

    if (SL[i] >= 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fL[i][nv];
      }
      state->press[i] = pL[i];

    }else if (SR[i] <= 0.0){

      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fR[i][nv];
      }
      state->press[i] = pR[i];

    }else{

      vL = VL[i]; uL = UL[i];
      vR = VR[i]; uR = UR[i];

  /* ----  define hll states  ----  */

      scrh = 1.0/(SR[i] - SL[i]);
      for (nv = 0; nv < NFLX; nv++){  
        Uhll[i][nv] =   SR[i]*uR[nv] - SL[i]*uL[nv] 
                      + fL[i][nv] - fR[i][nv];
        Uhll[i][nv] *= scrh;
  
        Fhll[nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                   + SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        Fhll[nv] *= scrh;
      }
      Uhll[i][MXn] += (pL[i] - pR[i])*scrh;
      Fhll[MXn] += (SR[i]*pL[i] - SL[i]*pR[i])*scrh;

      #if SHOCK_FLATTENING == MULTID   

      /* ---------------------------------------------
         Switch to HLL in proximity of strong shocks.
        --------------------------------------------- */

       if (CheckZone(i, FLAG_HLL) || CheckZone(i+1, FLAG_HLL)){
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
         continue;
       }
      #endif

   /* ---- define total pressure, vB in left and right states ---- */

      pl = EXPAND(vL[BX1]*vL[BX1], + vL[BX2]*vL[BX2], + vL[BX3]*vL[BX3]);
      pr = EXPAND(vR[BX1]*vR[BX1], + vR[BX2]*vR[BX2], + vR[BX3]*vR[BX3]);

      pl = vL[PRS] + 0.5*pl;  
      pr = vR[PRS] + 0.5*pr;

      vBl = EXPAND(vL[VX1]*vL[BX1], + vL[VX2]*vL[BX2], + vL[VX3]*vL[BX3]);
      vBr = EXPAND(vR[VX1]*vR[BX1], + vR[VX2]*vR[BX2], + vR[VX3]*vR[BX3]);

      vxl = vL[VXn]; 
      vxr = vR[VXn];

   /* ----  magnetic field ---- */

      EXPAND(Bxs = Uhll[i][BXn];  ,
             Bys = Uhll[i][BXt];  ,
             Bzs = Uhll[i][BXb];)

   /* ---- normal velocity vx  ----  */

      vxs = Uhll[i][MXn]/Uhll[i][RHO];
      ps  = Fhll[MXn] + Bxs*Bxs - Fhll[RHO]*vxs;
/*
        ps = vL[RHO]*(SL[i] - vxl)*(vxs - vxl) + pl - vL[BXn]*vL[BXn] + Bxs*Bxs; 
*/
      vBs = EXPAND(Uhll[i][BX1]*Uhll[i][MX1], + 
                   Uhll[i][BX2]*Uhll[i][MX2], + 
                   Uhll[i][BX3]*Uhll[i][MX3]);

      vBs /= Uhll[i][RHO];

      usl[RHO] = uL[RHO]*(SL[i] - vxl)/(SL[i] - vxs);
      usr[RHO] = uR[RHO]*(SR[i] - vxr)/(SR[i] - vxs);

      usl[ENG] = (uL[ENG]*(SL[i] - vxl) + 
                 ps*vxs - pl*vxl - Bxs*vBs + vL[BXn]*vBl)/(SL[i] - vxs);
      usr[ENG] = (uR[ENG]*(SR[i] - vxr) + 
                 ps*vxs - pr*vxr - Bxs*vBs + vR[BXn]*vBr)/(SR[i] - vxs);

      EXPAND(usl[MXn] = usl[RHO]*vxs;
             usr[MXn] = usr[RHO]*vxs;        ,

             usl[MXt] = (uL[MXt]*(SL[i] - vxl) - (Bxs*Bys - vL[BXn]*vL[BXt]))/(SL[i] - vxs);
             usr[MXt] = (uR[MXt]*(SR[i] - vxr) - (Bxs*Bys - vR[BXn]*vR[BXt]))/(SR[i] - vxs); ,

             usl[MXb] = (uL[MXb]*(SL[i] - vxl) - (Bxs*Bzs - vL[BXn]*vL[BXb]))/(SL[i] - vxs);
             usr[MXb] = (uR[MXb]*(SR[i] - vxr) - (Bxs*Bzs - vR[BXn]*vR[BXb]))/(SR[i] - vxs);)

      EXPAND(usl[BXn] = usr[BXn] = Bxs;   ,
             usl[BXt] = usr[BXt] = Bys;   ,
             usl[BXb] = usr[BXb] = Bzs;)

      #ifdef GLM_MHD
       usl[PSI_GLM] = usr[PSI_GLM] = vL[PSI_GLM];
      #endif

      if (vxs >= 0.0){
        for (nv = 0; nv < NFLX; nv++) {
          state->flux[i][nv] = fL[i][nv] + SL[i]*(usl[nv] - uL[nv]);
        }
        state->press[i] = pL[i];
      } else {
        for (nv = 0; nv < NFLX; nv++) {
          state->flux[i][nv] = fR[i][nv] + SR[i]*(usr[nv] - uR[nv]);
        }
        state->press[i] = pR[i];
      }
    }
  }

/* -----------------------------------------------------
               initialize source term
   ----------------------------------------------------- */

  #if MHD_FORMULATION == EIGHT_WAVES
/*
   ROE_DIVB_SOURCE (state, beg, end, grid);
*/
/*
   for (i = beg; i <= end; i++) {
     uR = state->uR[i]; uL = state->uL[i];
     scrh = 1.0 / (SR[i] - SL[i]);
     for (nv = 0; nv < NFLX; nv++) {
       Uhll[i][nv] = SR[i]*uR[nv] - SL[i]*uL[nv] +
                     fL[i][nv] - fR[i][nv];
       Uhll[i][nv] *= scrh;
     }
     Uhll[i][MXn] += (pL[i] - pR[i])*scrh;
   }
*/
   HLL_DIVB_SOURCE (state, Uhll, beg + 1, end, grid);
  #endif
}
#elif EOS == ISOTHERMAL 

/* ************************************************************* */
void HLLC_Solver (const State_1D *state, int beg, int end, 
                  real *cmax, Grid *grid)
/*
 *
 *
 *
 *************************************************************** */
{
  print1 ("! HLLC solver not implemented for Isothermal EOS\n");
  print1 ("! Use hll or hlld instead.\n");
  QUIT_PLUTO(1);
}

#endif /* end #if on EOS */


