#include"pluto.h"

static void PPM_COEFF (Grid *, real **, real **, real **, real **, real **);

/* ********************************************************************** */
void PPM_COEFF (Grid *grid, real **a, real **b, real **c, 
                real **d, real **e)
/* 
 *
 * PURPOSE
 *
 *   Compute PPM interpolation coefficients for arbitrary 
 *   mesh size. Coefficients are stored in the 2D arrays 
 *   a, b, c, d, e. The first index refers to the dimension,
 *   the second to the grid. 
 *   Coefficients for the steepening are computed in the 
 *   STEEPEN function.
 *
 *
 *
 ************************************************************************ */
{
  int    dim, i, beg, end;
  real scrh, *dx;

  for (dim = 0; dim < DIMENSIONS; dim++){

  /* ----------------------------------------------
      Initialize the coefficients everywhere to 
      the default values for uniform grids.
     ---------------------------------------------- */

    for (i = 0; i < NMAX_POINT; i++){ 
      a[dim][i] = 0.5;
      b[dim][i] = c[dim][i] = 1.0/6.0;
      d[dim][i] = e[dim][i] = 0.5;
    }
    
  /* -----------------------------------------------
      for non-uniform grids, the coefficients will
      depends on the grid spacing and are computed
      on a slightly smaller stencil. 
      Since initialization is done only once at the
      beginning, this should be avoided in Chombo
     ----------------------------------------------- */

    #ifndef CH_SPACEDIM
     beg = grid[dim].lbeg - grid[dim].nghost + 1;
     end = grid[dim].lend + grid[dim].nghost - 2;
     dx  = grid[dim].dx;
     for (i = beg; i <= end; i++) {
 
       scrh = 1.0/(dx[i - 1] + dx[i] + dx[i + 1] + dx[i + 2]);
 
       b[dim][i] = dx[i + 1]*scrh*(dx[i + 1] + dx[i + 2])/(dx[i] + 2.0*dx[i + 1]);
       c[dim][i] = dx[i]*scrh*(dx[i - 1] + dx[i])/(2.0*dx[i] + dx[i + 1]);

       a[dim][i] = dx[i]/(dx[i] + dx[i + 1]) + 
              + 2.0*(dx[i + 1]*c[dim][i] - dx[i]*b[dim][i])/(dx[i] + dx[i + 1]);

       scrh = dx[i]/(dx[i - 1] + dx[i] + dx[i + 1]);
       d[dim][i] = scrh*(2.0*dx[i - 1] + dx[i])/(dx[i + 1] + dx[i]);
       e[dim][i] = scrh*(2.0*dx[i + 1] + dx[i])/(dx[i - 1] + dx[i]);
     }
    #endif
  }
}

/* ********************************************************************** */
void RECONSTRUCT (const State_1D *state, int beg, int end, Grid *grid)
/*
 *
 *
 *  PURPOSE
 *    
 *    provide piece-wise parabolic reconstruction inside each 
 *    cell. Notice that wl and wr are left and right states
 *    with respect to the INTERFACE, while am and ap (later
 *    in this function) refer to the cell center, that is:
 *
 *                    vL-> <-vR
 *      |--------*--------|--------*--------|
 *       <-vm   (i)   vp->       (i+1)
 *
 *
 ************************************************************************ */
{
  int   i, nv, k;
  real  scrh1, scrh2, scrh3, scrh4;
  real  steepness[NVAR];
  real  **vm, **vp, **vc, *a, *b, *c, *d, *e;
  static real  **dv, **dv_lim;
  static real **aa, **bb, **cc, **dd, **ee;

/* --------------------------------------------------- 
     PPM stencil is +- 2 zones:

     |---*---|---*---|---*---|---*---|---*---|
        i-2     i-1      i      i+1     i+2

     Thus, only 3 ghost zones are necessary.
     However, if FLATTENING or STEEPENING are added, 
     one more zone is required and the number of 
     ghost zones becomes 4.

     We define the leftmost and rightmost indexes  
     in the domain as 

      beg - s = 0
      end + s = grid[DIR].np_tot - 1

     where s is the required stencil.
   ---------------------------------------------------- */

  if (dv == NULL){
    dv     = Array_2D(NMAX_POINT, NVAR, double);
    dv_lim = Array_2D(NMAX_POINT, NVAR, double);
    aa     = Array_2D(DIMENSIONS, NMAX_POINT, double);
    bb     = Array_2D(DIMENSIONS, NMAX_POINT, double);
    cc     = Array_2D(DIMENSIONS, NMAX_POINT, double);
    dd     = Array_2D(DIMENSIONS, NMAX_POINT, double);
    ee     = Array_2D(DIMENSIONS, NMAX_POINT, double);
    PPM_COEFF (grid, aa, bb, cc, dd, ee);
  }

  for (nv = NVAR; nv--; ) steepness[nv] = 2.0;

  vm = state->vm;
  vp = state->vp;
  vc = state->v;

  a = aa[DIR];
  b = bb[DIR];
  c = cc[DIR];
  d = dd[DIR];
  e = ee[DIR];

/* -------------------------------------------------------
         Compute undivided differences
   ------------------------------------------------------- */

  for (i = beg - 2; i <= end + 1; i++) {
  for (nv = 0; nv < NVAR; nv++) {
    dv[i][nv] = vc[i + 1][nv] - vc[i][nv];
  }}

/* ---------------------------------------------------------
       STEP #1:   single value edge interpolation
   --------------------------------------------------------- */

double dp, dm, dc, d2;
  for (i = beg-1; i <= end+1; i++){
    for (nv = 0; nv < NVAR; nv++) {
      dp = dv[i][nv]; dm = dv[i-1][nv];
      if (dp*dm > 0.0){
        d2 = 2.0*ABS_MIN(dp,dm);
        dc = d[i]*dp + e[i]*dm;
        dv_lim[i][nv] = ABS_MIN(d2,dc);
      }else{
        dv_lim[i][nv] = 0.0;
      }
    }
  }

/*
   for (i = beg - 1; i <= end + 1; i++) {
   for (nv = 0; nv < NVAR; nv++){
     scrh1 = dv[i][nv]*dv[i - 1][nv];
     if (scrh1 > 0.0) {
       scrh2 = d[i]*dv[i][nv] + e[i]*dv[i - 1][nv];
       scrh3 = steepness[nv]*MIN(fabs(dv[i][nv]), fabs(dv[i - 1][nv]));
       scrh3 = MIN(fabs(scrh2), scrh3);
       dv_lim[i][nv] = scrh3*dsign(scrh2);
     }else{
       dv_lim[i][nv] = 0.0;
     }
   }}
*/

  for (i = beg - 1; i <= end; i++) {
  for (nv = 0; nv < NVAR; nv++){
    vp[i][nv] = vc[i][nv] + a[i]*dv[i][nv]  
                + b[i]*dv_lim[i][nv] - c[i]*dv_lim[i + 1][nv];
    vm[i + 1][nv] = vp[i][nv];
  }}


/* -------------------------------------------------------
        S T E P     # 3:  Parabolic Limiter
   ------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    for (nv = 0; nv < NVAR; nv++){
      dp = vp[i][nv] - vc[i][nv];
      dm = vm[i][nv] - vc[i][nv];

      if (dp*dm >= 0.0) {
        vp[i][nv] = vm[i][nv] = vc[i][nv];
// dp = dm = 0.0;
      }else{
        if (fabs(dp) >= 2.0*fabs(dm)) {
          vp[i][nv] = vc[i][nv] - 2.0*dm;
//          dp = - 2.0*dm;
        } else if (fabs(dm) >= 2.0*fabs(dp)) {
          vm[i][nv] = vc[i][nv] - 2.0*dp;
//          dm = - 2.0*dp;
        }
      }
//      vp[i][nv] = vc[i][nv] + dp;
//      vm[i][nv] = vc[i][nv] + dm;
    }

    #if PHYSICS == RHD || PHYSICS == RMHD
     RELATIVISTIC_LIM (vc[i], vp[i], vm[i]);
    #endif
  }

/* --------------------------------------------------------
            S T E P   # 4 : Flattening
   -------------------------------------------------------- */

  #if SHOCK_FLATTENING != NO
   FLATTEN (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][B1] = state->vL[i][B1] = state->bn[i];
   }
  #endif
}

