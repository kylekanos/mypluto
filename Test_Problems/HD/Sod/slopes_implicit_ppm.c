
#include"param.h"
void PSLOPES_BL (real **a0, real **al, real **ar, int is, int ie,
                 struct GRID *GG, struct advance *method, int INCLUDE_ORIGIN);

void TRI_DIAG (real *a, real *b, real *c, real *d,
		 real *u, int npoint);

/*       PPM PARAMETERS        */


#if PHYSICS == HD
 #define   K0      0.1
 #define   ETA1    20.0
 #define   ETA2    0.05
 #define   EPS1    0.01

 #define   EPS2     0.33
 #define   OME1     0.75
 #define   OME2     10.0

 #define CONTACT_STEEPENING YES
#endif

#if (PHYSICS == MHD || PHYSICS == SRMHD)

 #define   K0      0.1
 #define   ETA1    20.0
 #define   ETA2    0.05
 #define   EPS1    0.01

 #define   EPS2     0.33
 #define   OME1     0.75
 #define   OME2     10.0

 #define CONTACT_STEEPENING NO
#endif

#if PHYSICS == RHD
 #define   K0      1.0
 #define   ETA1    5.0
 #define   ETA2    0.05
 #define   EPS1    0.1

 #define   EPS2     1.0
 #define   OME1     0.52
 #define   OME2     10.0
 
 #define CONTACT_STEEPENING YES

#endif


/*  Define extra variables when interpolation 
    is done in the norm */

#define INTERPOLATE_NORM   NO



#if INTERPOLATE_NORM == YES
  #define N_TOT_SET  2
  #define INORM      0
#else
  #define N_TOT_SET  1
#endif
/* ###################################################################### */
void RECONSTRUCT (real **wc, real **wl, real **wr, int is, int ie,
                  struct GRID *GXYZ, struct advance *method)
/*
 #
 #           PIECEWISE PARABOLIC METHOD
 #
 #
 ######################################################################## */
{
  int   l_shock[NMAX], j, nv, sj, npoint, ngh;
  int   n_var_set, n_tot_set,nvars;
  real  aj_12, delta2[NMAX];
  real  ald, ard, f_t[NMAX];
  real  del_a[NMAX], rho2[NMAX], da[NMAX];
  real  d2press[NMAX];
  real  dummy1, dummy2, dummy3, dummy4;
  real  eta_t, etaj;
  real  delta_m[NMAX], fj[NMAX];
  real  press[NMAX], vel[NMAX], min_press[NMAX];
  real  gammak0;
  real  steepness[NVAR + N_EXTRA];
  real  norm[NMAX];
  real  **a, **al, **ar;
  static real **we, **wel, **wer;
  struct GRID *GG;

real c1[NMAX], c2[NMAX], c3[NMAX], rhs[NMAX], qp[NMAX];


  GG = GXYZ + NSWEEP;  

  npoint  = GG->np_tot;
  ngh     = GG->nghost;
  gammak0 = g_gamma*K0;

  for (nv = 0; nv < NVAR + N_EXTRA; nv++) {
    steepness[nv] = 2.0;
  }

  a  = wc;
  al = wl;
  ar = wr;
/*
  EXPAND (steepness[VX1] = 1.0;  ,
          steepness[VX2] = 1.0;  ,
          steepness[VX3] = 1.0;)

  steepness[INORM] = 1.0;  
  steepness[ISIN]  = 1.0;  
  steepness[ICOS]  = 1.0;
*/        
#if BLONDIN_LUFKIN == YES
  if (NSWEEP == 0) {
    if (fabs (GG->xr[3]) < 1.e-12) {
      PSLOPES_BL (a, al, ar, is, ie, GG, method, 1);
    } else {
      PSLOPES_BL (a, al, ar, is, ie, GG, method, 0);
    }
    return;
  }
#endif

  for (j = is-3; j <= ie+3; j++) {
    *(press + j) = a[j][PRS];
    *(vel + j)   = a[j][V1];
  }

/* -------------------------------------------------------
           RESETTING SOME QUANTITIES...
   ------------------------------------------------------- */

  for (j = is-3; j <= ie+3; j++) f_t[j] = 0.0;

  for (j = is-2; j <= ie+2; j++) {
    d2press[j]   = press[j + 1] - press[j - 1];
    min_press[j] = MIN (press[j + 1], press[j - 1]);
  }

/* ---------------------------------------------------------
       CHECK TO SEE IF ZONE J IS INSIDE A SHOCK
   --------------------------------------------------------- */

  for (j = is-1; j <= ie+1; j++) {
    dummy1 = fabs (d2press[j]) / min_press[j];
    l_shock[j] = (vel[j - 1] > vel[j + 1]);
    l_shock[j] = (l_shock[j] && (dummy1 > EPS2));

    f_t[j] = 0.0;
    if (l_shock[j]) {
      dummy2  = press[j + 2] - press[j - 2];
      dummy1  = d2press[j]/dummy2 - OME1;
      dummy1 *= OME2;
      f_t[j]  = MAX(0.0, dummy1);
      f_t[j]  = MIN(1.0, f_t[j]);
    }
  }

  for (j = is; j <= ie; j++) {
    sj = -1;
    if (d2press[j] < 0.0)
      sj = 1;
/*
         dummy2 = f_t[j - 1];
         dummy3 = f_t[j];
         dummy4 = f_t[j + 1];
         dummy2 = MAX(dummy2,dummy3);
         dummy2 = MAX(dummy2,dummy4);
         fj[j]  = MAX(0.0,MIN(0.5,dummy2));
*/
    fj[j] = MAX (f_t[j], f_t[j + sj]);
  }

/*      Relativistic limiter     */

#if PHYSICS == RHD && INTERPOLATE_NORM == YES
  if (we == NULL){ 
    we  = matrix(npoint, 1);
    wer = matrix(npoint, 1);
    wel = matrix(npoint, 1);
  }

  for (j = is-3; j <= ie+3; j++) {
    dummy1 = EXPAND(a[j][VX1]*a[j][VX1], + a[j][VX2]*a[j][VX2], + a[j][VX3]*a[j][VX3]);
    dummy1 = sqrt(dummy1);
    we[j][INORM] = dummy1;
    norm[j]      = dummy1;    
    dummy1       = 1.0/(dummy1 + 1.e-14);
    EXPAND(wc[j][VX1] *= dummy1;  ,
           wc[j][VX2] *= dummy1;  ,
           wc[j][VX3] *= dummy1;)
  }
#endif

/* ---------------------------------------------------------
                 MAIN LOOP ON VARIABLES
   ---------------------------------------------------------  */

  for (n_var_set = 1; n_var_set <= N_TOT_SET ; n_var_set++){

    if (n_var_set == 1) {
      nvars = NVAR + N_EXTRA;
      a     = wc;
      al    = wl;
      ar    = wr;
    }else{
      nvars = 1;
      a     = we;
      al    = wel;
      ar    = wer;
    }

    for (nv = 0; nv < nvars; nv++) {


for (j = 0; j < ie+2; j++){
  c1[j] = c3[j] = 1.0;
  c2[j] = 4.0;
  rhs[j] = 3.0*(a[j+1][nv] + a[j][nv]);
  qp[j]  = a[j][nv];
}  

TRI_DIAG(c1 + is, c2 + is, c3 + is, rhs + is, qp + is, ie - is);

/* -------------------------------------------------
            RESETTING SOME QUANTITIES...
   ------------------------------------------------- */

      for (j = is-2; j <= ie+3; j++) {
        delta_m[j] = 0.0;
        da[j] = a[j][nv] - a[j - 1][nv];
      }

      for (j = is-2; j <= ie+2; j++) {
        del_a[j] = GG->c6[j]*da[j + 1] + GG->c7[j]*da[j];
      }

/* ---------------------------------------------------------
                            step #1
   --------------------------------------------------------- */

      for (j = is-2 ; j <= ie+2 ; j++) {
        dummy1 = da[j + 1]*da[j];
        if (dummy1 > 0.0) {
          dummy2 = steepness[nv]*MIN(fabs(da[j]), fabs(da[j + 1]));
          dummy2 = MIN(fabs(del_a[j]), dummy2);
          delta_m[j] = dummy2*dsign(del_a[j]);
        }
delta_m[j] = del_a[j];
      }

      for (j = is-2; j <= ie+1 ; j++) {
        aj_12 = a[j][nv] + GG->c1[j]*da[j + 1] + GG->c2[j]*
               (GG->c3[j]*da[j + 1] - GG->c4[j]*delta_m[j + 1] +
                GG->c5[j]*delta_m[j]);

//aj_12 = qp[j];
aj_12 = MIN(aj_12, MAX(a[j+1][nv], a[j][nv]));
aj_12 = MAX(aj_12, MIN(a[j+1][nv], a[j][nv]));

        ar[j][nv] = al[j + 1][nv] = aj_12;
      }


#if CONTACT_STEEPENING == YES
/* ---------------------------------------------------
          step # 2 :    CONTACT STEEPENING
   ---------------------------------------------------  */

      if (nv == DN && n_var_set == 1) {

/* -----------------------------------------------------
           APPLY ONLY TO DENSITY (nv = DN)
   ----------------------------------------------------- */

        for (j = is-2; j <= ie+2; j++) {
          delta2[j] = GG->c9[j] * da[j + 1] - GG->c10[j] * da[j];
          rho2[j] = a[j + 1][RHO] - a[j - 1][RHO];
        }

        for (j = is-1; j <= ie+1; j++) {
          dummy1 = gammak0*fabs(rho2[j]) /
                   MIN (a[j + 1][RHO], a[j - 1][RHO]);
          dummy2 = fabs (d2press[j]) / min_press[j];

/*     CHECK TO SEE IF ZONE J IS INSIDE A CONTACT DISC.  */

          if (dummy1 > dummy2) {
            dummy3 = -delta2[j + 1] * delta2[j - 1];
            dummy4 = fabs (rho2[j]) - EPS1 *
            MIN (fabs (a[j + 1][RHO]), fabs (a[j - 1][RHO]));


          /*    get ETA_TILDE */

            eta_t = 0.0;
            if ((dummy3 > 0.0) && (dummy4 > 0.0)) {
              eta_t = (delta2[j + 1] - delta2[j - 1]) / rho2[j];
              eta_t *= -GG->c8[j];
            }

  /*        get ETA       */

            dummy3 = MIN (ETA1 * (eta_t - ETA2), 1.0);
            etaj = MAX (0.0, dummy3);

            ald = a[j - 1][RHO] + 0.5 * delta_m[j - 1];
            ard = a[j + 1][RHO] - 0.5 * delta_m[j + 1];

            dummy3 = 1.0 - etaj;
            al[j][RHO] = al[j][RHO] * dummy3 + ald * etaj;
            ar[j][RHO] = ar[j][RHO] * dummy3 + ard * etaj;
          }
        }
      }
#endif

/* --------------------------------------------------------
                       S T E P   # 3
   -------------------------------------------------------- */

      for (j = is; j <= ie; j++) {
        if (l_shock[j]) {
          dummy1 = a[j][nv] * fj[j];
          dummy2 = 1.0 - fj[j];
          al[j][nv] = dummy1 + al[j][nv]*dummy2;
          ar[j][nv] = dummy1 + ar[j][nv]*dummy2;
        }
      }

/* -------------------------------------------------------
                  S T E P     # 4
   ------------------------------------------------------- */

      for (j = is; j <= ie; j++) {
        dummy1 = (ar[j][nv] - a[j][nv])*(a[j][nv] - al[j][nv]);
        if (dummy1 <= 0.0) {
          al[j][nv] = ar[j][nv] = a[j][nv];
        }
      }

      for (j = is; j <= ie; j++) {
        dummy2 = ar[j][nv] - al[j][nv];
        dummy3 = dummy2*(a[j][nv] - 0.5*(al[j][nv] + ar[j][nv]));
        dummy4 = dummy2*dummy2 / 6.0;

        if (dummy3 > dummy4)  al[j][nv] = 3.0*a[j][nv] - 2.0*ar[j][nv];
        if (-dummy4 > dummy3) ar[j][nv] = 3.0*a[j][nv] - 2.0*al[j][nv];
      }
    }
    for (j = is-1; j <= ie; j++) {
    for (nv = 0; nv < nvars; nv++) {
      al[j][nv] = ar[j][nv];
      ar[j][nv] = al[j + 1][nv];
    }}
  }


#if PHYSICS == RHD && INTERPOLATE_NORM == YES

/*  RELATIVISTIC LIMITER   */

  for (j = is-3; j <= ie+3; j++) {

    EXPAND(wc[j][VX1] *= norm[j];  ,
           wc[j][VX2] *= norm[j];  ,
           wc[j][VX3] *= norm[j];)
  }

  for (j = is-1; j <= ie; j++) {

    EXPAND(wl[j][VX1] *= wel[j][INORM];  ,
           wl[j][VX2] *= wel[j][INORM];  ,
           wl[j][VX3] *= wel[j][INORM];)

    EXPAND(wr[j][VX1] *= wer[j][INORM];  ,
           wr[j][VX2] *= wer[j][INORM];  ,
           wr[j][VX3] *= wer[j][INORM];)
  }
#endif

/*
  for (j = is-1; j <= ie; j++) {
  for (nv = 0; nv < NVAR; nv++) {
    wl[j][nv] = al[j][nv];
    wr[j][nv] = ar[j][nv];
  }}
*/
/* State symmetry check  
if (NSWEEP == 0){
for (nv = 0; nv<NVAR;nv++){ 
 sj = (nv == VX ? -1.0:1.0);
 if ( fabs(al[3][nv] - sj*ar[3][nv]) > 1.e-2 ||
      fabs(ar[2][nv] - sj*al[4][nv]) > 1.e-2       ){
   print ("no symm in var #%d  at %d \n",nv,j);
   for (j = 1; j <=6 ; j++){
     print ("%d  %d  %12.8e  %12.8e   %12.8e\n",nv, j,a[j][nv],al[j][nv],ar[j][nv]);
   }
   exit(1);
 }
}}
*/ 
  
#if PHYSICS == RHD22      /*  Check whether velocities are > 1 */
  for (j = is-1 ; j <= ie ; j++){

/* if the norm of the reconstructed velocity vector does not
   fall in between the point-wise norm to the left or 
   to the right, re-scale the velocity so that norm is equal 
   to the maximum between j and j+1 */

    dummy1 = EXPAND (al[j][VX1]*al[j][VX1],
                    +al[j][VX2]*al[j][VX2],
                    +al[j][VX3]*al[j][VX3]);

    dummy2 = EXPAND (ar[j][VX1]*ar[j][VX1],
                    +ar[j][VX2]*ar[j][VX2],
                    +ar[j][VX3]*ar[j][VX3]);

    dummy3 = MAX(norm[j], norm[j + 1]) + 1.e-12;

    if (dummy1 > dummy3) {
/*      
      WARNING(print ("Correcting left norm...\n");
              WHERE(j, GXYZ);)
*/      
      EXPAND (al[j][VX1] *= sqrt(dummy3/dummy1);  ,
              al[j][VX2] *= sqrt(dummy3/dummy1);  ,
              al[j][VX3] *= sqrt(dummy3/dummy1));
    }

    if (dummy2 > dummy3) {
/*      
      WARNING(print ("Correcting right norm...\n");
              WHERE(j, GXYZ);)
*/      
      EXPAND (ar[j][VX1] *= sqrt(dummy3/dummy2);  ,
              ar[j][VX2] *= sqrt(dummy3/dummy2);  ,
              ar[j][VX3] *= sqrt(dummy3/dummy2));

    }
  }
#endif
  return;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#if BLONDIN_LUFKIN == YES
/* ---------------------------------------------------------------------- */
void PSLOPES_BL (real **a0, real **al, real **ar, int is, int ie,
                 struct GRID *GG, struct advance *method, int INCLUDE_ORIGIN)


/* ---------------------------------------------------------------------
    PIECEWISE PARABOLIC METHOD in orthogonal curvilinear coordinates.
    It follows the general prescription given by Blondin & Lufkin (1993);

    Different switches might be required according to the combination of
    geometries and coordinate singularities; however this function requires
    that only one coordinate to be singular for a given geometry: 


     CARTESIAN    - 1,2,3D --> Regular  

     CYLINDRICAL  - 1D     --> r possible sing.
                  
                  - 2D     --> r possible sing. ; z cartesian
                  - 3D     --> r regular        ; z, phi are cartesian

     SPHERICAL    - 1D     --> r possible sing.
                  - 2D     --> r regular, theta possible sing
                  - 3D     --> r,theta,phi regular.


    Therefore the only allowed singularities are when 
  
     r-CYLINDRICAL   when  DIMENSIONS = 1,2
     r-SPHERICAL     when  DIMENSIONS = 1
     theta-SPHERICAL when  DIMENSIONS = 2
      
   --------------------------------------------------------------------- */
{
  real delta2[NMAX];
  real ald, ard, f_t[NMAX];
  real del_a[NMAX], rho2[NMAX], da[NMAX];
  real d2press[NMAX];
  real dummy1, dummy2, dummy3, dummy4;
  real eta_t, etaj;
  real delta_m[NMAX], fj[NMAX];
  real press[NMAX], vel[NMAX], min_press[NMAX];
  real gammak0;

  int l_shock[NMAX];
  int j, nv, nvt, sj;
  int npoint, ngh;

  real norm[NMAX];
  real **a, **aj_12;
  real *F, *G, *H, *a6, xR, xL, dx;
  real *AL_flat_R, *AR_flat_L, *A0_flat_L, *A0_flat_R;
  real Jx[NMAX], x;

  npoint = GG->np_tot;
  ngh = GG->nghost;
  gammak0 = g_gamma * K0;


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    Blondin & Lufkin  
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

  a = mbuffer[SLOPES_MBUFFER + 0];
  aj_12 = mbuffer[SLOPES_MBUFFER + 1];

  F = vbuffer[SLOPES_VBUFFER + 0];
  G = vbuffer[SLOPES_VBUFFER + 1];
  H = vbuffer[SLOPES_VBUFFER + 2];
  a6 = vbuffer[SLOPES_VBUFFER + 3];
  A0_flat_L = vbuffer[SLOPES_VBUFFER + 4];
  A0_flat_R = vbuffer[SLOPES_VBUFFER + 5];
  AR_flat_L = vbuffer[SLOPES_VBUFFER + 6];
  AL_flat_R = vbuffer[SLOPES_VBUFFER + 7];

  for (j = 0; j < npoint; j++) {
    dummy4 = 1.0;

    dx = GG->dx[j];
    x = GG->x[j];
#if GEOMETRY == SPHERICAL
    if (NSWEEP == 0)
      dummy1 = x * x + dx * dx / 12.0;
    if (NSWEEP == 1)
      dummy1 = sin (x) * (1.0 - dx * dx / 24.0 + dx * dx * dx * dx / 1920.);
#elif GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    if (NSWEEP == 0)
      dummy1 = x;
    if (NSWEEP == 1)
      dummy1 = 1.0;
#endif
    for (nv = 0; nv < NVAR; nv++) {
      a[j][nv] = a0[j][nv] * dummy1;
    }
  }

/* Define the flattening coefficient for the 
   parabola, and geometry dependent quantities   */

#if GEOMETRY == CARTESIAN
  for (j = 1; j < npoint - 1; j++) {
    A0_flat_L[j] = 3.0;
    AR_flat_L[j] = -2.0;

    A0_flat_R[j] = 3.0;
    AL_flat_R[j] = -2.0;
  }
#elif GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
  if (NSWEEP == 0) {
    for (j = 1; j < npoint - 1; j++) {
      xL = GG->xr[j - 1];
      dx = GG->dx[j];
      Jx[j - 1] = xL;
      F[j] = dx / (6.0 * xL + 3.0 * dx);

      A0_flat_L[j] = 6.0 / (2.0 - 3.0 * F[j]);
      AR_flat_L[j] = -(4.0 + 3.0 * F[j]) / (2.0 - 3.0 * F[j]);

      A0_flat_R[j] = 6.0 / (2.0 + 3.0 * F[j]);
      AL_flat_R[j] = -(4.0 - 3.0 * F[j]) / (2.0 + 3.0 * F[j]);
    }
  } else if (NSWEEP == 1) {
    for (j = 1; j < npoint - 1; j++) {
      Jx[j] = 1.0;

      A0_flat_L[j] = 3.0;
      AR_flat_L[j] = -2.0;

      A0_flat_R[j] = 3.0;
      AL_flat_R[j] = -2.0;
    }
  } else if (NSWEEP == 2) {
    print ("PPM phi- coord not set \n");
    exit (1);
  }
#elif GEOMETRY == SPHERICAL
  if (NSWEEP == 0) {
    for (j = 1; j < npoint - 1; j++) {
      xL = GG->xr[j - 1];
      dx = GG->dx[j];
      Jx[j - 1] = xL * xL;
      G[j] =
        (3.0 * xL * (xL + dx) + 0.9 * dx * dx) / (3.0 * xL * (xL + dx) +
                                                  dx * dx);
      H[j] = dx * (2.0 * xL + dx) / (6.0 * xL * (xL + dx) + 2.0 * dx * dx);

      A0_flat_L[j] = 6.0 / (3.0 - G[j] - 3.0 * H[j]);
      AR_flat_L[j] = -(3.0 + G[j] + 3.0 * H[j]) / (3.0 - G[j] - 3.0 * H[j]);

      A0_flat_R[j] = 6.0 / (3.0 - G[j] + 3.0 * H[j]);
      AL_flat_R[j] = -(3.0 + G[j] - 3.0 * H[j]) / (3.0 - G[j] + 3.0 * H[j]);

    }
  } else if (NSWEEP == 1) {
    xL = GG->xr[0];
    for (j = 1; j < npoint - 1; j++) {
      xR = GG->xr[j];
      dx = GG->dx[j];
      Jx[j] = sin (xR);
      dummy1 = cos (xL) - cos (xR);
      F[j] = (sin (xR) - sin (xL)) / dx;
      G[j] = 2.0 * dummy1 / (dx * dx) * (1.0 - 0.5 * dx / tan (0.5 * dx));

      A0_flat_L[j] = dummy1 / (cos (xL) - F[j] - G[j]);
      AR_flat_L[j] = (cos (xR) - F[j] - G[j]) / (cos (xL) - F[j] - G[j]);

      A0_flat_R[j] = -dummy1 / (cos (xR) - F[j] + G[j]);
      AL_flat_R[j] = (cos (xL) - F[j] + G[j]) / (cos (xR) - F[j] + G[j]);

      xL = xR;
    }
  } else if (NSWEEP == 2) {
    print ("PPM interpolation not set for Spherical phi - coord !\n");
    exit (1);
  }
#endif


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

  for (j = 0; j < npoint; j++) {
    *(press + j) = a0[j][PRS];
    *(vel + j) = a0[j][V1];
#if PHYSICS == RHD
    norm[j] = EXPAND (a0[j][VX1] * a0[j][VX1],
                      +a0[j][VX2] * a0[j][VX2], +a0[j][VX3] * a0[j][VX3]);
#endif
  }

/* -------------------------------------------------------
             RESETTING SOME QUANTITIES...
   ------------------------------------------------------- */

  for (j = 0; j < npoint; j++)
    f_t[j] = 0.0;

  for (j = 1; j < npoint - 1; j++) {
    d2press[j] = press[j + 1] - press[j - 1];
    min_press[j] = MIN (press[j + 1], press[j - 1]);
  }

/* ---------------------------------------------------------
          CHECK TO SEE IF ZONE J IS INSIDE A SHOCK
   --------------------------------------------------------- */

  for (j = 2; j < npoint - 2; j++) {
    dummy1 = fabs (d2press[j]) / min_press[j];
    l_shock[j] = (vel[j - 1] > vel[j + 1]);
    l_shock[j] = (l_shock[j] && (dummy1 > EPS2));

    f_t[j] = 0.0;
    if (l_shock[j]) {
      dummy2 = press[j + 2] - press[j - 2];
      dummy1 = d2press[j] / dummy2 - OME1;
      dummy1 *= OME2;
      f_t[j] = MAX (0.0, dummy1);
      f_t[j] = MIN (1.0, f_t[j]);
    }
  }

  for (j = 3; j < npoint - 3; j++) {
    sj = -1;
    if (d2press[j] < 0.0)
      sj = 1;
/*
         dummy2 = f_t[j - 1];
         dummy3 = f_t[j];
         dummy4 = f_t[j + 1];
         dummy2 = MAX(dummy2,dummy3);
         dummy2 = MAX(dummy2,dummy4);
         fj[j]  = MAX(0.0,MIN(0.5,dummy2));
*/
    fj[j] = MAX (f_t[j], f_t[j + sj]);
  }

/* ---------------------------------------------------------
                 MAIN LOOP ON VARIABLES
   ---------------------------------------------------------  */

  for (nv = 0; nv < NVAR; nv++) {

    nvt = MIN (nv, TRAC);

/* -------------------------------------------------
           RESETTING SOME QUANTITIES...
   ------------------------------------------------- */

    for (j = 1; j < npoint; j++) {
      delta_m[j] = 0.0;
      da[j] = a[j][nv] - a[j - 1][nv];
    }

    for (j = 1; j < npoint - 1; j++) {
      del_a[j] = GG->c6[j] * da[j + 1] + GG->c7[j] * da[j];
    }

/* ---------------------------------------------------------
                       step #1
   --------------------------------------------------------- */

    for (j = 1; j < npoint - 1; j++) {

      dummy1 = da[j + 1] * da[j];
      if (dummy1 > 0.0) {
        dummy2 = 2.0 * MIN (fabs (da[j]), fabs (da[j + 1]));
        dummy2 = MIN (fabs (del_a[j]), dummy2);
        delta_m[j] = dummy2 * dsign (del_a[j]);
      }

      delta_m[j] = del_a[j];
    }

    for (j = 1; j < npoint - 2; j++) {
      aj_12[j][nv] = a[j][nv] + GG->c1[j] * da[j + 1] + GG->c2[j] *
        (GG->c3[j] * da[j + 1] - GG->c4[j] * delta_m[j + 1] +
         GG->c5[j] * delta_m[j]);
    }

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    Blondin & Lufkin  
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    if (INCLUDE_ORIGIN == YES) {
      if (fabs (Jx[3]) < 1.e-18)
        Jx[3] = 1.0;
    }
    for (j = 1; j < npoint - 2; j++) {
      aj_12[j][nv] /= Jx[j];
    }

    if (INCLUDE_ORIGIN == 1) {
#if GEOMETRY == SPHERICAL
      if (NSWEEP == 0) {
        aj_12[3][nv] = 2.5 * a0[4][nv] - 1.5 * aj_12[4][nv];
        aj_12[3][VX1] = 0.0;
      } else if (NSWEEP == 1) {
        dx = GG->dx[4];
        aj_12[3][nv] = a0[4][nv] * (1.0 - cos (dx)) +
          aj_12[4][nv] * (G[4] + cos (dx) - F[4]);

        aj_12[3][nv] /= G[4] + 1.0 - F[4];
        aj_12[3][VX2] = 0.0;
      }
#endif

#if GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
      if (NSWEEP == 0) {
        aj_12[3][nv] = 2.0 * a0[4][nv] - aj_12[4][nv];
        aj_12[3][VX1] = 0.0;
      } else if (NSWEEP == 1) {
        aj_12[3][nv] = 1.5 * a0[4][nv] - 0.5 * aj_12[4][nv];
        aj_12[3][VX1] = 0.0;
      }
#endif
    }


    for (j = 2; j < npoint - 2; j++) {
      aj_12[j][nv] = MIN (aj_12[j][nv], MAX (a0[j][nv], a0[j + 1][nv]));
      aj_12[j][nv] = MAX (aj_12[j][nv], MIN (a0[j][nv], a0[j + 1][nv]));
    }

    for (j = 2; j < npoint - 2; j++) {
      ar[j][nv] = al[j + 1][nv] = aj_12[j][nv];
    }
/* --------------------------------------------------------
                    S T E P   # 3
   -------------------------------------------------------- */


    for (j = 3; j < npoint - 3; j++) {
      if (l_shock[j]) {
        dummy1 = a0[j][nv] * fj[j];
        dummy2 = 1.0 - fj[j];
        al[j][nv] = dummy1 + al[j][nv] * dummy2;
        ar[j][nv] = dummy1 + ar[j][nv] * dummy2;
      }
    }

/* -------------------------------------------------------
                      S T E P    # 4
   ------------------------------------------------------- */

    for (j = 4; j < npoint - 3; j++) {

      dummy1 = (ar[j][nv] - a0[j][nv]) * (a0[j][nv] - al[j][nv]);
      if (dummy1 <= 0.0) {
        if (j > 4)
          al[j][nv] = a0[j][nv];
        ar[j][nv] = a0[j][nv];
      }
    }

    /*    define  a6, geometry dependent coefficient   */

#if GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    if (NSWEEP == 0) {
      for (j = 3; j < npoint - 3; j++) {
        a6[j] = 6.0 * (a0[j][nv] - 0.5 * (al[j][nv] + ar[j][nv])
                       + 0.5 * (al[j][nv] - ar[j][nv]) * F[j]);
      }
    } else if (NSWEEP == 1) {
      for (j = 3; j < npoint - 3; j++) {
        a6[j] = 6.0 * (a0[j][nv] - 0.5 * (al[j][nv] + ar[j][nv]));
      }
    }
#elif GEOMETRY == SPHERICAL
    if (NSWEEP == 0) {
      for (j = 3; j < npoint - 3; j++) {
        a6[j] = 6.0 * (a0[j][nv] - 0.5 * (al[j][nv] + ar[j][nv])
                       + 0.5 * (al[j][nv] - ar[j][nv]) * H[j]) / G[j];
      }
    } else if (NSWEEP == 1) {
      xL = GG->xr[ngh - 1];
      for (j = 3; j < npoint - 3; j++) {
        xR = GG->xr[j];
        dummy1 = cos (xL) - cos (xR);
        a6[j] = (a0[j][nv] * dummy1 + ar[j][nv] * (cos (xR) - F[j])
                 - al[j][nv] * (cos (xL) - F[j])) / G[j];
        xL = xR;
      }
    }
#elif GEOMETRY == CARTESIAN
    for (j = 3; j < npoint - 3; j++) {
      a6[j] = 6.0 * (a0[j][nv] - 0.5 * (al[j][nv] + ar[j][nv]));
    }
#endif

/*       flat parabola         */

    for (j = 3 + INCLUDE_ORIGIN; j < npoint - 3; j++) {
      dummy2 = ar[j][nv] - al[j][nv];
      dummy3 = dummy2 * a6[j];
      dummy4 = dummy2 * dummy2;

      if (dummy3 > dummy4 && j > 4) {
        al[j][nv] = A0_flat_L[j] * a0[j][nv] + AR_flat_L[j] * ar[j][nv];
      }

      if (-dummy4 > dummy3) {
        ar[j][nv] = A0_flat_R[j] * a0[j][nv] + AL_flat_R[j] * al[j][nv];
      }
    }
  }

  for (j = 3; j < npoint - 3; j++) {
    for (nv = 0; nv < NVAR; nv++) {
      al[j][nv] = ar[j][nv];
      ar[j][nv] = al[j + 1][nv];
    }
#if PHYSICS == RHD              /*  Check whether velocities are > 1 */

    dummy1 = EXPAND (al[j][VX1] * al[j][VX1],
                     +al[j][VX2] * al[j][VX2], +al[j][VX3] * al[j][VX3]);

    dummy2 = EXPAND (ar[j][VX1] * ar[j][VX1],
                     +ar[j][VX2] * ar[j][VX2], +ar[j][VX3] * ar[j][VX3]);

    dummy3 = MAX (norm[j], norm[j + 1]) + 1.e-12;

    if (dummy1 > dummy3) {
      EXPAND (al[j][VX1] *= sqrt (dummy3 / dummy1);
              , al[j][VX2] *= sqrt (dummy3 / dummy1);
              , al[j][VX3] *= sqrt (dummy3 / dummy1));

    }
    if (dummy2 > dummy3) {

      EXPAND (ar[j][VX1] *= sqrt (dummy3 / dummy2);
              , ar[j][VX2] *= sqrt (dummy3 / dummy2);
              , ar[j][VX3] *= sqrt (dummy3 / dummy2));

    }
#endif
  }

/*
SHOW(al,3);
SHOW(ar,3);
print ("----------------- \n");
*/
  for (j = 3; j < npoint - 3; j++) {
    if (al[j][PRS] < 0.0 || ar[j][PRS] < 0.0) {
      print ("p < 0 in PSLOPES   ");
/*
      WHERE (j);
*/
      print ("pL   %12.6e  pr  %12.6e\n", al[j][PRS], ar[j][PRS]);
      for (j = 3; j < npoint - 3; j++) {
        print ("%d  %12.6e  %12.6e \n", j, al[j][PRS], ar[j][PRS]);
      }
      exit (1);

      al[j][PRS] = ar[j][PRS] =
        MAX (0.5 * (al[j - 1][PRS] + ar[j + 1][PRS]), 1.e-5);

    }
    if (al[j][RHO] < 0.0 || ar[j][RHO] < 0.0) {
      print ("rho < 0 in PSLOPES\n");
      print ("NSWEEP %d  j %d   al %12.6e      ar %12.6e\n",
              NSWEEP, j, al[j][RHO], ar[j][RHO]);
    }
  }
  return;
}
#endif

#undef   K0
#undef   ETA1
#undef   ETA2
#undef   EPS1

#undef   EPS2
#undef   OME1
#undef   OME2

#undef  CONTACT_STEEPENING 
#undef  INTERPOLATE_NORM

void TRI_DIAG (real *a, real *b, real *c, real *d,
		 real *u, int npoint)
/*

   SOLVE A TRI-DIAGONAL SYSTEM of THE FORM 

      A_j x_(j-1) + B_j x_j + C_j x_(j+1) = D_j

   WITH GIVEN (FIXED) BOUNDARY CONDITIONS in U[0], U[NPOINT-1]
*/
{
  int i;
  real *e, *f, den;


  e = vector (npoint);
  f = vector (npoint);

  e[0] = 0.0;
  f[0] = u[0];

  e[0] = 1.0;
  f[0] = 0.0;

  for (i = 1; i < npoint - 1; i++) {
    den = b[i] + a[i] * e[i - 1];
    e[i] = -c[i] / den;
    f[i] = (d[i] - a[i] * f[i - 1]) / den;
  }


  for (i = npoint - 2; i > 0; i--)
    u[i] = e[i] * u[i + 1] + f[i];

  free_vector (e);
  free_vector (f);

}
