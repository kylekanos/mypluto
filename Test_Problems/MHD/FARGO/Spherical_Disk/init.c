#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *  Setup a magnetized disk in 3D spherical coordinates.
 *  The configuration is similar to
 *
 *  "Turbulence and Steady Flows in 3D Global Stratified
 *   MHD Simulations of Accretion Disks"
 *
 *  Flock et al., ApJ (2011), 735:122
 *  
 *  The main difference is that a polodial magnetic field is
 *  used here and the boundary conditions have been simplified
 *  to outflow (or reflective) in the radial coordinates while
 *  a time-independent boundary is used in theta.
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double rnd, s, c, R, r, z;
  double H, cs, Bphi;

/*
  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;
    print1 (">>> suggested number of zones in x2:  N1 x r0 x %f\n",
            (DOM_XEND[JDIR]-DOM_XBEG[JDIR])/(DOM_XEND[IDIR]-DOM_XBEG[IDIR]));
    print1 (">>> suggested number of zones in x3:  N1 x r0 x %f\n",
            (DOM_XEND[KDIR]-DOM_XBEG[KDIR])/(DOM_XEND[IDIR]-DOM_XBEG[IDIR]));
  }
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);

*/
  rnd = 0.0;

  g_gamma = 1.0001;
  r   = x1;
  s   = sin(x2);  
  c   = cos(x2);

  R  = r*s;  
  z  = r*c;
  cs = g_inputParam[H_R]/sqrt(R);
  H  = g_inputParam[H_R]*R;

  us[RHO]    = exp((s-1.0)/(g_inputParam[H_R]*g_inputParam[H_R]))/(R*sqrt(R));
  us[VX1]    = us[VX2] = us[VX3] = 0.0;
  us[iVPHI] = 1.0/sqrt(r)*(1.0 - 2.5*g_inputParam[H_R]*g_inputParam[H_R]/s);

  us[PRS] = cs*cs*us[RHO];

  #if PHYSICS == MHD
   us[BX1] = us[BX2] = us[BX3] = 0.0;
   us[AX1] = us[AX2] = us[AX3] = 0.0;

   us[AX3]  = (sin(2.0*CONST_PI*R) - R*cos(2.0*CONST_PI*R))/R;
   us[AX3] *= 3.5e-4*exp(-pow(z/H,4.0))*exp(-pow((R-6.0)/2.0,4.0));
   if (fabs(z/H) > 2.5 || R < 2.5) us[AX3] = 0.0;
   us[PRS] *= 1.0 + 0.1*rnd; 
  #endif
  g_smallPressure = 1.e-9;
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int    i,j,k,nv;
  static int first_call = 1;
  double v[NVAR];
  double *r, *th, *dVr, *dVth, *dVphi, vol;
  double pm, kin;
  double pmtot, Tmtot;
  static double voltot;
  FILE *fp;

  r   = grid[IDIR].x;
  th  = grid[JDIR].x;
  dVr = grid[IDIR].dV; dVth = grid[JDIR].dV; dVphi = grid[KDIR].dV;

/* ---------------------------------------------------------
         compute total volume once at the beginning
   --------------------------------------------------------- */

  if (first_call){
    voltot = 0.0;
    DOM_LOOP(k,j,i) voltot += dVr[i]*dVth[j]*dVphi[k];
    #ifdef PARALLEL
     MPI_Allreduce (&voltot, &kin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Barrier (MPI_COMM_WORLD);
     voltot = kin;
    #endif
  }

/* --------------------------------------------------------
              compute volume averages 
   -------------------------------------------------------- */

  pmtot = Tmtot = 0.0;
  DOM_LOOP(k,j,i){
    vol = dVr[i]*dVth[j]*dVphi[k]/voltot;

    for (nv = NVAR; nv--;   ) v[nv] = d->Vc[nv][k][j][i];

    #if PHYSICS == MHD
     pm     = 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
     pmtot += pm*vol;
     Tmtot += v[BX1]*v[BX3]*vol;
    #endif
  }
 
  #ifdef PARALLEL
   MPI_Allreduce (&pmtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pmtot = vol;

   MPI_Allreduce (&Tmtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Tmtot = vol;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* -----------------------------------------------------
    Write "averages.dat" to disk. 
    The file is created as new if this is the very
    initial time step. Data is appended otherwise.
    Processor #0 does the actual writing.
   ----------------------------------------------------- */

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){
      fp = fopen("averages.dat","w");
      fprintf (fp,"#%s %12s  %12s  %12s\n", 
               "      t","dt","<B^2>/2","<Tm>");
    }else{
      if (tpos < 0.0){  /* obtain time coordinate of to last entry */
        char   sline[512];
        fp = fopen("averages.dat","r");
        while (fgets(sline, 512, fp))  {}
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen("averages.dat","a");
    }
    if (g_time > tpos){
      fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n",
               g_time, g_dt, pmtot, Tmtot);
    }
    fclose(fp);
  }
  first_call = 0; 
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double cs, qx, qy, qz, s, v[256];
  double *r, *th, *phi;
 
  r   = grid[IDIR].x;
  th  = grid[JDIR].x;
  phi = grid[KDIR].x;

  if (side == 0) {    /* -- check solution inside domain -- */
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= d->Vc[nv][k][j][IBEG];
        s = sin(th[j]);
        d->Vc[iVPHI][k][j][i] = 1.0/sqrt(r[i])*(1.0 - 2.5*g_inputParam[H_R]*g_inputParam[H_R]/s);

        if (d->Vc[iVR][k][j][i] > 0.0) d->Vc[iVR][k][j][i] = 0.0;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];  
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IBEG];  
      #endif
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= d->Vc[nv][k][j][IEND];
        s = sin(th[j]);
        d->Vc[iVPHI][k][j][i] = 1.0/sqrt(r[i])*(1.0 - 2.5*g_inputParam[H_R]*g_inputParam[H_R]/s);
        d->Vc[iVR][k][j][i]   = -d->Vc[iVR][k][j][2*IEND - i + 1];
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i)  d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IEND];
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i)  d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IEND];
      #endif
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        Init (v, r[i], th[j], phi[k]);
        cs = g_inputParam[H_R]/sqrt(r[i]*sin(th[j]));
        v[PRS] = v[RHO]*cs*cs; /* pressure is re-defined without perturbation */
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= v[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
   }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
      #endif
    }

  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        Init (v, r[i], th[j], phi[k]);
        cs = g_inputParam[H_R]/sqrt(r[i]*sin(th[j]));
        v[PRS] = v[RHO]*cs*cs; /* pressure is re-defined without perturbation */
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= v[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
      #endif
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = -1.0/(x1*x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  return 0.0;
}
#endif
