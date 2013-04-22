#include "pluto.h"

static real bmag;
static double Profile(double, int);
static void JETVAL (double *vjet);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int nv;
  double vjet[NVAR];

  us[RHO] = g_inputParam[ETA];
  us[PRS] = 1.0/g_gamma;

  bmag = sqrt(2.0*us[PRS]/g_inputParam[BETA]);

  EXPAND(us[VX1] = 0.0;  ,
         us[VX2] = 0.0;  ,
         us[VX3] = 0.0;)

  #if PHYSICS == MHD
   EXPAND(us[BX1] = 0.0;  ,
          us[BX2] = bmag; ,
          us[BX3] = 0.0;)

   us[AX1] = us[AX2] = 0.0;
   #if GEOMETRY == CARTESIAN
    us[AX3] = -x1*bmag;
   #elif GEOMETRY == CYLINDRICAL
    us[AX3] =  0.5*x1*bmag;
   #endif
  #endif

  if (x2 < -1.0){  /* -- if you need to put the jet inside at t=0 -- */
    JETVAL(vjet);
    for (nv = 0; nv < NVAR; nv++)
      us[nv] = us[nv] - (us[nv] - vjet[nv])*Profile(x1,nv);
  }
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 * 
 *
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *R, *Rp,  bxsout;
  double vjet[NVAR], vout[NVAR];

  if (side == X2_BEG){
    
    JETVAL(vjet);       /* -- beam/jet values -- */
    R  = grid[IDIR].xgc;  /* -- cylindrical radius -- */
 
    if (box->vpos == CENTER){  /* -- select cell-centered varaibles only -- */
      BOX_LOOP(box, k, j, i){  /* -- loop on boundary zones -- */
        for (nv = 0; nv < NVAR; nv++) vout[nv] = d->Vc[nv][k][2*JBEG-j-1][i];
        vout[VX2] *= -1.0;
        #if PHYSICS == MHD
         vout[BX1] *= -1.0;
        #endif
        for (nv = 0; nv < NVAR; nv++) /* -- smooth out the two solutions -- */
          d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(R[i],nv);
      }
    }else if (box->vpos == X1FACE){ /* -- select x1-staggered component -- */
      #ifdef STAGGERED_MHD
       Rp = grid[IDIR].A;           /* -- right interface area -- */
       BOX_LOOP(box, k, j, i){
         bxsout = -d->Vs[BX1s][k][2*JBEG - j - 1][i];
         d->Vs[BX1s][k][j][i] = bxsout*(1.0 - Profile(Rp[i],BX));
       }
     }
     #endif
  }
}

/* *********************************************** */
void JETVAL (double *vjet)
/* 
 *
 *  Assign jet values
 *
 ************************************************* */
{
  vjet[RHO] = 1.0;
  EXPAND( vjet[VX1] = 0.0;        ,
          vjet[VX2] = g_inputParam[MACH];  ,              
          vjet[VX3] = 0.0;)
  vjet[PRS] = 1.0/g_gamma;
  
  #if PHYSICS == MHD
   EXPAND( vjet[BX1] = 0.0;   ,
           vjet[BX2] = bmag;  ,              
           vjet[BX3] = 0.0;)
  #endif           

  #ifdef GLM_MHD
   vjet[PSI_GLM] = 0.0;
  #endif 
}
/* *********************************************** */
double Profile(double r, int nv)
/*
 *
 *  Set smooth profile between jet and ambient 
 *  in the r direction
 *
 ************************************************* */
{
  int m = 8;
  double r0=1.0;

  return 1.0/cosh(pow(r/r0,m));
}
