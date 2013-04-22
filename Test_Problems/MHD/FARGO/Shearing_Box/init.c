#include "pluto.h"

double sb_Omega = 1.e-3;
double sb_q     = 1.5;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double Bmag, pr0, rnd, dvy;
  double Lx,pi2x, pi2y, pi2z;

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;    
  }
  Lx   = g_domEnd[IDIR] - g_domBeg[IDIR];
  pi2x = 2.0*CONST_PI/(g_domEnd[IDIR] - g_domBeg[IDIR]);
  pi2y = 2.0*CONST_PI/(g_domEnd[JDIR] - g_domBeg[JDIR]);
  pi2z = 2.0*CONST_PI/(g_domEnd[KDIR] - g_domBeg[KDIR]);
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 2.0*sb_A*x1;
  dvy  = sin(pi2x*x1 + 0.20) + sin(2.0*pi2x*x1 - 0.37);
  dvy *= sin(pi2y*x2 + 0.13) + sin(2.0*pi2y*x2 + 0.04);
  dvy *= sin(pi2z*x3 + 0.56) + sin(2.0*pi2z*x3 + 0.62);

  us[VX2] += g_inputParam[EPS]*Lx*sb_Omega*dvy/8.0;
  us[VX3]  = 0.0;
 
  #if EOS == ISOTHERMAL
   g_isoSoundSpeed = g_inputParam[CS];
  #endif
  pr0 = g_isoSoundSpeed*g_isoSoundSpeed*us[RHO];

  Bmag   = sqrt(2.0*pr0/g_inputParam[BETA]);
  us[BX1] = 0.0;
  us[BX2] = 0.0;
  us[BX3] = Bmag*sin(pi2x*x1);
  
  us[AX1] = 0.0;
  us[AX2] = -Bmag*cos(pi2x*x1)/pi2x;
  us[AX3] = 0.0;
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int    i,j,k;
  static int first_call=1;
  double *dx, *dy, *dz;
  double *x, *y, *z;
  double pm, Mxy, Rxy, alpha, aM, aR, vol;
  FILE *fp;

  dx = grid[IDIR].dx; dy = grid[JDIR].dx; dz = grid[KDIR].dx;
   x = grid[IDIR].x;   y = grid[JDIR].x;   z = grid[KDIR].x;

  pm = Mxy = Rxy = alpha = 0.0;
  DOM_LOOP(k,j,i){
    pm += 0.5*(  d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i]
               + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i]
               + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]);

    aM = d->Vc[BX1][k][j][i]*d->Vc[BX2][k][j][i];
    aR = d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*
         (d->Vc[VX2][k][j][i] - 2.0*sb_A*x[i]);

    Mxy   += aM;
    Rxy   += aR;
    alpha += (aR - aM)/(g_isoSoundSpeed*g_isoSoundSpeed*d->Vc[RHO][k][j][i]);
  }

  vol  = dx[IBEG]*dy[JBEG]*dz[KBEG];
  vol /= g_domEnd[IDIR] - g_domBeg[IDIR];
  vol /= g_domEnd[JDIR] - g_domBeg[JDIR];
  vol /= g_domEnd[KDIR] - g_domBeg[KDIR];

  pm    *= vol;
  Mxy   *= vol;
  Rxy   *= vol;
  alpha *= vol;

  #ifdef PARALLEL
   MPI_Allreduce (&pm, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pm = vol;

   MPI_Allreduce (&Mxy, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Mxy = vol;

   MPI_Allreduce (&Rxy, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Rxy = vol;

   MPI_Allreduce (&alpha, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   alpha = vol;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){
      fp = fopen("averages.dat","w");
      fprintf (fp,"#   t\t\t  dt\t      <B^2/2>\t    <Bx*By>\t<rho*ux*duy>\t<alpha>\n");
      first_call = 0;
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
      fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", 
               g_time, g_dt, pm, Mxy, Rxy, alpha);
    }
    fclose(fp);
  }
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int i,j,k;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1.e-2) d->Vc[RHO][k][j][i] = 1.e-2;
    }
    return;
  }
}

/* ************************************************************** */
double FARGO_SetVelocity(double x1, double x2)
/*
 *
 * PURPOSE
 *
 *   Compute the shear angular velocity to be subtracted from 
 *   the HD or MHD equations.
 * 
 **************************************************************** */
{
  return 2.0*sb_A*x1;
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
  #ifdef FARGO
   g[IDIR] = 0.0;
   g[JDIR] = sb_q*sb_Omega*v[VX1];
   g[KDIR] = 0.0;
  #else
   g[IDIR] = 2.0*sb_Omega*sb_q*sb_Omega*x1;
   g[JDIR] = 0.0;
   g[KDIR] = 0.0;
  #endif
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
void BodyForcePotential(double *v, double *phi, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  *phi = 0.0;
}
#endif
