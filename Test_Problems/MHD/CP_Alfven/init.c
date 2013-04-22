#include "pluto.h"

static double AX_VEC (double x, double y, double z);
static double AY_VEC (double x, double y, double z);
static double AZ_VEC (double x, double y, double z);

static double ta = 0.0, tb = 0.0;
static double ca, sa, cg, sg;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 * Circularly polarized Alfven waves.
 * 
 * Reference paper:
 *
 *  "High-order conservative finite difference GLM-MHD 
 *   schemes for cell-centered MHD"
 * Mignone, Tzeferacos & Bodo JCP (2010)
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double tg, eps;
  double x, y, z, kx, ky,kz;
  double phi, cb, sb;
  double vx, vy, vz, Bx, By, Bz, rho, pr;
  double Lx, Ly, Lz;
  
  Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR];

  if (first_call == 1){
    D_EXPAND(      ,
      ta = Lx/Ly;  ,
      tb = Lx/Lz;)

    if (prank == 0) printf ("ta = %f; tb = %f\n",ta, tb);
    first_call = 0;
  }
  x = x1; y = x2; z = x3;

/* ------------------------------------------------------------ 
    This is a rotated version of a 1D setup where 
 
     rho = 1
     vx  = V0
     vy  = |eps|*cos(phi)
     vz  = |eps|*sin(phi)
     Bx  = 1 
     By  = eps*cos(phi)
     Bz  = eps*cos(phi)
     pr  = P0
 
    where:
   
     V0:  translation velocity in the 1D x direction (=g_inputParam[VEL0])
     psi: is the phase (kx*x in 1D and (kx*x + ky*y) in 2D);
     eps: is the wave amplitude (=g_inputParam[EPS]). 
          eps > 0 implies right going waves;
          eps < 0 implied left  going waves.

    With this normalization, the Alfven velocity VA = Bx/sqrt(rho)
    is always one.

    Rotations are specified by 
    
     ta = ky/kx  *and*  tb = kz/kx

    expressing the ratios between the y- and z- components of 
    the wave vector with the x component.
    We choose the wavelength in x-direction to be always 1, so
    that kx = 2*pi/lx = 2*pi, ky = 2*pi/ly, kz = 2*pi/lz so that
    tz = 1/ly, ty = 1/lz.
    By choosing the domain in the x-direction to be 1, the 
    extent in y and z should be adjusted so that Ly/ly = 1,2,3,4...
    (and similarly for Lz/lz) becomes an integer number to 
    ensure periodicity. If one wants 1 wave length in both 
    directions than 

     Lx = 1,  Ly = 1/ta, Lz = 1/tb
 
    If N wavelengths are specified in the y-direction than
  
     Lx = 1, Ly = N/ta,  Lz = N/tb

    1D is recovered by specifying tz = ty = 0.

    The final time step is one period and is found from 

                                
     VA = omega/|k| -->  1/T = sqrt(1 + ta^2 + tb^2)

   -------------------------------------------------------------- */     
  
  eps = g_inputParam[EPS];   /* wave amplitude */
  kx  = 2.0*CONST_PI/1.0;
  ky  = ta*kx;
  kz  = tb*kx;

  ca  = 1.0/sqrt(ta*ta + 1.0);
  cb  = 1.0/sqrt(tb*tb + 1.0);
  sa  = ta*ca;
  sb  = tb*cb;

  phi = D_EXPAND(kx*x, + ky*y, + kz*z);

/* -- define the 1-D solution -- */

  rho = 1.0;
  pr  = g_inputParam[PR0];
  vx  = g_inputParam[VEL0];  /* translational velocity */
  vy  = fabs(eps)*sin(phi);
  vz  = fabs(eps)*cos(phi);
  Bx  = 1.0;
  By  = eps*sqrt(rho)*sin(phi);
  Bz  = eps*sqrt(rho)*cos(phi);

/* -- rotate solution -- */

  us[RHO] = rho;
  us[PRS] = pr;

  tg = tb*ca;
  cg = 1.0/sqrt(tg*tg + 1.0); sg = cg*tg;

  us[VX1] =  vx*cg*ca - vy*sa - vz*sg*ca;
  us[VX2] =  vx*cg*sa + vy*ca - vz*sg*sa;
  us[VX3] =  vx*sg            + vz*cg;
  us[TR] = 0.0;

  us[BX1] =  Bx*cg*ca - By*sa - Bz*sg*ca;
  us[BX2] =  Bx*cg*sa + By*ca - Bz*sg*sa;
  us[BX3] =  Bx*sg            + Bz*cg;

  #ifdef STAGGERED_MHD
   us[AX1] = AX_VEC(x,y,z);
   us[AX2] = AY_VEC(x,y,z);
   us[AX3] = AZ_VEC(x,y,z);
  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int   i,j,k;
  double bmag, err[3], gerr[3], Vtot, dV;
  double ***Bx, ***By, ***Bz;
  double *dx, *dy, *dz;
  static double ***Bx0, ***By0, ***Bz0;
  FILE  *fp;
  
  Bx = d->Vc[BX1]; By = d->Vc[BX2]; Bz = d->Vc[BX3];
  dx = grid[IDIR].dx; dy = grid[JDIR].dx; dz = grid[KDIR].dx;

  if (Bx0 == NULL){
    Bx0 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    By0 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Bz0 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    print1 ("ANALYSIS: Initializing bmag\n");
    DOM_LOOP(k,j,i){
      Bx0[k][j][i] = Bx[k][j][i];
      By0[k][j][i] = By[k][j][i];
      Bz0[k][j][i] = Bz[k][j][i];
    }
    return;  
  }

  for (i = 0; i < 3; i++) err[i] = 0.0;
  DOM_LOOP(k,j,i){
    dV = D_EXPAND(dx[i], *dy[j], *dz[k]);
    err[0] += fabs(Bx[k][j][i] - Bx0[k][j][i])*dV;
    err[1] += fabs(By[k][j][i] - By0[k][j][i])*dV;
    err[2] += fabs(Bz[k][j][i] - Bz0[k][j][i])*dV;
  }

  Vtot = D_EXPAND(
                  (g_domEnd[IDIR] - g_domBeg[IDIR]),
                 *(g_domEnd[JDIR] - g_domBeg[JDIR]),
                 *(g_domEnd[KDIR] - g_domBeg[KDIR]));

  for (i = 0; i < 3; i++) err[i] /= Vtot;

  #ifdef PARALLEL 
   MPI_Allreduce (err, gerr, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   for (i = 0; i < 3; i++) err[i] = gerr[i];
   MPI_Barrier (MPI_COMM_WORLD);
  #endif
  
  if (prank == 0){
    double errtot;
    errtot = sqrt(err[0]*err[0] + err[1]*err[1] + err[2]*err[2]);
    fp  = fopen("bmag.dat","w");
    printf ("analysis: writing to disk...\n");
    fprintf (fp,"%d  %10.3e\n", (int)(1.0/dx[IBEG]),  errtot);
    fclose (fp);
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
{ }


double AX_VEC (double x, double y, double z)
{
  double xn, yn, zn, eps;
  double kx, ky, kz, phi;
  double kxn, kyn, kzn;
  double Axn, Ayn, Azn;

  eps = g_inputParam[EPS];   
  kx  = 2.0*CONST_PI/1.0;
  ky  = ta*kx;
  kz  = tb*kx;

  phi = D_EXPAND(kx*x, + ky*y, + kz*z);

  xn =  ca*cg*x + sa*cg*y + sg*z;
  yn =    -sa*x +    ca*y;
  zn = -ca*sg*x - sa*sg*y + cg*z;

  kxn =  ca*cg*kx + sa*cg*ky + sg*kz;
  kyn =    -sa*kx +    ca*ky;
  kzn = -ca*sg*kx - sa*sg*ky + cg*kz;

  Axn = 0.0;
  Ayn = eps*sin(phi)/kxn;
  Azn = eps*cos(phi)/kxn + yn;

  return(Axn*cg*ca - Ayn*sa - Azn*sg*ca);
}

double AY_VEC (double x, double y, double z)
{
  double xn, yn, zn, eps;
  double kx, ky, kz, phi;
  double kxn, kyn, kzn;
  double Axn, Ayn, Azn;

  eps = g_inputParam[EPS];   
  kx  = 2.0*CONST_PI/1.0;
  ky  = ta*kx;
  kz  = tb*kx;

  phi = D_EXPAND(kx*x, + ky*y, + kz*z);

  xn =  ca*cg*x + sa*cg*y + sg*z;
  yn =    -sa*x +    ca*y;
  zn = -ca*sg*x - sa*sg*y + cg*z;

  kxn =  ca*cg*kx + sa*cg*ky + sg*kz;
  kyn =    -sa*kx +    ca*ky;
  kzn = -ca*sg*kx - sa*sg*ky + cg*kz;

  Axn = 0.0;
  Ayn = eps*sin(phi)/kxn;
  Azn = eps*cos(phi)/kxn + yn;

  return(Axn*cg*sa + Ayn*ca - Azn*sg*sa);
}
double AZ_VEC (double x, double y, double z)
{
  double xn, yn, zn, eps;
  double kx, ky, kz, phi;
  double kxn, kyn, kzn;
  double Axn, Ayn, Azn;

  eps = g_inputParam[EPS];   
  kx  = 2.0*CONST_PI/1.0;
  ky  = ta*kx;
  kz  = tb*kx;

  phi = D_EXPAND(kx*x, + ky*y, + kz*z);

  xn =  ca*cg*x + sa*cg*y + sg*z;
  yn =    -sa*x +    ca*y;
  zn = -ca*sg*x - sa*sg*y + cg*z;

  kxn =  ca*cg*kx + sa*cg*ky + sg*kz;
  kyn =    -sa*kx +    ca*ky;
  kzn = -ca*sg*kx - sa*sg*ky + cg*kz;

  Axn = 0.0;
  Ayn = eps*sin(phi)/kxn;
  Azn = eps*cos(phi)/kxn + yn;

  return(Axn*sg             + Azn*cg);
}
