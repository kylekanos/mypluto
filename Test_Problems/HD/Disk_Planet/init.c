/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet problem.

  Disk-Planet setup.
 
  Reference paper:
   "A Conservative orbital advection scheme for simulations
    of magnetized shear flows with the PLUTO Code"
    Mignone et al, A&A (2012)
 
  -------------------------------------------------------------
   Independent of the coordinate system (polar/spherical), we
   adopt the following conventions:
 
    r = spherical radius
    R = cylindrical radius
    z = cylindrical height
   th = meridional angle

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

#define MIN_DENSITY 1e-8

static void NormalizeDensity (const Data *d, Grid *g);
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs;
  double scrh;

  #if EOS == IDEAL
   g_gamma = 1.01;
  #endif
  g_unitLength   = 5.2*CONST_au;
  g_unitDensity  = CONST_Msun/g_unitLength/g_unitLength/g_unitLength;
  g_unitVelocity = sqrt(CONST_G*g_inputParam[Mstar]*CONST_Msun/g_unitLength)/(2.*CONST_PI);

  #if ROTATING_FRAME == YES
   g_OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   g_OmegaZ *= 2.0*CONST_PI;
  #endif
  
  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   th = x2;
   R  = r*sin(th);
   z  = r*cos(th);
  #endif
  
  H      = 0.05*R;
  OmegaK = 2.0*CONST_PI/(R*sqrt(R));
  cs     = H*OmegaK;
  
  scrh   = (0.5*CONST_PI - th)*r/H;
  us[RHO] = 1.0/(R*sqrt(R))*exp(-0.5*scrh*scrh);
  
  us[VX1] = us[VX2] = us[VX3] = 0.0;

  us[iVPHI] = R*(OmegaK - g_OmegaZ);
  #if EOS == IDEAL
   us[PRS] = us[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
//   g_isoSoundSpeed = cs;
   g_isoSoundSpeed = CONST_PI*0.1;
  #endif

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
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, OmegaK, v[256];
  static int do_once = 1;
  
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  #if DIMENSIONS == 3
  if (side == 0){
    if (do_once){
      NormalizeDensity(d, grid);
      do_once = 0;
    }
  }
  #endif

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      for (nv = 0; nv < NVAR; nv++){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      }
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
       R = x1[i];
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      for (nv = 0; nv < NVAR; nv++){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      }
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
    }
  }
}

/* ************************************************************** */
void NormalizeDensity (const Data *d, Grid *grid)
/*
 *
 * Normalize density and pressure as   rho -> K*rho, where
 *
 *   K = M/(\sum rho*dV)
 *
 **************************************************************** */
{
  int   i, j, k;
  double *dV1, *dV2, *dV3;
  double dV, mass, gmass, mc;
        
  dV1 = grid[IDIR].dV; 
  dV2 = grid[JDIR].dV; 
  dV3 = grid[KDIR].dV;

  mass = 0.0;
  DOM_LOOP(k,j,i){
    dV    = dV1[i]*dV2[j]*dV3[k];
    mass += dV*d->Vc[RHO][k][j][i];
  }
                        
#ifdef PARALLEL
  gmass = 0.;
  MPI_Allreduce (&mass, &gmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  mass = gmass;
#endif
        
  mc  = 0.5*g_inputParam[Mdisk]*CONST_Msun;
  mc /= g_unitDensity*g_unitLength*g_unitLength*g_unitLength*mass;
  DOM_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= mc;
    #if EOS == IDEAL
     d->Vc[PRS][k][j][i] *= mc;
    #endif
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
  g[IDIR] = 0.0;
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
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi;

  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
   x  = R*cos(x2);
   y  = R*sin(x2);
  #elif (GEOMETRY == SPHERICAL)
   r  = x1;
   th = x2;
   R = r*sin(th);
   z = r*cos(th);
   x = r*sin(th)*cos(x3);
   y = r*sin(th)*sin(x3);
  #endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */

  #if ROTATING_FRAME == NO
   double OmegaZ;
   t = glob_time;
   if (ISTEP == 2) t += delta_t;
   OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   OmegaZ *= 2.0*CONST_PI;

   xp = cos(OmegaZ*t);
   yp = sin(OmegaZ*t);
  #else
   xp = 1.0/sqrt(2.0);  /* initial planet position */
   yp = 1.0/sqrt(2.0); 
  #endif

  d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z);
  rsm = 0.03*R;
  if (d > rsm) phiplanet = g_inputParam[Mplanet]/d;
  else phiplanet = g_inputParam[Mplanet]/d*(pow(d/rsm,4.)-2.*pow(d/rsm,3.)+2.*d/rsm);
  
  phi  = - 4.0*CONST_PI*CONST_PI/g_inputParam[Mstar];
  phi *= (g_inputParam[Mstar]/r + phiplanet*CONST_Mearth/CONST_Msun);

  return phi;
}
#endif
