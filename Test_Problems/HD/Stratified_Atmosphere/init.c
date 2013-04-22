#include "pluto.h"

static real acf, bcf, ccf; /*  polynomial coefficient for g when r < 1 */

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* ------------------------------------------------------------------

   Set initial conditions for an hydrostatic
   atmosphere in cylindrical coordinates.
   Gravity is given as 

     g = -1/R^2                for R > 1
     g = a*R + b*R^2 + c*R^3   for R < 1

    where R = sqrt(r*r + z*z) is the spherical 
    radius. The coefficients a, b and c are chosen
    to guarantee continuity of g, its first and second
    derivative (optionally).
    Density and pressure are tied by the isothermal 
    condition

      p = rho/alpha

    so that the hydrostatic condition is

      1/rho drho/dr = alpha*g

    with the normalization rho = 1 at R = 1.

  ---------------------------------------------------------------- */

  double rs, scrh;

/* ------------------------------------------
    with this choice g, g' and g''
    will be continuous at R = 1
   ----------------------------------------- */

  acf = -10.0;
  bcf =  15.0;
  ccf = -6.0;

/* -----------------------------------
     with this choice g and g' 
     will be continuous
   ---------------------------------- */

  acf = -3.0;
  bcf =  2.0;
  ccf = 0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL 
   rs = sqrt(x1*x1 + x2*x2);
  #elif GEOMETRY == SPHERICAL
   rs = sqrt(x1*x1);
  #endif

  if (rs > 1.0){
    us[RHO] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
  }else{
    scrh   = 0.5*acf*(rs*rs - 1.0) + 1.0/3.0*bcf*(rs*rs*rs - 1.0)
              + 0.25*ccf*(rs*rs*rs*rs - 1.0);
    us[RHO] = exp(scrh*g_inputParam[ALPHA]);
  }

  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = us[RHO]/g_inputParam[ALPHA];
  us[TRC] = 0.0;
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
  int   i, j, k;
  double *x1, *x2, *x3;
  double rs;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  if (side == X1_END) {

    X1_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }

  } else if (side == X2_END) {

    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }

  } else if (side == X3_END) {   /* Only Cartesian */

    X3_END_LOOP(k,j,i){  
      rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
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
  double gs, rs;
  double acf, bcf, ccf;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;
  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #endif

  if (rs > 1.0) gs = -1.0/rs/rs;
  else          gs = rs*(acf + rs*(bcf + rs*ccf));

  #if GEOMETRY == CARTESIAN
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = gs*x3/rs;
  #elif GEOMETRY == CYLINDRICAL
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = 0.0;
  #endif

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
  double rs, phi;
  double acf, bcf, ccf, C;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #endif

  C = (0.5*acf + bcf/3.0 + ccf*0.25);  /* integration constant to make phi continuous */
  if (rs > 1.0) phi = -1.0/rs;
  else          phi = -rs*rs*(0.5*acf + rs*(bcf/3.0 + rs*ccf*0.25)) + C - 1.0;

  return phi;
}
#endif
