#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* ------------------------------------------------------------
     Isentropic Vortex test problem.
     Reference: C.W. Shu, "High-order Finite Difference 
                and Finite Volume WENO Schemes and 
                Discontinuous Galerkin Methods for CFD"
                ICASE Report No.2001-11
                NASA/CR-2001-210865
   ------------------------------------------------------------ */
     
  double s, T, r2, xt, yt, eps;     
  double k0;

  g_gamma = 1.4;
  
  xt = x1 - 5.0;
  yt = x2 - 5.0;
   
  r2  = xt*xt + yt*yt;
  
  eps = 5.0;
  
  s   = 1.0;
  T   = 1.0 - (g_gamma - 1.0)*eps*eps/(8.0*g_gamma*CONST_PI*CONST_PI)*exp(1.0 - r2);
  
  k0     = eps/(2.0*CONST_PI)*exp(0.5*(1.0 - r2));
  us[RHO] = pow(T/s, 1.0/(g_gamma - 1.0));
  us[VX1] = 1.0 - k0*yt;
  us[VX2] = 1.0 + k0*xt;
  us[PRS] = T*us[RHO];
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
{ }

