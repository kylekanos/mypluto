#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *   Brief description: 
 *
 *   Magnetized accretion torus in 2D (Cyl., Sph.) or 3D (cartesian)
 *   in a radially stratified atmosphere. 
 *
 *   The gravitational potential can be initialized either as 
 *   Newtonian or pseudo-Newtonian (Paczynsky - Wiita) potential. 
 *
 *   
 *   The user defined parameters of this problem, as they appear 
 *   in pluto.ini, allows for control in the Torus shape and the 
 *   contrast of physical values: 
 *
 *   R_min        Minimum cylindrical radius for the Torus (inner rim)
 *   R_max        Radius of the Torus where pressure is maximum
 *   den_max      Maximum density of the torus (outer rim)
 *   den_cut      Minimum density to define the last contour of the magnetic vec. pot. 
 *   BETA         Plasma beta 
 *   D_Con        Density contrast between atmosphere and Torus
 *   T_Con        Temperature contrast between atmosphere and Torus
 *   R_0          "Gravitational" radius in P-W potential (for R_0 = 0 -> Newton) 
 *   R_Sphere     Radius of the sink region, must be greater than R_0 
 *  
 *  Written by Petros Tzeferacos (petros.tzeferacos@ph.unito.it)
 *              
 *  Last modified: 08/05/2010
 *
 *  References : Global Magnetohydrodynamical Simulations of Accretion Tori, 
 *               Hawley, J.F., ApJ (2000) 528 462     
 *
 *               The dynamical stability of differentially rotating discs with
 *               constant specific angular momentum
 *               Papaloizou, J.C.B, Pringle, J.E, MNRAS (1984) 208 721
 *
 *               The Acceleration Mechanism of Resistive Magnetohydrodynamic Jets
 *               Launched from Accretion Disks. 
 *               Kuwabara, T., Shibata, K., Kudoh, T., Matsumoto, R. ApJ (2005) 621 921   
 *
 *********************************************************************** */
{
  double r, R, z, A_phi, g_gamma1;
  double a,b,c, p0, rho0, kappa, l_k, l_k2, r_maxpr;
  double Bo, epsilon, beta, Vo, R_in, rho_max, a_max, Tt_max;
  double prt, pra, sch, phi, Temp_contr, rhoa_one, integral_g_spline;
  double rhot, rhoa, ETA;
  static double R_sphere;
  static double ag,bg,cg;
  double rhocut;
  
  g_gamma1 = g_gamma - 1.0;
  #if GEOMETRY == CARTESIAN && DIMENSIONS == 3
   r   = sqrt(x1*x1 + x2*x2);
   R   = sqrt(x1*x1 + x2*x2 + x3*x3);
   z   = x3;
   phi = atan2(x2,x1);  
  #elif GEOMETRY == CYLINDRICAL && DIMENSIONS == 2
   r = fabs(x1);
   z = x2;
   R = sqrt(r*r + z*z);
   phi = 0.0;  /* 2D */ 
  #elif GEOMETRY == SPHERICAL && DIMENSIONS == 2
   r = x1*sin(x2);
   z = x1*cos(x2);
   R = x1;
   phi = 0.0; /* 2D */ 
  #endif
  R_sphere = g_inputParam[R_Sphere];
  r_maxpr = g_inputParam[R_max];
  R_in = g_inputParam[R_min];
  l_k = sqrt(r_maxpr)*r_maxpr/(r_maxpr-g_inputParam[R_0]);
  l_k2 = l_k*l_k;
  Vo = l_k;
  rhocut = g_inputParam[den_cut];  

  /* *************************************
   * Solve for kappa and c assuming pra=0
   *  
   * ************************************* */  
 
  c =  -1.0/(R_in-g_inputParam[R_0]) + 0.5*l_k2/(R_in*R_in);
  kappa = ((g_gamma-1.)/g_gamma)*(c + 0.5*l_k2*(-1./r_maxpr/r_maxpr) 
                            + 1./(r_maxpr-g_inputParam[R_0]))/(pow(g_inputParam[den_max],g_gamma-1.));

  /* ************************************
   * Torus density + pressure
   * ************************************* */  

  a          = c + 1.0/(R-g_inputParam[R_0]) - 0.5*l_k2/(r*r);
  a          = MAX(a, 0.0);
  rhot       = pow(0.4*a/kappa,1.5);
  prt        = kappa*pow(rhot,g_gamma);
  ETA        = g_inputParam[D_Con];
  Temp_contr = g_inputParam[T_Con];
  Tt_max     = kappa*pow(g_inputParam[den_max],g_gamma-1.);

  /* ************************************
   * atmospheric density + pressure 
   * (force equilibrium Fg and \nabla P)
   * ************************************ */

  rhoa = ETA*g_inputParam[den_max]*exp( ( 1./(R-g_inputParam[R_0]) - 1./(r_maxpr-g_inputParam[R_0]) )/Temp_contr/Tt_max ); 
  pra  = rhoa*Tt_max*Temp_contr;
  rhoa_one = ETA*g_inputParam[den_max]*exp( ( 1./(R_sphere-g_inputParam[R_0]) - 1./(r_maxpr-g_inputParam[R_0]) )/Temp_contr/Tt_max );

rhoa = ETA*g_inputParam[den_max]*exp( ( 1./(R-g_inputParam[R_0]) - 1./(R_in-g_inputParam[R_0]) )/Temp_contr/Tt_max); 
pra  = rhoa*Temp_contr*Tt_max;

/*
rhoa = ETA*g_inputParam[den_max];
pra = rhoa*Tt_max*Temp_contr;
*/


  Bo = sqrt(2.*kappa*pow(g_inputParam[den_max],g_gamma)/g_inputParam[BETA])/g_inputParam[den_max];

  if (prt > pra && r > 2.0) {
    us[RHO] = rhot;
    us[PRS] = prt;
   
    b = Vo/r;

    #if GEOMETRY == CARTESIAN  
     us[VX1] = -b*sin(phi);
     us[VX2] =  b*cos(phi);
     us[VX3] = 0.0;
    #elif GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL
     us[VX1] = 0.0;
     us[VX2] = 0.0;
     us[VX3] = b;
    #endif    
    us[TR] = 1.0;
  }else{
    us[RHO] = rhoa;
    us[PRS] = pra;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
    us[VX3] = 0.0;
    us[TRC] = 0.0;
  }
  
  #if PHYSICS == MHD || PHYSICS == RMHD

   us[BX1] = 0.0;
   us[BX2] = 0.0;
   us[BX3] = 0.0;


   us[AX1] = 0.0;
   us[AX2] = 0.0;
   us[AX3] = 0.0;
   if (rhot > rhocut && r > 2.0){
      A_phi = Bo*(rhot - rhocut);
   }else{
      A_phi = 0.0;
   }
    #if GEOMETRY == CARTESIAN  
     us[AX1] = -A_phi*sin(phi);
     us[AX2] =  A_phi*cos(phi);
     us[AX3] = 0.0;
 
    #elif GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL
     us[AX1] = 0.0;
     us[AX2] = 0.0;
     us[AX3] = A_phi;
  #endif
  #endif
  /* ***************************************** *
   *  smooth atm profile where the spline is, 
   *  integration constant at R=R_sphere
   * ***************************************** */
  integral_g_spline = -0.5*ag*R_sphere*R_sphere - bg*R_sphere*R_sphere*R_sphere/3. - 0.25*cg*R_sphere*R_sphere*R_sphere*R_sphere ;	  

  #if GEOMETRY == CYLINDRICAL 
  if(R<=R_sphere){
    us[RHO] = rhoa_one*exp(integral_g_spline/Temp_contr/Tt_max);
    us[PRS] = rhoa_one*exp(integral_g_spline/Temp_contr/Tt_max)*Tt_max*Temp_contr;
  }
  #endif
  g_smallPressure = 1.e-8;
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
  int     i, j, k, nv;
  double  *x1, *x2, *x3;
  static int first_call=1;
  static double vR1[NVAR];
  double R;
  static double R_bound; 

  if (first_call){
    R_bound = g_inputParam[R_Sphere];
    Init(vR1, R_bound, 0.0, 0.0);
    first_call = 0;
  }
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
  
  if (side == 0){  /* -- set a density threshold and sink-- */
   
    DOM_LOOP(k,j,i){ 
      #if GEOMETRY == CYLINDRICAL && DIMENSIONS == 2
       R = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == CARTESIAN && DIMENSIONS == 3
       R = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == SPHERICAL && DIMENSIONS == 2
       R = x1[i];
      #endif

      if (R <= R_bound){   /* -- only inside R_sphere -- */
       d->Vc[RHO][k][j][i] = vR1[RHO];
       d->Vc[PRS][k][j][i] = vR1[PRS];

       d->Vc[VX1][k][j][i] = 0.0; 
       d->Vc[VX2][k][j][i] = 0.0; 
       d->Vc[VX3][k][j][i] = 0.0; 
       #if PHYSICS == MHD || PHYSICS == RMHD
/*
       d->Vc[BX1][k][j][i] = 0.0; 
       d->Vc[BX2][k][j][i] = 0.0; 
       d->Vc[BX3][k][j][i] = 0.0; 
*/
       #if NTRACER > 0
        d->Vc[TRC][k][j][i] = 0.0; 
       #endif
      #endif        
     }

 /* --------------------------------------------------
     set a threshold value for density and pressure
    -------------------------------------------------- */

     if (d->Vc[RHO][k][j][i] < 5.e-2*vR1[RHO]) d->Vc[RHO][k][j][i] = 5.e-2*vR1[RHO];
     if (d->Vc[PRS][k][j][i] < 1.e-3*vR1[PRS]) d->Vc[PRS][k][j][i] = 1.e-3*vR1[PRS];
     g_smallDensity = 5.e-2*vR1[RHO];
     g_smallPressure = 1.e-3*vR1[PRS];
     
    }
  }
  
  if (side == X1_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
        }
        d->Vc[VX1][k][j][i] = MIN(d->Vc[VX1][k][j][i],0.0); 
        #if GEOMETRY == CYLINDRICAL
         d->Vc[BX1][k][j][i] = 0.0;
         d->Vc[BX2][k][j][i] = 0.0;
         d->Vc[BX3][k][j][i] = 0.0;
        #endif
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
      #endif
    }
  }
}

/* ************************************************************** */
void USERDEF_BOUNDARY (const Data *d, int side, Grid *grid) 
/* 
 *
 *    Definition of the sink region, where density and pressure 
 *    are kept to their initial value and velocity is set to zero.
 *
 *
 *
 *
 **************************************************************** */
{
  int     i, j, k, nv;
  double  *x1, *x2, *x3;
  static int first_call=1;
  static double vR1[NVAR];
  double R;
  static double R_bound; 

return;
  if (first_call){
    R_bound = g_inputParam[R_Sphere];
    Init(vR1, R_bound, 0.0, 0.0);
    first_call = 0;
  }
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
  
  if (side == 0){  /* -- set a density threshold and sink-- */
   
    DOM_LOOP(k,j,i){ 
      #if GEOMETRY == CYLINDRICAL && DIMENSIONS == 2
       R = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == CARTESIAN && DIMENSIONS == 3
       R = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == SPHERICAL && DIMENSIONS == 2
       R = x1[i];
      #endif

      if (R <= R_bound){   /* -- only inside R_sphere -- */
       d->Vc[RHO][k][j][i] = vR1[RHO];
       d->Vc[PRS][k][j][i] = vR1[PRS];

       d->Vc[VX1][k][j][i] = 0.0; 
       d->Vc[VX2][k][j][i] = 0.0; 
       d->Vc[VX3][k][j][i] = 0.0; 
       #if PHYSICS == MHD || PHYSICS == RMHD
/*
       d->Vc[BX1][k][j][i] = 0.0; 
       d->Vc[BX2][k][j][i] = 0.0; 
       d->Vc[BX3][k][j][i] = 0.0; 
*/
       #if NTRACER > 0
        d->Vc[TR][k][j][i] = 0.0; 
       #endif
      #endif        
     }

 /* --------------------------------------------------
     set a threshold value for density and pressure
    -------------------------------------------------- */

     if (d->Vc[RHO][k][j][i] < 5.e-2*vR1[RHO]) d->Vc[RHO][k][j][i] = 5.e-2*vR1[RHO];
     if (d->Vc[PRS][k][j][i] < 1.e-3*vR1[PRS]) d->Vc[PRS][k][j][i] = 1.e-3*vR1[PRS];
     g_smallDensity = 5.e-2*vR1[RHO];
     g_smallPressure = 1.e-3*vR1[PRS];
     
    }
  }
  
  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      for (nv = 0; nv < NVAR; nv++){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
      }
      d->Vc[VX1][k][j][i] = MIN(d->Vc[VX1][k][j][i],0.0); 
      #if GEOMETRY == CYLINDRICAL
       d->Vc[BX1][k][j][i] = 0.0;
       d->Vc[BX2][k][j][i] = 0.0;
       d->Vc[BX3][k][j][i] = 0.0;
      #endif
    }
    #ifdef STAGGERED_MHD
     KTOT_LOOP(k) for (j = -1; j < NX2_TOT; j++) IBEG_LOOP(i)
       d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
    #endif
  }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  double r, R, R1, z, gR;
  static double R_sphere;
  static double ag,bg,cg;
  static int first_call=1;
    
  #if GEOMETRY == CARTESIAN && DIMENSIONS == 3
   r   = sqrt(x1*x1 + x2*x2);
   R   = sqrt(r*r + x3*x3);
   z   = x3;
  #elif GEOMETRY == CYLINDRICAL && DIMENSIONS == 2
   r = fabs(x1);
   z = x2;
   R = sqrt(r*r + z*z);
  #elif GEOMETRY == SPHERICAL && DIMENSIONS == 2
   r = x1*sin(x2);
   z = x1*cos(x2);
   R = x1;
  #endif
  R_sphere = g_inputParam[R_Sphere];
  R1 = R - g_inputParam[R_0];
  gR = -1.0/(R1*R1);
  /* compute spline coef for smoothing gravity from R_sphere to boundary*/
  if(first_call){
  
  cg = - 2./(pow(R_sphere,2.)*pow(R_sphere-g_inputParam[R_0],3.)) - 1./(pow(R_sphere,3.)*pow(R_sphere-g_inputParam[R_0],2.))
       - 3./(R_sphere*pow(R_sphere-g_inputParam[R_0],4.));
  bg = -3.*(1./pow(R_sphere - g_inputParam[R_0],4.) + cg*R_sphere);
  ag = 2./pow(R_sphere-g_inputParam[R_0],3.) - 2.*bg*R_sphere - 3.*cg*R_sphere*R_sphere;
  first_call = 0;
  }
   if (R <= R_sphere) gR=ag*R + bg*R*R+ cg*R*R*R; /* gR=aR+bR^2+cR^3, coef a,b,c are found by contin of gR, gR', gR'' at R=1 */
 
  #if GEOMETRY == CARTESIAN
   g[IDIR] = gR*x1/R;
   g[JDIR] = gR*x2/R;
   g[KDIR] = gR*x3/R;
   
  #elif GEOMETRY == CYLINDRICAL
     
   g[IDIR] = gR*x1/R;
   g[JDIR] = gR*x2/R;
   g[KDIR] = 0.0; 
   
  #elif GEOMETRY == SPHERICAL
     
   g[IDIR] = gR;
   g[JDIR] = 0.0;
   g[KDIR] = 0.0; 
   
  #endif
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
