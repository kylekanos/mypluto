/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Stellar wind test problem.

  Sets initial condition for a spherically symmetric radial wind blowing 
  from the origin of coordinates.
  Dimensions are chosen so that the spherical wind shell has radius 1,
  density 1 and velocity 1.
  Inside the shell, flow values are kept constant in time using the
  UserDefBoundary() function when \c side is equal to 0. 
  The input parameters that control the problem dynamics are
  
  \param CS_WIND:    sets the sound speed in the wind region;
  \param RHO_AMB:    sets the ambient density;
  \param CS_AMB:     sets the ambient sound speed;
  \param V_CSM:      sets the velocity of the star with respect 
                     to the background;

  \authors A. Mignone (mignone@ph.unito.it)\n
  
  \date   Sept 25, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*! 
 * Sets initial condition for the ambient medium.
 *
 *********************************************************************** */
{
  double rho, cs, R;

  rho = g_inputParam[RHO_AMB];
  cs  = g_inputParam[CS_AMB];
  us[RHO] = rho;
  us[PRS] = cs*cs*rho/g_gamma;

  #if GEOMETRY == CARTESIAN
   us[VX1] = 0.0;         
   us[VX2] = 0.0;         
   us[VX3] = -g_inputParam[V_CSM];       
  #elif GEOMETRY == CYLINDRICAL
   us[VX1] = 0.0;         
   us[VX2] = -g_inputParam[V_CSM];          
  #elif GEOMETRY == SPHERICAL
   us[VX1] =  g_inputParam[V_CSM]*cos(x2);
   us[VX2] = -g_inputParam[V_CSM]*sin(x2);
  #endif

  R = sqrt(x1*x1 + x2*x2);
  us[TRC] = (R <= 1.0 ? 1.0:0.0);
  g_smallPressure = 1.e-5;
}

/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/*
 *
 ****************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox * box, int side, Grid *grid) 
/*!
 * Sets inflow boundary condition at the top boundary (side == X2_END)
 * and the stellar wind region when side == 0.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double  r, r0, cs;
  double  Vwind  = 1.0, rho, vr;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  if (side == 0){
    r0 = 1.0;
    cs = g_inputParam[CS_WIND];
    TOT_LOOP(k,j,i){ 
      #if GEOMETRY == CARTESIAN
       r  = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
       if (r <= r0){  
         vr    = tanh(r/r0/0.1)*Vwind;
         rho   = Vwind*r0*r0/(vr*r*r);
         d->Vc[RHO][k][j][i] = rho;
         d->Vc[VX1][k][j][i] = Vwind*x1[i]/r;
         d->Vc[VX2][k][j][i] = Vwind*x2[j]/r;
         d->Vc[VX3][k][j][i] = Vwind*x3[k]/r;
         d->Vc[PRS][k][j][i] = cs*cs/g_gamma*pow(rho,g_gamma);
         d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
       }    
      #elif GEOMETRY == CYLINDRICAL
       r  = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
       if (r <= r0){  
         vr    = tanh(r/r0/0.1)*Vwind;
         rho   = Vwind*r0*r0/(vr*r*r);
         d->Vc[RHO][k][j][i] = rho;
         d->Vc[VX1][k][j][i] = Vwind*x1[i]/r;
         d->Vc[VX2][k][j][i] = Vwind*x2[j]/r;
         d->Vc[PRS][k][j][i] = cs*cs/g_gamma*pow(rho,g_gamma);
         d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
       }    
      #endif
    }
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
   X2_BEG_LOOP(k,j,i){ }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    
    cs  = g_inputParam[CS_AMB];
    rho = g_inputParam[RHO_AMB];
    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CYLINDRICAL
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = -g_inputParam[V_CSM];  
       d->Vc[RHO][k][j][i] =  rho;       
       d->Vc[PRS][k][j][i] =  cs*cs*rho/g_gamma; 
      #endif
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    X3_BEG_LOOP(k,j,i){}
  }
 
  if (side == X3_END){  /* -- X3_END boundary -- */
    
    cs  = g_inputParam[CS_AMB];
    rho = g_inputParam[RHO_AMB];
    X3_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = 0.0;
       d->Vc[VX3][k][j][i] = -g_inputParam[V_CSM];  
       d->Vc[RHO][k][j][i] =  rho;       
       d->Vc[PRS][k][j][i] =  cs*cs*rho/g_gamma; 
      #endif
    }
  }
}
