#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, theta, phi, B0;

/* ------------------------------------------------------------
    Set initial condition for the 2 or 3-D blast wave
    problem. 
    The 2-D versions are taken from 

     -test #1:  Balsara & Spicer (1999), JCP 149, 270 
  
    The 3-D version is taken from 

     -test #2, #3 (in axisymmetry):  Ziegler (2004), JCP 196,393 
            
    and from 

     -test #4,5: Gardiner & Stone (2008), JCP, 227
                 Mignone & Tzeferacos (2010), JCP, 229
     -test #6:   2D version of Mignone&Tzeferacos (2010), JCP, 229
     -test #7: like #5 but with GLM 

   ------------------------------------------------------------ */

  g_gamma = g_inputParam[GAMMA];
  r = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r);

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = g_inputParam[P_OUT];
  
  theta = g_inputParam[THETA]*CONST_PI/180.0;
  phi   =   g_inputParam[PHI]*CONST_PI/180.0;
  B0    = g_inputParam[BMAG];
 
  us[BX1] = B0*sin(theta)*cos(phi);
  us[BX2] = B0*sin(theta)*sin(phi);
  us[BX3] = B0*cos(theta);
  
  if (r <= g_inputParam[RADIUS]) us[PRS] = g_inputParam[P_IN];
  
  #if GEOMETRY == CARTESIAN
   us[AX1] = 0.0;
   us[AX2] =  us[BX3]*x1;
   us[AX3] = -us[BX2]*x1 + us[BX1]*x2;
  #elif GEOMETRY == CYLINDRICAL
   us[AX1] = us[AX2] = 0.0;
   us[AX3] = 0.5*us[BX2]*x1;
  #endif


  #if BACKGROUND_FIELD == YES
   us[BX1] = us[BX2] = us[BX3] =
   us[AX1] = us[AX2] = us[AX3] = 0.0;
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
  int   i, j, k, nv;
  double  x1, x2, x3;

  if (side == X1_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
        }
        d->Vc[VX1][k][j][i] *= -1.0;
        d->Vc[BX1][k][j][i] *= -1.0;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][2*IBEG-i-1];
      #endif
    }else if (box->vpos == X3FACE){
        #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][2*IBEG-i-1];
      #endif
    }
  }

  if (side == X2_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
         
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][k][2*JBEG - j - 1][i];
        }
        EXPAND(d->Vc[BX1][k][j][i] *= -1.0; ,
               d->Vc[VX2][k][j][i] *= -1.0; ,
               d->Vc[BX3][k][j][i] *= -1.0;)
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX1s][k][j][i] = -d->Vs[BX1s][k][2*JBEG-j-1][i];
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = -d->Vs[BX3s][k][2*JBEG-j-1][i];
      #endif
    }
  }

  if (side == X3_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP (box,k,j,i){
  
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][2*KBEG - k - 1][j][i];
        }
        d->Vc[VX3][k][j][i] *= -1.0;
        d->Vc[BX3][k][j][i] *= -1.0;
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX1s][k][j][i] = d->Vs[BX1s][2*KBEG-k-1][j][i];
      #endif
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][2*KBEG-k-1][j][i];
      #endif
    }
  }
}
#if BACKGROUND_FIELD == YES
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
  double theta, phi;

  theta = g_inputParam[THETA]*CONST_PI/180.0;
  phi   =   g_inputParam[PHI]*CONST_PI/180.0;
 
  B0[IDIR] = g_inputParam[BMAG]*sin(theta)*cos(phi);
  B0[JDIR] = g_inputParam[BMAG]*sin(theta)*sin(phi);
  B0[KDIR] = g_inputParam[BMAG]*cos(theta);
}
#endif
