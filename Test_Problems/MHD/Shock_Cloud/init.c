#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* -------------------------------------------------------
    Set initial condition for the 2D and 3D 
    shock-cloud problem.
    Reference paper:

     "The divB = 0 Constraint in SHock-Capturing MHD codes"
      G. Toth, JCP 161, 605-652 (2000)

     "A central constrained transport scheme for ideal
      magnetohydrodynamics" 
      U. Ziegler, JCP 196 (2004), 393
   -------------------------------------------------------- */

 double r, x0, y0, z0;

  x0 = 0.8;
  y0 = 0.5;
  z0 = 0.5;

  if (x1 < 0.6) {
    us[RHO] = 3.86859;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
    us[VX3] = 0.0;
    us[BX1] = 0.0;
    us[BX2] =  g_inputParam[B_POST];
    us[BX3] = -g_inputParam[B_POST];
    us[PRS] = 167.345;
  }else{
    us[RHO] = 1.0;
    us[VX1] = -11.2536;
    us[VX2] = 0.0;
    us[VX3] = 0.0;
    us[BX1] = 0.0;
    us[BX2] = g_inputParam[B_PRE];
    us[BX3] = g_inputParam[B_PRE];
    us[PRS] = 1.0;
  } 

  /*  ----  CLOUD  ----  */

  r = D_EXPAND(  (x1 - x0)*(x1 - x0) ,
               + (x2 - y0)*(x2 - y0) ,
               + (x3 - z0)*(x3 - z0));
                
  if (sqrt(r) < g_inputParam[RADIUS]) {
    us[RHO] = 10.0;
  }

/* no need for potential vector, 
   since CT_VEC_POT_INIT is set to NO

  us[AX1] = x3*us[BX2] - x2*us[BX3];
  us[AX2] = 0.0;
  us[AX3] = 0.0;
*/
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{ }


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
  int   i, j, k;

  if (side == X1_END){          /* -- select the boundary side -- */
    if (box->vpos == CENTER){   /* -- select the variable position -- */
      BOX_LOOP(box,k,j,i){      /* -- Loop over boundary zones -- */
        d->Vc[RHO][k][j][i] = 1.0;
        EXPAND(d->Vc[VX1][k][j][i] = -11.2536;  ,
               d->Vc[VX2][k][j][i] = 0.0;       ,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = 1.0;
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;        ,
               d->Vc[BX2][k][j][i] = g_inputParam[B_PRE]; ,
               d->Vc[BX3][k][j][i] = g_inputParam[B_PRE];)
      }
    }else if (box->vpos == X2FACE){  /* -- y staggered field -- */
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = g_inputParam[B_PRE];
      #endif
    }else if (box->vpos == X3FACE){  /* -- z staggered field -- */
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = g_inputParam[B_PRE];
      #endif
    }
  }
}


/* ************************************************************** */
void USERDEF_BOUNDARY (const Data *d, int side, Grid *grid) 
/* 
 *
 *
 **************************************************************** */
{
  int   i, j, k;

return;
  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] = 1.0;
      d->Vc[VX1][k][j][i] = -11.2536;
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[VX3][k][j][i] = 0.0;
      d->Vc[PRS][k][j][i] = 1.0;
      d->Vc[BX1][k][j][i] = 0.0;
      d->Vc[BX2][k][j][i] = g_inputParam[B_PRE];
      d->Vc[BX3][k][j][i] = g_inputParam[B_PRE];
 /*
      #ifdef STAGGERED_MHD
       D_EXPAND(                           ,
         d->Vs[BX2s][k][j][i] = g_inputParam[B_PRE]; ,
         d->Vs[BX3s][k][j][i] = g_inputParam[B_PRE];)
      #endif
 */
    }
    #ifdef STAGGERED_MHD
     D_EXPAND(                                                ;  ,
       KTOT_LOOP(k) for (j = -1; j < NX2_TOT; j++) IEND_LOOP(i)
         d->Vs[BX2s][k][j][i] = g_inputParam[B_PRE];                       ,

       for (k = -1; k < NX3_TOT; k++) JTOT_LOOP(j) IEND_LOOP(i)
         d->Vs[BX3s][k][j][i] = g_inputParam[B_PRE];)
    #endif
  }

}

