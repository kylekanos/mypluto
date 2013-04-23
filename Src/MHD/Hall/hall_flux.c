/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the resistive MHD flux.

  Compute the resistive fluxes for the induction and energy equations.
  In the induction equation, fluxes are computed by explicitly writing
  the curl operator in components. In Cartesian components, for instance,
  one has
  \f[
     \pd{\vec{B}}{t} = - \nabla\times{\vec{E}_{\rm res}} =
      \pd{}{x}\left(\begin{array}{c}
       0           \\ \noalign{\medskip}
       \eta_z J_z  \\ \noalign{\medskip}
     - \eta_y J_y \end{array}\right)
       +
      \pd{}{y}\left(\begin{array}{c}
     - \eta_z J_z  \\ \noalign{\medskip}
           0       \\ \noalign{\medskip}
       \eta_x J_x \end{array}\right)
        +
      \pd{}{z}\left(\begin{array}{c}
        \eta_y J_y \\ \noalign{\medskip}
       -\eta_x J_x \\ \noalign{\medskip}
            0    \end{array}\right)      
       \,,\qquad
     \left(\vec{E}_{\rm res} = \tens{\eta}\cdot\vec{J}\right)
  \f]
  where \f$\tens{\eta}\f$ is the resistive diagonal tensor and 
  \f$J = \nabla\times\vec{B}\f$ is the current density.
  The corresponding contribution to the energy equation is
  \f[
    \pd{E}{t} = -\nabla\cdot\Big[(\tens{\eta}\cdot\vec{J})\times\vec{B}\Big]
              = - \pd{}{x}\left(\eta_yJ_yB_z - \eta_zJ_zB_y\right)
                - \pd{}{y}\left(\eta_zJ_zB_x - \eta_xJ_xB_z\right)
                - \pd{}{z}\left(\eta_xJ_xB_y - \eta_yJ_yB_x\right)
  \f] 
  
  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics" \n
       Mignone et al, ApJS (2012) 198, 7M
       
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date   Nov 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void HallFlux (Data_Arr V, double **RF, double **dcoeff,
                    int beg, int end, Grid *grid)
/*! 
 *  Compute Hall fluxes for the induction and energy
 *  equations. Also, compute the diffusion coefficient.
 *  Called by either ParabolicFlux or ParabolicRHS.
 *
 * \param [in]   V      3D data array of primitive variables
 * \param [out]  RF     the resistive fluxes
 * \param [out]  dcoeff the diffusion coefficients evaluated at 
 *                      cell interfaces
 * \param [in]   beg    initial index of computation
 * \param [in]   end    final   index of computation
 * \param [in]  grid    pointer to an array of Grid structures.
 *
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double x1, x2, x3, scrh;
  double lHall, vi[NVAR];
  static double **J;

  if (J == NULL) J = ARRAY_2D(NMAX_POINT, 3, double);
  
  GetFullCurrent (V, J, grid);
 
  D_EXPAND(x1 = grid[IDIR].x[*g_i];  ,
           x2 = grid[JDIR].x[*g_j];  ,
           x3 = grid[KDIR].x[*g_k]; )

/* ----------------------------------------------- 
     Compute resistive flux for the induction
     and energy equations
   ----------------------------------------------- */
  
  if (g_dir == IDIR){
  
    j = *g_j; k = *g_k;    
    for (i = beg; i <= end; i++){

    /* -- interface value -- */

      x1 = grid[IDIR].xr[i];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);
      }
      lHall_Func (vi, x1, x2, x3, &lHall);

      EXPAND(                                       ,
             RF[i][BX2] = -lHall * (J[i][IDIR] * vi[BX2] - J[i][JDIR] * vi[BX1]);   ,
             RF[i][BX3] = -lHall * (J[i][IDIR] * vi[BX3] - J[i][KDIR] * vi[BX1]);)
/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[i] = scrh;
*/
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = |lHall|B0/dl^2
      (Whistler wave timescale)
     ---------------------------------------------- */
     
             dcoeff[i][MX2] = fabs(lHall)*vi[BX1];
    }

  }else if (g_dir == JDIR){
  
    i = *g_i; k = *g_k;    
    for (j = beg; j <= end; j++){

    /* -- interface value -- */

      x2 = grid[JDIR].xr[j];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);
      }
      lHall_Func (vi, x1, x2, x3, &lHall);

      EXPAND(RF[j][BX1] = -lHall * (J[j][JDIR] * vi[BX1] - J[j][IDIR] * vi[BX2]);   ,
                                                    ,
             RF[j][BX3] = -lHall * (J[j][JDIR] * vi[BX3] - J[j][KDIR] * vi[BX2]);)
/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[j] = scrh;
*/
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = |lHall|B0/dl^2
      (Whistler wave timescale)
     ---------------------------------------------- */

      dcoeff[j][MX2] = fabs(lHall)*vi[BX2];
    }

  }else if (g_dir == KDIR){
  
    i = *g_i; j = *g_j;    
    for (k = beg; k <= end; k++){

    /* -- interface value -- */

      x3 = grid[KDIR].xr[k];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k + 1][j][i]);
      }
      lHall_Func (vi, x1, x2, x3, &lHall);

      RF[k][BX1] = -lHall * (J[k][KDIR] * vi[BX1] - J[k][IDIR] * vi[BX3]);
      RF[k][BX2] = -lHall * (J[k][KDIR] * vi[BX2] - J[k][JDIR] * vi[BX3]);

/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[k] = scrh;
*/ 
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = |lHall|B0/dl^2
      (Whistler wave timescale)
     ---------------------------------------------- */

      dcoeff[k][MX2] = fabs(lHall)*vi[BX3];
    }
  }
}
