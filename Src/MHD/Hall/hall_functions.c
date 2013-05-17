/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the **FULL** electric current as needed by Hall
  
  \todo compute current in curvilinear coordinates

  \authors G. Lesur
           
  \date   April 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double hall_xh;
double hall_invdmin;

void ComputeCurrent(Data_Arr V, Grid *grid) {
    int  i, j, k;
    double dx_1, dy_1, dz_1;
    double *inv_dx,  *inv_dy,  *inv_dz;
    double ***Bx, ***By, ***Bz;
    double ***Jx, ***Jy, ***Jz;
    double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
    
    D_EXPAND(inv_dx  = grid[IDIR].inv_dx; inv_dxi = grid[IDIR].inv_dxi;  ,
             inv_dy  = grid[JDIR].inv_dx; inv_dyi = grid[JDIR].inv_dxi;  ,
             inv_dz  = grid[KDIR].inv_dx; inv_dzi = grid[KDIR].inv_dxi;)
    
    EXPAND(Bx = V[BX1]; ,
           By = V[BX2]; ,
           Bz = V[BX3]; )
    
    Jx = V[JX1];
    Jy = V[JX2];
    Jz = V[JX3];
    
    DOM_LOOP(k,j,i) {
        dx_1 = inv_dx[i];
        dy_1 = inv_dy[j];
        dz_1 = inv_dz[k];
        
        dxBy = 0.5*(By[k][j][i+1]-By[k][j][i-1])*dx_1;
        dxBz = 0.5*(Bz[k][j][i+1]-Bz[k][j][i-1])*dx_1;
        dyBz = 0.5*(Bz[k][j+1][i]-Bz[k][j-1][i])*dy_1;
        dyBx = 0.5*(Bx[k][j+1][i]-Bx[k][j-1][i])*dy_1;
        dzBx = 0.5*(Bx[k+1][j][i]-Bx[k-1][j][i])*dz_1;
        dzBy = 0.5*(By[k+1][j][i]-By[k-1][j][i])*dz_1;
        
        Jx[k][j][i] = dyBz - dzBy;
        Jy[k][j][i] = dzBx - dxBz;
        Jz[k][j][i] = dxBy - dyBx;
    }
    
}

/* ********************************************************************* */
void GetFullCurrent (Data_Arr V, double **curlB, Grid *grid)
/*!
 *
 * \param [in] V       A 3D array of cell-centered primitive quantities
 * \param [out] curlB  A 1D array containing two components of the current
 *                     in a the direction set by ::g_dir
 * \param [in] grid    a pointer to an array of Grid structures
 *   
 *********************************************************************** */
{
  int  i, j, k;
  double dx_1, dy_1, dz_1, r_1, s_1;
  double *inv_dx,  *inv_dy,  *inv_dz;
  double *inv_dxi, *inv_dyi, *inv_dzi;
  double *r, *th;
  double ***Bx, ***By, ***Bz;
  double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
  static double *one_dVr, *one_dmu;

  r  = grid[IDIR].x;
  th = grid[JDIR].x;

  if (one_dVr == NULL){
    one_dVr = ARRAY_1D(NX1_TOT, double);  /* -- intercell (i) and (i+1) volume -- */
    one_dmu = ARRAY_1D(NX2_TOT, double);  /* -- intercell (j) and (j+1) volume -- */
    for (i = 0; i < NX1_TOT - 1; i++){
      one_dVr[i] = r[i+1]*fabs(r[i + 1]) - r[i]*fabs(r[i]);
      one_dVr[i] = 2.0/one_dVr[i];
    }
    for (j = 0; j < NX2_TOT - 1; j++){
      one_dmu[j] =    1.0 - cos(th[j + 1]) 
                   - (1.0 - cos(th[j]))*(th[j] > 0.0 ? 1.0:-1.0);
      one_dmu[j] = 1.0/one_dmu[j];
    }
  }

  D_EXPAND(inv_dx  = grid[IDIR].inv_dx; inv_dxi = grid[IDIR].inv_dxi;  ,
           inv_dy  = grid[JDIR].inv_dx; inv_dyi = grid[JDIR].inv_dxi;  ,
           inv_dz  = grid[KDIR].inv_dx; inv_dzi = grid[KDIR].inv_dxi;)

  EXPAND(Bx = V[BX1]; ,
         By = V[BX2]; ,
         Bz = V[BX3]; )

  i = *g_i; j = *g_j; k = *g_k;
  dzBx = dzBy = dyBx = dyBz = dxBy = dxBz = 0.0;

/* ----------------------------------------------
      X1 g_dir:
   ---------------------------------------------- */

  if (g_dir == IDIR){

    D_EXPAND(            ,
      dy_1 = inv_dy[j];  ,
      dz_1 = inv_dz[k];
    )
    s_1 = 1.0/sin(th[j]);
    for (i = 0; i <= NX1_TOT - 2; i++){
      dx_1 = inv_dxi[i];

      #if GEOMETRY == CARTESIAN

/* Forget D_EXPAND and EXPAND MACROS for the moment */
        dxBy = (By[k][j][i+1] - By[k][j][i])*dx_1;
        dxBz = (Bz[k][j][i+1] - Bz[k][j][i])*dx_1;       
                                                       
        dyBx = 0.25*(  Bx[k][j+1][i]   - Bx[k][j-1][i]
                     + Bx[k][j+1][i+1] - Bx[k][j-1][i+1])*dy_1;
        dzBx = 0.25*(  Bx[k+1][j][i]   - Bx[k-1][j][i]
                     + Bx[k+1][j][i+1] - Bx[k-1][j][i+1])*dz_1;             
        dyBz = 0.25*(  Bz[k][j+1][i]   - Bz[k][j-1][i]
                     + Bz[k][j+1][i+1] - Bz[k][j-1][i+1])*dy_1;
        dzBy = 0.25*(  By[k+1][j][i]   - By[k-1][j][i]
                     + By[k+1][j][i+1] - By[k-1][j][i+1])*dz_1;

      #elif GEOMETRY == CYLINDRICAL
      	#error !HALL is not yet compatible with non-cartesian coordinates

     /* -------------------------------------------------------------------
         here dxBz = d(r Bphi)/(r dr) is computed in the following way:
         
                rp*Bp  -  |rm|*Bm
           2 -----------------------
                rp*rp - rm*|rm|
                
          which is obviously the desired derivative away from the axis 
          (where r = |r|) and it gives the correct behavior at the 
          axis, where Bphi --> -Bphi and therefore, 
          
           Bphi = 0  --> Bphi \sim alpha*r --> d(r Bphi)/(r dr) = 2*alpha 
           
          For efficiency purposes we save the geometrical factor 
          1/(rp*rp - rm*|rm|).
         ------------------------------------------------------------------- */

       EXPAND(                                                                ,
         dxBy = (By[k][j][i+1] - By[k][j][i])*dx_1;                           ,
         dxBz = (Bz[k][j][i+1]*r[i+1] - Bz[k][j][i]*fabs(r[i]))*one_dVr[i];) 

       D_EXPAND(                                                 ,
         dyBx = 0.25*( Bx[k][j+1][i]   - Bx[k][j-1][i]
                     + Bx[k][j+1][i+1] - Bx[k][j-1][i+1])*dy_1;  ,
         dzBx = 0.0;) /* -- axisymmetry -- */

      #elif GEOMETRY == POLAR
       	#error !HALL is not yet compatible with non-cartesian coordinates

       r_1  = 1.0/grid[IDIR].xr[i];

       EXPAND(                                             ,
         dxBy = 0.5*(By[k][j][i+1] + By[k][j][i])*r_1 
                  + (By[k][j][i+1] - By[k][j][i])*dx_1;    ,
         dxBz = (Bz[k][j][i+1] - Bz[k][j][i])*dx_1;     )            
 
       D_EXPAND(                                                    ,
         dyBx = 0.25*( Bx[k][j+1][i]   - Bx[k][j-1][i] 
                     + Bx[k][j+1][i+1] - Bx[k][j-1][i+1])*dy_1*r_1; ,
         dzBx = 0.25*( Bx[k+1][j][i]   - Bx[k-1][j][i]
                     + Bx[k+1][j][i+1] - Bx[k-1][j][i+1])*dz_1;     )

      #elif GEOMETRY == SPHERICAL
      	#error !HALL is not yet compatible with non-cartesian coordinates
       r_1 = 1.0/grid[IDIR].xr[i];

       EXPAND(                                           ,
         dxBy = 0.5*(By[k][j][i+1] + By[k][j][i])*r_1
                  + (By[k][j][i+1] - By[k][j][i])*dx_1;  ,
         dxBz = 0.5*(Bz[k][j][i+1] + Bz[k][j][i])*r_1
                  + (Bz[k][j][i+1] - Bz[k][j][i])*dx_1;)

       D_EXPAND(                                                        ,
         dyBx = 0.25*( Bx[k][j+1][i]   - Bx[k][j-1][i]
                     + Bx[k][j+1][i+1] - Bx[k][j-1][i+1])*dy_1*r_1;     ,
         dzBx = 0.25*( Bx[k+1][j][i]   - Bx[k-1][j][i]
                     + Bx[k+1][j][i+1] - Bx[k-1][j][i+1])*dz_1*r_1*s_1; )

      #endif
	  curlB[i][0] = dyBz - dzBy;
      curlB[i][1] = dzBx - dxBz; 
      curlB[i][2] = dxBy - dyBx;                          
    }

/* ----------------------------------------------
      X2 g_dir: 
   ---------------------------------------------- */

  }else if (g_dir == JDIR){

    D_EXPAND(dx_1 = inv_dx[i];  ,
                                ,
             dz_1 = inv_dz[k];)
    r_1  = 1.0/r[i];
    for (j = 0; j <= NX2_TOT - 2; j++){
      dy_1 = inv_dyi[j];

      #if GEOMETRY == CARTESIAN

       dyBx = (Bx[k][j+1][i] - Bx[k][j][i])*dy_1;  
       dyBz = (Bz[k][j+1][i] - Bz[k][j][i])*dy_1;
          
       dxBy = 0.25*(  By[k][j][i+1]   - By[k][j][i-1]
                    + By[k][j+1][i+1] - By[k][j+1][i-1])*dx_1;                                                      
       dzBy = 0.25*(  By[k+1][j][i]   - By[k-1][j][i]
                      + By[k+1][j+1][i] - By[k-1][j+1][i])*dz_1;
       dxBz = 0.25*(  Bz[k][j][i+1]   - Bz[k][j][i-1]
                    + Bz[k][j+1][i+1] - Bz[k][j+1][i-1])*dx_1;        
       dzBx = 0.25*(  Bx[k+1][j][i]   - Bx[k-1][j][i]
                      + Bx[k+1][j+1][i] - Bx[k-1][j+1][i])*dz_1;      
       

      #elif GEOMETRY == CYLINDRICAL
       	#error !HALL is not yet compatible with non-cartesian coordinates
       EXPAND(dyBx = (Bx[k][j+1][i] - Bx[k][j][i])*dy_1;   ,
                                                           ,
              dyBz = (Bz[k][j+1][i] - Bz[k][j][i])*dy_1;)

       dxBy = 0.25*(  By[k][j][i+1]   - By[k][j][i-1]
                    + By[k][j+1][i+1] - By[k][j+1][i-1])*dx_1; 

      #elif GEOMETRY == POLAR

       EXPAND(dyBx = (Bx[k][j+1][i] - Bx[k][j][i])*dy_1*r_1;  ,
                                                              , 
              dyBz = (Bz[k][j+1][i] - Bz[k][j][i])*dy_1*r_1;)

       D_EXPAND(
         dxBy = 0.25*( r[i+1]*(By[k][j][i+1] + By[k][j+1][i+1]) 
                     - r[i-1]*(By[k][j][i-1] + By[k][j+1][i-1]))*dx_1*r_1;  ,
                                                                            ,
         dzBy = 0.25*( By[k+1][j][i]   - By[k-1][j][i]
                     + By[k+1][j+1][i] - By[k-1][j+1][i])*dz_1;)

      #elif GEOMETRY == SPHERICAL
       	#error !HALL is not yet compatible with non-cartesian coordinates
       s_1  = 1.0/grid[JDIR].A[j];

     /* -------------------------------------------------------------------
         here dyBz \propto d(sin(theta) Bphi)/(sin(theta) dtheta) is 
         computed in the following way:
         
               sp*Bp - |sm|*Bm
           2 -----------------------       with s = sin(theta),
               mup - mum*sign(theta)            mu = 1 - cos(theta)
                
          which is obviously the desired derivative away from the axis 
          (where s = |s|) and it gives the correct behavior at the 
          axis, where Bphi --> -Bphi and therefore, 
          
           Bphi = 0  --> Bphi \sim alpha*s --> d(s Bphi)/(dmu) = alpha *(1+cp)
           
          For efficiency purposes we save the geometrical factor
          1/(mup - mum*sign(theta)).
         ------------------------------------------------------------------- */

       EXPAND(dyBx = (Bx[k][j+1][i] - Bx[k][j][i])*dy_1*r_1;   ,
                                                               ,
              dyBz = ( sin(th[j+1])*Bz[k][j + 1][i]
                    - fabs(sin(th[j]))*Bz[k][j][i])*r_1*one_dmu[j];)

       D_EXPAND(
         dxBy = 0.25*( r[i+1]*(By[k][j][i+1] + By[k][j+1][i+1])
                     - r[i-1]*(By[k][j][i-1] + By[k][j+1][i-1])
                                                              )*dx_1*r_1;    ,
                                                                             ,
         dzBy = 0.25*( By[k+1][j][i]   - By[k-1][j][i]
                     + By[k+1][j+1][i] - By[k-1][j+1][i])*dz_1*r_1*s_1; )
      #endif

      curlB[j][0] = dyBz - dzBy; 
      curlB[j][1] = dzBx - dxBz;
      curlB[j][2] = dxBy - dyBx; 
    }

/* ----------------------------------------------
      X3 g_dir: 
   ---------------------------------------------- */

  }else if (g_dir == KDIR) {

    dx_1 = inv_dx[i];
    dy_1 = inv_dy[j];
    r_1  = 1.0/r[i];
    for (k = 0; k <= NX3_TOT - 2; k++){
      dz_1 = inv_dzi[k];
      
      #if GEOMETRY == CARTESIAN

       dzBx = (Bx[k+1][j][i] - Bx[k][j][i])*dz_1;
       dzBy = (By[k+1][j][i] - By[k][j][i])*dz_1;
       dxBz = 0.25*(  Bz[k][j][i+1]   - Bz[k][j][i-1]
                    + Bz[k+1][j][i+1] - Bz[k+1][j][i-1])*dx_1;
       dyBz = 0.25*(  Bz[k][j+1][i]   - Bz[k][j-1][i]
                    + Bz[k+1][j+1][i] - Bz[k+1][j-1][i])*dy_1;
       dxBy = 0.25*(  By[k][j][i+1]   - By[k][j][i-1]
                    + By[k+1][j][i+1] - By[k+1][j][i-1])*dx_1;
       dyBx = 0.25*(  Bx[k][j+1][i]   - Bx[k][j-1][i]
                    + Bx[k+1][j+1][i] - Bx[k+1][j-1][i])*dy_1;
                    

      #elif GEOMETRY == POLAR
             	#error !HALL is not yet compatible with non-cartesian coordinates
       dxBz = 0.25*( Bz[k][j][i+1]   - Bz[k][j][i-1]
                   + Bz[k+1][j][i+1] - Bz[k+1][j][i-1])*dx_1;
       dyBz = 0.25*(  Bz[k][j+1][i]     - Bz[k][j-1][i]
                    + Bz[k+1][j+1][i] - Bz[k+1][j-1][i])*dy_1*r_1;
       dzBx = (Bx[k+1][j][i] - Bx[k][j][i])*dz_1;
       dzBy = (By[k+1][j][i] - By[k][j][i])*dz_1;

      #elif GEOMETRY == SPHERICAL
       	#error !HALL is not yet compatible with non-cartesian coordinates
       s_1  = 1.0/grid[JDIR].A[j];
       dxBz = 0.25*( r[i+1]*(Bz[k][j][i+1] + Bz[k+1][j][i+1])
                   - r[i-1]*(Bz[k][j][i-1] + Bz[k+1][j][i-1]))*dx_1*r_1;
       dyBz = 0.25*( sin(th[j+1])*(Bz[k][j+1][i] + Bz[k+1][j+1][i])
                   - sin(th[j-1])*(Bz[k][j-1][i] + Bz[k+1][j-1][i])
                                                            )*dy_1*r_1*s_1;
       dzBx = (Bx[k+1][j][i] - Bx[k][j][i])*dz_1*r_1*s_1;
       dzBy = (By[k+1][j][i] - By[k][j][i])*dz_1*r_1*s_1;

      #endif
      curlB[k][0] = dyBz - dzBy;              
      curlB[k][1] = dzBx - dxBz;
      curlB[k][2] = dxBy - dyBx;
    }
  }
}

