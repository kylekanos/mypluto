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

void ComputeCurrent(Data_Arr V, Grid *grid) {
    int  i, j, k;
    double dx_1, dy_1, dz_1;
    double *inv_dx,  *inv_dy,  *inv_dz;
    double ***Bx, ***By, ***Bz;
    double ***Jx, ***Jy, ***Jz;
    double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
    
    D_EXPAND(inv_dx  = grid[IDIR].inv_dx;  ,
             inv_dy  = grid[JDIR].inv_dx;  ,
             inv_dz  = grid[KDIR].inv_dx;)
    
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
        
        Jx[k][j][i] = (dyBz - dzBy);
        Jy[k][j][i] = (dzBx - dxBz);
        Jz[k][j][i] = (dxBy - dyBx);
    }
    
}

