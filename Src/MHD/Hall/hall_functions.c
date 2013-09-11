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
void ComputeJState(const Data *d, const Grid *grid, State_1D *state, int *in, int *i, int *j, int *k, int beg, int end) {
	/* Compute the current state->j associated to current sweep as a face-centered quantity */
	double dx_1, dy_1, dz_1;
    double *inv_dx,  *inv_dy,  *inv_dz;
	double ***Bx,  ***By,  ***Bz;
	double ***Bxs, ***Bys, ***Bzs;
	
	double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
	
	EXPAND(Bx = d->Vc[BX1]; ,
           By = d->Vc[BX2]; ,
           Bz = d->Vc[BX3]; )
           
    EXPAND(Bxs = d->Vs[BX1s]; ,
           Bys = d->Vs[BX2s]; ,
           Bzs = d->Vs[BX3s]; )   
           
    D_EXPAND(inv_dx  = grid[IDIR].inv_dx;  ,
             inv_dy  = grid[JDIR].inv_dx;  ,
             inv_dz  = grid[KDIR].inv_dx;)
                   
	if (g_dir == IDIR) {
		
		for(*in=beg ; *in <= end ; (*in)++) {
			dx_1 = inv_dx[*i];
        	dy_1 = inv_dy[*j];
        	dz_1 = inv_dz[*k];
        	
			dxBy = (By[*k][*j][*i+1]-By[*k][*j][*i])*dx_1;
			dxBz = (Bz[*k][*j][*i+1]-Bz[*k][*j][*i])*dx_1;
			
			dyBx = 0.5*(Bxs[*k][*j+1][*i]-Bxs[*k][*j-1][*i])*dy_1;
			dyBz = 0.25*(Bz[*k][*j+1][*i]+Bz[*k][*j+1][*i+1] - Bz[*k][*j-1][*i]-Bz[*k][*j-1][*i+1])*dy_1;
			
			dzBx = 0.5*(Bxs[*k+1][*j][*i]-Bxs[*k-1][*j][*i])*dz_1;
			dzBy = 0.25*(By[*k+1][*j][*i]+By[*k+1][*j][*i+1] - By[*k-1][*j][*i]-By[*k-1][*j][*i+1])*dz_1;
			
			state->j[*in][JX1] = dyBz-dzBy;
			state->j[*in][JX2] = dzBx-dxBz;
			state->j[*in][JX3] = dxBy-dyBx;

		} 
	}	
		
	else if (g_dir == JDIR) {
		for(*in=beg ; *in <= end ; (*in)++) {
			dx_1 = inv_dx[*i];
        	dy_1 = inv_dy[*j];
        	dz_1 = inv_dz[*k];
        	
			dxBy = 0.5*(Bys[*k][*j][*i+1]-Bys[*k][*j][*i-1])*dx_1;
			dxBz = 0.25*(Bz[*k][*j][*i+1]+Bz[*k][*j+1][*i+1] - Bz[*k][*j][*i-1]-Bz[*k][*j+1][*i-1])*dx_1;
			
			dyBx = (Bx[*k][*j+1][*i]-Bx[*k][*j][*i])*dy_1;
			dyBz = (Bz[*k][*j+1][*i]-Bz[*k][*j][*i])*dy_1;
			
			dzBx = 0.25*(Bx[*k+1][*j][*i]+Bx[*k+1][*j+1][*i] - Bx[*k-1][*j][*i]-Bx[*k-1][*j+1][*i])*dz_1;
			dzBy = 0.5*(Bys[*k+1][*j][*i]-Bys[*k-1][*j][*i])*dz_1;
			
			state->j[*in][JX1] = dyBz-dzBy;
			state->j[*in][JX2] = dzBx-dxBz;
			state->j[*in][JX3] = dxBy-dyBx;
			
		} 
	}
	
	else if (g_dir == KDIR) {
		for(*in=beg ; *in <= end ; (*in)++) {
			dx_1 = inv_dx[*i];
        	dy_1 = inv_dy[*j];
        	dz_1 = inv_dz[*k];
        	
			dxBy = 0.25*(By[*k][*j][*i+1]+By[*k+1][*j][*i+1] - By[*k][*j][*i-1]-By[*k+1][*j][*i-1])*dx_1;
			dxBz = 0.5*(Bzs[*k][*j][*i+1]-Bzs[*k][*j][*i-1])*dx_1;
						
			dyBx = 0.25*(Bx[*k][*j+1][*i]+Bx[*k+1][*j+1][*i] - Bx[*k][*j-1][*i]-Bx[*k+1][*j-1][*i])*dy_1;
			dyBz = 0.5*(Bzs[*k][*j+1][*i]-Bzs[*k][*j-1][*i])*dy_1;
			
			dzBx = (Bx[*k+1][*j][*i]-Bx[*k][*j][*i])*dz_1;
			dzBy = (By[*k+1][*j][*i]-By[*k][*j][*i])*dz_1;
			
			state->j[*in][JX1] = dyBz-dzBy;
			state->j[*in][JX2] = dzBx-dxBz;
			state->j[*in][JX3] = dxBy-dyBx;

		}
	}
	
}


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

