#include "pluto.h"

static double ***HallJx;
static double ***HallJy;
static double ***HallJz;

#ifdef SHEARINGBOX
/* ********************************************************************* */
void SB_CorrectCurrent (double t, Grid *grid)
/*!
 * Interpolate x and y J and properly correct leftmost
 * and rightmost cells to ensure conservation.
 *
 * \param [in,out] U data array containing cell-centered quantities
 * \param [in] t  the time step at which fluxes have  to be computed
 * \param [in] dt the time step being used in the integrator stage 
 * \param [in] grid pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dtdx;
  static double ***JxL, ***JxR;
  static double ***JyL, ***JyR;
  static double ***JzL, ***JzR;

  if (JxL == NULL){
    JxL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JxR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JyL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JyR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JzL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JzR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
  }

  if (grid[IDIR].lbound != 0){

  /* -- Store J components on the left -- */

    KTOT_LOOP(k) JTOT_LOOP(j){                                    
               JxL[k][j][0] = HallJx[k][j][IBEG - 1];
               JyL[k][j][0] = HallJy[k][j][IBEG - 1];
               JzL[k][j][0] = HallJz[k][j][IBEG - 1];
    }
  }

  if (grid[IDIR].rbound != 0){

  /* -- Store J components on the right -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
               JxR[k][j][0] = HallJx[k][j][IEND];
               JyR[k][j][0] = HallJy[k][j][IEND];
               JzR[k][j][0] = HallJz[k][j][IEND];
    }

  }



  #ifdef PARALLEL
   		ExchangeX (JxL[0][0], JxR[0][0], NX3_TOT*NX2_TOT, grid);
   		ExchangeX (JyL[0][0], JyR[0][0], NX3_TOT*NX2_TOT, grid); 
   		ExchangeX (JzL[0][0], JzR[0][0], NX3_TOT*NX2_TOT, grid); 
  #endif

/* --------------------------------------------------------------------- */
/*! \note 
    In parallel and if there's more than one processor in the
    x direction, we modify the values of FluxL and FluxR
    independently since the two boundary sides are handled by
    different processors.
    Conversely, when a processor owns both sides, we need to
    copy at least one of the two arrays before doing interpolation
    since their original values are lost when setting boundary
    conditions.                                                          */
/* --------------------------------------------------------------------- */

  if (grid[IDIR].lbound != 0){  /* ---- Left x-boundary ---- */

    RBox box;
    static double ***Jx1;
    static double ***Jy1;
    static double ***Jz1;

  /* -- set grid ranges of FluxR and exchange b.c. -- */

    box.ib = 0; box.ie = 0; 
    box.jb = 0; box.je = NX2_TOT-1; 
    box.kb = 0; box.ke = NX3_TOT-1;

  /* -- make a copy of FluxR to avoid overwriting one nproc_x = 1 -- */
	
	
       
    if (grid[IDIR].nproc == 1){
      if (Jx1 == NULL) {
      	  Jx1 = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
      	  Jy1 = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
      	  Jz1 = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
      }
      BOX_LOOP((&box), k, j, i) {
       		 Jx1[k][j][i] = JxR[k][j][i];
       		 Jy1[k][j][i] = JyR[k][j][i];
       		 Jz1[k][j][i] = JzR[k][j][i];
       }
       
    }
    else {
      Jx1 = JxR;
      Jy1 = JyR;
      Jz1 = JzR;
    }  
  /* -- set boundary conditions on Jx1, Jy1, Jz1 -- */

    SB_SetBoundaryVar(Jx1, &box, X1_BEG, t, grid);
    SB_SetBoundaryVar(Jy1, &box, X1_BEG, t, grid);
    SB_SetBoundaryVar(Jz1, &box, X1_BEG, t, grid);

    for (k = KBEG; k <= KEND; k++){
       for (j = JBEG; j <= JEND; j++) {
       	  HallJx[k][j][IBEG - 1] = 0.5*(JxL[k][j][0] + Jx1[k][j][0]);
       	  HallJy[k][j][IBEG - 1] = 0.5*(JyL[k][j][0] + Jy1[k][j][0]);
       	  HallJz[k][j][IBEG - 1] = 0.5*(JzL[k][j][0] + Jz1[k][j][0]);
      }
    }
  }   

  if (grid[IDIR].rbound != 0){  /* ---- Right x-boundary ---- */

  /* -- set grid ranges of JxR... and exchange b.c. -- */

    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;


    SB_SetBoundaryVar(JxL, &box, X1_END, t, grid);
    SB_SetBoundaryVar(JyL, &box, X1_END, t, grid);
    SB_SetBoundaryVar(JzL, &box, X1_END, t, grid);
    
    for (k = KBEG; k <= KEND; k++){
        for (j = JBEG; j <= JEND; j++) {
        	HallJx[k][j][IEND] = 0.5*(JxL[k][j][0] + JxR[k][j][0]);
        	HallJy[k][j][IEND] = 0.5*(JyL[k][j][0] + JyR[k][j][0]);
        	HallJz[k][j][IEND] = 0.5*(JzL[k][j][0] + JzR[k][j][0]);
        }
    }
  }
}
#endif





void ComputeJ(const Data *d, Grid *grid, double t) {
	int i,j,k;
	static int first_call = 1;
	double dx_1, dy_1, dz_1;
    double *inv_dx,  *inv_dy,  *inv_dz;
	double ***Bx,  ***By,  ***Bz;
	double ***Bxs, ***Bys, ***Bzs;
	double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
	
	if(first_call) {
		first_call=0;
		HallJx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);	
		HallJy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
		HallJz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);	
	}
	
	/* Compute the current state->j associated to current sweep as a face-centered quantity */

	
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
		for(k = 1 ; k < NX3_TOT-1 ; k++) {
			for(j = 1 ; j < NX2_TOT-1 ; j++) {
				for(i = 1 ; i < NX1_TOT-1 ; i++) {
					dx_1 = inv_dx[i];
        			dy_1 = inv_dy[j];
        			dz_1 = inv_dz[k];
        	
					dxBy = (By[k][j][i+1]-By[k][j][i])*dx_1;
					dxBz = (Bz[k][j][i+1]-Bz[k][j][i])*dx_1;
			
					dyBx = 0.5*(Bxs[k][j+1][i]-Bxs[k][j-1][i])*dy_1;
					dyBz = 0.25*(Bz[k][j+1][i]+Bz[k][j+1][i+1] - Bz[k][j-1][i]-Bz[k][j-1][i+1])*dy_1;
			
					dzBx = 0.5*(Bxs[k+1][j][i]-Bxs[k-1][j][i])*dz_1;
					dzBy = 0.25*(By[k+1][j][i]+By[k+1][j][i+1] - By[k-1][j][i]-By[k-1][j][i+1])*dz_1;
			
					HallJx[k][j][i] = dyBz-dzBy;
					HallJy[k][j][i] = dzBx-dxBz;
					HallJz[k][j][i] = dxBy-dyBx;
				}
			}
		}
#ifdef SHEARINGBOX
		SB_CorrectCurrent (t, grid);
#endif
	}
	
	if (g_dir == JDIR) {
		for(k = 1 ; k < NX3_TOT-1 ; k++) {
			for(j = 1 ; j < NX2_TOT-1 ; j++) {
				for(i = 1 ; i < NX1_TOT-1 ; i++) {
					dx_1 = inv_dx[i];
        			dy_1 = inv_dy[j];
        			dz_1 = inv_dz[k];
        	
					dxBy = 0.5*(Bys[k][j][i+1]-Bys[k][j][i-1])*dx_1;
					dxBz = 0.25*(Bz[k][j][i+1]+Bz[k][j+1][i+1] - Bz[k][j][i-1]-Bz[k][j+1][i-1])*dx_1;
			
					dyBx = (Bx[k][j+1][i]-Bx[k][j][i])*dy_1;
					dyBz = (Bz[k][j+1][i]-Bz[k][j][i])*dy_1;
			
					dzBx = 0.25*(Bx[k+1][j][i]+Bx[k+1][j+1][i] - Bx[k-1][j][i]-Bx[k-1][j+1][i])*dz_1;
					dzBy = 0.5*(Bys[k+1][j][i]-Bys[k-1][j][i])*dz_1;
			
					HallJx[k][j][i] = dyBz-dzBy;
					HallJy[k][j][i] = dzBx-dxBz;
					HallJz[k][j][i] = dxBy-dyBx;
				}
			}
		}
	}
	
	else if (g_dir == KDIR) {
		for(k = 1 ; k < NX3_TOT-1 ; k++) {
			for(j = 1 ; j < NX2_TOT-1 ; j++) {
				for(i = 1 ; i < NX1_TOT-1 ; i++) {
					dx_1 = inv_dx[i];
        			dy_1 = inv_dy[j];
        			dz_1 = inv_dz[k];
        	
					dxBy = 0.25*(By[k][j][i+1]+By[k+1][j][i+1] - By[k][j][i-1]-By[k+1][j][i-1])*dx_1;
					dxBz = 0.5*(Bzs[k][j][i+1]-Bzs[k][j][i-1])*dx_1;
						
					dyBx = 0.25*(Bx[k][j+1][i]+Bx[k+1][j+1][i] - Bx[k][j-1][i]-Bx[k+1][j-1][i])*dy_1;
					dyBz = 0.5*(Bzs[k][j+1][i]-Bzs[k][j-1][i])*dy_1;
			
					dzBx = (Bx[k+1][j][i]-Bx[k][j][i])*dz_1;
					dzBy = (By[k+1][j][i]-By[k][j][i])*dz_1;
			
					HallJx[k][j][i] = dyBz-dzBy;
					HallJy[k][j][i] = dzBx-dxBz;
					HallJz[k][j][i] = dxBy-dyBx;
				}
			}
		}
	}
	
}

			


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
void StoreJState(State_1D *state, int *in, int *i, int *j, int *k, int beg, int end) {
	/* Compute the current state->j associated to current sweep as a face-centered quantity */
     
    for(*in=beg ; *in <= end ; (*in)++) {		
		state->j[*in][JX1] = HallJx[*k][*j][*i];
		state->j[*in][JX2] = HallJy[*k][*j][*i];
		state->j[*in][JX3] = HallJz[*k][*j][*i];
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

