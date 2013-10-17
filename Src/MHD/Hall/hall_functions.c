#include "pluto.h"

static double ***HallJx;
static double ***HallJy;
static double ***HallJz;

#ifdef SHEARINGBOX

/* ********************************************************************* */
void Hall_SetCurrentRBox(RBox *center, RBox *x1face, RBox *x2face, RBox *x3face)
/* 
 *Specially adapted for current BC which are face centered but which do not
 *Fall into the usual face centered array layout.
 *
 *********************************************************************** */
{
  int s;

/* ---------------------------------------------------
            set X1_BEG grid index ranges
   --------------------------------------------------- */

  s = X1_BEG; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib = IBEG-1; center[s].ie =         0;
  center[s].jb =      0; center[s].je = NX2_TOT-1;
  center[s].kb =      0; center[s].ke = NX3_TOT-1;

  x1face[s] = x2face[s] = x3face[s] = center[s];

   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;

   x1face[s].ib--; 

/* ---------------------------------------------------
            set X1_END grid index ranges
   --------------------------------------------------- */
  
  s = X1_END; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib = IEND+1; center[s].ie = NX1_TOT-1;
  center[s].jb =      0; center[s].je = NX2_TOT-1;
  center[s].kb =      0; center[s].ke = NX3_TOT-1;

  x1face[s] = x2face[s] = x3face[s] = center[s];

  #ifndef CH_SPACEDIM /* -- useless for AMR -- */
   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;

  #endif

/* ---------------------------------------------------
            set X2_BEG grid index ranges
   --------------------------------------------------- */

  s = X2_BEG; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib =      0; center[s].ie = NX1_TOT-1;
  center[s].jb = JBEG-1; center[s].je =         0;
  center[s].kb =      0; center[s].ke = NX3_TOT-1;

  x1face[s] = x2face[s] = x3face[s] = center[s];

  #ifndef CH_SPACEDIM /* -- useless for AMR -- */
   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;


   x2face[s].jb--;
   
  #endif

/* ---------------------------------------------------
            set X2_END grid index ranges
   --------------------------------------------------- */
  
  s = X2_END; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib =      0; center[s].ie = NX1_TOT-1;
  center[s].jb = JEND+1; center[s].je = NX2_TOT-1;
  center[s].kb =      0; center[s].ke = NX3_TOT-1;

  x1face[s] = x2face[s] = x3face[s] = center[s];

  #ifndef CH_SPACEDIM /* -- useless for AMR -- */
   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;

  #endif

/* ---------------------------------------------------
            set X3_BEG grid index ranges
   --------------------------------------------------- */

  s = X3_BEG; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib =      0; center[s].ie = NX1_TOT-1;
  center[s].jb =      0; center[s].je = NX2_TOT-1;
  center[s].kb = KBEG-1; center[s].ke =         0;

  x1face[s] = x2face[s] = x3face[s] = center[s];

  #ifndef CH_SPACEDIM /* -- useless for AMR -- */
   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;


   x3face[s].kb--;
  #endif

/* ---------------------------------------------------
            set X3_END grid index ranges
   --------------------------------------------------- */
  
  s = X3_END; s -= X1_BEG;

  center[s].vpos = CENTER;

  center[s].ib =      0; center[s].ie = NX1_TOT-1; center[s].di = 1;
  center[s].jb =      0; center[s].je = NX2_TOT-1; center[s].dj = 1;
  center[s].kb = KEND+1; center[s].ke = NX3_TOT-1; center[s].dk = 1;

  x1face[s] = x2face[s] = x3face[s] = center[s];

  #ifndef CH_SPACEDIM /* -- useless for AMR -- */
   x1face[s].vpos = X1FACE;
   x2face[s].vpos = X2FACE;
   x3face[s].vpos = X3FACE;

  #endif
}


/* ********************************************************************* */
void Hall_SB_Boundary (double ***Jin, int side, Grid *grid) 
/*! 
 * Main wrapper function used to assign shearing-box boundary conditions
 * on flow variables.
 *
 * \param d    pointer to the PLUTO Data structure
 * \param side the side of the computational domain (X1_BEG or X1_END) 
 * \param grid pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *
 * \todo Check if sb_Ly and sb_vy need to be global.
 *********************************************************************** */
{
  int    i, j, k, nv;
  double t;
  RBox   box;

  t  = g_time;
  #if TIME_STEPPING == RK2
   if (g_intStage == 2) t = g_time + g_dt;
  #endif
/* -------------------------------------------------
                  X1 Beg Boundary
   ------------------------------------------------- */

  if (side == X1_BEG){
      box.ib =  0; box.ie = IBEG-1;
      box.jb = 0;  box.je = NX2_TOT-1;
      box.kb =  0; box.ke = NX3_TOT-1;
      SB_SetBoundaryVar(Jin, &box, side, t, grid);
  } /* -- END side X1_BEG -- */

/* -------------------------------------------------
                  X1 End Boundary
   ------------------------------------------------- */

  if (side == X1_END){

     box.ib = IEND+1; box.ie = NX1_TOT-1;
     box.jb =     0; box.je = NX2_TOT-1;
     box.kb =      0; box.ke = NX3_TOT-1;
     SB_SetBoundaryVar(Jin, &box, side, t, grid);
   
  }
}


/* ********************************************************************* */
void Hall_JBoundary (double ***Jin, Grid *grid)
/*!
 * Set boundary conditions on one or more sides of the computational
 * domain.
 *
 * \param [in,out] d  pointer to PLUTO Data structure containing the 
 *                    solution array (including centered and staggered
 *                    fields)
 * \param [in]   idim specifies on which side(s) of the computational
 *               domain boundary conditions must be set. Possible values
 *               are  
 *        - idim = IDIR   first dimension (x1)
 *        - idim = JDIR   second dimenson (x2)
 *        - idim = KDIR   third dimension (x3)
 *        - idim = ALL_DIR all dimensions
 * \param [in]  grid   pointer to an array of grid structures.
 ******************************************************************* */
{
  int  is, nv,idim;
  int  side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int  type[6], sbeg, send, vsign[NVAR];
  int  par_dim[3] = {0, 0, 0};
  static int first_call = 1;
  double ***q;
  static RBox center[8], x1face[8], x2face[8], x3face[8];

/* We only do periodic BCs in the x direction */

  idim = IDIR;
/* -----------------------------------------------------
     Set the boundary boxes on the six domain sides
   ----------------------------------------------------- */

  #ifndef CH_SPACEDIM
  if (first_call){
    Hall_SetCurrentRBox(center, x1face, x2face, x3face);
    first_call = 0;
  }
  #else /* -- with dynamic grids we need to re-define the RBox at each time -- */
   Hall_SetCurrentRBox(center, x1face, x2face, x3face);
  #endif

/* ---------------------------------------------------
    Check the number of processors in each direction
    Only perform transfers in x if required
   --------------------------------------------------- */

  D_EXPAND(par_dim[0] = grid[IDIR].nproc > 1;  ,
           par_dim[1] = 0;  ,
           par_dim[2] = 0;)

  
/* -------------------------------------
     Exchange data between processors 
   ------------------------------------- */
      
      
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
/*   NOT NEEDED
   if(g_dir==IDIR)   
   	AL_Exchange_dim ((char *)(Jin[0][0] - 1), par_dim, SZ_stagx); */
   if(g_dir==JDIR)
    AL_Exchange_dim ((char *)Jin[0][0]     , par_dim, SZ);
   if(g_dir==KDIR)
   	AL_Exchange_dim ((char *)Jin[0][0]     , par_dim, SZ);
   	
   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* ----------------------------------------------------------------
     When idim == ALL_DIR boundaries are imposed on ALL sides:
     a loop from sbeg = 0 to send = 2*DIMENSIONS - 1 is performed. 
     When idim = n, boundaries are imposed at the beginning and 
     the end of the i-th direction.
   ---------------------------------------------------------------- */ 

  if (idim == ALL_DIR) {
    sbeg = 0;
    send = 2*DIMENSIONS - 1;
  } else {
    sbeg = 2*idim;
    send = 2*idim + 1;
  }

/* --------------------------------------------------------
        Main loop on computational domain sides
   -------------------------------------------------------- */

  type[0] = grid[IDIR].lbound; type[1] = grid[IDIR].rbound;
  type[2] = grid[JDIR].lbound; type[3] = grid[JDIR].rbound;
  type[4] = grid[KDIR].lbound; type[5] = grid[KDIR].rbound;

  for (is = sbeg; is <= send; is++){


    if (type[is] == SHEARING) {  /* -- shearingbox boundary -- */

  /* ---------------------------------------------------------
      SHEARING-BOX boundary condition is implemented as

      1) apply periodic boundary conditions for all variables
        (except staggered BX)
      2) Perform spatial shift in the y-direction
     --------------------------------------------------------- */

       if (side[is] != X1_BEG && side[is] != X1_END){
         print1 ("! BOUNDARY: shearingbox can only be assigned at an X1 boundary\n");
         QUIT_PLUTO(1);
       }
       if (grid[IDIR].nproc == 1){
/*    NOT NEEDED
       	 if(g_dir==IDIR) 
       	 	PeriodicBound (Jin, x1face+is, side[is]);   */
       	 if(g_dir==JDIR) 
       	 	PeriodicBound (Jin, x2face+is, side[is]);
       	 if(g_dir==KDIR) 
       	 	PeriodicBound (Jin, x3face+is, side[is]);
       }
       Hall_SB_Boundary (Jin, side[is], grid);

    }
  }
}




/* ********************************************************************* */
void Hall_FixJ (double ***Jin, double t, Grid *grid)
/*!
 * Interpolate J and properly correct leftmost
 * and rightmost cells to ensure Shearing sheet BCs are satisfied.
 *
 * \param [in,out] U data array containing cell-centered quantities
 * \param [in] t  the time step at which fluxes have  to be computed
 * \param [in] dt the time step being used in the integrator stage 
 * \param [in] grid pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dtdx;
  static double ***JL, ***JR;

  if (JL == NULL){
    JL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    JR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
  }

  if (grid[IDIR].lbound != 0){

  /* -- Store J components on the left -- */

    KTOT_LOOP(k) JTOT_LOOP(j){                                    
               JL[k][j][0] = Jin[k][j][IBEG - 1];
    }
  }

  if (grid[IDIR].rbound != 0){

  /* -- Store J components on the right -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
               JR[k][j][0] = Jin[k][j][IEND];
    }

  }

  #ifdef PARALLEL
   		ExchangeX (JL[0][0], JR[0][0], NX3_TOT*NX2_TOT, grid); 
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
    static double ***J1;

  /* -- set grid ranges of FluxR and exchange b.c. -- */

    box.ib = 0; box.ie = 0; 
    box.jb = 0; box.je = NX2_TOT-1; 
    box.kb = 0; box.ke = NX3_TOT-1;

  /* -- make a copy of JR to avoid overwriting one nproc_x = 1 -- */
	
	
       
    if (grid[IDIR].nproc == 1){
      if (J1 == NULL) {
      	  J1 = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
      }
      BOX_LOOP((&box), k, j, i) {
       		 J1[k][j][i] = JR[k][j][i];
       }
       
    }
    else {
      J1 = JR;
    }  
  /* -- set boundary conditions on J1 -- */

    SB_SetBoundaryVar(J1, &box, X1_BEG, t, grid);

    for (k = 0; k < NX3_TOT; k++){
       for (j = 0; j < NX2_TOT; j++) {
       	  Jin[k][j][IBEG - 1] = 0.5*(JL[k][j][0] + J1[k][j][0]);
      }
    }
  }   

  if (grid[IDIR].rbound != 0){  /* ---- Right x-boundary ---- */

  /* -- set grid ranges of JxR... and exchange b.c. -- */

    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;


    SB_SetBoundaryVar(JL, &box, X1_END, t, grid);
    
    for (k = 0; k < NX3_TOT; k++){
       for (j = 0; j < NX2_TOT; j++) {
        	Jin[k][j][IEND] = 0.5*(JL[k][j][0] + JR[k][j][0]);
        }
    }
  }
}
#endif



#if HALL_MHD==RIEMANN

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
		// Symmetrize currents at x-cell faces
		Hall_FixJ (HallJx, t, grid);
		Hall_FixJ (HallJy, t, grid);
		Hall_FixJ (HallJz, t, grid);
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
#ifdef SHEARINGBOX
    Hall_JBoundary(HallJx, grid);
    Hall_JBoundary(HallJy, grid);
    Hall_JBoundary(HallJz, grid);
#endif
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
#ifdef SHEARINGBOX
    Hall_JBoundary(HallJx, grid);
    Hall_JBoundary(HallJy, grid);
    Hall_JBoundary(HallJz, grid); 
#endif
	}
}

		

void StoreJState(State_1D *state, int *in, int *i, int *j, int *k, int beg, int end) {
	/* Compute the current state->j associated to current sweep as a face-centered quantity */
     
    // Currents are automatically computed on the right side.
    for(*in=beg ; *in <= end ; (*in)++) {		
		state->j[*in][JX1] = HallJx[*k][*j][*i];
		state->j[*in][JX2] = HallJy[*k][*j][*i];
		state->j[*in][JX3] = HallJz[*k][*j][*i];
	}	
}

#endif	// HALL_MHD == RIEMANN

// Compute Hall EMFS required by Hall as a source term
void Hall_emf(EMF *emf, Data *d, Grid *grid) {
	static int first_call = 1;
	double dx_1, dy_1, dz_1;
    double *dx,  *dy,  *dz;
	double dxBy, dxBz, dyBx, dyBz, dzBx, dzBy;
	double v[NVAR];
	double lhall;
	
	double ***Bx,  ***By,  ***Bz;
	double ***Bxs, ***Bys, ***Bzs;
	
	static double ***Jx;
	static double ***Jy;
	static double ***Jz;
	static double ***Bxe;
	static double ***Bye;
	static double ***Bze;
	
	int i,j,k,n;
	
	EXPAND(Bx = d->Vc[BX1]; ,
           By = d->Vc[BX2]; ,
           Bz = d->Vc[BX3]; )
           
    EXPAND(Bxs = d->Vs[BX1s]; ,
           Bys = d->Vs[BX2s]; ,
           Bzs = d->Vs[BX3s]; )   
           
    D_EXPAND(dx  = grid[IDIR].dx;  ,
             dy  = grid[JDIR].dx;  ,
             dz  = grid[KDIR].dx;)
             
	
	if(first_call) {
		first_call=0;
		Jx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);	
		Jy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
		Jz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
		
		Bxe = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);	
		Bye = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
		Bze = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);	
	}
	
	// We start by computing the emf in the x direction (located on y-z edge)
	for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	dx_1=1.0/dx[i];
    	dy_1=1.0/(0.5*(dy[j]+dy[j+1]));
    	dz_1=1.0/(0.5*(dz[k]+dz[k+1]));
    	
    	dxBy = 0.25*(Bys[k][j][i+1]+Bys[k+1][j][i+1]-Bys[k][j][i-1]-Bys[k+1][j][i-1])*dx_1;
    	dxBz = 0.25*(Bzs[k][j][i+1]+Bzs[k][j+1][i+1]-Bzs[k][j][i-1]-Bzs[k][j+1][i-1])*dx_1;
    	dyBx = 0.5*(Bx[k][j+1][i]+Bx[k+1][j+1][i]-Bx[k][j][i]-Bx[k+1][j][i])*dy_1;
    	dzBx = 0.5*(Bx[k+1][j][i]+Bx[k+1][j+1][i]-Bx[k][j][i]-Bx[k][j+1][i])*dz_1;
    	
    	Jy[k][j][i]=dzBx-dxBz;
    	Jz[k][j][i]=dxBy-dyBx;
    	
    	Bye[k][j][i]=0.5*(Bys[k][j][i]+Bys[k+1][j][i]);
    	Bze[k][j][i]=0.5*(Bzs[k][j][i]+Bzs[k][j+1][i]);
    }}}
    // For y-z edge, there is no need to symmetrize the resuting fields
    
    for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	for(n=0 ; n < NVAR ; n++) 
    		v[n]=0.25*(d->Vc[n][k][j][i]+d->Vc[n][k+1][j][i]+d->Vc[n][k][j+1][i]+d->Vc[n][k+1][j+1][i]);
		lHall_Func(v, grid[IDIR].x[i], 0.5*(grid[JDIR].x[j]+grid[JDIR].x[j+1]), 0.5*(grid[KDIR].x[k]+grid[KDIR].x[k+1]), &lhall);
		
		emf->ex[k][j][i] += lhall*(Jy[k][j][i]*Bze[k][j][i]-Jz[k][j][i]*Bye[k][j][i]);
    }}}
    
    // EMF in the y direction (located on x-z edge)
	for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	dx_1=1.0/(0.5*(dx[i]+dx[i+1]));
    	dy_1=1.0/dy[j];
    	dz_1=1.0/(0.5*(dz[k]+dz[k+1]));
    	
    	dyBx = 0.25*(Bxs[k][j+1][i]+Bxs[k+1][j+1][i]-Bxs[k][j-1][i]-Bxs[k+1][j-1][i])*dy_1;
    	dyBz = 0.25*(Bzs[k][j+1][i]+Bzs[k][j+1][i+1]-Bzs[k][j-1][i]-Bzs[k][j-1][i+1])*dy_1;
    	dxBy = 0.5*(By[k][j][i+1]+By[k+1][j][i+1]-By[k][j][i]-By[k+1][j][i])*dx_1;
    	dzBy = 0.5*(By[k+1][j][i]+By[k+1][j][i+1]-By[k][j][i]-By[k][j][i+1])*dz_1;
    	
    	Jx[k][j][i]=dyBz-dzBy;
    	Jz[k][j][i]=dxBy-dyBx;
    	
    	Bxe[k][j][i]=0.5*(Bxs[k][j][i]+Bxs[k+1][j][i]);
    	Bze[k][j][i]=0.5*(Bzs[k][j][i]+Bzs[k][j][i+1]);
    }}}
    // For x-z edges, one has to symmetrize edge-centered fields
#ifdef SHEARINGBOX
    Hall_FixJ (Jx, g_time, grid);
    Hall_FixJ (Jz, g_time, grid);
    // Use similar procedure for edge-centered B, since J and Be are at the same location
    Hall_FixJ (Bxe, g_time, grid);
    Hall_FixJ (Bze, g_time, grid);
#endif    
    
    for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	for(n=0 ; n < NVAR ; n++) 
    		v[n]=0.25*(d->Vc[n][k][j][i]+d->Vc[n][k+1][j][i]+d->Vc[n][k][j][i+1]+d->Vc[n][k+1][j][i+1]);
		lHall_Func(v, 0.5*(grid[IDIR].x[i]+grid[IDIR].x[i+1]), grid[JDIR].x[j], 0.5*(grid[KDIR].x[k]+grid[KDIR].x[k+1]), &lhall);
		
		emf->ey[k][j][i] += lhall*(Jz[k][j][i]*Bxe[k][j][i]-Jx[k][j][i]*Bze[k][j][i]);
    }}}
    
    // EMF in the z direction (located on x-y edge)
	for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	dx_1=1.0/(0.5*(dx[i]+dx[i+1]));
    	dy_1=1.0/(0.5*(dy[j]+dy[j+1]));
    	dz_1=dz[k];
    	
    	dzBx = 0.25*(Bxs[k+1][j][i]+Bxs[k+1][j+1][i]-Bxs[k-1][j][i]-Bxs[k-1][j+1][i])*dz_1;
    	dzBy = 0.25*(Bys[k+1][j][i]+Bys[k+1][j][i+1]-Bys[k-1][j][i]-Bys[k-1][j][i+1])*dz_1;
    	dxBz = 0.5*(Bz[k][j][i+1]+Bz[k][j+1][i+1]-Bz[k][j][i]-Bz[k][j+1][i])*dx_1;
    	dyBz = 0.5*(Bz[k][j+1][i]+Bz[k][j+1][i+1]-Bz[k][j][i]-Bz[k][j][i+1])*dy_1;
    	
    	Jx[k][j][i]=dyBz-dzBy;
    	Jy[k][j][i]=dzBx-dxBz;
    	
    	Bxe[k][j][i]=0.5*(Bxs[k][j][i]+Bxs[k][j+1][i]);
    	Bye[k][j][i]=0.5*(Bys[k][j][i]+Bys[k][j][i+1]);
    }}}
    // For x-y edges, one has to symmetrize edge-centered fields
#ifdef SHEARINGBOX
    Hall_FixJ (Jx, g_time, grid);
    Hall_FixJ (Jy, g_time, grid);
    // Use similar procedure for edge-centered B, since J and Be are at the same location
    Hall_FixJ (Bxe, g_time, grid);
    Hall_FixJ (Bye, g_time, grid);
#endif    
    
    for (k = emf->kbeg; k <= emf->kend; k++){
    for (j = emf->jbeg; j <= emf->jend; j++){
    for (i = emf->ibeg; i <= emf->iend; i++){ 
    	for(n=0 ; n < NVAR ; n++) 
    		v[n]=0.25*(d->Vc[n][k][j][i]+d->Vc[n][k][j+1][i]+d->Vc[n][k][j][i+1]+d->Vc[n][k][j+1][i+1]);
		lHall_Func(v, 0.5*(grid[IDIR].x[i]+grid[IDIR].x[i+1]), 0.5*(grid[JDIR].x[j]+grid[JDIR].x[j+1]), grid[KDIR].x[k], &lhall);
		
		emf->ez[k][j][i] += lhall*(Jx[k][j][i]*Bye[k][j][i]-Jy[k][j][i]*Bxe[k][j][i]);
    }}}
	
}


