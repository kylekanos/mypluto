#include "pluto.h"


#if HALL_MHD = SOURCE
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
    Curr_FixJ (Jx, g_time, grid);
    Curr_FixJ (Jz, g_time, grid);
    // Use similar procedure for edge-centered B, since J and Be are at the same location
    Curr_FixJ (Bxe, g_time, grid);
    Curr_FixJ (Bze, g_time, grid);
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
    Curr_FixJ (Jx, g_time, grid);
    Curr_FixJ (Jy, g_time, grid);
    // Use similar procedure for edge-centered B, since J and Be are at the same location
    Curr_FixJ (Bxe, g_time, grid);
    Curr_FixJ (Bye, g_time, grid);
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

#endif		// HALL_MHD == SOURCE

