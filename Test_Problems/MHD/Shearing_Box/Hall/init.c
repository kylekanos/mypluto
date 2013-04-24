#include "pluto.h"

double sb_Omega = 1.0;
double sb_A     = -0.75;



/*********************************************/
/**
 Customized random number generator
 Allow one to have consistant random numbers
 generators on different architectures.
 
 Adapted from the snoopy code
 **/
/*********************************************/
double randm(void) {
	const int a	=	16807;
	const int m =	2147483647;
	static int in0 = 13763;
	int q;
	
	// When using mpi, this allows us to have different number series in each process...
    
#ifdef PARALLEL
	if(in0 == 13763) in0 += 2543 * prank;
#endif
    
	/* find random number  */
	q= (int) fmod((double) a * in0, m);
	in0=q;
	
	return((double)q/(double)m);
}


/* ************************************************************** */
void Init (double *v, double x1, double x2, double x3)
/* 
 *
 * 
 *
 *
 **************************************************************** */
{
  // Define some sort of initial condition   


   g_isoSoundSpeed  = g_inputParam[CS];
  
    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = g_inputParam[BZ0];
   v[VX1] = 0.2*(randm()-0.5);
   v[VX2] = 2.0*sb_A * x1;
   v[VX3] = 0.0;
    
   v[RHO] = exp( - x3*x3/(2.0*g_inputParam[CS]*g_inputParam[CS]));
    v[RHO]=1.0;
    
    
}

void WriteProfile(real ***d, Grid *grid) {
	// This routine prints a vertical profile in a text file
	
	static real *vertProfLoc;
	static real *vertProfGlob;
	static int first_call = 1;
	int i,j,k;
	int kglob;
	FILE *ht;
	
	if(first_call == 1) {
		
		first_call = 0;
		// Allocate 1 vertical Array
		vertProfLoc = (real *) malloc(sizeof(real) * grid[KDIR].np_int_glob);
		vertProfGlob = (real *) malloc(sizeof(real) * grid[KDIR].np_int_glob);
		
		// Should also init the profile file
		if(g_time==0.0) {
		if(prank==0) {
			ht = fopen("profile","w");
			fprintf(ht,"%08e\t",g_time);
			for(i = 0 ; i < grid[KDIR].np_int_glob ; i++ ) {
				fprintf(ht,"%08e\t",grid[KDIR].x[i+grid[KDIR].gbeg]);
			}
			fprintf(ht,"\n");
			fclose(ht);
		}
		}
	}
	
	// Set the array to 0
	
	for(i = 0 ; i < grid[KDIR].np_int_glob ; i++ ) {
		vertProfLoc[i] = 0.0;
	}

	// Compute the sum
	DOM_LOOP(k,j,i) {
		kglob = k - 2*grid[KDIR].nghost + grid[KDIR].beg;
		
		if(kglob >= grid[KDIR].np_int_glob) printf("rank %d overshoot\n",prank);
		if(kglob < 0) printf("rank %d undershoot\n",prank);
		vertProfLoc[kglob] += d[k][j][i];
		
	}
	
	// Reduce
	MPI_Reduce( vertProfLoc, vertProfGlob, grid[KDIR].np_int_glob, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
	
	// Divide by the number of points
	for(i = 0 ; i < grid[KDIR].np_int_glob ; i++ ) {
		vertProfGlob[i] = vertProfGlob[i] / ((double) grid[IDIR].np_int_glob * grid[JDIR].np_int_glob );
	}
	
	// Print to file
	if(prank==0) {
		ht = fopen("profile","a");
		fprintf(ht,"%08e\t",g_time);
		for(i = 0 ; i < grid[KDIR].np_int_glob ; i++ ) {
			fprintf(ht,"%08e\t",vertProfGlob[i]);
		}
		fprintf(ht,"\n");
		fclose(ht);
	}
	
	return;
}
	

//allocate a 3D array
double*** allocate3D(int l,int m,int n)
{
    double ***arr3D;
    double *array;
    int i,j,k;
    
    array = (double*) malloc( l*m*n*sizeof(double));
    
    arr3D = (double***)malloc(l * sizeof(double **));
    
    for(i=0;i<l;i++)
    {
        arr3D[i] = (double**)malloc(m * sizeof(double *));
        for(j=0;j<m;j++)
        {
            arr3D[i][j] = &array[m*n*i+n*j];
        }
    }
    
    return arr3D;
}

//deallocate a 3D array
void deallocate3D(double ***arr3D,int l,int m)
{
    int i,j;
    free(arr3D[0][0]);
    
    for(i=0;i<l;i++)
    {
    	free(arr3D[i]);
    }
    free(arr3D);
}



	
/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 **************************************************************** */
{
  int i,j,k;
  int ntot;
  FILE *ht;
  static int first_call = 1;
  
  real mass;
  real vxmax, vxmin, vymax, vymin, vzmax, vzmin;
  real bxmax, bxmin, bymax, bymin, bzmax, bzmin;
  real ek,em;
  real dmax, dmin;
  real bxby, dvxvy;
  real vxm, vym, vzm;
  real bxm, bym, bzm;
  real q0;
  real mflux;
  real txx,txy,txz,tyz,tyy,tzz, divv;
  real dxvx,dxvy,dxvz,dyvx,dyvy,dyvz,dzvx,dzvy,dzvz;
  real dx,dy,dz;
    
  static real ***temp1;
  static real ***temp2;
  static real ***temp3;
  static real ***temp4;
  static real ***vy;
  
    if(first_call == 1) {
        
        first_call = 0;
          // Reinitialize the timevar output file
        if(g_time==0.0) {
            if(prank==0) {
                ht = fopen("timevar","w");
                fprintf(ht,"t\t\t");
                fprintf(ht,"mass\t\t");
                fprintf(ht,"dmax\t\t");
                fprintf(ht,"dmin\t\t");
                fprintf(ht,"ek\t\t");
                fprintf(ht,"em\t\t");
                fprintf(ht,"bxby\t\t");
                fprintf(ht,"dvxvy\t\t");
                fprintf(ht,"vxmax\t\t");
                fprintf(ht,"vxmin\t\t");
                fprintf(ht,"vymax\t\t");
                fprintf(ht,"vymin\t\t");
                fprintf(ht,"vzmax\t\t");
                fprintf(ht,"vzmin\t\t");
                fprintf(ht,"bxmax\t\t");
                fprintf(ht,"bxmin\t\t");
                fprintf(ht,"bymax\t\t");
                fprintf(ht,"bymin\t\t");
                fprintf(ht,"bzmax\t\t");
                fprintf(ht,"bzmin\t\t");
                fprintf(ht,"vxm\t\t");
                fprintf(ht,"vym\t\t");
                fprintf(ht,"vzm\t\t");
                fprintf(ht,"bxm\t\t");
                fprintf(ht,"bym\t\t");
                fprintf(ht,"bzm\t\t");
                fprintf(ht,"vzd\t\t");
                fprintf(ht,"\n");
                fclose(ht);
            }}
        // Initialize the temporary array
        temp1 = allocate3D(NX3_TOT, NX2_TOT, NX1_TOT);
        temp2 = allocate3D(NX3_TOT, NX2_TOT, NX1_TOT);
        temp3 = allocate3D(NX3_TOT, NX2_TOT, NX1_TOT);
        temp4 = allocate3D(NX3_TOT, NX2_TOT, NX1_TOT);
        vy = allocate3D(NX3_TOT, NX2_TOT, NX1_TOT);
    
  }

  // Compute vy by removing the mean shear
  TOT_LOOP(k,j,i) {
      vy[k][j][i] = d->Vc[VX2][k][j][i] - 2.0* sb_A * grid[IDIR].x[i];
  }
    
    // Write some profiles  

  
  WriteProfile(d->Vc[RHO], grid);
  WriteProfile(d->Vc[VX1], grid);
  WriteProfile(vy        , grid);
  WriteProfile(d->Vc[VX3], grid);
    
  WriteProfile(d->Vc[BX1], grid);
  WriteProfile(d->Vc[BX2], grid);
  WriteProfile(d->Vc[BX3], grid);
    
    // Reynolds stress
  DOM_LOOP(k,j,i) {
      temp1[k][j][i] = d->Vc[RHO][k][j][i] * d->Vc[VX1][k][j][i] * vy[k][j][i];
  }
    
  WriteProfile(temp1, grid);
    
    // Maxwell stress
  DOM_LOOP(k,j,i) {
        temp1[k][j][i] = d->Vc[BX1][k][j][i] * d->Vc[BX2][k][j][i];
    }
  
    WriteProfile(temp1, grid);
    
    // total Magnetic pressure
    
    DOM_LOOP(k,j,i) {
        temp1[k][j][i] = 0.5*(d->Vc[BX1][k][j][i] * d->Vc[BX1][k][j][i] +
        					  d->Vc[BX2][k][j][i] * d->Vc[BX2][k][j][i] +
        					  d->Vc[BX3][k][j][i] * d->Vc[BX3][k][j][i] );
    }
  	WriteProfile(temp1, grid);
  	
  	// Current density
  	DOM_LOOP(k,j,i) {
  		// Jx=dyBz-dzBy
        temp1[k][j][i] = (d->Vc[BX3][k][j+1][i] - d->Vc[BX3][k][j-1][i]) / (grid[JDIR].x[j+1]-grid[JDIR].x[j-1]) -
        				 (d->Vc[BX2][k+1][j][i] - d->Vc[BX2][k-1][j][i]) / (grid[KDIR].x[k+1]-grid[KDIR].x[k-1]);
        // Jy=dzBx-dxBz				 
        temp2[k][j][i] = (d->Vc[BX1][k+1][j][i] - d->Vc[BX1][k-1][j][i]) / (grid[KDIR].x[k+1]-grid[KDIR].x[k-1]) -
        				 (d->Vc[BX3][k][j][i+1] - d->Vc[BX3][k][j][i-1]) / (grid[IDIR].x[i+1]-grid[IDIR].x[i-1]);
        
        // Jz=dxBy-dyBx				 
        temp3[k][j][i] = (d->Vc[BX2][k][j][i+1] - d->Vc[BX2][k][j][i-1]) / (grid[IDIR].x[i+1]-grid[IDIR].x[i-1]) -
        				 (d->Vc[BX1][k][j+1][i] - d->Vc[BX1][k][j-1][i]) / (grid[JDIR].x[j+1]-grid[JDIR].x[j-1]);
        				 
        //J2
        temp4[k][j][i] = temp1[k][j][i] * temp1[k][j][i] +
        				 temp2[k][j][i] * temp2[k][j][i] +
        				 temp3[k][j][i] * temp3[k][j][i];
  	}
  	
  	WriteProfile(temp4,grid);
    
    // Viscous stress dissipation and pressure work
    DOM_LOOP(k,j,i) {
        dx=(grid[IDIR].x[i+1]-grid[IDIR].x[i-1]);
        dy=(grid[JDIR].x[j+1]-grid[JDIR].x[j-1]);
        dz=(grid[KDIR].x[k+1]-grid[KDIR].x[k-1]);
        
        dxvx = (d->Vc[VX1][k][j][i+1] - d->Vc[VX1][k][j][i-1]) / dx;
        dxvy = (        vy[k][j][i+1] -         vy[k][j][i-1]) / dx;
        dxvz = (d->Vc[VX3][k][j][i+1] - d->Vc[VX3][k][j][i-1]) / dx;
        
        dyvx = (d->Vc[VX1][k][j+1][i] - d->Vc[VX1][k][j-1][i]) / dy;
        dyvy = (        vy[k][j+1][i] -         vy[k][j-1][i]) / dy;
        dyvz = (d->Vc[VX3][k][j+1][i] - d->Vc[VX3][k][j-1][i]) / dy;
        
        dzvx = (d->Vc[VX1][k+1][j][i] - d->Vc[VX1][k-1][j][i]) / dz;
        dzvy = (        vy[k+1][j][i] -         vy[k-1][j][i]) / dz;
        dzvz = (d->Vc[VX3][k+1][j][i] - d->Vc[VX3][k-1][j][i]) / dz;
        
        divv = dxvx+dyvy+dzvz;
    
        txx = 2.0*dxvx-2.0/3.0*divv;
        txy = dxvy+dyvx;
        txz = dxvz+dzvz;
        tyy = 2.0*dyvy-2.0/3.0*divv;
        tyz = dyvz+dzvy;
        tzz = 2.0*dzvz-2.0/3.0*divv;
        
        temp4[k][j][i] = 0.5*d->Vc[RHO][k][j][i]*(txx*txx + tyy*tyy + tzz*tzz + 2.0*txy*txy + 2.0*txz*txz + 2.0*tyz*tyz);
        
        temp1[k][j][i] = g_isoSoundSpeed*g_isoSoundSpeed*d->Vc[RHO][k][j][i] * divv;
        
    }
    WriteProfile(temp4,grid);
    WriteProfile(temp1,grid);
    
    // velocity dispersions (useful for dust settling)
    DOM_LOOP(k,j,i) {
        temp1[k][j][i] = d->Vc[VX1][k][j][i] * d->Vc[VX1][k][j][i];
        temp2[k][j][i] =         vy[k][j][i] *         vy[k][j][i];
        temp3[k][j][i] = d->Vc[VX3][k][j][i] * d->Vc[VX3][k][j][i];
    }
    WriteProfile(temp1,grid);
    WriteProfile(temp2,grid);
    WriteProfile(temp3,grid);
    
    // Work done by pressure forces
      
    
        
  // Timevar file
  
  vxmax = d->Vc[VX][KBEG][JBEG][IBEG];
  vymax = d->Vc[VY][KBEG][JBEG][IBEG];
  vxmin = d->Vc[VX][KBEG][JBEG][IBEG];
  vymin = d->Vc[VY][KBEG][JBEG][IBEG];
  
  bxmax = d->Vc[BX][KBEG][JBEG][IBEG];
  bymax = d->Vc[BY][KBEG][JBEG][IBEG];
  bxmin = d->Vc[BX][KBEG][JBEG][IBEG];
  bymin = d->Vc[BY][KBEG][JBEG][IBEG];
  
  dmax = d->Vc[DN][KBEG][JBEG][IBEG];
  dmin = d->Vc[DN][KBEG][JBEG][IBEG];

  mass = 0.0;
		
  bxby = 0.0;
  dvxvy = 0.0;
  
  ek=0.0;
  em=0.0;
    
  vxm = 0.0;
  vym = 0.0;
  vzm = 0.0;
  
  bxm = 0.0;
  bym = 0.0;
  bzm = 0.0;
  
  #if DIMENSIONS == 3
  vzmax = d->Vc[VZ][KBEG][JBEG][IBEG];
  vzmin = d->Vc[VZ][KBEG][JBEG][IBEG];
  bzmax = d->Vc[BZ][KBEG][JBEG][IBEG];
  bzmin = d->Vc[BZ][KBEG][JBEG][IBEG];
  
  #elif DIMENSIONS == 2
  vzmax = 0.0;
  vzmin = 0.0;
  bzmax = 0.0;
  bzmin = 0.0;
  #endif
  
  
  DOM_LOOP(k,j,i) {
	mass  += d->Vc[DN][k][j][i];
	bxby  += d->Vc[BX][k][j][i] * d->Vc[BY][k][j][i];
	dvxvy += d->Vc[DN][k][j][i] * d->Vc[VX][k][j][i] * vy[k][j][i];
      
    ek += 0.5*d->Vc[DN][k][j][i]*( d->Vc[VX][k][j][i] * d->Vc[VX][k][j][i] +
                                   vy[k][j][i]        * vy[k][j][i] +
                                  d->Vc[VZ][k][j][i] * d->Vc[VZ][k][j][i] );
    em += 0.5*( d->Vc[BX][k][j][i] * d->Vc[BX][k][j][i] +
                d->Vc[BY][k][j][i] * d->Vc[BY][k][j][i] +
                d->Vc[BZ][k][j][i] * d->Vc[BZ][k][j][i] );
	
	vxm += d->Vc[DN][k][j][i] * d->Vc[VX][k][j][i];
	vym += d->Vc[DN][k][j][i] * vy[k][j][i];
	
	bxm += d->Vc[BX][k][j][i];
	bym += d->Vc[BY][k][j][i];
	
	if( d->Vc[DN][k][j][i] > dmax) dmax = d->Vc[DN][k][j][i];
	if( d->Vc[DN][k][j][i] < dmin) dmin = d->Vc[DN][k][j][i];
	
	if( d->Vc[VX][k][j][i] > vxmax) vxmax = d->Vc[VX][k][j][i];
	if( d->Vc[VX][k][j][i] < vxmin) vxmin = d->Vc[VX][k][j][i];
	if( vy[k][j][i] > vymax) vymax = vy[k][j][i];
	if( vy[k][j][i] < vymin) vymin = vy[k][j][i];
	
	if( d->Vc[BX][k][j][i] > bxmax) bxmax = d->Vc[BX][k][j][i];
	if( d->Vc[BX][k][j][i] < bxmin) bxmin = d->Vc[BX][k][j][i];
	if( d->Vc[BY][k][j][i] > bymax) bymax = d->Vc[BY][k][j][i];
	if( d->Vc[BY][k][j][i] < bymin) bymin = d->Vc[BY][k][j][i];
	
	
	#if DIMENSIONS == 3
	vzm += d->Vc[DN][k][j][i] * d->Vc[VZ][k][j][i];
	bzm += d->Vc[BZ][k][j][i];
	
	if( d->Vc[VZ][k][j][i] > vzmax) vzmax = d->Vc[VZ][k][j][i];
	if( d->Vc[VZ][k][j][i] < vzmin) vzmin = d->Vc[VZ][k][j][i];
	
	if( d->Vc[BZ][k][j][i] > bzmax) bzmax = d->Vc[BZ][k][j][i];
	if( d->Vc[BZ][k][j][i] < bzmin) bzmin = d->Vc[BZ][k][j][i];
	
	#endif
	
	
  }
  
  /* Flux */
  mflux = 0.0;
  
  for( j = JBEG ; j <= JEND ; j++) {
	for( i = IBEG ; i <= IEND ; i++) {
		if(grid[KDIR].gend == grid[KDIR].end) {
			mflux += d->Vc[DN][KEND][j][i]*d->Vc[VZ][KEND][j][i];
		}
		if(grid[KDIR].gbeg == grid[KDIR].beg) {
			mflux -= d->Vc[DN][KBEG][j][i]*d->Vc[VZ][KBEG][j][i];
		}
	}
  }
		
		
  
  /* Normalization */
  
  #if DIMENSIONS == 2
  ntot = grid[IDIR].np_int_glob * grid[JDIR].np_int_glob ;
  #elif DIMENSIONS == 3
  ntot = grid[IDIR].np_int_glob * grid[JDIR].np_int_glob * grid[KDIR].np_int_glob ;
  #endif
  
  mass = mass / ntot;
  bxby = bxby / ntot;
  dvxvy = dvxvy / ntot;
    
    ek = ek / ntot;
    em = em / ntot;
  
  vxm = vxm / ntot;
  vym = vym / ntot;
  vzm = vzm / ntot;
  
  bxm = bxm / ntot;
  bym = bym / ntot;
  bzm = bzm / ntot;
  
  /* Flux normalization */
  #if DIMENSIONS == 2
  ntot = grid[IDIR].np_int_glob * grid[JDIR].np_int_glob ;
  #elif DIMENSIONS == 3
  ntot = grid[IDIR].np_int_glob * grid[JDIR].np_int_glob ;
  #endif
  
  mflux = mflux / (2.0*ntot);
  
  /* reduction */
  
  #ifdef PARALLEL
  
  MPI_Allreduce( &mass, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mass = q0;
  MPI_Allreduce( &bxby, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  bxby = q0;
  MPI_Allreduce( &dvxvy, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  dvxvy = q0;
  
  MPI_Allreduce( &ek, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ek = q0;
  MPI_Allreduce( &em, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  em = q0;

    
  MPI_Allreduce( &vxm, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  vxm = q0;
  MPI_Allreduce( &vym, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  vym = q0;
  MPI_Allreduce( &vzm, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  vzm = q0;
  
  MPI_Allreduce( &bxm, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  bxm = q0;
  MPI_Allreduce( &bym, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  bym = q0;
  MPI_Allreduce( &bzm, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  bzm = q0;
  
   MPI_Allreduce( &mflux, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mflux = q0;
  
  MPI_Allreduce( &dmax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  dmax = q0;
  MPI_Allreduce( &dmin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  dmin = q0;
  
  MPI_Allreduce( &vxmax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  vxmax = q0;
  MPI_Allreduce( &vxmin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  vxmin = q0;
  MPI_Allreduce( &vymax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  vymax = q0;
  MPI_Allreduce( &vymin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  vymin = q0;
  MPI_Allreduce( &vzmax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  vzmax = q0;
  MPI_Allreduce( &vzmin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  vzmin = q0;
  
  MPI_Allreduce( &bxmax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  bxmax = q0;
  MPI_Allreduce( &bxmin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  bxmin = q0;
  MPI_Allreduce( &bymax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  bymax = q0;
  MPI_Allreduce( &bymin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  bymin = q0;
  MPI_Allreduce( &bzmax, &q0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  bzmax = q0;
  MPI_Allreduce( &bzmin, &q0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  bzmin = q0;
  
  #endif	// PARALLEL
  
  
  if(prank==0) {
	ht = fopen("timevar","a");
	fprintf(ht,"%08e\t",g_time);
	fprintf(ht,"%08e\t",mass);
	fprintf(ht,"%08e\t",dmax);
	fprintf(ht,"%08e\t",dmin);
    fprintf(ht,"%08e\t",ek);
    fprintf(ht,"%08e\t",em);
	fprintf(ht,"%08e\t",bxby);
	fprintf(ht,"%08e\t",dvxvy);
	fprintf(ht,"%08e\t",vxmax);
	fprintf(ht,"%08e\t",vxmin);
	fprintf(ht,"%08e\t",vymax);
	fprintf(ht,"%08e\t",vymin);
	fprintf(ht,"%08e\t",vzmax);
	fprintf(ht,"%08e\t",vzmin);
	fprintf(ht,"%08e\t",bxmax);
	fprintf(ht,"%08e\t",bxmin);
	fprintf(ht,"%08e\t",bymax);
	fprintf(ht,"%08e\t",bymin);
	fprintf(ht,"%08e\t",bzmax);
	fprintf(ht,"%08e\t",bzmin);
	fprintf(ht,"%08e\t",vxm);
	fprintf(ht,"%08e\t",vym);
	fprintf(ht,"%08e\t",vzm);
	fprintf(ht,"%08e\t",bxm);
	fprintf(ht,"%08e\t",bym);
	fprintf(ht,"%08e\t",bzm);
	fprintf(ht,"%08e\t",mflux);

	fprintf(ht,"\n");
	fclose(ht);
   }
  
  
  return;
  

}


/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 * PURPOSE
 *
 *   Define user-defined boundary conditions.
 *
 *
 * ARGUMENTS
 * 
 *   d      the PLUTO data structure containing both cell
 *          centered primitive quantities (d->Vc) and 
 *          staggered (d->Vs) magnetic fields (when used) to 
 *          be filled.
 *
 *   side   specifies on which side boundary conditions need 
 *          to be assigned. side can assume the following 
 *          pre-definite values: X1_BEG, X1_END,
 *                               X2_BEG, X2_END, 
 *                               X3_BEG, X3_END
 *
 *                     When 
 *
 *                      side = Xn_BEG  -->  b.c. are assigned at the 
 *                                          beginning of the xn
 *                                          direction.
 *                      side = Xn_BEG  -->  b.c. are assigned at the 
 *                                          beginning of the xn
 *                                          direction.
 *
 *          The special value side == 0 is used to control
 *          a region inside the computational domain.
 *
 **************************************************************** */
{
  int   i, j, k, nv;
  real  *x1, *x2, *x3;
  if (side == 0) {    /* -- check solution inside domain -- */
  }


  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;


  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
  	if (box->vpos == CENTER) {
		BOX_LOOP(box,k,j,i) {
			d->Vc[DN][k][j][i] = d->Vc[DN][KBEG][j][i] * exp( - (x3[k]*x3[k]-x3[KBEG]*x3[KBEG])/(2.0*g_inputParam[CS]*g_inputParam[CS]));
			
			d->Vc[VX1][k][j][i] = d->Vc[VX1][k+1][j][i];
			d->Vc[VX2][k][j][i] = d->Vc[VX2][k+1][j][i];
			d->Vc[VX3][k][j][i] = -d->Vc[VX3][2*KBEG-k][j][i];

		}
  	}
	else if (box->vpos == X1FACE) {
		#ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) d->Vs[BX1s][k][j][i] = 0.0;
		#endif
	}
	else if (box->vpos == X2FACE) {
		#ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = 0.0;
		#endif
	}
  }

  
if (side == X3_END){  /* -- X3_END boundary -- */
  	if (box->vpos == CENTER) {
		BOX_LOOP(box,k,j,i) {
			d->Vc[DN][k][j][i] = d->Vc[DN][KEND][j][i] * exp( - (x3[k]*x3[k]-x3[KEND]*x3[KEND])/(2.0*g_inputParam[CS]*g_inputParam[CS]));
			
			d->Vc[VX1][k][j][i] = d->Vc[VX1][k-1][j][i];
			d->Vc[VX2][k][j][i] = d->Vc[VX2][k-1][j][i];
			d->Vc[VX3][k][j][i] = -d->Vc[VX3][2*KEND-k][j][i];

		
		}
  	}
	else if (box->vpos == X1FACE) {
		#ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) d->Vs[BX1s][k][j][i] = 0.0;
		#endif
	}
	else if (box->vpos == X2FACE) {
		#ifdef STAGGERED_MHD
			BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = 0.0;
		#endif
	}

  }


}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
#ifdef FARGO
    g[IDIR] = 0.0;
    g[JDIR] = -2.0*sb_A*v[VX1];
#else
    g[IDIR] = sb_Omega*sb_A*x1;
    g[JDIR] = 0.0;
#endif
    
    /* -- vertical gravity ? -- */
    
    
    g[KDIR] = 0.0;
    
    /*
    g[KDIR] = -sb_Omega*sb_Omega*x3;
    */ 
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
    return 0.0;
}
#endif

