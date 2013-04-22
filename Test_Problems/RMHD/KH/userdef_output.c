#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
  double ***T;

  T = GetUserVar("T");
  DOM_LOOP(k,j,i){
    T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
  }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

  SetDumpVar("pr",FLT_H5_OUTPUT, NO);

}





