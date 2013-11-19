/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Current MHD module header file.

  Contains prototypes for the module computing currents for nonideal MHD terms

  \authors G. Lesur
  \date   Nov 2013
*/
/* ///////////////////////////////////////////////////////////////////// */

void ComputeJ(const Data *d, Grid *grid, double t);
void StoreJState(State_1D *state, int *in, int *i, int *j, int *k, int beg, int end);
void Curr_FixJ (double ***Jin, double t, Grid *grid);