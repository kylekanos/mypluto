/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Hall MHD module header file.

  Contains prototypes for the resistive MHD module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
           G. Lesur
  \date   Apr 2013
*/
/* ///////////////////////////////////////////////////////////////////// */

//void HallFlux (Data_Arr, double **, double **, int, int, Grid *);
void ComputeJState(const Data *d, const Grid *grid, State_1D *state, int *in, int *i, int *j, int *k, int beg, int end);
void ComputeCurrent(Data_Arr, Grid *);
//void GetFullCurrent (Data_Arr, real **, Grid *);

void lHall_Func (real *, real, real, real, real *);
