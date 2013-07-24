/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Resistive MHD module header file.

  Contains prototypes for the resistive MHD module.

  \authors G. Lesur
  \date   Jul 24, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */

void AmbipolarFlux (Data_Arr, double **, double **, int, int, Grid *);

void GetCurrent (Data_Arr, real **, Grid *);

void AmbiETA_Func (real *, real, real, real, real *);
