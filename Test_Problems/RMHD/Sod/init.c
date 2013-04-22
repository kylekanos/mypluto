#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* ---------------------------------------------------
    Relativistic MHD shock tube problems.
    The tests considered are taken from
	
      Mignone & Bodo MNRAS(2006) 368,1040 [MB06]  
      Mignone et al. MNRAS(2009) 393,1141 [MUB09]
	
    Similar tests can be found in Balsara ApJS (2001) & 
    DelZanna et al., A&A 2003

    The different configurations refer to:

    Test    shock-tube
    ----    ----------

    01-02   1st in MB06
    03-04   2nd in MB06
    05-06   3rd in MB06
    07-08   4th in MB06
    09-10   2nd in MUB09
    11-12   4th in MUB09
    13-14   From Mignone et al.(2006) in Space Science Reviews
    	 
   --------------------------------------------------- */

  g_gamma = g_inputParam[GAMMA_EOS];

  if (x1 < 0.5){

    us[RHO] = g_inputParam[RHO_LEFT];
    us[VX1] = g_inputParam[VX_LEFT];
    us[VX2] = g_inputParam[VY_LEFT];
    us[VX3] = g_inputParam[VZ_LEFT];

    us[BX1] = g_inputParam[BX_CONST];
    us[BX2] = g_inputParam[BY_LEFT];
    us[BX3] = g_inputParam[BZ_LEFT];
    us[PRS] = g_inputParam[PR_LEFT];

  }else{

    us[RHO] = g_inputParam[RHO_RIGHT];
    us[VX1] = g_inputParam[VX_RIGHT];
    us[VX2] = g_inputParam[VY_RIGHT];
    us[VX3] = g_inputParam[VZ_RIGHT];
    
    us[BX1] = g_inputParam[BX_CONST];
    us[BX2] = g_inputParam[BY_RIGHT];
    us[BX3] = g_inputParam[BZ_RIGHT];
    
    us[PRS] = g_inputParam[PR_RIGHT];
 
  }


}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
}
