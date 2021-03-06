/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the HD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 4, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int ibeg, int iend)
/*!
 * Convert primitive variables to conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int  i, nv;
  double *v, *u;
  double rho, gmm1;

  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif
  for (i = ibeg; i <= iend; i++) {
  
    v = uprim[i];
    u = ucons[i];

    u[RHO] = rho = v[RHO];
    EXPAND (u[MX1] = rho*v[VX1];  ,
            u[MX2] = rho*v[VX2];  ,
            u[MX3] = rho*v[VX3];)

    #if EOS == IDEAL
     u[ENG] = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
     u[ENG] = 0.5*rho*u[ENG] + v[PRS]/gmm1;
    #endif

    for (nv = NFLX; nv < (NFLX + NSCL); nv++) u[nv] = rho*v[nv];
  }
}
/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
                unsigned char *flag)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]  ucons  array of conservative variables
 * \param [out] uprim  array of primitive variables
 * \param [in]  beg    starting index of computation
 * \param [in]  end    final index of computation
 * \param [out] flag   array of flags tagging zones where conversion
 *                     went wrong.
 * 
 * \return Return (0) if conversion was succesful in every zone 
 *         [ibeg,iend]. 
 *         Otherwise, return a non-zero integer number giving the bit 
 *         flag(s) turned on during the conversion process.
 *         In this case, flag contains the failure codes of those
 *         zones where where conversion did not go through.
 *
 *********************************************************************** */
{
  int  i, nv, status=0, use_energy;
  double tau, rho, gmm1;
  double kin, m2, rhog1;
  double *u, *v;

  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif

  for (i = ibeg; i <= iend; i++) {

    flag[i] = 0;
    u = ucons[i];
    v = uprim[i];

    m2  = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    
  /* -------------------------------------------
           Check density positivity 
     ------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim: negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= RHO_FAIL;
    }

    v[RHO] = rho = u[RHO];
    tau = 1.0/u[RHO];
    EXPAND(v[VX1] = u[MX1]*tau;  ,
           v[VX2] = u[MX2]*tau;  ,
           v[VX3] = u[MX3]*tau;)

  /* --------------------------------------------
       now try to recover pressure from 
       total energy or entropy
     -------------------------------------------- */

    #if EOS == IDEAL
   
     kin = 0.5*m2/u[RHO];
     use_energy = 1;     
     #if ENTROPY_SWITCH == YES
      #ifdef CH_SPACEDIM
       if (g_intStage > 0)   /* -- HOT FIX used with Chombo: avoid calling 
                              CheckZone when writing file to disk          -- */
      #endif
       if (CheckZone(i,FLAG_ENTROPY)){
         use_energy = 0;
         rhog1 = pow(rho, g_gamma - 1.0);
         v[PRS] = u[ENTR]*rhog1; 
         if (v[PRS] < 0.0){
           WARNING(
             print("! ConsToPrim: negative p(S) (%8.2e, %8.2e), ", v[PRS], u[ENTR]);
             Where (i, NULL);
           )
           v[PRS]   = g_smallPressure;
           flag[i] |= PRS_FAIL;
         }
         u[ENG] = v[PRS]/gmm1 + kin; /* -- redefine energy -- */
       }
     #endif  /* ENTROPY_SWITCH == YES */

     if (use_energy){
       if (u[ENG] < 0.0) {
         WARNING(
           print("! ConsToPrim: negative energy (%8.2e), ", u[ENG]);
           Where (i, NULL);
         )
         u[ENG]   = g_smallPressure/gmm1 + kin;
         flag[i] |= ENG_FAIL;
       }else{
         v[PRS] = gmm1*(u[ENG] - kin);
         if (v[PRS] < 0.0){
           WARNING(
             print("! ConsToPrim: negative p(E) (%8.2e), ", v[PRS]);
             Where (i, NULL);
           )
           v[PRS]   = g_smallPressure;
           flag[i] |= PRS_FAIL;
           u[ENG]   = v[PRS]/gmm1 + kin; /* -- redefine energy -- */
         }
       }
     }

    #endif  /* EOS == IDEAL */

  /* -- do remaning variables -- */
   
    for (nv = NFLX; nv < NVAR; nv++) v[nv] = u[nv]*tau;

    status |= flag[i];
  }

  return(status);
}
