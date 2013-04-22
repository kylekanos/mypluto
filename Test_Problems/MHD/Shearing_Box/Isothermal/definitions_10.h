#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  INCLUDE_BODY_FORCE      EXPLICIT
#define  INCLUDE_COOLING         NO
#define  INTERPOLATION           WENO5Z_FD
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     2

/* -- physics dependent declarations -- */

#define  EOS   ISOTHERMAL
#define  ENTROPY_SWITCH   NO
#define  MHD_FORMULATION   DIV_CLEANING
#define  INCLUDE_BACKGROUND_FIELD   NO
#define  RESISTIVE_MHD   NO
#define  THERMAL_CONDUCTION   NO
#define  VISCOSITY   NO

/* -- pointers to user-def parameters -- */

#define  CS                 0
#define  BETA               1

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  SHOCK_FLATTENING      NO
#define  USE_VECTOR_POTENTIAL  YES
#define  SAVE_VEC_POT          NO
