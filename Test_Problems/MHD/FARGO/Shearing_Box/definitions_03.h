#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     3

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  CS                 0
#define  BETA               1
#define  EPS                2

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               MC_LIM
#define  CT_EMF_AVERAGE        UCT_CONTACT
#define  CT_EN_CORRECTION      NO
#define  ASSIGN_VECTOR_POTENTIAL  NO
#define  UPDATE_VECTOR_POTENTIAL  NO
