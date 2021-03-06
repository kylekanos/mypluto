#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     5

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           EXPLICIT
#define  HALL_MHD                EXPLICIT
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               EXPLICIT
#define  ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  CS                 0
#define  BZ0                1
#define  NU                 2
#define  ETA                3
#define  LH                 4

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         NO
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      MULTID
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANLEER_LIM
#define  CT_EMF_AVERAGE        UCT_CONTACT
#define  CT_EN_CORRECTION      NO
#define  ASSIGN_VECTOR_POTENTIAL  NO
#define  UPDATE_VECTOR_POTENTIAL  NO
